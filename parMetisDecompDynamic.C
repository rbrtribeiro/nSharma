/*
 * parMetisDecompDynamic.cpp
 *
 *  Created on: Oct 3, 2014
 *      Author: rr
 */

#include "addToRunTimeSelectionTable.H"

extern "C" {
#   include "parmetis.h"
}

#include "parMetisDecompDynamic.H"
#include "parMetisDecomp.H"
#include "nSharma.H"
//#include "metisDecomp.H"
//#include "addToRunTimeSelectionTable.H"
//#include "floatScalar.H"
//#include "Time.H"
//#include "labelIOField.H"
#include "syncTools.H"
//#include "globalIndex.H"
#include "dynamicRefineExtFvMesh.H"

#include "fvMeshSubset.H"
#include "cellSet.H"
#include <vector>
#include <stdio.h>
using namespace Foam;

namespace Foam {
defineTypeNameAndDebug(parMetisDecompDynamic, 0);

addToRunTimeSelectionTable(decompositionMethod, parMetisDecompDynamic, dictionary);
}

Foam::parMetisDecompDynamic::parMetisDecompDynamic(
	const dictionary& decompositionDict) :
	Foam::parMetisDecomp(decompositionDict) {

}

Foam::parMetisDecompDynamic::parMetisDecompDynamic(
	const dictionary& decompositionDict, polyMesh& mesh, CommGraph* commGraph_,
	bool coarseMapping_, bool useCommGraph_) :
		Foam::parMetisDecomp(decompositionDict) {

	mesh_ = &mesh;
	commGraph = commGraph_;
	first = true;
	coarseMapping= coarseMapping_;
	useCommGraph= useCommGraph_;

}

/*
 * AMR proxy part I
 */
Foam::labelList Foam::parMetisDecompDynamic::decompose
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    CompactListList<label> coarseCellCells;
    calcCellCells
    (
        mesh,
        fineToCoarse,
        coarsePoints.size(),
        true,                       // use global cell labels
        coarseCellCells
    );

    // Decompose based on agglomerated points
    labelList coarseDistribution
    (
        decompose
        (
            coarseCellCells(),
            coarsePoints,
            coarseWeights,
			fineToCoarse
        )
    );

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}

/*
 * AMR proxy part II
 */
Foam::labelList Foam::parMetisDecompDynamic::decompose(
	const labelListList& globalCellCells, const pointField& cellCentres,
		const scalarField& cWeights, const labelList& fineToCoarse) {
	if (cellCentres.size() != globalCellCells.size()) {
		FatalErrorIn("parMetisDecomp::decompose(const labelListList&"
				", const pointField&, const scalarField&)")
				<< "Inconsistent number of cells (" << globalCellCells.size()
				<< ") and number of cell centres (" << cellCentres.size()
				<< ") or weights (" << cWeights.size() << ")." << exit(FatalError);
	}

	// Make Metis Distributed CSR (Compressed Storage Format) storage

	// Connections
	Field<int> adjncy;

	// Offsets into adjncy
	Field<int> xadj;

	calcCSR(globalCellCells, adjncy, xadj);

	// decomposition options. 0 = use defaults
	List<int> options(3, 0);

	// cell weights (so on the vertices of the dual)
	Field<int> cellWeights;

	// face weights (so on the edges of the dual)
	Field<int> faceWeights;

	// Check for externally provided cellweights and if so initialise weights
	scalar minWeights = gMin(cWeights);
	if (cWeights.size() > 0) {
		if (minWeights <= 0) {
			WarningIn("parMetisDecomp::decompose(const labelListList&"
					", const pointField&, const scalarField&)")
					<< "Illegal minimum weight " << minWeights << endl;
		}

		if (cWeights.size() != globalCellCells.size()) {
			FatalErrorIn("parMetisDecomp::decompose(const labelListList&"
					", const pointField&, const scalarField&)")
					<< "Number of cell weights " << cWeights.size()
					<< " does not equal number of cells " << globalCellCells.size()
					<< exit(FatalError);
		}

		// Convert to integers.
		cellWeights.setSize(cWeights.size());
		forAll(cellWeights, i)
		{
			cellWeights[i] = int(cWeights[i] / minWeights);
		}
	}

	// Do actual decomposition
	List<int> finalDecomp;

	if (useCommGraph){
		decomposeBasedCommGraph(xadj, adjncy, cellCentres, cellWeights, faceWeights, options,
			finalDecomp, fineToCoarse);
	} else {
		decompose(xadj, adjncy, cellCentres, cellWeights, faceWeights, options,
					finalDecomp);
	}

	// Copy back to labelList
	labelList decomp(finalDecomp.size());
	forAll(decomp, i)
	{
		decomp[i] = finalDecomp[i];
	}
	return decomp;
}

/*
 * static proxy
 */
Foam::labelList Foam::parMetisDecompDynamic::decompose(const polyMesh& mesh_,
		const pointField& cc, const scalarField& cWeights) {

	if (cc.size() != mesh_.nCells()) {
		FatalErrorIn("parMetisDecomp::decompose"
				"(const pointField&, const scalarField&)")
				<< "Can use this decomposition method only for the whole mesh"
				<< endl << "and supply one coordinate (cellCentre) for every cell."
				<< endl << "The number of coordinates " << cc.size() << endl
				<< "The number of cells in the mesh " << mesh_.nCells()
				<< exit(FatalError);
	}

	// Connections
	Field<int> adjncy;
	// Offsets into adjncy
	Field<int> xadj;

	Foam::calcDistributedCSR(mesh_, adjncy, xadj);

	// decomposition options. 0 = use defaults
	List<int> options(3, 0);

	// cell weights (so on the vertices of the dual)
	Field<int> cellWeights;

	// face weights (so on the edges of the dual)
	Field<int> faceWeights;

	// Check for externally provided cellweights and if so initialise weights
	scalar minWeights = gMin(cWeights);

	if (cWeights.size() > 0) {
		if (minWeights <= 0) {
			WarningIn("metisDecomp::decompose"
					"(const pointField&, const scalarField&)")
					<< "Illegal minimum weight " << minWeights << endl;
		}

		if (cWeights.size() != mesh_.nCells()) {
			FatalErrorIn("parMetisDecomp::decompose"
					"(const pointField&, const scalarField&)")
					<< "Number of cell weights " << cWeights.size()
					<< " does not equal number of cells " << mesh_.nCells()
					<< exit(FatalError);
		}

		// Convert to integers.
		cellWeights.setSize(cWeights.size());
		forAll(cellWeights, i)
		{
			cellWeights[i] = int(cWeights[i] / minWeights);
		}
	}

	// Do actual decomposition
	List<int> finalDecomp;

	Info << "useCommGraph " << useCommGraph << endl;

	if (useCommGraph){
		decomposeBasedCommGraph(xadj, adjncy, cc, cellWeights, faceWeights,
				options, finalDecomp,labelList() );

	} else {
		decompose(xadj, adjncy, cc, cellWeights, faceWeights, options,
					finalDecomp);
	}



	// Copy back to labelList
	labelList decomp(finalDecomp.size());
	forAll(decomp, i)
	{
		decomp[i] = finalDecomp[i];
	}
	return decomp;
}

//- Does prevention of 0 cell domains and calls parmetis.
Foam::label Foam::parMetisDecompDynamic::decompose(Field<int>& xadj,
		Field<int>& adjncy, const pointField& cellCentres,
		Field<int>& cellWeights, Field<int>& faceWeights,
		const List<int>& options, List<int>& finalDecomp) {

	Info << "Using parMetisDecompDynamic::decompose" << endl;

	// C style numbering
	int numFlag = 0;

	// Number of dimensions
	int nDims = 3;

	if (cellCentres.size() != xadj.size() - 1) {
		FatalErrorIn("parMetisDecomp::decompose(..)") << "cellCentres:"
				<< cellCentres.size() << " xadj:" << xadj.size()
				<< abort(FatalError);
	}

	// Get number of cells on all processors
	List<int> nLocalCells(Pstream::nProcs());
	nLocalCells[Pstream::myProcNo()] = xadj.size() - 1;
	Pstream::gatherList(nLocalCells);
	Pstream::scatterList(nLocalCells);

	// Get cell offsets.
	List<int> cellOffsets(Pstream::nProcs() + 1);
	int nGlobalCells = 0;
	forAll(nLocalCells, procI){
	cellOffsets[procI] = nGlobalCells;
	nGlobalCells += nLocalCells[procI];
}
	cellOffsets[Pstream::nProcs()] = nGlobalCells;

	// Convert pointField into the data type parMetis expects (float or double)
	Field<real_t> xyz(3 * cellCentres.size());
	int compI = 0;
	forAll(cellCentres, cellI){
	const point& cc = cellCentres[cellI];
	xyz[compI++] = float(cc.x());
	xyz[compI++] = float(cc.y());
	xyz[compI++] = float(cc.z());
}

// Make sure every domain has at least one cell
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// (Metis falls over with zero sized domains)
// Trickle cells from processors that have them up to those that
// don't.

// Number of cells to send to the next processor
// (is same as number of cells next processor has to receive)
	List<int> nSendCells(Pstream::nProcs(), 0);

	for (label procI = nLocalCells.size() - 1; procI >= 1; procI--) {
		if (nLocalCells[procI] - nSendCells[procI] < 1) {
			nSendCells[procI - 1] = nSendCells[procI] - nLocalCells[procI] + 1;
		}
	}

	// First receive (so increasing the sizes of all arrays)

	if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo() - 1] > 0) {
		// Receive cells from previous processor
		IPstream fromPrevProc(Pstream::blocking, Pstream::myProcNo() - 1);

		Field<int> prevXadj(fromPrevProc);
		Field<int> prevAdjncy(fromPrevProc);
		Field<real_t> prevXyz(fromPrevProc);
		Field<int> prevCellWeights(fromPrevProc);
		Field<int> prevFaceWeights(fromPrevProc);

		if (prevXadj.size() != nSendCells[Pstream::myProcNo() - 1]) {
			FatalErrorIn("parMetisDecomp::decompose(..)")
					<< "Expected from processor " << Pstream::myProcNo() - 1
					<< " connectivity for "
					<< nSendCells[Pstream::myProcNo() - 1]
					<< " nCells but only received " << prevXadj.size()
					<< abort(FatalError);
		}

		// Insert adjncy
		prepend(prevAdjncy, adjncy);
		// Adapt offsets and prepend xadj
		xadj += prevAdjncy.size();
		prepend(prevXadj, xadj);
		// Coords
		prepend(prevXyz, xyz);
		// Weights
		prepend(prevCellWeights, cellWeights);
		prepend(prevFaceWeights, faceWeights);
	}

	// Send to my next processor

	if (nSendCells[Pstream::myProcNo()] > 0) {
		// Send cells to next processor
		OPstream toNextProc(Pstream::blocking, Pstream::myProcNo() + 1);

		int nCells = nSendCells[Pstream::myProcNo()];
		int startCell = xadj.size() - 1 - nCells;
		int startFace = xadj[startCell];
		int nFaces = adjncy.size() - startFace;

		// Send for all cell data: last nCells elements
		// Send for all face data: last nFaces elements
		toNextProc << Field<int>::subField(xadj, nCells, startCell) - startFace
				<< Field<int>::subField(adjncy, nFaces, startFace)
				<< SubField<real_t>(xyz, nDims * nCells, nDims * startCell)
				<< (cellWeights.size() ?
						static_cast<const Field<int>&>(Field<int>::subField(
								cellWeights, nCells, startCell)) :
						Field<int>(0))
				<< (faceWeights.size() ?
						static_cast<const Field<int>&>(Field<int>::subField(
								faceWeights, nFaces, startFace)) :
						Field<int>(0));

		// Remove data that has been sent
		if (faceWeights.size()) {
			faceWeights.setSize(faceWeights.size() - nFaces);
		}
		if (cellWeights.size()) {
			cellWeights.setSize(cellWeights.size() - nCells);
		}
		xyz.setSize(xyz.size() - nDims * nCells);
		adjncy.setSize(adjncy.size() - nFaces);
		xadj.setSize(xadj.size() - nCells);
	}

	// Adapt number of cells
	forAll(nSendCells, procI){
	// Sent cells
	nLocalCells[procI] -= nSendCells[procI];

	if (procI >= 1)
	{
		// Received cells
		nLocalCells[procI] += nSendCells[procI-1];
	}
}
// Adapt cellOffsets
	nGlobalCells = 0;
	forAll(nLocalCells, procI){
	cellOffsets[procI] = nGlobalCells;
	nGlobalCells += nLocalCells[procI];
}

	if (nLocalCells[Pstream::myProcNo()] != (xadj.size() - 1)) {
		FatalErrorIn("parMetisDecomp::decompose(..)")
				<< "Have connectivity for " << xadj.size() - 1
				<< " cells but nLocalCells:" << nLocalCells[Pstream::myProcNo()]
				<< abort(FatalError);
	}

	// Weight info
	int wgtFlag = 0;
	int* vwgtPtr = NULL;
	int* adjwgtPtr = NULL;

	if (cellWeights.size()) {
		vwgtPtr = cellWeights.begin();
		wgtFlag += 2;       // Weights on vertices
	}
	if (faceWeights.size()) {
		adjwgtPtr = faceWeights.begin();
		wgtFlag += 1;       // Weights on edges
	}

	// Number of weights or balance constraints
	int nCon = 1;

	Field < real_t > processorWeights(Pstream::nProcs(), 1);
	Field<real_t> tpwgts(Pstream::nProcs(), 1. / nProcessors_);

	processorWeights = getProcessorWeights();

	tpwgts = processorWeights / sum(processorWeights);

	Info << "tpwgts " << tpwgts << endl;

	/*
	 * An array of size ncon that is used to specify the imbalance tolerance
	 * for each vertex weight, with 1 being perfect balance and nparts
	 * being perfect imbalance. A value of 1.05 for each of the ncon
	 * weights is recommended.
	 */

	// Imbalance tolerance
	Field<real_t> ubvec(nCon, parmetisImbalanceTolerance_);
	if (nProcessors_ == 1) {
		// If only one processor there is no imbalance.
		ubvec[0] = 1;
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	// output: cell -> processor addressing
	finalDecomp.setSize(nLocalCells[Pstream::myProcNo()]);

	// output: number of cut edges
	int edgeCut = 0;

	real_t itr(parmetisITR_);

//	List<int> localOptions(options);
//	localOptions[0]=1;
//	localOptions[1]=0;
//	localOptions[2]=15;
//	localOptions[3]=PARMETIS_PSR_COUPLED;

	/**
	 * http://dl.acm.org/citation.cfm?id=370498
	 */
	ParMETIS_V3_AdaptiveRepart(
			cellOffsets.begin(),    // vtxDist
			xadj.begin(), adjncy.begin(),
			vwgtPtr,                // vertexweights
			vwgtPtr, //idx_t * vsize, //?
			adjwgtPtr,              // edgeweights
			&wgtFlag, &numFlag, &nCon,
			&nProcessors_,          // nParts
			tpwgts.begin(), ubvec.begin(),
			&itr, //real_t * ipc2redist, //?
			const_cast<List<int>&>(options).begin(), &edgeCut,
			finalDecomp.begin(), &comm);




	// If we sent cells across make sure we undo it
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Receive back from next processor if I sent something
	if (nSendCells[Pstream::myProcNo()] > 0) {
		IPstream fromNextProc(Pstream::blocking, Pstream::myProcNo() + 1);

		List<int> nextFinalDecomp(fromNextProc);

		if (nextFinalDecomp.size() != nSendCells[Pstream::myProcNo()]) {
			FatalErrorIn("parMetisDecomp::decompose(..)")
					<< "Expected from processor " << Pstream::myProcNo() + 1
					<< " decomposition for " << nSendCells[Pstream::myProcNo()]
					<< " nCells but only received " << nextFinalDecomp.size()
					<< abort(FatalError);
		}

		append(nextFinalDecomp, finalDecomp);
	}

	// Send back to previous processor.
	if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo() - 1] > 0) {
		OPstream toPrevProc(Pstream::blocking, Pstream::myProcNo() - 1);

		int nToPrevious = nSendCells[Pstream::myProcNo() - 1];

		toPrevProc
				<< SubList<int>(finalDecomp, nToPrevious,
						finalDecomp.size() - nToPrevious);

		// Remove locally what has been sent
		finalDecomp.setSize(finalDecomp.size() - nToPrevious);
	}

	return edgeCut;
}

Foam::label topParentID(dynamicRefineExtFvMesh& m, label p) {
	label nextP = m.meshCutter().history().splitCells()[p].parent_;
	if (nextP < 0) {
		return p;
	} else {
		return topParentID(m, nextP);
	}
}

Foam::label Foam::parMetisDecompDynamic::decomposeBasedCommGraph(Field<int>& xadj,
	Field<int>& adjncy, const pointField& cellCentres, Field<int>& cellWeights,
		Field<int>& faceWeights, const List<int>& options, List<int>& finalDecomp,
		const labelList& fineToCoarse) {

	Info << "Using parMetisDecompDynamic::decomposeBasedCommGraph" << endl;


	int nparts = commGraph->getNumberOfNodes();
	List<int> finalsubDecomp;

	Field < float > processorWeights(nparts, 1);
	processorWeights = getCoarseProcessorWeights();

	subDecompose(xadj, adjncy, cellCentres, cellWeights, faceWeights, options,
			finalsubDecomp, nparts,  processorWeights, true);

	finalDecomp.setSize(cellCentres.size());


	for (int clusterID = 0; clusterID < nparts; clusterID++) {

		Info << "Processing cluster " << clusterID << endl;

		int sub_nparts = commGraph->getNumberOfProcsInNode(clusterID);
		// Connections
		Field<int> subAdjncy;

		// Offsets into adjncy
		Field<int> subXadj;
		scalarField cWeights;

		// decomposition options. 0 = use defaults
		List<int> options(3, 0);

		// cell weights (so on the vertices of the dual)
		Field<int> subCellWeights;

		// face weights (so on the edges of the dual)
		Field<int> subFaceWeights;

		if (coarseMapping) {



			std::map<int, std::vector<int> > coarseToFine;

			for (int i = 0; i < fineToCoarse.size(); i++) {

				std::map<int, std::vector<int> >::iterator it = coarseToFine.find(
						fineToCoarse[i]);
				if (it == coarseToFine.end()) {
					coarseToFine.insert(
							std::pair<int, std::vector<int> >(fineToCoarse[i],
									std::vector<int>()));
					coarseToFine[fineToCoarse[i]].push_back(i);

				} else {
					coarseToFine[fineToCoarse[i]].push_back(i);
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
			Info << "done coarseToFine" << endl;

			//Select the cells for each node
			labelHashSet cells;
			for (int f = 0; f < finalsubDecomp.size(); ++f) {

				if (finalsubDecomp[f] == clusterID) {
					std::vector<int>& childs = coarseToFine[f];
					for (int i = 0; i < childs.size(); ++i) {
						cells.insert(childs[i]);
					}
				}
			}


			fvMeshSubset submesh(*dynamic_cast<fvMesh*>(mesh_));
			submesh.setLargeCellSubset(cells, -1, true);
			fvMesh& clusterMesh = submesh.subMesh();


			dynamicRefineExtFvMesh& originalMesh =
					*dynamic_cast<dynamicRefineExtFvMesh*>(mesh_);


			const labelIOList& cellLevel = originalMesh.meshCutter().cellLevel();

			Map < label > coarseIDmap(100);

			labelList uniqueIndex(clusterMesh.nCells(), 0);

			label nCoarse = 0;

			forAll(clusterMesh.cells(), cellI)
			{
				int originalCell = submesh.cellMap()[cellI];

				if (cellLevel[originalCell] > 0) {
					uniqueIndex[cellI] = originalMesh.nCells()
							+ topParentID(originalMesh,
									originalMesh.meshCutter().history().parentIndex(
											originalCell));

				} else {
					uniqueIndex[cellI] = cellI;

				}

				if (coarseIDmap.insert(uniqueIndex[cellI], nCoarse)) {
					++nCoarse;
				}
			}



			labelList subfineToCoarse;
			pointField subcoarsePoints;
			scalarField subcoarseWeights;

			subfineToCoarse.resize(clusterMesh.nCells(), 0);
			subcoarsePoints.resize(nCoarse, vector::zero);// = pointField(nCoarse, vector::zero);
			subcoarseWeights.resize(nCoarse, 0.0);

			for (int cellI = 0; cellI < uniqueIndex.size(); ++cellI) {

				subfineToCoarse[cellI] = coarseIDmap[uniqueIndex[cellI]];

				// If 2D refinement (quadtree) is ever implemented, this '3'
				// should be set in general as the number of refinement
				// dimensions.
				int originalCell = submesh.cellMap()[cellI];


				label w = (1 << (3 * cellLevel[originalCell]));

				subcoarseWeights[subfineToCoarse[cellI]] += 1.0;

				//Seems fvSubSetMesh does nto construct C()
				subcoarsePoints[subfineToCoarse[cellI]] += originalMesh.C()[originalCell] / w;

			}

			std::map<int, std::vector<int> > subCoarseToFine;


			for (int i = 0; i < subfineToCoarse.size(); i++) {

				std::map<int, std::vector<int> >::iterator it = subCoarseToFine.find(
						subfineToCoarse[i]);
				if (it == subCoarseToFine.end()) {
					coarseToFine.insert(
							std::pair<int, std::vector<int> >(subfineToCoarse[i],
									std::vector<int>()));
					subCoarseToFine[subfineToCoarse[i]].push_back(i);

				} else {
					subCoarseToFine[subfineToCoarse[i]].push_back(i);
				}
			}


			CompactListList < label > coarseCellCells;
			calcCellCells(clusterMesh, subfineToCoarse, subcoarsePoints.size(),
					true,                       // use global cell labels
					coarseCellCells);

			const labelListList& globalCellCells = coarseCellCells();


			calcCSR(globalCellCells, subAdjncy, subXadj);

			cWeights = scalarField(subcoarseWeights.size(), 1);
			scalar minWeights = gMin(subcoarseWeights);

			if (cWeights.size() > 0) {
				if (minWeights <= 0) {
					WarningIn("metisDecomp::decompose"
							"(const pointField&, const scalarField&)")
							<< "Illegal minimum weight " << minWeights << endl;
				}

				if (cWeights.size() != globalCellCells.size()) {
					FatalErrorIn("parMetisDecomp::decompose"
							"(const pointField&, const scalarField&)")
							<< "Number of cell weights " << cWeights.size()
							<< " does not equal number of cells "
							<< globalCellCells.size() << exit(FatalError);
				}

				// Convert to integers.
				subCellWeights.setSize(cWeights.size());
				forAll(subCellWeights, i)
				{
					subCellWeights[i] = int(subcoarseWeights[i] / minWeights);
				}
			}

			Info << "subDecompose" << endl;

			List<int> finalsubDecomp;

			Field<float> processorWeights(sub_nparts, 1);
			processorWeights = getProcessorWeightsForNode(clusterID);

			subDecompose(subXadj, subAdjncy, subcoarsePoints,
					subCellWeights, subFaceWeights, options, finalsubDecomp, sub_nparts,
					processorWeights, false);

			forAll(finalsubDecomp, f)
			{


				int finalTargetRank = commGraph->getRankOfNodeIndex(clusterID, finalsubDecomp[f]);
				int originalCell = submesh.cellMap()[subCoarseToFine[f][0]];

				int originalCoarse = fineToCoarse[originalCell];


				finalDecomp[originalCoarse] = finalTargetRank;

			}


		} else {
			//Select the cells for each node
			labelHashSet cells;
			forAll(finalsubDecomp, f)
			{
				if (finalsubDecomp[f] == clusterID) {
					cells.insert(f);
				}
			}

			fvMeshSubset submesh(*dynamic_cast<fvMesh*>(mesh_));
			submesh.setLargeCellSubset(cells, -1, true);
			fvMesh& clusterMesh = submesh.subMesh();

			cWeights = scalarField(clusterMesh.cellCentres().size(), 1);

			Foam::calcDistributedCSR(clusterMesh, subAdjncy, subXadj);

			// Check for externally provided cellweights and if so initialise weights
			scalar minWeights = gMin(cWeights);

			if (cWeights.size() > 0) {
				if (minWeights <= 0) {
					WarningIn("metisDecomp::decompose"
							"(const pointField&, const scalarField&)")
							<< "Illegal minimum weight " << minWeights << endl;
				}

				if (cWeights.size() != clusterMesh.nCells()) {
					FatalErrorIn("parMetisDecomp::decompose"
							"(const pointField&, const scalarField&)")
							<< "Number of cell weights " << cWeights.size()
							<< " does not equal number of cells "
							<< clusterMesh.nCells() << exit(FatalError);
				}

				// Convert to integers.
				subCellWeights.setSize(cWeights.size());
				forAll(subCellWeights, i)
				{
					subCellWeights[i] = int(cWeights[i] / minWeights);
				}
			}

			// Do actual decomposition
			List<int> finalsubDecomp;

			Field<float> processorWeights(sub_nparts, 1);
			processorWeights = getProcessorWeightsForNode(clusterID);

			int error = subDecompose(subXadj, subAdjncy, clusterMesh.cellCentres(),
					cellWeights, faceWeights, options, finalsubDecomp, sub_nparts,
					processorWeights, false);

			if (error == -1)
				return error;

			finalDecomp.setSize(cellCentres.size());
			forAll(finalsubDecomp, f)
			{
				finalDecomp[submesh.cellMap()[f]] = commGraph->getRankOfNodeIndex(
						clusterID, finalsubDecomp[f]);
				//Pout << "(" << finalsubDecomp[f] << "," << clusterID* nparts + finalsubDecomp[f] << ")->" << submesh.cellMap()[f] << " ";
			}
		}
	}

}

//- Does prevention of 0 cell domains and calls parmetis.
Foam::label Foam::parMetisDecompDynamic::subDecompose(Field<int>& xadj,
	Field<int>& adjncy, const pointField& cellCentres, Field<int>& cellWeights,
	Field<int>& faceWeights, const List<int>& options, List<int>& finalDecomp,
	int nparts, Field<float>& processorWeights, bool clusterDecomp) {

	// C style numbering
	int numFlag = 0;

	// Number of dimensions
	int nDims = 3;
	int edgeCut = 0;

	if (cellCentres.size() != xadj.size() - 1) {
		FatalErrorIn("parMetisDecomp::decompose(..)") << "cellCentres:"
				<< cellCentres.size() << " xadj:" << xadj.size()
				<< abort(FatalError);
	}

	// Get number of cells on all processors
	List<int> nLocalCells(Pstream::nProcs());
	nLocalCells[Pstream::myProcNo()] = xadj.size() - 1;
	Pstream::gatherList(nLocalCells);
	Pstream::scatterList(nLocalCells);

        Info << "nLocalCells: " << nLocalCells << endl;

	std::vector<int> ranks;

	//exclude ranks with zero cells
	forAll(nLocalCells, procI)
	{
		if (nLocalCells[procI] > 0)
			ranks.push_back(procI);
	}

	int P = ranks.size();

	/**
	 * Create a MPI group with only non-zero cells ranks
	 * due to parmetis limitation
	 */

	// Get the rank and size in the original communicator
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the group of processes in MPI_COMM_WORLD
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	// Construct a group containing all of the prime ranks in world_group
	MPI_Group prime_group;
	MPI_Group_incl(world_group, ranks.size(), &(ranks[0]), &prime_group);

	// Create a new communicator based on the group
	MPI_Comm prime_comm;
	MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &prime_comm);

	int prime_rank = -1, prime_size = -1;
	// If this rank isn't in the new communicator, it will be
	// MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
	// MPI_Comm_size is erroneous
	if (MPI_COMM_NULL != prime_comm) {
		MPI_Comm_rank(prime_comm, &prime_rank);
		MPI_Comm_size(prime_comm, &prime_size);
	}

	List<int> nLocalCells2(P);
	std::vector<int> sub_ranks(Pstream::nProcs(), -1);

	int i = 0;
	forAll(nLocalCells, procI)
	{

		if (nLocalCells[procI] > 0) {
			sub_ranks[procI]= i;
			nLocalCells2[i++] = nLocalCells[procI];
		}

	}

	if (nLocalCells[Pstream::myProcNo()] > 0) {

		nLocalCells = nLocalCells2;

		// Get cell offsets.
		List<int> cellOffsets(P + 1);
		int nGlobalCells = 0;
		forAll(nLocalCells, procI)
		{
			cellOffsets[procI] = nGlobalCells;
			nGlobalCells += nLocalCells[procI];
		}
		cellOffsets[P] = nGlobalCells;

		// Convert pointField into the data type parMetis expects (float or double)
		Field < real_t > xyz(3 * cellCentres.size());
		int compI = 0;
		forAll(cellCentres, cellI)
		{
			const point& cc = cellCentres[cellI];
			xyz[compI++] = float(cc.x());
			xyz[compI++] = float(cc.y());
			xyz[compI++] = float(cc.z());
		}


		// Weight info
		int wgtFlag = 0;
		int* vwgtPtr = NULL;
		int* adjwgtPtr = NULL;

		if (cellWeights.size()) {
			vwgtPtr = cellWeights.begin();
			wgtFlag += 2;       // Weights on vertices
		}
		if (faceWeights.size()) {
			adjwgtPtr = faceWeights.begin();
			wgtFlag += 1;       // Weights on edges
		}

		// Number of weights or balance constraints
		int nCon = 1;

		Field < real_t > tpwgts(nparts, 1. / nparts);
		tpwgts = processorWeights / sum(processorWeights);

		//Pout << "tpwgts " << tpwgts << endl;

		/*
		 * An array of size ncon that is used to specify the imbalance tolerance
		 * for each vertex weight, with 1 being perfect balance and nparts
		 * being perfect imbalance. A value of 1.05 for each of the ncon
		 * weights is recommended.
		 */

		// Imbalance tolerance
		Field < real_t > ubvec(nCon, parmetisImbalanceTolerance_);
		if (nProcessors_ == 1) {
			// If only one processor there is no imbalance.
			ubvec[0] = 1;
		}


		// output: cell -> processor addressing
		List<int> currentDecomp;
		currentDecomp.setSize(nLocalCells[sub_ranks[Pstream::myProcNo()]]);



		/**
		 * Tell parmetis the current distribution of the cells
		 * This depends on wether we are decomposing for nodes or ranks inside nodes
		 * in the first iteration for node decomp we want the node partitions built from scratch
		 * for the rest, we tell parmetis how they are distributed
		 */

		forAll(currentDecomp, f)
		{

			if (clusterDecomp && !(nSharma::getInstance()->balanceEpisodeID() == 0)) {
				currentDecomp[f] = commGraph->getNodeID(Pstream::myProcNo());

			} else {
				currentDecomp[f] = sub_ranks[Pstream::myProcNo()];
			}
		}


		// output: number of cut edges

		real_t itr(parmetisITR_);

		List<int> coarseCommDecomp(currentDecomp);

		List<int> localOptions(options);
		localOptions[0] = 1;
		localOptions[1] = 1;
		localOptions[2] = 15;
		localOptions[3] = PARMETIS_PSR_UNCOUPLED;

                if (cellOffsets.begin() == NULL || xadj.begin() == NULL ||  adjncy.begin() == NULL){
                        Pout << "Something wrong in subCluster partitioning, skipping parMetis call" << endl;
                } else {

        		ParMETIS_V3_AdaptiveRepart(
				cellOffsets.begin(),    // vtxDist
				xadj.begin(), adjncy.begin(),
				vwgtPtr,                // vertexweights
				vwgtPtr, //idx_t * vsize, //?
				adjwgtPtr,              // edgeweights
				&wgtFlag, &numFlag, &nCon,
				&nparts,          // nParts
				tpwgts.begin(), ubvec.begin(),
				&itr, //real_t * ipc2redist, //?
				const_cast<List<int>&>(localOptions).begin(), &edgeCut,
				coarseCommDecomp.begin(), &prime_comm);


                }

		finalDecomp = coarseCommDecomp;

		MPI_Group_free(&world_group);
		MPI_Group_free(&prime_group);
		MPI_Comm_free(&prime_comm);

	}

	return edgeCut;
}


List<float>& Foam::parMetisDecompDynamic::getProcessorWeights() {
	return processorWeights_;
}

List<float> Foam::parMetisDecompDynamic::getCoarseProcessorWeights() {

	int nNodes = commGraph->getNumberOfNodes();
	List<float> coarse(nNodes, 0.0);


	for (int nodeID = 0; nodeID < nNodes; nodeID++){

		std::vector<int>& procs = commGraph->getProcsInNode(nodeID);


		for (int i = 0; i < procs.size(); i++){
			coarse[nodeID] += processorWeights_[procs[i]];
		}
	}

	return coarse;
}

List<float> Foam::parMetisDecompDynamic::getProcessorWeightsForNode(int node) {

	std::vector<int>& procs = commGraph->getProcsInNode(node);

	List<float> ws(procs.size(), 0.0);

	for (int i = 0; i < procs.size(); i++){
		ws[i] =  processorWeights_[procs[i]];
	}

	return ws;
}

void Foam::parMetisDecompDynamic::setProcessorWeights(
	List<float>& processorWeights) {

processorWeights_ = processorWeights;

}

void Foam::parMetisDecompDynamic::setParmetisImbalanceTolerance(float v) {

parmetisImbalanceTolerance_ = v;

}

void Foam::parMetisDecompDynamic::setParmetisITR(float v) {

parmetisITR_ = v;

}

