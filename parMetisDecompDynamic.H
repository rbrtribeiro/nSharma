/*---------------------------------------------------------------------------*
 Class
 Foam::parMetisDecompDynamic

 Description

 SourceFiles
 parMetisDecompDynamic.C

 \*---------------------------------------------------------------------------*/

#ifndef parMetisDecompDynamic_H
#define parMetisDecompDynamic_H

#include "parMetisDecomp.H"

#include "CommGraph.h"

//
//extern "C" {
//#   include "parmetis.h"
//}

namespace Foam {

/*---------------------------------------------------------------------------*\
                        Class parMetisDecompDynamic Declaration
 \*---------------------------------------------------------------------------*/

class parMetisDecompDynamic: public parMetisDecomp {
	// Private data

	//- Disallow default bitwise copy construct
	parMetisDecompDynamic(const parMetisDecompDynamic&);

	label decompose(Field<int>& xadj,
			Field<int>& adjncy, const pointField& cellCentres,
			Field<int>& cellWeights, Field<int>& faceWeights,
			const List<int>& options, List<int>& finalDecomp);

	label decomposeBasedCommGraph(Field<int>& xadj, Field<int>& adjncy,
			const pointField& cellCentres, Field<int>& cellWeights,
			Field<int>& faceWeights, const List<int>& options,
			List<int>& finalDecomp, const labelList& fineToCoarse);

	label subDecompose(Field<int>& xadj, Field<int>& adjncy,
			const pointField& cellCentres, Field<int>& cellWeights,
			Field<int>& faceWeights, const List<int>& options,
			List<int>& finalDecomp, int parts, Field<float>& processorWeights,
			bool clusterDecomp);

	List<float> processorWeights_;
	float parmetisImbalanceTolerance_;
	float parmetisITR_;
	bool first;

	polyMesh* mesh_;
	CommGraph* commGraph;
	bool coarseMapping;
	bool useCommGraph;

public:

	//- Runtime type information
	TypeName("parMetisDynamic");

	// Constructors

	//- Construct given the decomposition dictionary and mesh
	parMetisDecompDynamic(const dictionary& decompositionDict);

	parMetisDecompDynamic(const dictionary& decompositionDict, polyMesh& mesh,
			CommGraph* commGraph_, bool coarseMapping, bool useCommGraph);

	// Destructor

	~parMetisDecompDynamic() {
	}

	labelList decompose
	(
	    const polyMesh& mesh,
	    const labelList& fineToCoarse,
	    const pointField& coarsePoints,
	    const scalarField& coarseWeights
	);

	labelList decompose(const polyMesh& mesh, const pointField& cc,
			const scalarField& cWeights);
	labelList decompose(const labelListList& globalCellCells,
			const pointField& cellCentres, const scalarField& cWeights,
			const labelList& fineToCoarse);

	List<float>& getProcessorWeights();
	List<float> getCoarseProcessorWeights();
	List<float> getProcessorWeightsForNode(int node);

	void setProcessorWeights(List<float>& processorWeights);
	void setParmetisImbalanceTolerance(float v);
	void setParmetisITR(float v);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

#endif

// ************************************************************************* //
