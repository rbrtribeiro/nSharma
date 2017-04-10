/*--------------------------Foam::parMetisDecomp(decompositionDict, mesh) Foam::parMetisDecomp(decompositionDict, mesh) Foam::parMetisDecomp(decompositionDict, mesh) -------------------------------------------------*\
  =========                 |
 \\      /  F ield         | foam-extend: Open Source CFD
 \\    /   O peration     |
 \\  /    A nd           | For copyright notice see file Copyright
 \\/     M anipulation  |
 -------------------------------------------------------------------------------
 License
 This file is part of foam-extend.

 foam-extend is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation, either version 3 of the License, or (at your
 option) any later version.

 foam-extend is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

 \*---------------------------------------------------------------------------*/

#include "procedureProfPool.H"
#include "Partitioner.H"
#include "parMetisDecompDynamic.H"
#include "dynamicRefineExtFvMesh.H"
//#include "ptscotchDecomp.H"
#include "loadManager.H"

/*This parameter describes the ratio of inter-processor communication time compared to data redistribution
time. It should be set between 0.000001 and 1000000.0. If ITR is set high, a repartitioning
with a low edge-cut will be computed. If it is set low, a repartitioning that requires little data redistribution
will be computed. Good values for this parameter can be obtained by dividing inter-processor
communication time by data redistribution time. Otherwise, a value of 1000.0 is recommended.

TOL: 1.02, 1.05, 1.10
ITR: 0.001, 0.1, 1, 10

*/


using namespace Foam;

Partitioner::~Partitioner() {

}

ParMetisPartitioner::ParMetisPartitioner() :
		mesh(procedureProfPool::getMesh()), runTime(
				procedureProfPool::getRunTime()) {

	IOdictionary* decompositionDict = new IOdictionary(
			IOobject("decomposeParDict", runTime.system(), mesh,
					IOobject::MUST_READ, IOobject::NO_WRITE));

	decomposeEngine = new parMetisDecompDynamic(*decompositionDict, mesh);
	dynamic_cast<parMetisDecompDynamic*>(decomposeEngine)->
			setParmetisImbalanceTolerance(1.02);

	dynamic_cast<parMetisDecompDynamic*>(decomposeEngine)->
				setParmetisITR(10);

}

ParMetisPartitioner::~ParMetisPartitioner() {

}

void ParMetisPartitioner::decompose(List<float>& ratios, std::vector<labelList>& decomps) {

	Info << "New decomposition with ParMetisPartitioner" << endl;

	((parMetisDecompDynamic*) decomposeEngine)->setProcessorWeights(ratios);

	scalarField weights(mesh.cellCentres().size(), 1);

	parMetisDecompDynamic* partEngine = dynamic_cast<parMetisDecompDynamic*>(decomposeEngine);

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(0.001);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(0.1);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(1);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(10);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(0.001);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(0.1);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(1);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(10);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(0.001);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(0.1);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(1);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(10);
	decomps.push_back(decomposeEngine->decompose(mesh.cellCentres(), weights));

}

/*---------------------------------------------------------------------------*/

ParMetisRefinePartitioner::ParMetisRefinePartitioner() :
		ParMetisPartitioner() {

}

ParMetisRefinePartitioner::~ParMetisRefinePartitioner() {

}

void ParMetisRefinePartitioner::decompose(List<float>& ratios,
		std::vector<labelList>& decomps) {

	Info << "New decomposition with ParMetisRefinePartitioner" << endl;

	((parMetisDecompDynamic*) decomposeEngine)->setProcessorWeights(ratios);

	dynamicRefineExtFvMesh& extMesh = *(dynamicRefineExtFvMesh*) &mesh;


	/*RR: With dynamicRefineExtFvMesh and ParMetisRefinePartitioner,
	 * cell-to-coarse map is required and it's built in
	 * mesh.update(). If mesh.update() not triggered:
	 */
	if (!extMesh.changing())
		extMesh.buildCoarseMap();

	parMetisDecompDynamic* partEngine = dynamic_cast<parMetisDecompDynamic*>(decomposeEngine);

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(0.001);

		decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
				extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(0.1);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(1);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.02);
	partEngine->setParmetisITR(10);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(0.001);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(0.1);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(1);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.05);
	partEngine->setParmetisITR(10);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(0.001);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(0.1);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(1);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

	partEngine->setParmetisImbalanceTolerance(1.10);
	partEngine->setParmetisITR(10);

	decomps.push_back(decomposeEngine->decompose( extMesh.getLocalIndex(),
			extMesh.getCoarsePoints(), extMesh.getCoarseWeights()));

}

/*---------------------------------------------------------------------------*/
//
//PTscotchPartitioner::PTscotchPartitioner() :
//		mesh(procedureProfPool::getMesh()), runTime(
//				procedureProfPool::getRunTime()) {
//
//	IOdictionary* decompositionDict = new IOdictionary(
//			IOobject("decomposeParDict", runTime.system(), mesh,
//					IOobject::MUST_READ, IOobject::NO_WRITE));
//
//	decomposeEngine = new ptscotchDecomp(*decompositionDict);
//
//}
//
//PTscotchPartitioner::~PTscotchPartitioner() {
//
//}
//
//void PTscotchPartitioner::decompose(List<float>& ratios, labelList& newD) {
//
//
//	Info << "New decomposition with PTscotchPartitioner" << endl;
//
//	scalarField weights(mesh.cellCentres().size(), 1);
//	newD = decomposeEngine->decompose(mesh, mesh.cellCentres(), weights);
//
//}
