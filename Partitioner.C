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
#include "nSharma.H"

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

ParMetisPartitioner::ParMetisPartitioner(const nSharmaParameters& p, CommGraph* cG,
		bool coarseMapping) :
		mesh(procedureProfPool::getMesh()), runTime(
				procedureProfPool::getRunTime()) {

	IOdictionary* decompositionDict = new IOdictionary(
			IOobject("decomposeParDict", runTime.system(), mesh,
					IOobject::MUST_READ, IOobject::NO_WRITE));

	decomposeEngine = new parMetisDecompDynamic(*decompositionDict,
			procedureProfPool::getMesh(),
			cG, coarseMapping,
			p.useCommGraphPartitioning);

	Info << "Instancing ParMetisPartitioner, CommGraph support: " <<
			p.useCommGraphPartitioning << endl;

	dynamic_cast<parMetisDecompDynamic*>(decomposeEngine)->setParmetisImbalanceTolerance(
			1.02);

	dynamic_cast<parMetisDecompDynamic*>(decomposeEngine)->setParmetisITR(10);

	std::vector<float> tol;
	//tol.push_back(1.02);
	tol.push_back(1.075);
	//tol.push_back(1.10);

	std::vector<float> itr;
	itr.push_back(0.001);
	//itr.push_back( 0.1);
	//itr.push_back(1);
	//itr.push_back(500);

	int id = 0;
	for (int i = 0; i < tol.size(); i++) {
		for (int j = 0; j < itr.size(); j++) {

			votes.push_back(std::pair<int, int>(0, id));

			std::pair<float, float> p = std::pair<float, float>(tol[i], itr[j]);
			std::pair<int, std::pair<float, float> > v = std::pair<int,
					std::pair<float, float> >(id++, p);
			IDtoDecomp.insert(v);

		}
	}

	nMostVoted = votes.size();

}

ParMetisPartitioner::~ParMetisPartitioner() {

}

void ParMetisPartitioner::decompChosen(int id) {

	votes[id].first++;

}

//Select paramters that score above average
void ParMetisPartitioner::selectMostUsed(std::vector<int>& result) {

	std::vector < std::pair<int, int> > votesSorted;
	votesSorted = votes;
	std::sort(votesSorted.begin(), votesSorted.end());

	//each score is a balance episode
	int totalScore = 0;
	for (int j = 0; j < votes.size(); j++) {
		totalScore += votes[j].first;
		Info << j << ": " << votes[j].first << " " << "Id: " << votes[j].second
				<< " " << IDtoDecomp[votes[j].second].first << ","
				<< IDtoDecomp[votes[j].second].second << endl;
	}

	//Each 3 iterations reduce one, keep initial for the first two
	//keep a minimum of 2
	if (nMostVoted > 2 && totalScore > 2)
		//if ((totalScore % 3) == 0)
		nMostVoted--;

	std::vector<std::pair<int, int> >::reverse_iterator rit =
			votesSorted.rbegin();

	for (uint w = 0; w < nMostVoted; w++) {
		int id = (*rit).second;
		result.push_back(id);

		Info << "ID: " << id << " " << "Score: " << (*rit).first << " "
				<< IDtoDecomp[id].first << "," << IDtoDecomp[id].second << endl;

		rit++;
	}

}

void ParMetisPartitioner::smoothWeightCells(scalarField& cellWeights) {

//	float sumW = 0;
//	forAll(cellWeights, i)
//	{
//		sumW += cellWeights[i];
//	}
//
//	float meanW = sumW / cellWeights.size();
	float minW = min(cellWeights);
	float maxW = max(cellWeights);

	float meanMAXMIN = (maxW + minW) / static_cast<double>(2);

	if (meanMAXMIN == 0) {
		FatalErrorIn("parMetisDecomp::decompose(..)") << "something wrong:"
				<< abort(FatalError);
	}

	double sumsqr = 0;
	double stddeviation;

	sumsqr += powf(maxW - meanMAXMIN, 2);
	sumsqr += powf(minW - meanMAXMIN, 2);

	stddeviation = sqrtf(sumsqr);

	float RSD = (stddeviation / meanMAXMIN) * 100.f;

	Pout << "min cellWeights: " << minW << " max cellWeights: " << maxW
			<< " RSD: " << RSD << " stddeviation: " << stddeviation
			/*<< " meanW: " << meanW */<< " meanMAXMIN: " << meanMAXMIN << endl;

	if (RSD > 120) {
		float cap = meanMAXMIN + 0.5 * stddeviation;
		Pout << "Smoothing weights by : " << cap << endl;
		int count = 0;
		forAll(cellWeights, i)
		{
			if (cellWeights[i] > cap) {
				cellWeights[i] = cap;
				count++;
			}

		}
		Pout << "Smoothed : " << count << endl;
	}
}

void ParMetisPartitioner::decompose(List<float>& ratios,
		std::vector<std::pair<int, labelList> >& decomps) {

	Info << "New decomposition with ParMetisPartitioner" << endl;

	((parMetisDecompDynamic*) decomposeEngine)->setProcessorWeights(ratios);

	scalarField weights(mesh.cellCentres().size(), 1);

	parMetisDecompDynamic* partEngine =
			dynamic_cast<parMetisDecompDynamic*>(decomposeEngine);

	ADD_SECTION (parMetis);

	std::vector<int> toRequest;
	selectMostUsed(toRequest);

	float itr, tol;
	for (int k = 0; k < toRequest.size(); k++) {
		tol = IDtoDecomp[toRequest[k]].first;
		itr = IDtoDecomp[toRequest[k]].second;

		partEngine->setParmetisImbalanceTolerance(tol);
		partEngine->setParmetisITR(itr);
		decomps.push_back(
				std::pair<int, labelList>(toRequest[k],
						decomposeEngine->decompose(mesh, mesh.cellCentres(),
								weights)));
	}

	END_SECTION(parMetis);

}

/*---------------------------------------------------------------------------*/

ParMetisRefinePartitioner::ParMetisRefinePartitioner(const nSharmaParameters& p,
		CommGraph* cG) :
		ParMetisPartitioner(p, cG, true) {

}

ParMetisRefinePartitioner::~ParMetisRefinePartitioner() {

}

void ParMetisRefinePartitioner::decompose(List<float>& ratios,
		std::vector<std::pair<int, labelList> >& decomps) {

	Info << "New decomposition with ParMetisRefinePartitioner" << endl;

	((parMetisDecompDynamic*) decomposeEngine)->setProcessorWeights(ratios);

	dynamicRefineExtFvMesh& extMesh = *(dynamicRefineExtFvMesh*) &mesh;

	/*RR: With dynamicRefineExtFvMesh and ParMetisRefinePartitioner,
	 * cell-to-coarse map is required and it's built in
	 * mesh.update(). If mesh.update() not triggered:
	 */
	if (!extMesh.changing())
		extMesh.buildCoarseMap();

	parMetisDecompDynamic* partEngine =
			dynamic_cast<parMetisDecompDynamic*>(decomposeEngine);

	ADD_SECTION (parMetis);

	std::vector<int> toRequest;
	selectMostUsed(toRequest);

	//smoothWeightCells(extMesh.getCoarseWeights());

	float itr, tol;
	for (int k = 0; k < toRequest.size(); k++) {
		tol = IDtoDecomp[toRequest[k]].first;
		itr = IDtoDecomp[toRequest[k]].second;

		partEngine->setParmetisImbalanceTolerance(tol);
		partEngine->setParmetisITR(itr);
		decomps.push_back(
				std::pair<int, labelList>(toRequest[k],
						decomposeEngine->decompose(extMesh,
								extMesh.getLocalIndex(),
								extMesh.getCoarsePoints(),
								extMesh.getCoarseWeights())));
	}

	END_SECTION(parMetis);

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
