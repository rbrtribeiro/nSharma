/*---------------------------------------------------------------------------*\
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

#include "loadManager.H"
#include "MPIfiedProcedure.H"

#include <mpi.h>

#include "decompositionMethod.H"
#include "procedureProfPool.H"
//#include "parMetisDecompDynamic.H"
#include "CostLinearPerformanceModel.H"
#include "clockTime.H"
#include "dynamicFvMesh.H"

using namespace Foam;

class loadManager;

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <cassert>
#include <cfloat>

void linearRegression(const std::vector<float>& x, const std::vector<float>& y,
		float& m, float& b) {

	assert(x.size() == y.size());

    const size_t n    = x.size();
    const float s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const float s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const float s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const float s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    m    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

    b = (s_y/y.size())-m*(s_x/n);

}

float linearRegressionPredict(const float x, const float m, const float b){
	return m*x+b;
}

CostLinearPerformanceModel::~CostLinearPerformanceModel() {

}

CostLinearPerformanceModel::CostLinearPerformanceModel(loadManagerParameters& p) :
		LinearPerformanceModel(p) {

	redistributeTimes.reserve(10);
	cellsTransfered.reserve(10);
	lastTotalcellsTransfered=0;
	linearRegressionM=0;
	linearRegressionB=0;

}

void Foam::CostLinearPerformanceModel::PostLoadMigrationUpdate() {

	procedureProfPool::mapType& allInfo_ =
			procedureProfPool::getThePoolInstance()->map();

	if (redistributeTimes.size() < 10 && lastTotalcellsTransfered > 0) {



		float redistTime= float(allInfo_.find("local_redistribute")->second->totalTime());
		redistributeTimes.push_back(redistTime);

		//Info << "lastRedistributeTime " << redistTime << endl;
		//Info << "lastTotalcellsTransfered " << lastTotalcellsTransfered << endl;


		cellsTransfered.push_back(float(lastTotalcellsTransfered));



		Pout << "Base redistribute times: ";
		for (int i=0;  i < redistributeTimes.size(); i++){
			Pout << redistributeTimes[i] << " ";
		}
		Pout << endl;

		Pout << "Base cellsTransfered: ";
		for (int i=0;  i < cellsTransfered.size(); i++){
			Pout << cellsTransfered[i] << " ";
		}
		Pout << endl;





		if (redistributeTimes.size() > 2)
			linearRegression(cellsTransfered, redistributeTimes, linearRegressionM,
					linearRegressionB);

	}



}

bool Foam::CostLinearPerformanceModel::Balance(labelList& newD) {

	Info << "Using CostLinearPerformanceModel" << endl;

	calcWeightsDistAndTperCell();

	std::vector<labelList> rawDecompCandidates;
	partitioner->decompose(Weights, rawDecompCandidates);

	std::vector<decompDetails*> decompCandidates;
	for(std::vector<labelList>::iterator it = rawDecompCandidates.begin();
				it != rawDecompCandidates.end(); ++it) {

		decompCandidates.push_back(new decompDetails(*it));

	}

	ensureMinMoved(decompCandidates);
	bool doBalance = decompCandidates.size() > 0;

	if (!doBalance)
		Info << "\nBalance-episode-log: Very few cells required to move" << endl;

	if (doBalance) {

		calcCellsTransferedAndNonZero(decompCandidates);

		doBalance = decompCandidates.size() > 0;

		if (!doBalance){
			Info << "\nBalance-episode-log: Zero cells is non-sense, aborting" << endl;
			return doBalance;
		}

		label previousCells;
		previousCells = procedureProfPool::getMesh().nCells();

		//previous approximated computing time
		float T = previousCells * TperCell[Pstream::myProcNo()];

		float T__;
		MPI_Allreduce(&T, &T__, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
				T=T__;


		//TODO: HACK
		int expected_iterations = PARAMS_.simulation_iterations /*- loadManager::iterationCount()*/;
		expected_iterations = expected_iterations > 0 ? expected_iterations : 1;

		for(std::vector<decompDetails*>::iterator it = decompCandidates.begin();
				it != decompCandidates.end(); ++it) {

			decompDetails& D = *(*it);

			//new approximated computing time
			float Tstar = (D.totalReceived + D.count[Pstream::myProcNo()])
					* TperCell[Pstream::myProcNo()];

			float Tstar__;
			MPI_Allreduce(&Tstar, &Tstar__, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
			D.Tstar = Tstar__;

			//if first iterations, last migration cost is zero
			// TODO: do an initial guess?
			D.Tmigrt=linearRegressionPredict(D.totalcellsTransfered,
					linearRegressionM, linearRegressionB);

			if (D.Tmigrt <= 0){
				if (redistributeTimes.size() > 0 ){
					D.Tmigrt = redistributeTimes[0]/cellsTransfered[0] * D.totalcellsTransfered;
				}
				else D.Tmigrt = 0.0;
			}

			float Tmigrt__;
			MPI_Allreduce(&(D.Tmigrt), &Tmigrt__, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
			D.Tmigrt = Tmigrt__;

			Info << "Decomposition-candidate: "
						<< " = T:" << T
						<< " - Tstar:" << D.Tstar
						<< " + Tmigrt:" << D.Tmigrt
						<< ", totalcellsTransfered:" << D.totalcellsTransfered
						<< endl;

		}

		decompDetails* best=NULL;
		float bestTime=FLT_MAX;
		for(std::vector<decompDetails*>::iterator it = decompCandidates.begin();
						it != decompCandidates.end(); ++it) {

			decompDetails& D = *(*it);

			float time = D.Tstar + (D.Tmigrt/float(expected_iterations));
			if (time < bestTime){
				best = *it;
				bestTime = time;
			}
		}

		doBalance = best != NULL &&
				(T - (best->Tstar + (best->Tmigrt/float(expected_iterations)))) > (PARAMS_.minimal_gain * T);

		if (doBalance) lastTotalcellsTransfered = best->totalcellsTransfered;

		Pout << "previousCells " << previousCells <<
					" newCells " << best->totalReceived + best->count[Pstream::myProcNo()] <<
					 " T " << T <<
					 " Tstar " << best->Tstar <<
					 " Tmigrt " << best->Tmigrt <<
					 " totalSend " << best->totalSend <<
					 " totalReceived " << best->totalReceived
					 << endl;

		Info << "system totalcellsTransfered " << best->totalcellsTransfered << endl;

		if (Pstream::master()) {

			Info << "\nBalance-episode-log: " << doBalance
					<< " = T:" << T
					<< " - Tstar:" << best->Tstar
					<< " + Tmigrt:" << best->Tmigrt
					<< " / expected_iterations:" << expected_iterations
					<< " > "
					<< "PARAMS_.minimal_gain: " << PARAMS_.minimal_gain
					<< " * T:" << T
					<< endl;

		}

		newD = best->decomposition;

	}

	return doBalance;

}

