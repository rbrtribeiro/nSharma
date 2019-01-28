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




#include <stdlib.h>

#include "LinearLoadManagerModel.H"
#include <mpi.h>
#include "procedureProfPool.H"
#include "fvMeshDistribute.H"
#include "MPIfiedProcedure.H"
#include <math.h>
#include <cfloat>

using namespace Foam;


LinearLoadManagerModel::LinearLoadManagerModel(const nSharmaParameters& p,
		IOdictionary& dict, PerformanceModel* pm,
		CommGraph* cG) :
		PARAMS_(p), performanceModel(pm) {

	Weights.resize(Pstream::nProcs(), 1/float(Pstream::nProcs()));

	std::vector<float> tmp;
	tmp.resize(	Pstream::nProcs(), 1/float(Pstream::nProcs()));

	window_weights.reserve(p.window);

	float sum = 0;
	for (int i = p.window; i > 0; i--) {

		window_weights.push_back(log10(double(i)));

		sum +=window_weights.back();

		last_weights.push_back(tmp);
	}

	for (int i = 0; i < p.window; i++) {
		window_weights[i] = window_weights[i] / sum;

	}

	const word dynamicFvMeshTypeName(
			dict.lookup(
					"dynamicFvMesh"));

	if (dynamicFvMeshTypeName == "staticFvMesh")
		partitioner = new ParMetisPartitioner(p, cG);
	else {
		partitioner = new ParMetisRefinePartitioner(p,cG);
	}

	first = true;

}

LinearLoadManagerModel::~LinearLoadManagerModel() {

}

void Foam::LinearLoadManagerModel::print(std::vector<std::vector<double> > A) {
	int n = A.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n + 1; j++) {
			cout << A[i][j] << "\t";
			if (j == n - 1) {
				cout << "| ";
			}
		}
		cout << "\n";
	}
	cout << endl;
}

std::vector<double> Foam::LinearLoadManagerModel::gauss(
		std::vector<std::vector<double> > A) {
	int n = A.size();

	for (int i = 0; i < n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k < n + 1; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j < n + 1; j++) {
				if (i == j) {
					A[k][j] = 0;
				} else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	std::vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--) {
			A[k][n] -= A[k][i] * x[i];
		}
	}
	return x;
}

void Foam::LinearLoadManagerModel::PostLoadMigrationUpdate() {

}

void Foam::LinearLoadManagerModel::windowAverage(List<float> raw,
		std::vector<std::vector<float> >& last_,
		List<float>& r){

	label P = Pstream::nProcs();
	std::vector<float> std_raw;
	std_raw.resize(Pstream::nProcs(), 0.0f);
	for (int p = 0; p < P; p++) std_raw[p] = raw[p];


	//init window with first value
//	if (Foam::nSharma::balanceEpisodeID() == 0) {
//		for (uint w = 0; w < last_.size(); w++) {
//			last_[w] = std_raw;
//		}
//	} else {
	last_.push_back(std_raw);
//	}

	List<float> averaged;
	averaged.resize(Pstream::nProcs(), 0.0f);
	for (int p = 0; p < P; p++) {
		std::vector<std::vector<float> >::reverse_iterator rit =
				last_.rbegin();


		for (uint w = 0; w < Foam::nSharma::getInstance()->PARAMS().window; w++) {
//			Info << "p: " << p << " invw:" << w <<
//					" last_weights: " << (*rit)[p] << endl;

			averaged[p] += window_weights[w]*((*rit)[p]);
			rit++;
		}
	}

	r = averaged;

}

void Foam::LinearLoadManagerModel::calcWeightsDist() {

	performanceModel->update();


	if (Pstream::master()){


		label tCells;

		for (int i = 0; i < Pstream::nProcs(); i++) {
			tCells += nSharma::getInstance()->getCurrentCells(i);
		}

		List<float> tpercell = performanceModel->getTperCell();


		List<float> raw_weights;
		raw_weights.resize(Pstream::nProcs(), 0.0f);


		clockTime clock_;

		label P = Pstream::nProcs();
		std::vector < std::vector<double> > A(P + 1);

		for (int i = 0; i < P + 1; i++) {
			A[i].resize(P + 2);
			for (int j = 0; j < P + 2; j++) {
				A[i][j] = 0;

				if (i == j)
					A[i][j] = tpercell[i];

				if (i == P)
					A[i][j] = 1;
				if (j == P)
					A[i][j] = -1;
			}
		}

		A[P][P] = 0;
		A[P][P + 1] = tCells;

		//print (A);
		//gauss (A);

		std::vector<double> x(P);
		x = gauss(A);

		scalar elapsed = clock_.elapsedTime();
		Info << "Model_solving_time: " << elapsed << " seconds" << endl;

		// Print result
		//cout << "Result:\t";
		for (int i = 0; i < P; i++) {
			//cout << x[i] << " ";

			//divide calculated ideal solution by total number of cell to get
			//work distribution ratio for partitioner
			raw_weights[i] = x[i] / tCells;
			//pw[i] = 0.33;

		}



		if (first) {

			label P = Pstream::nProcs();
			std::vector<float> std_raw;
			std_raw.resize(Pstream::nProcs(), 0.0f);
			for (int p = 0; p < P; p++)
				std_raw[p] = raw_weights[p];

			for (uint w = 0; w < last_weights.size(); w++) {
				last_weights[w] = std_raw;
			}

			first = false;
		}



		Info << "raw-weights " << raw_weights << endl;
		Info << "previous-weights " << Weights << endl;
		windowAverage(raw_weights,last_weights, Weights);
		Info << "averaged-submitted-weights " << Weights << endl;

	}

	MPI_Bcast(&(Weights[0]), Pstream::nProcs(), MPI_FLOAT, Pstream::masterNo(),
				MPI_COMM_WORLD);


//	if (Pstream::master()){
//		Info << "PreviousEpisodeComputedCells " << endl;
//		for (int p = 0; p < Pstream::nProcs(); p++) {
//			Info << nSharma::getPreviousEpisodeComputedCells(p) << endl;
//		}
//	}


}

long Foam::LinearLoadManagerModel::getTotalCellsSend(Foam::labelList& decomp){

	long totalSend=0;
	for (int i = 0; i < Pstream::nProcs(); i++) {
			if (i != Pstream::myProcNo())
				totalSend += decomp[i];
		}

	return totalSend;
}

long Foam::LinearLoadManagerModel::getTotalCellsReceived(Foam::labelList& decomp){

	long totalReceived=0;
	for (int i = 0; i < Pstream::nProcs(); i++) {
				for (int j = 0; j < Pstream::nProcs(); j++) {
					if (i != j ){
						if (Pstream::myProcNo() == i)
							MPI_Send(&decomp[j], 1, MPI_INT, j, 0, MPI_COMM_WORLD);
						else if (Pstream::myProcNo() == j){
							int received=0;
							MPI_Recv(&received, 1, MPI_INT,i, 0,
												MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							totalReceived+=received;
						}
					}
				}
			}
	return totalReceived;
}

void Foam::LinearLoadManagerModel::ensureMinMoved(std::vector<decompDetails*>& decomps){

	return;

	for(std::vector<decompDetails*>::iterator it = decomps.begin();
			it != decomps.end();) {

		labelList decomp = (*it)->decomposition;
		Foam::labelList newDecompCounted = (*it)->count;

		for (int p = 0; p < Pstream::nProcs(); p++) {

			if (newDecompCounted[p] <= PARAMS_.minPercentageCellsMoved) {

				forAll(decomp, cellI){
					if ( decomp[cellI] == p)
						decomp[cellI] = Pstream::myProcNo();
				}

			}
		}

		newDecompCounted = fvMeshDistribute::countCells(decomp);

		bool doBalance = false;

		for (int p = 0; p < Pstream::nProcs(); p++) {

			if (p != Pstream::myProcNo() && newDecompCounted[p] > 0) {
				doBalance = true;
				break;
			}
		}

		reduceDoBalance(doBalance);

		if (!doBalance){
			delete *it;
			it = decomps.erase(it);
		}
		else {
			(*it)->count = newDecompCounted;
			it++;
		}
	}

}

void Foam::LinearLoadManagerModel::calcCellsTransferedAndNonZero(std::vector<decompDetails*>& decomps){

	for(std::vector<decompDetails*>::iterator it = decomps.begin();
			it != decomps.end(); ) {

				decompDetails& D = *(*it);

				D.totalSend = getTotalCellsSend(D.count);
				D.totalReceived = getTotalCellsReceived(D.count);

				long totalcellsTransfered = D.totalSend + D.totalReceived;

				long totalcellsTransfered__;
				MPI_Allreduce(&totalcellsTransfered, &totalcellsTransfered__, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
				totalcellsTransfered = totalcellsTransfered__;

				D.totalcellsTransfered = totalcellsTransfered;

				bool doBalance = ((D.totalReceived + D.count[Pstream::myProcNo()]) <= 1);

				reduceDoBalance(doBalance);

				if (doBalance) {
					delete *it;
					it = decomps.erase(it);
				} else {
					it++;
				}
	}

}

void Foam::LinearLoadManagerModel::reduceDoBalance(bool& doBalance){
	bool doBalance__;
	MPI_Allreduce(&doBalance, &doBalance__, 1, MPI_UNSIGNED_CHAR, MPI_LOR, MPI_COMM_WORLD);
	doBalance = doBalance__;

}

bool Foam::LinearLoadManagerModel::Balance(labelList& newD) {

	Info << "Using LinearLoadManagerModel" << endl;

	calcWeightsDist();

	std::vector<std::pair<int, labelList> > rawDecompCandidates;
	partitioner->decompose(Weights, rawDecompCandidates);

	std::vector<decompDetails*> decompCandidates;
	for(std::vector<std::pair<int, labelList> >::iterator it = rawDecompCandidates.begin();
				it != rawDecompCandidates.end(); ++it) {

		decompCandidates.push_back(new decompDetails((*it).second, (*it).first));

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


		List<float> TperCell = performanceModel->getTperCell();


		//previous approximated computing time
		float T = previousCells * TperCell[Pstream::myProcNo()];

		float T__;
		MPI_Allreduce(&T, &T__, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
				T=T__;

		int expected_iterations = PARAMS_.simulation_iterations -
				nSharma::getInstance()->iterationCount();
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

			float time = D.Tstar;
			if (time < bestTime){
				best = *it;
				bestTime = time;
			}
		}

		doBalance = best != NULL &&
				(T - (best->Tstar)) > (PARAMS_.minimal_gain * T);


		Pout << "previousCells " << previousCells <<
					" newCells " << best->totalReceived + best->count[Pstream::myProcNo()] <<
					 " T " << T <<
					 " Tstar " << best->Tstar <<
					 " totalSend " << best->totalSend <<
					 " totalReceived " << best->totalReceived
					 << endl;


		if (Pstream::master()) {

			Info << "\nBalance-episode-log: " << doBalance
					<< " = T:" << T
					<< " - Tstar:" << best->Tstar
					<< " > "
					<< "PARAMS_.minimal_gain: " << PARAMS_.minimal_gain
					<< " * T:" << T
					<< endl;

		}

		doBalance |= nSharma::getInstance()->forceNextBalance;

		if (doBalance){
			Info << "system-totalcellsTransfered " << best->totalcellsTransfered << endl;
		}


		newD = best->decomposition;

	}

	return doBalance;

}


