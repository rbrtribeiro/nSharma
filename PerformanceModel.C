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

#include "nSharma.H"
#include "MPIfiedProcedure.H"
#include "procedureProfPool.H"
#include "PerformanceModel.H"

#include "Pstream.H"

using namespace Foam;

TperCellPerformanceModel::TperCellPerformanceModel(int window) {

	TperCell.resize(Pstream::nProcs(), 1.0e-05);

	if (window > 1) {

		window_weights.reserve(window);

		std::vector<float> tmp2;
		tmp2.resize(Pstream::nProcs(), 1.0e-05);

		float sum = 0;
		for (int i = window; i > 0; i--) {

			window_weights.push_back(log10(double(i)));

			sum += window_weights.back();

			last_TperCell.push_back(tmp2);

		}

		for (int i = 0; i < window; i++) {
			window_weights[i] = window_weights[i] / sum;

		}

	}

	first = true;

}

TperCellPerformanceModel::~TperCellPerformanceModel() {

}

void TperCellPerformanceModel::windowAverage(List<float> raw,
		std::vector<std::vector<float> >& last_, List<float>& r) {

	if (Foam::nSharma::getInstance()->PARAMS().window <= 1) {
		r = raw;
		return;
	}

	label P = Pstream::nProcs();
	std::vector<float> std_raw;
	std_raw.resize(Pstream::nProcs(), 0.0f);
	for (int p = 0; p < P; p++)
		std_raw[p] = raw[p];

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
		std::vector<std::vector<float> >::reverse_iterator rit = last_.rbegin();

		for (uint w = 0; w < Foam::nSharma::getInstance()->PARAMS().window;
				w++) {
//			Info << "p: " << p << " invw:" << w <<
//					" last_weights: " << (*rit)[p] << endl;

			averaged[p] += window_weights[w] * ((*rit)[p]);
			rit++;
		}
	}

	r = averaged;

}

List<float> Foam::TperCellPerformanceModel::getTperCell() {
	return TperCell;
}

void Foam::TperCellPerformanceModel::update() {

	if (Pstream::master()) {

		List<float> raw_tpercell;
		raw_tpercell.resize(Pstream::nProcs(), 0.0f);

		std::list<Operation*> targetIDs;
		nSharma::loadDataType& loadData = nSharma::getInstance()->getLoadData();

		procedureProfPool::getThePoolInstance()->getOpsIDsOfType(targetIDs);

		unsigned int nOps = targetIDs.size();

		if (nOps <= 0)
			return;

		float v;

		for (int i = 0; i < Pstream::nProcs(); i++) {

			for (std::list<Operation*>::iterator it = targetIDs.begin();
					it != targetIDs.end(); ++it) {

				v = (loadData.at((*it)->id0()).at(i)).totalTime_;

				raw_tpercell[i] += v;
			}

			raw_tpercell[i] /=
					nSharma::getInstance()->getPreviousEpisodeComputedCells(i);

		}

		Info << "raw-tpercell " << raw_tpercell << endl;
		Info << "previous-tpercell " << TperCell << endl;

		int P = Pstream::nProcs();
		//init last_TperCell window with mean of the TperCell first episode
		if (first) {

			/**
			 *
			 */

//			float TperCell_mean = 0;
//			label P = Pstream::nProcs();
//			for (int i = 0; i < P; i++) {
//				TperCell_mean += raw_tpercell[i];
//			}
//
//			TperCell_mean /= P;
//
//			for (uint w = 0; w < last_TperCell.size(); w++) {
//				last_TperCell[w] = std::vector<float>(P, TperCell_mean);
//			}

			/**
			 *
			 */

			label P = Pstream::nProcs();
			std::vector<float> std_raw;
			std_raw.resize(Pstream::nProcs(), 0.0f);
			for (int p = 0; p < P; p++)
				std_raw[p] = raw_tpercell[p];

			for (uint w = 0; w < last_TperCell.size(); w++) {
				last_TperCell[w] = std_raw;
			}



			first = false;
		}

		windowAverage(raw_tpercell, last_TperCell, TperCell);

		Info << "averaged-TperCell " << TperCell << endl;

//		for (uint w = 0; w < last_TperCell.size(); w++) {
//			Info << "last_TperCell ";
//			for (int i = 0; i < P; i++) {
//				Info << last_TperCell[w][i] << " ";
//			}
//			Info << endl;
//		}

	}

	MPI_Bcast(&(TperCell[0]), Pstream::nProcs(), MPI_FLOAT, Pstream::masterNo(),
			MPI_COMM_WORLD);

}
