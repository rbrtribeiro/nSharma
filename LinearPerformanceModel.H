/*---------------------------------------------------------------------------*\
Class
 Foam::PerformanceModel

 Description

 SourceFiles
 PerformanceModel.C

 \*---------------------------------------------------------------------------*/
#ifndef LinearPerformanceModel_H
#define LinearPerformanceModel_H

#include "List.H"
#include "labelList.H"
#include "PerformanceModel.H"
#include "Partitioner.H"
#include "loadManager.H"

#include <cmath>
#include <vector>

namespace Foam {

class LinearPerformanceModel: public PerformanceModel {
protected:

	Partitioner* partitioner;

	List<float> Weights;
	std::vector<std::vector<float> > last_weights;

	List<float> TperCell;
	std::vector<std::vector<float> > last_TperCell;

	std::vector<float> window_weights;

	void print(std::vector< std::vector<double> > A);
	std::vector<double> gauss(std::vector< std::vector<double> > A);
	void calcWeightsDistAndTperCell();
	long getTotalCellsSend(Foam::labelList& decomp);
	long getTotalCellsReceived(Foam::labelList& decomp);
	void ensureMinMoved(std::vector<decompDetails*>& decomp);
	void reduceDoBalance(bool& doBalance);
	void calcCellsTransferedAndNonZero(std::vector<decompDetails*>& decomps);
	void windowAverage(List<float> raw,
			std::vector<std::vector<float> >& last_,
			List<float>& averaged);
public:
	LinearPerformanceModel(loadManagerParameters&);
	~LinearPerformanceModel();


	void PostLoadMigrationUpdate();
	bool Balance(labelList&);

};

}

#endif

