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
#include "fvMeshDistribute.H"

#include "PerformanceModel.H"
#include "Partitioner.H"
#include "nSharma.H"

#include <cmath>
#include <vector>

namespace Foam {

class LinearLoadManagerModel: public LoadManagerModel {
protected:

	class decompDetails {
	public:
		int index;
		labelList decomposition;
		labelList count;
		long totalReceived;
		long totalSend;
		long totalcellsTransfered;
		bool apply;
		float Tstar;
		float Tmigrt;

		decompDetails(labelList& newD, int index_) {
			index = index_;
			decomposition = newD;
			count = fvMeshDistribute::countCells(decomposition);
		}

	};

	const nSharmaParameters& PARAMS_;

	Partitioner* partitioner;
	PerformanceModel* performanceModel;

	List<float> Weights;
	std::vector<std::vector<float> > last_weights;

	std::vector<float> window_weights;

	bool first;

	void print(std::vector<std::vector<double> > A);

	std::vector<double> gauss(std::vector<std::vector<double> > A);

	void calcWeightsDist();

	long getTotalCellsSend(Foam::labelList& decomp);

	long getTotalCellsReceived(Foam::labelList& decomp);

	void ensureMinMoved(std::vector<decompDetails*>& decomp);

	void reduceDoBalance(bool& doBalance);

	void calcCellsTransferedAndNonZero(std::vector<decompDetails*>& decomps);

	void windowAverage(List<float> raw, std::vector<std::vector<float> >& last_,
			List<float>& averaged);

public:

	LinearLoadManagerModel(const nSharmaParameters&, IOdictionary& dict,
			PerformanceModel* pm, CommGraph* cG);

	~LinearLoadManagerModel();

	void PostLoadMigrationUpdate();

	bool Balance(labelList&);

};

}

#endif

