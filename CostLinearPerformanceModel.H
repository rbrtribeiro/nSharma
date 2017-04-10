/*---------------------------------------------------------------------------*\
Class
 Foam::PerformanceModel

 Description

 SourceFiles
 PerformanceModel.C

 \*---------------------------------------------------------------------------*/
#ifndef CostPerformanceModel_H
#define CostPerformanceModel_H

#include "LinearPerformanceModel.H"

#include "List.H"
#include "labelList.H"
#include "Partitioner.H"
#include <cmath>
#include <vector>

namespace Foam {


class CostLinearPerformanceModel: public LinearPerformanceModel {

public:
	CostLinearPerformanceModel(loadManagerParameters& p);
	~CostLinearPerformanceModel();

	scalar lastTotalcellsTransfered;
	std::vector<float> cellsTransfered;
	std::vector<float> redistributeTimes;

	float linearRegressionM;
	float linearRegressionB;

	bool Balance(labelList&);
	void PostLoadMigrationUpdate();

};

}

#endif

