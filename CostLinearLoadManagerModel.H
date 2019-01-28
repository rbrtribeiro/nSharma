/*---------------------------------------------------------------------------*\
Class
 Foam::PerformanceModel

 Description

 SourceFiles
 PerformanceModel.C

 \*---------------------------------------------------------------------------*/
#ifndef CostPerformanceModel_H
#define CostPerformanceModel_H

#include "LinearLoadManagerModel.H"

#include "List.H"
#include "labelList.H"
#include "Partitioner.H"
#include <cmath>
#include <vector>

namespace Foam {


class CostLinearLoadManagerModel: public LinearLoadManagerModel {

public:
	CostLinearLoadManagerModel(const nSharmaParameters& p,
			IOdictionary dict, PerformanceModel* pm,
			CommGraph* cG);
	~CostLinearLoadManagerModel();

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

