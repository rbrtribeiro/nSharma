/*
 * PowerManagerModel.h
 *
 *  Created on: Jul 13, 2017
 *      Author: rr
 */

#ifndef POWERMANAGERMODEL_H_
#define POWERMANAGERMODEL_H_

#include "PerformanceModel.H"
#include "CommGraph.h"
#include "PowerInterface.h"

namespace Foam {

typedef struct {
	int p;
	float time;
	float wcap;
	double* delta;
} time_constraint_data;

class PowerManagerModel {
public:
	virtual void assess(bool enabled)=0;

};

class NLoptPowerManagerModel: public PowerManagerModel {
public:

	PerformanceModel* performanceModel;
	CommGraph* commGraph;
	PowerInterface* powerInterface;
	List<float> tpercell;
	bool tpercell_set;

	NLoptPowerManagerModel(PerformanceModel* pm, CommGraph* cg,
			PowerInterface* pi);

	virtual ~NLoptPowerManagerModel();

	const List<float>& getTperCell();


	virtual void assess(bool enabled);

	virtual int solveNlopt(const std::vector<WattPair>& nodeBounds,
			const std::vector<int>& nodeCells,
			const std::vector<float>& nodeTperCells, std::vector<Watt>& nodeW,
			const void* arg1 = NULL, const void* arg2 = NULL) =0;

};

class upperBoundTPowerManagerModel: public NLoptPowerManagerModel {
public:

	upperBoundTPowerManagerModel(PerformanceModel* pm, CommGraph* cg,
			PowerInterface* pi);
	virtual ~upperBoundTPowerManagerModel();

	virtual int solveNlopt(const std::vector<WattPair>& nodeBounds,
			const std::vector<int>& nodeCells,
			const std::vector<float>& nodeTperCells, std::vector<Watt>& nodeW,
			const void* arg1 = NULL, const void* arg2 = NULL);

};

class MinTPowerManagerModel: public NLoptPowerManagerModel {
public:

	MinTPowerManagerModel(PerformanceModel* pm, CommGraph* cg,
			PowerInterface* pi);

	virtual ~MinTPowerManagerModel();

	virtual int solveNlopt(const std::vector<WattPair>& nodeBounds,
			const std::vector<int>& nodeCells,
			const std::vector<float>& nodeTperCells, std::vector<Watt>& nodeW,
			const void* arg1 = NULL, const void* arg2 = NULL);

};

class ScalarizationUpperTMinT: public NLoptPowerManagerModel {
public:

	ScalarizationUpperTMinT(PerformanceModel* pm, CommGraph* cg,
			PowerInterface* pi);

	virtual ~ScalarizationUpperTMinT();

	virtual int solveNlopt(const std::vector<WattPair>& nodeBounds,
			const std::vector<int>& nodeCells,
			const std::vector<float>& nodeTperCells, std::vector<Watt>& nodeW,
			const void* arg1 = NULL, const void* arg2 = NULL);

};

}

#endif /* POWERMANAGERMODEL_H_ */
