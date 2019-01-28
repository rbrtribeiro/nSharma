/*
 * PowerManagerModel.cpp
 *
 *  Created on: Jul 13, 2017
 *      Author: rr
 */

#include "PowerManagerModel.h"
#include "procedureProfPool.H"
#include "dynamicFvMesh.H"
#include <mpi.h>

#include "Pstream.H"

#include <nlopt/include/nlopt.h>

#include <cfloat>

using namespace Foam;

#define ALPHA 0.50
#define BETA 0.50

#define MY_NLOPT_TOL 1e-5
#define MY_NLOPT_EVALS 2500

ScalarizationUpperTMinT::ScalarizationUpperTMinT(PerformanceModel* pm,
		CommGraph* cg, PowerInterface* pi) :
		NLoptPowerManagerModel(pm, cg, pi) {

}

ScalarizationUpperTMinT::~ScalarizationUpperTMinT() {

}

//minimize alpha*Z + beta * SUM(W_p)
double objective_ScalarizationUpperTMinT(unsigned n, const double *x,
		double *grad, void *my_func_data) {

	time_constraint_data * data_timeC = (time_constraint_data*) my_func_data;

	/*
	 * SUM(W_p)
	 */

	double sumW = 0;
	for (int i = 0; i < n - 1; ++i) {

		sumW += x[i];

	}

	if (grad) {
		for (int i = 0; i < n - 1; ++i) {
			grad[i] = BETA;
		}

		grad[n - 1] = ALPHA;

	}

	return ALPHA * x[n - 1] + BETA * sumW;
}

//SUM(W_p) - W_cap <= 0
double constraint_ScalarizationUpperTMinT_MinTPowerManagerModel(unsigned n,
		const double *x, double *grad, void *data) {

	time_constraint_data * data_timeC = (time_constraint_data*) data;

	double sum = 0;
	for (int i = 0; i < n - 1; ++i) {

		sum += x[i];

	}

	if (grad) {
		for (int i = 0; i < n - 1; ++i) {

			grad[i] = 1;

		}

		grad[n - 1] = 0;

	}

	return sum - data_timeC[0].wcap;

}

// delta_p / W_p -T <= 0
double constraint_ScalarizationUpperTMinT_upperBoundTPowerManagerModel(
		unsigned n, const double *x, double *grad, void *data) {

	time_constraint_data * data_timeC = (time_constraint_data*) data;

	if (grad) {
		for (int i = 0; i < n - 1; ++i) {

			if (i == data_timeC->p) {
				grad[i] = -(data_timeC->delta[i] / (x[i] * x[i]));
			} else {
				grad[i] = 0;
			}

		}

		grad[n - 1] = 0;

	}

	return ((data_timeC->delta[data_timeC->p]) / x[data_timeC->p])
			- data_timeC->time;

}

// delta_p / W_p - z <= 0
double constraint_ScalarizationUpperTMinT_minmax(unsigned n, const double *x,
		double *grad, void *data) {

	time_constraint_data * data_timeC = (time_constraint_data*) data;

	if (grad) {
		for (int i = 0; i < n - 1; ++i) {
			if (i == data_timeC->p) {
				grad[i] = -(data_timeC->delta[i] / (x[i] * x[i]));
			} else {
				grad[i] = 0;
			}
		}
		grad[n - 1] = 0;
	}

	return ((data_timeC->delta[data_timeC->p]) / x[data_timeC->p]) - x[n - 1];

}

int ScalarizationUpperTMinT::solveNlopt(const std::vector<WattPair>& nodeBounds,
		const std::vector<int>& nodeCells,
		const std::vector<float>& nodeTperCells, std::vector<Watt>& nodeW,
		const void* arg1, const void* arg2) {

	float T  = *(float*)arg1;
	Watt wcap = *(Watt*) arg2;

	int nNodes = nodeBounds.size();

	std::vector<double> ub, lb;
	double gmax = 0;
	double gmin = FLT_MAX;
	for (int n = 0; n < nNodes; ++n) {

		ub.push_back(double(nodeBounds[n].upper));
		lb.push_back(double(nodeBounds[n].lower));

		if (double(nodeBounds[n].upper) > gmax)
			gmax = double(nodeBounds[n].upper);

		if (double(nodeBounds[n].lower) < gmin)
			gmin = double(nodeBounds[n].lower);

	}

	std::cout << "gmin: " << gmin << " gmax : " << gmax << std::endl;

	ub.push_back(gmax);
	lb.push_back(gmin);

	nlopt_opt opt;

	//NLOPT_GN_ISRES
	//NLOPT_LN_COBYLA
	//NLOPT_LD_SLSQP
	opt = nlopt_create(NLOPT_LN_COBYLA, nNodes + 1); /* algorithm and dimensionality */
	nlopt_set_upper_bounds(opt, &ub[0]);
	nlopt_set_lower_bounds(opt, &lb[0]);

	std::vector<time_constraint_data> data;
	std::vector<double> deltas;
	for (int i = 0; i < nNodes; ++i) {
		deltas.push_back(
				double(nodeCells[i] * nodeTperCells[i] * nodeBounds[i].upper));
	}

	data.resize(nNodes);
	for (int i = 0; i < nNodes; ++i) {

		data[i].p = i;
		data[i].wcap = wcap;
		data[i].time = T;
		data[i].delta = &deltas[0];

	}

	nlopt_set_min_objective(opt, objective_ScalarizationUpperTMinT, &data[0]);

	nlopt_add_inequality_constraint(opt,
			constraint_ScalarizationUpperTMinT_MinTPowerManagerModel, &data[0],
			MY_NLOPT_TOL);

	for (int i = 0; i < nNodes; ++i) {

		nlopt_add_inequality_constraint(opt,
				constraint_ScalarizationUpperTMinT_upperBoundTPowerManagerModel,
				&data[i], MY_NLOPT_TOL);

		nlopt_add_inequality_constraint(opt,
				constraint_ScalarizationUpperTMinT_minmax, &data[i], MY_NLOPT_TOL);
	}

	nlopt_set_xtol_rel(opt, MY_NLOPT_TOL);

	nlopt_set_maxeval(opt, MY_NLOPT_EVALS);

	/* some initial guess */
	std::vector<double> x(nNodes + 1, 0);
	for (int n = 0; n < nNodes; ++n) {
		x[n] = double(nodeBounds[n].upper);
	}

	x[nNodes] = double(nodeBounds[nNodes - 1].upper);

	double minf; /* the minimum objective value, upon return */

	int r = nlopt_optimize(opt, &x[0], &minf);

	for (int n = 0; n < nNodes; ++n) {
                nodeW[n] = Watt(round(x[n]));
		std::cout << "Estimated T: " << data[n].delta[n] / nodeW[n] << std::endl;

	}

        
	nlopt_destroy(opt);

	return r;
}
