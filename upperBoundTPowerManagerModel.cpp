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

using namespace Foam;

upperBoundTPowerManagerModel::upperBoundTPowerManagerModel(PerformanceModel* pm,
		CommGraph* cg, PowerInterface* pi) : NLoptPowerManagerModel(pm, cg, pi) {

}

upperBoundTPowerManagerModel::~upperBoundTPowerManagerModel() {

}


// minimize SUM(W_p)

double objective_upperBoundTPowerManagerModel(unsigned n, const double *x, double *grad,
		void *my_func_data) {

	double sum = 0;
	for (int i = 0; i < n; ++i) {

		sum += x[i];

	}

//	if (grad) {
//		for (int i = 0; i < n; ++i) {
//
//			grad[i] = 1;
//
//		}
//
//	}

	return sum;

}

// delta_p / W_p -T <= 0
double constraint_upperBoundTPowerManagerModel(unsigned n, const double *x, double *grad, void *data) {

	time_constraint_data * data_timeC = (time_constraint_data*) data;

//	if (grad) {
//		for (int i = 0; i < n; ++i) {
//			if (i == data_timeC->p) {
//				grad[i] = -(data_timeC->delta[i] / (x[i] * x[i]));
//			} else {
//				grad[i] = 0;
//			}
//
//		}
//
//	}

	return ((data_timeC->delta[data_timeC->p]) / x[data_timeC->p]) - data_timeC->time;

}

//SUM(W_p) - W_cap <= 0
double constraint_upperBoundTPowerManagerModel2(unsigned n, const double *x,
		double *grad, void *data) {

	time_constraint_data * data_timeC = (time_constraint_data*) data;

	double sum = 0;
	for (int i = 0; i < n; ++i) {

		sum += x[i];

	}

//	if (grad) {
//		for (int i = 0; i < n; ++i) {
//
//			grad[i] = 1;
//
//		}
//
//	}

	return sum - data_timeC[0].wcap;

}

int upperBoundTPowerManagerModel::solveNlopt(const std::vector<WattPair>& nodeBounds,
		const std::vector<int>& nodeCells,
		const std::vector<float>& nodeTperCells,
		std::vector<Watt>& nodeW,
		const void* arg1,
		const void* arg2) {

	float T  = *(float*)arg1;
	Watt wcap = *(Watt*) arg2;


	int nNodes = nodeBounds.size();

	std::vector<double> ub, lb;
	for (int n = 0; n < nNodes; ++n) {
		ub.push_back(double(nodeBounds[n].upper));
		lb.push_back(double(nodeBounds[n].lower));

	}

	nlopt_opt opt;
	//NLOPT_LD_SLSQP <- fails
	//NLOPT_GN_ISRES
	opt = nlopt_create(NLOPT_LN_COBYLA, nNodes); /* algorithm and dimensionality */
	nlopt_set_upper_bounds(opt, &ub[0]);
	nlopt_set_lower_bounds(opt, &lb[0]);

	// minimize SUM(W_p)
	nlopt_set_min_objective(opt,objective_upperBoundTPowerManagerModel, NULL);

	std::vector<time_constraint_data> data;
	std::vector<double> deltas;
	for (int i = 0; i < nNodes; ++i) {
			deltas.push_back(double(
					nodeCells[i] * nodeTperCells[i] * nodeBounds[i].upper));
		}

	data.resize(nNodes);
	for (int i = 0; i < nNodes; ++i) {

		data[i].p = i;
		data[i].time = T;
		data[i].wcap = wcap;
		data[i].delta = &deltas[0];
		/**
		 * results are much more smooth but more violent with equality.
		 */
		// delta_p / W_p -T <= 0
		nlopt_add_inequality_constraint(opt,
				constraint_upperBoundTPowerManagerModel, &data[i], 1e-6);
	}

	//SUM(W_p) - W_cap <= 0
	nlopt_add_inequality_constraint(opt,
			constraint_upperBoundTPowerManagerModel2, &data[0], 1e-6);


	nlopt_set_xtol_rel(opt, 1e-5);

	/* some initial guess */
	std::vector<double> x(nNodes, 0);
	for (int n = 0; n < nNodes; ++n) {
		x[n] = double(nodeBounds[n].upper);
	}

	double minf; /* the minimum objective value, upon return */

	int r = nlopt_optimize(opt, &x[0], &minf);

	for (int n = 0; n < nNodes; ++n) {
		nodeW[n] = Watt(x[n]);
	}

	nlopt_destroy(opt);

	return r;
}


