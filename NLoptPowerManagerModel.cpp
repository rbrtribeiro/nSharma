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
#include "nSharma.H"

#include "Pstream.H"

#include <nlopt/include/nlopt.h>

using namespace Foam;

NLoptPowerManagerModel::NLoptPowerManagerModel(PerformanceModel* pm,
		CommGraph* cg, PowerInterface* pi) :
		performanceModel(pm), commGraph(cg), powerInterface(pi) {

	tpercell_set = false;

}

NLoptPowerManagerModel::~NLoptPowerManagerModel() {

}

const List<float>& NLoptPowerManagerModel::getTperCell() {

	if (!tpercell_set) {

		tpercell = performanceModel->getTperCell();
		tpercell_set = true;
	}

	return tpercell;

}

void NLoptPowerManagerModel::assess(bool enabled) {

	DLM_PRINT_HEADER_START(Pstream::master())

	int nNodes = commGraph->getNumberOfNodes();

	std::vector<Watt> nodeW(nNodes, 0);
	label rankCells;
	std::vector<int> nCells(Pstream::nProcs(), 0);

	WattPair bounds; //Watt, force match below MPI call
	powerInterface->getHWWattLimits(bounds);

	std::vector<WattPair> ranksBounds;
	ranksBounds.resize(Pstream::nProcs());

	MPI_Gather(&bounds, 2, MPI_WATT, &ranksBounds[0], 2, MPI_WATT, 0,
			MPI_COMM_WORLD);

	if (enabled) {

		rankCells = procedureProfPool::getMesh().nCells();
		MPI_Gather(&rankCells, 1, MPI_WATT, &nCells[0], 1, MPI_WATT, 0,
				MPI_COMM_WORLD);

	}

	if (!tpercell_set && nSharma::getInstance()->PARAMS().enablePowerOptimize
			&& nSharma::getInstance()->iterationCount() > 1)
		performanceModel->update();

	if (enabled) {

		const List<float>& tpercell = getTperCell();
	}

	if (Pstream::master()) {

		std::vector<WattPair> nodeBounds;
		nodeBounds.resize(nNodes);

		for (int n = 0; n < nNodes; ++n) {
			nodeBounds[n] = ranksBounds[commGraph->getMasterRankOfNode(n)];
		}

		float cap = powerInterface->cap;
		Watt availW = 0;
		for (int n = 0; n < nNodes; ++n) {
			availW += nodeBounds[n].upper;
		}

		Watt wcap = round(availW * cap);

		for (int n = 0; n < nNodes; ++n) {

			nodeW[n] = round(nodeBounds[n].upper * cap);

		}

		if (enabled) {

			std::vector<float> predictedTime(Pstream::nProcs(), 0);
                        float sum =0;
                        float mean;
                        float max =0;
                        for (int i = 0; i < Pstream::nProcs(); ++i) {
                                sum += tpercell[i];
                                if (tpercell[i] > max) max = tpercell[i];
                        }

                        mean = sum/Pstream::nProcs();
                        //std::cout << "Using tpercell mean: " << mean << std::endl;
                        //std::cout << "Using tpercell max: " << max << std::endl;

			for (int i = 0; i < Pstream::nProcs(); ++i) {
				predictedTime[i] = nCells[i] * tpercell[i];
				//predictedTime[i] = nCells[i] * max;

			}

			std::vector<int> maxRankPerNode(nNodes, 0);

			float maxT = 0;
			for (int n = 0; n < nNodes; ++n) {
				std::vector<int> nodeRanks = commGraph->getProcsInNode(n);

				float maxExec = 0;
				int maxExecRank = 0;
				for (int nr = 0; nr < nodeRanks.size(); ++nr) {
					if (predictedTime[nodeRanks[nr]] > maxExec) {
						maxExec = predictedTime[nodeRanks[nr]];
						maxExecRank = nodeRanks[nr];
					}
				}

				maxRankPerNode[n] = maxExecRank;

				if (maxExec > maxT)
					maxT = maxExec;

			}

			std::vector<float> nodeTperCell(nNodes, 0);

			for (int n = 0; n < nNodes; ++n) {
				nodeTperCell[n] = tpercell[maxRankPerNode[n]];
				//nodeTperCell[n] = max;
			}

			std::vector<int> nodeCells(nNodes, 0);

			for (int n = 0; n < nNodes; ++n) {
				nodeCells[n] = nCells[maxRankPerNode[n]];
			}

			std::cout << "calling nlopt with arguments: " << std::endl;
			std::cout << "nodeBounds: ";
			for (int n = 0; n < nNodes; ++n) {
				std::cout << nodeBounds[n].lower << "," << nodeBounds[n].upper
						<< " ";

			}

			std::cout << std::endl << "nodeCells: ";
			for (int n = 0; n < nNodes; ++n) {
				std::cout << nodeCells[n] << " ";

			}
			std::cout << std::endl << "nodeTperCell: ";
			for (int n = 0; n < nNodes; ++n) {
				std::cout << nodeTperCell[n] << " ";

			}

			std::cout << std::endl << "maxT: " << maxT << std::endl;
			std::cout << "wcap: " << wcap << std::endl;

			int r = solveNlopt(nodeBounds, nodeCells, nodeTperCell, nodeW,
					&maxT, &wcap);

			std::cout << "solveNlopt return: " << r << std::endl;

		}

		std::cout << "powerManager result: W: ";
		Watt sum = 0;
		for (int n = 0; n < nNodes; ++n) {
			std::cout << nodeW[n] << " ";
			sum += nodeW[n];

		}
		std::cout << ", Sum: " << sum;
		std::cout << ", Wcap: " << wcap;

		std::cout << std::endl << std::endl;

	}

	if (enabled) {

		MPI_Bcast(&nodeW[0], nNodes, MPI_WATT, 0, MPI_COMM_WORLD);

		Watt myNodeWatt = nodeW[commGraph->getMyNodeID()];

		powerInterface->setWatts(myNodeWatt);
	}

	DLM_PRINT_HEADER_END(Pstream::master())

}
