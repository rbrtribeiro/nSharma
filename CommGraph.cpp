/*
 * CommGraph.cpp
 *
 *  Created on: Jun 2, 2017
 *      Author: rr
 */

#include "CommGraph.h"

#include "mpi.h"

namespace Foam {

CommGraph::CommGraph() {

	nodeIDcounter = 0;

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	//int half = world_size/2;
	//if (rank < half ) processor_name[0] = 'd';

//	if (rank % 2 == 0)
//		processor_name[0] = 'd';

	myNodeName = std::string(processor_name);

	if (Pstream::master()) {

		int nodeID = registerNode(std::string(processor_name));
		commGraph[nodeID].push_back(rank);

		for (int slave = Pstream::firstSlave(); slave <= Pstream::lastSlave();
				slave++) {

			char slave_processor_name[MPI_MAX_PROCESSOR_NAME];
			MPI_Recv(slave_processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
					slave, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int nodeID = registerNode(std::string(slave_processor_name));
			commGraph[nodeID].push_back(slave);

		}

	} else {

		MPI_Send(processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, 0,
				MPI_COMM_WORLD);
//
	}

	distributeGraph();

	for (cGraph_iterator it = commGraph.begin(); it != commGraph.end(); it++) {
		std::vector<int>& procs = it->second;
		for (std::vector<int>::iterator itProcs = procs.begin();
				itProcs != procs.end(); itProcs++) {
			rankToNodeID.insert(std::pair<int, int>(*itProcs, it->first));
		}
	}

	myNodeID = getNodeID(rank);

}

CommGraph::~CommGraph() {
	// TODO Auto-generated destructor stub
}

int CommGraph::getNumberOfNodes() {
	return commGraph.size();
}

int CommGraph::getNumberOfProcsInNode(int nodeID) {
	if (commGraph.find(nodeID) != commGraph.end()) {
		return commGraph[nodeID].size();
	} else
		return -1;
}

std::vector<int>& CommGraph::getProcsInNode(int nodeID) {
	if (commGraph.find(nodeID) != commGraph.end()) {
		return commGraph[nodeID];
	} else {
		Pout << "Error in getProcsInNode" << endl;
	}
}

int CommGraph::getRankOfNodeIndex(int nodeID, int index) {
	if (commGraph.find(nodeID) != commGraph.end()) {
		return commGraph[nodeID][index];
	} else {
		return -1;
	}
}

int CommGraph::getNodeID(int rank) {
	if (rankToNodeID.find(rank) != rankToNodeID.end()) {
		return rankToNodeID[rank];
	} else {
		return -1;
	}
}

int CommGraph::registerNode(std::string nodeName) {

	NodeToID_iterator it = nodeToID.find(nodeName);
	if (it == nodeToID.end()) {
		nodeToID.insert(std::pair<std::string, int>(nodeName, nodeIDcounter));
		commGraph.insert(
				std::pair<int, std::vector<int> >(nodeIDcounter,
						std::vector<int>()));
		return nodeIDcounter++;
	} else {
		return nodeToID[nodeName];
	}
}

void CommGraph::distributeGraph() {

	int nClusters;
	if (Pstream::master())
		nClusters = getNumberOfNodes();

	MPI_Bcast(&nClusters, 1, MPI_INT, Pstream::masterNo(), MPI_COMM_WORLD);

	for (int i = 0; i < nClusters; i++) {

		std::vector<int> value;
		int size;

		if (Pstream::master()) {
			value = commGraph[i];
			size = value.size();
		}
		MPI_Bcast(&size, 1, MPI_INT, Pstream::masterNo(), MPI_COMM_WORLD);

		if (!Pstream::master())
			value.resize(size);

		MPI_Bcast(&(value[0]), size, MPI_INT, Pstream::masterNo(),
				MPI_COMM_WORLD);
		if (!Pstream::master()) {
			commGraph.insert(
					std::pair<int, std::vector<int> >(i,
							std::vector<int>(value.begin(), value.end())));
		}

	}
}
bool CommGraph::isInDifferentNode(int rank, int neigh, int& neighNodeID) {
	if (getNodeID(rank) != getNodeID(neigh)) {
		neighNodeID = getNodeID(neigh);
		return true;
	} else
		return false;

}

int CommGraph::getMasterRankOfNode(int nodeID) {

	return commGraph[nodeID][0];

}

std::string CommGraph::getMyNodeName() {

	return myNodeName;

}

int CommGraph::getMyNodeID() {

	return myNodeID;

}

void CommGraph::print() {

	DLM_PRINT_HEADER_START(Pstream::master())

	Info << "Communication graph:" << endl;

	for (NodeToID_iterator it = nodeToID.begin(); it != nodeToID.end(); it++) {
		Info << "Node name: " << it->first << " ID: " << it->second << endl;
	}

	for (int i = 0; i < Pstream::nProcs(); i++) {
		if (Pstream::myProcNo() == i) {
			for (cGraph_iterator it = commGraph.begin(); it != commGraph.end();
					it++) {
				Pout << "NodeID: " << it->first << " Procs: ";
				std::vector<int>& procs = it->second;
				for (std::vector<int>::iterator itProcs = procs.begin();
						itProcs != procs.end(); itProcs++) {
					Pout << *itProcs << " ";
				}
				Pout << endl;
			}

		}
		MPI_Barrier (MPI_COMM_WORLD);
	}

	DLM_PRINT_HEADER_END(Pstream::master())
}

} //namespace foam
