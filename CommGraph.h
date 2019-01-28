/*
 * CommGraph.h
 *
 *  Created on: Jun 2, 2017
 *      Author: rr
 */

#ifndef COMMGRAPH_H_
#define COMMGRAPH_H_

#include <map>
#include <vector>
#include <string>
#include "Pstream.H"
#include "nSharma_common.h"


namespace Foam {

class CommGraph {

private:

	typedef std::map<int, std::vector<int> >::iterator cGraph_iterator;
	typedef std::map<std::string, int>::iterator NodeToID_iterator;

	long nodeIDcounter;
	std::map<int, std::vector<int> > commGraph;
	std::map<int, int> rankToNodeID;
	std::map<std::string, int> nodeToID;
	std::string myNodeName;
	int myNodeID;

	int registerNode(std::string nodeName);
	void distributeGraph();

public:
	CommGraph();
	virtual ~CommGraph();

	int getNumberOfNodes();
	int getNumberOfProcsInNode(int nodeID);
	int getRankOfNodeIndex(int nodeID, int index);
	std::vector<int>& getProcsInNode(int nodeID);
	int getNodeID(int rank);
	bool isInDifferentNode(int rank, int neigh, int& neighNodeID);
	void print();
	int getMasterRankOfNode(int nodeID);
	std::string getMyNodeName();
	int getMyNodeID();


};

}

#endif /* COMMGRAPH_H_ */
