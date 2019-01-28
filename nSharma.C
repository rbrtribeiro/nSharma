/*---------------------------------------------------------------------------*\
  =========                 |
 \\      /  F ield         | foam-extend: Open Source CFD
 \\    /   O peration     |
 \\  /    A nd           | For copyright notice see file Copyright
 \\/     M anipulation  |
 -------------------------------------------------------------------------------
 License
 This file is part of foam-extend.

 foam-extend is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation, either version 3 of the License, or (at your
 option) any later version.

 foam-extend is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

 \*---------------------------------------------------------------------------*/

#include "nSharma.H"
#include "CostLinearLoadManagerModel.H"
#include "procedureProfPool.H"
#include "MPIfiedProcedure.H"

#include "Pstream.H"
#include "UPstream.H"
#include "OPstream.H"
#include "dictionary.H"
#include "processorFvPatch.H"
#include "scalar.H"
#include "decompositionMethod.H"
#include "PstreamReduceOps.H"
#include "mapDistributePolyMesh.H"
#include "IOobjectList.H"
#include "dynamicRefineExtFvMesh.H"

#include <stddef.h>
#include <iomanip>
#include <sstream>
#include <map>
#include <cassert>
#include <numeric>
#include <vector>
#include <cfloat>
#include <cmath>
#include <list>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#include "mpi.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::nSharma* Foam::nSharma::instance_(NULL);

MPI_Datatype Foam::nSharma::MPI_LOADPROCEDURE;

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * * //

Foam::nSharma::nSharma(dictionary& decompositionDict, Time& t, fvMesh& m,
		nSharmaParameters p) :
		decompositionDict_(decompositionDict) {

	PARAMS_ = p;

	if (Pstream::master()) {
		worldEpisodeComputedCellCount.resize(Pstream::nProcs());
		worldCurrentCellCount.resize(Pstream::nProcs());
		busyPercentage_.resize(Pstream::nProcs());
		busyTime_.resize(Pstream::nProcs());
		idlePercentage_.resize(Pstream::nProcs());
		idleTime_.resize(Pstream::nProcs());
		totalIdleTime_.resize(Pstream::nProcs());
		totalBusyTime_.resize(Pstream::nProcs());

	}

	mesh_ = &m;
	time = &t;

	localEpisodeComputedCellCount = 0;

	IOobject io(dynamicFvMesh::defaultRegion,
			procedureProfPool::getRunTime().timeName(),
			procedureProfPool::getRunTime(), IOobject::MUST_READ);

	IOdictionary dict(
			IOobject("dynamicMeshDict", io.time().constant(),
					(io.name() == polyMesh::defaultRegion ? "" : io.name()),
					io.db(), IOobject::MUST_READ, IOobject::NO_WRITE, false));

	meshDict_ = dict;

	window_busyRSD.resize(WINDOW_BUSY_RSD_SIZE, 0.0);
	nSysBfaces_ = 0;
	nSysBfacesNodes_ = 0;

	for (int i = WINDOW_BUSY_RSD_SIZE; i > 0; i--)
		window_busyRSD_x.push_back(i);

	if (Pstream::parRun()) {

		commGraph = new CommGraph();

		commGraph->print();

		MPI_Barrier (MPI_COMM_WORLD);

		forceNextBalance = PARAMS().useCommGraphPartitioning;

		if (PARAMS().useCommGraphPartitioning)
			reconfigureTreeComm();

		printTreeComm();

		MPI_Barrier(MPI_COMM_WORLD);

	}

	iterationCount_ = 1;
	balanceEpisodeID_ = 0;
	lastBE = 0;
	currentIterationBalanced = false;

	IOdictionary nSharmaDict(
			IOobject("nSharmaDict", t.system(), m, IOobject::MUST_READ,
					IOobject::NO_WRITE));

	if (Pstream::parRun()) {

		performanceModel = new TperCellPerformanceModel(PARAMS().window);

		if (PARAMS().enableLoadBalance) {

			if (nSharmaDict.lookupOrAddDefault("includeMigrationCost", false)) {
				loadManager = new CostLinearLoadManagerModel(PARAMS(), dict,
						performanceModel, commGraph);
			} else
				loadManager = new LinearLoadManagerModel(PARAMS(), dict,
						performanceModel, commGraph);
		}

		std::string nodeName = commGraph->getMyNodeName();

		if ((nodeName.find("compute-6") != std::string::npos)
				|| (nodeName.find("compute-4") != std::string::npos)
				|| (nodeName.find("erform") != std::string::npos)
				|| (nodeName.find("002-1") != std::string::npos)) {
			powerInterface = new SeARCHPowerInterface(
					commGraph->getNodeID(Pstream::myProcNo()),
					commGraph->getMyNodeName(),
					std::string("/sys/devices/system/cpu"),
					commGraph->getMasterRankOfNode(commGraph->getMyNodeID())
							== Pstream::myProcNo(), PARAMS().powerCap);
		} else {

			powerInterface = new VoidPowerInterface(
					commGraph->getNodeID(Pstream::myProcNo()),
					commGraph->getMyNodeName(),
					commGraph->getMasterRankOfNode(commGraph->getMyNodeID())
							== Pstream::myProcNo(), PARAMS().powerCap);
		}

		powerManager = new ScalarizationUpperTMinT(performanceModel, commGraph,
				powerInterface);

	}

}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //

Foam::nSharma::~nSharma() {

	if (Pstream::parRun()) {

		delete commGraph;
		delete powerInterface;

	}

}

// * * * * * * * * * * * * * * * API Functions  * * * * * * * * * * * * * //

void Foam::nSharma::Init(Time& t, fvMesh& m) {

	nSharmaParameters& P = *(new nSharmaParameters());

	IOdictionary nSharmaDict(
			IOobject("nSharmaDict", t.system(), m, IOobject::MUST_READ,
					IOobject::NO_WRITE));

	IOdictionary controlDict(
			IOobject("controlDict", t.system(), m, IOobject::MUST_READ,
					IOobject::NO_WRITE));

	P.simulation_iterations = (controlDict.lookupOrAddDefault("endTime", 1.0)
			- controlDict.lookupOrAddDefault("startTime", 0.0))
			/ controlDict.lookupOrAddDefault("deltaT", 0.1);

	P.enableLoadBalance = nSharmaDict.lookupOrAddDefault("enableLoadBalance",
			false);

	P.enablePowerOptimize = nSharmaDict.lookupOrAddDefault(
			"enablePowerOptimize", false);

	if (P.enableLoadBalance && P.enablePowerOptimize) {
		std::cout << "nSharma FATAL: both load balance and power"
				" optimization enabled" << std::endl;

		FatalErrorIn("Foam::nSharma::Inits") << exit(FatalError);
	}

	string TWS = nSharmaDict.lookupOrAddDefault("targetWorkloadSection",
			string(PROF_INFO_MAIN_SECTION_NAME));

	strcpy(P.targetWorkloadSection, TWS.c_str());

	P.balancePeriod = nSharmaDict.lookupOrAddDefault("balancePeriod", 6);

	P.default_balance_period = P.balancePeriod;

	P.minPercentageCellsMoved = nSharmaDict.lookupOrAddDefault(
			"minPercentageCellsMoved", 0.025);

	P.powerCap = nSharmaDict.lookupOrAddDefault("powerCap", 1.0f);

	P.minimal_gain = nSharmaDict.lookupOrAddDefault("minimalGain", 0.05);

	P.window = P.default_balance_period * 3;
//if (!P.enablePowerOptimize)
//	P.window = nSharmaDict.lookupOrAddDefault("window", 5);

	P.useCommGraphPartitioning = nSharmaDict.lookupOrAddDefault(
			"useCommGraphPartitioning", false);

	procedureProfPool::initProfiling(t, m, P.targetWorkloadSection);

	if (Pstream::parRun()) {

		const int nitems = 5;
		int blocklengths[nitems] = { 1, 1, 1, 1, 1 };
		MPI_Datatype types[nitems] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,
				MPI_INT };
		MPI_Aint offsets[nitems];

		offsets[0] = offsetof(MPIfiedProcedure, calls_);
		offsets[1] = offsetof(MPIfiedProcedure, totalTime_);
		offsets[2] = offsetof(MPIfiedProcedure, childTime_);
		offsets[3] = offsetof(MPIfiedProcedure, id_);
		offsets[4] = offsetof(MPIfiedProcedure, opType_);

		MPI_Type_create_struct(nitems, blocklengths, offsets, types,
				&nSharma::MPI_LOADPROCEDURE);
		MPI_Type_commit(&nSharma::MPI_LOADPROCEDURE);

		int totalCells = procedureProfPool::thePool_->mesh_.nCells();

		int totalCells__;
		MPI_Allreduce(&totalCells, &totalCells__, 1, MPI_INT, MPI_SUM,
				MPI_COMM_WORLD);

		totalCells = totalCells__;

		P.minPercentageCellsMoved *= (totalCells / Pstream::nProcs());

	}

	if (!instance_) {

		instance_ = new nSharma(nSharmaDict, t, m, P);
		DLM_PRINT_HEADER_START(Pstream::master())
		Info << nSharmaDict;
		DLM_PRINT_HEADER_END(Pstream::master())
	}

}

void Foam::nSharma::Finalize() {

	Info << "Finishing load manager" << endl;
	Foam::nSharma::getInstance()->updateProfiling();

	Foam::nSharma::getInstance()->updateBusyIdleTime();

	Foam::nSharma::getInstance()->writeLoadData();

	delete nSharma::getInstance();

}

bool Foam::nSharma::Iteration() {

	bool balanced = nSharma::getInstance()->balance();

	nSharma::getInstance()->powerOptimize();

	return balanced;

}

void Foam::nSharma::IterationDone() {

	nSharma::getInstance()->iterationDone();

}

// * * * * * * * * * * * * * * * Access Functions  * * * * * * * * * * * * * //

void nSharma::printTreeComm() {
	Pout << Pstream::myProcNo() << " "
			<< Pstream::treeCommunication()[Pstream::myProcNo()].below()
			<< endl;
}

nSharma* Foam::nSharma::getInstance() {
	if (instance_ == NULL) {
		std::cout << "nSharma FATAL: instance_ not defined" << std::endl;

		FatalErrorIn("Foam::nSharma::getInstance") << exit(FatalError);
	}
	return instance_;
}

LoadManagerModel* Foam::nSharma::getLoadManagerInstance() {
	return nSharma::getInstance()->loadManager;
}

nSharma::loadDataType& Foam::nSharma::getLoadData() {
	return nSharma::getInstance()->loadData;
}

const Foam::nSharma_parameters& Foam::nSharma::PARAMS() {
	return PARAMS_;
}

Foam::nSharma_parameters& Foam::nSharma::editPARAMS() {
	return PARAMS_;
}

void Foam::nSharma::enablePool(bool v) {

	procedureProfPool::thePool_->setEnable(v);

}

// Get merging distance when matching face centres
scalar Foam::nSharma::getMergeDistance(const boundBox& bb) {
	scalar mergeTol = defaultMergeTol;

	scalar writeTol = Foam::pow(scalar(10.0),
			-scalar(IOstream::defaultPrecision()));

	Info << "Merge tolerance : " << mergeTol << nl << "Write tolerance : "
			<< writeTol << endl;

	scalar mergeDist = mergeTol * bb.mag();

	Info << "Overall meshes bounding box : " << bb << nl
			<< "Relative tolerance          : " << mergeTol << nl
			<< "Absolute matching distance  : " << mergeDist << nl << endl;

	return mergeDist;
}

void Foam::nSharma::printMeshData(Ostream& os, const polyMesh& mesh) {
	os << "Number of points:           " << mesh.points().size() << nl
			<< "          faces:            " << mesh.faces().size() << nl
			<< "          internal faces:   " << mesh.faceNeighbour().size()
			<< nl << "          cells:            " << mesh.cells().size() << nl
			<< "          boundary patches: " << mesh.boundaryMesh().size()
			<< nl << "          point zones:      " << mesh.pointZones().size()
			<< nl << "          face zones:       " << mesh.faceZones().size()
			<< nl << "          cell zones:       " << mesh.cellZones().size()
			<< nl;
}

MPIfiedProcedure& Foam::nSharma::getLoadDataElement(label opID, label proc) {
	return loadData.at(opID).at(proc);
}

void Foam::nSharma::setLoadDataProcBuffer(MPIfiedProcedure* infos, label proc) {

	label n = procedureProfPool::thePool_->map().size();
	loadData.resize(n);

	for (int i = 0; i < n; i++) {
		//do only once
		if (proc == Pstream::masterNo())
			loadData.at(i).resize(Pstream::nProcs());
		getLoadDataElement(i, proc) = infos[i];

	}

}

/**
 * master method
 * expensive method
 */
void Foam::nSharma::writeLoadData() {

//ADD_SECTION(writeLoadData);

//if (!Pstream::parRun())
//	return;

	if (!Pstream::master())
		return;

	label nProcs = Pstream::nProcs();

	procedureProfPool::mapType& allInfo_ = procedureProfPool::thePool_->map();
	loadDataType& loadData = getLoadData();

	DLM_PRINT_HEADER_START(Pstream::master())

	string lineBreak = "\n";
	int header_offset = 90;
	string line = SSTR(std::left << std::setw(header_offset) << "");
	string preline = SSTR(std::left << "");

	for (label itp = 0; itp < nProcs; itp++) {
		line += SSTR(std::right << std::setw(11) << itp);
	}

//line += SSTR(std::fixed << std::setw(9) << std::setprecision(2) << "RSD(%)");
//line += SSTR(std::fixed << std::setw(9) << std::setprecision(2) << "RMD(%)");
//line += SSTR(std::fixed << std::setw(9) << std::setprecision(2) << "MD(s)");

	line += lineBreak.c_str();

	cout << line;
	line = SSTR(std::left << "");

	std::map < label, string > mymap;

	for (procedureProfPool::mapConstIterator it = allInfo_.begin();
			it != allInfo_.end(); ++it) {

		if (it->second->id0() >= loadData.size())
			continue;

		const procedureProfInfo* info = it->second;
		preline = SSTR(std::left << "");

		while (info->parent().id() != info->id()) {
			info = &(info->parent());
			preline = SSTR(std::left << info->description() << " > ")+ preline;
		}

		//preline = SSTR(std::left << std::setw(3) << it->second->id() << " > ")+preline;

		line += SSTR(
				std::left << std::setw(header_offset) << preline + it->first);

		for (label itp = 0; itp < nProcs; itp++) {

			line +=

					SSTR(
							std::fixed << std::setw(11) << std::setprecision(3) << getLoadDataElement(it->second->id0(), itp).totalTime_);

		}

//		line +=
//				SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RSD_[it->second->id0()]);
//
//		line +=
//				SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RMD_[it->second->id0()]);
//
//		line +=
//				SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->MD_[it->second->id0()]);

		line += lineBreak.c_str();

		mymap[it->second->id()] = line;

		line = SSTR(std::left << "");
	}

	std::map<label, string>::iterator it;
	for (std::map<label, string>::iterator it = mymap.begin();
			it != mymap.end(); ++it)
		std::cout << it->second;

	/***************************************************************/

	preline = SSTR(std::left << "Number_of_cells");

	line += lineBreak.c_str();

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(3) << getCurrentCells(itp));

	}

//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RSD_[allInfo_.size()]);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RMD_[allInfo_.size()]);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->MD_[allInfo_.size()]);

	line += lineBreak.c_str();

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Busy_time");

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(3) << busyTime_[itp]);

	}

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Idle_time");

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(3) << idleTime_[itp]);

	}

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Idle(%)");

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(1) << idlePercentage_[itp]);

	}

	line += lineBreak.c_str();

	cout << line;

	line = SSTR(std::left << "");

	/***************************************************************/

	preline = SSTR(std::left << "Busy(%)");

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(1) << busyPercentage_[itp]);

	}

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Total_Busy");

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(3) << totalBusyTime_[itp]);

	}

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Total_Idle");

	line += SSTR(std::left << std::setw(header_offset) << preline);

	for (label itp = 0; itp < nProcs; itp++) {
		line +=
				SSTR(
						std::fixed << std::setw(11) << std::setprecision(3) << totalIdleTime_[itp]);

	}

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Busy_Relative_STD(%)");

	line += lineBreak.c_str();

	line += SSTR(std::left << std::setw(header_offset) << preline);

//label id_amul = allInfo_.find("op_b:AmulCore")->second->id0();

	line += SSTR(
			std::fixed << std::setw(9) << std::setprecision(2) << busyRSD_);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RMD_[allInfo_.size()]);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->MD_[allInfo_.size()]);

	cout << line;

	line = SSTR(std::left << "");

	/***************************************************************/

	preline = SSTR(std::left << "Idle_Relative_STD(%)");

	line += lineBreak.c_str();

	line += SSTR(std::left << std::setw(header_offset) << preline);

//label id_amul = allInfo_.find("op_b:AmulCore")->second->id0();

	line += SSTR(
			std::fixed << std::setw(9) << std::setprecision(2) << idleRSD_);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RMD_[allInfo_.size()]);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->MD_[allInfo_.size()]);

	line += lineBreak.c_str();

	cout << line;

	line = SSTR(std::left << "");

	/***************************************************************/

	preline = SSTR(std::left << "System_Proc_BFaces");

	line += lineBreak.c_str();

	line += SSTR(std::left << std::setw(header_offset) << preline);

//label id_amul = allInfo_.find("op_b:AmulCore")->second->id0();

	line += SSTR(
			std::fixed << std::setw(9) << std::setprecision(2) << nSysBfaces_);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RMD_[allInfo_.size()]);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->MD_[allInfo_.size()]);

	line += lineBreak.c_str();

	cout << line;

	line = SSTR(std::left << "");

	/***************************************************************/

	preline = SSTR(std::left << "System_Proc_BFaces_Nodes");

	line += lineBreak.c_str();

	line += SSTR(std::left << std::setw(header_offset) << preline);

//label id_amul = allInfo_.find("op_b:AmulCore")->second->id0();

	line +=
			SSTR(
					std::fixed << std::setw(9) << std::setprecision(2) << nSysBfacesNodes_);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->RMD_[allInfo_.size()]);
//
//	line +=
//			SSTR(std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->MD_[allInfo_.size()]);

	line += lineBreak.c_str();

	cout << line;

	line = SSTR(std::left << "");

	/***************************************************************/

	DLM_PRINT_HEADER_END(Pstream::master())

//END_SECTION(writeLoadData);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nSharma::reconfigureTreeComm() {

	int nNodes = commGraph->getNumberOfNodes();

	std::vector<int> rankMapper;
	for (label node = 0; node < nNodes; node++) {
		std::vector<int>& ranks = commGraph->getProcsInNode(node);
		for (int rank = 0; rank < ranks.size(); ++rank) {
			rankMapper.push_back(ranks[rank]);
		}

	}

	List < UPstream::commsStruct > t = Pstream::treeCommunication();

	Pout << "before: ";
	printTreeComm();

	for (label procID = 0; procID < Pstream::nProcs(); procID++) {
		UPstream::commsStruct& cs = t[procID];
		cs.commGraphRemap(rankMapper);
		Pstream::treeCommunication()[rankMapper[procID]] = cs;
	}

}

void Foam::nSharma::resetAllOperationTimes() {

	DLM_PRINT_HEADER_START(Pstream::master())

	procedureProfPool::mapType& allInfo_ = procedureProfPool::thePool_->map();

//Reset all ops from each process pool
	for (procedureProfPool::mapConstIterator it = allInfo_.begin();
			it != allInfo_.end(); ++it) {

		if (it->second->isOperation()) {
			it->second->totalTime_ = 0.0f;
			it->second->parent_.totalTime_ = 0.0f;
		}
	}

//reset target section. in the case of on-the-fly section it
// will set a new timer
	procedureProfPool::thePool_->getTargetWorkloadSection() =
			procedureProfPool::targetWorkloadSection(
					procedureProfPool::thePool_->getTargetWorkloadSection().name(),
					procedureProfPool::thePool_->getTargetWorkloadSection().section_ptr());

	if (allInfo_.find("local_redistribute") != allInfo_.end())
		allInfo_.find("local_redistribute")->second->totalTime_ = 0.0;

	resetComputedCellCount();

	Info << "Profilling data reset\n" << endl;
	DLM_PRINT_HEADER_END(Pstream::master())
}

void Foam::nSharma::updateRSD() {
	updateRSD_MIN_MAX();
//updateRSD_busyTime();
}

void Foam::nSharma::updateRSD_MIN_MAX() {

	label nProcs = Pstream::nProcs();

	double tmp = 0.0;
	double mean = 0.0;
	double min = FLT_MAX;
	double max = 0.0;

	for (label itp = 0; itp < nProcs; itp++) {
		if (busyTime_[itp] > max)
			max = busyTime_[itp];
		if (busyTime_[itp] < min)
			min = busyTime_[itp];
	}

	tmp += max;
	tmp += min;

	mean = tmp / static_cast<double>(2);

	if (mean == 0)
		return;

	double sumsqr = 0;
	double variance, stddeviation;

	sumsqr += powf(max - mean, 2);
	sumsqr += powf(min - mean, 2);

	variance = sumsqr;
	stddeviation = sqrtf(variance);

	busyRSD_ = (stddeviation / mean) * 100.f;

	{
		double tmp = 0.0;
		double mean = 0.0;
		double min = FLT_MAX;
		double max = 0.0;

		for (label itp = 0; itp < nProcs; itp++) {
			if (idleTime_[itp] > max)
				max = idleTime_[itp];
			if (idleTime_[itp] < min)
				min = idleTime_[itp];
		}

		tmp += max;
		tmp += min;

		mean = tmp / static_cast<double>(2);

		if (mean == 0)
			return;

		double sumsqr = 0;
		double variance, stddeviation;

		sumsqr += powf(max - mean, 2);
		sumsqr += powf(min - mean, 2);

		variance = sumsqr;
		stddeviation = sqrtf(variance);

		idleRSD_ = (stddeviation / mean) * 100.f;
	}

}

void Foam::nSharma::updateComputedCellCount() {

	localEpisodeComputedCellCount +=
			procedureProfPool::thePool_->mesh_.nCells();

}

// totalTime_ only
void Foam::nSharma::updateRSD_busyTime() {

	label nProcs = Pstream::nProcs();

	double tmp[2] = { 0, 0 };
	double mean[2];
	double max[2] = { 0, 0 };

	for (label itp = 0; itp < nProcs; itp++) {

		tmp[0] += busyTime_[itp];
		tmp[1] += idleTime_[itp];
	}

	mean[0] = tmp[0] / static_cast<double>(nProcs);
	mean[1] = tmp[1] / static_cast<double>(nProcs);

	if (mean[0] == 0 || mean[1] == 0)
		return;

	double sumsqr[2] = { 0, 0 };
	double variance[2], stddeviation[2];

	for (label itp = 0; itp < nProcs; itp++) {

		tmp[0] = busyTime_[itp];
		tmp[1] = idleTime_[itp];

		sumsqr[0] += powf(tmp[0] - mean[0], 2);
		sumsqr[1] += powf(tmp[1] - mean[1], 2);

		tmp[0] = fabs(tmp[0] - mean[0]);
		tmp[1] = fabs(tmp[1] - mean[1]);

		if (tmp[0] > max[0])
			max[0] = tmp[0];

		if (tmp[1] > max[1])
			max[1] = tmp[1];

	}

	variance[0] = sumsqr[0] / static_cast<double>(nProcs - 1);
	stddeviation[0] = sqrtf(variance[0]);

	variance[1] = sumsqr[1] / static_cast<double>(nProcs - 1);
	stddeviation[1] = sqrtf(variance[1]);

	busyRSD_ = (stddeviation[0] / mean[0]) * 100.f;
	idleRSD_ = (stddeviation[1] / mean[1]) * 100.f;

}

label Foam::nSharma::boilProfilingPool(MPIfiedProcedure*& infos) {

	procedureProfPool * pool = procedureProfPool::thePool_;

	procedureProfPool::mapType& allInfo_ = pool->map();
	label n = allInfo_.size();

	if (infos == NULL)
		infos = new MPIfiedProcedure[n];

//get all info and closed procedures
	for (procedureProfPool::mapConstIterator it = allInfo_.begin();
			it != allInfo_.end(); ++it) {

		it->second->MPIfy(infos[it->second->id0()]);

	}

//add those on-the-fly
	procedureProfStack& stack = pool->stack();
	procedureProfStack::const_iterator it = stack.begin();
	scalar oldElapsed = 0;

	do {
		//reference to info object which is already in the info map
		const procedureProfInfo &info = *(*it);

//		if (info.description() == "updateLM") {
//			//clockTime* t = new clockTime();
//			//stack.timers_[info.id()]=t;
//			Pout << info.parent().description() << endl;
//				}

		scalar elapsed = stack.timers_[info.id()]->elapsedTime();

		infos[info.id0()].calls_ += 1;
		infos[info.id0()].totalTime_ += elapsed;
		infos[info.id0()].childTime_ += oldElapsed;

		oldElapsed = elapsed;
		++it;
	} while (it != stack.end());

	return n;
}

void Foam::nSharma::updateBusyIdleTime() {

	ADD_SECTION(updateBusyIdleTime);

	DLM_PRINT_HEADER_START(Pstream::master())

	procedureProfPool * pool = procedureProfPool::thePool_;

	if (!Pstream::parRun())
		return;

	if (Pstream::master()) {

		scalar busy_idle[6];

		pool->getTargetWorkloadSection().getPercentages(busy_idle);

		busyPercentage_[0] = busy_idle[0];
		idlePercentage_[0] = busy_idle[1];
		busyTime_[0] = busy_idle[2];
		idleTime_[0] = busy_idle[3];
		totalBusyTime_[0] = pool->totalBusyTime_;
		totalIdleTime_[0] = pool->totalIdleTime_;

		for (int slave = Pstream::firstSlave(); slave <= Pstream::lastSlave();
				slave++) {

			MPI_Recv(busy_idle, 6, MPI_DOUBLE, /*Pstream::procID(slave)*/slave,
					0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			busyPercentage_[slave] = busy_idle[0];
			idlePercentage_[slave] = busy_idle[1];
			busyTime_[slave] = busy_idle[2];
			idleTime_[slave] = busy_idle[3];
			totalBusyTime_[slave] = busy_idle[4];
			totalIdleTime_[slave] = busy_idle[5];

		}

		updateRSD();

		window_busyRSD.push_back(busyRSD_);

	} else {

		scalar busy_idle[6];

		pool->getTargetWorkloadSection().getPercentages(busy_idle);

		busy_idle[4] = pool->totalBusyTime_;
		busy_idle[5] = pool->totalIdleTime_;

		MPI_Send(busy_idle, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

	}

	Info << "updateBusyIdleTime data updated\n" << endl;
	DLM_PRINT_HEADER_END(Pstream::master())
	END_SECTION(updateBusyIdleTime);

}

void Foam::nSharma::updateProfiling() {

	ADD_SECTION(updateProfiling);

	DLM_PRINT_HEADER_START(Pstream::master())

	long episodeCells = getComputedCellCount();
	long currentCells = procedureProfPool::thePool_->mesh_.nCells();

	if (Pstream::master()) {

		std::vector<long> episodeCellsV;
		episodeCellsV.resize(Pstream::nProcs());

		std::vector<long> currentCellsV;
		currentCellsV.resize(Pstream::nProcs());

		MPIfiedProcedure* infos = NULL;
		label nOps = boilProfilingPool(infos);

		setLoadDataProcBuffer(infos, Pstream::masterNo());

		episodeCellsV[0] = episodeCells;
		currentCellsV[0] = currentCells;

		if (Pstream::parRun()) {

			for (int slave = Pstream::firstSlave();
					slave <= Pstream::lastSlave(); slave++) {

				//MPIfiedProcedure* infos = managerInstance->getLoadDataProcBuffer(
				//		Pstream::procID(slave));

				MPI_Recv(infos, nOps, nSharma::MPI_LOADPROCEDURE,
				/*Pstream::procID(slave)*/slave, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				setLoadDataProcBuffer(infos,
				/*Pstream::procID(slave)*/slave);

				long ncells;
				MPI_Recv(&ncells, 1, MPI_LONG, /*Pstream::procID(slave)*/slave,
						0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				episodeCellsV[slave] = ncells;

				MPI_Recv(&ncells, 1, MPI_LONG, /*Pstream::procID(slave)*/slave,
						0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				currentCellsV[slave] = ncells;

			}

		}

		setPreviousEpisodeComputedCells(episodeCellsV);
		setCurrentCells(currentCellsV);

	} else {

		label n;
		MPIfiedProcedure* infos = NULL;

		n = boilProfilingPool(infos);

		MPI_Send(infos, n, nSharma::MPI_LOADPROCEDURE, 0, 0, MPI_COMM_WORLD);

		MPI_Send(&episodeCells, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&currentCells, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		delete[] infos;

	}

	if (Pstream::parRun()) {

		long nBFaces = 0;
		long nBFacesNodes = 0;
		forAll(procedureProfPool::getMesh().boundary(), patchi)
		{
			const fvPatch& p = procedureProfPool::getMesh().boundary()[patchi];
			if (isType < processorFvPatch > (p) && p.coupled()) {

				const processorFvPatch& procPatch =
						dynamic_cast<const processorFvPatch&>(p);

				nBFaces += procPatch.size();

				int neighNodeID;
				if (nSharma::getInstance()->commGraph->isInDifferentNode(
						procPatch.myProcNo(), procPatch.neighbProcNo(),
						neighNodeID)) {
					nBFacesNodes += procPatch.size();

				}

			}

		}


	MPI_Allreduce(&nBFaces, &(nSysBfaces_), 1, MPI_LONG, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&nBFacesNodes, &(nSysBfacesNodes_), 1, MPI_LONG, MPI_SUM,
			MPI_COMM_WORLD);

	}

//Info << "System boundary cells: " << old_sysBCells << " -> " << sysBCells
//		<< endl;

	Info << "updateProfiling data updated\n" << endl;
	DLM_PRINT_HEADER_END(Pstream::master())
	END_SECTION(updateProfiling);

}

void Foam::nSharma::updateBalanceEpisode() {

	ADD_SECTION(updateBalanceEpisode);

	int oP = PARAMS().balancePeriod;

	DLM_PRINT_HEADER_START(Pstream::master())

	if (Pstream::master()) {

		Info << "managerInstance->busyRSD_: " << busyRSD_ << endl;

		const size_t n = WINDOW_BUSY_RSD_SIZE;

		const float s_x = std::accumulate(window_busyRSD_x.begin(),
				window_busyRSD_x.end(), 0.0);

		const float s_y = std::accumulate(window_busyRSD.rbegin(),
				window_busyRSD.rbegin() + WINDOW_BUSY_RSD_SIZE, 0.0);

		const float s_xx = std::inner_product(window_busyRSD_x.begin(),
				window_busyRSD_x.end(), window_busyRSD_x.begin(), 0.0);

		const float s_xy = std::inner_product(window_busyRSD_x.begin(),
				window_busyRSD_x.end(), window_busyRSD.rbegin(), 0.0);

		float slope = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

		float rate_change = -slope + 50; //shift to positive
		rate_change = rate_change / 100.0 - 0.5; //normalize and back to negative
		rate_change *= PARAMS().default_balance_period * 2; //factor of current period

		float mag = -0.07 * busyRSD_ + 2.5; //empirical y=mx+b transformation

		float change = rate_change + mag;

		editPARAMS().balancePeriod += change;

		//Minimum trigger period associated with updateProfiling to ensure
		//that we never do a balance+reset whitout an updateProfiling before
		int min_max[2] = { int(PARAMS().default_balance_period / 2.0),
				PARAMS().default_balance_period * 6 };

		if (PARAMS().balancePeriod < min_max[0])
			editPARAMS().balancePeriod = min_max[0];
		else if (PARAMS().balancePeriod > min_max[1])
			editPARAMS().balancePeriod = min_max[1];

		Info << "rate_change " << rate_change << " mag " << mag << " change "
				<< change << endl;

		Info << "window_busyRSD";
		for (int i = 0; i < window_busyRSD.size(); i++) {
			Info << " " << window_busyRSD[i];
		}

		Info << endl;

	}

	MPI_Bcast(&(editPARAMS().balancePeriod), 1, MPI_INT, Pstream::masterNo(),
			MPI_COMM_WORLD);

	Info << "Updated balance period from " << oP << " to "
			<< PARAMS().balancePeriod << endl;

	DLM_PRINT_HEADER_END(Pstream::master())

	END_SECTION(updateBalanceEpisode);

}

void Foam::nSharma::redistribute(labelList finalDecomp) {

	fvMesh& mesh = procedureProfPool::thePool_->mesh_;

// TODO
//MPI_Barrier (MPI_COMM_WORLD);

// Global matching tolerance
//const scalar tolDim = getMergeDistance(mesh.globalData().bb());
	scalar tolDim = globalMeshData::matchTol_ * mesh.bounds().mag();

// Mesh distribution engine
	fvMeshDistribute distributor(mesh, tolDim);

// Do actual sending/receiving of mesh
	autoPtr < mapDistributePolyMesh > map = distributor.distribute(finalDecomp);

	const word dynamicFvMeshTypeName(getMeshDict().lookup("dynamicFvMesh"));

// distribute refinement history
	if (dynamicFvMeshTypeName == "dynamicRefineExtFvMesh") {
		dynamicRefineExtFvMesh& extMesh = *(dynamicRefineExtFvMesh*) &mesh;
		extMesh.distributeRefineCells(map);
	}

	correctBoundaries<scalar>();
	correctBoundaries<vector>();
	correctBoundaries<sphericalTensor>();
	correctBoundaries<symmTensor>();
	correctBoundaries<tensor>();

// Print a bit
//Pout << "After distribution mesh:" << endl;
//printMeshData(Pout, mesh);
//Pout << endl;

//nSharma::instance_->resetAllOperationTimes();

}

bool Foam::nSharma::loadBalanceCells() {

	DLM_PRINT_HEADER_START(Pstream::master())

	labelList finalDecomp;

	bool doBalance;

	if (Pstream::master()) {

		doBalance = busyRSD_ > 5;
	}

	doBalance |= forceNextBalance;

	MPI_Bcast(&doBalance, 1, MPI_UNSIGNED_CHAR, Pstream::masterNo(),
			MPI_COMM_WORLD);

	if (!doBalance)
		Info << "\nBalance-episode-log: Busy RSD too low" << endl;
	if (doBalance)
		doBalance = loadManager->Balance(finalDecomp);

	if (doBalance) {

		Info << "\nBalance-episode-" << balanceEpisodeID_ << "-in-iteration-"
				<< iterationCount_ << " Distributing... " << nl << endl;
		ADD_SECTION(redistribute);
		ADD_SECTION (local_redistribute);
		redistribute(finalDecomp);
		END_SECTION(local_redistribute);
		END_SECTION(redistribute);

		loadManager->PostLoadMigrationUpdate();

		Pout << "My cells after distribution: "
				<< procedureProfPool::getMesh().nCells() << endl;

	} else {
		Info << "\nBalance-episode-" << balanceEpisodeID_ << "-in-iteration-"
				<< iterationCount_ << " No_balance_needed." << nl << endl;
	}

	balanceEpisodeID_++;

	DLM_PRINT_HEADER_END(Pstream::master())

	return doBalance;

}

void Foam::nSharma::iterationDone() {

	if (!Pstream::parRun())
		return;

	MPI_Barrier (MPI_COMM_WORLD); // just to make prof times consistent

//One step ahead of Foam::nSharma::loadBalanceCells();
	incIterationCount();

	updateComputedCellCount();

	if (iterationCount_ > 1) {

		ADD_SECTION (updateLoadManager);

		updateBusyIdleTime();

		if (PARAMS().enableLoadBalance)
			updateBalanceEpisode();

		//Minimum trigger period associated with balance period to ensure
		//that we never do a balance+reset without an updateProfiling before
		if (iterationCount_ == 2 //force second iteration
				|| (iterationCount_
						% int(ceil(PARAMS().default_balance_period / 2.0)) == 0)) {

			updateProfiling();
			writeLoadData();

		}

		END_SECTION(updateLoadManager);

	}

	currentIterationBalanced = false;
	procedureProfPool::getMesh().topoChanging(false);

	MPI_Barrier(MPI_COMM_WORLD); // just to make prof times consistent

}

void Foam::nSharma::powerOptimize() {

	if (!Pstream::parRun())
		return;

	ADD_SECTION (optimizePower);

	if (iterationCount_ == 1)
		powerManager->assess(false);

	if (iterationCount_ == 2 //force second iteration
			|| (iterationCount_
					% int(ceil(PARAMS().default_balance_period / 2.0)) == 0)) {

		//only optmize after period to ensure a decent estimate of TperCell
		if (iterationCount_ >= PARAMS().default_balance_period * 3) {
			powerManager->assess(PARAMS().enablePowerOptimize);
		} else
			powerManager->assess(false);

		if (PARAMS().enablePowerOptimize)
			resetAllOperationTimes();
	}

	END_SECTION(optimizePower);

}

bool Foam::nSharma::balance() {

	if (!Pstream::parRun())
		return false;

	MPI_Barrier (MPI_COMM_WORLD); // just to make prof times consistent

	bool balanced = false;

	if (currentIterationBalanced)
		return false;

	if (iterationCount_ > 1)
		if (iterationCount_ == 2 //force second iteration
		|| ((lastBE + PARAMS().balancePeriod) <= iterationCount_)) {

			ADD_SECTION(loadBalanceCells);

			if (PARAMS().enableLoadBalance)
				balanced = loadBalanceCells();

			if (PARAMS().enableLoadBalance)
				resetAllOperationTimes();

			lastBE = iterationCount_;

			forceNextBalance = false;

			END_SECTION(loadBalanceCells);

		}

	currentIterationBalanced = true;

	MPI_Barrier(MPI_COMM_WORLD); // just to make prof times consistent

	return balanced;

}

void Foam::nSharma::moveCells(label n, label from, label to, fvMesh& mesh) {

	List<int> finalDecomp(mesh.nCells(), Pstream::myProcNo());

	if (Pstream::myProcNo() == from) {
		Info << mesh.cellCells()[n] << endl;
		//finalDecomp[n]= to;
	}

	fvMeshDistribute distributor(mesh, defaultMergeTol);

	Pout << "Wanted distribution:" << distributor.countCells(finalDecomp) << nl
			<< endl;

//	Info << "Distributing... " << nl << endl;
//		autoPtr<mapDistributePolyMesh> map = distributor.distribute(finalDecomp);

}

// ************************************************************************* //
