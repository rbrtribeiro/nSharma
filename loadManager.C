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

#include "loadManager.H"
#include "CostLinearPerformanceModel.H"
#include "Pstream.H"
#include "OPstream.H"
#include "procedureProfPool.H"
#include "dictionary.H"
#include "MPIfiedProcedure.H"
#include "loadManagerProf.H"

#include "mpi.h"
#include <stddef.h>
#include <stdio.h>
#include "scalar.H"

#include <iomanip>

#include <sstream>
#include <map>
#include <cassert>
#include <numeric>
#include <vector>
#include <cfloat>
#include <cmath>

#include "decompositionMethod.H"
#include "PstreamReduceOps.H"
//#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "IOobjectList.H"

//#include "topoMapper.H"

#include "dynamicRefineExtFvMesh.H"

//#include "parMetisDecompDynamic.H"
#include <list>
#include <string.h>
#include <time.h>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::loadManager* Foam::loadManager::managerInstance_(NULL);

MPI_Datatype Foam::loadManager::MPI_LOADPROCEDURE;

Foam::loadManagerParameters Foam::loadManager::PARAMS_;
PerformanceModel* Foam::loadManager::performanceModel;

uint Foam::loadManager::iterationCount_ = 1;
int Foam::loadManager::lastBE = 0;
int Foam::loadManager::balanceEpisodeID_ = 0;
bool Foam::loadManager::currentIterationBalanced = false;

// Tolerance (as fraction of the bounding box). Needs to be fairly lax since
// usually meshes get written with limited precision (6 digits)
const scalar Foam::loadManager::defaultMergeTol = 1E-6;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::loadManager::loadManager(dictionary& decompositionDict) :
		decompositionDict_(decompositionDict) {

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

	localEpisodeComputedCellCount=0;

	IOobject io(dynamicFvMesh::defaultRegion,
			procedureProfPool::getRunTime().timeName(),
			procedureProfPool::getRunTime(), IOobject::MUST_READ);

	IOdictionary dict(
			IOobject("dynamicMeshDict", io.time().constant(),
					(io.name() == polyMesh::defaultRegion ? "" : io.name()),
					io.db(), IOobject::MUST_READ, IOobject::NO_WRITE, false));

	meshDict_ = dict;

	window_busyRSD.resize(WINDOW_BUSY_RSD_SIZE, 0.0);

	for (int i = WINDOW_BUSY_RSD_SIZE; i > 0; i--)
		window_busyRSD_x.push_back(i);

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
void Foam::loadManager::finalize() {

	ADD_SECTION (updateLoadManager);

	Info << "Finishing load manager" << endl;
	Foam::loadManager::updateProfiling();

	Foam::loadManager::updateBusyIdleTime();

	Foam::loadManager::writeLoadData();

	END_SECTION(updateLoadManager);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

loadManager* Foam::loadManager::getManagerInstance() {
	return managerInstance_;
}

PerformanceModel* Foam::loadManager::getPerformanceModelInstance() {
	return loadManager::getManagerInstance()->performanceModel;
}

loadManager::loadDataType& Foam::loadManager::getLoadData() {
	return loadManager::getManagerInstance()->loadData;
}

const Foam::loadManager_parameters& Foam::loadManager::PARAMS() {
	return PARAMS_;
}

void Foam::loadManager::init(Time& t, fvMesh& m) {

	if (!Pstream::parRun())
		return;

	IOdictionary loadManagerDict(
			IOobject("loadManagerDict", t.system(), m, IOobject::MUST_READ,
					IOobject::NO_WRITE));

	IOdictionary controlDict(
			IOobject("controlDict", t.system(), m, IOobject::MUST_READ,
					IOobject::NO_WRITE));

	PARAMS_.simulation_iterations = (controlDict.lookupOrAddDefault("endTime",
			1.0) - controlDict.lookupOrAddDefault("startTime", 0.0))
			/ controlDict.lookupOrAddDefault("deltaT", 0.1);

	PARAMS_.enabled = loadManagerDict.lookupOrAddDefault("enabled", true);

	PARAMS_.targetWorkloadSection = loadManagerDict.lookupOrAddDefault(
			"targetWorkloadSection", string(PROF_INFO_MAIN_SECTION_NAME));

	PARAMS_.balancePeriod = loadManagerDict.lookupOrAddDefault("balancePeriod",
			6);

	PARAMS_.default_balance_period = PARAMS_.balancePeriod;

	PARAMS_.minPercentageCellsMoved = loadManagerDict.lookupOrAddDefault(
			"minPercentageCellsMoved", 0.025);

	PARAMS_.minimal_gain = loadManagerDict.lookupOrAddDefault("minimalGain",
			0.05);

	PARAMS_.window = loadManagerDict.lookupOrAddDefault("window", 5);

	procedureProfPool::initProfiling(t, m, PARAMS().targetWorkloadSection);

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
			&loadManager::MPI_LOADPROCEDURE);
	MPI_Type_commit(&loadManager::MPI_LOADPROCEDURE);

	int totalCells = procedureProfPool::thePool_->mesh_.nCells();

	int totalCells__;
	MPI_Allreduce(&totalCells, &totalCells__, 1, MPI_INT, MPI_SUM,
			MPI_COMM_WORLD);

	totalCells = totalCells__;

	PARAMS_.minPercentageCellsMoved *= (totalCells / Pstream::nProcs());

	if (/*Pstream::master() &&*/!managerInstance_) {

		managerInstance_ = new loadManager(loadManagerDict);
		DLM_PRINT_HEADER_START
		Info << loadManagerDict;
		DLM_PRINT_HEADER_END
	}

	if (loadManagerDict.lookupOrAddDefault("includeMigrationCost", false)) {
		performanceModel = new CostLinearPerformanceModel(PARAMS_);
	} else
		performanceModel = new LinearPerformanceModel(PARAMS_);

}

void Foam::loadManager::resetAllOperationTimes() {

	DLM_PRINT_HEADER_START

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

	ResetComputedCellCount();

	Info << "Profilling data reset\n" << endl;
	DLM_PRINT_HEADER_END
}

void Foam::loadManager::updateRSD() {
	updateRSD_MIN_MAX();
	//updateRSD_busyTime();
}

void Foam::loadManager::updateRSD_MIN_MAX() {

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

void Foam::loadManager::UpdateComputedCellCount() {

	loadManager::managerInstance_->localEpisodeComputedCellCount +=
			procedureProfPool::thePool_->mesh_.nCells();

}

// totalTime_ only
void Foam::loadManager::updateRSD_busyTime() {

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

label Foam::loadManager::BoilProfilingPool(MPIfiedProcedure*& infos) {

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

void Foam::loadManager::enablePool(bool v) {

	procedureProfPool::thePool_->setEnable(v);

}

// Get merging distance when matching face centres
scalar Foam::loadManager::getMergeDistance(const boundBox& bb) {
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

void Foam::loadManager::printMeshData(Ostream& os, const polyMesh& mesh) {
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

MPIfiedProcedure& Foam::loadManager::getLoadDataElement(label opID,
		label proc) {
	return loadData.at(opID).at(proc);
}

void Foam::loadManager::setLoadDataProcBuffer(MPIfiedProcedure* infos,
		label proc) {

	label n = procedureProfPool::thePool_->map().size();
	loadData.resize(n);

	for (int i = 0; i < n; i++) {
		//do only once
		if (proc == Pstream::masterNo())
			loadData.at(i).resize(Pstream::nProcs());
		getLoadDataElement(i, proc) = infos[i];

	}

}

void Foam::loadManager::updateBusyIdleTime() {

	ADD_SECTION(updateBusyIdleTime);

	DLM_PRINT_HEADER_START

	procedureProfPool * pool = procedureProfPool::thePool_;

	if (!Pstream::parRun())
		return;

	if (Pstream::master()) {

		loadManager *managerInstance = loadManager::managerInstance_;

		scalar busy_idle[6];

		pool->getTargetWorkloadSection().getPercentages(busy_idle);

		managerInstance_->busyPercentage_[0] = busy_idle[0];
		managerInstance_->idlePercentage_[0] = busy_idle[1];
		managerInstance_->busyTime_[0] = busy_idle[2];
		managerInstance_->idleTime_[0] = busy_idle[3];
		managerInstance_->totalBusyTime_[0] = pool->totalBusyTime_;
		managerInstance_->totalIdleTime_[0] = pool->totalIdleTime_;

		for (int slave = Pstream::firstSlave(); slave <= Pstream::lastSlave();
				slave++) {

			MPI_Recv(busy_idle, 6, MPI_DOUBLE, /*Pstream::procID(slave)*/slave,
					0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			managerInstance_->busyPercentage_[slave] = busy_idle[0];
			managerInstance_->idlePercentage_[slave] = busy_idle[1];
			managerInstance_->busyTime_[slave] = busy_idle[2];
			managerInstance_->idleTime_[slave] = busy_idle[3];
			managerInstance_->totalBusyTime_[slave] = busy_idle[4];
			managerInstance_->totalIdleTime_[slave] = busy_idle[5];

		}

		managerInstance->updateRSD();

		managerInstance->window_busyRSD.push_back(managerInstance->busyRSD_);

	} else {

		scalar busy_idle[6];

		pool->getTargetWorkloadSection().getPercentages(busy_idle);

		busy_idle[4] = pool->totalBusyTime_;
		busy_idle[5] = pool->totalIdleTime_;

		MPI_Send(busy_idle, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

	}

	Info << "updateBusyIdleTime data updated\n" << endl;
	DLM_PRINT_HEADER_END
	END_SECTION(updateBusyIdleTime);

}

void Foam::loadManager::updateProfiling() {

	ADD_SECTION(updateProfiling);

	DLM_PRINT_HEADER_START

	if (!Pstream::parRun())
		return;

	loadManager *managerInstance = loadManager::managerInstance_;

	long episodeCells = managerInstance->GetComputedCellCount();
	long currentCells = procedureProfPool::thePool_->mesh_.nCells();

	if (Pstream::master()) {

		std::vector<long> episodeCellsV;
		episodeCellsV.resize(Pstream::nProcs());

		std::vector<long> currentCellsV;
		currentCellsV.resize(Pstream::nProcs());

		MPIfiedProcedure* infos = NULL;
		label nOps = BoilProfilingPool(infos);

		managerInstance->setLoadDataProcBuffer(infos, Pstream::masterNo());

		episodeCellsV[0] = episodeCells;
		currentCellsV[0] = currentCells;

		for (int slave = Pstream::firstSlave(); slave <= Pstream::lastSlave();
				slave++) {

			//MPIfiedProcedure* infos = managerInstance->getLoadDataProcBuffer(
			//		Pstream::procID(slave));

			MPI_Recv(infos, nOps, loadManager::MPI_LOADPROCEDURE,
			/*Pstream::procID(slave)*/slave, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			managerInstance->setLoadDataProcBuffer(infos,
			/*Pstream::procID(slave)*/slave);

			long ncells;
			MPI_Recv(&ncells, 1, MPI_LONG, /*Pstream::procID(slave)*/slave, 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			episodeCellsV[slave] = ncells;

			MPI_Recv(&ncells, 1, MPI_LONG, /*Pstream::procID(slave)*/slave, 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			currentCellsV[slave] = ncells;

		}

		managerInstance->setPreviousEpisodeComputedCells(episodeCellsV);
		managerInstance->setCurrentCells(currentCellsV);

	} else {

		label n;
		MPIfiedProcedure* infos = NULL;

		n = BoilProfilingPool(infos);

		MPI_Send(infos, n, loadManager::MPI_LOADPROCEDURE, 0, 0,
				MPI_COMM_WORLD);

		MPI_Send(&episodeCells, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&currentCells, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		delete[] infos;

	}

	Info << "updateProfiling data updated\n" << endl;
	DLM_PRINT_HEADER_END
	END_SECTION(updateProfiling);

}

void Foam::loadManager::updateBalanceEpisode() {

	ADD_SECTION(updateBalanceEpisode);

	int oP = PARAMS_.balancePeriod;

	DLM_PRINT_HEADER_START

	if (Pstream::master()) {

		loadManager *managerInstance = loadManager::managerInstance_;
		Info << "managerInstance->busyRSD_: " << managerInstance->busyRSD_
				<< endl;

		const size_t n = WINDOW_BUSY_RSD_SIZE;

		const float s_x = std::accumulate(
				managerInstance->window_busyRSD_x.begin(),
				managerInstance->window_busyRSD_x.end(), 0.0);

		const float s_y = std::accumulate(
				managerInstance->window_busyRSD.rbegin(),
				managerInstance->window_busyRSD.rbegin() + WINDOW_BUSY_RSD_SIZE,
				0.0);

		const float s_xx = std::inner_product(
				managerInstance->window_busyRSD_x.begin(),
				managerInstance->window_busyRSD_x.end(),
				managerInstance->window_busyRSD_x.begin(), 0.0);

		const float s_xy = std::inner_product(
				managerInstance->window_busyRSD_x.begin(),
				managerInstance->window_busyRSD_x.end(),
				managerInstance->window_busyRSD.rbegin(), 0.0);

		float slope = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

		float rate_change = -slope + 50; //shift to positive
		rate_change = rate_change / 100.0 - 0.5; //normalize and back to negative
		rate_change *= PARAMS_.default_balance_period * 2; //factor of current period

		float mag = -0.07 * managerInstance->busyRSD_ + 2.5; //empirical y=mx+b transformation

		float change = rate_change + mag;

		PARAMS_.balancePeriod += change;

		//Minimum trigger period associated with updateProfiling to ensure
		//that we never do a balance+reset whitout an updateProfiling before
		int min_max[2] = { int(PARAMS_.default_balance_period / 2.0),
				PARAMS_.default_balance_period * 6 };

		if (PARAMS_.balancePeriod < min_max[0])
			PARAMS_.balancePeriod = min_max[0];
		else if (PARAMS_.balancePeriod > min_max[1])
			PARAMS_.balancePeriod = min_max[1];

		Info << "rate_change " << rate_change << " mag " << mag << " change "
				<< change << endl;

		Info << "window_busyRSD";
		for (int i = 0; i < managerInstance->window_busyRSD.size(); i++) {
			Info << " " << managerInstance->window_busyRSD[i];
		}

		Info << endl;

	}

	MPI_Bcast(&(PARAMS_.balancePeriod), 1, MPI_INT, Pstream::masterNo(),
			MPI_COMM_WORLD);

	Info << "Updated balance period from " << oP << " to "
			<< PARAMS_.balancePeriod << endl;

	DLM_PRINT_HEADER_END

	END_SECTION(updateBalanceEpisode);

}

void Foam::loadManager::redistribute(labelList finalDecomp) {

	fvMesh& mesh = procedureProfPool::thePool_->mesh_;

	// TODO
	//MPI_Barrier (MPI_COMM_WORLD);

	// Global matching tolerance
	const scalar tolDim = getMergeDistance(mesh.globalData().bb());
	//scalar tolDim = globalMeshData::matchTol_ * mesh.bounds().mag();

	// Mesh distribution engine
	fvMeshDistribute distributor(mesh, tolDim);

	// Do actual sending/receiving of mesh
	autoPtr < mapDistributePolyMesh > map = distributor.distribute(finalDecomp);

	const word dynamicFvMeshTypeName(
			loadManager::getManagerInstance()->getMeshDict().lookup(
					"dynamicFvMesh"));

	// distribute refinement history
	if (dynamicFvMeshTypeName == "dynamicRefineExtFvMesh") {
		dynamicRefineExtFvMesh& extMesh = *(dynamicRefineExtFvMesh*) &mesh;
		extMesh.distributeRefineCells(map);
	}

	// Print a bit
	//Pout << "After distribution mesh:" << endl;
	//printMeshData(Pout, mesh);
	//Pout << endl;

	//loadManager::managerInstance_->resetAllOperationTimes();

}

bool Foam::loadManager::loadBalanceCells() {

	DLM_PRINT_HEADER_START

	labelList finalDecomp;

	bool doBalance;

	if (Pstream::master()) {
		loadManager *managerInstance = loadManager::managerInstance_;
		doBalance = managerInstance->busyRSD_ > 5;
	}

	MPI_Bcast(&doBalance, 1, MPI_UNSIGNED_CHAR, Pstream::masterNo(),
			MPI_COMM_WORLD);

	if (!doBalance)
		Info << "\nBalance-episode-log: Busy RSD too low" << endl;
	if (doBalance)
		doBalance = performanceModel->Balance(finalDecomp);

	if (doBalance) {

		int old_nBcells = 0;
		forAll(procedureProfPool::getMesh().boundary(), patchi)
		{
			old_nBcells +=
					procedureProfPool::getMesh().boundary()[patchi].size();
		}

		Info << "\nBalance-episode-" << balanceEpisodeID_ << "-in-iteration-"
				<< loadManager::iterationCount_ << " Distributing... " << nl
				<< endl;
		ADD_SECTION(redistribute);
		ADD_SECTION (local_redistribute);
		redistribute(finalDecomp);
		END_SECTION(local_redistribute);
		END_SECTION(redistribute);

		performanceModel->PostLoadMigrationUpdate();

		int nBcells = 0;
		forAll(procedureProfPool::getMesh().boundary(), patchi)
		{
			nBcells += procedureProfPool::getMesh().boundary()[patchi].size();
		}

		int sysBCells;
		MPI_Allreduce(&nBcells, &sysBCells, 1, MPI_INT, MPI_SUM,
				MPI_COMM_WORLD);

		int old_sysBCells;
		MPI_Allreduce(&old_nBcells, &old_sysBCells, 1, MPI_INT, MPI_SUM,
				MPI_COMM_WORLD);

		Pout << "My cells after distribution: "
				<< procedureProfPool::getMesh().nCells() << endl;

		Pout << "Boundary cells: " << old_nBcells << " -> " << nBcells << endl;

		Info << "System boundary cells: " << old_sysBCells << " -> "
				<< sysBCells << endl;

	} else {
		Info << "\nBalance-episode-" << balanceEpisodeID_ << "-in-iteration-"
				<< loadManager::iterationCount_ << " No_balance_needed." << nl
				<< endl;
	}

	balanceEpisodeID_++;

	DLM_PRINT_HEADER_END

	return doBalance;

}

void Foam::loadManager::iterationDone() {

	MPI_Barrier (MPI_COMM_WORLD); // just to make prof times consistent

	//One step ahead of Foam::loadManager::loadBalanceCells();
	loadManager::incIterationCount();

	loadManager::UpdateComputedCellCount();

	if (loadManager::iterationCount_ > 1) {

		ADD_SECTION (updateLoadManager);

		Foam::loadManager::updateBusyIdleTime();

		if (loadManager::PARAMS().enabled)
				Foam::loadManager::updateBalanceEpisode();

		//Minimum trigger period associated with balance period to ensure
		//that we never do a balance+reset without an updateProfiling before
		if (loadManager::iterationCount_ == 2 //force second iteration
				|| (loadManager::iterationCount_
						% int(ceil(PARAMS_.default_balance_period / 2.0)) == 0)) {

			Foam::loadManager::updateProfiling();
			Foam::loadManager::writeLoadData();

		}

		END_SECTION(updateLoadManager);

	}

	loadManager::currentIterationBalanced = false;

	MPI_Barrier(MPI_COMM_WORLD); // just to make prof times consistent

}

bool Foam::loadManager::balance() {

	MPI_Barrier (MPI_COMM_WORLD); // just to make prof times consistent

	bool balanced = false;

	if (loadManager::currentIterationBalanced)
		return false;

	if (loadManager::iterationCount_ > 1)
		if (loadManager::iterationCount_ == 2 //force second iteration
				|| ((loadManager::lastBE + loadManager::PARAMS().balancePeriod)
						<= loadManager::iterationCount_)) {

			ADD_SECTION(loadBalanceCells);

			if (loadManager::PARAMS().enabled)
				balanced = Foam::loadManager::loadBalanceCells();

			// So that times per profiling episode match in logs with
			//and without nSharma
			loadManager::managerInstance_->resetAllOperationTimes();

			loadManager::lastBE = loadManager::iterationCount_;

			END_SECTION(loadBalanceCells);

		}

	loadManager::currentIterationBalanced = true;

	MPI_Barrier(MPI_COMM_WORLD); // just to make prof times consistent

	return balanced;

}

void Foam::loadManager::moveCells(label n, label from, label to, fvMesh& mesh) {

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

/**
 * master method
 * expensive method
 */
void Foam::loadManager::writeLoadData() {

//ADD_SECTION(writeLoadData);

	if (!Pstream::parRun())
		return;

	if (!Pstream::master())
		return;

	loadManager *managerInstance = loadManager::managerInstance_;

	label nProcs = Pstream::nProcs();

	procedureProfPool::mapType& allInfo_ = procedureProfPool::thePool_->map();
	loadManager::loadDataType& loadData = Foam::loadManager::getLoadData();

	DLM_PRINT_HEADER_START

	string lineBreak = "\n";
	int header_offset = 80;
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
							std::fixed << std::setw(11) << std::setprecision(3) << managerInstance->getLoadDataElement(it->second->id0(), itp).totalTime_);

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
						std::fixed << std::setw(11) << std::setprecision(3) << managerInstance->getCurrentCells(itp));

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
						std::fixed << std::setw(11) << std::setprecision(3) << managerInstance->busyTime_[itp]);

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
						std::fixed << std::setw(11) << std::setprecision(3) << managerInstance->idleTime_[itp]);

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
						std::fixed << std::setw(11) << std::setprecision(1) << managerInstance->idlePercentage_[itp]);

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
						std::fixed << std::setw(11) << std::setprecision(1) << managerInstance->busyPercentage_[itp]);

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
						std::fixed << std::setw(11) << std::setprecision(3) << managerInstance->totalBusyTime_[itp]);

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
						std::fixed << std::setw(11) << std::setprecision(3) << managerInstance->totalIdleTime_[itp]);

	}

	cout << line;

	line = SSTR(std::left << "");

	line += lineBreak.c_str();

	/***************************************************************/

	preline = SSTR(std::left << "Busy_Relative_STD(%)");

	line += lineBreak.c_str();

	line += SSTR(std::left << std::setw(header_offset) << preline);

	//label id_amul = allInfo_.find("op_b:AmulCore")->second->id0();

	line +=
			SSTR(
					std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->busyRSD_);
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

	line +=
			SSTR(
					std::fixed << std::setw(9) << std::setprecision(2) << managerInstance->idleRSD_);
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

	DLM_PRINT_HEADER_END

//END_SECTION(writeLoadData);

}

//// Debugging: compare two fields.
//void Foam::loadManager::compareFields(const scalar tolDim,
//		const volVectorField& a, const volVectorField& b) {
//	forAll(a, cellI){
//	if (mag(b[cellI] - a[cellI]) > tolDim) {
//		FatalErrorIn(
//				"compareFields"
//				"(const scalar, const volVectorField&, const volVectorField&)")
//		<< "Did not map volVectorField correctly:" << nl << "cell:"
//		<< cellI << " transfer b:" << b[cellI] << " real cc:"
//		<< a[cellI] << abort(FatalError);
//	}
//}
//forAll(a.boundaryField(), patchI)
//{
//	// We have real mesh cellcentre and
//	// mapped original cell centre.
//
//	const fvPatchVectorField& aBoundary = a.boundaryField()[patchI];
//
//	const fvPatchVectorField& bBoundary = b.boundaryField()[patchI];
//
//	if (!bBoundary.coupled()) {
//		forAll(aBoundary, i)
//		{
//			if (mag(aBoundary[i] - bBoundary[i]) > tolDim) {
//				FatalErrorIn("compareFields"
//						"(const scalar, const volVectorField&"
//						", const volVectorField&)")
//				<< "Did not map volVectorField correctly:" << endl
//				<< "patch:" << patchI << " patchFace:" << i
//				<< " cc:" << endl << "    real    :" << aBoundary[i]
//				<< endl << "    mapped  :" << bBoundary[i] << endl
//				<< abort(FatalError);
//			}
//		}
//	}
//}
//}

// ************************************************************************* //
