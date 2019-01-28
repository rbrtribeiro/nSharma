#ifndef nSharma_H
#define nSharma_H

#include "label.H"
#include "scalar.H"
#include "Time.H"
#include "boundBox.H"
#include "DynamicField.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"

#include "mpi.h"

#include <vector>
#include <map>

#include "procedureProfInfo.H"
#include "procedureProfStack.H"
#include "PerformanceModel.H"
#include "LoadManagerModel.H"
#include "PowerManagerModel.h"
#include "CommGraph.h"
#include "nSharma_common.h"
#include "PowerInterface.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define WINDOW_BUSY_RSD_SIZE 4

namespace Foam {

typedef struct nSharma_parameters {
	bool enableLoadBalance;
	bool enablePowerOptimize;
	char targetWorkloadSection[64];
	int balancePeriod;
	float minPercentageCellsMoved;
	uint window;
	uint simulation_iterations;
	float minimal_gain;
	uint default_balance_period;
	float powerCap;
	bool useCommGraphPartitioning;

} nSharmaParameters;

class fvMesh;
class polyMesh;
class MPIfiedProcedure;

/*---------------------------------------------------------------------------*\
                           Class procedureProfPool Declaration
 \*---------------------------------------------------------------------------*/
/**
 * Master process only instance
 */
class nSharma {
public:
	typedef std::vector<std::vector<MPIfiedProcedure> > loadDataType;

	template<class T>
	class resizable {
	public:
		typedef DynamicList<T, 0, 11, 10> ListType;
		typedef DynamicField<T, 0, 11, 10> FieldType;
	};

	/**
	 * this should be in procedureProfInfo. Avoids OpenFOAM MPI dependency
	 */
	static MPI_Datatype MPI_LOADPROCEDURE;

	bool forceNextBalance;

private:
	// Private data

	loadDataType loadData;
	std::vector<long> worldEpisodeComputedCellCount;
	std::vector<long> worldCurrentCellCount;
	long localEpisodeComputedCellCount;

	std::vector<scalar> busyPercentage_;
	std::vector<scalar> idlePercentage_;

	std::vector<scalar> busyTime_;
	std::vector<scalar> idleTime_;

	std::vector<scalar> totalBusyTime_;
	std::vector<scalar> totalIdleTime_;

	float busyRSD_;
	float idleRSD_;
	long nSysBfaces_;
	long nSysBfacesNodes_;

	std::vector<float> window_busyRSD;
	std::vector<float> window_busyRSD_x;


	// Tolerance (as fraction of the bounding box). Needs to be fairly lax since
	// usually meshes get written with limited precision (6 digits)
	const scalar defaultMergeTol = 1E-6;

	const dictionary& decompositionDict_;
	dictionary meshDict_;

	uint iterationCount_;
	int balanceEpisodeID_;
	int lastBE; //last balance episode iteration id
	bool currentIterationBalanced;

	fvMesh* mesh_;
	Time* time;
	CommGraph* commGraph;
	PerformanceModel* performanceModel;
	LoadManagerModel* loadManager;
	nSharmaParameters PARAMS_;

	PowerInterface* powerInterface;
	PowerManagerModel* powerManager;

	static nSharma* instance_;

public:

	// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

	nSharma(dictionary& decompositionDict_, Time& t, fvMesh& m,
			nSharmaParameters p);

	~nSharma();

	// * * * * * * * * * * * * * * * API  * * * * * * * * * * * * * //

	static void Init(Time& t, fvMesh& m);

	static bool Iteration();

	static void IterationDone();

	static void Finalize();

	// * * * * * * * * * * * * * * * get/set Functions  * * * * * * * * * * * * * //

	CommGraph* getCommGraph() {
		return commGraph;
	}

	inline uint iterationCount() {
		return iterationCount_;
	}

	inline uint balanceEpisodeID() {
		return balanceEpisodeID_;
	}

	inline void setPreviousEpisodeComputedCells(std::vector<long>& c) {

		nSharma::instance_->worldEpisodeComputedCellCount = c;

	}

	inline void setCurrentCells(std::vector<long>& c) {

		nSharma::instance_->worldCurrentCellCount = c;

	}

	inline int getCurrentCells(label n) {
		return nSharma::instance_->worldCurrentCellCount.at(n);
	}

	inline int getPreviousEpisodeComputedCells(label n) {
		return nSharma::instance_->worldEpisodeComputedCellCount.at(n);
	}

	inline void incIterationCount() {
		iterationCount_++;
	}

	inline void resetComputedCellCount() {

		nSharma::instance_->localEpisodeComputedCellCount = 0;

	}

	inline long getComputedCellCount() {

		return nSharma::instance_->localEpisodeComputedCellCount;

	}

	inline dictionary& getMeshDict() {
		return meshDict_;
	}

	MPIfiedProcedure& getLoadDataElement(label opID, label proc);

	const nSharmaParameters& PARAMS();

	nSharmaParameters& editPARAMS();

	static nSharma* getInstance();

	LoadManagerModel* getLoadManagerInstance();

	loadDataType& getLoadData();

	void enablePool(bool v);

	void setTargetWorkloadSectionName(string v);

	scalar getMergeDistance(const boundBox& bb);

private:

	// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

	template<class Type> static void correctBoundaries();

	void iterationDone();

	bool balance();

	void powerOptimize();

	void resetAllOperationTimes();

	void reconfigureTreeComm();

	void printTreeComm();

	void setLoadDataProcBuffer(MPIfiedProcedure* infos, label proc);

	label boilProfilingPool(MPIfiedProcedure*&);

	void updateProfiling();

	void updateBusyIdleTime();

	void updateBalanceEpisode();

	void updateRSD();

	void updateRSD_busyTime();

	void updateRSD_MIN_MAX();

	void writeLoadData();

	void printMeshData(Ostream& os, const polyMesh& mesh);

	void redistribute(labelList finalDecomp);

	bool loadBalanceCells();

	void moveCells(label n, label from, label to, fvMesh& mesh);

	void updateComputedCellCount();

};

}

template<class Type>
void Foam::nSharma::correctBoundaries() {
	typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

	HashTable<GeoField*> flds(
			nSharma::getInstance()->mesh_->objectRegistry::lookupClass<GeoField>());

forAllIter(typename HashTable<GeoField*>, flds, iter)
{
	GeoField& fld = *iter();

	//mimic "evaluate" but only for coupled patches (processor or cyclic)
	// and only for blocking or nonBlocking comms (no scheduled comms)
	if
	(
			Pstream::defaultCommsType == Pstream::blocking
			|| Pstream::defaultCommsType == Pstream::nonBlocking
	)
	{
		label nReq = Pstream::nRequests();

		forAll(fld.boundaryField(), patchi)
		{
			if(fld.boundaryField()[patchi].coupled())
			{
				fld.boundaryField()[patchi].initEvaluate
				(
						Pstream::defaultCommsType
				);
			}
		}

		// Block for any outstanding requests
		if
		(
				Pstream::parRun()
				&& Pstream::defaultCommsType == Pstream::nonBlocking
		)
		{
			Pstream::waitRequests(nReq);
		}

		forAll(fld.boundaryField(), patchi)
		{
			if(fld.boundaryField()[patchi].coupled())
			{
				fld.boundaryField()[patchi].evaluate
				(
						Pstream::defaultCommsType
				);
			}
		}
	}
	else
	{
		//Scheduled patch updates not supported
		FatalErrorIn
		(
				"dynamicRefineBalancedFvMeshTemplates::correctBoundaries"
		) << "Unsuported communications type "
		<< Pstream::commsTypeNames[Pstream::defaultCommsType]
		<< exit(FatalError);
	}
}
}

#endif

