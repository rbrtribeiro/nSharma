/*---------------------------------------------------------------------------*\
Class
 Foam::loadManager

 Description

 SourceFiles
 loadManager.C

 \*---------------------------------------------------------------------------*/
#ifndef loadManager_H
#define loadManager_H

#include "label.H"
#include "scalar.H"
#include "mpi.h"
#include <vector>
#include <map>

#include "procedureProfInfo.H"
#include "procedureProfStack.H"
//#include "MPIfiedProcedure.H"

//#include "fvMesh.H"
//#include "fvCFD.H"
#include "foamTime.H"
#include "boundBox.H"

#include "DynamicField.H"

#include "PerformanceModel.H"

#define WINDOW_BUSY_RSD_SIZE 4

#define DLM_PRINT_HEADER_START Info\
			<< "\n//****************************************** LOAD_MANAGER ******************************************//\n\n";
#define DLM_PRINT_HEADER_END Info\
			<< "\n//**************************************** END-LOAD_MANAGER ****************************************//\n\n";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

typedef struct loadManager_parameters {
	bool enabled; //1
	string targetWorkloadSection;
	int balancePeriod; //5
	float minPercentageCellsMoved; // 0.1 * mesh.nCells()
	uint window;
	uint simulation_iterations;
	float minimal_gain;
	uint default_balance_period;
} loadManagerParameters;

class fvMesh;
class polyMesh;
class MPIfiedProcedure;

/*---------------------------------------------------------------------------*\
                           Class procedureProfPool Declaration
 \*---------------------------------------------------------------------------*/
/**
 * Master process only instance
 */
class loadManager {
	// Private typedefs
public:
	typedef std::vector<std::vector<MPIfiedProcedure> > loadDataType;

	/**
	 * this should be in procedureProfInfo. Avoids OpenFOAM MPI dependency
	 */
	static MPI_Datatype MPI_LOADPROCEDURE;

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

	std::vector<float> window_busyRSD;
	std::vector<float> window_busyRSD_x;

	static loadManager *managerInstance_;
	static const scalar defaultMergeTol;

	template<class T>
	class resizable {
	public:
		typedef DynamicList<T, 0, 11, 10> ListType;
		typedef DynamicField<T, 0, 11, 10> FieldType;
	};

	static loadManagerParameters PARAMS_;

	const dictionary& decompositionDict_;
	dictionary meshDict_;

	static PerformanceModel* performanceModel;

	static uint iterationCount_;
	static int balanceEpisodeID_;
	static int lastBE; //last balance episode iteration id
	static bool currentIterationBalanced;

public:

	static inline uint iterationCount() {
		return iterationCount_;
	}

	static inline uint balanceEpisodeID() {
		return balanceEpisodeID_;
	}

	static void iterationDone();

	static bool balance();

	MPIfiedProcedure& getLoadDataElement(label opID, label proc);

	static const loadManagerParameters& PARAMS();

	void resetAllOperationTimes();

	//- Construct null
	loadManager(dictionary& decompositionDict_);

	//- Destructor
	~loadManager();

	// Member functions

	static loadManager* getManagerInstance();
	static PerformanceModel* getPerformanceModelInstance();
	static loadDataType& getLoadData();

	static void init(Time& t, fvMesh& m);

	static void finalize();

	void setLoadDataProcBuffer(MPIfiedProcedure* infos, label proc);

	static label BoilProfilingPool(MPIfiedProcedure*&);

	static void updateProfiling();

	static void updateBusyIdleTime();

	static void updateBalanceEpisode();

	void updateRSD();

	void updateRSD_busyTime();

	void updateRSD_MIN_MAX();

	static void writeLoadData();

	static void enablePool(bool v);

	static void setTargetWorkloadSectionName(string v);

	static scalar getMergeDistance(	const boundBox& bb);

	//static void compareFields(const scalar tolDim, const volVectorField& a,
	//		const volVectorField& b);

	static void printMeshData(Ostream& os, const polyMesh& mesh);

	static void redistribute(labelList finalDecomp);

	static bool loadBalanceCells();

	static void moveCells(label n, label from, label to, fvMesh& mesh);

	static void UpdateComputedCellCount();

	static inline void setPreviousEpisodeComputedCells(std::vector<long>& c) {

		loadManager::managerInstance_->worldEpisodeComputedCellCount = c;

	}

	static inline void setCurrentCells(std::vector<long>& c) {

		loadManager::managerInstance_->worldCurrentCellCount = c;

	}

	static inline int getCurrentCells(label n) {
		return loadManager::managerInstance_->worldCurrentCellCount.at(n);
	}

	static inline int getPreviousEpisodeComputedCells(label n) {
		return loadManager::managerInstance_->worldEpisodeComputedCellCount.at(
				n);
	}

	static inline void incIterationCount() {
		iterationCount_++;
	}

	static inline void ResetComputedCellCount() {

		loadManager::managerInstance_->localEpisodeComputedCellCount = 0;

	}

	static inline long GetComputedCellCount() {

		return loadManager::managerInstance_->localEpisodeComputedCellCount;

	}

	dictionary& getMeshDict() {
		return meshDict_;
	}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

