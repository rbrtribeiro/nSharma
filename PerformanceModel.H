/*---------------------------------------------------------------------------*\
Class
 Foam::PerformanceModel

 Description

 SourceFiles
 PerformanceModel.C

 \*---------------------------------------------------------------------------*/
#ifndef PerformanceModel_H
#define PerformanceModel_H

#include "List.H"
#include "labelList.H"

#include <cmath>
#include <vector>
#include "fvMeshDistribute.H"


namespace Foam {

class Partitioner;

struct loadManager_parameters;
typedef struct loadManager_parameters loadManagerParameters;

class PerformanceModel {
protected:
	class decompDetails {
	public:
		labelList decomposition;
		labelList count;
		long totalReceived;
		long totalSend;
		long totalcellsTransfered;
		bool apply;
		float Tstar;
		float Tmigrt;

		decompDetails(labelList& newD){
			decomposition = newD;
			count = fvMeshDistribute::countCells(decomposition);
		}

	};

public:
	PerformanceModel(loadManagerParameters&);
	virtual ~PerformanceModel()=0;

	bool virtual Balance(labelList&)=0;
	void virtual PostLoadMigrationUpdate()=0;

	loadManagerParameters& PARAMS_;


};

}

#endif

