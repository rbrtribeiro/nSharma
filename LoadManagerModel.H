/*---------------------------------------------------------------------------*\
Class
 Foam::LoadManagerModel

 Description

 SourceFiles
 LoadManagerModel.C

 \*---------------------------------------------------------------------------*/
#ifndef LoadManagerModel_H
#define LoadManagerModel_H

namespace Foam {

class Partitioner;

struct nSharma_parameters;
typedef struct nSharma_parameters nSharmaParameters;

class LoadManagerModel {
public:

	bool virtual Balance(labelList&)=0;
	void virtual PostLoadMigrationUpdate()=0;

};

}

#endif

