/*---------------------------------------------------------------------------*\
Class
 Foam::Partitioner

 Description

 SourceFiles
 Partitioner.C

 \*---------------------------------------------------------------------------*/
#ifndef Partitioner_H
#define Partitioner_H

#include "List.H"
#include "labelList.H"
#include "fvMesh.H"
#include "decompositionMethod.H"
#include "foamTime.H"
#include <vector>

namespace Foam {

class Partitioner {

public:

	 virtual ~Partitioner()=0;

	// TODO: A partitioner may not work with ratios.
	virtual void decompose(List<float>& ratios,std::vector<labelList>& decomps)=0;


};

class ParMetisPartitioner: public Partitioner {
protected:
	fvMesh& mesh;
	Time& runTime;
	decompositionMethod* decomposeEngine;

public:
	ParMetisPartitioner();
	~ParMetisPartitioner();

	void decompose(List<float>& ratios, std::vector<labelList>& decomps);

};

//class PTscotchPartitioner: public Partitioner {
//
//	fvMesh& mesh;
//	Time& runTime;
//	decompositionMethod* decomposeEngine;
//
//public:
//	PTscotchPartitioner();
//	~PTscotchPartitioner();
//
//	void decompose(List<float>& ratios, labelList&);
//
//};

class ParMetisRefinePartitioner: public ParMetisPartitioner {


public:
	ParMetisRefinePartitioner();
	~ParMetisRefinePartitioner();

	void decompose(List<float>& ratios, std::vector<labelList>& decomps);

};


}

#endif

