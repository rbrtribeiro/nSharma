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
#include "CommGraph.h"
#include "Time.H"
#include <vector>
#include <map>

namespace Foam {

struct nSharma_parameters;
typedef struct nSharma_parameters nSharmaParameters;

class Partitioner {

public:

	 virtual ~Partitioner()=0;

	// TODO: A partitioner may not work with ratios.
	virtual void decompose(List<float>& ratios,std::vector<std::pair<int,labelList> >& decomps)=0;
	virtual void decompChosen(int d)=0;


};

class ParMetisPartitioner: public Partitioner {
protected:
	fvMesh& mesh;
	Time& runTime;
	decompositionMethod* decomposeEngine;

	//<score, id>, score first enables for std::sort
	std::vector<std::pair<int,int> > votes;
	//convert id to itr,tol
	std::map<int, std::pair<float,float> > IDtoDecomp;
	//how many most scored/voted we select
	int nMostVoted;


	void selectMostUsed(std::vector<int>& result);

public:
	ParMetisPartitioner(
			const nSharmaParameters& p,
			CommGraph* cG,
			bool coarseMapping= false);

	~ParMetisPartitioner();

	void decompChosen(int d);
	void smoothWeightCells(scalarField& cellWeights);


	void decompose(List<float>& ratios, std::vector<std::pair<int,labelList> >& decomps);

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
	ParMetisRefinePartitioner(const nSharmaParameters& p, CommGraph* cG);
	~ParMetisRefinePartitioner();

	void decompose(List<float>& ratios, std::vector<std::pair<int, labelList> >& decomps);

};


}

#endif

