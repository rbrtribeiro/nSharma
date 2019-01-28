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
#include <vector>

namespace Foam {

struct nSharma_parameters;
typedef struct nSharma_parameters nSharmaParameters;

class PerformanceModel {
public:

	virtual void update()=0;
	virtual List<float> getTperCell()=0;

};


class TperCellPerformanceModel : public PerformanceModel {

private:
	std::vector<std::vector<float> > last_TperCell;
	std::vector<float> window_weights;
	bool first;

public:

	List<float> TperCell;


	TperCellPerformanceModel(int);
	virtual ~TperCellPerformanceModel();
	virtual void update();
	void windowAverage(List<float> raw,
			std::vector<std::vector<float> >& last_,
			List<float>& r);

	virtual List<float> getTperCell();

};

}

#endif

