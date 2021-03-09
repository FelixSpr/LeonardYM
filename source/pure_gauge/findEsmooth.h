#ifndef FINDESMOOTH_H_
#define FINDESMOOTH_H_

#include "actions/GaugeAction.h"
#include "PureGaugeUpdater.h"
#include "LatticeSweep.h"

namespace Update {

class findEsmooth {
public:
  findEsmooth();
	//findEsmooth(const findEsmooth&);
	//~findEsmooth();
 
  void execute(environment_t& environment, double Eint);
 
};
}

#endif /* FINDESMOOTH_H_ */