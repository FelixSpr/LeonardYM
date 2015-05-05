/*
 * Glueball.h
 *
 *  Created on: Oct 9, 2012
 *      Author: spiem_01
 */

#ifndef GLUEBALL_H_
#define GLUEBALL_H_
#include "LatticeSweep.h"

namespace Update {

class Glueball : public LatticeSweep {
public:
	Glueball();
	~Glueball();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* GLUEBALL_H_ */