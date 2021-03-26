#ifndef PLAQUETTE_H_
#define PLAQUETTE_H_

#include "LatticeSweep.h"

namespace Update {

class Plaquette: public Update::LatticeSweep {
public:
	Plaquette();
	~Plaquette();

	/**
	 * This function measure the plaquette on the lattice and it prints out its value.
	 * @param enviroment
	 * @param sweep
	 * @param n
	 */
	virtual void execute(environment_t& environment);
 
   

	static long_real_t temporalPlaquette(const extended_gauge_lattice_t& gaugeLinkConfiguration);
  long_real_t spacialPlaquette(environment_t& environment);
  long_real_t Plaquettevalue(environment_t& environment);
  long_real_t Plaquettestd(environment_t& environment,double average);
};

} /* namespace Update */
#endif /* PLAQUETTE_H_ */
