#ifndef PUREGAUGEOVERRELAXATIONCONST_H_
#define PUREGAUGEOVERRELAXATIONCONST_H_

#include "LatticeSweep.h"
#include "Environment.h"
#include "utils/RandomSeed.h"
#include "actions/GaugeAction.h"

namespace Update {

class PureGaugeOverrelaxationConst: public Update::LatticeSweep {
public:
	PureGaugeOverrelaxationConst();
	PureGaugeOverrelaxationConst(const PureGaugeOverrelaxationConst& copy);
	~PureGaugeOverrelaxationConst();

	/**
	 * This function performs a pure gauge update on the whole (fundamental) lattice, using the standard Wilson action
	 */
   
 	void execute(environment_t& environment);

	void updateLink(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* action);
   
	void executebeta(environment_t& environment, double Energylowerend, double a);

	void updateLinkbeta(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* action, double& actionE, double Energy, double deltaInterval, double a);

private:
#if NUMCOLORS > 2
	unsigned long int acceptance;
	unsigned long int nsteps;
#ifndef MULTITHREADING
	//The generator of random numbers
	random_generator_t randomGenerator;
	//The generator of random numbers, uniform distribution [0,1]
	random_uniform_generator_t randomUniform;
#endif
#ifdef MULTITHREADING
	//The generator of random numbers
	random_generator_t** randomGenerator;
	//The generator of random numbers, uniform distribution [0,1]
	random_uniform_generator_t** randomUniform;
#endif
#endif
};

} /* namespace Update */
#endif /* PUREGAUGEOVERRELAXATIONCONST_H_ */
