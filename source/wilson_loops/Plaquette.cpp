#include "Plaquette.h"
#include <iostream>
#include "io/GlobalOutput.h"

namespace Update {

Plaquette::Plaquette() : LatticeSweep() { }

Plaquette::~Plaquette() { }

void Plaquette::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0., plaquette_t = 0., plaquette_s = 0.;

#pragma omp parallel for reduction(+:plaquette, plaquette_t, plaquette_s)
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				real_t tmp = real(trace(environment.gaugeLinkConfiguration[site][mu]*environment.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(environment.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(environment.gaugeLinkConfiguration[site][nu])));
				plaquette += tmp;
				if (mu == 3 || nu == 3) {
					plaquette_t += tmp;
				}
				else {
					plaquette_s += tmp;
				}
			}
		}
	}
	reduceAllSum(plaquette);
	reduceAllSum(plaquette_s);
	reduceAllSum(plaquette_t);

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("plaquette");

		typedef extended_gauge_lattice_t::Layout Layout;
		std::cout << "Plaquette expectation value: " << plaquette/(6.*numberColors*Layout::globalVolume) << std::endl;
		output->write("plaquette", plaquette/(6.*numberColors*Layout::globalVolume));
		output->write("plaquette", plaquette_s/(3.*numberColors*Layout::globalVolume));
		output->write("plaquette", plaquette_t/(3.*numberColors*Layout::globalVolume));

		output->pop("plaquette");
	}

}

long_real_t Plaquette::temporalPlaquette(const extended_gauge_lattice_t& gaugeLinkConfiguration) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0.;
#pragma omp parallel for reduction(+:plaquette)
	for (int site = 0; site < gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 3; ++mu) {
			unsigned int nu = 3;
			plaquette += real(trace(gaugeLinkConfiguration[site][mu]*gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(gaugeLinkConfiguration[site][nu])));
		}
	}
	reduceAllSum(plaquette);
	if (isOutputProcess()) {
		typedef extended_gauge_lattice_t::Layout Layout;
		plaquette = plaquette/(3.*numberColors*Layout::globalVolume);
		std::cout << "Temporal Plaquette expectation value: " << plaquette << std::endl;
	}
	return plaquette;
}



long_real_t Plaquette::spacialPlaquette(environment_t& environment) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0., plaquette_t = 0., plaquette_s = 0.;

#pragma omp parallel for reduction(+:plaquette, plaquette_t, plaquette_s)
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				real_t tmp = real(trace(environment.gaugeLinkConfiguration[site][mu]*environment.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(environment.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(environment.gaugeLinkConfiguration[site][nu])));
				plaquette += tmp;
				if (mu == 3 || nu == 3) {
					plaquette_t += tmp;
				}
				else {
					plaquette_s += tmp;
				}
			}
		}
	}
	reduceAllSum(plaquette);
	reduceAllSum(plaquette_s);
	reduceAllSum(plaquette_t);
  typedef extended_gauge_lattice_t::Layout Layout;
  //plaquette_s=plaquette_s/(3.*numberColors*Layout::globalVolume);
  plaquette_s=plaquette_s/(3.*Layout::globalVolume);
	return plaquette_s;

}


long_real_t Plaquette::Plaquettevalue(environment_t& environment) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0., plaquette_t = 0., plaquette_s = 0.;

#pragma omp parallel for reduction(+:plaquette, plaquette_t, plaquette_s)
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				real_t tmp = real(trace(environment.gaugeLinkConfiguration[site][mu]*environment.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(environment.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(environment.gaugeLinkConfiguration[site][nu])));
				plaquette += tmp;
				if (mu == 3 || nu == 3) {
					plaquette_t += tmp;
				}
				else {
					plaquette_s += tmp;
				}
			}
		}
	}
	reduceAllSum(plaquette);
	reduceAllSum(plaquette_s);
	reduceAllSum(plaquette_t);

  typedef extended_gauge_lattice_t::Layout Layout;
  //plaquette = plaquette/(6.*numberColors*Layout::globalVolume);
  plaquette = plaquette/(6.*Layout::globalVolume);
	return plaquette;

}

long_real_t Plaquette::Plaquettestd(environment_t& environment,double average) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0., plaquette_t = 0., plaquette_s = 0., plaquettestandard = 0.;
  //long real_t plaquettestandard = 0.;
#pragma omp parallel for reduction(+:plaquette, plaquette_t, plaquette_s)
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				real_t tmp = real(trace(environment.gaugeLinkConfiguration[site][mu]*environment.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(environment.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(environment.gaugeLinkConfiguration[site][nu])));
        reduceAllSum(tmp);
        //std::cout << "plaquette at site: " << tmp << std::endl;
        plaquettestandard += pow(average- tmp,2);
				plaquette += tmp;
				if (mu == 3 || nu == 3) {
					plaquette_t += tmp;
				}
				else {
					plaquette_s += tmp;
				}
			}
		}
	}
  
	reduceAllSum(plaquette);
	reduceAllSum(plaquette_s);
	reduceAllSum(plaquette_t);

  typedef extended_gauge_lattice_t::Layout Layout;
  //plaquette = plaquette/(6.*numberColors*Layout::globalVolume);
  plaquette = plaquette/(6.*Layout::globalVolume);
  plaquettestandard = sqrt(plaquettestandard/(6.*Layout::globalVolume - 1.0));
	return plaquettestandard;

}


} /* namespace Update */
