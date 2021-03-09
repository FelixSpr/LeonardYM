#ifndef SIMULATION_H_
#define SIMULATION_H_
#include "Environment.h"
#include "LatticeSweep.h"
#include <list>

#include <iostream>
#include <iomanip>
#include <fstream>

namespace Update {

class Simulation {
public:
	//Create a Monte Carlo simulation for a given environment
	Simulation(const environment_t& _environment);
	~Simulation();

	//This function set the lattice, the environment and all other factors for the simulation
	void starter();
  void starterLLR(double Eint);

	//Execute the warm-up sweeps to thermalize the lattice
	void warmUp();
  void warmUpLLR(double a, double Energylowerend);

	//Run the main measurement-update sweeps
	void measurement();
  void measurementLLR(double a, double Energylowerend);

private:
	std::list<LatticeSweep*> listWarmUpSweeps;
	std::list<LatticeSweep*> listMeasurementSweeps;
	environment_t environment;
};

} /* namespace Update */
#endif /* SIMULATION_H_ */
