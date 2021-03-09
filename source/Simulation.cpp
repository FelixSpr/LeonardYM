#include "Simulation.h"
#include "starters/StartGaugeConfiguration.h"
#include "io/GlobalOutput.h"
#include <sys/time.h>

#include "actions/GaugeAction.h"
#include "pure_gauge/findEsmooth.h"
#include "pure_gauge/PureGaugeUpdaterConst.h"
#include "pure_gauge/PureGaugeOverrelaxation.h"

#include <iostream>
#include <iomanip>
#include <fstream>

namespace Update {

Simulation::Simulation(const environment_t& _environment) : environment(_environment) { }

Simulation::~Simulation() {
	std::list<LatticeSweep*>::iterator i;
	for (i = listWarmUpSweeps.begin(); i != listWarmUpSweeps.end(); ++i) {
		delete *i;
	}
	for (i = listMeasurementSweeps.begin(); i != listMeasurementSweeps.end(); ++i) {
		delete *i;
	}
}

void Simulation::starter() {
	if (isOutputProcess()) std::cout << "Starting the simulation ..." << std::endl;
	StartGaugeConfiguration* cl = StartGaugeConfiguration::getInstance(environment.configurations.get<std::string>("start"));
	cl->execute(environment);
	delete cl;
  //findEsmooth* finder = new findEsmooth();
  //finder->execute(environment);
}

void Simulation::starterLLR(double Eint) {
	if (isOutputProcess()) std::cout << "Starting the simulation ..." << std::endl;
	StartGaugeConfiguration* cl = StartGaugeConfiguration::getInstance(environment.configurations.get<std::string>("start"));
	cl->execute(environment);
	delete cl;
  findEsmooth* finder = new findEsmooth();
  finder->execute(environment,Eint);
}

void Simulation::warmUp() {
	if (isOutputProcess()) std::cout << "Loading the warm-up sweeps ..." << std::endl;
	std::vector<std::string> warmUpParameters = environment.configurations.get< std::vector< std::string > >("warm_up_sweeps");
	std::vector<std::string>::iterator i;
	for (i = warmUpParameters.begin(); i != warmUpParameters.end(); ++i) {
		listWarmUpSweeps.push_back(LatticeSweep::read(*i));
	}
	if (isOutputProcess()) std::cout << "Doing the warm-up sweeps ..." << std::endl;
	//Set the name for the output
	GlobalOutput* globalOutput = GlobalOutput::getInstance();
	globalOutput->setBaseName(environment.configurations.get<std::string>("output_name"));
	globalOutput->setBaseFolder(environment.configurations.get<std::string>("output_directory_measurements"));
	environment.measurement = false;
	unsigned int numberWarmUpSweeps = environment.configurations.get<unsigned int>("number_warm_up_sweeps");
	environment.sweep = 0;
	for (unsigned int i = 0; i < numberWarmUpSweeps; ++i) {
		timeval start, stop, result;
		gettimeofday(&start,NULL);
		std::list<LatticeSweep*>::iterator sweep;
		for (sweep = listWarmUpSweeps.begin(); sweep != listWarmUpSweeps.end(); ++sweep) {
			(*sweep)->call(environment);
		}
		gettimeofday(&stop,NULL);
		timersub(&stop,&start,&result);
		if (isOutputProcess()) std::cout << "Sweep cicle " << i << " done in: " << (double)result.tv_sec + result.tv_usec/1000000.0 << " sec" << std::endl;
		++environment.sweep;
	}
}

void Simulation::warmUpLLR(double a, double Energylowerend) {
	if (isOutputProcess()) std::cout << "Loading the warm-up sweeps ..." << std::endl;
  PureGaugeUpdaterConst* update = new PureGaugeUpdaterConst();
  GaugeAction* gaugeActionEoutput = GaugeAction::getInstance("StandardWilson",environment.configurations.get<real_t>("beta"));
	std::vector<std::string> warmUpParameters = environment.configurations.get< std::vector< std::string > >("warm_up_sweeps");
	std::vector<std::string>::iterator i;
	for (i = warmUpParameters.begin(); i != warmUpParameters.end(); ++i) {
		listWarmUpSweeps.push_back(LatticeSweep::read(*i));
	}
	if (isOutputProcess()) std::cout << "Doing the warm-up sweeps ..." << std::endl;
	//Set the name for the output
	GlobalOutput* globalOutput = GlobalOutput::getInstance();
	globalOutput->setBaseName(environment.configurations.get<std::string>("output_name"));
	globalOutput->setBaseFolder(environment.configurations.get<std::string>("output_directory_measurements"));
	environment.measurement = false;
	unsigned int numberWarmUpSweeps = environment.configurations.get<unsigned int>("number_warm_up_sweeps");
	environment.sweep = 0;
	for (unsigned int i = 0; i < numberWarmUpSweeps; ++i) {
		timeval start, stop, result;
		//gettimeofday(&start,NULL);
		//std::list<LatticeSweep*>::iterator sweep;
		//for (sweep = listWarmUpSweeps.begin(); sweep != listWarmUpSweeps.end(); ++sweep) {
		//	(*sweep)->call(environment);
		//}
		gettimeofday(&stop,NULL);
		timersub(&stop,&start,&result);
	//	if (isOutputProcess()) std::cout << "Sweep cicle " << i << " done in: " << (double)result.tv_sec + result.tv_usec/1000000.0 << " sec" << std::endl;
    update->executebeta(environment,a*environment.configurations.get<real_t>("beta"),Energylowerend);
		++environment.sweep;
	}
}

void Simulation::measurement() {
	if (isOutputProcess()) std::cout << "Loading the measurement sweeps ..." << std::endl;
	std::vector<std::string> measurementParameters = environment.configurations.get< std::vector< std::string > >("measurement_sweeps");
	std::vector<std::string>::iterator i;
	for (i = measurementParameters.begin(); i != measurementParameters.end(); ++i) {
		listMeasurementSweeps.push_back(LatticeSweep::read(*i));
	}
	if (isOutputProcess()) std::cout << "Doing the measurement sweeps ..." << std::endl;
	unsigned int numberMeasurementSweeps = environment.configurations.get<unsigned int>("number_measurement_sweeps");
	environment.sweep = 0;
	environment.measurement = true;
	//Set the name for the output
	GlobalOutput* globalOutput = GlobalOutput::getInstance();
	globalOutput->setBaseName(environment.configurations.get<std::string>("output_name"));
	globalOutput->setBaseFolder(environment.configurations.get<std::string>("output_directory_measurements"));
	for (unsigned int i = 0; i < numberMeasurementSweeps; ++i) {
		timeval start, stop, result;
		gettimeofday(&start,NULL);
		std::list<LatticeSweep*>::iterator sweep;
		for (sweep = listMeasurementSweeps.begin(); sweep != listMeasurementSweeps.end(); ++sweep) {
			(*sweep)->call(environment);
		}
		gettimeofday(&stop,NULL);
		timersub(&stop,&start,&result);
		if (isOutputProcess()) {
			std::cout << "Sweep cicle " << i << " done in: " << (double)result.tv_sec + result.tv_usec/1000000.0 << " sec" << std::endl;
			//Save the data
			globalOutput->print();
		}
		++environment.sweep;
	}
  //GaugeAction* gaugeActionEfind = GaugeAction::getInstance("StandardWilson",5.6);
	//long_real_t energyEfind = gaugeActionEfind->energy(environment);
  //double energyout = energyEfind;
  //std::cout << "Gaugeaction in simulation.cpp" << energyEfind << std::endl;
}

void Simulation::measurementLLR(double a, double Energylowerend) {
  std::ofstream ofs ("io_directory/energyvalues.txt", std::ofstream::out);

  //ofs << "lorem ipsum";

  //ofs.close();
	if (isOutputProcess()) std::cout << "Loading the measurement sweeps ..." << std::endl;
  //std::ofstream myfileenergyvalues;
  //myfileenergyvalues.open("energyvalues.txt");
  PureGaugeUpdaterConst* update = new PureGaugeUpdaterConst();
  PureGaugeOverrelaxation* updateover=new PureGaugeOverrelaxation();
  GaugeAction* gaugeActionEoutput = GaugeAction::getInstance("StandardWilson",environment.configurations.get<real_t>("beta"));
	//long_real_t energyEoutput = gaugeActionEoutput->energy(environment);
	//myfileenergyvalues.open("energyvalues.txt");
	std::vector<std::string> measurementParameters = environment.configurations.get< std::vector< std::string > >("measurement_sweeps");
	std::vector<std::string>::iterator i;
	for (i = measurementParameters.begin(); i != measurementParameters.end(); ++i) {
		listMeasurementSweeps.push_back(LatticeSweep::read(*i));
	}
	if (isOutputProcess()) std::cout << "Doing the measurement sweeps ..." << std::endl;
	unsigned int numberMeasurementSweeps = environment.configurations.get<unsigned int>("number_measurement_sweeps");
	environment.sweep = 0;
	environment.measurement = true;
	//Set the name for the output
	GlobalOutput* globalOutput = GlobalOutput::getInstance();
	globalOutput->setBaseName(environment.configurations.get<std::string>("output_name"));
	globalOutput->setBaseFolder(environment.configurations.get<std::string>("output_directory_measurements"));
	for (unsigned int i = 0; i < numberMeasurementSweeps; ++i) {
		timeval start, stop, result;
		gettimeofday(&start,NULL);
		//std::list<LatticeSweep*>::iterator sweep;
		//for (sweep = listMeasurementSweeps.begin(); sweep != listMeasurementSweeps.end(); ++sweep) {
		//	(*sweep)->call(environment);
		//}
		gettimeofday(&stop,NULL);
		timersub(&stop,&start,&result);
		if (isOutputProcess()) {
		//	std::cout << "Sweep cicle " << i << " done in: " << (double)result.tv_sec + result.tv_usec/1000000.0 << " sec" << std::endl;
			//Save the data
		//	globalOutput->print();
		}
    //std::cout << "energyvaluesfile: " << gaugeActionEoutput->energy(environment) << " \n";
    
    //myfileenergyvalues << "random stuff \n"; //<< gaugeActionEoutput->energy(environment) << " \n";
    update->executebeta(environment,a*environment.configurations.get<real_t>("beta"),Energylowerend);
    std::cout << "Energy after Heatbathconst: " << gaugeActionEoutput->energy(environment) << " \n";
    updateover->execute(environment);
    std::cout << "Energy after Overrelax: " << gaugeActionEoutput->energy(environment) << " \n";
    ofs << gaugeActionEoutput->energy(environment) << " \n";
		++environment.sweep;
	}
  //myfileenergyvalues.close();
  ofs.close();
  //GaugeAction* gaugeActionEfind = GaugeAction::getInstance("StandardWilson",5.6);
	//long_real_t energyEfind = gaugeActionEfind->energy(environment);
  //double energyout = energyEfind;
  //std::cout << "Gaugeaction in simulation.cpp" << energyEfind << std::endl;
}

} /* namespace Update */
