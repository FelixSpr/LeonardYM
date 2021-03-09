#include "findEsmooth.h"

namespace Update {


findEsmooth::findEsmooth() {
  //double Energy = -2500;
  //double delta = 200;
}

void findEsmooth::execute(environment_t & environment, double Eint) {
  
  
  bool Efound = false;
  double Energytarget = Eint;
  double delta = environment.configurations.get<double>("delta");
  //double energyref;
  double betaEfind = environment.configurations.get<double>("beta");
  double betaref = environment.configurations.get<double>("beta");
  PureGaugeUpdater* update = new PureGaugeUpdater();
  //update->executebeta(environment,betaEfind);
  GaugeAction* gaugeActionEfind = GaugeAction::getInstance("StandardWilson",betaref);
	long_real_t energyEfind = gaugeActionEfind->energy(environment);
  int counter = 0;
  std::cout << "GaugeEnergy::Energy value " << energyEfind << std::endl;
  if(energyEfind>=Energytarget && energyEfind <= (Energytarget+delta))
  {
    Efound = true;
    //betaEfind = betaref;
  }
  
  while(Efound == false && counter<4000)
  {
    update->executebeta(environment,betaEfind);
    //energy = action();
    energyEfind = gaugeActionEfind->energy(environment);
    //node0_printf("energy %.8g %.8g %.8g\n",
      //           energyref, beta, counter);
    //node0_printf("Eint = %.4g \n", Eint);
    //node0_printf("energy %.8g %.8g %.8g\n", energy, Eint, counter);
    //node0_printf("energyref %.8g %.8g %.8g\n", energyref, beta, counter);
    if(energyEfind>=Energytarget && energyEfind <= (Energytarget+delta))
    {
      Efound = true;
      //beta = betaref;
      std::cout << "Beta: " << betaEfind << " Counter: " << counter << std::endl;
    }
    else if(energyEfind>(Energytarget+delta))
    {
      betaEfind = betaEfind+0.1;
      std::cout << "Beta: " << betaEfind << " Counter: " << counter << std::endl;
    }
    else if(energyEfind<=(Energytarget+delta))
    {
      betaEfind = betaEfind-0.1;
      std::cout << "Beta: " << betaEfind << " Counter: " << counter << std::endl;
    }
    counter = counter +1;
  }
  std::cout << "Counter at the end: " << counter << std::endl; 
  
  
}

}