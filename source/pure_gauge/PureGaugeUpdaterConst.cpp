#include "PureGaugeUpdaterConst.h"
#include "MatrixTypedef.h"
#include "Checkerboard.h"
#include "utils/RandomSeed.h"

namespace Update {

#ifdef MULTITHREADING
PureGaugeUpdaterConst::PureGaugeUpdaterConst() {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif
#ifndef MULTITHREADING
PureGaugeUpdaterConst::PureGaugeUpdaterConst() : randomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif

#ifdef MULTITHREADING
PureGaugeUpdaterConst::PureGaugeUpdaterConst(const PureGaugeUpdaterConst& copy) : LatticeSweep(copy) {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif
#ifndef MULTITHREADING
PureGaugeUpdaterConst::PureGaugeUpdaterConst(const PureGaugeUpdaterConst& copy) : randomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif

#ifdef MULTITHREADING
PureGaugeUpdaterConst::~PureGaugeUpdaterConst() {
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomUniform[i];
	}
	delete[] randomGenerator;
	delete[] randomUniform;
}
#endif
#ifndef MULTITHREADING
PureGaugeUpdaterConst::~PureGaugeUpdaterConst() { }
#endif

void PureGaugeUpdaterConst::execute(environment_t & environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	real_t beta = environment.configurations.get<real_t>("beta");
  double delta = environment.configurations.get<real_t>("delta");
	//Get the gauge action
	GaugeAction* action = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<real_t>("beta"));
  //GaugeAction* action = GaugeAction::getInstance("StandardWilson",environment.configurations.get<real_t>("beta"));
  double actionE = action->energy(environment);
  //std::cout << "PureGaugeUpdaterConst did an Update " << std::endl;
	
#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif

#ifndef ENABLE_MPI
#ifdef MULTITHREADING
	for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#endif
#pragma omp parallel for
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,mu) == color) {
#endif
					this->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta, actionE,delta);
          //std::cout << "Always here? " << std::endl;
          //actionE=action->energy(environment);
          //std::cout << "Energy after Update on site and link: " << actionE << std::endl;
          //std::cout << "Energy after Update on site and link calc expl: " << action->energy(environment) << std::endl;
#ifdef MULTITHREADING
				}
#endif
      
			}
		}
   
		environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
	}
#endif
#endif

	//We suppose that always MPI+MTH
#ifdef ENABLE_MPI
	for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
		for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
			if (processor == Layout::this_processor) {
#pragma omp parallel for
				for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta, actionE,delta);
              std::cout << "Ever here? " << std::endl;
              //actionE=action->energy(environment);
              //std::cout << "Energy after Update on site and link: " << actionE << std::endl;
              //std::cout << "Energy after Update on site and link calc expl: " << action->energy(environment) << std::endl;
						}
              
					}
				}
			}
			else {
#pragma omp parallel for 
				for (int site = environment.gaugeLinkConfiguration.sharedsize; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta, actionE,delta);
              std::cout << "Ever here? " << std::endl;
              //actionE=action->energy(environment);
              //std::cout << "Energy after Update on site and link: " << actionE << std::endl;
              //std::cout << "Energy after Update on site and link calc expl: " << action->energy(environment) << std::endl;
						}
              
					}
				}
			}
		}
    
		environment.gaugeLinkConfiguration.updateHalo();
	}
  
#endif
  std::cout << "Energy after Update on site and link calc expl: " << action->energy(environment) << std::endl;
	environment.synchronize();
	delete action;
}



void PureGaugeUpdaterConst::executebeta(environment_t & environment,real_t beta, double Energylowerend) {
	typedef extended_gauge_lattice_t::Layout Layout;
  double delta = environment.configurations.get<real_t>("delta");
	//Get the gauge action
	GaugeAction* action = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<real_t>("beta"));
  double actionE = action->energy(environment);
	
#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif

#ifndef ENABLE_MPI
#ifdef MULTITHREADING
	for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#endif
#pragma omp parallel for
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,mu) == color) {
#endif
					this->updateLinkbeta(environment.gaugeLinkConfiguration, site, mu, action, beta,actionE,Energylowerend,delta);
          //std::cout << "Energy after Update on site and link: " << actionE << std::endl;
#ifdef MULTITHREADING
				}
#endif
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
	}
#endif
#endif

	//We suppose that always MPI+MTH
#ifdef ENABLE_MPI
	for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
		for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
			if (processor == Layout::this_processor) {
#pragma omp parallel for
				for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLinkbeta(environment.gaugeLinkConfiguration, site, mu, action, beta,actionE,Energylowerend,delta);
						}
					}
				}
			}
			else {
#pragma omp parallel for 
				for (int site = environment.gaugeLinkConfiguration.sharedsize; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLinkbeta(environment.gaugeLinkConfiguration, site, mu, action, beta,actionE,Energylowerend,delta);
						}
					}
				}
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
	}
#endif

	environment.synchronize();
	delete action;
}






#ifdef MULTITHREADING
real_t PureGaugeUpdaterConst::generate_radius(real_t b) {
	real_t x1 = log((*randomUniform[omp_get_thread_num()])()), x2 = log((*randomUniform[omp_get_thread_num()])()), x3 = pow(cos(2.*PI*(*randomUniform[omp_get_thread_num()])()), 2.);
	real_t s = 1. + b*(x1+x2*x3);
	real_t r = (*randomUniform[omp_get_thread_num()])();
	while ((1.+s-2.*r*r) < 0) {
		x1 = log((*randomUniform[omp_get_thread_num()])());
		x2 = log((*randomUniform[omp_get_thread_num()])());
		x3 = pow(cos(2.*PI*(*randomUniform[omp_get_thread_num()])()), 2.);
		s = 1. + b*(x1+x2*x3);
		r = (*randomUniform[omp_get_thread_num()])();
    //std::cout << "b = " << b << ", x1 = " << x1 << ", x2 = " << x2 << ", x3 = " << x3 << ", s = " << s << ", r = " << r << " --> " << (1.+s-2.*r*r) << std::endl;
	}
	return s;
}
#endif
#ifndef MULTITHREADING
real_t PureGaugeUpdaterConst::generate_radius(real_t b) {
	real_t x1 = log(randomUniform()), x2 = log(randomUniform()), x3 = pow(cos(2.*PI*randomUniform()), 2.);
	real_t s = 1. + b*(x1+x2*x3);
	real_t r = randomUniform();
	while ((1.+s-2.*r*r) < 0) {
		x1 = log(randomUniform());
		x2 = log(randomUniform());
		x3 = pow(cos(2.*PI*randomUniform()), 2.);
		s = 1. + b*(x1+x2*x3);
		r = randomUniform();
	}
	return s;
}
#endif

#ifdef MULTITHREADING
void PureGaugeUpdaterConst::generate_vector(real_t radius, real_t& u1, real_t& u2, real_t& u3) {
	real_t phi = 2.*PI*(*randomUniform[omp_get_thread_num()])();
	real_t theta = acos(2.*(*randomUniform[omp_get_thread_num()])()-1.);
	u1 = radius*sin(theta)*cos(phi);
	u2 = radius*sin(theta)*sin(phi);
	u3 = radius*cos(theta);
}
#endif
#ifndef MULTITHREADING
void PureGaugeUpdaterConst::generate_vector(real_t radius, real_t& u1, real_t& u2, real_t& u3) {
	real_t phi = 2.*PI*randomUniform();
	real_t theta = acos(2.*randomUniform()-1.);
	u1 = radius*sin(theta)*cos(phi);
	u2 = radius*sin(theta)*sin(phi);
	u3 = radius*cos(theta);
}
#endif


void PureGaugeUpdaterConst::updateLink(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* action, double beta, double& actionE, double delta) {
	GaugeGroup staple = action->staple(lattice, site, mu);
  
  double Energy = -3000.00;//-65000.0;
  
  
  //std::cout << "Energy: " << actionE << " lower border: " << Energy << " higher border: " << Energy+delta << std::endl;
  
#if NUMCOLORS > 2
	//take the plaquette
  
	GaugeGroup plaquette = lattice[site][mu]*(staple);
  GaugeGroup trialm = lattice[site][mu];
  GaugeGroup plaquettetrial = lattice[site][mu]*(staple);
  
  double difference;
  //double actionE = action->energy(environment);
  //double differencesum;
	for (unsigned int k = 0; k < numberColors-1; ++k) {
		for (int l = k+1; l < numberColors; ++l) {
			//Take the su2 subgroup matrix for the Cabibbo Marinari update
      
			real_t aeff = (imag(plaquette.at(k,k))-imag(plaquette.at(l,l)))/2.;
			real_t beff = (imag(plaquette.at(k,l))+imag(plaquette.at(l,k)))/2.;
			real_t ceff = (real(plaquette.at(k,l))-real(plaquette.at(l,k)))/2.;
			real_t deff = (real(plaquette.at(k,k))+real(plaquette.at(l,l)))/2.;
			real_t detStaple = aeff*aeff + beff*beff + ceff*ceff + deff*deff;
			real_t b = numberColors/(2.*beta*sqrt(detStaple));
			//Use the standard Kennedy-Pendleton algorithm for su2
			real_t u0, u1, u2, u3;
			u0 = generate_radius(b);
			generate_vector(sqrt(1.-u0*u0), u1, u2, u3);
			//Calculate the su2 update matrix
			matrix2x2_t subupdate;
			subupdate.at(0,0) = std::complex<real_t>(u0, u3);
			subupdate.at(0,1) = std::complex<real_t>(u2, u1);
			subupdate.at(1,0) = std::complex<real_t>(-u2, u1);
			subupdate.at(1,1) = std::complex<real_t>(u0, -u3);
			matrix2x2_t substaple;
			substaple.at(0,0) = std::complex<real_t>(deff, aeff);
			substaple.at(0,1) = std::complex<real_t>(ceff, beff);
			substaple.at(1,0) = std::complex<real_t>(-ceff, beff);
			substaple.at(1,1) = std::complex<real_t>(deff, -aeff);
			//extract the su2 update matrix, compute also one overrelaxation step
			subupdate = htrans(substaple)*htrans((subupdate)*(htrans(substaple)))*htrans(substaple)/pow(detStaple,3./2.);
			//Update the plaquette and the link
      std::complex<real_t> tmp3;
      std::complex<real_t> tmp4;
			for (int i = 0; i < numberColors; ++i) {
				std::complex<real_t> tmp1 = subupdate.at(0,0)*lattice[site][mu].at(k,i) + subupdate.at(0,1)*lattice[site][mu].at(l,i);
				std::complex<real_t> tmp2 = subupdate.at(1,0)*lattice[site][mu].at(k,i) + subupdate.at(1,1)*lattice[site][mu].at(l,i);
        //std::complex<real_t> save1 = lattice[site][mu].at(k,i);
        //std::complex<real_t> save2 = lattice[site][mu].at(l,i);
				//lattice[site][mu].at(k,i) = tmp1;
				//lattice[site][mu].at(l,i) = tmp2;
        trialm.at(k,i) = tmp1;
				trialm.at(l,i) = tmp2;
				tmp3 = subupdate.at(0,0)*plaquette.at(k,i) + subupdate.at(0,1)*plaquette.at(l,i);
				tmp4 = subupdate.at(1,0)*plaquette.at(k,i) + subupdate.at(1,1)*plaquette.at(l,i);
        //std::complex<real_t> save3 = plaquette.at(k,i);
        //std::complex<real_t> save4 = plaquette.at(l,i);
				plaquettetrial.at(k,i) = tmp3;
				plaquettetrial.at(l,i) = tmp4;
        //difference = action->deltaAction(lattice, trialm, staple, site, mu);
    
			}
      /*
      if ((actionE+difference) <= (Energy+delta) && (actionE+difference)>Energy)
      {
        actionE = actionE+difference;
        //std::cout << "accepted Link Update"<<  std::endl;
        for (int i = 0; i < numberColors; ++i) {
          lattice[site][mu].at(k,i) = trialm.at(k,i);
				  lattice[site][mu].at(l,i) = trialm.at(l,i);  
                                      
          std::complex<real_t> tmp3 = subupdate.at(0,0)*plaquette.at(k,i) + subupdate.at(0,1)*plaquette.at(l,i);
				  std::complex<real_t> tmp4 = subupdate.at(1,0)*plaquette.at(k,i) + subupdate.at(1,1)*plaquette.at(l,i);
				  plaquette.at(k,i) = tmp3;
				  plaquette.at(l,i) = tmp4;
        }
        
      }
      
      else
      {
        //std::cout << "rejected Link Update " << std::endl;
        for (int i = 0; i < numberColors; ++i) {
          trialm.at(k,i)=lattice[site][mu].at(k,i);
			    trialm.at(l,i)=lattice[site][mu].at(l,i);  
          
        }
        
      }
      */
		}
   
	}
 
  difference = action->deltaAction(lattice, trialm, staple, site, mu);
  if ((actionE+difference) <= (Energy+delta) && (actionE+difference)>Energy)
  {
    actionE = actionE+difference;
    //std::cout << "accepted Link Update"<<  std::endl;
    
    for (unsigned int k = 0; k < numberColors-1; ++k) {
      for (int l = k+1; l < numberColors; ++l) {
        for (int i = 0; i < numberColors; ++i) {
          
          
          lattice[site][mu].at(k,i) = trialm.at(k,i);
 	        lattice[site][mu].at(l,i) = trialm.at(l,i);  
                                      
 	        plaquette.at(k,i) = plaquettetrial.at(k,i);
 	        plaquette.at(l,i) = plaquettetrial.at(l,i);
          
        
        }
      }
    }
  }
      
  else
  {
    //std::cout << "rejected Link Update " << std::endl;
    for (unsigned int k = 0; k < numberColors-1; ++k) {
      for (int l = k+1; l < numberColors; ++l) {
        for (int i = 0; i < numberColors; ++i) {      
          trialm.at(k,i)=lattice[site][mu].at(k,i);
		      trialm.at(l,i)=lattice[site][mu].at(l,i);
          
          plaquettetrial.at(k,i) = plaquette.at(k,i);
 	        plaquettetrial.at(l,i) = plaquette.at(l,i); 
        }
      }
    }
  } 
   
   
   
   
#endif
#if NUMCOLORS == 2
	//compute the parameters for the kennedy algorithm
	real_t detStaple = abs(det(staple));
	real_t b = 1./(beta*sqrt(detStaple));
	real_t u0, u1, u2, u3;
	u0 = generate_radius(b);
	generate_vector(sqrt(1.-u0*u0), u1, u2, u3);
	//update the matrix
	lattice[site][mu].at(0,0) = std::complex<real_t>(u0, u3);
	lattice[site][mu].at(0,1) = std::complex<real_t>(u2, u1);
	lattice[site][mu].at(1,0) = std::complex<real_t>(-u2, u1);
	lattice[site][mu].at(1,1) = std::complex<real_t>(u0, -u3);
	//redefine the matrix U_New = R S^\dag, normalize the determinant
	lattice[site][mu] = (lattice[site][mu])*(htrans(staple))/sqrt(detStaple);
	//perform a single overrelaxation step
	lattice[site][mu] = (htrans(staple)*htrans(lattice[site][mu])*htrans(staple))/detStaple;
#endif
 // std::cout << "Energy after Update: " << actionE << std::endl;
}



void PureGaugeUpdaterConst::updateLinkbeta(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* action, double beta, double& actionE, double Energy, double delta) {
	GaugeGroup staple = action->staple(lattice, site, mu);
  
  //double Energy = -3000.00;//-65000.0;
  
  
  //std::cout << "Energy: " << actionE << " lower border: " << Energy << " higher border: " << Energy+delta << std::endl;
  
#if NUMCOLORS > 2
  //std::cout << "calc plaquette and staple"<< std::endl;
	//take the plaquette
	GaugeGroup plaquette = lattice[site][mu]*(staple);
  GaugeGroup trialm = lattice[site][mu];
  GaugeGroup latticeold = lattice[site][mu];
  GaugeGroup plaquettetrial = lattice[site][mu]*(staple);
  double difference;
  //double actionE = action->energy(environment);
  //double differencesum;
  //std::cout << "calc subupdate"<< std::endl;
	for (unsigned int k = 0; k < numberColors-1; ++k) {
		for (int l = k+1; l < numberColors; ++l) {
			//Take the su2 subgroup matrix for the Cabibbo Marinari update
      
			real_t aeff = (imag(plaquette.at(k,k))-imag(plaquette.at(l,l)))/2.;
			real_t beff = (imag(plaquette.at(k,l))+imag(plaquette.at(l,k)))/2.;
			real_t ceff = (real(plaquette.at(k,l))-real(plaquette.at(l,k)))/2.;
			real_t deff = (real(plaquette.at(k,k))+real(plaquette.at(l,l)))/2.;
			real_t detStaple = aeff*aeff + beff*beff + ceff*ceff + deff*deff;
			real_t b = numberColors/(2.*beta*sqrt(detStaple));
			//Use the standard Kennedy-Pendleton algorithm for su2
			real_t u0, u1, u2, u3;
      //std::cout << "generate radius" << std::endl;
			u0 = generate_radius(b);
      //std::cout << "generate vector" << std::endl;
			generate_vector(sqrt(1.-u0*u0), u1, u2, u3);
			//Calculate the su2 update matrix
			matrix2x2_t subupdate;
			subupdate.at(0,0) = std::complex<real_t>(u0, u3);
			subupdate.at(0,1) = std::complex<real_t>(u2, u1);
			subupdate.at(1,0) = std::complex<real_t>(-u2, u1);
			subupdate.at(1,1) = std::complex<real_t>(u0, -u3);
			matrix2x2_t substaple;
			substaple.at(0,0) = std::complex<real_t>(deff, aeff);
			substaple.at(0,1) = std::complex<real_t>(ceff, beff);
			substaple.at(1,0) = std::complex<real_t>(-ceff, beff);
			substaple.at(1,1) = std::complex<real_t>(deff, -aeff);
			//extract the su2 update matrix, compute also one overrelaxation step
			subupdate = htrans(substaple)*htrans((subupdate)*(htrans(substaple)))*htrans(substaple)/pow(detStaple,3./2.);
			//Update the plaquette and the link
      std::complex<real_t> tmp3;
      std::complex<real_t> tmp4;
			for (int i = 0; i < numberColors; ++i) {
				//std::complex<real_t> tmp1 = subupdate.at(0,0)*lattice[site][mu].at(k,i) + subupdate.at(0,1)*lattice[site][mu].at(l,i);
				//std::complex<real_t> tmp2 = subupdate.at(1,0)*lattice[site][mu].at(k,i) + subupdate.at(1,1)*lattice[site][mu].at(l,i);
        std::complex<real_t> tmp1 = subupdate.at(0,0)*trialm.at(k,i) + subupdate.at(0,1)*trialm.at(l,i);
				std::complex<real_t> tmp2 = subupdate.at(1,0)*trialm.at(k,i) + subupdate.at(1,1)*trialm.at(l,i);
        //std::complex<real_t> save1 = lattice[site][mu].at(k,i);
        //std::complex<real_t> save2 = lattice[site][mu].at(l,i);
				//lattice[site][mu].at(k,i) = tmp1;
				//lattice[site][mu].at(l,i) = tmp2;
        trialm.at(k,i) = tmp1;
				trialm.at(l,i) = tmp2;
				tmp3 = subupdate.at(0,0)*plaquette.at(k,i) + subupdate.at(0,1)*plaquette.at(l,i);
				tmp4 = subupdate.at(1,0)*plaquette.at(k,i) + subupdate.at(1,1)*plaquette.at(l,i);
        //std::complex<real_t> save3 = plaquette.at(k,i);
        //std::complex<real_t> save4 = plaquette.at(l,i);
				plaquettetrial.at(k,i) = tmp3;
				plaquettetrial.at(l,i) = tmp4;
        plaquette.at(k,i) = tmp3;
				plaquette.at(l,i) = tmp4;
        //difference = action->deltaAction(lattice, trialm, staple, site, mu);
    
			}
      /*
      if ((actionE+difference) <= (Energy+delta) && (actionE+difference)>Energy)
      {
        actionE = actionE+difference;
        //std::cout << "accepted Link Update"<<  std::endl;
        for (int i = 0; i < numberColors; ++i) {
          lattice[site][mu].at(k,i) = trialm.at(k,i);
				  lattice[site][mu].at(l,i) = trialm.at(l,i);  
                                      
          std::complex<real_t> tmp3 = subupdate.at(0,0)*plaquette.at(k,i) + subupdate.at(0,1)*plaquette.at(l,i);
				  std::complex<real_t> tmp4 = subupdate.at(1,0)*plaquette.at(k,i) + subupdate.at(1,1)*plaquette.at(l,i);
				  plaquette.at(k,i) = tmp3;
				  plaquette.at(l,i) = tmp4;
        }
        
      }
      
      else
      {
        //std::cout << "rejected Link Update " << std::endl;
        for (int i = 0; i < numberColors; ++i) {
          trialm.at(k,i)=lattice[site][mu].at(k,i);
			    trialm.at(l,i)=lattice[site][mu].at(l,i);  
          
        }
        
      }
      */
		}
   
	}
  //std::cout << "calc difference and change site and plaquette"<< std::endl;
  difference = action->deltaAction(lattice, trialm, staple, site, mu);
  if ((actionE+difference) <= (Energy+delta) && (actionE+difference)>Energy)
  {
    actionE = actionE+difference;
    //std::cout << "accepted Link Update"<<  std::endl;
    
    for (unsigned int k = 0; k < numberColors-1; ++k) {
      for (int l = k+1; l < numberColors; ++l) {
        for (int i = 0; i < numberColors; ++i) {
          
          
          lattice[site][mu].at(k,i) = trialm.at(k,i);
 	        lattice[site][mu].at(l,i) = trialm.at(l,i);  
                                      
 	        plaquette.at(k,i) = plaquettetrial.at(k,i);
 	        plaquette.at(l,i) = plaquettetrial.at(l,i);
          
        
        }
      }
    }
  }
      
  else
  {
    //std::cout << "rejected Link Update " << std::endl;
    for (unsigned int k = 0; k < numberColors-1; ++k) {
      for (int l = k+1; l < numberColors; ++l) {
        for (int i = 0; i < numberColors; ++i) {      
          trialm.at(k,i)=lattice[site][mu].at(k,i);
		      trialm.at(l,i)=lattice[site][mu].at(l,i);
          
          plaquettetrial.at(k,i) = plaquette.at(k,i);
 	        plaquettetrial.at(l,i) = plaquette.at(l,i); 
        }
      }
    }
  } 
   
   
   
   
#endif
#if NUMCOLORS == 2
	//compute the parameters for the kennedy algorithm
	real_t detStaple = abs(det(staple));
	real_t b = 1./(beta*sqrt(detStaple));
	real_t u0, u1, u2, u3;
	u0 = generate_radius(b);
	generate_vector(sqrt(1.-u0*u0), u1, u2, u3);
	//update the matrix
	lattice[site][mu].at(0,0) = std::complex<real_t>(u0, u3);
	lattice[site][mu].at(0,1) = std::complex<real_t>(u2, u1);
	lattice[site][mu].at(1,0) = std::complex<real_t>(-u2, u1);
	lattice[site][mu].at(1,1) = std::complex<real_t>(u0, -u3);
	//redefine the matrix U_New = R S^\dag, normalize the determinant
	lattice[site][mu] = (lattice[site][mu])*(htrans(staple))/sqrt(detStaple);
	//perform a single overrelaxation step
	lattice[site][mu] = (htrans(staple)*htrans(lattice[site][mu])*htrans(staple))/detStaple;
#endif
 // std::cout << "Energy after Update: " << actionE << std::endl;
}


} /* namespace Update */
