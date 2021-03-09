#include "Environment.h"
#include "io/StorageParameters.h"
#include "MatrixTypedef.h"
#include "utils/RandomSeed.h"
#include "utils/ToString.h"
#include "Simulation.h"
#include "LatticeSweep.h"
#include "io/GlobalOutput.h"
#include "MPILattice/ReducedStencil.h"
#include "MPILattice/StandardStencil.h"
#include "MPILattice/ExtendedStencil.h"
#include "MPILattice/LocalLayout.h"
#include "utils/LieGenerators.h"
#include "utils/ToString.h"
#include <iostream>
#include <fenv.h>

//#include "actions/GaugeEnergy.h"
#include "actions/GaugeAction.h"
#include <fstream>
//added by Felix
/*
#include "actions/GaugeEnergy.h"
#include "actions/GaugeAction.h"
void findEintsmooth(double Eint, double delta);
*/
//up to here

int Update::RandomSeed::counter = -1;
boost::mt19937 Update::RandomSeed::rng;
boost::uniform_int<> Update::RandomSeed::dist = boost::uniform_int<>(-10000000,10000000);

//MPI Datatype initilialization
#ifdef ENABLE_MPI
MPI_Datatype MpiType<short int>::type = MPI_SHORT;
MPI_Datatype MpiType<int>::type = MPI_INT;
MPI_Datatype MpiType<int[4]>::type = MPI_INT;
MPI_Datatype MpiType<Update::real_t>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::real_t[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup[6]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup[6]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalVector[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointVector[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup>::type = MPI_DOUBLE;
#ifdef ADJOINT
MPI_Datatype MpiType<Update::FermionicForceMatrix[4]>::type = MPI_DOUBLE;
#endif
#endif


namespace po = boost::program_options;


int main(int ac, char* av[]) {
#ifdef ENABLE_MPI
	//Initialize MPI if in MPI mode
	MPI_Init(&ac, &av);
#endif

	//The descriptor for the global options
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("list_sweeps", "produce the list of all the available sweeps")
		("configfile", po::value<std::string>(), "set the configuration file")
		("version", "Version message")

		//Action specifications
		("name_action", po::value<std::string>(), "the name of the gauge part of the action (\"StandardWilson/Improved\")")
		("dirac_operator", po::value<std::string>(), "the name of the dirac wilson operator that should be used (supported: DiracWilson/Improved/Overlap/ExactOverlap)")
		
		//Action and dirac operator options
		("beta", po::value<Update::real_t>(), "set the \\beta parameter of the simulations")
   
    ("Emin", po::value<Update::real_t>(), "set the min Energy of LLR")
    ("Emax", po::value<Update::real_t>(), "set the max Energy of LLR")
    ("delta", po::value<Update::real_t>(), "set the Energyinterval of LLR")
    ("ait", po::value<unsigned int>(), "set the number of iterations to find a in LLR")
    ("Njacknife", po::value<unsigned int>(), "set the number of Jackknife samples for LLR")
    
    
		("kappa", po::value<Update::real_t>(), "set the \\kappa parameter of the simulations")
		("mass", po::value<Update::real_t>(), "set the mass parameter of the overlap operator")
		("csw", po::value<Update::real_t>(), "The clover term coefficient")
		("stout_smearing_levels", po::value<int>(), "The levels for the stout smearing of the dirac operator")
		("stout_smearing_rho", po::value<Update::real_t>(), "The rho for the stout smearing of the dirac operator")
		
		//GRID parallelization and lattice options
		("glob_x", po::value<unsigned int>(), "The x lattice size")
		("glob_y", po::value<unsigned int>(), "The y lattice size")
		("glob_z", po::value<unsigned int>(), "The z lattice size")
		("glob_t", po::value<unsigned int>(), "The t lattice size")
		("pgrid_x", po::value<unsigned int>(), "The grid subdivision in the x direction")
		("pgrid_y", po::value<unsigned int>(), "The grid subdivision in the y direction")
		("pgrid_z", po::value<unsigned int>(), "The grid subdivision in the z direction")
		("pgrid_t", po::value<unsigned int>(), "The grid subdivision in the t direction")
		("number_threads", po::value<unsigned int>(), "The number of threads for openmp")
		("load_layout", "If the MPI layout should be loaded from the disk")
		("print_report_layout", "If the full report of the MPI layout should be printed")

		//Boundary conditions
		("boundary_conditions", po::value<std::string>(), "Boundary conditions to use: periodic (fermions), antiperiodic (fermions), spatialantiperiodic (fermion), open")

		//Input-output options
		("output_directory_configurations", po::value<std::string>(), "The directory for the output of the configurations")
		("output_directory_measurements", po::value<std::string>(), "The directory for the output of the measurements")
		("output_name", po::value<std::string>(), "The name for the beginning part of the file of the output")
		("input_directory_configurations", po::value<std::string>(), "The directory for the input of the configurations")
		("input_name", po::value<std::string>(), "The name for the beginning part of the file of the input")
		("input_number", po::value<unsigned int>(), "The number of the file of the input")
		("output_configuration_name", po::value<std::string>(), "The name of the output of the field")
		("output_offset", po::value<unsigned int>(), "The offset for the number of the output configurations")
		("format_name", po::value<std::string>(), "leonard_format/muenster_format for reading and writing configurations")
		("input_format_name", po::value<std::string>(), "leonard_format/muenster_format only for reading configurations")
		("output_format_name", po::value<std::string>(), "leonard_format/muenster_format only for writing configurations")
		("measurement_output_format", po::value<std::string>()->default_value("txt"), "output format for the measurements (xml/txt)")
		
		//Start, warm up and measurement specifications
		("start", po::value<std::string>(), "the start gauge configuration for the simulations (hotstart/coldstart/readstart)")
		("start_gauge_configuration_file", po::value<std::string>(), "The name of the beginning file with the gauge configuration for reading")
		("start_configuration_number", po::value<unsigned int>(), "The beginning number of the output configuration written (zero as default)")
		("number_warm_up_sweeps", po::value<unsigned int>(), "the number of warm-up sweeps")
		("number_measurement_sweeps", po::value<unsigned int>(), "the number of measurement sweeps")
		("warm_up_sweeps", po::value< std::string >(), "the vector of the warm up sweeps to do (example: {{PureGaugeCM,1,1},{Plaquette,1,1}} )")
		("measurement_sweeps", po::value< std::string >(), "the vector of the measurement sweeps to do (example: {{PureGaugeCM,1,1},{Plaquette,1,1}} )")
		
		//Scalar field options
		("adjoint_nf_scalars", po::value<unsigned int>()->default_value(0), "set the number of the adjoint scalar fields")
		("fundamental_nf_scalars", po::value<unsigned int>()->default_value(0), "set the number of the fundamental scalar fields")
		
		//HMC options
		("name_integrator", po::value<std::string>(), "the name of the type of integrator (\"second_order, omelyan, fourth_order, fourth_omelyan\")")
		("hmc_t_length", po::value<Update::real_t>(), "the length of a single HMC step (examples: 0.1, 0.05 ...)")
		("number_hmc_steps", po::value<std::string>(), "the vector of the numbers of HMC steps for a single trajectory (examples: 2, 7 ...)")
		
		//RHMC options
		("force_inverter_precision", po::value<Update::real_t>(), "The precision for the inverter in the force step")
		("metropolis_inverter_precision", po::value<Update::real_t>(), "The precision for the inverter in the metropolis step")
		("metropolis_inverter_max_steps", po::value<unsigned int>(),"maximum level of steps used by the inverters for computing the energy of the metropolis step")
		("force_inverter_max_steps", po::value<unsigned int>(),"maximum level of steps used by the inverters for computing the force")
		("number_pseudofermions", po::value<unsigned int>(), "the number of pseudofermions used")
		("number_force_levels", po::value<unsigned int>(), "the number of levels used for the force")
		("check_rational_approximations", po::value< std::string >(), "Set to true for checking the approximations in the beginning of the simulation")
		("twisted_mass_squared", po::value<Update::real_t>(), "The twisted mass squared used by twisted updater to regularize to RHMC")
		("number_twisted_correction_noise_vectors", po::value<unsigned int>(), "Number of noise vectors used to estimate the correction factor")
		("number_block_correction_noise_vectors", po::value<unsigned int>(), "Number of noise vectors used to estimate the correction factor")
		("theory_power_factor", po::value<Update::real_t>(), "The half of the number of fermion (1/4 for SUSY, 1 for two flavor, etc)")
		("twisted_breakup", po::value<unsigned int>(), "The value of the determinant breakup used in the correction step")
		("deflation_block_size", po::value<std::string>(), "The vector of the block size {lx,ly,lz,lt}")
		("preconditioner_recursions", po::value<unsigned int>(), "Number of preconditioner recursions used for evaluating the rational approximations")
		("preconditioner_precision", po::value<Update::real_t>(), "The precision for the inversion of the preconditioner used for evaluating the rational approximations")

		//Options for the ReadGaugeConfiguration sweep
		("read_start_number", po::value<unsigned int>(), "From which configurations start to read again for the analysis")
		("read_step", po::value<unsigned int>(), "The step in the reading analysis")

		//Overlap operator options
		("OverlapOperator::squareRootApproximation", po::value<std::string>(), "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("ExactOverlapOperator::squareRootApproximation", po::value<std::string>()->default_value(""), "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("ExactOverlapOperator::eigensolver::use_chebyshev", po::value<std::string>()->default_value("false"), "Use Chebyshev acceleration? (true/false)")
		("ExactOverlapOperator::eigensolver::chebyshev_left", po::value<double>()->default_value(0.2), "Left interval of the Chebyshev polynomial.")
		("ExactOverlapOperator::eigensolver::chebyshev_right", po::value<double>()->default_value(7.), "Right interval of the Chebyshev polynomial.")
		("ExactOverlapOperator::eigensolver::chebyshev_order", po::value<unsigned int>()->default_value(15), "Order of the Chebyshev acceleration. It must be an odd number")
		("ExactOverlapOperator::eigensolver::eigensolver_precision", po::value<double>()->default_value(0.000000001), "set the precision used by the eigensolver")
		("ExactOverlapOperator::eigensolver::number_extra_vectors", po::value<unsigned int>(), "Number of extra vectors for the Arnoldi algorithm used in the computation of the eigenvectors, increase this number to increase precision")
		("ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver", po::value<unsigned int>()->default_value(50), "Number of restarts for the implicitly restarted Arnoldi algorithm")
		("ExactOverlapOperator::eigensolver::number_eigenvalues", po::value<unsigned int>(), "Number of eigenvalues of the dirac wilson operator to be computed")

	    	//Options for the inverter
		("generic_inverter_max_steps", po::value<unsigned int>(),"maximum level of steps used by the inverters")
		("generic_inverter_precision", po::value<Update::real_t>(), "The precision for the inverter")
	;

	for (int level = 1; level < 4; ++level) {
		for (int i = 1; i < 33; ++i) {
			desc.add_options()((std::string("force_rational_fraction_")+Update::toString(i)+"_level_"+Update::toString(level)).c_str(),po::value<std::string>(),"Rational approximation of the force used by MultiStep integrator");
		}
		desc.add_options()(std::string("force_inverter_precision_level_"+Update::toString(level)).c_str(),po::value<Update::real_t>(),"Precision of the inverter of the force used by MultiStep integrator");
	}
	for (int i = 1; i < 33; ++i) {
		desc.add_options()((std::string("heatbath_rational_fraction_")+Update::toString(i)).c_str(),po::value<std::string>(),"the rational fraction approximation that should be used for the heatbath (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
		desc.add_options()((std::string("metropolis_rational_fraction_")+Update::toString(i)).c_str(),po::value<std::string>(),"the rational fraction approximation that should be used for the metropolis (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
		desc.add_options()((std::string("force_rational_fraction_")+Update::toString(i)).c_str(),po::value<std::string>(),"the rational fraction approximation that should be used for the force (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
	}

	Update::LatticeSweep::addParameters(desc);

	//Now we can parse the command line
	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	//If there is a configfile we parse also it
	if (vm.count("configfile")) {
		std::ifstream in(vm["configfile"].as<std::string>().c_str());
		po::store(po::parse_config_file(in, desc, true), vm);
		po::notify(vm);
	} else if (vm.count("help")) {
		//Print the help description of the parameters of the program
		std::cout << desc << std::endl;
		std::cout << "type --list_sweeps for having a list of all the possible sweeps" << std::endl;
		return 0;
	} else if (vm.count("version")) {
		std::cout << "LeonardQCD version: 1.0" << std::endl;
#ifdef ADJOINT
		std::cout << " with adjoint fermions" << std::endl;
#endif
#ifndef ADJOINT
		std::cout << " with fundamental fermions" << std::endl;
#endif
		std::cout << " with number of colors: " << Update::numberColors << std::endl;
#ifdef ENABLE_MPI
		std::cout << " with mpi enabled" << std::endl;
#endif
		return 0;
	} else if (vm.count("list_sweeps")) {
		Update::LatticeSweep::printSweepsName();
		return 0;
	}

	//Initialized up/down stencil
	Lattice::StandardStencil::initializeNeighbourSites();
	Lattice::ExtendedStencil::initializeNeighbourSites();
	Lattice::ReducedStencil::initializeNeighbourSites();

	//Initialize lattice layout
#ifndef ENABLE_MPI
	Lattice::LocalLayout::pgrid_t = 1;
	Lattice::LocalLayout::pgrid_x = 1;
	Lattice::LocalLayout::pgrid_y = 1;
	Lattice::LocalLayout::pgrid_z = 1;

	Lattice::LocalLayout::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::LocalLayout::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::LocalLayout::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::LocalLayout::glob_z = vm["glob_z"].as<unsigned int>();

	Lattice::LocalLayout::initialize();

	if (isOutputProcess()) std::cout << "Lattice size (x,y,z,t): (" << Lattice::LocalLayout::glob_x << "," << Lattice::LocalLayout::glob_y << "," << Lattice::LocalLayout::glob_z << "," << Lattice::LocalLayout::glob_t << ")" << std::endl;
#endif
#ifdef ENABLE_MPI
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_t = vm["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_x = vm["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_y = vm["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_z = vm["pgrid_z"].as<unsigned int>();
	
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_t = vm["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_x = vm["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_y = vm["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_z = vm["pgrid_z"].as<unsigned int>();

	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_t = vm["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_x = vm["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_y = vm["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_z = vm["pgrid_z"].as<unsigned int>();


	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_z = vm["glob_z"].as<unsigned int>();
	
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_z = vm["glob_z"].as<unsigned int>();

	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_z = vm["glob_z"].as<unsigned int>();

	if (vm.count("load_layout")) {
		std::string basename = vm["output_directory_configurations"].as<std::string>();
		Lattice::MpiLayout<Lattice::ExtendedStencil>::load(basename+"extended_layout");
		Lattice::MpiLayout<Lattice::StandardStencil>::load(basename+"standard_layout");
		Lattice::MpiLayout<Lattice::ReducedStencil>::load(basename+"reduced_layout");
	} else {
		Lattice::MpiLayout<Lattice::ExtendedStencil>::initialize();
		Lattice::MpiLayout<Lattice::StandardStencil>::initialize();
		Lattice::MpiLayout<Lattice::ReducedStencil>::initialize();
		
		std::string basename = vm["output_directory_configurations"].as<std::string>();
		Lattice::MpiLayout<Lattice::ExtendedStencil>::save(basename+"extended_layout");
		Lattice::MpiLayout<Lattice::StandardStencil>::save(basename+"standard_layout");
		Lattice::MpiLayout<Lattice::ReducedStencil>::save(basename+"reduced_layout");
	}

	if (vm.count("print_report_layout")) {
		if (isOutputProcess()) Lattice::MpiLayout<Lattice::ExtendedStencil>::printReport();
		if (isOutputProcess()) Lattice::MpiLayout<Lattice::ReducedStencil>::printReport();
	}

	if (isOutputProcess()) std::cout << "Lattice size (x,y,z,t): (" << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_x << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_y << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_z << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_t << ")" << std::endl;
	if (isOutputProcess()) std::cout << "Mpi grid (px,py,pz,pt): (" << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_x << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_y << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_z << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_t << ")" << std::endl;
#endif

	//Initialize the enviroment of the program
	Update::environment_t* environment = new Update::environment_t(vm);

	//Set the number of threads
	if (vm.count("number_threads")) {
#ifdef MULTITHREADING
		omp_set_num_threads(vm["number_threads"].as<unsigned int>());
		int num_threads = omp_get_max_threads();
		if (isOutputProcess()) std::cout << "Number of threads: " << num_threads << std::endl;
#endif
	}

	//Set the output to format
	Update::GlobalOutput* output = Update::GlobalOutput::getInstance();
	output->setFormat(vm["measurement_output_format"].as<std::string>());

	//Finally create the simulation
	Update::Simulation simulation(*environment);
	//Start and run warmup
	simulation.starterLLR(-660.0);
	simulation.warmUpLLR(1.0,-660.0);
	//Perform the measurements
	simulation.measurementLLR(1.0,-660.0);
   
  //double betaread = environment->configurations.get<std::string>("beta");
  std::cout << "See if this works: " << environment->configurations.get<Update::real_t>("delta") << std::endl;
  output->print();
	output->destroy();

	//destroy the environment
	delete environment;
  
  environment = new Update::environment_t(vm);
  
  double Emax = environment->configurations.get<Update::real_t>("Emax");
  double Emin = environment->configurations.get<Update::real_t>("Emin");
  double delta = environment->configurations.get<Update::real_t>("delta");
  int ait = environment->configurations.get<unsigned int>("ait");
  int Njacknife = environment->configurations.get<unsigned int>("Njacknife");
  int trajecs = environment->configurations.get<unsigned int>("number_measurement_sweeps");
  
  
  double x0[(int)((Emax-Emin)/delta)];  //lower end of energy interval
	double a[(int)((Emax-Emin)/delta)]; 
	double a_i[Njacknife]; 
	double a_i_new;
	int RobMarchcount; //number of iteration of the rob march algorithm
	double Reweightexpect; //Reweighted expectationvalue of the energy 
	double meassurement[trajecs]; //smh does not work with meassurement[trajecs]
	//double varianz;
  //bool Einterval = false; 
  
  //FILE *fpavalue;           
                  
  int Eint = 0;
  int jcount;
  int acounter;
  int k;
  
  //std::ofstream avaluesoutput ("io_directory/avalues.txt", std::ofstream::out);
  std::ofstream ofs ("io_directory/avalues.txt", std::ofstream::out);

  
/*
  std::ifstream infile("io_directory/energyvalues.txt");
  k=0;
  std::string stuff;
  while (std::getline (infile, stuff)) {
    // Output the text from the file
    std::cout << stuff << std::endl;
    meassurement[k]= std::stod(stuff);
    std::cout << meassurement[k] << std::endl;
    k = k + 1;
  }
  std:: cout <<"k: "<< k << std::endl;
  */
/*
        while (infile >> stuff)
        {
          std::cout <<  stuff << std::endl;
          //Reweightexpect = Reweightexpect + meassurement[k];
          k = k+1;
        }
        std:: cout <<"k: "<< k << std::endl;
        */
  //meassurement = std::malloc(trajecs * sizeof(double));
  /*
  for(Eint=0;Eint<=((Emax-Emin)/delta);Eint++)
  {
    x0[Eint] = (double)(Eint)*delta + Emin;
    a[Eint] = 0.0;
    for(jcount = 0;jcount<Njacknife;jcount++)
    {
      a_i[jcount] = 2.0;
		  a_i_new = 2.0;
        
      RobMarchcount = 0;
      environment = new Update::environment_t(vm);
      simulation.starterLLR(x0[Eint]);
      
      
      for(acounter=0;acounter<ait;acounter++)
      {
        a_i[jcount] = a_i_new;
        
        
        
        simulation.warmUpLLR(1.0,x0[Eint]);
          
        Reweightexpect=0;
				//varianz=0;
        
        simulation.measurementLLR(a_i[jcount],x0[Eint]);
        
        
        std::ifstream infile("io_directory/energyvalues.txt");
        k=0;
        std::string stuff;
        while (std::getline (infile, stuff)) {
        
        std::cout << stuff << std::endl;
        meassurement[k]= std::stod(stuff);
        Reweightexpect = Reweightexpect + meassurement[k];
        std::cout << meassurement[k] << std::endl;
        k = k + 1;
        }
        std:: cout <<"k: "<< k << std::endl;
        
        
        //for(k = 0;k<trajecs;k++)
				//{
          
					//meassurement[k] = action();
          //node0_printf("action = %.4g off \n", meassurement[k]);
					//Reweightexpect = Reweightexpect + meassurement[k];
          //updateconst_e(x0[Eint],delta,a_i[jcount]);
				//}
        Reweightexpect = Reweightexpect/trajecs;
        
        
        Reweightexpect = Reweightexpect - x0[Eint] - 0.5*delta;
        
        //node0_printf("Reweightexpect = %.4g off \n", Reweightexpect);
        
        if(RobMarchcount<80)
				{
					a_i_new = a_i[jcount] + 12/(delta*delta)*Reweightexpect;
           //node0_printf("a = %.4g off \n", a_i_new);
           //avaluesoutput << "a = " << a_i_new << " off \n";
           ofs << a_i_new << " \n";
           std::cout << "a = " << a_i_new << " off" <<std::endl;
				}
				else
				{
					a_i_new = a_i[jcount] + 12/(delta*delta*(RobMarchcount+1-80))*Reweightexpect;
           //node0_printf("a = %.4g \n", a_i_new);
           //avaluesoutput << "a = " << a_i_new << " \n";
           ofs << a_i_new << " \n";
           std::cout << "a = " << a_i_new <<std::endl;
				}
        RobMarchcount = RobMarchcount + 1;
        std::cout << "RobMarchcount: " << RobMarchcount << std::endl;
        
      }
      
      a[Eint] = a[Eint]+a_i_new/Njacknife;
      
    }
    //node0_printf("a = %.4g \n", a[Eint]);
    //printf("%d \n",1.0);
    //fprintf(fpavalue,"Enter a sentence:\n");
    //fprintf(fpavalue, "%d \n",1.0);
  }
  //avaluesoutput.close();
  ofs.close();
  */
 /*
  Update::GaugeAction* gaugeActionEfind = Update::GaugeAction::getInstance("StandardWilson",5.6);
	Update::long_real_t energyEfind = gaugeActionEfind->energy(*environment);
  double energyout = energyEfind;
  std::cout << "Gaugeaction " << energyEfind << std::endl;
  */
	//Print and finalize the output to file
	

#ifdef ENABLE_MPI
	Lattice::MpiLayout<Lattice::ExtendedStencil>::destroy();
	Lattice::MpiLayout<Lattice::StandardStencil>::destroy();
	Lattice::MpiLayout<Lattice::ReducedStencil>::destroy();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	return 0;
}

