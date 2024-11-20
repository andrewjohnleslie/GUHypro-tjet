/*!
 *  \brief      GP optimisation test case
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
 */
 /* This file is part of HyPro.
 *
 * HyPro is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HyPro is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>. */

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
#include <fstream>
#include <backward/strstream>

#include <stdlib.h>
//#include <new.h>    // For the new-handler
#include <math.h>   // fabs()
#include <string>

#include <list>

// Include header file of genetic programming system.
#include "gp.h"
#include "gpconfig.h"

#include <MyPopulation.h>
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Phase.h"

#include <solvers/EffIsenDuct.h>
#include <core/systemModule.h>
#include <combustion/Combustion.h>
#include <injectors/InjectionPhi.h>
#include <AdaptedThroatInlet.h>
#include <SupersonicInlet.h>
#include <Wedge.h>
#include <Rocket.h>
#include <injectors/InjectionPlate.h>
#include <ParMixer.h>
#include <AdaptedNozzle.h>
#include <core/MultiModeMachP0.h>
#include <AdaptedInlet.h>
#include <NeutralLink.h>
#include <core/FeedBack.h>
#include <TermModule.h>
#include <solvers/Friction.h>
#include <examplewindow.h>
#include <nonLinear.h>
#include <GPTypedGeneSet.h>
#include <ModuleAdditionGene.h>
#include <FeedbackDefinitionGene.h>
#include <ParameterDefinitionGene.h>
#include <CrossAreaConstraint.h>
#include <ComplexityPenalty.h>

#include <fenv.h>

#ifdef MATLAB
#include "mex.h"
#endif
using namespace hypro;

// The TeX-file
ofstream tout;
int printTexStyle=0;

// Define configuration parameters and the neccessary array to
// read/write the configuration to a file. If you need more variables,
// just add them below and insert an entry in the configArray.
GPVariables cfg;
const char *InfoFileName="data";
struct GPConfigVarInformation configArray[]=
{
  {"PopulationSize", DATAINT, &cfg.PopulationSize},
  {"NumberOfGenerations", DATAINT, &cfg.NumberOfGenerations},
  {"CreationType", DATAINT, &cfg.CreationType},
  {"CrossoverProbability", DATADOUBLE, &cfg.CrossoverProbability},
  {"CreationProbability", DATADOUBLE, &cfg.CreationProbability},
  {"MaximumDepthForCreation", DATAINT, &cfg.MaximumDepthForCreation},
  {"MaximumDepthForCrossover", DATAINT, &cfg.MaximumDepthForCrossover},
  {"SelectionType", DATAINT, &cfg.SelectionType},
  {"TournamentSize", DATAINT, &cfg.TournamentSize},
  {"DemeticGrouping", DATAINT, &cfg.DemeticGrouping},
  {"DemeSize", DATAINT, &cfg.DemeSize},
  {"DemeticMigProbability", DATADOUBLE, &cfg.DemeticMigProbability},
  {"SwapMutationProbability", DATADOUBLE, &cfg.SwapMutationProbability},
  {"ShrinkMutationProbability", DATADOUBLE, &cfg.ShrinkMutationProbability},
  {"SteadyState", DATAINT, &cfg.SteadyState},
  {"AddBestToNewPopulation", DATAINT, &cfg.AddBestToNewPopulation},
  {"InfoFileName", DATASTRING, &InfoFileName},
  {"", DATAINT, NULL}
};

///Input directory
/**
 * - Change to . if input files are in the running directory
 * - append desired subfolder to choose between subsonic test case and supersonic test case
 */
std::string inDir = QUOTE(GPRescource);
//const char* inDir = inDirStr.c_str();

// This function is the divide with closure. Basically if you divide
// anything by zero you get an error so we have to stop this
// process. We check for a very small denominator and return a very
// high value.
inline double divide (double x, double y)
{
  if (fabs (y)<1e-6)
    {
      if (x*y<0.0)
	return -1e6;
      else
	return 1e6;
    }
  else
    return x/y;
}

///Pointer to nodes used by node set modules
static Node* nodes[10];

/// Create functions and terminal set
/**
 * @param[out] adfNs created Node set
 * @param freeStream reference to free stream node
 */
void createNodeSet (GPAdfNodeSet& nodeSet, Node& freeStream)
{
	double foot =  0.3048; //conversion feet to meters
	double lbm = 0.45359; //conversion libre to kg

	// Reserve space for the node sets
	nodeSet.reserveSpace (1);

	// Now define the function and terminal set for each ADF and place
	// function/terminal sets into overall ADF container
	GPTypedGeneSet& ns0=*new GPTypedGeneSet (7);
	nodeSet.put (0, ns0);

	//Define nodes
	double Aprim = 11.25 - 8.24;
	//Nodes                 N1_    I0    I1       0     1       2     3    4    5      6     7   N2_
	const int numNodes = 12;
	double A[numNodes] =   {27.0, 27.0, 0.25*27, 11.25, 11.25, 11.25, 0.0, 0, 22.5, 22.5, 22.5, 95.0};
	/*name = {'Pre Intake','Intake','Throat','Pinch Point','Primary','Mixer End',...
  		    'Dummy','Rocket Outlet','End Diffuser','Injection','End Chamber','Nozzle'};*/
	double foot2 = std::pow(foot,2.0);
	A[3] = A[3] - Aprim;
	A[7] = Aprim;

	freeStream.A(A[0]*foot2);
	for(int i=0; i<8; i++){
		nodes[i] = new Node(freeStream);
		nodes[i]->A(A[i+3]*foot2);
	}
	Node* nozzExit = new Node(freeStream);
	nodes[8] = nozzExit;
	Node* dummy = new Node(freeStream);
	dummy->Name_ = "dummy";
	nodes[9] = dummy;
	nozzExit->A(A[numNodes-1]*foot2);
	nozzExit->Amax_ = A[numNodes-1]*foot2;

	// diffuser
	EffIsenDuct* diffuser = new EffIsenDuct("Duct",dummy,nodes[1],2);
	diffuser->etap_ = 0.98;
	ns0.putNode (*diffuser);
	diffuser->N2().Name_ = "Duct-out";
	ModuleAdditionGene* diffgene = new ModuleAdditionGene(*diffuser);
	diffgene->addParameter(
			new ModuleAdditionGene::DefinedPar<Node>(0.5*diffuser->N2().A(), 1.5*diffuser->N2().A(), &Node::A));
	ns0.putGene(diffgene);

	//Inlet
	AdaptedThroatInlet* inlet = new AdaptedThroatInlet("Inlet",dummy,nodes[0],1);
	inlet->N1(freeStream);
	inlet->nodes_[0].get()->A(A[1]*foot2);
	inlet->nodes_[0]->Amax_ = inlet->nodes_[0]->A();
	inlet->nodes_[1]->A(A[2]*foot2);
	inlet->nodes_[1]->Amax_ = inlet->nodes_[1]->A();
	SupersonicInlet* conv = (SupersonicInlet*)inlet->modules_.back().get();
	conv->etap_ = 0.98;
	//conv->etaT_ = 0.9;
	inlet->verbosity_ =0;
	ns0.putNode (*inlet);
	inlet->N2().Name_ = "Inlet-out";
	ModuleAdditionGene* ingene = new ModuleAdditionGene(*inlet);
	ingene->addParameter(
				new ModuleAdditionGene::DefinedPar<Node>(0.1,
						7.0, &Node::Amax_, 0));
	ingene->addParameter(
			new ModuleAdditionGene::DefinedPar<Node>(0.5*inlet->nodes_[1]->Amax_,
					3.0*inlet->nodes_[1]->Amax_, &Node::Amax_, 1));
	ingene->addParameter(
					new ModuleAdditionGene::DefinedPar<Node>(0.5*inlet->N2().A(),
							3.0*inlet->N2().A(), &Node::A));
	ns0.putGene(ingene);

	//Rocket
	InjectionPlate* inj = new InjectionPlate("RocketInj",dummy,nodes[4],3);
	std::vector<double> Xr1(4,0.0);
	Xr1[0] = 2.0/3.0;
	Xr1[1] = 1.0/3.0;
	inj->Y_ = Xr1;
	inj->M_ = 0.001;
	inj->mfrMax_ = 216*lbm;
	inj->unreduce();
	inj->T_ = 20;
	ns0.putNode (*inj);
	inj->N2().Name_ = "RocketInj-out";
	ns0.putGene(new ModuleAdditionGene(*inj));

	//mixer
	ParMixer* mixer = new ParMixer("Mixer",dummy,nodes[2],dummy,4);
	mixer->verbosity_ = 0;
	mixer->eff_ = 0.67;
	ns0.putNode (*mixer);
	mixer->N2().Name_ = "Mixer-out";
	ns0.putGene(new ModuleAdditionGene(*mixer));

	//Injection
	InjectionPhi<Friction>* injection = new InjectionPhi<Friction>("injection",dummy,nodes[6],5);
	injection->PhiMax_ = 1;
	injection->Phi_ = 1;
	std::vector<double> Xf(4,0.0);
	Xf[0] = 1;
	injection->Xf_ = Xf;
	injection->Tf_ = 20;
	injection->Rcoeff_[0] = 2;
	injection->Rcoeff_[1] = 1;
	injection->Cf_ = 0.01;
	injection->L_ = 1.0;
	injection->D_ = 0.1;
	injection->N2().Name_ = "Injec-out";
	ns0.putNode (*injection);
	ModuleAdditionGene* injgene = new ModuleAdditionGene(*injection);
	injgene->addParameter(
			new ModuleAdditionGene::DefinedPar<InjectionPhi<Friction>>(0.1, 1.0, &InjectionPhi<Friction>::PhiMax_));
	ns0.putGene(injgene);


	//Combustor
	Combustion<Friction>* combustor = new Combustion<Friction>("combustor",dummy,nodes[7],6);
	combustor->Cf_ = 0.01;
	combustor->L_ = 1.0;
	combustor->D_ = 0.1;
	combustor->Rcoeff_[0] = 2;
	combustor->Rcoeff_[1] = 1;
	combustor->Pcoeff_[3] = 2;
	combustor->N2().Name_ = "comb-out";
	ns0.putNode (*combustor);
	ModuleAdditionGene* gene = new ModuleAdditionGene(*combustor);
	gene->addParameter(
			new ModuleAdditionGene::DefinedPar<Combustion<Friction>>(0.5, 2.0, &Combustion<Friction>::L_));
	ns0.putGene(gene);

	// Nozzle
	AdaptedNozzle<EffIsenDuct>* nozzle = new AdaptedNozzle<EffIsenDuct>("nozzle",dummy,nozzExit,7);
	nozzle->choked_ = true;
	nozzle->etap_ = 0.98;
	nozzle->etaT_ = 0.9;
	nozzle->freeStream_ = &freeStream;
	nozzle->N2().Name_ = "nozzle-out";
	ns0.putNode (*nozzle);
	ns0.putGene(new ModuleAdditionGene(*nozzle));

	ns0.putGene(new FeedbackDefinitionGene(0, 8));
	ns0.putGene(new FeedbackDefinitionGene(1, 9));
	ns0.putGene(new FeedbackDefinitionGene(2, 10));

	for(int i=0; i<10; i++){ //TODO add up to 11 not 10
		double p = i/10.0;
		ns0.putGene(new ParameterDefinitionGene(p, i+11));
	}
}

///Delete all nodes used by the Node set modules
void deleteNodes(){
	for(int i=0; i<10; i++){
		delete(nodes[i]);
	}
}


///Error handler
/**
 * @deprecated this function is not used anymore
 */
void newHandler ()
{
	cerr << "\nFatal error: Out of memory." << endl;
	exit (1);
}

///Create the free stream node
/**
 * Read from the file freeStream.dat and create the freeStream node accordingly.
 * @return Return a pointer to the freeStream node
 */
Node createFreeStream(){
	ostrstream strInFile;
	strInFile  << inDir << "/freeStream.dat" << ends;
//	strInFile  << "/home/mark/git/GUHypro/test/GP/resources/subsonic/freeStream.dat" << ends;
//	std::cout << "\n" << strInFile.str() << "\n";
	std::ifstream fil(strInFile.str());
	if(fil.fail()){
		throw std::runtime_error("Error: cannot open file freeStream.dat");
	}

	double p;
	double T;
	double Mach;

	while(!fil.eof()){
		std::string property;
		fil >> property;
		if(strcmp(property.c_str(),"pressure=")==0){
			fil >> p;
		}else if(strcmp(property.c_str(),"Temperature=")==0){
			fil >> T;
		}else if(strcmp(property.c_str(),"Mach=")==0){
			fil >> Mach;
		}else{
			throw std::runtime_error("Error: unrecognized property " + property);
		}
	}

	std::cout << "\np is: " << p << "\nT is: " << T << "\nMach is: " << Mach << "\n";


    Cantera::compositionMap initialComposition;
    initialComposition.emplace("O2", 0.209);
    initialComposition.emplace("N2", 1 - 0.209);

//	Node freeStream("/usr/share/cantera/data/gri30_highT.yaml", "gri30");
	Node freeStream("/home/andrew/GUHypro-tjet/test/GP/airNASA9.yaml", "airNASA9");

	freeStream.setTPX(T,p,initialComposition);
	freeStream.setU(Mach*freeStream.geta());
	return freeStream;
}

///Delete freeStream node and gas model within it.
void deleteFreeStream(Node& freeStream){
	delete(&freeStream);
}

///Error handler
void fpe_handler(int sig){
	if(sig == SIGFPE) {
		std::cerr << "Error: Floating point exception!" << std::endl;
	}else if (sig == SIGSEGV){
		std::cerr << "Error: Segmentation fault!" << std::endl;
	}else{
		std::cerr << "Error signal " << sig << std::endl;
	}
	exit(1);
}

///Inspect results for post processing
/**
 * @param popID generation number of the population to inspect
 * @param gpID GP ID number of the GP individual to inspect
 * @param best whether the GP to investigate is the best of its population. If True gpID is ignored
 * @param rmId ID of the module to be removed from the system module. If negative nothing is removed
 * @param argc number of generic arguments
 * @param argv generic arguments
 * @return pointer to the inspected GP individual
 */
OptimSystem* inspect(int popID, int gpID, bool best, int rmId, int argc, char *argv[]){

	cout << cfg << endl;

	Node freeStream = createFreeStream();
	Node exit(freeStream);

	// Create the adf function/terminal set and print it out.
	GPAdfNodeSet nodeSet;
	createNodeSet (nodeSet, freeStream);
	cout << nodeSet << endl;

	// Create a population with this configuration
	ostrstream inpFile;
	inpFile << "Generation-" << popID << ".obj" << std::ends;
	cout << "Loading file " << inpFile.str() << "..." << endl;
	ifstream loadPop (inpFile.str());
	MyPopulation* pop = new MyPopulation(freeStream,exit);
	char* errMsg = pop->load(loadPop);
	if(errMsg){
		std::cerr << errMsg << std::endl;
		std::exit(1);
	}
	pop->setNodeSets (nodeSet);
	cout << "Ok." << endl;

	OptimSystem* sys;
	if(best){
		sys = (OptimSystem*)pop->NthGP(pop->bestOfPopulation);
	}else{
		sys = (OptimSystem*)pop->NthGP(gpID);
	}
	sys->N1(*pop->freeStream_);
	sys->N2(*pop->exit_);
	sys->verbosity_ = 1;

	std::cout << "Fitness=" << sys->getFitness() << " " << sys->message_ << std::endl;
	std::cout << "Calculating " << *(GP*)sys << std::endl;
	sys->evaluate();

	if(rmId>=0){
		std::list<OptimSystem::ModulePtr>::iterator module(sys->modules_.begin());
		for(int i=0; i<rmId; i++){
			module++;
		}
		std::cout << "Removing module: " << (*module)->Name_ << std::endl;
		propModule* outM = (*module)->N2().downstream();
		propModule* inM = (*module)->N1().upstream();
		sys->remove(*module);
		outM->N1(inM->N2());
		sys->calculate();
	}

	ostrstream outFile;
	if(best){
		outFile << "Gen" << popID << "-GPbest.xml" << std::ends;
	}else{
		outFile << "Gen" << popID << "-GP" << gpID << ".xml" << std::ends;
	}
//	std::ofstream fil(outFile.str());
//	((systemModule*)sys)->save(fil);
//	fil.close();

#ifndef MATLAB
	std::cout << "Fitness=" << sys->getFitness() << " " << sys->message_ << std::endl;
	sys->printOut(std::cout);
	std::cout << "Thrust =" << sys->thrust() << "\nSystem Cross Area =" << sys->crossArea() << std::endl;

	sys->thrustPath();

	sys->createTree(argc,argv);
	sys->show(argc,argv);
#endif

	deleteNodes();

#ifdef MATLAB
	return new OptimSystem(*sys);
#else
	deleteFreeStream(freeStream);

	return NULL;
#endif
}

///Run GP optimization
/**
 * @param outLev level of output verbosity. If 2 save details at each generation.
 * 	If 1 print advancement of optimisation. If 0 silent.\n
 * @param append whether to append new result to preexisting generation data
 */
MyPopulation* run(int outLev, bool append){
	// Open the main output file for data and statistics file. First set
	// up names for data file. We use also a TeX-file where the
	// equations are written to in TeX-style. Very nice to look at!
	// Remember we should delete the string from the stream, well just a
	// few bytes
	ofstream fout;
	ofstream bout;
	if(outLev>0){
		ostrstream strOutFile, strStatFile, strTeXFile;
		strOutFile  << InfoFileName << ".dat" << ends;
		strStatFile << InfoFileName << ".stc" << ends;
		strTeXFile  << InfoFileName << ".tex" << ends;

		if(append){
			fout.open(strOutFile.str(), ios::app);
			bout.open(strStatFile.str(), ios::app);
			tout.open (strTeXFile.str(), ios::app);
		}else{
			fout.open(strOutFile.str());
			bout.open(strStatFile.str());
			tout.open (strTeXFile.str(), ios::out);
			tout << endl
					<< "\\documentstyle[a4]{article}" << endl
					<< "\\begin{document}" << endl;

			// Print the configuration to the files just opened
			fout << cfg << endl;
			if(outLev>0){
				cout << cfg << endl;
			}
			tout << "\\begin{verbatim}\n" << cfg << "\\end{verbatim}\n" << endl;
		}
	}

	Node freeStream = createFreeStream();
	Node exit(freeStream);

	// Create the adf function/terminal set and print it out.
	GPAdfNodeSet nodeSet;
	createNodeSet (nodeSet, freeStream);
	if(outLev>0){
		cout << nodeSet << endl;
		fout << nodeSet << endl;
	}

	// Create a population with this configuration
	MyPopulation* pop;
	int gen0 = 0;
	if(!append){
		if(outLev>0){
			cout << "Creating initial population ..." << endl;
		}
		pop = new MyPopulation (cfg, nodeSet, freeStream, exit);
		pop->GPtype_ = std::unique_ptr<OptimSystem>(new ComplexityPenalty());
		pop->create ();
		if(outLev>0){
			cout << "Ok." << endl;
			pop->createGenerationReport (1, 0, fout, bout);
		}
		if(outLev>1){
			ofstream detout("Generation-0.dat");
			ofstream saving("Generation-0.obj");

			pop->printOn(detout);
			pop->save(saving);
		}
	}else{
		// Create a population with this configuration
		gen0 = 20;
		ostrstream inpFile;
		inpFile << "Generation-" << gen0 << ".obj" << std::ends;
		if(outLev>0){
			cout << "Loading file " << inpFile.str() << "..." << endl;
		}
		ifstream loadPop (inpFile.str());
		pop = new MyPopulation (freeStream, exit);
		char* errMsg = pop->load(loadPop);
		if(errMsg){
			std::cerr << errMsg << std::endl;
			std::exit(1);
		}
		pop->setNodeSets (nodeSet);
		if(outLev>0){
			cout << "Ok." << endl;
		}
	}

	if(outLev>0){
		// Print the best of generation to the LaTeX-file.
		printTexStyle=1;

		tout << *pop->NthGP (pop->bestOfPopulation);
		printTexStyle=0;
	}

	// This next for statement is the actual genetic programming system
	// which is in essence just repeated reproduction and crossover loop
	// through all the generations .....
	MyPopulation* newPop=NULL;
	for (int i=1; i<=cfg.NumberOfGenerations; i++)
	{
		int gen = i + gen0;
		// Create a new generation from the old one by applying the
		// genetic operators
		if (!cfg.SteadyState)
			newPop=new MyPopulation (cfg, nodeSet, freeStream, exit);
		pop->generate (*newPop);

		// Delete the old generation and make the new the old one
		if (!cfg.SteadyState)
		{
			delete pop;
			pop=newPop;
		}

		if(outLev>0){
			// Print the best of generation to the LaTeX-file.
			printTexStyle=1;
			tout << "Generation " << gen << ", fitness "
					<< pop->NthGP (pop->bestOfPopulation)->getFitness()
					<< endl;
			tout << *pop->NthGP (pop->bestOfPopulation);
			printTexStyle=0;

			// Create a report of this generation and how well it is doing
			pop->createGenerationReport (0, gen, fout, bout);
		}

		if(outLev>1){
			ostrstream detFile;
			detFile << "Generation-" << gen << ".dat" << std::ends;
			ofstream detout1(detFile.str());
			pop->printOn(detout1);

			ostrstream saveFile;
			saveFile << "Generation-" << gen << ".obj" << std::ends;
			ofstream saving1(saveFile.str());
			pop->save(saving1);
		}
	}

	if(outLev>0){
		// TeX-file: end of document
		tout << endl
				<< "\\end{document}"
				<< endl;
		tout.close ();

		cout << "\nResults are in "
				<< InfoFileName << ".dat,"
				<< InfoFileName << ".tex,"
				<< InfoFileName << ".stc." << endl;
	}

//	OptimSystem* GPINSPECT;
	inspect(2,1,true,0,0,0);
//	OptimSystem* inspect(int popID, int gpID, bool best, int rmId, int argc, char *argv[]){

	deleteFreeStream(freeStream);
	deleteNodes();
	return pop;
}

///Initialise the GP algorithm classes
/**
 * @param prevSeed if True the seed is loaded from the seed.dat file.
 * 	If False it is generated with time.
 */
void setup(bool prevSeed, std::string seedFile = "seed.dat"){
	nonLinear<systemModule>::maxIter_ = 100;
	systemModule::resetLoopLimits_ = true;
	// We set up a new-handler, because we might need a lot of memory,
	// and we don't know it's there.
	set_new_handler (newHandler);

	// Read seed
	long seed;
	if(prevSeed){
		ifstream fin;
		fin.open(seedFile.c_str(), ios::in);
		fin >> seed;
		fin.close();
	}else{
		seed = time (NULL);
	}

	//Store seed to seed file
	ofstream fout;
	fout.open("seed.dat", ios::out);
	fout << seed;
	fout.close();

	GPInit (0, seed);
	// Init HyPro
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	signal(SIGFPE,(*fpe_handler));
	signal(SIGSEGV,(*fpe_handler));
	thermoState::warning_ = false;
	MyPopulation::registerOptimClasses();

	// Declare the GP Variables, set defaults and read configuration
	// file.  The defaults will be overwritten by the configuration file
	// when read.  If it doesn't exist, the defaults will be written to
	// the file.
	ostrstream strInFile;
	strInFile  << inDir << "/GP.ini" << ends;
//	strInFile  << "~/git/GUHypro/build/test/GP/GP.ini" << ends;
//	std::cout << "\n" << strInFile.str() << "\n";
	GPConfiguration config (cout, strInFile.str(), configArray);
}
