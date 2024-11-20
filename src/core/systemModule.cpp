/*!
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

#include "systemModule.h"
#include "nonLinear.h"
#include "SystemGraph.h"
#include "AdaptedNozzle.h"
#include "solvers/isentropicDuct.h"
#include "Inlet.h"
#include "Archive.h"

namespace hypro {
	bool systemModule::resetLoopLimits_ = false;

	systemModule::systemModule(std::string name, Node *N1, Node *N2)
			:
			Collection(name, N1, N2),
			offError_{1.0}
	{

	}

	systemModule::systemModule(const systemModule &mod)
			:
			Collection(mod),
			offError_{mod.offError_}
{
	}

/*systemModule::systemModule(std::istream& is, const NodeMap& nodeMap, const ModuleMap& moduleMap)
:
		Collection(is,nodeMap,moduleMap){
}*/


	systemModule::~systemModule() {
		// TODO Auto-generated destructor stub
	}

	double systemModule::calculate() {

		double rChok;
		int rOD;


		// First, find all modules that are running off design and index them in a list


		std::list<ModulePtr>::iterator ODstart = modules_.begin();
		std::list<ModulePtr> offdesignlist;
		std::list<ModulePtr>::iterator OffPtrit = ODstart;
		ModulePtr OffPtr;


		while (OffPtr != *(modules_.end())) {//not end

			const void *par[1] = {&OffPtrit};
			OffPtr = offdesigncheck(par);

			if (OffPtr != *(ODstart) && OffPtr != *(modules_.end()))  {
				offdesignlist.push_front(OffPtr);
				OffPtrit = offdesignlist.begin();
				for (std::list<ModulePtr>::iterator it = modules_.begin(); it != modules_.end(); it++)
					{

					if ((*(*(OffPtrit))).Name_ == (*it)->Name_){

						OffPtrit =  ++it;
					}

				}

			}

		}

		offdesignlist.reverse();
		std::list<ModulePtr> offdesignlistMaster (offdesignlist);


		// Next create a list of pairs of corresponding off design modules

		OffDesignPairedList PairedList;

		if (offdesignlist.size() ){

			const unsigned int paircounter = offdesignlist.size()/2;
			rOD = offdesignlist.size();

			for (std::size_t loop2 = 0; loop2 < paircounter; loop2++){


				std::list<ModulePtr>::iterator Start = offdesignlist.begin();
				offdesignlist.reverse();
				std::list<ModulePtr>::iterator End = offdesignlist.begin();
				offdesignlist.reverse();
				if (verbosity_ > 0) {
					std::cout<< "Start Points to: "<< (*Start)->Name_	<< "\n";
					std::cout<< "End Points to: "<< (*End)->Name_	<< "\n";
				}

				PairedList.push_back(std::pair<ModulePtr,ModulePtr>((*Start),(*End)));
				offdesignlist.erase(End);
				offdesignlist.erase(Start);
				if (verbosity_ > 0) {
					std::cout<< "The first (start) module is " <<(*PairedList.begin()).first->Name_<< "\n";
					std::cout<< "The second (end) module is " <<(*PairedList.begin()).second->Name_<< "\n";
				}

			}

		}

		PairedList_ = PairedList;


		//////////////////////////////////////////////////////////////////////////
		// Original calculate code with choking feedback loop

		std::list<ModulePtr>::iterator it = modules_.begin();
		while (it != modules_.end()) {
			if (verbosity_ > 0) {
				std::cout << "Calculating module: " << (*it)->Name_ << "\n";
			}


			rChok = (*it)->calculate();

			if (verbosity_ > 0) {
				(*it)->N1().printOut(std::cout, true);
				(*it)->N2().printOut(std::cout, true); // Print node by node output during run through
			}
			//std::cout << "rChoke = " << rChok << std::endl;	// 
			if (rChok < 0) { // might not be a zero error situation
				if (verbosity_ > 0) {
					std::cout << "The flow is choked." << "\n";
				}
				//std::cout << "TESTING - SYST_MOD.CPP - LINE 156 " << std::endl;	// TRIGGERED
				std::list<ModulePtr>::iterator Start = std::find(modules_.begin(), modules_.end(),
																 (ModulePtr)((*it)->chokingFeedback_));
				//std::cout << "TESTING - SYST_MOD.CPP - LINE 159 " << std::endl;	// TRIGGERED
				std::list<ModulePtr>::iterator End = it;
				End++;
				const void *par[2] = {&Start, &End};
				//std::cout << "TESTING - SYST_MOD.CPP - LINE 163 " << std::endl;	// TRIGGERED

				double xMin = 1e-5; //the lowest value considered in reduction is not 0 to avoid error in calculation
				//std::cout << "choking feedback = " << chokingFeedback_ << std::endl;	// 0

				if ((*it)->chokingFeedback_ == 0) {
					throw std::runtime_error("Error: flow choked in module: " + (*it)->Name_ +
											 ", that has not choking feedback defined");
				}
				//std::cout << "TESTING - SYST_MOD.CPP - LINE 171 " << std::endl;	// NOT TRIGGERED
				//Make sure the full reduction is enough to avoid choking
				if (resetLoopLimits_) {
					(*it)->chokingFeedback_->unreduce(); //set reduction parameter to maximum
				}
				double xMax = (*it)->chokingFeedback_->reduce(); //get maximum value of reduction parameter
				double r_xMin = findChok(xMin, par);
				if (r_xMin < 0) {
					throw std::runtime_error("Error: flow always choked in module: " + (*it)->Name_ +
											 ", not possible to reduce " + (*Start)->Name_ + " lower than 0");
				}

				nonLinear<systemModule>::SolveFalsePos(&systemModule::findChok, xMin, xMax, *this, par);

				/*Reduce the found result by a small amount to make sure it actually avoid choking
                 * Indeed The SolveFalsePos find a solution that gives almost 0
                 * but do not guarantees that it is greater than 0*/
				double xred = (*it)->chokingFeedback_->reduce();
				double r_fin = -1;
				const int maxTry = 10;
				int i = 0;
				while (r_fin < 0 && i < maxTry) {
					xred *= 1 - 1e-5;
					r_fin = findChok(xred, par);
					i++;

				}
				if (r_fin < 0) {
					throw std::runtime_error("Error: Final solution found still produce choking.");
				}
            } else {
				//std::cout << "choking feedback = " << chokingFeedback_ << std::endl;	//
				it++;
		    }
        

		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Off Design Calculation Comes here


//		if (PairedList.size()){
//
////			printOut(std::cout);
////
////			std::cout << "\n" << "Now Calling offdesigncalculatelist() " << "\n" << "\n";
////
////			offdesignlistcalculate(PairedList);
//		}

		const int pair = 2;


		for (std::size_t loop = 0; loop < 3; loop++){

			offdesignlist = offdesignlistMaster;
			rOD = offdesignlist.size();

			while (rOD != 0)
			{

				if (verbosity_ > 0) {
					std::cout << "\n" << "Off design system module count has found: " << rOD << " modules running off design point" <<  "\n" << "\n";
					for (std::list<ModulePtr>::iterator ii = offdesignlist.begin(); ii != offdesignlist.end(); ii++)
					{
						 std::cout << (*ii)->Name_ << "\n";
					}
				}


				std::list<ModulePtr>::iterator OffStart = offdesignlist.begin();
//				std::advance(OffStart, ((rOD/pair) - 1));
				offdesignlist.reverse();
				std::list<ModulePtr>::iterator OffEnd = offdesignlist.begin();
//				std::advance(OffEnd, ((rOD/pair) - 1));
				offdesignlist.reverse();


				const void *ODpar[2] = {&OffStart, &OffEnd};
				offdesigncalculate(ODpar);

				std::list<ModulePtr>::iterator itt;
				for (std::list<ModulePtr>::iterator it = modules_.begin(); it != modules_.end(); it++){

					if ((*(OffEnd))->Name_ == (*it)->Name_){
						itt =  it;

					}

				}

				while (itt != modules_.end()) {

//					std::cout << "Propogating downstream " << "\n";
					rChok = (*itt)->calculate();
					itt++;
				}


				offdesignlist.erase(OffEnd);
				offdesignlist.erase(OffStart);
				rOD = offdesignlist.size();


			}

		}



        return rChok;
    }




    double systemModule::findChok(double &x, const void *par[]) {

		std::list<ModulePtr>::iterator Start = *(std::list<ModulePtr>::iterator *) (par[0]);
		std::list<ModulePtr>::iterator End = *(std::list<ModulePtr>::iterator *) (par[1]);

		double rChok = 1.0;

		if (verbosity_ > 0) {
			std::cout << "Reducing module: " << (*Start)->Name_ << " x=" << x << "\n";
		}
		(*Start)->reduce(x);

		for (std::list<ModulePtr>::iterator it = Start; it != End; it++) {
			if (rChok < 0.0) {
				throw std::runtime_error("Module choked within reduce loop.");
				//TODO maybe this should be handled with recursive calling of main system iteration function.
				//TODO this probably influence the solution of GP optimisation, but probably it is right!
			}
			if (verbosity_ > 0) {
				std::cout << "Calculating module: " << (*it)->Name_ << "\n";
			}
			rChok = (*it)->calculate();
		}
		if (verbosity_ > 0) {
			std::cout << " r=" << rChok << "\n";
		}
		return rChok;
	}

//


    double systemModule::offdesigncalculate(const void *par[]){

       	double rOD  = 1.0;
       	offError_ = 1.0;

       	std::list<ModulePtr>::iterator OffDesignModuleBegin = *(std::list<ModulePtr>::iterator *) (par[0]);
   		std::list<ModulePtr>::iterator OffDesignModuleEnd   = *(std::list<ModulePtr>::iterator *) (par[1]);
   		std::list<ModulePtr>::iterator OffDesignIterationStart;
   		std::list<ModulePtr>::iterator OffDesignIterationEnd;



   		for (std::list<ModulePtr>::iterator it = modules_.begin(); it != modules_.end(); it++) {
   			if ((*(*(OffDesignModuleBegin))).Name_ == (*it)->Name_){
   				OffDesignIterationStart =  it;
   			}

   			if ((*(*(OffDesignModuleEnd))).Name_ == (*it)->Name_){
   				OffDesignIterationEnd =  ++it;

   			}
   		}

//   		std::cout << "\n" << "Beginning Iteration Module Iterator is: " << (*OffDesignIterationStart)->Name_;
//   		std::cout << "\n" << "End Iteration Module Iterator is: " << (*OffDesignIterationEnd)->Name_;
//   		std::cout << "\n" << "Off Design Start Control Module is: " << (*OffDesignModuleBegin)->Name_;
//   		std::cout << "\n" << "Off Design End Control Module is: " << (*OffDesignModuleEnd)->Name_;


   		double T = ((*OffDesignModuleEnd))->getT04();
   		(*(*OffDesignModuleBegin)).setT04(T);
   		double first = 0.0;


   		while (offError_ >offToll_){

   			double T04start = ((*OffDesignModuleBegin))->getT04();

   			for (std::list<ModulePtr>::iterator it = OffDesignIterationStart; it != OffDesignIterationEnd; it++) {

   				int test = (*it)->isoffdesign();

   				if (test == 0) {
   					(*it)->calculate();
   					if (verbosity_ > 0) {
   						std::cout << "\n" << "Off Design Calculating module: " << (*it)->Name_;	// COMMENTED
   						std::cout << "\n" << "Calculating module: " << (*it)->Name_ << "\n";		// COMMENTED
   						(*it)->N1().printOut(std::cout, true);					// COMMENTED
   						(*it)->N2().printOut(std::cout, true);					// COMMENTED
   					}
   				}
   				else {
   					(*it)->calculate();
   					if (verbosity_ > 0) {
   						std::cout << "\n" << "Off Design Calculating module: " << (*it)->Name_ ;	// COMMENTED	
   						std::cout << "\n" << "Calculating module: " << (*it)->Name_ << "\n";		// COMMENTED
   						(*it)->N1().printOut(std::cout, true);					// COMMENTED
   						(*it)->N2().printOut(std::cout, true);					// COMMENTED
   					}
   				}

   			}

   			double T04end = ((*OffDesignModuleEnd))->getT04();
   			offError_ = (T04end - T04start);

   			if (verbosity_ > 0) {
   				if (!first){
   					std::cout << "\n";
   					fmt::print(std::cout, "{:<12} {:<12} {:<12}\n", "T04 Start ", "T04 End", "Error");
   					first +=1.0;
   				}
   				std::cout << "\n";
   				fmt::print(std::cout, "{:<12.5f} {:<12.5f} {:<12.9f}\n", T04start, T04end, offError_);
   			}

   			double newTemp = (T04end + T04start)/2;
   			(*(*OffDesignModuleBegin)).setT04(newTemp);

   			offError_ = std::abs(offError_);
   		}

   		return rOD;


       }



	std::shared_ptr<propModule> systemModule::offdesigncheck(const void *par[]) {

		std::list<ModulePtr>::iterator OffStart = *(std::list<ModulePtr>::iterator *) (par[0]);

		ModulePtr fOffPtr;
		ModulePtr DefPtr = *(modules_.end());

		for (std::list<ModulePtr>::iterator it = OffStart; it != modules_.end(); it++) {
            if ((*it)->isoffdesign() == 0 && it != modules_.end()) {
            	fOffPtr = (*it);
            	return fOffPtr;
            }
        }

        return DefPtr;
    }



    void systemModule::offdesignlistcalculate(OffDesignPairedList PairedList){


		int layers = PairedList.size();

		std::list<ModulePtr>::iterator ODstart;
		std::list<ModulePtr>::iterator ODend;
		OffDesignPairedList::iterator itpl = PairedList.begin();

		for(int i = 0; i < 2; i++) {

			for (OffDesignPairedList::iterator itpl = PairedList.begin(); itpl != PairedList.end(); itpl++){



					for (std::list<ModulePtr>::iterator it = modules_.begin(); it != modules_.end(); it++) {
							if ((*itpl).first->Name_ == (*it)->Name_){
								ODstart =  it;
							}

							if ((*itpl).second->Name_ == (*it)->Name_){
								ODend =  it;

							}
					}

				const void *par[2] = {&ODstart, &ODend};

//				std::cout << "Calling offdesignTemperatureSeek()" << "\n";
				nonLinear<systemModule>::funSolve(&systemModule::offdesignTemperatureSeek, Tmin_, Tmax_, *this, par);



				for (std::list<ModulePtr>::iterator it = ODend; it != modules_.end(); it++) {

						(*it)->calculate();
				}
			}

		}

	}


    double systemModule::offdesignTemperatureSeek(double &T, const void *par[]){

    	std::list<ModulePtr>::iterator OffDesignModuleBegin = *(std::list<ModulePtr>::iterator *) (par[0]);
		std::list<ModulePtr>::iterator OffDesignModuleEnd   = *(std::list<ModulePtr>::iterator *) (par[1]);
		std::list<ModulePtr>::iterator OffDesignIterationStart;
		std::list<ModulePtr>::iterator OffDesignIterationEnd;


		// Finds the two components running off design in their order as designated in modules_ and creates a list to iterate through.

		for (std::list<ModulePtr>::iterator it = modules_.begin(); it != modules_.end(); it++) {
			if ((*(*(OffDesignModuleBegin))).Name_ == (*it)->Name_){
				OffDesignIterationStart =  it;
			}

			if ((*(*(OffDesignModuleEnd))).Name_ == (*it)->Name_){
				OffDesignIterationEnd =  ++it;

			}
		}

//		std::cout << "\n" << "Off Design Start Control Module is: " << (*OffDesignModuleBegin)->Name_;
//		std::cout << "\n" << "Off Design End Control Module is: " << (*OffDesignModuleEnd)->Name_;
//		std::cout << "\n";

		(*(*OffDesignModuleBegin)).setT04(T);

		double T04start = ((*OffDesignModuleBegin))->getT04();

		for (std::list<ModulePtr>::iterator it = OffDesignIterationStart; it != OffDesignIterationEnd; it++) {


			for (OffDesignPairedList::iterator itpl = PairedList_.begin(); itpl != PairedList_.end(); itpl++){

//				std::cout << "Looping" << "\n";


//					if ((*it)->Name_ == itpl->first->Name_ && (*it)->Name_ != (*(OffDesignIterationStart))->Name_) {
//							//Use find to search the PairedList for a thing of that name
//						std::list<ModulePtr>::iterator ODstartInterrupt;
//						std::list<ModulePtr>::iterator ODendInterrupt;
//
//						std::cout << "Inside" << "\n";
//
//
//						for (std::list<ModulePtr>::iterator it2 = modules_.begin(); it2 != modules_.end(); it2++) {
//
//								if ((*itpl).first->Name_ == (*it2)->Name_){
//									ODstartInterrupt =  it2;
//									std::cout << "\n" << "ODstartInterrupt: " << (*it2)->Name_ << "\n";
//
//								}
//
//								if ((*itpl).second->Name_ == (*it2)->Name_){
//									ODendInterrupt =  it2;
//
//								}
//						}
//
//
//
//						const void *par[2] = {&ODstartInterrupt, &ODendInterrupt};
//						nonLinear<systemModule>::funSolve(&systemModule::offdesignTemperatureSeek, Tmin_, Tmax_, *this, par);
//
//					}
				}


			(*it)->calculate();

		}



		double T04end = ((*OffDesignModuleEnd))->getT04();
		offError_ = (T04end - T04start);


		return offError_;


	}





	double systemModule::thrust() const {

	  //double T = N2_->A() * N2_->rho() * std::pow(N2_->getU(), 2.0) - N1_->A() * N1_->rho() * std::pow(N1_->getU(), 2.0) +
		//		   N2_->A() * (N2_->getPress() - N1_->getPress());
		double T = N2_->A() * N2_->rho() * std::pow(N2_->getU(), 2.0) +
				   N2_->A() * (N2_->getPress() - N1_->getPress());
		
		double D = 0.0;
		for (std::list<ModulePtr>::const_iterator it = modules_.begin(); it != modules_.end(); it++) {
			D = D + (*it)->drag();
		}

		return T - D;
	}

	void systemModule::thrustPath() {
		Node N2tp(*N2_);
		N2tp.Amax_ = std::numeric_limits<double>::infinity();
		N2tp.Amin_ = 0;
		AdaptedNozzle<isentropicDuct> nozzle("nozzle", N2_, &N2tp);
		nozzle.choked_ = true;
		nozzle.freeStream_ = N1_;
		N2tp.Name_ = "thrustPath-out";
		for (TreeIterator it = --end(); !it.isPastBegin(); it--) {
			nozzle.N1(it->N2());
			try {
				nozzle.calculate();
			} catch (std::runtime_error &e) {
				std::cout << it->Name_ << " " << "Not defined" << "\n";
				continue;
			}
			double T = N2tp.A() * N2tp.rho() * std::pow(N2tp.getU(), 2.0) +
					   N2tp.A() * (N2tp.getPress() - N1_->getPress());

			double D = 0.0;
			for (TreeIterator it1 = TreeIterator(&*it, *this); !it1.isPastBegin(); it1--) {
				D = D + it1->drag();
			}
			std::cout << it->Name_ << " " << T - D << "\n";
		}
	}

	Glib::RefPtr<ModuleGraph> systemModule::draw() {
		graph_ = SystemGraph::create(this);
		return graph_;
	}


	double systemModule::crossArea()const{
		// Calculate engine cross section area
		double S = 0.0;
		for(std::list<ModulePtr>::const_iterator it=modules_.begin(); it!=modules_.end(); it++){
			if((*it)->drag()!=0.0){
				if(!(*it)->isCollection()){
					throw std::logic_error("Error: cross area can be calculated only for collections.");
				}
				Inlet* mod = (Inlet*)((*it).get());
				S += mod->nodes_[0]->Amax_;
			}
		}
		return S;
	}

	void systemModule::serialize(Archive& ar)const{

		ar.put("Thrust",this->thrust());
		ar.put("MFR", (this->N2().mfr() - this->N1().mfr()));

		Collection::serialize(ar);

	}

	void systemModule::deleteAll(systemModule* sys){
		//delete(&sys->N1().specieGas());
		delete(&sys->N1());
		delete(&sys->N2());
		delete(sys);
	}
}

