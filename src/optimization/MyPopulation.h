/*!
 *  \brief      implementation of GP Population for HyPro
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

#ifndef MYPOPULATION_H_
#define MYPOPULATION_H_

#include <gp.h>
#include <OptimSystem.h>

namespace hypro {
class MyPopulation: public GPPopulation {
	///No input constructor only to be used for class registration
	MyPopulation():
		GPPopulation(),
		freeStream_(NULL),
		exit_(NULL),
		GPtype_(new OptimSystem()){}

public:
	///GP classes IDs, used during serialisation
	enum ClassID{
			MyPopulationID = GPUserID,
			OptimSystemID  = GPUserID + 1,
			MyGeneID       = GPUserID + 2,
			ModuleAdditionGeneID = GPUserID + 3,
			FeedbackDefinitionGeneID = GPUserID + 4,
			GPTypedNodeSetID = GPUserID + 5,
			ParameterDefinitionGeneID = GPUserID + 6,
			CrossAreaConstraintID = GPUserID + 7,
			ComplexityPenaltyID = GPUserID + 8
		};

	///The selected GP type to be used in this population
	std::unique_ptr<OptimSystem> GPtype_;

	Node* freeStream_; ///<Reference to the HyPro system free stream node (i.e. node N1 of class systemModule)
	Node* exit_; ///<Reference to the HyPro system exit node (i.e. node N2 of class systemModule)

	///Constructor with with no GP inputs but only free stream and exit
	MyPopulation(Node& freeStream, Node& exit):
			GPPopulation(),
			freeStream_(&freeStream),
			exit_(&exit){}

	///Constructor
	/**
	 * @param GPVar_ GP set up parameters
	 * @param adfNs_ node sets used during the GP optimization
	 * @param freeStream Reference to the HyPro system free stream node (i.e. node N1 of class systemModule)
	 * @param exit Reference to the HyPro system exit node (i.e. node N2 of class systemModule)
	 */
	MyPopulation (GPVariables& GPVar_, GPAdfNodeSet& adfNs_, Node& freeStream, Node& exit);

	///Copy constructor
	MyPopulation (MyPopulation& gpo);

	///Duplicate the MyPopulation object
	virtual GPObject& duplicate () { return *(new MyPopulation(*this)); }

	///Check if the GP genome is valid
	/**
	 * Don't check for ultimate diversity as it takes a long time.
	 * Accept every created GP
	 */
	virtual int checkForValidCreation (GP&) { return 1; }

	///Instantiates a permanent GP object
	/**
	 * This function passes also a reference of freeStream_ and exit_ to the GP object
	 * @param numOfGenes total number of genes contained in the GP object
	 * @return pointer to the created GP object
	 */
	virtual GP* createGP (int numOfGenes);

	///Prints out the complete population
	virtual void printOn (ostream& os);

	///Returns the class ID. Used during serialisation
	virtual int isA () { return MyPopulationID; }

	///Instantiate permanent object of class MyPopulation
	virtual GPObject* createObject(Node& freeStream, Node& exit) { return new MyPopulation(freeStream,exit); }

	///Register class IDs for serialisation purposes
	static void registerOptimClasses ();

	///Load population from input file
	virtual char* load (istream& is);
};
}
#endif /* MYPOPULATION_H_ */
