/*!
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
 */
 /* This file is part of HyPro.
 *
 * HyPro is free software: you redistribute it and/or modify
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

#ifndef PROPMODULE_H_
#define PROPMODULE_H_

#include <Node.h>
#include <memory>
#include <goocanvasmm.h>
#include <map>
#include <execinfo.h>
#include <unistd.h>
#include <iostream>
#include <typedefs.h>
#include <fstream>

#include <MyNode.h>
// #include <Node.h>
// #include "typedefs.h"

// include headers that implement a archive in simple text format
//#include <boost/archive/xml_oarchive.hpp>
//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/serialization/export.hpp>

namespace hypro {
	class Collection;
	class ModuleGraph;
	class Archive;

///Base class for all the propulsive modules.
	class propModule : public MyNode {
		friend class Collection;

		///ID of the next module that will be created
		/**It is assigned during creation of the module and then it increased.
         * It starts from 1, since 0 is reserved to identify NULL pointers.
         */
		static unsigned nextID_;

		///Class used to create objects of subclasses given the class type name.
		static class Populator {
		public:
			Populator();

			~Populator();
		} pop_;

		///ID of the chokingFeedback module
		/**It is used only during unserialization.
         * Indeed a two step serialization is needed because the feedback module could
         * not be defined during the first step.
         */
		unsigned chokingFeedbackID_;

	protected:
		static const double pi_; ///<Pi value
		static const double toll_; ///<Value of tolerance used in many numerical calculations
		static const int maxIter_ = 1000; ///<Maximum number of iterations for numerical calculations
		static const double g0_; ///<Gravity acceleration for Earth [m.s-2]
		static const std::string defName_; ///<Default name for new modules

		///Pointer to the graph representing the module
		Glib::RefPtr<ModuleGraph> graph_;

		Node *N1_; ///<Input node
		Node *N2_; ///<Output node

		///Map between nodes and their IDs. Used in serialization
		typedef std::map<unsigned, std::shared_ptr<Node> > NodeMap;

		///Map between modules and their IDs. Used in serialization
		typedef std::map<unsigned, std::shared_ptr<propModule> > ModuleMap;

		///The factory function used to create an object in unserialize
		typedef std::shared_ptr<propModule> (*Factory)(std::istream &,
													   const NodeMap &nodeMap,
													   const ModuleMap &moduleMap);

		///Stores a pointer to the creator function for each object type
		/**This map is populated in the .cpp file of each class.
         * The map is then called in the unserialize method.
         */
		static std::map<std::string, Factory> typeCreators_;

		static std::string indentation_; ///<Indentation, modified during serialization

		///Construct from input stream and external maps
		/**This constructor is protected since it is only meant to be called by
         * the Create method.
         * The object constructed from stream does not include information about GP node
         */
		propModule(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

	public:

		unsigned ID_; ///<Unique identifier
		///Verbosity during module calculation.
		/**- `0` Silent calculation
         * - `1` Information about convergence and iteration process printed
         *
         * Other verbosity levels can be added in future versions
         */
		int verbosity_;
		static double chokingMach_; ///<Define the Mach number at which the choking is detected.
		static double chokingToll_; ///<Tolerance for detecting choking
		std::shared_ptr<propModule> chokingFeedback_; ///<Points to the module where the choking feedback is directed

		propModule();

		///Constructor with name
		/**@param[in] N1 input node
         * @param[in] N2 output node
         * @param[in] name name of the module
         *
         * The name of the module will be initialized to the default one.
         */
		propModule(std::string name, Node *N1, Node *N2, int ID = defID_,
				   int args = 1);  ///<To be used in case of module attached to a constant Node (i.e. free stream)
		propModule(const propModule &mod); ///<copy constructor
		virtual ~propModule(); ///<Destructor

		///Assign N1_
		virtual void N1(Node &N1);

		///Get N1_
		Node &N1() const;

		///Assign N2_
		virtual void N2(Node &N2);

		///Get N2_
		Node &N2() const;

		///Get input nodes
		virtual std::vector<const Node *> inNode() const;

		virtual void inNode(std::size_t, Node& N);

		///calculate the module, gives as output a parameter negative if module is over choked
		/**This function has to be implemented in base classes
         * depending on the physics of the module.
         */
		virtual double calculate()=0;

		///Determines whether the module is running off design. Will remain as 1.0 if there is no offdesign considerations or if there are but operating at design point
		/**This function has to be implemented in base classes that contain off design capabilities
		 * depending on the physics of the module.
		 */

		virtual double isoffdesign ();


		///
		/*
		*/
		virtual void setT04(double& T);
		virtual double getT04();

		///Bring back the reduced parameter to its default value.
		/**Need to be reimplemented in classes that can receive choking feedbacks.
         * Does nothing by default.
         */
		virtual void unreduce();

		///Return true if the module has been reduced as consequence of a choking feedback.
		/**Need to be reimplemented in classes that can receive choking feedbacks.
         * Returns false by default.
         */
		virtual bool isreduced() const;

		///Returns true if the module is choked
		virtual bool ischoked() const;

		///Reduce the module as action following to a choking feedback
		/**
         * @param[in] x the value at which the reduce parameter has to be set.
         * It has to be between 1e-5 and the value it reaches when method unreduce() is called.
         *
         * Need to be reimplemented in classes that can receive choking feedbacks.
         * Throws an error by default.
         */
		virtual void reduce(const double &x);

		///return the value of the parameter reduced in method reduce().
		/**Need to be reimplemented in classes that can receive choking feedbacks.
         * Throws an error by default.
         */
		virtual double reduce() const;

		///Returns true if the module is possible to reduce the module
		/**The module is reduceable if it is possible to call the method
         * reduce without throwing any error.
         * All the subclasses that reimplement the reduce methods must
         * reimplement also this method to return true.
         *
         * @return true if reduceable
         */
		virtual bool isreduceable() const;

		virtual bool canChoke()const;

		virtual bool inputFreeStream()const;

		///return the drag produced by the module, such as the intake drag.
		/**Needs to be reimplemented by subclasses with drag.
         * It returns 0 by default.
         */
		virtual double drag() const;

		///Return the mass flow rate of propellant linked with the module.
		/**Normally only the first module that introduce the propellant into the engine
         * has to implement this method, in order to avoid counting the amount of
         * propellant consumed multiple times.
         *
         * @return a list containing the mass flow rate for each chemical species.
         */
		virtual std::vector<double> propellantMfr() const;

		///Draws a graph representing the module
		virtual Glib::RefPtr<ModuleGraph> draw();

		///Returns a pointer to the graph representing the module
		Glib::RefPtr<ModuleGraph> graph() const;

		///Returns true if the module is a collection of modules.
		/**It returns false by default and
         * it is reimplemented in class ::Collection to return true.
         */
		virtual bool isCollection() const;

		///Name of the module type
		/**It is used for serialization.
         */
		virtual std::string typeName() const =0;

		///Write the module to stream
		/**
         * @param[out] os output stream where the object is serialized.
         */
		virtual void serialize(Archive& ar)const;

		virtual void unserialize(const Archive& ar);


		///Create a new module from a stream
		/**The stream must be created previously using method serialize().
         * @param[in] istr is the input stream
         * @param[in] nodeMap is the map of the nodes external of the module being created
         * @param[in] moduleMap is the map of the modules external of the module being created
         */
		static std::shared_ptr<propModule> unserialize(std::istream &istr,
													   const NodeMap &nodeMap,
													   const ModuleMap &moduleMap);

//	friend class boost::serialization::access;

		// When the class Archive corresponds to an output archive, the
		// & operator is defined similar to <<.  Likewise, when the class Archive
		// is a type of input archive the & operator is defined similar to >>.
		template<class Archive>
		void save(Archive &ar, const unsigned int version) const;

		template<class Archive>
		void load(Archive &ar, const unsigned int version);

		void saveJson(std::ostream& os)const;
		void saveJson(std::string fname)const;
		static std::shared_ptr<propModule> loadJson(std::istream& is);
		static std::shared_ptr<propModule> loadJson(std::string fname);
	};

#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#define propModuleSERIAL(TYPENAME) \
	virtual std::string typeName()const{ \
		return QUOTE(TYPENAME); \
	} \
	/**Creates a modules of this class from stream.*/ \
	/**This function is used by unserialize. */ \
	static std::shared_ptr<propModule> Create(std::istream& is, \
			const NodeMap& nodeMap, \
			const ModuleMap& moduleMap){ \
		try{\
			return std::shared_ptr<propModule>(new TYPENAME(is,nodeMap,moduleMap)); \
		}catch(std::exception& e){\
			void *array[10];\
			size_t size;\
			size = backtrace(array, 10);\
			/*print out all the frames to stderr*/ \
			backtrace_symbols_fd(array, size, STDERR_FILENO);\
			throw std::runtime_error("\nError: problem in unserialization for type " + std::string(QUOTE(TYPENAME)));\
		}\
	} \
	TYPENAME(std::istream& is, \
			const NodeMap& nodeMap, \
			const ModuleMap& moduleMap); \
	static std::shared_ptr<propModule> Create(){ \
		return std::shared_ptr<propModule>(new TYPENAME()); \
	}

#define propModuleSERIALtemplate(TYPENAME) \
	virtual std::string typeName()const; \
	/**Creates a modules of this class from stream.*/ \
	/**This function is used by unserialize. */ \
	static std::shared_ptr<propModule> Create(){ \
		return std::shared_ptr<propModule>(new TYPENAME()); \
	} \
	TYPENAME(std::istream& is, \
			const propModule::NodeMap& nodeMap, \
			const propModule::ModuleMap& moduleMap);
}
#endif /* PROPMODULE_H_ */
