/*!
 *  \brief      GP gene implementation for HyPro
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

#ifndef MYGENE_H_
#define MYGENE_H_

#include <gp.h>
#include <core/Node.h>
#include <OptimSystem.h>
#include <gtkmm.h>
#include <examplewindow.h>
#include <glibmm/refptr.h>

namespace hypro {
class propModule;
class OptimSystem;
class GPTypedGeneSet;

class MyGene : public GPGene
{
protected:
	///Maximum number of attempts in creating a child gene of the correct type
	static const int createAttempt_ = 100;

public:
	///Type of the node
	enum Type {
		MODULE = 0, //!< MODULE
		FEEDBACK = 1, //!< FEEDBACK
		PARAMETER = 2 //!< PARAMETER
	};

	int begin(Type type)const;
	int end(Type type)const;

	friend OptimSystem;

	int verbosity_; ///<Verbosity level

	MyGene():
		GPGene(),
		verbosity_(0){}

	///Constructor based on GP node
	MyGene (GPNode& gpo);

	///Copy constructor
	MyGene (const MyGene& gpo);

	virtual MyGene& operator = (const MyGene& gpo);

	virtual ~MyGene();

	///Reimplements the createChild function
	virtual GPGene* createChild (GPNode& gpo);

	///Grow GPs according to the given parameters
	virtual void create (enum GPCreationType ctype, int allowableDepth,
			GPTypedGeneSet& ns)=0;

	virtual void create (enum GPCreationType ctype, int allowableDepth,
				GPNodeSet& ns);

	///Print a Gene.
	virtual void printOn (ostream& os);

	///Print a Gene, this for now is the only implemented style
	virtual void printMathStyle (ostream& os, int lastPrecedence=0);

	///Print a Gene, this style is for now implemented as a copy of the MathStyle
	void printTeXStyle (ostream& os, int lastPrecedence=0);

	///Returns the nth child of the gene
	MyGene* NthMyChild (int n) {
		return (MyGene*) GPContainer::Nth (n); }

	///Returns the nth child of the gene
	const MyGene* NthMyChild (int n)const {
		return (MyGene*) GPContainer::Nth (n); }

	///Shadow of GPGene functions that should be constant
	virtual int isFunction()const{
		return GPGene::isFunction();
	}

	///Shadow of GPGene functions that should be constant
	virtual int isTerminal()const{
		return GPGene::isTerminal();
	}

	///Create the GP systemModel object starting from the GP tree
	/**
	 * This function is called by the root GP gene that in turns call this same function
	 * for its child and so on.
	 * So the GP tree is spanned recursively.
	 * @param sys reference to the system module being created.
	 * @param parent pointer to the module created by the parent Gene
	 */
	virtual void createSystem(OptimSystem& sys, systemModule::ModulePtr parent=NULL)const=0;

	///Returns the Gene node value. That is the underlying propModule object
	virtual propModule& getNode()const;

	///Returns the Gene type
	virtual MyGene::Type type()const=0;

	///Return the type of the node
	virtual MyNode::Type nodeType()const=0;

	virtual std::vector<Type> argsType()const=0;

	///Returns the value of the gene
	virtual int value()const;

	///Return the name of the gene
	virtual std::string name()const;

	///Check that the type rules are satisfied
	void checkValidity()const;

	///Returns the type ID of this class. Used for serialisation
	virtual int isA ();

	MyGene** findNthNode (MyGene** rootPtr, int findFunction,
				      int &iLengthCount, MyGene::Type type);

	MyGene** choose (MyGene** rootPtr, MyGene::Type type);

	MyGene** choose (MyGene** rootPtr){return (MyGene**)GPGene::choose((GPGene**)rootPtr);}

	int length (MyGene::Type type);

	void reserveSpace (int numObjects);

	void resolveNodeValues (GPTypedGeneSet& ns);

	virtual void save (ostream& os);

	///Create GUI tree representing the sub-tree starting from this gene
	/**
	 * @param m_refTreeModel Reference to the GUI tree object
	 * @param m_Columns Reference to the column model used by the GUI tree
	 * @param parent parent row of this tree node
	 * @param ID reference to the current column ID
	 */
	void createTree(Glib::RefPtr<Gtk::TreeStore> m_refTreeModel,
			const ExampleWindow::ModelColumns& m_Columns,
			Gtk::TreeModel::Row& parent, int& ID)const;

	///Create GUI tree representing the tree starting from this gene. This gene is considered the root of the whole tree.
	/**
	 * @param Reference to the GUI tree object
	 * @param m_Columns Reference to the column model used by the GUI tree
	 */
	void createTree(Glib::RefPtr<Gtk::TreeStore> m_refTreeModel,
			const ExampleWindow::ModelColumns& m_Columns)const;
};
}
#endif /* MYGENE_H_ */
