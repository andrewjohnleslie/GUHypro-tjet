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

#include "MyGene.h"
#include <core/propModule.h>
#include <ParMixer.h>
#include <nonLinear.h>
#include <MyPopulation.h>
#include <ModuleAdditionGene.h>
#include <GPTypedGeneSet.h>

namespace hypro {
int printTexStyle=0;

MyGene::MyGene (GPNode& gpo)
:
	GPGene(gpo),
	verbosity_(0){
}

MyGene::MyGene(const MyGene& gpo)
:
		GPGene(gpo),
		verbosity_(gpo.verbosity_){
}

MyGene& MyGene::operator= (const MyGene& gpo){
	verbosity_ = gpo.verbosity_;
	node = gpo.node;
	return *this;
}

MyGene::~MyGene(){
}

void MyGene::createSystem(OptimSystem& sys, systemModule::ModulePtr parent)const{
	throw std::logic_error("Error: not implemented");
}

propModule& MyGene::getNode()const{
	return *(propModule*)node;
}

void MyGene::printMathStyle (ostream& os, int lastPrecedence){
	MyNode& node_ = *(MyNode*)node;
	// Function or terminal?
	if (isFunction()){
		os << node_.Name_ << "(";
		for(int i=0; i<containerSize(); i++){
			NthMyChild(i)->printMathStyle(os);
			if(i<containerSize()-1){
				os << ", ";
			}
			os << ")";
		}

		// Print the terminal
		if (isTerminal()) {
			os << node_.Name_;
		}
	}
}

// Print out a gene in LaTeX-style. Don't be confused, I don't make a
// difference whether this gene is the main program or an ADF, I
// assume the internal structure is correct.
	void MyGene::printTeXStyle(ostream &os, int lastPrecedence) {
		printMathStyle(os);
	}


// Print a Gene.
	void MyGene::printOn(ostream &os) {
		if (printTexStyle)
			printTeXStyle(os);
		else
			printMathStyle(os);
	}

//Reimplements the createChild function
GPGene* MyGene::createChild (GPNode& gpo){
//	return new MyGene (gpo);
	throw std::logic_error("Error: not implemented");
}

void MyGene::checkValidity()const{
	std::vector<Type> argsT(argsType());

	for(int i=0;i<containerSize(); i++){
		NthMyChild(i)->checkValidity();
		if(NthMyChild(i)->type() != argsT[i]){
			std::stringstream err;
			err << "Error: Type of child " << NthMyChild(i)->getNode().Name_
					<< " is wrong for Gene " << getNode().Name_ ;
			throw std::logic_error(err.str());
		}
	}
}

int MyGene::isA() {
		return MyPopulation::MyGeneID;
	}

void MyGene::createTree(Glib::RefPtr<Gtk::TreeStore> m_refTreeModel,
		const ExampleWindow::ModelColumns& m_Columns,
		Gtk::TreeModel::Row& parent, int& ID)const{
	Gtk::TreeModel::Row row = *(m_refTreeModel->append(parent.children()));
	row[m_Columns.m_col_name] = name();

		for (int i = 0; i < containerSize(); i++) {
			NthMyChild(i)->createTree(m_refTreeModel, m_Columns, row, ID);
		}
		row[m_Columns.m_col_id] = ID;
		ID++;
	}
void MyGene::createTree(Glib::RefPtr<Gtk::TreeStore> m_refTreeModel,
		const ExampleWindow::ModelColumns& m_Columns)const{
	Gtk::TreeModel::Row row = *(m_refTreeModel->append());
	row[m_Columns.m_col_name] = name();

	int ID = 0;
	for(int i=0; i<containerSize(); i++){
		NthMyChild(i)->createTree(m_refTreeModel,m_Columns,row,ID);
	}
	row[m_Columns.m_col_id] = ID;
}

MyGene** MyGene::choose (MyGene** rootPtr, MyGene::Type type){
	MyGene **pg=NULL;

	// Calculate the length of the subtree starting at the given gene
	int totalLength=(**rootPtr).length (type);

	if (totalLength==0) return NULL;

	// loop 10 times
	for (int i=0; i<10; i++){
		// Calculate a random number between 1..totalLength
		int iLengthCount = (GPrand() % totalLength) + 1;

		// Find gene with this value
		pg = findNthNode (rootPtr, 0, iLengthCount, type);

#if GPINTERNALCHECK
		if (!pg)
			GPExitSystem ("GPGene::choose", "Didn't find tree node");
#endif

		// If this pointer points to a function rather than terminal
		// return this value else keep going around the loop
		if ((**pg).isFunction ())
		  return pg;
	}

	// If after 10 loops still don't have function well return terminal
	return pg;
}

MyGene** MyGene::findNthNode (MyGene** rootPtr, int findFunction,
			      int &iLengthCount, MyGene::Type typ)
{
	if(*rootPtr != this){
		throw std::logic_error("Error: rootPtr must point to this object.");
	}

	// if we are looking for any node, or if we are looking just for
	// function nodes and this happens to be one, then ...
	if (type()==typ &&
			(!findFunction || (findFunction && containerSize()>0)))
	{
	  // ... decrement the length counter and return if it is zero
	  if (--iLengthCount<=0){
		  return rootPtr;
	  }

	}

	// Do the same for all children
	for (int n=0; n<containerSize(); n++){
		// Get pointer to pointer of n-th child
		MyGene** childPtr=(MyGene**) getPointerAddress (n);

		// If child exists
		if (*childPtr){
			// Recursive call for the children
			MyGene** found = (**childPtr).findNthNode (childPtr, findFunction, iLengthCount, typ);

			// If we found it, return immediately
			if (found)
				return found;
		}
	}

	// We didn't find the correct one so far
	return NULL;
}

int MyGene::length (MyGene::Type typ){
	int lengthSoFar;
	if(type() == typ){
		lengthSoFar = 1;   // Thats me!
	}else{
		lengthSoFar = 0;   // Not counting if another type
	}

	// Do same for all children, if there are any, and add up length
	MyGene* current;
	for (int n=0; n<containerSize(); n++)
	if ((current=(MyGene*)NthChild (n)))
	  lengthSoFar += current->length (typ);

	// Return the length
	return lengthSoFar;
}

void MyGene::create (enum GPCreationType ctype, int allowableDepth,
				GPNodeSet& ns){
	create(ctype, allowableDepth, (GPTypedGeneSet&)ns);
}

int MyGene::value()const{
	if(node){
		return getNode().value();
	}else{
		return nodeValue;
	}
}

std::string MyGene::name()const{
	return getNode().Name_;
}

int MyGene::begin(Type type)const{
	int i = -1;
	do{
		i++;
		if(i >= containerSize()){
			break;
		}
	}while(NthMyChild(i)->type() != type);
	return i;
}

int MyGene::end(Type type)const{
	int i = containerSize();
	do{
		i--;
		if(i < 0){
			break;
		}
	}while(NthMyChild(i)->type() != type);
	i++;

	return i;
}

void MyGene::reserveSpace (int numObjects)
{
  // We can't alloc space, if this is 0.  This is not an error!  A
  // terminal node, for example, doesn't have children.
  if (numObjects<=0)
    {
      contSize=0;
      return;
    }

  // Save the container size in object
  contSize=numObjects;

  // Alloc array of pointers and set them to NULL
  container=new GPObject*[numObjects];
  for (int i=0; i<numObjects; i++)
    container[i]=NULL;
}

void MyGene::save (ostream& os){
	os << value() << ' ';
	GPContainer::save (os);
}

void MyGene::resolveNodeValues (GPTypedGeneSet& ns){
  // Transform the node value to an address to the appropriate node
  *this = *ns.searchForGene (value());

  // Same for all children
  MyGene* current;
  for (int n=0; n<containerSize(); n++)
    if ((current=NthMyChild (n)))
      current->resolveNodeValues (ns);
}
}

