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

#include "OptimSystem.h"
#include <ConvNozzle.h>
#include <GPTypedGeneSet.h>
#include <ModuleAdditionGene.h>
#include <MyGene.h>
#include <MyPopulation.h>
#include <fenv.h>
#include <fstream>

namespace hypro
{
// The TeX-file
ofstream tout;
int OptimSystem::printTexStyle = 0;
bool OptimSystem::checkValidity_ = false;

OptimSystem::OptimSystem() : systemModule(defName_, (Node *)NULL, NULL), GP() {}

OptimSystem::OptimSystem(std::string name, Node *N1, Node *N2, int genes)
		: systemModule(name, N1, N2), GP(genes) {}

OptimSystem::OptimSystem(const OptimSystem &module)
		: systemModule(module), GP(module), message_(module.message_) {}

OptimSystem::~OptimSystem()
{
	// TODO Auto-generated destructor stub
}

GPObject &OptimSystem::duplicate() { return *((GP *)(new OptimSystem(*this))); }

GPGene *OptimSystem::createGene(GPNode &gpo)
{
	return new ModuleAdditionGene(gpo);
}

MyGene *OptimSystem::NthMyGene(int n) { return (MyGene *)GPContainer::Nth(n); }

MyGene *OptimSystem::NthMyGene(int n) const
{
	return (MyGene *)GPContainer::Nth(n);
}

void OptimSystem::evaluate()
{
	try
	{
		nodes_.clear();
		modules_.clear();

		NthMyGene(0)->createSystem(*this);
		N2_ = nodes_.back().get();

		unreduce();
		calculate();

		/**For now we define tinhe fitness function as the specific consumption.
* Note the it satisfy the requirement the the smaller the better.
*/
		if (N2_->M() < 1 && N2_->getPress() != N1_->getPress())
		{
			NodePtr N2(new Node(*N2_));
			nodes_.push_back(N2);
			modules_.back()->N2(*N2);
			N2_->Amax_ = std::numeric_limits<double>::infinity();
			N2_->Amin_ = 0.0;
			ModulePtr nozzle(
					new ConvNozzle("AutoNozzle", &modules_.back()->N2(), N2_));
			ConvNozzle *nozzleC = (ConvNozzle *)nozzle.get();
			nozzleC->freeStream_ = N1_;
			N2_->Name_ = "Nozz-out";

			add(nozzle);

			calculate();
		}

		calcFitness();
	}
	catch (std::runtime_error &e)
	{
		stdFitness = std::numeric_limits<double>::infinity();
		message_ = e.what();
	}
	catch (std::logic_error &e)
	{
		std::cerr << e.what() << std::endl
							<< std::endl
							<< "Logic error during evaluation of " << *(GP *)this;
		exit(1);
	}
	fitnessValid = 1;
}

void OptimSystem::calcFitness()
{
	std::vector<double> Xmfr = propellantMfr();
	double mfr = 0.0;
	for (std::size_t i = 0; i < Xmfr.size(); i++)
	{
		mfr += Xmfr[i];
	}
	double T = thrust();

	if (T > 0 && mfr > 0)
	{
		stdFitness = mfr * 9.81 / thrust();
		message_ = "";
	}
	else
	{
		stdFitness = std::numeric_limits<double>::infinity();
		message_ = "The thrust or the mfr is negative.";
	}
}

// Print a GP. If we want a LaTeX-output, we must provide for the
// equation environment, otherwise we simply call the print function
// of our base class.
void OptimSystem::printOn(ostream &os)
{
	// If we use LaTeX-style, we provide here for the right equation
	// overhead used for LaTeX.
	if (printTexStyle)
	{
		tout << "\\begin{eqnarray}" << endl;

		// Print all ADF's, if any
		GPGene *current;
		for (int n = 0; n < containerSize(); n++)
		{
			if (n != 0)
				os << "\\\\" << endl;
			os << "f_" << n + 1 << " & = & ";
			if ((current = NthGene(n)))
				os << *current;
			else
				os << " NONE";
			os << "\\nonumber ";
		}
		tout << endl
				 << "\\end{eqnarray}" << endl
				 << endl;
	}
	else
		// Just call the print function of our base class to do the
		// standard job.
		GP::printOn(os);
}

void OptimSystem::printOn() { printOn(std::cout); }

GPContainer &OptimSystem::cross(GPContainer *parents, int maxdepthforcrossover)
{
#if GPINTERNALCHECK
	// We are conservative: Only two sexes allowed
	if (parents->containerSize() != 2)
		GPExitSystem("GP::cross", "Only two parents allowed for crossover");
#endif

	// Get pointers of mum and dad from container
	OptimSystem &dad = *(OptimSystem *)(GP *)parents->Nth(0);
	OptimSystem &mum = *(OptimSystem *)(GP *)parents->Nth(1);

#if GPINTERNALCHECK
	// GP's certainly must have the same number of trees
	if (dad.containerSize() != mum.containerSize())
		GPExitSystem("GP::cross", "Mum and Dad must have same number of trees");
	if (dad.containerSize() == 0)
		GPExitSystem("GP::cross", "Parents contain no trees");
#endif

	// Work out which adf branch we are going to cut from
	int randTree = GPrand() % dad.containerSize();

	// Get the adresses of the pointers to the root genes of the ADF
	// branch we are going to cut from
	MyGene **rootGene1 = (MyGene **)dad.getPointerAddress(randTree);
	MyGene **rootGene2 = (MyGene **)mum.getPointerAddress(randTree);

#if GPINTERNALCHECK
	if (!*rootGene1 || !*rootGene2)
		GPExitSystem("GP::cross", "Genetic tree of Mum or Dad is NULL");
#endif

	// Loop around this as long as we find two points on the trees so
	// that the maxdepthforcrossover is not exceeded
	int maxDepth1, maxDepth2;
	do
	{
		// Determine the cut points by choosing a node within mum and
		// dad
		MyGene **cutPoint2 = NULL;
		MyGene **cutPoint1 = NULL;
		for (int i = 0; i < 100; i++)
		{
			cutPoint1 = (**rootGene1).choose(rootGene1);
			cutPoint2 = (**rootGene2).choose(rootGene2, (*cutPoint1)->type());

			if (cutPoint2)
				break;
		}
		if (!cutPoint2)
			throw std::runtime_error(
					"Error: Could not find two genes of the same type in mum and dad.");

		// Swap the whole subtrees.  Easy, isn't it? And so fast...
		MyGene *tmp = *cutPoint1;
		*cutPoint1 = *cutPoint2;
		*cutPoint2 = tmp;

		// Here the maximum depth of the new trees and only that trees
		// is calculated as other trees we assume to be under the
		// maximum depth of crossover
		maxDepth1 = (**rootGene1).depth();
		maxDepth2 = (**rootGene2).depth();

		// Make sure that maximum depth is not too high.  If so, swap
		// the subtrees back
		if (maxDepth1 > maxdepthforcrossover || maxDepth2 > maxdepthforcrossover)
		{
			tmp = *cutPoint1;
			*cutPoint1 = *cutPoint2;
			*cutPoint2 = tmp;
		}
	} while (maxDepth1 > maxdepthforcrossover ||
					 maxDepth2 > maxdepthforcrossover);

	// After crossover, the fitness of the GP is no longer valid, so we
	// set the corresponding flag.  The length and depth has to be
	// recalculated as well.
	dad.fitnessValid = 0;
	mum.fitnessValid = 0;
	dad.calcLength();
	dad.calcDepth();
	mum.calcLength();
	mum.calcDepth();

	dad.checkValidity();
	mum.checkValidity();

	// We return the same container, so we don't have to allocate a new
	// one and delete the parents container
	return *parents;
}

void OptimSystem::swapMutation(GPAdfNodeSet &adfNs)
{
#if GPINTERNALCHECK
	if (containerSize() == 0)
		GPExitSystem("GP::swapMutation", "GP contains no trees");
#endif

	// Point to a particular adf genetic tree and get that node set for
	// that tree.
	int randtree = GPrand() % containerSize();
	GPTypedGeneSet &gs = *(GPTypedGeneSet *)adfNs.NthNodeSet(randtree);

	GPGene *rootGene = NthGene(randtree);
	if (rootGene)
	{
		// Select a gene on that branch, e.g. get the address of the
		// pointer that points to the gene
		MyGene **g =
				(MyGene **)(rootGene->choose((GPGene **)getPointerAddress(randtree)));

		// First of all we have to work out how many arguments this
		// function (if it is one) has.  There are two ways: Using
		// containerSize() of the Gene or have a look at the node.  The
		// user could have modified the number of children of a gene, so
		// we don't use the node information
		int args = (*g)->containerSize();

		// Try a certain number of times to find a node with different
		// node identification value than before.
		for (int i = 0; i < swapAttempts_; i++)
		{
			// Choose a node at random from node set with the same
			// number of arguments than the gene we chose
			MyGene *gene = gs.chooseGeneWithArgs(args, (*g)->type());
			if (gene)
			{
				/* Check that the node has a different identification value,
* and that it is of the same type of the gene to be swept.*/
				if ((gene->value() != (*g)->value()) &&
						(gene->argsType() == (*g)->argsType()))
				{
					for (int j = 0; j < args; j++)
					{
						gene->put(j, *(*g)->NthMyChild(j));
					}
					*g = gene;
					// TODO there might be memory leak here because the old gene is not
					// deleted
					break;
				}
			}
		}
	}
	checkValidity();
}

void OptimSystem::shrinkMutation()
{
#if GPINTERNALCHECK
	if (containerSize() == 0)
		GPExitSystem("GP::swapMutation", "GP contains no trees");
#endif

	// Select a particular adf genetic tree by random
	int randtree = GPrand() % containerSize();

	// Get root gene
	GPGene *rootGene = NthGene(randtree);
	if (rootGene)
	{
		// Select a function gene on that branch, e.g. get the address
		// of the pointer that points to the gene.
		GPGene **rootGenePtr = (GPGene **)getPointerAddress(randtree);
		MyGene **g = (MyGene **)rootGene->chooseFunctionNode(rootGenePtr);

		// If function node exists (it may happen that there is no
		// function node at all)
		if (g)
		{
			// Try a certain number of times to find a node with different
			// node identification value than before.
			for (int i = 0; i < swapAttempts_; i++)
			{
				// Choose one subtree (or child) of the chosen function gene
				int subTree = GPrand() % (**g).containerSize();
				MyGene *child = (**g).NthMyChild(subTree);

				/*Check that the selected subtree is compatible in type with the
* parent gene to be removed.
*/
				if (child->type() == (*g)->type())
				{
					// Set the pointer the parent uses to point to the child to
					// NULL, so that the child won't be deleted if the parent is
					// killed.
					GPGene **childPtr = (GPGene **)(**g).getPointerAddress(subTree);
					*childPtr = NULL;

					// Delete parent.
					delete *g;

					// Put the child on the position of the former parent
					*g = child;

					// Recalculate length and depth
					calcLength();
					calcDepth();

					break;
				}
			}
		}
	}
	checkValidity();
}

void OptimSystem::mutate(GPVariables &GPVar, GPAdfNodeSet &ns)
{
	if (GPRandomPercent(GPVar.SwapMutationProbability))
	{
		swapMutation(ns);
		fitnessValid = 0;
	}

	if (GPRandomPercent(GPVar.ShrinkMutationProbability))
	{
		shrinkMutation();
		fitnessValid = 0;
	}
}

void OptimSystem::create(enum GPCreationType ctype, int allowableDepth,
												 GPAdfNodeSet &adfNs)
{
#if GPINTERNALCHECK
	// Check argument
	if (ctype != GPGrow && ctype != GPVariable)
		GPExitSystem("GP::create", "Argument ctype must be GPGrow or GPVariable");
#endif

	// As the first node is always a function that is created here in
	// this routine, decrease allowableDepth
	allowableDepth--;

	// loop through each adf
	for (int adf = 0; adf < containerSize(); adf++)
	{
		// Set the node set pointer to the correct node set (most
		// interesting work is when these are different for each ADF)
		GPTypedGeneSet &ns = (GPTypedGeneSet &)*adfNs.NthNodeSet(adf);

		// Choose a function from function set which complies with
		// Genetic Programming book or code.  Create a new gene with that
		// function.
		MyGene &g = ns.chooseFunction(MyGene::MODULE);
		//      GPGene& g=*createGene (tempfunc);

		// Create tree structure
		g.create(ctype, allowableDepth, ns);

		// Now we put the child into the container.  Nice, isn't it?
		put(adf, g);
	}

	// Calculate length and depth
	calcLength();
	calcDepth();
}

int OptimSystem::isA() { return MyPopulation::OptimSystemID; }

// Load operation
char *OptimSystem::load(istream &is)
{
	// Load variables
	is >> fitnessValid;
	is >> stdFitness;
	if (is.fail())
	{
		std::string s;
		is.clear();
		is >> s;
		if (s == "inf")
		{
			stdFitness = std::numeric_limits<double>::infinity();
		}
	}

	// Load container
	char *errMsg = GPContainer::load(is);
	if (errMsg)
		return errMsg;

	// We didn't save the length and depth as we can calculate them now.
	// Anyway, it's safer.
	calcLength();
	calcDepth();

	// Return NULL
	return errMsg;
}

int OptimSystem::createTree(int argc, char *argv[])
{
	Glib::RefPtr<Gtk::Application> app =
			Gtk::Application::create(argc, argv, "org.gtkmm.example");

	ExampleWindow window;
	NthMyGene(0)->createTree(window.m_refTreeModel, window.m_Columns);

	// Shows the window and returns when it is closed.
	fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	int out = app->run(window);
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	return out;
}

void OptimSystem::checkValidity() const
{
	if (checkValidity_)
	{
		NthMyGene(0)->checkValidity();
	}
}

void OptimSystem::resolveNodeValues(GPAdfNodeSet &adfNs)
{
	MyGene *current;
	for (int n = 0; n < containerSize(); n++)
		if ((current = NthMyGene(n)))
			current->resolveNodeValues(*(GPTypedGeneSet *)adfNs.NthNodeSet(n));
}
}