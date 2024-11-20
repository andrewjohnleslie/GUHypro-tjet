/*!
 *  \brief      GP gene set for typed GP
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2017 Alessandro Mogavero
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

#ifndef TYPEDGENESET_H_
#define TYPEDGENESET_H_

#include <gp.h>
#include <map>
#include <MyGene.h>

namespace hypro {
class GPTypedGeneSet : public GPNodeSet
{
	std::list<std::unique_ptr<MyGene>> genes_;
	std::map<const MyGene::Type, std::vector<MyGene*>> funcByType_;
	std::map<const MyGene::Type, std::vector<MyGene*>> termByType_;

public:
	GPTypedGeneSet ();
	GPTypedGeneSet(int n);
	GPTypedGeneSet (const GPTypedGeneSet& gpo);

	virtual GPObject& duplicate () { return *(new GPTypedGeneSet(*this)); }

	virtual void putGene (MyGene* gpo);

	MyGene& chooseTerminal(MyGene::Type type);

	MyGene& chooseFunction(MyGene::Type type);

	MyGene* chooseGeneWithArgs (int args, MyGene::Type type);

	virtual int isA ();
//  virtual char* load (istream& is); TODO add implementation

	virtual GPObject* createObject() { return new GPTypedGeneSet; }

	virtual char* load(std::istream&){
		throw std::logic_error("Error: not implemented method.");
	}

	virtual void save(std::ostream&){
		throw std::logic_error("Error: not implemented method.");
	}

	virtual void printOn(std::ostream& os){
		GPNodeSet::printOn(os);
	}

	MyGene* searchForGene (int value);
};
}
#endif /* TYPEDGENESET_H_ */
