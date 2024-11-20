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

#ifndef PROPMODULE3_H_
#define PROPMODULE3_H_

#include "solvers/balanceMach.h"

namespace hypro {
	class ModuleGraph3;

///Propulsive module with 3 nodes
/**The model has two inputs nodes `N1_` and `N3_`
 */
	class PropModule3 : public balanceMach {
	protected:
		Node *N3_; ///<Third node. It is the second input node beside node `N1_`

	public:
		PropModule3();

		PropModule3(std::string name, Node *N1, Node *N2, Node *N3, int ID = defID_);

		PropModule3(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

		PropModule3(const PropModule3 &mod);

		virtual ~PropModule3();

		///Assign `N3_`
		virtual void N3(Node &N3);

		///Get `N3_`
		Node &N3() const;

		virtual std::vector<const Node *> inNode() const;

		virtual void inNode(std::size_t i, Node& N);

		virtual Glib::RefPtr<ModuleGraph> draw();

		virtual void serialize(Archive& ar) const;
		virtual void unserialize(const Archive& ar);

		template<class Archive>
		void save(Archive &ar, const unsigned int version) const;

		template<class Archive>
		void load(Archive &ar, const unsigned int version);
	};
}

#endif /* PROPMODULE3_H_ */
