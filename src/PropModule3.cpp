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

#include "Archive.h"
#include "PropModule3.h"
#include "Graphics/ModuleGraph3.h"

namespace hypro {

	PropModule3::PropModule3(std::string name, Node *N1, Node *N2, Node *N3, int ID)
			:
			balanceMach(name, N1, N2, ID, 3),
			N3_(N3) {
	}

	PropModule3::PropModule3()
			:
			balanceMach(),
			N3_(NULL) {
	}

	PropModule3::PropModule3(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			balanceMach(is, nodeMap, moduleMap) {
		unsigned N3id;
		is >> N3id;
		N3_ = nodeMap.at(N3id).get();
	}

	PropModule3::PropModule3(const PropModule3 &mod)
			:
			balanceMach(mod),
			N3_(mod.N3_) {
	}

	PropModule3::~PropModule3() {
		// TODO Auto-generated destructor stub
	}

	void PropModule3::N3(Node &N3) {
		N3_ = &N3;
		N3.downstream_ = this;
	}

	Node &PropModule3::N3() const {
		return *N3_;
	}

	Glib::RefPtr<ModuleGraph> PropModule3::draw() {
		graph_ = ModuleGraph3::create(this);
		return graph_;
	}

	std::vector<const Node *> PropModule3::inNode() const {
		std::vector<const Node *> out = propModule::inNode();
		out.push_back(N3_);
		return out;
	}

	void PropModule3::serialize(Archive& ar) const{
		propModule::serialize(ar);

		ar.put("N3", N3_->ID_);
	}

	void PropModule3::unserialize(const Archive& ar) {
		propModule::unserialize(ar);

		N3(*ar.getNodeRef("N3"));
	}

	void PropModule3::inNode(std::size_t i, Node& N){
		propModule::inNode(i, N);
		if(i == 1){
			N3(N);
		}
	}
}
