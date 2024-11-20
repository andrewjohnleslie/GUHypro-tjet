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
#include "Node.h"
#include "NodeGraph.h"
#include "typedefs.h"
#include <core/propModule.h>
#include <InfNanFloat.h>

namespace hypro {
	std::string Node::defName_ = "Node";
	unsigned Node::nextID_ = 0;

	Node::Node(const std::string &infile, std::string id_)
			:
			thermoKineticState(infile, id_),
			graph_(),
			A_(0),
			Amin_(0),
			Amax_(std::numeric_limits<double>::infinity()),
			Name_(defName_),
			ID_(nextID_),
			upstream_(NULL),
			downstream_(NULL) {
		nextID_++;
	}

	Node::Node() :
			thermoKineticState(),
			graph_(),
			A_(0),
			Amin_(0),
			Amax_(std::numeric_limits<double>::infinity()),
			Name_(defName_),
			ID_(nextID_),
			upstream_(NULL),
			downstream_(NULL) {
		nextID_++;
	}

	Node::Node(std::istream &is, const std::string &infile, std::string Name, unsigned ID, std::string id_)
			:
			thermoKineticState(is, infile, id_),
			Name_(Name),
			ID_(ID),
			upstream_(NULL),
			downstream_(NULL) {
		is >> A_;

		try {
			is >> Amin_;
		} catch (std::ios_base::failure &e) {
			is.clear();
			std::string tmp;
			is >> tmp;
			if (std::strcmp(tmp.c_str(), "inf") == 0) {
				Amin_ = std::numeric_limits<double>::infinity();
			} else {
				throw e;
			}
		}

		try {
			is >> Amax_;
		} catch (std::ios_base::failure &e) {
			is.clear();
			std::string tmp;
			is >> tmp;
			if (std::strcmp(tmp.c_str(), "inf") == 0) {
				Amax_ = std::numeric_limits<double>::infinity();
			} else {
				throw e;
			}
		}
	}

	Node::Node(const Node &N)
			:
			thermoKineticState(N),
			graph_(),
			Amin_(N.Amin_),
			Amax_(N.Amax_),
			Name_(N.Name_),
			ID_(nextID_),
			upstream_(NULL),
			downstream_(NULL) {
		nextID_++;
		A(N.A());
	}

	Node::~Node() {
		// TODO Auto-generated destructor stub
	}

	Node &Node::operator=(Node N) {
		thermoKineticState::operator=(N);

		graph_ = Glib::RefPtr<NodeGraph>();
		Amin_ = N.Amin_;
		Amax_ = N.Amax_;
		A(N.A());
		Name_ = N.Name_;

		return *this;
	}

	double Node::mfr() const {
		return rho() * getU() * A_;
	}

	double Node::Nfr() const {
		//return mfr()/gas_.W();
		return mfr() / W();
	}

//Mass flow rate for each chemical specie
	vector<double> Node::mfrX() const {
		//Foam::List<double> Mfr(X_);
		std::vector<double> Mfr(X());
		//Mfr = Nfr()*Mfr;
		std::transform(Mfr.begin(), Mfr.end(), Mfr.begin(),
					   std::bind1st(std::multiplies<double>(), Nfr()));

		//for(int i=0; i<Mfr.size(); i++){
		//	Mfr[i] = specieGas_.speciesData()[i].W()*Mfr[i];
		//}
		vector<double> Molecular_weights(WX());
		for (unsigned i = 0; i < Mfr.size(); i++) {
			Mfr[i] = Molecular_weights[i] * Mfr[i];
		}

		return Mfr;
	}

//Mole flow rate for each chemical specie
	vector<double> Node::NfrX() const {
		//Foam::List<double> nfr(X_);
		//nfr = Nfr()*nfr;
		std::vector<double> nfr(X());
		std::transform(nfr.begin(), nfr.end(), nfr.begin(),
					   std::bind1st(std::multiplies<double>(), Nfr()));
		return nfr;
	}

	void Node::A(const double &A) {
		if (A > Amax_ || A < Amin_) {
			std::ostringstream er;
			er << "Error: A = " << A << " is out of range <" << Amin_ << "-" << Amax_ << "> at Node: " << Name_;
			throw std::range_error(er.str());
		}

		A_ = A;
	}

	double Node::A() const {
		return A_;
	}

	Glib::RefPtr<NodeGraph> Node::draw(double x, double y) {
		if (graph_) {
			throw std::logic_error("Error: Trying to re-draw an already drawn node.");
		}
		graph_ = NodeGraph::create(*this, x, y);
	return graph_;
	}

void Node::undraw(){
	if(graph_){
		graph_->remove();
		graph_ = Glib::RefPtr<NodeGraph>();
	}
}


	Glib::RefPtr<NodeGraph> Node::graph() const {
		if (graph_) {
			return graph_;
		} else {
			std::ostringstream er;
			er << "Error: Node " << Name_ << " has not been drawn yet.";
			throw std::logic_error(er.str());
		}
	}


//	void Node::printOut(std::ostream &os, bool label) const {
//		if (label) {
////			fmt::format(os, "{}", "Node");
//			fmt::print(os, "{:<20} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n", "Node", "A", "p", "T", "T0", "U", "M", "mfr", "I*A", "p0", "H");
//		}
//		fmt::print(os, "{:<20} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.3f} {:<12.2f} {:<12} {:<12} {:<12}\n", Name_.c_str(), A(), getPress(), getTemp(), T0(), getU(), M(), mfr(), I()*A(), p0(), H());
//	}

	void Node::printOut(std::ostream &os, bool label){
		if (label) {
//			fmt::format(os, "{}", "Node");
			fmt::print(os, "{:<20} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<14} {:<12} {:<12} {:<12} {:<12}\n", "Node", "A", "p", "p0", "T", "T0", "U", "M", "mfr", "Entropy", "H", "gamma", "Cp");
		}
		fmt::print(os, "{:<20} {:<12.4f} {:<12.2f} {:<12.2f} {:<12.6f} {:<12.2f} {:<12.5f} {:<12.5f} {:<14.10f} {:<12} {:<12.3f} {:<12} {:<12}\n", Name_.c_str(), A(), getPress(), getp0(), getTemp(), getT0(), getU(), M(), mfr(), s(),  H(), gamma(), Cp());
	}

propModule* Node::downstream()const{
	return downstream_;
}

propModule* Node::upstream()const{
	return upstream_;
}

void Node::serialize(Archive& ar) const{
	thermoKineticState::serialize(ar);

	ar.put("ID", ID_);
	ar.put("Name", Name_);
	ar.put("MassFlowRate", mfr());
	ar.put("Area", A_);
	ar.put("Area min", Amin_);
	ar.put("Area max", Amax_);
}

void Node::unserialize(const Archive& ar){
	thermoKineticState::unserialize(ar);

	ID_ = ar.get<int>("ID");
	Name_ = ar.get<std::string>("Name");
	Amin_ = ar.getDouble("Area min");
	Amax_ = ar.getDouble("Area max");
	A(ar.getDouble("Area"));
}

std::shared_ptr<Node> Node::unserialize(std::istream& istr,
										const std::string& infile, std::string id_)
{
	istr.exceptions ( std::istream::failbit | std::istream::badbit );

	std::string tName;
	istr >> tName;
	if(std::strcmp(tName.c_str(),"Node")!=0){
		throw std::runtime_error("Error: Expected Node found: " + tName);
	}

	unsigned ID;
	std::string Name;
	char quote;
	istr >> ID >> quote >> Name;
	while(Name.back()!='"' && !istr.eof()){
		std::string name;
		istr >> name;
		Name = Name + " " + name;
	}

	Name.pop_back();

	return std::shared_ptr<Node>(new Node(istr,infile,Name,ID, id_));
}
}
