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

#include "utilities/Archive.h"
#include "propModule.h"
#include "ModuleGraph.h"
#include "Mixer.h"
#include "solvers/isentropicDuct.h"
#include "combustion/Combustion.h"
#include "injectors/InjectionPhi.h"
#include "ParMixer.h"
#include "NeutralLink.h"
#include "Rocket.h"
#include "Wedge.h"
#include "AdaptedInlet.h"
#include "AdaptedThroatInlet.h"
#include "chokedConv.h"
#include "ClosedInlet.h"
#include "Inlet.h"
#include "SupersonicInlet.h"
#include "AdaptedNozzle.h"
#include "ConvNozzle.h"
#include "systemModule.h"
#include "injectors/InjectionPlate.h"
#include "injectors/InjectionPlatePressure.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
#include "combustion/EquilCombustion.h"
#include "combustion/CombustionReactor.h"

//#include <boost/serialization/shared_ptr.hpp>

namespace hypro {
	const double propModule::toll_ = 1e-8;
	const double propModule::pi_ = 3.14159265359;
	double propModule::chokingMach_ = 0.9;
	double propModule::chokingToll_ = 0.01;
	const double propModule::g0_ = 9.81;
	const std::string propModule::defName_ = "Module";

	std::map<std::string, propModule::Factory> propModule::typeCreators_ = std::map<std::string, propModule::Factory>();
	propModule::Populator propModule::pop_ = propModule::Populator();
	unsigned propModule::nextID_ = 1;
	std::string propModule::indentation_ = "";

	propModule::Populator::Populator() {
//		typeCreators_["Combustion"] = &Combustion<balanceMach>::Create;
//		typeCreators_["CombFriction"] = &Combustion<Friction>::Create;
//		typeCreators_["InjectionPhi"] = &InjectionPhi<balanceMach>::Create;
//		typeCreators_["InjPhiFriction"] = &InjectionPhi<Friction>::Create;
//		typeCreators_["InjectionPlate"] = &InjectionPlate::Create;
//		typeCreators_["InjectionPlatePressure"] = &InjectionPlatePressure::Create;
//		typeCreators_["isentropicDuct"] = &isentropicDuct::Create;
//		typeCreators_["EffIsenDuct"] = &EffIsenDuct::Create;
//		typeCreators_["Mixer"] = &Mixer::Create;
//		typeCreators_["ParMixer"] = &ParMixer::Create;
//		typeCreators_["NeutralLink"] = &NeutralLink::Create;
//		typeCreators_["Rocket"] = &Rocket::Create;
//		typeCreators_["Wedge"] = &Wedge::Create;
//		typeCreators_["AdaptedInlet"] = &AdaptedInlet::Create;
//		typeCreators_["AdaptedThroatInlet"] = &AdaptedThroatInlet::Create;
//		typeCreators_["chokedConv"] = &chokedConv::Create;
//		typeCreators_["ClosedInlet"] = &ClosedInlet::Create;
//		typeCreators_["Inlet"] = &Inlet::Create;
//		typeCreators_["SupersonicInlet"] = &SupersonicInlet::Create;
//		typeCreators_["AdaptedNozzle"] = &AdaptedNozzle<isentropicDuct>::Create;
//		typeCreators_["AdaptedNozzEff"] = &AdaptedNozzle<EffIsenDuct>::Create;
//		typeCreators_["ConvNozzle"] = &ConvNozzle::Create;
//		typeCreators_["systemModule"] = &systemModule::Create;
//		typeCreators_["EquilCombustion"] = &EquilCombustion<balanceMach>::Create;
//		typeCreators_["EquilCombustionFriction"] = &EquilCombustion<Friction>::Create;
//       typeCreators_["separatedFlowCombustor"] = &separatedFlowCombustor::Create;
//        typeCreators_["CombustionReactor"] = &CombustionReactor::Create;

	}

	propModule::Populator::~Populator() {
	}

	propModule::propModule() {

	}

	propModule::propModule(std::string name, Node *N1, Node *N2, int ID, int args)
			:
			MyNode(ID, name, MODULE, args),
			graph_(),
			N1_(NULL),
			N2_(NULL),
			ID_(nextID_),
			verbosity_(0),
			chokingFeedback_(0),
			chokingFeedbackID_(0) {
		nextID_++;
		this->N1(*N1);
		this->N2(*N2);
	}


	propModule::propModule(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
			:
			MyNode(defID_, defName_, MODULE, 1),
			graph_() {
		char quote;
		is >> ID_ >> quote >> Name_;
		while (Name_.back() != '"' && !is.eof()) {
			std::string name;
			is >> name;
			Name_ = Name_ + " " + name;
		}
		Name_.pop_back();

		unsigned N1id, N2id, cFeed;
		is >> N1id >> N2id >> verbosity_;
		try {
			N1_ = nodeMap.at(N1id).get();
		} catch (std::out_of_range &e) {
			std::ostringstream er;
			er << "Error: Node ID " << N1id << " not found.";
			throw std::range_error(er.str());
		}
		try {
			N2_ = nodeMap.at(N2id).get();
		} catch (std::out_of_range &e) {
			std::ostringstream er;
			er << "Error: Node ID " << N2id << " not found.";
			throw std::range_error(er.str());
		}

		try {
			is >> cFeed;
			chokingFeedbackID_ = cFeed;
		} catch (std::ios_base::failure &e) {
			is.clear();
			std::string tmp;
			is >> tmp;
			if (std::strcmp(tmp.c_str(), "NULL") == 0) {
				chokingFeedbackID_ = 0;
			} else {
				throw e;
			}
		}
	}

	propModule::propModule(const propModule &mod)
			:
			MyNode(mod),
			graph_(),
			N1_(mod.N1_),
			N2_(mod.N2_),
			ID_(nextID_),
			verbosity_(mod.verbosity_),
			chokingFeedback_(mod.chokingFeedback_),
			chokingFeedbackID_(mod.chokingFeedbackID_) {
		nextID_++;
	}

	propModule::~propModule() {
		// TODO Auto-generated destructor stub
	}

	void propModule::N1(Node &N1) {
		N1_ = &N1;
		if (N1_ != NULL) {
			N1.downstream_ = this;
		}
	}

	Node &propModule::N1() const {
		return *N1_;
	}

	void propModule::N2(Node &N2) {
		N2_ = &N2;
		if (N2_ != NULL) {
			N2.upstream_ = this;
		}
	}

	Node &propModule::N2() const {
		return *N2_;
	}

	double propModule::isoffdesign() {
		return 1.0; // To be redefined in modules with offdesign features, otherwise left at 1.0 to indicate no off design considerations
	}

	double propModule::getT04() {
		throw std::runtime_error("Error: getT04 method undefined.");
		return -1;
	}

	void propModule::setT04(double& T) {
		throw std::runtime_error("Error: setT04 method undefined.");
	}


	void propModule::unreduce() {
	}

	bool propModule::isreduced() const {
		return false;
	}

	bool propModule::ischoked() const {
		return (N2_->M() > (chokingMach_ - chokingToll_)) && (N2_->M() < (chokingMach_ + chokingToll_));
	}

	void propModule::reduce(const double &par) {
		throw std::runtime_error("Error: reduce method undefined.");
	}

	double propModule::reduce() const {
		throw std::runtime_error("Error: reduce method undefined.");
		return -1;
	}

	bool propModule::isreduceable() const {
		return false;
	}

/**return the drag produced by the module, such as the intake drag.
 * it is 0 by default
 */
	double propModule::drag() const {
		return 0.0;
	}

/**Returns the mfr of propellants consumed by the module.
 * the mass flow rate is given for each chemical specie.
 */
	std::vector<double> propModule::propellantMfr() const {
		return std::vector<double>(N1_->WX().size(), 0.0);
	}

	Glib::RefPtr<ModuleGraph> propModule::draw() {
		graph_ = ModuleGraph::create(this);
		return graph_;
	}

	Glib::RefPtr<ModuleGraph> propModule::graph() const {
		return graph_;
	}

/**Returns true if the module is a collection of modules */
	bool propModule::isCollection() const {
		return false;
	}

	std::vector<const Node *> propModule::inNode() const {
		return std::vector<const Node *>(1, N1_);
	}

void propModule::inNode(std::size_t i, Node& N){
	if(i >= inNode().size()){
		throw std::logic_error("Error: index of input node out of range.");
	}
	if(i == 0){
		N1(N);
	}
}

void propModule::serialize(Archive& ar)const{
	ar.put("Type", typeName());
	ar.put("Name", Name_);
	ar.put("ID", ID_);
	ar.put("N1", N1_->ID_);
	ar.put("N2", N2_->ID_);
	ar.put("verbosity", verbosity_);
	if(chokingFeedback_.get()){
		ar.put("Choking Feedback", chokingFeedback_->ID_);
	}
}

void propModule::unserialize(const Archive& ar){
	ID_ = ar.get<int>("ID");
	Name_ = ar.get<std::string>("Name");
	N1(*ar.getNodeRef("N1"));
	N2(*ar.getNodeRef("N2"));
	verbosity_ = ar.get<int>("verbosity");
	try{
		 ar.getModuleRef("Choking Feedback", chokingFeedback_);
	}catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::property_tree::ptree_bad_path> >& e){
		chokingFeedback_ = NULL;
	}
}

std::shared_ptr<propModule> propModule::unserialize(std::istream& istr,
		const NodeMap& nodeMap,
		const ModuleMap& moduleMap){
	istr.exceptions ( std::istream::failbit | std::istream::badbit );

	std::string tName;
	istr >> tName;

	Factory tf;
	try{
		tf = typeCreators_.at(tName);
	}catch(std::out_of_range& e){
    std::cout << "in this propmodule place" << endl;
    throw std::range_error("Error: Unrecognized module type: " + tName);
	}
  std::cout << "In this propmodule place 2" << endl;
	return (*tf)(istr,nodeMap,moduleMap);
}
  

bool propModule::canChoke()const{
	return true;
}

bool propModule::inputFreeStream()const{
	return false;
}

void propModule::saveJson(std::ostream& os)const{
	Archive ar;
	ar.put("Node 1", &N1());
	ar.put("Model", this);
	ar.put("Node 2", &N2());

	ar.write_jason(os);
}

std::shared_ptr<propModule> propModule::loadJson(std::istream& is){
	Archive::initialize();

	Archive ar;
	ar.read_jason(is);

	ar.getNode("N1");
	ar.getNode("N2");
	auto out = ar.getPropModule("Model");

	Archive::finalize();

	return out;
}

void propModule::saveJson(std::string fname)const{
	std::ofstream ofs(fname);
	saveJson(ofs);
	ofs.close();
}

std::shared_ptr<propModule> propModule::loadJson(std::string fname){
	std::ifstream ifs(fname);
	if(ifs.fail()){
		throw std::runtime_error("Error: cannot open file " + fname);
	}

	auto out = loadJson(ifs);
	ifs.close();
	return out;
}
}
