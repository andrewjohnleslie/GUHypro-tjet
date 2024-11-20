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
#include "solvers/isentropicDuct.h"
#include "Mixer.h"
#include "simpleSplitter.h"
#include "inletModels/chokedConv.h"
#include "combustion/Combustion.h"
#include "injectors/InjectionPhi.h"
#include "ParMixer.h"
#include "NeutralLink.h"
#include "Rocket.h"
#include "Wedge.h"
#include "inletModels/AdaptedInlet.h"
#include "inletModels/AdaptedThroatInlet.h"
#include "inletModels/ClosedInlet.h"
#include "inletModels/Inlet.h"
#include "inletModels/SupersonicInlet.h"
#include "nozzleModels/AdaptedNozzle.h"
#include "nozzleModels/ConvNozzle.h"
#include "core/systemModule.h"
#include "injectors/InjectionPlate.h"
#include "solvers/EffIsenDuct.h"
#include "solvers/Friction.h"
#include "combustion/EffComb.h"
#include <sstream>

namespace hypro {
Archive::ModuleMap Archive::modules_ = Archive::ModuleMap();
Archive::NodeMap Archive::nodes_ = Archive::NodeMap();
std::vector<std::pair<Collection::ModulePtr&, unsigned> > Archive::moduleBuffer_ = std::vector<std::pair<Collection::ModulePtr&, unsigned> >();

std::map<std::string, Archive::Factory> Archive::typeCreators_ = std::map<std::string, Archive::Factory>();

Archive::Populator Archive::pop_ = Archive::Populator();
Archive::Populator::Populator(){
	typeCreators_["Combustion"] = &Combustion<balanceMach>::Create;
	typeCreators_["CombFriction"] = &Combustion<Friction>::Create;
	typeCreators_["InjectionPhi"] = &InjectionPhi<balanceMach>::Create;
	typeCreators_["InjPhiFriction"] = &InjectionPhi<Friction>::Create;
	typeCreators_["InjectionPlate"] = &InjectionPlate::Create;
	typeCreators_["isentropicDuct"] = &isentropicDuct::Create;
	typeCreators_["EffIsenDuct"] = &EffIsenDuct::Create;
	typeCreators_["Mixer"] = &Mixer::Create;
	typeCreators_["simpleSplitter"] =&simpleSplitter::Create;
	typeCreators_["ParMixer"] = &ParMixer::Create;
	typeCreators_["NeutralLink"] = &NeutralLink::Create;
	typeCreators_["Rocket"] = &Rocket::Create;
	typeCreators_["Wedge"] = &Wedge::Create;
	typeCreators_["AdaptedInlet"] = &AdaptedInlet::Create;
	typeCreators_["AdaptedThroatInlet"] = &AdaptedThroatInlet::Create;
	typeCreators_["chokedConv"] = &chokedConv::Create;
	typeCreators_["ClosedInlet"] = &ClosedInlet::Create;
	typeCreators_["Inlet"] = &Inlet::Create;
	typeCreators_["SupersonicInlet"] = &SupersonicInlet::Create;
	typeCreators_["AdaptedNozzle"] = &AdaptedNozzle<isentropicDuct>::Create;
	typeCreators_["AdaptedNozzEff"] = &AdaptedNozzle<EffIsenDuct>::Create;
	typeCreators_["ConvNozzle"] = &ConvNozzle::Create;
	typeCreators_["systemModule"] = &systemModule::Create;
	typeCreators_["Collection"] = &Collection::Create;
	typeCreators_["EffCombFriction"] = &EffComb<Friction>::Create;
}

Archive::Populator::~Populator(){
}

Archive::Archive():
			boost::property_tree::ptree(){
}

Archive::Archive(std::string s):
	boost::property_tree::ptree(s){
}

Archive::~Archive() {
	// TODO Auto-generated destructor stub
}

template<class T>
void Archive::operator& (std::pair<std::string, T> p){
	Archive a;
	p.last.save(a);
	push_back(boost::property_tree::ptree::value_type(p.first, a));
}

template<>
void Archive::operator&<double> (std::pair<std::string, double> p){
	std::ostringstream strs;
	strs << p.second;

	push_back(boost::property_tree::ptree::value_type(p.first, Archive(strs.str())));
}

template<class T>
void Archive::operator& (std::pair<std::string, std::vector<T> > p){
	Archive a;
	for (const auto it=p.second.begin(); it!=p.second.end(); it++){
		Archive ac;
		it->save(ac);
		a.push_back(boost::property_tree::ptree::value_type("", ac));
	}

	push_back(boost::property_tree::ptree::value_type(p.first, a));
}

template<class T>
Archive& Archive::operator<< (std::pair<std::string, T > p){
	std::ostringstream strs;
	strs << p.second;

	push_back(boost::property_tree::ptree::value_type(p.first, Archive(strs.str())));

	return *this;
}

template Archive& Archive::operator<< <std::string>(std::pair<std::string, std::string > p);

Archive& Archive::put(std::string name, const propModule* p){
	Archive a;
	p->serialize(a);
	push_back(boost::property_tree::ptree::value_type(name, a));

	return *this;
}

Archive& Archive::put (std::string name, propModule* p){
	const propModule* p_cnst = p;
	return put(name, p_cnst);
}

Archive& Archive::put (std::string name, Collection::ModulePtr p){
	const propModule* p_cnst = p.get();
	return put(name, p_cnst);
}

Archive& Archive::put (std::string name, const std::list<Collection::ModulePtr>& p){
	Archive a;
	for (auto it=p.begin(); it!=p.end(); ++it){
		a.put("", *it);
	}

	auto it = push_back(boost::property_tree::ptree::value_type(name, a));

	return (Archive&)it->second;
}

Archive& Archive::put (std::string name, const std::vector<double>& p){
	Archive a;
	int i = 0;
	for (auto it=p.begin(); it!=p.end(); ++it){
		std::stringstream ss;
		ss << i;
		a.put(ss.str(), *it);
		i++;
	}

	auto it = push_back(value_type(name, a));

	return (Archive&)it->second;
}

std::vector<double> Archive::getVector (std::string name)const{
	std::vector<double> out;

	int i = 0;
	const Archive& ar = get_child(name);
	for(auto it=ar.begin(); it!=ar.end(); it++){
		out.push_back(std::stod(it->second.data()));
		i++;
	}

	return out;
}

Archive& Archive::put(std::string name, const Node* p){
	Archive a;
	p->serialize(a);
	push_back(value_type(name, a));

	return *this;
}

Archive& Archive::put (std::string name, Node* p){
	const Node* p_cnst = p;
	return put(name, p_cnst);
}

Archive& Archive::put (std::string name, Collection::NodePtr p){
	const Node* p_cnst = p.get();
	return put(name, p_cnst);
}

Archive& Archive::put (std::string name, const std::vector<Collection::NodePtr>& p){
	Archive a;
	for (auto it=p.begin(); it!=p.end(); ++it){
		a.put("", *it);
	}

	auto it = push_back(value_type(name, a));

	return (Archive&)it->second;
}

Archive& Archive::put (std::string name, const std::array<std::vector<double>,2>& p){
	Archive a;
	for(auto it = p[0].begin(), it1 = p[1].begin(); it != p[0].end() && it1 != p[1].end(); ++it, ++it1){
		a.put(std::to_string(*it), *it1);
	}

	auto it = push_back(value_type(name, a));

	return (Archive&)it->second;
}

void Archive::write_jason(std::ostream& os)const{
	const boost::property_tree::ptree& pt = static_cast<const boost::property_tree::ptree&>(*this);
	boost::property_tree::json_parser::write_json(os, pt);
}

void Archive::read_jason(std::istream& in){
	boost::property_tree::ptree& pt = static_cast<boost::property_tree::ptree&>(*this);
	boost::property_tree::json_parser::read_json(in, pt);
}

Collection::ModulePtr Archive::getPropModule(std::string name)const{
	const Archive& ar = get_child(name);

	std::string typeName = ar.get<std::string>("Type");
	Collection::ModulePtr out;
	try{
		out = typeCreators_.at(typeName)();
	}catch(const std::out_of_range& oor){
		throw std::ios_base::failure("Error: Module type " + typeName + " not recognised.");
	}

	out->unserialize(ar);

	if(modules_.find(out->ID_) != modules_.end()){
		throw std::ios_base::failure("Error: Module ID " + std::to_string(out->ID_) + " already exists.");
	}
	modules_[out->ID_] = out;
	return out;
}

void Archive::getModuleRef(std::string name, Collection::ModulePtr& out)const{
	unsigned id = get<unsigned>(name);

	try{
		out = modules_.at(id);
	}catch(const std::out_of_range& oor){
		moduleBuffer_.push_back(std::pair<Collection::ModulePtr&, unsigned>(out, id));
	}
}

std::list<Collection::ModulePtr> Archive::getModuleList(std::string name)const{
	std::list<Collection::ModulePtr> out;

	const Archive& ar = get_child(name);
	for(auto it=ar.begin(); it!=ar.end(); it++){
		Archive a;
		a.push_back(value_type("Module", it->second));
		out.push_back(a.getPropModule("Module"));
	}

	return out;
}

const Archive& Archive::get_child(const std::string& path) const{
	return static_cast<const Archive&>(boost::property_tree::ptree::get_child(path));
}

Node* Archive::getNode(std::string name)const{
	Node* out = new Node();

	const Archive& ar = get_child(name);
	out->unserialize(ar);

	if(nodes_.find(out->ID_) != nodes_.end()){
		throw std::ios_base::failure("Error: Node ID " + std::to_string(out->ID_) + " already exists.");
	}
	nodes_[out->ID_] = out;

	return out;
}

Node* Archive::getNodeRef(std::string name)const{
	unsigned id = get<unsigned>(name);

	Node* out;
	try{
		out = nodes_.at(id);
	}catch(const std::out_of_range& oor){
		throw std::ios_base::failure("Error: Node ID " + std::to_string(id) + " does not exists.");
	}

	return out;
}

std::vector<Collection::NodePtr> Archive::getNodeVector(std::string name)const{
	std::vector<Collection::NodePtr> out;

	const Archive& ar = get_child(name);
	for(auto it=ar.begin(); it!=ar.end(); it++){
		Archive a;
		a.push_back(value_type("Node", it->second));
		out.push_back(Collection::NodePtr(a.getNode("Node")));
	}

	return out;
}

double Archive::getDouble(std::string name)const{
	try{
		return get<double>(name);
	}catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::property_tree::ptree_bad_data> >& e){
		std::string out = get<std::string>(name);
		if (out == "inf"){
			return std::numeric_limits<double>::infinity();
		}else{
			throw std::ios_base::failure("Error: " + out + " is not a number.");
		}
	}
}

std::array<std::vector<double>,2> Archive::getLookup(std::string name)const{
	std::array<std::vector<double>,2> out;

	const Archive& ar = get_child(name);
	for(auto it=ar.begin(); it!=ar.end(); it++){
		out[0].push_back(std::stod(it->first));
		out[1].push_back(std::stod(it->second.data()));
	}

	return out;
}

void Archive::finalize(){
	for(auto it=moduleBuffer_.begin(); it!=moduleBuffer_.end(); it++){
		try{
			it->first = modules_.at(it->second);
		}catch(const std::out_of_range& oor){
			throw std::ios_base::failure("Error: Module ID " + std::to_string(it->second) + " does not exist.");
		}
	}
}

void Archive::initialize(){
	modules_.clear();
	nodes_.clear();
	moduleBuffer_.clear();
}
}