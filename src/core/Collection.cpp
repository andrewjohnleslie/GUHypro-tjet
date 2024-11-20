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
#include "Collection.h"
#include "solvers/isentropicDuct.h"
#include "AdaptedInlet.h"
#include "AdaptedThroatInlet.h"
#include "NeutralLink.h"
#include "AdaptedNozzle.h"
#include "ClosedInlet.h"
#include "solvers/EffIsenDuct.h"
#include "inletModels/ScramjetInlet.h"
#include "inletModels/SupersonicInlet.h"
#include "injectors/InjectionPhi.h"
#include "Friction.h"
#include "combustion/EffComb.h"
#include "combustion/EquilCombustion.h"


#include <gtkmm.h>
#include <goocanvasmm.h>
#include <cairomm.h>
#include <NodeGraph.h>
#include <CollGraph.h>

#include <fenv.h>
//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/shared_ptr.hpp>
//#include <boost/serialization/list.hpp>

namespace hypro{
Collection::Collection(std::string name,Node* N1,Node* N2,int ID)
:
	Collection(name,N1,N2,ID,1){
}

Collection::Collection(std::string name,Node* N1,Node* N2,int ID, int args)
:
	propModule(name,N1,N2,ID,args),
	nodes_(),
	mech_()
{
	N2point_ = 0;
}

Collection::Collection():
	Collection(defName_, NULL, NULL, 0, 0){
}

Collection::Collection(std::istream& is, const NodeMap& nodeMap, const ModuleMap& moduleMap)
:
			propModule(is,nodeMap,moduleMap),
			nodes_(),
			mech_() {
	N2point_ = 0;

	NodeMap nodeMapInt;
	ModuleMap moduleMapInt;
	nodeMapInt[N1_->ID_] = nodeMap.at(N1_->ID_);
	nodeMapInt[N2_->ID_] = nodeMap.at(N2_->ID_);

	std::string tag;
	is >> tag;
	if(std::strcmp(tag.c_str(),"nodes")!=0){
		throw std::runtime_error("Error: expected nodes found: " + tag);
	}
	is >> tag;
	while(!is.eof()){
		try{
			NodePtr N = Node::unserialize(is,N1_->XbyNames()[1]); //TODO THIS IS DEFINITELY WRONG
			if(nodeMapInt.count(N->ID_)){
				std::cout << "Warning: Node ID " << N->ID_ << " already used. A new node will not be created." << std::endl;
				nodes_.push_back(nodeMapInt.at(N->ID_));
			}else{
				nodes_.push_back(N);
				nodeMapInt[N->ID_] = N;
			}
		}catch(std::runtime_error& e){
			break;
		}
	}

	is >> tag;
	if(std::strcmp(tag.c_str(),"modules")!=0){
		throw std::runtime_error("Error: expected modules found: " + tag);
	}
	is >> tag;
	while(!is.eof()){
		std::streampos pos = is.tellg();
		is >> tag;
		if(std::strcmp(tag.c_str(),"}")==0){
			break;
		}else{
			is.seekg(pos);
		}

		ModulePtr M = propModule::unserialize(is,nodeMapInt,moduleMapInt);
		add(M);
		moduleMapInt[M->ID_] = M;
	}
	
	//Set choking Feedback and connectivity
	for(std::list<ModulePtr>::iterator it=modules_.begin(); it!=modules_.end(); it++){
		if(!(*it)->chokingFeedbackID_){
			(*it)->chokingFeedback_ = NULL;
		}else{
			try{
				(*it)->chokingFeedback_ = moduleMapInt.at((*it)->chokingFeedbackID_);
			}catch(std::out_of_range& e){
				std::ostringstream er;
				er << "Error: Module ID " << (*it)->chokingFeedbackID_ << " not found.";
				throw std::range_error(er.str());
			}
		}
	}
}

Collection::Collection(const Collection& mod)
:
		propModule(mod),
		N1point_(mod.N1point_),
		N2point_(mod.N2point_),
		modules_(mod.modules_),
		nodes_(mod.nodes_),
		mech_(mod.mech_)

{
}

GPObject& Collection::duplicate (){
	Collection* newMod(new Collection(*this));

	duplicate(newMod);

	return *newMod;
}

void Collection::duplicate (Collection* newModPtr){
	Collection& newMod(*newModPtr);

	newMod.clear();

	Node* nodes[nodes_.size()];
	Node* newNodes[nodes_.size()];
	int i = 0;
	for(std::vector<NodePtr>::const_iterator it = nodes_.begin(); it != nodes_.end(); it++){
		NodePtr newN(new Node(**it));
		newMod.nodes_.push_back(newN);
		nodes[i] = it->get();
		newNodes[i] = newN.get();
		i++;
	}

	int j = 0;
	for(std::list<ModulePtr>::const_iterator it = modules_.begin(); it != modules_.end(); it++){
		ModulePtr module((propModule*)&(*it)->duplicate());
		newMod.modules_.push_back(module);

		//Connect modules of the new collection to the new nodes
		if (&module->N1()!=N1_) {
			//NodePtr p((NodePtr)&module->N1());
			Node** I = std::find(nodes, nodes + nodes_.size(), &module->N1());
			module->N1(**(I - nodes + newNodes));
		}
		if (&module->N2()!=N2_) {
			Node** I = std::find(nodes, nodes + nodes_.size(), &module->N2());
			module->N2(**(I - nodes + newNodes));
		}

		newMod.updateNpoint(module);
		j++;
	}
}

Collection::~Collection() {
	// TODO Auto-generated destructor stub
}

void Collection::updateNpoint(const ModulePtr& module){
	if(module->N1_==N1_){
		N1point_.push_back(&(module->N1_));
	}
	if(module->N2_==N2_){
		N2point_ = &(module->N2_);
	}
}

void Collection::add(ModulePtr module){
	modules_.push_back(module);
	updateNpoint(module);
}

template <class Module>
void Collection::add(std::string name, std::vector<int> nodes){
	std::vector<Node*> nod;
	for(unsigned i=0;i<nod.size();i++)
	{
		if(nodes[i]==-1){
			nod[i] = N1_;
		}else if(nodes[0]==-2){
			nod[i] = N2_;
		}else{
			nod[i] = nodes_[nodes[i]].get();
		}
	}
	ModulePtr mod(new Module(name,nod[0],nod[1]));

	add(mod);
}

template <class Module>
void Collection::addshaft(std::string name, std::vector<int> shaftID){
	std::vector<mechanicalLink*> mech;
	for(unsigned i=0;i<mech.size();i++)
	{
		if(shaftID[i]==-1){
			mech[i] = N1_;
		}else if(shaftID[0]==-2){
			mech[i] = N2_;
		}else{
			mech[i] = mech_[shaftID[i]].get();
		}
	}
	ModulePtr mod(new Module(name,mech[0],mech[1]));
	add(mod);
}

template void Collection::add<isentropicDuct>(std::string name, std::vector<int> nodes);


void Collection::exchange(const ModulePtr& prevMod, const ModulePtr& newMod){
	std::list<ModulePtr>::iterator I = std::find(modules_.begin(), modules_.end(), prevMod);
	*I = newMod;

	/*Remove old pointer to N1point_ vector, this is needed because updateNpoint
	 * always add a new element and does not delete anything
	 */
	if(&prevMod->N1()==N1_){
		N1point_.remove(&prevMod->N1_);
	}
	updateNpoint(newMod);
}

template <class Module>
const Collection::ModulePtr Collection::exchange(const ModulePtr& prevMod, std::string name){
	ModulePtr newMod;
	newMod = ModulePtr(new Module(name,&prevMod->N1(),&prevMod->N2()));
	exchange(prevMod, newMod);

	return newMod;
}
template const Collection::ModulePtr Collection::exchange<AdaptedInlet>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<AdaptedThroatInlet>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<NeutralLink>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<isentropicDuct>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<AdaptedNozzle<isentropicDuct> >(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<AdaptedNozzle<EffIsenDuct> >(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<EffIsenDuct>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<ClosedInlet>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<InjectionPhi<Friction> >(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<EffComb<Friction> >(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<EquilCombustion<Friction> >(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<SupersonicInlet>(const ModulePtr& prevMod, std::string name);
template const Collection::ModulePtr Collection::exchange<ScramjetInlet>(const ModulePtr& prevMod, std::string name);

void Collection::remove(const ModulePtr& module){
	modules_.remove(module);
	if(module->N1_==N1_){
		N1point_.remove(&module->N1_);
	}
	if(&(module->N2_)==N2point_){
		N2point_ = 0;
	}
}

void Collection::initNodes(int n,Node N){
	nodes_ = std::vector<NodePtr>(n,NodePtr());
	for(int i=0;i<n;i++){
		nodes_[i] = NodePtr(new Node(N));
	}
}


void Collection::initShaft(int n){
	mech_ = std::vector<MechPtr>(n,MechPtr());
	for(int i=0;i<n;i++){

		mech_[i] = MechPtr(new mechanicalLink());

	}
}


double Collection::calculate(){
	double rChok;
	for(std::list<ModulePtr>::iterator it = modules_.begin(); it != modules_.end(); it++){
		if(verbosity_>0){
			std::cout << "Calculating module: " << (*it)->Name_ << std::endl;
		}
		rChok = (*it)->calculate();
		if(rChok<0){
			break;
		}
	}
	return rChok;
}

void Collection::N1(Node& N1){
	propModule::N1(N1);
	for(std::list<Node**>::const_iterator it=N1point_.begin(); it!=N1point_.end(); it++){
		**it = N1_;
	}
}

void Collection::N2(Node& N2){
	propModule::N2(N2);
	if(N2point_!=0){
		*N2point_ = N2_;
	}
}

vector<double> Collection::propellantMfr()const{
	vector<double> Mfr(N1_->X().size(),0.0);
	vector<double> MfrI(Mfr);
	for(std::list<ModulePtr>::const_iterator it = modules_.begin(); it != modules_.end(); it++){
		MfrI = (*it)->propellantMfr();
		Mfr = Mfr + MfrI;
	}
	return Mfr;
}

void Collection::unreduce(){
	for(std::list<ModulePtr>::const_iterator it = modules_.begin(); it != modules_.end(); it++){
		(*it)->unreduce();
	}
}

void Collection::clear(){
	modules_.clear();
	nodes_.clear();
	mech_.clear();
	N1point_.clear();
}

void Collection::printOut(std::ostream& os)const{
	os << Name_ << ":" << std::endl;
	N1_->printOut(os,true);

	for(std::vector<NodePtr>::const_iterator it=nodes_.begin(); it!=nodes_.end(); it++){
		(*it)->printOut(os,false);
	}
	N2_->printOut(os,false);
}

/* Gives detailed output of the collection.
 * For each node it gives in this order: A(), p_, p0_, T_, T0_, M()
 */
Collection::OutType Collection::out()const{
	OutType out;

	for(std::list<ModulePtr>::const_reverse_iterator it=modules_.rbegin(); it!=modules_.rend(); it++){
		std::vector<double> outTmp;
		outTmp.push_back((*it)->N2_->A());
		outTmp.push_back((*it)->N2_->Amax_);
		outTmp.push_back((*it)->N2_->Amin_);
		outTmp.push_back((*it)->N2_->getPress());
		outTmp.push_back((*it)->N2_->getTemp());
		outTmp.push_back((*it)->N2_->getU());
		for(std::size_t i=0; i<(*it)->N2_->X().size(); i++){
			outTmp.push_back((*it)->N2_->X()[i]);
		}
		out.push_back(std::make_pair((*it)->N2_->Name_,outTmp));

		if((*it)->isCollection()){
			std::vector< std::pair< std::string, std::vector<double> > > out1(
					((Collection*)it->get())->out() );
			out.insert(out.end(),++out1.rbegin(),--out1.rend());
		}
	}
	std::vector<double> outTmp;
	outTmp.push_back(modules_.front()->N1_->A());
	outTmp.push_back(modules_.front()->N1_->Amax_);
	outTmp.push_back(modules_.front()->N1_->Amin_);
	outTmp.push_back(modules_.front()->N1_->getPress());
	outTmp.push_back(modules_.front()->N1_->getTemp());
	outTmp.push_back(modules_.front()->N1_->getU());
	for(std::size_t i=0; i<modules_.front()->N1_->X().size(); i++){
		outTmp.push_back(modules_.front()->N1_->X()[i]);
	}
	out.push_back(std::make_pair(modules_.front()->N1_->Name_,outTmp));

	OutType out1;
	out1.resize(out.size());
	std::reverse_copy(out.begin(), out.end(), out1.begin());

	return out1;
}



bool Collection::isCollection()const{
	return true;
}

void Collection::show(int argc, char* argv[]){
	Gtk::Main app(&argc, &argv);
	Goocanvas::init();

	std::string filename = "image.ps";
	double width = 2000;
	double height = 1000;
	Cairo::RefPtr<Cairo::PsSurface> surface =
			Cairo::PsSurface::create(filename, width, height);

	Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);

	Gtk::Window win;
	Goocanvas::Canvas m_canvas;

	win.set_title(Name_);
	m_canvas.set_size_request(1000, 480);
	m_canvas.set_bounds(0, 0, 2000, 1000);
	Glib::RefPtr<Goocanvas::Item> root = m_canvas.get_root_item();
	Glib::RefPtr<NodeGraph> Node2 = N2().draw(1900.0,800.0);
	root->add_child(Node2);

	root->add_child(draw());

	Gtk::ScrolledWindow* sw = Gtk::manage(new Gtk::ScrolledWindow());
	sw->add(m_canvas);
	sw->set_min_content_width(640);
	sw->set_min_content_height(480);
	Glib::RefPtr<Gtk::Adjustment> hor = sw->get_hadjustment();
	hor->set_value(hor->get_upper());
	Glib::RefPtr<Gtk::Adjustment> vert = sw->get_vadjustment();
	vert->set_value(vert->get_upper());
	win.add(*sw);
	win.show_all_children();

	fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	Gtk::Main::run(win);
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	cr->save();
	m_canvas.render(cr,1.0);
	cr->stroke();
	cr->restore();
	cr->show_page();
}

Glib::RefPtr<ModuleGraph> Collection::draw(){
	graph_ = CollGraph::create(this);
	return graph_;
}

double Collection::drag()const{
	double D = 0.0;
	for(std::list<ModulePtr>::const_iterator it = modules_.begin(); it != modules_.end(); it++){
		D = D + (*it)->drag();
	}

	return D;
}

//------------Tree Iterator------------------------

Collection::TreeIterator::TreeIterator(propModule* x, Collection& cont)
:
		p_(x),
		container_(cont),
		isPastEnd_(false),
		isPastBegin_(false),
		branch_(){
	if(!p_){
		throw std::logic_error("Error: A tree iterator cannot be a NULL pointer.");
	}
}

Collection::TreeIterator::TreeIterator(const TreeIterator& mit)
:
		p_(mit.p_),
		container_(mit.container_),
		isPastEnd_(mit.isPastEnd_),
		isPastBegin_(mit.isPastBegin_),
		branch_(mit.branch_){
}

Collection::TreeIterator& Collection::TreeIterator::operator=(const TreeIterator& rhs){
	if(&container_!=&rhs.container_){
		throw std::logic_error("Error: Assignment between two TreeIterators of two different containers.");
	}
	p_ = rhs.p_;
	isPastEnd_ = rhs.isPastEnd_;
	isPastBegin_ = rhs.isPastBegin_;
	branch_ = rhs.branch_;
	return *this;
}

Collection::TreeIterator& Collection::TreeIterator::operator++(){
	if(isPastBegin_){
		isPastBegin_ = false;
	}else if(isPastEnd_){
		throw std::logic_error("Error: Tried to increment a past the end TreeIterator.");
	}else{
		if(p_->N2().downstream()){
			p_ = p_->N2().downstream();
		}else{
			isPastEnd_ = true;
		}
	}
	return *this;
}

Collection::TreeIterator Collection::TreeIterator::operator++(int){
	TreeIterator tmp(*this);
	operator++();
	return tmp;
}

Collection::TreeIterator& Collection::TreeIterator::operator--(){
	if(isPastEnd_){
		isPastEnd_ = false;
	}else if(isPastBegin_){
		throw std::logic_error("Error: Tried to decrement a past the begin TreeIterator.");
	}else{
		if(branch_.size()==0){
			if(p_->inNode().size()>1){
				const std::vector<const Node*>& inNode(p_->inNode());
				for(std::vector<const Node*>::const_iterator it=inNode.begin(); it!=inNode.end(); it++){
					branch_.push_back(TreeIterator((*it)->upstream(),container_));
				}
				p_ = branch_.front().p_;
			}else{
				if(p_->N1().upstream()){
					p_ = p_->N1().upstream();
				}else{
					isPastBegin_ = true;
				}
			}
		}else{
			branch_.front()--;
			if(branch_.front().isPastBegin_){
				if(branch_.size()>1){
					branch_.pop_front();
				}else{
					isPastBegin_ = true;
				}
			}
			p_ = branch_.front().p_;
		}
	}
	return *this;
}

Collection::TreeIterator Collection::TreeIterator::operator--(int){
	TreeIterator tmp(*this);
	operator--();
	return tmp;
}

bool Collection::TreeIterator::operator==(const TreeIterator& rhs){
	return p_==rhs.p_ && isPastEnd_==rhs.isPastEnd_;
}

bool Collection::TreeIterator::operator!=(const TreeIterator& rhs){
	return !operator==(rhs);
}

propModule& Collection::TreeIterator::operator*(){
	if(isPastEnd_){
		throw std::logic_error("Error: Tried to dereference a past the end TreeIterator.");
	}
	return *p_;
}

propModule* Collection::TreeIterator::operator->(){
	return &(**this);
}

Collection::TreeIterator& Collection::TreeIterator::operator+=(const int& rhs){
	if(rhs>=0){
		int i = 0;
		while(i<rhs && !isPastEnd_){
			operator++();
			i++;
		}
	}else{
		std::vector<TreeIterator> branch;
		int i = 0;
		while(i>rhs && !isPastBegin_){
			operator--();
			i--;
		}
	}
	return *this;
}

bool Collection::TreeIterator::isPastBegin()const{
	return isPastBegin_;
}

Collection::TreeIterator& Collection::TreeIterator::operator-=(const int& rhs){
	return operator+=(-rhs);
}

Collection::TreeIterator Collection::begin(){
	return Collection::TreeIterator(modules_.front().get(),*this);
}

Collection::TreeIterator Collection::end(){
	Collection::TreeIterator it(modules_.back().get(),*this);
	it.isPastEnd_ = true;
	return it;
}

Collection::ModulePtr Collection::TreeIterator::getModulePtr()const{
	ModulePtr mod;
	for(std::list<ModulePtr>::iterator it=container_.modules_.begin(); it!=container_.modules_.end(); it++){
		if(it->get()==p_){
			mod = *it;
			break;
		}
	}
	return mod;
}

void Collection::serialize(Archive& ar)const{
	propModule::serialize(ar);

	ar.put("Nodes", nodes_);
	ar.put("Modules", modules_);
}

void Collection::unserialize(const Archive& ar){
	propModule::unserialize(ar);

	nodes_ = ar.getNodeVector("Nodes");
	modules_ = ar.getModuleList("Modules");
}

bool Collection::canChoke()const{
	for(std::list<ModulePtr>::const_iterator it=modules_.cbegin(); it!=modules_.cend(); it++){
		if((*it)->canChoke()){
			return true;
		}
	}

	return false;
}
}

//BOOST_CLASS_EXPORT_IMPLEMENT(Collection);


