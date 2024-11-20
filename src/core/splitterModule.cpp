//
// Created by robert on 24/11/17.
//

#include "Archive.h"
#include "splitterModule.h"

namespace hypro {

    splitterModule::splitterModule(std::string name, Node *N1, Node *N2, Node *N3, int ID)
            :
            propModule(name, N1, N2, ID, 3),
            N3_(NULL) {
        this->N3(*N3);
    }

    splitterModule::splitterModule()
            :
            propModule(),
            N3_(NULL) {
    }

    splitterModule::splitterModule(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
            :
            propModule(is, nodeMap, moduleMap) {
        unsigned N3id;
        is >> N3id;
        N3_ = nodeMap.at(N3id).get();
    }

    splitterModule::splitterModule(const splitterModule &mod)
            :
            propModule(mod),
            N3_(mod.N3_) {
    }

    splitterModule::~splitterModule() {
        // TODO Auto-generated destructor stub
    }


    void splitterModule::N3(Node &N3) {
        N3_ = &N3;
        N3.upstream_ = this;
    }

    Node &splitterModule::N3() const {
        return *N3_;
    }

//    Glib::RefPtr<ModuleGraph> PropModule3::draw() {
//        graph_ = ModuleGraph3::create(this);
//        return graph_;
//    }

//    std::vector<const Node *> splitterModule::inNode() const {
//        std::vector<const Node *> out = propModule::inNode();
//        out.push_back(N3_);
//        return out;
//    }


    std::string splitterModule::typeName() const {

        std::string name = "splitterModule";
        return name;
    }

    void splitterModule::serialize(Archive& ar) const{
		propModule::serialize(ar);

		ar.put("N3", N3_->ID_);
	}

	void splitterModule::unserialize(const Archive& ar) {
		propModule::unserialize(ar);

		N3(*ar.getNodeRef("N3"));
	}
}
