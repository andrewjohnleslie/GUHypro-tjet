//
// Created by robert on 24/11/17.
//

#include "mixerModule.h"

namespace hypro {

    mixerModule::mixerModule(std::string name, Node *N1, Node *N2, Node *N3, int ID)
            :
            propModule(name, N1, N2, ID, 3),
            N3_(NULL) {
        this->N3(*N3);
    }

    mixerModule::mixerModule()
            :
            propModule(),
            N3_(NULL) {
    }

    mixerModule::mixerModule(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
            :
            propModule(is, nodeMap, moduleMap) {
        unsigned N3id;
        is >> N3id;
        N3_ = nodeMap.at(N3id).get();
    }

    mixerModule::mixerModule(const mixerModule &mod)
            :
            propModule(mod),
            N3_(mod.N3_) {
    }

    mixerModule::~mixerModule() {
        // TODO Auto-generated destructor stub
    }


    void mixerModule::N3(Node &N3) {
        N3_ = &N3;
        N3.downstream_ = this;
    }

    Node &mixerModule::N3() const {
        return *N3_;
    }

//    Glib::RefPtr<ModuleGraph> PropModule3::draw() {
//        graph_ = ModuleGraph3::create(this);
//        return graph_;
//    }

    std::vector<const Node *> mixerModule::inNode() const {
        std::vector<const Node *> out = propModule::inNode();
        out.push_back(N3_);
        return out;
    }


    // void mixerModule::serialize(std::ostream &os) const {
    //     mixerModule::serialize(os);

    //     os << " " << N3_->ID_;
    // }

    // template<class Archive>
    // void mixerModule::save(Archive &ar, const unsigned int version) const {
    //     boost::serialization::void_cast_register<mixerModule, propModule>();
    //     ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(propModule);

    //     Node *n3 = &N3();
    //     ar & boost::serialization::make_nvp("N3", n3);
    // }

    // template<class Archive>
    // void mixerModule::load(Archive &ar, const unsigned int version) {
    //     ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(propModule);

    //     Node *n3 = NULL;
    //     ar & boost::serialization::make_nvp("N3", n3);
    //     N3(*n3);
    // }

    // template void mixerModule::save<boost::archive::xml_oarchive>(
    //         boost::archive::xml_oarchive &ar, const unsigned int version) const;

    // template void mixerModule::load<boost::archive::xml_iarchive>(
    //         boost::archive::xml_iarchive &ar, const unsigned int version);
}
