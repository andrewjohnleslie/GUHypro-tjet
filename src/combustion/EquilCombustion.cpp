//
// Created by robert on 11/04/17.
//

#include "EquilCombustion.h"
#include "solvers/balanceMach.h"
#include "solvers/Friction.h"

namespace hypro {
    template<class Delta>
    EquilCombustion<Delta>::EquilCombustion(std::string name, Node *N1, Node *N2):
            Delta(name, N1, N2, this->defID_, 2) {
        // this->argType_[1] = this->FEEDBACK;
    }


    template<class Delta>
    EquilCombustion<Delta>::EquilCombustion():
            Delta(this->defName_, NULL, NULL, this->defID_, 2) {
    }


    template<class Delta>
    EquilCombustion<Delta>::EquilCombustion(std::istream &is,
                                            const propModule::NodeMap &nodeMap, const propModule::ModuleMap &moduleMap)
            :
            Delta(is, nodeMap, moduleMap) {
    }


    template<class Delta>
    EquilCombustion<Delta>::EquilCombustion(const EquilCombustion &mod)
            :
            Delta(mod) {
    }

    template<class Delta>
    EquilCombustion<Delta>::~EquilCombustion() {
    }

    template<class Delta>
    void EquilCombustion<Delta>::changeComposition() {
//        this->N2_->equilibrium = true;
        this->N2_->setEquilibrium(true);
        this->N2_->X(this->N1_->X());
    }

    template<class Delta>
    double EquilCombustion<Delta>::calculate() {
        return Delta::calculate();
    }


    // template<class Delta>
    // void EquilCombustion<Delta>::serialize(std::ostream &os) const {
    //     Delta::serialize(os);
    // }

    // template<class Delta>
    // template<class Archive>
    // void EquilCombustion<Delta>::serialize(Archive &ar, const unsigned int version) {
    //     boost::serialization::void_cast_register<EquilCombustion, Delta>();
    //     ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Delta);
    // }


    template<>
    std::string EquilCombustion<balanceMach>::typeName() const {
        return "CombustionEquil";
    }

    template<>
    std::string EquilCombustion<Friction>::typeName() const {
        return "CombustionEquilFriction";
    }


//Template implementation
    template
    class EquilCombustion<balanceMach>;

    template
    class EquilCombustion<Friction>;
}
// BOOST_CLASS_EXPORT_IMPLEMENT(hypro::EquilCombustion<hypro::balanceMach>);
// BOOST_CLASS_EXPORT_IMPLEMENT(hypro::EquilCombustion<hypro::Friction>);