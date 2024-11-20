//
// Created by robert on 11/04/17.
//

#ifndef EQUILCOMBUSTION_H_
#define EQUILCOMBUSTION_H_

#include <solvers/balanceMach.h>

namespace hypro {
    class Friction;


    template<class Delta>
    class EquilCombustion : public Delta {
    protected:
        void changeComposition();

    public:
        EquilCombustion();

        EquilCombustion(std::string name, Node *N1, Node *N2);

        EquilCombustion(const EquilCombustion &mod);

        virtual ~EquilCombustion();

        double calculate();

        propModuleSERIALtemplate(EquilCombustion);

        // void serialize(std::ostream &os) const;

        // template<class Archive>
        // void serialize(Archive &ar, const unsigned int version);

    };
}

// BOOST_CLASS_EXPORT_KEY(hypro::EquilCombustion<hypro::balanceMach>);
// BOOST_CLASS_EXPORT_KEY(hypro::EquilCombustion<hypro::Friction>);

#endif /* EQUILCOMBUSTION_H_ */
