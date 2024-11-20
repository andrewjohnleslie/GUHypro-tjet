//
// Created by robert on 24/11/17.
//

#ifndef HYPRO_MIXERMODULE_H
#define HYPRO_MIXERMODULE_H

#include "core/propModule.h"

namespace hypro {
    // TODO Modified ModuleGraph3 for generic mixerModule class

    class mixerModule : public propModule {
    protected:
        Node *N3_; ///< Third node. It is the second input node beside node 'N1_'

    public:
        mixerModule();

        mixerModule(std::string name, Node *N1, Node *N2, Node *N3, int ID = defID_);

        mixerModule(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

        mixerModule(const mixerModule &mod);

        virtual ~mixerModule();

        /// Set N3_
        virtual void N3(Node &N3);

        /// Get N3_
        Node &N3() const;

        virtual std::vector<const Node *> inNode() const;

//        virtual Glib::RefPtr<ModuleGraph> draw();

        // virtual void serialize(std::ostream &os) const;

        // template<class Archive>
        // void save(Archive &ar, const unsigned int version) const;

        // template<class Archive>
        // void load(Archive &ar, const unsigned int version);

        // BOOST_SERIALIZATION_SPLIT_MEMBER()
    };

}
#endif //HYPRO_MIXERMODULE_H
