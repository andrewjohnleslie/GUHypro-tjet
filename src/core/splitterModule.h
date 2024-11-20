//
// Created by robert on 24/11/17.
//

#ifndef HYPRO_SPLITTERMODULE_H
#define HYPRO_SPLITTERMODULE_H

#include "core/propModule.h"

namespace hypro {
    class splitterModule : public propModule {
    protected:
        Node *N3_; ///< Third node. It is the second exit node besides 'N2_'

    public:
        splitterModule();

        splitterModule(std::string name, Node *N1, Node *N2, Node *N3, int ID = defID_);

        splitterModule(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

        splitterModule(const splitterModule &mod);

        virtual ~splitterModule();

        /// Set N3_
        virtual void N3(Node &N3);

        /// Get N3_
        Node &N3() const;

        // TODO outNode needed?
//        virtual std::vector<const Node *> inNode() const;

        virtual std::string typeName() const;

        virtual void serialize(Archive& ar) const;
		virtual void unserialize(const Archive& ar);

		template<class Archive>
		void save(Archive &ar, const unsigned int version) const;

		template<class Archive>
		void load(Archive &ar, const unsigned int version);
    };
}

#endif //HYPRO_SPLITTERMODULE_H
