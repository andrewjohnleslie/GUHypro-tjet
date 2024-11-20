//
// Created by robert on 24/11/17.
//

#ifndef HYPRO_SIMPLESPLITTER_H
#define HYPRO_SIMPLESPLITTER_H

#include "splitterModule.h"

namespace hypro {
    class simpleSplitter : public splitterModule {

    public:
        simpleSplitter();

        simpleSplitter(std::string name, Node *N1, Node *N2, Node *N3, int ID = defID_);

        // simpleSplitter(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap);

        simpleSplitter(const simpleSplitter &mod);

        virtual ~simpleSplitter();

        double calculate();

        propModuleSERIAL(simpleSplitter);
        virtual void serialize(Archive& ar) const;
	    virtual void unserialize(const Archive& ar);


    };
}

#endif //HYPRO_SIMPLESPLITTER_H
