//
// Created by robert on 24/11/17.
//
#include "Archive.h"
#include "simpleSplitter.h"

namespace hypro{

    simpleSplitter::simpleSplitter(std::string name, Node *N1, Node *N2, Node *N3, int ID)
            :
            splitterModule(name, N1, N2, N3, ID) {

    }

    simpleSplitter::simpleSplitter():
            splitterModule() {

    }

    simpleSplitter::simpleSplitter(std::istream &is, const NodeMap &nodeMap, const ModuleMap &moduleMap)
            :
            splitterModule(is, nodeMap, moduleMap) {

    }

    simpleSplitter::simpleSplitter(const simpleSplitter &mod)
            :
            splitterModule(mod){
    }

    simpleSplitter::~simpleSplitter() {

    }

    double simpleSplitter::calculate(){
//       std::cout << "Areas, N1_: " << N1_->A() << ", N2_: " << N2_->A() << ", N3_: " << N3_->A() << std::endl;
//        std::cout << N1_->report() << std::endl;
        N2_->setTPX(N1_->getTemp(), N1_->getPress(), N1_->X());
        N2_->setU(N1_->getU());
//        std::cout << N2_->report() << std::endl;
        N3_->setTPX(N1_->getTemp(), N1_->getPress(), N1_->X());
        N3_->setU(N1_->getU());
//        std::cout << N3_->report() << std::endl;
//        std::cout << "MFR, N1_: " << N1_->mfr() << ", N2_: " << N2_->mfr() << ", N3_: " << N3_->mfr() << std::endl;
//        std::cout << "U, N1_: " << N1_->getU() << ", N2_: " << N2_->getU() << ", N3_: " << N3_->getU() << std::endl;
        return 0;
    }

    void simpleSplitter::serialize(Archive& ar) const{
            splitterModule::serialize(ar);
    }

    void simpleSplitter::unserialize(const Archive& ar){
            splitterModule::unserialize(ar);
    }

}
