//
// Created by robert on 30/04/18.
//

#ifndef HYPRO_RETURNDATA_H
#define HYPRO_RETURNDATA_H

namespace hypro {
    struct outputData {
        double thrust;
        vector<double> Mfr;
        Collection::OutType out;
        double dynamic_pressure;
        vector<double> emissions;
        double total_flow;
        systemModule *sys;
        double kH;
        double Aratio;
        double T04;
    };
}
#endif //HYPRO_RETURNDATA_H
