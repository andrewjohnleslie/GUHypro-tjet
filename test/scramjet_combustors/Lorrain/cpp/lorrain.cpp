//
// Created by robert on 13/11/17.
//

#include "../lorrain.h"

//int main() {
int run_lorrain(int mode, double combeff_par, double split) {
/* Modes:
 *      Two mode parameters, geometry and combustion.
 *      Geometry:
 *          1: Geometry is injector - isolator - combustor,
 *          2: Geometry is injector - combustor
 *      Combustion:
 *          1: Combustion is perfect heat addition
 *          2: Combustion is chemical equilibrium
 *          3: Combustion is split perfect heat addition
 *          4: Combustion is split chemical equilibrium
 */

    /*  Select the mode
     *     i               EffComb         EquilComb       Split EffComb   Split EquilComb  SplitReactor
     *     Isolator        0               1               2               3                4
     *
     */
    int i = mode;

    // Initial Conditions
    double T = 779.456;
    double U = 2602.0;
    double p = 42010;
    Cantera::compositionMap inletComposition;
    inletComposition.emplace("O2", 0.227);
    inletComposition.emplace("N2", 0.75);
    inletComposition.emplace("H2", 0.023);

    // Initiate and set some variables
    int combustionMode;
    int numNodes;
    std::vector<double> A;
    double combustorFriction;
    double combustorLength;
    int combustorNode1;
    double isolatorFriction = 0;
    double isolatorLength = 0;
    double combeff;
    int variableArea = 1;
    double reactFrac;

    // Calculate the global equivalence ratio, used by injectionPhi
    double phi = ((1.0 * 0.4661) / (0.258 * 35.889)) / (1.0 / 8.0);

    if (i == 0) {
        // Isolator, EffComb
        numNodes = 3;
        A = {0.013, 0.013, 0.013};

        combustorFriction = 0.0035; // 0 // equivalent to balanceMach
        combustorLength = 0.25;
        combustorNode1 = 1;
        isolatorFriction = 0.0035; // 0 // equivalent to balanceMach
        isolatorLength = 0.15;
        combustionMode = 0;
        combeff = combeff_par;
    } else if (i == 1) {
        // Isolator, EquilComb
        if (variableArea == 1) {
            numNodes = 5;
            A = {0.089, 0.089, 0.089, 0.089, 0.1048};
        } else {
            numNodes = 4;
            A = {0.089, 0.089, 0.089, 0.089};
        }
        combustorFriction = 0.003695;
        combustorLength = 0.126;
        combustorNode1 = 1;
        isolatorFriction = 0.0036048;
        isolatorLength = 0.23;
        combustionMode = 1;
    } else if (i == 2) {
        // Isolator, Split EffComb
        reactFrac = split;
        if (variableArea == 1) {
            numNodes = 9;
            A = {0.089, 0.089, 0.089, 0.089 - 0.089 * reactFrac, 0.089 * reactFrac, 0.089 - 0.089 * reactFrac,
                 0.089 * reactFrac, 0.089, 0.1048};
        } else {
            numNodes = 8;
            //      N1,  [0]    [1]   [2]                       [3]              [4]                      [5]
            A = {0.089, 0.089, 0.089, 0.089 - 0.089 * reactFrac, 0.089 * reactFrac, 0.089 - 0.089 * reactFrac,
                 0.089 * reactFrac, 2};
        }
//        for (auto &massFlow : A){
//            std::cout << massFlow << ", ";
//        }
        combustorFriction = 0.003695;
        combustorLength = 0.126;
        combustorNode1 = 1;
        isolatorFriction = 0.0036048;
        isolatorLength = 0.23;
        combustionMode = 2;
        combeff = combeff_par;
    } else if (i == 3) {
        // Isolator Split EquilComb
        reactFrac = split;
        if (variableArea == 1) {
            numNodes = 9;
            A = {0.089, 0.089, 0.089, 0.089 - 0.089*reactFrac, 0.089*reactFrac, 0.089 - 0.089*reactFrac, 0.089*reactFrac, 0.089, 0.1048};
        } else {
            numNodes = 8;
            //      N1,  [0]    [1]   [2]                       [3]              [4]                      [5]
            A = {0.089, 0.089, 0.089, 0.089 - 0.089 * reactFrac, 0.089 * reactFrac, 0.089 - 0.089 * reactFrac,
                 0.089 * reactFrac, 0.089};
        }
//        for (auto &massFlow : A){
//            std::cout << massFlow << ", ";
//        }
        combustorFriction = 0.003695;
        combustorLength = 0.126;
        combustorNode1 = 1;
        isolatorFriction = 0.0036048;
        isolatorLength = 0.23;
        combustionMode = 3;
    } else if (i == 4) {
        // Isolator, Split Reactor
        numNodes = 8;
        reactFrac = split;
        //      N1,  [0]    [1]   [2]                       [3]              [4]                      [5]
        A = {0.089, 0.089, 0.089, 0.089 - 0.089*reactFrac, 0.089*reactFrac, 0.089 - 0.089*reactFrac, 0.089*reactFrac, 0.089};
        for (auto &massFlow : A){
            std::cout << massFlow << ", ";
        }
        combustorFriction = 0.003695;
        combustorLength = 0.126;
        combustorNode1 = 1;
        isolatorFriction = 0.0036048;
        isolatorLength = 0.23;
        combustionMode = 3;
    } else {
        std::cout << "Error" << std::endl;
    }

    std::cout << i << ", " << variableArea << ", " << combeff << ", " << reactFrac << ", ";
    lorrain(numNodes, A, combustorFriction, combustorLength, combustorNode1, combustionMode, isolatorFriction,
            isolatorLength, p, T, U, inletComposition, phi, combeff, variableArea);
//    std::cout << "\ni = " << i << ", variableDuct = " << variableArea << ", comfeff = " << combeff << std::endl;
}

int main(){
    std::vector<double> engineModes = {2};
    std::vector<double> combEffs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
//    std::vector<double> splits = {0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

//    for (auto &engineMode : engineModes ) {
//        for (auto &combEff : combEffs ){
//            for (auto & split : splits)
//            run_lorrain(engineMode, combEff, split);
//        }
//    }

    // for (auto &engineMode : engineModes ) {
    //     run_lorrain(engineMode, 0.4, 0);
    // }
//
//    for (auto &engineMode : engineModes ) {
//        for (auto & split : splits){
//            run_lorrain(engineMode, 0, split);
//        }
//    }
//
   for (auto &engineMode : engineModes){
       for (auto &combEff : combEffs) {
           run_lorrain(engineMode, combEff, 0.116);
       }
   }

    return 0;
}