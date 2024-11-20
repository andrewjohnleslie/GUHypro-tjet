//
// Created by robert on 7/6/17.
//

#include <boost/python.hpp>
#include "Hyperion.h"

boost::python::tuple runHyperion(double p, double T){
    using namespace hypro;
    double Thrust;
    std::vector<double> Mfr;
    Collection::OutType out;
    double dynamic_pressure;
    vector<double> emissions;
    double total_flow;
    std::cout << "Here we are" << std::endl;
   // systemModule *sys = Hyperion(101325, 230, 2, 1, 1, 2, Thrust, Mfr, out, dynamic_pressure, emissions, total_flow);
    systemModule *sys = Hyperion(p, T, 7, 1, 2, 1, Thrust, Mfr, out, dynamic_pressure, emissions, total_flow);
    // p = 11664, T = 216.65
    systemModule::deleteAll(sys);
    return boost::python::make_tuple(Thrust, Mfr, out, dynamic_pressure, emissions, total_flow);
}

BOOST_PYTHON_MODULE(libHyperion_python)
{
    //using namespace boost::python;
    //boost::python::def("runHyperion", runHyperion, boost::python::args("p","T"));
    boost::python::def("runHyperion", runHyperion);
}