#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include "Hyperion.h"

BOOST_PYTHON_MODULE(libpypro)
{
  boost::python::class_<std::map<std::string, double> >("compositionMap")
  .def(boost::python::map_indexing_suite<std::map<std::string, double> >() )
;

void    (Node::*setA)(const double &)              = &Node::A;
double  (Node::*getA)()   const          = &Node::A;

void (Node::*setTPXmap)(double, double, std::map<std::string, double>) = &Node::setTPX;
  boost::python::class_<Node>("Node", boost::python::init<const std::string &, const std::string & >())
  .def(boost::python::init<const Node&>())
  .def("mfr", &Node::mfr)
  .def("A", setA)
  .def("A", getA)
  .def("setTPX", setTPXmap)
  .def("report", &Node::report)
;

    struct propwrap : propModule, boost::python::wrapper<propModule> {
        propwrap() {
            propModule(std::string name, Node *N1, Node *N2, int ID = defID_,
                    int args = 1);
        }

        double calculate() {
            return this->get_override("calculate")();
        }

        std::string typeName() const{
            return this->get_override("typeName")();
        }
    };

  boost::python::class_<propwrap, boost::noncopyable>("propModule", boost::python::no_init)
          .def(boost::python::init<std::string, Node , Node , boost::python::optional< int , int>>())
          .def("calculate", boost::python::pure_virtual(&propModule::calculate))
          .def("typeName", boost::python::pure_virtual(&propModule::typeName));

}

