#include "MAT.h"
#include <boost/python/numpy.hpp>
#include <iostream>

namespace p = boost::python;
using namespace boost::python;
namespace np = boost::python::numpy;


BOOST_PYTHON_MODULE(Matrix){
    class_<VEC>("VEC", init<int>())
        .def(init<int>())
        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / self)
        .def(self += self)
        .def(self -= self)
        .def(self *= self)
        .def(self /= self)
        .def(self + double())
        .def(self - double())
        .def(self * double())
        .def(self / double())
        .def(double() + self)
        .def(double() - self)
        .def(double() * self)
        .def(double() / self)
        .def("Print", &VEC::print)
    ;
    class_ <MAT>("MAT", init<int, int>())
        .def(init<int, int>())
        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self / self)
        .def(self += self)
        .def(self -= self)
        .def(self *= self)
        .def(self /= self)
        .def(self += double())
        .def(self -= double())
        .def(self *= double())
        .def(self /= double())
        .def(self + double())
        .def(self - double())
        .def(self * double())
        .def(self / double())
        .def(double() + self)
        .def(double() - self)
        .def(double() * self)
        .def(double() / self)
        .def("Print", &MAT::print)
    ;
}
