#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <TriangMesh.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(SWE_FVM, m) {
    // Optional docstring
    m.doc() = "My library";

    py::class_<Topology>(m, "Topology")
        .def(py::init<
                NodeTag,
                const EdgeTagArray&,
                const EdgeTagArray&,
                const TriangTagArray&,
                const TriangTagArray&,
                const TriangTagArray&
                >())
        .def_static("create", &Topology::create,
                "numNodes"_a.noconvert(),
                "edgeNodes"_a.noconvert(),
                "edgeElements"_a.noconvert(),
                "elementNodes"_a.noconvert(),
                "elementEdges"_a.noconvert(),
                "elementNeighbours"_a.noconvert())
        .def("num_nodes", &Topology::NumNodes)
        .def("num_edges", &Topology::NumEdges)
        .def("num_elements", &Topology::NumTriangles);
}

