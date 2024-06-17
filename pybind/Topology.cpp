#include <pybind11/pybind11.h>
#include <TriangMesh.h>

namespace py = pybind11;

PYBIND11_MODULE(SWE_FVM, m) {
    // Optional docstring
    m.doc() = "My library";

    py::class_<Topology>(m, "Topology")
        .def(py::init<
                const Eigen::Ref<NodeTagArray>&,
                const Eigen::Ref<EdgeTagArray>&,
                const Eigen::Ref<EdgeTagArray>&,
                const Eigen::Ref<TriangTagArray>&,
                const Eigen::Ref<TriangTagArray>&,
                const Eigen::Ref<TriangTagArray>&
                >())
        .def("num_nodes", &Topology::NumNodes)
        .def("num_edges", &Topology::NumEdges)
        .def("num_elements", &Topology::NumTriangles);
}

