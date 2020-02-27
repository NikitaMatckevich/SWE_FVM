#pragma once
#include "TriangMesh.hpp" // see this for details
namespace Eigen {
  template<typename T, size_t N>
  struct ArrayHasher {
    using arg_type = Array<T, N, 1>;
    size_t operator()(arg_type const & a) const {
      std::hash<T> hasher;
      size_t res = 0;
      for (size_t i = 0; i < N; i++) {
        res = res * 31 + hasher(a[i]);
      }
      return res;
    }
  };
  template<typename T, size_t N>
  struct ArrayEqualTo {
    using arg_type = Array<T, N, 1>;
    bool operator()(arg_type const & l, arg_type const & r) const {
      Eigen::Array<bool, N, 1> b = l == r;
      for (int k = 0; k < l.SizeAtCompileTime; k++) {
        if (b[k] == false) return false;
      }
      return true;
    }
  };
}
using edgeHasher  = Eigen::ArrayHasher<node_tag, 2>;
using edgeEqualTo = Eigen::ArrayEqualTo<node_tag, 2>;

//Unstructured triangular mesh
class UnstructTriangMesh : public TriangMesh {
public:
  //Parse .msh files
  int read_from_gmsh(std::string const & filename) {
    //Create input stream
    std::ifstream msh(filename);

    //Check file consistency
    if (!msh.is_open()) {
      std::cerr << "Can't open mesh file " << filename << std::endl;
      return -1;
    }
    std::string line;
    std::getline(msh, line);
    if (line.find("$MeshFormat") != 0) {
      std::cerr << "File " << filename << " is not in MSH (ver 4) format"
        << std::endl;
      return -1;
    }
    double version;
    int file_type, data_size;
    msh >> version >> file_type >> data_size;
    if (version < 4.0) {
      std::cerr << "Mesh format in " << filename << " is too old " << version
        << std::endl;
      return -1;
    }
    if (file_type != 0) {
      std::cerr << "Mesh format in " << filename << " is binary. Not supported."
        << std::endl;
      return -1;
    }

    //Parse boundary types from $Entities section
    while (std::getline(msh, line) && line.find("$Entities") != 0);
    if (msh.eof()) {
      std::cerr << "EOF while looking for $Entities section in " << filename
        << std::endl;
      return -1;
    }
    int num_pts, num_curves;
    msh >> num_pts >> num_curves;
    for (node_tag i = 0; i < num_pts; i++) {
      std::getline(msh, line);
    }
    std::unordered_map<node_tag, int> boundary_types;
    for (node_tag i = 0; i < num_curves; i++) {
      int curve_tag;
      msh >> curve_tag;
      double min_x, min_y, min_z, max_x, max_y, max_z;
      msh >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z;
      int boundary_type;
      msh >> boundary_type;
      boundary_types[curve_tag] = boundary_type;
    }

    int num_entity_blocks;

    //Skip till $Nodes section
    while (std::getline(msh, line) && line.find("$Nodes") != 0);
    if (msh.eof()) {
      std::cerr << "EOF while looking for $Nodes section in " << filename
        << std::endl;
      return -1;
    }
    //  int num_pts;
    msh >> num_entity_blocks >> num_pts;
    if (num_pts >= std::numeric_limits<node_tag>::max()) {
      std::cerr << "node_tag is too small to handle " << filename << std::endl;
      return -1;
    }
    p_.resize(2, num_pts);
    node_tag node_id = 0;
    //Read
    for (node_tag i = 0; i < num_entity_blocks; i++) {
      int entity_tag, dim, parametric, nodes_in_block;
      msh >> entity_tag >> dim >> parametric >> nodes_in_block;
      for (node_tag j = 0; j < nodes_in_block; j++) {
        node_tag tag_node;
        msh >> tag_node;
        if (tag_node != ++node_id) {
          std::cerr << "Warning: nodes in mesh " << filename
            << " are not consequent " << std::endl;
        }
        if (tag_node > num_pts) {
          std::cerr << "Node tag overflow in " << filename << std::endl;
          return -1;
        }
        msh >> p_(0, tag_node - 1) >> p_(1, tag_node - 1);
        double z;
        msh >> z;
      }
    }

    //Skip till $Elements section
    while (std::getline(msh, line) && line.find("$Elements") != 0);
    if (msh.eof()) {
      std::cerr << "EOF while looking for $Elements section in " << filename
        << std::endl;
      return -1;
    }
    int num_elements;
    msh >> num_entity_blocks >> num_elements;
    if (num_elements >= std::numeric_limits<node_tag>::max()) {
      std::cerr << "node_tag is too small to handle " << filename << std::endl;
      return -1;
    }
    //Resizing is not appropriate (not all the elements are triangles)
    //Reserve memory for triangles
    triang_points_.resize(num_elements, NoChange);
    triang_edges_.resize(num_elements, NoChange);
    //Reserve memory for edges
    int num_edges = num_pts + num_elements - 1;
    if (num_edges >= std::numeric_limits<node_tag>::max()) {
      std::cerr << "node_tag is too small to handle " << filename << std::endl;
      return -1;
    }
    edge_points_.resize(num_edges, NoChange);
    edge_triangs_.resize(num_edges, NoChange);

    node_tag cur_triang = 0, cur_edge = 0;

    std::unordered_map<edge_tag, node_tag, edgeHasher, edgeEqualTo> edge_indices;
    for (node_tag i = 0; i < num_entity_blocks; i++) {
      int tag_entity, dim, type, num_elements;
      msh >> tag_entity >> dim >> type >> num_elements;
      if (type == 2) {// 3-node triangle in MSH format
        for (node_tag j = 0; j < num_elements; j++) {
          triang_tag pts, edges;
          node_tag tag_element;
          msh >> tag_element >> pts[0] >> pts[1] >> pts[2];
          for (int k = 0; k < 3; k++) pts[k]--;
          for (int k = 0; k < 3; k++) {
            edge_tag e = { pts[k], pts[(k + 1) % 3] };
            if (e[0] > e[1]) std::swap(e[0], e[1]);
            node_tag edge_idx;
            auto edge_iter = edge_indices.find(e);
            if (edge_iter == edge_indices.end()) {
              //new edge
              edge_indices[e] = edge_idx = cur_edge;
              edge_points_.row(cur_edge) = e;
              edge_triangs_.row(cur_edge++) << cur_triang, -1;
            }
            else {
              //known edge
              edge_idx = (*edge_iter).second;
              edge_indices.erase(edge_iter);
              edge_triangs_(edge_idx, 1) = cur_triang;
            }
            edges[k] = edge_idx;
          }
          triang_points_.row(cur_triang) = pts;
          triang_edges_.row(cur_triang++) = edges;
        }
      }
      else {
        //Ignore everything else. Skip num_elements lines + 1 for cur line
        for (node_tag j = 0; j < num_elements + 1; j++) {
          msh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
      }
    }
    //Now edge_indices contains only boundary edges, we can use it somehow
    edge_points_.conservativeResize(cur_edge, NoChange);
    edge_triangs_.conservativeResize(cur_edge, NoChange);
    triang_points_.conservativeResize(cur_triang, NoChange);
    triang_edges_.conservativeResize(cur_triang, NoChange);
    triang_triangs_.conservativeResize(cur_triang, NoChange);
    for (node_tag i = 0; i < num_triangles(); i++) {
      triang_tag triangs;
      for (int k = 0; k < 3; k++) {
        if (i == edge_triangs_(triang_edges_(i, k), 0))
          triangs[k] = edge_triangs_(triang_edges_(i, k), 1);
        else
          triangs[k] = edge_triangs_(triang_edges_(i, k), 0);
      }
      triang_triangs_.row(i) = triangs;
    }

    return 0;
  }
};
