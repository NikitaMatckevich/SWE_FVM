#pragma once
#include "TriangMesh.hpp" // see this for details

//Structured triangular mesh
//Structure is as follows:
//   b = bottom triangle, see function triang_b
//   ...
//    ____________  _____________
//   |\           /|\           /|
//   | \         / | \    t    / |
//   |  \       /  |  \       /  |
//   |   \     /   |   \     /   |
//   |    \   /    |    \   /    |
//   |_____\ /_____|  l  \ /  r  |
//   |     /|\ t  /|     / \     | 
//   |    / | \  / |    /   \    |
//   |   /  |l \/ r|   /     \   |
//   |  /   |  /\  |  /       \  |
//   | /    | /b \ | /    b    \ |
//   |/_____|/____\|/___________\|


class StructTriangMesh : public TriangMesh {
public:
  StructTriangMesh(int Ni, int Nj, double h)
    : Ni_(Ni), Nj_(Nj)
    , num_verts_((Ni + 1) * (Nj + 1))
    , num_sqrs_(Ni      *  Nj) {
    int num_triangs = 4 * num_sqrs_;
    triang_points_.resize(num_triangs, NoChange);
    triang_edges_.resize(num_triangs, NoChange);
    triang_triangs_.resize(num_triangs, NoChange);
    int num_pts = num_verts_ + num_sqrs_;
    p_.resize(NoChange, num_pts);
    int num_edges = num_triangs + num_pts - 1;
    edge_points_.resize(num_edges, NoChange);
    edge_triangs_.resize(num_edges, NoChange);
    // fill p_ array with square's vertices
    for (int j = 0; j < Nj_ + 1; j++) {
      for (int i = 0; i < Ni_ + 1; i++) {
        p_(0, j*(Ni_ + 1) + i) = i * h;
        p_(1, j*(Ni_ + 1) + i) = j * h;
      }
    }
    // fill p_ array with square's centers
    for (int j = 0; j < Nj_; j++) {
      for (int i = 0; i < Ni_; i++) {
        p_(0, num_verts_ + j * Ni_ + i) = (i + 0.5) * h;
        p_(1, num_verts_ + j * Ni_ + i) = (j + 0.5) * h;
      }
    }
    // determine triangles and edges
    for (int j = 0; j < Nj_; j++) {
      for (int i = 0; i < Ni_; i++) {
        triang_b(i, j); //-bottom triangle of this square
        triang_r(i, j); //-right  triangle of this square
        triang_t(i, j); //-top    triangle of this square
        triang_l(i, j); //-left   triangle of this square
      }
    }
  }
  int Ni() const { return Ni_; }
  int Nj() const { return Nj_; }
private:
  node_tag curr_e_ = 0;
  int Ni_, Nj_, num_verts_, num_sqrs_;

  void triang_b(int i, int j) {
    node_tag it = 4 * (j * Ni_ + i) + 0;
    node_tag it0 = (j == 0) ? SOLID_WALL :
      4 * ((j - 1) * Ni_ + i) + 2;
    node_tag it1 = 4 * (j      * Ni_ + i) + 1;
    node_tag it2 = 4 * (j      * Ni_ + i) + 3;
    triang_triangs_.row(it) << it0, it1, it2;
    node_tag ip0 = j * (Ni_ + 1) + i;
    node_tag ip1 = j * (Ni_ + 1) + i + 1;
    node_tag ip2 = num_verts_ + j * Ni_ + i;
    triang_points_.row(it) << ip0, ip1, ip2;
    if (j == 0) {
      edge_points_.row(curr_e_) << ip0, ip1;
      edge_triangs_.row(curr_e_++) << it , it0;
    }
    edge_points_.row(curr_e_) << ip1, ip2 ;
    edge_triangs_.row(curr_e_++) << it , it1 ;
    edge_points_.row(curr_e_) << ip0, ip2 ;
    edge_triangs_.row(curr_e_++) << it , it2 ;
    if (j == 0) {
      // we've already put all 3 edges in mesh
      triang_edges_.row(it) << curr_e_ - 3, curr_e_ - 2, curr_e_ - 1;
      return;
    }
    if (j == 1) {
      // take bottom edge from previous step
      // (it is actually top edge of some top
      // triangle in row j - 1)
      triang_edges_.row(it) << curr_e_ - 7 * Ni_ + i + 2, curr_e_ - 2, curr_e_ - 1;
      return;
    }
    // in general we do same thing as
    // in case (j == 1), but as we
    // have put more bottom edges in mesh
    // then in case (j == 1), we should
    // now decrement curr_e_ in other way
    triang_edges_.row(it) << curr_e_ - 6 * Ni_ + 1, curr_e_ - 2, curr_e_ - 1;
  }
  void triang_r(int i, int j) {
    node_tag it = 4 * (j * Ni_ + i) + 1;
    node_tag it0 = (i == Ni_ - 1) ? SOLID_WALL :
      4 * (j * Ni_ + i + 1) + 3;
    node_tag it1 = 4 * (j * Ni_ + i) + 2;
    node_tag it2 = 4 * (j * Ni_ + i) + 0;
    triang_triangs_.row(it)  << it0, it1, it2 ;
    node_tag ip0 = j * (Ni_ + 1) + i + 1;
    node_tag ip1 = (j + 1) * (Ni_ + 1) + i + 1;
    node_tag ip2 = num_verts_ + j * Ni_ + i;
    triang_points_.row(it) << ip0, ip1, ip2 ;
    edge_points_ .row(curr_e_) << ip0, ip1 ;
    edge_triangs_.row(curr_e_++) << it , it0 ;
    edge_points_ .row(curr_e_) << ip1, ip2 ;
    edge_triangs_.row(curr_e_++) << it , it1 ;
    triang_edges_.row(it) << curr_e_ - 2, curr_e_ - 1, curr_e_ - 4;
  }
  void triang_t(int i, int j) {
    node_tag it = 4 * (j * Ni_ + i) + 2;
    node_tag it0 = (j == Nj_ - 1) ? SOLID_WALL :
      4 * ((j + 1) * Ni_ + i) + 0;
    node_tag it1 = 4 * (j * Ni_ + i) + 3;
    node_tag it2 = 4 * (j * Ni_ + i) + 1;
    triang_triangs_.row(it)  << it0, it1, it2 ;
    node_tag ip0 = (j + 1) * (Ni_ + 1) + i + 1;
    node_tag ip1 = (j + 1) * (Ni_ + 1) + i;
    node_tag ip2 = num_verts_ + j * Ni_ + i;
    triang_points_.row(it) << ip0, ip1, ip2 ;
    edge_points_.row(curr_e_) << ip1, ip0 ;
    edge_triangs_.row(curr_e_++) << it , it0 ;
    edge_points_.row(curr_e_) << ip1, ip2 ;
    edge_triangs_.row(curr_e_++) << it , it1 ;
    triang_edges_.row(it) << curr_e_ - 2, curr_e_ - 1, curr_e_ - 3;
  }
  void triang_l(int i, int j) {
    node_tag it = 4 * (j * Ni_ + i) + 3;
    node_tag it0 = (i == 0) ? SOLID_WALL :
      4 * (j * Ni_ + i - 1) + 1;
    node_tag it1 = 4 * (j * Ni_ + i) + 0;
    node_tag it2 = 4 * (j * Ni_ + i) + 2;
    triang_triangs_.row(it) << it0, it1, it2;
    node_tag ip0 = (j + 1) * (Ni_ + 1) + i;
    node_tag ip1 = j * (Ni_ + 1) + i;
    node_tag ip2 = num_verts_ + j * Ni_ + i;
    triang_points_.row(it) << ip0, ip1, ip2;
    if (i == 0) {
      // we add left edge of this square-box
      edge_points_.row(curr_e_) << ip1, ip0;
      edge_triangs_.row(curr_e_++) << it , it0;
      // so now our left edge is exactly
      // curr_e_-1
      triang_edges_.row(it) << curr_e_ - 1 , curr_e_ - 6, curr_e_ - 2;
      return;
    }
    if (i == 1) {
      // but in this case we already have edge
      // with this coordinates
      // it's right edge of square {0, j}
      // so we should decrement curr_e_ in other way
      if (j == 0)
        triang_edges_.row(it) << curr_e_ - 12, curr_e_ - 5, curr_e_ - 1;
      else
        triang_edges_.row(it) << curr_e_ - 11, curr_e_ - 5, curr_e_ - 1;
      return;
    }
    // in general we do same thing as
    // in case (i == 1), but as the 
    // previous square had number {i-1, j}
    // and i-1 != 0, we now should decrement
    // curr_e_ on smaller value then in case
    // (i == 1)
    if (j == 0)
      triang_edges_.row(it) << curr_e_ - 11, curr_e_ - 5, curr_e_ - 1;
    else
      triang_edges_.row(it) << curr_e_ - 10, curr_e_ - 5, curr_e_ - 1;
  }
};