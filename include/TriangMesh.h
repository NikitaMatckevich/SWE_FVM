#pragma once
#include <Exceptions.h>
#include <PointOperations.h>

using NodeTag = Idx; 
using EdgeTag   = Eigen::Array<Idx, 1, 2, Eigen::RowMajor>;
using TriangTag = Eigen::Array<Idx, 1, 3, Eigen::RowMajor>;

using NodeTagArray   = Eigen::Array<Idx, 1, Eigen::Dynamic, Eigen::RowMajor>;
using EdgeTagArray   = Eigen::Array<Idx, Eigen::Dynamic, 2, Eigen::RowMajor>;
using TriangTagArray = Eigen::Array<Idx, Eigen::Dynamic, 3, Eigen::RowMajor>;

struct Topology {
 
    static Topology create(
            NodeTag numNodes,
            const EdgeTagArray& edgeNodes,
            const EdgeTagArray& edgeElements,
            const TriangTagArray& elementNodes,
            const TriangTagArray& elementEdges,
            const TriangTagArray& elementNeighbours)
    {
        return Topology(numNodes, edgeNodes, edgeElements, elementNodes, elementEdges, elementNeighbours);
    }

    Topology(
        NodeTag numNodes,
        const EdgeTagArray& edgeNodes,
        const EdgeTagArray& edgeElements,
        const TriangTagArray& elementNodes,
        const TriangTagArray& elementEdges,
        const TriangTagArray& elementNeighbours)
    : m_numNodes(numNodes)
    , m_edge_points(edgeNodes)
    , m_edge_triangs(edgeElements)
    , m_triang_points(elementNodes)
    , m_triang_edges(elementEdges)
    , m_triang_triangs(elementNeighbours) {}

    inline NodeTag NumNodes() const { return m_numNodes; }
    inline NodeTag NumEdges() const { return m_edge_points->col(0).size(); }
    inline NodeTag NumTriangles() const { return m_triang_points->col(0).size(); }
    
    EdgeTag EdgePoints(NodeTag i) const;
    EdgeTag EdgeTriangs(NodeTag i) const;
    bool IsEdgeBoundary(NodeTag i) const;
    
    TriangTag TriangPoints(NodeTag i) const;
    TriangTag TriangEdges(NodeTag i) const;
    TriangTag TriangTriangs(NodeTag i) const;
    bool IsTriangleBoundary(NodeTag i) const;
    
private:
    const NodeTag  m_numNodes; // number of mesh nodes
    const EdgeTagArray* m_edge_points;  // ends of each edge
    const EdgeTagArray* m_edge_triangs; // triangle pairs that have 1 edge in common
                                        // some special types are used for boundary
                                        // edge's "ghost cells"
  
    const TriangTagArray* m_triang_points;  // points of each triangle
    const TriangTagArray* m_triang_edges;   // edges of each triangle
    const TriangTagArray* m_triang_triangs; // triangles of each triangle
};

