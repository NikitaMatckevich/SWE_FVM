import gmsh

def parse_msh_file(filename):
    nodes = {}
    triangles = {}
    edges = {}
    boundary_edges = {}
    
    gmsh.initialize()
    gmsh.open(filename)
    
    # Get nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    for tag, coord in zip(node_tags, node_coords):
        nodes[tag] = coord
    
    # Get elements (triangles and edges)
    element_tags, element_types, element_nodes = gmsh.model.mesh.getElements()
    for tag, type_code, node_ids in zip(element_tags, element_types, element_nodes):
        if type_code == 2:  # Triangle
            triangles[tag] = node_ids
            for i in range(3):
                edge = tuple(sorted((node_ids[i], node_ids[(i+1)%3])))
                if edge not in edges:
                    edges[edge] = []
                edges[edge].append(tag)
        elif type_code == 1:  # Edge
            edge = tuple(node_ids)
            edges[edge] = tag

            # Check if edge is on boundary
            boundary_tag = gmsh.model.mesh.getPhysicalGroupsForEntity(1, tag)
            if boundary_tag:
                if boundary_tag[0] not in boundary_edges:
                    boundary_edges[boundary_tag[0]] = []
                boundary_edges[boundary_tag[0]].append(tag)

    gmsh.finalize()
    
    return nodes, triangles, edges, boundary_edges

def write_topology_file(triangles, edges, boundary_edges, filename):
    with open(filename, 'w') as f:
        # Write edges
        for edge, triangles_adjacent in edges.items():
            node1, node2 = edge
            if isinstance(triangles_adjacent, list):
                for triangle_id in triangles_adjacent:
                    f.write(f"{node1} {node2} {triangle_id} {-1}\n")  # -1 indicates boundary edge
            else:
                physical_group_id = boundary_edges[triangles_adjacent]
                f.write(f"{node1} {node2} {triangles_adjacent} {-physical_group_id}\n")

        # Write triangles
        for triangle_id, node_ids in triangles.items():
            edge_ids = [edges[tuple(sorted((node_ids[i], node_ids[(i+1)%3])))] for i in range(3)]
            f.write(f"{node_ids[0]} {node_ids[1]} {node_ids[2]} {edge_ids[0]} {edge_ids[1]} {edge_ids[2]} ")
            for i in range(3):
                adjacent_triangle = [tri for tri in edges[(node_ids[i], node_ids[(i+1)%3])]]
                if len(adjacent_triangle) == 1:
                    f.write(f"{adjacent_triangle[0]} ")
                else:
                    f.write(f"-1 ")
            f.write("\n")

def write_geometry_file(nodes, filename):
    with open(filename, 'w') as f:
        for node_id, coordinates in nodes.items():
            f.write(f"{node_id} {coordinates[0]} {coordinates[1]}\n")

# Example usage
nodes, triangles, edges, boundary_edges = parse_msh_file("basic.msh")
write_topology_file(triangles, edges, boundary_edges, "topology.txt")
write_geometry_file(nodes, "geometry.txt")
