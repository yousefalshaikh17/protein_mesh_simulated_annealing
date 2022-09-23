def write_ffea_output(output_file, nvoxel, tet_array, points, faces, original_ids, comment=""):
    """
    write the ffea input files
    Args:
        nvoxel (int): number of voxels
        tet_array ():
        points ():
        faces
        original_ids
        comment (str): any user comment
    """
    write_ffea_elements(output_file, nvoxel, tet_array, comment)
    write_ffea_nodes(output_file, points, comment)
    write_ffea_faces(output_file, faces, original_ids, comment)

def write_ffea_elements(output_file, nvoxel, tet_array, comment=""):
    """
    write ffea elements file .ele
    First line: <# of tetrahedra> <nodes per tetrahedron> <# of attributes>
    Remaining lines list of # of tetrahedra:
    <tetrahedron #> <node> <node> <node> <node> ... [attributes]
    Args:
        output_file (pathlib.Path): root name of file, will have .node added
        nvoxel (int): number of voxels
        tet_array ():
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.ele"), "w") as ele:
        ele.write(f"{nvoxel*5} 4 0\n")
        for i, tet in enumerate(tet_array):
            ele.write(f"{i+1} {tet[0]+1} {tet[1]+1} {tet[2]+1} {tet[3]+1}\n")
        ele.write(comment)

def write_ffea_nodes(output_file, points, comment=""):
    """
    write ffea .node file
    First line: <# of points> <dimension (must be 3)> <# of attributes> <# of boundary markers (0 or 1)>
    Remaining lines list # of points:
    <point #> <x> <y> <z>
    Args:
        output_file (pathlib.Path): root name of file, will have .node added
        points ():
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.node"), "w") as node:
        node.write(f"{len(points)} 3 0 0\n")
        for i, point in enumerate(points):
            node.write(f"{i+1} {point[0]} {point[1]} {point[2]}\n")
        node.write(comment)

def write_ffea_faces(output_file, faces, original_ids, comment=""):
    """
    write .face file
    First line: <# of faces> <boundary marker (0 or 1)>
    Remaining lines list of # of faces:
    <face #> <node> <node> <node> [boundary marker]
    Args:
        output_file (pathlib.Path): root name of file, will have .face added
        faces
        original_ids
        comment (str): any user comment
    """
    with open(output_file.with_suffix(".1.face"), "w") as face_file:
        face_file.write(f"{len(faces)} 1\n")
        for i, face in enumerate(faces):
            face_file.write(
                f"{i+1} {original_ids[face[1]]+1} {original_ids[face[2]]+1} {original_ids[face[3]]+1} -1\n")
        face_file.write(comment)