import bpy
import bmesh
import re
import numpy as np
from . import _math, read_Files
from .Chem_Data import ELEMENTS_DEFAULT, BONDS_DEFAULT

'----------------------------------------------------------------------------------'
language = 1 if "zh_HAN" in bpy.context.preferences.view.language else 0
version = bpy.app.version
if language == 0:
    bsdf = "Principled BSDF"
else:
    if version[:2] == (4,0) or version[:2] == (4,1):
        bsdf = "原理化BSDF"
    else:
        bsdf = "原理化 BSDF"
'##################################################################################'
'##########################   Objects and Materials  ##############################'
# This function creates mesh objects from verts, edges and faces.
def create_object(name, coll, verts, edges, faces):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, edges, faces)
    object = bpy.data.objects.new(name, mesh)
    coll.objects.link(object)
    return object

# This function creates objects from mesh data.
def create_obj_from_mesh(coll, mesh):
    obj = bpy.data.objects.new(mesh.name, mesh)
    coll.objects.link(obj)
    return obj

def get_object(name):
    target_obj = None
    for obj in bpy.data.objects:
        if name in obj.name:
            target_obj = obj
            break
    return target_obj

def get_modifiers(obj, name):
    target_modifier = []
    for modifier in obj.modifiers:
        if name in modifier.name:
            target_modifier = modifier
            break
    return target_modifier

def join_objects(obj_list):
    object = obj_list[0]
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    for obj in obj_list: obj.select_set(True)
    bpy.ops.object.join()
    return object

def remove_doubles(object):
    bpy.context.view_layer.objects.active = object
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.remove_doubles(threshold=0.001)
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

def dissolve_flat_edges(object):
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(object.data)
    bpy.ops.mesh.select_all(action='DESELECT')
    for edge in bm.edges:
        if len(edge.link_faces) == 2:
            if edge.calc_face_angle() <= 0.01:
                edge.select_set(True)
    bpy.ops.mesh.dissolve_edges()
    bmesh.update_edit_mesh(object.data)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')

# This function creates materials upon exist elements(ATOMS) and materials.
# Elements here is a list of ATOMS
def create_atom_materials(Elements):
    for material in bpy.data.materials:
        if material.name in Elements:
            Elements.remove(material.name)
    mats = []
    for key in list(ELEMENTS_DEFAULT.keys()):
        if key in Elements:
            mat = bpy.data.materials.new(name=key)
            mat.use_nodes = True
            _BSDF = mat.node_tree.nodes[bsdf]
            _BSDF.inputs['Base Color'].default_value = ELEMENTS_DEFAULT[key][3]  # Color
            _BSDF.inputs['Roughness'].default_value = 0.2   # Roughness
            mats.append(mat)
    return mats




'----------------------------------------------------------------------------------'
'##################################################################################'
'###############################   Attributes  ####################################'
# Add an attribute to the given object on the given domain.
types = {
    'FLOAT': 'value',
    'INT': 'value',
    'BOOLEAN': 'value',
    'QUATERNION': 'value',
    'FLOAT_COLOR': 'color',
    'FLOAT_VECTOR': 'vector',
}
def get_attirbute(obj, attr_name, datatype, domain):
    bpy.ops.object.mode_set(mode='OBJECT')
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    attr = obj.data.attributes.get(attr_name)
    count = len(bm.verts) if domain == 'VERT' else len(bm.edges)
    if datatype == 'FLOAT_COLOR':
        m = 4
    elif datatype == 'FLOAT_VECTOR':
        m = 3
    else:
        m = 1
    tot_count = count*m
    seq = [0]*tot_count
    attr.data.foreach_get(types[datatype], seq)
    attr.data.update()
    
    new_seq = []
    if datatype == 'FLOAT_COLOR':
        for i in range(count):
            new_seq.append(seq[m*i])
            new_seq.append(seq[m*i+1])
            new_seq.append(seq[m*i+2])
            new_seq.append(seq[m*i+3])
    elif datatype == 'FLOAT_VECTOR':
        for i in range(count):
            new_seq.append(seq[m*i])
            new_seq.append(seq[m*i+1])
            new_seq.append(seq[m*i+2])
    else:
        new_seq = seq
    bm.free()
    return new_seq

def add_attribute(obj, attr_name, datatype, domain, values):
    attr = obj.data.attributes.get(attr_name)
    if not attr:
        attr = obj.data.attributes.new(attr_name, datatype, domain)
    attr.data.foreach_set(types[datatype], values)
    return attr

# add atoms and bonds radii attributes
def add_base_radii_attr(obj, Atoms, Bonds, Bond_Orders):
    
    atom_ids = [ELEMENTS_DEFAULT[element][0] for element in Atoms]
    atom_radii_B = [ELEMENTS_DEFAULT[element][4]/2 for element in Atoms]
    atom_radii_S = [ELEMENTS_DEFAULT['Bond'][4]*1.8 for element in Atoms]
    atom_radii_W = [0]*len(Atoms)
    atom_radii_F = [ELEMENTS_DEFAULT[element][7] if ELEMENTS_DEFAULT[element][7] else 1.5 for element in Atoms]
    atom_colors = []
    for element in Atoms:
        for value in ELEMENTS_DEFAULT[element][3]:
            atom_colors.append(value)
    bond_radii_B = [ELEMENTS_DEFAULT['Bond'][4]]*len(Bond_Orders)
    for i, order in enumerate(Bond_Orders): bond_radii_B[i] = ELEMENTS_DEFAULT['Bond'][4] if order==1 else 0.8*ELEMENTS_DEFAULT['Bond'][4]
    bond_radii_S = [ELEMENTS_DEFAULT['Bond'][4]*1.8]*len(Bond_Orders)
    bond_radii_W = [ELEMENTS_DEFAULT['Bond'][5]]*len(Bond_Orders)
    bond_aids_1 = [ELEMENTS_DEFAULT[Atoms[bond[0]]][0] for bond in Bonds]
    bond_aids_2 = [ELEMENTS_DEFAULT[Atoms[bond[1]]][0] for bond in Bonds]

    attributes = (
        {'name': 'bond_order',    'type':'INT',      'domain':'EDGE',   'values':Bond_Orders},
        {'name': 'bond_aids_1',   'type':'INT',      'domain':'EDGE',   'values':bond_aids_1},
        {'name': 'bond_aids_2',   'type':'INT',      'domain':'EDGE',   'values':bond_aids_2},
        {'name': 'bond_radii_B',  'type':'FLOAT',    'domain':'EDGE',   'values':bond_radii_B},
        {'name': 'bond_radii_S',  'type':'FLOAT',    'domain':'EDGE',   'values':bond_radii_S},
        {'name': 'bond_radii_W',  'type':'FLOAT',    'domain':'EDGE',   'values':bond_radii_W},
        {'name': 'atom_id',       'type':'INT',      'domain':'POINT',  'values':atom_ids},
        {'name': 'radii_B',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_B},
        {'name': 'radii_S',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_S},
        {'name': 'radii_W',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_W},
        {'name': 'radii_F',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_F},
        {'name': 'colour',   'type':'FLOAT_COLOR',    'domain':'POINT',  'values':atom_colors},
    )
    for attr in attributes:
        add_attribute(obj, attr['name'], attr['type'], attr['domain'], attr['values'])

def get_bond_key(atom_A, atom_B):
    if ELEMENTS_DEFAULT[atom_A][0] <= ELEMENTS_DEFAULT[atom_B][0]:
        bond_key = f'{atom_A},{atom_B}'
    else:
        bond_key = f'{atom_B},{atom_A}'
    if bond_key not in BONDS_DEFAULT: bond_key = 'Default'
    return bond_key

# set atoms
def set_sel_atoms_attr(obj, element):
    sel_vert_idxs = get_sel_idxs(obj, 'VERT')
    atom_ids = get_attirbute(obj, 'atom_id', 'INT', 'VERT')
    radii_B = get_attirbute(obj, 'radii_B', 'FLOAT', 'VERT')
    colors = get_attirbute(obj, 'colour', 'FLOAT_COLOR', 'VERT')
    for idx in sel_vert_idxs:
        atom_ids[idx] = ELEMENTS_DEFAULT[element][0]
        radii_B[idx] = ELEMENTS_DEFAULT[element][4]/2
        colors[4*idx],colors[4*idx+1],colors[4*idx+2],colors[4*idx+3] = ELEMENTS_DEFAULT[element][3]
    add_attribute(obj, 'atom_id', 'INT', 'VERT', atom_ids)
    add_attribute(obj, 'radii_B', 'FLOAT', 'VERT', radii_B)
    add_attribute(obj, 'colour', 'FLOAT_COLOR', 'VERT', colors)

    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    ordi_num_layer = bm.verts.layers.int.get('atom_id')
    for vert in bm.verts:
        if vert.select and len(vert.link_edges) == 1:
            edge = vert.link_edges[0]
            if edge.verts[0].index != vert.index:
                link_vert = edge.verts[0]
            else:
                link_vert = edge.verts[1]
            atom_id1 = vert[ordi_num_layer]
            atom_id2 = link_vert[ordi_num_layer]
            atom1 = list(ELEMENTS_DEFAULT.keys())[atom_id1]
            atom2 = list(ELEMENTS_DEFAULT.keys())[atom_id2]
            key = get_bond_key(atom1, atom2)
            dir = _math.normalize(np.array(vert.co)-np.array(link_vert.co))
            vert.co = np.array(link_vert.co) + dir * (BONDS_DEFAULT[key][2]+BONDS_DEFAULT[key][3])/2
    bmesh.update_edit_mesh(obj.data)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')

# set bonds
def set_sel_bonds_attr(obj, bond_order):
    sel_edge_idxs = get_sel_idxs(obj, 'EDGE')
    Bond_Orders = get_attirbute(obj, 'bond_order', 'INT', 'EDGE')
    if bond_order >= 0:
        for idx in sel_edge_idxs: Bond_Orders[idx] = bond_order
    bond_radii = [1]*len(Bond_Orders)
    bond_radii_B = [ELEMENTS_DEFAULT['Bond'][4]]*len(Bond_Orders)
    for i, order in enumerate(Bond_Orders): bond_radii_B[i] = ELEMENTS_DEFAULT['Bond'][4] if order==1 else 0.8*ELEMENTS_DEFAULT['Bond'][4]
    bond_radii_S = [ELEMENTS_DEFAULT['Bond'][4]*1.8]*len(Bond_Orders)
    bond_radii_W = [ELEMENTS_DEFAULT['Bond'][5]]*len(Bond_Orders)
    add_attribute(obj, 'bond_order', 'INT', 'EDGE', Bond_Orders)
    add_attribute(obj, 'bond_radii_B', 'FLOAT', 'EDGE', bond_radii_B)
    add_attribute(obj, 'bond_radii_S', 'FLOAT', 'EDGE', bond_radii_S)
    add_attribute(obj, 'bond_radii_W', 'FLOAT', 'EDGE', bond_radii_W)

def set_sel_bonds_aromatic(obj, is_aromatic):
    sel_edge_idxs = get_sel_idxs(obj, 'EDGE')
    sel_vert_idxs = get_sel_idxs(obj, 'VERT')
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    coord = np.array((0.0,0.0,0.0))
    for idx in sel_vert_idxs:
        coord += bm.verts[idx].co
    avg_coord = coord/len(sel_vert_idxs)
    bond_vertvec_layer = bm.edges.layers.float_vector.get('bond_vertvec')
    for idx in sel_edge_idxs:
        edge = bm.edges[idx]
        vec1 = edge[bond_vertvec_layer]
        vec2 = avg_coord - np.array(edge.verts[0].co)
        if np.dot(vec1,vec2)<0:
            edge[bond_vertvec_layer] *= -1
    bmesh.update_edit_mesh(obj.data)
    bm.free()
    Aromatic_Edges = get_attirbute(obj, 'is_aromatic_edge', 'INT', 'EDGE')
    for idx in sel_edge_idxs: Aromatic_Edges[idx] = is_aromatic
    add_attribute(obj, 'is_aromatic_edge', 'INT', 'EDGE', Aromatic_Edges)
    

def get_link_vert(vert, edge):
    if edge.verts[0].index == vert.index:
        return edge.verts[1]
    else:
        return edge.verts[0]

def get_layer(bm, datatype, domain, name):
    if datatype == 'INT':
        if domain == 'VERT':
            attr_layer = bm.verts.layers.int.get(name)
            if not attr_layer:
                attr_layer = bm.verts.layers.int.new(name)
        elif domain == 'EDGE':
            attr_layer = bm.edges.layers.int.get(name)
            if not attr_layer:
                attr_layer = bm.edges.layers.int.new(name)
    if datatype == 'FLOAT':
        if domain == 'VERT':
            attr_layer = bm.verts.layers.float.get(name)
            if not attr_layer:
                attr_layer = bm.verts.layers.float.new(name)
        elif domain == 'EDGE':
            attr_layer = bm.edges.layers.float.get(name)
            if not attr_layer:
                attr_layer = bm.edges.layers.float.new(name)
    if datatype == 'FLOAT_VECTOR':
        if domain == 'VERT':
            attr_layer = bm.verts.layers.float_vector.get(name)
            if not attr_layer:
                attr_layer = bm.verts.layers.float_vector.new(name)
        elif domain == 'EDGE':
            attr_layer = bm.edges.layers.float_vector.get(name)
            if not attr_layer:
                attr_layer = bm.edges.layers.float_vector.new(name)
    return attr_layer

def is_aromatic(scaffold):
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bm = bmesh.from_edit_mesh(scaffold.data)
    bm.verts.ensure_lookup_table()
    edge_counts = len(bm.edges)
    for vert in bm.verts:
        if not vert.link_edges:
            vert.select_set(False)
    bpy.ops.mesh.edge_face_add()
    bond_order_layer = bm.edges.layers.int.get('bond_order')
    bond_vertvec_layer = bm.edges.layers.float_vector.get('bond_vertvec')
    atom_id_layer = bm.verts.layers.int.get('atom_id')
    is_aromatic_edge = get_layer(bm, 'INT','EDGE','is_aromatic_edge')
    ring_id_layer = get_layer(bm, 'INT','EDGE','aromatic_ring_id')
    aromatic_ring_id = 0
    bm.faces.ensure_lookup_table()
    for face in bm.faces:
        avg_coord = avg_coords(face.verts)
        is_aromatic_ring = False
        # is_aromatic judgement
        hexaring_FLAG = True
        
        if len(face.verts) == 5:
            orders_1, orders_2 = 0,0
            for edge in face.edges:
                if edge[bond_order_layer] == 1: orders_1 += 1
                if edge[bond_order_layer] == 2: orders_2 += 1
            if orders_1 == 3 and orders_2 == 2:
                is_aromatic_ring = True
            
        if len(face.verts) == 6:
            orders_1, orders_2, aromatic_edge = 0,0,0
            for edge in face.edges:
                if edge[bond_order_layer] == 1: orders_1 += 1
                if edge[bond_order_layer] == 2: orders_2 += 1
                if edge[is_aromatic_edge] == 1: aromatic_edge += 1
            if orders_1 == 3 and orders_2 == 3:
                is_aromatic_ring = True
            elif orders_1 > 3 and orders_2+aromatic_edge >= 3:
                is_aromatic_ring = True
            for vert in face.verts:
                if vert[atom_id_layer] not in [6,7]:
                    hexaring_FLAG = False
                    break

        if len(face.verts):
            aromatic_edge = 0
            for edge in face.edges:
                if edge[is_aromatic_edge] == 1: aromatic_edge += 1
            if aromatic_edge == len(face.edges):
                is_aromatic_ring = True

        if is_aromatic_ring and hexaring_FLAG:
            aromatic_ring_id += 1
            for edge in face.edges: edge[is_aromatic_edge] = 1
            # for vert in face.verts: vert[is_aromatic_vert] = 1
            for edge in face.edges: edge[ring_id_layer] = aromatic_ring_id
            for edge in face.edges:
                vec1 = edge[bond_vertvec_layer]
                vec2 = avg_coord - np.array(edge.verts[0].co)
                if np.dot(vec1,vec2)<0:
                    edge[bond_vertvec_layer] *= -1
    
    bpy.ops.mesh.delete(type='ONLY_FACE')
    bm.edges.ensure_lookup_table()
    for edge in bm.edges:
        if edge.index >= edge_counts:
            edge.select_set(True)
            bpy.ops.mesh.delete(type='EDGE')
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')


'----------------------------------------------------------------------------------'
'##################################################################################'
'################################   Selection  ####################################'
def select_all(obj, domain):
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type=domain)
    bpy.ops.mesh.select_all(action='SELECT')

def deselect_all(obj):
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

def get_sel_idxs(obj, domain):
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    sel_idxs = []
    if domain == 'VERT':
        sel_idxs = [vert.index for vert in bm.verts if vert.select]
    elif domain == 'EDGE':
        sel_idxs = [edge.index for edge in bm.edges if edge.select]
    bmesh.update_edit_mesh(obj.data)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')
    return sel_idxs

def select_verts(obj, elements):   # 'C', 'H'
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    atom_id = get_layer(bm, 'INT', 'VERT', 'atom_id')
    Elements = list(ELEMENTS_DEFAULT.keys())
    for vert in bm.verts:
        if Elements[vert[atom_id]] in elements:
            vert.select_set(True)

def select_edges(obj, bond_syms):   # ['C','H']
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    atom_id = get_layer(bm, 'INT', 'VERT', 'atom_id')
    Elements = list(ELEMENTS_DEFAULT.keys())
    for edge in bm.edges:
        atom1 = Elements[edge.verts[0][atom_id]]
        atom2 = Elements[edge.verts[1][atom_id]]
        if [atom1,atom2] in bond_syms or [atom2,atom1] in bond_syms:
            edge.select_set(True)

def select_bond_orders(obj, bond_orders):
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bond_order = get_layer(bm, 'INT', 'EDGE', 'bond_order')
    for edge in bm.edges:
        if re.sub('[\W]','',f'{edge[bond_order]}') in bond_orders:
            edge.select_set(True)
    




'----------------------------------------------------------------------------------'
'##################################################################################'
'################################   Calculation  ##################################'
# This function measures the length of an edge in geometry
def calc_length():
    # in Edit Mode (select an edge)
    if bpy.context.mode == 'EDIT_MESH':
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)
        locations = []
        for edge in bm.edges:
            if edge.select:
                locations.append(edge.verts[0].co)
                locations.append(edge.verts[1].co)
                break
        if len(locations) == 0: return "N.A."

        pos1 = np.array(locations[0])
        pos2 = np.array(locations[1])
        dist = _math.get_length(pos1-pos2)
        dist = str(float('%.3f' % dist)) + " Å"   # retain three digits behind decimal point
        return dist
    
def calc_angle():
    # in Edit Mode (selct two edges)
    if bpy.context.mode == 'EDIT_MESH':
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)
        bm.edges.ensure_lookup_table()
        locations = []
        for edge in bm.edges:
            if edge.select:
                locations.append(edge.verts[0].co)
                locations.append(edge.verts[1].co)
                if len(locations) == 4: break
        if len(locations) < 4: return "N.A."

        pos1 = np.array(locations[0])
        pos2 = np.array(locations[1])
        pos3 = np.array(locations[2])
        pos4 = np.array(locations[3])
        epsilon=0.000001
        if _math.get_length(pos2-pos3)<epsilon or _math.get_length(pos1-pos4)<epsilon:
            vec1 = _math.normalize(pos2-pos1)
            vec2 = _math.normalize(pos3-pos4)
        else:
            vec1 = _math.normalize(pos1-pos2)
            vec2 = _math.normalize(pos3-pos4)
        rad = np.arccos(np.dot(vec1,vec2))
        angle = np.rad2deg(rad)
        angle = str(float('%.2f' % angle)) + "°"   # retain two digits behind decimal point
        return angle


def avg_dummy_xyz(obj):
    # in EDIT mode
    bm = bmesh.from_edit_mesh(obj.data)
    sum_xyz = np.array((0.0, 0.0, 0.0))
    counts = 0
    for vert in bm.verts:
        if vert.select == True:
            sum_xyz += vert.co
            counts += 1
    avg_xyz = sum_xyz / counts
    return avg_xyz

def avg_coords(vertices):
    sum_xyz = np.array((0.0, 0.0, 0.0))
    for vert in vertices:
        sum_xyz += np.array(vert.co)
    avg_xyz = sum_xyz / len(vertices)
    return avg_xyz

# get heading, pitch, bank vectors of a vert
def vert_hpb(vert, rotate_angle):
    vecs = []
    for edge in vert.link_edges:
        link_vert = get_link_vert(vert, edge)
        vecs.append(_math.normalize(np.array(vert.co) - np.array(link_vert.co)))
    bank = np.array((0.0, 0.0, 0.0))
    for vec in vecs: bank += vec
    bank = _math.normalize(bank)
    # 此处有待优化
    # -------------
    if len(vecs)>=3:
        if abs(np.dot(np.cross(vecs[0],vecs[1]),vecs[2])) <= 0.0025:   # coplanar vecs
            bank = np.cross(vecs[0],vecs[1])
    # -------------
    if len(vecs) == 0:
        bank = np.array((1, 0, 0))
        pitch = np.array((0, 1, 0))
    elif len(vecs) == 1:
        pitch = _math.vertical_vec(np.array((0,0,0)), bank)
    elif len(vecs) >= 2:
        pitch = _math.normalize(vecs[1]-vecs[0])
    pitch = _math.rotate_vec(pitch, bank, rotate_angle)
    heading = np.cross(bank, pitch)
    return (bank, pitch, heading)

# calculate vertical vector for all bonds
def bond_vertical_dir(scaffold):
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    bond_vertic_vecs = []
    for edge in bm.edges:
        links_1 = edge.verts[0].link_edges
        links_2 = edge.verts[1].link_edges
        if len(links_1) > 1: 
            try:
                pos1,pos2,pos3 = list(set([tuple(np.array(links_1[0].verts[0].co)),
                                           tuple(np.array(links_1[0].verts[1].co)),
                                           tuple(np.array(links_1[1].verts[0].co)),
                                           tuple(np.array(links_1[1].verts[1].co)),]))
            except:
                pos1,pos2,pos3 = list(set([tuple(np.array(links_1[0].verts[0].co)),
                                           tuple(np.array(links_1[0].verts[1].co)),
                                           tuple(np.array(links_1[2].verts[0].co)),
                                           tuple(np.array(links_1[2].verts[1].co)),]))
        elif len(links_2) > 1:
            try:
                pos1,pos2,pos3 = list(set([tuple(np.array(links_2[0].verts[0].co)),
                                           tuple(np.array(links_2[0].verts[1].co)),
                                           tuple(np.array(links_2[1].verts[0].co)),
                                           tuple(np.array(links_2[1].verts[1].co)),]))
            except:
                pos1,pos2,pos3 = list(set([tuple(np.array(links_2[0].verts[0].co)),
                                           tuple(np.array(links_2[0].verts[1].co)),
                                           tuple(np.array(links_2[2].verts[0].co)),
                                           tuple(np.array(links_2[2].verts[1].co)),]))
        else:
            pos1,pos2,pos3 = (edge.verts[0].co, edge.verts[1].co, np.array(edge.verts[0].co)+np.array((1.,0.,0.)))
        normal = _math.get_normal(pos1, pos2, pos3)
        vertic_vec = _math.normalize(np.cross((np.array(edge.verts[0].co)-np.array(edge.verts[1].co)), normal))
        for component in vertic_vec: bond_vertic_vecs.append(component)
    return bond_vertic_vecs


'----------------------------------------------------------------------------------'
'##################################################################################'
'###################################   Others  ####################################'
# The mesh object should have attribute 'atom_id'
def atom_type_symbols(scaffold):
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    ordi_num_layer = bm.verts.layers.int.get('atom_id')
    atom_type_symbols = []
    for vert in bm.verts:
        ordi = vert[ordi_num_layer]
        sym = list(ELEMENTS_DEFAULT.keys())[ordi]
        atom_type_symbols.append(sym)
    bm.free()
    return atom_type_symbols


def atom_orient(scaffold):
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    atom_orientations = []
    for vert in bm.verts:
        atom_dir = np.array((0.0,0.0,0.0))
        for edge in vert.link_edges:
            link_vert = get_link_vert(vert, edge)
            # calculate the orientation of the axis of each ball.
            atom_dir += np.array(link_vert.co) - np.array(vert.co)
        for component in _math.normalize(atom_dir): atom_orientations.append(component)
    bm.free()
    return atom_orientations


def bond_type_symbols(scaffold):
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    ordi_num_layer = bm.verts.layers.int.get('atom_id')
    bond_type_symbols = []
    atom_orientations = []
    for vert in bm.verts:
        ordi_A = vert[ordi_num_layer]
        sym_A = list(ELEMENTS_DEFAULT.keys())[ordi_A]
        atom_dir = np.array((0.0,0.0,0.0))
        for edge in vert.link_edges:
            if edge.verts[0].index != vert.index:
                link_vert = edge.verts[0]
            else:
                link_vert = edge.verts[1]
            ordi_B = int(link_vert[ordi_num_layer])
            sym_B = list(ELEMENTS_DEFAULT.keys())[ordi_B]
            if (sym_A, sym_B) not in bond_type_symbols: bond_type_symbols.append((sym_A, sym_B))
            if (sym_B, sym_A) not in bond_type_symbols: bond_type_symbols.append((sym_B, sym_A))
            
            # calculate the orientation of the axis of each ball.
            atom_dir += np.array(link_vert.co) - np.array(vert.co)
        for component in _math.normalize(atom_dir): atom_orientations.append(component)
    bm.free()
    return (bond_type_symbols, atom_orientations)


def struct_to_block(scaffold):
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    atoms_num = len(bm.verts)
    bonds_num = len(bm.edges)
    ordi_num_layer = bm.verts.layers.int.get('atom_id')
    bond_order_layer = bm.edges.layers.int.get('bond_order')
    Atoms, Coords, Bonds, Bond_Orders = [],[],[],[]
    for vert in bm.verts:
        ordi_A = vert[ordi_num_layer]
        sym_A = list(ELEMENTS_DEFAULT.keys())[ordi_A]
        Atoms.append(sym_A)
        coord = ["{:.4f}".format(vert.co[0]),"{:.4f}".format(vert.co[1]),"{:.4f}".format(vert.co[2])]
        Coords.append(coord)
    for edge in bm.edges:
        bvert_id_1 = edge.verts[0].index+1
        bvert_id_2 = edge.verts[1].index+1
        Bonds.append((bvert_id_1, bvert_id_2))
        Bond_Orders.append(edge[bond_order_layer])
    Header_Block = '\n'+'    RDKIT    3D\n'+'\n'
    Counts_Line = f' {atoms_num} {bonds_num} 0 0 0 0 0 0 0 0999 V2000\n'
    Atom_Block, Bond_Block = '',''
    for atom, coord in zip(Atoms, Coords):
        for fract in coord:
            space = 9-len(fract)
            Atom_Block += ' '*(space+1)+f'{fract}'
        Atom_Block += f' {atom}   0 0 0 0 0 0 0 0 0 0 0 0\n'
    for bond, bond_order in zip(Bonds, Bond_Orders):
        for id in bond:
            if id > 9:
                Bond_Block += f' {id}'
            else:
                Bond_Block += f'  {id}'        
        Bond_Block += f'  {bond_order}  0\n'
    End_Block = 'M  END'
    Block = Header_Block + Counts_Line + Atom_Block + Bond_Block + End_Block
    return Block

def struct_optimize(scaffold, addHs):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    block = struct_to_block(scaffold)
    m = Chem.MolFromMolBlock(block)
    if addHs: m = Chem.AddHs(m)
    AllChem.MMFFOptimizeMolecule(m)
    new_block = Chem.MolToMolBlock(m)
    file = new_block.split('\n')

    molecules = read_Files.read_MOL(file, 'mol', 'small_mol', np.array((0.0,0.0,0.0)), 1)
    Atoms, Coords, Bonds, Bond_Orders = molecules[0]
    coll = scaffold.users_collection[0]
    Coords = [coord+np.array((0.001,0.001,0.001)) for coord in Coords]   # deviation from origin point
    new_struct = create_object('new_struct', coll, Coords, Bonds, [])
    new_struct['Type'] = 'small_mol'
    new_struct['Elements'] = list(set(Atoms))
    create_atom_materials(['H'])
    add_base_radii_attr(new_struct, Atoms, Bonds, Bond_Orders)
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.context.view_layer.objects.active = new_struct

    # remove previous scaffold mesh_data
    obj_name = scaffold.name
    mesh_name = scaffold.data.name
    bpy.data.meshes.remove(scaffold.data)
    new_struct.name = obj_name
    new_struct.data.name = mesh_name
    return new_struct