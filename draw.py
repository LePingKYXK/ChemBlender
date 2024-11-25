import bpy
import bmesh
import numpy as np
from math import sin,tan,pi
from . import _mesh, _math, _shade
from .Chem_Data import ELEMENTS_DEFAULT, BONDS_DEFAULT, FUNCTIONAL_GROUPS

def branches_dir(vert, branches, deflection_angle):
    bank, pitch, heading = _mesh.vert_hpb(vert, rotate_angle=0)
    banks = []
    pitches = []
    if branches == 1:
        bank = _math.rotate_vec(bank, heading, deflection_angle*pi/180)
        pitch = np.cross(heading, bank)
        banks.append(bank)
        pitches.append(pitch)
    elif branches == 2:
        bank_1 = _math.rotate_vec(bank, heading, deflection_angle*pi/180)
        bank_2 = _math.rotate_vec(bank_1, bank, pi)
        pitch_1 = np.cross(heading, bank_1)
        pitch_2 = np.cross(heading, bank_2)
        banks.append(bank_1)
        banks.append(bank_2)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
    elif branches == 3:
        bank_1 = _math.normalize(bank+2.824*pitch)
        bank_2 = _math.rotate_vec(bank_1, bank, pi*2/3)
        bank_3 = _math.rotate_vec(bank_1, bank, pi*4/3)
        pitch_1 = np.cross(heading, bank_1)
        pitch_2 = _math.rotate_vec(pitch_1, bank, pi*2/3)
        pitch_3 = _math.rotate_vec(pitch_1, bank, pi*4/3)
        banks.append(bank_1)
        banks.append(bank_2)
        banks.append(bank_3)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
        pitches.append(pitch_3)
    elif branches == 4:
        bank_1 = pitch
        bank_2 = _math.rotate_vec(bank_1, bank, pi*2/3)
        bank_3 = _math.rotate_vec(bank_1, bank, pi*4/3)
        pitch_1 = np.cross(bank, bank_1)
        pitch_2 = np.cross(bank, bank_2)
        pitch_3 = np.cross(bank, bank_3)
        banks.append(bank_1)
        banks.append(bank_2)
        banks.append(bank_3)
        banks.append(bank)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
        pitches.append(pitch_3)
        pitches.append(pitch)
    elif branches == 5:
        bank_1 = pitch
        bank_2 = _math.rotate_vec(bank_1, bank, pi/2)
        bank_3 = _math.rotate_vec(bank_1, bank, pi)
        bank_4 = _math.rotate_vec(bank_1, bank, pi*3/2)
        pitch_1 = np.cross(bank, bank_1)
        pitch_2 = np.cross(bank, bank_2)
        pitch_3 = np.cross(bank, bank_3)
        pitch_4 = np.cross(bank, bank_4)
        banks.append(bank_1)
        banks.append(bank_2)
        banks.append(bank_3)
        banks.append(bank_4)
        banks.append(bank)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
        pitches.append(pitch_3)
        pitches.append(pitch_4)
        pitches.append(pitch)
    return (banks, pitches)


def para_branches(mol_scaffold, new_atom, bond_order, branches, calc_SA, default_angle, rotate_angle, deflection_angle):
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(mol_scaffold.data)
    ordi_num_layer = bm.verts.layers.int.get('atom_id')
    bond_order_layer = bm.edges.layers.int.get('bond_order')
    origin_atoms, origin_coords, new_atoms, new_coords, bond_orders = [],[],[],[],[]
    for vert in bm.verts:
        if vert.select == True:
            atom_id = vert[ordi_num_layer]
            origin_atom = list(ELEMENTS_DEFAULT.keys())[atom_id]
            bond_key = _mesh.get_bond_key(origin_atom, new_atom)
            bond_length = (BONDS_DEFAULT[bond_key][2]+BONDS_DEFAULT[bond_key][3])/2
            bank, pitch, heading = _mesh.vert_hpb(vert, rotate_angle=0)
            saturation = 0
            for edge in vert.link_edges: saturation += edge[bond_order_layer]
            if calc_SA: branches = abs(ELEMENTS_DEFAULT[origin_atom][8])-saturation
            if default_angle:
                rotate_angle = 90
                deflection_angle = 54.735 if branches == 2 else 0
                if atom_id in [8,16] and branches == 1:   # O,S atom
                    rotate_angle = 0
                    deflection_angle = 60
            banks = branches_dir(vert, branches, deflection_angle)[0]
            for vec in banks:
                vec = _math.rotate_vec(vec, bank, rotate_angle*pi/180)
                if atom_id == 7: # N atom
                    axis = pitch if branches == 1 else heading
                    vec = _math.rotate_vec(vec, axis, 15*pi/180)
                new_coord = np.array(vert.co)+vec*bond_length
                origin_atoms.append(origin_atom)
                origin_coords.append(vert.co)
                new_atoms.append(new_atom)
                bond_orders.append(bond_order)
                new_coords.append(new_coord)
    bpy.ops.mesh.select_all(action='DESELECT')
    return (origin_atoms, origin_coords, new_atoms, new_coords, bond_orders)


def add_atoms(mol_scaffold, coll, origin_atoms, origin_coords, new_atoms, new_coords, bond_orders):
    verts, edges, faces= [],[],[]
    Atoms = []
    count=0
    for origin_co, new_co, origin_atom, new_atom in zip(origin_coords, new_coords,origin_atoms, new_atoms):
        verts.append(origin_co)
        verts.append(new_co)
        Atoms.append(origin_atom)
        Atoms.append(new_atom)
        edges.append((count,count+1))
        count += 2
    if origin_atoms:
        new_struct = _mesh.create_object('new_struct', coll, verts, edges, faces)
        _mesh.add_base_radii_attr(new_struct, Atoms, edges, bond_orders)
        new_scaffold = _mesh.join_objects([mol_scaffold, new_struct])
        _mesh.remove_doubles(new_scaffold)
        atom_orientations = _mesh.atom_orient(new_scaffold)
        bond_vertic_vecs = _mesh.bond_vertical_dir(new_scaffold)
        _mesh.add_attribute(new_scaffold, 'atom_orient', 'FLOAT_VECTOR', 'POINT', atom_orientations)
        _mesh.add_attribute(new_scaffold, 'bond_vertvec', 'FLOAT_VECTOR', 'EDGE', bond_vertic_vecs)
        _mesh.is_aromatic(new_scaffold)


def add_functional_groups(mol_scaffold, coords, fg_name, n_carbon):
    coll = bpy.context.collection
    if fg_name == 'n-Chain':
        Atoms = ['C']*n_carbon
        Bonds = [(i,i+1) for i in range(n_carbon-1)]
        Bond_Orders = [1]*(n_carbon-1)
    else:
        Atoms = FUNCTIONAL_GROUPS[fg_name]['Atoms']
        Bonds = FUNCTIONAL_GROUPS[fg_name]['Bonds']
        Bond_Orders = FUNCTIONAL_GROUPS[fg_name]['Bond_Orders']
    new_struct = _mesh.create_object('new_molecule', coll, coords, Bonds, [])
    _mesh.add_base_radii_attr(new_struct, Atoms, Bonds, Bond_Orders)
    _mesh.create_atom_materials(set(Atoms))
    _shade.create_mat_hetero('Hetero')
    new_scaffold = _mesh.join_objects([mol_scaffold, new_struct])
    _mesh.remove_doubles(new_scaffold)
    atom_orientations = _mesh.atom_orient(mol_scaffold)
    bond_vertic_vecs = _mesh.bond_vertical_dir(mol_scaffold)
    _mesh.add_attribute(mol_scaffold, 'atom_orient', 'FLOAT_VECTOR', 'POINT', atom_orientations)
    _mesh.add_attribute(mol_scaffold, 'bond_vertvec', 'FLOAT_VECTOR', 'EDGE', bond_vertic_vecs)
    _mesh.is_aromatic(mol_scaffold)

def draw_cell_edges(filename, coll_mol, length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma, space_group):
    fract_corners_xyz = [[0,0,0], [1,0,0], [1,1,0], [0,1,0], [0,0,1], [1,0,1], [1,1,1], [0,1,1]]
    cartn_corners = [_math.fract_to_cartn(corner, length_a, length_b, length_c, 
                                            angle_alpha, angle_beta, angle_gamma) for corner in fract_corners_xyz]
    cell_edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    mesh = bpy.data.meshes.new('Cell_Edges_'+filename)
    new_obj = bpy.data.objects.new(mesh.name, mesh)
    coll_mol.objects.link(new_obj)
    mesh.from_pydata(cartn_corners,cell_edges,[])
    bpy.context.view_layer.objects.active = new_obj
    new_obj.select_set(True)

    # store basic infomation of crystal cell in customic properties of Cell_Edges object
    new_obj['Type'] = 'unit_cell'
    new_obj['cell lengths'] = f'{length_a},{length_b},{length_c}'
    new_obj['cell angles'] = f'{angle_alpha},{angle_beta},{angle_gamma}'
    new_obj['space group'] = space_group
    bpy.ops.object.convert(target='CURVE')
    bpy.context.object.data.bevel_depth = 0.01
    bpy.ops.object.select_all(action='DESELECT')
    return new_obj


def add_CC():
    cursor = bpy.context.scene.cursor.location
    cursor = np.array(cursor)
    verts = [cursor-np.array((0.586,0.0,0.0)), cursor+np.array((0.586,0.0,0.8288))]
    edges = [(0,1)]
    faces = []
    coll = bpy.context.collection
    new_molecule = _mesh.create_object('new_molecule', coll, verts, edges, faces)
    Atoms, Bonds, Bond_Orders = ['C','C'], [(0,1)], [1]
    _mesh.add_base_radii_attr(new_molecule, Atoms, Bonds, Bond_Orders)
    _mesh.create_atom_materials({'C'})
    _shade.create_mat_hetero('Hetero')
    new_molecule.select_set(True)
    return new_molecule

def functional_group_coords(vert, bank, pitch, fg_name, rotate_angle, n_carbon):
    origin = np.array(vert.co)
    if fg_name == 'Benzyl':   # bond_length_CC: 1.4355
        C1 = origin
        C2 = C1 + bank * 0.71775 + pitch * 1.24318
        C3 = C2 + bank * 1.4355
        C4 = C1 + bank * 2.871
        C6 = C1 + bank * 0.71775 - pitch * 1.24318
        C5 = C6 + bank * 1.4355
        coords = [C1, C2, C3, C4, C5, C6]
    elif fg_name == 'Carboxyl':   # bond_length_CO: 1.3125
        C1 = origin
        O1 = C1 + bank * 0.65625 + pitch * 1.13666
        O2 = C1 + bank * 0.65625 - pitch * 1.13666
        coords = [C1, O1, O2]
    elif fg_name == 'Hexatomic Ring':
        C1 = origin
        bank_1 = _math.normalize(bank+2.824*pitch)
        bank_2 = _math.rotate_vec(bank_1, bank, pi*2/3)
        bank_3 = _math.rotate_vec(bank_1, bank, pi*4/3)
        C2 = C1 + bank_2 * 1.4355
        C3 = C2 + bank * 1.4355
        C4 = C3 + bank_3 * 1.4355
        C5 = C4 - bank_2 * 1.4355
        C6 = C5 - bank * 1.4355
        coords = [C1, C2, C3, C4, C5, C6]
    elif fg_name == 'n-Chain':
        C1 = origin
        bank_1 = _math.normalize(bank+2.824*pitch)
        banks = [bank_1, bank]
        coords = [C1]
        for i in range(n_carbon-1):
            new_C = origin + banks[i%2] * 1.4355
            coords.append(new_C)
            origin = new_C
        origin = C1
    coords = [_math.rotate_vec(coord-origin, bank, rotate_angle*pi/180)+origin for coord in coords]
    return coords