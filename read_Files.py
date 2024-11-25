import numpy as np
import re
from .Chem_Data import ELEMENTS_DEFAULT, BONDS_DEFAULT, SPACEGROUP_DEFAULTS
from . import _math, _mesh

# generate BONDS list from pure ATOMS list(.xyz and some .pdb/.cif format),
# judging whether a bond is formed based on the distance between atoms.
def add_BONDS(ATOMS, COORDS, factor):
    BONDS, BOND_ORDERS = [],[]
    for i, atom1 in enumerate(ATOMS):
        for j, atom2 in enumerate(ATOMS):
            if i < j:
                key = f"{atom1},{atom2}" if ELEMENTS_DEFAULT[atom1][0] <= ELEMENTS_DEFAULT[atom2][0] else f"{atom2},{atom1}"
                if key not in BONDS_DEFAULT: key = "Default"
                coord1 = np.array(COORDS[i])
                coord2 = np.array(COORDS[j])
                distance = np.linalg.norm(coord1-coord2)
                if distance <= BONDS_DEFAULT[key][3]*factor:
                    BONDS.append((i,j))
                    BOND_ORDERS.append(1)
    return (BONDS, BOND_ORDERS)

def read_MOL(file, filetype, mol_type, cursor, bond_length_f):
    ATOMS, COORDS = [],[]
    BONDS, BOND_ORDERS = [],[]
    molecules = []
    cursor = np.array(cursor)
    ################################################################################
    if mol_type == 'small_mol':
        # MOL: MDL Molfile
        if filetype == 'mol':
            tot_atoms = int(file[3][0:3])
            tot_bonds = int(file[3][3:6])
            for line in file[4: 4+tot_atoms]:
                ATOMS.append(line.split()[3])
                coord = np.array([float(value) for value in line.split()[0:3]])
                COORDS.append(cursor+coord)
            ATOMS = [atom if atom in ELEMENTS_DEFAULT else 'C' for atom in ATOMS]
            for line in file[4+tot_atoms: 4+tot_atoms+tot_bonds]:
                BONDS.append((int(line[0:3])-1, int(line[3:6])-1))
                BOND_ORDERS.append(int(line[6:9]))
            
        ################################################################################
        # SDF: structure-data format (V2000 or V3000)
        if filetype == 'sdf':
            subfiles = []
            start = 0
            count = 0
            # split multiple molecules into a subfiles list
            for line in file:
                count += 1
                if '$$$$' in line:
                    subfile = [line for line in file[start:start+count]]
                    start = count
                    subfiles.append(subfile)
                    continue
            # get the version of the sdf file
            file_version = file[3].split()[-1]
            molecules = []
            coord_shift = np.array([0.0, 0.0, 0.0])
            for file in subfiles:
                ATOMS, COORDS = [],[]
                BONDS, BOND_ORDERS, ATOM_LINKS = [],[],[]
                if file_version == 'V2000':
                    tot_atoms = int(file[3][0:3])
                    tot_bonds = int(file[3][3:6])
                    for line in file[4: 4+tot_atoms]:
                        ATOMS.append(line.split()[3])
                        coord = np.array([float(value) for value in line.split()[0:3]])
                        COORDS.append(cursor+coord+coord_shift)
                    for line in file[4+tot_atoms: 4+tot_atoms+tot_bonds]:
                        BONDS.append((int(line[0:3])-1, int(line[3:6])-1))
                        BOND_ORDERS.append(int(line[6:9]))

                elif file_version == 'V3000':
                    tot_atoms = int(file[5].split()[3])
                    tot_bonds = int(file[5].split()[4])
                    for line in file[7: 7+tot_atoms]:
                        ATOMS.append(line.split()[3])
                        coord = np.array([float(value) for value in line.split()[4:7]])
                        COORDS.append(cursor+coord+coord_shift)
                    for line in file[9+tot_atoms: 9+tot_atoms+tot_bonds]:
                        BONDS.append((int(line.split()[4])-1, int(line.split()[5])-1))
                        BOND_ORDERS.append(int(line.split()[3]))

                ATOMS = [atom if atom in ELEMENTS_DEFAULT else 'C' for atom in ATOMS]
                mol = [ATOMS, COORDS, BONDS, BOND_ORDERS]
                molecules.append(mol)
                coord_shift += np.array([10.0, 0.0, 0.0])

        ################################################################################
        # XYZ: a chemical file format storing the Cartesian coordinates of atoms.
        # The units are generally in angstroms.
        # Note: Only atomic info is present in the xyz file, thus bond order info
        # cannot be restored, and here all defaults to 1.
        # A 'Stick' or 'Space Filling' type is considered more suitable for this format.
        elif filetype == 'xyz':
            tot_atoms = int(re.sub('[\W_a-z_A-Z]+','',file[0]))
            count = 0
            for line in file[1:]:
                count += 1
                if line.split()[0][0] != '#':
                    break
            for line in file[count:count+tot_atoms]:
                text = line.split()[0] if line.split()[0].isalpha() else line.split()[1]
                coord_values = line.split()[1:4] if line.split()[0].isalpha() else line.split()[2:5]
                if text in ELEMENTS_DEFAULT:
                    ATOMS.append(text)
                    coord = np.array([float(value) for value in coord_values])
                    COORDS.append(cursor+coord)
            ATOMS = [atom if atom in ELEMENTS_DEFAULT else 'C' for atom in ATOMS]
            tot_atoms = len(ATOMS)
            BONDS, BOND_ORDERS = add_BONDS(ATOMS, COORDS, bond_length_f)
        
        ################################################################################
        # PDB: protein data bank format, a standard for files containing atomic coordinates.
        # if there is a CONECT block, then use it to get BONDS; if not, create one by add_BONDS.
        # Similarly, all bond orders default to 1 here.
        elif filetype == 'pdb':
            tot_atoms, tot_bonds = 0, 0
            HETATM, CONECT = [], []
            for line in file:
                if line.split() == []:
                    pass
                else:
                    if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                        HETATM.append(line)
                    #else:
                    #    continue
                    if line.split()[0] == 'CONECT':
                        CONECT.append(line)
                    else:
                        continue
            
            for line in HETATM:
                if line[12:13] == " " or line[12:13].isdigit() == True:
                    short_name = line[13:14]
                    if short_name.isdigit() == True:
                        short_name = line[14:15]
                    if line[14:15].islower() == True:
                        short_name = short_name + line[14:15]
                elif line[12:13].isupper() == True:
                    short_name = line[12:13]
                    if line[13:14].isalpha() == True:
                        short_name = short_name + line[13:14]
                else:
                    short_name = 'Default'
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:55])
                coord = np.array((x,y,z))
                ATOMS.append(short_name.lower().capitalize())
                COORDS.append(cursor+coord)
            ATOMS = [atom if atom in ELEMENTS_DEFAULT else 'C' for atom in ATOMS]
            tot_atoms = len(ATOMS)
            
            # a double for-loop is carried out, looping over each pair of atoms.
            # This may take some time. The same thing may happen to add_BONDS function.
            CONECT = list(set(CONECT))
            for line in CONECT:
                start_idx = int(line.split()[1])-1
                end_idxes = [int(v)-1 for v in line.split()[2:]]
                for end_idx in end_idxes:
                    if end_idx > start_idx:
                        BONDS.append((start_idx, end_idx))
                        BOND_ORDERS.append(1)
            
            if BONDS == []:
                BONDS, BOND_ORDERS = add_BONDS(ATOMS, COORDS, bond_length_f)
            
        ################################################################################
        elif filetype == 'json':
            atoms_dict = file['PC_Compounds'][0]['atoms']
            bonds_dict = file['PC_Compounds'][0]['bonds']
            coords_list = file['PC_Compounds'][0]['coords'][0]
            x_list = coords_list['conformers'][0]['x']
            y_list = coords_list['conformers'][0]['y']
            z_list = coords_list['conformers'][0]['z']

            elem_nums = atoms_dict['element']
            bonds_aid1 = bonds_dict['aid1']
            bonds_aid2 = bonds_dict['aid2']
            bonds_order = bonds_dict['order']
            
            # import atoms information: elements and coordinates
            for i,num in enumerate(elem_nums):
                atom = list(ELEMENTS_DEFAULT)[num]
                x = x_list[i]
                y = y_list[i]
                z = z_list[i]
                ATOMS.append(atom)
                COORDS.append(np.array((x,y,z))+cursor)
            ATOMS = [atom if atom in ELEMENTS_DEFAULT else 'C' for atom in ATOMS]
            tot_atoms = len(ATOMS)

            # import bonds information: bond_atoms_id and bond_orders
            for aid1, aid2, order in zip(bonds_aid1, bonds_aid2, bonds_order):
                BONDS.append((aid1-1,aid2-1))
                BOND_ORDERS.append(order)

            
        ################################################################################
        # a standard CIF(Crystallographic Information File) format file of small molecules
        elif filetype == 'cif':
            from gemmi import cif
            block = file.sole_block()
            # elements and coordinates
            atom_symbols = ["_chem_comp_atom.type_symbol", "_atom_site.label_atom_id", "_atom_site_type_symbol"]
            atom_Cartn_xs = ["_chem_comp_atom.pdbx_model_Cartn_x_ideal", "_atom_site.cartn_x", "_atom_site_fract_x"]
            atom_Cartn_ys = ["_chem_comp_atom.pdbx_model_Cartn_y_ideal", "_atom_site.cartn_y", "_atom_site_fract_y"]
            atom_Cartn_zs = ["_chem_comp_atom.pdbx_model_Cartn_z_ideal", "_atom_site.cartn_z", "_atom_site_fract_z"]

            symbol = [block.find_values(sym) for sym in atom_symbols if block.find_values(sym)]
            cartnx = [block.find_values(x) for x in atom_Cartn_xs if block.find_values(x)]
            cartny = [block.find_values(y) for y in atom_Cartn_ys if block.find_values(y)]
            cartnz = [block.find_values(z) for z in atom_Cartn_zs if block.find_values(z)]

            # bond aids and bond orders
            _bond_aids_1 = [cif.as_string(string) for string in block.find_values("_chem_comp_bond.atom_id_1")]
            _bond_aids_2 = [cif.as_string(string) for string in block.find_values("_chem_comp_bond.atom_id_2")]
            if not _bond_aids_1: _bond_aids_1 = [cif.as_int(value) for value in block.find_values("_struct_conn.ptnr1_atom_site_id")]
            if not _bond_aids_2: _bond_aids_2 = [cif.as_int(value) for value in block.find_values("_struct_conn.ptnr2_atom_site_id")]
            _bond_orders = [cif.as_string(string) for string in block.find_values("_chem_comp_bond.value_order")]

            try:
                _atom_symbol = [string for string in symbol[0]]
                _atom_Cartn_x = [cif.as_number(value) for value in cartnx[0]]
                _atom_Cartn_y = [cif.as_number(value) for value in cartny[0]]
                _atom_Cartn_z = [cif.as_number(value) for value in cartnz[0]]

                ATOMS = _atom_symbol
                ATOMS = [atom if atom in ELEMENTS_DEFAULT else 'C' for atom in ATOMS]
                tot_atoms = len(ATOMS)
                COORDS = [np.array((x,y,z))+cursor for x,y,z in zip(_atom_Cartn_x, _atom_Cartn_y, _atom_Cartn_z)]
                
                if _bond_aids_1:
                    if type(_bond_aids_1[0]) is str:
                        _element_id = [cif.as_string(string) for string in block.find_values("_chem_comp_atom.atom_id")]
                        _atom_id = [cif.as_int(value) for value in block.find_values("_chem_comp_atom.pdbx_ordinal")]
                        atom_ordinal_dic = {key:value for key,value in zip(_element_id, _atom_id)}
                        BONDS = [(atom_ordinal_dic[aid1]-1, atom_ordinal_dic[aid2]-1) for aid1, aid2 in zip(_bond_aids_1, _bond_aids_2)]
                    else:
                        BONDS = [(aid1-1, aid2-1) for aid1, aid2 in zip(_bond_aids_1, _bond_aids_2)]
                else:
                    BONDS, BOND_ORDERS = add_BONDS(ATOMS, COORDS, bond_length_f)

                if _bond_orders:
                    bond_order_dic = {'SING':1, 'DOUB':2, 'TRIP':3}
                    BOND_ORDERS = [bond_order_dic[order] for order in _bond_orders]
                else:
                    BOND_ORDERS = [1] * len(BONDS)

            except Exception as e:
                print('Sorry %s' % e)
            
        # else: other filetypes...
        ################################################################################
        if filetype != 'sdf': molecules.append([ATOMS, COORDS, BONDS, BOND_ORDERS])
        return molecules


def read_CIF(filepath, bond_length_f, boundary):
    from gemmi import cif
    doc = cif.read_file(filepath)
    block = doc.sole_block()

    length_a = cif.as_number(block.find_value("_cell_length_a"))
    length_b = cif.as_number(block.find_value("_cell_length_b"))
    length_c = cif.as_number(block.find_value("_cell_length_c"))
    angle_alpha = cif.as_number(block.find_value("_cell_angle_alpha"))
    angle_beta = cif.as_number(block.find_value("_cell_angle_beta"))
    angle_gamma = cif.as_number(block.find_value("_cell_angle_gamma"))
    space_groups = ["_symmetry_space_group_name_H-M","_space_group_name_Hall","_space_group_name_H-M_alt"]
    space_group = [block.find_value(x) for x in space_groups if block.find_value(x)][0]
    space_group = re.sub(" ", "", space_group)

    atom_sites_fx = []
    atom_sites_fy = []
    atom_sites_fz = []
    symop_operations = []
    atom_site_labels = []
    atom_type_symbols = []

    try:
        atom_sites_fx = [cif.as_number(value) for value in block.find_values("_atom_site_fract_x")]
        atom_sites_fy = [cif.as_number(value) for value in block.find_values("_atom_site_fract_y")]
        atom_sites_fz = [cif.as_number(value) for value in block.find_values("_atom_site_fract_z")]
        atom_site_labels = [cif.as_string(string) for string in block.find_values("_atom_site_label")]
        atom_type_symbols = [cif.as_string(string) for string in block.find_values("_atom_site_type_symbol")]
        symop_operations = [cif.as_string(string) for string in block.find_values("_space_group_symop_operation_xyz")]
        if not symop_operations:
            symop_operations = [cif.as_string(string) for string in block.find_values("_symmetry_equiv_pos_as_xyz")]
        if not symop_operations:
            symop_operations = SPACEGROUP_DEFAULTS[eval(space_group).capitalize()][3]
        
    except Exception as e:
        print(e)

    ATOMS = []
    FRACT_COORDS = []
    CARTN_COORDS = []
    atom_type_symbols = atom_type_symbols if atom_type_symbols else atom_site_labels # Asymmetric unit. 
    fx, fy, fz = atom_sites_fx, atom_sites_fy, atom_sites_fz
    atom_sites_fracts = [[fx[i], fy[i], fz[i]] for i in range(len(atom_type_symbols))]
    atom_sites_dic = {}
    for symbol, fract_xyz in zip(atom_type_symbols, atom_sites_fracts):
        atom_sym = re.sub('[\W0-9]', '', symbol)
        # Symmetry equivalent positions (non-normalized)
        sym_fracts = _math.fract_symop(fract_xyz, symop_operations)
        # After the coordinates translation, some of the atoms can overlap
        # remove duplicates (same position)
        sym_fracts = list(set([tuple(fract.tolist()) for fract in sym_fracts]))
        # coordinates normalization (e.g. 1.15 becomes 0.15, -0.25 becomes 0.75, etc)
        sym_fracts = _math.fracts_normalize(sym_fracts, boundary)
        
        if atom_sym in atom_sites_dic:
            for fract in sym_fracts:
                atom_sites_dic[atom_sym].append(fract)
        else:
            atom_sites_dic.update({atom_sym: sym_fracts})
    
    for key, value in atom_sites_dic.items():
        for fract in value:
            ATOMS.append(key)
            FRACT_COORDS.append(fract)
            CARTN_COORDS.append(_math.fract_to_cartn(fract, length_a, length_b, length_c, 
                                angle_alpha, angle_beta, angle_gamma))
    print(FRACT_COORDS)
    print(CARTN_COORDS)        
    BONDS, BOND_ORDERS = add_BONDS(ATOMS, CARTN_COORDS, bond_length_f)

    cell_lengths = [length_a, length_b, length_c]
    cell_angles = [angle_alpha, angle_beta, angle_gamma]

    molecules = [[ATOMS, CARTN_COORDS, BONDS, BOND_ORDERS],
                 [cell_lengths, cell_angles, space_group, symop_operations],
                 ]
    return molecules

    
def read_PDB(filepath, filename, coll):
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb
    file = pdb.PDBFile.read(filepath)
    biomol = pdb.get_structure(file, extra_fields = ['b_factor', 'charge', 'occupancy', 'atom_id'], include_bonds = True)
    if not biomol.bonds:
        biomol.bonds = struc.connect_via_distances(biomol[0], inter_residue=True)
    verts = biomol[0].coord
    edges = biomol.bonds.as_array()[:,[0,1]]
    Bond_Orders = [1]*len(edges)
    bond_radii = [1]*len(edges)
    object = _mesh.create_object(filename, coll, verts, edges, faces=[])
    elements = [e.capitalize() for e in list(set(biomol.element))]
    atom_ids = [ELEMENTS_DEFAULT[element.capitalize()][0] for element in biomol.element]
    object['Type'] = 'biomacro'
    object['Elements'] = elements

    atom_radii_B = [ELEMENTS_DEFAULT[element.capitalize()][4]/2 for element in biomol.element]
    atom_radii_S = [ELEMENTS_DEFAULT['Bond'][4]*1.8 for element in biomol.element]
    atom_radii_W = [0]*len(biomol.element)
    atom_radii_F = [ELEMENTS_DEFAULT[element.capitalize()][7] if ELEMENTS_DEFAULT[element.capitalize()][7] else 2 for element in biomol.element]
    atom_colors = []
    for element in biomol.element:
        for value in ELEMENTS_DEFAULT[element.capitalize()][3]:
            atom_colors.append(value)

    bond_radii_B = [ELEMENTS_DEFAULT['Bond'][4]]*len(Bond_Orders)
    for i, order in enumerate(Bond_Orders): bond_radii_B[i] = ELEMENTS_DEFAULT['Bond'][4] if order==1 else 0.8*ELEMENTS_DEFAULT['Bond'][4]
    bond_radii_S = [ELEMENTS_DEFAULT['Bond'][4]*1.8]*len(Bond_Orders)
    bond_radii_W = [ELEMENTS_DEFAULT['Bond'][5]]*len(Bond_Orders)
    bond_aids_1 = [ELEMENTS_DEFAULT[biomol.element[bond[0]].capitalize()][0] for bond in edges]
    bond_aids_2 = [ELEMENTS_DEFAULT[biomol.element[bond[1]].capitalize()][0] for bond in edges]
    attributes = (
        {'name': 'bond_order',    'type':'INT',      'domain':'EDGE',   'values':Bond_Orders},
        {'name': 'bond_aids_1',   'type':'INT',      'domain':'EDGE',   'values':bond_aids_1},
        {'name': 'bond_aids_2',   'type':'INT',      'domain':'EDGE',   'values':bond_aids_2},
        {'name': 'bond_radii_B',  'type':'FLOAT',    'domain':'EDGE',   'values':bond_radii_B},
        {'name': 'bond_radii_S',  'type':'FLOAT',    'domain':'EDGE',   'values':bond_radii_S},
        {'name': 'bond_radii_W',  'type':'FLOAT',    'domain':'EDGE',   'values':bond_radii_W},
        {'name': 'radii_B',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_B},
        {'name': 'radii_S',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_S},
        {'name': 'radii_W',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_W},
        {'name': 'radii_F',       'type':'FLOAT',    'domain':'POINT',  'values':atom_radii_F},
        {'name': 'colour',   'type':'FLOAT_COLOR',   'domain':'POINT',  'values':atom_colors},
    )

    # add attributes
    for attr in attributes:
        _mesh.add_attribute(object, attr['name'], attr['type'], attr['domain'], attr['values'])
    _mesh.add_attribute(object, 'bond_order', 'INT', 'EDGE', Bond_Orders)
    _mesh.add_attribute(object, 'bond_radius', 'FLOAT', 'EDGE', bond_radii)
    return (biomol, object)