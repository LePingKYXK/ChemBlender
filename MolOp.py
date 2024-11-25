import bpy
import os
import sys
import json
import numpy as np
from bpy.props import(StringProperty,
                      FloatProperty,
                      EnumProperty,
                      BoolProperty,
                      IntProperty,)

from . import _mesh, _node, _math, _shade, draw, read_Files, node_assets
from .Chem_Data import ELEMENTS_DEFAULT, lipophobicity, atom_names, atom_charge

# get molecule information from a file.
def get_mol_info():
    from gemmi import cif
    mytool = bpy.context.scene.my_tool
    filepath = mytool.datafile
    if os.path.exists(filepath):
        filename = os.path.basename(filepath).split('.')[0]
        filetype = os.path.basename(filepath).split('.')[-1]
        if filetype == 'json':
            context = open(filepath,'r')
            file = json.loads(context.read())
        elif filetype == 'cif':
            file = cif.read_file(filepath)
            try:
                length_a = file.sole_block().find_value("_cell_length_a")
            except:
                length_a = 0
        elif filetype in ['xyz','sdf','mol','pdb']:
            with open(filepath) as in_file: file = in_file.readlines()
        else:
            file = ''
        
        if filetype in ['xyz','sdf','mol','json']:
            mol_type = 'small_mol'
        elif filetype == 'cif':
            mol_type = 'crystal' if length_a else 'small_mol'
        elif filetype in ['pdb', 'pdbx', 'cif']:
            prefix = file[0].split(' ')[0].split('_')[0]
            mol_type = 'biomacro' if prefix == 'HEADER' or prefix == 'data' else 'small_mol'
        elif filetype in ['png','jpg','jpeg']:
            mol_type = 'image'
        else:
            mol_type = ''
        return (file, filename, filetype, mol_type)
    else:
        return None

def mol_scaffold_build(filename, mol_type, coll_mol, molecule):
    Atoms, Coords, Bonds, Bond_Orders = molecule
    print(Atoms)
    print(Coords)
    print(Bonds)
    print(Bond_Orders)
    _shade.create_mat_hetero('Hetero')
    _mesh.create_atom_materials(list(set(Atoms)))
    mol_scaffold = _mesh.create_object(filename, coll_mol, Coords, Bonds, [])
    mol_scaffold['Type'] = mol_type
    mol_scaffold['Elements'] = list(set(Atoms))
    # add attributes
    _mesh.add_base_radii_attr(mol_scaffold, Atoms, Bonds, Bond_Orders)

    bpy.context.view_layer.objects.active = mol_scaffold
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.remove_doubles()
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    if mol_type in ['small_mol','crystal']:
        atom_orientations = _mesh.atom_orient(mol_scaffold)
        bond_vertic_vecs = _mesh.bond_vertical_dir(mol_scaffold)
        _mesh.add_attribute(mol_scaffold, 'atom_orient', 'FLOAT_VECTOR', 'POINT', atom_orientations)
        _mesh.add_attribute(mol_scaffold, 'bond_vertvec', 'FLOAT_VECTOR', 'EDGE', bond_vertic_vecs)
        _mesh.is_aromatic(mol_scaffold)
    return mol_scaffold


class MESH_OT_Build(bpy.types.Operator):
    '''Operator for 'Build 3D Molecule' button in Panel'''
    bl_idname = 'mol3d.structure_build'
    bl_label = "创建结构" if 'zh_HAN' in bpy.context.preferences.view.language else "Build Structure"
    bl_options = {'REGISTER','UNDO','GRAB_CURSOR'}

    # small molecules
    struct_type: EnumProperty(
        name = '',
        default = 'B',
        items = [('B', '球棍模型', ''),
                 ('S', '棒状模型', ''),
                 ('W', '键线式', ''),
                 ('F', '空间填充模型', '')]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('B', 'Ball and Stick', ''),
                 ('S', 'Sticks', ''),
                 ('W', 'Wire Frame', ''),
                 ('F', 'Space Filling', '')]
    ) # type:ignore
    
    # biomolecules
    bio_style: EnumProperty(
        name = '',
        default = 'C',
        items = [('B', '球棍模型', ''),
                 ('F', '空间填充模型', ''),
                 ('C', '卡通模型', ''),
                 ('R', '管道模型', ''),]
                 # ('S', '表面模型', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('B', 'Ball and Stick', ''),
                 ('F', 'Space Filling', ''),
                 ('C', 'Cartoon', ''),
                 ('R', 'Ribbon', ''),]
                 # ('S', 'Surface', ''),]
    ) # type:ignore
    hide_H: BoolProperty(name='Hide Hydrogens', default=False) # type: ignore
    atom_scale: FloatProperty(name='', default = 1.0, min = 0.0, soft_max = 2.5) # type: ignore
    bond_scale: FloatProperty(name='', default = 1.0, min = 0.0, soft_max = 2.5) # type: ignore
    atom_subdiv: IntProperty(name='', default = 3, min = 0, soft_max = 5) # type: ignore
    bond_subdiv: IntProperty(name='', default = 12, min = 3, soft_max = 24) # type: ignore

    # crystal structure
    count_a: IntProperty(name='', default=1, min=1, soft_max=20) # type: ignore
    count_b: IntProperty(name='', default=1, min=1, soft_max=20) # type: ignore
    count_c: IntProperty(name='', default=1, min=1, soft_max=20) # type: ignore
    filter: StringProperty(name='', default='') # type: ignore
    grow_by_fract: FloatProperty(name='', default=0.0, min=0.0, max=1.0) #type: ignore
    grow_direct: IntProperty(name='', default=0, min=0, soft_max=6) #type: ignore
    bond_length_f: FloatProperty(name='', default=1.0, min=0.0, soft_max=2.0) #type: ignore
    boundary: FloatProperty(name='', default=0.0, min=0.0, soft_max=0.5) #type: ignore
    

    def draw(self, context):
        mol_type = get_mol_info()[3]
        if mol_type == 'small_mol':
            layout = self.layout
            row = layout.row()
            text = "结构类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Struct Type:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'struct_type')
            row = layout.row()
            text = "键长倍增系数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Length Factor:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bond_length_f')
            
        elif mol_type == 'crystal':
            layout = self.layout
            row = layout.row()
            text = "晶胞重复数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Cell Cycles:"
            row.label(text=text)
            row.prop(self, 'count_a')
            row.prop(self, 'count_b')
            row.prop(self, 'count_c')
            row = layout.row()
            text = "过滤原子：" if 'zh_HAN' in bpy.context.preferences.view.language else "Filter Atom:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'filter')
            row = layout.row()
            text = "按位置分数扩展:" if 'zh_HAN' in bpy.context.preferences.view.language else "Grow by Fract:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'grow_by_fract')
            row = layout.row()
            text = "直接扩展:" if 'zh_HAN' in bpy.context.preferences.view.language else "Grow Direct:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'grow_direct')
            row = layout.row()
            text = "键长倍增系数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Length Factor:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bond_length_f')
            row = layout.row()
            text = "边界连键:" if 'zh_HAN' in bpy.context.preferences.view.language else "Boundary Bonding:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'boundary')
        
        elif mol_type == 'biomacro':
            layout = self.layout
            row = layout.row()
            text = "结构类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Biomol Style:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bio_style')
    
    def create_Ball_Stick(self, filename, mol_scaffold):
        # add geometry node modifier
        GN_mol = _node.create_GeoNode(mol_scaffold, 'Struct_'+filename, 'GN_'+filename)  
        struct_type = self.bio_style if mol_scaffold['Type'] == 'biomacro' else self.struct_type
        atom_syms = mol_scaffold['Elements']
        _node.Ball_Stick_nodes(GN_mol.node_group, atom_scale=self.atom_scale, bond_scale=self.bond_scale, atom_subdiv=self.atom_subdiv, 
                               bond_subdiv=self.bond_subdiv, atom_syms=atom_syms, type=struct_type, hide_H=self.hide_H, 
                               color_type='H', dash_count=5, dash_radius=0.1, category='C2', display_mode='M2', hide_segs='')

    def execute(self,context):
        mytool = bpy.context.scene.my_tool
        filepath = mytool.datafile

        filter_text = self.filter
        for sep in ',;/| ': filter_text = filter_text.replace(sep, ' ')
        filters = filter_text.split()
        filters = [filter.capitalize() for filter in filters]

        if get_mol_info():
            mol_file, filename, filetype, mol_type = get_mol_info()
            coll_mol = bpy.data.collections.new('Struct_'+filename)
            bpy.context.scene.collection.children.link(coll_mol)
            cursor = bpy.context.scene.cursor.location

            if mol_type == 'small_mol':
                molecules = read_Files.read_MOL(mol_file, filetype, mol_type, cursor, self.bond_length_f)
                if filetype == 'sdf':
                    for molecule in molecules:
                        mol_scaffold = mol_scaffold_build(filename, mol_type, coll_mol, molecule)
                        self.create_Ball_Stick(filename, mol_scaffold)
                else:
                    mol_scaffold = mol_scaffold_build(filename, mol_type, coll_mol, molecules[0])
                    self.create_Ball_Stick(filename, mol_scaffold)
            
            elif mol_type == 'crystal':
                crys_molecules = read_Files.read_CIF(mytool.datafile, self.bond_length_f, self.boundary)
                crys_scaffold = mol_scaffold_build('Crystal_'+filename, mol_type, coll_mol, crys_molecules[0])
                length_a, length_b, length_c = crys_molecules[1][0]
                angle_alpha, angle_beta, angle_gamma = crys_molecules[1][1]
                space_group = crys_molecules[1][2]
                cell_edges = draw.draw_cell_edges(filename, coll_mol, length_a, length_b, length_c,
                                                  angle_alpha, angle_beta, angle_gamma, space_group)

                vec_a = _math.fract_to_cartn((1,0,0), length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                vec_b = _math.fract_to_cartn((0,1,0), length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                vec_c = _math.fract_to_cartn((0,0,1), length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                edge_vecs = (vec_a, vec_b, vec_c)
                cell_cycles = (self.count_a, self.count_b, self.count_c)
                _node.expand_scaffold(crys_scaffold, cell_edges, edge_vecs, self.grow_direct, self.grow_by_fract)  # convert_to_mesh
                _node.unit_cycles(crys_scaffold, filename, 'Crys_Cycles_'+filename, cell_cycles, edge_vecs)
                _node.crys_filter(crys_scaffold, filename, 'Crys_Filter_'+filename, filters)
                self.create_Ball_Stick(filename, crys_scaffold)
            
            elif mol_type == 'biomacro':
                biomol, bioscaffold = read_Files.read_PDB(filepath, 'Bio_'+filename, coll_mol)
                MyBiomol = BioMol(biomol, bioscaffold, calc_ss=True)
                chain_ids, count_res, charge, lipophobicity = MyBiomol.set_attributes()
                _shade.create_mat_hetero('Hetero')
                _shade.create_mat_id(chain_ids, 'ChainID', 'chain_id')
                _shade.create_mat_id([1,2,3], 'Cartoon', 'sec_struct')
                _shade.create_mat_rainbow([0,count_res])
                _shade.create_mat_charge(charge)
                _shade.create_mat_hydrophobic(lipophobicity)
                _mesh.create_atom_materials(list(set(MyBiomol.element)))

                if self.bio_style in ['B','F']:
                    self.create_Ball_Stick(filename, bioscaffold)
                else:
                    GN_biomol = _node.create_GeoNode(bioscaffold, 'Struct_'+filename, 'GN_'+filename)
                    group = GN_biomol.node_group
                    if self.bio_style == 'C':
                        _node.Cartoon_nodes(group, sharp=True, width=1.8, thickness=0.5, coil_radius=0.3, 
                                            resolution=6, subdivision=0, color_type='H', category='C2',
                                            display_mode='M2', hide_segs='')
                    elif self.bio_style == 'R':
                        _node.Ribbon_nodes(group, radius=1.6, resolution=4, vertices=6, subdivision=2, color_type='H', 
                                           category='C2', display_mode='M2', hide_segs='')
            elif mol_type == 'image':
                import torch
                import cv2
                from molscribe import MolScribe
                from rdkit import Chem
                path = sys.executable[:-14]+'\\lib\\site-packages\\checkpoints\\swin_base_char_aux_1m.pth'
                model = MolScribe(path, device=torch.device('cpu'))
                img_file = filepath
                output = model.predict_image_file(img_file)
                block = output['smiles']+output['molfile']
                m = Chem.MolFromMolBlock(block)
                basename = os.path.basename(filepath)
                mol_filepath = filepath.removesuffix(basename)+filename+'.mol'
                Chem.MolToMolFile(m, mol_filepath)
                with open(mol_filepath) as mol_file: file = mol_file.readlines()
                molecules = read_Files.read_MOL(file, 'mol', 'small_mol', cursor, self.bond_length_f)
                mol_scaffold = mol_scaffold_build(filename, 'small_mol', coll_mol, molecules[0])
                self.create_Ball_Stick(filename, mol_scaffold)
                    
            return{'FINISHED'}
        else:
            return{'CANCELLED'}
        


class BioMol():
    def __init__(self, biomol, bioscaffold, calc_ss):
        self.biomol = biomol
        self.bioscaffold = bioscaffold
        self.calc_ss = calc_ss

        self.element = [e.capitalize() for e in biomol.element]
        self.bonds_array = biomol.bonds.as_array()
        self.bond_types = self.bonds_array[:, 2]
        self.attr_res_id = biomol.res_id
        self.attr_atom_id = biomol.atom_id
        self.attr_b_factor = biomol.b_factor
        self.attr_occupancy = biomol.occupancy
        self.count = len(self.element)

    def comp_secondary_structure(self):
        from biotite.structure import annotate_sse, spread_residue_wise
        conv_sse_char_int = {'a': 1, 'b': 2, 'c': 3, '': 0}  # 'a' is alpha-helix, 'b' is beta-sheet, 'c' is loop
        char_sse = annotate_sse(self.biomol[0])
        int_sse = np.array([conv_sse_char_int[char] for char in char_sse], dtype=int)
        atom_sse = spread_residue_wise(self.biomol, int_sse)
        return atom_sse

    def set_attributes(self):
        import biotite.structure as struc

        def attr_chain_id():
            return np.searchsorted(np.unique(self.biomol.chain_id), self.biomol.chain_id)

        def attr_entity_id():
            try:
                return self.biomol.entity_id
            except Exception:
                return [-1]*len(self.biomol[0])

        def attr_ordi_num():
            return [ELEMENTS_DEFAULT[element.capitalize()][0] for element in self.element]

        def attr_lipophobicity():
            lipo = np.array(list(map(
                lambda x, y: lipophobicity.get(x,{"0":0}).get(y,0),
                self.biomol.res_name, self.biomol.atom_name
            )))
            return lipo
        
        def attr_charge():
            charge = np.array(list(map(
                lambda x, y: atom_charge.get(x,{"0":0}).get(y,0),
                self.biomol.res_name, self.biomol.atom_name
            )))
            return charge

        def attr_atom_name():
            atom_name = np.array(list(map(
                lambda x: atom_names.get(x,-1),
                self.biomol.atom_name
            )))
            return atom_name

        def attr_is_alpha():
            return np.isin(self.biomol.atom_name, 'CA')

        def attr_is_nucleic():
            return struc.filter_nucleotides(self.biomol)
        
        def attr_is_peptide():
            aa = struc.filter_amino_acids(self.biomol)
            con_aa = struc.filter_canonical_amino_acids(self.biomol)
            return aa | con_aa
        
        def attr_is_hetero():
            return self.biomol.hetero
        
        def attr_is_carb():
            return struc.filter_carbohydrates(self.biomol)
        
        def attr_is_solvent():
            return struc.filter_solvent(self.biomol)

        def attr_is_backbone():
            backbone_atom_names = [
                'N', 'C', 'CA', 'O',                     # peptide backbone atoms
                "P", "O5'", "C5'", "C4'", "C3'", "O3'",  # 'continuous' nucleic backbone atoms
                "O1P", "OP1", "O2P", "OP2",              # alternative names for phosphate O's
                "O4'", "C1'", "C2'", "O2'"               # remaining ribose atoms
            ]
            is_backbone = np.logical_and(
                np.isin(self.biomol.atom_name, backbone_atom_names),
                np.logical_not(struc.filter_solvent(self.biomol))
            )

            return is_backbone

        def attr_sec_struct():
            if hasattr(self.biomol, "sec_strcut"):
                return self.biomol.sec_struct
            if self.calc_ss:
                return self.comp_secondary_structure()

        attributes = (
            {'name': 'res_id',        'type':'INT',      'domain':'POINT',   'values':self.attr_res_id},
            {'name': 'chain_id',      'type':'INT',      'domain':'POINT',   'values':attr_chain_id()},
            {'name': 'entity_id',     'type':'INT',      'domain':'POINT',   'values':attr_entity_id()},
            {'name': 'sec_struct',    'type':'INT',      'domain':'POINT',   'values':attr_sec_struct()},
            {'name': 'atom_name',     'type':'INT',      'domain':'POINT',   'values':attr_atom_name()},
            {'name': 'atom_id',       'type':'INT',      'domain':'POINT',   'values':attr_ordi_num()},

            {'name': 'b_factor',      'type':'FLOAT',    'domain':'POINT',   'values':self.attr_b_factor},
            {'name': 'occupancy',     'type':'FLOAT',    'domain':'POINT',   'values':self.attr_occupancy},
            {'name': 'lipophobicity', 'type':'FLOAT',    'domain':'POINT',   'values':attr_lipophobicity()},
            {'name': 'charge',        'type':'FLOAT',    'domain':'POINT',   'values':attr_charge()},
            
            {'name': 'is_alpha',      'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_alpha()},
            {'name': 'is_nucleic',    'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_nucleic()},
            {'name': 'is_peptide',    'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_peptide()},
            {'name': 'is_hetero',     'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_hetero()},
            {'name': 'is_carb',       'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_carb()},
            {'name': 'is_solvent',    'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_solvent()},
            {'name': 'is_backbone',   'type':'BOOLEAN',  'domain':'POINT',   'values':attr_is_backbone()},
        )

        # assign the attibutes to the object
        _mesh.add_attribute(self.bioscaffold, 'bond_type', 'INT', 'EDGE', self.bond_types)

        for attr in attributes:
            try:
                _mesh.add_attribute(self.bioscaffold, attr['name'], attr['type'], attr['domain'], attr['values'])
            except Exception as e:
                print(e)

        res_ids = [res_id for i,res_id in enumerate(self.attr_res_id) if not attr_is_solvent()[i] and not attr_is_hetero()[i]]
        
        return (attr_chain_id(), max(res_ids), attr_charge(), attr_lipophobicity())