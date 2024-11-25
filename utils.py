import bpy, bmesh
import re
from bpy.types import Operator
from bpy.props import FloatProperty, IntProperty, BoolProperty, EnumProperty, StringProperty

from . import _mesh, _node, _math, _shade, draw
from .MolOp import get_mol_info, mol_scaffold_build
from .read_Files import read_CIF, add_BONDS
from .Chem_Data import ELEMENTS_DEFAULT, metals, SPACEGROUP_DEFAULTS
from . import version as v
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
def simplify_text(text):
    for sep in ',;/| ': text = text.replace(sep, ' ')
    elements = text.split()
    elements = [e.capitalize() if e.isalpha() else e for e in elements]
    elements = ['bond' if e == 'Bond' else e for e in elements]
    if 'all' in text.lower() or '*' in text:
        elements += ['atom', 'bond']
    if 'atom' in text.lower() or 'ball' in text.lower():
        elements += ['atom']
    if 'bond' in text.lower() or 'ball' in text.lower():
        elements += ['bond']
    # e.g. ['C','H']
    bond_syms = [[x.split('-')[0].capitalize(),x.split('-')[1].capitalize()] for x in elements if '-' in x]
    return (list(set(elements)), bond_syms)

def add_legends(obj, elements, show_legend):
    bpy.ops.object.mode_set(mode='OBJECT')
    FLAG = False
    for coll in bpy.data.collections:
        if coll.name == 'LEGENDS':
            FLAG = True
            break
    if not FLAG:
        coll = bpy.data.collections.new('LEGENDS')
        bpy.context.scene.collection.children.link(coll)
    else:
        coll = bpy.data.collections['LEGENDS']
    coll.hide_viewport = False
    coll.hide_render = False

    if show_legend:
        for i,element in enumerate(elements):
            FLAG = False
            for obj in coll.objects:
                if obj.name == 'Legend_'+element:
                    FLAG = True
                    break
            if not FLAG:
                bpy.ops.mesh.primitive_uv_sphere_add(location=(3, -(i+1)*3-3, 0), radius = ELEMENTS_DEFAULT[element][4]/2)
                legend = bpy.context.object
                legend.name = 'Legend_'+element
                legend.data.materials.append(bpy.data.materials[element])
                legend.users_collection[0].objects.unlink(legend)
                coll.objects.link(legend)
                
                bpy.ops.object.text_add(location=(5, -(i+1)*3-3.6, 0), radius = 2)
                text = bpy.context.object
                text.data.body = element
                text.name = 'Text_'+element
                text.data.materials.append(bpy.data.materials[element])
                coll.objects.link(text)
                bpy.context.collection.objects.unlink(text)
                bpy.ops.object.select_all(action='DESELECT')
        
    else:      
        if FLAG:
            bpy.data.collections['LEGENDS'].hide_viewport=True
            bpy.data.collections['LEGENDS'].hide_render=True
# ----------------------------------------------------------------------------------------
# #########################################################################################################
# #########################################################################################################
# Mol3D Structure Building
class SmartDispalyButton(Operator):
    bl_idname = "mol3d.mol_display"
    bl_label = "智能显示" if 'zh_HAN' in bpy.context.preferences.view.language else "Smart Display"
    bl_options = {'REGISTER','UNDO'}

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

    color_type: EnumProperty(
        name = '',
        default = 'H',
        items = [('H', '默认', ''),
                 ('C', '自定义', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('H', 'Hetero', ''),
                 ('C', 'Customize', ''),]
    ) # type:ignore
    
    text = "隐藏氢原子" if 'zh_HAN' in bpy.context.preferences.view.language else "Hide Hydrogens"
    hide_H: BoolProperty(name=text, default=False) # type: ignore
    text = "芳香环显示" if 'zh_HAN' in bpy.context.preferences.view.language else "Aromatic Ring"
    aromatic: BoolProperty(name=text, default=False) # type: ignore
    atom_scale: FloatProperty(name='', default = 1.0, min = 0.0, soft_max = 2.5) # type: ignore
    bond_scale: FloatProperty(name='', default = 1.0, min = 0.0, soft_max = 2.5) # type: ignore
    atom_subdiv: IntProperty(name='', default = 3, min = 1, soft_max = 5) # type: ignore
    bond_subdiv: IntProperty(name='', default = 12, min = 3, soft_max = 24) # type: ignore

    count_a: IntProperty(name = '', default = 1, min = 1, soft_max = 20) # type: ignore
    count_b: IntProperty(name = '', default = 1, min = 1, soft_max = 20) # type: ignore
    count_c: IntProperty(name = '', default = 1, min = 1, soft_max = 20) # type: ignore
    text = "启用过滤" if 'zh_HAN' in bpy.context.preferences.view.language else "Enable Filtering"
    filter_toggle: BoolProperty(name = text, default = True,) # type: ignore

    # hide_backbone: BoolProperty(name = 'Hide Backbone', default = False) # type: ignore
    hide_segs: StringProperty(name = '', default='') # type: ignore
    bio_style: EnumProperty(
        name = '',
        default = 'C',
        items = [('B', '球棍模型', ''),
                 ('F', '空间填充模型', ''),
                 ('C', '卡通模型', ''),
                 ('R', '管道模型', ''),
                 ('S', '表面模型', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('B', 'Ball and Stick', ''),
                 ('F', 'Space Filling', ''),
                 ('C', 'Cartoon', ''),
                 ('R', 'Ribbon', ''),
                 ('S', 'Surface', ''),]
    ) # type:ignore
    color_style: EnumProperty(
        name = '',
        default = 'H',
        items = [('H', '异质', ''),
                 ('C', '按链', ''),
                 ('R', '光谱', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('H', 'Hetero', ''),
                 ('C', 'Chain', ''),
                 ('R', 'Rainbow', ''),]
    ) # type:ignore
    surf_color_style: EnumProperty(
        name = '',
        default = 'H',
        items = [('H', '异质', ''),
                 ('C', '按链', ''),
                 ('R', '光谱', ''),]
                 #('E', '电荷', ''),
                 #('L', '疏水性', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('H', 'Hetero', ''),
                 ('C', 'Chain', ''),
                 ('R', 'Rainbow', ''),]
                 #('E', 'Electrostatic', ''),
                 #('L', 'Hydrophobic', ''),]
    ) # type:ignore

    category: EnumProperty(
        name = '',
        default = 'C2',
        items = [('C0', '仅蛋白', ''),
                 ('C1', '仅核酸', ''),
                 ('C2', '两者', '')]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('C0', 'Only Protein', ''),
                 ('C1', 'Only Nucleic', ''),
                 ('C2', 'Protein and Nucleic', '')]
    ) # type:ignore
    
    display_mode: EnumProperty(
        name = '',
        default = 'M2',
        items = [('M0', '仅主链', ''),
                 ('M1', '仅侧基', ''),
                 ('M2', '显示全部', '')]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('M0', 'Only Backbone', ''),
                 ('M1', 'Only Residues', ''),
                 ('M2', 'Show All', '')]
    ) # type:ignore

    surface_mode: EnumProperty(
        name = '',
        default = 'E',
        items = [('A', '溶剂可及表面', ''),
                 ('E', '溶剂排除表面', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('A', 'SAS', ''),
                 ('E', 'SEA', '')]
    ) # type:ignore
    
    text = "移除溶剂" if 'zh_HAN' in bpy.context.preferences.view.language else "Remove Solvent"
    remove_solvent: BoolProperty(name = text, default = True,) # type: ignore
    radius: FloatProperty(name = '', default = 1.6, min = 0.0, soft_max = 5.0) # type: ignore
    vertices: IntProperty(name = '', default = 12, min = 3, soft_max = 24) # type: ignore

    text = "锐利/平滑" if 'zh_HAN' in bpy.context.preferences.view.language else "Sharp / Smooth"
    sharp: BoolProperty(name = text, default = True) # type: ignore
    width: FloatProperty(name = '', default = 1.8, min = 0.0, soft_max = 5.0) # type: ignore
    thickness: FloatProperty(name = '', default = 0.5, min = 0.0, soft_max = 2.0) # type: ignore
    loop_radius: FloatProperty(name = '', default = 0.3, min = 0.0, soft_max = 2.0) # type: ignore
    resolution: IntProperty(name = '', default = 6, min =0, soft_max = 12) # type: ignore
    subdivision: IntProperty(name = '', default = 0, min = 0, soft_max = 4) # type: ignore

    probe_radius: FloatProperty(name='', default = 1.4, min = 0.0, soft_max = 2.5) # type: ignore
    surface_res: IntProperty(name = '', default = 1, min = 1, max = 5) # type: ignore
    smooth_iteration: IntProperty(name = '', default = 3, min = 0, soft_max = 5) # type: ignore
    blur_iteration: IntProperty(name = '', default = 0, min = 0, soft_max = 5) # type: ignore

    def draw(self, context):
        layout = self.layout
        try:
            mol_type = self.ao['Type']
        except:
            mol_type = get_mol_info()[3]
        if mol_type == 'small_mol':
            row = layout.row()
            text = "结构类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Struct Type:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'struct_type')
            row = layout.row()
            text = "着色类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Coloring Type:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'color_type')
            layout.prop(self, 'hide_H')
            layout.prop(self, 'aromatic')
        elif mol_type == 'crystal' or mol_type == 'polyhedra':
            layout = self.layout
            row = layout.row()
            text = "晶胞重复数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Cell Cycles:"
            row.label(text = text)
            row.prop(self, 'count_a')
            row.prop(self, 'count_b')
            row.prop(self, 'count_c')
            layout.prop(self, 'filter_toggle')
            row = layout.row()
            text = "着色类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Coloring Type:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'color_type')
        elif mol_type == 'biomacro':
            row = layout.row()
            text = "结构类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Biomol Style:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bio_style')
            if self.bio_style == 'S':
                row = layout.row()
                text = "表面计算方式:" if 'zh_HAN' in bpy.context.preferences.view.language else "Surface Mode:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'surface_mode')

            row = layout.row()
            text = "显示类别:" if 'zh_HAN' in bpy.context.preferences.view.language else "Show Category:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'category')
            row = layout.row()
            text = "显示模式:" if 'zh_HAN' in bpy.context.preferences.view.language else "Show Mode:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'display_mode')

            row = layout.row()
            text = "着色类型:" if 'zh_HAN' in bpy.context.preferences.view.language else "Coloring Type:"
            row.label(text=text)
            row.scale_x = 1.8
            if self.bio_style != 'S':
                row.prop(self, 'color_style')
            else:
                row.prop(self, 'surf_color_style')
            row = layout.row()
            text = "隐藏片段:" if 'zh_HAN' in bpy.context.preferences.view.language else "Hide Segments:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'hide_segs')
            row = layout.row()
            layout.prop(self, 'remove_solvent')

            if self.bio_style == 'R':
                row = layout.row()
                text = "半径:" if 'zh_HAN' in bpy.context.preferences.view.language else "Radius:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'radius')
                row = layout.row()
                text = "径向顶点数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Vertices:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'vertices')
            elif self.bio_style == 'C':
                layout.prop(self, 'sharp')
                row = layout.row()
                text = "宽度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Width:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'width')
                row = layout.row()
                text = "厚度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Thickness:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'thickness')
                row = layout.row()
                text = "管道半径:" if 'zh_HAN' in bpy.context.preferences.view.language else "Loop Radius:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'loop_radius')
            elif self.bio_style == 'S':
                row = layout.row()
                text = "滚球半径:" if 'zh_HAN' in bpy.context.preferences.view.language else "Probe Radius:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'probe_radius')
                row = layout.row()
                text = "分辨率:" if 'zh_HAN' in bpy.context.preferences.view.language else "Resolution:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'surface_res')
                row = layout.row()
                text = "平滑迭代次数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Smooth Iteration:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'smooth_iteration')
                row = layout.row()
                text = "模糊迭代次数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Blur Iteration:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'blur_iteration')


            if self.bio_style in ['R','C']:
                row = layout.row()
                text = "主链分辨率:" if 'zh_HAN' in bpy.context.preferences.view.language else "Resolution:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'resolution')
                row = layout.row()
                text = "表面细分:" if 'zh_HAN' in bpy.context.preferences.view.language else "Subdivision:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'subdivision')

        if mol_type != 'biomacro' or self.bio_style in ['B','F']:
            row = layout.row()
            text = "原子半径缩放:" if 'zh_HAN' in bpy.context.preferences.view.language else "Atom Scale:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'atom_scale')
            row = layout.row()
            text = "键半径缩放:" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Scale:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bond_scale')
            row = layout.row()
            text = "原子细分:" if 'zh_HAN' in bpy.context.preferences.view.language else "Atom Subdiv:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'atom_subdiv')
            row = layout.row()
            text = "键细分:" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Subdiv:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bond_subdiv')


    def execute(self, context):
        ao = bpy.context.object
        self.ao = ao
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'polyhedra':
                ao = bpy.data.objects['Crystal_'+ao.name[10:]]
                type = 'crystal'
            bpy.ops.object.mode_set(mode='OBJECT')
            atom_syms = ao['Elements']
            dash_count = 5
            dash_radius = 0.1
            if type == 'small_mol':
                # add attributes for new created structure
                new_atom_orientations = _mesh.atom_orient(ao)
                new_bond_vertic_vecs = _mesh.bond_vertical_dir(ao)
                _mesh.add_attribute(ao, 'atom_orient', 'FLOAT_VECTOR', 'POINT', new_atom_orientations)
                _mesh.add_attribute(ao, 'bond_vertvec', 'FLOAT_VECTOR', 'EDGE', new_bond_vertic_vecs)
                _mesh.is_aromatic(ao)
                geonode = _mesh.get_modifiers(ao, 'Struct_')
                group = geonode.node_group
                
                # toggle display styles
                _node.Ball_Stick_nodes(group, self.atom_scale, self.bond_scale, self.atom_subdiv, self.bond_subdiv, atom_syms, self.struct_type, self.hide_H, 
                                       self.color_type, dash_count, dash_radius, self.category, self.display_mode, self.hide_segs)
                bonds = _node.get_node(group, 'Bonds')
                sepa_bonds = _node.get_node(group, 'Separate Bonds')
                aromatic = _node.get_node(group, 'Aromatic')
                join_aromatic = _node.get_node(group, 'Join')
                bonds.inputs[5].default_value = True if self.aromatic else False
                sepa_bonds.inputs[3].default_value = True if self.aromatic else False
                if self.aromatic:
                    _node.nodes_link(group, aromatic, 0, join_aromatic, 0)

            elif type == 'crystal':
                geonode = _mesh.get_modifiers(ao, 'Crys_Cycles_')
                group = geonode.node_group
                node_a = _node.get_node(group, 'count_a')
                node_b = _node.get_node(group, 'count_b')
                node_c = _node.get_node(group, 'count_c')
                node_a.integer = self.count_a
                node_b.integer = self.count_b
                node_c.integer = self.count_c

                geonode = _mesh.get_modifiers(ao, 'Crys_Filter_')
                if geonode:
                    if self.filter_toggle:
                        geonode.show_viewport = True
                        geonode.show_render = True
                    else:
                        geonode.show_viewport = False
                        geonode.show_render = False

                geonode = _mesh.get_modifiers(ao, 'Struct_')
                group = geonode.node_group
            
                # toggle display styles
                _node.Ball_Stick_nodes(group, self.atom_scale, self.bond_scale, self.atom_subdiv, self.bond_subdiv, atom_syms, self.struct_type, self.hide_H, 
                                       self.color_type, dash_count, dash_radius, self.category, self.display_mode, self.hide_segs)
            elif type == 'biomacro':
                geonode = _mesh.get_modifiers(ao, 'Struct_')
                group = geonode.node_group
                if self.bio_style in ['B','F']:
                    
                    _node.Ball_Stick_nodes(group, self.atom_scale, self.bond_scale, self.atom_subdiv, self.bond_subdiv, atom_syms, self.bio_style, self.hide_H, 
                                           self.color_style, dash_count, dash_radius, self.category, self.display_mode, self.hide_segs)
                    node = _node.get_node(group, 'Remove Solvents')
                    if self.remove_solvent:
                        node.inputs[v.sw_sockets_in['switch_1']].default_value = True
                elif self.bio_style == 'C':
                    _node.Cartoon_nodes(group, self.sharp, self.width, self.thickness, self.loop_radius, 
                                        self.resolution, self.subdivision, self.color_style, self.category, self.display_mode, self.hide_segs)
                elif self.bio_style == 'R':
                    _node.Ribbon_nodes(group, self.radius, self.resolution, self.vertices, self.subdivision, 
                                       self.color_style, self.category, self.display_mode, self.hide_segs)
                elif self.bio_style == 'S':
                    _node.Surface_nodes(group, self.probe_radius, self.surface_mode, self.surface_res, self.smooth_iteration, self.surf_color_style, self.blur_iteration, 
                                        self.category, self.display_mode, self.hide_segs)

            bpy.ops.object.select_all(action='DESELECT')
            return {'FINISHED'}

# #########################################################################################################
# #########################################################################################################
# Mol3D Tools
class SelectButton(Operator):
    bl_idname = "mol3d.select"
    bl_label = "选择" if 'zh_HAN' in bpy.context.preferences.view.language else "Select"
    bl_options = {'REGISTER','UNDO'}

    def execute(self, context):
        mytool  = context.scene.my_tool
        text = mytool.select_text
        text, bond_syms = simplify_text(text)
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type in ['small_mol','crystal','biomacro']:
                _mesh.deselect_all(ao)
                _mesh.select_verts(ao, text)
                _mesh.select_edges(ao, bond_syms)
                _mesh.select_bond_orders(ao, text)
                if 'atom' in text: _mesh.select_all(ao, 'VERT')
                if 'bond' in text: _mesh.select_all(ao, 'EDGE')
                if 'all' in text or '*' in text: _mesh.select_all(ao, 'VERT')

            return {'FINISHED'}


class ScaleButton(Operator):
    bl_idname = "mol3d.scale"
    bl_label = "缩放" if 'zh_HAN' in bpy.context.preferences.view.language else "Scale"
    bl_options = {'REGISTER','UNDO'}

    text = "缩放比例:" if 'zh_HAN' in bpy.context.preferences.view.language else "Scale:"
    scale: FloatProperty(name=text, default=1.0, min=0.0, soft_max=10.0) # type:ignore
    text = "绝对尺寸" if 'zh_HAN' in bpy.context.preferences.view.language else "Absolute Dimensions"
    abs_dimensions: BoolProperty(name = text, default = False) # type:ignore

    def execute(self,context):
        mytool  = context.scene.my_tool
        text = mytool.select_text
        text, bond_syms = simplify_text(text)
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type in ['small_mol','crystal','biomacro']:
                bpy.ops.object.mode_set(mode='OBJECT')
                if 'atom' in text: text = list(set(text+ao['Elements']))
                atom_id = _mesh.get_attirbute(ao, 'atom_id', 'INT', 'VERT')
                bond_aids_1 = _mesh.get_attirbute(ao, 'bond_aids_1', 'INT', 'EDGE')
                bond_aids_2 = _mesh.get_attirbute(ao, 'bond_aids_2', 'INT', 'EDGE')
                bond_orders = _mesh.get_attirbute(ao, 'bond_order', 'INT', 'EDGE')
                Elements = list(ELEMENTS_DEFAULT.keys())
                for radii in ['radii_B', 'radii_S', 'radii_F']:
                    radii_list = _mesh.get_attirbute(ao, radii, 'FLOAT', 'VERT')
                    new_radii_list = []
                    for i,radius in enumerate(radii_list):
                        element = Elements[atom_id[i]]
                        if element not in text:
                            new_radii_list.append(radius)
                        else:
                            if self.abs_dimensions:
                                new_radii_list.append(self.scale)
                            else:
                                new_radii_list.append(radius*self.scale)
                    _mesh.add_attribute(ao, radii, 'FLOAT', 'VERT', new_radii_list)

                for bond_radii in ['bond_radii_B', 'bond_radii_S', 'bond_radii_W']:
                    bond_radii_list = _mesh.get_attirbute(ao, bond_radii, 'FLOAT', 'EDGE')
                    new_bond_radii_list = []
                    for i,bond_radius in enumerate(bond_radii_list):
                        atom1 = Elements[bond_aids_1[i]]
                        atom2 = Elements[bond_aids_2[i]]
                        bond_order = bond_orders[i]
                        if [atom1,atom2] in bond_syms or [atom2,atom1] in bond_syms or f'{bond_order}' in text or 'bond' in text:
                            if self.abs_dimensions:
                                new_bond_radii_list.append(self.scale)
                            else:
                                new_bond_radii_list.append(bond_radius*self.scale)
                        else:
                            new_bond_radii_list.append(bond_radius)
                    _mesh.add_attribute(ao, bond_radii, 'FLOAT', 'EDGE', new_bond_radii_list)

            # bpy.ops.object.select_all(action='DESELECT')
            return {'FINISHED'}


class DefaultScaleButton(Operator):
    bl_idname = "mol3d.default_size"
    bl_label = "Default Size"
    bl_options = {'REGISTER','UNDO'}

    def execute(self,context):
        mytool  = context.scene.my_tool
        text = mytool.select_text
        text, bond_syms = simplify_text(text)
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type in ['small_mol','crystal','biomacro']:
                bpy.ops.object.mode_set(mode='OBJECT')
                if 'atom' in text: text = list(set(text+ao['Elements']))
                atom_id = _mesh.get_attirbute(ao, 'atom_id', 'INT', 'VERT')
                bond_aids_1 = _mesh.get_attirbute(ao, 'bond_aids_1', 'INT', 'EDGE')
                bond_aids_2 = _mesh.get_attirbute(ao, 'bond_aids_2', 'INT', 'EDGE')
                bond_orders = _mesh.get_attirbute(ao, 'bond_order', 'INT', 'EDGE')
                Elements = list(ELEMENTS_DEFAULT.keys())
                for radii in ['radii_B', 'radii_S', 'radii_F']:
                    radii_list = _mesh.get_attirbute(ao, radii, 'FLOAT', 'VERT')
                    new_radii_list = []
                    for i,radius in enumerate(radii_list):
                        element = Elements[atom_id[i]]
                        if element not in text:
                            new_radii_list.append(radius)
                        else:
                            if radii == 'radii_B':
                                new_radii_list.append(ELEMENTS_DEFAULT[element][4]/2)
                            elif radii == 'radii_S':
                                new_radii_list.append(ELEMENTS_DEFAULT['Bond'][4]*1.8)
                            else:
                                new_radii_list.append(ELEMENTS_DEFAULT[element][7])
                    _mesh.add_attribute(ao, radii, 'FLOAT', 'VERT', new_radii_list)

                for bond_radii in ['bond_radii_B', 'bond_radii_S', 'bond_radii_W']:
                    bond_radii_list = _mesh.get_attirbute(ao, bond_radii, 'FLOAT', 'EDGE')
                    new_bond_radii_list = []
                    for i,bond_radius in enumerate(bond_radii_list):
                        atom1 = Elements[bond_aids_1[i]]
                        atom2 = Elements[bond_aids_2[i]]
                        bond_order = bond_orders[i]
                        if [atom1,atom2] in bond_syms or [atom2,atom1] in bond_syms or f'{bond_order}' in text or 'bond' in text:
                            if bond_radii == 'bond_radii_B':
                                multiply = 1 if bond_order == 1 else 0.8
                                new_bond_radii_list.append(ELEMENTS_DEFAULT['Bond'][4]*multiply)
                            elif bond_radii == 'bond_radii_S':
                                new_bond_radii_list.append(ELEMENTS_DEFAULT['Bond'][4]*1.8)
                            else:
                                new_bond_radii_list.append(ELEMENTS_DEFAULT['Bond'][5])
                        else:
                            new_bond_radii_list.append(bond_radius)
                    _mesh.add_attribute(ao, bond_radii, 'FLOAT', 'EDGE', new_bond_radii_list)
            bpy.ops.object.select_all(action='SELECT')
            bpy.ops.object.select_all(action='DESELECT')
            return {'FINISHED'}


class AddLegendsButton(Operator):
    bl_idname = 'mol3d.add_legends'
    bl_label = 'Show Legends'
    bl_options = {'REGISTER','UNDO'}

    text = "显示原子标注" if 'zh_HAN' in bpy.context.preferences.view.language else "Show Atoms Legend"
    show_legend: BoolProperty(name = text, default = True) # type: ignore

    def draw(self,context):
        layout = self.layout
        layout.prop(self, 'show_legend')

    def execute(self, context):
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                elements = ao['Elements']
            except:
                elements = None
            if elements:
                add_legends(ao, elements, self.show_legend)
                return {'FINISHED'}     


class ConnectDisplayButton(Operator):
    bl_idname = "mol3d.showconnect"
    bl_label = "Connect Display Toggle"
    bl_options = {'REGISTER','UNDO'}

    def execute(self, context):
        mytool = context.scene.my_tool
        text = mytool.select_text
        text, bond_syms = simplify_text(text)
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type in ['small_mol','crystal','biomacro']:
                _mesh.deselect_all(ao)
                _mesh.select_edges(ao, bond_syms)
                _mesh.select_bond_orders(ao, text)
                if 'bond' in text: _mesh.select_all(ao, 'EDGE')
                sel_edge_idxs = _mesh.get_sel_idxs(ao, 'EDGE')
                Bond_Orders = _mesh.get_attirbute(ao, 'bond_order', 'INT', 'EDGE')
                for idx in sel_edge_idxs: Bond_Orders[idx] *= -1
                _mesh.add_attribute(ao, 'bond_order', 'INT', 'EDGE', Bond_Orders)
                
            return {'FINISHED'}

class DistanceButton(Operator):
    bl_idname = "mol3d.button_distance"
    bl_label = "Measure Distance"
    bl_description = "Measure the distance between two atoms"

    def execute(self, context):
        mytool = context.scene.my_tool
        ao = bpy.context.object
        if not ao or bpy.context.mode=='OBJECT':
            return {'CANCELLED'}
        else:
            dist = _mesh.calc_length()
            # Put the distance into the string of the output field.
            mytool.distance = dist
            return {'FINISHED'}
        
class AngleButton(Operator):
    bl_idname = "mol3d.button_angle"
    bl_label = "Measure Angle"
    bl_description = "Measure the angle of two edges"

    def execute(self, context):
        mytool = context.scene.my_tool
        ao = bpy.context.object
        if not ao or bpy.context.mode=='OBJECT':
            return {'CANCELLED'}
        else:
            angle = _mesh.calc_angle()
            # Put the angle into the string of the output field.
            mytool.angle = angle
            return {'FINISHED'}
        

# #########################################################################################################
# #########################################################################################################
# #########################################################################################################
# #########################################################################################################
# Crystal3D Tools
class PolyhedraBuildButton(Operator):
    bl_idname = "mol3d.polyhedra_build"
    bl_label = "创建多面体" if 'zh_HAN' in bpy.context.preferences.view.language else "Build Polyhedra"
    bl_options = {'REGISTER','UNDO'}

    build_mode: EnumProperty(
        name = 'Type',
        default = '0',
        items = [('0','默认',''),
                 ('1','自定义',''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('0','Default',''),
                 ('1','Customize',''),]
    ) # type:ignore

    center: StringProperty(name = "",) # type: ignore
    ligand: StringProperty(name = "",) # type: ignore
    text = "默认键长" if 'zh_HAN' in bpy.context.preferences.view.language else "R_Defaults"
    r_default: BoolProperty(name = text, default = True) # type: ignore
    bond_length_f: FloatProperty(name='', default=1.0, min=0.0, soft_max=2.0) #type: ignore
    homo: BoolProperty(name = "Homoligands", default = True) # type: ignore
    RMin: FloatProperty(name="RMin", default=0.0, min = 0.0,) # type: ignore
    RMax: FloatProperty(name="RMax", default=0.0, min = 0.0,) # type: ignore

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'build_mode')
        if self.build_mode == '1':
            row = layout.row(align=True)
            col = row.column(align=True)
            col.label(text="")
            col.prop(self, 'r_default')
            
            col = row.column(align=True)
            text = "中心原子" if 'zh_HAN' in bpy.context.preferences.view.language else "Center"
            col.label(text=text)
            col.prop(self, 'center')
            col = row.column(align=True)
            text = "配位原子" if 'zh_HAN' in bpy.context.preferences.view.language else "Ligand"
            col.label(text=text)
            col.prop(self, 'ligand')

            if not self.r_default:
                row = layout.row(align=True)
                row.prop(self, 'homo')
                row.prop(self, 'RMin')
                row.prop(self, 'RMax')
        row = layout.row()
        text = "键长倍增系数" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Length Factor:"
        row.label(text=text)
        row.scale_x = 1.5
        row.prop(self, 'bond_length_f')
                
    def set_polyhedra_base_attr(self, polyhedra):
        atom_ids = _mesh.get_attirbute(polyhedra, 'atom_id', 'INT', 'VERT')
        Atoms = [list(ELEMENTS_DEFAULT.keys())[id] for id in atom_ids]
        Bonds = []
        bm = bmesh.new()
        bm.from_mesh(polyhedra.data)
        for edge in bm.edges:
            Bonds.append([edge.verts[0].index,edge.verts[1].index])
        Bond_Orders = [1]*len(Bonds)
        _mesh.add_base_radii_attr(polyhedra, Atoms, Bonds, Bond_Orders)

    def execute(self, context):
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            bpy.ops.object.mode_set(mode='OBJECT')
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'crystal':
                sign = ao.name.find('_') + 1
                crys_name = ao.name[sign:].split('.')[0]
                coll = bpy.data.collections['Struct_'+crys_name]
                if self.build_mode == '0':
                    centers = ao['Elements']
                    centers = [e for e in centers if e in metals]
                    ligands=[element for element in ao['Elements'] if element in ['O', 'N', 'P', 'S']]
                else:
                    centers = simplify_text(self.center)[0]
                    ligands = simplify_text(self.ligand)[0]
                if 'atom' in centers: centers = ao['Elements']
                if 'atom' in ligands: ligands = ao['Elements']

                # create materials for polyhedras
                for center in centers:
                    FLAG = False    # ctr_center material is not existing
                    for material in bpy.data.materials:
                        if material.name == 'ctr_'+center:
                            FLAG = True
                            break
                    if FLAG:
                        ctr_mat = bpy.data.materials['ctr_'+center]
                    else:
                        ctr_mat = bpy.data.materials.new(name='ctr_'+center)
                    ctr_mat.use_nodes = True
                    try:
                        _BSDF = ctr_mat.node_tree.nodes["Principled BSDF"]
                    except:
                        _BSDF = ctr_mat.node_tree.nodes["原理化BSDF"]
                    _BSDF.inputs[0].default_value = ELEMENTS_DEFAULT[center][3]
                    _BSDF.inputs[2].default_value = 0.2  # Roughness
                    _BSDF.inputs[4].default_value = 0.2  # Alpha

                # add 'poly_aid' attribute to crys_scaffold object and convert 
                _node.add_attr_poly_aid(ao, centers, ligands, self.bond_length_f, self.RMax, self.RMin, self.r_default)
                
                # create polyhedra object from crys_scaffold mesh data
                polyhedra = _mesh.create_obj_from_mesh(coll, ao.data)
                polyhedra.name = 'Polyhedra_'+crys_name
                polyhedra['Type'] = 'polyhedra'
                _node.crys_polyhedra(polyhedra, crys_name, coll, centers, ligands, self.r_default, self.RMax, self.RMin, self.bond_length_f)
                self.set_polyhedra_base_attr(polyhedra)
                _mesh.dissolve_flat_edges(polyhedra)

                # add Unit Cycles geometry node_group
                GN_cycles = polyhedra.modifiers.new(crys_name+'_Polyhedra_Cycles', 'NODES')
                GN_cycles.node_group = bpy.data.node_groups[ao.modifiers[0].node_group.name]
                # add polyhedra edges geometry node_group
                _node.crys_polyhedra_edges(polyhedra, crys_name, bond_subdiv=12, radii=0.6)
                
            bpy.ops.object.select_all(action='DESELECT')
            return {'FINISHED'}


class PolyhedraEditButton(Operator):
    bl_idname = "mol3d.polyhedra_edit"
    bl_label = "多面体编辑" if 'zh_HAN' in bpy.context.preferences.view.language else "Edit Polyhedra"
    bl_options = {'REGISTER','UNDO'}

    edge_subdiv: IntProperty(name='', default=12, min=3, soft_max=24) # type: ignore
    edge_radius: FloatProperty(name='', default=0.6, min=0.0, soft_max=2.0) # type: ignore
    center_scale: FloatProperty(name='', default=1.0, min=0.0, soft_max=5.0) # type: ignore
    ligand_scale: FloatProperty(name='', default=1.0, min=0.0, soft_max=5.0) # type: ignore
    text = "显示多面体骨架" if 'zh_HAN' in bpy.context.preferences.view.language else "Show Polyhedra Skeleton"
    show_skeleton: BoolProperty(name=text, default=True) # type: ignore
    transparent: FloatProperty(name='', default=0.2, min=0.0, max=1.0) # type: ignore

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        text = "棱边半径:" if 'zh_HAN' in bpy.context.preferences.view.language else "Edge Radius:"
        row.label(text=text)
        row.prop(self, 'edge_radius')
        row = layout.row(align=True)
        text = "边细分数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Edge Subdiv:"
        row.label(text=text)
        row.prop(self, 'edge_subdiv')

        row = layout.row(align=True)
        text = "中心原子缩放:" if 'zh_HAN' in bpy.context.preferences.view.language else "Center Scale:"
        row.label(text=text)
        row.prop(self, 'center_scale')
        row = layout.row(align=True)
        text = "配位原子缩放:" if 'zh_HAN' in bpy.context.preferences.view.language else "Ligand Scale:"
        row.label(text=text)
        row.prop(self, 'ligand_scale')

        layout.prop(self, 'show_skeleton')

    def execute(self, cointext):
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'polyhedra':
                geonode = _mesh.get_modifiers(ao, 'Polyhedra_Edges_')
                bonds = _node.get_node(geonode.node_group, 'Bonds')
                bonds.inputs[3].default_value = self.edge_subdiv
                bonds.inputs[2].default_value = self.edge_radius

                sign = ao.name.find('_') + 1
                crys_scaffold = bpy.data.objects['Crystal_'+ao.name[sign:]]
                geonode = _mesh.get_modifiers(crys_scaffold, 'Struct_')
                atoms = _node.get_node(geonode.node_group, 'Atoms')
                
                atoms.inputs[5].default_value = self.ligand_scale
                atoms.inputs[4].default_value = self.center_scale
                bonds = _node.get_node(geonode.node_group, 'Bonds')
                bonds.inputs[4].default_value = True if self.show_skeleton else False

            return {'FINISHED'}


# calculate an average fract_xyz for selected vertices in the crystal scaffold
class AvgFractButton(Operator):
    bl_idname = "mol3d.button_avgfract"
    bl_label = "Calc. Average Fract"
    bl_options = {'REGISTER','UNDO'}

    def execute(self, context):
        mytool = context.scene.my_tool
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'crystal':
                sign = ao.name.find('_')+1
                crys_name = ao.name[sign:].split('.')[0]
                try:
                    cell_edges = _mesh.get_object('Cell_Edges_'+crys_name)
                    length_a, length_b, length_c = eval(cell_edges['cell lengths'])
                    angle_alpha, angle_beta, angle_gamma = eval(cell_edges['cell angles'])
                except Exception as e:
                    print(e)
                
                if ao.mode == 'EDIT':
                    avg_xyz = _mesh.avg_dummy_xyz(ao)
                    avg_fract = _math.cartn_to_fract(avg_xyz, length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                    fract = (float('%.3f' % avg_fract[0]), float('%.3f' % avg_fract[1]), float('%.3f' % avg_fract[2]))
                    fract_str = f'{fract[0]}, {fract[1]}, {fract[2]}'
                    mytool.avgfract = fract_str
            return {'FINISHED'}

class DummyButton(Operator):
    bl_idname = "mol3d.button_dummy"
    bl_label = "Create Dummy"
    bl_options = {'REGISTER','UNDO'}

    def execute(self, context):
        mytool = context.scene.my_tool
        filepath = mytool.datafile
        avg_fract = mytool.avgfract
        if avg_fract:
            if avg_fract != 'x, y, z':
                fract_xyz = eval(avg_fract)

        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'crystal':
                sign = ao.name.find('_')+1
                crys_name = ao.name[sign:].split('.')[0]
                try:
                    cell_edges = _mesh.get_object('Cell_Edges_'+crys_name)
                    length_a, length_b, length_c = eval(cell_edges['cell lengths'])
                    angle_alpha, angle_beta, angle_gamma = eval(cell_edges['cell angles'])
                    space_group = cell_edges['space group']
                    space_group = re.sub(" ", "", space_group)
                    space_group = eval(space_group).capitalize()
                    if filepath:
                        symop_operations = read_CIF(filepath, bond_length_f=1.0, boundary=0.0)[1][3]
                    else:
                        symop_operations = SPACEGROUP_DEFAULTS[space_group][3]
                except Exception as e:
                    print(e)
                
                sym_fracts = _math.fract_symop(fract_xyz, symop_operations)
                sym_fracts = list(set([tuple(fract.tolist()) for fract in sym_fracts]))
                sym_fracts = _math.fracts_normalize(sym_fracts, boundary=0)
                sym_cartns = [_math.fract_to_cartn(fract, length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma) for fract in sym_fracts]

                bpy.ops.object.mode_set(mode='OBJECT')
                mesh = bpy.data.meshes.new('Dummy')
                dummy_obj = bpy.data.objects.new(mesh.name, mesh)
                bpy.context.collection.objects.link(dummy_obj)
                mesh.from_pydata(sym_cartns,[],[])
                bpy.ops.object.select_all(action='DESELECT')
                ao.select_set(True)
                dummy_obj.select_set(True)
                bpy.ops.object.join()

                ao = bpy.context.object
                bm = bmesh.new()
                bm.from_mesh(ao.data)
                ordi_num_layer = bm.verts.layers.int.get('atom_id')
                for vert in bm.verts:
                    if vert[ordi_num_layer] == 0:
                        vert[ordi_num_layer] = -1
                bm.to_mesh(ao.data)

                bpy.ops.object.mode_set(mode='EDIT')
                FLAG = False
                for mat in bpy.data.materials:
                    if mat.name == 'Dummy':
                        FLAG = True
                        break
                if not FLAG: _mesh.create_atom_materials('Dummy')

            return {'FINISHED'}

class CrysSectionButton(Operator):
    bl_idname = "mol3d.crystal_section"
    bl_label = "晶面显示" if 'zh_HAN' in bpy.context.preferences.view.language else "Crystal Section Display"
    bl_options = {'REGISTER','UNDO'}

    h: IntProperty(name = '', default = 1) # type: ignore
    k: IntProperty(name = '', default = 1) # type: ignore
    l: IntProperty(name = '', default = 1) # type: ignore
    bottom_section: FloatProperty(name='', default=0.0)  # type: ignore
    top_section: FloatProperty(name='', default=10.0)  # type: ignore

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        text = "晶面指数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Miller Indices:"
        row.label(text = text)
        row.scale_x = 0.685
        row.prop(self, 'h')
        row.prop(self, 'k')
        row.prop(self, 'l')
        row = layout.row()
        text = "底部截距:" if 'zh_HAN' in bpy.context.preferences.view.language else "Bottom Intercept:"
        row.label(text = text)
        row.prop(self, 'bottom_section')
        row = layout.row()
        text = "顶部截距:" if 'zh_HAN' in bpy.context.preferences.view.language else "Top Intercept:"
        row.label(text = text)
        row.prop(self, 'top_section')

    def execute(self, context):
        ao = context.object
        if not ao:
            return {'CANCELLED'}
        else:
            bpy.ops.object.mode_set(mode='OBJECT')
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'crystal':
                sign = ao.name.find('_') + 1
                crys_name = ao.name[sign:].split('.')[0]
                hkl = [self.h, self.k, self.l]
                _node.crys_face_section(ao, crys_name, self.bottom_section, self.top_section, hkl)
        
        return {'FINISHED'}
    

# ------------------ Here is Crystal Structure Deformation related --------------------
class UpdateDeformButton(Operator):
    bl_idname = "mol3d.deformation_add"
    bl_label = "变形参数" if 'zh_HAN' in bpy.context.preferences.view.language else "Defomation Setting"
    bl_options = {'REGISTER','UNDO'}

    roll_angle: FloatProperty(name='', default = 360, soft_min = 0, soft_max = 360) # type: ignore
    align_axis: EnumProperty(
        name='',
        default = 'y',
        items = [('x','X',''),
                 ('y','Y',''),
                 ('z','Z','')]) # type: ignore
    transform1_x: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform1_y: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform1_z: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform2_x: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform2_y: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform2_z: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform3_x: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform3_y: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    transform3_z: FloatProperty(name = '', default = 1.0, soft_min = 0.0, soft_max = 2.5) # type: ignore
    text = "切换" if 'zh_HAN' in bpy.context.preferences.view.language else "Toggle"
    toggle1: BoolProperty(name=text, default = False) # type: ignore
    toggle2: BoolProperty(name=text, default = False) # type: ignore
    text = "逆向" if 'zh_HAN' in bpy.context.preferences.view.language else "Clockwise"
    clockwise: BoolProperty(name=text, default = False) # type: ignore
    
    noise_dir: EnumProperty(
        name='',
        default = 'z',
        items = [('x','X',''),
                 ('y','Y',''),
                 ('z','Z','')]) # type: ignore
    wavelength: FloatProperty(name = '', default = 15.0, soft_min = 0.0, soft_max = 50.0) # type: ignore
    noise_scale: FloatProperty(name = '', default = 0.05, soft_min = -0.5, soft_max = 0.5) # type: ignore
    amplitude: FloatProperty(name = '', default = 1.5, soft_min = 0.0, soft_max = 25.0) # type: ignore
    phase_angle: FloatProperty(name = '', default = 0.0, soft_min = -20.0, soft_max = 20.0) # type: ignore


    def draw(self, context):
        mytool = context.scene.my_tool
        deform_mode = mytool.deformation
        layout = self.layout
        if deform_mode == 'Roll':
            row = layout.row()
            text = "卷曲角度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Rolling Angle:"
            row.label(text = text)
            row.scale_x = 2.2
            row.prop(self, 'roll_angle', slider = True)
            row = layout.row()
            text = "卷曲轴向:" if 'zh_HAN' in bpy.context.preferences.view.language else "Rolling Axis:"
            row.label(text = text)
            row.scale_x = 1.0
            row.prop(self, 'align_axis')
            row.prop(self, 'toggle1')
            row = layout.row(align=True)
            text = "变换:" if 'zh_HAN' in bpy.context.preferences.view.language else "Transform:"
            row.label(text = text)
            row.scale_x = 0.685
            row.prop(self, 'transform1_x')
            row.prop(self, 'transform1_y')
            row.prop(self, 'transform1_z')
            row = layout.row()
            row.label(text= '')
            row.scale_x = 2.3
            row.prop(self, 'clockwise')
        elif deform_mode == 'Wave':
            row = layout.row()
            text = "波长:" if 'zh_HAN' in bpy.context.preferences.view.language else "Wavelength:"
            row.label(text = text)
            row.scale_x = 2.2
            row.prop(self, 'wavelength')
            row = layout.row()
            text = "幅度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Amplitude:"
            row.label(text = text)
            row.scale_x = 2.2
            row.prop(self, 'amplitude')
            row = layout.row()
            text = "相位:" if 'zh_HAN' in bpy.context.preferences.view.language else "Phase Angle:"
            row.label(text = text)
            row.scale_x = 2.2
            row.prop(self, 'phase_angle')
            row = layout.row()
            text = "波动朝向:" if 'zh_HAN' in bpy.context.preferences.view.language else "Waving Axis:"
            row.label(text = text)
            row.scale_x = 1.0
            row.prop(self, 'align_axis')
            row.prop(self, 'toggle2')
        elif deform_mode == 'Noise':
            row = layout.row()
            text = "噪波比例:" if 'zh_HAN' in bpy.context.preferences.view.language else "Noise Scale:"
            row.label(text = text)
            row.scale_x = 2.2
            row.prop(self, 'noise_scale')
            row = layout.row()
            text = "幅度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Amplitude:"
            row.label(text = text)
            row.scale_x = 2.2
            row.prop(self, 'amplitude')
            row = layout.row(align=True)
            text = "变换:" if 'zh_HAN' in bpy.context.preferences.view.language else "Transform:"
            row.label(text = text)
            row.scale_x = 0.685
            row.prop(self, 'transform3_x')
            row.prop(self, 'transform3_y')
            row.prop(self, 'transform3_z')
            row = layout.row()
            text = "起伏方向:" if 'zh_HAN' in bpy.context.preferences.view.language else "Noise Direction:"
            row.label(text = text)
            row.scale_x = 1.0
            row.prop(self, 'noise_dir')

    def execute(self, context):
        mytool = context.scene.my_tool
        deform_mode = mytool.deformation
        transform1 = (self.transform1_x, self.transform1_y, self.transform1_z)
        transform2 = (self.transform2_x, self.transform2_y, self.transform2_z)
        transform3 = (self.transform3_x, self.transform3_y, self.transform3_z)
        if bpy.context.object:
            ao = bpy.context.object    # ao means active object
            MyDeformation = _node.CrystalDeformation(ao, deform_mode)
            GN_deform = MyDeformation.updateDeform()
            if deform_mode == 'Roll':
                MyDeformation.roll_deformation(GN_deform, self.roll_angle, self.clockwise, self.align_axis, transform1, self.toggle1)
            elif deform_mode == 'Wave':
                MyDeformation.wave_deformation(GN_deform, self.wavelength, self.amplitude, self.phase_angle, self.align_axis, transform2, self.toggle2)
            elif deform_mode == 'Noise':
                MyDeformation.noise_deformation(GN_deform, self.noise_scale, self.amplitude, self.noise_dir, transform3)
            
            return {'FINISHED'}

class RemoveDeformButton(Operator):
    bl_idname = "mol3d.deformation_remove"
    bl_label = "Remove Defomation"
    bl_options = {'REGISTER','UNDO'}
        
    def execute(self, context):
        # clear previous data-blocks
        mytool  = context.scene.my_tool
        deform_mode = mytool.deformation
        for node_group in bpy.data.node_groups:
            if deform_mode in node_group.name:
                bpy.data.node_groups.remove(node_group)
        for obj in bpy.data.objects:
            for modifier in obj.modifiers:
                if deform_mode in modifier.name:
                    obj.modifiers.remove(modifier)
        return {'FINISHED'}


# #########################################################################################################
# #########################################################################################################
# Custom Drawing
class SetAtomsButton(Operator):
    bl_idname = 'mol3d.set_atoms'
    bl_label = "原子设置" if 'zh_HAN' in bpy.context.preferences.view.language else "Atom Setting"
    bl_options = {'REGISTER','UNDO'}

    element: StringProperty(name = '', default = 'C') # type: ignore
    rotate_angle: FloatProperty(name = '', default = 0.0, min=-180, max=180) # type: ignore
    H_rot_angle: FloatProperty(name = '', default = 90.0, min=-180, max=180) # type: ignore
    deflection_angle: FloatProperty(name = '', default = 0.0, min=-180, max=180) # type: ignore
    functional_group: EnumProperty(
        name = '',
        default = 'F0',
        items = [('F0', '元素', ''),
                 ('F1', '苯基', ''),
                 ('F2', '羧基', ''),
                 ('F3', '六元环', ''),
                 ('F4', '脂链', ''),]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('F0', 'Element', ''),
                 ('F1', 'Benzyl', ''),
                 ('F2', 'Carboxyl', ''),
                 ('F3', 'Hexatomic Ring', ''),
                 ('F4', 'n-Chain', ''),]
    ) # type:ignore
    n_carbon: IntProperty(name = '', default = 4, min=2) # type: ignore
    count: IntProperty(name = '', default = 0, min=0, max=1) # type: ignore
    text = "补齐氢原子" if 'zh_HAN' in bpy.context.preferences.view.language else "With Hydrogen"
    with_H: BoolProperty(name=text, default = False) # type: ignore
    text = "默认角度" if 'zh_HAN' in bpy.context.preferences.view.language else "Default Angle"
    default_angle: BoolProperty(name=text, default = True) # type: ignore

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        text = "官能团:" if 'zh_HAN' in bpy.context.preferences.view.language else "Functional Group:"
        row.label(text=text)
        row.scale_x = 1.8
        row.prop(self, 'functional_group')
        if self.functional_group == 'F0':
            row = layout.row(align=True)
            text = "元素名称:" if 'zh_HAN' in bpy.context.preferences.view.language else "Element Name:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'element')
        else:
            row = layout.row(align=True)
            text = "旋转角度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Rotate Angle:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'rotate_angle')
            if self.functional_group == 'F4':
                row = layout.row(align=True)
                text = "碳原子数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Carbon Num:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'n_carbon')
            
            row = layout.row()
            row.prop(self,'with_H')
            if self.with_H:
                row.scale_x = 2.0
                row.prop(self,'default_angle')
                if not self.default_angle:
                    row = layout.row(align=True)
                    text = "氢原子偏转:" if 'zh_HAN' in bpy.context.preferences.view.language else "Deflection Angle:"
                    row.label(text=text)
                    row.scale_x = 1.8
                    row.prop(self, 'deflection_angle')
                    row = layout.row(align=True)
                    text = "氢原子旋转:" if 'zh_HAN' in bpy.context.preferences.view.language else "Hydrogen Rotate:"
                    row.label(text=text)
                    row.scale_x = 1.8
                    row.prop(self, 'H_rot_angle')
    
    def draw_struct(self, obj, fg_name):
        _mesh.set_sel_atoms_attr(obj, 'C')
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        vert_counts = len(bm.verts)
        vert_ids = []
        for vert in bm.verts:
            if vert.select:
                vert_ids.append(vert.index)
                bank, pitch, heading = _mesh.vert_hpb(vert,0)
                coords = draw.functional_group_coords(vert, bank, pitch, fg_name, self.rotate_angle, self.n_carbon)
                draw.add_functional_groups(obj, coords, fg_name, self.n_carbon)
        if self.with_H:
            bpy.ops.object.mode_set(mode='EDIT')
            bm = bmesh.from_edit_mesh(obj.data)
            for vert in bm.verts:
                if vert.index >= vert_counts-1 or vert.index in vert_ids:
                    vert.select_set(True)
            lists = draw.para_branches(obj, 'H', bond_order=1, branches=1, calc_SA=True, default_angle = self.default_angle,
                                       rotate_angle=self.H_rot_angle, deflection_angle=self.deflection_angle)
            origin_atoms, origin_coords, new_atoms, new_coords, bond_orders = lists
            coll = obj.users_collection[0]
            draw.add_atoms(obj, coll, origin_atoms, origin_coords, new_atoms, new_coords, bond_orders)

    def execute(self, context):
        ao = bpy.context.object
        element = self.element.capitalize()
        if ao:
            try:
                type = ao['Type']
            except:
                type = None
            if type in ['small_mol','crystal','biomacro']:
                if self.functional_group == 'F0':
                    _mesh.set_sel_atoms_attr(ao, element)
                    _shade.create_mat_hetero('Hetero')
                    _shade.create_material(element)
                    ao['Elements'] = list(set(ao['Elements']+[element]))
                elif self.functional_group == 'F1':
                    self.draw_struct(ao, 'Benzyl')
                elif self.functional_group == 'F2':
                    self.draw_struct(ao, 'Carboxyl')
                elif self.functional_group == 'F3':
                    self.draw_struct(ao, 'Hexatomic Ring')
                elif self.functional_group == 'F4':
                    self.draw_struct(ao, 'n-Chain')
                
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}
        else:
            return {'CANCELLED'}


class SetBondsButton(Operator):
    bl_idname = 'mol3d.set_bonds'
    bl_label = "键设置" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Setting"
    bl_options = {'REGISTER','UNDO'}

    text = "键级" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Order"
    bond_order: IntProperty(name =text, default = -1, min=-1, max=3) # type: ignore
    text = "芳香环" if 'zh_HAN' in bpy.context.preferences.view.language else "Aromatic Ring"
    is_aromatic: BoolProperty(name=text, default=False) # type: ignore

    def execute(self, context):
        ao = bpy.context.object
        if ao:
            _mesh.set_sel_bonds_attr(ao, self.bond_order)
            _mesh.set_sel_bonds_aromatic(ao, self.is_aromatic)
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}
        else:
            return {'CANCELLED'}


class AddCellButton(Operator):
    bl_idname = 'mol3d.add_unit_cell'
    bl_label = "晶胞参数设置" if 'zh_HAN' in bpy.context.preferences.view.language else "Set Unit Cell Parameters"
    bl_options = {'REGISTER','UNDO'}

    length_a: FloatProperty(name ='', default = 5.0, min=0.0) # type: ignore
    length_b: FloatProperty(name ='', default = 5.0, min=0.0) # type: ignore
    length_c: FloatProperty(name ='', default = 5.0, min=0.0) # type: ignore
    angle_alpha: FloatProperty(name ='', default = 90, min=0.0, max=180) # type: ignore
    angle_beta: FloatProperty(name ='', default = 90, min=0.0, max=180) # type: ignore
    angle_gamma: FloatProperty(name ='', default = 90, min=0.0, max=180) # type: ignore
    space_group: StringProperty(name = '', default = 'P1') # type: ignore
    crys_name: StringProperty(name = '', default = '') # type: ignore

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        text = "晶胞轴长:" if 'zh_HAN' in bpy.context.preferences.view.language else "Cell Lengths:"
        row.label(text = text)
        row.scale_x = 0.61
        row.prop(self, 'length_a')
        row.prop(self, 'length_b')
        row.prop(self, 'length_c')
        row = layout.row(align=True)
        text = "晶胞轴角:" if 'zh_HAN' in bpy.context.preferences.view.language else "Cell Angles:"
        row.label(text = text)
        row.scale_x = 0.61
        row.prop(self, 'angle_alpha')
        row.prop(self, 'angle_beta')
        row.prop(self, 'angle_gamma')
        row = layout.row(align=True)
        text = "空间群:" if 'zh_HAN' in bpy.context.preferences.view.language else "Space Group:"
        row.label(text = text)
        row.scale_x = 1.8
        row.prop(self, 'space_group')
        row = layout.row(align=True)
        text = "晶体名称:" if 'zh_HAN' in bpy.context.preferences.view.language else "Crystal Name:"
        row.label(text = text)
        row.scale_x = 1.8
        row.prop(self, 'crys_name')

    def execute(self, context):
        coll_mol = bpy.data.collections.new('Struct_'+self.crys_name)
        context.scene.collection.children.link(coll_mol)
        space_group = self.space_group
        draw.draw_cell_edges(self.crys_name, coll_mol, self.length_a, self.length_b, self.length_c, self.angle_alpha, self.angle_beta, self.angle_gamma, space_group)
        return {'FINISHED'}

        
class AddStructButton(Operator):
    bl_idname = 'mol3d.add_struct'
    bl_label = "添加新的结构" if 'zh_HAN' in bpy.context.preferences.view.language else "Add New Structure"
    bl_options = {'REGISTER','UNDO'}

    new_atom_1: StringProperty(name ='', default = 'H') # type: ignore
    new_atom_2: StringProperty(name ='', default = '') # type: ignore
    new_atom_3: StringProperty(name ='', default = '') # type: ignore
    new_atom_4: StringProperty(name ='', default = '') # type: ignore
    new_atom_5: StringProperty(name ='', default = '') # type: ignore
    bond_order: IntProperty(name ='', default = 1, min=0, max=3) # type: ignore
    rotate_angle: FloatProperty(name = '', default = 0.0, min=-180, max=180) # type: ignore
    deflection_angle: FloatProperty(name = '', default = 0.0, min=-180, max=180) # type: ignore
    branches: IntProperty(name ='', default = 1, min=1, max=5) # type: ignore
    text = "计算饱和度" if 'zh_HAN' in bpy.context.preferences.view.language else "Calc. Saturation"
    calc_SA: BoolProperty(name=text, default = True) # type: ignore
    text = "默认角度" if 'zh_HAN' in bpy.context.preferences.view.language else "Default Angle"
    default_angle: BoolProperty(name=text, default = True) # type: ignore
    fract_x1: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_y1: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_z1: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_x2: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_y2: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_z2: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_x3: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_y3: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_z3: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_x4: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_y4: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_z4: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_x5: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_y5: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    fract_z5: FloatProperty(name ='', default = 0.0, min=0.0, max=1.0) # type: ignore
    grow_by_fract: FloatProperty(name='', default=0.0, min=0.0, max=1.0) #type: ignore
    grow_direct: IntProperty(name='', default=0, min=0, soft_max=6) #type: ignore
    bond_length_f: FloatProperty(name='', default=1.0, min=0.0, soft_max=2.0) #type: ignore
    boundary: FloatProperty(name='', default=0.0, min=0.0, soft_max=0.5) #type: ignore

    def condition(self):
        ao = bpy.context.object
        try:
            mode = bpy.context.object.mode
        except:
            mode = None
        try:
            type = ao['Type']
        except:
            type = None
        return (mode, type)

    def draw(self, context):
        mode, type = self.condition()
        layout = self.layout
        if mode == 'EDIT' or type in ['small_mol','biomacro']:
            row = layout.row()
            text = "新的原子:" if 'zh_HAN' in bpy.context.preferences.view.language else "New Atom:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'new_atom_1')
            row = layout.row()
            text = "键级:" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Order:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bond_order')

            row = layout.row()
            text = "数量:" if 'zh_HAN' in bpy.context.preferences.view.language else "Branches:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'branches')
            row = layout.row()
            row.prop(self, 'calc_SA')
            row.scale_x = 1.8
            row.prop(self, 'default_angle')

            if not self.default_angle:
                row = layout.row()
                text = "旋转角度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Rotate Angle:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'rotate_angle')
                row = layout.row()
                text = "偏转角度:" if 'zh_HAN' in bpy.context.preferences.view.language else "Deflection Angle:"
                row.label(text=text)
                row.scale_x = 1.8
                row.prop(self, 'deflection_angle')

        elif mode == 'OBJECT' and type in ['crystal', 'unit_cell']:
            row = layout.row()
            col = row.column(align=True)
            text = "新的原子:" if 'zh_HAN' in bpy.context.preferences.view.language else "New Atom:"
            col.label(text=text)
            col.scale_x = 0.4
            col = row.column(align=True)
            text = "位置分数:" if 'zh_HAN' in bpy.context.preferences.view.language else "Cell Fracts:"
            col.label(text = text)
            
            row = layout.row()
            col = row.column(align=True)
            col.prop(self, 'new_atom_1')
            col.scale_x = 1.2
            col = row.column(align=True)
            row = col.row(align=True)
            row.prop(self, 'fract_x1')
            row.prop(self, 'fract_y1')
            row.prop(self, 'fract_z1')

            row = layout.row()
            col = row.column(align=True)
            col.prop(self, 'new_atom_2')
            col.scale_x = 1.2
            col = row.column(align=True)
            row = col.row(align=True)
            row.prop(self, 'fract_x2')
            row.prop(self, 'fract_y2')
            row.prop(self, 'fract_z2')

            row = layout.row()
            col = row.column(align=True)
            col.prop(self, 'new_atom_3')
            col.scale_x = 1.2
            col = row.column(align=True)
            row = col.row(align=True)
            row.prop(self, 'fract_x3')
            row.prop(self, 'fract_y3')
            row.prop(self, 'fract_z3')
            
            row = layout.row()
            col = row.column(align=True)
            col.prop(self, 'new_atom_4')
            col.scale_x = 1.2
            col = row.column(align=True)
            row = col.row(align=True)
            row.prop(self, 'fract_x4')
            row.prop(self, 'fract_y4')
            row.prop(self, 'fract_z4')
            '''
            row = layout.row()
            col = row.column(align=True)
            col.prop(self, 'new_atom_5')
            col.scale_x = 1.2
            col = row.column(align=True)
            row = col.row(align=True)
            row.prop(self, 'fract_x5')
            row.prop(self, 'fract_y5')
            row.prop(self, 'fract_z5')
            '''

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
            text = "键长倍增系数" if 'zh_HAN' in bpy.context.preferences.view.language else "Bond Length Factor:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'bond_length_f')
            row = layout.row()
            text = "边界连键:" if 'zh_HAN' in bpy.context.preferences.view.language else "Boundary Bonding:"
            row.label(text=text)
            row.scale_x = 1.8
            row.prop(self, 'boundary')


    def execute(self, context):
        mode, type = self.condition()
        if type not in ['unit_cell','crystal'] and mode != 'EDIT':
            new_molecule = draw.add_CC()
            GN_mol = _node.create_GeoNode(new_molecule, 'Struct_'+new_molecule.name, 'GN_'+new_molecule.name)
            group = GN_mol.node_group
            _node.Ball_Stick_nodes(group, atom_scale=1.0, bond_scale=1.0, atom_subdiv=3, bond_subdiv=12, atom_syms = ['C'], type='B', hide_H=False, 
                                   color_type='H', dash_count=5, dash_radius=0.1, category='C2', display_mode='M2', hide_segs='')
            new_molecule['Type'] = 'small_mol'
            new_molecule['Elements'] = ['C']
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}
        else:
            ao = context.object
            coll = ao.users_collection[0]
            new_atom_1 = self.new_atom_1.capitalize()
            new_atom_2 = self.new_atom_2.capitalize()
            new_atom_3 = self.new_atom_3.capitalize()
            new_atom_4 = self.new_atom_4.capitalize()
            new_atom_5 = self.new_atom_5.capitalize()
            if mode == 'EDIT':
                bpy.ops.object.mode_set(mode='EDIT')
                try:
                    lists = draw.para_branches(ao, new_atom_1, self.bond_order, self.branches, self.calc_SA, 
                                               self.default_angle, self.rotate_angle, self.deflection_angle)
                    origin_atoms, origin_coords, new_atoms, new_coords, bond_orders = lists
                    draw.add_atoms(ao, coll, origin_atoms, origin_coords, new_atoms, new_coords, bond_orders)
                    _shade.create_mat_hetero('Hetero')
                    _shade.create_material(new_atom_1)
                except Exception as e:
                    print(e)
                ao['Elements'] = list(set(ao['Elements']+[new_atom_1]))
                bpy.ops.object.mode_set(mode='EDIT')
            elif mode == 'OBJECT' and type == 'unit_cell':
                new_crystal_name = 'Crystal_'+ao.name[11:]
                length_a, length_b, length_c = eval(ao['cell lengths'])
                angle_alpha, angle_beta, angle_gamma = eval(ao['cell angles'])
                vec_a = _math.fract_to_cartn((1,0,0), length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                vec_b = _math.fract_to_cartn((0,1,0), length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                vec_c = _math.fract_to_cartn((0,0,1), length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)
                cell_cycles = (1, 1, 1)
                edge_vecs = (vec_a, vec_b, vec_c)
                filters = ''
                space_group = ao['space group']
                space_group = re.sub(" ", "", space_group).capitalize()
                symop_operations = SPACEGROUP_DEFAULTS[space_group][3]

                fract_xyzs = [(self.fract_x1, self.fract_y1, self.fract_z1),
                                (self.fract_x2, self.fract_y2, self.fract_z2),
                                (self.fract_x3, self.fract_y3, self.fract_z3),
                                (self.fract_x4, self.fract_y4, self.fract_z4),
                                (self.fract_x5, self.fract_y5, self.fract_z5),]
                new_atoms = [new_atom_1, new_atom_2, new_atom_3, new_atom_4, new_atom_5]
                Atoms, Coords = [],[]
                for new_atom,fract_xyz in zip(new_atoms,fract_xyzs):
                    sym_fracts = _math.fract_symop(fract_xyz, symop_operations)
                    sym_fracts = _math.fracts_normalize(sym_fracts, self.boundary)
                    atoms = [new_atom]*len(sym_fracts)
                    coords = [_math.fract_to_cartn(fract, length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma) for fract in sym_fracts]
                    for atom, coord in zip(atoms,coords):
                        if atom:
                            Atoms.append(atom)
                            Coords.append(coord)
                try:
                    Bonds, Bond_Orders = add_BONDS(Atoms, Coords, self.bond_length_f)
                    molecule = [Atoms, Coords, Bonds, Bond_Orders]
                    # create new crystal scaffold
                    new_crystal = mol_scaffold_build(new_crystal_name, 'crystal', coll, molecule)
                    atom_syms = new_crystal['Elements']
                    # add geometry node group
                    _node.expand_scaffold(new_crystal, ao, edge_vecs, self.grow_direct, self.grow_by_fract)
                    _node.unit_cycles(new_crystal, new_crystal_name, 'Crys_Cycles_'+ao.name[11:], cell_cycles, edge_vecs)
                    _node.crys_filter(new_crystal, new_crystal_name, 'Crys_Filter_'+new_crystal_name, filters)
                    GN_mol = _node.create_GeoNode(new_crystal, 'Struct_'+new_crystal_name, 'GN_'+new_crystal_name)
                    _node.Ball_Stick_nodes(GN_mol.node_group, atom_scale=1.0, bond_scale=1.0, atom_subdiv=3, bond_subdiv=12, atom_syms=atom_syms, type='B', hide_H=False,
                                        color_type='H', dash_count=5, dash_radius=0.1, category='C2', display_mode='M2', hide_segs='')
                except Exception as e:
                    print(e)
                
            return {'FINISHED'}
        
class MMFFOptimizeButton(Operator):
    bl_idname = 'mol3d.mmff_optimize'
    bl_label = "Conformation Optimize"
    bl_options = {'REGISTER','UNDO'}

    def execute(self, context):
        ao = context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'small_mol':
                _mesh.select_verts(ao,'H')
                bpy.ops.mesh.delete(type='VERT')
                bpy.ops.object.mode_set(mode='OBJECT')
                new_scaffold = _mesh.struct_optimize(ao, addHs=False)
                coll = new_scaffold.users_collection[0]
                
                bpy.ops.object.mode_set(mode='EDIT')
                try:
                    lists = draw.para_branches(new_scaffold, 'H', bond_order=1, branches=1, calc_SA=True, 
                                               default_angle=True, rotate_angle=0, deflection_angle=0)
                    origin_atoms, origin_coords, new_atoms, new_coords, bond_orders = lists
                    draw.add_atoms(new_scaffold, coll, origin_atoms, origin_coords, new_atoms, new_coords, bond_orders)
                    _shade.create_mat_hetero('Hetero')
                    _shade.create_material('H')
                except Exception as e:
                    print(e)
                new_scaffold['Elements'] = list(set(new_scaffold['Elements']+['H']))
                bpy.ops.object.mode_set(mode='OBJECT')
                
                
                new_atom_orientations = _mesh.atom_orient(new_scaffold)
                print(new_atom_orientations)
                new_bond_vertic_vecs = _mesh.bond_vertical_dir(new_scaffold)
                _mesh.add_attribute(new_scaffold, 'atom_orient', 'FLOAT_VECTOR', 'POINT', new_atom_orientations)
                _mesh.add_attribute(new_scaffold, 'bond_vertvec', 'FLOAT_VECTOR', 'EDGE', new_bond_vertic_vecs)
                _mesh.is_aromatic(new_scaffold)

                
                filename = new_scaffold.name
                GN_mol = _node.create_GeoNode(new_scaffold, 'Struct_'+filename, 'GN_'+filename)  
                atom_syms = new_scaffold['Elements']
                _node.Ball_Stick_nodes(GN_mol.node_group, atom_scale=1, bond_scale=1, atom_subdiv=3, bond_subdiv=12,
                                       atom_syms=atom_syms, type='B', hide_H=False, color_type='H', dash_count=5, 
                                       dash_radius=0.1, category='C2', display_mode='M2', hide_segs='')
                
            return {'FINISHED'}
        
class MMFFUpdateButton(Operator):
    bl_idname = 'mol3d.mmff_update'
    bl_label = "Conformation Update"
    bl_options = {'REGISTER','UNDO'}

    def execute(self, context):
        ao = context.object
        if not ao:
            return {'CANCELLED'}
        else:
            try:
                type = ao['Type']
            except:
                type = None
            if type == 'small_mol':
                _mesh.select_verts(ao,'H')
                bpy.ops.mesh.delete(type='VERT')
                bpy.ops.object.mode_set(mode='OBJECT')
                new_scaffold = _mesh.struct_optimize(ao, addHs=True)
                new_atom_orientations = _mesh.atom_orient(new_scaffold)
                new_bond_vertic_vecs = _mesh.bond_vertical_dir(new_scaffold)
                _mesh.add_attribute(new_scaffold, 'atom_orient', 'FLOAT_VECTOR', 'POINT', new_atom_orientations)
                _mesh.add_attribute(new_scaffold, 'bond_vertvec', 'FLOAT_VECTOR', 'EDGE', new_bond_vertic_vecs)
                _mesh.is_aromatic(new_scaffold)

                filename = new_scaffold.name
                GN_mol = _node.create_GeoNode(new_scaffold, 'Struct_'+filename, 'GN_'+filename)  
                atom_syms = new_scaffold['Elements']
                _node.Ball_Stick_nodes(GN_mol.node_group, atom_scale=1, bond_scale=1, atom_subdiv=3, bond_subdiv=12,
                                       atom_syms=atom_syms, type='B', hide_H=False, color_type='H', dash_count=5, 
                                       dash_radius=0.1, category='C2', display_mode='M2', hide_segs='')
            return {'FINISHED'}