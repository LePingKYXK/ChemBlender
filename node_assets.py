import bpy
from bl_ui import node_add_menu
from .Chem_Data import ELEMENTS_DEFAULT
from . import _node
from . import version as v
from ._node import (add_node, get_node, add_attr_node, add_math_node, add_vec_math_node, add_compare_node, add_mat_node,
                    add_sample_index_node, add_store_attr_node, delete_mismatch, mesh_to_curve, sepa_by_attr, 
                    set_node_values, set_node_operation, set_node_datatype, nodes_link)

class MOL3D_MT_STRUCT(bpy.types.Menu):
    bl_idname = 'MOL3D_MT_STRUCT'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        # layout.operator_context = "INVOKE_DEFAULT"
        MOL3D_MT_NODE_names = ["Mol3D_MT_Atoms",
                               "Mol3D_MT_Bonds", 
                               "Mol3D_MT_Sepa_Bonds",
                               "Mol3D_MT_Dashline",
                               "Mol3D_MT_Aromatic_Bonds",
                               "Mol3D_MT_Nucleic_Restruct",
                               "Mol3D_MT_Peptide_Restruct",
                               "Mol3D_MT_Alpha_Helix",
                               "Mol3D_MT_Beta_Sheet",
                               "Mol3D_MT_Loop_Coil",
                               "Mol3D_MT_Nucleic",
                               ]
        for node_name in MOL3D_MT_NODE_names: node_add_menu.add_node_type(layout, node_name)
      
class MOL3D_MT_ATTRIBUTE(bpy.types.Menu):
    bl_idname = 'MOL3D_MT_ATTRIBUTE'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        # layout.operator_context = "INVOKE_DEFAULT"
        MOL3D_MT_NODE_names = ["Mol3D_MT_BioAttr_Preset",
                               "Mol3D_MT_BioAttr_Reset",]
        for node_name in MOL3D_MT_NODE_names: node_add_menu.add_node_type(layout, node_name)


class MOL3D_MT_NODES(bpy.types.Menu):
    bl_idname = "MOL3D_MT_NODES"
    bl_label = "Mol3D Nodes"

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"
        layout.menu('MOL3D_MT_STRUCT', text='Structure')
        layout.menu('MOL3D_MT_ATTRIBUTE', text='Attribute')


def Mol3D_node_menu(self, context):
    layout = self.layout
    layout.menu("MOL3D_MT_NODES", text='Mol3D Nodes')


# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
class Mol3D_MT_Atoms(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_ATOMS'
    bl_label='Atoms'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Scaffold', 'INPUT', 'NodeSocketGeometry')   # inputs
        self.add_new_socket('Atom Radii', 'INPUT', 'NodeSocketString')
        self.add_new_socket('Atom Scale', 'INPUT', 'NodeSocketFloat')
        self.add_new_socket('Atom Subdiv', 'INPUT', 'NodeSocketInt')
        self.add_new_socket('Center Scale', 'INPUT', 'NodeSocketFloat')
        self.add_new_socket('Ligand Scale', 'INPUT', 'NodeSocketFloat')
        self.add_new_socket('Atoms', 'OUTPUT', 'NodeSocketGeometry')   # outputs
        # defalut values
        self.inputs[1].default_value = 'radii_B'
        self.inputs[2].default_value = 1.0
        self.inputs[3].default_value = 3
        self.inputs[4].default_value = 1.0
        self.inputs[5].default_value = 1.0

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a-600,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+600,b), '')
        _instance = add_node(group, 'GeometryNodeInstanceOnPoints', (a+200,b), '')
        _realize = add_node(group, 'GeometryNodeRealizeInstances', (a+400,b), '')
        _uv_sphere = add_node(group, 'GeometryNodeMeshUVSphere', (a,b-200), '')
        _align_euler2vec = add_node(group, 'FunctionNodeAlignEulerToVector',(a,b-400), '')
        _align_euler2vec.axis = 'Z'
        _attr_atom_orient = add_attr_node(group, (a-200,b-500), '', [(0,'atom_orient')], 'FLOAT_VECTOR')
        _attr_radii = add_attr_node(group, (a-200,b-700), '', [], 'FLOAT')
        _multiply = add_math_node(group, (a,b-600), '', [], 'MULTIPLY')
        _multiply_segs = add_math_node(group, (a-200,b), '', [(1,8)], 'MULTIPLY')
        _multiply_rings = add_math_node(group, (a-200, b-200), '', [(1,5)], 'MULTIPLY')

        _attr_poly_aid = add_attr_node(group, (a-400,b-900), '', [(0,'poly_aid')], 'FLOAT')
        _not_equal = add_compare_node(group, (a-200,b-900),'', [(1,0.0)], 'FLOAT', 'NOT_EQUAL')
        _equal = add_compare_node(group, (a-200,b-1100),'', [(1,1.0)], 'FLOAT', 'EQUAL')
        _switch_1 = add_node(group, 'GeometryNodeSwitch', (a,b-800), '')
        _switch_1.input_type = 'FLOAT'
        set_node_values(_switch_1,[(v.sw_sockets_in['floatA'],1.0)])
        _switch_2 = add_node(group, 'GeometryNodeSwitch', (a,b-1000), '')
        _switch_2.input_type = 'FLOAT'
        _multiply_2 = add_math_node(group, (a+200,b-600), '', [], 'MULTIPLY')
        nodes_link(group, _attr_poly_aid, v.na_sockets_out['float'], _not_equal, 0)
        nodes_link(group, _attr_poly_aid, v.na_sockets_out['float'], _equal, 0)
        nodes_link(group, _not_equal, 0, _switch_1, 0)
        nodes_link(group, _equal, 0, _switch_2, 0)
        nodes_link(group, _switch_2, 0, _switch_1, v.sw_sockets_in['floatB'])
        nodes_link(group, _switch_1, 0, _multiply_2, 1)

        nodes_link(group, _input, 0, _instance, 0)
        nodes_link(group, _input, 1, _attr_radii, 0)
        nodes_link(group, _input, 2, _multiply, 0)
        nodes_link(group, _input, 3, _multiply_rings, 0)
        nodes_link(group, _input, 3, _multiply_segs, 0)
        nodes_link(group, _input, 4, _switch_2, v.sw_sockets_in['floatB'])
        nodes_link(group, _input, 5, _switch_2, v.sw_sockets_in['floatA'])
        nodes_link(group, _multiply_segs, 0, _uv_sphere, 0)
        nodes_link(group, _multiply_rings, 0, _uv_sphere, 1)
        nodes_link(group, _instance, 0, _realize, 0)
        nodes_link(group, _realize, 0, _output, 0)
        nodes_link(group, _uv_sphere, 0, _instance, 2)
        nodes_link(group, _attr_atom_orient, 0, _align_euler2vec, 2)
        nodes_link(group, _align_euler2vec, 0, _instance, 5)
        nodes_link(group, _attr_radii, v.na_sockets_out['float'], _multiply, 1)
        nodes_link(group, _multiply, 0, _multiply_2, 0)
        nodes_link(group, _multiply_2, 0, _instance, 6)


class Mol3D_MT_Bonds(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_BONDS'
    bl_label='Bonds'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Scaffold', 'INPUT', 'NodeSocketGeometry')   # inputs
        self.add_new_socket('Bond Radii', 'INPUT', 'NodeSocketString') # socket 1
        self.add_new_socket('Bond Scale', 'INPUT', 'NodeSocketFloat') # socket 2
        self.add_new_socket('Bond Subdiv', 'INPUT', 'NodeSocketInt') # socket 3
        self.add_new_socket('Polyhedra Skeleton', 'INPUT', 'NodeSocketBool') # socket 4
        self.add_new_socket('Show Aromatic', 'INPUT', 'NodeSocketBool') # socket 5
        self.add_new_socket('Bonds', 'OUTPUT', 'NodeSocketGeometry')   # outputs
        # defalut values
        self.inputs[1].default_value = 'bond_radii_B'
        self.inputs[2].default_value = 1.0
        self.inputs[3].default_value = 12
        self.inputs[4].default_value = True
        self.inputs[5].default_value = False

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+3200,b), '')
        # split edges and hide polyhedra skeletons
        _split_edges = add_node(group, 'GeometryNodeSplitEdges', (a+200,b), '')
        _sepa_bonds = add_node(group, 'GeometryNodeSeparateGeometry', (a+800,b), '')
        _sepa_bonds.domain = 'EDGE'
        _attr_bond_order = add_attr_node(group, (a+600,b-200), '', [(0,'bond_order')], 'INT')
        _greater = add_compare_node(group, (a+800,b-200),'', [(3,0)], 'INT', 'GREATER_THAN')
        nodes_link(group, _split_edges, 0, _sepa_bonds, 0)
        nodes_link(group, _attr_bond_order, v.na_sockets_out['int'], _greater, 2)
        nodes_link(group, _greater, 0, _sepa_bonds, 1)

        # preset attributes for bonds coloring
        _store_attr_dir, _store_attr_length = preset_attr_bonds_color(group, (a+200,b), _split_edges)
        
        _subdivide = add_node(group, 'GeometryNodeSubdivideMesh', (a+1600,b+200), '')
        _delete_geo = add_node(group, 'GeometryNodeDeleteGeometry', (a+600,b), '')
        _delete_geo.domain = 'EDGE'
        _attr_poly_aid = add_attr_node(group, (a-200, b+200), '', [(0,'poly_aid')], 'FLOAT')
        _equal_1 = add_compare_node(group, (a,b+600), '', [(1,1.25)], 'FLOAT', 'EQUAL')  # 1.25 and 1.75 depends on attribute 'poly_aid' values
        _equal_2 = add_compare_node(group, (a,b+400), '', [(1,1.75)], 'FLOAT', 'EQUAL')
        _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a+200,b+200), 'Bonds in Poly')
        _boolean_AND.operation = 'AND'
        _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (a+200,b+400), 'Bonds in Poly')
        _boolean_OR.operation = 'OR'
        _switch = add_node(group, 'GeometryNodeSwitch', (a+400,b), '')
        _switch.input_type = 'BOOLEAN'
        nodes_link(group, _attr_poly_aid, v.na_sockets_out['float'], _equal_1, 0)
        nodes_link(group, _attr_poly_aid, v.na_sockets_out['float'], _equal_2, 0)
        nodes_link(group, _equal_1, 0, _boolean_AND, 0)
        nodes_link(group, _equal_2, 0, _boolean_AND, 1)
        nodes_link(group, _equal_1, 0, _boolean_OR, 0)
        nodes_link(group, _equal_2, 0, _boolean_OR, 1)
        nodes_link(group, _boolean_OR, 0, _switch, v.sw_sockets_in['boolA'])
        nodes_link(group, _boolean_AND, 0, _switch, v.sw_sockets_in['boolB'])
        nodes_link(group, _switch, v.sw_sockets_out['bool'], _delete_geo, 1)

        _mesh2curve = add_node(group, 'GeometryNodeMeshToCurve', (a+1600,b), '')
        _set_crv_radius = add_node(group, 'GeometryNodeSetCurveRadius', (a+1800,b), '')
        _curve2mesh = add_node(group, 'GeometryNodeCurveToMesh', (a+2000,b), '')
        set_node_values(_curve2mesh, [(2,True)])
        _curve_circle = add_node(group, 'GeometryNodeCurvePrimitiveCircle', (a+1800,b-200), 'Curve Circle')
        _attr_bond_radii = add_attr_node(group, (a+1600,b-200), '', [], 'FLOAT')
        _attr_statistic = add_node(group, 'GeometryNodeAttributeStatistic', (a+1600,b-400), '')
        _attr_is_aromatic = add_attr_node(group, (a+1600,b-800), '', [(0,'is_aromatic_edge')], 'BOOLEAN')
        _switch_1 = add_node(group, 'GeometryNodeSwitch', (a+1800,b-800), '')
        _switch_1.input_type = 'FLOAT'
        _switch_2 = add_node(group, 'GeometryNodeSwitch', (a+1800,b-600), '')
        _switch_2.input_type = 'FLOAT'

        nodes_link(group, _sepa_bonds, 0, _store_attr_dir, 0)
        nodes_link(group, _store_attr_length, 0, _subdivide, 0)
        nodes_link(group, _subdivide, 0, _delete_geo, 0)
        nodes_link(group, _delete_geo, 0, _mesh2curve, 0)
        nodes_link(group, _mesh2curve, 0, _set_crv_radius, 0)
        nodes_link(group, _set_crv_radius, 0, _curve2mesh, 0)
        nodes_link(group, _attr_bond_radii, v.na_sockets_out['float'], _switch_1, v.sw_sockets_in['floatA'])
        nodes_link(group, _attr_bond_radii, v.na_sockets_out['float'], _switch_2, v.sw_sockets_in['floatA'])
        nodes_link(group, _mesh2curve, 0, _attr_statistic, 0)
        nodes_link(group, _attr_bond_radii, v.na_sockets_out['float'], _attr_statistic, 2)
        nodes_link(group, _attr_is_aromatic, v.na_sockets_out['bool'], _switch_1, 0)
        nodes_link(group, _attr_statistic, v.na_sockets_out['int'], _switch_1, v.sw_sockets_in['floatB'])
        nodes_link(group, _switch_1, 0, _switch_2, v.sw_sockets_in['floatB'])
        nodes_link(group, _switch_2, 0, _set_crv_radius, 2)
        nodes_link(group, _curve_circle, 0, _curve2mesh, 1)
        # link parameters
        nodes_link(group, _input, 0, _split_edges, 0)
        nodes_link(group, _input, 1, _attr_bond_radii, 0)
        nodes_link(group, _input, 2, _curve_circle, 4)
        nodes_link(group, _input, 3, _curve_circle, 0)
        nodes_link(group, _input, 4, _switch, 0)
        nodes_link(group, _input, 5, _switch_2, 0)

        _sepa_segs, _join_geo = add_attr_bonds_color(group, (a+2000,b), _split_edges)
        nodes_link(group, _curve2mesh, 0, _sepa_segs, 0)
        nodes_link(group, _join_geo, 0, _output, 0)


class Mol3D_MT_Sepa_Bonds(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_SEPA_BONDS'
    bl_label='Separate Bonds'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
    
    def sepa_bonds_I(self, location):
        group = self.node_tree
        a,b = location
        _sepa_bonds_I = add_node(group, 'GeometryNodeSeparateGeometry', (a,b), '')
        _sepa_bonds_I.domain = 'FACE'
        _attr_bond_order = add_attr_node(group, (a-600,b), '', [(0,'bond_order')], 'INT')
        _attr_is_aromatic = add_attr_node(group, (a-600,b-400), '', [(0,'is_aromatic_edge')], 'BOOLEAN')
        _equal = add_compare_node(group, (a-400,b),'', [(3,1)], 'INT', 'EQUAL')
        _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (a-400,b-200), '')
        _boolean_OR.operation = 'OR'
        _switch = add_node(group, 'GeometryNodeSwitch', (a-200,b), '')
        _switch.input_type = 'BOOLEAN'
        nodes_link(group, _attr_bond_order, v.na_sockets_out['int'], _equal, 2)
        nodes_link(group, _equal, 0, _boolean_OR, 0)
        nodes_link(group, _attr_is_aromatic, v.na_sockets_out['bool'], _boolean_OR, 1)
        nodes_link(group, _equal, 0, _switch, v.sw_sockets_in['boolA'])
        nodes_link(group, _boolean_OR, 0, _switch, v.sw_sockets_in['boolB'])
        nodes_link(group, _switch, v.sw_sockets_out['bool'], _sepa_bonds_I, 1)
        return (_sepa_bonds_I, _switch)

    def sepa_bonds_II(self, location):
        group = self.node_tree
        a,b = location
        _sepa_bonds_II = add_node(group, 'GeometryNodeSeparateGeometry', (a,b), '')
        _sepa_bonds_II.domain = 'FACE'
        _attr_bond_order = add_attr_node(group, (a-600,b), '', [(0,'bond_order')], 'INT')
        _attr_is_aromatic = add_attr_node(group, (a-600,b-200), '', [(0,'is_aromatic_edge')], 'BOOLEAN')
        _equal = add_compare_node(group, (a-400,b),'', [(3,2)], 'INT', 'EQUAL')
        _boolean_NOT = add_node(group, 'FunctionNodeBooleanMath', (a-400,b-200), '')
        _boolean_NOT.operation = 'NOT'
        _switch = add_node(group, 'GeometryNodeSwitch', (a-200,b), '')
        _switch.input_type = 'BOOLEAN'
        _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a-200,b-200), '')
        _boolean_AND.operation = 'AND'
        nodes_link(group, _attr_bond_order, v.na_sockets_out['int'], _equal, 2)
        nodes_link(group, _equal, 0, _boolean_AND, 0)
        nodes_link(group, _attr_is_aromatic, v.na_sockets_out['bool'], _boolean_NOT, 0)
        nodes_link(group, _boolean_NOT, 0, _boolean_AND, 1)
        nodes_link(group, _equal, 0, _switch, v.sw_sockets_in['boolA'])
        nodes_link(group, _boolean_AND, 0, _switch, v.sw_sockets_in['boolB'])
        nodes_link(group, _switch, v.sw_sockets_out['bool'], _sepa_bonds_II, 1)
        return (_sepa_bonds_II, _switch)
    
    def sepa_bonds_III(self, location):
        group = self.node_tree
        a,b = location
        sepa_bonds_III = add_node(group, 'GeometryNodeSeparateGeometry', (a,b), '')
        sepa_bonds_III.domain = 'FACE'
        _attr_bond_order = add_attr_node(group, (a-400,b), '', [(0,'bond_order')], 'INT')
        _equal = add_compare_node(group, (a-200,b),'', [(3,3)], 'INT', 'EQUAL')
        nodes_link(group, _attr_bond_order, v.na_sockets_out['int'], _equal, 2)
        nodes_link(group, _equal, 0, sepa_bonds_III, 1)
        return sepa_bonds_III

    def bond_offset(self, location):
        group = self.node_tree
        a,b = location
        _set_pos = add_node(group, 'GeometryNodeSetPosition', (a,b), '')
        _scale = add_vec_math_node(group, (a,b-200), '', [], 'SCALE')
        nodes_link(group, _scale, 0, _set_pos, 3)
        return (_set_pos, _scale)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Bonds', 'INPUT', 'NodeSocketGeometry')   # inputs 0
        self.add_new_socket('Double Offset', 'INPUT', 'NodeSocketFloat') # socket 1
        self.add_new_socket('Triple Offset', 'INPUT', 'NodeSocketFloat') # socket 2
        self.add_new_socket('Show Aromatic', 'INPUT', 'NodeSocketBool') # socket 3
        self.add_new_socket('Bonds I', 'OUTPUT', 'NodeSocketGeometry')   # outputs 0
        self.add_new_socket('Bonds II', 'OUTPUT', 'NodeSocketGeometry')   # outputs 1
        self.add_new_socket('Bonds III', 'OUTPUT', 'NodeSocketGeometry')   # outputs 2
        # defalut values
        self.inputs[1].default_value = 0.12
        self.inputs[2].default_value = 0.2
        self.inputs[3].default_value = False

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+800,b), '')
        _attr_bond_order = add_attr_node(group, (a,b+200), '', [(0,'bond_order')], 'INT')
        _attr_bond_vertvec = add_attr_node(group, (a,b+200), '', [(0,'bond_vertvec')], 'FLOAT_VECTOR')
        _sepa_bonds_I, _switch_I = self.sepa_bonds_I((a+200,b+200))
        _sepa_bonds_II, _switch_II = self.sepa_bonds_II((a+200,b+600))
        _sepa_bonds_III = self.sepa_bonds_III((a+200,b-200))
        nodes_link(group, _input, 0, _sepa_bonds_I, 0)
        nodes_link(group, _input, 0, _sepa_bonds_II, 0)
        nodes_link(group, _input, 0, _sepa_bonds_III, 0)

        _set_pos_IIA, _scale_IIA = self.bond_offset((a+400,b+600))
        _set_pos_IIB, _scale_IIB = self.bond_offset((a+600,b+600))
        _set_pos_IIIA, _scale_IIIA = self.bond_offset((a+400,b-200))
        _set_pos_IIIB, _scale_IIIB = self.bond_offset((a+600,b-200))
        _multiply_IIB = add_math_node(group, (a+400,b+200), '', [(1,-1)], 'MULTIPLY')
        _multiply_IIIB = add_math_node(group, (a+400,b), '', [(1,-1)], 'MULTIPLY')
        nodes_link(group, _attr_bond_vertvec, 0, _scale_IIA, 0)
        nodes_link(group, _attr_bond_vertvec, 0, _scale_IIB, 0)
        nodes_link(group, _attr_bond_vertvec, 0, _scale_IIIA, 0)
        nodes_link(group, _attr_bond_vertvec, 0, _scale_IIIB, 0)
        nodes_link(group, _sepa_bonds_II, 0, _set_pos_IIA, 0)
        nodes_link(group, _sepa_bonds_II, 0, _set_pos_IIB, 0)
        nodes_link(group, _sepa_bonds_III, 0, _set_pos_IIIA, 0)
        nodes_link(group, _sepa_bonds_III, 0, _set_pos_IIIB, 0)
        nodes_link(group, _multiply_IIB, 0, _scale_IIB, 3)
        nodes_link(group, _multiply_IIIB, 0, _scale_IIIB, 3)

        nodes_link(group, _input, 1, _scale_IIA, 3)
        nodes_link(group, _input, 1, _multiply_IIB, 0)
        nodes_link(group, _input, 2, _scale_IIIA, 3)
        nodes_link(group, _input, 2, _multiply_IIIB, 0)
        nodes_link(group, _input, 3, _switch_I, 0)
        nodes_link(group, _input, 3, _switch_II, 0)

        _join_II = add_node(group, 'GeometryNodeJoinGeometry', (a+800,b+400), 'Bonds II')
        _join_III = add_node(group, 'GeometryNodeJoinGeometry', (a+800,b-200), 'Bonds III')
        nodes_link(group, _sepa_bonds_I, 0, _output, 0)
        nodes_link(group, _set_pos_IIA, 0, _join_II, 0)
        nodes_link(group, _set_pos_IIB, 0, _join_II, 0)
        nodes_link(group, _sepa_bonds_III, 0, _join_III, 0)
        nodes_link(group, _set_pos_IIIA, 0, _join_III, 0)
        nodes_link(group, _set_pos_IIIB, 0, _join_III, 0)
        nodes_link(group, _join_II, 0, _output, 1)
        nodes_link(group, _join_III, 0, _output, 2)


class Mol3D_MT_Dashline(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_DASHED'
    bl_label='Dashed Bonds'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Scaffold', 'INPUT', 'NodeSocketGeometry')   # inputs
        self.add_new_socket('Radius', 'INPUT', 'NodeSocketFloat')
        self.add_new_socket('Dash Counts', 'INPUT', 'NodeSocketInt')
        self.add_new_socket('Dash Bonds', 'OUTPUT', 'NodeSocketGeometry')   # outputs
        
        # defalut values
        self.inputs[1].default_value = 0.1
        self.inputs[2].default_value = 7

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+200,b), '')
        # separate dashed line: bond_order = 0
        _split_edges = add_node(group, 'GeometryNodeSplitEdges', (a+200, b), '')
        _sepa_dash = add_node(group, 'GeometryNodeSeparateGeometry', (a+800,b), '')
        _sepa_dash.domain = 'EDGE'
        _attr_bond_order = add_attr_node(group, (a+600,b-200), '', [(0,'bond_order')], 'INT')
        _equal = add_compare_node(group, (a+800,b-200),'', [(3,0)], 'INT', 'EQUAL')
        nodes_link(group, _input, 0, _split_edges, 0)
        nodes_link(group, _split_edges, 0, _sepa_dash, 0)
        nodes_link(group, _attr_bond_order, v.na_sockets_out['int'], _equal, 2)
        nodes_link(group, _equal, 0, _sepa_dash, 1)

        # set attributes: dir, pos1, length
        _store_attr_dir, _store_attr_length = preset_attr_bonds_color(group, (a+200,b), _split_edges)

        a += 1600
        # instance cylinder on points
        _mesh2curve = add_node(group, 'GeometryNodeMeshToCurve', (a,b), '')
        _resample = add_node(group, 'GeometryNodeResampleCurve', (a+200,b), '')
        _add = add_math_node(group, (a,b-100), '', [(1,2)], 'ADD')
        nodes_link(group, _sepa_dash, 0, _store_attr_dir, 0)
        nodes_link(group, _store_attr_length, 0, _mesh2curve, 0)
        nodes_link(group, _mesh2curve, 0, _resample, 0)
        nodes_link(group, _input, 2, _add, 0)
        nodes_link(group, _add, 0, _resample, 2)

        _instance_on_points = add_node(group, 'GeometryNodeInstanceOnPoints', (a+400,b), '')
        _realize_instance = add_node(group, 'GeometryNodeRealizeInstances', (a+600,b), '')
        _attr_dir = add_attr_node(group, (a,b-200), '', [(0,'dir')], 'FLOAT_VECTOR')
        _align_euler2vec = add_node(group, 'FunctionNodeAlignEulerToVector',(a+200,b-200), '')
        _align_euler2vec.axis = 'Z'
        _attr_length = add_attr_node(group, (a-200,b-400), '', [(0,'length')], 'FLOAT')
        _add = add_math_node(group, (a-200,b-600), '', [(1,4)], 'ADD')
        _multiply = add_math_node(group, (a-200,b-600), '', [(1,0.5)], 'MULTIPLY')
        _divide = add_math_node(group, (a,b-400), '', [], 'DIVIDE')
        _combine_xyz = add_node(group, 'ShaderNodeCombineXYZ', (a+200,b-500), '')
        set_node_values(_combine_xyz,[(0,1),(1,1)])
        nodes_link(group, _resample, 0, _instance_on_points, 0)
        nodes_link(group, _instance_on_points, 0, _realize_instance, 0)
        nodes_link(group, _input, 2, _add, 0)
        nodes_link(group, _add, 0, _multiply, 0)
        nodes_link(group, _multiply, 0, _divide, 1)
        nodes_link(group, _attr_length, v.na_sockets_out['float'], _divide, 0)
        nodes_link(group, _divide, 0, _combine_xyz, 2)
        nodes_link(group, _combine_xyz, 0, _instance_on_points, 6)
        nodes_link(group, _attr_dir, v.na_sockets_out['vec'], _align_euler2vec, 2)
        nodes_link(group, _align_euler2vec, 0, _instance_on_points, 5)

        _cylinder = add_node(group, 'GeometryNodeMeshCylinder', (a+200,b+350), '')
        set_node_values(_cylinder, [(1,2),(4,1)])
        _spline_parameter = add_node(group, 'GeometryNodeSplineParameter', (a-200,b+300), '')
        _greater = add_compare_node(group, (a,b+500), '', [], 'FLOAT', 'GREATER_THAN')
        _less = add_compare_node(group, (a,b+300), '', [(1,1.0)], 'FLOAT', 'LESS_THAN')
        _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a+200,b+550), '')
        _boolean_AND.operation = 'AND'
        nodes_link(group, _input, 1, _cylinder, 3)
        nodes_link(group, _spline_parameter, 0, _greater, 0)
        nodes_link(group, _spline_parameter, 0, _less, 0)
        nodes_link(group, _greater, 0, _boolean_AND, 0)
        nodes_link(group, _less, 0, _boolean_AND, 1)
        nodes_link(group, _boolean_AND, 0, _instance_on_points, 1)
        nodes_link(group, _cylinder, 0, _instance_on_points, 2)

        _sepa_segs, _join_geo = add_attr_bonds_color(group, (a+400,b), _split_edges)
        _output.location = (a+1800,b)
        nodes_link(group, _realize_instance, 0, _sepa_segs, 0)
        nodes_link(group, _join_geo, 0, _output, 0)


def preset_attr_bonds_color(group, location, in_node):
    a,b = location
    _sample_nearest_I = add_node(group, 'GeometryNodeSampleNearest', (a+200,b-400), '')
    _sample_idx_I = add_sample_index_node(group, (a+400,b-400), '', 'FLOAT', 'POINT')
    _sample_nearest_II = add_node(group, 'GeometryNodeSampleNearest', (a+200,b-600), '')
    _sample_idx_II = add_sample_index_node(group, (a+400,b-600), '', 'FLOAT', 'POINT')
    _attr_atom_id = add_attr_node(group, (a,b-600), '', [(0,'atom_id')], 'INT')
    nodes_link(group, in_node, 0, _sample_nearest_I, 0)
    nodes_link(group, in_node, 0, _sample_nearest_II, 0)
    nodes_link(group, in_node, 0, _sample_idx_I, 0)
    nodes_link(group, in_node, 0, _sample_idx_II, 0)
    nodes_link(group, _attr_atom_id, v.na_sockets_out['int'], _sample_idx_I, v.si_sockets_in['int'])
    nodes_link(group, _attr_atom_id, v.na_sockets_out['int'], _sample_idx_II, v.si_sockets_in['int'])
    nodes_link(group, _sample_nearest_I, 0, _sample_idx_I, v.si_sockets_in['index'])
    nodes_link(group, _sample_nearest_II, 0, _sample_idx_II, v.si_sockets_in['index'])

    _edge_vertices = add_node(group, 'GeometryNodeInputMeshEdgeVertices', (a,b-800), '')
    _less_equal = add_compare_node(group, (a+600,b-500), '', [], 'INT', 'LESS_EQUAL')
    _switch_I = add_node(group, 'GeometryNodeSwitch', (a+600,b-700), '')
    _switch_I.input_type = 'VECTOR'
    _switch_II = add_node(group, 'GeometryNodeSwitch', (a+600,b-900), '')
    _switch_II.input_type = 'VECTOR'
    nodes_link(group, _edge_vertices, 2, _sample_nearest_I, 1)
    nodes_link(group, _edge_vertices, 3, _sample_nearest_II, 1)
    nodes_link(group, _edge_vertices, 2, _switch_I, v.sw_sockets_in['vecA'])
    nodes_link(group, _edge_vertices, 3, _switch_I, v.sw_sockets_in['vecB'])
    nodes_link(group, _edge_vertices, 2, _switch_II, v.sw_sockets_in['vecB'])
    nodes_link(group, _edge_vertices, 3, _switch_II, v.sw_sockets_in['vecA'])
    nodes_link(group, _sample_idx_I, 0, _less_equal, 2)
    nodes_link(group, _sample_idx_II, 0, _less_equal, 3)
    nodes_link(group, _less_equal, 0, _switch_I, 0)
    nodes_link(group, _less_equal, 0, _switch_II, 0)

    _subtract = add_vec_math_node(group, (a+800, b-500), '', [], 'SUBTRACT')
    _normalize = add_vec_math_node(group, (a+800,b-300), '', [], 'NORMALIZE')
    _distance = add_vec_math_node(group, (a+800,b-900), '', [], 'DISTANCE')
    _multiply = add_math_node(group, (a+1000,b-900), '', [(1,0.5)], 'MULTIPLY')
    nodes_link(group, _switch_I, v.sw_sockets_out['vec'], _subtract, 0)
    nodes_link(group, _switch_II, v.sw_sockets_out['vec'], _subtract, 1)
    nodes_link(group, _subtract, 0, _normalize, 0)
    nodes_link(group, _switch_I, v.sw_sockets_out['vec'], _distance, 0)
    nodes_link(group, _switch_II, v.sw_sockets_out['vec'], _distance, 1)
    nodes_link(group, _distance, 1, _multiply, 0)
    
    _store_attr_dir = add_store_attr_node(group, (a+800,b), 'dir', 'FLOAT_VECTOR', 'POINT', 'dir')
    _store_attr_pos1 = add_store_attr_node(group, (a+1000,b), 'pos1', 'FLOAT_VECTOR', 'POINT', 'pos1')
    _store_attr_pos2 = add_store_attr_node(group, (a+1000,b-250), 'pos1', 'FLOAT_VECTOR', 'POINT', 'pos2')
    _store_attr_length = add_store_attr_node(group, (a+1200,b), 'length', 'FLOAT', 'POINT', 'length')
    nodes_link(group, _store_attr_dir, 0, _store_attr_pos1, 0)
    nodes_link(group, _store_attr_pos1, 0, _store_attr_pos2, 0)
    nodes_link(group, _store_attr_pos2, 0, _store_attr_length, 0)
    nodes_link(group, _switch_I, v.sw_sockets_out['vec'], _store_attr_pos1, v.sa_sockets_in['vec'])
    nodes_link(group, _switch_II, v.sw_sockets_out['vec'], _store_attr_pos2, v.sa_sockets_in['vec'])
    nodes_link(group, _normalize, 0, _store_attr_dir, v.sa_sockets_in['vec'])
    nodes_link(group, _multiply, 0, _store_attr_length, v.sa_sockets_in['float'])

    return (_store_attr_dir, _store_attr_length)

def add_attr_bonds_color(group, location, in_node):
    # set color attribute
    a,b = location
    _sepa_segs = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b), '')
    _sepa_segs.domain = 'FACE'
    _store_attr_color_I = add_store_attr_node(group, (a+800,b), '', 'FLOAT_COLOR', 'FACE', 'colour')
    _store_attr_color_II = add_store_attr_node(group, (a+800,b+300), '', 'FLOAT_COLOR', 'FACE', 'colour')
    _join_geo = add_node(group, 'GeometryNodeJoinGeometry', (a+1000,b+150), '')
    
    nodes_link(group, _sepa_segs, 0, _store_attr_color_I, 0)
    nodes_link(group, _sepa_segs, 1, _store_attr_color_II, 0)
    nodes_link(group, _store_attr_color_I, 0, _join_geo, 0)
    nodes_link(group, _store_attr_color_II, 0, _join_geo, 0)

    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a,b+400), '')
    _attr_pos1 = add_attr_node(group, (a,b+200), '', [(0,'pos1')], 'FLOAT_VECTOR')
    _attr_length = add_attr_node(group, (a+200,b+200), '', [(0,'length')], 'FLOAT')
    _distance = add_vec_math_node(group, (a+200,b+400), '', [], 'DISTANCE')
    _less = add_compare_node(group, (a+400,b+200), '', [], 'FLOAT', 'LESS_THAN')
    _boolean_NOT = add_node(group, 'FunctionNodeBooleanMath', (a+600,b+200), '')
    _boolean_NOT.operation = 'NOT'
    nodes_link(group, _input_pos, 0, _distance, 0)
    nodes_link(group, _attr_pos1, v.na_sockets_out['vec'], _distance, 1)
    nodes_link(group, _distance, 1, _less, 0)
    nodes_link(group, _attr_length, v.na_sockets_out['float'], _less, 1)
    nodes_link(group, _less, 0, _boolean_NOT, 0)
    nodes_link(group, _less, 0, _sepa_segs, 1)
    nodes_link(group, _less, 0, _store_attr_color_I, 1)
    nodes_link(group, _boolean_NOT, 0, _store_attr_color_II, 1)

    _attr_pos1 = add_attr_node(group, (a+200,b-400), '', [(0,'pos1')], 'FLOAT_VECTOR')
    _attr_pos2 = add_attr_node(group, (a+200,b-200), '', [(0,'pos2')], 'FLOAT_VECTOR')
    _attr_color = add_attr_node(group, (a+400,b-600), '', [(0,'colour')], 'FLOAT_COLOR')
    _sample_nearest_I = add_node(group, 'GeometryNodeSampleNearest', (a+400,b-400), '')
    _sample_idx_I = add_sample_index_node(group, (a+600,b-500), '', 'FLOAT_COLOR', 'POINT')
    _sample_nearest_II = add_node(group, 'GeometryNodeSampleNearest', (a+400,b-200), '')
    _sample_idx_II = add_sample_index_node(group, (a+600,b-200), '', 'FLOAT_COLOR', 'POINT')
    nodes_link(group, in_node, 0, _sample_nearest_I, 0)
    nodes_link(group, in_node, 0, _sample_nearest_II, 0)
    nodes_link(group, in_node, 0, _sample_idx_I, 0)
    nodes_link(group, in_node, 0, _sample_idx_II, 0)
    nodes_link(group, _attr_pos1, 0, _sample_nearest_I, 1)
    nodes_link(group, _attr_pos2, 0, _sample_nearest_II, 1)
    nodes_link(group, _sample_nearest_I, 0, _sample_idx_I, v.si_sockets_in['index'])
    nodes_link(group, _sample_nearest_II, 0, _sample_idx_II, v.si_sockets_in['index'])
    nodes_link(group, _attr_color, v.na_sockets_out['color'], _sample_idx_I, v.si_sockets_in['color'])
    nodes_link(group, _attr_color, v.na_sockets_out['color'], _sample_idx_II, v.si_sockets_in['color'])
    nodes_link(group, _sample_idx_I, v.si_sockets_out['color'], _store_attr_color_I, v.sa_sockets_in['color'])
    nodes_link(group, _sample_idx_II, v.si_sockets_out['color'], _store_attr_color_II, v.sa_sockets_in['color'])

    return (_sepa_segs, _join_geo)


class Mol3D_MT_Aromatic_Bonds(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_AROMATIC'
    bl_label='Aromatic Bonds'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Scaffold', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Offset', 'INPUT', 'NodeSocketFloat')   # socket 1
        self.add_new_socket('Radius', 'INPUT', 'NodeSocketFloat')   # socket 2
        self.add_new_socket('Dash Counts', 'INPUT', 'NodeSocketInt')   # socket 3
        self.add_new_socket('Length Scale', 'INPUT', 'NodeSocketFloat')   # socket 4
        self.add_new_socket('Dash Bonds', 'OUTPUT', 'NodeSocketGeometry')   # outputs
        
        # defalut values
        self.inputs[1].default_value = 0.25
        self.inputs[2].default_value = 0.05
        self.inputs[3].default_value = 3
        self.inputs[4].default_value = 1.05

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+1400,b), '')

        _sepa_aromatic = add_node(group, 'GeometryNodeSeparateGeometry', (a+200,b), '')
        _sepa_aromatic.domain = 'EDGE'
        _split_edges = add_node(group, 'GeometryNodeSplitEdges', (a+400, b), '')
        _set_pos = add_node(group, 'GeometryNodeSetPosition', (a+600,b), '')
        _scale_elements = add_node(group, 'GeometryNodeScaleElements', (a+800,b), '')
        _scale_elements.domain = 'EDGE'
        _store_attr_bond_order = add_store_attr_node(group, (a+1000,b), '', 'INT', 'EDGE', 'bond_order')
        _dashed_bonds = add_node(group, 'Mol3D_MT_Dashline', (a+1200,b), '')
        nodes_link(group, _input, 0, _sepa_aromatic, 0)
        nodes_link(group, _sepa_aromatic, 0, _split_edges, 0)
        nodes_link(group, _split_edges, 0, _set_pos, 0)
        nodes_link(group, _set_pos, 0, _scale_elements, 0)
        nodes_link(group, _scale_elements, 0, _store_attr_bond_order, 0)
        nodes_link(group, _store_attr_bond_order, 0, _dashed_bonds, 0)
        nodes_link(group, _dashed_bonds, 0, _output, 0)

        _attr_is_aromatic = add_attr_node(group, (a,b-200), '', [(0,'is_aromatic_edge')], 'BOOLEAN')
        _attr_bond_vertvec = add_attr_node(group, (a,b-400), '', [(0,'bond_vertvec')], 'FLOAT_VECTOR')
        _scale = add_vec_math_node(group, (a+200,b-200), '', [], 'SCALE')
        _subtract = add_math_node(group, (a+200,b-400), '', [], 'SUBTRACT')
        nodes_link(group, _attr_is_aromatic, v.na_sockets_out['bool'], _sepa_aromatic, 1)
        nodes_link(group, _attr_bond_vertvec, 0, _scale, 0)
        nodes_link(group, _scale, 0, _set_pos, 3)
        nodes_link(group, _subtract, 0, _scale_elements, 2)
        nodes_link(group, _input, 1, _scale, 3)
        nodes_link(group, _input, 1, _subtract, 1)
        nodes_link(group, _input, 2, _dashed_bonds, 1)
        nodes_link(group, _input, 3, _dashed_bonds, 2)
        nodes_link(group, _input, 4, _subtract, 0)
        
        
class Mol3D_MT_Nucleic_Restruct(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_NUCLEIC_MESHLINE'
    bl_label='Nucleic Meshline'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Nucleic', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Cutoff', 'INPUT', 'NodeSocketFloat')   # inputs socket 1
        self.add_new_socket('Points', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 0
        self.add_new_socket('Line', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1

        # defalut values
        self.inputs[1].default_value = 7.7
        
        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+1400,b), '')
        _sepa_nucleic = sepa_by_attr(group, (a+200,b), 'POINT', 'sec_struct', 'INT', 0, 'EQUAL')
        _attr_atom_name = add_attr_node(group, (a+200,b+200), '', [(0,'atom_name')], 'INT')
        _equal_1 = add_compare_node(group, (a+400,b+100),'', [(3,57)], 'INT', 'EQUAL')   # C3 ring
        _equal_2 = add_compare_node(group, (a+400,b+300),'', [(3,55)], 'INT', 'EQUAL')   # C4 ring
        _sepa_geo_1 = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b+100), '')
        _sepa_geo_2 = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b+300), '')
        _sample_idx_1 = add_sample_index_node(group, (a+800,b+100), '', 'FLOAT_VECTOR', 'POINT')
        _sample_idx_2 = add_sample_index_node(group, (a+800,b+350), '', 'FLOAT_VECTOR', 'POINT')
        _input_pos = add_node(group, 'GeometryNodeInputPosition', (a+600,b+400), '')
        _input_idx = add_node(group, 'GeometryNodeInputIndex', (a+600,b+500), '')
        _vec_mix_1 = add_node(group, 'ShaderNodeMix', (a+1000, b+100), '')
        set_node_datatype(_vec_mix_1, 'VECTOR')
        _vec_mix_2 = add_node(group, 'ShaderNodeMix', (a+1000, b+350), '')
        set_node_datatype(_vec_mix_2, 'VECTOR')
        _sepa_C1 = sepa_by_attr(group, (a+1200,b), 'POINT', 'atom_name', 'INT', 61, 'EQUAL')   # C1 ring # connection to the base
        _sample_idx = add_sample_index_node(group, (a+1400,b+300), '', 'FLOAT_VECTOR', 'POINT')
        _domain_size = add_node(group, 'GeometryNodeAttributeDomainSize', (a+1400,b), '')
        _mesh_line = add_node(group, 'GeometryNodeMeshLine', (a+1600,b), '')
        _set_pos = add_node(group, 'GeometryNodeSetPosition', (a+1800,b), '')
        _delete_mismatch, _greater_than = delete_mismatch(group, (a+2000,b), max_length=7.7)
        nodes_link(group, _input, 0, _sepa_nucleic, 0)
        nodes_link(group, _input, 1, _greater_than, 1)
        nodes_link(group, _attr_atom_name, v.na_sockets_out['int'], _equal_1, 2)
        nodes_link(group, _attr_atom_name, v.na_sockets_out['int'], _equal_2, 2)
        nodes_link(group, _equal_1, 0, _sepa_geo_1, 1)
        nodes_link(group, _equal_2, 0, _sepa_geo_2, 1)
        nodes_link(group, _sepa_nucleic, 0, _sepa_geo_1, 0)
        nodes_link(group, _sepa_nucleic, 0, _sepa_geo_2, 0)
        nodes_link(group, _input_pos, 0, _sample_idx_1, v.si_sockets_in['vec'])
        nodes_link(group, _input_pos, 0, _sample_idx_2, v.si_sockets_in['vec'])
        nodes_link(group, _input_pos, 0, _vec_mix_2, 5)
        nodes_link(group, _input_idx, 0, _sample_idx_1, v.si_sockets_in['index'])
        nodes_link(group, _input_idx, 0, _sample_idx_2, v.si_sockets_in['index'])
        nodes_link(group, _sepa_geo_1, 0, _sample_idx_1, 0)
        nodes_link(group, _sepa_geo_2, 0, _sample_idx_2, 0)
        nodes_link(group, _sample_idx_1, v.si_sockets_out['vec'], _vec_mix_1, 4)
        nodes_link(group, _sample_idx_2, v.si_sockets_out['vec'], _vec_mix_1, 5)
        nodes_link(group, _vec_mix_1, 1, _vec_mix_2, 4)
        nodes_link(group, _vec_mix_2, 1, _sample_idx, v.si_sockets_in['vec'])
        nodes_link(group, _input_idx, 0, _sample_idx, v.si_sockets_in['index'])
        nodes_link(group, _sepa_nucleic, 0, _sepa_C1, 0)
        nodes_link(group, _sepa_C1, 0, _sample_idx, 0)
        nodes_link(group, _sepa_C1, 0, _domain_size, 0)
        nodes_link(group, _domain_size, 0, _mesh_line, 0)
        nodes_link(group, _mesh_line, 0, _set_pos, 0)
        nodes_link(group, _sample_idx, v.si_sockets_out['vec'], _set_pos, 2)
        nodes_link(group, _set_pos, 0, _delete_mismatch, 0)
        nodes_link(group, _sepa_C1, 0, _output, 0)
        nodes_link(group, _delete_mismatch, 0, _output, 1)


class Mol3D_MT_Peptide_Restruct(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_PEPTIDE_MESHLINE'
    bl_label='Peptide Meshline'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Peptide', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Cutoff', 'INPUT', 'NodeSocketFloat')   # inputs socket 1
        self.add_new_socket('Points', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 0
        self.add_new_socket('Line', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1

        # defalut values
        self.inputs[1].default_value = 7.5
        
        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+1200,b), '')
        _sepa_peptide = add_node(group, 'GeometryNodeSeparateGeometry', (a+200,b), '')
        _attr_is_alpha = add_attr_node(group, (a+200,b-200), '', [(0,'is_alpha')], 'BOOLEAN')
        _input_pos = add_node(group, 'GeometryNodeInputPosition', (a+200,b+200), '')
        _input_idx = add_node(group, 'GeometryNodeInputIndex', (a+200,b+300), '')
        _sample_idx = add_sample_index_node(group, (a+400,b+300), '', 'FLOAT_VECTOR', 'POINT')
        _domain_size = add_node(group, 'GeometryNodeAttributeDomainSize', (a+400,b), '')
        _mesh_line = add_node(group, 'GeometryNodeMeshLine', (a+600,b), '')
        _set_pos = add_node(group, 'GeometryNodeSetPosition', (a+800,b), '')
        _delete_mismatch, _greater_than = delete_mismatch(group, (a+1000,b), max_length=7.5)
        nodes_link(group, _input, 0, _sepa_peptide, 0)
        nodes_link(group, _input, 1, _greater_than, 1)
        nodes_link(group, _attr_is_alpha, v.na_sockets_out['bool'], _sepa_peptide, 1)
        nodes_link(group, _sepa_peptide, 0, _sample_idx, 0)
        nodes_link(group, _sepa_peptide, 0, _domain_size, 0)
        nodes_link(group, _input_pos, 0, _sample_idx, v.si_sockets_in['vec'])
        nodes_link(group, _input_idx, 0, _sample_idx, v.si_sockets_in['index'])
        nodes_link(group, _domain_size, 0, _mesh_line, 0)
        nodes_link(group, _mesh_line, 0, _set_pos, 0)
        nodes_link(group, _sample_idx, v.si_sockets_out['vec'], _set_pos, 2)
        nodes_link(group, _set_pos, 0, _delete_mismatch, 0)
        nodes_link(group, _sepa_peptide, 0, _output, 0)
        nodes_link(group, _delete_mismatch, 0, _output, 1)

attributes = (
    {'name': 'res_id',        'type':'INT',      'domain':'POINT'},
    {'name': 'atom_id',       'type':'INT',      'domain':'POINT'},
    {'name': 'chain_id',      'type':'INT',      'domain':'POINT'},
    #{'name': 'entity_id',     'type':'INT',      'domain':'POINT'},
    {'name': 'sec_struct',    'type':'INT',      'domain':'POINT'},
    {'name': 'scale_beta',    'type':'FLOAT',    'domain':'POINT'},
    {'name': 'coil_end',      'type':'FLOAT',    'domain':'POINT'},

    {'name': 'guide_X',      'type':'FLOAT_VECTOR',   'domain':'POINT'},
    {'name': 'guide_Z',      'type':'FLOAT_VECTOR',   'domain':'POINT'},
    
    {'name': 'index',         'type':'FLOAT',    'domain':'POINT'},
    {'name': 'b_factor',      'type':'FLOAT',    'domain':'POINT'},
    #{'name': 'occupancy',     'type':'FLOAT',    'domain':'POINT'},
    
    {'name': 'is_alpha',      'type':'BOOLEAN',  'domain':'POINT'},
    {'name': 'is_nucleic',    'type':'BOOLEAN',  'domain':'POINT'},
    {'name': 'is_peptide',    'type':'BOOLEAN',  'domain':'POINT'},
    {'name': 'is_hetero',     'type':'BOOLEAN',  'domain':'POINT'},
    {'name': 'is_carb',       'type':'BOOLEAN',  'domain':'POINT'},
    {'name': 'is_backbone',   'type':'BOOLEAN',  'domain':'POINT'},

    {'name': 'bond_order',    'type':'INT',      'domain':'EDGE'},
    {'name': 'bond_aids_1',   'type':'INT',      'domain':'EDGE'},
    {'name': 'bond_aids_2',   'type':'INT',      'domain':'EDGE'},
    {'name': 'bond_radii_B',  'type':'FLOAT',    'domain':'EDGE'},
    {'name': 'bond_radii_S',  'type':'FLOAT',    'domain':'EDGE'},
    {'name': 'bond_radii_W',  'type':'FLOAT',    'domain':'EDGE'},
    {'name': 'radii_B',       'type':'FLOAT',    'domain':'POINT'},
    {'name': 'radii_S',       'type':'FLOAT',    'domain':'POINT'},
    {'name': 'radii_W',       'type':'FLOAT',    'domain':'POINT'},
    {'name': 'radii_F',       'type':'FLOAT',    'domain':'POINT'},
    {'name': 'colour',        'type':'FLOAT_COLOR',   'domain':'POINT'},
)

class Mol3D_MT_BioAttr_Reset(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_ATTRIBUTE_RESET'
    bl_label='Bioattribute Reset'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('Points', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Line', 'INPUT', 'NodeSocketGeometry')   # inputs socket 1
        self.add_new_socket('BioScaffold', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 0
    
        # defalut values
        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+200,b), '')

        for i, attr in enumerate(attributes):
            _sample_idx = _node.add_attr_trans_nodes(group, attr['name'], attr['type'], attr['domain'], (a+200*i,b), f"Attr_{'%03d' % i}")[0]
            nodes_link(group, _input, 0, _sample_idx, 0)
        for i in range(len(attributes)-1):
            nodes_link(group, get_node(group, f"Attr_{'%03d' % i}"), 0, get_node(group, f"Attr_{'%03d' % (i+1)}"), 0)
        nodes_link(group, _input, 1, get_node(group, 'Attr_000'), 0)
        _last_attr = get_node(group, f"Attr_{'%03d' % (len(attributes)-1)}")
        
        a = a+200*len(attributes)+200
        _store_new_guide_Z = _node.guide_Z_revised(group, (a,b))
        _store_angle_XZ = _node.angle_XZ(group, (a+200,b))
        _output.location = (a+400,b)
        nodes_link(group, _last_attr, 0, _store_new_guide_Z, 0)
        nodes_link(group, _store_new_guide_Z, 0, _store_angle_XZ, 0)
        nodes_link(group, _store_angle_XZ, 0, _output, 0)


class Mol3D_MT_BioAttr_Preset(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_ATTRIBUTE_PRESET'
    bl_label='Bioattribute Preset'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('BioScaffold', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Nucleic', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 0
        self.add_new_socket('Peptide', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1
    
        # defalut values
        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a-200,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+1200,b), '')

        _attr_res_id = add_attr_node(group, (a,b+150), '', [(0,'res_id')], 'INT')
        _res_idx = add_store_attr_node(group, (a,b), '', 'INT', 'POINT', 'index')
        _guide_X, _guide_Z = _node.add_guide_vector_nodes(group, (a+200,b), is_nucleic=False)
        _sepa_alpha_C = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b), '')
        _attr_is_alpha = add_attr_node(group, (a+400, b-200), '', [(0,'is_alpha')], 'BOOLEAN')
        _sample_index = add_node(group, 'GeometryNodeSampleIndex', (a+400, b+200), '')
        set_node_datatype(_sample_index, 'FLOAT_VECTOR')
        _scale_beta = _node.add_scale_beta_nodes(group, (a+800,b))
        _coil_end = _node.add_coil_end_nodes(group, (a+1000,b+500))
        _guide_Z_inverse = _node.add_guide_Z_inverse_nodes(group, (a+1200,b))
        
        nodes_link(group, _input, 0, _res_idx, 0)
        nodes_link(group, _attr_res_id, v.na_sockets_out['int'], _res_idx, v.sa_sockets_in['int'])
        nodes_link(group, _res_idx, 0, _guide_X, 0)
        nodes_link(group, _guide_Z, 0, _sepa_alpha_C, 0)
        nodes_link(group, _sepa_alpha_C, 0, _scale_beta, 0)
        nodes_link(group, _attr_is_alpha, v.na_sockets_out['bool'], _sepa_alpha_C, 1)
        nodes_link(group, _scale_beta, 0, _coil_end, 0)
        nodes_link(group, _coil_end, 0, _guide_Z_inverse, 0)
        nodes_link(group, _guide_Z, 0, _output, 0)
        nodes_link(group, _guide_Z_inverse, 0, _output, 1)


class Mol3D_MT_Alpha_Helix(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_ALPHA_HELIX'
    bl_label='Alpha Helix'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('BioScaffold', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Resolution', 'INPUT', 'NodeSocketInt')   # inputs socket 1
        self.add_new_socket('Width', 'INPUT', 'NodeSocketFloat')   # inputs socket 2
        self.add_new_socket('Thickness', 'INPUT', 'NodeSocketFloat')   # inputs socket 3
        self.add_new_socket('Sharp/Smooth', 'INPUT', 'NodeSocketBool')   # inputs socket 4
        self.add_new_socket('Alpha Helix', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1
    
        # defalut values
        self.inputs[1].default_value = 6
        self.inputs[2].default_value = 2
        self.inputs[3].default_value = 0.5
        self.inputs[4].default_value = False

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+2000,b), '')
        _sepa_alpha_helix = sepa_by_attr(group, (a+200,b), 'POINT', 'sec_struct', 'INT', 1, 'EQUAL')
        _mesh2curve_alpha, _curve_res, _resample_curve_alpha = mesh_to_curve(group, (a+400,b), resolution=6)
        _curve2mesh_alpha = add_node(group, 'GeometryNodeCurveToMesh', (a+1400,b), '')
        set_node_values(_curve2mesh_alpha, [(2,True)])    # Fill Caps True
        _instance_alpha, _sample_index_alpha = _node.cross_instance(group, (a+1800,b-200))
        _switch = add_node(group, 'GeometryNodeSwitch', (a+400,b-200), '')
        _switch.input_type = 'INT'
        set_node_values(_switch, [(v.sw_sockets_in['intB'],4)])
        _multiply = add_math_node(group, (a+400,b-400), '', [(1,8)], 'MULTIPLY')
        _curve_circle, _cross_section = _node.cross_section_transform(group, False, 2, 0.5, (a+600,b-200), 6)
        _combine_xyz = add_node(group, 'ShaderNodeCombineXYZ', (a+600,b-400), '')
        set_node_values(_combine_xyz, [(2,1.0)])
        _radians = add_math_node(group, (a+600,b-600), '', [], 'RADIANS')
        _attr_angle_XZ = add_attr_node(group, (a+400,b-600), '', [(0,'angle_XZ')], 'FLOAT')
        _rotate_euler = _node.rotate_euler(group, angle=45*3.1416/180, location=(a+800,b-600))
        _reset_pos_alpha = add_node(group, 'GeometryNodeSetPosition', (a+1800,b), '')
        
        nodes_link(group, _input, 0, _sepa_alpha_helix, 0)
        nodes_link(group, _input, 1, _multiply, 0)
        nodes_link(group, _input, 1, _curve_res, 2)
        nodes_link(group, _input, 2, _combine_xyz, 0)
        nodes_link(group, _input, 3, _combine_xyz, 1)
        nodes_link(group, _input, 4, _switch, 0)
        nodes_link(group, _attr_angle_XZ, v.na_sockets_out['float'], _radians, 0)
        nodes_link(group, _sepa_alpha_helix, 0, _mesh2curve_alpha, 0)
        nodes_link(group, _resample_curve_alpha, 0, _curve2mesh_alpha, 0)
        nodes_link(group, _combine_xyz, 0, _cross_section, 3)
        nodes_link(group, _resample_curve_alpha, 0, _instance_alpha, 0)
        nodes_link(group, _multiply, 0, _switch, v.sw_sockets_in['intA'])
        nodes_link(group, _switch, v.sw_sockets_out['int'], _curve_circle, 0)
        nodes_link(group, _curve_circle, 0, _curve2mesh_alpha, 1)
        nodes_link(group, _cross_section, 0, _instance_alpha, 2)
        nodes_link(group, _radians, 0, _rotate_euler, 3)
        nodes_link(group, _rotate_euler, 0, _instance_alpha, 5)
        nodes_link(group, _curve2mesh_alpha, 0, _reset_pos_alpha, 0)
        nodes_link(group, _sample_index_alpha, v.si_sockets_out['vec'], _reset_pos_alpha, 2)
        nodes_link(group, _reset_pos_alpha, 0, _output, 0)


class Mol3D_MT_Beta_Sheet(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_BETA_SHEET'
    bl_label='Beta Sheet'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('BioScaffold', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Resolution', 'INPUT', 'NodeSocketInt')   # inputs socket 1
        self.add_new_socket('Width', 'INPUT', 'NodeSocketFloat')   # inputs socket 2
        self.add_new_socket('Thickness', 'INPUT', 'NodeSocketFloat')   # inputs socket 3
        self.add_new_socket('Arrow Multiply', 'INPUT', 'NodeSocketFloat')   # inputs socket 4
        self.add_new_socket('Sharp/Smooth', 'INPUT', 'NodeSocketBool')   # inputs socket 5
        self.add_new_socket('Rotate Angle', 'INPUT', 'NodeSocketFloat')   # inputs socket 6
        self.add_new_socket('Arrow Offset', 'INPUT', 'NodeSocketFloat')   # inputs socket 7
        self.add_new_socket('Beta Sheet', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1
    
        # defalut values
        self.inputs[1].default_value = 6
        self.inputs[2].default_value = 1.8
        self.inputs[3].default_value = 0.5
        self.inputs[4].default_value = 2.0
        self.inputs[5].default_value = False
        self.inputs[6].default_value = 0
        self.inputs[7].default_value = 0.0

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a-200,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+2000,b), '')
        _sepa_beta_sheet = sepa_by_attr(group, (a,b), 'POINT', 'sec_struct', 'INT', 2, 'EQUAL')
        _mesh2curve_beta, _curve_res, _resample_curve_beta = mesh_to_curve(group, (a+200,b), 6)
        _curve2mesh_beta = add_node(group, 'GeometryNodeCurveToMesh', (a+1200, b), '')
        set_node_values(_curve2mesh_beta, [(2,True)])    # Fill Caps True
        _set_pos_beta, _mix_vec, _equal = _node.beta_arrow_bottom_shift(group, 0.0, (a+1000,b-200), 6, 2.0)
        _switch_set_pos = add_node(group, 'GeometryNodeSwitch', (a+1200,b+200), '')
        _subtract = add_math_node(group, (a+200,b-600), '', [(0,1)], 'SUBTRACT')
        _divide = add_math_node(group, (a+200,b-800), '', [(0,1)], 'DIVIDE')
        _instance_beta, _sample_index_beta = _node.cross_instance(group, (a+1200, b-200))
        _switch = add_node(group, 'GeometryNodeSwitch', (a+400,b-600), '')
        _switch.input_type = 'INT'
        set_node_values(_switch, [(v.sw_sockets_in['intB'],4)])
        _multiply = add_math_node(group, (a+400,b-800), '', [(1,8)], 'MULTIPLY')
        _curve_circle_beta, _cross_section_beta = _node.cross_section_transform(group, False, 2, 0.5, (a+600,b-600), 6)
        _combine_xyz = add_node(group, 'ShaderNodeCombineXYZ', (a+600,b-800), '')
        set_node_values(_combine_xyz, [(2,1.0)])
        _radians = add_math_node(group, (a+400,b-1000), '', [], 'RADIANS')
        _rotate_euler_beta = _node.rotate_euler(group, angle=0*3.1416/180, location=(a+600,b-1000))
        _reset_pos_beta = add_node(group, 'GeometryNodeSetPosition', (a+1400, b), '')

        # scale_beta adjustment
        _attr_scale_beta = add_attr_node(group, (a+400,b-1500), '', [(0,'scale_beta')], 'FLOAT')
        _equal_1 = add_compare_node(group, (a+600,b-1400), 'Equal', [(1,1.0)], 'FLOAT', 'EQUAL')   # greater than critical
        _multi_add = add_math_node(group, (a+600,b-1600), '', [(1,2.0),(2,0.01)], 'MULTIPLY_ADD')   # multiply = 2.0
        _switch_1 = add_node(group, 'GeometryNodeSwitch', (a+800,b-1500), '')
        _switch_1.input_type = 'FLOAT'
        set_node_values(_switch_1, [(v.sw_sockets_in['floatB'],1.0)])
        _combine = add_node(group, 'ShaderNodeCombineXYZ', (a+1000,b-1500), '')
        set_node_values(_combine, [(0,1.0),(1,1.0)])
        _less_equal = add_compare_node(group, (a+800,b-1300), '', [(1,1.0)], 'FLOAT', 'LESS_EQUAL')
        _input_int = add_node(group, 'FunctionNodeInputInt', (a+1000,b-1300), '')
        _input_int.integer = 1
        _switch_2 = add_node(group, 'GeometryNodeSwitch', (a+1200,b-1300), '')
        _switch_2.input_type = 'VECTOR'
        
        nodes_link(group, _input, 0, _sepa_beta_sheet, 0)
        nodes_link(group, _input, 1, _multiply, 0)
        nodes_link(group, _input, 1, _curve_res, 2)
        nodes_link(group, _input, 1, _divide, 1)
        nodes_link(group, _input, 2, _combine_xyz, 0)
        nodes_link(group, _input, 3, _combine_xyz, 1)
        nodes_link(group, _input, 4, _multi_add, 1)
        nodes_link(group, _input, 4, _less_equal, 0)
        nodes_link(group, _input, 5, _switch, 0)
        nodes_link(group, _input, 6, _radians, 0)
        nodes_link(group, _input, 7, _mix_vec, 0)
        nodes_link(group, _sepa_beta_sheet, 0, _mesh2curve_beta, 0)
        nodes_link(group, _resample_curve_beta, 0, _curve2mesh_beta, 0)
        nodes_link(group, _resample_curve_beta, 0, _set_pos_beta, 0)
        nodes_link(group, _resample_curve_beta, 0, _switch_set_pos, v.sw_sockets_in['geoB'])
        nodes_link(group, _set_pos_beta, 0, _switch_set_pos, v.sw_sockets_in['geoA'])
        nodes_link(group, _divide, 0, _subtract, 1)
        nodes_link(group, _subtract, 0, _equal, 1)
        nodes_link(group, _combine_xyz, 0, _cross_section_beta, 3)
        nodes_link(group, _switch_set_pos, v.sw_sockets_out['geo'], _instance_beta, 0)
        nodes_link(group, _multiply, 0, _switch, v.sw_sockets_in['intA'])
        nodes_link(group, _switch, v.sw_sockets_out['int'], _curve_circle_beta, 0)
        nodes_link(group, _curve_circle_beta, 0, _curve2mesh_beta, 1)
        nodes_link(group, _cross_section_beta, 0, _instance_beta, 2)
        nodes_link(group, _radians, 0, _rotate_euler_beta, 3)
        nodes_link(group, _rotate_euler_beta, 0, _instance_beta, 5)
        nodes_link(group, _curve2mesh_beta, 0, _reset_pos_beta, 0)
        nodes_link(group, _sample_index_beta, v.si_sockets_out['vec'], _reset_pos_beta, 2)

        nodes_link(group, _attr_scale_beta, v.na_sockets_out['float'], _equal_1, 0)
        nodes_link(group, _attr_scale_beta, v.na_sockets_out['float'], _multi_add, 0)
        nodes_link(group, _equal_1, 0, _switch_1, 0)
        nodes_link(group, _multi_add, 0, _switch_1, v.sw_sockets_in['floatA'])
        nodes_link(group, _switch_1, 0, _combine, 1)
        nodes_link(group, _switch_1, 0, _combine, 2)
        nodes_link(group, _less_equal, 0, _switch_2, 0)
        nodes_link(group, _less_equal, 0, _switch_set_pos, v.sw_sockets_in['switch_1'])
        nodes_link(group, _input_int, 0, _switch_2, v.sw_sockets_in['vecB'])
        nodes_link(group, _combine, 0, _switch_2, v.sw_sockets_in['vecA'])
        nodes_link(group, _switch_2, v.sw_sockets_out['vec'], _instance_beta, 6)
        nodes_link(group, _reset_pos_beta, 0, _output, 0)



class Mol3D_MT_Loop_Coil(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_LOOPS'
    bl_label='Loop Coil'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('BioScaffold', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Resolution', 'INPUT', 'NodeSocketInt')   # inputs socket 1
        self.add_new_socket('Radius', 'INPUT', 'NodeSocketFloat')   # inputs socket 2
        self.add_new_socket('Loop Coil', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1
    
        # defalut values
        self.inputs[1].default_value = 8
        self.inputs[2].default_value = 0.3

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a-200,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+1400,b), '')
        _sepa_coil = add_node(group, 'GeometryNodeSeparateGeometry', (a,b), 'Coil')
        _attr_sec_struct = add_attr_node(group, (a-200,b-200), '', [(0,'sec_struct')], 'INT')
        _coil_sec_equal = add_compare_node(group, (a,b-200), '', [(3, 3)], 'INT', 'EQUAL')
        _attr_coil_end = add_attr_node(group, (a+200,b-600), '', [(0,'coil_end')], 'FLOAT')
        _coil_end_equal = add_compare_node(group, (a+200,b-400), '', [(1, 1)], 'FLOAT', 'EQUAL')
        _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (a+200,b-200), '')
        set_node_operation(_boolean_OR, 'OR')
        _store_attr_sec = add_node(group, 'GeometryNodeStoreNamedAttribute', (a+200,b+200), '')
        _store_attr_sec.data_type = 'INT'
        set_node_values(_store_attr_sec, [(2,'sec_struct'),(v.sa_sockets_in['int'],3)])
        _mesh2curve_coil, _curve_res, _resample_curve_coil = mesh_to_curve(group, (a+200,b), 6)
        _curve2mesh_coil = add_node(group, 'GeometryNodeCurveToMesh', (a+1200,b), '')
        set_node_values(_curve2mesh_coil, [(2,True)])    # Fill Caps True
        _curve_circle_coil = add_node(group, 'GeometryNodeCurvePrimitiveCircle', (a+1000, b-200), '')
        set_node_values(_curve_circle_coil, [(0,8),(4,0.3)])
        nodes_link(group, _input, 0, _sepa_coil, 0)
        nodes_link(group, _input, 1, _curve_res, 2)
        nodes_link(group, _input, 1, _curve_circle_coil, 0)
        nodes_link(group, _input, 2, _curve_circle_coil, 4)
        nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _coil_sec_equal, 2)
        nodes_link(group, _attr_coil_end, v.na_sockets_out['float'], _coil_end_equal, 0)
        nodes_link(group, _coil_sec_equal, 0, _boolean_OR, 0)
        nodes_link(group, _coil_end_equal, 0, _boolean_OR, 1)
        nodes_link(group, _boolean_OR, 0, _sepa_coil, 1)
        nodes_link(group, _sepa_coil, 0, _store_attr_sec, 0)
        nodes_link(group, _store_attr_sec, 0, _mesh2curve_coil, 0)
        nodes_link(group, _resample_curve_coil, 0, _curve2mesh_coil, 0)
        nodes_link(group, _curve_circle_coil, 0, _curve2mesh_coil, 1)
        nodes_link(group, _curve2mesh_coil, 0, _output, 0)


class Mol3D_MT_Nucleic(bpy.types.GeometryNodeCustomGroup):
    bl_name='MOL3D_MT_NUCLEIC'
    bl_label='Nucleic'
    bl_icon = 'NONE'

    def init(self, context):
        self.createNodetree(self.name)

    def add_new_socket(self, socket_name, in_out, socket_type):
        self.node_tree.interface.new_socket(socket_name, in_out=in_out, socket_type=socket_type)
  
    def createNodetree(self, name):
        self.node_tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
        group = self.node_tree
        self.add_new_socket('BioScaffold', 'INPUT', 'NodeSocketGeometry')   # inputs socket 0
        self.add_new_socket('Resolution', 'INPUT', 'NodeSocketInt')   # inputs socket 1
        self.add_new_socket('Radius', 'INPUT', 'NodeSocketFloat')   # inputs socket 2
        self.add_new_socket('Loop Coil', 'OUTPUT', 'NodeSocketGeometry')   # outputs socket 1
    
        # defalut values
        self.inputs[1].default_value = 8
        self.inputs[2].default_value = 0.3

        a,b=(0,0)
        _input = add_node(group, 'NodeGroupInput', (a-200,b), '')
        _output = add_node(group, 'NodeGroupOutput', (a+1400,b), '')
        _sepa_nucleic = sepa_by_attr(group, (a,b), 'POINT', 'sec_struct', 'INT', 0, 'EQUAL')
        _mesh2curve, _curve_res, _resample_curve = mesh_to_curve(group, (a+200,b), 8)
        _curve2mesh_nucleic = add_node(group, 'GeometryNodeCurveToMesh', (a+1200,b), '')
        set_node_values(_curve2mesh_nucleic, [(2,True)])    # Fill Caps True
        _curve_circle = add_node(group, 'GeometryNodeCurvePrimitiveCircle', (a+800, b-200), '')
        set_node_values(_curve_circle, [(0,8),(4,0.3)])
        nodes_link(group, _input, 0, _sepa_nucleic, 0)
        nodes_link(group, _input, 1, _curve_res, 2)
        nodes_link(group, _input, 1, _curve_circle, 0)
        nodes_link(group, _input, 2, _curve_circle, 4)
        nodes_link(group, _sepa_nucleic, 0, _mesh2curve, 0)
        nodes_link(group, _resample_curve, 0, _curve2mesh_nucleic, 0)
        nodes_link(group, _curve_circle, 0, _curve2mesh_nucleic, 1)
        nodes_link(group, _curve2mesh_nucleic, 0, _output, 0)


def Mol3D_Atoms(group, location, _input):
    a,b=location
    _instance = add_node(group, 'GeometryNodeInstanceOnPoints', (a+200,b), '')
    _realize = add_node(group, 'GeometryNodeRealizeInstances', (a+400,b), '')
    _uv_sphere = add_node(group, 'GeometryNodeMeshUVSphere', (a,b-200), '')
    _align_euler2vec = add_node(group, 'FunctionNodeAlignEulerToVector',(a,b-400), '')
    _align_euler2vec.axis = 'Z'
    _attr_atom_orient = add_attr_node(group, (a-200,b-500), '', [(0,'atom_orient')], 'FLOAT_VECTOR')
    _attr_radii = add_attr_node(group, (a-200,b-700), '', [], 'FLOAT')
    _multiply = add_math_node(group, (a,b-600), '', [], 'MULTIPLY')
    _multiply_segs = add_math_node(group, (a-200,b), '', [(1,8)], 'MULTIPLY')
    _multiply_rings = add_math_node(group, (a-200, b-200), '', [(1,5)], 'MULTIPLY')

    _attr_poly_aid = add_attr_node(group, (a-400,b-900), '', [(0,'poly_aid')], 'FLOAT')
    _not_equal = add_compare_node(group, (a-200,b-900),'', [(1,0.0)], 'FLOAT', 'NOT_EQUAL')
    _equal = add_compare_node(group, (a-200,b-1100),'', [(1,1.0)], 'FLOAT', 'EQUAL')
    _switch_1 = add_node(group, 'GeometryNodeSwitch', (a,b-800), '')
    _switch_1.input_type = 'FLOAT'
    set_node_values(_switch_1,[(2,1.0)])
    _switch_2 = add_node(group, 'GeometryNodeSwitch', (a,b-1000), '')
    _switch_2.input_type = 'FLOAT'
    _multiply_2 = add_math_node(group, (a+200,b-600), '', [], 'MULTIPLY')
    nodes_link(group, _attr_poly_aid, 1, _not_equal, 0)
    nodes_link(group, _attr_poly_aid, 1, _equal, 0)
    nodes_link(group, _not_equal, 0, _switch_1, 0)
    nodes_link(group, _equal, 0, _switch_2, 0)
    nodes_link(group, _switch_2, 0, _switch_1, 3)
    nodes_link(group, _switch_1, 0, _multiply_2, 1)

    nodes_link(group, _input, 0, _instance, 0)
    nodes_link(group, _multiply_segs, 0, _uv_sphere, 0)
    nodes_link(group, _multiply_rings, 0, _uv_sphere, 1)
    nodes_link(group, _instance, 0, _realize, 0)
    nodes_link(group, _uv_sphere, 0, _instance, 2)
    nodes_link(group, _attr_atom_orient, 0, _align_euler2vec, 2)
    nodes_link(group, _align_euler2vec, 0, _instance, 5)
    nodes_link(group, _attr_radii, 1, _multiply, 1)
    nodes_link(group, _multiply, 0, _multiply_2, 0)
    nodes_link(group, _multiply_2, 0, _instance, 6)
    return _realize