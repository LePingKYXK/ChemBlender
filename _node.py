import bpy
import bmesh
from .Chem_Data import ELEMENTS_DEFAULT, BONDS_DEFAULT
from . import _node, _mesh, _math
from . import version as v

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
def Ball_Stick_nodes(group, atom_scale, bond_scale, atom_subdiv, bond_subdiv, atom_syms, type, hide_H, 
                     color_type, dash_count, dash_radius, category, display_mode, hide_segs):
    a,b = (0,0)
    _input, _output = set_io_nodes(group, (a,b), (a+1400,b))
    delete_ios(group)

    # hide hydrogen switch
    _switch_H = add_node(group, 'GeometryNodeSwitch', (a+200,b), 'Remove Hydrogens')
    _delete_H = add_node(group, 'GeometryNodeDeleteGeometry', (a+200,b-200), '')
    _attr_atom_id = add_attr_node(group, (a+200,b-600), '', [(0,'atom_id')], 'INT')
    _equal = add_compare_node(group, (a+200,b-400),'', [(3,1)], 'INT', 'EQUAL')
    nodes_link(group, _attr_atom_id, v.na_sockets_out['int'], _equal, 2)
    nodes_link(group, _equal, 0, _delete_H, 1)
    nodes_link(group, _input, 0, _switch_H, v.sw_sockets_in['geoA'])
    nodes_link(group, _input, 0, _delete_H, 0)
    nodes_link(group, _delete_H, 0, _switch_H, v.sw_sockets_in['geoB'])
    if hide_H: set_node_values(_switch_H, [(v.sw_sockets_in['switch_1'],True)])

    # remove solvents
    _switch_sol = add_node(group, 'GeometryNodeSwitch', (a+400,b), 'Remove Solvents')
    _delete_sol = add_node(group, 'GeometryNodeDeleteGeometry', (a+400,b-200), '')
    _attr_is_sol = add_attr_node(group, (a+400,b-400), '', [(0,'is_solvent')], 'BOOLEAN')
    nodes_link(group, _attr_is_sol, v.na_sockets_out['bool'], _delete_sol, 1)
    nodes_link(group, _switch_H, v.sw_sockets_out['geo'], _switch_sol, v.sw_sockets_in['geoA'])
    nodes_link(group, _delete_sol, 0, _switch_sol, v.sw_sockets_in['geoB'])
    nodes_link(group, _switch_H, v.sw_sockets_out['geo'], _delete_sol, 0)
    remove_solvent = False
    if remove_solvent: set_node_values(_switch_sol, [(1,True)])

    # dispaly proteins and/or nucleics/ligands
    _sepa_peptide = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b), 'Peptides')
    _attr_is_peptide = add_attr_node(group, (a+600,b-200), '', [(0,'is_peptide')], 'BOOLEAN')
    nodes_link(group, _attr_is_peptide, v.na_sockets_out['bool'], _sepa_peptide, 1)
    nodes_link(group, _switch_sol, v.sw_sockets_out['geo'], _sepa_peptide, 0)
    Outputs = {'C0':[_sepa_peptide,0], 'C1':[_sepa_peptide,1], 'C2':[_switch_sol,v.sw_sockets_out['geo']]}

    # display backbones and/or residues
    _sepa_backbone = add_node(group, 'GeometryNodeSeparateGeometry', (a+1000,b), 'Backbone')
    _attr_is_alpha = add_attr_node(group, (a+1000,b-300), '', [(0,'is_alpha')], 'BOOLEAN')
    _attr_is_backbone = add_attr_node(group, (a+800, b-300), '', [(0,'is_backbone')], 'BOOLEAN')
    _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (a+1000,b-150), '')
    set_node_operation(_boolean_OR, 'OR')
    _boolean_NOT = add_node(group, 'FunctionNodeBooleanMath', (a+800,b-150), '')
    set_node_operation(_boolean_NOT, 'NOT')
    nodes_link(group, _attr_is_alpha, v.na_sockets_out['bool'], _boolean_OR, 1)
    nodes_link(group, _attr_is_backbone, v.na_sockets_out['bool'], _boolean_NOT, 0)
    nodes_link(group, _boolean_NOT, 0, _boolean_OR, 0)
    socket_out = 0
    bio_style = 'B'
    if bio_style in ['B','F']:
        if display_mode == 'M0':  # only backbone
            nodes_link(group, _attr_is_backbone, v.na_sockets_out['bool'], _sepa_backbone, 1)
        elif display_mode == 'M1':  # only residues
            nodes_link(group, _boolean_OR, 0, _sepa_backbone, 1)
    else:
        if display_mode == 'M0':  # only backbone
            nodes_link(group, _attr_is_backbone, v.na_sockets_out['bool'], _sepa_backbone, 1)
        elif display_mode == 'M1':  # only residues
            socket_out = 1
    nodes_link(group, Outputs[category][0], Outputs[category][1], _sepa_backbone, 0)

    a += 800
    Atom_Radii = {'B':'radii_B', 'S':'radii_S', 'W':'radii_W', 'F':'radii_F'}
    Bond_Radii = {'B':'bond_radii_B', 'S':'bond_radii_S', 'W':'bond_radii_W', 'F':'bond_radii_F'}
    # create Atoms
    _atoms = add_node(group, 'Mol3D_MT_Atoms',(a+400,b),'Atoms')
    set_node_values(_atoms, [(1,Atom_Radii[type]),(2,atom_scale),(3,atom_subdiv)])
    # dashed lines
    _dashed_bonds = add_node(group, 'Mol3D_MT_Dashline', (a+400,b+150), 'Dashed Bonds')
    set_node_values(_dashed_bonds, [(1,dash_radius),(2,dash_count)])
    # create Bonds
    _bonds = add_node(group, 'Mol3D_MT_Bonds', (a+400,b+350), 'Bonds')
    set_node_values(_bonds, [(1,Bond_Radii[type]),(2,bond_scale),(3,bond_subdiv)])
    # separate Bonds
    _sepa_bonds = add_node(group, 'Mol3D_MT_Sepa_Bonds', (a+600,b+350), 'Separate Bonds')
    _join_bonds = add_node(group, 'GeometryNodeJoinGeometry', (a+800,b+350), 'All bonds')
    # aromatic Bonds
    _aromatic = add_node(group, 'Mol3D_MT_Aromatic_Bonds', (a+400,b-200), 'Aromatic')    
    # a+=2000
    _join_aromatic = add_node(group, 'GeometryNodeJoinGeometry', (a+600,b-200), 'Join')
    # join
    _join_all = add_node(group, 'GeometryNodeJoinGeometry', (a+1000,b), 'All Struct')
    nodes_link(group, _sepa_backbone, socket_out, _atoms, 0)
    nodes_link(group, _sepa_backbone, socket_out, _dashed_bonds, 0)
    nodes_link(group, _sepa_backbone, socket_out, _bonds, 0)
    nodes_link(group, _sepa_backbone, socket_out, _aromatic, 0)
    nodes_link(group, _bonds, 0, _sepa_bonds, 0)
    nodes_link(group, _sepa_bonds, 2, _join_bonds, 0)
    nodes_link(group, _sepa_bonds, 1, _join_bonds, 0)
    nodes_link(group, _sepa_bonds, 0, _join_bonds, 0)
    nodes_link(group, _join_aromatic, 0, _join_all, 0)
    nodes_link(group, _atoms, 0, _join_all, 0)
    nodes_link(group, _dashed_bonds, 0, _join_all, 0)
    if type in ['B','W']:
        nodes_link(group, _join_bonds, 0, _join_all, 0)
    else:
        nodes_link(group, _bonds, 0, _join_all, 0)
    
    # set material and smooth
    _store_attr_sharp_face = add_store_attr_node(group, (a+1200,b), '', 'BOOLEAN', 'FACE', 'sharp_face')
    _edge_angle = add_node(group, 'GeometryNodeInputMeshEdgeAngle', (a+1200,b-400), '')
    _greater = add_compare_node(group, (a+1200,b-200), '', [(1,1.500)], 'FLOAT', 'GREATER_THAN')
    nodes_link(group, _join_all, 0, _store_attr_sharp_face, 0)
    nodes_link(group, _edge_angle, 0, _greater, 0)
    nodes_link(group, _greater, 0, _store_attr_sharp_face, v.sa_sockets_in['bool'])

    mat = bpy.data.materials['Hetero']
    _hetero = add_mat_node(group, (a+1600,b+200), 'Hetero', mat)
    _attr_color = add_attr_node(group, (a+1400,b+200), '', [(0,'colour')], 'FLOAT_COLOR')
    for i,element in enumerate(atom_syms):
        _mat_element = add_mat_node(group, (a+1600,b-150*i), f"Mat_{'%03d' % i}", bpy.data.materials[element])
        _input_color = add_node(group, 'FunctionNodeInputColor', (a+1400,b), '')
        try:
            _input_color.color = ELEMENTS_DEFAULT[element][3]
        except:
            _input_color.value = ELEMENTS_DEFAULT[element][3]
        _equal = add_compare_node(group, (a+1400,b), '', [], 'RGBA', 'EQUAL')
        nodes_link(group, _attr_color, v.na_sockets_out['color'], _equal, 6)
        nodes_link(group, _input_color, 0, _equal, 7)
        nodes_link(group, _equal, 0, _mat_element, 1)
    for i in range(len(atom_syms)-1):
        nodes_link(group, get_node(group, f"Mat_{'%03d' % i}"), 0, get_node(group, f"Mat_{'%03d' % (i+1)}"), 0)

    _last_mat = get_node(group, f"Mat_{'%03d' % (len(atom_syms)-1)}")
    _output.location = (a+1800,b)
    nodes_link(group, _store_attr_sharp_face, 0, _hetero, 0)
    nodes_link(group, _store_attr_sharp_face, 0, get_node(group, 'Mat_000'), 0)
    if color_type == 'H':
        nodes_link(group, _hetero, 0, _output, 0)
    else:
        nodes_link(group, _last_mat, 0, _output, 0)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# create a new GeometryNode modifier for given object
def create_GeoNode(obj, GN_name, group_name, ):
    bpy.context.view_layer.objects.active = obj
    geonode = obj.modifiers.new(GN_name, 'NODES')
    bpy.ops.node.new_geometry_node_group_assign()
    geonode.node_group.name = group_name
    return geonode

# set locations of Group Input and Group Output nodes
def set_io_nodes(group, in_location, out_location):
    try:
        _input = group.nodes['Group Input']
    except:
        _input = group.nodes['组输入']
    _input.location = in_location
    try:
        _output = group.nodes['Group Output']
    except:
        _output = group.nodes['组输出']
    _output.location = out_location
    _input.label = 'in_node'
    _output.label = 'out_node'
    return (_input, _output)

# add a common geometry node in given node_group
def add_node(group, name, location, label):
    node = group.nodes.new(name)
    node.location = location
    node.label = label
    return node

# find the node in node_group through node label
def get_node(group, node_label):
    target_node = []
    for node in group.nodes:
        if node.label == node_label:
            target_node = node
    return target_node

# create links between two nodes
def nodes_link(group, node_a, socket_out, node_b, socket_in):
    link = group.links.new
    link(node_a.outputs[socket_out], node_b.inputs[socket_in])


# ----------------------------------------------------------------------------
# -------------------------- add single node ---------------------------------
# set socket value of the node inputs
def set_node_values(node, socket_values):
    for socket_value in socket_values:
        node.inputs[socket_value[0]].default_value = socket_value[1]

def set_node_datatype(node, datatype):
    node.data_type = datatype

def set_node_operation(node, operation):
    node.operation = operation

# add a Set Material node
def add_mat_node(group, location, label, material):
    node = add_node(group, 'GeometryNodeSetMaterial', location, label)
    node.inputs[2].default_value = material
    return node

# add a Compare node
def add_compare_node(group, location, label, socket_values, datatype, operation):
    node = add_node(group, 'FunctionNodeCompare', location, label)
    set_node_values(node, socket_values)
    set_node_datatype(node, datatype)
    set_node_operation(node, operation)
    return node

# add a Named Attribute node
def add_attr_node(group, location, label, socket_values, datatype):
    node = add_node(group, 'GeometryNodeInputNamedAttribute', location, label)
    set_node_values(node, socket_values)
    set_node_datatype(node, datatype)
    return node

# add a Store Named Attribute node
def add_store_attr_node(group, location, label, data_type, domain, attr_name):
    node = add_node(group, 'GeometryNodeStoreNamedAttribute', location, label)
    node.data_type = data_type
    node.domain = domain
    set_node_values(node, [(2,attr_name)])
    return node

# add a Sample Index node
def add_sample_index_node(group, location, label, data_type, domain):
    node = add_node(group, 'GeometryNodeSampleIndex', location, label)
    node.data_type = data_type
    node.domain = domain
    return node

# add a Math node
def add_math_node(group, location, label, socket_values, operation):
    node = add_node(group, 'ShaderNodeMath', location, label)
    set_node_values(node, socket_values)
    set_node_operation(node, operation)
    return node

# add a Vector Math node
def add_vec_math_node(group, location, label, socket_values, operation):
    node = add_node(group, 'ShaderNodeVectorMath', location, label)
    set_node_values(node, socket_values)
    set_node_operation(node, operation)
    return node

# for removing free atoms
def remove_solvents(group, _input, location):
    link = group.links.new
    a,b = location
    _named_attr = add_attr_node(group, (a-200,b-100), '',[(0,'is_solvent')], 'BOOLEAN')
    _delete_geometry = add_node(group, 'GeometryNodeDeleteGeometry', location, 'Remove')

    nodes_link(group, _input, 0, _delete_geometry, 0)
    nodes_link(group, _named_attr, v.na_sockets_out['bool'], _delete_geometry, 1)

    return _delete_geometry
    # -----------------------

# 
def display(group, location, _input, bio_style, category, display_mode):
    a,b = location
    _remove = remove_solvents(group, _input, (a, b))
    _sepa_geo = add_node(group, 'GeometryNodeSeparateGeometry', (a+200,b), '')
    attr_name = 'is_peptide' if category=='C0' else 'is_nucleic'
    _attr = add_attr_node(group, (a,b-200), '', [(0,attr_name)], 'BOOLEAN')
    nodes_link(group, _remove, 0, _sepa_geo, 0)
    nodes_link(group, _attr, v.na_sockets_out['bool'], _sepa_geo, 1)
    _remove = _remove if category=='C2' else _sepa_geo
    # hide_backbone
    a += 400
    _sepa_backbone, socket_out = add_hide_backbone_nodes(group, (a,b), bio_style, display_mode, _remove)
    _attr_res_id = add_attr_node(group, (a+200,b+150), '', [(0,'res_id')], 'INT')
    _index = add_store_attr_node(group, (a+200,b), '', 'INT', 'POINT', 'index')
    nodes_link(group, _sepa_backbone, socket_out, _index, 0)
    nodes_link(group, _attr_res_id, v.na_sockets_out['int'], _index, v.sa_sockets_in['int'])
    a += 200
    return (_index, a)


# delete all old nodes except _input and _output
def delete_ios(group): 
    try:
        for node in group.nodes:
            if node.label not in ['in_node', 'out_node']:
                group.nodes.remove(node)
    except Exception as e:
        print(e)


# separate geometry through named attribute
def sepa_by_attr(group, location, sepa_domain, attr_name, datatype, value, operation):
    link = group.links.new
    a,b = location
    _sepa_geo = add_node(group, 'GeometryNodeSeparateGeometry', location, '')
    _sepa_geo.domain = sepa_domain
    _attr = add_attr_node(group, (a,b-400), '', [(0,attr_name)], datatype)
    type = {
       'VECTOR':'vec',
       'FLOAT':'float',
       'COLOR':'color',
       'BOOLEAN':'bool',
       'INT':'int',
       'QUATERNION':'quater', 
    }
    _equal = add_compare_node(group, (a,b-200),'', [(v.cp_sockets_in[type[datatype]+'B'],value)], datatype, operation)
    link(_attr.outputs[v.na_sockets_out[type[datatype]]], _equal.inputs[v.cp_sockets_in[type[datatype]+'A']])
    link(_equal.outputs[0], _sepa_geo.inputs[1])
    return _sepa_geo

# select smooth
def sharp_shade(group, location, domain, critical):
    link = group.links.new
    a,b = location
    _store_attr = add_node(group, 'GeometryNodeStoreNamedAttribute', location, '')
    _edge_angle = add_node(group, 'GeometryNodeInputMeshEdgeAngle', (a,b-400), '')
    _store_attr.data_type = 'BOOLEAN'
    _store_attr.domain = domain
    sharp = 'sharp_face' if domain=='FACE' else 'sharp_edge'
    set_node_values(_store_attr, [(2,sharp)])
    _compare = add_compare_node(group, (a,b-200), 'Greater', [(1,critical)], 'FLOAT', 'GREATER_THAN')   # greater than critical
    link(_edge_angle.outputs[0], _compare.inputs[0])
    link(_compare.outputs[0], _store_attr.inputs[6])
    return _store_attr


def add_attr_trans_nodes(group, attr_name, data_type, domain, location, label):
    a,b = location
    # domain = 'POINT'
    _store_attr = add_store_attr_node(group, location, label, data_type, domain, attr_name)
    _sample_idx = add_sample_index_node(group, (a,b+240), '', data_type, domain)
    _input_idx =  add_node(group, 'GeometryNodeInputIndex', (a,b+300), '')
    _input_attr = add_node(group, 'GeometryNodeInputNamedAttribute', (a,b+450), '')
    _input_attr.data_type = data_type
    set_node_values(_input_attr, [(0,attr_name)])
    type = {
       'FLOAT_VECTOR':'vec',
       'FLOAT':'float',
       'FLOAT_COLOR':'color',
       'BOOLEAN':'bool',
       'INT':'int',
       'QUATERNION':'quater', 
    }
    nodes_link(group, _input_attr, v.na_sockets_out[type[data_type]], _sample_idx, v.si_sockets_in[type[data_type]])
    nodes_link(group, _input_idx, 0, _sample_idx, v.si_sockets_in['index'])
    nodes_link(group, _sample_idx, v.si_sockets_out[type[data_type]], _store_attr, v.sa_sockets_in[type[data_type]])
    return (_sample_idx, _store_attr)



# ----------------------------------------------------------------------------
def add_guide_vector_nodes(group, location, is_nucleic):
    a,b = location
    link = group.links.new
    _store_attr_guide_X = add_store_attr_node(group, location, '', 'FLOAT_VECTOR', 'POINT', 'guide_X')
    _store_attr_guide_Z = add_store_attr_node(group, (a+200,b), '', 'FLOAT_VECTOR', 'POINT', 'guide_Z')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-1000,b-500), '')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a-1000, b-600), '')
    if not is_nucleic:
        ids = [-1,1,2,1]
    else:
        ids = [-3,-5,1,0]
        _add_0 = add_math_node(group, (a-1200,b-600), '', [(1,1)], 'ADD')
        _field_idx_0 = add_node(group, 'GeometryNodeFieldAtIndex', (a-1200,b-800), '')
        set_node_datatype(_field_idx_0, 'INT')
        _equal = add_compare_node(group, (a-1200,b-1000), '', [(3,62)], 'INT', 'EQUAL')
        _attr_atom_name = add_attr_node(group, (a-1400,b-800), '', [(0,'atom_name')], 'INT')
        _switch = add_node(group, 'GeometryNodeSwitch', (a-1000,b-750), '')
        _switch.input_type = 'INT'
        set_node_values(_switch, [(4,7),(5,4)])
        link(_input_idx.outputs[0], _add_0.inputs[0])
        link(_add_0.outputs[0], _field_idx_0.inputs[0])
        link(_attr_atom_name.outputs[4], _field_idx_0.inputs[2])
        link(_field_idx_0.outputs[1], _equal.inputs[2])
        link(_equal.outputs[0], _switch.inputs[0])
    
    _subtract1 = add_vec_math_node(group, (a-400,b-200), '', [], 'SUBTRACT')
    _subtract2 = add_vec_math_node(group, (a-400,b-600), '', [], 'SUBTRACT')
    _normalize1 = add_vec_math_node(group, (a-200,b-200), '', [], 'NORMALIZE')
    _normalize2 = add_vec_math_node(group, (a-200,b-600), '', [], 'NORMALIZE')
    for i in range(4):
        _add = add_math_node(group, (a-800,b-200*(i+1)), '', [(1,ids[i])], 'ADD')
        _field_idx = add_node(group, 'GeometryNodeFieldAtIndex', (a-600,b-200*(i+1)), '')
        set_node_datatype(_field_idx, 'FLOAT_VECTOR')
        link(_input_idx.outputs[0], _add.inputs[0])
        link(_add.outputs[0], _field_idx.inputs[0])
        nodes_link(group, _input_pos, 0, _field_idx, v.fi_sockets_in['vec'])
        if i < 2:
            nodes_link(group, _field_idx, v.fi_sockets_out['vec'], _subtract1, i)
        else:
            nodes_link(group, _field_idx, v.fi_sockets_out['vec'], _subtract2, i-2)
        if is_nucleic and i==2: nodes_link(group, _switch, v.sw_sockets_out['switch_1'], _add, 1)

    link(_subtract1.outputs[0], _normalize1.inputs[0])
    link(_subtract2.outputs[0], _normalize2.inputs[0])
    nodes_link(group, _normalize1, 0, _store_attr_guide_X, v.sa_sockets_in['vec'])
    nodes_link(group, _normalize2, 0, _store_attr_guide_Z, v.sa_sockets_in['vec'])
    link(_store_attr_guide_X.outputs[0], _store_attr_guide_Z.inputs[0])

    return (_store_attr_guide_X, _store_attr_guide_Z)


def add_guide_Z_inverse_nodes(group, location,):
    a,b = location
    link = group.links.new
    _store_attr_guide_Z = add_store_attr_node(group, location, '', 'FLOAT_VECTOR', 'POINT', 'guide_Z')
    _attr_sec_struct = add_attr_node(group, (a-400,b-200), '', [(0,'sec_struct')], 'INT')
    _equal = add_compare_node(group, (a-200,b-200), '', [(3,2)], 'INT', 'EQUAL')
    _switch = add_node(group, 'GeometryNodeSwitch', (a,b-200), '')
    _switch.input_type = 'VECTOR'
    _attr_guide_Z = add_attr_node(group, (a-400,b-350), '', [(0,'guide_Z')], 'FLOAT_VECTOR')
    _scale = add_vec_math_node(group, (a-200,b-400), '', [], 'SCALE')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-600,b-350), '')
    _modulo = add_math_node(group, (a-600,b-500), '', [(1,2)], 'FLOORED_MODULO')
    _subtract = add_math_node(group, (a-400,b-500), '', [(1,0.5)], 'SUBTRACT')
    _multiply = add_math_node(group, (a-200,b-500), '', [(1,2)], 'MULTIPLY')
    
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _equal, v.cp_sockets_in['intA'])
    nodes_link(group, _attr_guide_Z, 0, _switch, v.sw_sockets_in['vecA'])
    nodes_link(group, _attr_guide_Z, 0, _scale, 0)
    nodes_link(group, _scale, 0, _switch, v.sw_sockets_in['vecB'])
    nodes_link(group, _equal, 0, _switch, 0)
    nodes_link(group, _input_idx, 0, _modulo, 0)
    nodes_link(group, _modulo, 0, _subtract, 0)
    nodes_link(group, _subtract, 0, _multiply, 0)
    nodes_link(group, _multiply, 0, _scale, 3)
    nodes_link(group, _switch, v.sw_sockets_out['vec'], _store_attr_guide_Z, v.sa_sockets_in['vec'])
    return _store_attr_guide_Z

def guide_Z_revised(group, location):
    a,b = location
    _store_attr_guide_Z = add_store_attr_node(group, location, '', 'FLOAT_VECTOR', 'POINT', 'guide_Z')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-800,b+200), '')
    _attr_sec_struct = add_attr_node(group, (a-600,b+200), '', [(0,'sec_struct')], 'INT')
    _equal_1 = add_compare_node(group, (a-400,b+200), '', [(3,1)], 'INT', 'EQUAL')
    _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a-200,b+200), '')
    set_node_operation(_boolean_AND, 'AND')
    _add = add_math_node(group, (a-600,b+400), '', [(1,1)], 'ADD')
    _field_idx_1 = add_node(group, 'GeometryNodeFieldAtIndex', (a-400,b+400), '')
    set_node_datatype(_field_idx_1, 'INT')
    _equal_2 = add_compare_node(group, (a-200,b+400), '', [(3,3)], 'INT', 'EQUAL')
    _subtract = add_math_node(group, (a-600,b), '', [(1,1)], 'SUBTRACT')
    _field_idx_2 = add_node(group, 'GeometryNodeFieldAtIndex', (a-400,b), '')
    set_node_datatype(_field_idx_2, 'FLOAT_VECTOR')
    _switch = add_node(group, 'GeometryNodeSwitch', (a-200,b), '')
    _switch.input_type = 'VECTOR'
    _attr_guide_Z = add_attr_node(group, (a-600,b-200), '', [(0,'guide_Z')], 'FLOAT_VECTOR')

    nodes_link(group, _input_idx, 0, _add, 0)
    nodes_link(group, _input_idx, 0, _subtract, 0)
    nodes_link(group, _add, 0, _field_idx_1, 0)
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _field_idx_1, v.fi_sockets_in['int'])
    nodes_link(group, _field_idx_1, v.fi_sockets_out['int'], _equal_2, v.cp_sockets_in['intA'])
    nodes_link(group, _equal_2, 0, _boolean_AND, 0)
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _equal_1, v.cp_sockets_in['intA'])
    nodes_link(group, _equal_1, 0, _boolean_AND, 1)
    nodes_link(group, _subtract, 0, _field_idx_2, 0)
    nodes_link(group, _attr_guide_Z, 0, _field_idx_2, v.fi_sockets_in['vec'])
    nodes_link(group, _field_idx_2, v.fi_sockets_out['vec'], _switch, v.sw_sockets_in['vecB'])
    nodes_link(group, _attr_guide_Z, 0, _switch, v.sw_sockets_in['vecA'])
    nodes_link(group, _boolean_AND, 0, _switch, 0)
    nodes_link(group, _switch, v.sw_sockets_out['vec'], _store_attr_guide_Z, v.sa_sockets_in['vec'])

    return _store_attr_guide_Z

def angle_XZ(group, location):
    a,b = location
    _store_attr_angle_XZ = add_store_attr_node(group, location, '', 'FLOAT', 'POINT', 'angle_XZ')
    _attr_guide_X = add_attr_node(group, (a-800,b), '', [(0,'guide_X')], 'FLOAT_VECTOR')
    _attr_guide_Z = add_attr_node(group, (a-800,b-200), '', [(0,'guide_Z')], 'FLOAT_VECTOR')
    _dot_product = add_vec_math_node(group, (a-600,b), '', [], 'DOT_PRODUCT')
    _arccos = add_math_node(group, (a-400,b), '', [], 'ARCCOSINE')
    _rad2deg = add_math_node(group, (a-200,b), '', [], 'DEGREES')
    _subtract = add_math_node(group, (a-200,b-200), '', [(1,109.47)], 'SUBTRACT')

    nodes_link(group, _attr_guide_X, 0, _dot_product, 0)
    nodes_link(group, _attr_guide_Z, 0, _dot_product, 1)
    nodes_link(group, _dot_product, 1, _arccos, 0)
    nodes_link(group, _arccos, 0, _rad2deg, 0)
    nodes_link(group, _rad2deg, 0, _subtract, 0)
    nodes_link(group, _subtract, 0, _store_attr_angle_XZ, v.sa_sockets_in['float'])

    return _store_attr_angle_XZ



# attribute 'scale_beta' is used to determine the arrow shape in beta-sheets  
def add_scale_beta_nodes(group, location):
    a,b = location
    link = group.links.new
    _store_attr_scale_beta = add_store_attr_node(group, location, '', 'FLOAT', 'POINT', 'scale_beta')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-800,b+300), '')
    _add = add_math_node(group, (a-600,b+400), '', [(1,1)], 'ADD')
    _attr_sec_struct = add_attr_node(group, (a-600,b+200), '', [(0,'sec_struct')], 'INT')
    _field_idx = add_node(group, 'GeometryNodeFieldAtIndex', (a-400,b+400), '')
    set_node_datatype(_field_idx, 'INT')
    _equal = add_compare_node(group, (a-400,b+200), '', [(3,2)], 'INT', 'EQUAL')
    _not_equal = add_compare_node(group, (a-200,b+400), '', [(3,2)], 'INT', 'NOT_EQUAL')
    _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a-200,b+200), '')
    set_node_operation(_boolean_AND, 'AND')
    _boolean_NOT = add_node(group, 'FunctionNodeBooleanMath', (a,b+200), '')
    set_node_operation(_boolean_NOT, 'NOT')

    link(_input_idx.outputs[0], _add.inputs[0])
    link(_add.outputs[0], _field_idx.inputs[0])
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _field_idx, v.fi_sockets_in['int'])
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _equal, v.cp_sockets_in['intA'])
    nodes_link(group, _field_idx, v.fi_sockets_out['int'], _not_equal, v.cp_sockets_in['intA'])
    link(_not_equal.outputs[0], _boolean_AND.inputs[0])
    link(_equal.outputs[0], _boolean_AND.inputs[1])
    link(_boolean_AND.outputs[0], _boolean_NOT.inputs[0])
    nodes_link(group, _boolean_NOT, 0, _store_attr_scale_beta, v.sa_sockets_in['float'])
    return _store_attr_scale_beta


# attribute 'coil_end' is used to form completed loop
def add_coil_end_nodes(group, location):
    a,b = location
    link = group.links.new
    _store_attr_coil_end = add_store_attr_node(group, location, '', 'FLOAT', 'POINT', 'coil_end')

    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-1000,b+400), '')
    _add = add_math_node(group, (a-800,b+600), '', [(1,1)], 'ADD')
    _attr_sec_struct = add_attr_node(group, (a-800,b+400), '', [(0,'sec_struct')], 'INT')
    _subtract = add_math_node(group, (a-800,b+200), '', [(1,-1)], 'ADD')
    _field_idx1 = add_node(group, 'GeometryNodeFieldAtIndex', (a-600,b+600), '')
    set_node_datatype(_field_idx1, 'INT')
    _not_equal = add_compare_node(group, (a-600,b+400), '', [(3,3)], 'INT', 'NOT_EQUAL')
    _field_idx2 = add_node(group, 'GeometryNodeFieldAtIndex', (a-600,b+200), '')
    set_node_datatype(_field_idx2, 'INT')
    _equal_1 = add_compare_node(group, (a-400,b+600), '', [(3,3)], 'INT', 'EQUAL')
    _equal_2 = add_compare_node(group, (a-400,b+200), '', [(3,3)], 'INT', 'EQUAL')
    _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a-400,b+400), '')
    set_node_operation(_boolean_AND, 'AND')
    _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (a-200,b+200), '')
    set_node_operation(_boolean_OR, 'OR')

    link(_input_idx.outputs[0], _add.inputs[0])
    link(_input_idx.outputs[0], _subtract.inputs[0])
    link(_add.outputs[0], _field_idx1.inputs[0])
    link(_subtract.outputs[0], _field_idx2.inputs[0])
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _field_idx1, v.fi_sockets_in['int'])
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _field_idx2, v.fi_sockets_in['int'])
    nodes_link(group, _attr_sec_struct, v.na_sockets_out['int'], _not_equal, v.cp_sockets_in['intA'])
    nodes_link(group, _field_idx1, v.fi_sockets_out['int'], _equal_1, v.cp_sockets_in['intA'])
    nodes_link(group, _field_idx2, v.fi_sockets_out['int'], _equal_2, v.cp_sockets_in['intA']) 
    link(_not_equal.outputs[0], _boolean_AND.inputs[0])
    link(_equal_1.outputs[0], _boolean_AND.inputs[1])
    link(_boolean_AND.outputs[0], _boolean_OR.inputs[0])
    link(_equal_2.outputs[0], _boolean_OR.inputs[1])
    nodes_link(group, _boolean_OR, 0, _store_attr_coil_end, v.sa_sockets_in['float'])
    return _store_attr_coil_end


def delete_mismatch(group, location, max_length):
    a,b = location
    link = group.links.new
    _delete_geometry = add_node(group, 'GeometryNodeDeleteGeometry', location, '')
    _delete_geometry.domain ='EDGE'
    _edge_verts = add_node(group, 'GeometryNodeInputMeshEdgeVertices', (a-600,b-200), '')
    _distance = add_vec_math_node(group, (a-400,b-200), '', [], 'DISTANCE')
    # 6.5 here is the critical distance
    _greater_than = add_compare_node(group, (a-200,b-200), '', [(1,max_length)], 'FLOAT', 'GREATER_THAN')

    link(_edge_verts.outputs[2], _distance.inputs[0])
    link(_edge_verts.outputs[3], _distance.inputs[1])
    link(_distance.outputs[1], _greater_than.inputs[0])
    link(_greater_than.outputs[0], _delete_geometry.inputs[1])
    return (_delete_geometry, _greater_than)


def mesh_to_curve(group, location, resolution):
    link = group.links.new
    a,b = location
    # convert mesh to curve
    _mesh2curve = add_node(group, 'GeometryNodeMeshToCurve', (a, b), '')
    # set Bezier type
    _curve_type = add_node(group, 'GeometryNodeCurveSplineType', (a+200, b), '')
    _curve_type.spline_type = 'BEZIER'
    # set curve handles
    _curve_handles = add_node(group, 'GeometryNodeCurveSetHandles', (a+400, b), '')
    # set curve resolutionsocket_out
    _curve_res = add_node(group, 'GeometryNodeSetSplineResolution', (a+600, b), 'spline_res')
    set_node_values(_curve_res, [(2,resolution)])   # spline Resolution
    # resample curve evaluated
    _resample_curve = add_node(group, 'GeometryNodeResampleCurve', (a+800, b), '')
    _resample_curve.mode = 'EVALUATED'
   
    link(_mesh2curve.outputs[0], _curve_type.inputs[0])
    link(_curve_type.outputs[0], _curve_handles.inputs[0])
    link(_curve_handles.outputs[0], _curve_res.inputs[0])
    link(_curve_res.outputs[0], _resample_curve.inputs[0])
    return (_mesh2curve, _curve_res, _resample_curve)


def separate_geometry(group, sec_struct_name, domain, attr_name, sec_struct_id, location):
    link = group.links.new
    a,b = location
    _sepa_sec_struct = add_node(group, 'GeometryNodeSeparateGeometry', (a,b), sec_struct_name)
    _sepa_sec_struct.domain = domain
    _attr_sec_struct = add_attr_node(group, (a, b-400), '', [(0, attr_name)], 'INT')
    _compare = add_compare_node(group, (a,b-200), '', [(3, sec_struct_id)], 'INT', 'EQUAL')
    link(_attr_sec_struct.outputs[4], _compare.inputs[2])
    link(_compare.outputs[0], _sepa_sec_struct.inputs[1])
    return _sepa_sec_struct

def add_hide_segs_nodes(group, attr_name, seg_idx_min, seg_idx_max, location):
    link = group.links.new
    a,b = location
    _sepa_geometry = add_node(group, 'GeometryNodeSeparateGeometry', location, '')
    _sepa_geometry.domain = 'EDGE'
    _attr = add_attr_node(group, (a, b-800), '', [(0,attr_name)], 'INT')
    _greater_equal = add_compare_node(group, (a,b-600), '', [(3, seg_idx_min)], 'INT', 'GREATER_EQUAL')
    _less_equal = add_compare_node(group, (a,b-400), '', [(3, seg_idx_max)], 'INT', 'LESS_EQUAL')
    _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a,b-200), '')
    set_node_operation(_boolean_AND, 'AND')
    link(_attr.outputs[4], _greater_equal.inputs[2])
    link(_attr.outputs[4], _less_equal.inputs[2])
    link(_greater_equal.outputs[0], _boolean_AND.inputs[0])
    link(_less_equal.outputs[0], _boolean_AND.inputs[1])
    link(_boolean_AND.outputs[0], _sepa_geometry.inputs[1])
    return _sepa_geometry

def hide_segs_by_input_texts(group, location, hide_segs, _input):
    link = group.links.new
    a,b = location
    _geometry = _input
    hide_list = hide_segs.split(',')
    socket_out = 0
    for str in hide_list:
        segs = str[1:].split('-') if str.find('~') == 0 else str.split('-')
        if len(segs) == 1:
            try:
                seg_idx = int(segs[0])
                _sepa_geometry = separate_geometry(group, '', 'EDGE','res_id', seg_idx, (a,b))
                a += 200
                link(_geometry.outputs[socket_out], _sepa_geometry.inputs[0])
                socket_out = 0 if str.find('~') == 0 else 1
                _geometry = _sepa_geometry
            except Exception as e:
                print(e)
        elif len(segs) == 2:
            try:
                seg_idx_min = int(segs[0])
                seg_idx_max = int(segs[1])
                if seg_idx_min > seg_idx_max:
                    min = seg_idx_max
                    seg_idx_max = seg_idx_min
                    seg_idx_min = min
                _sepa_geometry = add_hide_segs_nodes(group, 'res_id', seg_idx_min, seg_idx_max, (a,b))
                a += 200
                link(_geometry.outputs[socket_out], _sepa_geometry.inputs[0])
                socket_out = 0 if str.find('~') == 0 else 1
                _geometry = _sepa_geometry
            except Exception as e:
                print(e)
    return (_geometry, socket_out, a)


def add_hide_backbone_nodes(group, location, bio_style, display_mode, _input):
    link = group.links.new
    a,b = location
    _sepa_backbone = add_node(group, 'GeometryNodeSeparateGeometry', location, 'backbone')
    _attr_is_alpha = add_attr_node(group, (a, b-400), '', [(0,'is_alpha')], 'BOOLEAN')
    _attr_is_backbone = add_attr_node(group, (a-200, b-400), '', [(0,'is_backbone')], 'BOOLEAN')
    _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (a,b-200), '')
    set_node_operation(_boolean_OR, 'OR')
    _boolean_NOT = add_node(group, 'FunctionNodeBooleanMath', (a-200,b-200), '')
    set_node_operation(_boolean_NOT, 'NOT')
    nodes_link(group, _attr_is_alpha, v.na_sockets_out['bool'], _boolean_OR, 1)
    nodes_link(group, _attr_is_backbone, v.na_sockets_out['bool'], _boolean_NOT, 0)
    nodes_link(group, _boolean_NOT, 0, _boolean_OR, 0)
    socket_out = 0
    if bio_style in ['B','F']:
        if display_mode == 'M0':  # only backbone
            nodes_link(group, _attr_is_backbone, v.na_sockets_out['bool'], _sepa_backbone, 1)
        elif display_mode == 'M1':  # only residues
            nodes_link(group, _boolean_OR, 0, _sepa_backbone, 1)
    else:
        if display_mode == 'M0':  # only backbone
            nodes_link(group, _attr_is_backbone, v.na_sockets_out['bool'], _sepa_backbone, 1)
        elif display_mode == 'M1':  # only residues
            socket_out = 1

    link(_input.outputs[0], _sepa_backbone.inputs[0])
    return (_sepa_backbone, socket_out)


# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# crystal scaffold growing direct
def expand_scaffold(scaffold, cell_edges, edge_vecs, grow_iter, grow_fract):
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    if bm.edges:
        domain = 'EDGE'
    else:
        domain = 'POINT'

    GN_expand = create_GeoNode(scaffold, 'GN_expand', 'Ex_scaffold')
    group = GN_expand.node_group
    link = group.links.new
    a,b = (0,0)
    _input, _output = set_io_nodes(group, (a-400,b), (a+2400,b))

    va, vb, vc = edge_vecs
    neighbors = [va,-va,vb,-vb,vc,-vc,
                 va+vb,va-vb,-va+vb,-va-vb,

                 va+vb+vc,vb+vc,-va+vb+vc,
                 va+vc,-va+vc,
                 va-vb+vc,-vb+vc,-va-vb+vc,

                 va+vb-vc,vb-vc,-va+vb-vc,
                 va-vc,-va-vc,
                 va-vb-vc,-vb-vc,-va-vb-vc]

    _join_neighbors = add_node(group, 'GeometryNodeJoinGeometry', (a-200,b-200), '')
    for i,neighbor in enumerate(neighbors):
        _transform = add_node(group, 'GeometryNodeTransform', (a-400,b-400*(i+1)), '')
        set_node_values(_transform, [(1,neighbor)])
        link(_input.outputs[0], _transform.inputs[0])
        link(_transform.outputs[0], _join_neighbors.inputs[0])
    link(_join_neighbors.outputs[0], _output.inputs[0])

    # grow by fract
    _join_geo = add_node(group, 'GeometryNodeJoinGeometry', (a,b), '')
    _sepa_ex = add_node(group, 'GeometryNodeSeparateGeometry', (a,b-200), '')
    _sepa_ex.domain = domain
    _store_attr_0 = add_store_attr_node(group, (a+200,b), '', 'FLOAT', 'POINT', 'origin')
    set_node_values(_store_attr_0, [(v.sa_sockets_in['float'],1.0)])
    _store_attr_1 = add_store_attr_node(group, (a+200,b-200), '', 'FLOAT', 'POINT', 'origin')
    set_node_values(_store_attr_1, [(v.sa_sockets_in['float'],2.0)])
    nodes_link(group, _input, 0, _join_geo, 0)
    nodes_link(group, _sepa_ex, 1, _join_geo, 0)
    nodes_link(group, _join_neighbors, 0, _sepa_ex, 0)
    nodes_link(group, _join_geo, 0, _store_attr_0, 0)
    nodes_link(group, _sepa_ex, 0, _store_attr_1, 0)

    _obj_cell_edges = add_node(group, 'GeometryNodeObjectInfo', (a-1600,b+200), '')
    set_node_values(_obj_cell_edges,[(0,cell_edges)])
    _bound_box = add_node(group, 'GeometryNodeBoundBox', (a-1400,b+200), '')
    _subtract_1 = add_vec_math_node(group, (a-1200, b+400), '', [], 'SUBTRACT')
    _scale_1 = add_vec_math_node(group, (a-1200, b+200), '', [], 'SCALE')
    set_node_values(_scale_1,[(3,0.5)])
    _scale_2 = add_vec_math_node(group, (a-1200, b), '', [], 'SCALE')
    set_node_values(_scale_2,[(3,-1)])
    _subtract_2 = add_vec_math_node(group, (a-800, b+200), '', [], 'SUBTRACT')
    _transform_1 = add_node(group, 'GeometryNodeTransform', (a-1000,b), '')
    _transform_2 = add_node(group, 'GeometryNodeTransform', (a-800,b), '')
    set_node_values(_transform_2,[(3,[grow_fract*2+1]*3)])
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a-1000, b+400), '')
    _raycast = add_node(group, 'GeometryNodeRaycast', (a-600,b+600), '')
    _dot_product = add_vec_math_node(group, (a-400, b+200), '', [], 'DOT_PRODUCT')
    _less = add_compare_node(group, (a-200,b+200), '', [], 'FLOAT', 'LESS_THAN')
    nodes_link(group, _obj_cell_edges, 3, _bound_box, 0)
    nodes_link(group, _bound_box, 0, _transform_1, 0)
    nodes_link(group, _transform_1, 0, _transform_2, 0)
    nodes_link(group, _transform_2, 0, _raycast, 0)
    nodes_link(group, _bound_box, 1, _subtract_1, 1)
    nodes_link(group, _bound_box, 2, _subtract_1, 0)
    nodes_link(group, _subtract_1, 0, _scale_1, 0)
    nodes_link(group, _scale_1, 0, _scale_2, 0)
    nodes_link(group, _scale_2, 0, _transform_1, 1)
    nodes_link(group, _scale_1, 0, _transform_2, 1)
    nodes_link(group, _scale_1, 0, _subtract_2, 0)
    nodes_link(group, _input_pos, 0, _subtract_2, 1)
    nodes_link(group, _subtract_2, 0, _raycast, v.rc_sockets_in['dir'])
    nodes_link(group, _subtract_2, 0, _dot_product, 1)
    nodes_link(group, _raycast, 2, _dot_product, 0)
    nodes_link(group, _dot_product, 1, _less, 0)
    nodes_link(group, _less, 0, _sepa_ex, 1)

    # grow direct
    _repeat_input = add_node(group, 'GeometryNodeRepeatInput', (a+400,b), '')
    set_node_values(_repeat_input, [(0,grow_iter)])
    _repeat_output = add_node(group, 'GeometryNodeRepeatOutput', (a+2000,b), '')
    _repeat_input.pair_with_output(_repeat_output)
    _merge_by_distance = add_node(group,'GeometryNodeMergeByDistance', (a+2200,b), '')
    
    link(_store_attr_0.outputs[0], _repeat_input.inputs[1])

    _join_geometry = add_node(group, 'GeometryNodeJoinGeometry', (a+600,b), '')
    _merge = add_node(group, 'GeometryNodeMergeByDistance', (a+800,b), '')
    _attr_origin = add_attr_node(group, (a+800,b-200), '', [(0,'origin')], 'FLOAT')
    _edge_vertices = add_node(group, 'GeometryNodeInputMeshEdgeVertices', (a+800,b-400), '')
    _sample_idx_1 = add_sample_index_node(group, (a+1000,b-100), '', 'FLOAT', 'POINT')
    _sample_idx_2 = add_sample_index_node(group, (a+1000,b-400), '', 'FLOAT', 'POINT')
    _multiply = add_math_node(group, (a+1200,b-100), '', [(1,0.5)], 'MULTIPLY')
    _add = add_math_node(group, (a+1200,b-400), '', [], 'ADD')
    _store_attr_2 = add_store_attr_node(group, (a+1400,b), '', 'FLOAT', 'EDGE', 'expand')
    _sepa_expand = sepa_by_attr(group, (a+1600,b), 'EDGE', 'expand', 'FLOAT', 2.0, 'NOT_EQUAL')
    _store_attr_3 = add_store_attr_node(group, (a+1800,b), '', 'FLOAT', 'POINT', 'origin')

    link(_repeat_input.outputs[0], _join_geometry.inputs[0])
    link(_store_attr_1.outputs[0], _join_geometry.inputs[0])
    link(_join_geometry.outputs[0], _merge.inputs[0])
    link(_merge.outputs[0], _sample_idx_1.inputs[0])
    link(_merge.outputs[0], _sample_idx_2.inputs[0])
    link(_merge.outputs[0], _store_attr_2.inputs[0])
    link(_attr_origin.outputs[v.na_sockets_out['float']], _sample_idx_1.inputs[v.si_sockets_in['float']])
    link(_attr_origin.outputs[v.na_sockets_out['float']], _sample_idx_2.inputs[v.si_sockets_in['float']])
    link(_edge_vertices.outputs[0], _sample_idx_1.inputs[v.si_sockets_in['index']])
    link(_edge_vertices.outputs[1], _sample_idx_2.inputs[v.si_sockets_in['index']])
    link(_sample_idx_1.outputs[0], _add.inputs[0])
    link(_sample_idx_2.outputs[0], _add.inputs[1])
    link(_add.outputs[0], _multiply.inputs[0])
    link(_multiply.outputs[0], _store_attr_2.inputs[v.sa_sockets_in['float']])
    link(_store_attr_2.outputs[0], _sepa_expand.inputs[0])
    link(_sepa_expand.outputs[0], _store_attr_3.inputs[0])
    link(_store_attr_3.outputs[0], _repeat_output.inputs[0])
    link(_repeat_output.outputs[0], _merge_by_distance.inputs[0])
    link(_merge_by_distance.outputs[0], _output.inputs[0])

    bpy.ops.object.modifier_apply(modifier='GN_expand')


# duplicate the scaffold in cycles l*m*n
def unit_cycles(unit, filename, GN_name, cell_cycles, edge_vecs):
    GN_cycles = create_GeoNode(unit, GN_name, 'GN_Unit_Cycles')
    group = GN_cycles.node_group
    link = group.links.new
    _input, _output = set_io_nodes(group, (-400,100), (400,400))

    # cell cycles in three directions: a, b and c.
    labels = ['a', 'b', 'c']
    for i in range(3):
        _vector = add_node(group,'FunctionNodeInputVector', (-1000, 400*(i+1)), 'vec_'+labels[i])
        _vector.vector = edge_vecs[i]
        _integer = add_node(group,'FunctionNodeInputInt', (-1000, 400*i+270), 'count_'+labels[i])
        _integer.integer = cell_cycles[i]
        _vector_math = add_vec_math_node(group, (-800, 400*(i+1)), '', [], 'MULTIPLY')
        _math = add_math_node(group, (-1000, 400*i+180), '', [(1,1)], 'SUBTRACT')
        _curve_line = add_node(group,'GeometryNodeCurvePrimitiveLine', (-600, 400*(i+1)), '')
        set_node_values(_curve_line, [(1, edge_vecs[i]*(cell_cycles[i]-1))])
        _curve_resample = add_node(group,'GeometryNodeResampleCurve', (-400,400*(i+1)), '')
        set_node_values(_curve_resample, [(2, cell_cycles[i])])
        _instance = add_node(group,'GeometryNodeInstanceOnPoints', (-200, 400*(i+1)), f'Instance_{i+1}')

        link(_vector.outputs[0], _vector_math.inputs[0])
        link(_integer.outputs[0], _math.inputs[0])
        link(_integer.outputs[0], _curve_resample.inputs[2])
        link(_math.outputs[0], _vector_math.inputs[1])
        link(_vector_math.outputs[0], _curve_line.inputs[1])
        link(_curve_line.outputs[0], _curve_resample.inputs[0])
        link(_curve_resample.outputs[0], _instance.inputs[0])
    
    link(_input.outputs[0], get_node(group, 'Instance_1').inputs[2])
    link(get_node(group, 'Instance_1').outputs[0], get_node(group, 'Instance_2').inputs[2])
    link(get_node(group, 'Instance_2').outputs[0], get_node(group, 'Instance_3').inputs[2])

    # merge overlapped points
    _realize_instances = add_node(group,'GeometryNodeRealizeInstances', (0, 400), '')
    _merge_by_distance = add_node(group,'GeometryNodeMergeByDistance', (200, 400), '')
    
    link(get_node(group, 'Instance_3').outputs[0], _realize_instances.inputs[0])
    link(_realize_instances.outputs[0], _merge_by_distance.inputs[0])
    link(_merge_by_distance.outputs[0], _output.inputs[0])

# do not show designed atoms in crystal structure 
def crys_filter(scaffold, filename, GN_name, filters):
    if filters:
        GN_filter = create_GeoNode(scaffold, GN_name, 'GN_Crys_Filter')
        group = GN_filter.node_group
        link =group.links.new
        _input, _output = set_io_nodes(group, (0, -100), (200*len(filters)+200, 0))
        _named_attr = add_attr_node(group, (0, -200), '', [(0,'atom_id')], 'FLOAT')

        for j, filter in enumerate(filters):
            _delete_geometry = add_node(group, 'GeometryNodeDeleteGeometry', (200*j+200, 0), f"Delete_{'%03d' % (j+1)}")
            _compare = add_compare_node(group,(200*j+200, -200), '', [(1, ELEMENTS_DEFAULT[filter][0])], 'FLOAT', 'EQUAL')
            link(_compare.outputs[0], _delete_geometry.inputs[1])
            link(_named_attr.outputs[v.na_sockets_out['float']], _compare.inputs[0])
        
        if filters:
            link(_input.outputs[0], get_node(group,'Delete_001').inputs[0])
            
            if len(filters)>1:
                link(get_node(group,'Delete_001').outputs[0], get_node(group,'Delete_002').inputs[0])
                for j in range(len(filters)-2):
                    link(get_node(group,f"Delete_{'%03d' % (j+1)}").outputs[0], get_node(group,f"Delete_{'%03d' % (j+2)}").inputs[0])
                link(get_node(group,f"Delete_{'%03d' % (len(filters)-1)}").outputs[0], _output.inputs[0])
            else:
                link(get_node(group,'Delete_001').outputs[0], _output.inputs[0])
        else:
            link(_input.outputs[0], _output.inputs[0])

# crystal face section
def crys_face_section(scaffold, filename, bottom, top, hkl):
    GN_crystalFace = create_GeoNode(scaffold, "Crys_Section_"+filename, 'GN_Cross_Section')
    bpy.ops.object.modifier_move_to_index(modifier=GN_crystalFace.name, index=1)
    group = GN_crystalFace.node_group
    _input, _output = set_io_nodes(group, (0,0), (400,0))
    _delete_geo = add_node(group, 'GeometryNodeDeleteGeometry', (200,0), '')
    nodes_link(group, _input, 0, _delete_geo, 0)
    nodes_link(group, _delete_geo, 0, _output, 0)

    _hkl = add_node(group, 'FunctionNodeInputVector', (-1200,0), 'hkl')
    _hkl.vector = hkl
    _normalize = add_vec_math_node(group, (-1000,0), '', [], 'NORMALIZE')
    _scale_bottom = add_vec_math_node(group, (-800,0), '', [(3,bottom)], 'SCALE')
    _scale_top = add_vec_math_node(group, (-800,-150), '', [(3,top)], 'SCALE')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (-800,100), '')
    _subtract_I = add_vec_math_node(group, (-600,0), '', [], 'SUBTRACT')
    _subtract_II = add_vec_math_node(group, (-600,-150), '', [], 'SUBTRACT')
    _dot_product_I = add_vec_math_node(group, (-400,0), '', [], 'DOT_PRODUCT')
    _dot_product_II = add_vec_math_node(group, (-400,-150), '', [], 'DOT_PRODUCT')
    _less_than = add_compare_node(group, (-200,0), '', [], 'FLOAT', 'LESS_THAN')
    _greater_than = add_compare_node(group, (-200,-150), '', [], 'FLOAT', 'GREATER_THAN')
    _boolean_OR = add_node(group, 'FunctionNodeBooleanMath', (0,-100), '')
    set_node_operation(_boolean_OR, 'OR')
    nodes_link(group, _hkl, 0, _normalize, 0)
    nodes_link(group, _normalize, 0, _scale_bottom, 0)
    nodes_link(group, _normalize, 0, _scale_top, 0)
    nodes_link(group, _normalize, 0, _dot_product_I, 1)
    nodes_link(group, _normalize, 0, _dot_product_II, 1)
    nodes_link(group, _input_pos, 0, _subtract_I, 0)
    nodes_link(group, _input_pos, 0, _subtract_II, 0)
    nodes_link(group, _scale_bottom, 0, _subtract_I, 1)
    nodes_link(group, _scale_top, 0, _subtract_II, 1)
    nodes_link(group, _subtract_I, 0, _dot_product_I, 0)
    nodes_link(group, _subtract_II, 0, _dot_product_II, 0)
    nodes_link(group, _dot_product_I, 1, _less_than, v.cp_sockets_in['floatA'])
    nodes_link(group, _dot_product_II, 1, _greater_than, v.cp_sockets_in['floatA'])
    nodes_link(group, _less_than, 0, _boolean_OR, 0)
    nodes_link(group, _greater_than, 0, _boolean_OR, 1)
    nodes_link(group, _boolean_OR, 0, _delete_geo, 1)



# create polyhedra for crystal structure
def crys_polyhedra(scaffold, filename, collection, centers, ligands, r_default, RMax, RMin, bond_length_f):
    GN_polyhedra = create_GeoNode(scaffold, "Crys_Polyhedra_"+filename, filename+' Polyhedra')
    group = GN_polyhedra.node_group
    _input, _output = set_io_nodes(group, (0, 0), (3400, 0))
    _join_hetero_centers = add_node(group,'GeometryNodeJoinGeometry', (2600, 0), 'Hetero Centers')
    
    for i, center in enumerate(centers):
        ctr_atom_id = ELEMENTS_DEFAULT[center][0]
        _join_hetero_ligands = add_node(group,'GeometryNodeJoinGeometry', (2200, -700*i*len(ligands)), 'Hetero Ligands')
        material = bpy.data.materials['ctr_'+center] if ctr_atom_id > 0 else bpy.data.materials['Dummy']
        _set_material = add_mat_node(group, (2400, -700*i*len(ligands)), '', material)
        for j, ligand in enumerate(ligands):
            ligand_atom_id = ELEMENTS_DEFAULT[ligand][0]
            if ctr_atom_id < ligand_atom_id:
                key = f'{center},{ligand}'
            elif ctr_atom_id > ligand_atom_id:
                key = f'{ligand},{center}'
            else:
                key = 'Vac'
            if key not in BONDS_DEFAULT: key = 'Default'
            if r_default:
                RMax = BONDS_DEFAULT[key][3]*bond_length_f
                RMin = 0.0

            b = -700*j-700*i*len(ligands)
            _sepa_center = sepa_by_attr(group, (200,b), 'POINT', 'atom_id', 'INT', ctr_atom_id, 'EQUAL')
            _domain_size = add_node(group, 'GeometryNodeAttributeDomainSize', (400,b-100), '')
            _repeat_input = add_node(group, 'GeometryNodeRepeatInput', (600,b), '')
            _repeat_output = add_node(group, 'GeometryNodeRepeatOutput', (800,b), '')
            _repeat_output.repeat_items.new('INT', 'index')
            _repeat_input.pair_with_output(_repeat_output)
            nodes_link(group, _input, 0, _sepa_center, 0)
            nodes_link(group, _sepa_center, 0, _domain_size, 0)
            nodes_link(group, _domain_size, 0, _repeat_input, 0)

            _sepa_ctr_pt = add_node(group, 'GeometryNodeSeparateGeometry', (800,b+200), 'Ctr_pts')
            _equal = add_compare_node(group, (800,b), '', [], 'INT', 'EQUAL')
            _input_idx = add_node(group, 'GeometryNodeInputIndex', (800,b-200), '')
            _geo_proximity = add_node(group, 'GeometryNodeProximity', (1000,b+200), '')
            _geo_proximity.target_element = 'POINTS'
            _less_equal = add_compare_node(group, (1200,b+200), '', [(1,RMax)], 'FLOAT', 'LESS_EQUAL')
            _not_equal = add_compare_node(group, (1000,b), '', [(1,0.0)], 'FLOAT', 'NOT_EQUAL')
            _boolean_AND_1 = add_node(group, 'FunctionNodeBooleanMath', (1200,b), '')
            _boolean_AND_1.operation = 'AND'
            _boolean_AND_2 = add_node(group, 'FunctionNodeBooleanMath', (1400,b), '')
            _boolean_AND_2.operation = 'AND'
            _attr_atom_id = add_attr_node(group, (1200,b-200), '', [(0,'atom_id')], 'INT')
            _equal_2 = add_compare_node(group, (1400,b-200), '', [(3,ligand_atom_id)], 'INT', 'EQUAL')
            nodes_link(group, _sepa_center, 0, _sepa_ctr_pt, 0)
            nodes_link(group, _input_idx, 0, _equal, 2)
            nodes_link(group, _equal, 0, _sepa_ctr_pt, 1)
            nodes_link(group, _repeat_input, 1, _equal, 3)
            nodes_link(group, _sepa_ctr_pt, 0, _geo_proximity, 0)
            nodes_link(group, _geo_proximity, 1, _less_equal, 0)
            nodes_link(group, _geo_proximity, 1, _not_equal, 0)
            nodes_link(group, _less_equal, 0, _boolean_AND_1, 0)
            nodes_link(group, _not_equal, 0, _boolean_AND_1, 1)
            nodes_link(group, _boolean_AND_1, 0, _boolean_AND_2, 0)
            nodes_link(group, _attr_atom_id, 4, _equal_2, 2)
            nodes_link(group, _equal_2, 0, _boolean_AND_2, 1)

            _sepa_ligand_pt = add_node(group, 'GeometryNodeSeparateGeometry', (1600,b-200), 'Ligand_pts')
            _geo_input = add_node(group, 'NodeGroupInput', (1600,b), '')
            _convex_hull = add_node(group,'GeometryNodeConvexHull', (1800,b-200) , '')
            _join_geo = add_node(group,'GeometryNodeJoinGeometry', (1800, b-350), '')
            _repeat_output.location = (2000,b)
            _add = add_math_node(group, (1000,b-200), '', [(1,1)], 'ADD')
            nodes_link(group, _geo_input, 0, _sepa_ligand_pt, 0)
            nodes_link(group, _boolean_AND_2, 0, _sepa_ligand_pt, 1)
            nodes_link(group, _sepa_ligand_pt, 0, _convex_hull, 0)
            nodes_link(group, _convex_hull, 0, _join_geo, 0)
            nodes_link(group, _repeat_input, 0, _join_geo, 0)
            nodes_link(group, _join_geo, 0, _repeat_output, 0)
            nodes_link(group, _add, 0, _repeat_output, 1)
            nodes_link(group, _repeat_input, 1, _add, 0)
            
            nodes_link(group, _repeat_output, 0, _join_hetero_ligands, 0)
        nodes_link(group, _join_hetero_ligands, 0, _set_material, 0)
        nodes_link(group, _set_material, 0, _join_hetero_centers, 0)
    
    
    _store_attr_atom_id = add_store_attr_node(group, (3000,0), '', 'INT', 'POINT', 'atom_id')
    _sample_idx = add_sample_index_node(group, (3000,200), '', 'INT', 'POINT')
    _attr_atom_id = add_attr_node(group, (2400,400), '', [(0,'atom_id')], 'INT')
    _geo_input = add_node(group, 'NodeGroupInput', (2600,400), '')
    _sample_nearest = add_node(group, 'GeometryNodeSampleNearest', (2600,200), '')
    _domain_size = add_node(group, 'GeometryNodeAttributeDomainSize', (3200,0), '')
    _switch = add_node(group, 'GeometryNodeSwitch', (3200,200), '')
    nodes_link(group, _join_hetero_centers, 0, _store_attr_atom_id, 0)
    nodes_link(group, _store_attr_atom_id, 0, _domain_size, 0)
    nodes_link(group, _domain_size, 0, _switch, 1)
    nodes_link(group, _geo_input, 0, _switch, 14)
    nodes_link(group, _store_attr_atom_id, 0, _switch, 15)
    nodes_link(group, _switch, 6, _output, 0)
    nodes_link(group, _geo_input, 0, _sample_nearest, 0)
    nodes_link(group, _geo_input, 0, _sample_idx, 0)
    nodes_link(group, _sample_nearest, 0, _sample_idx, 7)
    nodes_link(group, _attr_atom_id, 4, _sample_idx, 2)
    nodes_link(group, _sample_idx, 1, _store_attr_atom_id, 7)
    bpy.ops.object.modifier_apply(modifier="Crys_Polyhedra_"+filename, single_user=True)

# poly_atom_id is used to create polyhedra edges
def add_attr_poly_aid(obj, centers, ligands, bond_length_f, RMax, RMin, r_default):
    GN_polyattr = create_GeoNode(obj, "poly_aid_attr", 'Poly aid')
    bpy.ops.object.modifier_move_to_index(modifier="poly_aid_attr", index=0)
    group = GN_polyattr.node_group
    link = group.links.new
    _input, _output = set_io_nodes(group, (0, 0), (1600, 0))

    for i, center in enumerate(centers):
        ctr_atom_id = ELEMENTS_DEFAULT[center][0]
        for j, ligand in enumerate(ligands):
            ligand_atom_id = ELEMENTS_DEFAULT[ligand][0]
            num = i*len(ligands)+j
            b = -500*num

            if ctr_atom_id < ligand_atom_id:
                key = f'{center},{ligand}'
            elif ctr_atom_id > ligand_atom_id:
                key = f'{ligand},{center}'
            else:
                key = 'Vac'
            if key not in BONDS_DEFAULT: key = 'Default'
            if r_default:
                RMax = BONDS_DEFAULT[key][3]*bond_length_f
                RMin = 0.0

            _store_attr_poly_aid_1 = add_store_attr_node(group, (200,b), '', 'FLOAT', 'POINT', 'poly_aid')
            set_node_values(_store_attr_poly_aid_1, [(4,1.0)])
            _store_attr_poly_aid_1.label = f"center id_{'%03d' % num}"
            _store_attr_poly_aid_2 = add_store_attr_node(group, (1200,b), '', 'FLOAT', 'POINT', 'poly_aid')
            set_node_values(_store_attr_poly_aid_2, [(4,2.0)])
            _store_attr_poly_aid_2.label = f"ligand id_{'%03d' % num}"

            _attr_atom_id = add_attr_node(group, (0,b-200), '', [(0,'atom_id')], 'INT')
            _equal_1 = add_compare_node(group, (200,b-200), '', [(3,ctr_atom_id)], 'INT', 'EQUAL')
            _sepa_geo = add_node(group, 'GeometryNodeSeparateGeometry', (400,b-100), '')
            _geo_proximity = add_node(group, 'GeometryNodeProximity', (600,b-100), '')
            _geo_proximity.target_element = 'POINTS'
            _less_equal = add_compare_node(group, (800,b-100), '', [(1,RMax)], 'FLOAT', 'LESS_EQUAL')
            _equal_2 = add_compare_node(group, (800,b-300), '', [(3,ligand_atom_id)], 'INT', 'EQUAL')
            _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (1000,b), '')
            _boolean_AND.operation = 'AND'
            nodes_link(group, _store_attr_poly_aid_1, 0, _store_attr_poly_aid_2, 0)
            nodes_link(group, _store_attr_poly_aid_1, 0, _sepa_geo, 0)
            nodes_link(group, _sepa_geo, 0, _geo_proximity, 0)
            nodes_link(group, _attr_atom_id, 4, _equal_1, 2)
            nodes_link(group, _attr_atom_id, 4, _equal_2, 2)
            nodes_link(group, _equal_1, 0, _store_attr_poly_aid_1, 1)
            nodes_link(group, _equal_1, 0, _sepa_geo, 1)
            nodes_link(group, _geo_proximity, 1, _less_equal, 0)
            nodes_link(group, _less_equal, 0, _boolean_AND, 0)
            nodes_link(group, _equal_2, 0, _boolean_AND, 1)
            nodes_link(group, _boolean_AND, 0, _store_attr_poly_aid_2, 1)

    if len(centers)*len(ligands) == 0:
        nodes_link(group, _input, 0, _output, 0)
    else:
        nodes_link(group, _input, 0, get_node(group, 'center id_000'), 0)
        for k in range(len(centers)*len(ligands)-1):
            nodes_link(group, get_node(group, f"ligand id_{'%03d' % k}"), 0, get_node(group, f"center id_{'%03d' % (k+1)}"), 0)
        nodes_link(group, get_node(group, f"ligand id_{'%03d' % (len(centers)*len(ligands)-1)}"), 0, _output, 0)
    bpy.ops.object.modifier_apply(modifier="poly_aid_attr")

# create polyhedra edges
def crys_polyhedra_edges(polyhedra, filename, bond_subdiv, radii):
    GN_polyhedra = create_GeoNode(polyhedra, "Crys_Polyhedra_Edges_"+filename, filename+' Polyhedra Edges')
    group = GN_polyhedra.node_group
    link = group.links.new
    _input, _output = set_io_nodes(group, (0,100), (1000,100))
    _delete = add_node(group, 'GeometryNodeDeleteGeometry', (200,0), '')
    _delete.mode = 'ONLY_FACE'
    _bonds = add_node(group, 'Mol3D_MT_Bonds', (400,0), 'Bonds')
    set_node_values(_bonds, [(1,'bond_radii_B'),(2,radii),(3,bond_subdiv)])
    mat = bpy.data.materials['Hetero']
    _set_material = add_mat_node(group, (600,0), '', mat)
    _join_geo = add_node(group,'GeometryNodeJoinGeometry', (800,100), '')
    nodes_link(group, _input, 0, _delete, 0)
    nodes_link(group, _delete, 0, _bonds, 0)
    nodes_link(group, _bonds, 0, _set_material, 0)
    nodes_link(group, _set_material, 0, _join_geo, 0)
    nodes_link(group, _input, 0, _join_geo, 0)
    nodes_link(group, _join_geo, 0, _output, 0)



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

class CrystalDeformation():
    def __init__(self, ao, deform_mode):
        self.mode = deform_mode
        self.obj = ao

    # determine if the deformation modifier exists 
    def is_deform_exist(self, object):
        FLAG = False
        for modifier in object.modifiers:
            if modifier.name.split('_')[-1] == self.mode:
                FLAG = True
                break
        return FLAG
    
    def updateDeform(self):
        prefix = self.obj.name.split('_')[0]
        sign = self.obj.name.find('_') + 1
        crys_name = self.obj.name[sign:].split('.')[0]
        Crys_Name = crys_name + '_Crys_' + self.mode
        Poly_Name = crys_name + '_Polyhedra_' + self.mode
        ng_name = crys_name + ' ' + self.mode
        exist_name = crys_name + '_' + self.mode

        FLAG = self.is_deform_exist(self.obj)
        if not FLAG:
            # create judgement
            if prefix == 'Crystal':
                GN_name = Crys_Name
                GN_deform = _node.create_GeoNode(self.obj, GN_name, ng_name)
                bpy.ops.object.modifier_move_up(modifier=GN_name)
            elif prefix == 'Polyhedra':
                GN_name = Poly_Name
                GN_deform = self.obj.modifiers.new(GN_name, 'NODES')
                bpy.ops.node.new_geometry_node_group_assign()
                crys_obj = _mesh.get_object('Crystal_'+crys_name)
                crys_FLAG = self.is_deform_exist(crys_obj)
                bpy.ops.object.modifier_move_up(modifier=GN_name)

                if crys_FLAG:  # crys_obj already has this deformation modifier
                    GN_deform.node_group = bpy.data.node_groups[ng_name]
                else:  # crys_obj doesn't have this deformation modifier
                    GN_deform.node_group.name = ng_name
                    GN_name = Crys_Name
                    bpy.context.view_layer.objects.active = crys_obj
                    crys_GN_deform = crys_obj.modifiers.new(GN_name, 'NODES')
                    bpy.ops.node.new_geometry_node_group_assign()
                    crys_GN_deform.node_group = bpy.data.node_groups[ng_name]
                    bpy.ops.object.modifier_move_up(modifier=GN_name)
            else:
                GN_name = self.obj.name + '_' + self.mode
                GN_deform = _node.create_GeoNode(self.obj, GN_name, self.obj.name + ' ' + self.mode)
        else:
            if prefix == 'Crystal':
                GN_deform = self.obj.modifiers[Crys_Name]
            elif prefix == 'Polyhedra':
                GN_deform = self.obj.modifiers[Poly_Name]
            else:
                GN_deform = self.obj.modifiers[exist_name]

        return GN_deform


    def roll_deformation(self, GN_deform, roll_angle, clockwise, axis, transform, toggle):
        type = self.obj['Type']
        sign = self.obj.name.find('_') + 1
        crys_name = self.obj.name[sign:].split('.')[0]
        # get cell_cycles and edge_vecs from Crystal Unit Cycles Modifier
        if type == 'crystal' or type == 'polyhedra':
            Crystal = _mesh.get_object('Crystal_'+crys_name)
            GN_cycles = Crystal.modifiers['Crys_Cycles_'+crys_name]
            group = GN_cycles.node_group

            node_a = get_node(group, 'count_a')
            node_b = get_node(group, 'count_b')
            node_c = get_node(group, 'count_c')
            cell_cycles = (node_a.integer, node_b.integer, node_c.integer)

            Cell_Edges = _mesh.get_object('Cell_Edges_'+crys_name)
            length_a, length_b, length_c = eval(Cell_Edges['cell lengths'])
            angle_alpha, angle_beta, angle_gamma = eval(Cell_Edges['cell angles'])

            unit_x = abs(_math.fract_to_cartn((1,0,0),length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)[0])
            unit_y = abs(_math.fract_to_cartn((0,1,0),length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)[1])
            unit_z = abs(_math.fract_to_cartn((0,0,1),length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma)[2])

            Lx = cell_cycles[0] * unit_x
            Ly = cell_cycles[1] * unit_y
            Lz = cell_cycles[2] * unit_z  
        else:
            Lx = self.obj.bound_box[4][0] - self.obj.bound_box[0][0]
            Ly = self.obj.bound_box[2][1] - self.obj.bound_box[0][1]
            Lz = self.obj.bound_box[1][2] - self.obj.bound_box[0][2]
            
        R_yz = Lx/(roll_angle*3.1416/180)
        R_xy = Lz/(roll_angle*3.1416/180)
        R_zx = Ly/(roll_angle*3.1416/180)  
            
        if axis == 'x':
            R = R_zx if toggle else R_xy
            H = Lz if toggle else Ly
            u,v,w = (1,0,2) if toggle else (2,0,1)
        elif axis == 'y':
            R = R_xy if toggle else R_yz
            H = Lx if toggle else Lz
            u,v,w = (2,1,0) if toggle else (0,1,2)
        else:
            R = R_yz if toggle else R_zx
            H = Ly if toggle else Lx
            u,v,w = (0,2,1) if toggle else (1,2,0)

        group = GN_deform.node_group
        link = group.links.new
        _input, _output = set_io_nodes(group, (800, 0), (1400, 0))

        # delete all old nodes except _input and _output
        delete_ios(GN_deform.node_group)

        _transform = add_node(group,'GeometryNodeTransform', (1000, 100), '')
        set_node_values(_transform, [(3, transform)])
  
        _set_pos = add_node(group,'GeometryNodeSetPosition', (1200, 0), '')
        _input_pos = add_node(group,'GeometryNodeInputPosition', (0, -200), '')

        _separate_xyz = add_node(group,'ShaderNodeSeparateXYZ', (200, -200), '')
        _combine_xyz = add_node(group,'ShaderNodeCombineXYZ', (1000, -200), '')

        _multiply1 = add_math_node(group, (400, -200), '', [(1,1/R)], 'MULTIPLY')
        _add = add_math_node(group, (400, -600), '', [(1,R)], 'ADD')
        _sin = add_math_node(group, (600, -200), '', [], 'SINE')
        _multiply2 = add_math_node(group, (800, -200), '', [], 'MULTIPLY')
        _cosine = add_math_node(group, (600, -400), '', [], 'COSINE')
        _multiply3 = add_math_node(group, (800, -400), '', [], 'MULTIPLY')

        if clockwise:
            multi_value = 1
            add_value = 0
            subs_value = R
        else:
            multi_value = -1
            add_value = H
            subs_value = R + H

        _subtract = add_math_node(group, (1000, -400), '', [(1,subs_value)], 'SUBTRACT')
        _multiply4 = add_math_node(group, (400, -400), '', [(1,multi_value)], 'MULTIPLY')
        _multiply5 = add_math_node(group, (1200, -400), '', [(1,multi_value)], 'MULTIPLY')
        _add1 = add_math_node(group, (600, -600), '', [(1,add_value)], 'ADD')
        
        link(_input_pos.outputs[0], _separate_xyz.inputs[0])
        link(_separate_xyz.outputs[u], _multiply1.inputs[0])
        link(_multiply1.outputs[0], _sin.inputs[0])
        link(_sin.outputs[0], _multiply2.inputs[0])
        
        link(_separate_xyz.outputs[w], _multiply4.inputs[0])
        link(_multiply4.outputs[0], _add.inputs[0])
        link(_add.outputs[0], _add1.inputs[0])
        link(_add1.outputs[0], _multiply2.inputs[1])
        link(_multiply2.outputs[0], _combine_xyz.inputs[u])
        link(_multiply1.outputs[0], _cosine.inputs[0])
        link(_cosine.outputs[0], _multiply3.inputs[0])

        link(_add1.outputs[0], _multiply3.inputs[1])
        link(_multiply3.outputs[0], _subtract.inputs[0])
        link(_subtract.outputs[0], _multiply5.inputs[0])
        link(_multiply5.outputs[0], _combine_xyz.inputs[w])

        link(_input.outputs[0], _transform.inputs[0])
        link(_transform.outputs[0], _set_pos.inputs[0])
        link(_set_pos.outputs[0], _output.inputs[0])
        link(_separate_xyz.outputs[v], _combine_xyz.inputs[v])
        
        # no position change if rolling angle is zero
        if roll_angle != 0:
            link(_combine_xyz.outputs[0], _set_pos.inputs[2])


    def wave_deformation(self, GN_deform, wavelength, amplitude, phase_angle, axis, transform, toggle):
        group = GN_deform.node_group
        link = group.links.new
        _input, _output = set_io_nodes(group, (800, 0), (1400, 0))

        # delete all old nodes except _input and _output
        for node in group.nodes[2:]: group.nodes.remove(node)
        
        _transform = add_node(group,'GeometryNodeTransform', (1000, 150), '')
        set_node_values(_transform, [(3, transform)])
        _set_pos = add_node(group,'GeometryNodeSetPosition', (1200, 0), '')
        _input_pos = add_node(group,'GeometryNodeInputPosition', (-400, -200), '')
        _separate_xyz = add_node(group,'ShaderNodeSeparateXYZ', (-200, -200), '')
        _combine_xyz = add_node(group,'ShaderNodeCombineXYZ', (1000, -200), '')
        value1 = 2*3.1416/wavelength if wavelength != 0 else 2*3.1416/0.0000001
        _multiply1 = add_math_node(group, (0, -400), '', [(1,value1)], 'MULTIPLY')
        value = phase_angle*wavelength*2*3.1416/180
        _add = add_math_node(group, (200, -400), '', [(1,value)], 'ADD')
        _sin = add_math_node(group, (400, -400), '', [], 'SINE')
        _multiply2 = add_math_node(group, (600, -400), '', [(1,amplitude)], 'MULTIPLY')
        _add1 = add_math_node(group, (800, -200), '', [], 'ADD')

        if axis == 'x':
            u,v,w = (0,2,1) if toggle else (0,1,2)
        elif axis == 'y':
            u,v,w = (1,0,2) if toggle else (1,2,0)
        else:
            u,v,w = (2,1,0) if toggle else (2,0,1)

        link(_input_pos.outputs[0], _separate_xyz.inputs[0])
        link(_separate_xyz.outputs[u], _multiply1.inputs[0])
        link(_multiply1.outputs[0], _add.inputs[0])
        link(_add.outputs[0], _sin.inputs[0])
        link(_sin.outputs[0], _multiply2.inputs[0])
        link(_multiply2.outputs[0], _add1.inputs[0])
        link(_separate_xyz.outputs[w], _add1.inputs[1])
        link(_separate_xyz.outputs[u], _combine_xyz.inputs[u])
        link(_separate_xyz.outputs[v], _combine_xyz.inputs[v])
        link(_add1.outputs[0], _combine_xyz.inputs[w])
        link(_input.outputs[0], _transform.inputs[0])
        link(_transform.outputs[0], _set_pos.inputs[0])
        link(_set_pos.outputs[0], _output.inputs[0])

        # no position change if wavelength is zero
        if wavelength != 0:
            link(_combine_xyz.outputs[0], _set_pos.inputs[2])

    def noise_deformation(self, GN_deform, nscale, amplitude, axis, transform):
        group = GN_deform.node_group
        link = group.links.new
        a,b=(0,0)
        _input, _output = set_io_nodes(group, (a, b), (a+600, b))

        # delete all old nodes except _input and _output
        for node in group.nodes[2:]: group.nodes.remove(node)
        
        _transform = add_node(group,'GeometryNodeTransform', (a+200, b+150), '')
        set_node_values(_transform, [(3, transform)])
        _set_pos = add_node(group,'GeometryNodeSetPosition', (a+400, b), '')
    
        _combine_xyz = add_node(group,'ShaderNodeCombineXYZ', (a+200, b-200), '')
        _multiply = add_math_node(group, (a, b-400), '', [(1,amplitude)], 'MULTIPLY')
        _subtract = add_math_node(group, (a-200, b-400), '', [], 'SUBTRACT')
        _noise_texture = add_node(group,'ShaderNodeTexNoise', (a-400, b-400), '')
        set_node_values(_noise_texture, [(2,nscale),(3,0)])
        
        if axis == 'x':
            socket = 0
        elif axis == 'y':
            socket = 1
        else:
            socket = 2
        
        link(_noise_texture.outputs[0], _subtract.inputs[0])
        link(_subtract.outputs[0], _multiply.inputs[0])
        link(_multiply.outputs[0], _combine_xyz.inputs[socket])
        link(_input.outputs[0], _transform.inputs[0])
        link(_transform.outputs[0], _set_pos.inputs[0])
        link(_combine_xyz.outputs[0], _set_pos.inputs[3])
        link(_set_pos.outputs[0], _output.inputs[0])



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
)
#####################################################################################

def peptide_meshline(group, location):
    a,b = location
    _sepa_alpha_C = add_node(group, 'GeometryNodeSeparateGeometry', (a,b), 'alpha_C')
    _named_attr = add_attr_node(group, (a, b-200), '', [(0,'is_alpha')], 'BOOLEAN')
    _sample_index = add_node(group, 'GeometryNodeSampleIndex', (a+400, b+150), '')
    set_node_datatype(_sample_index, 'FLOAT_VECTOR')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a+200, b), '')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a+200, b-100), '')
    _domain_size = add_node(group, 'GeometryNodeAttributeDomainSize', (a+200, b-200), '')
    _mesh_line = add_node(group, 'GeometryNodeMeshLine', (a+400, b-100), '')
    _set_pos = add_node(group, 'GeometryNodeSetPosition', (a+600, b), '')
    _delete_mismatch, _greater_than = delete_mismatch(group, (a+600, b-250), max_length=7.5)
    nodes_link(group, _named_attr, 3, _sepa_alpha_C, 1)
    nodes_link(group, _sepa_alpha_C, 0, _sample_index, 0)
    nodes_link(group, _sepa_alpha_C, 0, _domain_size, 0)
    nodes_link(group, _input_pos, 0, _sample_index, 3)
    nodes_link(group, _input_idx, 0, _sample_index, 7)
    nodes_link(group, _domain_size, 0, _mesh_line, 0)
    nodes_link(group, _mesh_line, 0, _set_pos, 0)
    nodes_link(group, _sample_index, 2, _set_pos, 2)
    nodes_link(group, _set_pos, 0, _delete_mismatch, 0)
    return (_sepa_alpha_C, _delete_mismatch)
    
def nucleic_meshline(group, location):
    a,b = location
    _sepa_nucleic = separate_geometry(group, sec_struct_name='nucleic', domain='POINT', 
                                         attr_name='sec_struct', sec_struct_id=0, location=(a,b))
    _sepa_C1 = add_node(group, 'GeometryNodeSeparateGeometry', (a+400,b-200), 'Nucleic_C1')
    _attr_atom_name = add_attr_node(group, (a+200, b-400), '', [(0,'atom_name')], 'INT')
    _equal = add_compare_node(group, (a+400,b-400), '', [(3,61)], 'INT', 'EQUAL')
    _sample_index = add_node(group, 'GeometryNodeSampleIndex', (a+800, b+150), '')
    set_node_datatype(_sample_index, 'FLOAT_VECTOR')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a+600, b-50), '')
    _domain_size = add_node(group, 'GeometryNodeAttributeDomainSize', (a+600, b-200), '')
    _mesh_line = add_node(group, 'GeometryNodeMeshLine', (a+800, b-100), '')
    _set_pos = add_node(group, 'GeometryNodeSetPosition', (a+1000, b), '')
    _delete_mismatch, _greater_than = delete_mismatch(group, (a+1000,b-350), max_length=7.7)
    nodes_link(group, _sepa_nucleic, 0, _sepa_C1, 0)
    nodes_link(group, _attr_atom_name, 4, _equal, 2)
    nodes_link(group, _equal, 0, _sepa_C1, 1)
    nodes_link(group, _sepa_C1, 0, _sample_index, 0)
    nodes_link(group, _sepa_C1, 0, _domain_size, 0)
    nodes_link(group, _input_idx, 0, _sample_index, 7)
    nodes_link(group, _domain_size, 0, _mesh_line, 0)
    nodes_link(group, _mesh_line, 0, _set_pos, 0)
    nodes_link(group, _sample_index, 2, _set_pos, 2)
    nodes_link(group, _set_pos, 0, _delete_mismatch, 0)

    _attr_atom_name = add_attr_node(group, (a+200, b+300), '', [(0,'atom_name')], 'INT')
    _equal_1 = add_compare_node(group, (a+400,b+400), '', [(3,55)], 'INT', 'EQUAL')
    _equal_2 = add_compare_node(group, (a+400,b+200), '', [(3,57)], 'INT', 'EQUAL')
    _sepa_geo_1 = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b+400), '')
    _sepa_geo_2 = add_node(group, 'GeometryNodeSeparateGeometry', (a+600,b+200), '')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a+400, b+650), '')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a+400, b+550), '')
    _sample_idx_1 = add_node(group, 'GeometryNodeSampleIndex', (a+600, b+900), '')
    set_node_datatype(_sample_idx_1, 'FLOAT_VECTOR')
    _sample_idx_2 = add_node(group, 'GeometryNodeSampleIndex', (a+600, b+650), '')
    set_node_datatype(_sample_idx_2, 'FLOAT_VECTOR')
    _vec_mix_1 = add_node(group, 'ShaderNodeMix', (a+800, b+900), '')
    set_node_datatype(_vec_mix_1, 'VECTOR')
    _vec_mix_2 = add_node(group, 'ShaderNodeMix', (a+800, b+650), '')
    set_node_datatype(_vec_mix_2, 'VECTOR')
    nodes_link(group, _attr_atom_name, 4, _equal_1, 2)
    nodes_link(group, _attr_atom_name, 4, _equal_2, 2)
    nodes_link(group, _equal_1, 0, _sepa_geo_1, 1)
    nodes_link(group, _equal_2, 0, _sepa_geo_2, 1)
    nodes_link(group, _sepa_nucleic, 0, _sepa_geo_1, 0)
    nodes_link(group, _sepa_nucleic, 0, _sepa_geo_2, 0)
    nodes_link(group, _input_pos, 0, _sample_idx_1, 3)
    nodes_link(group, _input_pos, 0, _sample_idx_2, 3)
    nodes_link(group, _input_pos, 0, _vec_mix_2, 5)
    nodes_link(group, _input_idx, 0, _sample_idx_1, 7)
    nodes_link(group, _input_idx, 0, _sample_idx_2, 7)
    nodes_link(group, _sepa_geo_1, 0, _sample_idx_1, 0)
    nodes_link(group, _sepa_geo_2, 0, _sample_idx_2, 0)
    nodes_link(group, _sample_idx_1, 2, _vec_mix_1, 4)
    nodes_link(group, _sample_idx_2, 2, _vec_mix_1, 5)
    nodes_link(group, _vec_mix_1, 1, _vec_mix_2, 4)
    nodes_link(group, _vec_mix_2, 1, _sample_index, 3)
    return (_sepa_nucleic, _sepa_C1, _delete_mismatch)


# rotate cross-section shape using 'guide_X' and 'guide_Z' vectors
def rotate_euler(group, angle, location):
    link = group.links.new
    a, b = location
    _align_euler2vec_1 = add_node(group, 'FunctionNodeAlignEulerToVector',(a,b), '')
    _align_euler2vec_1.axis = 'Z'
    _attr_guide_Z = add_attr_node(group, (a, b-200), '', [(0,'guide_Z')], 'FLOAT_VECTOR')
    _align_euler2vec_2 = add_node(group, 'FunctionNodeAlignEulerToVector',(a+200,b), '')
    _align_euler2vec_2.axis = 'X'
    _align_euler2vec_2.pivot_axis = 'Z'
    _attr_guide_X = add_attr_node(group, (a+200, b-200), '', [(0,'guide_X')], 'FLOAT_VECTOR')
    _rotate_euler = add_node(group, 'FunctionNodeRotateEuler',(a+400,b), '')
    if v.version == (4,0,2):
        _rotate_euler.type = 'AXIS_ANGLE'
    else:
        _rotate_euler.rotation_type = 'AXIS_ANGLE'
    set_node_values(_rotate_euler, [(3, angle)])
    link(_attr_guide_Z.outputs[0], _align_euler2vec_1.inputs[2])
    link(_attr_guide_X.outputs[0], _align_euler2vec_2.inputs[2])
    link(_attr_guide_X.outputs[0], _rotate_euler.inputs[2])
    link(_align_euler2vec_1.outputs[0], _align_euler2vec_2.inputs[0])
    link(_align_euler2vec_2.outputs[0], _rotate_euler.inputs[0])
    return _rotate_euler

def cross_section_transform(group, sharp, width, thickness, location, resolution):
    link = group.links.new
    a,b = location
    _curve_circle = add_node(group, 'GeometryNodeCurvePrimitiveCircle', (a, b), '')
    circle_res = 4 if sharp else 8*resolution
    set_node_values(_curve_circle, [(0,circle_res)])
    _transform1 = add_node(group, 'GeometryNodeTransform', (a+200, b), '')
    set_node_values(_transform1, [(2,(0.0, 0.0, 3.1416/4))])
    _transform2 = add_node(group, 'GeometryNodeTransform', (a+400, b), 'Cross_Section')    # Width and Thickness
    set_node_values(_transform2, [(2,(0, 3.1416/2, 0.0)), (3,(width, thickness, 1.0))])
    link(_curve_circle.outputs[0], _transform1.inputs[0])
    link(_transform1.outputs[0], _transform2.inputs[0])
    return (_curve_circle, _transform2)

def cross_instance(group, location):
    link = group.links.new
    a,b = location
    _instance = add_node(group, 'GeometryNodeInstanceOnPoints', (a, b), '')
    _realize_instance = add_node(group, 'GeometryNodeRealizeInstances', (a, b-300), '')
    _sample_index = add_node(group, 'GeometryNodeSampleIndex', (a, b-400), '')
    set_node_datatype(_sample_index, 'FLOAT_VECTOR')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a-200, b-400), '')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-200, b-500), '')
    link(_instance.outputs[0], _realize_instance.inputs[0])
    link(_realize_instance.outputs[0], _sample_index.inputs[0])
    nodes_link(group, _input_pos, 0, _sample_index, v.si_sockets_in['vec'])
    nodes_link(group, _input_idx, 0, _sample_index, v.si_sockets_in['index'])
    return (_instance, _sample_index)


def beta_arrow_bottom_shift(group, factor, location, resolution, multiply):
    link = group.links.new
    a,b = location
    _set_pos = add_node(group, 'GeometryNodeSetPosition', (a, b), '')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a-800, b-150), '')
    _input_idx = add_node(group, 'GeometryNodeInputIndex', (a-800, b), '')
    _add = add_math_node(group, (a-600,b), '', [(1,1.0)], 'ADD')
    _field_idx = add_node(group, 'GeometryNodeFieldAtIndex', (a-400,b), '')
    set_node_datatype(_field_idx, 'FLOAT_VECTOR')
    _mix_vec = add_node(group, 'ShaderNodeMix', (a-200,b), '')
    set_node_values(_mix_vec, [(0,factor)])
    set_node_datatype(_mix_vec, 'VECTOR')
    _attr_scale_beta = add_attr_node(group, (a-600, b-200), '', [(0,'scale_beta')], 'FLOAT')
    _field_idx2 = add_node(group, 'GeometryNodeFieldAtIndex', (a-400,b-200), '')
    set_node_datatype(_field_idx2, 'FLOAT')
    critical = 1-1/resolution if resolution != 0 else 0
    _equal = add_compare_node(group, (a-200,b-200), '', [(1, critical)], 'FLOAT', 'EQUAL')

    link(_input_idx.outputs[0], _add.inputs[0])
    link(_add.outputs[0], _field_idx.inputs[0])
    link(_add.outputs[0], _field_idx2.inputs[0])
    nodes_link(group, _field_idx, v.fi_sockets_out['vec'], _mix_vec, 4)
    nodes_link(group, _input_pos, 0, _field_idx, v.fi_sockets_in['vec'])
    link(_input_pos.outputs[0], _mix_vec.inputs[5])
    if multiply != 1: link(_mix_vec.outputs[1], _set_pos.inputs[2])
    link(_attr_scale_beta.outputs[1], _field_idx2.inputs[1])
    link(_field_idx2.outputs[0], _equal.inputs[0])
    link(_equal.outputs[0], _set_pos.inputs[1])

    return (_set_pos, _mix_vec, _equal)

#####################################################################################
def Ribbon_nodes(group, radius, resolution, vertices, subdivision, color_type, category, display_mode, hide_segs):
    a, b = (0, 0)
    _input, _output = set_io_nodes(group, (a-600, b), (a+2000,b))
    delete_ios(group)
    _geometry, a = display(group, (a-400,b), _input, 'R', category, display_mode)
    
    _nucleic_meshline = add_node(group, 'Mol3D_MT_Nucleic_Restruct', (a+200,b), '')
    _peptide_meshline = add_node(group, 'Mol3D_MT_Peptide_Restruct', (a+200,b-200), '')
    _join_1 = add_node(group, 'GeometryNodeJoinGeometry', (a+400, b), '')
    _join_2 = add_node(group, 'GeometryNodeJoinGeometry', (a+400, b-200), '')
    _bioattr_reset = add_node(group, 'Mol3D_MT_BioAttr_Reset', (a+600,b), '')
    nodes_link(group, _geometry, 0, _nucleic_meshline, 0)
    nodes_link(group, _geometry, 0, _peptide_meshline, 0)
    nodes_link(group, _nucleic_meshline, 0, _join_1, 0)
    nodes_link(group, _peptide_meshline, 0, _join_1, 0)
    nodes_link(group, _nucleic_meshline, 1, _join_2, 0)
    nodes_link(group, _peptide_meshline, 1, _join_2, 0)
    nodes_link(group, _join_1, 0, _bioattr_reset, 0)
    nodes_link(group, _join_2, 0, _bioattr_reset, 1)
    
    _geometry, socket_out, a = hide_segs_by_input_texts(group, (a+800,b), hide_segs, _bioattr_reset)
    _mesh2curve, _curve_res, _resample_curve = mesh_to_curve(group, (a+200,b), resolution)
    _curve2mesh = add_node(group, 'GeometryNodeCurveToMesh', (a+1200, b), '')
    set_node_values(_curve2mesh, [(2,True)])    # Fill Caps True
    _curve_circle = add_node(group, 'GeometryNodeCurvePrimitiveCircle', (a+1000, b-200), 'circle')
    set_node_values(_curve_circle, [(0,vertices),(4,radius)])   # circle Resolution
    _subdiv_surface = add_node(group, 'GeometryNodeSubdivisionSurface', (a+1400, b), '')
    set_node_values(_subdiv_surface, [(1,subdivision)])   # subdivision level
    _shade_smooth = add_node(group, 'GeometryNodeSetShadeSmooth', (a+1600, b), '')
    if color_type == 'C':
        mat = bpy.data.materials['ChainID']
    elif color_type == 'H':
        mat = bpy.data.materials['ChainID']
    elif color_type == 'R':
        mat = bpy.data.materials['Rainbow']
    _set_material = add_mat_node(group, (a+1800, b), '', mat)
    _output.location = (a+2000,b)

    nodes_link(group, _geometry, socket_out, _mesh2curve, 0)
    nodes_link(group, _resample_curve, 0, _curve2mesh, 0)
    nodes_link(group, _curve_circle, 0, _curve2mesh, 1)    
    nodes_link(group, _curve2mesh, 0, _subdiv_surface, 0)
    nodes_link(group, _subdiv_surface, 0, _shade_smooth, 0)
    nodes_link(group, _shade_smooth, 0, _set_material, 0)
    nodes_link(group, _set_material, 0, _output, 0)
    


#####################################################################################
def Cartoon_nodes(group, sharp, width, thickness, coil_radius, resolution, subdivision,
                  color_type, category, display_mode, hide_segs):
    link = group.links.new
    a, b = (0, 0)
    _input, _output = set_io_nodes(group, (a-1600, b), (a+2400,b))
    delete_ios(group)
    # set attributes: guide_X, guide_Z and scale_beta
    _geometry, a = display(group, (a-1400,b), _input, 'C', category, display_mode)
    
    _bioattr_preset = add_node(group, 'Mol3D_MT_BioAttr_Preset', (a,b), '')
    _nucleic_meshline = add_node(group, 'Mol3D_MT_Nucleic_Restruct', (a+200,b), '')
    _peptide_meshline = add_node(group, 'Mol3D_MT_Peptide_Restruct', (a+200,b-200), '')
    _join_1 = add_node(group, 'GeometryNodeJoinGeometry', (a+400, b), '')
    _join_2 = add_node(group, 'GeometryNodeJoinGeometry', (a+400, b-200), '')
    _bioattr_reset = add_node(group, 'Mol3D_MT_BioAttr_Reset', (a+600,b), '')
    nodes_link(group, _geometry, 0, _bioattr_preset, 0)
    nodes_link(group, _bioattr_preset, 0, _nucleic_meshline, 0)
    nodes_link(group, _bioattr_preset, 1, _peptide_meshline, 0)
    nodes_link(group, _nucleic_meshline, 0, _join_1, 0)
    nodes_link(group, _peptide_meshline, 0, _join_1, 0)
    nodes_link(group, _nucleic_meshline, 1, _join_2, 0)
    nodes_link(group, _peptide_meshline, 1, _join_2, 0)
    nodes_link(group, _join_1, 0, _bioattr_reset, 0)
    nodes_link(group, _join_2, 0, _bioattr_reset, 1)
    
    # nodes for structure
    _geometry, socket_out, a = hide_segs_by_input_texts(group, (a+800,b), hide_segs, _bioattr_reset)
    # alpha_helix
    _alpha_helix = add_node(group, 'Mol3D_MT_Alpha_Helix', (a+200,b), '')
    set_node_values(_alpha_helix, [(1,resolution),(2,width),(3,thickness),(4,sharp)])
    # beta_sheet
    _beta_sheet = add_node(group, 'Mol3D_MT_Beta_Sheet', (a+200,b-250), '')
    set_node_values(_beta_sheet, [(1,resolution),(2,width),(3,thickness),(5,sharp)])
    # loop_coil
    _loop_coil = add_node(group, 'Mol3D_MT_Loop_Coil', (a+200,b+200), '')
    set_node_values(_loop_coil, [(1,resolution),(2,coil_radius)])
    # nucleic
    _nucleic = add_node(group, 'Mol3D_MT_Nucleic', (a+200,b+400), '')
    set_node_values(_nucleic, [(1,resolution),(2,coil_radius)])

    nodes_link(group, _geometry, socket_out, _alpha_helix, 0)
    nodes_link(group, _geometry, socket_out, _beta_sheet, 0)
    nodes_link(group, _geometry, socket_out, _loop_coil, 0)
    nodes_link(group, _geometry, socket_out, _nucleic, 0)

    # ------------------------------------------------------------------------------------------
    _join_geometry = add_node(group, 'GeometryNodeJoinGeometry', (a+400, b), '')
    _edge_angle = add_node(group, 'GeometryNodeInputMeshEdgeAngle', (a+600, b-400), '')
    _store_attr = add_node(group, 'GeometryNodeStoreNamedAttribute', (a+600, b), '')
    _store_attr.data_type = 'BOOLEAN'
    _store_attr.domain = 'EDGE'
    set_node_values(_store_attr, [(2,'sharp_edge')])
    critical = 1.05 if subdivision == 0 else 2*1.05
    _compare2 = add_compare_node(group, (a+600,b-200), 'Greater', [(1,critical)], 'FLOAT', 'GREATER_THAN')   # greater than critical
    _subdiv_surf = add_node(group, 'GeometryNodeSubdivisionSurface', (a+800, b), 'Subdivision')   # subdivision level
    set_node_values(_subdiv_surf, [(1,subdivision)])

    if color_type == 'C':
        mat = bpy.data.materials['ChainID']
    elif color_type == 'H':
        mat = bpy.data.materials['Cartoon']
    elif color_type == 'R':
        mat = bpy.data.materials['Rainbow']
    _mat_cartoon = add_mat_node(group, (a+1000,b), '', mat)
    _output.location = (a+1200,b)

    link(_nucleic.outputs[0], _join_geometry.inputs[0])
    link(_alpha_helix.outputs[0], _join_geometry.inputs[0])
    link(_beta_sheet.outputs[0], _join_geometry.inputs[0])
    link(_loop_coil.outputs[0], _join_geometry.inputs[0])
    link(_join_geometry.outputs[0], _store_attr.inputs[0])
    link(_edge_angle.outputs[0], _compare2.inputs[0])
    nodes_link(group, _compare2, 0, _store_attr, v.sa_sockets_in['bool'])
    link(_store_attr.outputs[0], _subdiv_surf.inputs[0])
    link(_subdiv_surf.outputs[0], _mat_cartoon.inputs[0])
    link(_mat_cartoon.outputs[0], _output.inputs[0])


def Surface_nodes(group, probe_radius, surface_mode, surface_res, smooth_iter, color_type, blur_iteration,
                  category, display_mode, hide_segs):
    link = group.links.new
    a, b = (0, 0)
    _input, _output = set_io_nodes(group, (a,b), (a+200,b))
    delete_ios(group)
    # set attributes: guide_X, guide_Z and scale_beta
    _backbone, a = display(group, (a+200,b), _input, 'C', category, display_mode)
    # nodes for structure
    _geometry, socket_out, a = hide_segs_by_input_texts(group, (a+200,b), hide_segs, _backbone)

    _instance = add_node(group,'GeometryNodeInstanceOnPoints', (a+180,b), '')
    _ico_sphere = add_node(group, 'GeometryNodeMeshIcoSphere', (a,b-50), '')
    set_node_values(_ico_sphere, [(1,3)])
    _attr_radii_vdW = add_attr_node(group, (a-200,b-200), '', [(0,'radii_F')], 'FLOAT')
    _add = add_math_node(group, (a,b-200), '', [(1,probe_radius)], 'ADD')
    _realize_instance = add_node(group,'GeometryNodeRealizeInstances', (a+360,b), '')
    _mesh2vol = add_node(group, 'GeometryNodeMeshToVolume', (a+560,b), '')
    _mesh2vol.resolution_mode = 'VOXEL_SIZE'
    voxel_size = max(1-surface_res*0.2, 0.5/surface_res) if surface_mode == 'A' else 0.3
    set_node_values(_mesh2vol, [(2,voxel_size)])
    _vol2mesh_SAS = add_node(group, 'GeometryNodeVolumeToMesh', (a+800,b), '')
    _vol2mesh_SAS.resolution_mode = 'VOXEL_SIZE'
    set_node_values(_vol2mesh_SAS, [(1,voxel_size*1.1)])
    nodes_link(group, _geometry, socket_out, _instance, 0)
    nodes_link(group, _ico_sphere, 0, _instance, 2)
    nodes_link(group, _attr_radii_vdW, v.na_sockets_out['float'], _add, 0)
    nodes_link(group, _add, 0, _instance, 6)
    nodes_link(group, _instance, 0, _realize_instance, 0)
    nodes_link(group, _realize_instance, 0, _mesh2vol, 0)
    nodes_link(group, _mesh2vol, 0, _vol2mesh_SAS, 0)

    _sample_nearest = add_node(group, 'GeometryNodeSampleNearest', (a+1000,b+300), '')
    _sample_idx_I = add_sample_index_node(group, (a+1200,b+250), '', 'FLOAT_VECTOR', 'POINT')
    _sample_idx_II = add_sample_index_node(group, (a+1200,b+500), '', 'FLOAT_VECTOR', 'POINT')
    _input_pos = add_node(group, 'GeometryNodeInputPosition', (a+1000,b+500), '')
    _input_normal = add_node(group, 'GeometryNodeInputNormal', (a+1000,b+100), '')
    _boundbox = add_node(group, 'GeometryNodeBoundBox', (a+1200,b), '')
    _subtract = add_vec_math_node(group, (a+1400,b+700), '', [], 'SUBTRACT')
    _dot_product = add_vec_math_node(group, (a+1400,b+500), '', [], 'DOT_PRODUCT')
    _distance = add_vec_math_node(group, (a+1400,b+250), '', [], 'DISTANCE')
    _less = add_compare_node(group, (a+1600,b+500), '', [], 'FLOAT', 'LESS_THAN')
    _greater_equal = add_compare_node(group, (a+1600,b+250), '', [(v.cp_sockets_in['floatB'],probe_radius)], 'FLOAT', 'GREATER_EQUAL')
    _boolean_AND = add_node(group, 'FunctionNodeBooleanMath', (a+1800,b+400), '')
    _boolean_AND.operation = 'AND'
    _volume_cube = add_node(group, 'GeometryNodeVolumeCube', (a+2000,b+50), '')
    _subtract_1 = add_vec_math_node(group, (a+1400,b-150), '', [], 'SUBTRACT')
    _scale = add_vec_math_node(group, (a+1600,b-150), '', [(3,surface_res)], 'SCALE')
    _separate_XYZ = add_node(group, 'ShaderNodeSeparateXYZ', (a+1800,b-150), '')
    nodes_link(group, _vol2mesh_SAS, 0, _sample_nearest, 0)
    nodes_link(group, _vol2mesh_SAS, 0, _sample_idx_I, 0)
    nodes_link(group, _vol2mesh_SAS, 0, _sample_idx_II, 0)
    nodes_link(group, _vol2mesh_SAS, 0, _boundbox, 0)
    nodes_link(group, _input_pos, 0, _sample_idx_I, v.si_sockets_in['vec'])
    nodes_link(group, _input_pos, 0, _subtract, 0)
    nodes_link(group, _input_pos, 0, _distance, 0)
    nodes_link(group, _sample_nearest, 0, _sample_idx_I, v.si_sockets_in['index'])
    nodes_link(group, _input_normal, 0, _sample_idx_II, v.si_sockets_in['vec'])
    nodes_link(group, _sample_nearest, 0, _sample_idx_II, v.si_sockets_in['index'])
    nodes_link(group, _sample_idx_I, v.si_sockets_out['vec'], _subtract, 1)
    nodes_link(group, _sample_idx_I, v.si_sockets_out['vec'], _distance, 1)
    nodes_link(group, _subtract, 0, _dot_product, 0)
    nodes_link(group, _sample_idx_II, v.si_sockets_out['vec'], _dot_product, 1)
    nodes_link(group, _dot_product, 1, _less, 0)
    nodes_link(group, _distance, 1, _greater_equal, 0)
    nodes_link(group, _less, 0, _boolean_AND, 0)
    nodes_link(group, _greater_equal, 0, _boolean_AND, 1)
    nodes_link(group, _boolean_AND, 0, _volume_cube, 0)
    nodes_link(group, _boundbox, 2, _subtract_1, 0)
    nodes_link(group, _boundbox, 1, _subtract_1, 1)
    nodes_link(group, _subtract_1, 0, _scale, 0)
    nodes_link(group, _scale, 0, _separate_XYZ, 0)
    nodes_link(group, _boundbox, 1, _volume_cube, 2)
    nodes_link(group, _boundbox, 2, _volume_cube, 3)
    nodes_link(group, _separate_XYZ, 0, _volume_cube, 4)
    nodes_link(group, _separate_XYZ, 1, _volume_cube, 5)
    nodes_link(group, _separate_XYZ, 2, _volume_cube, 6)

    _vol2mesh_SES = add_node(group, 'GeometryNodeVolumeToMesh', (a+2200,b), '')
    _repeat_input = add_node(group, 'GeometryNodeRepeatInput', (a+2400,b), '')
    set_node_values(_repeat_input, [(0,smooth_iter)])   # repeat iterations
    _repeat_output = add_node(group, 'GeometryNodeRepeatOutput', (a+2800,b), '')
    _repeat_input.pair_with_output(_repeat_output)
    _set_pos = add_node(group,'GeometryNodeSetPosition', (a+2600,b), '')
    _edge_vertices = add_node(group, 'GeometryNodeInputMeshEdgeVertices', (a+2100,b-220), '')
    _add = add_vec_math_node(group, (a+2300,b-220), '', [], 'ADD')
    _scale = add_vec_math_node(group, (a+2500,b-220), '', [(3,0.5)], 'SCALE')
    _shade_smooth = add_node(group, 'GeometryNodeSetShadeSmooth', (a+3000, b), '')
    _output.location = (a+3200,0)
    nodes_link(group, _volume_cube, 0, _vol2mesh_SES, 0)
    _vol2mesh = _vol2mesh_SAS if surface_mode == 'A' else _vol2mesh_SES
    nodes_link(group, _vol2mesh, 0, _repeat_input, 1)
    nodes_link(group, _repeat_input, 0, _set_pos, 0)
    nodes_link(group, _set_pos, 0, _repeat_output, 0)
    nodes_link(group, _edge_vertices, 2, _add, 0)
    nodes_link(group, _edge_vertices, 3, _add, 1)
    nodes_link(group, _add, 0, _scale, 0)
    nodes_link(group, _scale, 0, _set_pos, 2)
    nodes_link(group, _repeat_output, 0, _shade_smooth, 0)

    a += 3000
    _sample_nearest = add_node(group, 'GeometryNodeSampleNearest', (a,b+200), '')
    _sample_idx = add_sample_index_node(group, (a+200,b), '', 'FLOAT_COLOR', 'POINT')
    if color_type == 'C':
        attr_name = 'chain_id'
        datatype = 'INT'
        socket_type = 'int'
        mat = bpy.data.materials['ChainID']
    elif color_type == 'H':
        attr_name = 'colour'
        datatype = 'FLOAT_COLOR'
        socket_type = 'color'
        mat = bpy.data.materials['Hetero']
    elif color_type == 'R':
        attr_name = 'index'
        datatype = 'INT'
        socket_type = 'int'
        mat = bpy.data.materials['Rainbow']
    elif color_type == 'E':
        attr_name = 'charge'
        datatype = 'FLOAT'
        socket_type = 'float'
        mat = bpy.data.materials['Charge']
    elif color_type == 'L':
        attr_name = 'lipophobicity'
        datatype = 'FLOAT'
        socket_type = 'float'
        mat = bpy.data.materials['Hydrophobic']

    _attr_color = add_attr_node(group, (a,b-200), '', [(0,attr_name)], datatype)
    _blur_attr = add_node(group, 'GeometryNodeBlurAttribute', (a+400,b), '')
    set_node_datatype(_blur_attr, 'FLOAT_COLOR')
    set_node_values(_blur_attr, [(v.ba_sockets_in['iter'],blur_iteration)])
    _store_attr_color = add_store_attr_node(group, (a+600,b), '', 'FLOAT_COLOR', 'POINT', attr_name)
    _material = add_mat_node(group, (a+800,b), '', mat)
    _output.location = (a+1000,b)

    nodes_link(group, _backbone, 0, _sample_nearest, 0)
    nodes_link(group, _backbone, 0, _sample_idx, 0)
    nodes_link(group, _sample_nearest, 0, _sample_idx, v.si_sockets_in['index'])
    nodes_link(group, _attr_color, v.na_sockets_out[socket_type], _sample_idx, v.si_sockets_in['color'])
    nodes_link(group, _sample_idx, v.si_sockets_out['color'], _blur_attr, v.ba_sockets_in['color'])
    nodes_link(group, _blur_attr, v.ba_sockets_out['color'], _store_attr_color, v.sa_sockets_in['color'])
    nodes_link(group, _store_attr_color, 0, _material, 0)
    nodes_link(group, _shade_smooth, 0, _store_attr_color, 0)
    nodes_link(group, _material, 0, _output, 0)

