import bpy
import colorsys
import numpy as np
from . import _node
from .Chem_Data import ELEMENTS_DEFAULT

language = 1 if "zh_HAN" in bpy.context.preferences.view.language else 0
version = bpy.app.version
if language == 0:
    bsdf = "Principled BSDF"
else:
    if version[:2] == (4,0) or version[:2] == (4,1):
        bsdf = "原理化BSDF"
    else:
        bsdf = "原理化 BSDF"


def create_material(mat_name):
    FLAG = False
    for material in bpy.data.materials:
        if material.name == mat_name:
            mat = material
            FLAG = True
            break
    if not FLAG:
        mat = bpy.data.materials.new(name=mat_name)
        mat.use_nodes = True
        if mat_name in list(ELEMENTS_DEFAULT.keys()):
            _BSDF = mat.node_tree.nodes[bsdf]
            _BSDF.inputs['Base Color'].default_value = ELEMENTS_DEFAULT[mat_name][3]  # Color
            _BSDF.inputs['Roughness'].default_value = 0.2   # Roughness
    return mat

def create_mat_hetero(mat_name):
    FLAG = False
    for material in bpy.data.materials:
        if material.name == mat_name:
            mat = material
            FLAG = True
            break
    if not FLAG:
        mat = bpy.data.materials.new(name=mat_name)
        mat.use_nodes = True
        _BSDF = mat.node_tree.nodes[bsdf]
        _BSDF.inputs[2].default_value = 0.2   # Roughness
        _attr_color = _node.add_node(mat.node_tree, 'ShaderNodeAttribute', (-300,300), '')
        _attr_color.attribute_name = 'colour'
        link = mat.node_tree.links.new
        link(_attr_color.outputs[0], _BSDF.inputs[0])
    return mat

# Cartoon material
def cartoon_mat(color_helix, color_sheet, color_coil):
    cartoon = create_material('Cartoon')
    cartoon.use_nodes = True
    _principled_bsdf = cartoon.node_tree.nodes['Principled BSDF']
    _color_ramp = _node.add_node(cartoon.node_tree, 'ShaderNodeValToRGB', (-300,300), '')
    _color_ramp.color_ramp.interpolation = 'CONSTANT'
    color_list = (color_helix, color_sheet, color_coil)
    _color_ramp.color_ramp.elements.new(0.667)
    for i in range(3):
        _color_ramp.color_ramp.elements[i].position = i*0.333
        _color_ramp.color_ramp.elements[i].color = color_list[i]
    
    _map_range = _node.add_node(cartoon.node_tree, 'ShaderNodeMapRange', (-500,300), 'Range')
    _node.set_node_values(_map_range, [(1,1),(2,3)])
    _attr_sec_struct = _node.add_node(cartoon.node_tree, 'ShaderNodeAttribute', (-700,300),'')
    _attr_sec_struct.attribute_name = 'sec_struct'

    link = cartoon.node_tree.links.new
    link(_attr_sec_struct.outputs[2], _map_range.inputs[0])
    link(_map_range.outputs[0], _color_ramp.inputs[0])
    link(_color_ramp.outputs[0], _principled_bsdf.inputs[0])
    return cartoon


def color_by_id(id_list):
    ids = np.unique(id_list)
    num_colors = len(ids)
    hues = [i/(num_colors+1) for i in range(num_colors)]
    # convert HSL to RGB
    colors = [colorsys.hsv_to_rgb(hue, 0.8, 0.8) for hue in hues]
    colors = [(r,g,b,1) for (r,g,b) in colors]
    return colors

def create_mat_id(id_list, mat_name, id_name):
    mat = create_material(mat_name)
    mat.use_nodes = True
    _BSDF = mat.node_tree.nodes[bsdf]
    _color_ramp = _node.add_node(mat.node_tree, 'ShaderNodeValToRGB', (-300,300), '')
    _color_ramp.color_ramp.interpolation = 'CONSTANT'

    colors = color_by_id(id_list)
    if len(colors) > 2:
        for i in range(len(colors)-2):
            _color_ramp.color_ramp.elements.new(0)
    for i in range(len(colors)):
        _color_ramp.color_ramp.elements[i].position = i/len(colors)
        _color_ramp.color_ramp.elements[i].color = colors[i]
    _map_range = _node.add_node(mat.node_tree, 'ShaderNodeMapRange', (-500,300), 'Range')
    _node.set_node_values(_map_range, [(1,min(id_list)),(2,max(id_list))])
    _attr = _node.add_node(mat.node_tree, 'ShaderNodeAttribute', (-700,300),'')
    _attr.attribute_name = id_name

    link = mat.node_tree.links.new
    link(_attr.outputs[2], _map_range.inputs[0])
    link(_map_range.outputs[0], _color_ramp.inputs[0])
    link(_color_ramp.outputs[0], _BSDF.inputs[0])
    return mat


def create_mat_rainbow(id_list):
    mat = create_material('Rainbow')
    mat.use_nodes = True
    _BSDF = mat.node_tree.nodes[bsdf]
    _color_ramp = _node.add_node(mat.node_tree, 'ShaderNodeValToRGB', (-300,300), '')
    _color_ramp.color_ramp.color_mode = 'HSV'
    _color_ramp.color_ramp.hue_interpolation = 'FAR'
    _color_ramp.color_ramp.elements[0].color = (1.0, 0.1, 0.1, 1.0)
    _color_ramp.color_ramp.elements[1].color = (0.115, 0.0, 1.0, 1.0)
    _map_range = _node.add_node(mat.node_tree, 'ShaderNodeMapRange', (-500,300), 'Range')
    _node.set_node_values(_map_range, [(1,min(id_list)),(2,max(id_list))])
    _attr = _node.add_node(mat.node_tree, 'ShaderNodeAttribute', (-700,300),'')
    _attr.attribute_name = 'index'

    link = mat.node_tree.links.new
    link(_attr.outputs[2], _map_range.inputs[0])
    link(_map_range.outputs[0], _color_ramp.inputs[0])
    link(_color_ramp.outputs[0], _BSDF.inputs[0])
    return mat


# Eletrostatic material
def create_mat_charge(charge):
    mat = create_material('Charge')
    mat.use_nodes = True
    _BSDF = mat.node_tree.nodes[bsdf]
    _color_ramp = _node.add_node(mat.node_tree, 'ShaderNodeValToRGB', (-300,300), '')
    _color_ramp.color_ramp.elements.new(0)
    _color_ramp.color_ramp.elements[0].color = (1.0, 0.0, 0.0, 1.0)
    _color_ramp.color_ramp.elements[0].position = 0.0
    _color_ramp.color_ramp.elements[1].color = (1.0, 1.0, 1.0, 1.0)
    _color_ramp.color_ramp.elements[1].position = 0.5
    _color_ramp.color_ramp.elements[2].color = (0.0, 0.0, 1.0, 1.0)
    _color_ramp.color_ramp.elements[2].position = 1.0
    _map_range = _node.add_node(mat.node_tree, 'ShaderNodeMapRange', (-500,300), 'Range')
    _node.set_node_values(_map_range, [(1,min(charge)),(2,max(charge))])
    # _node.set_node_values(_map_range, [(1,-1),(2,1)])
    _attr = _node.add_node(mat.node_tree, 'ShaderNodeAttribute', (-700,300),'')
    _attr.attribute_name = 'charge'

    link = mat.node_tree.links.new
    link(_attr.outputs[2], _map_range.inputs[0])
    link(_map_range.outputs[0], _color_ramp.inputs[0])
    link(_color_ramp.outputs[0], _BSDF.inputs[0])
    return mat



# Hydrophobic material
def create_mat_hydrophobic(lipophobicity):
    mat = create_material('Hydrophobic')
    mat.use_nodes = True
    _BSDF = mat.node_tree.nodes[bsdf]
    _color_ramp = _node.add_node(mat.node_tree, 'ShaderNodeValToRGB', (-300,300), '')
    _color_ramp.color_ramp.elements.new(0)
    _color_ramp.color_ramp.elements[0].color = (0.0, 0.5, 1.0, 1.0)
    _color_ramp.color_ramp.elements[0].position = 0.0
    _color_ramp.color_ramp.elements[1].color = (1.0, 1.0, 1.0, 1.0)
    _color_ramp.color_ramp.elements[1].position = 0.5
    _color_ramp.color_ramp.elements[2].color = (1.0, 0.5, 0.0, 1.0)
    _color_ramp.color_ramp.elements[2].position = 1.0
    _map_range = _node.add_node(mat.node_tree, 'ShaderNodeMapRange', (-500,300), 'Range')
    _node.set_node_values(_map_range, [(1,min(lipophobicity)),(2,max(lipophobicity))])
    # _node.set_node_values(_map_range, [(1,-1),(2,1)])
    _attr = _node.add_node(mat.node_tree, 'ShaderNodeAttribute', (-700,300),'')
    _attr.attribute_name = 'lipophobicity'

    link = mat.node_tree.links.new
    link(_attr.outputs[2], _map_range.inputs[0])
    link(_map_range.outputs[0], _color_ramp.inputs[0])
    link(_color_ramp.outputs[0], _BSDF.inputs[0])
    return mat