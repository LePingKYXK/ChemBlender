import bpy
from bpy.props import StringProperty, EnumProperty

# ------------------------------------------------------------------------------------
# This class lists properties shown in Panel but not in Operator.
class Mol3DSettings(bpy.types.PropertyGroup):
    datafile: StringProperty(
        name = '',
        default = '请选择一个分子文件' if 'zh_HAN' in bpy.context.preferences.view.language else 'Choose Your Molecule File',
        subtype = 'FILE_PATH'
        ) # type: ignore
    
    select_text: StringProperty(
        name = '',
        default = "输入文本" if 'zh_HAN' in bpy.context.preferences.view.language else "Input Text Here",
        description = "Text of Element Name for Scaling"
    ) # type:ignore

    distance: StringProperty(
        name = '',
        default = "距离 (Å)" if 'zh_HAN' in bpy.context.preferences.view.language else "Distance (Å)",
        description = "Length of edge in Angstrom"
    ) # type:ignore

    angle: StringProperty(
        name = '',
        default = "角度 (°)" if 'zh_HAN' in bpy.context.preferences.view.language else "Angle (°)",
        description = "Angle of two edges in degree"
    ) # type:ignore

    avgfract: StringProperty(name='', default='x, y, z') # type:ignore
    deformation: EnumProperty(
        name = '',
        default = 'Roll',
        items = [('Wave','波浪',''),
                 ('Roll','卷曲',''),
                 ('Noise','噪波','')]
                if 'zh_HAN' in bpy.context.preferences.view.language else 
                [('Wave','Wave',''),
                 ('Roll','Roll',''),
                 ('Noise','Noise','')]
    ) # type:ignore


# This class to draw layout of the panel
class MOL3D_PT_Build(bpy.types.Panel):
    bl_label = "分子结构创建" if 'zh_HAN' in bpy.context.preferences.view.language else "Mol3D Structure Building"
    bl_idname = "MOL3D_PT_BUILD"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "ChemBlender"
        
    def draw(self, context):
        mytool = context.scene.my_tool
        layout = self.layout
        layout.prop(mytool, "datafile")
        
        layout.row()
        text = "创建结构" if 'zh_HAN' in bpy.context.preferences.view.language else "Build Structure"
        layout.operator("mol3d.structure_build", text=text, icon="GREASEPENCIL")
        text = "智能显示" if 'zh_HAN' in bpy.context.preferences.view.language else "Smart Display"
        layout.operator("mol3d.mol_display", text=text, icon="OUTLINER_OB_LIGHT")
        

class MOL3D_PT_TOOLS(bpy.types.Panel):
    bl_label = "分子工具" if 'zh_HAN' in bpy.context.preferences.view.language else "Mol3D Tools"
    bl_idname = "MOL3D_PT_UTILITY"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "ChemBlender"
    bl_options = {'DEFAULT_CLOSED'}
        
    def draw(self, context):
        layout = self.layout
        mytool = context.scene.my_tool

        row = layout.row(align=True)
        row.prop(mytool, "select_text")
        row.scale_x = 0.75
        text = "选择" if 'zh_HAN' in bpy.context.preferences.view.language else "Select"
        row.operator('mol3d.select', text = text)
        row = layout.row(align=True)
        text = "缩放" if 'zh_HAN' in bpy.context.preferences.view.language else "Scale"
        row.operator("mol3d.scale", text = text)
        row.scale_x = 1.5
        text = "默认尺寸" if 'zh_HAN' in bpy.context.preferences.view.language else "Default Size"
        row.operator("mol3d.default_size", text = text)
        row = layout.row(align=True)
        text = "显示标注" if 'zh_HAN' in bpy.context.preferences.view.language else 'Show Legends'
        row.operator("mol3d.add_legends", text = text)
        
        # text = "显示连接/取消显示" if 'zh_HAN' in bpy.context.preferences.view.language else "Connect Display Toggle"
        # layout.operator("mol3d.showconnect", text=text, icon="MOD_LENGTH")
        box = layout.box()
        row = box.row(align=True)
        text = "测量长度" if 'zh_HAN' in bpy.context.preferences.view.language else "Length Calc."
        row.operator("mol3d.button_distance", text=text)
        row.prop(mytool,"distance")
        row = box.row(align=True)
        text = "测量角度" if 'zh_HAN' in bpy.context.preferences.view.language else "Angle Calc"
        row.operator("mol3d.button_angle", text=text)
        row.prop(mytool,"angle")
        

class CRYSTAL_PT_TOOLS(bpy.types.Panel):
    bl_label = "晶体工具" if 'zh_HAN' in bpy.context.preferences.view.language else "Crystal3D Tools"
    bl_idname = "MOL3D_PT_CRYSTAL"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "ChemBlender"
    bl_options = {'DEFAULT_CLOSED'}
        
    def draw(self, context):
        layout = self.layout
        mytool = context.scene.my_tool
        row = layout.row(align=True)
        text = "创建多面体" if 'zh_HAN' in bpy.context.preferences.view.language else "Build Polyhedra"
        row.operator("mol3d.polyhedra_build", text=text)
        row.scale_x = 0.4
        text = "编辑" if 'zh_HAN' in bpy.context.preferences.view.language else "Edit"
        row.operator("mol3d.polyhedra_edit", text=text)

        row = layout.row()
        text = "剖面显示" if 'zh_HAN' in bpy.context.preferences.view.language else "Cross Section"
        row.operator("mol3d.crystal_section", text=text)

        box = layout.box()
        row = box.row(align=True)
        text = "平均位置分数" if 'zh_HAN' in bpy.context.preferences.view.language else "AVG. Fracts"
        row.operator("mol3d.button_avgfract", text=text)
        row.scale_x = 0.65
        row.prop(mytool, "avgfract")
        row = box.row(align=True)
        text = "创建Dummy" if 'zh_HAN' in bpy.context.preferences.view.language else "Create Dummy"
        row.operator("mol3d.button_dummy", text=text, icon="ALIASED")
        # text = "结构变形：" if 'zh_HAN' in bpy.context.preferences.view.language else "Structure Deformation:"
        # layout.label(text = text)
        box = layout.box()
        row = box.row()
        text = "结构变形：" if 'zh_HAN' in bpy.context.preferences.view.language else "Deform Mode:"
        row.label(text = text)
        row.scale_x = 1.2
        row.prop(mytool, "deformation")
        row = box.row()
        text = "更新" if 'zh_HAN' in bpy.context.preferences.view.language else "Update"
        row.operator("mol3d.deformation_add", text = text)
        text = "移除" if 'zh_HAN' in bpy.context.preferences.view.language else "Remove"
        row.operator("mol3d.deformation_remove", text = text)

class CUSTOM_PT_DRAW(bpy.types.Panel):
    bl_label = "自定义绘制" if 'zh_HAN' in bpy.context.preferences.view.language else "Custom Drawing"
    bl_idname = "MOL3D_PT_CUSTOM_DRAW"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "ChemBlender"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        text = "设置原子" if 'zh_HAN' in bpy.context.preferences.view.language else "Set Atoms"
        row.operator("mol3d.set_atoms", text = text)
        text = "设置键" if 'zh_HAN' in bpy.context.preferences.view.language else "Set Bonds"
        row.operator("mol3d.set_bonds", text = text)

        row = layout.row()
        text = "添加新结构" if 'zh_HAN' in bpy.context.preferences.view.language else "New Struct"
        row.operator("mol3d.add_struct", text = text)
        text = "添加晶胞" if 'zh_HAN' in bpy.context.preferences.view.language else "New Cell"
        row.operator("mol3d.add_unit_cell", text = text)

        row = layout.row()
        text = "构象优化" if 'zh_HAN' in bpy.context.preferences.view.language else "MolOptimize"
        row.operator("mol3d.mmff_optimize", text = text)
        text = "构象更新" if 'zh_HAN' in bpy.context.preferences.view.language else "MolUpdate"
        row.operator("mol3d.mmff_update", text = text)

        