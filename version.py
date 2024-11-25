import bpy
version = bpy.app.version
if version[:2] == (4,0):
    na_sockets_out = {  # named attribute
        'vec':0, 'float':1, 'color':2, 'bool':3, 'int':4, 'quater':5
    }
    si_sockets_in = {  # sample index
        'float':1, 'int':2, 'vec':3, 'color':4, 'bool':5, 'quater':6, 'index':7
    }
    si_sockets_out = {  # sample index
        'float':0, 'int':1, 'vec':2, 'color':3, 'bool':4, 'quater':5
    }
    sa_sockets_in = {  # store attribute
        'vec':3, 'float':4, 'color':5, 'bool':6, 'int':7, 'quater':8
    }
    cp_sockets_in = {  # compare
        'floatA':0, 'floatB':1, 'intA':2, 'intB':3, 'vecA':4, 'vecB':5, 'colorA':6, 'colorB':7, 'strA':8, 'strB':9
    }
    # from float to color: switch;   from object to material: switch_1
    sw_sockets_in = {  # switch
        'switch':0, 'switch_1':1, 'floatA':2, 'floatB':3, 'intA':4, 'intB':5, 'boolA':6, 'boolB':7, 'vecA':8, 'vecB':9,
        'colorA':10, 'colorB':11, 'strA':12, 'strB':13, 'geoA':14, 'geoB':15, 'objA':16, 'objB':17,
        'collA':18, 'collB':19, 'textA':20, 'textB':21, 'matA':22, 'matB':23, 'imageA':24, 'imageB':25,
        'rotA':26, 'rotB':27
    }
    sw_sockets_out = {  # switch
        'float':0, 'int':1, 'bool':2, 'vec':3, 'color':4, 'str':5, 'geo':6, 'obj':7, 'coll':8, 'text':9, 'mat':10, 'image':11, 'rot':12
    }
    rc_sockets_in = {  # raycast
        'dir':8,
    }
    fi_sockets_in = {  # field at index
        'float':1, 'int':2, 'vec':3, 'color':4, 'bool':5, 'quater':6
    }
    fi_sockets_out = {  # field at index
        'float':0, 'int':1, 'vec':2, 'color':3, 'bool':4, 'quater':5
    }
    ba_sockets_in = {  # blur attribute
        'float':0, 'int':1, 'vec':2, 'color':3, 'iter':4, 
    }
    ba_sockets_out = {  # blur attribute
        'float':0, 'int':1, 'vec':2, 'color':3 
    }
else:
    na_sockets_out = {  # named attribute
        'vec':0, 'float':0, 'color':0, 'bool':0, 'int':0, 'quater':0
    }
    si_sockets_in = {  # sample index
        'float':1, 'int':1, 'vec':1, 'color':1, 'bool':1, 'quater':1, 'index':2
    }
    si_sockets_out = {  # sample index
        'float':0, 'int':0, 'vec':0, 'color':0, 'bool':0, 'quater':0
    }
    sa_sockets_in = {  # store attribute
        'vec':3, 'float':3, 'color':3, 'bool':3, 'int':3, 'quater':3
    }
    cp_sockets_in = {  # compare
        'floatA':0, 'floatB':1, 'intA':2, 'intB':3, 'vecA':4, 'vecB':5, 'colorA':6, 'colorB':7, 'strA':8, 'strB':9
    }
    # from float to color: switch;   from object to material: switch_1
    sw_sockets_in = {  # switch
        'switch':0, 'switch_1':0, 'floatA':1, 'floatB':2, 'intA':1, 'intB':2, 'boolA':1, 'boolB':2, 'vecA':1, 'vecB':2,
        'colorA':1, 'colorB':2, 'strA':1, 'strB':2, 'geoA':1, 'geoB':2, 'objA':1, 'objB':2,
        'collA':1, 'collB':2, 'textA':1, 'textB':2, 'matA':1, 'matB':2, 'imageA':1, 'imageB':2,
        'rotA':1, 'rotB':2
    }
    sw_sockets_out = {
        'float':0, 'int':0, 'bool':0, 'vec':0, 'color':0, 'str':0, 'geo':0, 'obj':0, 'coll':0, 'text':0, 'mat':0, 'image':0, 'rot':0
    }
    rc_sockets_in = {  # raycast
        'dir':3,
    }
    fi_sockets_in = {  # field at index
        'float':1, 'int':1, 'vec':1, 'color':1, 'bool':1, 'quater':1
    }
    fi_sockets_out = {  # field at index
        'float':0, 'int':0, 'vec':0, 'color':0, 'bool':0, 'quater':0
    }
    fd_sockets_in = {  # field on domain
        'float':0, 'int':0, 'vec':0, 'color':0, 'bool':0, 'quater':0
    }
    fd_sockets_out = {  # field on domian
        'float':0, 'int':0, 'vec':0, 'color':0, 'bool':0, 'quater':0
    }
    ba_sockets_in = {  # blur attribute
        'float':0, 'int':0, 'vec':0, 'color':0, 'iter':1, 
    }
    ba_sockets_out = {  # blur attribute
        'float':0, 'int':0, 'vec':0, 'color':0 
    }