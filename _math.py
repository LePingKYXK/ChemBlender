import numpy as np
from math import cos,acos,sin,asin,atan2,sqrt
from mathutils import Vector, Matrix

# ------------------------------------------------------------------------------------
def get_length(vec):
    length = np.linalg.norm(vec)
    return length

def normalize(vec):
    return vec / get_length(vec) if get_length(vec) else vec

# rotate a vector around an axis
def rotate_vec(vec, axis, angle):
    u = normalize(axis)
    x, y, z = u[0], u[1], u[2]
    cos = np.cos(angle)
    sin = np.sin(angle)
    R = np.array([[cos+x*x*(1-cos), x*y*(1-cos)-z*sin, x*z*(1-cos)+y*sin],
                  [x*y*(1-cos)+z*sin, cos+y*y*(1-cos), y*z*(1-cos)-x*sin],
                  [x*z*(1-cos)-y*sin, y*z*(1-cos)+x*sin, cos+z*z*(1-cos)],])
    return np.dot(R, vec)

# calculate a vertical vector of the vec
def vertical_vec(pos1, pos2):
    v = normalize(pos2-pos1)
    if v[0] == 0 and v[1] == 0:
        vertical = np.array((1,0,0))
    else:
        vertical = normalize((-v[1], v[0], 0))
    return vertical

# calculate normal vector from three points in a plane
def get_normal(pos1, pos2, pos3):
    x1,y1,z1 = pos1
    x2,y2,z2 = pos2
    x3,y3,z3 = pos3
    a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
    b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
    c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
    return normalize((a,b,c))

# get the foot point from the known point to a line
def getFootPoint(pos,line_p1,line_p2):
    x0, y0, z0 = pos
    x1, y1, z1 = line_p1
    x2, y2, z2 = line_p2
    if (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 != 0:
        k = -((x1-x0)*(x2-x1) + (y1-y0)*(y2-y1) + (z1-z0)*(z2-z1))/\
            ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    else:
        k = 0
    xn = k*(x2-x1)+x1
    yn = k*(y2-y1)+y1
    zn = k*(z2-z1)+z1
    return (xn,yn,zn)

# This function calculates euler angle when rotate one vector to another.
def Euler_angle(vec0, vec1):    
    theta = acos(np.dot(vec0, vec1))
    crs = np.cross(vec0, vec1)
    n = crs / np.linalg.norm(crs)
    n_invert = np.array(([0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]))
    R = np.eye(3) + sin(theta)*n_invert + np.dot(n_invert,n_invert)*(1-cos(theta))
    alpha = atan2(R[1][0],R[0][0])
    beta = asin(-R[2][0])
    gamma = atan2(R[2][1], R[2][2])
    if abs(np.dot(vec1,vec0)) == 1.0:
        return vec0
    else:
        return (gamma,beta,alpha)

# A Simple One. If using this, vec1 should be an Vector data.
def simple_Euler_angle(vec1):
    vec0 = Vector([0.0, 0.0, 1.0])
    # Angle with respect to the z-axis
    angle = vec1.angle(vec0, 0)
    # Cross-product between vec1 and vec0. It is the vector of rotation.
    axis = vec0.cross(vec1)
    # Calculate Euler angles
    euler = Matrix.Rotation(angle,4,axis).to_euler()
    return euler


# ----------------------------------------------------------------------------------------

# get symmetry operation component from a string
# ref. International Tables for Crystallography (2006). Vol. A, Section 11.1.1, p.810.
# 11.1. Point coordinates, symmetry operations and their symbols
# BY W.FISCHER AND E.KOCH
def convert_symop_xyz_to_vec(string):  # eg. string = "3/4-x"
    rotate_comp = []
    for symbol in ['x','y','z']:
        i = string.find(symbol)
        value = 0.0 if i == -1 else -1.0 if string[i-1]=='-' else 1.0
        rotate_comp.append(value)
    i = string.find("/")
    transl_comp = 0.0 if i == -1 else float(string[max(i-2,0):i])/float(string[i+1:i+2])
    return(rotate_comp, transl_comp)
    
# get symmetry operation matrices (W, w) from strings. W is the rotation part, and w is the
# translation part.
def symop_xyz_to_matrix(strings):    # eg. strings = "-x+y, 1/2+y, -z-1/2"
    rotate_matrix = [convert_symop_xyz_to_vec(str)[0] for str in strings.split(",")]
    transl_matrix = [convert_symop_xyz_to_vec(str)[1] for str in strings.split(",")]
    return (rotate_matrix, transl_matrix)

# from one fract_xyz to multi equiv fracts in one cell through symmetry operations.
def fract_symop(fract_xyz, symop_operations):
    sym_fracts = []
    for symop in symop_operations:
        sym_matrix = symop_xyz_to_matrix(symop)
        npS = np.array(sym_matrix[0])
        npT = np.array(sym_matrix[1])
        npA = np.array(fract_xyz)
        npB = np.matmul(npS,npA) + npT    # Applying symmetry and translation.
        sym_fracts.append(npB)
    return sym_fracts

def fract_to_cartn(fract_xyz,a,b,c,alpha,beta,gamma):
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)
    v = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    M = np.array([[a, b*cos(gamma), c*cos(beta)],
                [0, b*sin(gamma), c*v],
                [0, 0, c*sqrt(pow(sin(beta),2)-v*v)]])
    cartn_xyz = np.matmul(M, np.array(fract_xyz))
    return cartn_xyz

def cartn_to_fract(cartn_xyz,a,b,c,alpha,beta,gamma):
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)
    v = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    M = np.array([[a, b*cos(gamma), c*cos(beta)],
                    [0, b*sin(gamma), c*v],
                    [0, 0, c*sqrt(pow(sin(beta),2)-v*v)]])
    inv_M = np.linalg.inv(M)
    fract_xyz = np.matmul(inv_M, np.array(cartn_xyz))
    return fract_xyz

def fracts_normalize(sym_fracts, boundary):
    # coordinates normalization (e.g. 1.15 becomes 0.15, -0.25 becomes 0.75, etc)
    for i in range(len(sym_fracts)):
        if sym_fracts[i][0]<-boundary: sym_fracts[i] = np.array(sym_fracts[i]) + np.array((1,0,0))
        if sym_fracts[i][0]>1+boundary: sym_fracts[i] = np.array(sym_fracts[i]) - np.array((1,0,0))
        if sym_fracts[i][1]<-boundary: sym_fracts[i] = np.array(sym_fracts[i]) + np.array((0,1,0))
        if sym_fracts[i][1]>1+boundary: sym_fracts[i] = np.array(sym_fracts[i]) - np.array((0,1,0))
        if sym_fracts[i][2]<-boundary: sym_fracts[i] = np.array(sym_fracts[i]) + np.array((0,0,1))
        if sym_fracts[i][2]>1+boundary: sym_fracts[i] = np.array(sym_fracts[i]) - np.array((0,0,1))
    
    # remove duplicates (same position)
    sym_fracts = list(set([tuple(fract) for fract in sym_fracts]))
    new_sym_fracts = sym_fracts
    # boundary = 0.0  # the value is between 0.0 and 1.0
    for fract in sym_fracts:
        if fract[0] <= boundary and fract[0] >= -boundary:
            new_fract = (fract[0]+1, fract[1], fract[2])
            new_sym_fracts.append(new_fract)
        if fract[1] <= boundary and fract[0] >= -boundary:
            new_fract = (fract[0], fract[1]+1, fract[2])
            new_sym_fracts.append(new_fract)
        if fract[2] <= boundary and fract[0] >= -boundary:
            new_fract = (fract[0], fract[1], fract[2]+1)
            new_sym_fracts.append(new_fract)
    
    for fract in sym_fracts:
        if fract[0] >= 1-boundary and fract[0] <= 1+boundary:
            new_fract = (fract[0]-1, fract[1], fract[2])
            new_sym_fracts.append(new_fract)
        if fract[1] >= 1-boundary and fract[0] <= 1+boundary:
            new_fract = (fract[0], fract[1]-1, fract[2])
            new_sym_fracts.append(new_fract)
        if fract[2] >= 1-boundary and fract[0] <= 1+boundary:
            new_fract = (fract[0], fract[1], fract[2]-1)
            new_sym_fracts.append(new_fract)
    
    sym_fracts = list(set([tuple(fract) for fract in new_sym_fracts]))
    print(sym_fracts)
    return sym_fracts


    