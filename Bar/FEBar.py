#!/usr/bin/env python3
"""
@mainpage FEBar

This is the main page for FEBar

@section Goals

A simple 2D bar finite element

@section Todo

@todo
 - 19.04.11 - everything


@section Teste



"""

import sys  # os, codecs , fileinput
# import string
# import math
import numpy as np

from numpy.linalg import norm, solve
from math import radians #, pi, cos, sin, sqrt,


# classes locais -------------
from Data.class_data   import *
from Util.utilities import align_header
from math import atan2

from Data.class_results import *


    

def compute_K_bar(A,E,L ):
    
    
        
    K_elem_bar = np.matrix([[1,-1],[-1,1]])
    
    K_elem_bar *= (A*E/L)
   

    return K_elem_bar


def rotation_matrix_bar(cs,ss ):
    
    T = np.matrix([[ cs,  ss,    0,   0],   
                    [0,   0,   cs,  ss ]])
    
    return T


'''
 main code for the bar elements

'''

def FEBar(data):
    
    
    print '* * * * '
    print '\n**********************************\n'
    print 'computing static finite element solution using a bar element'
    print '\n**********************************\n'
    
    
    
    
    total_dof = data.ndof* data.n_nodes
    
    print 'loop in the elements'
    
    
    K_gl = np.zeros((total_dof,total_dof))
    
    for i in range (0,data.n_elem):
        # print i
        
        i_prop = data.elements[i].prop
        n1 = data.elements[i].nodes[0]
        n2 = data.elements[i].nodes[1]
        
        x1 = data.nodes[n1-1].coord
        x2 = data.nodes[n2-1].coord
        
        dx = x2 - x1
        
        Le = norm(dx)
        
        Area = data.sections[i_prop-1].A
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E = data.materials[i_mat-1].E
        
        
        dx_norm = dx/Le
        ctheta = dx_norm[0]
        stheta = dx_norm[1]
        
        K_elem = compute_K_bar(Area, E, Le)
        T = rotation_matrix_bar(ctheta, stheta)
        
        K_e_gl = np.dot(T.transpose(), np.dot(K_elem, T ))
        
        '''
        print dx, Le, Area, E, ctheta, stheta
        print T
        print K_e_gl
        '''
        
        K11 = K_e_gl[0:2,0:2]
        K12 = K_e_gl[0:2,2:4]
        K21 = K_e_gl[2:4,0:2]
        K22 = K_e_gl[2:4,2:4]
        
        '''
        print 'K11'
        print K11
        print K12
        print K21
        print K22
        '''
       
        i1_i   = (n1-1)*data.ndof ;  
        i1_f   = (n1-1)*data.ndof + 1;  
        i2_i   = (n2-1)*data.ndof ; 
        i2_f   = (n2-1)*data.ndof + 1;
         
        K_gl[i1_i:i1_f+1 , i1_i:i1_f+1 ] += K11
        K_gl[i1_i:i1_f+1 , i2_i:i2_f+1 ] += K12
        K_gl[i2_i:i2_f+1 , i1_i:i1_f+1 ] += K21
        K_gl[i2_i:i2_f+1 , i2_i:i2_f+1 ] += K22
        
        #-- end of elements loop
        
    #print K_gl    
    
    # apply boundary cc
    print 'apply boundary condition'
    for bc in data.bconditions:
        
        dofbc =  (bc.id -1)*data.ndof + bc.dir -1;
        
        K_gl[0:total_dof, dofbc ] = 0.
        K_gl[dofbc,0:total_dof  ] = 0.
        K_gl[dofbc,dofbc  ] = 1
        
        
        
        
    # force vector
    print 'assemble the force vector'
    f_gl = np.zeros((total_dof,1))
    
    for force in data.forces:
        
        doff = (force.id -1)*data.ndof + force.dir -1;
        
        f_gl [doff,0]=force.A
        
        
    #print K_gl  
    #print f_gl  
    
    print 'solve the linear system, K.q = f'
    u_gl =  solve(K_gl, f_gl)
    
    # print u_gl
    
    #-------------------------
    # pos-processamento
    u_nodal  = np.zeros((data.n_nodes,2))   
    
    header = 'i'.rjust(4) + '|'
    header+= align_header('u_1',13)
    header+= align_header('u_2',13)
    print( header) 
          
    for i in range(0,data.n_nodes):
    
        gdl_1 = i*data.ndof + 0;  
        gdl_2 = i*data.ndof + 1;  
        
        u_nodal[i,0] = u_gl[gdl_1,0];
        u_nodal[i,1] = u_gl[gdl_2,0];
        
        
        outline = '%3d |' % i
        outline+= '%12.8f |' % u_nodal[i,0]
        outline+= '%12.8f |' % u_nodal[i,1]
    
        print (outline)
    
    
    
    
    
    Results = C_Results()
    
    Results.u_global = u_nodal
    
    Results.save_to_gmsh(data.path,data.filebase,data.nodes,data.elements)
        
    return


