#!/usr/bin/env python3
"""
@mainpage FEBeam


@section Goals


@section Todo

@todo


@section Teste


@package PyEolica

@author carlos eduardo de souza
@version 18.09.19  start

@attention: esse arquivo ainda esta com muita cara de c++! Como usar comandos em python?
@note
    Esse programa deve se manter sem interface grafica.
    Usar sempre o conceito: rodar sozinho!!!
    Posteriomente, devo criar uma interface que chame as rotinas atuais.

@note
Para documentacao, ver o site <a href="http://notemagnet.blogspot.com.br/2009/10/using-doxypy-for-python-code.html">Using doxypy for Python code documentation </a>

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


    

def compute_K_beam2d(A,E,I,L):
    
    
    K_B = E * A / L
    K_V = E * I / L**3
    
        
    K_beam2d = np.matrix([ [K_B  ,        0,          0,  -K_B,        0,           0 ], 
                           [0    ,   12*K_V,    6*L*K_V,     0,  -12*K_V,     6*L*K_V ],
                           [0    ,  6*L*K_V, 4*L**2*K_V,     0, -6*L*K_V,  2*L**2*K_V ],
                           [-K_B ,        0,          0,   K_B,        0,           0 ],
                           [0    ,  -12*K_V,   -6*L*K_V,     0,   12*K_V,    -6*L*K_V ],
                           [0    ,  6*L*K_V, 2*L**2*K_V,     0, -6*L*K_V,  4*L**2*K_V ]])
    
    
    
   

    return K_beam2d


def rotation_matrix_beam2d(cs,ss ):
    
    T = np.matrix([[ cs,  ss,    0,   0,   0,  0 ],  
                   [-ss,  cs,    0,   0,   0,  0 ],
                   [  0,   0,    1,   0,   0,  0 ],
                   [  0,   0,    0,  cs,  ss,  0 ],
                   [  0,   0,    0, -ss,  cs,  0 ],
                   [  0,   0,    0,   0,   0,  1 ]])
    
    return T


'''
 main code for the bar elements

'''

def FEBeam2d(data):
    
    
    print ('* * * * ')
    print ('\n**********************************\n')
    print ('computing static finite element     ') 
    print ('solution using a 2d beam element    ')
    print ('\n**********************************\n')
    
    
    
    
    total_dof = data.ndof* data.n_nodes
    
    print ('loop on  elements')
    
    
    K_gl = np.zeros((total_dof,total_dof))
    
    for i in range (0,data.n_elem):
        # print i
        
        n1 = data.elements[i].nodes[0]
        n2 = data.elements[i].nodes[1]
        
        x1 = data.nodes[n1-1].coord
        x2 = data.nodes[n2-1].coord
        
        dx = x2 - x1
        
        Le = norm(dx)
        
        # geometry
        i_prop = data.elements[i].prop
        Area = data.sections[i_prop-1].A
        Ixx = data.sections[i_prop-1].Ixx
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E = data.materials[i_mat-1].E
        
        
        dx_norm = dx/Le
        ctheta = dx_norm[0]
        stheta = dx_norm[1]
        
        
        
        K_elem = compute_K_beam2d(Area, E, Ixx, Le)
        T = rotation_matrix_beam2d(ctheta, stheta)
        
        
        data.elements[i].Le = Le
        data.elements[i].T = T
        
        K_e_gl = np.dot(T.transpose(), np.dot(K_elem, T ))
        
        
        print (dx, Le, Area, E, ctheta, stheta)
        print (T)
        print (K_e_gl)
        
        
        K11 = K_e_gl[0:3,0:3]
        K12 = K_e_gl[0:3,3:6]
        K21 = K_e_gl[3:6,0:3]
        K22 = K_e_gl[3:6,3:6]
        
        '''
        print 'K11'
        print K11
        print K12
        print K21
        print K22
        '''
       
        i1_i   = (n1-1)*data.ndof ;  
        i1_f   = (n1-1)*data.ndof + 2;  
        i2_i   = (n2-1)*data.ndof ; 
        i2_f   = (n2-1)*data.ndof + 2;
         
        K_gl[i1_i:i1_f+1 , i1_i:i1_f+1 ] += K11
        K_gl[i1_i:i1_f+1 , i2_i:i2_f+1 ] += K12
        K_gl[i2_i:i2_f+1 , i1_i:i1_f+1 ] += K21
        K_gl[i2_i:i2_f+1 , i2_i:i2_f+1 ] += K22
        
        #-- end of elements loop
        
    #print K_gl    
    
    #--- apply boundary cc
    print ('apply boundary condition')
    for bc in data.bconditions:
        
        dofbc =  (bc.id -1)*data.ndof + bc.dir -1;
        
        K_gl[0:total_dof, dofbc ] = 0.
        K_gl[dofbc,0:total_dof  ] = 0.
        K_gl[dofbc,dofbc  ] = 1
        
        
        
        
    #--- force vector
    print ('assemble the force vector')
    f_gl = np.zeros((total_dof,1))
    
    for force in data.forces:
        
        doff = (force.id -1)*data.ndof + force.dir -1;
        
        f_gl [doff,0]=force.A
        
        
    #print K_gl  
    #print f_gl  
    
    print ('solve the linear system, K.q = f')
    u_gl =  solve(K_gl, f_gl)
    
    # print u_gl
    
    #-------------------------
    #--- pos-processing
    u_nodal  = np.zeros((data.n_nodes,2))   
    
    
    print(' Displacement results:')
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
    print(' Strain and stress results (need to correct it!):')
    header = 'i'.rjust(4) + '|'
    header+= align_header('Le',13)
    header+= align_header('Lef',13)
    header+= align_header('epsilon',13)
    header+= align_header('sigma',13)
    print( header) 
    for e in data.elements:
         
        # geometry
        i_prop = e.prop
        Area = data.sections[i_prop-1].A
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E = data.materials[i_mat-1].E    
          
        n1 = e.nodes[0]
        n2 = e.nodes[1]
        
        x1 = data.nodes[n1-1].coord
        x2 = data.nodes[n2-1].coord
        
        x1f = x1 + [u_nodal[n1-1,0], u_nodal[n1-1,1], 0]
        x2f = x2 + [u_nodal[n2-1,0], u_nodal[n2-1,1], 0]
        dxf = x2f - x1f
        Lef = norm(dxf)
        
        epsilon = (Lef-Le)/Le
        sigma = E * epsilon
        
        Results.eps.append(epsilon)
        
        Results.S.append(sigma/1e6)
        
         
        outline = '%3d |' % e.id
        outline+= '%12.8f |' % Le
        outline+= '%12.8f |' % Lef
        outline+= '%12.8f |' % epsilon
        outline+= '%12.8f |' % (sigma/1.e6)
        
    
        print (outline)
    
    
    Results.u_global = u_nodal
    
    Results.save_to_gmsh(data.path,data.filebase,data.nodes,data.elements)
        
    return


