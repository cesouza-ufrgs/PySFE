#!/usr/bin/env python3
"""
@mainpage FEBar

This is the main page of FEBar

@section Goals

Bar finite elements - 

@section Todo

@todo
 - 18.09.18 - everything, planning



@package PyFE

@author carlos eduardo de souza
@date 19.04.10
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


    

def compute_K_bar(A,E,L ):
        
    K_elem_bar = np.matrix([[1,-1],[-1,1]])
    
    K_elem_bar *= (A*E/L)

    return K_elem_bar


def rotation_matrix_bar(cs,ss ):
    
    T = np.matrix([[ cs,  ss,    0,   0],   
                    [0,   0,   cs,  ss ]])
    
    return T

def rotation_matrix_bar3d(vl):
    
    ''' full matrix:
    T = np.matrix([[ vl[0],  vm[0],  vn[0],   0,   0,   0],   
                   [ vl[1],  vm[1],  vn[1],   0,   0,   0],   
                   [ vl[2],  vm[2],  vn[2],   0,   0,   0],   
                   [      0,     0,      0, vl[0],  vm[0],  vn[0]],   
                   [      0,     0,      0, vl[1],  vm[1],  vn[1]],   
                   [      0,     0,      0, vl[2],  vm[2],  vn[2]]])
    '''
    T = np.matrix([[ vl[0],  vl[1],  vl[2],   0,   0,   0],  
                   [      0,     0,      0, vl[0],  vl[1],  vl[2]]])
   
    return T


def rotate_Kbar_3d(Ke,de):
    
    
    
    return 


'''
 computation of the 2d bar element problem

'''

def FEBar2d(data):
    
    
    print ('* * * * ')
    print ('\n**********************************\n')
    print ('computing static finite element solution using a bar element')
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
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E = data.materials[i_mat-1].E
        
        
        dx_norm = dx/Le
        ctheta = dx_norm[0]
        stheta = dx_norm[1]
        
        
        
        K_elem = compute_K_bar(Area, E, Le)
        T = rotation_matrix_bar(ctheta, stheta)
        
        
        data.elements[i].Le = Le
        data.elements[i].T = T
        
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
    print(' Strain and stress results:')
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
        
        Results.epsilon.append(epsilon)
        
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


def FEBar3d(data):
    
    
    print ('* * * * ')
    print ('\n**********************************\n')
    print ('computing static finite element solution using a bar element')
    print ('\n**********************************\n')
    
    
    
    
    total_dof = data.ndof* data.n_nodes
    
    print ('loop on  elements')
    
    
    K_gl = np.zeros((total_dof,total_dof))
    
    for e in data.elements:
        # print i
        
        n1 = e.nodes[0]
        n2 = e.nodes[1]
        
        x1 = data.nodes[n1-1].coord
        x2 = data.nodes[n2-1].coord
        
        d = x2 - x1
        
        Le = norm(d)
        
        # norm
        de = d/Le
        
        # geometry
        i_prop = e.prop
        Area = data.sections[i_prop-1].A
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E     = data.materials[i_mat-1].E
        
                
                
        K_elem = compute_K_bar(Area, E, Le)
        
        T = rotation_matrix_bar3d(de)
                
        e.Le = Le
        e.T = T
        
        K_e_gl = np.dot(T.transpose(), np.dot(K_elem, T ))
        
        '''
        print dx, Le, Area, E, ctheta, stheta
        print T
        print K_e_gl
        '''
        
        K11 = K_e_gl[0:data.ndof          ,0:data.ndof]
        K12 = K_e_gl[0:data.ndof          ,data.ndof:2*data.ndof]
        K21 = K_e_gl[data.ndof:2*data.ndof,0:data.ndof]
        K22 = K_e_gl[data.ndof:2*data.ndof,data.ndof:2*data.ndof]
        
        '''
        print ('K11')
        print (K11)
        print (K12)
        '''
       
        i1_i   = (n1-1)*data.ndof ;  
        i1_f   = (n1  )*data.ndof ;  
        i2_i   = (n2-1)*data.ndof ; 
        i2_f   = (n2  )*data.ndof ;
         
        K_gl[i1_i:i1_f , i1_i:i1_f ] += K11
        K_gl[i1_i:i1_f , i2_i:i2_f ] += K12
        K_gl[i2_i:i2_f , i1_i:i1_f ] += K21
        K_gl[i2_i:i2_f , i2_i:i2_f ] += K22
        
        #-- end of elements loop
        
    #print K_gl    
    
    #--- apply boundary cc
    print ('apply boundary condition')
    for bc in data.bconditions:
        
        dofbc =  (bc.id -1)*data.ndof + bc.dir -1;
        
        K_gl[0:total_dof, dofbc ] = 0.
        K_gl[dofbc      ,0:total_dof  ] = 0.
        K_gl[dofbc      ,dofbc  ] = 1
        
        
        
        
    #--- force vector
    print ('assemble the force vector')
    f_gl = np.zeros((total_dof,1))
    
    for force in data.forces:
        
        doff = (force.id - 1 )*data.ndof + force.dir -1;
        
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
    print(' Strain and stress results:')
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
        
        Results.epsilon.append(epsilon)
        
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
