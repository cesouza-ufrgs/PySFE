"""
@mainpage FEPlane3


@section Goals


@section Todo

@todo


@section Teste


@package PyFE

@author carlos eduardo de souza
@version 19.05.08

@attention: 
@note
   

"""

import sys  # os, codecs , fileinput
# import string
# import math
import numpy as np

from numpy.linalg import norm, solve, det,inv
from math import radians #, pi, cos, sin, sqrt,


# classes locais -------------
from Data.class_data   import *
from Util.utilities import align_header
from math import atan2

from Data.class_results import *



def compute_Jac_plane3(x1,x2,x3):

    print('computing Jacobian')
    print(x1,x2,x3)
    dphi1_dxi1 = -1.;
    dphi2_dxi1 =  1.; 
    dphi3_dxi1 =  0.;
    dphi1_dxi2 = -1.;
    dphi2_dxi2 =  0.;
    dphi3_dxi2 =  1.;
    
    Jac=np.zeros((2,2));
       
    Jac[0,0] = dphi1_dxi1 * x1[0] + dphi2_dxi1*x2[0] + dphi3_dxi1*x3[0]
    Jac[0,1] = dphi1_dxi1 * x1[1] + dphi2_dxi1*x2[1] + dphi3_dxi1*x3[1]
      
    Jac[1,0] = dphi1_dxi2 * x1[0] + dphi2_dxi2*x2[0] + dphi3_dxi2*x3[0]
    Jac[1,1] = dphi1_dxi2 * x1[1] + dphi2_dxi2*x2[1] + dphi3_dxi2*x3[1]          
    
    return Jac
    
def compute_B_plane3(invJ):

    dphi1_dxi1 = -1.;
    dphi2_dxi1 =  1.; 
    dphi3_dxi1 =  0.;
    dphi1_dxi2 = -1.;
    dphi2_dxi2 =  0.;
    dphi3_dxi2 =  1.;
   
    B = np.zeros((3,6))
    
    B[0,0] = dphi1_dxi1 * invJ[0,0] + dphi1_dxi2 * invJ[0,1];
    B[0,2] = dphi2_dxi1 * invJ[0,0] + dphi2_dxi2 * invJ[0,1];
    B[0,4] = dphi3_dxi1 * invJ[0,0] + dphi3_dxi2 * invJ[0,1];      
    
    B[1,1] = dphi1_dxi1 * invJ[1,0] + dphi1_dxi2 * invJ[1,1];
    B[1,3] = dphi2_dxi1 * invJ[1,0] + dphi2_dxi2 * invJ[1,1];
    B[1,5] = dphi3_dxi1 * invJ[1,0] + dphi3_dxi2 * invJ[1,1];
    
    B[2,0] = B[1,1]
    B[2,1] = B[0,0]
    B[2,2] = B[1,3]
    B[2,3] = B[0,2]
    B[2,4] = B[1,5]
    B[2,5] = B[0,4]   
    
    return B


def compute_K_plane3(x1,x2,x3,h,E,nu):
    
    dphi1_dxi1 = -1.;
    dphi2_dxi1 =  1.; 
    dphi3_dxi1 =  0.;
    dphi1_dxi2 = -1.;
    dphi2_dxi2 =  0.;
    dphi3_dxi2 =  1.;
    
    n_gauss = 3
    w =  [2/6, 1/6, 1/6]                       

    K_e = np.zeros((6,6));
   
    C = np.matrix([ [ 1.,   nu,   0.],
                    [ nu,   1.,   0. ],
                    [ 0.,   0.,   (1.-nu)/2.]])
   
    C *= E /(1.-nu*nu) 
              
    Jac = compute_Jac_plane3(x1, x2, x3)
    invJ = inv(Jac)
    detJ = det(Jac)
    
    print('Jacobian')
    print(Jac)
    print(invJ)
    
    B = compute_B_plane3(invJ)
    
    CB = np.dot(C,  B)
    BT = B.transpose()
   
    for i_g in range(0,n_gauss):
   
        #print('cb\n',CB)
        #print('bt\n',BT)
        K_e += h * w[i_g] * detJ * np.dot(BT,CB) 
        
    print('K_e computed')
    
    return K_e
    

'''

 main code for the plane3 element
 
 It actually computes a mesh with only plane3 elements, here. 

'''

def FEPlane3(data):
    
    
    print ('* * * * ')
    print ('\n*************************************\n')
    print ('computing static finite element        ') 
    print ('solution using a 2d - 3n plane element ')
    print ('\n*************************************\n')
    
    total_dof = data.ndof * data.n_nodes
    
    print ('loop on  elements')
    
    K_gl = np.zeros((total_dof,total_dof))
    
    for i in range (0,data.n_elem):
        print('element', i)
        
        n1 = data.elements[i].nodes[0]
        n2 = data.elements[i].nodes[1]
        n3 = data.elements[i].nodes[2]
        
        print(n1,n2,n3)
        
        x1 = data.nodes[n1-1].coord
        x2 = data.nodes[n2-1].coord
        x3 = data.nodes[n3-1].coord
        
        # geometry
        i_prop = data.elements[i].prop
        #Area = data.sections[i_prop-1].A
        h = data.sections[i_prop-1].h
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E  = data.materials[i_mat-1].E
        nu = data.materials[i_mat-1].nu
        
        K_e_gl = compute_K_plane3(x1,x2,x3, h, E, nu)
        
        print (x1,x2,x3,h,E,nu)
        # print (K_e_gl)      
        
        #---- ------------------------------------     
        # global stiffness matrix superposition
        # hee, the symmetry could be used to simplify the proccess 
        
        K11 = K_e_gl[0:2,0:2]
        K12 = K_e_gl[0:2,2:4]
        K13 = K_e_gl[0:2,4:6]
        
        K22 = K_e_gl[2:4,2:4]
        K23 = K_e_gl[2:4,4:6]        
        
        K33 = K_e_gl[4:6,4:6]        
    
        # o numero de graus de liberdade.
        i1_i  = (n1-1)*data.ndof + 0;  
        i1_f  = (n1-1)*data.ndof + 1;  
        i2_i  = (n2-1)*data.ndof + 0; 
        i2_f  = (n2-1)*data.ndof + 1;
        i3_i  = (n3-1)*data.ndof + 0; 
        i3_f  = (n3-1)*data.ndof + 1;
       
        # 1 1 
        K_gl[i1_i:i1_f+1,i1_i:i1_f+1] += K11;
        K_gl[i1_i:i1_f+1,i2_i:i2_f+1] += K12;
        K_gl[i1_i:i1_f+1,i3_i:i3_f+1] += K13;

        # 2 1
        K_gl[i2_i:i2_f+1,i1_i:i1_f+1] += K12.transpose();
        K_gl[i2_i:i2_f+1,i2_i:i2_f+1] += K22;
        K_gl[i2_i:i2_f+1,i3_i:i3_f+1] += K23;
   
        # 2 1
        K_gl[i3_i:i3_f+1,i1_i:i1_f+1] += K13.transpose();
        K_gl[i3_i:i3_f+1,i2_i:i2_f+1] += K23.transpose();
        K_gl[i3_i:i3_f+1,i3_i:i3_f+1] += K33;
               
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
    
    
    
    #epsilon = np.zeros((data.n_elem,3))
    #sigma   = np.zeros((data.n_elem,3))
    #s_vm    = np.zeros((data.n_elem,1))
    
    Results = C_Results()
    print(' Strain and stress results (need to correct it!):')
    header = 'i'.rjust(4) + '|'
    header+= align_header('Le',13)
    header+= align_header('Lef',13)
    header+= align_header('epsilon',13)
    header+= align_header('sigma',13)
    print( header) 
    for i in range (0,data.n_elem):
        print('element', i)
        
        n1 = data.elements[i].nodes[0]
        n2 = data.elements[i].nodes[1]
        n3 = data.elements[i].nodes[2]
        
        print(n1,n2,n3)
        
        x1 = data.nodes[n1-1].coord
        x2 = data.nodes[n2-1].coord
        x3 = data.nodes[n3-1].coord
        
        
        print(n1,n2,n3)
        print(x1,x2,x3)
        # geometry
        i_prop = data.elements[i].prop
        h = data.sections[i_prop-1].h
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E  = data.materials[i_mat-1].E    
        nu = data.materials[i_mat-1].nu
        
        
        C = np.matrix([ [ 1.,   nu,   0.],
                        [ nu,   1.,   0. ],
                        [ 0.,   0.,   (1.-nu)/2.]])
         
        C *= E /(1.-nu*nu) 
        
        Jac = compute_Jac_plane3(x1, x2, x3)
        
        print('Jacobian')
        print(Jac)
        invJ = inv(Jac)
        detJ = det(Jac)
      
        print('Jacobian')
        print(invJ)
      
        B = compute_B_plane3(invJ)
        
        # u element
        i1_i  = (n1-1)*data.ndof + 0;  
        i1_f  = (n1-1)*data.ndof + 1;  
        i2_i  = (n2-1)*data.ndof + 0; 
        i2_f  = (n2-1)*data.ndof + 1;
        i3_i  = (n3-1)*data.ndof + 0; 
        i3_f  = (n3-1)*data.ndof + 1;
       
        u_e = np.zeros((3*data.ndof,1))
        
        u_e [0:2,0] = u_gl[i1_i:i1_f+1,0];
        u_e [2:4,0] = u_gl[i2_i:i2_f+1,0];
        u_e [4:6,0] = u_gl[i3_i:i3_f+1,0];
        
        epsilon = np.dot(B,u_e)
        
        sigma = np.dot(C , epsilon)
        
        Results.eps.append(epsilon)
        
        Results.S.append(sigma/1e6)
        
         
        outline = '%3d |' % data.elements[i].id
        for ep in epsilon:
            outline+= '%12.8f |' % ep
        for si in sigma:
            outline+= '%12.8f |' % (si/1.e6)
        
    
        print (outline)
    
    
    Results.u_global = u_nodal
    
    Results.save_to_gmsh(data.path,data.filebase,data.nodes,data.elements,data.elemtype)
        
    return


