"""
@mainpage FEPlane4


@section Goals


@section Todo

@todo


@section Teste


@package PyFE

@author carlos eduardo de souza
@version 20.05.03

@attention: 
@note
   

"""

import sys  # os, codecs , fileinput
# import string
# import math
import numpy as np

from numpy.linalg import norm, solve, det,inv
from math import sqrt #, pi, cos, sin, sqrt,


# classes locais -------------
from Data.class_data   import *
from Util.utilities import align_header, ruler
from math import atan2

from Data.class_results import *
 

"""
Define matrices globally, to avoid reallocation
"""
X = np.zeros((4,2))
dXdx = np.zeros((2,2))
phi  = np.zeros((4,1))
dphidxi = np.zeros((4,2))
u_e = np.zeros((8,1))
weight =  [1., 1., 1., 1.]     
B = np.zeros((3,8))
pg = sqrt(3)/3.
strain6 = np.zeros((6,1))
stress6 = np.zeros((6,1))

xi = np.array([[-pg, -pg],
               [ pg, -pg],
               [ pg,  pg],
               [-pg,  pg]])
   
n_gauss = 4
     
   
def compute_gauss(r,s):
    
    phi[0] = .25*(1.-r)*(1.-s)
    phi[1] = .25*(1.+r)*(1.-s)
    phi[2] = .25*(1.+r)*(1.+s)    
    phi[3] = .25*(1.-r)*(1.+s)

    dphidxi[0,0] = -.25*(1-s)
    dphidxi[1,0] =  .25*(1-s)
    dphidxi[2,0] =  .25*(1+s)
    dphidxi[3,0] = -.25*(1+s)
    
    dphidxi[0,1] = -.25*(1-r)
    dphidxi[1,1] = -.25*(1+r)
    dphidxi[2,1] =  .25*(1+r)
    dphidxi[3,1] =  .25*(1-r)

    
def update_B_plane4(dXdx):

    dxdX = np.linalg.inv(dXdx)    

    B[0,0] = dphidxi[0,0] * dxdX[0,0] + dphidxi[0,1] * dxdX[0,1]
    B[0,2] = dphidxi[1,0] * dxdX[0,0] + dphidxi[1,1] * dxdX[0,1]
    B[0,4] = dphidxi[2,0] * dxdX[0,0] + dphidxi[2,1] * dxdX[0,1]
    B[0,6] = dphidxi[3,0] * dxdX[0,0] + dphidxi[3,1] * dxdX[0,1]
    
    B[1,1] = dphidxi[0,0] * dxdX[1,0] + dphidxi[0,1] * dxdX[1,1]
    B[1,3] = dphidxi[1,0] * dxdX[1,0] + dphidxi[1,1] * dxdX[1,1]
    B[1,5] = dphidxi[2,0] * dxdX[1,0] + dphidxi[2,1] * dxdX[1,1]
    B[1,7] = dphidxi[3,0] * dxdX[1,0] + dphidxi[3,1] * dxdX[1,1]    
    
    B[2,0] = B[1,1]
    B[2,1] = B[0,0]
    B[2,2] = B[1,3]
    B[2,3] = B[0,2]
    B[2,4] = B[1,5]
    B[2,5] = B[0,4]   
    B[2,6] = B[1,7]   
    B[2,7] = B[0,6]     
    

def compute_K_plane4(X,e,E,nu):
    """
    Computes the stiffness matrix for the plane4 element

    Taylor element, bi-linear
    """
    # gaus points coordinates

    K_e = np.zeros((8,8))
   
    C = np.matrix([ [ 1.,   nu,   0.],
                    [ nu,   1.,   0. ],
                    [ 0.,   0.,   (1.-nu)/2.]])
   
    C *= E /(1.-nu*nu) 
                           
    
   
    for i_g in range(0,n_gauss):
    
        # print('gauss point')
        #  print(xi[i_g,:])
        # update phi and dphidxi
        compute_gauss(xi[i_g,0],xi[i_g,1]  )

        # print(dphidxi)

        # compute dXdx - deformation gradient
        dXdx[0,0] = np.dot(dphidxi[:,0], X[:,0])
        dXdx[0,1] = np.dot(dphidxi[:,1], X[:,0])
        dXdx[1,0] = np.dot(dphidxi[:,0], X[:,1])
        dXdx[1,1] = np.dot(dphidxi[:,1], X[:,1])

        # jacobian
        Jac = np.linalg.det(dXdx)
 
        update_B_plane4(dXdx)

        CB = np.dot(C,  B)
        BT = B.transpose()
        
        K_e += e * weight[i_g] * Jac * np.dot(BT,CB) 
        
    # print('K_e computed\n', K_e)
    
    return K_e
#---------------------------------------------------------------    
#---------------------------------------------------------------
'''

 main code for the plane3 element
 
 It actually computes a mesh with only plane3 elements, here. 

'''

def FEPlane4(data):
    
    ruler()
    print ('computing static finite element        ') 
    print ('solution using a 2d - 4n plane element ')
    ruler()
    
    total_dof = data.ndof * data.n_nodes
    
    # print ('loop on  elements')
    
    K_gl = np.zeros((total_dof,total_dof))

 
    for i in range (0,data.n_elem):
        # print('element', i)
        
        n1 = data.elements[i].nodes[0]
        n2 = data.elements[i].nodes[1]
        n3 = data.elements[i].nodes[2]
        n4 = data.elements[i].nodes[3]
          
        X[0,0:2] =data.nodes[n1-1].coord[0:2]
        X[1,0:2] =data.nodes[n2-1].coord[0:2]
        X[2,0:2] =data.nodes[n3-1].coord[0:2]
        X[3,0:2] =data.nodes[n4-1].coord[0:2]
        
        # geometry
        i_prop = data.elements[i].prop
        #Area = data.sections[i_prop-1].A
        e = data.sections[i_prop-1].e
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E  = data.materials[i_mat-1].E
        nu = data.materials[i_mat-1].nu
        
        #---- ------------------------------------  
        K_e_gl = compute_K_plane4(X, e, E, nu)        
        
        #---- ------------------------------------     
        # global stiffness matrix superposition
        # hee, the symmetry could be used to simplify the proccess 
        
        K11 = K_e_gl[0:2,0:2]
        K12 = K_e_gl[0:2,2:4]
        K13 = K_e_gl[0:2,4:6]
        K14 = K_e_gl[0:2,6:8]
        
        K22 = K_e_gl[2:4,2:4]
        K23 = K_e_gl[2:4,4:6]
        K24 = K_e_gl[2:4,6:8]       
        
        K33 = K_e_gl[4:6,4:6] 
        K34 = K_e_gl[4:6,6:8]  

        K44 = K_e_gl[6:8,6:8]          
    
        # o numero de graus de liberdade.
        i1_i  = (n1-1)*data.ndof + 0;  
        i1_f  = (n1-1)*data.ndof + 1;  
        i2_i  = (n2-1)*data.ndof + 0; 
        i2_f  = (n2-1)*data.ndof + 1;
        i3_i  = (n3-1)*data.ndof + 0; 
        i3_f  = (n3-1)*data.ndof + 1;
        i4_i  = (n4-1)*data.ndof + 0; 
        i4_f  = (n4-1)*data.ndof + 1;
       
        # 1 1 
        K_gl[i1_i:i1_f+1,i1_i:i1_f+1] += K11;
        K_gl[i1_i:i1_f+1,i2_i:i2_f+1] += K12;
        K_gl[i1_i:i1_f+1,i3_i:i3_f+1] += K13;
        K_gl[i1_i:i1_f+1,i4_i:i4_f+1] += K14;

        # 2 1
        K_gl[i2_i:i2_f+1,i1_i:i1_f+1] += K12.transpose();
        K_gl[i2_i:i2_f+1,i2_i:i2_f+1] += K22;
        K_gl[i2_i:i2_f+1,i3_i:i3_f+1] += K23;
        K_gl[i2_i:i2_f+1,i4_i:i4_f+1] += K24;
   
        # 2 1
        K_gl[i3_i:i3_f+1,i1_i:i1_f+1] += K13.transpose();
        K_gl[i3_i:i3_f+1,i2_i:i2_f+1] += K23.transpose();
        K_gl[i3_i:i3_f+1,i3_i:i3_f+1] += K33;
        K_gl[i3_i:i3_f+1,i4_i:i4_f+1] += K34;    

        # 2 2
        K_gl[i4_i:i4_f+1,i1_i:i1_f+1] += K14.transpose();
        K_gl[i4_i:i4_f+1,i2_i:i2_f+1] += K24.transpose();
        K_gl[i4_i:i4_f+1,i3_i:i3_f+1] += K34.transpose();
        K_gl[i4_i:i4_f+1,i4_i:i4_f+1] += K44;
               
        #-- end of elements loop
        
    #print K_gl    
    #--- apply boundary cc
    # print ('apply boundary condition')
    for bc in data.bconditions:
        
        dofbc =  (bc.id -1)*data.ndof + bc.dir -1;
        
        K_gl[0:total_dof, dofbc ] = 0.
        K_gl[dofbc,0:total_dof  ] = 0.
        K_gl[dofbc,dofbc  ] = 1
        
    #--- force vector
    # print ('assemble the force vector')
    f_gl = np.zeros((total_dof,1))
    
    for force in data.forces:
        
        doff = (force.id -1)*data.ndof + force.dir -1;
        
        f_gl [doff,0]=force.A
        
        
    #print K_gl  
    #print f_gl  
    
    # print ('solve the linear system, K.q = f')
    u_gl =  solve(K_gl, f_gl)
    
    # print(u_gl)

    # print('Solution obtained!')
    
    #-------------------------
    #--- pos-processing
    u_nodal  = np.zeros((data.n_nodes,2))   
    strain6 = np.zeros((data.n_elem,6))
    stress6 = np.zeros((data.n_elem,6))
    
    
    # print(' Displacement results:')
    header = 'i'.rjust(4) + '|'
    header+= align_header('u_1',13)
    header+= align_header('u_2',13)
 
          
    for i in range(0,data.n_nodes):
    
        gdl_1 = i*data.ndof + 0;  
        gdl_2 = i*data.ndof + 1;  
        
        u_nodal[i,0] = u_gl[gdl_1,0];
        u_nodal[i,1] = u_gl[gdl_2,0];        
        
        outline = '%3d |' % i
        outline+= '%12.8f |' % u_nodal[i,0]
        outline+= '%12.8f |' % u_nodal[i,1]
    
        # print (outline)
    
    
    
    #epsilon = np.zeros((data.n_elem,3))
    #sigma   = np.zeros((data.n_elem,3))
    #s_vm    = np.zeros((data.n_elem,1))
    
    print(' Strain and stress results (need to correct it!):')


    Results = C_Results()
    Results.strain = []
    Results.stress = []
    header = 'i'.rjust(4) + '|'
    header+= align_header('Le',13)
    header+= align_header('Lef',13)
    header+= align_header('epsilon',13)
    header+= align_header('sigma',13)
    # print( header) 
    for i in range (0,data.n_elem):
        # print('element', i)
         
        n1 = data.elements[i].nodes[0]
        n2 = data.elements[i].nodes[1]
        n3 = data.elements[i].nodes[2]
        n4 = data.elements[i].nodes[3]
          
        X[0,0:2] =data.nodes[n1-1].coord[0:2]
        X[1,0:2] =data.nodes[n2-1].coord[0:2]
        X[2,0:2] =data.nodes[n3-1].coord[0:2]
        X[3,0:2] =data.nodes[n4-1].coord[0:2]
        
        # geometry
        i_prop = data.elements[i].prop
        #Area = data.sections[i_prop-1].A
        e = data.sections[i_prop-1].e
        
        # material
        i_mat = data.sections[i_prop-1].imat
        E  = data.materials[i_mat-1].E
        nu = data.materials[i_mat-1].nu  
   
        C = np.matrix([ [ 1.,   nu,   0.],
                        [ nu,   1.,   0. ],
                        [ 0.,   0.,   (1.-nu)/2.]])
    
        C *= E /(1.-nu*nu) 
                        
         
        compute_gauss(0.,0.)


        # compute dXdx - deformation gradient
        dXdx[0,0] = np.dot(dphidxi[:,0], X[:,0])
        dXdx[0,1] = np.dot(dphidxi[:,1], X[:,0])
        dXdx[1,0] = np.dot(dphidxi[:,0], X[:,1])
        dXdx[1,1] = np.dot(dphidxi[:,1], X[:,1])

        update_B_plane4(dXdx)
        
        # u element
        i1_i  = (n1-1)*data.ndof + 0;  
        i1_f  = (n1-1)*data.ndof + 1;  
        i2_i  = (n2-1)*data.ndof + 0; 
        i2_f  = (n2-1)*data.ndof + 1;
        i3_i  = (n3-1)*data.ndof + 0; 
        i3_f  = (n3-1)*data.ndof + 1;
        i4_i  = (n4-1)*data.ndof + 0; 
        i4_f  = (n4-1)*data.ndof + 1;
       
        
        u_e [0:2,0] = u_gl[i1_i:i1_f+1,0];
        u_e [2:4,0] = u_gl[i2_i:i2_f+1,0];
        u_e [4:6,0] = u_gl[i3_i:i3_f+1,0];
        u_e [6:8,0] = u_gl[i4_i:i4_f+1,0];
        

        strain3 = np.dot(B,u_e)

        stress3 = np.dot(C , strain3)

        strain33 = -nu/E*(stress3[0]+stress3[1])


        strain6[i,0] = strain3[0]   # xx
        strain6[i,1] = strain3[1]   # yy
        strain6[i,2] = strain33     # zz 
        strain6[i,3] = strain3[2]   # xy
        strain6[i,4] = 0.           # yz
        strain6[i,5] = 0.           # xz  

        stress6[i,0] = stress3[0]   # xx
        stress6[i,1] = stress3[1]   # yy
        stress6[i,2] = 0.0          # zz 
        stress6[i,3] = stress3[2]   # xy
        stress6[i,4] = 0.           # yz
        stress6[i,5] = 0.           # xz  


        # print(epsilon)
        # print(sigma/1e6)
        
        Results.strain = strain6
        
        Results.stress = stress6
        
         
    
    
    print(Results.stress) 
    
    Results.u_global = u_nodal
    
    Results.save_to_gmsh(data.path,data.filebase,data.nodes,data.elements,data.elemtype)
        
    return


