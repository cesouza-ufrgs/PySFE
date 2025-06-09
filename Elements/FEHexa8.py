"""
@mainpage FEHexa8


@section Goals


@section Todo

@todo


@section Teste


@package PyFE

@author carlos eduardo de souza
@version 20.10.05

@attention: 
@note
   

"""

import sys  # os, codecs , fileinput
# import string
# import math
import numpy as np
from tabulate import tabulate

from numpy.linalg import norm, solve, det,inv
from math import radians, sqrt #, pi, cos, sin, sqrt,


# classes locais -------------
from Data.class_data   import *
from Util.utilities import align_header
from math import atan2

from Data.class_results import *


class class_Phi: 
   
    def __init__(self,nfunc, ndim):  
        self.phi  = np.zeros((nfunc))
        self.dphi = np.zeros((nfunc,ndim))
        self.dphinc  = np.zeros((ndim,ndim))

#-----------------------
#--- global info ------
# defined here to avoid resizing

n_nodes_h8 = 8

n_int_h8 = 8
           
G = (1./sqrt(3.));      
int_points = np.matrix([  [-G, -G,  -G],   
                          [-G,  G,  -G],
                          [ G,  G,  -G],
                          [ G, -G,  -G],
                          [-G, -G,   G],            
                          [-G,  G,   G],
                          [ G,  G,   G],
                          [ G, -G,   G]])

W_G8 =  [1., 1., 1., 1.,1., 1., 1., 1.]  

# each row: integration point
# each col: value of function (1 function per node.)

# loop on gauss - precompute all functions
PHI  = []
for i_g in range(0,8):
   
    print('i_g',i_g)
    xi1 = int_points[i_g,0];
    xi2 = int_points[i_g,1];
    xi3 = int_points[i_g,2];
   
    Phi = class_Phi(8,3) 
  
    # each col is the node position
   
    Phi.phi[0] = (1-xi1)*(1-xi2)*(1+xi3)/8         
    Phi.phi[1] = (1-xi1)*(1-xi2)*(1-xi3)/8
    Phi.phi[2] = (1-xi1)*(1+xi2)*(1-xi3)/8
    Phi.phi[3] = (1-xi1)*(1+xi2)*(1+xi3)/8
    Phi.phi[4] = (1+xi1)*(1-xi2)*(1+xi3)/8
    Phi.phi[5] = (1+xi1)*(1-xi2)*(1-xi3)/8
    Phi.phi[6] = (1+xi1)*(1+xi2)*(1-xi3)/8
    Phi.phi[7] = (1+xi1)*(1+xi2)*(1+xi3)/8 
    
    
    Phi.dphi[0,0] =  -(1-xi2)*(1+xi3)/8 
    Phi.dphi[1,0] =  -(1-xi2)*(1-xi3)/8 
    Phi.dphi[2,0] =  -(1+xi2)*(1-xi3)/8 
    Phi.dphi[3,0] =  -(1+xi2)*(1+xi3)/8 
    Phi.dphi[4,0] =   (1-xi2)*(1+xi3)/8  
    Phi.dphi[5,0] =   (1-xi2)*(1-xi3)/8  
    Phi.dphi[6,0] =   (1+xi2)*(1-xi3)/8  
    Phi.dphi[7,0] =   (1+xi2)*(1+xi3)/8              
    
    Phi.dphi[0,1] =   -(1-xi1)*(1+xi3)/8
    Phi.dphi[1,1] =   -(1-xi1)*(1-xi3)/8
    Phi.dphi[2,1] =    (1-xi1)*(1-xi3)/8 
    Phi.dphi[3,1] =    (1-xi1)*(1+xi3)/8 
    Phi.dphi[4,1] =   -(1+xi1)*(1+xi3)/8
    Phi.dphi[5,1] =   -(1+xi1)*(1-xi3)/8
    Phi.dphi[6,1] =    (1+xi1)*(1-xi3)/8 
    Phi.dphi[7,1] =    (1+xi1)*(1+xi3)/8 
   
    Phi.dphi[0,2] =   (1-xi1)*(1-xi2)/8
    Phi.dphi[1,2] =  -(1-xi1)*(1-xi2)/8
    Phi.dphi[2,2] =  -(1-xi1)*(1+xi2)/8
    Phi.dphi[3,2] =   (1-xi1)*(1+xi2)/8
    Phi.dphi[4,2] =   (1+xi1)*(1-xi2)/8
    Phi.dphi[5,2] =  -(1+xi1)*(1-xi2)/8
    Phi.dphi[6,2] =  -(1+xi1)*(1+xi2)/8
    Phi.dphi[7,2] =   (1+xi1)*(1+xi2)/8  
    
    
    Phi.dphinc[0,0] = xi1
    Phi.dphinc[1,1] = xi2
    Phi.dphinc[2,2] = xi3
    
    PHI.append(Phi)
    
    
#--- non conforming - bubble
   
wC =  [2.]  
PhiC = class_Phi(8,3) 
   
PhiC.dphi[0,0]=  -1/8 
PhiC.dphi[1,0]=  -1/8 
PhiC.dphi[2,0]=  -1/8 
PhiC.dphi[3,0]=  -1/8 
PhiC.dphi[4,0]=   1/8  
PhiC.dphi[5,0]=   1/8  
PhiC.dphi[6,0]=   1/8  
PhiC.dphi[7,0]=   1/8           

PhiC.dphi[0,1]=   -1/8
PhiC.dphi[1,1]=   -1/8
PhiC.dphi[2,1]=    1/8 
PhiC.dphi[3,1]=    1/8 
PhiC.dphi[4,1]=   -1/8
PhiC.dphi[5,1]=   -1/8
PhiC.dphi[6,1]=    1/8 
PhiC.dphi[7,1]=    1/8 

PhiC.dphi[0,2] =   1/8
PhiC.dphi[1,2] =  -1/8
PhiC.dphi[2,2] =  -1/8
PhiC.dphi[3,2] =   1/8
PhiC.dphi[4,2] =   1/8
PhiC.dphi[5,2] =  -1/8
PhiC.dphi[6,2] =  -1/8
PhiC.dphi[7,2] =   1/8 
PHIC  = []
PHIC.append(PhiC)   



#--- 
# matrices sizes

Ke33 =np.zeros((11*3,11*3)) 
Ke = np.zeros((8*3,8*3)); 

    
def compute_B_hexa8(dphi_dX,dphinc_dX):
   
    B = np.zeros((2*3,11*3))
    
      
    for I in range(0,8) : 
        col = 3*I
        B[0,col  ] =  dphi_dX[I,0]  
        B[1,col+1] =  dphi_dX[I,1]   
        B[2,col+2] =  dphi_dX[I,2] 
        
        B[3,col+1] =  dphi_dX[I,2]  
        B[3,col+2] =  dphi_dX[I,1] 
                
        B[4,col+0] =  dphi_dX[I,2]  
        B[4,col+2] =  dphi_dX[I,0] 
        
        B[5,col+0] =  dphi_dX[I,1]  
        B[5,col+1] =  dphi_dX[I,0]  
        
    for I in range(0,3) :
        col = 24+3*I
        
        B[0,col  ] =  dphinc_dX[I,0]  
        B[1,col+1] =  dphinc_dX[I,1]   
        B[2,col+2] =  dphinc_dX[I,2] 
        
        B[3,col+1] =  dphinc_dX[I,2]  
        B[3,col+2] =  dphinc_dX[I,1] 
                
        B[4,col+0] =  dphinc_dX[I,2]  
        B[4,col+2] =  dphinc_dX[I,0] 
        
        B[5,col+0] =  dphinc_dX[I,1]  
        B[5,col+1] =  dphinc_dX[I,0] 
    
    
    return B


'''
@brief Computation of

'''
def compute_K_hexa8(ni,nodes,E,nu,):
    
    #print ('compute K hexa8')
 
    #print(ni)
 
    X = np.zeros((3,8))
     
    for i_n in range(0,8):
       #print('i_n',i_n,ni[i_n-1])
       X[0,i_n] =  nodes[ni[i_n]-1].coord[0]
       X[1,i_n] =  nodes[ni[i_n]-1].coord[1]
       X[2,i_n] =  nodes[ni[i_n]-1].coord[2]
        
    #print(X)
   
    C = np.matrix([ [ 1.,   nu,  nu,  0.,   0.,   0. ],
                    [ nu,   1.,  nu,  0.,   0.,   0. ],
                    [ nu,   nu,  1.,  0.,   0.,   0. ],
                    [ 0.,   0.,   0., (1.-nu)/2.,         0.,   0. ],
                    [ 0.,   0.,   0.,         0., (1.-nu)/2.,   0. ],
                    [ 0.,   0.,   0.,         0.,         0., (1.-nu)/2.]])
   
    C *= E /(1.-nu*nu) 
              
    J0 = np.dot(X, PHIC[0].dphi)
    invJ0 = inv(J0)
    detJ0 = det(J0)
    
    #print(PHIC[0].dphi)
    
           
    #J = np.zeros((3,3))
    global Ke33
    
    Ke33.fill(0.)
    volume_e = 0.;
    # print(tabulate (Ke33, floatfmt=".3f", tablefmt="grid"))
    
    for i_g in range(0,n_int_h8):
        
        
        # (3 x 8 ) (8 , 3)
        J = np.dot(X, PHI[i_g].dphi)
        detJ = det(J);
        invJ = inv(J);  
        #print(PHI[i_g].dphi)
        
        #print('Jacobian')
        #print(J)
        #print(invJ)
        #print(detJ)   
   
      
        dphi_dX = np.dot(PHI[i_g].dphi, invJ);    
        #print(dphi_dX)      
   
        # non conforming dPhinc_dX
        
        dphinc_dX = np.dot(PHI[i_g].phinc, invJ0); 
   
        B = compute_B_hexa8(dphi_dX,dphinc_dX)
        
        #print('B\n',B.transpose())
        #print('bt\n',BT)
        Ke33 += W_G8[i_g] * detJ * np.dot( B.transpose(),np.dot(C,B) ) 
        volume_e += detJ
        
    #print(tabulate (Ke33, floatfmt=".3f", tablefmt="grid"))
    # print('K_e computed', volume_e)
    
    ## Condensacao estatica de 33x33 para 24x24
    KA = Ke33[ 0:24, 0:24]
    KB = Ke33[ 0:24,24:33]
    KC = Ke33[24:33,24:33]
    
    global Ke 
    Ke = KA -  np.dot( KB , np.dot( inv(KC),  KB.transpose() )  )
    
'''
@brief Computation of

'''
def compute_stress_hexa8(ni,nodes,E,nu,u_global,Results):
   
   
    u_e = np.zeros((24,1))
   
    for i_n in range(0,8):
       #print('i_n',i_n,ni[i_n-1])
       for i_g in range(0,3):
           posgl = (ni[i_n]-1)*3+i_g
           u_e[i_n*3+i_g] =  u_global[posgl]
       
        
    X = np.zeros((3,8))
          
    
    for i_n in range(0,8):
       #print('i_n',i_n,ni[i_n-1])
       X[0,i_n] =  nodes[ni[i_n]-1].coord[0]
       X[1,i_n] =  nodes[ni[i_n]-1].coord[1]
       X[2,i_n] =  nodes[ni[i_n]-1].coord[2]
        
    #print(X)
   
    C = np.matrix([ [ 1.,   nu,  nu,  0.,   0.,   0. ],
                    [ nu,   1.,  nu,  0.,   0.,   0. ],
                    [ nu,   nu,  1.,  0.,   0.,   0. ],
                    [ 0.,   0.,   0., (1.-nu)/2.,         0.,   0. ],
                    [ 0.,   0.,   0.,         0., (1.-nu)/2.,   0. ],
                    [ 0.,   0.,   0.,         0.,         0., (1.-nu)/2.]])
   
    C *= E /(1.-nu*nu) 
    
    dphinc_dX = np.zeros((6,9)); 
    #print(PHIC[0].dphi)
           
    
    # print(tabulate (Ke33, floatfmt=".3f", tablefmt="grid"))
    epsilon = np.zeros((6,1))
    sigma   = np.zeros((6,1))
    
    for i_g in range(0,n_int_h8):
                
        # (3 x 8 ) (8 , 3)
        J = np.dot(X, PHI[i_g].dphi)
        detJ = det(J);
        invJ = inv(J);  
        #print(PHI[i_g].dphi)
        
        #print('Jacobian')
        #print(J)
        #print(invJ)
        #print(detJ)   
   
      
        dphi_dX = np.dot(PHI[i_g].dphi, invJ);    
        #print(dphi_dX)      
   
        # non conforming dPhinc_dX
        
   
        B = compute_B_hexa8(dphi_dX,dphinc_dX)
        
        #print('B\n',B.transpose())
        #print('bt\n',BT)
        
        epsilon_g = np.dot(B[0:7,0:24],u_e)
        
        #print(epsilon_g)
        
        sigma_g = np.dot(C , epsilon_g)
        
        
        epsilon += epsilon_g/8
        sigma += sigma_g/8
        
        
   
        
    Results.eps.append(epsilon)
        
    Results.S.append(sigma)     
   

'''

 main code for the plane3 element
 
 It actually computes a mesh with only plane3 elements, here. 

'''

def FEHexa8(data):
    
    
    print ('* * * * ')
    print ('\n*************************************\n')
    print ('computing static finite element        ') 
    print ('solution using a 3d - 8n hexa element ')
    print ('\n*************************************\n')
    
    total_dof = data.ndof * data.n_nodes
    
    print ('loop on  elements', total_dof)
    
    K_gl = np.zeros((total_dof,total_dof))
    
    
    for ie in range (0,data.n_elem):
        
        # print('element', ie)
                
        # geometry
        i_mat = data.elements[ie].prop
        
        # material
        E  = data.materials[i_mat-1].E
        nu = data.materials[i_mat-1].nu
        
        global Ke
        
        compute_K_hexa8(data.elements[ie].ni, data.nodes, E, nu)
        
          
        #print(tabulate(Ke[0:13,0:13], tablefmt="plain") )
        
        #---- ------------------------------------     
        # global stiffness matrix superposition
        # hee, the symmetry could be used to simplify the proccess 
        for I in range(0,8):
            
            ie1=data.ndof*I
            ie3=data.ndof*I+3
            
            inode=data.elements[ie].ni[I]-1
            
            # o numero de graus de liberdade.
            i_1  = data.ndof*inode  ;  
            i_3  = data.ndof*inode + 3;  
                        
            #print('  i ',data.elements[ie].ni[I],ie1,ie3, i_1,i_3)
            
            for J in range(0,8):
               
               je1=data.ndof*J
               je3=data.ndof*J+3
            
               jnode=data.elements[ie].ni[J]-1
           
               j_1  = data.ndof*jnode ;  
               j_3  = data.ndof*jnode + 3;  
              
               #print('j ',data.elements[ie].ni[J],je1,je3, j_1,j_3)
               # 1 1 
               K_gl[i_1:i_3 , j_1:j_3 ] +=  Ke[ie1:ie3,je1:je3];
             
                  
        #-- end of elements loop
        
    
    #print(tabulate(K_gl[0:13,0:13]))
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
    
    Results = C_Results()
    
    u_nodal  = np.zeros((data.n_nodes,3))   
   
    for i in range(0,data.n_nodes):
            
        u_nodal[i,0:3] = u_gl[i*data.ndof :i*data.ndof+3,0];
       
    
    
    for ie in range (0,data.n_elem):
        
        # print('element', ie)
                
        # geometry
        i_mat = data.elements[ie].prop
        
        # material
        E  = data.materials[i_mat-1].E
        nu = data.materials[i_mat-1].nu
        
        
        compute_stress_hexa8(data.elements[ie].ni,data.nodes,E,nu,u_gl,Results)
        
        
    #epsilon = np.zeros((data.n_elem,3))
    #sigma   = np.zeros((data.n_elem,3))
    #s_vm    = np.zeros((data.n_elem,1))
    
    
    Results.u_global = u_nodal
    
    
    Results.save_to_gmsh(data.path,data.filebase,data.nodes,data.elements,data.elemtype)
        
    return


