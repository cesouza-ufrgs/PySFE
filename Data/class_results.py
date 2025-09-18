 

 
import sys, os, codecs , fileinput

from Data.class_data import *

# casse Resultados - parametros calculados e que devem ser salvos. 
class C_Results():
    '''
     
    @class C_Resultados
    
    @brief Classe para organizar resultados da analise  e salvar em arquivo.
    
    
    '''
    
    def __init__(self):
        
        self.p_q = 0
        
        self.energia = 0
        
        self.u_gl = []
        
        self.stress = []
        
        self.strain = []

        self.mainstress = []

        self.maindir1 = []
        self.maindir2 = []
        self.maindir3 = []
        
        return
        
    def save_to_gmsh(self,path,filebase,nodes,elements,elemtype):  
        '''
        funcao para salvar parametros calculados em um arquivo.    
        '''
       
        n_n = len(nodes)
        n_e = len(elements)
        
        file_out_gmsh = path + '/' + filebase + '_gmsh.msh'
        
        
        print(file_out_gmsh)
        
        fout = open(file_out_gmsh, 'w')    
    
        fout.write('$MeshFormat\n')    
        fout.write('2.2 0 4\n')    
        fout.write('$EndMeshFormat\n')  
            
        if elemtype == 'bar2d':
            self.save_to_gmsh_bar2d(fout,n_n,nodes,n_e,elements)
        elif   elemtype == 'plane4':    
            self.save_to_gmsh_plane4(fout,n_n,nodes,n_e,elements)      
        
        # ----------------------------------------      
        fout.close()    
    
        return
    
    
    
    def save_to_gmsh_bar2d(self,fout,n_n,nodes,n_e,elements):  
    
    
        #--- write nodes
        fout.write('$Nodes\n')  
        
        #linha = str("%5d  \n" % n_n)        
        fout.write(str("%5d  \n" % n_n))
        
        for no in nodes:
            
            linha = str("%5d  " % no.id)
            linha += str('%2.8f '  % no.coord[0])         
            linha += str('%2.8f '  % no.coord[1])
            linha += str('%2.8f '  % no.coord[2])
        
            fout.write(linha+'\n')
       
        fout.write('$EndNodes\n')  
        
        #--- write elements
        
        fout.write('$Elements\n')  
            
        fout.write(str("%5d  \n" % n_e))
        
        for el in elements:
            
            linha = str("%5d  1 2 0 0 " % el.id)                
        
            linha += str('%5d '  % el.nodes[0])         
            linha += str('%5d '  % el.nodes[1])
        
            fout.write(linha+'\n')
       
        fout.write('$EndElements\n')  
        
        #--- write displacements
        '''
        $NodeData
          numStringTags(ASCII int)
          stringTag(string) ...
          numRealTags(ASCII int)
          realTag(ASCII double) ...
          numIntegerTags(ASCII int)
          integerTag(ASCII int) ...
          nodeTag(size_t) value(double) ...
          ...
        $EndNodeData
        '''
        fout.write('$NodeData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"u\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('3\n')
        fout.write(str("%5d  \n" % n_n)) # number of nodes with results
        
        
        
        for i in range(0, n_n):
            
            linha = str('%5d  ' % (i+1))
            linha += str('%12.8f '  % self.u_global[i,0])         
            linha += str('%12.8f '  % self.u_global[i,1])      
            linha += str('%12.8f '  % 0.0)
            fout.write(linha+'\n')
      
        fout.write('$EndNodeData\n')  


    def save_to_gmsh_plane4(self,fout,n_n,nodes,n_e,elements):  
    
    
        #--- write nodes
        fout.write('$Nodes\n')  
        
        #linha = str("%5d  \n" % n_n)        
        fout.write(str("%5d  \n" % n_n))
        
        for no in nodes:
            
            linha = str("%5d  " % no.id)
            linha += str('%2.8f '  % no.coord[0])         
            linha += str('%2.8f '  % no.coord[1])
            linha += str('%2.8f '  % no.coord[2])
        
            fout.write(linha+'\n')
       
        fout.write('$EndNodes\n')  
        
        #--- write elements
        
        fout.write('$Elements\n')  
            
        fout.write(str("%5d  \n" % n_e))
        
        for el in elements:
            
            linha = str("%5d  3 2 1 1" % el.id)                
        
            linha += str('%5d '  % el.nodes[0])         
            linha += str('%5d '  % el.nodes[1])       
            linha += str('%5d '  % el.nodes[2])       
            linha += str('%5d '  % el.nodes[3])
        
            fout.write(linha+'\n')
       
        fout.write('$EndElements\n')  
        
        #--- write displacements
        '''
        $NodeData
          numStringTags(ASCII int)
          stringTag(string) ...
          numRealTags(ASCII int)
          realTag(ASCII double) ...
          numIntegerTags(ASCII int)
          integerTag(ASCII int) ...
          nodeTag(size_t) value(double) ...
          ...
        $EndNodeData
        '''
        fout.write('$NodeData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"u\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('3\n')
        fout.write(str("%5d  \n" % n_n)) # number of nodes with results
        
        
        
        for i in range(0, n_n):
            
            linha = str('%5d  ' % (i+1))
            linha += str('%12.8f '  % self.u_global[i,0])         
            linha += str('%12.8f '  % self.u_global[i,1])      
            linha += str('%12.8f '  % 0.0)
            fout.write(linha+'\n')
      
        fout.write('$EndNodeData\n')  

        '''
        $ElementData
        numStringTags(ASCII int)
        stringTag(string) ...
        numRealTags(ASCII int)
        realTag(ASCII double) ...
        numIntegerTags(ASCII int)
        integerTag(ASCII int) ...
        elementTag(int) value(double) ...
        ...
        $EndElementData
        '''

        fout.write('$ElementData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"Strain\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('9\n')
        fout.write(str("%5d  \n" % n_e)) # number of elements with results
        

        for i in range(0, n_e):
            
            linha = str('%5d  ' % (elements[i].id))
            linha += str('%15.9f '  % self.strain[i][0])   #xx    
            linha += str('%15.9f '  % self.strain[i][3])   #xy   
            linha += str('%15.9f '  % self.strain[i][5])   #xz
            linha += str('%15.9f '  % self.strain[i][3])   #yx
            linha += str('%15.9f '  % self.strain[i][1])   #yy
            linha += str('%15.9f '  % self.strain[i][4])   #yz
            linha += str('%15.9f '  % self.strain[i][5])   #zx
            linha += str('%15.9f '  % self.strain[i][4])   #zy
            linha += str('%15.9f '  % self.strain[i][2])   #zz 
            fout.write(linha+'\n')
      
        fout.write('$EndElementData\n')  


        fout.write('$ElementData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"S\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('9\n')
        fout.write(str("%5d  \n" % n_e)) # number of elements with results
        

        for i in range(0, n_e):
            
            linha = str('%5d  ' % (elements[i].id))
            linha += str('%15.9f '  % self.stress[i][0])   #xx        
            linha += str('%15.9f '  % self.stress[i][3])   #xy       
            linha += str('%15.9f '  % self.stress[i][5])   #xz   
            linha += str('%15.9f '  % self.stress[i][3])   #yx   
            linha += str('%15.9f '  % self.stress[i][1])   #yy   
            linha += str('%15.9f '  % self.stress[i][4])   #yz  
            linha += str('%15.9f '  % self.stress[i][5])   #zx   
            linha += str('%15.9f '  % self.stress[i][4])   #zy   
            linha += str('%15.9f '  % self.stress[i][2])   #zz 
            fout.write(linha+'\n')
      
        fout.write('$EndElementData\n')  


        fout.write('$ElementData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"vec1\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('3\n')
        fout.write(str("%5d  \n" % n_e)) # number of elements with results
        

        for i in range(0, n_e):
            
            linha = str('%5d  ' % (elements[i].id))
            linha += str('%15.9f '  % self.mainvec1[i][0])   #xx        
            linha += str('%15.9f '  % self.mainvec1[i][1])   #xy       
            linha += str('%15.9f '  % self.mainvec1[i][2])   #xz    
            fout.write(linha+'\n')
      
        fout.write('$EndElementData\n')  


        fout.write('$ElementData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"vec2\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('3\n')
        fout.write(str("%5d  \n" % n_e)) # number of elements with results
        

        for i in range(0, n_e):
            
            linha = str('%5d  ' % (elements[i].id))
            linha += str('%15.9f '  % self.mainvec2[i][0])   #xx        
            linha += str('%15.9f '  % self.mainvec2[i][1])   #xy       
            linha += str('%15.9f '  % self.mainvec2[i][2])   #xz    
            fout.write(linha+'\n')
      
        fout.write('$EndElementData\n')  
        fout.write('$ElementData\n')
        fout.write('1\n')     # number-of-string-tags 
        fout.write('\"vec3\"\n')   # name
        fout.write('1\n')     # number-of-real-tags
        fout.write('0.0\n')    # real - time value, por ex.
        fout.write('3\n')    # number of integer tags
        fout.write('1\n')
        fout.write('3\n')
        fout.write(str("%5d  \n" % n_e)) # number of elements with results
        

        for i in range(0, n_e):
            
            linha = str('%5d  ' % (elements[i].id))
            linha += str('%15.9f '  % self.mainvec3[i][0])   #xx        
            linha += str('%15.9f '  % self.mainvec3[i][1])   #xy       
            linha += str('%15.9f '  % self.mainvec3[i][2])   #xz    
            fout.write(linha+'\n')
      
        fout.write('$EndElementData\n')  