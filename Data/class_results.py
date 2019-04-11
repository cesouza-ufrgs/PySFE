 

 
import sys, os, codecs , fileinput

from class_data import *

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
        
        self.S = []
        
        
        return
        
    def save_to_gmsh(self,path,filebase,nodes,elements):  
        '''
        funcao para salvar parametros calculados em um arquivo.    
        '''
    
        print 'arquivo input: ', path, filebase
        #fi = open(arquivo_input, 'r');
       
        n_n = len(nodes)
        n_e = len(elements)
        
        file_out_gmsh = path + '/' + filebase + '_gmsh.msh'
        
        
        print file_out_gmsh
        
        fout = open(file_out_gmsh, 'w')
    
    
        fout.write('$MeshFormat\n')    
        fout.write('2.2 0 4\n')    
        fout.write('$EndMeshFormat\n')  
            
      
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
        
        # ----------------------------------------      
        fout.close()    
    
        return
    
    
    
    
    