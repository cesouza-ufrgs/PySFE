#!/usr/bin/env python3
"""
@mainpage PyFE

This is the main page of PyFE - a Python code for Finite Element teaching. 

@section Goals

initially, only bar and beam models

@section Todo

@todo
 - 19.04.10 - everythinh


@section Teste


@package PyEolica

@author carlos eduardo de souza
@date 19.04.10
@version 19.04.10
@version   19.04.10
           
           
                     
            
@attention: esse arquivo ainda esta com muita cara de c++! Como usar comandos em python?
@note 
    Esse programa deve se manter sem interface grafica. 
    Usar sempre o conceito: rodar sozinho!!!
    Posteriomente, devo criar uma interface que chame as rotinas atuais.
    
@note 
Para documentacao, ver o site <a href="http://notemagnet.blogspot.com.br/2009/10/using-doxypy-for-python-code.html">Using doxypy for Python code documentation </a>    

"""
 
import sys # os, codecs , fileinput
# import string
# import math
import numpy as np
from numpy import linalg as LA

from math import pi,atan,log


# classes locais -------------
from Data.class_data   import *
#from Cargas.calcula_cargas   import *
#from BEM.calcula_bem   import *
# from ClassesGlobais.class_resultados import *
from Graphics.graphics import *
from Util.utilities import *


from Bar.FEBar import *
#------------------------------------------------------------------------------|
#------------------------------------------------------------------------------|




 
#--- calcula vento--------------------------------------------------------------
 
def calcula_vento(dados):
    '''
    @summary  calcula velocidade de vento necessaria para atender a potencia 
              requerida e dimensoes
              
              
    @version  14.02.22 
    
    @todo    Tudo.
 
    @param dados: class C_Dados input data structure
    '''
     
     
    return



    
#--- main ----------------------------------------------------------------------
if __name__ == '__main__':
    
    '''
    @brief program entry. The specific code is selected from the inputfile.
    @date 19.04.11
       
    '''
    print ('\n\n\n')   
    tela_regua()    
    print ('$                   PyFE                      $')  
    print ('$                  v 19.04                    $') 
    tela_regua()     
    print ('$  a simplified finite element code           $')
    tela_regua()
    #print sys.argv
    
    ## number of command line arguments
    n_arg =  len(sys.argv)  
        
    
    if n_arg > 1:
        
        #print sys.argv[1:]
        #print sys.argv[1]        
        ## name of input data file.
        inputfile = sys.argv[1]  
                
        tela_regua()
        print ('A FAZER: \n')
        print (' 19.04.10 - everything')
        #desenha_grafico()                   
        tela_regua()
        
        
        # input data 
        data = C_Data()
        data.read_file(inputfile)
        
        
        if data.analysis == 'static':
            
            if 'bar' in data.elemtype:
                FEBar(data)
        
        
        
        
    else:
        print ('\n**********************************\n')
        print ('O programa PyEolica exige a definicao de um arquivo de entrada de dados')
        print ('Execute fazendo: >>./PyEolica.py arquivo.dat')
        print ('\n saindo \n ')
        print ('\n**********************************\n')
            

    print ('\n**********************************\n')
