"""
@file class_dados

@brief Definition of class and commands for input data parsing.

"""
 
 
import  os, math
import numpy as np
from Util.utilities import align_header, ruler

def extrai_comando(linha):
    '''
    extracts a command from a line of the input data.
    
    @param linha: string, line read from the input data file.
    '''

 
    pos_virgula = -1
    comando = ' '
    virg = ","
    
    if virg in linha:
                    
        pos_virgula = linha.index(virg)
        
        
        
        comando = linha[0:pos_virgula]


        print ('comando: '+ comando)    

    return comando

#-----------------------------------------------------------
'''
@brief - extrai um comando de uma linha. Eh a primeira palavra antes do @

@note: essa funcao deve ser passada para um outro input_file/biblioteca/?

@attention: deve existir uma maneira mais 'pythonica' de fazer isso!
'''
def separa_linha(linha):

 
    pos_virgula = -1
    comando = ' '
    parametro = ' '
    virg = ","
       
    
    if not '$' in linha[0] or '#' in linha[0]:
    
        if virg in linha:
                        
                        
                        
            pos_virgula = linha.index(virg)        
            
            
            comando = linha[0:pos_virgula]
            parametro = linha[pos_virgula+1:-1]
            
            comando.strip()
            parametro.lstrip()
            parametro.rstrip()
    
            #print 'comando:'+ comando + '*param:*' + parametro +'*' 
            #print len(parametro)
            
            if len(parametro)<2:
                print (' ATENCAO: \n') 
                print ('  erro no parametro')
                print (' ver linha do comando ' + comando)
                print (' ATENCAO!!!\n' )
                parametro = "  0   "
    

    return comando,parametro


'''--------------
'''
def le_comando(linha):


    virg = ","
       
    comando = []
    
    
    if not ('$' in linha[0] or '#' in linha[0] or '%' in linha[0]):
    
        #print linha
        if virg in linha:
                        
            cmd = linha.split(virg)
            for c in cmd:
                #print c
                comando.append(c.strip())
        else:
            comando.append(linha.strip())

            #print 'separando o comando'
            #print 'antes ', cmd
            #print 'depois', comando
            
    
    return comando


class C_Nodes():
    
    def __init__(self):
        self.id = 0
        
        self.coord = []

class C_Elements():
    
    def __init__(self):
        self.id = 0
        
        self.nodes = []
        
        self.prop = []
        

class C_Materials():
    
    def __init__(self):
        self.id = 0
        
        self.type = 'isotropic'
        
        self.E = 0.
        
        self.nu = 0.       
        
        self.G = 0.

        self.rho = 0.       

        self.sig_e = 0.       

        
class C_Sections():
    
    def __init__(self):
        self.id = 0
        
        self.type = 'bar'
        
        self.imat = 0
        
        self.A = 0.
        
        self.Ixx = 0.       
        
        self.Iyy = 0.

        self.J0 = 0.    

        self.e = 0.   
  

class C_Forces():
    
    def __init__(self):

        self.id = 0
        
        self.dir = 0
        
        self.A = 0.  


#------------------------------------------------------------------------------- 
 
class C_Data():
    '''
    @class C_Data 
    
    @brief Input data  
    
    Aqui estao todos os tipos de dados, entao eh importante tomar cuidado com 
    o uso.
    
    @todo 
    - listar todas as variaveis
    - 
    
    '''
    
    def __init__(self):
        '''
        Init function
        '''
        
        
        ## defines the type of analysis
        self.analysis = 'static'
        
        # options= modal, response  
        
        ## defines if various graphs are to be plotted during the simulation
        self.plotfactor = 0.0
        
        self.ndof = 2
        
        self.elemtype = 'none'        
        
        #----------------------------------------------
        # contantes
        self.grav = 9.81  # aceleracao da gravidade
        
        
        self.n_nodes = 0
        self.nodes = []
        self.elements = []
        self.materials = []
        self.sections = []
        self.section2d = []
        
        self.forces = []
        self.bconditions = []
        
        #----------------------------------------------
        # outras variaveis sao definidas na leitura e outros locais
        # WARNING: preciso listar todas aqui! 14.03.25
    
        return
    
    
    def v_ang(self):
        '''
        Returns the angular velocity at the tip
        '''
        v = self.v_ave * self.r        
        
        return v
    
    
    def read_file(self,input_file):
        '''
        @function read_file
        @brief Abre o input_file e read_file os dados la contidos.
        
        
        @param input_file - caminho para o input_file de dados.
        '''
        ruler()
        print ('        Parsing input file')
        
        # separa o nome do input_file
        (diretorio, nome) = os.path.split(input_file)
                
        # verifica a ext do input_file        
        (filebase,ext)  = os.path.splitext(nome)   
        
                       
        
        print ('path:          '+ diretorio)
        print ('file name:     '+ nome)
        print ('file only:     '+ filebase )
        print ('file extension:'+ ext  )            
      
        self.input_file = input_file
        self.path       = diretorio
        self.filename   = nome
        self.filebase   = filebase
        self.ext        = ext
       
      
        if nome =='':
            return
        
        
        print ('opening  file  for reading:\n' +  input_file )
        
          
        #--- inicia a leitura do input_file de dados
        # em 13.11.22
        #print "Name of the file: ", fo.name

        
        
        #--- reading looping 
        with open(input_file) as fileobject:  
            
            read_nodes = False
            read_elements = False
            read_materials = False
            read_sections = False
            read_section2d = False
            read_section3d = False
            read_bconditions = False
            read_forces = False

            current_element = None
            
         
            for in_line in fileobject:  
    
                #line = fo.readline()
                # print "Read Line: %s" % (in_line)    
                                            
                
                   
                #---   read nodes
                if read_nodes:
                    #print 'reading nodes'
                    
                    cmds = le_comando(in_line)
                                                
                    if "end" in cmds[0]:
                        read_nodes = False
                    else:                             
                        node = C_Nodes()                       
                             
                        node.id =   int(cmds[0])
                        node.coord = np.array([float(cmds[1]), float(cmds[2]) , float(cmds[3])   ])
                        
                        self.nodes.append(node)  
                             
                #-----    read elements
                if read_elements:
                    #print 'reading elements'
                    


                    cmds = le_comando(in_line)
                        
                        
                    if "end" in cmds[0]:
                        read_elements = False
                    else:                         
                        elem = C_Elements()                       
                        
                        elem.prop =   1     
                        elem.id   =   int(cmds[0])

                        n_nodes = len(cmds)-1
                                                
                        for i in range(1,n_nodes+1):
                            elem.nodes.append(int(cmds[i]))
                        
                        self.elements.append(elem)  
                             
                #-----    read materials
                if read_materials:
                    #print 'reading elements'
                    
                    cmds = le_comando(in_line)
                        
                    if "end" in cmds[0]:
                        read_materials = False
                    else:                         
                        mat = C_Materials()                       
                             
                        mat.id     = int(cmds[0])
                        mat.E      = float(cmds[1])
                        mat.nu     = float(cmds[2])
                        mat.rho    = float(cmds[3])
                        mat.sig_e  = float(cmds[4])                       
                        
                        
                        self.materials.append(mat)  
                                       
                #-----    read sections
                if read_sections:
                    #print 'reading elements'
                    
                    cmds = le_comando(in_line)
                        
                    if "end" in cmds[0]:
                        read_sections = False
                    else:                         
                        sec = C_Sections()                       
                             
                        sec.id    = int(cmds[0])
                        sec.imat  = int(cmds[1])
                        sec.A     = float(cmds[2])
                        sec.Ixx   = float(cmds[3])
                        sec.Iyy   = float(cmds[4])
                        if len(cmds)>5:
                            sec.J0    = float(cmds[5])                       
                                                
                        self.sections.append(sec)                           
                #-----    read sections
                if read_section2d:

                    # print('reading section2d for plane elements ')
                    
                    cmds = le_comando(in_line)
                        
                    if "end" in cmds[0]:
                        read_section2d = False
                    else:                         
                        sec = C_Sections()                       
                             
                        sec.id    = int(cmds[0])
                        sec.imat  = int(cmds[1])
                        sec.e     = float(cmds[2])                      
                                                
                        self.sections.append(sec)  
                        
                #-----    read forces
                if read_forces:
                    
                    cmds = le_comando(in_line)
                        
                    if "end" in cmds[0]:
                        read_forces = False
                    else:                       
                        F = C_Forces()                         
                             
                        F.id    = int(cmds[0])
                        F.dir   = int(cmds[1])
                        F.A     = float(cmds[2])
                                                                        
                        self.forces.append(F)  
                        

                #-----    read boundary cond
                if read_bconditions:
                    
                    cmds = le_comando(in_line)
                        
                    if "end" in cmds[0]:
                        read_bconditions = False
                    else:                       
                        bc = C_Forces()                         
                             
                        bc.id    = int(cmds[0])
                        bc.dir   = int(cmds[1])
                                                                        
                        self.bconditions.append(bc)      
                        
                        
                                            
                
                #--- separa a in_line em 'card' e 'parametros'
                card_line = le_comando(in_line)
                
                #print 'card' , card_line
                
                if len(card_line)>0:                
                    card = card_line[0]
                else:                
                    card = 'nada'
                      
                if "analysis" in card: 
                    self.analysis = card_line[1]
                    
                if "plotfactor" in card:                     
                    self.plotfactor = float(card_line[1])
                    
                if "ndof" in card:                     
                    self.ndof = int(card_line[1])
                    
                if "elemtype" in card:                     
                    self.elemtype = card_line[1]
                       
                #---- -FE -----        
                elif  card == "nodes" :  
                    
                    read_nodes = True  
                    
                elif  card == "elements" :  
                    
                    read_elements = True
                    
                elif  card == "materials" :  
                    
                    read_materials = True   
         
                elif  card == "sections" :  
                    
                    read_sections = True   

                elif  card == "section2d" :  
                    
                    read_section2d = True  
           
           
                elif  card == "forces" :  
                    
                    read_forces = True   
                    
                elif  card == "bconditions" :  
                    
                    read_bconditions = True
                 
                    
                    
        #--- verifying stuff            
        #mensagens:
        
        print ('Reading completed.')
        ruler()
        
        
        self.n_nodes = len(self.nodes)
        self.n_elem = len(self.elements)
        self.n_mat = len(self.materials)
        self.n_sec = len(self.sections)
        self.n_forces = len(self.forces)
        self.n_bconditions = len(self.bconditions)
        
        show_log = False
        
        if show_log:
            print ('Log:')
            #--------------            
            print ('n_nodes' , self.n_nodes)
            for n in self.nodes:
                print (n.id , n.coord)
                
                
            print ('n_elem' , self.n_elem)
            for e in self.elements:
                print (e.id , e.prop, e.nodes)                
                
                
            print ('n_mat' , self.n_mat)
            for m in self.materials:
                print (m.id , m.E, m.rho)
            
                
            print ('n_mat' ,  self.n_mat)
            for m in self.materials:
                print (m.id , m.E, m.rho)
                
            print ('n_sec' ,  self.n_sec)
            for s in self.sections:
                print (s.id , s.A, s.Ixx, s.e)

        
        #self.v_ave = self.v_des / 1.4
        #--------------

        print ('           end of reading')
        ruler()
        return
