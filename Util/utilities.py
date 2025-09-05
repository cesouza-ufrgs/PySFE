


from math import pi, pow, sqrt, atan, cos, sin, acos

def tela_regua():
    
    #print '************************************'
        
    print ('_______________________________________________')
    return

def calc_v_res(v_a,v_b):
    
    
    return sqrt(v_a*v_a + v_b*v_b)


def align_header(text,leng):
    
    
    return text.center(leng)+ '|' 




def format_header(header,data,text,value):
    
    header +=  align_header(text,9)
    
    data  += ' %7.5f |' % value
    
    return (header, data)