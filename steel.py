import numpy as np

def raynor(Es,fy,fsu,esh,esu,dels,C1,Ey):
    es  = np.linspace(0,esu,int(esu/dels+1))
    fs  = es*0
    ey  = fy/Es;
    fsh = fy + (esh-ey)*Ey
    
    for i in range(len(es)):
        if es[i] < ey:
            fs[i] = Es*es[i]
            
        if es[i] >= ey and es[i] <= esh:
            fs[i] = fy+(es[i]-ey)*Ey

        if es[i] > esh:
            fs[i] = fsu-(fsu-fsh)*(((esu-es[i])/(esu-esh))**C1)
            
    return es, fs
            
def king(Es,fy,fsu,esh,esu,dels):
    r = esu - esh
    m = ((fsu/fy)*((30*r+1)**2)-60*r-1)/(15*(r**2))
    es = np.linspace(0,esu,int(esu/dels+1))
    fs  = es*0
    ey = fy/Es
    
    for i in range(len(es)):
        if es[i] < ey:
            fs[i] = Es*es[i]
        if es[i] >= ey and es[i] <= esh:
            fs[i] = fy
        if es[i] > esh:
            fs[i] = ((m*(es[i]-esh)+2)/(60*(es[i]-esh)+2) + (es[i]-esh)*(60-m)/(2*((30*r+1)**2)))*fy
            
    return es, fs