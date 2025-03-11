import math 
import random 
import numpy as np
import pandas as pd 
from sympy import symbols,diff, Matrix,cos,sin,DotProduct, simplify
from sympy import Inverse

#Constants 


def R():
    u,v=symbols("u v")
    return np.array([cos(v)*sin(u),sin(v)*sin(u),cos(u)])

def PartialDerivative(): 
    U,V=symbols("u v")
    r=R()
    F=r[0]
    G=r[1]
    H=r[2]
    dRdu=np.array([diff(F,U), diff(G,U), diff(H,U)])
    dRdv=np.array([diff(F,V), diff(G,V), diff(H,V)])

    print(f"dR/du={dRdu}",end='\n')
    print(f"dR/dv={dRdv}",end="\n")
    return np.array([dRdu,dRdv])


def SecondPartialDerivative(): 
    U,V=symbols("u v")
    dR=PartialDerivative()

    dRdu=dR[0]
    dRdv=dR[1]
    
    d2Rdu2=np.array([diff(dRdu[0],U), diff(dRdu[1],U),diff(dRdu[2],U)])
    d2Rdudv=np.array([diff(dRdu[0],V), diff(dRdu[1],V),diff(dRdu[2],V)])
    d2Rdvdu=np.array([diff(dRdv[0],U), diff(dRdv[1],U),diff(dRdv[2],U)])
    d2Rdv2=np.array([diff(dRdv[0],V), diff(dRdv[1],V),diff(dRdv[2],V)])

    print("\n")
    print(f"d2R/du2={d2Rdu2}",end='\n')
    print(f"d2R/dudv={d2Rdudv}",end="\n")
    print(f"d2R/dvdu={d2Rdvdu}",end="\n")
    print(f"d2R/dv={d2Rdv2}",end="\n")
    print('\n')
    return np.array([d2Rdu2,d2Rdudv,d2Rdvdu,d2Rdv2])

def MetricTensor():
    U,V=symbols('u v')
    R=PartialDerivative()
    drdu=R[0]
    drdv=R[1]

    g = Matrix([[simplify(np.dot(drdu,drdu)),simplify(np.dot(drdu,drdv))],
                  [simplify(np.dot(drdv,drdu)),simplify(np.dot(drdv,drdv))]])
    print("Metric Tensor: ", end='\n')
    print(g,end="\n") 

    return g

def ChristoffelSymbols():
    U,V=symbols('u v')
    dR=PartialDerivative()
    d2R=SecondPartialDerivative()

    #First Derivatives 
    drdu=dR[0]
    drdv=dR[1]

    #Second Derivatives
    d2Rdu2=d2R[0]
    d2Rdudv=d2R[1]
    d2Rdvdu=d2R[2]
    d2Rdv2=d2R[3]
    

    g=MetricTensor()
    gInv=g.inv()

    gamma_uu_u=simplify(np.dot(d2Rdu2,drdu)*gInv[0]+np.dot(d2Rdu2,drdv)*gInv[2])
    gamma_uu_v=simplify(np.dot(d2Rdu2,drdu)*gInv[1]+np.dot(d2Rdu2,drdv)*gInv[3])

    gamma_uv_u=simplify(np.dot(d2Rdudv,drdu)*gInv[0]+np.dot(d2Rdudv,drdv)*gInv[2])
    gamma_uv_v=simplify(np.dot(d2Rdudv,drdu)*gInv[1]+np.dot(d2Rdudv,drdv)*gInv[3])

    gamma_vu_u=simplify(np.dot(d2Rdvdu,drdu)*gInv[0]+np.dot(d2Rdvdu,drdv)*gInv[2])
    gamma_vu_v=simplify(np.dot(d2Rdvdu,drdu)*gInv[1]+np.dot(d2Rdvdu,drdv)*gInv[3])

    gamma_vv_u=simplify(np.dot(d2Rdv2,drdu)*gInv[0]+np.dot(d2Rdv2,drdv)*gInv[2])
    gamma_vv_v=simplify(np.dot(d2Rdv2,drdu)*gInv[1]+np.dot(d2Rdv2,drdv)*gInv[3])

    print("gInv=",gInv,end="\n\n")
    
    CH_Symbols=Matrix([gamma_uu_u,gamma_uu_v,gamma_uv_u,gamma_uv_v,gamma_vu_u,gamma_vu_v,gamma_vv_u,gamma_vv_v])
    print("CH_Symbols",CH_Symbols,end="\n\n")

def GeodesicEquations():
    uτ=2
    vτ=math.pi/4
    dudτ=0.5
    dvdτ=0.5
    dτ=0.1

    U,V=symbols('u v')
    dR=PartialDerivative()
    d2R=SecondPartialDerivative()
    CH_Symbols=ChristoffelSymbols()
    

    d2udτ2=-((CH_Symbols[0]*(dudτ)*(dudτ))
            +(CH_Symbols[2]*(dudτ)*(dvdτ))
            +(CH_Symbols[4]*(dudτ)*(dvdτ))
            +(CH_Symbols[6]*(dvdτ)*(dvdτ)))
    d2vdτ2=-(CH_Symbols[1]*(dudτ)*(dudτ)
            +(CH_Symbols[3]*(dudτ)*(dvdτ))
            +(CH_Symbols[5]*(dudτ)*(dvdτ))
            +(CH_Symbols[7]*(dvdτ)*(dvdτ)))
        
    #Inplement the Runge-Kutta Method of Order 4 to Solve The Geodesic Equations
    


  
def main():
    pass
    
main()

'''
Excess Code: 

MetricTensor()
    print('\n')
    ChristoffelSymbols()
'''