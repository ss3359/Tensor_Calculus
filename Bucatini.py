import random 
import numpy as np  
import pandas as pd 
import math 
import pygame
from sympy import symbols,diff,Matrix,cos,sin


# Constants 
c = 3*(10**8) # speed of light m/s
G = 6.674*(10**-11) # gravatational constant m^3/(kg s^2)



#Classes/Stucts 
class MATRIX:
    def __init__(self,r,c): 
        self.rows=r
        self.columns=c
        self.m=np.zeros((self.rows,self.columns))
    def __add__(self,other): 
        if(self.rows==other.rows and self.columns==other.columns): 
            result=Matrix(self.rows,self.columns)
            for i in range(self.rows):
                for j in range(self.columns):
                    result.m[i][j]= self.m[i][j]+other.m[i][j]
            return result
        else:
            print("The Dimensions Do Not Equal! ")

    def __mul__(self,other): 
        if(self.columns==other.rows):    
            result=Matrix(self.rows,other.columns)
            for i in range(result.rows):
                for j in range(result.columns):
                    c_ij=0
                    for k in range(self.columns): 
                        c_ij+=self.m[i][k]*other.m[k][j]
                    result.m[i][j]= c_ij
    
            return result
        else:
            print(f"The Dimensions Do Not Equal! ")
            print(f"Self.Columns: {self.columns}")
            print(f"Other.Rows: {other.rows}")

    def __sub__(self,other): 
        if(self.rows==other.rows and self.columns==other.columns): 
            result=Matrix(self.rows,self.columns)
            for i in range(self.rows):
                for j in range(self.columns):
                    result.m[i][j]= self.m[i][j]-self[i][j]
            return result
        else:
            print("The Dimensions Do Not Equal! ")

    def mul_scalar(self,num): 
        result=Matrix(self.rows,self.columns)
        for i in range(result.rows):
            for j in range(result.columns):
                c_ij=0
                for k in range(self.columns): 
                    c_ij=num*self.m[i][k]
                result.m[i][j]= c_ij
        return result 
    def Input_Values(self): 
        for i in range(0,self.rows): 
            for j in range(0,self.columns):
                self.m[i][j]=float(input())
    def Print_Matrix(self): 
        print("The Matrix Is: ")
        for i in range(0,self.rows): 
            for j in range(0,self.columns):
                print(self.m[i][j], end="\t")
            print()
    def Determinant(self):
        if(self.columns==1 and self.rows==1):
            return self.m[0][0]
        elif(self.columns==2 and self.rows==2): 
            return self.m[0][0]*self.m[1][1]-self.m[0][1]*self.m[1][0]
        det=0
        for i in range(self.columns): 
            SUB_MATRIX=Matrix(self.rows-1,self.columns-1)
            for j in range(1,self.rows):
                sub_col=0
                for k in range(self.columns):
                    if(k==i):
                        continue
                    SUB_MATRIX.m[j-1][sub_col]=self.m[j][k]
                    sub_col+=1
            det+=((-1)**i)*self.m[0][i]*SUB_MATRIX.Determinant()
        return det      

    def RowReduce(self): 
        Result=Matrix(self.rows,self.columns)
        pass


class Body:
    def __init__(self,t,r,theta,rho,m,x,y,z):
        self.t=t
        self.r=r
        self.theta=theta
        self.rho=rho
        self.mass=m
        self.x=x
        self.y=y
        self.z=z
        
    def Christosffel_Symbol(self): 
        gamma_tt_r=((G*self.mass)/(self.r**2))*(1-((2*G*self.mass)/(c*c*self.r)))
        gamma_rr_r=-((G*self.mass)/(self.r**2))*(1-((2*G*self.mass)/(c*c*self.r)))**-1
        gamma_phiphi_r=-self.r*(1-((2*G*self.mass)/(c*c*self.r)))
        
        return [gamma_tt_r,gamma_rr_r,gamma_phiphi_r]
        
    def SchwarschildRadius(self):
        return 2*G*self.mass/(c*c)
    def SchwartzchildMetric(self):
        pass
    def UpdatePosition(self): 
        pass


class Manifold:
    pass
class Sphere: 
    def __init__(self,x,y,z):
        self.rho=np.sqrt(x**2 + y**2 + z**2)
        self.phi=np.arccos(1/self.rho)
        self.theta=np.arctan(y/x)

    def Jacobian(self,x,y,z): 
        drho_dx=x/self.rho
        drho_dy=y/self.rho
        drho_dz=z/self.rho
        print(f"({drho_dx}, {drho_dy}, {drho_dz})")
        
class Saddle: 
    def __init__(self,x,y,z): 
        self.u=x 
        self.v=y
        self.w=x*y
    def Jacobian(self,x,y,z): 
        pass
    
class Cylinder: 
    def __init__(self,x,y,z):
        self.r=np.sqrt(x**2 + y**2)
        self.theta=np.arctan(y/x)
        self.z=z
    def Jacobian(self,x,y,z):
        drdx=x/self.r
        drdy=y/self.r
        drdz=0

       
        

#Functions Not Belonging to classes
def R(): 
    u,v =symbols('u v')
    return Matrix([
        [u*cos(v)],
        [u*sin(v)],
        [u]])     #This is a column vector
def PartialR_u():
    
    u,v =symbols('u v')
    f = R()
    return diff(f,u)
def PartialR_v():
    u,v =symbols('u v')
    f = R()
    return diff(f,v)
def MetricTensor(): 
     u,v =symbols('u v')
     R_u=PartialR_u()
     R_v=PartialR_v()
     return Matrix([[R_u.T*R_u, R_u.T*R_v],
            [R_v.T*R_u,R_v.T*R_v]])

def main(): 
    x=3
    y=4
    z=5
    coord_sph=Sphere(x,y,z)
    coord_sph.Jacobian(x,y,z)
  

main()






































'''
Excess Code
    Sun=Body(0,1.496e11, 0,0,1.989e30) # Sun at Earths Orbit

    print(f"Christoffell Values: {Sun.Christosffel_Symbol()}")
    print(f"Schwartzchild Radius (Sun): {Sun.SchwarschildRadius()}")

      r=2
    c=2
    A=Matrix(r,c)
    B=Matrix(r,c)
    RESULT=Matrix(r,c)
    RESULT_2=Matrix(r,c)
    print("Enter Values for A")
    A.Input_Values()
    print("Enter Values For B")
    B.Input_Values()

    print()
    RESULT=A*B
    RESULT_2=A+B
    RESULT.Print_Matrix()
    print()
    RESULT_2.Print_Matrix()
    
       # r=3
    # c=3
    # A=Matrix(r,c)
    # print("Enter Values for A")
    # A.Input_Values()

    # print(f"Determinant of A is {A.Determinant()}")

   print(f"DR/du ={PartialR_u()}")
   print(f"DR/dv = {PartialR_v()}")
   print(f"Metric Tensor: \n {MetricTensor()}")

   class Tensor:
    def __init__(self,A):
        self.A=A
'''