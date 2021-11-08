# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 11:24:13 2021

@author: dell
"""

import numpy as np
import sympy as sym
##g is your metric,lv means vector of your variables
#g是你的度规
c,R,rho=sym.symbols('c R rho')
t,r,the,phi=sym.symbols("t r theta phi")
lv=sym.Matrix([t,r,the,phi])
#define some functions of your variable
#lv是你的变量
a=sym.Function("a")(t)
#an test metric
#g=sym.Matrix([[-c**2*(1-a/r),0,0,0],[0,1/(1-a/r),0,0],[0,0,r**2,0],[0,0,0,(r*sym.sin(the))**2]])
g=sym.Matrix([[-c**2,0,0,0],[0,a**2/(1-rho*r**2),0,0],[0,0,a**2*r**2,0],[0,0,0,(a*r*sym.sin(the))**2]])
g_inv=g.inv()

longt=len((g[:,0]))
gamma=sym.zeros(longt,longt**2)
##calculate Christoffel Symbols \Gamma^{i}_{jk}
#print Christoffel Symbols that is not 0
#输出不为0的克氏符号,第一个数字就是上标
for i in range(len(g[:,0])):
   for j in range(len(g[:,0])):
       for k in range(len(g[:,0])):
           h=0
           for lam in range(len(g[:,0])):                      
              gam=h+0.5*(g_inv[i,lam])*(sym.diff(g[lam,k],lv[j])+sym.diff(g[j,lam],lv[k])-sym.diff(g[j,k],lv[lam]))
              h=gam
           if(gam!=0):
            print(i,j,k) 
            gam=sym.simplify(gam)
            gamma[i,longt*j+k]=gam
            sym.print_latex(gam)

            
#calculate Riemanian curvature tensor R_{abcd}
#计算不为零的黎曼曲率张量
Rem=sym.zeros(longt**2,longt**2)
Rem_up=sym.zeros(longt**2,longt**2)
K=0
for i in range(len(g[:,0])):
   for j in range(len(g[:,0])):
       for k in range(len(g[:,0])):
           for l in range(len(g[:,0])):
               h=0
               Re=0
               
               for n in range(len(g[:,0])):
                  for m in range(len(g[:,0])):
                      h=h+g[n,m]*(gamma[n,l*longt+i]*gamma[m,j*longt+k]-gamma[n,k*longt+i]*gamma[m,j*longt+l])
               Re=1/2*(sym.diff(sym.diff(g[i,l],lv[j]),lv[k])+sym.diff(sym.diff(g[j,k],lv[i]),lv[l])-sym.diff(sym.diff(g[j,l],lv[i]),lv[k])-sym.diff(sym.diff(g[i,k],lv[j]),lv[l]))+h
               if (Re!=0):
                   #print riemanian curvature tensor that is not 0
                   
                   print(i,j,k,l)
                   Re=sym.simplify(Re)
                   Rem[i*longt+j,k*longt+l]=Re
                   Rem_up[i*longt+j,k*longt+l]=Re*g_inv[i,i]*g_inv[j,j]*g_inv[k,k]**g_inv[l,l]
                   #some dirty methods to calculate Kretschmann scalar
                   K=Re*Re*g_inv[i,i]*g_inv[j,j]*g_inv[k,k]*g_inv[l,l]+K
                   K=sym.simplify(K)
                   sym.print_latex(Re)
                 
#print Kretschmann scalar after "ok"
#计算Kretschmann标量
print("ok")  
K=sym.simplify(K)              
sym.print_latex(K)
#妈的，手写作业有毒