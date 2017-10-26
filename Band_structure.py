#-----------------------------------------------------------------------------------------------------------------------------
# UPDATE LOG
#-----------------------------------------------------------------------------------------------------------------------------

# 10/24 Just fixed the eigenvalue matrix: values are now about the right order of magnitude
#
#
#
#
#
#
#
#
#



#-----------------------------------------------------------------------------------------------------------------------------
# MODULE IMPORT PROCEDURE
#-----------------------------------------------------------------------------------------------------------------------------

# -*- coding: utf-8 -*-
from decimal import*        #don't mess with any of this, or stuff will not work properly
from visual import*
from visual.graph import *
import numpy
from numpy.linalg import det
from sympy import*

#-----------------------------------------------------------------------------------------------------------------------------
# A BUNCH OF POORLY ORGANIZED VARIABLES
#-----------------------------------------------------------------------------------------------------------------------------

G1 = vector(10e10,0,0) # G1 designates the first reciprocal lattice vector

U = 10e-19   # U1 designates the first fourier component of the expansion of the periodic potential function U(r)

h_bar = 1.05e-34 # h_bar designates the reduced Planck constant

m = 9.11e-31 # m designates the mass of an electron in the lattice

c0 = (h_bar**2)/(2*m) # c0 is h_bar^2 over 2m

BZ_boundary = int(0.5*G1.mag) # self explanitory

k = 0.5*10e10 # k designates the wavevector

T = 3 # T designates the number of fourier coefficients used to approximate the eigenvalues.  It also gives the # of Eigenvalues

a = [] # a is the temp array variable that will later be turned into the matrix that contains the set of equations used to solve for E

#-----------------------------------------------------------------------------------------------------------------------------
# SOLVING FOR EIGENVALUES
#-----------------------------------------------------------------------------------------------------------------------------

for i in range (T):                 # this whole thing is extremely convoluted, don't bother too much with it
    a.append([])
    for j in range (T):
        a[i].append(0)
 
for i in range (T):                 # this is the addition of the wavefunction coefficients to the matrix

    if(i==0):                       # special case for i=0
       a[i][i]= c0*(k-i*G1.mag)**2-E
    elif(i==1):
         a[i][i]= c0*(k+i*G1.mag)**2-E
    elif(i%2==0):                   # special case for i=even
       a[i][i]= c0*(k+(1-i)*G1.mag)**2-E
    else:                           # special case for i=odd
       a[i][i]= c0*(k+(i-1)*G1.mag)**2-E

for i in range (T):                 # this is the addition of the potential energy coefficients to the matrix

    if(i==0):                       # special case for i=0
        a[i][1]=U
        if(T>2): 
            a[i][2]=U
    elif(i==1):                     # special case for i=1
        a[i][i-1]=U
        if(T>3):
           a[i][i+2]=U
    elif(i+2<T):                    # general case for i=n
        a[i][i+2]=U
    if(i-2>=0):
        a[i][i-2]=U

central_equation_set= Matrix(a) # this is the matrix that contains all the set of equations
central_equation_set_det = central_equation_set.det() # this takes the determinant of the matrix

print central_equation_set
print ""
print central_equation_set_det

#-----------------------------------------------------------------------------------------------------------------------------
# PLOTTING OF FREE ELECTRON APPROXIMATION
#-----------------------------------------------------------------------------------------------------------------------------

energyplot = gdisplay(x=0,y=200,width=600,height=1500,xmin=0,xmax=40e10,ymin=0,ymax=1e-16,title="Band Structure",xtitle="Wavevector",ytitle="Energy")

free_electron1 = gcurve(color=color.white)                  # graphing of the free electron polts about a few RL points
free_electron2 = gcurve(color=color.white)
free_electron3 = gcurve(color=color.white)
free_electron4 = gcurve(color=color.white)

k = -2e10

while(k<=4*BZ_boundary and k>=-4*BZ_boundary):
    
     free_electron1.plot(pos=(k-G1.mag,(c0)*(k**2)))
     free_electron2.plot(pos=(k,(c0)*(k**2)))
     free_electron3.plot(pos=(k+G1.mag,(c0)*(k**2)))
     free_electron4.plot(pos=(k+2*G1.mag,(c0)*(k**2)))
     k = k+1e7
     
#-----------------------------------------------------------------------------------------------------------------------------
# PLOTTING OF BAND STRUCTURE
#-----------------------------------------------------------------------------------------------------------------------------

valence_band = gcurve(color=color.green)        
conduction_band = gcurve(color=color.red)

k = 0
while(k>=0 and k<=2*BZ_boundary):
    
    Eigenvalue1 = 0.5*(c0*(k-G1.mag)**2+c0*k**2)-(0.25*(c0*(k-G1.mag)**2-(c0*k**2))**2+U**2)**(0.5)
    valence_band.plot(pos=(k,Eigenvalue1))
    k = k+1e7

k = 0
while(k>=0 and k<=2*BZ_boundary):
 
    Eigenvalue2 = 0.5*(c0*(k-G1.mag)**2+c0*k**2)+((0.25*(c0*(k-G1.mag)**2-c0*k**2)**2)+(U**2))**(0.5)
    conduction_band.plot(pos=(k,Eigenvalue2))
    k = k +1e7

valence_band = gcurve(color=color.green)        #this looks redundant, but it actually is necessary for a reason I don't remember
conduction_band = gcurve(color=color.red)

k = G1.mag
while(k>=G1.mag and k<=2*BZ_boundary+G1.mag):
    
    Eigenvalue1 = 0.5*(c0*(k-G1.mag-G1.mag)**2+c0*(k-G1.mag)**2)-(0.25*(c0*(k-G1.mag-G1.mag)**2-(c0*(k-G1.mag)**2))**2+U**2)**(0.5)
    valence_band.plot(pos=(k,Eigenvalue1))
    k = k+1e7

k = G1.mag
while(k>=G1.mag and k<=2*BZ_boundary+G1.mag):
 
    Eigenvalue2 = 0.5*(c0*(k-G1.mag-G1.mag)**2+c0*(k-G1.mag)**2)+((0.25*(c0*(k-G1.mag-G1.mag)**2-c0*(k-G1.mag)**2)**2)+(U**2))**(0.5)
    conduction_band.plot(pos=(k,Eigenvalue2))
    k = k +1e7

