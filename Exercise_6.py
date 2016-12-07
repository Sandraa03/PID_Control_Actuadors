
# coding: utf-8

# In[1]:

#CONTROL KP/KD (x³)-(2*x²)+((1.00002+(0.02*kd))*x)-(0.02*kd)=0
#Import library
from math import sqrt, acos, cos, pi
import cmath
import numpy as np
import pylab as pl
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.transforms as mtransforms
import matplotlib.text as mtext
import matplotlib.pyplot as polar


import _tkinter


# In[3]:

#Definición función para utilizar una raíz cubica negativa, python de entrada no permite calcular raices cubicas negativas.
def s3(a):
    if a<0:
        return -abs(a)**(1/3.)
    else: return a**(1/3.)


# In[49]:

#Declaración de variables del control PID
thau=0.001
gamma=20.
kp=10.
kd=np.linspace(1., 5., num=10)

#Polinomio de 3er grado para el control, ax**3+bx**2+cx+d=0
#Definimos las componentes

a = 1.
b = -2.
c = 1+gamma*kp*thau**2+gamma*thau*kd
d = -gamma*thau*kd


# In[50]:

#Definición del método de Cardano para python para resolver polinomios de 3er grado, 
#establecemos Discriminante (D) y decidimos en el If donde entramos para calcular.

def Cardano(a,b,c,d):
        
    p = c-b**2/3.
    q = (1/27.)*(27*d-9*b*c+2*b**3)
    
    D = (p/3.)**3+(q/2.)**2
    
    if D>0:
        A = s3(-q/2.+sqrt(D))
        B = s3(-q/2.-sqrt(D))
        y1 = A + B
        y2 = -(1/2.)*y1 + 1j*(sqrt(3)/2.)*(A-B)
        y3 = -(1/2.)*y1 - 1j*(sqrt(3)/2.)*(A-B)
    if D==0:
        if q==0: 
            if p==0:
                y1=0
                y2=0
                y3=0
        else:
            y1=A+B
            y2=-0.5*y1
            y3=y2
    if D<0:
        a=acos(sqrt((q**2/4.)/(-p**3/27.)))
        if q<=0:
            y1=2*sqrt(-p/3.)*cos(a/3.)
            y2=2*sqrt(-p/3.)*cos((a+2*pi)/3.)
            y3=2*sqrt(-p/3.)*cos((a+4*pi)/3.)
        else:
            y1=-2*sqrt(-p/3.)*cos(a/3.)
            y2=-2*sqrt(-p/3.)*cos((a+2*pi)/3.)
            y3=-2*sqrt(-p/3.)*cos((a+4*pi)/3.)

    x1 = y1 - b/3.
    x2 = y2 - b/3.
    x3 = y3 - b/3.
    
    return x1,x2,x3


# In[51]:

#Creamos arrays para almacenar los datos x1,x2,x3 en función de Kd
x1 = np.zeros(kd.size,dtype='cfloat')
x2 = np.zeros(kd.size,dtype='cfloat')
x3 = np.zeros(kd.size,dtype='cfloat')

for i in range(kd.size):
    a = 1.
    b = -2.
    c = 1+gamma*kp*thau**2+gamma*thau*i
    d = -gamma*thau*i
    x1[i],x2[i],x3[i] = Cardano(a,b,c,d)


# In[52]:

#Plot de las distintas raizes para cada Kd, x1,x2,x3
plt.figure(0)

plt.subplot(311)
plt.title('x1')
plt.plot(x1.real,x1.imag, 'bo')
#plt.xlim(-0.025,0.25)
plt.grid(True)

plt.subplot(312)
plt.title('x2')
plt.plot(x2.real,x2.imag, 'ro')
plt.ylabel('Imaginary')
#plt.xlim(0.994,1.0005)
#plt.ylim(-0.016,0.017)
plt.grid(True)

plt.subplot(313)
plt.title('x3')
plt.plot(x3.real,x3.imag, 'go')
plt.xlabel('Real')
#plt.xlim(0.75,1.01)
#plt.ylim(-0.016,0.016)
plt.grid(True)

plt.show()


# In[54]:

x1


# In[55]:

x2


# In[56]:

x3


# In[57]:

#Preparamos y definimos la secuencia y los plots

get_ipython().magic(u'matplotlib notebook')

import matplotlib.pyplot as plt
import numpy as np

from __future__ import print_function
from ipywidgets import interact, interactive, fixed
import ipywidgets as widgets

def plotSequence(y):
    n = np.linspace(0, y.size, y.size)
    plt.scatter(n, y)
    plt.plot([n, n], [np.zeros(n.size), y], color='gray', linestyle="--")
    return

def circlePlot(a):
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1)
    circ = plt.Circle((0, 0), 1, color='g', fill=False)
    ax.add_patch(circ)
    plt.show()
    for x in range(len(a)):
        plt.plot([0,a[x].real],[0,a[x].imag],'r-',label='python')
        plt.plot(a[x].real,a[x].imag,'ro',label='python')
    limit=np.max(np.ceil(np.absolute(a))) # set limits for axis
    plt.xlim((-limit,limit))
    plt.ylim((-limit,limit))
    plt.axis('equal')
    plt.ylabel('Imaginary')
    plt.xlabel('Real')
    plt.show()


# In[60]:

#Ploteamos las raízes (reales e imaginarias) en función de kd establecida
for i in range(kd.size):
    lmbda = np.array([x1.real[i]+x1.imag[i]*1j,x2.real[i]+x2.imag[i]*1j ,x3.real[i]+x3.imag[i]*1j])
    circlePlot(lmbda)


# Consider initial conditions as:
# 
# $\begin{matrix} 
# y[0] = 1\\
# y[1] = 1\\
# y[2] = 1\\
# \end{matrix}$
# 

# In[115]:

##Creación y plot de la respuesta del sistema
#Me falta comprender los argumentos de la función, 
#entiendo que parece un linspace de valor a valor con i incrementos pero no entiendo el porque. (Consultar en clase)

def PDControlThirdOrderCar(kp,kd):
    # start with the robot at 3m from beacon
    d0 = 1
    d1 = 1
    d2 = 1
    
    # This initializes the sequence
    n = 2000
    d = np.zeros(n)
    for i in range(n):
        if i == 0:  # first initial condition
            d[i] = d0
        if i == 1: # second initial condition
            d[i] = d1
        if i == 2: # third initial condition
            d[i] = d2
        if i > 2:
            d[i] = 2*d[i-1] - (1 + (gamma*kp*thau**2) + (thau*kd*gamma))*d[i-2] + (thau*kd*gamma)*d[i-3]
    # Plot the sequence
    plt.figure()
    plotSequence(d)
    return 

interact(PDControlThirdOrderCar, kp=(0,100,0.1), kd=(0,3,0.1))


# In[ ]:



