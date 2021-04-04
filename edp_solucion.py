# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

### DATOS
orden=10
L=50e-6
C0=1000
t0=0.36
T=300   
A=0.0002 
D=2e-11
F= 96485
I0=30e-3 
R=8.314
k0=2.1e-4
k=3.8e-2
Kd=2*k*R*T*(1-t0)/F
### para el 14
colour='b'
def concentracion4(x,t,n,L):
    ### Funcion de la concentracion del ejercicio 4
    c=C0+(1-t0)*I0*(x-L/2.)/(F*A*D)
    suma=0
    for i in range(n):
        suma+=(np.cos((2*i+1)*np.pi*x/L)/(2*i+1)**2)*np.exp(-(2*i+1)**2*D*np.pi**2*t/L**2)
    c+=suma*4*(1-t0)*I0*L/(F*A*D*np.pi**2)
    return c

def concentracion12(x,t,n,w):
    ### Funcion de la concentracion del ejercicio 12
    c=C0+(1-t0)*I0*(L/2-x*np.cos(w*t))/(A*F*D)+(1-t0)*I0*L*(np.cos(w*t)-1)/(A*F*D*2)
    suma1=0
    suma2=0
    for i in range(n):
        suma1+=(np.cos((2*i+1)*np.pi*x/L)/(2*i+1)**2)*np.exp(-(2*i+1)**2*D*np.pi**2*t/L**2)
        suma2+=w*np.cos((2*i+1)*np.pi*x/L)*(D*np.sin(w*t)*((2*i+1)*np.pi/L)**2-w*(np.cos(w*t)-np.exp(-(2*i+1)**2*D*np.pi**2*t/L**2)))/(w**2+(D**2)*((2*i+1)*np.pi/L)**4)
    c+=4*(1-t0)*I0*L*(suma2-suma1)/(A*F*D*np.pi**2)
    return c

def deltav(C_0,C_L):
    num=np.ones(len(t))+np.sqrt(np.ones(len(t))+4*C_L*(F*k0)**2/(C0*(I0/A)**2))
    den=np.ones(len(t))+np.sqrt(np.ones(len(t))+4*C_0*(F*k0)**2/(C0*(I0/A)**2))
    DeltaV=2.*R*T*np.log(num/den)/F+Kd*np.log(C_0/C_L)/k+I0*L/(A*k)
    return DeltaV

def C(x,t,terms):
    aux = (1-t0)*I0/(D*F*A)
    sumatory = 0
    for m in range(terms):
        npi = (2*m+1)*np.pi
        sumatory += (2/npi)**2*aux*L*np.cos(npi*x/L)*np.exp(-(npi/L)**2*D*t)
    C = C0 + aux*(x-L/2) + sumatory
    return C

# -- Gráficas C(x,t), para un array dado de valores de L -----------------------

def plot_Cxt(L1):

    t_max = 20
    t1 = np.linspace(0,int(t_max),100000)

    x1 = np.linspace(0,L1,10000)
    lj = round(L1*10**6)/1000.
    print(f"L = {lj} mm")

    plt.figure(1)
    plt.clf()

    for i in range(len(x1)):
        if i%round(len(x1)/10.) == 0 or i in [0,len(x1)-1]:
            print(f"\t> i: {i}; x: {x1[i]}")
            x_i = round(x1[i]*1000)/1000.
            
            plt.title(f"L = {lj} mm")
            plt.ylabel(f"C(x,t)")
            plt.xlabel(f"t")
            plt.plot(t,C(x1[i],t, orden), label=f"C(X_i={i},t)")

    plt.legend(loc=4)
    plt.show()

# ------------------------------------------------------------------------------


def concentracion8_9(t):
    C_0=C0-(1-t0)*2*I0*np.sqrt(t)/(A*F*np.sqrt(np.pi*D))
    C_L=C0+(1-t0)*2*I0*np.sqrt(t)/(A*F*np.sqrt(np.pi*D))
    return C_0,C_L

def concentracion15(x,t,w):
    expo=(1-t0)*I0*np.exp(x*np.sqrt(w/(2*D)))*(np.sqrt(D/(2*w)))/(A*F*D)
    A1=np.cos(x*np.sqrt(w/(2*D)))*np.cos(w*t)
    A2=np.cos(x*np.sqrt(w/(2*D)))*np.sin(w*t)
    A3=np.sin(x*np.sqrt(w/(2*D)))*np.sin(w*t)
    A4=np.sin(x*np.sqrt(w/(2*D)))*np.cos(w*t)
    sol=expo*(A1+A2+A3-A4)
    return sol

def Zb_15(w,orden):
    Zr=(1-t0)*L*Kd/(A*F*D*k*C0)-L/(k*A)
    Zi=0
    for i in range(orden):
        Zr+=-8*(1-t0)*L*w**2*Kd/((C0*k*A*F*D*(2*i+1)**2*np.pi**2)*(((2*i+1)*np.pi/L)**4*D**2+w**2))
        Zi+=8*(1-t0)*w*Kd/(C0*k*L*A*F*((2*i+1)*np.pi/L)**4*D**2+w**2)
    return Zr,Zi

tau=np.pi*D*(((F*C0*A)/((1-t0)*I0))**2)/4
#print(tau)

entrada=float(input("Que ejercicio quieres representar? "))



t=np.linspace(0,30,10000)       ### mas o menos a partir de t=17, C_0 da numeros negativos
x_muestra = np.linspace(0,L,1000)

if entrada==6:
    C_01=concentracion4(0,t,orden,L/2.)
    C_L1=concentracion4(L/2.,t,orden,L/2.)
    DeltaV1=deltav(C_01,C_L1)
    C_02=concentracion4(0,t,orden,L)
    C_L2=concentracion4(L,t,orden,L)
    DeltaV2=deltav(C_02,C_L2)
    C_03=concentracion4(0,t,orden,2*L)
    C_L3=concentracion4(2*L,t,orden,2*L)
    DeltaV3=deltav(C_03,C_L3)
    C_04=concentracion4(0,t,orden,10*L)
    C_L4=concentracion4(10*L,t,orden,10*L)
    DeltaV4=deltav(C_04,C_L4)
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(t,DeltaV1,colour)
    axs[0, 0].grid(True, linestyle = '--')
    axs[0, 0].set_title('L= ' +str(L/2.))
    axs[0, 0].set_ylabel(r'$\Delta$V [V]')
    axs[0, 1].plot(t,DeltaV2, colour)
    axs[0, 1].set_title('L= ' +str(L))
    axs[0, 1].grid(True, linestyle = '--')
    axs[1, 0].plot(t,DeltaV3, colour)
    axs[1, 0].set_title('L= ' +str(2*L))
    axs[1, 0].set_ylabel(r'$\Delta$V [V]')
    axs[1, 0].grid(True, linestyle = '--')
    axs[1, 0].set_xlabel('Tiempo [s]')
    axs[1, 1].plot(t,DeltaV4, colour)
    axs[1, 1].set_title('L= ' +str(10*L))
    axs[1, 1].set_xlabel('Tiempo [s]')
    axs[1, 1].grid(True, linestyle = '--')
    plt.show()
elif entrada==7:
    plot_Cxt(L)
elif entrada==11:
    Cini,Cfin=concentracion8_9(t)
    #plt.plot(t,Cini)
    #plt.plot(t,Cfin)
    vectau=np.ones(10)*tau
    vecaux=np.linspace(-0.6,0.2,10)
    difpot=deltav(Cini,Cfin)
    plt.plot(t,difpot,colour)
    plt.plot(vectau,vecaux,'--r')
    plt.grid(True, linestyle = '--')
    plt.title('$\Delta$V vs t')
    plt.ylabel(r'$\Delta$V [V]')
    plt.xlabel('Tiempo [s]')
    plt.show()
elif entrada==14:
    I0=200e-6       ### Aqui cambia la intensidad          ### frecuencias a probar: 0.001 Hz, 1 Hz, 1000 Hz, 100000Hz
    punto=0
    C12_1=concentracion12(punto,t,orden,0.001)
    C12_2=concentracion12(punto,t,orden,1)
    C12_3=concentracion12(punto,t,orden,1000)
    C12_4=concentracion12(punto,t,orden,100000)
    ### hago esto para que podamos representar el que nos pete, aunque solo pide en 0
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(t,C12_1,colour)
    axs[0, 0].grid(True, linestyle = '--')
    axs[0, 0].set_title('w= ' +str(0.001) + 'Hz')
    axs[0, 0].set_ylabel('Concentracion [mol/m**2]')
    axs[0, 1].plot(t,C12_2,colour)
    axs[0, 1].set_title('w= ' +str(1)+ 'Hz')
    axs[0, 1].grid(True, linestyle = '--')
    axs[1, 0].plot(t,C12_3,colour)
    axs[1, 0].set_title('w= ' +str(100)+ 'Hz')
    axs[1, 0].set_ylabel('Concentracion [mol/m**2]')
    axs[1, 0].grid(True, linestyle = '--')
    axs[1, 0].set_xlabel('Tiempo [s]')
    axs[1, 1].plot(t,C12_4,colour)
    axs[1, 1].set_title('w= ' +str(100000)+ 'Hz')
    axs[1, 1].set_xlabel('Tiempo [s]')
    axs[1, 1].grid(True, linestyle = '--')
    plt.show()
elif entrada==16:
    w=1000
    punto=0
    C15_1=concentracion15(punto,t,0.001)
    C15_2=concentracion15(punto,t,1)
    C15_3=concentracion15(punto,t,1000)
    C15_4=concentracion15(punto,t,100000)
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(t,C15_1,colour)
    axs[0, 0].grid(True, linestyle = '--')
    axs[0, 0].set_title('w= ' +str(0.001) + 'Hz')
    axs[0, 0].set_ylabel('Concentracion [mol/m**2]')
    axs[0, 1].plot(t,C15_2,colour)
    axs[0, 1].set_title('w= ' +str(1)+ 'Hz')
    axs[0, 1].grid(True, linestyle = '--')
    axs[1, 0].plot(t,C15_3,colour)
    axs[1, 0].set_title('w= ' +str(1000)+ 'Hz')
    axs[1, 0].set_ylabel('Concentracion [mol/m**2]')
    axs[1, 0].grid(True, linestyle = '--')
    axs[1, 0].set_xlabel('Tiempo [s]')
    axs[1, 1].plot(t,C15_4,colour)
    axs[1, 1].set_title('w= ' +str(100000)+ 'Hz')
    axs[1, 1].set_xlabel('Tiempo [s]')
    axs[1, 1].grid(True, linestyle = '--')
    plt.show()
elif entrada==19:
    #frec=np.array([0.001,0.01,0.1,1,10,100,1000])
    frec=np.linspace(0.001,1000,1000000)
    Zr,Zi=Zb_15(frec,orden)
    print(Zr)
    print(Zi)
    plt.plot(Zr,Zi,colour)
    #plt.plot(Zr,Zi,'og')
    plt.xlabel('Re(Zb)')
    plt.ylabel('Im(Zb)')
    plt.title('Zb(w) [V/A]')
    plt.grid()
    plt.show()
else:
    print("Dicho ejercicio no esta resuelto aqui.")