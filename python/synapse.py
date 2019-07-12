import numpy as np
from scipy.integrate import odeint
import pylab as plt

def alpha_n(vr,v):
    return (0.1-0.01*(v-vr))/(np.exp(1.0-0.1*(v-vr))-1.0)
def alpha_m(vr,v):
    return (2.5-0.1*(v-vr))/(np.exp(2.5-0.1*(v-vr))-1.0)
def alpha_h(vr,v):
    return 0.07 * np.exp((vr-v)/20.0)
def beta_n(vr,v):
    return 0.125*np.exp((vr-v)/80.0)
def beta_m(vr,v):
    return  4.0 * np.exp((vr-v)/18.0)
def beta_h(vr,v):
    return 1.0 / ( 1.0 + np.exp( 3.0 -0.1*(v-vr) ))

def odesys(x, t, C1,C2, Iext1, Iext2, gNa, gK, gL, vNa, vK, vL, vr, vth, g21):
    v1,n1,m1,h1, v2,n2,m2,h2 , s21 = x
    dv1dt= (Iext1 - gNa *m1*m1*m1*h1*(v1-vNa) - gK *n1*n1*n1*n1*(v1-vK) - gL*(v1-vL))/C1
    dn1dt= alpha_n(vr,v1)*(1.0-n1) -beta_n(vr,v1)*n1
    dm1dt= alpha_m(vr,v1)*(1.0-m1) -beta_m(vr,v1)*m1
    dh1dt= alpha_h(vr,v1)*(1.0-h1) -beta_h(vr,v1)*h1

    ds21dt= (.5*(1.0+np.tanh(5.0*(v1-vth)))*(1.0-s21)-.2*s21)
    dv2dt= (Iext2 -g21*s21*(v2-20) - gNa *m2*m2*m2*h2*(v2-vNa) - gK *n2*n2*n2*n2*(v2-vK) - gL*(v2-vL) )/C2
    dn2dt= alpha_n(vr,v2)*(1.0-n2) -beta_n(vr,v2)*n2
    dm2dt= alpha_m(vr,v2)*(1.0-m2) -beta_m(vr,v2)*m2
    dh2dt= alpha_h(vr,v2)*(1.0-h2) -beta_h(vr,v2)*h2

    dxdt = [dv1dt, dn1dt , dm1dt , dh1dt, dv2dt, dn2dt , dm2dt , dh2dt, ds21dt]
    return dxdt

vr = -65.0
C1=1.0
C2=1.0
Iext1=8.0
Iext2=8.0
g21=-.00
gNa=120.0
gK= 36.0
gL= 0.3
vNa= 50.0
vK= -77.0 
vL= -54.4
v0 = -65.0
n0 = 0.31768111085600037
m0 = 5.2934208791209428E-002
h0 = 0.59610933633356666  
x0 = [v0, n0, m0, h0, v0, n0, m0, h0,0.0]
vth=-40.0
t = np.linspace(0, 100, 1001)
sol = odeint(odesys, x0, t, args=(C1,C2, Iext1, Iext2, gNa, gK, gL, vNa, vK, vL, vr, vth, g21))

# print(sol[-1,:])

plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(t, sol[:, 0], 'black', label='v1(t)')
plt.plot(t, sol[:, 4], 'gray', label='v2(t)')
plt.legend(loc='upper right')
plt.xlabel('t')
plt.grid()
plt.subplot(2, 1, 2)
plt.plot(t, sol[:, 1], 'r', label='n1(t)')
plt.plot(t, sol[:, 2], 'g', label='m1(t)')
plt.plot(t, sol[:, 3], 'b', label='h1(t)')
plt.legend(loc='upper right')
plt.xlabel('t')
plt.grid()
# plt.savefig('order_param_r.eps', format='eps', dpi=300)
plt.savefig('2neuron.png', format='png', dpi=300)
plt.show()