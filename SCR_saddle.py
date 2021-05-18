import numpy as np

## saddle design parameters
A = 320
b = 24
H = 0
L = 1491
P = 96.7
Q = 116451/2
R = 18.5
ts = 0.298
S = 12984
y = 50038
E = 0.85

K1 = 0.335
K2 = 1.171
K3 = 0.319
K4 = 0.88
K5 = 0.401
K6 = 0.013
K7 = 0.76
K8 = 0.603
Kstiff = 3.14

## stress at saddles - K1 (tension); K2(compression); Kstiff(exceed S)
s1a = Q*A/(Kstiff*(R**2)*ts)
slb = 1-A/L+(((R**2)-(H**2))/(2*A*L))
slc = 1+((4*H)/(3*L))
sl_saddles = s1a*(1-(slb/slc))

## stress at midspan
s1d = Q*L/(4*np.pi*(R**2)*ts)
sle = 1+2*(((R**2)-(H**2))/L**2)
sld = 1 + (4*H)/(3*L)
slf = 4*A/L
sl_midspan = s1d*(sle/sld-slf)

## stress due to internal pressure
int_pressure = (P*R)/(2*ts)

## tensional stress
tension = sl_midspan + int_pressure

## tangential stress
tang = (K3*Q)/(R*ts)*((L-2*A)/(L+(4*H)/4))

## circumferential stress
circ = -(K7*Q)/(ts*(b+1.56*((R*ts)**(0.5))))

print("at saddles: ", sl_saddles)
print("at midspan: ", sl_midspan)
print("internal pressure stress: ", int_pressure)
print("tension sum: ", tension, " must below lower than: ", S)
print("tangential stress: ", tang)
print("circumferential stress ", circ)
