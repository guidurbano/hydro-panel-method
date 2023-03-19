# This file is a part of the Hidrodinâmica II – PNV3413 project (Pro. Alexandre Simos)
# Autor: Guilherme Donatoni Urbano [guilherme.urbano@usp.br] NUSP: 9835985
# All Rights Reserved

import numpy as np
import math as m
import matplotlib.pyplot as plt
import cmath as cm

ro = 1025
viscosidade = 1.2 * 10 ** -6
g = 9.81

Dc = 17.5
Dp = 12.5
Lp = 50
b = 52.5
d = 27.5
KG = 25
Kxx = 30
Kyy = 35
Kzz = 35

B_heave = 2.25 * 10 ** 6
B_roll = 9.13 * 10 ** 8
B_pitch = 1.82 * 10 ** 9

M = ro * (2 * m.pi * Dp ** 2 * Lp / 4 + 4 * m.pi * Dc ** 2 * d / 4)
print(f"M: {M}")
V = (2 * m.pi * Dp ** 2 * Lp / 4 + 4 * m.pi * Dc ** 2 * d / 4)

Vc = m.pi * Dc ** 2 / 4 * d
Vp = m.pi * Dp ** 2 / 4 * Lp

Ac = m.pi * Dc ** 2 / 4
Ap = m.pi * Dp ** 2 / 4

Awl = 4 * m.pi * Dc ** 2 / 4

print(Vp)
# item a) #############################################################
# surge
a11 = 4 * ro * Vc

# sway
a22 = 2 * ro * Vp + 4 * ro * Vc

# heave
print("HEAVE")

a33 = 2 * ro * m.pi * Dp ** 2 / 4 * Lp
print(f"m33: {a33}")

c33 = ro * g * Awl
w3 = m.sqrt(c33 / (M + a33))

T3 = 2 * m.pi / w3
print(f"T3: {T3}")

zeta = B_heave / (2 * (M + a33) * w3)

T3_amort = T3 / m.sqrt(1 - zeta ** 2)
print(f"T3_amort: {T3_amort}")
print("")

# roll
print("ROLL")

I44 = M * Kxx ** 2

a44 = 2 * (Vp * ro * ((b / 2) ** 2 + (KG - Dp / 2) ** 2)) + 4 * Ac * ro * d ** 3 / 3
print(f"a44: {a44}")

KB = (4 * Vc * d / 2 + 2 * Vp * Dp / 2) / V
print(f"KB: {KB}")

Syy = 4 * (m.pi * Dc ** 4 / 64 + Ac * (b / 2) ** 2)
print(f"Syy: {Syy}")

GMt = KB + (Syy / V) - KG
print(f"GMt: {GMt}")

c44 = ro * g * V * GMt

w4 = m.sqrt(c44 / (I44 + a44))

T4 = 2 * m.pi / w4
print(f"T4: {T4}")

zeta = B_roll / (2 * (I44 + a44) * w3)

T4_amort = T4 / m.sqrt(1 - zeta ** 2)
print(f"T4_amort: {T4_amort}")
print("")

# pitch
print("")
print("PITCH")

I55 = M * Kyy ** 2

a55 = 2 * (ro * Ap * Lp ** 3 / 3) + 4 * (ro * Ac * d ** 3 / 3)  # subtrair KG das componentes??
print(f"a55: {a55}")

Sxx = 4 * (m.pi * Dc ** 4 / 64 + Ac * ((Lp + Dc) / 2) ** 2)
print(f"Sxx: {Sxx}")

GMl = KB + (Sxx / V) - KG
print(f"GMl: {GMl}")

c55 = ro * g * V * GMl

w5 = m.sqrt(c55 / (I55 + a55))

T5 = 2 * m.pi / w5
print(f"T5: {T5}")

zeta = B_pitch / (2 * (I55 + a55) * w5)

T5_amort = T5 / m.sqrt(1 - zeta ** 2)
print(f"T5_amort: {T5_amort}")
print("")

I66 = M * Kzz ** 2

a66 = 4 * (Vp * ro * ((b / 2) ** 2 + (Lp / 2 + Dc / 2) ** 2)) + 2 * Ap * ro * Lp ** 3 / 3

"""ka=0,35 (coluna) e ka=0,25 (pontoon)
Coluna: KC(0)=0,35 e KC(-d)=0,12
Pontoon: KC(0)=0,50 e KC(-d)=0,16
"""
""""""
# item b) #############################################################
C_M = 2
C_D = 1
A = 1
T = 7

w = 2 * m.pi / T
k = w ** 2 / g

y_c = (b / 2, b / 2, -b / 2, -b / 2)
y_p = (b / 2, -b / 2)
z_p = -d + Dp / 2

tempo = np.linspace(0, 20, 100)
lista_F_inercial = []
lista_F_viscoso = []
lista_F_total = []

for t in tempo:  # Para cada instante de tempo
    F_inercial = 0
    F_viscoso = 0
    for y in y_c:  # Para cada coluna
        F_inercial += ro * Ac * C_M * g * A * (1 - m.exp(-k * d)) * m.cos(k * y - w * t)
        F_viscoso += -(1 / 4) * ro * Dc * C_D * g * A ** 2 * (1 - m.exp(-2 * k * d)) * m.sin(k * y - w * t) * abs(
            m.sin(k * y - w * t))
    for y in y_p:
        F_inercial += Lp * ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) * m.cos(k * y - w * t)
        F_viscoso += -0.5 * Lp * ro * Dp * C_D * w ** 2 * A ** 2 * (m.exp(2 * k * z_p)) * m.sin(k * y - w * t) * abs(
            m.sin(k * y - w * t))
    lista_F_inercial.append(F_inercial)
    lista_F_viscoso.append(F_viscoso)
    lista_F_total.append(F_inercial + F_viscoso)

"""

plt.title(f'Força de Excitação (Sway)')
plt.ylabel('Força (N)')
plt.xlabel('Tempo (s)')
plt.grid(True)
plt.plot(tempo, lista_F_inercial, color='red', label="Inercial")
plt.plot(tempo, lista_F_viscoso, color='blue', label="Viscoso")
plt.plot(tempo, lista_F_total, color='black', label="Total")
plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
plt.show()
"""
# item c) #############################################################
print("item c)\n")
print(f"Força de inércia máxima: {max(lista_F_inercial)}")
print(f"Força de viscosa máxima: {max(lista_F_viscoso)}")

print(f"Arrasto / Inercia = {max(lista_F_viscoso) / max(lista_F_inercial)}")

# item d) #############################################################
"""
tempo = np.linspace(0, 20, 100)
Periodo = np.linspace(4, 30, 27)  # Periodo = np.arange(4, 30.5, 0.5)

y_c = (b/2, b/2, -b/2, -b/2)
y_p = (b/2, -b/2)
z_p = -d + Dp/2

F2 = []
F3 = []
F4 = []

RAO2 = []
RAO3 = []
RAO4 = []
w_lista=[]

print("item d)")
for T in Periodo: # Para cada periodo
    w = 2 * m.pi / T
    w_lista.append(w)
    k = w ** 2 / g
    lista_F_inercial2 = []
    lista_F_inercial3 = []
    lista_F_inercial4 = []
    for t in tempo:  # Para cada instante de tempo
        F_inercial2 = 0
        F_inercial3 = 0
        F_inercial4 = 0
        for y in y_c:  # Para cada coluna
            F_inercial2 += ro * Ac * C_M * g * A * (1 - m.exp(-k * d)) * m.cos(k * y - w * t)
            F_inercial3 += - ro * g * A * Ac * m.exp(-k * d) * m.sin(k * y - w * t)
            F_inercial4 += ((-(k * KG + 1) * m.exp(-d * k) - d * k + KG * k + 1) / k**2) *ro * Ac * C_M * w**2 * A * m.cos(k * y - w * t) + (
                - y * ro * g * Ac * m.exp(-k * d) * A * m.sin(k * y - w * t))
        for y in y_p:
            F_inercial2 += Lp * ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) * m.cos(k * y - w * t)
            F_inercial3 += Lp * ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) * m.sin(k * y - w * t)
            F_inercial4 += (z_p - d + KG) * Lp * ro * Ap * C_M * w**2 * A * m.exp(k * z_p) * m.cos(k * y - w * t) + (
               y * Lp * ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) * m.sin(k * y - w * t))

        lista_F_inercial2.append(F_inercial2)
        lista_F_inercial3.append(F_inercial3)
        lista_F_inercial4.append(F_inercial4)

    F2.append(max(lista_F_inercial2))
    F3.append(max(lista_F_inercial3))
    F4.append(max(lista_F_inercial4))
    RAO2.append(max(lista_F_inercial2) / (w**2 * (M + a22)))
    RAO3.append(abs(max(lista_F_inercial3) / (-w**2 * (M + a33) + w*B_heave*1j + c33)))
    RAO4.append(abs(max(lista_F_inercial4) / (-w**2 * (I44 + a44) + w*B_roll*1j + c44))*180/m.pi)

    #if T in (5, 10, 15, 20, 30):
       # print(f"T = {T}: RAO = {max(lista_F_inercial2) / (w**2 * (M + a22))}")

print(F3)
"""
"""
plt.title(f'Amplitude de Força (Sway)')
plt.ylabel('Força (N)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, F2, color='red')
plt.show()

plt.title(f'RAO (Sway)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO2, color='blue')
plt.show()
"""
"""
plt.title(f'Amplitude de Força (Heave)')
plt.ylabel('Força (N)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, F3, color='red')
plt.show()

plt.title(f'RAO (Heave)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO3, color='blue')
plt.show()


plt.title(f'Amplitude de Momento (Roll)') # Pra roll fiz Periodo 80s
plt.ylabel('Momento (N.m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, F4, color='red')
plt.show()

plt.title(f'RAO (Roll)') 
plt.ylabel('Amplitude (°/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO4, color='blue')
plt.show()
"""

# item f) #############################################################

delta_T = 0.1

tempo = np.linspace(0, 10, 100)
Periodo = np.arange(4, 60, delta_T)

# Pra achar os valores tabelados
Excel = np.arange(4, 20.5, 0.5)
Excel = np.append(Excel, np.arange(21, 30.5, 1))

# Coordenadas
x_c = ((Lp + Dc) / 2, -(Lp + Dc) / 2)
y_c = (b / 2, b / 2, -b / 2, -b / 2)
y_p = (b / 2, -b / 2)
z_p = -d + Dp / 2

# Angulo de incidencia (beta)
s = m.sin(m.pi / 4)

# JONSWAP
Tp = 14
wp = 2 * m.pi / Tp
Hs = 6.5
gama = 3.3
jonswap = []
jonswap_heave = []
m0_heave = 0
m2_heave = 0

# Ponto P
alfa = m.atan(10 / 35)
Vetor_P = m.sqrt(35 ** 2 + 10 ** 2)
H = 10

# RAO
RAO1 = []
RAO2 = []
RAO3 = []
RAO4 = []
RAO4l = []  # em relacao a Lagua
RAO5 = []
RAO5l = []  # em relacao a Lagua
RAO6 = []

RAO3f = []
RAO4f = []
RAO5f = []

lista_k = []
lista_w = []

print("item e)")
for T in Periodo:  # Para cada periodo
    w = 2 * m.pi / T
    k = w ** 2 / g
    lista_k.append(k)
    lista_w.append(w)
    if w < wp:
        sigma = 0.07
    else:
        sigma = 0.09
    A_jon = m.exp(- (w / wp - 1) ** 2 / (2 * sigma ** 2))

    jon = 320 * Hs ** 2 / Tp ** 4 / w ** 5 * m.exp(-1950 / Tp ** 4 / w ** 4) * gama ** A_jon

    lista_F_inercial1 = []
    lista_F_inercial2 = []
    lista_F_inercial3 = []
    lista_F_inercial4 = []
    lista_F_inercial4l = []
    lista_F_inercial5 = []
    lista_F_inercial5l = []
    lista_F_inercial6 = []

    for t in tempo:  # Para cada instante de tempo
        F_inercial1 = 0
        F_inercial2 = 0
        F_inercial3 = 0
        F_inercial4 = 0
        F_inercial4l = 0
        F_inercial5 = 0
        F_inercial5l = 0
        F_inercial6 = 0

        for x in ((Lp + Dc) / 2, -(Lp + Dc) / 2):  # para cada coluna
            for y in (b / 2, -b / 2):
                F_inercial1 += s * ro * Ac * C_M * g * A * (1 - m.exp(-k * d)) * m.cos(s * k * (x + y) - w * t)
                F_inercial2 += s * ro * Ac * C_M * g * A * (1 - m.exp(-k * d)) * m.cos(s * k * (x + y) - w * t)
                F_inercial3 += - ro * g * A * Ac * m.exp(-k * d) * m.sin(s * k * (x + y) - w * t)
                F_inercial4 += - s * ro * Ac * C_M * A * w ** 2 * m.cos(s * k * (x + y) - w * t) * (
                        (k * KG + 1) * m.exp(-k * d) + k * d - k * KG - 1) / k ** 2 - y * ro * g * Ac * A * m.exp(
                    -k * d) * m.sin(s * k * (x + y) - w * t)
                F_inercial5 += - s * ro * Ac * C_M * A * w ** 2 * m.cos(s * k * (x + y) - w * t) * (
                        (k * KG + 1) * m.exp(-k * d) + k * d - k * KG - 1) / k ** 2 - x * ro * g * Ac * A * m.exp(
                    -k * d) * m.sin(s * k * (x + y) - w * t)
                F_inercial6 += s * y * ro * Ac * C_M * g * A * (1 - m.exp(-k * d)) * m.cos(
                    s * k * (x + y) - w * t) + s * x * ro * Ac * C_M * g * A * (1 - m.exp(-k * d)) * m.cos(
                    s * k * (x + y) - w * t)

                # Em relacao a origem

                F_inercial4l += - s * ro * Ac * C_M * A * g * m.cos(s * k * (x + y) - w * t) * (
                        d * m.exp(-k * d) + (m.exp(-k * d) - 1) / k) - y * ro * g * Ac * A * m.exp(
                    -k * d) * m.sin(s * k * (x + y) - w * t)
                F_inercial5l += - s * ro * Ac * C_M * A * g * m.cos(s * k * (x + y) - w * t) * (
                        d * m.exp(-k * d) + (m.exp(-k * d) - 1) / k) - x * ro * g * Ac * A * m.exp(
                    -k * d) * m.sin(s * k * (x + y) - w * t)

        for y in (b / 2, -b / 2):  # Para cada pontoon

            F_inercial2 += ro * Ap * C_M * g * A * m.exp(k * z_p) * (
                    m.sin(s * k * (y + Lp / 2) - w * t) - m.sin(s * k * (y - Lp / 2) - w * t))

            F_inercial3 += 1 / s * ro * Ap * C_M * g * A * m.exp(k * z_p) * (
                    -m.cos(s * k * (y + Lp / 2) - w * t) + m.cos(s * k * (y - Lp / 2) - w * t))

            F_inercial4 += (z_p + d - KG) * ro * Ap * C_M * g * A * m.exp(k * z_p) * (
                    m.sin(s * k * (y + Lp / 2) - w * t) - m.sin(
                s * k * (y - Lp / 2) - w * t)) + 1 / s * y * ro * Ap * C_M * g * A * m.exp(k * z_p) * (
                                   -m.cos(s * k * (y + Lp / 2) - w * t) + m.cos(s * k * (y - Lp / 2) - w * t))
            F_inercial5 += ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) / k ** 2 * (
                    -s * k * Lp * m.cos(s * k * (y + Lp / 2) - w * t) - s * k * Lp * m.cos(s * k * (y - Lp / 2) - w * t)
                    + 2 * m.sin(s * k * (y + Lp / 2) - w * t) - 2 * m.sin(s * k * (y - Lp / 2) - w * t))

            F_inercial6 += s * ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) / k ** 2 * (
                    2 * m.cos(s * k * (y + Lp / 2) - w * t) - 2 * m.cos(s * k * (y - Lp / 2) - w * t)
                    + s * k * Lp * m.sin(s * k * (y + Lp / 2) - w * t) + s * k * Lp * m.sin(
                s * k * (y - Lp / 2) - w * t))

            # Em relacao a origem
            F_inercial4l += z_p * ro * Ap * C_M * g * A * m.exp(k * z_p) * (
                    m.sin(s * k * (y + Lp / 2) - w * t) - m.sin(
                s * k * (y - Lp / 2) - w * t)) + 1 / s * y * ro * Ap * C_M * g * A * m.exp(k * z_p) * (
                                    -m.cos(s * k * (y + Lp / 2) - w * t) + m.cos(s * k * (y - Lp / 2) - w * t))
            F_inercial5l += ro * Ap * C_M * w ** 2 * A * m.exp(k * z_p) / k ** 2 * (
                    -s * k * Lp * m.cos(s * k * (y + Lp / 2) - w * t) - s * k * Lp * m.cos(s * k * (y - Lp / 2) - w * t)
                    + 2 * m.sin(s * k * (y + Lp / 2) - w * t) - 2 * m.sin(s * k * (y - Lp / 2) - w * t))

        lista_F_inercial1.append(F_inercial1)
        lista_F_inercial2.append(F_inercial2)
        lista_F_inercial3.append(F_inercial3)
        lista_F_inercial4.append(F_inercial4)
        lista_F_inercial4l.append(F_inercial4l)
        lista_F_inercial5.append(F_inercial5)
        lista_F_inercial5l.append(F_inercial5l)
        lista_F_inercial6.append(F_inercial6)

    RAO1_ = max(lista_F_inercial1) / (w ** 2 * (M + a11))

    RAO2_ = max(lista_F_inercial2) / (w ** 2 * (M + a22))

    RAO3_m = abs(max(lista_F_inercial3) / (-w ** 2 * (M + a33) + w * B_heave * 1j + c33))
    RAO3_f = cm.phase((max(lista_F_inercial3) / (-w ** 2 * (M + a33) + w * B_heave * 1j + c33)))

    RAO4_ = abs(max(lista_F_inercial4) / (-w ** 2 * (I44 + a44) + w * B_roll * 1j + c44)) * 180 / m.pi
    RAO4l_m = abs(max(lista_F_inercial4l) / (-w ** 2 * (I44 + a44) + w * B_roll * 1j + c44))
    RAO4l_f = cm.phase((max(lista_F_inercial4l) / (-w ** 2 * (I44 + a44) + w * B_roll * 1j + c44)))

    RAO5_ = abs(max(lista_F_inercial5) / (-w ** 2 * (I55 + a55) + w * B_pitch * 1j + c55)) * 180 / m.pi
    RAO5l_m = abs(max(lista_F_inercial5l) / (-w ** 2 * (I55 + a55) + w * B_pitch * 1j + c55))
    RAO5l_f = cm.phase((max(lista_F_inercial5l) / (-w ** 2 * (I55 + a55) + w * B_pitch * 1j + c55)))

    RAO6_ = abs(max(lista_F_inercial6) / (-w ** 2 * (I66 + a66))) * 180 / m.pi

    jonswap.append(jon)

    jon_heave = jon * RAO3_m ** 2
    jonswap_heave.append(jon_heave)

    #if T in Excel:
     #   print(round(RAO6_, 6))

    RAO1.append(RAO1_)
    RAO2.append(RAO2_)
    RAO3.append(RAO3_m)
    RAO4.append(RAO4_)
    RAO4l.append(RAO4l_m)
    RAO5.append(RAO5_)
    RAO5l.append(RAO5l_m)
    RAO6.append(RAO6_)

    RAO3f.append(RAO3_f)
    RAO4f.append(RAO4l_f)
    RAO5f.append(RAO5l_f)


lista_w.append(w)
for o in range(0, len(Periodo)):
    m0_heave += jonswap_heave[o] * abs(lista_w[o+1] - lista_w[o])
    m2_heave += jonswap_heave[o] * lista_w[o] ** 2 * abs(lista_w[o+1] - lista_w[o])

print(f"Amplitude significativa: {2 * m.sqrt(m0_heave)}")
print(f"Período ascendente: {2 * m.pi * m.sqrt(m0_heave / m2_heave)}")

RAOP = []
RAOP2 = []
jonswap_P = []
m0_P = 0

for i in range(0, len(Periodo)):
    RAOP_ = []
    RAOP2_ = []
    for j in range(0, len(tempo)):
        rao3P = (RAO3[i] * m.sin(RAO3f[i] - w * j))
        rao4P = (m.sin(RAO4l[i] + alfa) * Vetor_P - H) * m.sin(RAO4f[i] + m.pi / 2 - w * j)
        rao5P = (m.sin(RAO5l[i] + alfa) * Vetor_P - H) * m.sin(RAO5f[i] + m.pi / 2 - w * j)
        Elevacao = -A*m.sin(s*lista_k[i]*(35+35) - tempo[j] * (lista_k[i]*g)**0.5)

        RAOP_.append(rao3P + rao4P + rao5P)
        RAOP2_.append(rao3P + rao4P + rao5P + Elevacao)

    a = max(RAOP_)
    b = min(RAOP2_)
    print(round(b, 6))

    jon_P = jonswap[i] * b ** 2

    m0_P += jon_P * abs(lista_w[i+1] - lista_w[i])

    RAOP.append(a)
    RAOP2.append(b)
    jonswap_P.append(jon_P)

print(f"Probabilidade de Slamming: {(m.exp(-H**2/(2*2*m0_P**0.5)))*100}")

"""
plt.title(f'RAO (Surge)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO1, color='blue')
plt.show()
"""
"""
plt.title(f'RAO (Sway)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO2, color='blue')
plt.show()
"""
"""
plt.title(f'RAO (Heave)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO3_m, color='blue')
plt.show()
"""
"""
plt.title(f'RAO (Roll)') 
plt.ylabel('Amplitude (°/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO4, color='blue')
plt.show()
"""
"""
plt.title(f'RAO (Pitch)')
plt.ylabel('Amplitude (°/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO5, color='blue')
plt.show()
"""
"""
plt.title(f'RAO (Yaw)')
plt.ylabel('Amplitude (°/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO6, color='blue')
plt.show()
"""
"""
plt.title(f'RAO - Origem (Roll)')
plt.ylabel('Amplitude (°/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO4l_m, color='blue')
plt.show()
"""
"""
plt.title(f'RAO - Origem (Pitch)')
plt.ylabel('Amplitude (°/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAO5l_m, color='blue')
plt.show()
"""
"""
plt.title(f'JONSWAP')
plt.ylabel('Espectro de energia (m².s)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, jonswap, color='orange')
plt.show()
"""
"""
plt.title(f'Espectro de heave')
plt.ylabel('Espectro de energia (m².s)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, jonswap_heave, color='orange')
plt.show()
"""
"""
plt.title(f'RAO do ponto P (Vertical)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAOP, color='blue')
plt.show()
"""

plt.title(f'RAO Relativo do ponto P (Vertical)')
plt.ylabel('Amplitude (m/m)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, RAOP2, color='blue')
plt.show()
"""
plt.title(f'Espectro de movimento relativo de P')
plt.ylabel('Espectro de energia (m².s)')
plt.xlabel('Período (s)')
plt.grid(True)
plt.plot(Periodo, jonswap_P, color='orange')
plt.show()
"""
"""
plt.title(f'Séries temporais das Forças (T = 10s)')
plt.ylabel('Amplitude da Força (N ou Nm)')
plt.xlabel('Tempo (s)')
plt.grid(True)
plt.plot(tempo, lista_F_inercial3, label="Heave")
plt.plot(tempo, lista_F_inercial4l, label="Roll")
plt.plot(tempo, lista_F_inercial5l, label="Pitch")
plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
plt.show()
"""