# -*- coding: utf-8 -*-
# This file is a part of the Hidrodinâmica II – PNV3413 project.
# Método dos Painéis para escoamentos potenciais sem superfície livre 2D
# Autor: Guilherme Donatoni Urbano [guilherme.urbano@usp.br] NUSP: 9835985
# All Rights Reserved

import numpy as np
import math as m
import matplotlib.pyplot as plt
import painel

# ===================== SETUP ========================
EP = 2  # EP1 = 1;  EP2 m11 = 2; EP2 m22 = 3; EP3 = 4
# ====================================================

# ================================= PARAMETROS DO ESCOAMENTO =======================================
N = painel.N  # Numero de paineis
ang_ataque = 0  # ângulo de ataque em graus 0, 5, 10, 15
ro = 1025  # Densidade
pressao0 = 0  # Pressao ao longe
U_longe = {1: [1, 0],
           2: [1, 0],
           3: [0, 1],
           4: [1 * m.cos(m.radians(ang_ataque)), 1 * m.sin(m.radians(ang_ataque))]}  # Velocidade ao longe em [x, y]
# ==================================================================================================

def calcular_coord_ref_local(x, y, X_inicio, Y_inicio, beta):
    """

    :param x: coordenadaX no ref. global
    :param y: coordenadaY no ref. global
    :param X_inicio: coordenadaX no ref. local
    :param Y_inicio: coordenadaX no ref. local
    :param beta: angulo entre ref. local e ref. global
    :return: coordenadas do ponto P no referencial local
    """

    X = (x - X_inicio) * m.cos(beta) + (y - Y_inicio) * m.sin(beta)
    Y = - (x - X_inicio) * m.sin(beta) + (y - Y_inicio) * m.cos(beta)

    return X, Y


def calcular_vel_normal_induzida(X, Y, b_j, i, j, EP):
    """

    :param X: coordenadaX no ref. local
    :param Y: coordendaY no ref. local
    :param b_j: tamanho do painel
    :param i: numero do painel do ponto de colocacao
    :param j: numero do painel da fonte
    :return: velocidade normal induzida no ponto X e Y
    """


    if EP == 4:
        U = (1 / (4 * m.pi)) * m.log(((X - b_j) ** 2 + Y ** 2) / (
                X ** 2 + Y ** 2))
    else:
        U = (1 / (4 * m.pi)) * m.log((X ** 2 + Y ** 2) / (
                (X - b_j) ** 2 + Y ** 2))
    if i != j:
        V = (1 / (2 * m.pi)) * (m.atan(np.float64(X) / Y) - m.atan(np.float64(X - b_j) / Y))

    else:
        V = - 1 / 2

    return U, V


def calcular_vel_global(U, V, beta):
    """

    :param U: velocidade normal ind. em X
    :param V: velocidade normal ind. em Y
    :param beta: angulo entre ref. local e global
    :return: velocidade normal induzida no ponto X e Y (ref. global)
    """
    u = U * m.cos(beta) - V * m.sin(beta)
    v = U * m.sin(beta) + V * m.cos(beta)

    return u, v


def grafico_Cp(grid, matriz_Cp):
    """

    :param grid: malha com os paineis
    :param matriz_Cp: valores de Cp para pontos de colocacao
    :return: grafico Cp x teta
    """
    if EP == 2:
        matriz_teta = []

        for painel_i in grid.values():
            matriz_teta.append(painel_i.teta)
        matriz_teta.sort()

        plt.figure()
        plt.plot(matriz_teta, matriz_Cp, marker='o')
        plt.title('Cp x Teta (°) (a/b = 3.0)')
        plt.xlabel('Teta (°)')
        plt.ylabel('Cp')
        plt.grid(True)
        plt.show()
        plt.savefig("EP1 - Cp x Teta")

    if EP == 4:
        x = []

        for painel_i in grid.values():
            if painel_i.Y_inicio >= 0:
                x.append(painel_i.X_inicio)
        matriz_Cp = matriz_Cp[0:16]
        plt.figure()
        plt.plot(x, matriz_Cp, marker='o')
        plt.title('Cp x (x/c)')
        plt.xlabel('x/c')
        plt.ylabel('Cp')
        plt.grid(True)
        plt.savefig("EP4 - Cp x Posição")


def grafico_Vel(grid, matriz_Vel_total):
    """

    :param grid: malha com os paineis
    :param matriz_Vel_total:
    :return: grafico da velocidade induzida + ao longe
    """
    x = []
    y = []
    Dist_Vel = []
    x_coloq = []
    plt.figure()

    for painel_i in grid.values():
        if EP != 4:
            plt.arrow(painel_i.colocacaoX, painel_i.colocacaoY,
                      matriz_Vel_total[painel_i.numero - 1][0] + U_longe[EP][0],
                      matriz_Vel_total[painel_i.numero - 1][1] + U_longe[EP][1],
                      head_width=0.07, width=0.001, color='black', length_includes_head=True
                      )
        else:
            plt.arrow(painel_i.colocacaoX, painel_i.colocacaoY,
                      matriz_Vel_total[painel_i.numero - 1][0] + U_longe[EP][0],
                      matriz_Vel_total[painel_i.numero - 1][1] + U_longe[EP][1],
                      head_width=0.02, width=0.001, color='black', length_includes_head=True
                      )
            Dist_Vel.append((matriz_Vel_total[painel_i.numero - 1][2] + (m.sqrt(U_longe[EP][0] + U_longe[EP][1])) ** 2) ** 2)  # Seria dividido por 1

        x_coloq.append(painel_i.colocacaoX)

        x.append(painel_i.X_inicio)
        y.append(painel_i.Y_inicio)

    x.append(painel_i.X_fim)
    y.append(painel_i.Y_fim)

    plt.axis('scaled')
    if EP == 4:

        plt.subplot(2, 1, 1)
        plt.plot(x_coloq, Dist_Vel)
        plt.title('Velocidades')
        plt.xlim(0, 1)
        plt.ylim(0, 2)
        plt.ylabel("(v/U)²")

        plt.subplot(2, 1, 2)
        for painel_i in grid.values():
            plt.arrow(painel_i.colocacaoX, painel_i.colocacaoY,
                      matriz_Vel_total[painel_i.numero - 1][0] + U_longe[EP][0],
                      matriz_Vel_total[painel_i.numero - 1][1] + U_longe[EP][1],
                      head_width=0.02, width=0.001, color='black', length_includes_head=True
                      )
        plt.plot(x, y)
        plt.xlim(0, 2.1)
        plt.ylim(-0.4, 0.4)
        plt.xlabel("x/c")

    else:
        plt.xlim(-4.1, 4.1)
        plt.ylim(-4.1, 4.1)
        plt.title('Velocidades')

    plt.plot(x, y)


    plt.grid(True)
    plt.savefig(f"EP{EP} - Velocidades")


def forca(grid, matriz_pressao):
    """

    :param grid: malha de paineis
    :param matriz_pressao: Matriz com os valores de pressao para cada painel
    :return: força resultante por unidade de comprimento
    """
    F1 = 0
    F2 = 0
    for painel_i in grid.values(): F1 += matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_1 * painel_i.b_j
    for painel_i in grid.values(): F2 += matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_2 * painel_i.b_j

    return (print(
        f"A força resultante por un/comprimento é: F = {round(F1)} i + {round(F2)} j ou |F| = {round(m.sqrt(F1 ** 2 + F2 ** 2))} \n "
        f"O Coeficiente de sustentação (lift) é de C_L = {F2 / (1 / 2 * ro * U_longe_total ** 2)}"))


def calcular_massa_adicional(grid, matriz_sigma, matriz_X, matriz_Y, matriz_Vel_total):
    """

    :param grid: malha com os paineis
    :param matriz_sigma: intensidade sigma de cada painel
    :param matriz_X: coord. local X de i em j
    :param matriz_Y: coord. local Y de i em j
    :param matriz_Vel_total:  velocidade de cada painel
    :return: calcula a massa adicional
    """

    massa_adicional = 0

    for i in grid.values():
        gradiente_n = np.dot([matriz_Vel_total[i.numero - 1][0], matriz_Vel_total[i.numero - 1][1]],
                             [i.versor_normal_1, i.versor_normal_2])  # gradiente de fi vezes normal
        fi = 0

        for j in grid.values():
            sigma = matriz_sigma[j.numero - 1]  # intensidade do painel j
            X = matriz_X[i.numero - 1][j.numero - 1]  # coord. local X de i em j
            Y = matriz_Y[i.numero - 1][j.numero - 1]  # coord. local Y de i em j
            a = 0  # referencial local do painel ser em um dos vértices
            b = j.b_j  # referencial local do painel ser em um dos vértices
            fi += (sigma / (4 * m.pi)) * (
                    - ((X - b) * m.log((X - b) ** 2 + Y ** 2) - 2 * (X - b) + 2 * Y * m.atan2((X - b), Y)) + (
                    (X - a) * m.log((X - a) ** 2 + Y ** 2) - 2 * (X - a) + 2 * Y * m.atan2((X - a), Y)))

        massa_adicional += fi * gradiente_n * i.b_j  # b_j é dSb da fórmula

    return massa_adicional * ro  # multiplica pela densidade


###############################################################
###############   ROTINA MÉTODO DOS PAINEIS    ################
###############################################################

# Grid com todos os paineis
grid = {i: painel.Painel(i, EP) for i in range(1, N + 1)}

# Montagem do sistema matricial
A = []  # Matriz de influencia

matriz_u_globais = []  # Matriz de vel. globais em x (u_ij) --- Usadas para calculos das velocidades
matriz_v_globais = []  # Matriz de vel. globais em y (v_ij) --- Usadas para calculos das velocidades

matriz_X = []  # --- Usadas no calculo de massa adicional
matriz_Y = []  # --- Usadas no calculo de massa adicional

for painel_i in grid.values():  # Ponto P
    A_i = []  # Linha da matriz

    matriz_u_globais_i = []
    matriz_v_globais_i = []

    matriz_X_i = []
    matriz_Y_i = []

    for painel_j in grid.values():  # Todos os paineis
        X, Y = calcular_coord_ref_local(painel_i.colocacaoX, painel_i.colocacaoY,
                                        painel_j.X_inicio, painel_j.Y_inicio, painel_j.beta)

        matriz_X_i.append(X)
        matriz_Y_i.append(Y)

        if EP != 4:  # Velocidades induzidas de Vórtices são o contrário de Fonte
            U, V = calcular_vel_normal_induzida(X, Y, painel_j.b_j, painel_i.numero, painel_j.numero, EP)
        else:
            V, U = calcular_vel_normal_induzida(X, Y, painel_j.b_j, painel_i.numero, painel_j.numero, EP)

        u, v = calcular_vel_global(U, V, painel_j.beta)

        matriz_u_globais_i.append(u)
        matriz_v_globais_i.append(v)

        if painel_i.numero != painel_j.numero:  # Se i =\= j
            A_ij = u * painel_i.versor_normal_1 + v * painel_i.versor_normal_2  # Elemento da matriz

        else:  # Se i == j
            A_ij = - 1 / 2

        if EP == 4 and painel_i.numero == N:  # Condicao de Kutta
            # intensidades no primeiro e ultimo paineis são opostas
            if painel_j.numero == 1 or painel_j.numero == N:
                A_ij = 1
            else:
                A_ij = 0

        A_i.append(A_ij)

    A.append(A_i)
    matriz_u_globais.append(matriz_u_globais_i)
    matriz_v_globais.append(matriz_v_globais_i)
    matriz_X.append(matriz_X_i)
    matriz_Y.append(matriz_Y_i)

# Matriz U x n de dimensao "Nx1"
matriz_Uxn = []
for painel_i in grid.values():

    if EP == 4 and painel_i.numero == N:  # Condicao de Kutta
        matriz_Uxn.append(0)
    else:
        matriz_Uxn.append(
            - np.dot(U_longe[EP], [painel_i.versor_normal_1, painel_i.versor_normal_2]))

# Resolvendo o sistema, obtemos sigma

matriz_sigma = np.linalg.solve(A, matriz_Uxn)

# Calculo das velocidades potenciais nos pontos de colocacao dos paineis
matriz_Vel_total = []  # Matriz possui N termos: (Vel_x, Vel_y, |Vel|)

for i in range(0, N):
    Vel_x = 0
    Vel_y = 0
    for j in range(0, N):
        Vel_x += matriz_sigma[j] * matriz_u_globais[i][j]
        Vel_y += matriz_sigma[j] * matriz_v_globais[i][j]
    matriz_Vel_total.append([Vel_x, Vel_y, m.sqrt(Vel_x ** 2 + Vel_y ** 2)])

# Calculo da pressao no ponto de colocacao e e do adimensional Cp
if EP in (1, 2, 4):
    matriz_pressao = []  # Matriz com os valores de pressao para cada painel
    matriz_Cp = []  # Matriz com os valores de Cp para cada painel

    U_longe_total = m.sqrt((U_longe[EP][0]) ** 2 + (U_longe[EP][1] ** 2))  # modulo da vel. ao longe

    for i in range(0, N):
        pressao = - ro * U_longe_total * matriz_Vel_total[i][0] - (1 / 2 * ro * matriz_Vel_total[i][2] ** 2)
        matriz_pressao.append(pressao)
        matriz_Cp.append((pressao - pressao0) / (1 / 2 * ro * U_longe_total ** 2))  # adimensionalizar
    grafico_Cp(grid, matriz_Cp)  # plotar gráfico

    # Calculo da forca resultante
    forca(grid, matriz_pressao)  # forca total deve ser zero

if EP in (2, 3):
    massa_adicional = calcular_massa_adicional(grid, matriz_sigma, matriz_X, matriz_Y, matriz_Vel_total)

grafico_Vel(grid, matriz_Vel_total)  # grafico da vel. induzida + ao longe

