# -*- coding: utf-8 -*-
# This file is a part of the Hidrodinâmica II – PNV3413 project.
# TRABALHO INDIVIDUAL: MÉTODO DOS PAINÉIS
# Autor: Guilherme Donatoni Urbano [guilherme.urbano@usp.br] NUSP: 9835985
# All Rights Reserved

import numpy as np
import math as m
import matplotlib.pyplot as plt

# ==================================== SETUP ===========================================================================
EX = 6  # Selecione o item do trabalho que quer executar (EX = 1,2,3,4,5,6)
# ======================================================================================================================

# ============================== CONSTANTES DO PROBLEMA ================================================================
b = 1  # Semi eixo vertical [m]
ro = 1025  # Densidade [kg/m3]
g = -9.81  # Aceleracao da gravidade [m/s2]
pressao0 = 0  # Pressao ao longe
# =============================== VARIAVEIS DO PROBLEMA  ===============================================================
# N: Numero de paineis
# a: Semi eixo horizontal [m]
# U_longe: Velocidade horizontal do escoamento ao longe [m/s]
# h: Distancia do fundo [m]

# NOTE: Esses parametros são determinados na sessão EXERCICIOS.
# É possivel usar a classe MP() do Metodo dos Paineis para testar variações isoladas do problema.
# ======================================================================================================================

class Painel(object):
    """
    Definicao de cada variavel
    importante dos paineis

    OBS: Não foi empregado a modelagem de corpos simétricos
    """

    def __init__(self, i, a, h, N):
        self.numero = i  # Numero do painel (Numeracao no sentido anti-horario)

        self.X_inicio = a * m.cos(2 * m.pi * (i - 1) / N)  # CoordenadaX do ponto inicial
        self.Y_inicio = b * m.sin(2 * m.pi * (i - 1) / N)  # CoordenadaY do ponto inicial
        self.X_fim = a * m.cos(2 * m.pi * i / N)  # CoordenadaX do ponto final
        self.Y_fim = b * m.sin(2 * m.pi * i / N)  # CoordenadaY do ponto final

        self.colocacaoX = (self.X_inicio + self.X_fim) / 2  # CoordenadaX do ponto medio
        self.colocacaoY = (self.Y_inicio + self.Y_fim) / 2  # CoordenadaY do ponto medio

        self.beta = m.atan2(self.Y_fim - self.Y_inicio, self.X_fim - self.X_inicio)  # Angulo do painel

        self.b_j = m.sqrt((self.X_inicio - self.X_fim) ** 2 + (self.Y_inicio - self.Y_fim) ** 2)  # Tamanho do painel

        self.versor_normal_1 = - m.sin(self.beta)  # Versor normal calculado no ponto de colocação do mesmo (1)
        self.versor_normal_2 = m.cos(self.beta)  # Versor normal calculado no ponto de colocação do mesmo (2)

        # Determinamos o angulo em graus do ponto de colocação (sentido anti-horário):
        if i <= N / 2:
            self.teta = m.degrees(m.atan2(self.colocacaoY, self.colocacaoX))
        else:
            self.teta = 360 + m.degrees(
                m.atan2(self.colocacaoY, self.colocacaoX))  # Angulo em graus dos pontos de colocacao

        # Os efeitos da parede sao como uma imagem do corpo espelhada
        # Os parâmetros geometricos da reflexão do corpo sofrem modificacoes no eixo vertical
        if i > N:
            self.Y_inicio -= (2 * h + 2 * b)
            self.Y_fim -= (2 * h + 2 * b)
            self.colocacaoY -= (2 * h + 2 * b)


def calcular_coord_ref_local(x, y, X_inicio, Y_inicio, beta):
    """
    Calcula as coordenadas de um painel no referencial de outro painel

    :param x: coordenadaX no ref. global
    :param y: coordenadaY no ref. global
    :param X_inicio: coordenadaX no ref. local
    :param Y_inicio: coordenadaX no ref. local
    :param beta: angulo entre ref. local e ref. global
    :return: coordenadas (do ponto P) no referencial local
    """

    X = (x - X_inicio) * m.cos(beta) + (y - Y_inicio) * m.sin(beta)
    Y = - (x - X_inicio) * m.sin(beta) + (y - Y_inicio) * m.cos(beta)

    return X, Y


def calcular_vel_normal_induzida(X, Y, b_j, i, j):
    """
    Calcula a velocidade normal induzida por um painel pela
    distribuicao de fontes sobre outro painel

    :param X: coordenadaX no ref. local
    :param Y: coordendaY no ref. local
    :param b_j: tamanho do painel
    :param i: numero do painel do ponto de colocacao
    :param j: numero do painel da fonte
    :return: velocidade normal induzida no ponto X e Y
    """

    U = (1 / (4 * m.pi)) * m.log((X ** 2 + Y ** 2) / ((X - b_j) ** 2 + Y ** 2))

    if i != j:
        V = (1 / (2 * m.pi)) * (m.atan(np.float64(X) / Y) - m.atan(np.float64(X - b_j) / Y))
    else:  # Contribuicao do painel no proprio painel
        V = - 1 / 2

    return U, V


def calcular_vel_global(U, V, beta):
    """
    Transforma as componentes de velocidade do ref. local para ref. global

    :param U: velocidade normal ind. em X
    :param V: velocidade normal ind. em Y
    :param beta: angulo entre ref. local e global
    :return: velocidade normal induzida no ponto X e Y (ref. global)
    """
    u = U * m.cos(beta) - V * m.sin(beta)
    v = U * m.sin(beta) + V * m.cos(beta)

    return u, v


class MP(object):
    """
    Cria uma instancia de solução do metodo dos paineis 2D
    sem superficie livre para o problema de 'squat'

    """

    def __init__(self, N, a, U_longe, h):

        self.U_longe = U_longe
        self.h = h
        self.N = N
        self.a = a
        self.grid = {i: Painel(i, a, h, N) for i in range(1, 2 * N + 1)}  # Criacao do grid
        self.matriz_sigma, self.matriz_Vel_total = self.cacula_sigma_e_VelTotal()  # intensidades, vel.induzidas
        self.matriz_pressao, self.matriz_Cp = self.Cp()  # pressoes, Cp

    def cacula_sigma_e_VelTotal(self):
        """
        Resolve o sistema característico do problema encontrando as intensidades das fontes
        de cada painel. Tambem calcula a velocidade total induzida no ponto de colocação

        :return: lista com intensidades de fontes dos paineis

        """

        # Montagem do sistema matricial
        A = []  # Matriz de influencia

        matriz_u_globais = []  # Matriz de vel. globais em x (u_ij) --- Usadas para calculos das velocidades
        matriz_v_globais = []  # Matriz de vel. globais em y (v_ij) --- Usadas para calculos das velocidades

        for painel_i in self.grid.values():  # painel_i recebe as contribuicoes das fontes
            A_i = []  # Linha da matriz

            matriz_u_globais_i = []
            matriz_v_globais_i = []

            for painel_j in self.grid.values():  # painel_j induz as contribuicoes nos painel_i
                X, Y = calcular_coord_ref_local(painel_i.colocacaoX, painel_i.colocacaoY,
                                                painel_j.X_inicio, painel_j.Y_inicio, painel_j.beta)

                U, V = calcular_vel_normal_induzida(X, Y, painel_j.b_j, painel_i.numero, painel_j.numero)

                u, v = calcular_vel_global(U, V, painel_j.beta)

                matriz_u_globais_i.append(u)
                matriz_v_globais_i.append(v)

                if painel_i.numero != painel_j.numero:  # Se i =\= j
                    A_ij = u * painel_i.versor_normal_1 + v * painel_i.versor_normal_2  # Elemento da matriz

                else:  # Se i == j
                    A_ij = - 1 / 2

                A_i.append(A_ij)

            A.append(A_i)
            matriz_u_globais.append(matriz_u_globais_i)
            matriz_v_globais.append(matriz_v_globais_i)

        # Matriz - [U_longe x vetor normal]
        matriz_Uxn = []
        for painel_i in self.grid.values():
            matriz_Uxn.append(- np.dot([self.U_longe, 0], [painel_i.versor_normal_1, painel_i.versor_normal_2]))

        # Resolvendo o sistema, obtemos a intensidade dos paineis
        matriz_sigma = np.linalg.solve(A, matriz_Uxn)

        # Calculo das velocidades induzidas nos pontos de colocacao dos paineis
        matriz_Vel_total = []  # Matriz possui N termos: (Vel_x, Vel_y, |Vel|)

        for i in range(0, 2 * self.N):
            Vel_x = 0
            Vel_y = 0
            for j in range(0, 2 * self.N):
                Vel_x += matriz_sigma[j] * matriz_u_globais[i][j]
                Vel_y += matriz_sigma[j] * matriz_v_globais[i][j]
            matriz_Vel_total.append([Vel_x, Vel_y, m.sqrt(Vel_x ** 2 + Vel_y ** 2)])

        return matriz_sigma, matriz_Vel_total

    def Cp(self):
        """
        Retorna os valores de pressao e os coeficientes de pressao (Cp) sobre
        o corpo submerso e a reflexão

        :return: matriz com os valores dos coeficientes de pressao

        """

        # Calculo da pressao no ponto de colocacao e e do adimensional Cp
        matriz_pressao = []  # Matriz com os valores de pressao para cada painel
        matriz_Cp = []  # Matriz com os valores de Cp para cada painel

        for i in range(0, 2 * self.N):
            pressao = - ro * self.U_longe * self.matriz_Vel_total[i][0] - (
                    1 / 2 * ro * self.matriz_Vel_total[i][2] ** 2)
            matriz_pressao.append(pressao)
            matriz_Cp.append((pressao - pressao0) / (1 / 2 * ro * self.U_longe ** 2))  # adimensionalizar

        return matriz_pressao, matriz_Cp


# ================================= FUNÇÕES ============================================================================

def graficos_grid_e_Vel(grid, matriz_Vel_total, N, h, a, U_longe):
    """
    Plota o grid utilizado com a numeracao dos paineis
    (as velocidades nos paineis também se descomentados)

    :param grid: malha com os paineis
    :param matriz_Vel_total: matriz com as velocidades induzidas
    :param N: numero de paineis
    :param h: distancia
    :param a: semi eixo horizontal
    :return: grafico da velocidade induzida + ao longe

    """
    x = []  # Armazena a coord. x
    y = []  # Armazena a coord. y
    controle_X = []  # Armazena os pontos de controle (coord. X)
    controle_Y = []  # Armazena os pontos de controle (coord. Y)

    for painel_i in grid.values():
        x.append(painel_i.X_inicio)
        y.append(painel_i.Y_inicio)

        controle_X.append(painel_i.colocacaoX)
        controle_Y.append(painel_i.colocacaoY)

        ####### Para colocar setas com velocidades nos pontos de colocacao, descomentar:
        # plt.arrow(painel_i.colocacaoX, painel_i.colocacaoY,
        # matriz_Vel_total[painel_i.numero - 1][0] + U_longe,
        # matriz_Vel_total[painel_i.numero - 1][1],
        # head_width=0.07, width=0.001, color='black', length_includes_head=True
        # )

        # Para "fechar" o grid com o ultimo ponto
        if painel_i.numero == N:
            x.append(a)
            y.append(0)
        if painel_i.numero == 2 * N:
            x.append(a)
            y.append(-(2 * h + 2 * b))

        plt.scatter(controle_X, controle_Y, s=7, color='red')

        if painel_i.numero <= N:
            plt.text(painel_i.colocacaoX, painel_i.colocacaoY, painel_i.numero, **dict(size=7.5, color='red'))

    plt.plot(x[0:N + 1], y[0: N + 1])
    plt.plot(x[N + 1:2 * N + 2], y[N + 1:2 * N + 2])

    plt.axis('scaled')

    # Para só o corpo usar limites de eixoX: (-a - 0.1 , a + 0.1)
    # Para a imagem usar limites eixoX: (-5.1, -2.9)
    # Para ver corpo e imagem, comentar os dois limites de eixos
    plt.xlim(-a - 0.1, a + 0.1)
    plt.ylim(-1.1, 1.1)

    plt.title(f'Grid (N = {N}, h = {h}m)')
    plt.grid(True)
    plt.show()


def grafico_Cp(Instancias, N):
    """
    Plota o gráfico do coeficiente de pressão sobre o corpo

    :param Instancias: instancias com as resoluções dos metodos dos paineis
    :param N: numero de paineis
    :return: grafico Cp x teta

    """

    for MP in Instancias.values():
        matriz_teta = []

        for painel_i in MP.grid.values():
            if painel_i.numero <= N:
                matriz_teta.append(painel_i.teta)
        matriz_teta.sort()

        plt.plot(matriz_teta, MP.matriz_Cp[0:N], marker='o', label=f"h ={MP.h}")
    plt.title(f'Cp x Teta (°) (a/b = {a / b})')
    plt.xlabel('Teta (°)')
    plt.ylabel('Cp')
    plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
    plt.grid(True)
    plt.show()


def convergencia_forca(Instancias, N_max):
    """
    Calcula as forças hidrodinamicas (horizontal e vertical) para diferentes numero de paineis
    e informa a convergencia numerica

    :param Instancias: instancias com as resoluções dos metodos dos paineis
    :return: grafico com forças (horizontal e vertical) de cada uma das instancias
    """

    Fx = []
    Fy = []
    paineis = []

    for MP in Instancias.values():
        F1 = 0
        F2 = 0
        for painel_i in MP.grid.values():
            if painel_i.numero <= MP.N:
                F1 += MP.matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_1 * painel_i.b_j
                F2 += MP.matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_2 * painel_i.b_j
        Fx.append(F1)
        Fy.append(F2)
        paineis.append(MP.N)

    plt.subplot(2, 1, 1)
    plt.title(f'Força x Nº de Painéis')
    plt.ylabel('Força Horizontal [N]')
    plt.grid(True)
    plt.plot(paineis, Fx, color='red')

    plt.subplot(2, 1, 2)
    plt.plot(paineis, Fy, color='blue')
    plt.xlabel('Nº de Painéis')
    plt.ylabel('Força Vertical [N]')
    plt.grid(True)

    plt.show()

    # Analise da convergencia numerica com erro relativo de 0,01%
    print("-----------Convergência Numérica-------------")
    print("")

    # Para a forca horizontal:
    for i in range(1, int((N_max - 2) / 2)):
        if (abs(Fx[i + 1] - Fx[i])) / Fx[i] < 0.00001:
            print(f"Para Força Horizontal com erro de 0,001% em: N = {2 * i + 4} \n"
                  f"(Fx = {Fx[i + 1]} [N])")
            break

    print("")

    # Para a forca vertical:
    for i in range(1, int((N_max - 2) / 2)):
        if (Fy[i + 1] - Fy[i]) / Fy[i] < 0.00001:
            print(f"Para Força Vertical com erro de 0,001% em: N = {2 * i + 4} \n"
                  f"(Fy = {Fy[i + 1]} [N])")
            break

    print("")

    # Para o modulo da forca:
    for i in range(1, int((N_max - 2) / 2)):
        if abs(m.sqrt(Fy[i + 1] ** 2 + Fx[i + 1] ** 2) - m.sqrt(Fy[i] ** 2 + Fx[i] ** 2)) / m.sqrt(
                Fy[i] ** 2 + Fx[i] ** 2) < 0.00001:
            print(f"Para Força Total com erro de 0,001% em: N = {2 * i + 4} \n"
                  f"(|F| = {m.sqrt(Fy[i + 1] ** 2 + Fx[i + 1] ** 2)} [N])")
            N_numerico = 2 * i + 4
            break

    print("--------------------------------------------")
    print("-----------Comparação com Analítico-------------")
    # Comparacao com o valor analitico (módulo) com erro relativo de de 1%

    for i in range(1, int((N_max - 2) / 2)):
        MP_analitico = Instancias[2 * i + 4]
        F = (m.pi / 2) * ((MP_analitico.a) ** 4 / (MP_analitico.a + MP_analitico.h) ** 3) * ro * MP_analitico.U_longe

        if (abs((m.sqrt(Fx[i] ** 2 + Fy[i] ** 2) - F)) / F) < 0.01:
            print(
                f"Diferença percentual para N = {N_numerico}: {round(abs(F - m.sqrt(Fy[int((N_numerico - 4) / 2)] ** 2 + Fx[int((N_numerico - 4) / 2)] ** 2)) * 100 / F, 2)}% \n"
                f"Para atingir Valor Analítico com erro de 1% em: N = {2 * i + 2} \n"
                f"(|F| = {F} [N])")
            break
    print("--------------------------------------------")

    return N_numerico


def forca(Instancias, caso):
    """
    Plota os graficos das forças em funcao da distancia e em funcao da velocidade
    de avanço

    :param Instancias: instancias com as resoluções dos metodos dos paineis
    :param caso: variar distancia (caso 1) e variar velocidade (caso 2)
    :return: grafico com forças (horizontal e vertical) quando varia-se  a
            distancia (caso 1) e a velocidade (caso 2)
    """

    if caso == 1:  # força x distancia
        for i in range(0, len(Instancias)):
            Fx = []
            Fy = []
            distancia = []
            for MP in Instancias[i]:
                F1 = 0
                F2 = 0
                for painel_i in MP.grid.values():
                    if painel_i.numero <= MP.N:
                        F1 += MP.matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_1 * painel_i.b_j
                        F2 += MP.matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_2 * painel_i.b_j
                Fx.append(F1)
                Fy.append(F2)
                distancia.append(MP.h)

            plt.plot(distancia, Fx, label=f"Força Horizontal (a/b ={MP.a})")
            plt.plot(distancia, Fy, label=f"Força Vertical (a/b ={MP.a})")

        plt.title(f'Força x Distância')
        plt.xlabel('h [m]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
        plt.show()

    if caso == 2:  # força x velocidade de avanço
        for i in range(0, len(Instancias)):
            Fx = []
            Fy = []
            U = []
            for MP in Instancias[i]:
                F1 = 0
                F2 = 0
                for painel_i in MP.grid.values():
                    if painel_i.numero <= MP.N:
                        F1 += MP.matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_1 * painel_i.b_j
                        F2 += MP.matriz_pressao[painel_i.numero - 1] * painel_i.versor_normal_2 * painel_i.b_j
                Fx.append(F1)
                Fy.append(F2)
                U.append(MP.U_longe)

            plt.plot(U, Fx, label=f"Força Horizontal (a/b ={MP.a})")
            plt.plot(U, Fy, label=f"Força Vertical (a/b ={MP.a})")

        plt.title(f'Força x Velocidade de avanço')
        plt.xlabel('U [m/s]')
        plt.ylabel('Força [N]')
        plt.grid(True)
        plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
        plt.show()


def mapa_velocidades(MP):
    """
    Plota o mapa da intensidades de velocidades do escoamento

    :param MP: instancia do metodo dos paineis
    :return: mapa com intensidade
    """

    Pontos = 100  # Numero de pontos fora do corpo que serao plotados (relatorio está 200)

    Velocidade = np.zeros((Pontos, Pontos))  # Armazena as velocidades dos pontos

    x = np.linspace(-5, 5, num=Pontos)
    y = np.linspace(-2 * b, 2 * b, num=Pontos)
    A, B = np.meshgrid(x, y)

    for i in range(0, Pontos):
        for j in range(0, Pontos):

            if ((x[j] ** 2 / MP.a ** 2) + (y[i] ** 2 / b ** 2)) > 1:  # Se nao pertence a elipse
                Vel_x = MP.U_longe  # Todos os pontos tem no minimo a vel. ao longe
                Vel_y = 0

                for painel_j in MP.grid.values():  # painel_j induz as contribuicoes nos pontos P
                    X, Y = calcular_coord_ref_local(x[j], y[i], painel_j.X_inicio, painel_j.Y_inicio, painel_j.beta)

                    U = (1 / (4 * m.pi)) * m.log((X ** 2 + Y ** 2) / ((X - painel_j.b_j) ** 2 + Y ** 2))

                    V = (1 / (2 * m.pi)) * (m.atan(np.float64(X) / Y) - m.atan(np.float64(X - painel_j.b_j) / Y))

                    u, v = calcular_vel_global(U, V, painel_j.beta)

                    Vel_x += MP.matriz_sigma[painel_j.numero - 1] * u
                    Vel_y += MP.matriz_sigma[painel_j.numero - 1] * v

                Velocidade[i][j] = m.sqrt(Vel_x ** 2 + Vel_y ** 2)

    levels = np.linspace(0, 2.5, 26)  # Valores de altura no plano x-y
    plt.contourf(A, B, Velocidade, levels=levels, cmap='seismic')

    plt.title(f"Mapa de velocidades (h = {MP.h}m)")
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.axis('scaled')
    plt.colorbar()
    plt.show()


# ======================================================================================================================


# ================================================ EXERCICIOS ==========================================================
# Pode-se alterar os parametros de cada exercicio


if EX == 1:
    N = 64
    h = 1
    a = 1
    U_longe = 1

    # Plota os paineis com numeracao para o corpo, para plotar reflexão, seguir instrucoes da funcao
    MP = MP(N=N, a=a, U_longe=1, h=h)
    graficos_grid_e_Vel(MP.grid, MP.matriz_Vel_total, N, h, a, U_longe)

if EX == 2:
    N = 32
    U_longe = 1
    h = 1
    a = (1, 3)  # a/b = 1.0 e a/b = 3.0
    Instancias = {a: MP(N=N, a=a, U_longe=1, h=h) for a in a}

    for MP in Instancias.values():
        print(f"---------------------------[a/b = {MP.a}]------------------------------------")
        print(f"Nº Intensidade [m²/s] Módulo de velocidade [m/s]")

        for i in range(1, 2 * N + 1):
            # A velocidade total do escoamento é dada pela soma da velocidade induzida (em modulo)
            # e da velocidade no infinito
            print("%s %s %s" % (i,
                                MP.matriz_sigma[i - 1],
                                m.sqrt((MP.matriz_Vel_total[i - 1][0] + U_longe) ** 2 +
                                       MP.matriz_Vel_total[i - 1][1] ** 2)))
        print("----------------------------------------------------------------------------")

if EX == 3:
    N = 64
    U_longe = 1
    h = (1, 2)  # Distancias h = 1m e h = 2m
    for a in (1.0, 1.5, 2.0, 2.5, 3.0):  # Um grafico para cada a/b com duas curvas de h
        Instancias = {h: MP(N=N, a=a, U_longe=U_longe, h=h) for h in h}
        grafico_Cp(Instancias, N)  # plota grafico do coeficiente de pressao

if EX == 4:
    U_longe = 1
    h = 2
    N_max = 128  # Numeros de paineis que analisaremos na convergencia numerica
    N = []
    for i in range(4, N_max, 2): N.append(i)  # Passo de 2 paineis
    Instancias = {N: MP(N=N, a=1, U_longe=U_longe, h=h) for N in N}

    # Plota o grafico com as forças horizontais e verticais e chega a convergencia numerica de N
    N_numerico = convergencia_forca(Instancias, N_max)

if EX == 5:
    N = 48  # Foi definido no EX4 que deve-se usar essa quantidade de paineis
    a = (1, 3)

    # Caso 1: aumento da distancia
    h = np.linspace(1.0, 20.0, num=20)

    Instancias = []

    for comp in a:
        Inst = []
        for distancia in h:
            MP1 = MP(N=N, a=comp, U_longe=1, h=distancia)
            Inst.append(MP1)
        Instancias.append(Inst)
    forca(Instancias, caso=1)

    # Caso 2: aumento da velocidade de avanço
    U_longe = np.linspace(0.5, 4.0, num=14)

    Instancias = []

    for comp in a:
        Inst = []
        for U in U_longe:
            MP1 = MP(N=N, a=comp, U_longe=U, h=2)
            Inst.append(MP1)
        Instancias.append(Inst)
    forca(Instancias, caso=2)

if EX == 6:
    N = 48  # Foi definido no EX4 que deve-se usar essa quantidade de paineis
    a = 2
    U_longe = 1
    h = (1.5, 2.0, 4.0, 30.0)

    Instancias = {h: MP(N=N, a=a, U_longe=U_longe, h=h) for h in h}
    for MP in Instancias.values():
        mapa_velocidades(MP)  # Para alterar a resolucao, alterar "Pontos" na funcao

# ==================================================================================================
