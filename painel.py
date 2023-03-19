import math as m

# =========== PARAMETROS GEOMÉTRICOS =================
N = 64  # Numero de paineis
R = 1  # Raio da circunf
b = 1  # Semi eixo vertical
a = 3  # Semi eixo horizontal
c = 1  # Corda
r_espessura = 0.12  # Razão de espessura t/c


# ====================================================

class Painel(object):
    """

    Definicao de cada variavel
    primaria dos paineis
    """

    def __init__(self, i, EP):
        self.numero = i  # Numero do painel (sentido anti-horario)

        if EP == 1:
            self.X_inicio = R * m.cos(2 * m.pi * (i - 1) / N)  # CoordenadaX do ponto inicial
            self.Y_inicio = R * m.sin(2 * m.pi * (i - 1) / N)  # CoordenadaY do ponto inicial
            self.X_fim = R * m.cos(2 * m.pi * i / N)  # CoordenadaX do ponto final
            self.Y_fim = R * m.sin(2 * m.pi * i / N)  # CoordenadaY do ponto final

        elif EP == 2 or EP == 3:
            self.X_inicio = a * m.cos(2 * m.pi * (i - 1) / N)  # CoordenadaX do ponto inicial
            self.Y_inicio = b * m.sin(2 * m.pi * (i - 1) / N)  # CoordenadaY do ponto inicial
            self.X_fim = a * m.cos(2 * m.pi * i / N)  # CoordenadaX do ponto final
            self.Y_fim = b * m.sin(2 * m.pi * i / N)  # CoordenadaY do ponto final

        elif EP == 4:

            if i <= N / 2:  # Metade de cima
                self.X_inicio = c - (c * (i - 1) / (N / 2))  # CoordenadaX do ponto inicial
                self.Y_inicio = (r_espessura / 0.2) * (0.2969 * m.sqrt(self.X_inicio) - 0.126 * self.X_inicio
                                                       - 0.3516 * self.X_inicio ** 2 + 0.2843 * self.X_inicio ** 3 - 0.1015 * self.X_inicio ** 4)  # CoordenadaY do ponto inicial
                if i == 1:  # forçando a fechar
                    self.Y_inicio = 0
                self.X_fim = c - (c * i / (N / 2))  # CoordenadaX do ponto final
                self.Y_fim = (r_espessura / 0.2) * (0.2969 * m.sqrt(self.X_fim) - 0.126 * self.X_fim
                                                      - 0.3516 * self.X_fim ** 2 + 0.2843 * self.X_fim ** 3 - 0.1015 * self.X_fim ** 4)  # CoordenadaY do ponto final
            else:  # Metade de baixo
                self.X_inicio = c * (i - (N / 2 + 1)) / (N / 2)  # CoordenadaX do ponto inicial
                self.Y_inicio = -(r_espessura / 0.2) * (0.2969 * m.sqrt(self.X_inicio) - 0.126 * self.X_inicio
                                                        - 0.3516 * self.X_inicio ** 2 + 0.2843 * self.X_inicio ** 3 - 0.1015 * self.X_inicio ** 4)  # CoordenadaY do ponto inicial
                self.X_fim = c * (i - N / 2) / (N / 2)  # CoordenadaX do ponto final
                self.Y_fim = -(r_espessura / 0.2) * (0.2969 * m.sqrt(self.X_fim) - 0.126 * self.X_fim
                                                     - 0.3516 * self.X_fim ** 2 + 0.2843 * self.X_fim ** 3 - 0.1015 * self.X_fim ** 4)  # CoordenadaY do ponto final
                if i == N:  # forçando a fechar
                    self.Y_fim = 0

        self.colocacaoX = (self.X_inicio + self.X_fim) / 2  # CoordenadaX do ponto medio
        self.colocacaoY = (self.Y_inicio + self.Y_fim) / 2  # Coordenada Y do ponto medio

        self.beta = m.atan2(self.Y_fim - self.Y_inicio, self.X_fim - self.X_inicio)  # Angulo do painel

        # if i <= N/4:  # 1 quadrante
        #    beta1 += 0
        # elif i <= N/2 and i > N/4:  # 2 quadrante
        #    beta1 += 2 * m.pi
        # elif i > N/2 and i <= 3 * N/4:  # 3 quadrante
        #    beta1 += 2 * m.pi
        # elif i > 3 * N/4:  # 4 quadrante
        #    beta1 += 0
        # self.beta = beta1

        self.b_j = m.sqrt((self.X_inicio - self.X_fim) ** 2 + (self.Y_inicio - self.Y_fim) ** 2)  # Tamanho do painel

        self.versor_normal_1 = - m.sin(self.beta)  # Versor normal calculado no ponto de colocação do mesmo (1)
        self.versor_normal_2 = m.cos(self.beta)  # Versor normal calculado no ponto de colocação do mesmo (2)

        # self.versor_normal_1 = - (self.Y_fim-self.Y_inicio) / self.b_j
        # self.versor_normal_2 = (self.X_fim-self.X_inicio) / self.b_j

        # Apenas para plotar o grafico:
        if i <= N / 2:
            self.teta = m.degrees(m.atan2(self.colocacaoY, self.colocacaoX))  # Angulo em graus dos pontos de colocacao
        else:
            self.teta = 360 + m.degrees(
                m.atan2(self.colocacaoY, self.colocacaoX))  # Angulo em graus dos pontos de colocacao
