import numpy as np
import matplotlib.pyplot as plt 
import math
from sympy import *
from sympy.solvers.ode.systems import *
from matplotlib.patches import Arc, RegularPolygon


class Retangulo:
    """
    Essa classe representa seções transversais retangulares e suas principais características.
    """

    def __init__(self, b, h):
        """
        :b (float): tamanho da base do retângulo em metros.
        :h (float): tamanho da altura do retângulo em metros.

        Método que inicializa objetos do tipo Retangulo.
        """

        self.b = b
        self.h = h
        self.I = self.momento_inercia()

    def momento_inercia(self):
        """
        :output I (float): momento de inércia. 

        Método responsável pelo cálculo do momento de inercia da seção transversal, neste caso, um retângulo.
        """

        I = (self.b*self.h**3)/12

        return I

    def esquema(self, estilo='--', cor='black', cor_preench='black'):
        """
        :estilo (str): estilo a ser utilizado na borda da seção transversal.
        :cor (str): cor a ser utilizada na borda da seção transversal.
        :cor_preench (str): cor a ser utilizada no preenchimento da seção transversal.

        Método responsável pela representação gráfica da seção transversal retangular.
        """

        x = np.arange(-self.b/2, self.b/2, 0.0001)
        y = np.arange(-self.h/2, self.h/2, 0.0001)

        ax = plt.subplot(111)

        plt.plot(x, -self.h/2+0*x, estilo, c=cor)
        plt.plot(x, self.h/2+0*x, estilo, c=cor)
        plt.plot(self.b/2+0*y, y, estilo, c=cor)
        plt.plot(-self.b/2+0*y, y, estilo, c=cor)

        plt.fill_between(x, -self.h/2, self.h/2, color=cor_preench, alpha=0.2, hatch=".")

        plt.axis('equal')
        plt.xticks([-self.b/2, self.b/2])
        plt.yticks([-self.h/2, self.h/2])

        ax.spines['left'].set_position(('data', 0))
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_color('none')

        plt.show()


class Circulo:
    """
    Essa classe representa seções transversais circulares e suas principais características.
    """

    def __init__(self, r):
        """
        :r (float): tamanho do raio do círculo em metros.

        Método que inicializa objetos do tipo Circulo.
        """

        self.r = r 
        self.I = self.momento_inercia()

    def momento_inercia(self):
        """
        :output I (float): momento de inércia. 

        Método responsável pelo cálculo do momento de inercia da seção transversal, neste caso, um círculo.
        """

        I = np.pi*self.r**4/4

        return I

    def esquema(self, estilo='--', cor='black', cor_preench='black'):
        """
        :estilo (str): estilo a ser utilizado na borda da seção transversal.
        :cor (str): cor a ser utilizada na borda da seção transversal.
        :cor_preench (str): cor a ser utilizada no preenchimento da seção transversal.

        Método responsável pela representação gráfica da seção transversal circular.
        """

        theta = np.arange(0, 2*np.pi, 0.0001)
        x = self.r*np.cos(theta)
        y = self.r*np.sin(theta)

        ax = plt.subplot(111)

        plt.plot(x, y, estilo, c=cor)

        plt.fill_between(x, y, -y, color=cor_preench, alpha=0.2, hatch=".")

        plt.axis('equal')
        plt.xticks([-self.r, self.r])
        plt.yticks([-self.r, self.r])

        ax.spines['left'].set_position(('data', 0))
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_color('none')

        plt.show()


class Fixação:
    """
    Essa classe representa fixações e suas principais características.
    """

    def __init__(self, x, y, fixa_x=False, fixa_y=False, fixa_rot=False):
        """
        :x (float): posição da fixação no eixo x.
        :y (float): posição da fixação no eixo y.
        :fixa_x (bool): indica se a fixação permite ou não a movimentação no eixo x.
        :fixa_y (bool): indica se a fixação permite ou não a movimentação no eixo y.
        :fixa_rot (bool): indica se a fixação permite ou não a rotação em torno do eixo z.

        Método que inicializa objetos do tipo Fixação.
        """

        self.x = x
        self.y = y
        self.fixa_x = fixa_x
        self.fixa_y = fixa_y
        self.fixa_rot = fixa_rot

        self.posição = np.array([x, y])

        self.mag_x = None
        self.mag_y = None

    def esquema(self, ax, tamanho = 8, cor = 'red'):
        """
        :ax (object): eixo para plotar a fixação.
        :tamanho (float): tamanho do marcador da fixação.
        :cor (str): cor a ser utilizada na fixação.

        Método responsável pela representação gráfica da fixação.
        """

        if self.fixa_x and self.fixa_y and self.fixa_rot:

            ax.plot(self.x, self.y, 'x', markersize=tamanho, c=cor)

        elif self.fixa_x and self.fixa_y:

            ax.plot(self.x, self.y, '.', c=cor)

        elif self.fixa_x:

            ax.plot(self.x, self.y, marker=4, markersize=tamanho, c=cor)

        elif self.fixa_y:

            ax.plot(self.x, self.y, marker=6, markersize=tamanho, c=cor)


class CarregamentoPontual:
    """
    Essa classe representa carregamentos pontuais e suas principais características.
    """

    def __init__(self, x, y, mag_x, mag_y):
        """
        :x (float): posição da aplicação da carga no eixo x.
        :y (float): posição da aplicação da carga no eixo y.
        :mag_x (float): magnitude da força no direção do eixo x.
        :mag_y (float): magnitude da força no direção do eixo y.

        Método que inicializa objetos do tipo CarregamentoPontual.
        """

        self.x = float(x)
        self.y = float(y)
        self.mag_x = float(mag_x)
        self.mag_y = float(mag_y)

        self.posição = np.array([x, y])
        self.vetor = np.array([self.mag_x, self.mag_y])

    def esquema(self, ax, cor = 'red', tamanho = 0.004, raio=1):
        """
        :ax (object): eixo para plotar o carregamento.
        :cor (str): cor a ser utilizada no vetor.
        :tamanho (float): tamanho da grossura do vetor.

        Método responsável pela representação gráfica do carregamento.
        """

        if self.vetor.any() != 0:

            ax.quiver(self.x, self.y, self.vetor[0], self.vetor[1], pivot='tip', color=cor, width=tamanho, ls='-')


class Momento:
    """
    Essa classe representa momentos em torno do eixo z e suas principais características.
    """

    def __init__(self, x, y, mag_z):
        """
        :x (float): posição da aplicação do momento no eixo x.
        :y (float): posição da aplicação do momento no eixo y.
        :mag_z (float): magnitude do momento no direção do eixo z.

        Método que inicializa objetos do tipo Momento.
        """

        self.x = float(x)
        self.y = float(y)
        self.mag_z = float(mag_z)

        self.posição = np.array([x, y])

    def esquema(self, ax, cor = 'red', raio = 1):
        """
        :ax (object): eixo para plotar o momento.
        :cor (str): cor a ser utilizada na representação do momento.
        :raio (float): raio da representação do momento.

        Método responsável pela representação gráfica do momento.
        """

        arc = Arc([self.x , self.y], theta1=40, theta2=320, color=cor, height=raio, width=raio)
        
        ax.add_patch(arc)

        if self.mag_z >= 0:

            x_seta = self.x + (raio/2)*np.cos(np.deg2rad(320)) 
            y_seta = self.y + (raio/2)*np.sin(np.deg2rad(320))

            ax.add_patch(RegularPolygon((x_seta, y_seta), 3, raio/9, rad(320), color=cor))

        else: 

            x_seta = self.x + (raio/2)*np.cos(np.deg2rad(40)) 
            y_seta = self.y + (raio/2)*np.sin(np.deg2rad(40))

            ax.add_patch(RegularPolygon((x_seta, y_seta), 3, raio/9, rad(50+60), color=cor))


class EsforçoCortante:
    """
    Essa classe representa esforços cortantes e suas principais características.
    """

    def __init__(self, func, intervalo, expr):
        """
        :func (function): função que mapeia uma posição na barra a um esforço cortante.
        :intervalo (array): intervalo de posições da barra que o esforço atua.
        :expr (object): expressão do esforço cortante.

        Método que inicializa objetos do tipo EsforçoCortante.
        """

        self.func = func
        self.intervalo = intervalo
        self.expr = expr

    def esquema(self, ax, cor = 'red'):
        """
        :ax (object): eixo para plotar o esforço cortante.

        Método responsável pela representação gráfica do esforço cortante.
        """
        
        Q = self.func(self.intervalo)

        ax.plot(self.intervalo, Q, c=cor, label='Esforço Cortante')


class MomentoFletor:
    """
    Essa classe representa momentos fletores e suas principais características.
    """

    def __init__(self, func, intervalo, expr):
        """
        :func (function): função que mapeia uma posição na barra a um momento fletor.
        :intervalo (array): intervalo de posições da barra que o momento atua.
        :expr (object): expressão do momento.

        Método que inicializa objetos do tipo MomentoFletor.
        """

        self.func = func
        self.intervalo = intervalo
        self.expr = expr

    def esquema(self, ax, cor = 'blue'):
        """
        :ax (object): eixo para plotar o momento fletor.

        Método responsável pela representação gráfica do momento fletor.
        """

        M = self.func(self.intervalo)

        ax.plot(self.intervalo, M, c=cor, label='Momento Fletor')


class Deflexão:
    """
    Essa classe representa deflexões e suas principais características.
    """

    def __init__(self, func, intervalo):
        """
        :func (function): função que mapeia uma posição na barra a uma deflexão.
        :intervalo (array): intervalo de posições da barra que a deflexão atua.

        Método que inicializa objetos do tipo Deflexão.
        """

        self.func = func
        self.intervalo = intervalo

    def esquema(self, ax, estilo='--', cor = 'black'):
        """
        :ax (object): eixo para plotar a delfexão.

        Método responsável pela representação gráfica da deflexão.
        """

        d = self.func(self.intervalo)

        ax.plot(self.intervalo, d, estilo, c=cor)
        

class Barra:
    """
    Essa classe representa barras e suas principais características.
    """

    def __init__(self, x_inicio, y_inicio, inclinação, comprimento, seção, E):
        """
        :x_inicio (float): posição de inicio da barra no eixo x.
        :y_inicio (float): posição de inicio da barra no eixo y.
        :inclinação (float): inclinação da barra em radianos.
        :comprimento (float): comprimento da barra em metros.
        :seção (object): objeto de seção transversal.
        :E (float): módulo de elasticidade.

        Método que inicializa objetos do tipo Barra.
        """

        self.x_inicio = x_inicio
        self.y_inicio = y_inicio

        self.inicio = np.array([x_inicio, y_inicio])
        
        self.inclinação = inclinação
        self.comprimento = comprimento 
        self.seção = seção
        self.E = E

        self.I_seção = self.seção.I

        self.fixações = []
        self.carregamentos = []
        self.esforços_externos = []
        self.esforços_internos = []
        self.deflexões = []

        self.t = np.arange(0, self.comprimento, 0.0001)

        self.x = self.x_inicio + self.t*np.cos(self.inclinação)
        self.y = self.y_inicio + self.t*np.sin(self.inclinação)

        self.x_final = self.x_inicio + self.comprimento*np.cos(self.inclinação)
        self.y_final = self.y_inicio + self.comprimento*np.sin(self.inclinação)

        self.nos = [(self.x_inicio, self.y_inicio), (self.x_final, self.y_final)]

    def adicionar_fixação(self, fixação):
        """
        :fixação (object): objeto de fixação a ser adicionado à estrutura.

        Método que adiciona fixações à barra.
        """
        
        try:

            inclinação_fixação = (fixação.y - self.y_inicio)/(fixação.x - self.x_inicio)
        
        except:

            inclinação_fixação = 0

        condição_inclinação = math.isclose(inclinação_fixação, np.tan(self.inclinação), rel_tol = 0.001)

        condição_intervalo = (min(self.x) <= fixação.x <= max(self.x)) and (min(self.y) <= fixação.y <= max(self.y))

        condição_bordas = (((math.isclose(fixação.x, min(self.x), rel_tol = 0.001)) and (math.isclose(fixação.y, min(self.y), rel_tol = 0.001))) or 
                            (math.isclose(fixação.x, max(self.x), rel_tol = 0.001) and (math.isclose(fixação.y, max(self.y), rel_tol = 0.001))))

        if condição_inclinação and (condição_intervalo or condição_bordas):

            self.fixações += [fixação]
            self.nos += [(fixação.x, fixação.y)]

        elif condição_bordas:

            self.fixações += [fixação]
            self.nos += [(fixação.x, fixação.y)]

        else:

            raise ValueError('Fixação não está contida no objeto alvo.')

    def adicionar_carregamento(self, carregamento):
        """
        :carregamento (object): objeto de carregamento a ser adicionado à estrutura.

        Método que adiciona carregamentos pontuais à barra.
        """
        
        try:

            inclinação_fixação = (carregamento.y - self.y_inicio)/(carregamento.x - self.x_inicio)
        
        except:

            inclinação_fixação = 0

        condição_inclinação = math.isclose(inclinação_fixação, np.tan(self.inclinação), rel_tol = 0.001)

        condição_intervalo = (min(self.x) <= carregamento.x <= max(self.x)) and (min(self.y) <= carregamento.y <= max(self.y))

        condição_bordas = (((math.isclose(carregamento.x, min(self.x), rel_tol = 0.001)) and (math.isclose(carregamento.y, min(self.y), rel_tol = 0.001))) or 
                            (math.isclose(carregamento.x, max(self.x), rel_tol = 0.001) and (math.isclose(carregamento.y, max(self.y), rel_tol = 0.001))))

        if condição_inclinação and (condição_intervalo or condição_bordas):

            self.carregamentos += [carregamento]
            self.nos += [(carregamento.x, carregamento.y)]

        else:

            raise ValueError('Carregamento não está contida no objeto alvo.')
    
    def esquema(self, estilo='-', cor='black'):
        """
        :estilo (str): estilo do traço que representa a barra.
        :cor (str): cor a ser utilizada na representação da barra.

        Método responsável pela representação gráfica da barra.
        """

        plt.figure('Esquema da Barra')
        ax = plt.subplot(111)

        plt.plot(self.x, self.y, estilo, c=cor)

        ax.spines['left'].set_position(('data', 0))
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_color('none')

        if self.esforços_externos == []:

            plt.axis('equal')
            plt.axis('off')

            for carregamento in self.carregamentos:

                carregamento.esquema(ax)

            for fixação in self.fixações:

                fixação.esquema(ax)

        elif self.esforços_internos == []:

            plt.axis('equal')
            plt.axis('off')

            for carregamento in self.carregamentos:

                carregamento.esquema(ax)

            for esforço_externo in self.esforços_externos:
                
                esforço_externo.esquema(ax, cor='blue', raio=self.comprimento/15)        

        elif self.deflexões == []:

            plt.axis('off')
            plt.tight_layout()
            
            for fixação in self.fixações:

                fixação.esquema(ax)

            for esforço_interno in self.esforços_internos:
                
                esforço_interno.esquema(ax)        

            
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            
            plt.legend(by_label.values(), by_label.keys(), frameon=False)

        else:
            
            plt.axis('off')
            
            for fixação in self.fixações:

                fixação.esquema(ax)

            for deflexão in self.deflexões:
                
                deflexão.esquema(ax)            

        plt.show()

    def solução_esforços_externos(self):
        """
        Método responsável pela solução dos esforços externos aplicados nas fixações.

        :output resultados (dict): dicionário com as soluções para os esforços externos aplicados na estrutura. 
        """

        incognitas_x = []
        incognitas_y = []
        incognitas_M = []

        soma_incognitas_x = ''
        soma_incognitas_y = ''
        soma_incognitas_M = ''

        carregamentos_x = []
        carregamentos_y = []
        carregamentos_M = []

        for i, fixação in enumerate(self.fixações):

            if fixação.fixa_x:

                incognitas_x += ['x_' + str(i)]
                soma_incognitas_x += '+x_' + str(i)

                soma_incognitas_M += '+x_' + str(i) + '*' + str((fixação.posição - self.inicio)[1])

            if fixação.fixa_y:
                
                incognitas_y += ['y_' + str(i)]
                soma_incognitas_y += '+y_' + str(i)

                soma_incognitas_M += '-y_' + str(i) + '*' + str((fixação.posição - self.inicio)[0])

            if fixação.fixa_rot:

                incognitas_M += ['M_' + str(i)]
                soma_incognitas_M += '-M_' + str(i)

        for carregamento in self.carregamentos:

                carregamentos_x += [carregamento.mag_x]
                carregamentos_y += [carregamento.mag_y]
                carregamentos_M += [-float(np.cross(carregamento.posição - self.inicio, carregamento.vetor))]

        soma_carregamentos_x = sum(carregamentos_x)
        soma_carregamentos_y = sum(carregamentos_y)
        soma_carregamentos_M = sum(carregamentos_M)

        equilibrio_x = soma_incognitas_x + '+(' + str(soma_carregamentos_x) + ')'
        equilibrio_y = soma_incognitas_y + '+(' + str(soma_carregamentos_y) + ')'
        equilibrio_M = soma_incognitas_M + '+(' + str(soma_carregamentos_M) + ')'

        expr_equilibrio_x = sympify(equilibrio_x)
        expr_equilibrio_y = sympify(equilibrio_y)
        expr_equilibrio_M = sympify(equilibrio_M)

        incognitas = incognitas_x + incognitas_y + incognitas_M
        incognitas.sort()

        solução = solve_poly_system([expr_equilibrio_x, expr_equilibrio_y, expr_equilibrio_M])

        resultados = dict(zip(incognitas, solução[0]))

        for variavel, resultado in zip(incognitas, solução[0]):

            i = int(variavel.split('_')[1])
            nome_variavel = variavel.split('_')[0]

            x = self.fixações[i].x
            y = self.fixações[i].y

            mag_x = 0 
            mag_y = 0
            mag_z = 0

            if nome_variavel == 'x':

                mag_x = resultado

                self.esforços_externos += [CarregamentoPontual(x, y, mag_x, mag_y)]

            elif nome_variavel == 'y':

                mag_y = resultado

                self.esforços_externos += [CarregamentoPontual(x, y, mag_x, mag_y)]

            elif nome_variavel == 'M':

                mag_z = resultado

                self.esforços_externos += [Momento(x, y, mag_z)]

        return resultados

    def solução_esforços_internos(self):
        """
        Método responsável pela solução dos esforços internos aplicados nas fixações.
        """

        if self.esforços_externos == []:

            self.solução_esforços_externos()

        nos_ordenados = sorted(list(set(self.nos)), key=lambda tupla: (tupla[0], tupla[1])) 

        inicio = 0
        fim = 0

        for i in range(len(nos_ordenados)):

            if i != 0:

                fim += np.linalg.norm(np.asarray(nos_ordenados[i])-np.asarray(nos_ordenados[i-1]))

                t = np.arange(inicio, fim, 0.0001)
                x = ['x' for _ in range(len(t))]

                inicio = fim

                metade_x = (nos_ordenados[i][0] + nos_ordenados[i-1][0]) / 2
                metade_y = (nos_ordenados[i][1] + nos_ordenados[i-1][1]) / 2

                esforços_esquerda = [esforço for esforço in self.esforços_externos if (esforço.x <= metade_x and esforço.y <= metade_y)]
                esforços_esquerda += [carregamento for carregamento in self.carregamentos if (carregamento.x <= metade_x and carregamento.y <= metade_y)] 

                soma_esforços_cortantes = ''
                soma_momento_fletor = ''

                for esforço in esforços_esquerda:
                    
                    try:
                        
                        soma_esforços_cortantes = soma_esforços_cortantes + '+' + str(esforço.mag_y)
                        soma_momento_fletor = soma_momento_fletor + '-' + str(esforço.mag_y) + '*(x-' + str(esforço.x) + ')'

                    except:

                        soma_momento_fletor = soma_momento_fletor + '+' + str(esforço.mag_z) 
                
                equilibrio_y = soma_esforços_cortantes + '+x'
                equilibrio_M = soma_momento_fletor 

                expr_equilibrio_y = utilities.lambdify('x', sympify(equilibrio_y))
                expr_equilibrio_M = utilities.lambdify('x', sympify(equilibrio_M))
        
                self.esforços_internos += [EsforçoCortante(expr_equilibrio_y, t, sympify(equilibrio_y))]
                self.esforços_internos += [MomentoFletor(expr_equilibrio_M, t, sympify(equilibrio_M))]

    def solução_deflexões(self):
        """
        Método responsável pela solução das deflexões na barra.
        """

        if self.esforços_externos == []:

            self.solução_esforços_externos()

        if self.esforços_internos == []:

            self.solução_esforços_internos()

        momentos_fletores = [esforço_interno for esforço_interno in self.esforços_internos if isinstance(esforço_interno, MomentoFletor)]

        deflexões = [Function('f_'+str(i)) for i in range(len(momentos_fletores))]

        resultados = []
        expr_momentos_fletores = []
        condições_iniciais_equações = []

        flag_inicio = True 

        for i, esforço_interno in enumerate(momentos_fletores):

            f = deflexões[i]

            self.I_seção = 0.1

            expr_momento_fletor = Derivative(f('x'), ('x',2)) + (esforço_interno.expr/(self.E*self.I_seção))
            
            expr_momentos_fletores += [Eq(0, expr_momento_fletor)]

        resultados_equações = dsolve(expr_momentos_fletores)

        resultados_expr = [sympify(resultado.rhs) for resultado in resultados_equações]

        for i, esforço_interno in enumerate(momentos_fletores):

            f = resultados_expr[i]

            for fixação in self.fixações:

                if math.isclose(min(esforço_interno.intervalo), fixação.x, abs_tol=0.1) or math.isclose(max(esforço_interno.intervalo), fixação.x, abs_tol=0.1) or min(esforço_interno.intervalo) <= fixação.x <= max(esforço_interno.intervalo):

                    if fixação.fixa_y:

                        condições_iniciais_equações += [Eq(f.subs('x', fixação.x), 0)]
                    
                    if fixação.fixa_rot:

                        condições_iniciais_equações += [Eq(diff(f, 'x').subs('x',fixação.x), 0)]
            
            if not flag_inicio:
                
                f_anterior = resultados_expr[i-1]

                condições_iniciais_equações += [Eq(f.subs('x', min(esforço_interno.intervalo)), f_anterior.subs('x', min(esforço_interno.intervalo)))]

                condições_iniciais_equações += [Eq(diff(f, 'x').subs('x', min(esforço_interno.intervalo)), diff(f_anterior, 'x').subs('x', min(esforço_interno.intervalo)))]

            flag_inicio = False

        constantes = solve(condições_iniciais_equações)

        resultados_expr = [i.subs(constantes) for i in resultados_expr]
        resultados = [utilities.lambdify('x', sympify(resultado)) for resultado in resultados_expr]

        for i, resultado in enumerate(resultados): 
            
            intervalo = momentos_fletores[i].intervalo

            self.deflexões += [Deflexão(resultado, intervalo)]