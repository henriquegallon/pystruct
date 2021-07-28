# Importando as classes e métodos do arquivo estrutura.py
from bibstruct import *

# Definindo uma seção transversal retângular
ret = Retangulo(0.2, 0.1)

# Definindo uma barra de 5 metros com seção transversal retângular
barra = Barra(0, 0, 0, 5, ret)

# Instânciando a fixação na extremidade da viga
fix_1 = Fixação(0, 0, fixa_x=True, fixa_y=True, fixa_rot=True)

# Instânciando um carregamento pontual de 2KN na extremidade da viga
f_1 = CarregamentoPontual(5, 0, 0, -2000)

# Adicionando as fixações e carregamentos à barra
barra.adicionar_fixação(fix_1)
barra.adicionar_carregamento(f_1)

# Representação da estrutura construída
barra.esquema()

# Solucionando os esforços externos
barra.solução_esforços_externos()

# Representação da estrutura com os esforços externos
barra.esquema()

# Representação da estrutura com os esforços internos
barra.solução_esforços_internos()