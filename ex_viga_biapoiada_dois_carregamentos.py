# Importando as classes e métodos do arquivo estrutura.py
from bibstruct import *

# Definindo uma seção transversal retângular
ret = Retangulo(0.2, 0.1)

# Definindo uma barra de 5 metros com seção transversal retângular
barra = Barra(0, 0, 0, 5, ret, 2*10**11)

# Instânciando as fixações nas extremidades da viga
fix_1 = Fixação(0, 0, fixa_x=True, fixa_y=True)
fix_2 = Fixação(5, 0, fixa_y=True)

# Instânciando um carregamento pontual de 2KN na metade da viga
f_1 = CarregamentoPontual(2, 0, 0, -2000)
f_2 = CarregamentoPontual(3, 0, 0, 2000)

# Adicionando as fixações e carregamentos à barra
barra.adicionar_fixação(fix_1)
barra.adicionar_fixação(fix_2)
barra.adicionar_carregamento(f_1)
barra.adicionar_carregamento(f_2)

# Representação da estrutura construída
barra.esquema()

# Solucionando os esforços externos
barra.solução_esforços_externos()

# Representação da estrutura com os esforços externos
barra.esquema()

# Solucionando os esforços internos
barra.solução_esforços_internos()

# Representação da estrutura com os esforços internos
barra.esquema()

# Solucionando as deflexões
barra.solução_deflexões()

# Representação da estrutura com as deflexões
barra.esquema()