#!/bin/sh

#  sim.sh
#  
#
#  Created by Luis Enrique Quintanar Cortés on 13/01/15.
#
#!/bin/bash

red=${red:-triangular}
modelo=${modelo:-Ising-antiferro}
algoritmo=${algoritmo:-WL2D}
L=${L:-18}
h=${h:-0}
f=${f:-0.0001}
caminantes=${caminantes:-30}
max_bloques=${max_bloques:-no}
frec_animacion=${frec_animacion:-no}
numero_de_simulacion=${numero_de_simulacion:-no}

echo $caminantes
echo $max_bloques
echo "Frecuencia de animación = $frec_animacion"

./simulacion $red $modelo $algoritmo $L $h $f $caminantes $max_bloques  $frec_animacion $numero_de_simulacion

# primer argumento: tipo de red: 1-cuad 2-triang
# Segundo argumento: Tipo de modelo: 1-Ising    2-Antiferro
# Tercer argumento: Tipo de algoritmo: 1-Metrop 2-WL1D  3-WL2D
# Cuarto argumento: Longitud de la red por lado
# Quinto argumento: Campo externo
# Sexto argumento: Límite de factor de modificación
# Séptimo argumento: Número de caminantes
# Octavo argumento: Máximo número de bloques de sweeps
# Noveno argumento: Frecuencia de graficado. En cantidad de series de sweeps
# Décimo argumento: Número de repetición de la simulación

