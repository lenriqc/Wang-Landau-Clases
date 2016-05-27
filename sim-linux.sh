#!/bin/bash

#  sim.sh
#  
#
#  Created by Luis Enrique Quintanar Cortés on 13/01/15.
#
#!/bin/bash
#./simulacion 2 2 2 20 0 1
./simulacion 2 2 3 50 2 0.0001 10 2 1

# primer argumento: tipo de red: 1-cuad 2-triang
# Segundo argumento: Tipo de modelo: 1-Ising    2-Antiferro
# Tercer argumento: Tipo de algoritmo: 1-Metrop 2-WL1D  3-WL2D
# Cuarto argumento: Longitud de la red por lado
# Quinto argumento: Campo externo
# Sexto argumento: Límite de factor de modificación
# Séptimo argumento: número de caminantes
# Octavo argumento: 1-animar-histograma 2-mostrar observables
# Noveno argumento: Número de sweeps por cada histograma

#2 tipo de red