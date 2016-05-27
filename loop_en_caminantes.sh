#!/bin/bash

#  loop_en_L_simulacion.sh
#  WL-CLASES
#
#  Created by Luis Enrique Quintanar Cortés on 19/02/16.
#
#!/bin/bash
export max_bloques=9000
for k in {1,8,10,16,32,50,64,128,200,256}
do
export caminantes=$k
for l in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
do
export numero_de_simulacion=$l
{ echo "$k caminantes----simulación $l" ; time ./sim.sh ; } 2>> tiempos_red_cuadrada_bloques_$max_bloques_repeticiones_de_caminantes.txt
done
done