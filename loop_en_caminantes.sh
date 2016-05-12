#!/bin/bash

#  loop_en_L_simulacion.sh
#  WL-CLASES
#
#  Created by Luis Enrique Quintanar Cortés on 19/02/16.
#
#!/bin/bash

for k in {1,8,10,16,32,50,64,128,200,256}
do
export caminantes=$k
export max_bloques=40000
{ time ./sim.sh ; } 2>> tiempos_red_cuadrada_bloques_$max_bloques.txt
done