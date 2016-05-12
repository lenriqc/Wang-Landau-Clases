#!/bin/bash

#  loop_en_L_simulacion.sh
#  WL-CLASES
#
#  Created by Luis Enrique Quintanar Cortés on 19/02/16.
#
#!/bin/bash

for k in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
do
export numero_de_simulacion=$k
{ time ./sim.sh ; } 2>> tiempos_red_cuadrada_serie_para_promedios.txt
done