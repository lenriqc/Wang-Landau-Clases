#!/bin/bash

#  loop_en_L_simulacion.sh
#  WL-CLASES
#
#  Created by Luis Enrique Quintanar Cortés on 19/02/16.
#
#!/bin/bash

for k in {20,50,100}
do
export L=$k
{ time ./sim.sh ; } 2>> tiempos_red_cuadrada_variando_L.txt
done
