#!/bin/tcsh

convert +append -gravity center PAU_CMD_LePhare.png GAMA_PAU_nz.png GAMA_CMD.png fig1.png
convert -append -gravity center fig1.png PAU_CMD_Cigale.png PAU_GAMA_zCMD.png 
rm fig1.png

