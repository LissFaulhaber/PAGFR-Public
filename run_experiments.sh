#!/bin/bash

echo "--------------ModelBellmans-------------------">>modelBell_50.txt
for instancia in $(find data_new/ -name "*_50.txt" | sort -t '/' -k2.5)
do
	bin/EVCSLP $instancia 2 >> modelBell_50.txt 
done
