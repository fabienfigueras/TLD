#!/bin/bash

Echo "PBS Abacus Creation"

Str1="python3 ./PBS-v118.py 0.224 43 800 "
Str2=" 200 N 0 6 N 0.0001 N 0 N G1 1 Y"
Dist=100

echo "#!/bin/bash" >./PBS_Abacus_Creation.bsh

while [ $Dist -le 1600 ]
do

Cmd="$Str1$(printf "%d" $Dist)$Str2"
echo $Cmd
echo $Cmd >>./PBS_Abacus_Creation.bsh
echo "sleep 5" >>./PBS_Abacus_Creation.bsh
echo "wait" >>./PBS_Abacus_Creation.bsh

Dist=$[$Dist+100]
done
