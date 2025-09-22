#!/bin/bash
if [ ! $(basename $(pwd)) = "validation" ]; then
   echo "ERROR : clean_result.sh must be launched in its own directory 'validation'"
   exit 0
fi
#----------------------------------------
PWD0=/u/q/dg765/amitex_fftp/resultats
cd $PWD0
PWD1=$(pwd)

if [ "$PWD0" = "$PWD1" ]
then
echo "Cleaning $PWD1"
for fich in ./*
do
  rm "$fich"/"$fich"*.std  2> /dev/null
  rm "$fich"/"$fich"*.mstd 2> /dev/null
  rm "$fich"/"$fich"*.zstd 2> /dev/null
  rm "$fich"/"$fich"*.log  2> /dev/null
  rm "$fich"/"$fich"*.xml  2> /dev/null
  rm "$fich"/"$fich"*.vtk  2> /dev/null
done
fi

