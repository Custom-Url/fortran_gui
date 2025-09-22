#!/bin/bash
if [ ! $(basename $(pwd)) = "cas_tests_PF" ]; then
   echo "ERROR : clean_result_PF.sh must be launched in its own directory 'cas_tests_PF'"
   exit 0
fi
#----------------------------------------
PWD0=/u/q/dg765/amitex_fftp/resultats_PF
cd $PWD0
PWD1=$(pwd)

if [ "$PWD0" = "$PWD1" ]
then
echo "Cleaning $PWD1"
for fich in ./*
do
  rm -rf * 
done
touch empty_file #provisoire (git refuse de suivre un repertoire vide)
fi

