#!/bin/bash
  cd /u/q/dg765/amitex_fftp/resultats

for fich in ./*
do
  mv "$fich"/"$fich".std "$fich"/reference 2> /dev/null
  mv "$fich"/"$fich".log "$fich"/reference 2> /dev/null
done

