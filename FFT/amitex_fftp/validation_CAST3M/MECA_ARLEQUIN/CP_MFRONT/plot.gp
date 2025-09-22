
# COMPARAISON CUB8RI / FFT - arlequin 3
plot 'TYPE_3_1.txt_cub8ri' u 3:6, 'res3/res.std' u 25:4
pause -1 "'Return' for next simulation"

# COMPARAISON CUB8 / CUB8RI / FFT - arlequin 9
plot 'TYPE_9_1.txt_cub8ri' u 3:6,'TYPE_9_1.txt' u 3:6, 'res9/res.std' u 25:4
pause -1 "'Return' for next simulation"

# COMPARAISON arlequin 9 / 10 (CUB8 / FFT)  
plot 'TYPE_9_1.txt' u 3:6, 'TYPE_10_1.txt' u 3:6, 'res9/res.std' u 25:4, 'res10/res.std' u 25:4
pause -1 "'Return'"




