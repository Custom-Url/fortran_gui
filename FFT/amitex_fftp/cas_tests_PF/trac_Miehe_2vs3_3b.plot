plot "../resultats_PF/SENT3/SENT.std" u 9:3,"../resultats_PF/SENT3b/SENT.std" u 9:3,"../resultats_PF/SENT2/SENT.std" u 9:3
pause -1 "'Return' for next simulation"

plot "../resultats_PF/SENT_VARgclc3/SENT_FIXgc_VARlc_2.std" u 8:2,"../resultats_PF/SENT_VARgclc3b/SENT_FIXgc_VARlc_2.std" u 8:2,"../resultats_PF/SENT_VARgclc2/SENT_FIXgc_VARlc_2.std" u 8:2
pause -1 "'Return' for next simulation"
plot "../resultats_PF/SENT_VARgclc3/SENT_FIXlc_VARgc_2.std" u 8:2,"../resultats_PF/SENT_VARgclc3b/SENT_FIXlc_VARgc_2.std" u 8:2,"../resultats_PF/SENT_VARgclc2/SENT_FIXlc_VARgc_2.std" u 8:2

pause -1 "'Return' for next simulation"
plot "../resultats_PF/shearCrack3/SC_xy.std" u 11:5,"../resultats_PF/shearCrack3b/SC_xy.std" u 11:5,"../resultats_PF/shearCrack2/SC_xy.std" u 11:5
pause -1 "'Return' for next simulation"
plot "../resultats_PF/shearCrack3/SC_xz.std" u 12:6,"../resultats_PF/shearCrack3b/SC_xz.std" u 12:6,"../resultats_PF/shearCrack2/SC_xz.std" u 12:6
pause -1 "'Return' for next simulation"
plot "../resultats_PF/shearCrack3/SC_yz.std" u 13:7,"../resultats_PF/shearCrack3b/SC_yz.std" u 13:7, "../resultats_PF/shearCrack2/SC_yz.std" u 13:7
pause -1 "'Return' for next simulation"

plot "../resultats_PF/DENT3/DENTasym.std" u 9:3,"../resultats_PF/DENT3b/DENTasym.std" u 9:3,"../resultats_PF/DENT2/DENTasym.std" u 9:3
pause -1 "'Return' for next simulation"
plot "../resultats_PF/DENT3/DENTsym.std" u 9:3,"../resultats_PF/DENT3b/DENTsym.std" u 9:3,"../resultats_PF/DENT2/DENTsym.std" u 9:3
pause -1 "'Return' for next simulation"

plot "../resultats_PF/bimat3/bimat.std" u 8:2,"../resultats_PF/bimat3b/bimat.std" u 8:2,"../resultats_PF/bimat2/bimat.std" u 8:2
pause -1 "Last result"


