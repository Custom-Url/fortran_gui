C***********************************************************************
C* Calcul de la "cinematique" en grandes transformations sur un pas    *
C***********************************************************************
C* Remarque importante :                                               *
C* ---------------------                                               *
C* Dans les sous-programmes ci-dessous, les gradients de la transforma-*
C* tion sont associes a des matrices (3x3) stockees sous forme de      *
C* vecteurs a 9 composantes selon chaque colonne de la matrice :       *
C*  Vec(1) <=> MAT(1,1)    Vec(2) <=> MAT(2,1)    Vec(3) <=> MAT(3,1)  *
C*  Vec(4) <=> MAT(1,2)    Vec(5) <=> MAT(2,2)    Vec(6) <=> MAT(3,2)  *
C*  Vec(7) <=> MAT(1,3)    Vec(8) <=> MAT(2,3)    Vec(9) <=> MAT(3,3)  *
C***********************************************************************

C-----------------------------------------------------------------------
C- Increment du gradient de la transformation : DGRDF
C- DGRDF = GRDF1 - GRDF0  suppose constant sur le pas de temps
C-----------------------------------------------------------------------
      SUBROUTINE UGD_CINE0 (GRDF0,GRDF1, DGRDF)

      include 'IMPLICIT.INC'

      DIMENSION GRDF0(*),GRDF1(*), DGRDF(*)

      DGRDF(1) = GRDF1(1) - GRDF0(1)
      DGRDF(2) = GRDF1(2) - GRDF0(2)
      DGRDF(3) = GRDF1(3) - GRDF0(3)
      DGRDF(4) = GRDF1(4) - GRDF0(4)
      DGRDF(5) = GRDF1(5) - GRDF0(5)
      DGRDF(6) = GRDF1(6) - GRDF0(6)
      DGRDF(7) = GRDF1(7) - GRDF0(7)
      DGRDF(8) = GRDF1(8) - GRDF0(8)
      DGRDF(9) = GRDF1(9) - GRDF0(9)

      RETURN
      END

C-----------------------------------------------------------------------
C- Vitesse du gradient de la transformation : VGRDF
C- VGRDF = DGRDF / DTPS  constante sur le pas de temps
C-----------------------------------------------------------------------
      SUBROUTINE UGD_CINE1 (DGRDF,DTPS, VGRDF)

      include 'IMPLICIT.INC'
      include 'CMATH.INC'

      DIMENSION DGRDF(*), VGRDF(*)

      UNsDT = cUn / DTPS

      VGRDF(1) = DGRDF(1) * UNsDT
      VGRDF(2) = DGRDF(2) * UNsDT
      VGRDF(3) = DGRDF(3) * UNsDT
      VGRDF(4) = DGRDF(4) * UNsDT
      VGRDF(5) = DGRDF(5) * UNsDT
      VGRDF(6) = DGRDF(6) * UNsDT
      VGRDF(7) = DGRDF(7) * UNsDT
      VGRDF(8) = DGRDF(8) * UNsDT
      VGRDF(9) = DGRDF(9) * UNsDT

      RETURN
      END

C-----------------------------------------------------------------------
C- Gradient de la transformation au cours du pas de temps : GRDFI
C- GRDFI = (1-ALPI).GRDF0 + ALPI.GRDF1 avec 0 <= ALPI <= 1
C-----------------------------------------------------------------------
      SUBROUTINE UGD_CINE2 (GRDF0,GRDF1,DGRDF, ALPI, GRDFI)

      include 'IMPLICIT.INC'
      include 'CMATH.INC'

      DIMENSION GRDF0(*),GRDF1(*),DGRDF(*), GRDFI(*)

      IF (ALPI.GE.cUn) THEN
        GRDFI(1) = GRDF1(1)
        GRDFI(2) = GRDF1(2)
        GRDFI(3) = GRDF1(3)
        GRDFI(4) = GRDF1(4)
        GRDFI(5) = GRDF1(5)
        GRDFI(6) = GRDF1(6)
        GRDFI(7) = GRDF1(7)
        GRDFI(8) = GRDF1(8)
        GRDFI(9) = GRDF1(9)
      ELSE IF (ALPI.LE.cZero) THEN
        GRDFI(1) = GRDF0(1)
        GRDFI(2) = GRDF0(2)
        GRDFI(3) = GRDF0(3)
        GRDFI(4) = GRDF0(4)
        GRDFI(5) = GRDF0(5)
        GRDFI(6) = GRDF0(6)
        GRDFI(7) = GRDF0(7)
        GRDFI(8) = GRDF0(8)
        GRDFI(9) = GRDF0(9)
      ELSE IF (ALPI.EQ.c1s2) THEN
        GRDFI(1) = c1s2 * (GRDF1(1) + GRDF0(1))
        GRDFI(2) = c1s2 * (GRDF1(2) + GRDF0(2))
        GRDFI(3) = c1s2 * (GRDF1(3) + GRDF0(3))
        GRDFI(4) = c1s2 * (GRDF1(4) + GRDF0(4))
        GRDFI(5) = c1s2 * (GRDF1(5) + GRDF0(5))
        GRDFI(6) = c1s2 * (GRDF1(6) + GRDF0(6))
        GRDFI(7) = c1s2 * (GRDF1(7) + GRDF0(7))
        GRDFI(8) = c1s2 * (GRDF1(8) + GRDF0(8))
        GRDFI(9) = c1s2 * (GRDF1(9) + GRDF0(9))
      ELSE
        GRDFI(1) = GRDF0(1) + ALPI * DGRDF(1)
        GRDFI(2) = GRDF0(2) + ALPI * DGRDF(2)
        GRDFI(3) = GRDF0(3) + ALPI * DGRDF(3)
        GRDFI(4) = GRDF0(4) + ALPI * DGRDF(4)
        GRDFI(5) = GRDF0(5) + ALPI * DGRDF(5)
        GRDFI(6) = GRDF0(6) + ALPI * DGRDF(6)
        GRDFI(7) = GRDF0(7) + ALPI * DGRDF(7)
        GRDFI(8) = GRDF0(8) + ALPI * DGRDF(8)
        GRDFI(9) = GRDF0(9) + ALPI * DGRDF(9)
      ENDIF

      RETURN
      END

C-----------------------------------------------------------------------
C- Gradient de la transformation au cours du pas de temps : GRDFI
C- GRDFI = (1-ALPI).GRDF0 + ALPI.GRDF1 avec 0 <= ALPI <= 1
C-       = GRDF0 + ALPI.DGRDF
C- Taux de la transformation au cours du pas de temps : DLI
C- DLI = DGRDF . Inv(GRDFI)
C-----------------------------------------------------------------------
      SUBROUTINE UGD_CINE3 (GRDF0,GRDF1,DGRDF, ALPI, GRDFI, DLI)

      include 'IMPLICIT.INC'
      include 'CMATH.INC'

      DIMENSION GRDF0(*),GRDF1(*),DGRDF(*), GRDFI(*),DLI(*)
      DIMENSION dum1(9)

C- Calcul de GRDFI pour l'instant ALPI
      CALL UGD_CINE2(GRDF0,GRDF1,DGRDF, ALPI, GRDFI)
C- Inverse de GRDFI
      CALL INV(GRDFI, dum1,detFI)
c      CALL INV(GRDFI, dum1)
      IF (detFI.EQ.0) STOP 'Stop dans le calcul de 1/GRDFI'
C- Calcul L : LI = DGRDF.Inverse(GRDFI)
      CALL PROD33(DGRDF,dum1, DLI)

      RETURN
      END

