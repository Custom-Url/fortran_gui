C***********************************************************************
C* Calcul de la "cinematique" en petites deformations sur un pas    *
C***********************************************************************
C* Remarque importante :                                               *
C* ---------------------                                               *
C* Dans les sous-programmes ci-dessous, les tenseurs                   *
C* sont associes a des matrices (3x3) stockees sous forme de           *
C* vecteurs a 9 composantes selon chaque colonne de la matrice :       *
C*  Vec(1) <=> MAT(1,1)    Vec(2) <=> MAT(2,1)    Vec(3) <=> MAT(3,1)  *
C*  Vec(4) <=> MAT(1,2)    Vec(5) <=> MAT(2,2)    Vec(6) <=> MAT(3,2)  *
C*  Vec(7) <=> MAT(1,3)    Vec(8) <=> MAT(2,3)    Vec(9) <=> MAT(3,3)  *
C***********************************************************************

C-----------------------------------------------------------------------
C- Increment du tenseur de deformation : DEPS
C- DEPS = EPS1 - EPS0  suppose constant sur le pas de temps
C-----------------------------------------------------------------------
      SUBROUTINE UHPP_CINE0 (EPS0,EPS1, DEPS)

      INCLUDE 'IMPLICIT.INC'

      DIMENSION EPS0(*),EPS1(*), DEPS(*)

      DEPS(1) = EPS1(1) - EPS0(1)
      DEPS(2) = EPS1(2) - EPS0(2)
      DEPS(3) = EPS1(3) - EPS0(3)
      DEPS(4) = EPS1(4) - EPS0(4)
      DEPS(5) = EPS1(5) - EPS0(5)
      DEPS(6) = EPS1(6) - EPS0(6)
      DEPS(7) = EPS1(7) - EPS0(7)
      DEPS(8) = EPS1(8) - EPS0(8)
      DEPS(9) = EPS1(9) - EPS0(9)

      RETURN
      END

C-----------------------------------------------------------------------
C- Vitesse de deformation : DEPSDT
C- DEPSDT = DEPS / DTPS  constante sur le pas de temps
C-----------------------------------------------------------------------
      SUBROUTINE UHPP_CINE1 (DEPS,DTPS, DEPSDT)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION DEPS(*), DEPSDT(*)

      UNsDT = cUn / DTPS

      DEPSDT(1) = DEPS(1) * UNsDT
      DEPSDT(2) = DEPS(2) * UNsDT
      DEPSDT(3) = DEPS(3) * UNsDT
      DEPSDT(4) = DEPS(4) * UNsDT
      DEPSDT(5) = DEPS(5) * UNsDT
      DEPSDT(6) = DEPS(6) * UNsDT
      DEPSDT(7) = DEPS(7) * UNsDT
      DEPSDT(8) = DEPS(8) * UNsDT
      DEPSDT(9) = DEPS(9) * UNsDT

      RETURN
      END

C-----------------------------------------------------------------------
C- Tenseur de deformation au cours du pas de temps : DEFOI
C- DEFOI = (1-ALPI).EPS0 + ALPI.EPS1 avec 0 <= ALPI <= 1
C-----------------------------------------------------------------------
      SUBROUTINE UHPP_CINE2 (EPS0,EPS1,DEPS, ALPI, DEFOI)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION EPS0(*),EPS1(*),DEPS(*), DEFOI(*)

      IF (ALPI.GE.cUn) THEN
        DEFOI(1) = EPS1(1)
        DEFOI(2) = EPS1(2)
        DEFOI(3) = EPS1(3)
        DEFOI(4) = EPS1(4)
        DEFOI(5) = EPS1(5)
        DEFOI(6) = EPS1(6)
        DEFOI(7) = EPS1(7)
        DEFOI(8) = EPS1(8)
        DEFOI(9) = EPS1(9)
      ELSE IF (ALPI.LE.cZero) THEN
        DEFOI(1) = EPS0(1)
        DEFOI(2) = EPS0(2)
        DEFOI(3) = EPS0(3)
        DEFOI(4) = EPS0(4)
        DEFOI(5) = EPS0(5)
        DEFOI(6) = EPS0(6)
        DEFOI(7) = EPS0(7)
        DEFOI(8) = EPS0(8)
        DEFOI(9) = EPS0(9)
      ELSE IF (ALPI.EQ.c1s2) THEN
        DEFOI(1) = c1s2 * (EPS1(1) + EPS0(1))
        DEFOI(2) = c1s2 * (EPS1(2) + EPS0(2))
        DEFOI(3) = c1s2 * (EPS1(3) + EPS0(3))
        DEFOI(4) = c1s2 * (EPS1(4) + EPS0(4))
        DEFOI(5) = c1s2 * (EPS1(5) + EPS0(5))
        DEFOI(6) = c1s2 * (EPS1(6) + EPS0(6))
        DEFOI(7) = c1s2 * (EPS1(7) + EPS0(7))
        DEFOI(8) = c1s2 * (EPS1(8) + EPS0(8))
        DEFOI(9) = c1s2 * (EPS1(9) + EPS0(9))
      ELSE
        DEFOI(1) = EPS0(1) + ALPI * DEPS(1)
        DEFOI(2) = EPS0(2) + ALPI * DEPS(2)
        DEFOI(3) = EPS0(3) + ALPI * DEPS(3)
        DEFOI(4) = EPS0(4) + ALPI * DEPS(4)
        DEFOI(5) = EPS0(5) + ALPI * DEPS(5)
        DEFOI(6) = EPS0(6) + ALPI * DEPS(6)
        DEFOI(7) = EPS0(7) + ALPI * DEPS(7)
        DEFOI(8) = EPS0(8) + ALPI * DEPS(8)
        DEFOI(9) = EPS0(9) + ALPI * DEPS(9)
      ENDIF

      RETURN
      END

C-----------------------------------------------------------------------
C- Tenseur des deformations au cours du pas de temps : DEFOI
C- DEFOI = (1-ALPI).EPS0 + ALPI.EPS1 avec 0 <= ALPI <= 1
C-       = EPS0 + ALPI.DEPS
C-----------------------------------------------------------------------
      SUBROUTINE UHPP_CINE3 (EPS0,EPS1,DEPS, ALPI, DEFOI)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION EPS0(*),EPS1(*),DEPS(*), DEFOI(*)

C- Calcul de DEFOI pour l'instant ALPI
      CALL UHPP_CINE2(EPS0,EPS1,DEPS, ALPI, DEFOI)

      RETURN
      END

