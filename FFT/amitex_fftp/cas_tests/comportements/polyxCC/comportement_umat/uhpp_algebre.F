C***********************************************************************
C* Calcul du produit matriciel C=A.B (tenseurs d'ordre 2)              *
C***********************************************************************
      SUBROUTINE PROD (A,B, C)

      INCLUDE 'IMPLICIT.INC'

      DIMENSION A(*),B(*), C(*)

      C(1) = A(1)*B(1) + A(4)*B(2) + A(7)*B(3)
      C(2) = A(2)*B(1) + A(5)*B(2) + A(8)*B(3)
      C(3) = A(3)*B(1) + A(6)*B(2) + A(9)*B(3)
      C(4) = A(1)*B(4) + A(4)*B(5) + A(7)*B(6)
      C(5) = A(2)*B(4) + A(5)*B(5) + A(8)*B(6)
      C(6) = A(3)*B(4) + A(6)*B(5) + A(9)*B(6)
      C(7) = A(1)*B(7) + A(4)*B(8) + A(7)*B(9)
      C(8) = A(2)*B(7) + A(5)*B(8) + A(8)*B(9)
      C(9) = A(3)*B(7) + A(6)*B(8) + A(9)*B(9)

      RETURN
      END

C***********************************************************************
C* Calcul du produit matriciel C=A^T.B (tenseurs d'ordre 2)              *
C***********************************************************************
      SUBROUTINE PRODT (A,B, C)

      INCLUDE 'IMPLICIT.INC'

      DIMENSION A(*),B(*), C(*)

      C(1) = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      C(2) = A(4)*B(1) + A(5)*B(2) + A(6)*B(3)
      C(3) = A(7)*B(1) + A(8)*B(2) + A(9)*B(3)
      C(4) = A(1)*B(4) + A(2)*B(5) + A(3)*B(6)
      C(5) = A(4)*B(4) + A(5)*B(5) + A(6)*B(6)
      C(6) = A(7)*B(4) + A(8)*B(5) + A(9)*B(6)
      C(7) = A(1)*B(7) + A(2)*B(8) + A(3)*B(9)
      C(8) = A(4)*B(7) + A(5)*B(8) + A(6)*B(9)
      C(9) = A(7)*B(7) + A(8)*B(8) + A(9)*B(9)

      RETURN
      END

C***********************************************************************
C  Calcul du produit symetrique     C= 1/2((Btrans)*(Atrans)+A*B)
C***********************************************************************
      subroutine PRODS(A,B,C)
      INTEGER I,J,K
      REAL*8  A(3,3),B(3,3),C(3,3),D(3,3)
c
      DO J=1,3
        DO I=1,3
          C(I,J)=0.0D0
          D(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
            D(I,J)=D(I,J)+B(K,I)*A(J,K)
          ENDDO
          C(I,J)=0.5D0*(C(I,J)+D(I,J))
        ENDDO
      ENDDO
      END


C***********************************************************************
C* Calcul du determinant (DET) de la matrice A                         *
C***********************************************************************
      FUNCTION DET (A)

      INCLUDE 'IMPLICIT.INC'

      DIMENSION A(*)

      A11 = A(1)
      A21 = A(2)
      A31 = A(3)
      A12 = A(4)
      A22 = A(5)
      A32 = A(6)
      A13 = A(7)
      A23 = A(8)
      A33 = A(9)
      W1 = A22*A33 - A32*A23
      W2 = A31*A23 - A21*A33
      W3 = A21*A32 - A31*A22
      DET = 1. / (A11*W1 + A12*W2 + A13*W3)

      RETURN
      END

C***********************************************************************
C* Calcul du determinant (ADET) et de l'inverse (AINV) de la matrice A *
C***********************************************************************
      SUBROUTINE INV (A, AINV,ADET)

      INCLUDE 'IMPLICIT.INC'

      DIMENSION A(*), AINV(*)

      A11 = A(1)
      A21 = A(2)
      A31 = A(3)
      A12 = A(4)
      A22 = A(5)
      A32 = A(6)
      A13 = A(7)
      A23 = A(8)
      A33 = A(9)
      W1 = A22*A33 - A32*A23
      W2 = A31*A23 - A21*A33
      W3 = A21*A32 - A31*A22
      ADET = 1. / (A11*W1 + A12*W2 + A13*W3)
      AINV(1) = ADET * W1
      AINV(2) = ADET * W2
      AINV(3) = ADET * W3
      AINV(4) = ADET * (A32*A13 - A12*A33)
      AINV(5) = ADET * (A11*A33 - A31*A13)
      AINV(6) = ADET * (A31*A12 - A11*A32)
      AINV(7) = ADET * (A12*A23 - A22*A13)
      AINV(8) = ADET * (A21*A13 - A11*A23)
      AINV(9) = ADET * (A11*A22 - A21*A12)

      RETURN
      END

C-----------------------------------------------------------------------
C- Decomposition en parties symetrique et antisymetrique
C-----------------------------------------------------------------------
      SUBROUTINE SYMANT (TENS, SYME,ANTI)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION TENS(*), SYME(*),ANTI(*)

C- Tenseur SYMEtrique : SYME = (TENS+transpose(TENS))/2
      Z12 = c1s2 * (TENS(2) + TENS(4))
      Z13 = c1s2 * (TENS(3) + TENS(7))
      Z23 = c1s2 * (TENS(6) + TENS(8))
      SYME(1) = TENS(1)
      SYME(2) = Z12
      SYME(3) = Z13
      SYME(4) = Z12
      SYME(5) = TENS(5)
      SYME(6) = Z23
      SYME(7) = Z13
      SYME(8) = Z23
      SYME(9) = TENS(9)

C- Tenseur ANTIsymetrique : ANTI = (TENS-transpose(TENS))/2
      Z12 = c1s2 * (TENS(4) - TENS(2))
      Z13 = c1s2 * (TENS(7) - TENS(3))
      Z23 = c1s2 * (TENS(8) - TENS(6))
      ANTI(1) = cZero
      ANTI(2) = -Z12
      ANTI(3) = -Z13
      ANTI(4) = Z12
      ANTI(5) = cZero
      ANTI(6) = -Z23
      ANTI(7) = Z13
      ANTI(8) = Z23
      ANTI(9) = cZero

      RETURN
      END
      
C-----------------------------------------------------------------------
C- Calcul de la partie elastique du gradient de la transformation FE
C- a partir de F et FP
C-----------------------------------------------------------------------
      SUBROUTINE FEFFP (F, FP, FE)

      INCLUDE 'IMPLICIT.INC'

      DIMENSION F(*), FP(*),FE(*),FPINV(9)
      call INV(FP,FPINV,DETFP)
      call PROD(F,FPINV,FE)
      END


C************************************************
C  Calcule EEL a partir de FE
C***********************************************
      SUBROUTINE EELFE(FE,EEL1)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION FE(*),EEL1(*)

      call PRODT(FE,FE,EEL1)
      EEL1(1)=EEL1(1)-cUn
      EEL1(5)=EEL1(5)-cUn
      EEL1(9)=EEL1(9)-cUn
      do i=1,9
       EEL1(i)=cUn/cDeux*EEL1(i)
      enddo
c
      END
  

C************************************************
C  Calcule contrainte dans la base micro
C***********************************************

      SUBROUTINE ELASTIC(DEF,C11,C12,C44,SIG)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION DEF(*),SIG(*)
c
      SIG(1)=C11*DEF(1)+C12*(DEF(5)+DEF(9))
      SIG(5)=C11*DEF(5)+C12*(DEF(1)+DEF(9))
      SIG(9)=C11*DEF(9)+C12*(DEF(1)+DEF(5))
      SIG(2)=cDeux*C44*DEF(2)
      SIG(3)=cDeux*C44*DEF(3)
      SIG(4)=cDeux*C44*DEF(4)
      SIG(6)=cDeux*C44*DEF(6)
      SIG(7)=cDeux*C44*DEF(7)
      SIG(8)=cDeux*C44*DEF(8)
      END

C************************************************
C ROTATECC(TCH,R,T)  T=R*TCH*(Rtrans)
C************************************************
      SUBROUTINE ROTATECC(TCH,R,T)
      INTEGER I,J,K,L
      REAL*8  S,RJL,A(3,3),T(3,3),R(3,3),TCH(3,3)
C
      DO 210 J=1,3
         DO 200 I=1,3
            S=0.0
            DO 110 L=1,3
            RJL=R(J,L)
               DO 100 K=1,3
                  S=S+R(I,K)*RJL*TCH(K,L)
100           CONTINUE
110         CONTINUE
            A(I,J)=S
200     CONTINUE
210   CONTINUE
C
      DO I=1,3
       DO J=1,3
        T(I,J)=A(I,J)
       ENDDO
      ENDDO
      END 

c***************************************
C ROTINVCC(T,R,TCH)  TCH=(Rtrans)*T*R
c***************************************
      SUBROUTINE ROTINVCC(T,R,TCH)
      INTEGER  I,J,K,L
      REAL*8   S,RLJ,A(3,3),T(3,3),R(3,3),TCH(3,3)
C
      DO 210 J=1,3
         DO 200 I=1,3
            S=0.0
            DO 110 L=1,3
               RLJ=R(L,J)
               DO 100 K=1,3
                  S=S+R(K,I)*RLJ*T(K,L)
100           CONTINUE
110        CONTINUE
            A(I,J)=S
200     CONTINUE
210   CONTINUE
C
      DO I=1,3
       DO J=1,3
        TCH(I,J)=A(I,J)
       ENDDO
      ENDDO
      END 

c***************************************
C Initialisation de la matrice de Taylor nimis
c***************************************
       SUBROUTINE CGNS(GNS)
       
      INCLUDE 'IMPLICIT.INC'
       
       PARAMETER(MAXS=12,MAXS3=36)
       DIMENSION XNSYS(MAXS3),XMSYS(MAXS3),GNS(3,3,MAXS)
       INCLUDE 'nsms.inc'

c-------- cree les nismis 
c-------- 
      DO IS=1,MAXS
       DO I=1,3
        DO J=1,3
         GNS(I,J,IS)= (XMSYS(IS+(I-1)*MAXS)*XNSYS(IS+(J-1)*MAXS))   
        END DO
       END DO         
      END DO
  
      END

c***************************************
C Initialisation de la matrice de Taylor nimis symetrique
c***************************************
       SUBROUTINE CGNSSYM(GNSSYM)
       
      INCLUDE 'IMPLICIT.INC'
       
       PARAMETER(MAXS=12,MAXS3=36)
       DIMENSION XNSYS(MAXS3),XMSYS(MAXS3),GNSSYM(3,3,MAXS)
       INCLUDE 'nsms.inc'

c-------- cree les nismis 
c-------- 
      DO IS=1,MAXS
       DO I=1,3
        DO J=1,3
         GNSSYM(I,J,IS)= 0.5D0*
     & ((XMSYS(IS+(I-1)*MAXS)*XNSYS(IS+(J-1)*MAXS)) + 
     & (XNSYS(IS+(I-1)*MAXS)*XMSYS(IS+(J-1)*MAXS)))
        END DO
       END DO         
      END DO
  
      END

c***************************************
c CISRED(TAU,DSB,TAUSB) calcule la cission reduite sur tous les systemes
c***************************************
      SUBROUTINE  CISRED(TAU,DSB,TAUSB)
      INTEGER   MAXS,ISY,I,J
      PARAMETER (MAXS=12)
      REAL*8 TAUSB(MAXS)
      REAL*8 TAU(3,3), DSB(3,3,MAXS)
C
      DO 300 ISY=1,MAXS
         TAUSB(ISY)=0.0
         DO 200 I=1,3
            DO 100 J=1,3
               TAUSB(ISY)=TAUSB(ISY)+DSB(I,J,ISY)*TAU(I,J)
100        CONTINUE
200       CONTINUE
300   CONTINUE
      END 
      
c***************************************
c VECTOMAT(STRAN,DEFOR) Transforme un vecteur à 6 composantes en matrice 3x3
c***************************************
      SUBROUTINE  VECTOMAT(STRAN,DEFOR)
      REAL*8 STRAN(6),DEFOR(3,3)
      DEFOR(1,1)=STRAN(1)
      DEFOR(2,2)=STRAN(2)
      DEFOR(3,3)=STRAN(3)
      DEFOR(2,1)=STRAN(4)
      DEFOR(3,1)=STRAN(5)
      DEFOR(3,2)=STRAN(6)
      DEFOR(1,2)=DEFOR(2,1)
      DEFOR(1,3)=DEFOR(3,1)
      DEFOR(2,3)=DEFOR(3,2)
      END
c***************************************
c MATTOVEC(DEFOR,STRAN) Transforme une matrice sym 3x3 en vecteur à 6 composantes 
c***************************************
      SUBROUTINE  MATTOVEC(DEFOR,STRAN)
      REAL*8 STRAN(6),DEFOR(3,3)
      STRAN(1)=DEFOR(1,1)
      STRAN(2)=DEFOR(2,2)
      STRAN(3)=DEFOR(3,3)
      STRAN(4)=DEFOR(2,1)
      STRAN(5)=DEFOR(3,1)
      STRAN(6)=DEFOR(3,2)
      END

