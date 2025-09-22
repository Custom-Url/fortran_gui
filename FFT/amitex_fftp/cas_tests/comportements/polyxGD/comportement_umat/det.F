C DET       SOURCE    CHAT      05/01/12    22:44:04     5004           
C DETMAT3   SOURCE    MAGN      00/05/16    21:19:14     3851           
C***********************************************************************
C*                                                                     *
C* PROJET      : Opérateur LIMI                                        *
C* NOM         : DETMAT3, fonction                                               *
C* DESCRIPTION : Cette fonction calcule le déterminant d'une matrice   *
C*               carrée de rang 3                                      *
C* LANGAGE     : Esope                                                 *
C* AUTEUR      : Guillaume VENCO - DRN/DMT/SEMT/LTMF                   *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C* APPELES     : Aucun                                                 *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C* ENTREES     : A   : -Matrice carrée 3*3                             *
C*                                                                     * 
C* SORTIES     : DET : -Déterminant de la matrice A                    *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C* VERSION     : 20/04/2000                                            *
C* CREATION    :                                                       *
C*                                                                     *
C***********************************************************************
       
      FUNCTION DETMAT3(A)
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C*       *********************************
C*       *** Déclaration des variables ***
C*       *********************************

        REAL*8 A(3,3),DETMAT3

C*       *********************************
C*       ******** Calcul de DETMAT3  *********
C*       *********************************

        DETMAT3= A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)
     &      +A(2,1)*A(3,2)*A(1,3)-A(3,1)*A(2,2)*A(1,3)
     &      -A(2,1)*A(1,2)*A(3,3)-A(1,1)*A(3,2)*A(2,3) 

        RETURN
        END

C*       *********************************
C*       ****  Fin Sous-Programme  *******
C*       *********************************
 
 
 
