c=======================================================================
C=  Prepare les donnes pour l'appel a la procedure d'integration numerique
C=  Notamment, certaines variables definies dans les fichiers .dgibi
C=  ne sont pas integrees au moyen d'une equation differentielle.
C= Ces variables ne doivent donc pas faire partie des variables internes
C= definies pour la routine rk43 ou rk21
C=======================================================================
C=                                                                     =
C=  Entrees :                                                          =
C=  ---------                                                          =
C=  nLdC     Numero du modele de comportement (idem que dans UMAT)     =
C=  iCALC    Indicateur du mode de calcul (cf. IFOMOD)                 =
C=  TIME     Instant initial "t"                                       =
C=  DTIME    Pas de temps "Dt"                                         =
C=  TEMP     Temperature T(t) au debut du pas de temps "t"             =
C=  DTEMP    Increment de temperature "DT" sur le pas de temps "Dt"    =
C=  STRAN    Tenseur des deformations a l'instant "t"                  =
C=  DSTRAN   Increment de deformations entre l'instant "t" et "t+Dt"   =
C=  PROPS    Proprietes materielles a l'instant "t+Dt"                 =
C=  NPROPS   Nombre de coefficients materiau (PROPS(1:NPROPS))         =
C=  PREDEF   Parametres externes a l'instant "t"                       =
C=  DPRED    Variation de ces parametres sur le pas de temps "Dt"      =
C=  STATEV   Variables internes du modele au debut du pas de temps "t" =
C=  NSTATEV  Nombre de variables internes (STATEV(1:NSTATV))           =
C=  STRESS   Contraintes de Cauchy au debut du pas de temps "t"        =
C=  NTENS    Nombre de composantes de contraintes (STRESS(1:NTENS))    =
C=                                                                     =
C=  Sorties :                                                          =
C=  ---------                                                          =
C=  STATEV   Variables internes calculees a la fin du pas de temps     =
C=  STRESS   Contraintes de Cauchy calculees a la fin du pas de temps  =
C=  KINC     Indicateur d'erreur (si KINC different de 1)              =
C=                                                                     =
C=  Remarque importante :                                              =
C=  ---------------------                                              =
C=  Le vecteur PROPS(NPROPS) contient les proprietes materielles (cal- =
C=  culees) en FIN de pas de temps (a l'instant "t+Dt").               =
C=======================================================================

      SUBROUTINE UMATHPP ( nLdC, iCALC, TIME,DTIME, TEMP,DTEMP,
     &                    STRAN,DSTRAN, PROPS,NPROPS, PREDEF,DPRED,
     &                    STATEV,NSTATV, STRESS,NTENS, KINC )

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      DIMENSION STRAN(*),DSTRAN(*),PROPS(*),PREDEF(*),DPRED(*),
     &          STATEV(*),STRESS(*)

C  -----------------------
C  Coefficients materiau : Vecteur PROPS(1:NPROPS)
C  -----------------------
C  
C  PROPS ( 1) : 'YOUN'   PROPS ( 2) : 'NU  '   PROPS(3) : '   '
C
C  --------------------
C  Variables internes : Vecteur STATEV(1:NSTATEV)
C  --------------------
C  
C  STATEV( 1) : '    '
C  
C  STATEV( 2) : '    '   STATEV( 3) : '    '   STATEV( 4) : '    '
C  
C  STATEV( 5) : '    '   STATEV( 6) : '    '   STATEV( 7) : '    '
C
C  ---------------------
C  Parametres externes : Vecteurs PREDEF(1:*) et DPRED(1:*)
C  ---------------------
C  

C  Pas de temps minimal
C  (eviter un pas de temps nul en cas de convergence forcee)
      PARAMETER ( cPasMin = 1.E-15 )

C  Variables locales utilisees pour l'integration
c  Coefficients materiaux : Tous les coefficients du .dgibi
c  Variables internes : 12gammap, 12rho, 6EPSP
c  Variables auxiliaires:6stresses,TAUMAX
      PARAMETER (nCoe=27,nChgt=2,nVInt=42,nVAux=10)
      DIMENSION Coef(nCoe),Chgt(nChgt),dChgt(nChgt),
     &          VInt(nVInt),dVIntdt(nVInt),VAux(nVAux)
      PARAMETER (MAXS=12)
      DIMENSION SIGMA(9),RELAS(9),E0CRIS(9),E1CRIS(9),DELTAEPS(9),
     &DEV(9),ETEMP(9),EPLUS(9),DEFOR(9),DDEFOR(9)
      
      LOGICAL check
      
      include 'ParamRK.inc'
C-- 1 --
C  Coef contient les "coefficients" utiles pour l'integration :
C  S'il n'y a aucun coefficient necessaire a l'integration, 
C  ecrire "nCoe=0" et "DIMENSION Coef(1)"
C-- 2 --
C  Definition des valeurs des parametres de chargement utiles a l'integration :
C  - Chgt(-)  : valeur des parametres a l'instant t
C  - dChgt(-) : increment de Chgt(-) sur le pas de temps DTIME
C  La valeur des parametres de chargement a l'instant t+x*DTIME vaut :
C   Chgt_en_x(-) = Chgt(-) + x * dChgt(-) avec 0 <= x <= 1.
C   * Chgt(1) = temperature (=TEMP)
C  S'il n'y a aucun parametre de chargement necessaire a l'integration, 
C  ecrire "nChgt=0" et "DIMENSION Chgt(1)"
C-- 3 --
C  VInt contient uniquement les variables internes a integrer.
C-- 4 --
C  VAux contient les variables auxiliaires utiles a l'integration.
C  Par exemple, les contraintes a chaque point de calcul et, en particulier,
C  en fin de pas, suite a l'integration.
C  S'il n'y a aucune variable auxiliaire, ecrire "nVaux=0" et 
C  "DIMENSION Vaux(1)"

C  Initialisation des calculs
C  "Coefficients" materiaux :
      do i=1,nCoe
      Coef(i) = PROPS(i)
      enddo
      BURG = coef(11)
      RHI  = coef(16)
      rhir = coef(20)
C  Parametres de "chargement"
      Chgt(1)  = TEMP
      dChgt(1) = DTEMP
C  Recopie des variables internes a integrer
      do i=1,MAXS
      VInt(i) = STATEV(i)
      VInt(MAXS+i)=STATEV(MAXS+i)
      Vint(2*MAXS+i)=STATEV(2*MAXS+i)
      enddo
      do i=1,6
      VInt(3*MAXS+i)=STATEV(3*MAXS+i)
      enddo      
c-------------------------------------------------------------------
c--------initialisation au tout debut du calcul
c
c     les densites de disloc et les densites de défauts d'irradiation
c     ----------------------
      DO 704 I=1,MAXS
       if (VInt(MAXS+I).eq.0.0D0) then
        VInt(MAXS+I)=RHI*BURG**2
c        if (i.eq.5) VInt(MAXS+I)=10.0D0*VInt(MAXS+I)
       endif
       if (VInt(2*MAXS+I).eq.0.0D0) VInt(2*MAXS+I)=rhir*BURG**2
704   CONTINUE 

c ----calcul de la matrice de passage repere cristallo>repere materiau
       c1=dcos(coef(nCoe-2))
       s1=dsin(coef(nCoe-2))
       c2=dcos(coef(nCoe))
       s2=dsin(coef(nCoe))
       cg=dcos(coef(nCoe-1))
       sg=dsin(coef(nCoe-1))
c
       RELAS(1)=c1*c2-s1*s2*cg
       RELAS(5)=-s1*s2+c1*c2*cg
       RELAS(9)=cg
       RELAS(4)=-c1*s2-s1*c2*cg
       RELAS(7)=s1*sg
       RELAS(8)=-c1*sg
       RELAS(2)=s1*c2+c1*s2*cg
       RELAS(3)=s2*sg
       RELAS(6)=c2*sg
c       
c------Fin initialisation 
c      DO i=1,NTENS
c       write(6,*) 'STRESS(',i,')=',STRESS(i)
c      ENDDO
c
c     TEST SUR LE PAS DE DEFORMATION AVEC REDECOUPAGE AUTOMATIQUE DU PAS DE TEMPS
C     SI PAS DE DEFORMATION TROP GRAND

c     PASSAGE VECTEUR A MATRICE.
      CALL VECTOMAT(STRAN,DEFOR)
      CALL VECTOMAT(DSTRAN,DDEFOR)
      DO I=1,9
       IF ((I.NE.1).AND.(I.NE.5).AND.(I.NE.9)) THEN
       DEFOR(I)=0.5D0*DEFOR(I)
       DDEFOR(I)=0.5D0*DDEFOR(I)
       ENDIF
      ENDDO
c
      npart=0
1975  test=0.0D0
      DO I=1,9
        DELTAEPS(I)=(DDEFOR(I))/2.D0**npart
        ETEMP(I)=DEFOR(I)
        EPLUS(I)=ETEMP(I)+DELTAEPS(I)
      ENDDO
      traceps=DELTAEPS(1)+DELTAEPS(5)+DELTAEPS(9)
      do i=1,9
       IF ((I.NE.1).AND.(I.NE.5).AND.(I.NE.9)) THEN
        DEV(I)=DELTAEPS(I)
       else
        DEV(I)=DELTAEPS(I)-1.0D0/3.0D0*traceps
       endif
       test = test+DEV(I)**2
      enddo
      test=(2.0D0/3.0D0*test)**0.5
c      write(6,*) 'test=',test
      if (test.ge.pasF) then
c       write(6,*) 'DELTAEPS trop grand, redecoupage du pas de temps,',
c     &' npart=',npart, ' test=',test
       npart=npart+1
c       write(6,*) 'npart=',npart
       goto 1975
      endif
C  Increment de temps (s)
       DTIMEPART = DTIME/cDeux**npart
       PasDt = MAX (DTIMEPART,cPasMin)

C***********************************************************************
c                                                        *
c  SOUS ITERATIONS SI LES PAS DE DEFORMATIONS PROPOSES   *
c        PAR CAST3M SONT TROP IMPORTANTS                 *
c                                                        *
c*********************************************************
c
      DO 1976 ibou=1,2**npart
c___________________________________________
c
c      PASSAGE DANS LE REPERE CRISTALLIN 
c___________________________________________
c
      call ROTINVCC(ETEMP,RELAS,E0CRIS)
      call ROTINVCC(EPLUS,RELAS,E1CRIS)
c      write(6,*) 'E0CRIS='
c      write(6,*) E0CRIS(1),E0CRIS(4),E0CRIS(7)
c      write(6,*) E0CRIS(2),E0CRIS(5),E0CRIS(8)
c      write(6,*) E0CRIS(3),E0CRIS(6),E0CRIS(9)
c      write(6,*) 'E1CRIS='
c      write(6,*) E1CRIS(1),E1CRIS(4),E1CRIS(7)
c      write(6,*) E1CRIS(2),E1CRIS(5),E1CRIS(8)
c      write(6,*) E1CRIS(3),E1CRIS(6),E1CRIS(9)
c-------------------------------------------------------------------
c
c      write(6,*) 'entree dans UGD_RK43'
C  Integration sur le pas de temps
      CALL UHPP_RK43(nLdC, PasDt, VInt,dVIntdt,nVInt, Coef,nCoe,
     &         Chgt,dChgt,nChgt,E0CRIS,E1CRIS, VAux,nVAux, KINC,check)
c
       if (check) then
        npart=npart+1
        write(6,*) 'erreur trop grande dans RK, npart=',npart
        goto 1975
       endif
      do i=1,9
       ETEMP(i)=ETEMP(i)+DELTAEPS(i)
       EPLUS(i)=EPLUS(i)+DELTAEPS(i)
      enddo
1976  continue


c      write(6,*) 'sortie de UGD_RK43'
C  Mise a jour des contraintes a la fin du pas de temps
C  --> Si calculees lors de l'integration (sous-programme UGD_RKMOD)
c___________________________________________
c
c      PASSAGE DANS LE REPERE MACRO 
c___________________________________________
c
       CALL VECTOMAT(VAux,SIGMA)
       call ROTATECC(SIGMA,RELAS,SIGMA)
       CALL MATTOVEC(SIGMA,STRESS)
c      DO i=1,NTENS
c       write(6,*) 'STRESS(',i,')=',STRESS(i)
c      ENDDO

C  Mise a jour de variables internes integrees en fin de pas
      do i=1,3*MAXS
      STATEV(i) = VInt(i) 
      enddo
      do i=1,6
      STATEV(3*MAXS+i) = VInt(3*MAXS+i)
      enddo
c    Mise a jour du TauMax
      STATEV(3*MAXS+6+1)=VAux(6+1)
C    Mise a jour du glissement plastique cumule sur tous les syst
      STATEV(3*MAXS+9)=cZero
      do i=1,MAXS
       STATEV(3*MAXS+9)=STATEV(3*MAXS+9)+dabs(STATEV(i))
      enddo
C    Mise a jour de la densite de disloc cumulee sur tous les syst
       STATEV(3*MAXS+10)=cZero
       do i=1,MAXS
        STATEV(3*MAXS+10)=STATEV(3*MAXS+10)+STATEV(MAXS+i)
       enddo
c    Mise a jour du tauTL, Longueur de segment vis et Tau_LR
c    pour systeme principal
        STATEV(3*MAXS+11)=VAux(7+1)
        STATEV(3*MAXS+12)=VAux(7+2)
        STATEV(3*MAXS+13)=VAux(7+3)
       
c      DO i=1,NSTATV
c       write(6,*) 'STATEV(',i,')=',STATEV(i)
c      ENDDO
c       write(6,*) 'sortie ugd_umat'      

      
C  Mise a jour de variables internes calculees en fin de pas
C  --> Variables internes qui n'ont pas a etre integrees
C Par ex. :   STATEV(1) = Vaux(1) / Young - VInt(1)

      RETURN
      END

