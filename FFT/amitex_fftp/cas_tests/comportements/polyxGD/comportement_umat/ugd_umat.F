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
C=  DFGRD0   Gradient de la transformation a l'instant "t"             =
C=  DFGRD1   Gradient de la transformation a l'instant "t+Dt"          =
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

      SUBROUTINE UMATGD ( nLdC, iCALC, TIME,DTIME, TEMP,DTEMP,
     &                    DFGRD0,DFGRD1, PROPS,NPROPS, PREDEF,DPRED,
     &                    STATEV,NSTATV, STRESS,NTENS, KINC )

      include 'IMPLICIT.INC'
      include 'CMATH.INC'

      DIMENSION DFGRD0(*),DFGRD1(*),PROPS(*),PREDEF(*),DPRED(*),
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
c  Variables internes : 12gammap, 12rho, 9FP
c  Variables auxiliaires:6stresses,TAUMAX
      PARAMETER (nCoe=28,nChgt=2,nVInt=81,nVAux=31)
      DIMENSION Coef(nCoe),Chgt(nChgt),dChgt(nChgt),
     &          VInt(nVInt),dVIntdt(nVInt),VAux(nVAux)
      PARAMETER (MAXS=24)
      DIMENSION SIGMA(9),RELAS(9),F0CRIS(9),F1CRIS(9),
     &DELTAF(9),FTEMP(9),FPLUS(9),
     &GPINI(MAXS),ETOT(9),EDEV(9)
      
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
C  VAux(1) jusqu'à VAux(6) : tenseur des contraintes
C  VAux(7) TAUMAX
C  VAux(8) jusqu'à VAux(31) : incrément de glissement plastique par système

C  Initialisation des calculs
C  "Coefficients" materiaux :
      do i=1,nCoe
      Coef(i) = PROPS(i)
      enddo
      BURG = coef(13)
      RHI  = coef(18)
      rhir = coef(22)
C  Parametres de "chargement"
      Chgt(1)  = TEMP
      dChgt(1) = DTEMP
C  Recopie des variables internes a integrer
      do i=1,MAXS
      VInt(i) = STATEV(i)
      VInt(MAXS+i)=STATEV(MAXS+i)
      Vint(2*MAXS+i)=STATEV(2*MAXS+i)
      enddo
      do i=1,9
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
       endif
       if (VInt(2*MAXS+I).eq.0.0D0) VInt(2*MAXS+I)=rhir*BURG**2
704   CONTINUE 

c      Partie plastique du gradient de la transformation, FP
c     dans la base micro
      j=3*MAXS
      if ((VInt(j+1).eq.cZero).and.(VInt(j+2).eq.cZero).and.
     & (VInt(j+3).eq.cZero).and.(VInt(j+4).eq.cZero).and.
     & (VInt(j+5).eq.cZero).and.(VInt(j+6).eq.cZero).and.
     & (VInt(j+7).eq.cZero).and.(VInt(j+8).eq.cZero).and.
     & (VInt(j+9).eq.cZero)) then
       VInt(j+1)=cUn
       VInt(j+5)=cUn
       VInt(j+9)=cUn
      endif

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
c------- taux de glissement plastique sur chaque système au pas de temps précédent
      DO i=1,MAXS
       VAUX(7+i)=STATEV(3*MAXS+9+i)
      ENDDO

c     TEST SUR LE PAS DE F AVEC REDECOUPAGE AUTOMATIQUE DU PAS DE TEMPS
C     SI PAS DE F TROP GRAND
      npart=0
1975  test=0.0D0
      DO I=1,9
        DELTAF(I)=(DFGRD1(I)-DFGRD0(I))/2.D0**npart
	test = test+DELTAF(I)**2
	FTEMP(I)=DFGRD0(I)
	FPLUS(I)=FTEMP(I)+DELTAF(I)
      ENDDO
      test=(2.0D0/3.0D0*test)**0.5
      if (test.ge.pasF) then
c       write(6,*) 'DELTAF trop grand, redecoupage du pas de temps,',
c     &' npart=',npart, ' test=',test
       npart=npart+1
c       write(6,*) 'npart=',npart
       goto 1975
      endif
C  Increment de temps (s)
       DTIMEPART = DTIME/cDeux**npart
       PasDt = MAX (DTIMEPART,cPasMin)
       DO I=1,MAXS
        GPINI(i)= STATEV(i)
       ENDDO

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
      call ROTINVCC(FTEMP,RELAS,F0CRIS)
      call ROTINVCC(FPLUS,RELAS,F1CRIS)
c-------------------------------------------------------------------
c
c      write(6,*) 'entree dans UGD_RK43'
C  Integration sur le pas de temps
      CALL UGD_RK43(nLdC, PasDt, VInt,dVIntdt,nVInt, Coef,nCoe,
     &         Chgt,dChgt,nChgt,F0CRIS,F1CRIS, VAux,nVAux, KINC,check)
c      CALL UGD_RK21(nLdC, PasDt, VInt,dVIntdt,nVInt, Coef,nCoe,
c     &         Chgt,dChgt,nChgt,F0CRIS,F1CRIS, VAux,nVAux, KINC,check)
c
       if (check) then
        npart=npart+1
        write(6,*) 'erreur trop grande dans RK, npart=',npart
        goto 1975
       endif
      do i=1,9
       FTEMP(i)=FTEMP(i)+DELTAF(i)
       FPLUS(i)=FPLUS(i)+DELTAF(i)
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
       SIGMA(1)=VAux(1)
       SIGMA(5)=VAux(2)
       SIGMA(9)=VAux(3)
       SIGMA(4)=VAux(4)
       SIGMA(7)=VAux(5)
       SIGMA(8)=VAux(6)
       SIGMA(2)=SIGMA(4)
       SIGMA(3)=SIGMA(7)
       SIGMA(6)=SIGMA(8)
      call ROTATECC(SIGMA,RELAS,SIGMA)
      STRESS(1) = SIGMA(1)
      STRESS(2) = SIGMA(5)
      STRESS(3) = SIGMA(9)
      STRESS(4) = SIGMA(4)
      STRESS(5) = SIGMA(7)
      STRESS(6) = SIGMA(8)
c      DO i=1,NTENS
c       write(6,*) 'STRESS(',i,')=',STRESS(i)
c      ENDDO

C  Mise a jour de variables internes integrees en fin de pas
      do i=1,3*MAXS
      STATEV(i) = VInt(i) 
      enddo
      do i=1,9
      STATEV(3*MAXS+i) = VInt(3*MAXS+i)
      enddo
c    Mise à jour du TauMax
      STATEV(4*MAXS+9+1)=VAux(6+1)
c    Mise à jour de la déformation equivalente au sens de von Mises, EPSE
       call EELFE(DFGRD1,ETOT)
       do i=1,9
        EDEV(i)=ETOT(i)
       enddo
       EDEV(1)=ETOT(1)-c1s3*(ETOT(1)+ETOT(2)+ETOT(3))
       EDEV(2)=ETOT(2)-c1s3*(ETOT(1)+ETOT(2)+ETOT(3))
       EDEV(3)=ETOT(3)-c1s3*(ETOT(1)+ETOT(2)+ETOT(3))
       ESCAL=cZero
       do i=1,9
        ESCAL = ESCAL+EDEV(i)**cDeux
       enddo
       STATEV(4*MAXS+9+2)=(cDeux*c1s3*ESCAL)**0.5      
       
      DO i=1,MAXS
       STATEV(3*MAXS+9+i)=(VInt(i)-GPINI(i))/DTIME
      ENDDO
      
C  Mise a jour de variables internes calculees en fin de pas
C  --> Variables internes qui n'ont pas a etre integrees
C Par ex. :   STATEV(1) = Vaux(1) / Young - VInt(1)

      RETURN
      END

