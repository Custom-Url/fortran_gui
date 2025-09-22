C=======================================================================
C=  TOUS MODELES DE COMPORTEMENT :  sousprogramme UHPP_RKMOD            =
C=======================================================================
C=  Calcul de l'evolution des grandeurs du MODELE a l'instant "x"      =
C=                                                                     =
C=  Entrees :                                                          =
C=  ---------                                                          =
C=  nLdC     Numero du modele de comportement (idem que dans UMAT)     =
C=  x        Instant de calcul = t+x*Dt                                =
C=  y        Valeurs des "variables" du modele a l'instant "x"         =
C=  ny       Nombre de "variables" du modele a integrer                =
C=  coef     "Coefficients" du modele (parametres materiau + autres)   =
C=  ncoef    Nombre de "coefficients" du modele                        =
C=  paext    Parametres de chargement externes a l'instant "x"         =
C=  dpaextdx Vitesse de variation des parametres externes (constantes) =
C=  npaext   Nombre de parametres de chargement externes               =
C=  epsx     Tenseur des deformations a l'instant "x"                  =
C=  depsdt   Vitesse de deformation a l'instant "x"                    =
C=  nvaux    Nombre de grandeurs "auxiliaires"                         =
C=                                                                     =
C=  Sorties :                                                          =
C=  ---------                                                          =
C=  dydx     Derivees des "variables" a l'instant de calcul "x"        =
C=  vaux     Grandeurs "auxiliaires" mises a jour si necessaire        =
C=  iok      Indicateur de problemes (different de 1 dans ce cas)      =
C=           ( Cas particulier : iok = 999 pour forcer le decoupage    =
C=             souspas courant, iok remis a 1 par la suite )           =
C=======================================================================
C= Les variables "epsx" et "depsdt" sont des tenseurs 3x3 stockes sous =
C= forme de vecteurs a 9 composantes : T(i,j) -> T(k) avec k=i+3*(j-1).=
C=======================================================================

      SUBROUTINE UHPP_RKMOD (nLdC, x, y,dydx,ny, coef,ncoef,
     &                      paext,dpaextdx,npaext, epsx,depsdt,
     &                      vaux,nvaux, iok)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'

      PARAMETER (MAXS=12,MAXS3=36)
      DIMENSION y(*),dydx(*), coef(*), paext(*),dpaextdx(*),
     &          epsx(*),depsdt(*), vaux(*)
      DIMENSION EPSP(9),DEPSP(3,3),EEL(9),r(MAXS),
     &GNS(3,3,MAXS),GNSSYM(3,3,MAXS),TAUS(MAXS),SIGMA(9),
     &ASU(MAXS,MAXS),XNSYS(MAXS3),XMSYS(MAXS3),
     &DELTASU(MAXS,MAXS),rirr(MAXS)
     
      INCLUDE 'nsms.inc'
      INCLUDE 'inter.inc'

C=======================================================================
C=  EXEMPLE :
C=======================================================================
c      IF ( nLdC .EQ. 80 ) THEN

c--------initialisation 

c    ATTENTION ISOTROPIE ELASTIQUE !!!!!
c
      YOUNG = coef(1)
      petinu= coef(2)
      Temper=coef(24)
      YT = YOUNG - 0.04589998*TEMPER*1000.
c      C11=coef(3)
c      C12=coef(4)
c      C44=coef(5)
c      pmu = C44-0.017*Temper*1000.
      C11 = YT*(cUn-petinu)/(cUn-cDeux*petinu)/(cUn+petinu)
      C12 = YT*petinu/(cUn-cDeux*petinu)/(cUn+petinu)
      C44 = YT/cDeux/(cUn+petinu)
      pmu=C44
      tau0=coef(6)
      tauF=coef(7)
      Dg=coef(8)
      Dobs=coef(9)
      granh=coef(10)
      burg=coef(11)
      gkself=coef(12)
      gksgkf=coef(13)
      yAT=coef(14)
      deltag0=coef(15)
      rhoini=coef(16)
      rhom=coef(17)
      gp0=coef(18)
      pn=coef(19)
      rhir=coef(20)
      airr=coef(21)
      xsi=coef(22)
      gpmoyen=coef(23)
      petitk=8.617D-5
      dg0kt=deltag0/petitk/Temper
      petitlc = 1.4D-9*Temper**CDeux
      gkf=gkself/min(cUn+(gksgkf-cUn)/300.d0*TEMPER,gksgkf)
      
c      write(6,*) pmu,tau0,tauF,Dg,Dobs,granh,burg,gkself,
c     &gksgkf,yAT,deltag0,rhoini,rhom,gp0,pn,rhirr,airr,xsi,gpmoyen
   
c      Deformations plastiques
c
       call VECTOMAT(y(3*MAXS+1),EPSP)
       DO I=1,9
       IF ((I.NE.1).AND.(I.NE.5).AND.(I.NE.9)) THEN
       EPSP(I)=0.5D0*EPSP(I)
       ENDIF
       ENDDO
c      write(6,*) 'EPSP='
c      write(6,*) EPSP(1),EPSP(4),EPSP(7)
c      write(6,*) EPSP(2),EPSP(5),EPSP(8)
c      write(6,*) EPSP(3),EPSP(6),EPSP(9)

c     Somme des rhos et initialisation des variables r(is) et rirr(is)
      DO I=1,MAXS
       r(I)=y(MAXS+I)
       rirr(I)=y(2*MAXS+I)
      ENDDO
c

c   Calcul de la deformation elastique
      call UHPP_CINE0(EPSP,EPSX,EEL)
c      write(6,*) 'EEL='
c      write(6,*) EEL(1),EEL(4),EEL(7)
c      write(6,*) EEL(2),EEL(5),EEL(8)
c      write(6,*) EEL(3),EEL(6),EEL(9)

c   Calcul de la contrainte
      call ELASTIC(EEL,C11,C12,C44,SIGMA)
c      write(6,*) 'SIGMA='
c      write(6,*) SIGMA(1),SIGMA(4),SIGMA(7)
c      write(6,*) SIGMA(2),SIGMA(5),SIGMA(8)
c      write(6,*) SIGMA(3),SIGMA(6),SIGMA(9)

c   Calcul de la cission resolue sur chaque system s
      call CGNS(GNS)
      call CISRED(SIGMA,GNS,TAUS)

c----------------------------------------------------
c
c     Boucle sur les systemes de glissement
c----------------------------------------------------
c
      VAux(7)=cZero
      DO 7005 IS=1,MAXS

c        Calcul de rhof(is)
       rhof=0.0D0
       DO JS=1,MAXS
        rhof=rhof+DELTASU(IS,JS)*r(JS)
       ENDDO
       rhof=rhof+rirr(IS)

c      Calcul de alphaat(is)
         alphaat=0.0D0
         DO JS=1,MAXS
          alphaat=alphaat+
     &deltasu(IS,JS)*asu(is,js)*r(js)/(rhof)
       ENDDO
       alphaat=alphaat+airr*rirr(is)/(rhof)
       alphaat=alphaat**0.5

c      Calcul du rayon de courbure approximatif de la portion non-vis
c      des dislocations
c----------------------------------------------------
       petitlsv=burg/(rhof)**0.5
       fact = rhom*burg*granh*petitlsv/gpmoyen
       deltaG = min(deltaG0,petitk*Temper*log(fact))
       taueff = tau0*(CUn-deltaG/deltag0)**CDeux
       if (taueff.lt.(pmu*burg/CDeux*1.D-6)) then
        Rs=1.D6
       else
        Rs=pmu*burg/CDeux/(taueff)
       endif
       testmin = min((rhof**0.5)/burg,
     &((Dobs+CDeux*Rs)*(rhof))/burg**cDeux)
       GLS=CUn/testmin- Dobs
       test = CDeux*alphaat*Rs+petitlc
       if (GLS.gt.test) then
        tauTLs=0.D0
        petitlsv = GLS-CDeux*alphaat*Rs
       else
        tauTLs=alphaat*pmu*burg*(CUn/GLS-CUn/
     &(CDeux*alphaat*Rs+petitlc))
        petitlsv = petitlc
       endif

c     Calcul de l'increment de deformation viscoplastique sur chaque systeme
c----------------------------------------------------------------------------
c
       if (VAux(7).lt.dabs(TAUS(IS))) VAux(7)=dabs(TAUS(IS))
       tauLD=pmu*(asu(IS,IS)*r(IS))**0.5
       tauc = tauF+(tauLD**cDeux+tauTLs**cDeux)**0.5
       taueff = max(0.0D0,dabs(taus(IS))-tauc)
       if (taueff.le.0.0D0) then
        dydx(IS)=0.0D0
        dydx(MAXS+IS)=0.0D0
        dydx(2*MAXS+IS)=0.0D0
       else
c      pour prendre en compte l'augmentation rapide de rho
c      sur le système actif pendant le stade 0
        dgnuc = rhom*burg*granh*petitlsv*
     &dexp(-dg0kt*(cUn-(taueff/tau0)**0.5))*sign(cUn,TAUS(IS))
c        if (IS.EQ.5) dgnuc=10.0D0*dgnuc
c      attention, modif egalement de la loi d'ecoulement
        dgmig = sign(cUn,TAUS(IS))*gp0*(dabs(taus(IS))/
     &tauc)**pn
        dydx(IS) = cUn/dgnuc+cUn/dgmig
        dydx(IS) = cUn/dydx(IS)

        ys= cUn/(cUn/yAT+cDeux*cPi*taueff/pmu/burg)

c     Evolution de la densite de dislocations sur chaque systeme
c-----------------------------------------------------------------
       GLIBPM = burg/Dg + (CUn-taueff/tau0)*((asu(IS,IS)*r(IS))**0.5/
     &gkself + alphaat*GLS*rhof/burg/gkf)
       dydx(MAXS+IS)=(GLIBPM-ys*r(IS)/BURG)*DABS(dydx(IS))

c     Evolution de la densite de défauts d'irradiation sur chaque systeme
c-----------------------------------------------------------------
       dydx(2*MAXS+IS)=-xsi*rirr(iS)*DABS(dydx(IS))

        endif

       if (IS.EQ.5) then
       VAux(7+1)=tauTLs
       VAux(7+2)=petitlsv
       VAux(7+3)=tauLD
       endif

7005  CONTINUE
c
c        Loi d'evolution de EPSP
       DO i=1,3
        do j=1,3
        DEPSP(i,j)=0.0D0
        enddo
       ENDDO

      call CGNSSYM(GNSSYM)

      DO IS=1,MAXS
       DO I=1,3
        DO J=1,3
         DEPSP(I,J)=DEPSP(I,J)+GNSSYM(I,J,IS)*dydx(IS)
        ENDDO
       ENDDO
      ENDDO
      call MATTOVEC(DEPSP,dydx(3*MAXS+1))
      dydx(3*MAXS+4)=cDeux*dydx(3*MAXS+4)
      dydx(3*MAXS+5)=cDeux*dydx(3*MAXS+5)
      dydx(3*MAXS+6)=cDeux*dydx(3*MAXS+6)
c
c      write(6,*) 'depsp='
c      write(6,*) dydx(3*MAXS+1),dydx(3*MAXS+4),dydx(3*MAXS+5)
c      write(6,*) dydx(3*MAXS+4),dydx(3*MAXS+2),dydx(3*MAXS+6)
c      write(6,*) dydx(3*MAXS+5),dydx(3*MAXS+6),dydx(3*MAXS+3)


c   Mise à jour des Variables auxiliaires
C  Calcul des contraintes a l'instant courant :
C  -> A stocker dans "vaux" si besoin
C       vaux(1) = ...

      VAux(1)=SIGMA(1)
      VAux(2)=SIGMA(5)
      VAux(3)=SIGMA(9)
      VAux(4)=SIGMA(4)
      VAux(5)=SIGMA(7)
      VAux(6)=SIGMA(8)
c      do i=1,7
c       write(6,*) 'VAux(',i,')=',VAux(i)
c      enddo
C  Calcul des differentes derivees :
C       IF ( SEQ .GT. seuil ) THEN
C         vpcum =  vitesse_plastique_fct_de_vaux
C         dydx( 1) = vpcum
C         ...
C       ELSE
C         dydx( 1) = cZero
C         ...
C       ENDIF
        RETURN


C=======================================================================
C=======================================================================
c      ELSE
c        iok = -1
c      ENDIF

      RETURN
      END

