C=======================================================================
C=  TOUS MODELES DE COMPORTEMENT :  sousprogramme UGD_RKMOD            =
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
C=  grdfx    Gradient de la transformation a l'instant "x"             =
C=  grdlx    Vitesse  de la transformation a l'instant "x"             =
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
C= Les variables "grdfx" et "grdlx" sont des tenseurs 3x3 stockes sous =
C= forme de vecteurs a 9 composantes : T(i,j) -> T(k) avec k=i+3*(j-1).=
C=======================================================================

      SUBROUTINE UGD_RKMOD (nLdC, x, y,dydx,ny, coef,ncoef,
     &     paext,dpaextdx,npaext, grdfx,grdlx,
     &     vaux,nvaux, iok)

      include 'IMPLICIT.INC'
      include 'CMATH.INC'

      PARAMETER (MAXS=24,MAXS3=72)
      DIMENSION y(*),dydx(*), coef(*), paext(*),dpaextdx(*),
     &     grdfx(*),grdlx(*), vaux(*)
      DIMENSION FP(9),FE(9),EEL(9),r(MAXS),PIOLA2(9),
     &     DEEPLUSI(9),TAUB(9),GNS(3,3,MAXS),TAUS(MAXS),SIGMA(9),
     &     ASU(MAXS,MAXS),XNSYS(MAXS3),XMSYS(MAXS3),
     &     GNSFP(9),DELTASU(MAXS,MAXS),rirr(MAXS)
      
      INCLUDE 'nsms.inc'
      INCLUDE 'inter.inc'
      
      IF ( nLdC .EQ. 91 ) THEN
         STOP
      endif
C=======================================================================
C     =  EXEMPLE :
C=======================================================================
c     IF ( nLdC .EQ. 90 ) THEN

c--------initialisation 

c     ATTENTION ISOTROPIE ELASTIQUE !!!!!
c     
      YOUNG = coef(1)
      petinu= coef(2)
      Temper=coef(ncoef-3)
      YT = YOUNG - 0.04589998*TEMPER*1000.
c     C11=coef(3)
c     C12=coef(4)
c     C44=coef(5)
c     pmu = C44-0.017*Temper*1000.
      C11 = YT*(cUn-petinu)/(cUn-cDeux*petinu)/(cUn+petinu)
      C12 = YT*petinu/(cUn-cDeux*petinu)/(cUn+petinu)
      C44 = YT/cDeux/(cUn+petinu)
      pmu=C44
      tau01=coef(6)
      tau02T=coef(7)
      tau02AT = coef(8)
      tauF=coef(9)
      Dg=coef(10)
      Dobs=coef(11)
      granh=coef(12)
      burg=coef(13)
      gkself=coef(14)
      gksgkf=coef(15)
      yAT=coef(16)
      deltag0=coef(17)
      rhoini=coef(18)
      rhom=coef(19)
      gp0=coef(20)
      pn=coef(21)
      rhir=coef(22)
      airr=coef(23)
      xsi=coef(24)
      petitk=8.617D-5
      dg0kt=deltag0/petitk/Temper
      petitlc = 300.d0*burg*(Temper/400.D0)**CDeux
      gkf=gkself/min(cUn+(gksgkf-cUn)/300.d0*TEMPER,gksgkf)
      
c     Partie plastique du gradient de la transformation, FP
c     
      do i=1,9
         FP(i)=y(3*MAXS+i)
      enddo

c     Somme des rhos et initialisation des variables r(is) et rirr(is)
      DO I=1,MAXS
         r(I)=y(MAXS+I)
         rirr(I)=y(2*MAXS+I)
      ENDDO
c     

c     Calcul de la deformation elastique de GREEN-LAGRANGE
      call FEFFP(grdfx,FP,FE)
      call EELFE(FE,EEL)

c     Calcul de la contrainte de Piola-Kirchoff 2
      call ELASTIC99(EEL,C11,C12,C44,PIOLA2)

c     Calcul de la cission resolue sur chaque system s
      call PRODT33(FE,FE,DEEPLUSI)
      call PROD33(DEEPLUSI,PIOLA2,TAUB)
      call CGNS(GNS)
      call CISRED(TAUB,GNS,TAUS)
c----------------------------------------------------
c     
c     Boucle sur les systemes de glissement
c----------------------------------------------------
c     
      VAux(7)=cZero
      DO 7005 IS=1,MAXS

c     Calcul de rhof(is)
         rhof=0.0D0
         DO JS=1,MAXS
            rhof=rhof+DELTASU(IS,JS)*r(JS)
         ENDDO
         rhof=rhof+rirr(IS)

c     Calcul de alphaat(is)
         alphaat=0.0D0
         DO JS=1,MAXS
            alphaat=alphaat+
     &           deltasu(IS,JS)*asu(is,js)*r(js)/(rhof)
         ENDDO
         alphaat=alphaat+airr*rirr(is)/(rhof)
         alphaat=alphaat**0.5

c     Calcul du rayon de courbure approximatif de la portion non-vis
c     des dislocations
c----------------------------------------------------
         if (IS.lt.13) then
            tau0 = tau01
         else if (TAUS(IS).ge.0.0D0) then
            tau0 = tau02T
         else
            tau0 = tau02AT
         endif
         petitlsv=burg/(rhof)**0.5
         fact = rhom*burg*granh*petitlsv/max(dabs(VAUX(7+IS)),1.D-8)
         deltaG = min(deltaG0,petitk*Temper*log(fact))
         taueff = tau0*(CUn-deltaG/deltag0)**CDeux
         if (taueff.lt.(pmu*burg/CDeux*1.D-6)) then
            Rs=1.D6
         else
            Rs=pmu*burg/CDeux/(taueff)
         endif
         testmin = min((rhof**0.5)/burg,
     &        ((Dobs+CDeux*Rs)*(rhof))/burg**cDeux)
         GLS=CUn/testmin- Dobs
         test = CDeux*alphaat*Rs+petitlc
         if (GLS.gt.test) then
            tauTLs=0.D0
            petitlsv = GLS-CDeux*alphaat*Rs
         else
            tauTLs=alphaat*pmu*burg*(CUn/GLS-CUn/
     &           (CDeux*alphaat*Rs+petitlc))
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
            dgnuc = rhom*burg*granh*petitlsv*
     &           dexp(-dg0kt*(cUn-(taueff/tau0)**0.5))
     &           *sign(cUn,TAUS(IS))
            dgmig = sign(cUn,TAUS(IS))*gp0*(dabs(taus(IS))/
     &           tauc)**pn
            dydx(IS) = cUn/dgnuc+cUn/dgmig
            dydx(IS) = cUn/dydx(IS)

            ys= cUn/(cUn/yAT+cDeux*cPi*taueff/pmu/burg)

c     Evolution de la densite de dislocations sur chaque systeme
c-----------------------------------------------------------------
            GLIBPM = burg/Dg + (CUn-taueff/tau0)*((asu(IS,IS)
     &           *r(IS))**0.5/gkself + alphaat*GLS*rhof/burg/gkf)
            dydx(MAXS+IS)=(GLIBPM-ys*r(IS)/BURG)*DABS(dydx(IS))

c     Evolution de la densite de défauts d'irradiation sur chaque systeme
c-----------------------------------------------------------------
            dydx(2*MAXS+IS)=-xsi*rirr(iS)*DABS(dydx(IS))

         endif

 7005 CONTINUE
c     
c     Loi d'evolution de FP
      do i=1,9
         dydx(3*MAXS+i)=cZero
      enddo
      DO IS=1,MAXS
         call PROD33(GNS(1,1,IS),FP,GNSFP)
         do i=1,9
            dydx(3*MAXS+i)=dydx(3*MAXS+i)+dydx(IS)*GNSFP(i)
         enddo
      ENDDO

c     
c     Mise à jour des Variables auxiliaires
C     Calcul des contraintes a l'instant courant :
C     -> A stocker dans "vaux" si besoin
C     vaux(1) = ...

      call ROTATECC(PIOLA2,FE,SIGMA)
      VAux(1)=SIGMA(1)
      VAux(2)=SIGMA(5)
      VAux(3)=SIGMA(9)
      VAux(4)=SIGMA(4)
      VAux(5)=SIGMA(7)
      VAux(6)=SIGMA(8)

      RETURN
      END

