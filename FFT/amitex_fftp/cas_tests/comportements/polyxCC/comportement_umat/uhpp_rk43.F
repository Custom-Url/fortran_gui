C=======================================================================
C=  TOUS MODELES DE COMPORTEMENT :  sousprogramme UHPP_RK43             =
C=======================================================================
C=  Sousprogramme d'integration d'un modele de comportement par la me- =
C=  thode de RUNGE-KUTTA d'ordre 4 avec determination du pas de temps  =
C=  (estimation de l'erreur par une methode de RUNGE-KUTTA d'ordre 3)  =
C=  La methode d'ordre 4 est la methode de RUNGE-KUTTA "classique",    =
C=  basee sur la formule de Simpson. La methode d'ordre 3 emboitee,    =
C=  utilisee pour l'estimation de l'erreur, est de type RK-Felhberg.   =
C=                                                                     =
C=  Entrees :                                                          =
C=  ---------                                                          =
C=  nLdC     Numero du modele de comportement (utilise dans UGD_RKMod) =
C=  pas      Duree du pas d'integration                                =
C=  ny       Nombre de grandeurs a integrer                            =
C=  coef     "Coefficients" du modele (parametres materiau + autres)   =
C=  ncoef    Nombre de "coefficients" du modele                        =
C=  paext    Parametres de chargement externes au debut du pas         =
C=  dpaext   Taux de variation des parametres externes sur le pas      =
C=  npaext   Nombre de parametres de chargement externes               =
C=  epsini   Tenseur des deformations au debut du pas                  =
C=  epsfin   Tenseur des deformations a la fin du pas                  =
C=  nvaux    Nombre de grandeurs auxiliaires a calculer                =
C=                                                                     =
C=  Sorties :                                                          =
C=  ---------                                                          =
C=  y        Valeurs des "variables" du modele a la fin du pas         =
C=  dydx     Derivees des "variables" au dernier instant de calcul     =
C=  vaux     Valeurs de "variables" auxiliaires a la fin du pas        =
C=  iok      Indicateur de probleme(s)                                 =
C=======================================================================
C=  Note : La variable iok sert a indiquer qu'il y a eu un probleme au =
C=  ------ cours de l'integration. Dans ce cas, iok prend une valeur   =
C=         differente de 1. Il est inutile de l'initialiser a 1.       =
C=         (Voir la signification de KINC lors de l'appel a UMAT.)     =
C=======================================================================
C=  CE SOUSPROGRAMME EST GENERIQUE ET INDEPENDANT DU MODELE A INTEGRER =
C=======================================================================

      SUBROUTINE UHPP_RK43 (nLdC, pas, y,dydx,ny, coef,ncoef,
     &                     paext,dpaext,npaext, epsini,epsfin,
     &                     vaux,nvaux, iok,check)

      INCLUDE 'IMPLICIT.INC'
      INCLUDE 'CMATH.INC'
      INCLUDE 'ParamRK.inc'

      DIMENSION y(*),dydx(*), coef(*), paext(*),dpaext(*),
     &          epsini(*),epsfin(*), vaux(*)

      DIMENSION deps(9), epsx(9),depsdt(9)
      DIMENSION ye(nVloc), yloc(nVloc),
     &          dydx1(nVloc),dydx2(nVloc),dydx3(nVloc),dydx4(nVloc)
      DIMENSION yext(nPext), dyextdx(nPext)
      LOGICAL check

      check=.false.

CD   IF (nVloc.LT.ny) THEN
CD     WRITE(*,*) 'Valeur de nVloc a revoir'
CD     iok = -237
CD     RETURN
CD   ENDIF
CD   IF (nPext.LT.npaext) THEN
CD     WRITE(*,*) 'Valeur de nPext a revoir'
CD     iok = -238
CD     RETURN
CD   ENDIF

      x = cZero
      h = pas
      alp_x = cZero
      unspas = cUn / pas

C Calcul de l'increment du tenseur de deformation :
      CALL UHPP_CINE0(epsini,epsfin, deps)

C Initialisation des parametres de chargement et leur vitesse :
      DO i = 1, npaext
        yext(i)    = paext(i)
        dyextdx(i) = dpaext(i) * unspas
      ENDDO

C  "Tolerance" sur les grandeurs integrees :
      DO i = 1, ny
        ye(i) = MAX(toly,ABS(y(i)))
      ENDDO

C  Calcul de k1 = dydx1 en debut de pas avant entree dans boucle 10
      CALL UHPP_CINE3(epsini,epsfin,deps,alp_x,epsx)
      CALL UHPP_CINE1(deps,pas, depsdt)
      CALL UHPP_RKMod(nLdC, x, y,dydx1,ny, coef,ncoef,
     &               yext,dyextdx,npaext, epsx,depsdt,
     &               vaux,nvaux, iok)
      IF (iok.NE.1) THEN
        write(6,*) 'iok=',iok
        iok = -643
        RETURN
      ENDIF

C  Boucle 10 tant que : 1) la fin du pas est atteinte
C                    ou 2) nombre maximal d'iterations atteint
      Iter = 0
CD    IterOK = 0
CD    Iter999 = 0
C  ------------
 10   CONTINUE
        Iter = Iter + 1
C  IterMax atteint ==> ERREUR
        IF (Iter.GT.IterMax) THEN
          iok = -666
          RETURN
        ENDIF
C  Calcul de k1 = dydx1 :
C    - soit calcul en debut de pas avant entree dans boucle 10
C    - soit recuperation du sous-pas precedent (dydx1 <- dydx)
C  Calcul de termes a l'instant x+h/2
        hs2 = c1s2 * h
        xhs2 = x + hs2
        alp_x = xhs2 * unspas
        CALL UHPP_CINE3(epsini,epsfin,deps,alp_x, epsx)
        CALL UHPP_CINE1(epsx,pas,depsdt)
        DO i = 1, npaext
          yext(i) = paext(i) + alp_x * dpaext(i)
        ENDDO
C  Calcul de k2 = dydx2
        DO i = 1, ny
          yloc(i) = y(i) + hs2 * dydx1(i)
        ENDDO
        CALL UHPP_RKMod(nLdC, xhs2, yloc,dydx2,ny, coef,ncoef,
     &                 yext,dyextdx,npaext, epsx,depsdt,
     &                 vaux,nvaux, iok)
        IF (iok.NE.1) THEN
          IF (iok.EQ.999) THEN
            iok = 1
CD          Iter999 = Iter999+1
            h = hopt999 * h
            GOTO 10
          ENDIF
          RETURN
        ENDIF
C  Calcul de k3 = dydx3
        DO i = 1, ny
          yloc(i) = y(i) + hs2 * dydx2(i)
        ENDDO
        CALL UHPP_RKMod(nLdC, xhs2, yloc,dydx3,ny, coef,ncoef,
     &                 yext,dyextdx,npaext, epsx,depsdt,
     &                 vaux,nvaux, iok)
        IF (iok.NE.1) THEN
          IF (iok.EQ.999) THEN
            iok = 1
CD          Iter999 = Iter999+1
            h = hopt999 * h
            GOTO 10
          ENDIF
          RETURN
        ENDIF
C  Calcul de termes a l'instant x+h
        xh = x + h
        alp_x = xh * unspas
        CALL UHPP_CINE3(epsini,epsfin,deps,alp_x, epsx)
        CALL UHPP_CINE1(epsx,pas, depsdt)
        DO i = 1, npaext
          yext(i) = paext(i) + alp_x * dpaext(i)
        ENDDO
C  Calcul de k4 = dydx4
        DO i = 1, ny
          yloc(i) = y(i) + h * dydx3(i)
        ENDDO
        CALL UHPP_RKMod(nLdC, xh, yloc,dydx4,ny, coef,ncoef,
     &                 yext,dyextdx,npaext, epsx,depsdt,
     &                 vaux,nvaux, iok)
        IF (iok.NE.1) THEN
          IF (iok.EQ.999) THEN
            iok = 1
CD          Iter999 = Iter999+1
            h = hopt999 * h
            GOTO 10
          ENDIF
          RETURN
        ENDIF
C  Calcul de y en fin de sous-pas (estimation)
        hs6 = c1s6 * h
        DO i = 1, ny
          yloc(i) = y(i) + hs6 * ( dydx1(i)+cDeux*(dydx2(i)+dydx3(i))
     &                            +dydx4(i) )
        ENDDO
C  Calcul de dydx = k5 pour la methode emboitee (estimation de l'erreur)
        CALL UHPP_RKMod(nLdC, xh, yloc,dydx,ny, coef,ncoef,
     &                 yext,dyextdx,npaext, epsx,depsdt,
     &                 vaux,nvaux, iok)
        IF (iok.NE.1) THEN
          IF (iok.EQ.999) THEN
            iok = 1
CD          Iter999 = Iter999+1
            h = hopt999 * h
            GOTO 10
          ENDIF
          RETURN
        ENDIF
C  Calcul de l'erreur locale et du sous-pas optimal
        err = cZero
        imax=0
        DO i = 1, ny
c          err = MAX(err,ABS(dydx4(i)-dydx(i))/ye(i))
           IF (err.lt.(ABS(dydx4(i)-dydx(i))/ye(i))) then
              err=ABS(dydx4(i)-dydx(i))/ye(i)
              imax =i
           endif
        ENDDO
        err = hs6 * err
c       write(6,*) '    '
c       write(6,*) 'Iter=',Iter, ' err=',err,' imax=', imax
c       write(6,*) '    '
C  Suivant la valeur de l'erreur :
C     1) stockage resultats et calcul du sous-pas suivant
C  ou 2) diminution du sous-pas de temps et recalcul
        IF (err.LE.eps) THEN
CD        IterOK = IterOK+1
          x = xh
          DO i = 1, ny
            y(i) = yloc(i)
          ENDDO
C  Fin du pas atteinte : retour au sousprogramme appelant
          IF (x.GE.pas) THEN
CD          WRITE(*,*) 'URK43 (',pas,') = ',IterOK,Iter999,' / ',Iter
            RETURN
          ENDIF
C  Estimation du sous-pas d'integration suivant
          hopt = h*MIN(hoptsec*(eps/err)**hc43pui,hoptmax)
C  Preparation des valeurs pour le sous-pas suivant
C  Recuperation de k1 = dydx1 = dydx (!)
          DO i = 1, ny
            dydx1(i) = dydx(i)
          ENDDO
c	ELSE IF (err.gt.(1.D5*eps)) then
c	 check=.true.
c	 return
        ELSE
C  Decoupage du sous-pas d'integration courant
          hopt = h*MAX(hoptsec*(eps/err)**hc43pui,hoptmin)
        ENDIF
        h = MIN(pas-x,hopt)
      GOTO 10
C  ------------

      RETURN
      END

