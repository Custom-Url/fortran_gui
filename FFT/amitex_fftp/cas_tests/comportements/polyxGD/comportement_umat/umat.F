C UMAT      SOURCE    CHAT      11/03/16    21:32:49     6902         
C-----------------------------------------------------------------------
C>
C>    DESCRIPTION FONCTIONNELLE :
C>    -------------------------
C>
C>    Version du module externe UMAT integree au code pour eviter tout
C>    probleme d'edition de lien et pour permettre le passage de tests
C>    de non regression.
C>
C>    En utilisation normale, ce module doit etre surcharge par celui
C>    de meme nom fourni par l'utilisateur.
C>
C>    Si cette version de UMAT est appelee, cela signifie qu'il n'y a
C>    pas eu surcharge par le module utilisateur
C>    ATTENTION : si l'utilisateur oublie de surcharger et affecte a sa
C>    loi le meme numero que l'un de ceux pris en charge ci-apres, il y
C>    a risque d'erreur non controlee a l'execution.
C>
C>    CONVENTION : on autorise la programmation de l'utilisateur a se
C>    servir de KINC comme code retour, de la meme maniere que l'on
C>    s'en sert ci-dessous.\n Regles a respecter :
C>       1. Pas d'initialisation superflue de KINC en entrant dans UMAT.
C>        KINC est initialise a 1 par WKUMA0, avant l'appel a UMAT.
C>       2. En cas d'erreur, KINC est affecte d'une valeur differente de 1
C>
C>    Dans ce dernier cas, l'utilisateur a la possibilite de completer
C>     le message d'erreur affiche (de numero 965) en surchargeant le
C>     sousprogramme UMAERR (umaerr.eso).
C>
C>    DESCRIPTION DE L'INTERFACE :
C>    --------------------------
C>
C>    \param[in,out] STRESS: REAL*8(NTENS), tenseur des contraintes
C>                       - En entree : tenseur des contraintes a t0
C>                       - En sortie : tenseur des contraintes a t0+dt
C>                       .
C>                   Dans le cas des grandes transformations, il
C>                   s'agit des contraintes de Cauchy ("actuelles").
C>    \param[in,out] STATEV: REAL*8(*), variables internes
C>                    - En entree : variables internes a t0
C>                    - En sortie : variables internes a t0+dt
C>
C>    \param[out] DDSDDE: REAL*8(NTENS,NTENS), matrice jacobienne du
C>                      modele (matrice de Hooke tangente) a t0+dt
C>                      Sortie facultative, non exploitee par CAST3M
C>                      pour l'instant
C>
C>    \param[in,out] SSE: REAL*8, energie de deformation elastique
C>    \param[in,out]         SPD: REAL*8, dissipation plastique
C>    \param[in,out]         SCD: REAL*8, dissipation visqueuse
C>                      Valeurs a t0 en entree, a t0+dt en sortie
C>                      Entrees/sorties facultatives, non exploitees par
C>                      CAST3M pour l'instant
C>
C>    \param[out]    RPL: REAL*8, puissance calorifique volumique degagee
C>                      par le travail mecanique, a t0+dt
C>    \param[out]         DDSDDT: REAL*8(NTENS),derivee du tenseur des contraintes
C>                      par rapport a la temperature, a t0+dt
C>    \param[out]         DRPLDE: REAL*8(NTENS), derivees de RPL par rapport aux
C>                      composantes du tenseur des deformations, a t0+dt
C>    \param[out]         DRPLDT: REAL*8, derivee de RPL par rapport a la
C>                      temperature, a t0+dt
C>                      Sorties facultatives, non exploitees par CAST3M
C>                      pour l'instant
C>                      Ces sorties sont prevues pour un couplage fort
C>                      entre thermique et mecanique
C>
C>    \param[in]      STRAN: REAL*8(NTENS), tenseur des deformations
C>                      totales a t0
C>    \param[in]       DSTRAN: REAL*8(NTENS), tenseur des increments de
C>                      deformation totale par rapport a l'etat
C>                      de reference a t0
C>
C>    \param[in]     TIME: REAL*8(2), TIME(1) = 0
C>                                 TIME(2) = t0
C>    \param[in]         DTIME: REAL*8,    DTIME   = dt
C>                      t0 : precedent instant d'equilibre atteint
C>                      dt : nouveau pas de temps propose par PASAPAS
C>                      pour atteindre l'equilibre avec l'increment de
C>                      deformation totale impose (DSTRAN)
C>                      Remarque : les notions attachees a ces entrees
C>                      sont differentes lorsque l'on adhere a ABAQUS
C>                      (dans ce cas la variable TIME(1) a son utilite)
C>
C>    \param[in]     TEMP: REAL*8, temperature a t0
C>    \param[in]         DTEMP: REAL*8, increment de temperature a t0+dt
C>    \param[in]         PREDEF: REAL*8(*), vecteur des parametres externes de
C>                      la loi de comportement, valeurs a t0
C>    \param[in]         DPRED: REAL*8(*), increments des parametres externes
C>                      a t0+dt
C>
C>    \param[in]     CMNAME: CHARACTER*16, identifiant de la loi de
C>                      comportement.\n
C>                      On conserve le type 'chaine de caracteres' pour
C>                      l'identifiant de la loi, afin de preserver la
C>                      compatibilite avec ABAQUS.
C>                      Dans le cas d'une adherence a CAST3M, la loi est
C>                      identifiee par le numero qui lui a ete attribue:
C>                      argument 'NUME_LOI' de l'operateur MODE.
C>                      Par convention, ce numero est encode dans les 4
C>                      derniers caracteres de la chaine, et doit etre
C>                      recupere dans UMAT par une instruction du type
C>                      K4ILOI = CMNAME(13:16)
C>                      avec K4ILOI variable locale de type CHARACTER*4
C>
C>    \param[in]   NDI: INTEGER definissant le type de calcul CAST3M
C>                    - =  2 : tridimensionnel (3D)
C>                    - =  0 : axisymetrie (2D)
C>                    - = -1 : deformations planes (2D)
C>                    - = -2 : contraintes planes (2D)
C>                    - = -3 : deformations planes generalisees (2D)
C>                    - =  3 a 11 : unidimensionnel plan (massif 1D)
C>                    - = 12 a 14 : unidimensionnel axisymetrique (massif 1D)
C>                    - = 15 : unidimensionnel spherique (massif 1D)
C>
C>    \param[in]     NSHR: INTEGER, nombre de composantes extradiagonales
C>                      du tenseur des contraintes
C>                      Entree non active dans le cas d'une adherence
C>                      a CAST3M
C>
C>    \param[in]     NTENS: INTEGER, nombre de composantes du tenseur des
C>                      contraintes
C>    \param[in]         NSTATV: INTEGER, nombre de variables internes
C>
C>    \param[in]     PROPS: REAL*8(NPROPS), vecteur des constantes
C>                      de materiau.\n
C>                      N.B. Il peut s'agir de coefficients constants
C>                      intervenant dans des formules analytiques
C>                      decrivant les evolutions des proprietes du
C>                      materiau en fonction de la temperature et/ou
C>                      d'autres parametres externes.
C>                      OU il peut s'agir plus simplement de proprietes
C>                      de materiau constantes.\n
C>    \param[in]         NPROPS: INTEGER, nombre de constantes de materiau
C>
C>    \param[in]     COORDS: REAL*8(3), coordonnees cartesiennes du point
C>                      d'integration courant
C>
C>    \param[in]      DROT: REAL*8(3,3)
C>                      - Dans le cas d'une adherence a CAST3M\n
C>                      Matrice de passage du repere local de l'element
C>                      fini massif au repere general du maillage
C>                      DROT(l,c) : composante l du vecteur c de la base
C>                      locale dans la base generale
C>                      Entree non active pour les elements finis non
C>                      massifs
C>                      - Dans le cas d'une adherence a ABAQUS\n
C>                      Matrice d'increments de rotation
C>                      Cette matrice decrit la rotation sur le pas de
C>                      temps de la base dans laquelle sont exprimes
C>                      les tenseurs de contraintes et de deformations
C>
C>    \param[out]    PNEWDT: REAL*8, rapport entre le nouveau pas de temps
C>                      suggere et le pas de temps donne en entree
C>                      (NEWDT = PNEWDT * DTIME)
C>
C>    \param[in]      CELENT: REAL*8, longueur caracteristique de l'element
C>                      Dans le cas d'une adherence a CAST3M, cette
C>                      longueur caracteristique est determinee par
C>                      LOCARA comme la distance maximale entre deux
C>                      noeuds de l'element.
C>
C>    \param[in]     DFGRD0: REAL*8(3,3), tenseur gradient de deplacement
C>                      a t0
C>    \param[in]         DFGRD1: REAL*8(3,3), tenseur gradient de deplacement
C>                      a t0+dt
C>
C>    \param[in]     NOEL: INTEGER, numero de l'element courant
C>    \param[in]         NPT: INTEGER, numero du point d'integration courant
C>
C>    \param[in]   LAYER: INTEGER, numero de couche pour des coques
C>                      composites ou des solides multi-couches
C>    \param[in]   KSPT: INTEGER, numero de section dans la couche
C>                      courante
C>                      Entrees non actives dans le cas d'une adherence
C>                      a CAST3M
C>
C>    \param[in]     KSTEP: INTEGER
C>    \param[in]     KINC: INTEGER
C>                      Entrees n'ayant pas de sens dans le cas
C>                      d'une adherence a CAST3M
C>
C-----------------------------------------------------------------------
  

      SUBROUTINE UMATBCCGD ( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     &                  RPL, DDSDDT, DRPLDE, DRPLDT,
     &                  STRAN, DSTRAN, TIME, DTIME,
     &                  TEMP, DTEMP, PREDEF, DPRED,
     &                  CMNAME, NDI, NSHR, NTENS, NSTATV,
     &                  PROPS, NPROPS, COORDS,
     &                  DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,
     &                  NOEL, NPT, LAYER, KSPT, KSTEP, KINC )

      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL*8(A-H,O-Z)

C     Arguments de l'interface
C
      CHARACTER*16  CMNAME
C
      INTEGER       NDI, NSHR, NTENS, NSTATV, NPROPS,
     &              NOEL, NPT, LAYER, KSPT, KSTEP, KINC
C
      REAL*8        STRESS(NTENS), STATEV(*),
     &              DDSDDE(NTENS,NTENS),
     &              SSE, SPD, SCD,
     &              RPL, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT,
     &              STRAN(NTENS), DSTRAN(NTENS),
     &              TIME(2), DTIME,
     &              TEMP, DTEMP, PREDEF(*), DPRED(*),
     &              PROPS(NPROPS),
     &              COORDS(3),
     &              DROT(3,3),
     &              PNEWDT,
     &              CELENT,
     &              DFGRD0(3,3), DFGRD1(3,3)

      CHARACTER*4   K4ILOI

      
c      K4ILOI = CMNAME(13:16)
      KINC = 1

c      IF (K4ILOI.EQ.'  90') THEN
         nLdC = 90
         CALL UMATGD(  nLdC,    NDI, TIME(2),  DTIME,   TEMP, DTEMP,
     &        DFGRD0, DFGRD1,   PROPS, NPROPS, PREDEF, DPRED,
     &        STATEV, NSTATV,  STRESS,  NTENS,   KINC)
c      ELSEIF (K4ILOI.EQ.'  91') THEN
c         nLdC = 91
c         CALL UMATGD(  nLdC,    NDI, TIME(2),  DTIME,   TEMP, DTEMP,
c     &        DFGRD0, DFGRD1,   PROPS, NPROPS, PREDEF, DPRED,
c     &        STATEV, NSTATV,  STRESS,  NTENS,   KINC)

C
C=======================================================================
C     MODELE CONTRAINTE IMPOSEE
C=======================================================================
C
c      ELSE IF (K4ILOI.EQ.'   0') THEN
         
c         call contrainte_imposee(STRESS,PROPS)

C
C=======================================================================
C     MODELE ELASTIQUE LINEAIRE ISOTROPE
C=======================================================================
C
c      ELSE IF (K4ILOI.EQ.'   1') THEN
        
c        call elasiso(STRESS,STRAN+DSTRAN,PROPS)

C
C=======================================================================
C     MODELE THERMOELASTIQUE LINEAIRE ISOTROPE
C=======================================================================
C     
c      ELSEIF (K4ILOI.EQ.'   2') THEN
         
c         call thermoelasiso(STRESS,STRAN+DSTRAN,PROPS,TEMP)

C     
C=======================================================================
C     MODELE VISCOELASTIQUE LINEAIRE MAXWELL GENERALISE POUR k ET mu
C=======================================================================
C     
c      ELSE IF (K4ILOI.EQ.'  40') THEN
         
c         call viscoelas_maxwell(STRESS, STATEV, DSTRAN, DTIME, 
c     &        NTENS, NSTATV, PROPS, NPROPS)         
         
C=======================================================================
C     Modele externe non reconnu
C=======================================================================

c      ELSE

c         KINC = -2

c      ENDIF

      RETURN
      END

