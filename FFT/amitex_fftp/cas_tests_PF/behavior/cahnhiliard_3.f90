!=======================================================================
!!                      FONCTION CAHNHILLIARD_3
!!
!! MODELE DE CHAMP DE PHASE CAHN-HILLIARD (AVEC COUPLAGE MECA) 
!!                     INTERFACE "3" (NON-LOCAL)
!!       ASSOCIE A user_CahnHilliard3/non_local_user_mod.f90
!!
!!
!! J. Boisse (equations), L. Gelebart (amitex, non-local framework)
!!----------------------------------------------------------------------
!!
!!----------------------------------------------------------------------
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) increment tenseur des deformations
!!
!! \param[in]  PROPS   (tableau reel(30)) :
!!
!!        PROPS(1:9)            9 coefficients elastiques-orthotrope phase principale / matrice
!!        PROPS(10:15)          6 coordonnees vecteurs de bases du grain, e1 = PROPS(10:12) ; e2 = PROPS(13:15)
!!        PROPS(16:24)          9 coefficients elastiques-orthotrope champ de phase "1"
!!        PROPS(25:30)          6 coefficients tenseur deformation residuel / "bain strain" champ de phase "1"
!!        PROPS(31)             gradient coefficient 
!!                              "alpha" in user_CahnHilliard3/non_local_user_mod.f90
!!        PROPS(32)             kinetic coefficient
!!                             "Kinetic" in user_CahnHilliard3/non_local_user_mod.f90
!!
!! \param[out] STRESS  (tableau reel(6)) tenseur des contraintes
!!
!! \Intvar[in][out]   STATEV    (tableau reel(2)) variables internes, 
!!
!!        STATEV(1)   (out) force motrice df/dEta pour champ de phase (sans la contribution gradient)
!!                          -> Nloc(1)%Var(:,1) dans user_CahnHilliard3/non_local_user_mod.f90
!!        STATEV(2)   (in)  valeur champ de phase mise a jour depuis user_CahnHilliard3/non_local_user_mod.f90
!!                          -> GNloc(1)%Var(:,1) dans user_CahnHilliard3/non_local_user_mod.f90
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================

SUBROUTINE CahnHilliard_3( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC )


  implicit none

!     UMAT INTERFACE ARGUMENTS                            (DO NOT TOUCH)
!-----------------------------------------------------------------------
      CHARACTER*16 CMNAME
!
      INTEGER*8    NDI, NSHR, NTENS, NSTATV, NPROPS,&
                   NOEL, NPT, LAYER, KSPT, KSTEP, KINC
!
      REAL*8       STRESS(NTENS), STATEV(NSTATV),&
                   DDSDDE(NTENS,NTENS),&
                   SSE, SPD, SCD,&
                   RPL, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT,&
                   STRAN(NTENS), DSTRAN(NTENS),&
                   TIME(2), DTIME,&
                   TEMP, DTEMP, PREDEF(*), DPRED(*),&
                   PROPS(NPROPS),&
                   COORDS(3),&
                   DROT(3,3),&
                   PNEWDT,&
                   CELENT,&
                   DFGRD0(3,3), DFGRD1(3,3)

!                                                         USER VARIABLES
!-----------------------------------------------------------------------

  double precision,dimension(6)                 :: Def
  double precision,dimension(6)                 :: ER
  double precision,dimension(6)                 :: Def_elast
  double precision,dimension(6)                 :: S
  double precision,dimension(6)                 :: dER_dEta
  double precision,dimension(9)                 :: CE
  double precision,dimension(9)                 :: dCE_dEta

  double precision                              :: a
  double precision                              :: f_homo
  double precision                              :: f_elast

  integer,dimension(6),parameter                :: voigt_strain_coefficients=(/1,1,1,2,2,2/)

  double precision,dimension(3)                 :: e1,e2,e3
  double precision,dimension(3)                 :: v1,v2,v3

  double precision,dimension(6)                 :: def1
  double precision,dimension(6)                 :: sig1

  !=====================================================================
  !                                                             MECHANIC

  ! update total strain
  ! here = elastic strain + bain strain (ER)
  Def = STRAN+DSTRAN

  ! update the elastic matrix as a fonction of phase fields
  ! associated with phase field STATEV(2)
  CE = PROPS(1:9) + STATEV(2)*(PROPS(16:24) - PROPS(1:9))

  ! self deformation intrinsic to phase tranformation 
  ! associated with phase field STATEV(2)
  ER = STATEV(2)*PROPS(25:30)

  !! rotation matrix defined by vectors [e1,e2,e3] of local base expressed in the global base
  !! (e1,e2,e3) base should be orhtonormal

  ! Calcul de e1 (e1 normalise) 
  e1 = PROPS(10:12)/(norm2(PROPS(10:12)))

  ! Calcul de e2 (deuxieme vecteur de la base locale) e2 = e2 - (e2.e1)e1, normalise
  e2 = PROPS(13:15) - (dot_product(PROPS(13:15),e1) * e1) 
  e2 = e2 / norm2(e2)

  ! Calcul de e3 (troisieme vecteur de la base local)
  e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
  e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
  e3(3) = e1(1)*e2(2) - e1(2)*e2(1)

  ! transform total strain from global (def) to local base (def1)
  call rotate_symetric_voigt_vector(def/dble(voigt_strain_coefficients), def1, e1, e2, e3)
  def1 = def1*dble(voigt_strain_coefficients)

  ! compute elastic strain in the local base
  Def_elast = def1 - ER

  ! compute stress in the local base => sig1
  call voigt_product_orthotropic(CE,Def_elast,sig1)

  !! rotation matrix defined by vectors [v1,v2,v3] of global base expressed in the local base
  ! [v1,v2,v3] = transposition of rotation matrix [e1,e2,e3]
  v1(1) = e1(1)
  v1(2) = e2(1)
  v1(3) = e3(1)
  v2(1) = e1(2)
  v2(2) = e2(2)
  v2(3) = e3(2)
  v3(1) = e1(3)
  v3(2) = e2(3)
  v3(3) = e3(3)

  ! transform stress from local (sig1) to global base (STRESS)
  call rotate_symetric_voigt_vector(sig1, STRESS, v1, v2, v3)

  !=====================================================================
  !                                                          PHASE FIELD

  !> compute total driving force for phase field implementation 
  !> in local base (more convenient)

  ! dER_dEta = dER/dETA_i associated with phase field
  dER_dEta = PROPS(25:30)

  ! dCE_dEta = dCE/dEta_i associated with phase field 
  dCE_dEta = (PROPS(16:24) - PROPS(1:9))

  ! compute elastic contribution for driving force
  ! partial derivative of elastic energy 
  ! here from [Kochmann2016] p.97 expres.(31) 
  ! elastic energy = 0.5*(Def-ER({Etas}))*CE({Etas}):(Def-ER({Etas}))
  f_elast = 100._8 ! ! elastic constant term for driving force / should be equal to 1.
  !f_elast = 0._8 !> to switch off elastic effects 
  STATEV(1) = -sum(dER_dEta*sig1)
  call voigt_product_orthotropic(dCE_dEta,def_elast,S)
  STATEV(1) = STATEV(1) + sum(Def_elast*S)
  call voigt_product_orthotropic(CE,dER_dEta,S)
  STATEV(1) = STATEV(1) - sum(Def_elast*S)
  STATEV(1) = STATEV(1)*0.5_8*f_elast

  ! homogeneous contribution for driving force df/dEta
  ! partial derivative of homogeneous energy
  ! f(eta) = a*eta*eta*(1 - eta)**2
  a = 16._8
  f_homo = 1._8
  STATEV(1) =  STATEV(1) + f_homo * 4._8 * a * STATEV(2)*(STATEV(2) - 0.5_8) &
                                  *(STATEV(2) - 1._8)

end subroutine CahnHilliard_3
