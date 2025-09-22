!=======================================================================
!!                      FONCTION MARTENSITE_KOCHMANN2016_polygrain
!!
!!----------------------------------------------------------------------
!! modele de champ de phase pour transformations martensitique 2D
!! application du modele 2D presente dans 
!! [Kochmann2016] J. Kochmann et al. / comput. Methods Appl. Mech. Engrg. 305 (2016) 89-110
!! resultats p. 102 Fig. 6
!!
!!----------------------------------------------------------------------
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) incremenet tenseur des deformations
!!
!! \param[in]  PROPS   (tableau reel(30)) :
!!
!!        PROPS(1:9)            9 coefficients elastiques-orthotrope phase principale / matrice
!!        PROPS(10:15)          6 coordonnees vecteurs de bases du grain, e1 = PROPS(10:12) ; e2 = PROPS(13:15)
!!        PROPS(16:24)          9 coefficients elastiques-orthotrope champ de phase "1"
!!        PROPS(25:30)          6 coefficients tenseur deformation residuel / "bain strain" champ de phase "1"
!!        ...
!!        PROPS({1:9} + i*15)   9 coefficients elastiques-orthotrope champ de phase "i"
!!        PROPS({10:15}+i*15)   6 coefficients tenseur deformation residuel / "bain strain" champ de phase "i"
!!
!! \param[out] STRESS  (tableau reel(6)) tenseur des contraintes
!!
!! \Intvar[in][out]   STATEV    (tableau reel(4)) variables internes, 
!!                              1 -> mise a jour depuis resolutionPF_mod.f90
!!                              2-4 -> envoyees dans resolutionPF_mod.f90
!!
!!        STATEV(1 + (i-1)*4)   valeur champ de phase "i" mise a jour depuis resolutionPF_mod.f90
!!        STATEV(2 + (i-1)*4)   force motrice df/dEta_i pour champ de phase "i" (sans la contribution gradient)
!!        STATEV(3 + (i-1)*4)   coefficient energie d interface / gradient pour champ de phase "i"
!!        STATEV(4 + (i-1)*4)   coefficient cinetique pour champ de phase "i"
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================

SUBROUTINE martensite_Kochmann2016_polygrain( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC )


  implicit none

!     Arguments de l'interface umat
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
!-----------------------------------------------------------------------

  double precision,dimension(6)                 :: Def
  double precision,dimension(6)                 :: ER
  double precision,dimension(6)                 :: Def_elast
  double precision,dimension(6)                 :: S
  double precision,dimension(6)                 :: dER_dEta
  double precision,dimension(9)                 :: CE
  double precision,dimension(9)                 :: dCE_dEta
  integer                                       :: i

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
  ! associated with phase field "i" STATEV(i)
  !CE = PROPS(1:9) + STATEV(1)*(PROPS(16:24) - PROPS(1:9))
  ! here from [Kochmann2016] p.99 expres.(43)
  CE = PROPS(1:9) + STATEV(1)*STATEV(1)*(PROPS(16:24) - PROPS(1:9))

  ! self deformation intrinsic to phase tranformation 
  ! associated with phase field "i" STATEV(i)
  !ER = STATEV(1)*PROPS(31:36)
  ! compute bain / residual strain in case of martensite transformations
  ! here from [Kochmann2016] p.99 expres.(43) 
  ER(1) = STATEV(1)*STATEV(1)*PROPS(25)
  ER(2) = ER(1) 
  ER(3) = 0._8
  ER(4) = STATEV(1)*PROPS(28)*2._8
  ER(5) = 0._8
  ER(6) = 0._8

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
  ! here from [Kochmann2016] p.97 expres.(31) 
  ! elastic energy = 0.5*(Def-ER({Etas}))*CE({Etas}):(Def-ER({Etas}))
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

  !> compute total driving force for phase field "i" implementation 
  !> in local base (more convenient)

  ! dER_dEta = dER/dETA_i associated with phase field "i"
  !dER_dEta = PROPS(25:30)
  dER_dEta(1) = 2._8*STATEV(1)*PROPS(25)
  dER_dEta(2) = dER_dEta(1)
  dER_dEta(3) = 0._8
  dER_dEta(4) = PROPS(28)*2._8
  dER_dEta(5) = 0._8
  dER_dEta(6) = 0._8

  ! dCE_dEta = dCE/dEta_i associated with phase field "i"
  !dCE_dEta = (PROPS(16:24) - PROPS(1:9))
  dCE_dEta = 2._8*STATEV(1)*(PROPS(16:24) - PROPS(1:9))

  ! compute elastic contribution for driving force
  ! partial derivative of elastic energy 
  ! here from [Kochmann2016] p.97 expres.(31) 
  ! elastic energy = 0.5*(Def-ER({Etas}))*CE({Etas}):(Def-ER({Etas}))
  f_elast = 1._8 ! ! elastic constant term for driving force / should be equal to 1.
  !f_elast = 0._8 !> to switch off elastic effects 
  STATEV(2) = -sum(dER_dEta*sig1)
  call voigt_product_orthotropic(dCE_dEta,def_elast,S)
  STATEV(2) = STATEV(2) + sum(Def_elast*S)
  call voigt_product_orthotropic(CE,dER_dEta,S)
  STATEV(2) = STATEV(2) - sum(Def_elast*S)
  STATEV(2) = STATEV(2)*0.5_8*f_elast

  ! homogeneous contribution for driving force df/dEta
  ! partial derivative of homogeneous energy
  ! here from [Kochmann2016] p.99 expres. (44)
  ! with f(eta) = (3*a*a - 1 2*eta*eta)(1 - eta *eta)**2/(3*a*a - 1)
  a = 0.35_8           !<  ! from [Kochmann2016] p.100 table 1
  !f_homo = 0.00049042_8  !<  ! consideration p.103 part 5.3 is not respected 
                               !<  ! I think there is an error in article
  f_homo = 0.049042_8     !<  ! Psy_H0 in [Kochmann2016] p.100 table 1 is kept
  STATEV(2) =  STATEV(2) - f_homo * 12._8*STATEV(1)*(1._8 - STATEV(1)*STATEV(1)) &
                                  *(a*a - STATEV(1)*STATEV(1)) &
                                  /(3*a*a - 1._8)


  !> other phase field coefficients sent to resolutionPF_mod.f90

  ! gradient coefficient 
  ! "alpha" in resolutionPF_mod.f90
  ! here Psy_G0 in [Kochmann2016] p.95 expres. 20 or p.96 expres. 28
  !STATEV(3) = 1.0904_8  !< ! Psy_G0 in [Kochmann2016] p.100 table 1
  STATEV(3) = 0.010904_8 !< ! Psy_G0 in [Kochmann2016] p.100 table 1 with consideration p.103 part 5.3
                              !< ! if only this parameter smaller => interfaces thiner

  ! kinetic coefficient
  ! "Kinetic" in resolutionPF_mod.f90
  ! here mobility m_0 in [Kochmann2016] p.95 expres. 20
  ! TODO => mobility m_0 can be activated only when Stress is above a certain Threshold value 
  !      => see [Kochmann2016] part 3.3 p.31
  !      => thise value is not given in the article
  STATEV(4) = 5._8     !<  ! m_0 in [Kochmann2016] p.100 table 1

end subroutine martensite_Kochmann2016_polygrain

