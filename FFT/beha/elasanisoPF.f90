!=======================================================================
!>           FONCTION ELASANISO for Phase-Field
!!-----------------------------------------------------------------------
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (tableau reel(6)) tenseur des deformations
!! \param[in]  DSTRAN  (tableau reel(6)) inccremenet tenseur des deformations
!! \param[in]  PROPS   (tableau reel(2)) coefficients (lambda, mu)
!!
!! \param[inout]  STATEV (tableau reel(2)) damage and history fields
!!
!! \param[out] STRESS    (tableau reel(6)) tenseur des contraintes
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!! NOTE: this subroutine requires another subroutine "rs" for computing the eigen*
!=================================================================================
SUBROUTINE elasanisoPF( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                        RPL, DDSDDT, DRPLDE, DRPLDT,&
                        STRAN, DSTRAN, TIME, DTIME,&
                        TEMP, DTEMP, PREDEF, DPRED,&
                        CMNAME, NDI, NSHR, NTENS, NSTATV,&
                        PROPS, NPROPS, COORDS,&
                        DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                        NOEL, NPT, LAYER, KSPT, KSTEP, KINC )


  implicit none

!     Arguments de l'interface umat
!--------------------------------------------------------
      CHARACTER*16  CMNAME
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
!--------------------------------------------------------

  double precision,dimension(6)    :: Def, sig1, def1
  double precision,dimension(3)    :: Def_sw
  double precision,dimension(3)    :: e1,e2,e3,v1,v2,v3
  double precision                 :: psiP
  integer                          :: i
  integer ( kind = 4 )             :: ierr
  !get the current strain
  Def = STRAN+DSTRAN
   
  !thermal expansion - orthotropic
  Def_sw(1:3) = PROPS(10:12) * DTEMP

  !Calcul de e1 (e1 normalisé) 
  e1 = PROPS(13:15)/(norm2(PROPS(13:15)))

  !Calcul de e2 (deuxième vecteur local) e2 = e2 - (e2.e1)e1, normalisé
  e2 = PROPS(16:18) - (dot_product(PROPS(16:18),e1) * e1) 
  e2 = e2 / norm2(e2)

  !Calcul de e3 (troisieme vecteur local)
  e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
  e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
  e3(3) = e1(1)*e2(2) - e1(2)*e2(1)

  !Passage de la deformation dans la base locale (def1)
  def1(1) = (           Def(1)*e1(1) + 0.5_8*Def(4)*e1(2) + 0.5_8*Def(5)*e1(3))* e1(1) &
          + (0.5_8*Def(4)*e1(1) +            Def(2)*e1(2) + 0.5_8*Def(6)*e1(3))* e1(2) & 
          + (0.5_8*Def(5)*e1(1) + 0.5_8*Def(6)*e1(2) +            Def(3)*e1(3))* e1(3)

  def1(2) = (           Def(1)*e2(1) + 0.5_8*Def(4)*e2(2) + 0.5_8*Def(5)*e2(3))* e2(1) &
          + (0.5_8*Def(4)*e2(1) +            Def(2)*e2(2) + 0.5_8*Def(6)*e2(3))* e2(2) & 
          + (0.5_8*Def(5)*e2(1) + 0.5_8*Def(6)*e2(2) +            Def(3)*e2(3))* e2(3)

  def1(3) = (           Def(1)*e3(1) + 0.5_8*Def(4)*e3(2) + 0.5_8*Def(5)*e3(3))* e3(1) &
          + (0.5_8*Def(4)*e3(1) +            Def(2)*e3(2) + 0.5_8*Def(6)*e3(3))* e3(2) & 
          + (0.5_8*Def(5)*e3(1) + 0.5_8*Def(6)*e3(2) +            Def(3)*e3(3))* e3(3)

  def1(4) = (           Def(1)*e1(1) + 0.5_8*Def(4)*e1(2) + 0.5_8*Def(5)*e1(3))* e2(1) &
          + (0.5_8*Def(4)*e1(1) +            Def(2)*e1(2) + 0.5_8*Def(6)*e1(3))* e2(2) & 
          + (0.5_8*Def(5)*e1(1) + 0.5_8*Def(6)*e1(2) +            Def(3)*e1(3))* e2(3) &
          + (           Def(1)*e2(1) + 0.5_8*Def(4)*e2(2) + 0.5_8*Def(5)*e2(3))* e1(1) &
          + (0.5_8*Def(4)*e2(1) +            Def(2)*e2(2) + 0.5_8*Def(6)*e2(3))* e1(2) & 
          + (0.5_8*Def(5)*e2(1) + 0.5_8*Def(6)*e2(2) +            Def(3)*e2(3))* e1(3)

  def1(5) = (           Def(1)*e1(1) + 0.5_8*Def(4)*e1(2) + 0.5_8*Def(5)*e1(3))* e3(1) &
          + (0.5_8*Def(4)*e1(1) +            Def(2)*e1(2) + 0.5_8*Def(6)*e1(3))* e3(2) & 
          + (0.5_8*Def(5)*e1(1) + 0.5_8*Def(6)*e1(2) +            Def(3)*e1(3))* e3(3) &
          + (           Def(1)*e3(1) + 0.5_8*Def(4)*e3(2) + 0.5_8*Def(5)*e3(3))* e1(1) &
          + (0.5_8*Def(4)*e3(1) +            Def(2)*e3(2) + 0.5_8*Def(6)*e3(3))* e1(2) & 
          + (0.5_8*Def(5)*e3(1) + 0.5_8*Def(6)*e3(2) +            Def(3)*e3(3))* e1(3)

  def1(6) = (           Def(1)*e2(1) + 0.5_8*Def(4)*e2(2) + 0.5_8*Def(5)*e2(3))* e3(1) &
          + (0.5_8*Def(4)*e2(1) +            Def(2)*e2(2) + 0.5_8*Def(6)*e2(3))* e3(2) & 
          + (0.5_8*Def(5)*e2(1) + 0.5_8*Def(6)*e2(2) +            Def(3)*e2(3))* e3(3) &
          + (           Def(1)*e3(1) + 0.5_8*Def(4)*e3(2) + 0.5_8*Def(5)*e3(3))* e2(1) &
          + (0.5_8*Def(4)*e3(1) +            Def(2)*e3(2) + 0.5_8*Def(6)*e3(3))* e2(2) & 
          + (0.5_8*Def(5)*e3(1) + 0.5_8*Def(6)*e3(2) +            Def(3)*e3(3))* e2(3)

  !include thermal expansion to total strain (in local coords)
  def1(1:3) = def1(1:3) - Def_sw


  !compute the stress
  sig1(1) = PROPS(1)*def1(1) + PROPS(2)*def1(2) + PROPS(3)*def1(3) 
  sig1(2) = PROPS(2)*def1(1) + PROPS(4)*def1(2) + PROPS(5)*def1(3) 
  sig1(3) = PROPS(3)*def1(1) + PROPS(5)*def1(2) + PROPS(6)*def1(3) 
  sig1(4) = PROPS(7)*def1(4) 
  sig1(5) = PROPS(8)*def1(5) 
  sig1(6) = PROPS(9)*def1(6) 

  !Retour de la contrainte dans la base macrosopique
  v1(1) = e1(1)
  v1(2) = e2(1)
  v1(3) = e3(1)

  v2(1) = e1(2)
  v2(2) = e2(2)
  v2(3) = e3(2)

  v3(1) = e1(3)
  v3(2) = e2(3)
  v3(3) = e3(3)

  STRESS(1) =  (sig1(1)*v1(1) + sig1(4)*v1(2) + sig1(5)*v1(3))* v1(1) &
          + (sig1(4)*v1(1) + sig1(2)*v1(2) + sig1(6)*v1(3))* v1(2) & 
          + (sig1(5)*v1(1) + sig1(6)*v1(2) + sig1(3)*v1(3))* v1(3)

  STRESS(2) =  (sig1(1)*v2(1) + sig1(4)*v2(2) + sig1(5)*v2(3))* v2(1) &
          + (sig1(4)*v2(1) + sig1(2)*v2(2) + sig1(6)*v2(3))* v2(2) & 
          + (sig1(5)*v2(1) + sig1(6)*v2(2) + sig1(3)*v2(3))* v2(3)

  STRESS(3) =  (sig1(1)*v3(1) + sig1(4)*v3(2) + sig1(5)*v3(3))* v3(1) &
          + (sig1(4)*v3(1) + sig1(2)*v3(2) + sig1(6)*v3(3))* v3(2) & 
          + (sig1(5)*v3(1) + sig1(6)*v3(2) + sig1(3)*v3(3))* v3(3)

  STRESS(4) =  (sig1(1)*v1(1) + sig1(4)*v1(2) + sig1(5)*v1(3))* v2(1) &
          + (sig1(4)*v1(1) + sig1(2)*v1(2) + sig1(6)*v1(3))* v2(2) & 
          + (sig1(5)*v1(1) + sig1(6)*v1(2) + sig1(3)*v1(3))* v2(3) &
          + (sig1(1)*v2(1) + sig1(4)*v2(2) + sig1(5)*v2(3))* v1(1) &
          + (sig1(4)*v2(1) + sig1(2)*v2(2) + sig1(6)*v2(3))* v1(2) & 
          + (sig1(5)*v2(1) + sig1(6)*v2(2) + sig1(3)*v2(3))* v1(3)

  STRESS(5) =  (sig1(1)*v1(1) + sig1(4)*v1(2) + sig1(5)*v1(3))* v3(1) &
          + (sig1(4)*v1(1) + sig1(2)*v1(2) + sig1(6)*v1(3))* v3(2) & 
          + (sig1(5)*v1(1) + sig1(6)*v1(2) + sig1(3)*v1(3))* v3(3) &
          + (sig1(1)*v3(1) + sig1(4)*v3(2) + sig1(5)*v3(3))* v1(1) &
          + (sig1(4)*v3(1) + sig1(2)*v3(2) + sig1(6)*v3(3))* v1(2) & 
          + (sig1(5)*v3(1) + sig1(6)*v3(2) + sig1(3)*v3(3))* v1(3)
 
  STRESS(6) =  (sig1(1)*v2(1) + sig1(4)*v2(2) + sig1(5)*v2(3))* v3(1) &
          + (sig1(4)*v2(1) + sig1(2)*v2(2) + sig1(6)*v2(3))* v3(2) & 
          + (sig1(5)*v2(1) + sig1(6)*v2(2) + sig1(3)*v2(3))* v3(3) & 
          + (sig1(1)*v3(1) + sig1(4)*v3(2) + sig1(5)*v3(3))* v2(1) &
          + (sig1(4)*v3(1) + sig1(2)*v3(2) + sig1(6)*v3(3))* v2(2) & 
          + (sig1(5)*v3(1) + sig1(6)*v3(2) + sig1(3)*v3(3))* v2(3)

  STRESS(4)=STRESS(4)/2._8
  STRESS(5)=STRESS(5)/2._8
  STRESS(6)=STRESS(6)/2._8


  !compute the history field
  psiP = 0.

  !update the history field
  STATEV(2) = psiP
 !TODO
 !   1. in mat.xml, change param values
 !                  for Fibres, "Coeff_composite" --> need to be changed
 !   2. in src/user_pfczm/amitex_user....f90, check the index for gc, lc for fibres
 !   3. in umat (this file), PARAM(i), check the indices to be consistent with mat.xml
 
end subroutine elasanisoPF

