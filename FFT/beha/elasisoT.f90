!=======================================================================
!>           ELASISO FUNCTION for Phase-Field
!!-----------------------------------------------------------------------
!!
!!-----------------------------------------------------------------------
!!
!! \param[in]  STRAN   (real array(6)) strain tensor
!! \param[in]  DSTRAN  (real array(6)) strain tensor increment
!! \param[in]  PROPS   (real array(2)) coefficients (lambda, mu)
!!
!! \param[inout]  STATEV (real array(2)) damage and history fields
!!
!! \param[out] STRESS    (real array(6)) stress tensor
!!
!! CAST3M notation (Voigt, ordre 11 22 33 12 13 23)
!! NOTE: this subroutine requires another subroutine "rs" for computing the eigen*
!=================================================================================
SUBROUTINE elasisoT( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
                       RPL, DDSDDT, DRPLDE, DRPLDT,&
                       STRAN, DSTRAN, TIME, DTIME,&
                       TEMP, DTEMP, PREDEF, DPRED,&
                       CMNAME, NDI, NSHR, NTENS, NSTATV,&
                       PROPS, NPROPS, COORDS,&
                       DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,&
                       NOEL, NPT, LAYER, KSPT, KSTEP, KINC )

  implicit none

  !     UMAT interface arguments
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

  double precision,dimension(6)    :: Def
  double precision                 :: Def_sw
  double precision                 :: Tr_Def,Mu2, Lambda, Mu
  integer                          :: i

  ! Material Properties
  Lambda = PROPS(1)
  Mu = PROPS(2)
  Mu2 = 2._8 * Mu

  !get the current strain
  Def = STRAN+DSTRAN
   
  !thermal expansion (assume constant CTE)
  Def_sw = PROPS(3) * DTEMP !/ 3
  STATEV = STATEV + Def_sw   !accumulated thermal strain

  Def(1:3) = Def(1:3) - STATEV

  !calculate stress 
  Tr_Def = Def(1) + Def(2) + Def(3)
  do i = 1,3
     STRESS(i) = Lambda * Tr_Def + Mu2 * Def(i)
  end do
  do i = 4,6
     STRESS(i) = Mu * Def(i)
  end do

  !print *, 'TEMP:', TEMP+DTEMP, 'eps_th:', STATEV
end subroutine elasisoT

