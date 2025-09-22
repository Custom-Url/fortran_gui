!=======================================================================
!!                      FONCTION graingrowth_10grains_3
!!
!!----------------------------------------------------------------------
!! modele de champ de phase pour ###
!!
!!----------------------------------------------------------------------
!! \param[in]  PROPS   (tableau reel(20)) :
!!
!!        PROPS(1:2)          2 coeff (alpha,kinetic for grain 1)
!!        PROPS(3:4)          2 coeff (alpha,kinetic for grain 2)
!!        .
!!        .
!!        PROPS(19:20)          2 coeff (alpha,kinetic for grain 10)
!!
!! \param[out] STRESS  (tableau reel(6)) tenseur des contraintes ici = 1
!!
!! \Intvar[in][out]   STATEV    (tableau reel(20)) variables internes, 
!!
!!        STATEV(1:10)  (out) force motrice df/dEta_i pour champ de phase "i" (sans la contribution gradient)
!!        STATEV(11:20)  (in)  valeur champ de phase "i" mise a jour depuis user_graingrowth3
!!
!! Notation CAST3M (Voigt, ordre 11 22 33 12 13 23)
!=======================================================================

SUBROUTINE graingrowth_10grains_3( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,&
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

integer                                       :: i

double precision                              :: a, b, c
double precision                              :: f_homo


double precision                              :: sum_eta_squared

!=====================================================================
!                                                             MECHANIC

STRESS=1.

!=====================================================================
!                                                          PHASE FIELD


! Ginzburg Landau Type
! f(eta) = a*eta*eta*(1 - eta)**2
a = 0.14_8
b = 12.42_8
c = 12.28_8
f_homo = 1._8

! for multiple eta variables
! not optimized
sum_eta_squared = 0
do i=1,10
  sum_eta_squared = sum_eta_squared + STATEV(i+10)*STATEV(i+10)
end do

do i=1,10

  ! for independant eta variables
  !STATEV(i) =  f_homo * STATEV(i+10)*(a + STATEV(i+10)*STATEV(i+10)*(-b &
  !                                                 + c*STATEV(i+10)*STATEV(i+10)))

  ! for multiple eta variables with interactions
  STATEV(i) =  f_homo * STATEV(i+10)*(a - b*STATEV(i+10)*STATEV(i+10) &
                                       + c*sum_eta_squared*sum_eta_squared)

end do

end subroutine graingrowth_10grains_3
