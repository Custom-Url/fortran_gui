!===============================================================================
!
!       MODULE AMITEX_USER_MOD (minimalist white module in src)
!>
!>   - Declaration and initialization of USER defined variables 
!>   - User-defined procedures uses in resolution_user_mod or standard_user_mod 
!>
!!
!!
!===============================================================================
module amitex_user_mod

  use ISO_FORTRAN_ENV

  use decomp_2d

  use amitex_mod
  use error_mod
  use simu_mod,  only : SIMU_AMITEX,SIMU,& ! type
                        init_simu_command_line,nullify_pointers_simu, associate_pointers_simu ! functions
  use param_algo_mod, only : user_param  ! pointer 

  implicit none

  private

  public :: init_user_variables
  
  public :: SIMUloc  
  
  type(SIMU_AMITEX), dimension(1),target :: SIMUloc

contains


!!------------------------------------------------------------------------------
!>                                         ALLOCATE AND INITIALIZE PUBLIC FIELDS 

subroutine init_user_variables()

  implicit none
  
  character(len=200) :: version = "user_biscale2_0.0.0"  
  character(len=1000):: cmd_file         !< Command file name (from user_param%p_string(1))
  
  !!  3 variables inutiles dans init_simu_command_line (car presence cmd_file)
  integer            :: coeff2print =0   !< indice du coeff a sortir (ligne de commande) 
  integer            :: varint2print=0   !< indice de varint a sortir (ligne de commande) 
  logical            :: help=.false.     !< test pour avoir l'aide et sortir du code

  integer            :: i
  integer            :: err = 0

  !  Load local simulation in SIMUloc from cmd_file
  !------------------------------------------------
  cmd_file = user_param%p_string(1)
  call init_simu_command_line(SIMUloc(1),version, help,coeff2print,varint2print,&
                 trim(cmd_file),Igrid=2)
                 
  !rename simulation for output
  SIMUloc(1)%simu_name="LOCAL"

  !  Check loading (times and discretization)
  !-------------------------------------------------
  if (size(SIMU(1)%load) /=  size(SIMUloc(1)%load)) then 
    err = 1 
  else
    do i = 1, size(SIMU(1)%load)
      if (SIMU(1)%load(i)%Nincr /= SIMUloc(1)%load(i)%Nincr) err =1
      if (maxval(abs(SIMU(1)%load(i)%time(1:) - SIMUloc(1)%load(i)%time(1:)) / SIMU(1)%load(i)%time(1:)) > 1e-6) err = 1
      if (SIMU(1)%load(i)%Discretization /= SIMUloc(1)%load(i)%Discretization) err = 1
    end do
  end if


  if (err == 1) then
     call amitex_abort("ERROR : user_biscale2/init_user_variables : loading time discretization (local /= global)", 2, 0)
  end if

  ! Come back to SIMU(1)
  !-------------------------------------------------
  call nullify_pointers_simu()
  call associate_pointers_simu(SIMU,1)


end subroutine init_user_variables

end module amitex_user_mod
