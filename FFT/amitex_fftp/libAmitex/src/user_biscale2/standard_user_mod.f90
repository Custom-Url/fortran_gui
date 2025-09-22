!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE STANDARD_USER_MOD (minimalist white module in src)
!>      
!>
!!
!!
module standard_user_mod

  use ISO_FORTRAN_ENV

  use decomp_2d
  use mpi

  use amitex_mod
  use error_mod

  use amitex_user_mod, only : SIMUloc ! variables
  use param_algo_mod,  only : algo_param ! pointers
  use loading_mod,     only : LOCALLOAD,& ! types
                              local_load, load, extract  ! pointers
  use simu_mod, only        : SIMU,& ! variables
                              nullify_pointers_simu, associate_pointers_simu ! functions
  use NL_base_mod, only  : unpas_NL_base, log_final_times, log_after_unpas,& ! functions
                              times_b ! pointers
  use field_mod, only       : times_f !pointers
  use material_mod, only    : times_m, mattotP !pointers
  use green_mod, only       : times_g !pointers
  use io2_amitex_mod, only  : times_io,&  !pointers
                              writeVTK !functions
  use sortie_std_mod, only  : times_s,& !pointers
                              sortie_std ! functions
  

  implicit none

  private 

  public :: before_unpas_user, after_unpas_user

contains

!===========================================================================================================
!!==========================================================================================================
!> before_unpas_user : write what must be done just BEFORE unpas (using the standard algorithm)
!!                       
!!                       
!!==========================================================================================================
subroutine before_unpas_user(load_n,load_incr, ind_tps)
  
  implicit none

  integer,intent(in)     :: load_n,load_incr, ind_tps

 ! trick to avoid gcc-warning
 integer                 :: bidon
 bidon = load_n + load_incr + ind_tps

end subroutine before_unpas_user


!===========================================================================================================
!!==========================================================================================================
!> after_unpas_user : write what must be done just AFTER unpas (using the standard algorithm)
!!                       
!!                       
!!==========================================================================================================
subroutine after_unpas_user(load_n,load_incr, ind_tps)
  
  implicit none

  integer,intent(in)       :: load_n,load_incr, ind_tps
  
  integer                  :: I0=1, i, ierror
  
  integer                  :: nIT,nIt0=0, nIttot=0   
  type(LOCALLOAD)          :: locload        
  
  real(mytype)             :: t_fft0,t_ifft0,t_behavior0,t_green0,&
                              t_crit0,t_wvtk0,t_wstd0,t_unpas0
  logical                  :: testLaminate, test_load_interruption
  real(mytype),dimension(algo_param%nTensDef) :: defMoy, sigMoy      !<  def et sig moyens fin d'increment
  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy, FluxDMoy !<  grad et flux moyens fin d'increment
                                                                     !<  pas utilisees / sorties de sortie_std

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
  !---------------------------------------------------------------------
  !                                                      INITIALISATIONS
  !

  nIt = 0

  !---------------------------------------------------------------------
  !                                          RECUPERATION CHARGEMENT LOC
  !                                 (seul locload%t1 nous interesse ici)
  !
  locload = local_load ! V0 : chargement GLOB = LOC
   
  !---------------------------------------------------------------------
  !                                           TRANSFERT SIMU GLOB -> LOC
  call nullify_pointers_simu()
  call associate_pointers_simu(SIMUloc,I0)


  !---------------------------------------------------------------------
  !                                          INITIALISATIONS AVANT UNPAS
  !                                            (voir resolution_NL_base)
  !

  t_fft0      = times_f%fft
  t_ifft0     = times_f%ifft
  t_behavior0 = times_m%behavior
  t_green0    = times_g%apply
  t_crit0     = times_g%crit
  t_wvtk0     = times_io%wvtk
  t_wstd0     = times_s%wtot
  t_unpas0    = times_b%unpas
  
  testLaminate = .false. !TODO : mettre au claire sortie laminate/composite...
  test_load_interruption = .false.
  
  !---------------------------------------------------------------------
  !                                                      LANCEMENT UNPAS

  call unpas_NL_base(load_n,load_incr,ind_tps,nIt,locload)
                
  !---------------------------------------------------------------------
  !                                                         APRES  UNPAS
  !
  !                                            (voir resolution_NL_base)
  !

  nittot = nittot + nit
  
  !--- SORTIES STD
  if(load_incr .eq. 0) nIt0 = nIt
  if(load_incr .eq. 1) nIt = nIt + nIt0

  if(load_incr/=0) then  ! Pas de sortie std pour les pas de temps fictifs
     call sortie_std(local_load%t_load(0),ind_tps,nIt,&
                    sigMoy,defMoy,FluxDMoy, GradQDMoy)
  end if
  
  if(load_incr .eq. 1) nIt = nIt - nIt0     

  !--- SORTIES VTK
  if(load_incr/=0) then
  if(extract%tpsVTK(ind_tps) .OR. test_load_interruption)then 
       call writeVTK(ind_tps)
  end if
  end if
  
  !--- LOG APRES UN PAS
  call log_after_unpas(load_incr,ind_tps, nIt,&
               t_fft0,t_ifft0,t_behavior0,t_green0,t_crit0,t_wvtk0,t_wstd0,t_unpas0,&
               testLaminate)

  !--- LOG FIN DE CALCUL
  if (load_n == size(load) .AND. load_incr == load(load_n)%Nincr) then
      call log_final_times(nittot)
  end if 
  
  !---  MAJ CONTRAINTE DE CAUCHY, DES VARIABLES INTERNES
  do i=1,size(MattotP)
       MattotP(i)%VarInt0=MattotP(i)%VarInt
  end do
  Sig0 = Sig

  !---------------------------------------------------------------------
  !                                           TRANSFERT SIMU LOC -> GLOB
  call nullify_pointers_simu()
  call associate_pointers_simu(SIMU,I0)


end subroutine after_unpas_user


end module standard_user_mod



