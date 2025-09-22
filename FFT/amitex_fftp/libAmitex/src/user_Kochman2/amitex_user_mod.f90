!===============================================================================
!
!       MODULE AMITEX_USER_MOD : 
!>
!>      Module for the USER defined procedure :
!>
!>                  resolution_user() in resolution_user_mod.f90
!>                  OR
!>                  before_unpas_user, after_unpas_user in standard_user_mod.f90
!!
!!
!!     TODO : EXTENSION TO MULTIPLE FIELDS
!!
!===============================================================================
module amitex_user_mod


  use ISO_FORTRAN_ENV

! MPI AND 2DECOMP MODULES
!------------------------
  use MPI             
  use decomp_2d
  use decomp_2d_fft

! ALL AMITEX MODULES
!-----------------------
  use algo_functions_mod  ! complete list of amitex modules
  use amitex_mod
  !!use amitex_user_mod     ! this file
  use error_mod
  use field_mod
  use green_mod
  use io2_amitex_mod
  use io_amitex_mod    
  use linear_mod
  use loading_mod
  use material_mod
  use non_local_mod
  use param_algo_mod
  !!use resolution_mod       these modules use amitex_user_mod
  !!use resolution_user_mod  
  use sortie_std_mod
  !!use standard_user_mod  
  
#ifdef OPENMP
  use omp_lib
#endif

  implicit none

  private

! Declare the procedure used in amitex_fftp before the resolution (REQUIRED) 
! Initialize global user-variables

  public :: init_user_variables

! Declare global user-variables as public

  public :: Eta,Eta0,df_dEta                       
  public :: EtaF, FREQ2                          
  public :: alpha,Kinetic,Neta   
                       
! Declare global procedures as public (used in resolution_user_mod or standard_user_mod)

  public :: updateIntVar_JB, get_Eta, get_df_dEta,&
            get_interface_energy_coef, get_kinetic_coef

!!------------------------------------------------------------------------------
!>                                                                 PUBLIC FIELDS 

  !> real space (ntot,ncomp)
  real(mytype),allocatable,dimension(:,:)        :: Eta,Eta0,df_dEta    

  !> Fourier space ((nx,ny,nz,ncomp)
  complex(mytype),allocatable,dimension(:,:,:,:) :: EtaF    

  real(mytype),allocatable,dimension(:,:,:)      :: FREQ2               


!!------------------------------------------------------------------------------
!>                                                                  COEFFICIENTS 
 
  real(mytype),allocatable,dimension(:)          :: alpha         !< cefficient d energy d'interface / gradient
  real(mytype),allocatable,dimension(:)          :: Kinetic       !< coefficient cinetique

!!------------------------------------------------------------------------------
!>                                                              NUMBER OF FIELDS 
 
  integer                                        :: Neta 





contains


!!------------------------------------------------------------------------------
!>                                         ALLOCATE AND INITIALIZE PUBLIC FIELDS 

subroutine init_user_variables()

  implicit none
  integer                                     :: alloc_stat           !< erreur lors d'une allocation memoire
  integer                                     :: i
  character(len=200)                          :: tmp_char2,tmp_char1,tmp_char   !< variable caractere temporaire

  !! get Neta number of phase fields from UMAT : /4 (voir martensite_Kochmann2016.f90 : 4 varint / champ)
  do i=1,size(mattotP)
    !number of internal variables
    Neta=(MattotP(i)%NvarInt)/4
  end do
  if (nrank==0) write(OUTPUT_UNIT,*) "Number of phases Neta = ",Neta
  if(Neta .NE. 1) call amitex_abort("Neta different from 1 NOT YET IMPLEMENTED (init_user_variables)",2)
  
  !! allocate and initialize FIELDS
  allocate(Eta(xsize(1)*xsize(2)*xsize(3),Neta),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 1)",2)

  allocate(Eta0(xsize(1)*xsize(2)*xsize(3),Neta),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 2)",2)

  allocate(df_dEta(xsize(1)*xsize(2)*xsize(3),Neta),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 3)",2)

  allocate(EtaF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),Neta),&
           stat=alloc_stat)   
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 4)",2)

  allocate(FREQ2(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 4)",2)

  Eta=0.
  Eta0=0.
  df_dEta=0.
  EtaF=0.
  FREQ2 = FREQ(:,:,:,1)**2 + FREQ(:,:,:,2)**2 + FREQ(:,:,:,3)**2


  !! update interface and Kinetic parameters associated with each phase field
  allocate(alpha(1:Neta),stat=alloc_stat)
  allocate(Kinetic(1:Neta),stat=alloc_stat)

  !! Initialization Eta,Eta0 : get Eta0 from UMAT initialisation
  call get_Eta(eta0,Neta)
  eta = eta0

  !! Sortie vtk des champs de phase initiaux
  write(tmp_char1,"(I4)") 0
  do i=1,Neta
    write(tmp_char2,"(I4)") i
    tmp_char="_eta"//trim(adjustl(tmp_char2))//"_initial_"//trim(adjustl(tmp_char1))
    call print_field_vtk(eta0(:,i),trim(fic_vtk)//trim(tmp_char)//".vtk","eta"//trim(adjustl(tmp_char2)))
  end do


  !! mise a jour des variables internes dans l'UMAT
  !! ==============================================
  call updateIntVar_JB(eta,Neta)  !eta->VarInt(end); VarInt->VarInt0


end subroutine init_user_variables




!!==========================================================================================================
!!==========================================================================================================
!                                 subroutines for resolutionPF_JB
!!==========================================================================================================
!!==========================================================================================================

!!==========================================================================================================
!!==========================================================================================================
!> updateIntVar_JB(eta,Neta) : mise a jour de la variable interne UMAT
!!                             associee aux valeurs des champ de phases
!!                             TODO : compatibilite avec plus de 1 champ a verifier/optimiser
!!==========================================================================================================
subroutine updateIntVar_JB(eta,Neta)
  !!update the internal variables and assign the phase field into one of the internal variables
  ! param [in]
  !       - eta(:,i) : phase field indice i
  implicit none
  integer,intent(in) :: Neta
  real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),Neta),intent(in) :: eta
  integer(INT64) :: i,j,k,l,m, nVarInt, n
#ifdef OPENMP
    integer  :: nbthread 

!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif

  do i=1,size(mattotP)
    !number of internal variables
    nVarInt=MattotP(i)%NvarInt
    !indice minimum de zone
    l=1
    do j=1,size(MattotP(i)%zone(:,1))
#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m)
        !$OMP DO
#endif
      do k=l,MattotP(i)%zone(j,1)
        !linear index of voxel position
        m = MattotP(i)%pos(k)
        do n=1,Neta
          !update internal variable by phase field variable (eta(n)->VarInt(end-2*n+1))
          MattotP(i)%VarInt(nVarInt+4*(n-Neta)-3,k) = eta(m,n)
        end do
        !update internal variable at t-1 (VarInt->VarInt0)
        MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
      end do
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
      l=MattotP(i)%zone(j,1)+1
    end do
  end do
end subroutine updateIntVar_JB

!!==========================================================================================================
!!==========================================================================================================
!> get_Eta(Eta,Neta) :  recuperation dans l UMAT de Eta
!!                              TODO : compatibilite avec plus de 1 champ a verifier/optimiser
!!==========================================================================================================
subroutine get_Eta(Eta,Neta)
  ! param [out]
  !       - Eta(:,i) : phase field indice i
  implicit none
  integer,intent(in) :: Neta
  real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),Neta),intent(out) :: Eta
  integer(INT64) :: i,j,k,l,m, nVarInt, n
#ifdef OPENMP
    integer  :: nbthread
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif

  !get Eta(i) by internal variables
  !------------------------------------------------
  do i=1,size(mattotP)
    !number of internal variables
    nVarInt=MattotP(i)%NvarInt
    !indice minimum de zone
    l=1
    do j=1,size(MattotP(i)%zone(:,1))
#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m)
        !$OMP DO
#endif
      do k=l,MattotP(i)%zone(j,1)
        !linear index of voxel position
        m = MattotP(i)%pos(k)
        do n=1,Neta
          Eta(m,n) = MattotP(i)%VarInt0(nVarInt+4*(n-Neta)-3,k)
          !Eta(m,n) = MattotP(i)%VarInt(nVarInt+4*(n-Neta)-3,k)
        end do
        !update internal variable at t-1 (VarInt->VarInt0)
        MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
      end do
      !write(*,*) sum(df_dEta(:,1))
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
      l=MattotP(i)%zone(j,1)+1
    end do
  end do

end subroutine get_Eta

!!==========================================================================================================
!!==========================================================================================================
!> get_df_dEta(df_dETA,Neta) :  recuperation dans l UMAT de la force motrice
!!                              associee a chacun des champs de phase
!!                              TODO : compatibilite avec plus de 1 champ a verifier/optimiser
!!==========================================================================================================
subroutine get_df_dEta(df_dETA,Neta)
  ! param [out]
  !       - df_dETA(:,i) : phase field indice i
  implicit none
  integer,intent(in) :: Neta
  real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),Neta),intent(out) :: df_dEta
  integer(INT64) :: i,j,k,l,m, nVarInt, n
#ifdef OPENMP
    integer  :: nbthread
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif

  !get df/dEta(i) by internal variables
  !------------------------------------------------
  do i=1,size(mattotP)
    !number of internal variables
    nVarInt=MattotP(i)%NvarInt
    !indice minimum de zone
    l=1
    do j=1,size(MattotP(i)%zone(:,1))
#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m)
        !$OMP DO
#endif
      do k=l,MattotP(i)%zone(j,1)
        !linear index of voxel position
        m = MattotP(i)%pos(k)
        do n=1,Neta
          df_dEta(m,n) = MattotP(i)%VarInt(nVarInt+4*(n-Neta)-2,k)
        end do
        !update internal variable at t-1 (VarInt->VarInt0)
        MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
      end do
      !write(*,*) sum(df_dEta(:,1))
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
      l=MattotP(i)%zone(j,1)+1
    end do
  end do

end subroutine get_df_dEta

!!==========================================================================================================
!!==========================================================================================================
!> get_interface_energy_coef(alpha,Neta) : recuperation dans l UMAT du coefficient d energie d interface
!!                                         / terme de gradient associee a chacun des champs de phase
!!                                         TODO : peut etre generalise a un tenseur
!!                                         TODO : peut varier en fonction de la position / temps
!!                                         TODO : compatibilite avec plus de 1 champ a verifier/optimiser
!!==========================================================================================================
subroutine get_interface_energy_coef(alpha,Neta)
  ! param [out]
  !       - alpha : scalar
  integer,intent(in) :: Neta
  real(mytype),dimension(Neta),intent(out) :: alpha
  integer(INT64) :: i,j,k,l,m, nVarInt, n
#ifdef OPENMP
    integer  :: nbthread
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif

  !get df/dEta(i) by internal variables
  !------------------------------------------------
  do i=1,size(mattotP)
    !number of internal variables
    nVarInt=MattotP(i)%NvarInt
    !indice minimum de zone
    l=1
    do j=1,size(MattotP(i)%zone(:,1))
#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m)
        !$OMP DO
#endif
      do k=l,MattotP(i)%zone(j,1)
        !linear index of voxel position
        m = MattotP(i)%pos(k)
        do n=1,Neta
          alpha(n) = MattotP(i)%VarInt(nVarInt+4*(n-Neta)-1,k)
        end do
        !update internal variable at t-1 (VarInt->VarInt0)
        MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
      end do
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
      l=MattotP(i)%zone(j,1)+1
    end do
  end do

end subroutine get_interface_energy_coef

!!==========================================================================================================
!!==========================================================================================================
!> get_kinetic_coef(L,Neta) : recuperation dans l UMAT du coefficient cinetique
!!                            associe a chacun des champs de phase
!!                            TODO : peut varier en fonction de la position / temps
!!                            TODO : compatibilite avec plus de 1 champ a verifier/optimiser
!!==========================================================================================================
subroutine get_kinetic_coef(Kinetic,Neta)
!! to comment

  implicit none
  integer,intent(in) :: Neta
  real(mytype),dimension(Neta),intent(out) :: Kinetic
  integer(INT64) :: i,j,k,l,m, nVarInt, n
#ifdef OPENMP
    integer  :: nbthread
!$OMP PARALLEL
  nbthread = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif

  !get df/dEta(i) by internal variables
  !------------------------------------------------
  do i=1,size(mattotP)
    !number of internal variables
    nVarInt=MattotP(i)%NvarInt
    !indice minimum de zone
    l=1
    do j=1,size(MattotP(i)%zone(:,1))
#ifdef OPENMP
        !$OMP PARALLEL NUM_THREADS(nbthread) private(m)
        !$OMP DO
#endif
      do k=l,MattotP(i)%zone(j,1)
        !linear index of voxel position
        m = MattotP(i)%pos(k)
        do n=1,Neta
          Kinetic(n) = MattotP(i)%VarInt(nVarInt+4*(n-Neta),k)
        end do
        !update internal variable at t-1 (VarInt->VarInt0)
        MattotP(i)%VarInt0(:,k) = MattotP(i)%VarInt(:,k)
      end do
#ifdef OPENMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
      l=MattotP(i)%zone(j,1)+1
    end do
  end do

end subroutine get_kinetic_coef









end module amitex_user_mod
