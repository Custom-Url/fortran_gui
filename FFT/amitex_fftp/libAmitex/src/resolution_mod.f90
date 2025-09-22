!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE RESOLUTION_MOD : 
!>
!> Solving the problem using an FFT-based method
!!
!!  
!!  Subroutines
!!
!! - resolution_NL_base : loop over time steps
!!                        call unpas_NL_base
!!                        but also (for user defined 'couplings') : 
!!                               - Nloc_call
!!                               OR
!!                               - before(after)_unpas_user
!!
!!
!!
!
module resolution_mod

  use ISO_FORTRAN_ENV

  use MPI             ! Loading modules to gain access to functions
  use decomp_2d, only : mytype, nrank

  use io2_amitex_mod
  use algo_functions_mod
  use param_algo_mod
  use loading_mod
  use material_mod
  use green_mod
  use error_mod
  use amitex_mod
  use field_mod
  use non_local_mod
  use sortie_std_mod
  use NL_base_mod
  use standard_user_mod
  use user_functions_mod

    private

  !> "Public" pointers to components of SIMU_AMITEX
  ! nothing

  !> "Public" variables used during initialization
  ! nothing

  !> Public types (for definition of SIMU_AMITEX)
  ! nothing

  !> Public functions
  public :: resolution_NL_base
 
contains

!===========================================================================================================
!!==========================================================================================================
!> resolution_NL_base: Basic scheme for solving a nonlinear problem
!!                     using Fast Fourier Transforms (FFT)
!!                     This procedure contains the resolution algorithm
!!                     and handles both mechanics and diffusion
!!
!! NOTATIONS (xml input and fields) Sigma and Def (HPP format):
!!                 CAST3M (Voigt notation + order 11 22 33 12 13 23)
!!                 PK1 and grad(u) (GDEF): order 11 22 33 12 13 23 21 31 32 without Voigt notation
!!
!!                             Variables in MATERIAL_MOD
!!--------------------------------------------------------
!!   LambdaMu0   Coefficients [Lambda0, mu0] of the reference medium
!!   C0          Stiffness matrix of the reference material
!!   k0D         Conductivity coefficient (diffusion)
!!               size (nVarD)
!!
!!   nmateriaux  Number of materials
!!
!!   VCinfo      Info related to laminate algorithm
!!
!!           Parallel fields accessible in AMITEX_MOD
!!--------------------------------------------------------
!! Description of fields in Fourier space:
!!             Z-pencils, (nx/2+1, ny, nz, 3, nVarD or NTensDef)
!!             3D indices: (fft_start(1):fft_end(1), fft_start(2):fft_end(2),
!!                         fft_start(3):fft_end(3), xx)
!!
!! Description of fields in real space:
!!             X-pencils, (nx, ny, nz, 3, nVarD) or (nx, ny, nz, NTensDef)
!!             1D index: (1:xsize(1)*xsize(2)*xsize(3), NTensDef or 3, nVarD)
!!
!!   Sig      Cauchy stress field
!!   PK1      First Piola-Kirchhoff stress field for large deformations
!!   Def      Deformation field in HPP and grad(u) in GDEF format
!!   Def0     Deformation field at previous time step
!!   Sig0     Cauchy stress field at previous time step
!!   SigF     Cauchy stress in HPP or Piola-Kirchhoff stress in GD,
!!            Fourier space
!!   DefF     Deformation in HPP or grad(u) in GD, Fourier space
!!            Z-pencil (nx/2+1, ny, nz, nTensDef)
!!   FluxD0   Flux vector field (:, 3, nVarD) (diffusion)
!!   FluxD    Flux vector field (:, 3, nVarD) (diffusion)
!!   GradQD0  Gradient field of the diffused variable (:, 3, nVarD) (diffusion)
!!   GradQD   Gradient field of the diffused variable (:, 3, nVarD) (diffusion)
!!   FluxDF   Flux vector field in Fourier space (diffusion)
!!   GradQDF  Gradient field in Fourier space
!!
!! Tables for convergence acceleration
!!                 for mechanics (ntot, NTens, 4)
!!                 for diffusion (ntot, 3*nVarD, 4)
!!   ACT3_R   Residuals for mechanics convergence acceleration
!!   ACT3_U   Increments in deformation fields for mechanics convergence
!!   ACT3_RD  Same as ACT3_R for diffusion
!!   ACT3_UD  Same as ACT3_U for diffusion
!!
!!                        Other variables in AMITEX_MOD
!!--------------------------------------------------------
!!   Flog       Logical units for log file output
!!   fic_vtk    Root name of the vtk or (mz)std output files
!!
!!                                Variables in GREEN_MOD
!!--------------------------------------------------------
!!   FREQ     Frequency array, Fourier space
!!            Z-pencil (nx/2+1, ny, nz, 3)
!!            3D indices: (fft_start(1):fft_end(1), fft_start(2):fft_end(2),
!!                        fft_start(3):fft_end(3), 3)
!!
!!                              Variables in LOADING_MOD
!!--------------------------------------------------------
!!   local_load%t1         Imposed mechanical loading at each iteration
!!   local_loadD%t1        Imposed diffusion loading at each iteration
!!
!!
!!                              Variables in NL_BASE_MOD
!!--------------------------------------------------------
!!   crit_b                Convergence criteria for resolution_NL_base
!!   times_b               Elapsed computation times for resolution_NL_base
!!
!!------------------------------------------------------------
!!                             OTHER DERIVED-TYPE VARIABLES
!!
!!      MattotP         Array of MATERIAL structures   (material_mod.f90)
!!            -> description of the material
!!
!!      load            Array of LOADING structures    (loading_mod.f90)
!!            -> description of the loading
!!
!!      initValExt      INITLOADEXT structure          (loading_mod.f90)
!!            -> description of the initial external loading
!!
!!      extract         PARAMEXTRACT structure         (loading_mod.f90)
!!            -> description of vtk outputs
!!
!!      algo_param      PARAM_ALGO structure           (param_algo_mod.f90)
!!            -> description of algorithm parameters
!!
!!      grid            GRID_DIM structure             (green_mod.f90)
!!            -> description of the FFT grid
!!
!!      VCinfo          VOXCOMP_INFO structure         (material_mod.f90)
!!            -> tracking information for composite voxel algorithms
!!
!!===============================================================================
subroutine resolution_NL_base()

        
!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================

  implicit none
!!------------------------------------------------------------------------------
!>                                                     PARAMETERS "PARALLELISM"

  integer :: ierror                                      !< error related to MPI function

!!------------------------------------------------------------------------------
!>                                                       PARAMETERS "LOADING"

  integer                                     :: ind_tps !< index of the current time step
  integer                                     :: load_n, load_incr
                                                         !< number of the partial loading
                                                         !< increment number of the loading

  logical      :: test_load_interruption  !< test for interruption of the current loading
  real(mytype) :: alpha_interruption      !< parameter to adjust loading in case of interruption

!!------------------------------------------------------------------------------
!>                                           PARAMETERS FOR ALGORITHM MONITORING

  integer               :: nIt                    !< number of iterations in one loading step
  integer               :: nIt0                   !< number of iterations in the "fictitious" time step
  integer               :: nitTot                 !< total number of iterations

!!------------------------------------------------------------------------------
!>                                  MONITORING PARAMETERS FOR LAMINATE ALGORITHM

  logical                       :: testLaminate,testLaminateloc !< true if composite laminate voxels are present
                                                                !< only useful for outputting laminate info

  !< laminate algorithm performance tracking (material_mod)
  !<               VCinfo%nSub_lam         : number of time step subdivisions
  !<               VCinfo%nIncr_lam        : number of sub time steps per subdivision
  !<               VCinfo%nIt_lam          : number of iterations per call to the laminate algorithm
  !<               VCinfo%nCall_lam        : number of calls to the laminate algorithm
  !<               VCinfo%nSub_lam_pas
  !<               VCinfo%nIncr_lam_pas
  !<               VCinfo%nIt_lam_pas
  !<               VCinfo%nCall_lam_pas

!!------------------------------------------------------------------------------
!>                                            PARAMETERS OF THE SOLVING ALGORITHM

  real(mytype),dimension(algo_param%nTensDef) :: defMoy, sigMoy      !< average stress, strain... at the end of increment
  real(mytype),dimension(algo_param%nTensDef) :: defMoy0, sigMoy0    !< average stress, strain... at the start of increment
                                                                     !   useful in user_load_interruption_average
  real(mytype),dimension(algo_param%nTensDef) :: defMoy00, sigMoy00  !< average stress, strain... at the start of partial loading

  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy, FluxDMoy     !< average gradient and flux at the end of increment
  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy0, FluxDMoy0   !< average gradient and flux at the start of increment
  real(mytype),dimension(3,algo_param%nVarD)  :: gradQDMoy00, FluxDMoy00 !< average gradient and flux at the start of partial loading

!!------------------------------------------------------------------------------
!>                                     PARAMETERS FOR EXECUTION TIME MONITORING

  !! double precision type instead of real_mytype (=> required by MPI_WTIME)
  double precision :: t1
  double precision :: t_fft0, t_ifft0, t_behavior0, t_green0, t_crit0, t_wvtk0,&
                      t_wstd0, t_unpas0
  !< t1          : reference time before calling functions
  !! t_xxx0      : cumulative time at the start of the step

!!------------------------------------------------------------------------------
!>                                                                        MISCELLANEOUS
  
  integer      :: i,j                      !< loop index 
  logical      :: tmp_bool
  real(mytype) :: tmp_real

!! Variables for restart test
!  character(len=250)   :: totofile
!  integer              :: totofid

!!==============================================================================
!!                                                               INITIALISATIONS
!!==============================================================================
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  times_b%startalgo = MPI_WTIME()

  nIttot= 0
  nIt0 = 0

  if (algo_param%Mechanics) then
    defMoy=0                 !< average stresses and strains
    sigMoy=0
    defMoy0 = 0
    sigMoy0 = 0
    defMoy00 = 0
    sigMoy00 = 0
  end if 
  if (algo_param%Diffusion) then  
    GradQDMoy = 0
    FluxDMoy = 0
    GradQDMoy0 = 0
    FluxDMoy0 = 0
    GradQDMoy00 = 0
    FluxDMoy00 = 0
  end if
  ind_Tps= 1               !< index of the first step

  t1 = 0                   !< exectuion time 

  testlaminateloc = .false.  !< test for prescence of "laminate" material
  do i=1,size(mattotP)
   if (trim(mattotP(i)%lawname) == "laminate") testlaminateloc = .true.
  end do
  call MPI_Allreduce(testlaminateloc, testlaminate, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)    

  !***************************************************************************************************
  !                                                                    START LOOP OVER PARTIAL LOADING
  BOUCLE_CHARGEMENTS:&
  do load_n=1,size(load)
     test_load_interruption = .false.
     alpha_interruption = 1.
     !************************************************************************************************
     !                                                 START LOOP OVER INCREMENTS OF A PARTIAL LOADING
     BOUCLE_INCREMENTS:&
     do load_incr = 0,load(load_n)%NIncr

        !================================================================================
        !                                                     TIME MEASUREMENT REFERENCES
        !

        t_fft0      = times_f%fft
        t_ifft0     = times_f%ifft
        t_behavior0 = times_m%behavior
        t_green0    = times_g%apply
        t_crit0     = times_g%crit
        t_wvtk0     = times_io%wvtk
        t_wstd0     = times_s%wtot
        t_unpas0    = times_b%unpas
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
        times_b%startstep = MPI_WTIME()
 
        !================================================================================
        !                                             UPDATE OF THE LOADING local_load(D)
        !
        !                             actualisation of local_load(D)%t0, t1, dt and t_load
        !

        if(load_incr == 0) ind_tps=ind_tps-1  ! we do not increment ind_tps for the fictitious time step
        
        call update_local_load(load_n,load_incr,defMoy00, sigMoy00, gradQDMoy, FluxDMoy)

        !================================================================================
        !                                                           CALL TO UNPAS_NL_BASE
        !
        !
 
        if (user_param%test .AND. trim(user_param%algo) == "standard") then
           call before_unpas_user(load_n,load_incr,ind_tps)
        end if

        if (algo_param%Nloc) then !---- Non-local mechanics explicit resolution
           call Nloc_call(load_n,load_incr,ind_tps,"befor",algo_param%nloc_explicit)
        end if
        
!print *,"BEFORE UNPAS_NL_BASE", nrank
        call unpas_NL_base(load_n,load_incr,ind_tps,nIt)

        if (user_param%test .AND. trim(user_param%algo) == "standard") then
           call after_unpas_user(load_n,load_incr,ind_tps)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierror) !necessary for non-local computation ?

        if (algo_param%Nloc) then !---- Non-local mehcanics explicit resolution
!print *,"BEFORE NLOC_CALL", nrank
           call Nloc_call(load_n,load_incr,ind_tps,"after",algo_param%nloc_explicit)
        end if

!print *,"BEFORE USER LOAD INTERRUPTION", nrank
        !================================================================================
        !                                                          USER LOAD INTERRUPTION 

        !--TODO : plan to retrieve averages from unpas_NL_base

        if (algo_param%Mechanics .AND. load(load_n)%user_interruption_flag) then
        if(algo_param%HPP) then
           call field_Mean(SIG,grid%ntot,algo_param%nTensDef,sigMoy)
        else
           call field_Mean(PK1,grid%ntot,algo_param%nTensDef,sigMoy)
        end if
        call field_Mean(Def,grid%ntot,algo_param%nTensDef,defMoy)
        
        print *, '>>> user_load_interruption_average called at step:', load_n

        !--- Stop criterion based on average stress-strain
        call user_load_interruption_average(&
                          tmp_bool,tmp_real,&
                          load_n,sigMoy0,sigMoy,defMoy0,defMoy,&
                          local_load%t0,local_load%t1,local_load%t_load,&
                          load(load_n)%user_interruption)
        test_load_interruption = test_load_interruption .OR. tmp_bool
        alpha_interruption = min(alpha_interruption, tmp_real)


        !--- Stop criterion based on internal variables
        do i=1,size(MattotP)
             call user_load_interruption_internal_variables(&
                          tmp_bool,tmp_real,&
                          load_n,MattotP(i)%VarInt0,MattotP(i)%VarInt,MattotP(i)%numM,&
                          load(load_n)%user_interruption) 
             test_load_interruption = test_load_interruption .OR. tmp_bool
             alpha_interruption = min(alpha_interruption, tmp_real)
        end do
        end if
        
!print *,"BEFORE INTERRUPTION ALPHA_INTERRUPTION", nrank
        !================================================================================
        !                                                 INTERRUPTION ALPHA_INTERRUPTION
        ! 
        !         we adjust SIG, DEF, mattot%Varint, local_load%t1 and local_load%t_load(0) 
        !
        !                    + UPDATE load()%time : 
        !                    current load: we freeze the time of the following increments
        !                    subsequent loads: we shift the times
        !
        ! TODO : to review or limit to an interruption without restart
 
        if (algo_param%Mechanics .AND. load(load_n)%user_interruption_flag) then       
        if (test_load_interruption .and. load_incr .NE. 0 .and. alpha_interruption .ne. 1.) then
           if (alpha_interruption > 1.) then
             call amitex_abort("Load Interruption (in user_functions_mod) with alpha > 1 ",2)
           else
             Sig  = Sig0 + alpha_interruption * (Sig - Sig0)
             Def  = Def0 + alpha_interruption * (Def - Def0)
             if (.not. algo_param%HPP) call SigToPK1(Sig,Def,PK1)
             do i=1,size(MattotP)
               MattotP(i)%VarInt = MattotP(i)%VarInt0 &
                                 + alpha_interruption*(MattotP(i)%VarInt-MattotP(i)%VarInt0) 
             end do
             local_load%t1(algo_param%nTensDef:) = local_load%t0(algo_param%nTensDef:) +&
               alpha_interruption*(local_load%t1(algo_param%nTensDef:)-local_load%t0(algo_param%nTensDef:))
             local_load%t_load(0) = local_load%t_load(-1) +&
                                    alpha_interruption*(local_load%t_load(0) - local_load%t_load(-1))
             
             ! Shift the time of subsequent loads
             do i=load_n+1,size(load)
                load(i)%time=load(i)%time - (load(load_n)%time(load(load_n)%Nincr) - local_load%t_load(0))
             end do 
             ! Freeze the rest of the current load (time)
             load(load_n)%time(load_incr:) = load(load_n)%time(load_incr)
 
           end if
        end if  
        end if

!------------------------------------------------------------------------------
!RESTART TEST: WRITE / READ ON NPROC FILES
!
!SAVE
!      write(totofile,*) nrank
!
!      open(newu0nit=totofid,file="/tmp/titi"//trim(adjustl(totofile)),form='UNFORMATTED',&
!             action='WRITE',status='REPLACE',position='REWIND',&
!             access='sequential')
!      write(totofid) FluxD,GradQD,MattotP
!      close(totofid)
!      call mpi_barrier(mpi_comm_world,ierror)
!
!CORRUPT THE SAVED VARIABLES
!      FluxD=0./0.
!      GradQD=0./0.
!      MattotP(1)%LibName=""
!      MattotP(1)%LibNameK=""
!
!READ BACK THE SAVED VARIABLES
!      open(newunit=totofid,file="/tmp/titi"//trim(adjustl(totofile)),form='UNFORMATTED',&
!             action='READ',position='REWIND',&
!             access='sequential')
!      read(totofid) FluxD,GradQD,MattotP
!      close(totofid)
!      call mpi_barrier(mpi_comm_world,ierror)
!-------------------------------------------------------------------------------------

!print *,"BEFORE STANDARD OUTPUT FILE", nrank
        !================================================================================
        !                                                   STANDARD OUTPUT FILE
        !
        ! average stresses and strains, standard deviations
        ! + check average values over the cell
        !
        ! \TODO Separate standard output and evaluation of averages
        !
        !

        !> for first loading step, we account for the iterations of the fictitious time step 
        if(load_incr .eq. 0) nIt0 = nIt
        if(load_incr .eq. 1) nIt = nIt + nIt0

        if(load_incr/=0) then  ! No standard output for fictitious time steps
            call sortie_std(local_load%t_load(0),ind_tps,nIt,&
                            sigMoy,defMoy,FluxDMoy, GradQDMoy)
        end if

        !> we substract nIt0 above for a first loading step
        if(load_incr .eq. 1) nIt = nIt - nIt0 
         
!print *,"BEFORE VTK OUTPUT FILE", nrank
        !================================================================================
        !                                                                 VTK OUTPUT FILE
        
        if(load_incr/=0) then ! No VTK output for ficticious time steps
           ! output on request or upon interupption
           if(extract%tpsVTK(ind_tps) .OR. test_load_interruption)then 

              call writeVTK(ind_tps)

           end if
        end if
        
!print *,"BEFORE CAUCHY STRESS UPDATE", nrank
        !================================================================================
        !                                  UPDATE OF CAUCHY STRESS AND INTERNAL VARIABLES
        !
        
        if (algo_param%Mechanics) then
          do i=1,size(MattotP)
             MattotP(i)%VarInt0=MattotP(i)%VarInt
          end do
          !Sig0 = Sig
          ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises
          do i=1,6
          do j= 1,size(Sig0,dim=1)
             Sig0(j,i) = Sig(j,i)
          end do
          end do
          sigMoy0 = sigMoy
          defMoy0 = defMoy
        end if 
        if (algo_param%Diffusion) then
          FluxD0 = FluxD
          GradQDMoy0 = GradQDMoy
          FluxDMoy0  = FluxDMoy 
        end if 

        !================================================================================
        !                       UPDATE OF AVERAGE VALUES AT THE BEGINNING OF PARTIAL LOAD
        !        
                   
        if (algo_param%Mechanics) then                  
        if(load_incr==load(load_n)%NIncr .OR. test_load_interruption)then  
            sigMoy00 = SigMoy
            defMoy00 = DefMoy  
        end if
        end if
        if (algo_param%Diffusion) then
        if(load_incr==load(load_n)%NIncr)then  
            GradQDMoy00 = GradQDMoy
            FluxDMoy00  = FluxDMoy 
        end if
        end if

        !================================================================================
        !                                                         UPDATE ITERATION COUNTS

        nitTot = nitTot + nIt                                           ! Increment total number of iterations

        VCinfo%nCall_lam_pas = VCinfo%nCall_lam_pas * (nIt + 1)         ! Increment number of calls to umatLaminate for this step
        VCinfo%nCall_lam     = VCinfo%nCall_lam + VCinfo%nCall_lam_pas ! Increment total number of calls to umatLaminate
        VCinfo%nIt_lam       = VCinfo%nIt_lam + VCinfo%nIt_lam_pas     ! Increment total number of iterations for the laminate algorithm
        VCinfo%nSub_lam      = VCinfo%nSub_lam + VCinfo%nSub_lam_pas   ! Increment total number of times it was
                                                                        ! necessary to subdivide the time step in laminate
        VCinfo%nIncr_lam     = VCinfo%nIncr_lam + VCinfo%nIncr_lam_pas ! Increment total number of sub time steps used for the laminate algorithm


        !================================================================================
        !                                                                      LOG OUTPUT

        call log_after_unpas(load_incr, ind_tps, nIt, &
                  t_fft0, t_ifft0, t_behavior0, t_green0, t_crit0, t_wvtk0, t_wstd0, t_unpas0, &
                  testLaminate)

        call check_amitex_abort(0) !necessary???

        !================================================================================
        !                                                                  LOAD INCREMENT
        !                                                          USER LOAD INTERRUPTION

        ind_tps = ind_tps + 1

        if (test_load_interruption .AND. load_incr .NE. 0) exit BOUCLE_INCREMENTS

!print *,"END INCREMENT LOOP", ind_tps
     end do BOUCLE_INCREMENTS
     !                                               END OF LOOP OVER INCREMENTS OF A PARTIAL LOADING
     !************************************************************************************************
  end do BOUCLE_CHARGEMENTS
  !                                                            END OF LOOP OVER PARTIAL LOADINGS
  !***************************************************************************************************

  !================================================================================
  !                                                                FINAL LOG OUTPUT

  call log_final_times(nitTot)

end subroutine resolution_NL_base

!==============================================================================
!      SUBROUTINE UPDATE_LOCAL_LOAD
!
!>  RETRIEVAL OF CURRENT LOAD (load_n, load_incr)
!!                                      -> local_load(D)%t0, t1 and %dt
!!                                      -> UPDATE local_load%t_load
!!  UPDATE local_load%t_load :   t(-2)=t(-1), t(-1)=t(0), t(0)=load(load_nb)%time(incr)
!!
!! REMINDER (see loading_mod.f90)
!! MECHANICS
!! local_load%t1(1:algo_param%nTensDef) :                       control variables
!! local_load%t1(algo_param%nTensDef+1:2*algo_param%nTensDef) : values associated with control
!! local_load%t1(2*algo_param%nTensDef+1) :                     temperature
!! local_load%t1(nb_param next indices) :                       external parameters (if any)
!! local_load%t1(27 next indices) :                             gradgradU components (if any)
!!
!! DIFFUSION
!! local_loadD%t1(1:3*algo_param%nVarD) :                        control variables
!! local_loadD%t1(3*algo_param%nVarD+1:6*algo_param%nVarD) :     values associated with control
!! local_loadD%t1(6*algo_param%nVarD+1) :                        temperature
!! local_loadD%t1(6*algo_param%nVarD+2:) :                       external parameters (if any)
!!
!!
!!  \param[in]    load_n, partial loading number
!!  \param[in]    load_incr, index of the calculation step
!!  \param[in]    defMoy     | average quantities at beginning of step
!!  \param[in]    sigMoy     |    -> useful in case of "hold" at change
!!  \param[in]    gradQDMoy  |       of partial loading
!!  \param[in]    FluxDMoy   |
!!
!------------------------------------------------------------------------------
subroutine update_local_load(load_n,load_incr,defMoy, sigMoy, gradQDMoy, FluxDMoy)

  implicit none

  integer, intent(in)                                     :: load_n, load_incr
                                                                       !< partial loading number
                                                                       !< loading increment number
                                                                       !< (load_incr = 0 : extra increment at the beginning of loading)

real(mytype), dimension(algo_param%nTensDef), intent(in) :: defMoy, sigMoy
                                                                       !< Average strain and stress at the beginning of the partial loading

real(mytype), dimension(3,algo_param%nVarD), intent(in)  :: gradQDMoy, FluxDMoy
                                                                       !< Useful in case of hold (e.g., creep)


  if (algo_param%Mechanics) then
     local_load%t0 = local_load%t1
     call get_current_loading(load_n,load_incr,defMoy, sigMoy, local_load%t_load,&
                              local_load%t0,local_load%t1)
     local_load%dt = local_load%t1 - local_load%t0
  end if

  if (algo_param%Diffusion) then
     local_loadD%t0 = local_loadD%t1
     call get_current_loadingD(load_n,load_incr,gradQDMoy, FluxDMoy, local_loadD%t_load,&
                               local_loadD%t0,local_loadD%t1)
     local_loadD%dt = local_loadD%t1 - local_loadD%t0
  end if 

end subroutine update_local_load
!============================================================================== 




end module resolution_mod
