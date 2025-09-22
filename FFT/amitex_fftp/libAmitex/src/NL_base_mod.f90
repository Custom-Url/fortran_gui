!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE NL_BASE_MOD 
!>
!> SOLVE 'ONE STEP' WITH THE BASIC SCHEME (WITH CONVERGENCE ACCELARATION)
!!
!!  
!!  Subroutines
!! - unpas_NL_base : solves one step (unpas) with the basic Scheme (Moulinec-Suquet)
!!                   + CAST3M convergence acceleration (ACT3)
!!
!! - log_after_unpas, log_final_times, write_crit_log : 
!!                    Write in "log" file and/or standard output
!
module NL_base_mod

  use ISO_FORTRAN_ENV

  use MPI             ! load modules to access MPO functions
  use decomp_2d, only : mytype, nrank

  use io_amitex_mod,  only: write_stdout0 ! function
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
!  use user_functions_mod

  private
   
  !> Public pointers to components of SIMU_AMITEX
  public ::  crit_b, times_b

  !> Public variables used during initialisation
  ! nothing

  !> Public types (for SIMU_AMITEX definition)
  public ::  CRIT_BASE, TIMES_BASE

  !> Public functions
  public :: unpas_NL_base, log_after_unpas, log_final_times, write_crit_log

  !!---------------------------------------------------------------
  !>                                   TYPES CRIT_BASE & TIMES_BASE
  !>                     gather criteria used by resolution_NL_base
  !> 
  type CRIT_BASE
                                                        !< MECHANICAL
    real(mytype)                 :: eq                  !< equilibrium criterion
    real(mytype)                 :: SigMoy,DefMoy       !< means criterion 
    real(mytype)                 :: Cptb                !< compbaitibilty criterion
  
                                                                  !< DIFFUSION
    real(mytype),allocatable,dimension(:)  :: eqD                 !< equilibrium criterion
    real(mytype),allocatable,dimension(:)  :: FluxDMoy,GradQDMoy  !< means criterion
    real(mytype),allocatable,dimension(:)  :: CptbD               !< compatibility criterion

  end type CRIT_BASE
  
  type TIMES_BASE
     double precision       :: startalgo    = 0  
     double precision       :: startstep    = 0  
     double precision       :: unpas    = 0       
  end type TIMES_BASE

  type(CRIT_BASE), pointer       :: crit_b 
  type(TIMES_BASE), pointer      :: times_b 
  
contains

!===========================================================================================================
!!==========================================================================================================
!> unpas_NL_base: Basic scheme for solving a nonlinear problem
!!                using Fast Fourier Transforms
!!                This procedure contains the resolution algorithm
!!                and handles both mechanics and diffusion
!!
!!        Solves the problem over a single step, called by resolution_NL_base
!!        which iterates over successive steps
!!
!! NOTATIONS (XML input and fields) Sigma and Def (HPP):
!!        CAST3M (Voigt + order 11 22 33 12 13 23)
!!        PK1 and grad(u) (GDEF): order 11 22 33 12 13 23 21 31 32 with no notation
!!
!!
!!        LIST OF MODIFIED VARIABLES (To be completed!)
!!--------------------------------------------------------
!!   VCinfo%nIt_lam_pas, nSub_lam_pas, nIncr_lam_pas, nCall_lam_pas
!!   local_load%t_load, UPDATE of the time triplet [t-2, t-1, t]
!!             used by resolution_NL_base in sortie_std(t_load(0)...)
!!             and in the ALPHA_INTERRUPTION section
!!
!!                     Variables in MATERIAL_MOD
!!--------------------------------------------------------
!!   LambdaMu0   Coefficients [Lambda0, mu0] of the reference medium
!!   C0          Stiffness matrix of the reference material
!!   k0D         Conductivity coefficient (diffusion)
!!               size (nVarD)
!!
!!   nmateriaux  Number of materials
!!
!!   VCinfo      Info on composite model integration
!!
!!         Parallel fields accessible in AMITEX_MOD
!!--------------------------------------------------------
!! Description of fields in Fourier space:
!!        Z pencils, (nx/2+1, ny, nz, 3, nVarD or NTensDef)
!!        3D indices: (fft_start(1):fft_end(1), fft_start(2):fft_end(2), fft_start(3):fft_end(3), xx)
!!
!! Description of fields in real space:
!!        X pencils, (nx, ny, nz, 3, nVarD) or (nx, ny, nz, NTensDef)
!!        1D index: (1:xsize(1)*xsize(2)*xsize(3), NTensDef or 3, nVarD)
!!
!!   Sig      Cauchy stress field
!!   PK1      Piola-Kirchhoff stress field under large deformation
!!   Def      Strain field in HPP and grad(u) in GD
!!   Def0     Previous time step strain field in HPP and grad(u) in GD
!!   Sig0     Previous time step Cauchy stress field
!!   Def_star Additional strain field such that: Def = E + Def_per + Def_star
!!   SigF     Cauchy stress in HPP or Piola-Kirchhoff stress in GD,
!!            in Fourier space
!!   DefF     Strain in HPP or grad(u) in GD, in Fourier space
!!            Z pencils (nx/2+1, ny, nz, nTensDef)
!!   FluxD0   Vector field of fluxes (:,3,nVarD) (diffusion)
!!   FluxD    Vector field of fluxes (:,3,nVarD) (diffusion)
!!   GradQD0  Gradient field of diffusion variable (:,3,nVarD) (diffusion)
!!   GradQD   Gradient field of diffusion variable (:,3,nVarD) (diffusion)
!!   FluxDF   Flux vector field in Fourier space (diffusion)
!!   GradQDF  Gradient field in Fourier space
!!
!! Convergence acceleration arrays
!!        For mechanics (ntot, Ntens, 4)
!!        For diffusion (ntot, 3*nVarD, 4)
!!   ACT3_R   Residuals for convergence acceleration
!!   ACT3_U   Increments of strain fields for convergence acceleration
!!   ACT3_RD  Same as ACT3_R for diffusion
!!   ACT3_UD  Same as ACT3_U for diffusion
!!
!!                  Other variables in AMITEX_MOD
!!--------------------------------------------------------
!!   Flog       Logical units for opening the .log file
!!   fic_vtk    Root name for vtk or .(mz)std output files
!!
!!                        Variables in GREEN_MOD
!!--------------------------------------------------------
!!   FREQ     Frequency array, Fourier space
!!            Z-pencils (nx/2+1, ny, nz, 3)
!!            3D indices: (fft_start(1):fft_end(1), fft_start(2):fft_end(2), fft_start(3):fft_end(3), 3)
!!
!!                         Variables in LOADING_MOD
!!--------------------------------------------------------
!!   local_load%t1 and t0 and dt   mechanical loading, temperature, and external parameters imposed at each iteration
!!                                 %t1 is used to impose mechanical loading
!!                                 %t0 and %dt are used in behavior to impose loading in T and external parameters
!!
!!   local_loadD%t1 and dt         "diffusion" loading imposed at each iteration
!!
!!                         Variables in NL_BASE_MOD
!!--------------------------------------------------------
!!   crit_b                criteria used in resolution_NL_base
!!   times_b               elapsed times for resolution_NL_base
!!
!!------------------------------------------------------------
!!
!!   \param[in]    ind_tps   (time step index)
!!   \param[in]    load_n, load_incr
!!                 partial loading number and calculation step index
!!   \param[out]   nIt
!!                 number of iterations performed
!!   \param[in]    locload   (optional)   LOCALLOAD type variable. If present:
!!   \param[in]    locloadD  (optional)   overrides local_loading / local_loadingD
!!                                      after get_current_loading
!!
!!------------------------------------------------------------------------------
!!                         OTHER DERIVED-TYPE VARIABLES
!!
!!      MattotP         array of MATERIAL structures   (material_mod.f90)
!!            -> material description
!!
!!      load            array of LOADING structures    (loading_mod.f90)
!!            -> loading description
!!
!!      initValExt      INITLOADEXT structure          (loading_mod.f90)
!!            -> initialization description of external loading
!!
!!      extract         PARAMEXTRACT structure         (loading_mod.f90)
!!            -> VTK output description
!!
!!      algo_param      PARAM_ALGO structure           (param_algo_mod.f90)
!!            -> algorithmic parameters description
!!
!!      grid            GRID_DIM structure             (green_mod.f90)
!!            -> grid description
!!
!!
!!===============================================================================


subroutine unpas_NL_base(load_n,load_incr,ind_tps,nIt,locload,locloadD)

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================

  implicit none
!!------------------------------------------------------------------------------
!>                                                     "PARALLELISM" PARAMETERS

  integer                          :: ierror             !< error related to MPI funcitons
  
!!------------------------------------------------------------------------------
!>                                                   "GRID DIMENSION" PARAMETERS

  integer                          :: ni2              !< dimension divided by 2


!!------------------------------------------------------------------------------
!>                                               "REFERENCE MATERIAL" PARAMETERS

  real(mytype)                     :: lambda, mu

!!------------------------------------------------------------------------------
!>                                                         "LOADING" PARAMETERS"

  integer,intent(in)                                     :: ind_tps    !< index of the current time step 
  real(mytype),dimension(algo_param%nTensDef)            :: Eimp       !< imposed average strain
  real(mytype),dimension(3*algo_param%nVarD)             :: EimpD      !< imposed average gradQ
  integer,intent(in)                                     :: load_n, load_incr
                                                                       !< partial load number 
                                                                       !< increment number of the loading
                                                                       !<(load_incr=0 : extra increment at the start of loading)
  integer                                                :: nb_stress_pilot,nb_flux_pilot 
                                                                       !< number of directions pilotd in stress
                                     
!!------------------------------------------------------------------------------
!>                                                 OPTIONAL "LOADING" PARAMETERS

  type(LOCALLOAD), intent(in), optional                  :: locload, locloadD

!!------------------------------------------------------------------------------
!>                                         CONVERGENCE ACCELERATION (IF ENABLED)
!>                                                 parallel array (ntot,Ntens,4)
!>                                             parallel arrayss (ntot,3*nVarD,4)
!>      increments of strain fields and residuals from the 4 previous time steps

  integer                           :: iact3           !< counter for ACV (convergence acceleration)
  logical                           :: lact3           !< flag indicating if ACV is enabled
  logical                           :: acv_test, acv_testD        
                                                                  !< act3 output : 
                                                                  !< prevents storing a non-accelerated field

!!------------------------------------------------------------------------------
!>                                               FORCED CONVERGENCE (IF ENABLED)
!>                     
  logical                           :: testCVFOR        !< flag to activate forced convergence 
  integer                           :: countCVFOR       !< counter for forced convergence 

!!------------------------------------------------------------------------------
!>                                                 ALGORITHM TRACKING PARAMETERS

  integer,intent(out)               :: nIt                    !< number of iterations in one loading step

!!------------------------------------------------------------------------------
!>                                        LAMINATE ALGORITHM TRACKING PARAMETERS

  !< performance tacking of the laminate algorithm
  !<               nSub_laminate         : number of time step subdivisions 
  !<               nIncr_laminate        : number if sub-time steos per subdivision 
  !<               nIt_laminate          : number of iterations per laminate algorithm call
  !<               nCall_laminate        : number of calls to the laminate algorithm
  integer(kind=8)               :: nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp !< per behavior call

!!------------------------------------------------------------------------------
!>                                               RESOLUTION ALGORITHM PARAMETERS

  logical                                                  :: test_out            !< exit condition for while loop


!!------------------------------------------------------------------------------
!>                                                          ABAQUS COMPATIBILITY
  integer(kind=8)                    :: KSTEP 

!!------------------------------------------------------------------------------
!>                                       PARAMETERS FOR TRACKING EXECUTION TIMES

!! double precision type instead of real_mytype (=> required by MPI_WTIME)
  double precision                   :: t_ini,t1,t10

  !< t1          : timestep before entering function blocks
  !! t10         : cumulative time taken by the last 10 iterations  
  
!!------------------------------------------------------------------------------
!>                                                                 MISCELLANEOUS

  integer               :: i,j,k,l,m !< loop indices 
  integer               :: nVarD     !< number of diffusion variables
  logical               :: testNaN1,testNaN2,testdef

  integer               :: TPRED_deb !< index in local_load%t1 for T 
  integer               :: TPRED_fin !< index in local_load%t1 for the last param_ext

!!==============================================================================
!!                                                                    INIT TIMES
!!==============================================================================
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  t_ini = MPI_WTIME() 
  t10 = t_ini  !< timestamp before start of while loop for cumulative timing

!!==============================================================================
!!                                                           AUXILIARY VARIABLES
!!==============================================================================

  TPRED_deb = 2*algo_param%nTensDef+1 
  TPRED_fin = 2*algo_param%nTensDef+1+initValExt%nb_param    

!================================================================================ 
!                                                                 INITIALISATIONS
!

  lambda=0; mu=0
  nVarD = algo_param%nVarD
  if (algo_param%Mechanics) then
    lambda = Matref%LambdaMu0(1)    !< reference material
    mu = Matref%LambdaMu0(2)
    Eimp=0
  end if
  if (algo_param%Diffusion) then  
    EimpD=0
  end if
  nSub_laminate_temp=0
  nIncr_laminate_temp=0
  nIt_laminate_temp=0

       !================================================================================
       !                                                         FORCED CONVERGENCE LOOP
       !
       !  Forces convergence if testCVFOR becomes .false. after exiting the while loop
       ! 
       testCVFOR  = .true.
       countCVFOR = 0
       do while (testCVFOR)

        VCinfo%nSub_lam_pas          = 0    !< resetting laminate algorithm 
        VCinfo%nIncr_lam_pas         = 0    !  performance indicators
        VCinfo%nIt_lam_pas           = 0
        VCinfo%nCall_lam_pas         = 0
        !----------------------------------------------------------------------
        !                               MANAGEMENT OF ADDITIONAL HALF_TIME STEP

        !! load_incr=0 corresponds to an intermediate load between
        !! the final load imposed at load_n-1 and 
        !! the first iteration of load step load_n 
        !! the ratio of this intermediate is defined in loading_mod
        KSTEP = ind_tps       ! Note: due to the fictitous time step, 
                              ! KSTEP is repeated twice at each start of a new load
                     
        if(nrank==0)then
          if(load_incr == 0)then
              write(OUTPUT_UNIT,"(/,A)")"=================================================="
              write(OUTPUT_UNIT,"(A,I8,A)")"Pseudo time step :",ind_tps,".5"
              write(Flog,"(/,A)")"=================================================="
              write(Flog,"(A,I8,A)")"Psuedo time step :",ind_tps,".5"
          else
              write(OUTPUT_UNIT,"(/,A)")"=================================================="
              write(OUTPUT_UNIT,"(A,I8)")"Time step :",ind_tps
              write(Flog,"(/,A)")"=================================================="
              write(Flog,"(A,I8)")"Time steps :",ind_tps
          end if 
        end if

        !----------------------------------------------------------------------      
        !                                            IF A LOCAL_LOAD IS PRESENT      
        !

        if (algo_param%Mechanics .AND. present(locload)) local_load = locload
        if (algo_param%Diffusion .AND. present(locloadD)) local_loadD = locloadD

        !----------------------------------------------------------------------      
        !                                                         LOAD ANALYSIS       
        !

        if (algo_param%Mechanics) nb_stress_pilot = count(local_load%t1(1:algo_param%nTensDef) == 1)
        if (algo_param%Diffusion) nb_flux_pilot = count(local_loadD%t1(1:3*nVarD) == 1)  

        !----------------------------------------------------------------------
        !                                           PROPOSAL FOR A STRAIN FIELD
        !
        ! with forced convergence:
        !    First initialisation same as without forced convergence
        !    THEN
        !    if algo_param%init_cvfor=='default' : 
        !       initialisation same as without forced convergence.
        !    if algo_param%init_cvfor=='last' or 'best': 
        !       recover the DEF field chosen at the end of the previous forced convergence
        !    'default' seems to give the best results
        !
        ! wihtout forced convergence: the test below is always true
        !
        if (algo_param%init_cvfor=='default' .OR. countCVFOR == 0) then
        if (algo_param%Mechanics) then
            if (algo_param%init_def /= 'None') then !-- may be set to 'None' in a user function to propose a different initialisation 
            if (algo_param%HPP_nsym) then
                ! consider full displacement gradient in  
                ! HPP calculation: initialise Def_nsym
                call initDef(local_load%t1(1:2*algo_param%nTensDef),&
                         nb_stress_pilot,local_load%t_load,load_incr,load_n,Matref%C0,Def0,Def,&
                         Def_nsym0,Def_nsym); !load_incr,load_n  for testing special cases 0 and 1
            else 
                call initDef(local_load%t1(1:2*algo_param%nTensDef),&
                         nb_stress_pilot,local_load%t_load,load_incr,load_n,Matref%C0,Def0,Def);

            end if
            else
                call write_stdout0("STRAINS (DEF,DEF0) NOT INITIALIZED IN UNPAS_NL_BASE")
            end if
        end if

        if (algo_param%Diffusion) then
            call initDefD(local_loadD%t1(1:6*algo_param%nVarD),&
                         nb_flux_pilot,local_loadD%t_load,load_incr,load_n,Matref%K0D,GradQD0,GradQD);
        end if
        end if

        !----------------------------------------------------------------------
        !                               INITIALISATION OF CRITERIA AND COUNTERS

        if (algo_param%Diffusion) then
            crit_b%eqD=0_mytype
            crit_b%CptbD=0_mytype
            crit_b%FluxDMoy=0_mytype
            crit_b%GradQDMoy=0_mytype  
        end if
        if (algo_param%Mechanics) then
            crit_b%eq=0_mytype
            crit_b%Cptb=0_mytype
            crit_b%DefMoy=0_mytype
            crit_b%SigMoy=0_mytype
        end if

        nIt=-1  ! number of iterations for the current load step
                ! becomes 0 when entering the while loop
                ! then the algorithm forces at least one iteration
        iact3=0
        lact3=.false.
        acv_test = .true.
        acv_testD = .true.

        !----------------------------------------------------------------------
        !                                                          DefF=FFT(Def)
        !                  format DefF: (pinceaux-Z, global size (nx/2+1)*ny*nz)
        !                                             is normalised DefF by ntot

        if (algo_param%Mechanics) then
           call field_fft(Def,DefF,algo_param%nTensDef,1)
           DefF=DefF/ real(grid%ntot,mytype)
        end if

        if (algo_param%Diffusion) then
           call field_fft(gradQD,GradQDF,3,algo_param%nvarD)
           GradQDF=GradQDF / real(grid%ntot,mytype)
        end if

       ! To test eva;uation of deformation criterion on a simple field
       !  call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
       !        load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
       !        crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)


        !----------------------------------------------------------------------
        !                                   FOR THE 1st LOAD STEP (load_incr=0)
        !                PROPOSAL FOR A DEFORMATION FIELD (SUITE): ADD gragradU
        !                                                          ADD Def_star
        !
        !                           gradgradU must be aff AFTER evaluating DefF
        !
        if (algo_param%init_cvfor=='default' .OR. countCVFOR == 0) then
        if (algo_param%Mechanics) then
                ! Case with imposed gradgradU
            ! Def = def + (gradgradU) (_sym si HPP)
            if (n_gradgradU==27 .AND. load_incr==0) then ! not very elegant test <=> test if gradgradU imposed
               call add_gradgradU(Def,algo_param%nTensDef,&
                          local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
            end if
            if (add_def_star .AND. load_incr==0) then 
               Def = Def + Def_star
            end if            
         end if
         end if

        !----------------------------------------------------------------------
        !                                                             BEHAVIOUR

        if (algo_param%Mechanics) then !----- MECHANICS
             if (algo_param%HPP_nsym)  then
                ! take into account full displacement gradient in 
                ! HPP calculation
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,VCinfo%nCall_lam_pas,&
                           Def_nsym=Def_nsym,Def_nsym0=Def_nsym0)
             else
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,VCinfo%nCall_lam_pas)
 
             end if
           if(.not. algo_param%HPP) then
              call SigToPK1(Sig,Def,PK1)
           end if
           if (algo_param%Nloc .and. .false. ) then !---- None-local mechanics, implicit resoltuion (TO BE REVIEWED)
              !TODO: review the possibility of introducing implicit treatment
              call Nloc_call(load_n,load_incr,ind_tps,"after", algo_param%nloc_implicit)
           end if
        end if

        if (algo_param%Diffusion) then !----- DIFFUSION
           call behaviorD(FluxD,FluxD0,GradQD,GradQD0,local_loadD%t0(6*algo_param%nVarD+1:),&
                          local_loadD%dt(6*algo_param%nVarD+1:),local_loadD%t_load,KSTEP)
        end if

        ! Increment counters after behaviour call
        !       at thjis point, nCall_laminate contains number of composite voxels "laminate"
        !       for this load step (in case this number changes between steps) 
        VCinfo%nSub_lam_pas   = VCinfo%nSub_lam_pas  + nSub_laminate_temp
        VCinfo%nIncr_lam_pas  = VCinfo%nIncr_lam_pas + nIncr_laminate_temp
        VCinfo%nIt_lam_pas    = VCinfo%nIt_lam_pas   + nIt_laminate_temp

        !================================================================================
        !                                                                    BOUCLE WHILE
        do
          nIt=nIt+1   ! actualisation du compteur de boucle, 
                      ! au 1er passage avant tests des criteres on a i=0

!interet de cette initialisation des criteres moyens?
!         if (algo_param%Mechanics) then
!              crit_b%DefMoy=0
!              crit_b%SigMoy=0  
!         end if
!         if (algo_param%Diffusion) then
!              crit_b%FluxDMoy=0
!              crit_b%GradQDMoy=0      
!         end if
           !----------------------------------------------------------------------
           !                                FFT DE LA CONTRAINTE : SigF = FFT(Sig)
           !              format SigF: (pinceaux-Z, taille globale (nx/2+1)*ny*nz)
          

           if(algo_param%Mechanics) then  !-------MECANIQUE
           if(algo_param%HPP) then !>--HPP
              call field_fft(Sig,SigF,algo_param%nTensDef,1)
           else                    !>---GRANDES DEF
              call field_fft(PK1,SigF,algo_param%nTensDef,1)
              call mpi_barrier(mpi_comm_world,ierror)
              ! correction dimensions paires
              if (algo_param%CorrPair) then
              if(modulo(grid%ny,2) == 0) then
                 ni2 = grid%ny/2 +1
                 if(ni2>= fft_start(2) .AND. ni2<= fft_end(2) ) then 
                    SigF(fft_start(1):fft_end(1),ni2,fft_start(3):fft_end(3),1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
              if(modulo(grid%nz,2) == 0) then
                 ni2 = grid%nz/2 +1
                 if(ni2>= fft_start(3) .AND. ni2<= fft_end(3) ) then 
                    SigF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),ni2,1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
              end if
           end if 
           end if

           if(algo_param%Diffusion) then  !-------DIFFUSION
              call field_fft(FluxD,FluxDF,3,algo_param%nvarD)
           end if

!print *, "AV ACV",nrank,shape(ACT3_U),shape(Def)
           !----------------------------------------------------------------------
           !                                 STOCKAGE ACCCELERATION DE CONVERGENCE
           ! TODO : dissocier les deux acv meca et diffusion

           ! on ne stocke pas s'il n'y a pas eu ACV (acv_test = .false.)
           if(algo_param%acc_CV .and. acv_test .and. acv_testD) then
              iact3 = iact3+1
              if (iact3==5) then 
                 iact3 = 1
              end if
              if (algo_param%Mechanics) then
                if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : stockage de gradient de u pour l'acceleration
                    !ACT3_U(:,:,iact3) = Def_nsym(:,:)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i = 1,algo_param%nTensDef
                    do j = 1, size(ACT3_U,dim=1)
                       ACT3_U(j,i,iact3) = Def_nsym(j,i)
                    end do
                    end do

                else
                    !ACT3_U(:,:,iact3) = Def(:,:)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i = 1, algo_param%nTensDef
                    do j = 1, size(ACT3_U,dim=1)
                       ACT3_U(j,i,iact3) = Def(j,i)
                    end do
                    end do
                end if
              end if
              if (algo_param%Diffusion) then
                  do j = 1,nVarD
                  do i = 1,3
                  do l = 1, size(ACT3_UD,dim=1)
                     k = i+3*(j-1)
                     ACT3_UD(l,k,iact3) = GradQD(l,i,j)           
                  end do         
                  end do
                  end do
              end if
           end if
!print *, "AP ACV",nrank

           !----------------------------------------------------------------------
           !                                 EVALUATION DES CRITERES EN CONTRAINTE
           !                                                          + VERIF. NaN
           
           if(algo_param%Mechanics) then !----------MECANIQUE
           if(algo_param%HPP) then
              call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                   load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                   crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,0_INT64)
              !A DECOMMENTER pour critere en def. a chaque pas
              ! call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
              !     load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
              !     crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)
!if (nrank==0) print *, "CRITEQ                                   =   ", crit_b%eq
           else
              call eval_criteres_GD(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                   load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                   crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,0_INT64)
              !
              !block pour tester l'utilisation du criteres HPP
              !------------------------------------------------
              !block
              !complex(mytype),allocatable,dimension(:,:,:,:)  ::SigF2
              !allocate(SigF2(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:6))
              !SIGF2=(SigF(:,:,:,1:6)+ SigF(:,:,:,(/1,2,3,7,8,9/)))/2.
              !call eval_criteres(SigF2,DefF(:,:,:,1:6),&
              !     FREQ,(/ local_load%t1(1:6), local_load%t1(10:15) /),&
              !     load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
              !     crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,0_INT64)
              !end block
           end if
           end if

           if(algo_param%Diffusion) then !----------DIFFUSION
              call eval_criteresD(FluxDF,GradQDF,FREQ,local_loadD%t1(1:6*nVarD),&
                   crit_b%eqD,crit_b%FluxDMoy,crit_b%GradQDMoy,crit_b%CptbD,0_INT64)        
           end if 

                                         !----------TESTS CRITERES NAN
           testNaN1=.false.
           testNaN2=.false.
           if (algo_param%Mechanics) then
                  testNaN1=isNAN(crit_b%eq)
                  testNaN2=isNAN(crit_b%SigMoy)
           end if
           if (algo_param%Diffusion) then
                do j = 1,nVarD
                        testNaN1 = testNaN1 .or. isNaN(crit_b%eqD(j))
                        testNaN2 = testNaN2 .or. isNaN(crit_b%FluxDmoy(j))
                end do
           end if
           if(testNaN1) then
              if(nrank==0) then
                 write (Flog,"(A,I8)")"-------- ERROR : EQUILIBRIUM CRITERION NaN----- iteration ",nIT
                 if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "crit_eq = ", crit_b%eq
                 if (algo_param%Diffusion) write (Flog,"(A,(E15.8))") "crit_eqD = ", crit_b%eqD
              end if
              call amitex_abort("Equilibrium criterion is  NaN",2,0)
           end if
           if(testNaN2) then
              if(nrank==0) then
                 write (Flog,"(A,I8)")&
                 "-------- ERROR : AVARAGE STRESS (or FLUX) CRITERION NaN----- iteration ",nIT
                 if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "critSigMoy = ", crit_b%SigMoy
                 if (algo_param%Diffusion) write (Flog,"(A,(E15.8))") "critFluxDmoy) = ", crit_b%FluxDmoy
              end if
              call amitex_abort("Average stress (or flux) criterion is NaN",2,0)
           end if

!print *, "AV STOCKAGE CV FOR",nrank
           !----------------------------------------------------------------------
           !                                           STOCKAGE CONVERGENCE FORCEE
           !
           ! nIT >1 est INDISPENSABLE dans le cas algo_param%init_cvfor='best'
           !        pour eviter de tomber sur des cycles d'iterations de CV forcee identiques
           if (algo_param%CVfor .AND. (nIt > 1) ) then
           if (crit_b%eq < CVFORsauv%critmin) then
              CVFORsauv%critmin=crit_b%eq
              CVFORsauv%Def=Def
              if (nrank==0) write (Flog,"(A,(E15.8))") "crit_eq NON-CONV= ", crit_b%eq
              if (nrank==0) write (OUTPUT_UNIT,"(A,(E15.8))") "crit_eq NON-CONV= ", crit_b%eq
           end if
           end if

           !----------------------------------------------------------------------
           !                                  SORTIE CONDITIONNELLE BOUCLE INFINIE
           !
           !    criteres Contrainte (Flux) Macro et Equilibre < tol
           !    Si demande (Nitermin=1 et/ou Nitermin_acv=1) :
           !      minimum 1 itération (= application de Green)
           !      minimum 1 itération après accélération de convergence (= appli. Green)
           !
           !     => application de Green assure la compatibilité 
           !    (on verifie néanmoins le critere de compatibilité avant de sortie)
           !
           !                                  AVANT DE SORTIR: 

           !! EVALUATION DU TEST DE SORTIE
           test_out=.true.

           if (algo_param%NiterMin_acv==1) then 
              test_out = (.NOT. lact3)
           end if

           if (algo_param%NiterMin==1) then
              test_out = test_out .AND. (nIt .ge. 1) 
           end if

           if (algo_param%Mechanics) then 
                test_out = test_out .AND. (crit_b%eq .LT. algo_param%tol_criteq) .AND. &
                           (crit_b%SigMoy .LT. algo_param%tol_critsig)
           end if

           if (algo_param%Diffusion) then
                  test_out = test_out .AND. &
                        (max(maxval(crit_b%eqD),maxval(crit_b%FluxDMoy)) .LT. algo_param%tol_criteq)
           end if

           !! SORTIE
           if (test_out) then
              !>-----------EVALUATION DES CRITERES EN DEFORMATION (OU GRADQ)

              if(algo_param%Mechanics) then  !-------MECANIQUE
              if(algo_param%HPP) then
                 call eval_criteres(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                      load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                      crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)
              else
                 call eval_criteres_GD(SigF,DefF,FREQ,local_load%t1(1:2*algo_param%nTensDef),&
                      load(load_n)%DirStress_flag, load(load_n)%DirStress2,& 
                      crit_b%eq,crit_b%SigMoy,crit_b%DefMoy,crit_b%Cptb,grid%ntot)
              end if
              end if
              if(algo_param%Diffusion) then !-------DIFFUSION
                 call eval_criteresD(FluxDF,GradQDF,FREQ,local_loadD%t1(1:6*nVarD),&
                   crit_b%eqD,crit_b%FluxDMoy,crit_b%GradQDMoy,crit_b%CptbD,grid%ntot)
              end if

!print *, "AP EVAL CRITERES DEF",nrank
              !>-----------TESTS NAN 
              testNaN1=.false.
              testNaN2=.false.
              if (algo_param%Mechanics) then
                     testNaN1=isNAN(crit_b%Cptb)
                     testNaN2=isNAN(crit_b%DefMoy)
              end if
              if (algo_param%Diffusion) then
                     do j = 1,nVarD
                        testNaN1 = testNaN1 .or. isNaN(crit_b%CptbD(j))
                        testNaN2 = testNaN2 .or. isNaN(crit_b%GradQDmoy(j))
                     end do
              end if
              if(testNaN1) then
                 if(nrank==0) then
                    write (Flog,"(A)")"-------- ERROR: COMPATIBLILTY CRTIERION IS  NaN-----"
                    if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "critCptb = ", crit_b%eq
                    if (algo_param%Diffusion)write (Flog,"(A,(E15.8))") "critCptbD = ", crit_b%eqD
                 end if
                 call amitex_abort("Compatibilty criterion is NaN",2,0)
              end if
              if(testNaN2) then
                 if(nrank==0) then
                    write (Flog,"(A)")&
                         "-------- ERROR: AVERAGE STRAIN (or GRADQD) IS NaN-----"
                    if (algo_param%Mechanics) write (Flog,"(A,E15.8)") "critDefMoy = ", crit_b%DefMoy
                    if (algo_param%Diffusion) write (Flog,"(A,(E15.8))") "critFluxDmoy) = ", crit_b%GradQDmoy
                 end if
                 call amitex_abort("Average strain (or GRADQD) criterion is NaN",2,0)
              end if

!print *, "AV SORTIE SUR CRITERES EN DEF",nrank
              !>-----------SORTIE SUR CRITERES EN DEF (ET GRADQD) (compatibilite & moy) 
              testdef = .true.
              if (algo_param%Mechanics) then
                  testdef = max(crit_b%DefMoy,crit_b%Cptb) .LT. algo_param%tol_critdef
              end if
              if (algo_param%Diffusion) then
                 testdef = testdef & 
                     .and. (max(maxval(crit_b%GradQDMoy),maxval(crit_b%CptbD)) .LT. algo_param%tol_critdef)
              end if
!print *, "testdef",testdef
              if(testdef) then
                 !! si ce critere est respecte, on sort de la boucle
                 !! Avant de sortir, en grandes transformations et si CorrPair=.true.
                 !! et si la resolution en y ou en z est paire, on fait une FFT inverse 
                 !! pour avoir le tenseur de
                 !! Piola-Kirchhoff sans oscillation et on recupere le tenseur des
                 !! contraintes de Cauchy 

                   if (algo_param%Mechanics .AND. algo_param%CorrPair) then
                   if((.NOT. algo_param%HPP) .AND. (modulo(grid%ny,2)==0 .OR. modulo(grid%nz,2)==0)) then
                       call field_ifft(PK1,SigF,algo_param%nTensDef,1)
                       PK1=PK1/ real(grid%ntot,mytype)
                       call PK1ToSig(Sig,Def,PK1)
                   end if 
                   end if
                 testCVFOR=.false.  ! Sortie sur criteres respectes : test convergence forcee => .FALSE.
                 !! on reinitialise le critere min  de la CV forcee pour une iteration de CVfor  
                 CVFORsauv%critmin = 1.e100_mytype
                 exit
              end if
           end if
           
!print *, "BEFORE EVALUATING CONDITIONAL EXIT ERROR",nrank
           !----------------------------------------------------------------------
           !                                             CONDITIONAL EXIT ON ERROR 
           !                maximum number of iterations reached without forced CV

           if(nIt .gt. algo_param%nIterMax .AND. (.not. algo_param%CVfor)) then
                   !! If too many iterations are done, stop the computation
              call amitex_abort("Maximum number of iterations reached, computation stopped.",1,0)
              t1 = MPI_Wtime()

              !! Before exiting the progra :
              !! output the VTK images of the fields from the last completed iteration
              if (.not. algo_param%CVfor  .or. (algo_param%CVfor .and. (nit .eq. algo_param%Nit_cvfor))) then
                 if(nrank == 0) write(Flog,"(A)")"Writing VTK of the last converged time stepe."
                 call writeVTK(ind_tps)
              end if

              if (nrank == 0) then
                 write(OUTPUT_UNIT,"(A)") "Maximum number of iterations reached, computation stopped."
                 write(OUTPUT_UNIT,"(/,A,/)")"-----------------------------"
                 write (Flog,"(A)") "Maximum number of iterations reached, computation stopped."
                 write (Flog,"(A,/)")"-----------------------------"
              end if
              exit
           end if
           
!print *, "BEFORE EVALUATING CONDITIONAL EXIT FORCED CONVERGENCE",nrank
           !----------------------------------------------------------------------
           !                                CONDITIONAL EXIT ON FORCED CONVERGENCE 
           !                                            1ST ITERATION OF FORCED CV
           !                                            criterion n_terMax reached
           !                         NOTE: no forced convergence on the first step
           ! 
           if(algo_param%CVfor .AND. (load_incr .ne. 0) .AND. &
              (countCVfor == 0) .AND. (nIt .gt. algo_param%nIterMax)) then              
              testCVFOR=.true.
              exit
           end if

           !----------------------------------------------------------------------
           !                                CONDITIONAL EXIT ON FORCED CONVERGENCE
           !                                         OTHER ITERATIONS OF FORCED CV
           !                                           criterion nit_cvfor reached
           !                         NOTE: no forced convergence on the first step
           !  
           if (algo_param%CVfor .and. (nit .eq. algo_param%Nit_cvfor) .and. (load_incr .ne. 0)) then
              if (countCVFOR .eq. (algo_param%Ncvfor-1) ) then
                if(nrank == 0) write(Flog,"(A)")"Writing vtk of the last forced convergence step."
                call writeVTK(ind_tps)
              end if
              testCVFOR=.true.
              exit
           end if
           
!print *, "AV MAJ Eimp",nrank
           !----------------------------------------------------------------------
           !                            MISE A JOUR DE LA DEFORMATION IMPOSEE Eimp

           if(algo_param%Mechanics) then !--MECANIQUE
           if(nb_stress_pilot>0) then 
              call updateEimp(Eimp, local_load%t1(1:2*algo_param%nTensDef),nb_stress_pilot,&
                load(load_n)%DirStress_flag, load(load_n)%DirStress_cauchy, &
                load(load_n)%DirStress,load(load_n)%DirStress2, &
                Matref%C0, SigF, DefF, grid%ntot,fft_start)
           else
              Eimp = local_load%t1(algo_param%nTensDef+1:2*algo_param%nTensDef);
           end if
           end if

           if(algo_param%Diffusion) then !--DIFFUSION
           if(nb_flux_pilot>0) then
              call updateEimpD(EimpD, local_loadD%t1(1:6*algo_param%nVarD),&
                   Matref%K0D, FluxDF, GradQDF, grid%ntot,fft_start)
           else
              EimpD = local_loadD%t1(3*nVarD+1:6*nVarD);
           end if
           end if
           
!print *, "AV GREEN",nrank
           !----------------------------------------------------------------------
           !                                                     APPLICATION GREEN  
           ! format dans Fourier: pinceaux-Z, taille globale ((nx/2)*ny*nz)
           ! ATTENTION : on traite ici les cas des grilles paires
           ! ATTENTION : en simple precision la FFT doit etre realisee sur un champs normalise
           !            (sinon, risque de depassement du reel max) 
           !            sortie de apply_green non normalisee
           ! REMARQUE : SIGF est utilise ici en entree/sortie de appli_Green
           ! 

        if(algo_param%Mechanics) then !----- MECANIQUE
           !! DefF = -\Gamma_0 * (SigF - C_0:DefF)          !modif intel/oneapi (:,:,:,1) -> (i,j,k,1)
!print *, "AV DefF = -\Gamma_0 * (SigF - C_0:DefF)"         !   sinon segfault aleatoires
           do k=fft_start(3),fft_end(3)
           do j=fft_start(2),fft_end(2)
           do i=fft_start(1),fft_end(1)     
           SigF(i,j,k,1) = SigF(i,j,k,1) / real(grid%ntot,mytype) &
                -  lambda*(DefF(i,j,k,1)+DefF(i,j,k,2)+DefF(i,j,k,3))- 2._mytype*mu*DefF(i,j,k,1)
           SigF(i,j,k,2) = SigF(i,j,k,2) / real(grid%ntot,mytype) &
                -  lambda*(DefF(i,j,k,1)+DefF(i,j,k,2)+DefF(i,j,k,3))- 2._mytype*mu*DefF(i,j,k,2)
           SigF(i,j,k,3) = SigF(i,j,k,3) / real(grid%ntot,mytype) &
                -  lambda*(DefF(i,j,k,1)+DefF(i,j,k,2)+DefF(i,j,k,3))- 2._mytype*mu*DefF(i,j,k,3)
           end do
           end do
           end do

           if(algo_param%HPP) then
              do k=fft_start(3),fft_end(3)
              do j=fft_start(2),fft_end(2)
              do i=fft_start(1),fft_end(1)     
              SigF(i,j,k,4:algo_param%nTensDef) = SigF(i,j,k,4:algo_param%nTensDef) &
                   / real(grid%ntot,mytype) - mu*DefF(i,j,k,4:algo_param%nTensDef)
              end do
              end do
              end do
              !! Ici SigF = FFT(sigma - C0:epsilon)
              if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : SigF = sym(DefF_nsym)
                call apply_Green(SigF,FREQ,Matref%LambdaMu0,DefF_nsym)
              else
                call apply_Green(SigF,FREQ,Matref%LambdaMu0)
              end if
           else
              if(algo_param%C0sym) then ! -> version symetrique : C0:gradU_sym
                 do k=fft_start(3),fft_end(3)
                 do j=fft_start(2),fft_end(2)
                 do i=fft_start(1),fft_end(1)     
                 SigF(i,j,k,4:algo_param%nTensDef) = SigF(i,j,k,4:algo_param%nTensDef) &
                    / real(grid%ntot,mytype) &
                    - mu*( DefF(i,j,k,4:algo_param%nTensDef) + DefF(i,j,k,(/7,8,9,4,5,6/)) ) 
                 end do
                 end do
                 end do 
              else                      ! -> version non symetrique : C0:gradU
                 do k=fft_start(3),fft_end(3)
                 do j=fft_start(2),fft_end(2)
                 do i=fft_start(1),fft_end(1)     
                 SigF(i,j,k,4:algo_param%nTensDef) = SigF(i,j,k,4:algo_param%nTensDef) &
                    / real(grid%ntot,mytype) - 2._mytype*mu*DefF(i,j,k,4:algo_param%nTensDef)
                 end do
                 end do
                 end do 
              end if
              !! Ici SigF = FFT(Pi - C0:grad(u))
              call apply_Green_GD(SigF,FREQ,Matref%LambdaMu0)
           end if
           !! En sortie SigF = FFT(-Gamma_0 * Tau)
           !DefF = SigF
           ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
!print *, "AV DefF=SigF"
           do l=1,algo_param%nTensDef
           do k=fft_start(3),fft_end(3)
           do j=fft_start(2),fft_end(2)
           do i=fft_start(1),fft_end(1)     
               DefF(i,j,k,l) = SigF(i,j,k,l)
           end do
           end do
           end do
           end do
!print *, "AP DefF=SigF"
           
           !! On impose la valeur moyenne Eimp a Def en modifiant DefF(0)
           !! si la frequence 0 est dans le processus
           if ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
              DefF(1,1,1,1:algo_param%nTensDef) = cmplx(Eimp,kind=mytype)
           end if
           !! Si HPP et prise en compte complete du gradient du deplacement
           !! On impose la valeur moyenne du gradient du deplacement a Eimp 
           !!  en modifiant DefF_nsym(0) si la frequence 0 est dans le processus
           if (algo_param%HPP_nsym) then
                if ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
                    DefF_nsym(1,1,1,1) = cmplx(Eimp(1),kind=mytype)
                    DefF_nsym(1,1,1,2) = cmplx(Eimp(2),kind=mytype)
                    DefF_nsym(1,1,1,3) = cmplx(Eimp(3),kind=mytype)
                    DefF_nsym(1,1,1,4) = cmplx(Eimp(4),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,5) = cmplx(Eimp(5),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,6) = cmplx(Eimp(6),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,7) = cmplx(Eimp(4),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,8) = cmplx(Eimp(5),kind=mytype)/2._mytype
                    DefF_nsym(1,1,1,9) = cmplx(Eimp(6),kind=mytype)/2._mytype
                end if
           end if

           !! On annule les fréquences extremes si ny et nz sont pairs
           if(.NOT. algo_param%HPP .AND. algo_param%CorrPair) then 
              if(modulo(grid%ny,2) == 0) then
                 ni2 = grid%ny/2 +1
                 if(ni2>= fft_start(2) .AND. ni2<= fft_end(2) ) then 
                    DefF(fft_start(1):fft_end(1),ni2,fft_start(3):fft_end(3),1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
              if(modulo(grid%nz,2) == 0) then
                 ni2 = grid%nz/2 +1
                 if(ni2>= fft_start(3) .AND. ni2<= fft_end(3) ) then 
                    DefF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),ni2,1:9)=&
                         cmplx(0,kind=mytype)
                 end if
              end if
           end if

        end if !-----FIN MECANIQUE-------------------------------------

        if(algo_param%Diffusion) then !----- DIFFUSION-----------------

           !! DefF = -\Gamma_0 * (SigF - C_0:DefF)
           do i=1,nVarD
           do j=1,3
           do m=fft_start(3),fft_end(3)   ! modif intel15/oneapi2022 (:,:,:,j,i) -> (k,l,m,j,i)
           do l=fft_start(2),fft_end(2)
           do k=fft_start(1),fft_end(1)     

                FluxDF(k,l,m,j,i) =  FluxDF(k,l,m,j,i) / real(grid%ntot,mytype) &
                                 +  Matref%K0D(i)*GradQDF(k,l,m,j,i)
           end do
           end do
           end do
           end do
           end do

           !! Ici SigF = FFT(sigma - C0:epsilon)
           call apply_GreenD(FluxDF, FREQ,Matref%K0D)

           !! En sortie SigF = FFT(-Gamma_0 * Tau)
           GradQDF = FluxDF

           !! On impose la valeur moyenne EimpD a GradQ en modifiant GradQDF(0)
           !! si la frequence 0 est dans le processus
           if ((fft_start(1) .eq. 1) .AND. (fft_start(2) .eq. 1) .AND. (fft_start(3) .eq. 1)) then
              do i = 1,nVarD
                   GradQDF(1,1,1,1:3,i) = cmplx(EimpD(1+3*(i-1):3*i),kind=mytype)
              end do
           end if

        end if !-----FIN DIFFUSION

           !----------------------------------------------------------------------
           !                                                      Def = iFFT(DefF)
           !                                         here Def stores (E + Def_per)
           ! retour espace reel => format (pinceaux-X, taille globale, nx,ny,nz)      

           if(algo_param%Mechanics) then !----- MECANIQUE
               call field_ifft(Def,DefF,algo_param%nTensDef,1)
               if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP 
                    call field_ifft(Def_nsym,DefF_nsym,9,1)
               end if
           end if 

           if(algo_param%Diffusion) then !----- DIFFUSION
               call field_ifft(GradQD,GradQDF,3,algo_param%nvarD)
           end if 

           !----------------------------------------------------------------------
           !                           ON AJOUTE LE GRADIENT DE DEFORMATION IMPOSE  
           !                                                           ET DEF_STAR    
           if(algo_param%Mechanics) then !----- MECANIQUE
               !Def = Def + Defimp_hetero
               if (n_gradgradU==27) then
                  call add_gradgradU(Def,algo_param%nTensDef,&
                          local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
               end if
               if (add_def_star) Def = Def + Def_star
           end if 

           !----------------------------------------------------------------------
           !                          APPLICATION DE L'ACCELERATION DE CONVERGENCE
           ! lact3 : flag pour signaler que l'ACV vient d'être utilisée 
           ! TODO : dissocier les deux acv

           ! on stocke le residu
           if (algo_param%acc_CV .and. acv_test .and. acv_testD) then 
              if (algo_param%Mechanics) then
                if (algo_param%HPP_nsym) then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP : residu calcule sur le gradient de u
                    !ACT3_R(:,:,iact3) = Def_nsym(:,:) - ACT3_U(:,:,iact3)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i=1,algo_param%nTensDef
                    do j=1,size(ACT3_R,dim=1)
                        ACT3_R(j,i,iact3) = Def_nsym(j,i) - ACT3_U(j,i,iact3)
                    end do
                    end do
                else 
                    !ACT3_R(:,:,iact3) = Def(:,:) - ACT3_U(:,:,iact3)
                    ! modif intel15/oneapi2022 pour cas test fissure_1proc sur marquises machine batch 28p (ok sur machine interactive 14p)
                    do i=1,algo_param%nTensDef
                    do j=1,size(ACT3_R,dim=1)
                        ACT3_R(j,i,iact3) = Def(j,i) - ACT3_U(j,i,iact3)
                    end do
                    end do
                end if
              end if
              if (algo_param%Diffusion) then
                  do j = 1,nVarD
                  do i = 1,3
                  do l=1,size(ACT3_RD,dim=1)
                     k = i+3*(j-1)
                     ACT3_RD(l,k,iact3) = GradQD(l,i,j) - ACT3_UD(l,k,iact3)                     
                  end do
                  end do
                  end do
              end if
           end if

           ! on rebascule le test ACV a .true.
           if (.not. acv_test) acv_test = .true.
           if (.not. acv_testD) acv_testD = .true.

           ! on applique l'ACV
           lact3 = .false.
           if (algo_param%acc_CV) then 
              if(nIt>=4 .AND. modulo(nIt-4,algo_param%modACV)==0) then 
                 lact3 = .true.
                 if (algo_param%Mechanics) then
                    if (algo_param%HPP_nsym) then     
                         ! prise en compte du gradient complet du deplacement dans un 
                         ! calcul HPP : ACV
                       call act3(ACT3_U,ACT3_R,Def_nsym,GradQD,acv_test)
                       call field_fft(Def_nsym,DefF_nsym,9,1)
                       DefF_nsym=DefF_nsym/real(grid%ntot,mytype)

                       ! Calcul de Def = sym(Def_nsym) correspondant 
                       do j=1,size(Def,dim=1)              ! modif intel/oneapi (:,1) -> (j,1)
                       Def(j,1) = Def_nsym(j,1)
                       Def(j,2) = Def_nsym(j,2)
                       Def(j,3) = Def_nsym(j,3)
                       Def(j,4) = Def_nsym(j,4) + Def_nsym(j,7)
                       Def(j,5) = Def_nsym(j,5) + Def_nsym(j,8)
                       Def(j,6) = Def_nsym(j,6) + Def_nsym(j,9)   
                       end do
                       call field_fft(Def,DefF,algo_param%nTensDef,1)
                       DefF=DefF/ real(grid%ntot,mytype)                                      
                    else
                       call act3(ACT3_U,ACT3_R,Def,GradQD,acv_test)
                       if (n_gradgradU==27) then !si gradgradU, on retire la partie gradient pour le calcul de defF
                          call add_gradgradU(Def,algo_param%nTensDef,&
                           -local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
                       end if
                       if (add_def_star) Def = Def - Def_star        
                                      
                       call field_fft(Def,DefF,algo_param%nTensDef,1)
                       
                       DefF=DefF/ real(grid%ntot,mytype) 
                       if (n_gradgradU==27) then
                          call add_gradgradU(Def,algo_param%nTensDef,&
                           local_load%t1(TPRED_fin+1:TPRED_fin+27),load(load_n)%gevolution,&
                           grid%dx,grid%dy,grid%dz,grid%nx,grid%ny,grid%nz)
                       end if       
                       if (add_def_star) Def = Def + Def_star                                                         
                    end if
                 end if
                 if (algo_param%Diffusion) then
                   call act3(ACT3_UD,ACT3_RD,Def,GradQD,acv_testD)
                   call field_fft(GradQD,GradQDF,3,algo_param%nvarD)
                   GradQDF=GradQDF/ real(grid%ntot,mytype)
                 end if
              end if
           end if

           !----------------------------------------------------------------------
           !                                            EVALUATION DU COMPORTEMENT

           if (algo_param%Mechanics) then !-----MECANIQUE
             if (algo_param%HPP_nsym)  then
                ! prise en compte du gradient complet du deplacement dans un 
                ! calcul HPP 
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,&
                           Def_nsym=Def_nsym,Def_nsym0=Def_nsym0)
             else
                call behavior(SIG,SIG0,DEF,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                           local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                           nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp)
             end if
             if(.not. algo_param%HPP) then
                call SigToPK1(Sig,Def,PK1)
             end if
             if (algo_param%Nloc .and. .false.) then !---- Mecanique non locale resolution implicite (A REVOIR)
                ! TODO revoir possiblite d'introduire de l'implicite
                call Nloc_call(load_n,load_incr,ind_tps,"after",algo_param%nloc_implicit)
             end if
           end if

           if (algo_param%Diffusion) then !-----DIFFUSION
             call behaviorD(FluxD,FluxD0,GradQD,GradQD0,local_loadD%t0(6*algo_param%nVarD+1:),&
                           local_loadD%dt(6*algo_param%nVarD+1:),local_loadD%t_load,KSTEP)
           end if

           VCinfo%nSub_lam_pas  = VCinfo%nSub_lam_pas  + nSub_laminate_temp
           VCinfo%nIncr_lam_pas = VCinfo%nIncr_lam_pas + nIncr_laminate_temp
           VCinfo%nIt_lam_pas   = VCinfo%nIt_lam_pas   + nIt_laminate_temp

           !----------------------------------------------------------------------
           !                          DISPLAY to fichier.log (EVERY 10 ITERATIONS)

           
           !nIt=nIt+1;      ! incrementat  number of iterations 
             ! UNCOMMENT to output criteria at each iteration to .log 
             !  if (nrank==0) call write_crit_log(Flog)
             !  if (nrank==0) write (Flog,"(/,A,/)")"-----------------------------"
             ! UNCOMMENT to output nIt, crit_b%eq at each iteration to .log 
             ! if(nrank==0)then
             !     write (Flog,"(I8,E15.8)") nIt,crit_b%eq
             ! end if

!-- deformation energy output
!block
!   real(mytype) :: energy_loc,energy
  
!   energy_loc = sum(SIG*DEF)
!   call MPI_Allreduce(energy_loc,energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

!   if(nrank==0)then
!     if (nIT .eq. 0 .and.load_incr==0) then
!       open(unit=78941, file=trim(fic_vtk)//"_energy.log",form="formatted",status="replace", action="write")
!     end if
!     write (78941,"(E15.8,I8)") Energy,nIT
!   end if

!end block


           if((nIt .gt. 0)  .and. modulo(nIt,10)==0) then  ! display every  10 iterations
              call MPI_BARRIER(MPI_COMM_WORLD,ierror) !needed?
              t1 = MPI_WTIME()

              if(nrank==0)then
                 write (Flog,"(/,A,/)")"-----------------------------"
                 if(load_incr==0) then
                    write (Flog,"(A,I8,A)")"ficticious loading step ",ind_tps,".5"
                 else
                    write (Flog,"(A,I8)")"loading step ",ind_tps
                 end if
                 write (Flog,"(I8,A)") nIt, " iterations"
                 !sortie fichier log
                 call write_crit_log(Flog)
                 !sortie ecran
                 call write_crit_log(OUTPUT_UNIT)
                 write (Flog,"(A,E15.8)") "Total time after initialisation :",t1-times_b%startalgo
                 write (Flog,"(A,E15.8)") "Time for the last 10 iterations :",t1-t10
                 write (Flog, "(A,/,A)") "Temps :",&
                   "FFT_contrainte IFFT_deformation     umat        apply_green   eval_criteres  ecriture"
                 !! A afficher en sortie
                 write(OUTPUT_UNIT,"(/,A,/)")"-----------------------------"
                 if(load_incr==0) then
                    write(OUTPUT_UNIT,"(A,I8,A)")"ficticous loading step",ind_tps,".5"
                 else
                    write(OUTPUT_UNIT,"(A,I8)")"loading step ",ind_tps
                 end if
                 write(OUTPUT_UNIT,"(I8,A)") nIt, " iterations"
                 write (Flog,"(6E15.8)") times_f%fft, times_f%ifft, times_m%behavior, times_g%apply, times_g%crit, times_io%wvtk
              end if
              t10 = MPI_WTIME()
           end if

        end do
        !                                                                  END WHILE LOOP
        !================================================================================
        !                                                              FORCED CONVERGENCE 
        !                                              
        !
        !------------FORCED UPDATE OF INTERNAL VARIABLES        
        !
!print *,"AP WHILE - BEFORE CVFOR"
        if (testCVFOR) then

          !--Recompute the behaviour to evaluate internal variables 
          !  using the best solution in strain.
          if (algo_param%HPP_nsym)  then
             ! take into account the full displacement gradient in 
             ! HPP calculation
             call behavior(SIG,SIG0,CVFORsauv%Def,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                        local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                        nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp,&
                        Def_nsym=Def_nsym,Def_nsym0=Def_nsym0)
          else
             call behavior(SIG,SIG0,CVFORsauv%Def,DEF0,local_load%t0(TPRED_deb:TPRED_fin),&
                        local_load%dt(TPRED_deb:TPRED_fin),local_load%t_load,KSTEP,&
                        nSub_laminate_temp,nIncr_laminate_temp,nIt_laminate_temp)
          end if

          !--Update internal variables and stress at beginning of time step
          do i=1,size(MattotP)
             MattotP(i)%VarInt0(1:MattotP(i)%nVarInt,:)=MattotP(i)%VarInt(1:MattotP(i)%nVarInt,:)
          end do  
          Sig0=Sig

         !-- Initialisation choice for the next forced CV iteration CV forcee
         !   pour le cas init_cvfor='best'
         if (algo_param%init_cvfor=='best') then
             Def = CVFORsauv%Def
         end if

         !-- Update counter
         if (nrank==0) write (Flog,"(A,I4,A,E15.8)") &
                       "countCVFOR= ", countCVFOR, " / crit eq = ", CVFORsauv%critmin
         if (nrank==0) write (OUTPUT_UNIT,"(A,I4,A,E15.8)") &
                       "countCVFOR= ", countCVFOR," / crit eq = ", CVFORsauv%critmin
         CVFORsauv%critmin=1.e100_mytype
         countCVFOR = countCVFOR + 1 
        end if

        !
        !------------EXIT PROGRAM TEST ON FORCED CONVRGENCE 
        !
        if (countCVFOR .eq. algo_param%Ncvfor) then 
           if (nrank == 0) then
              write (OUTPUT_UNIT,"(A,I8,A)") "More than ",algo_param%Ncvfor," forced convergences, computation stopped."
              write (OUTPUT_UNIT,"(/,A,/)")"-----------------------------"
              write (Flog,"(A,I8,A)") "More than ",algo_param%Ncvfor," forced convergences, computation stopped."
              write (Flog,"(/,A,/)")"-----------------------------"
           end if
           call amitex_abort("Maximim number of forced convergences reached",2,0)
        end if

        end do
        !                                                     END FORCED CONVERGENCE LOOP
        !================================================================================

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  times_b%unpas = times_b%unpas + MPI_WTIME()- t_ini

  !print *,"END UNPAS_NL_BASE"
end subroutine unpas_NL_base



!==================================================================================
!==================================================================================
!                         SUBROUTINE LOG_AFTER_UNPAS
!
!> Write informations after a loading step in resolution_NL_base
!!
!!
!==================================================================================
subroutine log_after_unpas(load_incr,ind_tps, nIt,&
               t_fft0,t_ifft0,t_behavior0,t_green0,t_crit0,t_wvtk0,t_wstd0,t_unpas0,&
               testLaminate,msg)
  
  implicit none
  
  integer, intent(in)      :: load_incr,ind_tps, nIt
  real(mytype), intent(in) :: t_fft0,t_ifft0,t_behavior0,t_green0,&
                              t_crit0,t_wvtk0,t_wstd0,t_unpas0
  logical, intent(in)      :: testLaminate
  character(len=*),optional,intent(in) :: msg
  integer                  :: ierror  ! mpi error
  double precision         :: t1      ! time measurement


  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  t1 = MPI_WTIME()

  if (nrank==0) then
    if (present(msg)) write(OUTPUT_UNIT,"(A)") msg
    write(OUTPUT_UNIT,"(/,A)")     "------------------------------------------ ("//trim(simu_name)//")"
    write(OUTPUT_UNIT,"(A,I8)")"Number of iterations    :",nIt
    write(OUTPUT_UNIT,"(A,/,/)")"=================================================="
    !! In the log
    write(Flog,"(/,A)")"--------------------------------------------------"
    if(load_incr==0) then
            write(Flog,"(A,I8,A)")"End of ficticious loading step :",ind_tps,".5"
    else
            write(Flog,"(A,I8)")"End of loading step :",ind_tps
    end if
    write(Flog,"(A,I8)")"Number of iterations    :",nIt 
    call write_crit_log(Flog)

    !! Info umatLaminate
    if (testLaminate) then
       write(Flog,"(A)") "UmatLaminate : "
       write(Flog,"(A,F12.4)")"          - Average number of iterations :",&
                   real(VCinfo%nIt_lam_pas)/real(VCinfo%nCall_lam_pas)
       !! Info éventuels pilotages effectués pour umatLaminate
       if (VCinfo%nSub_lam_pas > 0) then
          write(Flog,"(A,I12,A,I12,A,F12.4,A)")"          - Number of time step subdivisions & 
            & /Number of calls  :",VCinfo%nSub_lam_pas,&
              " / ",VCinfo%nCall_lam_pas,"  (",100*real(VCinfo%nSub_lam_pas)/real(VCinfo%nCall_lam_pas),"%)"
          write(Flog,"(A,F12.4)")"          - Average number of sub-time steps used for &
            & subdivisons      :", real(VCinfo%nIncr_lam_pas)/real(VCinfo%nSub_lam_pas)
       end if
    end if
    
    write(Flog,"(A,E15.8)") "Total time             :",t1-times_b%startalgo
    write(Flog, "(A,/,A)") &
      "FFT             IFFT           umat           apply_green    eval_crit      write VTK      write STD",&
      "Total times :"
    t1 = MPI_Wtime()
    write(Flog,"(7E15.8)") times_f%fft, times_f%ifft, times_m%behavior, times_g%apply,&
                           times_g%crit, times_io%wvtk, times_s%wtot
    write(Flog,"(A,/,7E15.8)")"Loading step times:",times_f%fft - t_fft0, &
         times_f%ifft - t_ifft0, times_m%behavior - t_behavior0, times_g%apply - t_green0,&
         times_g%crit - t_crit0, times_io%wvtk - t_wvtk0, times_s%wtot - t_wstd0
    write(Flog,"(A,E15.8,A,E15.8,A)")"Time in loading step time (and in unpas_NL_base) :",&
         t1-times_b%startstep,"(",times_b%unpas - t_unpas0,")"
    write(Flog,"(A)")"=================================================="
    write(Flog,"(A,/,/)")"=================================================="
    
  end if

end subroutine log_after_unpas

!==================================================================================
!==================================================================================
!                         SUBROUTINE LOG_FINAL_TIMES
!
!> Write final times informations after resolution_NL_base
!!
!!
!==================================================================================
subroutine log_final_times(nitTot)

  implicit none
 
  integer, intent(in)      :: nitTot
  integer                  :: ierror  ! mpi error
  double precision         :: t1      ! time measurement
  

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  t1 = MPI_WTIME()
  
  if (nrank==0) then
     write(Flog,"(/,A)")     "========================================== ("//trim(simu_name)//")"
     write(Flog,"(A)")       "Algorithm times"
     write(Flog,"(A)")       "---------------"    
     write(Flog,"(A,E15.8)") "Cumulated time for  behavior      : ", times_m%behavior
     write(Flog,"(A,E15.8)") "Cumulated time for  eval_crit     : ", times_g%crit
     write(Flog,"(A,E15.8)") "Cumulated time for  apply_green   : ", times_g%apply
     write(Flog,"(A,E15.8)") "Cumulated time for  FFT           : ", times_f%fft
     write(Flog,"(A,E15.8)") "Cumulated time for  iFFT          : ", times_f%ifft
     write(Flog,"(A,E15.8)") "Cumulated time for  write VTK     : ", times_io%wvtk
     write(Flog,"(A,E15.8)") "Cumulated time for  write STD     : ", times_s%wtot
     write(Flog,"(A)")       "-----------------------------------------------------"    
     write(Flog,"(A,E15.8)") "Sum                               : ",&
          times_m%behavior + times_g%crit + times_g%apply + times_f%fft &
          + times_f%ifft + times_io%wvtk + times_s%wtot
     write(Flog,"(A,E15.8)") "Cumulated time in unpas_NL_base   : ", times_b%unpas
     write(Flog,"(A,E15.8)") "Total time in resolution_NL_base  : ", t1-times_b%startalgo

     write(Flog,"(A,/, I12)") "total number of iterations" ,  nitTot ! used in  validation/comparaison.f90

     write(OUTPUT_UNIT,"(/,A)")     "========================================== ("//trim(simu_name)//")"
     write(OUTPUT_UNIT,"(/,A)")     "Algorithm times"
     write(OUTPUT_UNIT,"(A)")       "---------------"    
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  behavior      : ", times_m%behavior
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  eval_crit     : ", times_g%crit
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  apply_green   : ", times_g%apply
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  FFT           : ", times_f%fft
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  iFFT          : ", times_f%ifft
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  write VTK     : ", times_io%wvtk
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time for  write STD     : ", times_s%wtot
     write(OUTPUT_UNIT,"(A)")       "-----------------------------------------------------"    
     write(OUTPUT_UNIT,"(A,E15.8)") "Sum                               : ",&
          times_m%behavior + times_g%crit + times_g%apply + times_f%fft &
          + times_f%ifft + times_io%wvtk + times_s%wtot
     write(OUTPUT_UNIT,"(A,E15.8)") "Cumulated time in unpas_NL_base   : ", times_b%unpas
     write(OUTPUT_UNIT,"(A,E15.8)") "Total time in resolution_NL_base  : ", t1-times_b%startalgo
     
     write(OUTPUT_UNIT,"(/,A,/, I12,/)") "Total number of iterations :" ,  nitTot
  end if

end subroutine log_final_times

!==============================================================================
!      SUBROUTINE WRITE_CRIT_LOG
!
!>  WRITE CRITERIA TO THE LOG FILE
!!      if diffusion is active, display the max values for each variable
!!      otherwise, display 0
!!
!!  crit : type(CRIT_BASE) variable from this module (NL_base_mod.f90)
!!
!!  \param[in]    FU, logical unit (Flog ou OUTPUT_UNIT)
!!
!------------------------------------------------------------------------------
subroutine write_crit_log(FU)

   implicit none 

   integer, intent(in)                  :: FU

   real(mytype)                         :: A,B,C,D
   character(len=5),parameter           :: FMT_real="E15.8"             !<real number output format
!   character(len=6),parameter          :: FMT_real="E25.18"            !<very precise real number output format


!> If diffusion is not activated, the associated criteria are displayed as 0.
   if (algo_param%Diffusion) then   
      A = maxval(crit_b%eqD)
      B = maxval(crit_b%FluxDMoy)
      C = maxval(crit_b%CptbD) 
      D = maxval(crit_b%GradQDMoy)
   else  
      A = 0.
      B = 0.
      C = 0. 
      D = 0.
   end if

   !> wrtie to file
   write(FU,"(A)")     "------------------------------------------ ("//trim(simu_name)//")"                               
   write(FU,"(2(A,"//FMT_real//"))") &
                      "equilibrium criteria                (meca,diff) :", crit_b%eq,' , ', A
   write(FU,"(2(A,"//FMT_real//"))") &
                      "mean stress criteria                (meca,diff) :", crit_b%SigMoy,' , ', B
   write(FU,"(2(A,"//FMT_real//"))") &
                      "compatibility criteria              (meca,diff) :", crit_b%Cptb,' , ', C
   write(FU,"(2(A,"//FMT_real//"))") &
                      "mean strain criteria                (meca,diff) :", crit_b%DefMoy,' , ', D

end subroutine write_crit_log
!==============================================================================


end module NL_base_mod
