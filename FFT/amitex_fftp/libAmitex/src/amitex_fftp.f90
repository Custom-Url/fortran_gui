!===============================================================================
!>                  MAIN PROGRAM OF AMITEX_FFTP CODE
!!
!! THE MAIN PROGRAM HANDLES:
!!      1/ Allocating and initializing, from input files,
!!         all necessary data for calculation in SIMU(1)
!!      2/ Associating pointers corresponding to fields in SIMU(1)
!!      3/ Launching the calculation (e.g., resolution_NL_base)
!!         using these pointers
!!
!! NOTATIONS (applicable to input code and fields):
!!         Sigma and Def (HPP): CAST3M (Voigt + order 11 22 33 12 13 23)
!!         Pi and grad(u) (GDEF) (order 11 22 33 12 13 23 21 31 32 without notation)
!!
!! PROGRAMMING RULES: use types defined in 2decomp
!!         remove double precision => real(mytype),
!!             except for double precision output from the MPI library
!!             (e.g., time output from MPI_WTIME())
!!         in MPI functions, replace MPI types:
!!             MPI_DOUBLE_PRECISION => real_mytype
!!         creating complex numbers: 1 + i = cmplx(1,1,kind=mytype)
!!         pay attention to constants => pi = 3.14.._mytype
!!         float => real(nx,mytype)
!!         When passing Fourier transforms DefF and SigF to other functions
!!         do not forget to specify that arrays go from fft_start to fft_end
!!         otherwise Fortran allocates them from 1 to fft_end-fft_start + 1
!!
!! WARNING: in single precision, FFT must be performed on a normalized field
!!         (otherwise, risk of exceeding maximum real value)
!!
!!==============================================================================
!!
!! USAGE OF EXECUTABLE:
!!      amitex_fftp -help
!!
!!==============================================================================

program amitex_fftp



  use ISO_FORTRAN_ENV

  use MPI             
  use decomp_2d,      only : mytype, nrank, xsize, decomp_2d_finalize
  use decomp_2d_fft,  only : decomp_2d_fft_finalize

  use io_amitex_mod,  only : write_stdout0 ! function
  use io2_amitex_mod, only : print_coeff_varint_vtk ! function
  use green_mod,      only : FREQ ! pointer
  use param_algo_mod, only : algo_param, user_param ! pointers
  use material_mod ,  only : VCInfo, testComposite  ! pointers
  use error_mod,      only : amitex_abort !function
  use resolution_mod, only : resolution_NL_base ! function
  use amitex_mod,     only : Sig, Flog, &       ! variables
                             AMITEXenv, LOGO    ! pointers
  use field_mod,      only : field_isNaN_real,& ! function 
                             ntotFP             ! pointer
  use non_local_mod,  only : init_user_nloc_variables_call ! function
  use simu_mod,       only : SIMU, & ! variable
                             nullify_pointers_simu, init_simu_command_line,&
                             associate_pointers_simu !functions

 !! USERS
  use amitex_user_mod,     only : init_user_variables  ! function
  use resolution_user_mod, only : resolution_user      ! function

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================
  implicit none

!!------------------------------------------------------------------------------
!>                                                                    SIMULATION  
  integer                                           :: choix_simu=1 
                                        !< 1 : 1 standard simulation 
                                        !! 2 : (dev) 2 simu simulations, the 2nd is performed

!!------------------------------------------------------------------------------
!>                                                                  PARALLELISM 

  integer :: ierror                     !< error related to MPI functions


!!------------------------------------------------------------------------------
!>                                                          EXECUTION TRACKING 

  double precision :: t0,t1,t2          !< double precision instead of real_mytype required by MPI_WTIME
                                        !! t0 : start time of the algorithm
                                        !! t1 : time to measure duration of specific functions
                                        !! t2 : start time of the program (for initialisation duration)

!!------------------------------------------------------------------------------
!>                                                                 VERIFICATIONS

  logical :: resolution=.true.          !< launch the resolution
  logical :: check                      !< verifications before resolution
  integer :: coeff2print =0             !< index of coefficient to output (command line) 
  integer :: varint2print=0             !< index of varint ot output (command line) 
  logical :: help=.false.               !< check to display help and exit the code

!!------------------------------------------------------------------------------
!>                                                                        MISC

  integer            :: alloc_stat      !< error during memory allocation
  character(len=200) :: version         !< version number written in .log
  logical            :: bool0 

  t0 = 0 ! avoid gcc-warning (unused var)
  
!!==============================================================================
!>                                                            MPI INITIALISATION 
!>                                                        + nrank INITIALISATION
!>                            -> necessary as long as 2decomp is not initialised
!!==============================================================================
  call MPI_INIT(ierror)
  t2 = MPI_WTIME()
  call MPI_COMM_RANK(mpi_comm_world,nrank,ierror)

  ! write logo on stdout
  call write_stdout0(achar(10)//LOGO//achar(10))

  ! write version on stdout
  version="10.0.0"
  call write_stdout0("VERSION "//trim(version))

!!==============================================================================
!!                                                     SIMULATION INITIALISATION 
!!==============================================================================

  if (choix_simu == 1) allocate(SIMU(1))
  if (choix_simu == 2) allocate(SIMU(2))
  
  call nullify_pointers_simu()  
  call init_simu_command_line(SIMU(1),version, help,coeff2print,varint2print, Igrid = 1)

  ! end program if  "amitex_fftp -help" is used
  if (help) then
     goto 666
  end if

  if (choix_simu == 2) then
     call nullify_pointers_simu()
     call init_simu_command_line(SIMU(2),version, help,coeff2print,varint2print,&
                 "cmd_files/bilayer_2mat_2zone_xmla.in", Igrid = 1)
  end if 
  
!!==============================================================================
!!                                                             SIMULATION CHOICE
!!==============================================================================
  call nullify_pointers_simu()

  select case(choix_simu)
  case(1)
     call associate_pointers_simu(SIMU,1)
  case(2)
     call associate_pointers_simu(SIMU,2)   
  end select
  
!!==============================================================================
!!                                         CHECK BEFORE STARTING THE CALCULATION
!!                                TODO : create a check_simu function in simu_mod
!!==============================================================================
! option to disable all checks 
  check=.true.

  if (check) then

  ! check NaN Frequency
  bool0=.false.
  call field_isNaN_real(Freq,3*ntotFP,bool0)
  if (bool0) call amitex_abort("Presence of NaN in Freq (amitex)",2)

  ! check "Composite Voxels" AND "None-Local Models"
  if (testComposite .AND. algo_param%Nloc) then 
      call amitex_abort("Simultaneous use of non-local models &
              & and composite voxels is under development (amitex)",0,0)
  end if

  ! check "User (full algo or functions)" AND "Non-Local Models"
  if (user_param%test .AND. algo_param%Nloc) then 
      call amitex_abort("Simultaneous use of non-local models &
              & and a <User> node in the algorithm file (.xml) (amitex)",2)
  end if

  end if
  call write_stdout0("After a few check") 
  
!!==============================================================================
!!                                                         VTK OUTPUT FOR CHECKS 
!!                                  (option -coeff2print i and/or varin2print i)
!!                                 Sig is used to construct the coefficent field    
!!                                               but the calculation is not used
!!==============================================================================
  ! In pure diffusion, Sig is not allocated => allocate a field with one component
  if(.not. allocated(SIMU(1)%Sig)) then    
        allocate(SIMU(1)%Sig(xsize(1)*xsize(2)*xsize(3),1),stat=alloc_stat)
        if(alloc_stat /=0) call amitex_abort("insufficient avalable memory (amitex13)",2)
        Sig => SIMU(1)%Sig
  end if

  call print_coeff_varint_vtk(coeff2print,varint2print)
  
  if ((coeff2print > 0) .OR. (varint2print > 0)) then
    resolution=.false.
  end if 

!!==============================================================================
!!                                                           USER INITIALISATION
!!==============================================================================

if (user_param%test) then
  call write_stdout0("Before init_user_variable") 
  call init_user_variables()
  call write_stdout0("After init_user_variable") 
end if

!!==============================================================================
!!                                                      NLOC USER INITIALISATION
!!      allocation and initialisation of additional variables for user_nloc1,2,3
!!                                                                  if necessary
!!==============================================================================

if (algo_param%Nloc) then
    call write_stdout0("Before init_user_nloc_variables_call") 
    call init_user_nloc_variables_call()
    call write_stdout0("After init_user_nloc_variables_call") 
end if

!!==============================================================================
!!                                                                     USER ALGO        
!!==============================================================================

if (resolution .AND. user_param%test .AND. trim(user_param%algo) == "user") then

  !> total initialisation time   
  call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
  t0 = MPI_WTIME()
  if(nrank==0) then
      write(Flog,fmt="(A,E15.8)") "Duration TOTAL INITIALISATION (s) : ", t0 - t2
      write(OUTPUT_UNIT,fmt="(A,E15.8)") "Duration TOTAL INITIALISATION (s) : ", t0 - t2
  end if

  !> Launch of USER algorithm 
  call resolution_user()

end if


!!==============================================================================
!!                                                            STANDARD ALGORITHM
!!==============================================================================
!!
!!         FOR THE CALCULATION: THE ALGORITHM MAY USE, IN ADDITION TO THE PASSED 
!!          PARAMATERS, THE FOLLOWING STRUCTURES ASSIGNED DURING INITIALISATION:
!!      
!!      MattotP         array of MATERIAL strucutes     (material_mod.f90)
!!              -> material description 
!!
!!      load            array of LOADING structures     (loading_mod.f90))
!!              -> loading description      
!!
!!      initValExt      INITLOADEXT strucute            (loading_mod.f90)
!!              -> description of external load insitialisation 
!!
!!      extract         PARAMEXTRACT structure          (loading_mod.f90)
!!              -> description of VTK outputs
!!
!!      algo_param      PARM_ALGO structures            (param_algo_mod.f90)
!!              -> description of algorithm parameters     
!!
!!      grid            GRID_DIM structure              (green_mod.f90)
!!              -> grid description         
!!
!!==============================================================================
call write_stdout0("Just before standard algorithm") 
if(resolution .AND. .NOT. (user_param%test .AND. trim(user_param%algo) == "user" )) then

  !> total initialisation time   
  call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
  t0 = MPI_WTIME()
  if(nrank==0) then
      write(Flog,fmt="(A,E15.8)") "Duration TOTAL INITAIALIZATION (s) : ", t0 - t2
      write(OUTPUT_UNIT,fmt="(A,E15.8)") "Duration TOTAL INITAIALIZATION (s) : ", t0 - t2
  end if

  !> Launch of algorithm      
  if (algo_param%scheme_type == "basic_scheme") then
     call resolution_NL_base()
  end if

  if(nrank==0) then
     if (testComposite) then
        write(Flog,"(A)") "------------------------------"
        write(Flog,"(A)") "UmatLaminate "
        write(Flog,"(A,F12.4)") "          - average number of iterations per voxel (laminate algo): " ,&
                    real(VCinfo%nIt_lam)/real(VCinfo%nCall_lam)
        write(Flog,"(A,I12,A,I12,A,F12.4,A)") &
          "          - Number of time step subdivisions / Number of calls: ",VCinfo%nSub_lam,&
          " / ",VCinfo%nCall_lam,"  (",100*real(VCinfo%nSub_lam)/real(VCinfo%nCall_lam),"%)"
        if (VCinfo%nSub_lam >0) write(Flog,"(A,F12.4)") &
          "          - Average number of subdivisions in PilotageLaminate: ",&
          real(VCinfo%nIncr_lam)/real(VCinfo%nSub_lam)
        write(Flog,"(A,/)") "------------------------------"
     end if
  end if

end if


!!==============================================================================
!!                                            DEALLOCATIONS AND END OF PROGRAMME
!!==============================================================================
  call write_stdout0("Deallocations and end of program")

  if(allocated(SIMU)) deallocate(SIMU)

  if(allocated(AMITEXenv%path)) deallocate(AMITEXenv%path)


  call decomp_2d_fft_finalize
  call decomp_2d_finalize 

  call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
  t1 = MPI_WTIME()
  if (nrank==0) then
    write(Flog,"(A,E15.8)") 'Total time (excluding initialisation) =', t1-t0
    write(OUTPUT_UNIT,"(A,E15.8)") 'Total time (excluding initialisation) =', t1-t0
    close(Flog)
  end if
666  call MPI_FINALIZE(ierror)
end program amitex_fftp

!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
