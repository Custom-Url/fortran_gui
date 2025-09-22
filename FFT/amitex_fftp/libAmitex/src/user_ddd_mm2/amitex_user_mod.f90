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
  
  use dd_mod

#ifdef OPENMP
  use omp_lib
#endif

  implicit none

  private

! Declare the procedure used in amitex_fftp before the resolution (REQUIRED) 
! Initialize global user-variables

  public :: init_user_variables

! Declare global user-variables as public

  public :: DEIGS, key_ascii,idt,timestep_ratio,ratio,nprocDD
  public :: tmdc0,tmdc1,tmdc2,tmdc3,tmdc4,tmdc5,tmdc6,tcomm

                       
! Declare global procedures as public (used in resolution_user_mod or standard_user_mod)

  public :: write_header_bin

!!------------------------------------------------------------------------------
!>                                                                 PUBLIC FIELDS 

real(mytype),allocatable,dimension(:,:)    :: DEigs 
                                             !< EigenStrain  - WARNING pas (ntot,6)!!
                                             !6 components of the straintensor: xx,yy,zz,xy,xz,yz

!!------------------------------------------------------------------------------
!>                                                               OTHER VARIABLES 

logical             :: key_ascii=.false.     !< =.false. pour ecriture fichier restart 
integer             :: timestep_ratio        !< ratio between microMegas et amitex time step, should be >=1
real(mytype)        :: ratio                 !< ratio between mM grid and amitex grid
integer             :: idt                   !< indice boucle construction des eigenstrains 
integer             :: nprocDD               !< number of procs requested to run the DD code mM
real(mytype)        :: tmdc0,tmdc1,tmdc2,tmdc3,tmdc4,tmdc5,tmdc6,tcomm
                                             !< temps d'execution

contains


!!------------------------------------------------------------------------------
!>                                         ALLOCATE AND INITIALIZE PUBLIC FIELDS 

subroutine init_user_variables()

  implicit none
  
  integer         :: alloc_stat, ierror,i
  integer         :: restart_key             !< key to restart calculation
  real(mytype)    :: deltat,time0,time1      !< time step calculated from tfinal and nincr defined in the loading file
                                             !< here it is read again from an other file
  character(len=1000) :: mm_exe              !< path to the executable tu run the DD code mM                                       

  allocate(Deigs(6,xsize(1)*xsize(2)*xsize(3)),stat=alloc_stat)  !! WARNING pas (ntot,6)!!
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (init_user_variables 1)",2)

  DEigs = 0._mytype
  
  tmdc0 = 0._mytype
  tmdc1 = 0._mytype
  tmdc2 = 0._mytype          
  tmdc3 = 0._mytype
  tmdc4 = 0._mytype
  tmdc5 = 0._mytype
  tmdc6 = 0._mytype
  tcomm = 0._mytype
  tcomm1 = 0._mytype
  tcomm2 = 0._mytype


  !!2****************************************************INITIALISATION MM
  !restart_key, if restart_key different from 1 we set the restart_key equal to 
  restart_key=int(user_param%P_real(2))  
  if (restart_key /= 1 ) restart_key=0

  
  kkdebug=int(user_param%p_real(4))

  timestep_ratio=int(user_param%p_real(5)) 
  time0 = load(1)%time(1)
  time1 = load(1)%time(2)
  deltat=time1-time0
  print * , " The time increment in Amitex is : ", deltat," and "
  print *, " the update of long range stress field through the FFT solution is done every ",timestep_ratio," steps" 

  nprocDD=int(user_param%p_real(6))
  mm_exe=user_param%p_string(1)

  key_cputime=int(user_param%p_real(7))
  
  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "before open mM. rank=",nrank  
  end if
  
  ! Opening communication with the DD code microMegas 

  call open_mM(nprocDD,mm_exe)

  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after open mM. rank=",nrank  
  end if


  ! At present we assume our material to be homogeneous and isotropic     
  ! send elastic constants to the DD code microMegas
  call send_matrix_of_elasticity(MattotP(1)%Coeff(1,1),MattotP(1)%Coeff(2,1))

  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after send elastic constants. rank=",nrank  
  end if

  
  call send_restart_key(restart_key)
  
  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after send_restart_key. rank=",nrank  
  end if

  ! Received from the DD code microMegas some key parameters  
  call get_global_data(ratio)
  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after get global data. rank=",nrank  
  end if
  
  call send_global_data(ratio)
  
  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after send global data. rank=",nrank  
  end if
  !!2*********************************************FIN INITIALISATION MM

  !!3****************************************************CALCUL EIGENSTRAIN PAS 0

  ! Calcul eigenstrain from the initial dislocation configuration
  !

  if (restart_key/=1) call eigenstrain_calculation(DEigs,ratio,1,0)
  
  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after eigenstrain calculation at step 0. rank=",nrank  
  end if

  !transmet timestep_ratio a mM pour savoir s'il doit appeler ou non amitex
  call send_timestep(deltat,timestep_ratio)
  
  if (kkdebug > 0 ) then
    call mpi_barrier(mpi_comm_world,ierror)
    if (nrank==0) print *, "after time step exchange. rank=",nrank  
  end if

  !!3****************************************************FIN CALCUL EIGENSTRAIN PAS 0
  
  !!4**************************************************** TRANFERT EIGENSTRAIN VAR.INT MM PAS 0
            
  ! transfer to MattoP%VarInt0 (before un_pas the eigenstrain transfer is done to VarInt0  )
  !!>MattotP(1)%Varint(:,:) variables internes dimensions: (nVarInt,NPosition)
  do i=1,size(MattotP)
     MattotP(i)%VarInt0(1:MattotP(i)%nVarInt,:)=MattotP(i)%VarInt0(1:MattotP(i)%nVarInt,:)+DEigs(1:MattotP(i)%nVarInt,:)
  end do 
        
  !!4*********************************************FIN  TRANFERT EIGENSTRAIN VAR.INT MM PAS 0
  

end subroutine init_user_variables


subroutine write_header_bin(nomFic, nx,ny,nz)

  implicit none
  character(len=*), intent(in)  :: nomFic
  integer, intent(in)           :: nx,ny,nz
  integer                       :: Fvtk
  integer                       :: ios 
  integer(kind=8)               :: ntot

  ntot = nx
  ntot = ntot*ny
  ntot = ntot*nz

!ouverture du fichier de sortie
    open(newunit=Fvtk, file=nomFic,form="formatted", &
                              status="replace", action="write",iostat= ios)
      if ( ios /= 0 ) then
        print *," Probleme a l'ouverture (fichier: ",nomFic,")",ios
        stop 
      end if
    
    write(Fvtk,"(I12)")ntot
    write(Fvtk,"(A)")"float"
  close(Fvtk)


end subroutine write_header_bin



end module amitex_user_mod
