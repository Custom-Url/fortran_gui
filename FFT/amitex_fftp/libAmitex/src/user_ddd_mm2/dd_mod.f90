!==============================================================================
!       MODULE DD_MOD : 
!> Definitions, functions and subroutines used to communicate from amitex to DD
!> (at present microMegas code) and vice-versa 
!!
!! structure 
!! DD code : receives the elastic constants from amitex, generates the dislocation 
!! microstructure, sends the information about the areas swept by the dislocations 
!! to amitex, send the integration points to amitex (segment centers in micromes),
!! receives the stress at the integration points , moves the dislocation segments
!! send again the information about the swept areas etc ...
!!
!! AMITEX- : send the elastic constants to the DD code, receives information about 
!! the swept areas and reconstructs the eigenstrains, calculate the stress in 
!! the simulation box, interpolate the stress at the integration points received
!! from DD (segment centers for microMegas), send the stress values at the 
!! integration points to the DD code, receive again the information about the swept 
!! ares, etc ...  
!! 
!!
!! subroutines 
!! open_mM :  start the communication with the DD code microMegas           
!!                        
!! send_matrix_of_elasticity : send the elastic constants to the DD code microMegas  
!!      REMARQUE : at present micromegas can deal only with homogeneous and isotropic materials           
!! initParam_umat :       
!! deallocate_MattotP :   
!!
!! REMARQUE 1
!!  
        
        
module dd_mod

  use MPI 
  use error_mod            
  use decomp_2d
  use green_mod  
  
  implicit none
  real(mytype), PARAMETER  :: Pi = 3.1415927,M_CUTOFF = 1.75
  integer                  :: mMintercomm
  integer                  :: fesmpierr,fesmpistatus(MPI_STATUS_SIZE) 
  integer                  :: fesmpirecvbufsize=0,fesmpirecvbufpos,fesmpisendbufsize=0,fesmpisendbufpos
  character, pointer       :: fesmpirecvbuf(:),fesmpisendbuf(:)
  integer                  :: ierror,kkdebug  
  real(mytype)             :: b,h,avalue,box(3)
  integer                  :: flagfreesurf      !< flag 1 if there are free surfaces
  real(mytype)             :: tcomm1,tcomm2
  integer             :: key_cputime           !< key for evaluation of the time spent in the different subroutine of the main loop
                                               !< 1=true, 0 or other values = false
  


  
!***********************************************************************
!Dist2psData: structure storing the position of a point with respect to a segment
!//
!//    l2: square of the length of the segment
!//    dn2: ( (P-C).n )^2
!//    dt2: ( (P-C).t )^2
!//    sdt2: sign( (P-C).t )
!//
!//      legend:
!//        P: point
!//        S0, S1: segment extremities
!//        C: segment center
!//        n: unit vector normal to the segment
!//        t: unit vector colinear to the segment oriented from S0 to S1
!//
!//                            P
!//                           x
!//                           ^
!//                           |
!//                           |dn2
!//                           |             ^
!//       S0           C      v     S1     n|
!//      x------------+------------x        ---->
!//                   <------->                t
!//                      dt2
!//      <------------------------->
!//                  l2
!//
  type Dist2psData
    real(mytype)               :: l2 
    real(mytype)               :: dn2 
    real(mytype)               :: dt2 
    real(mytype)               :: sdt2 !sign of dt2           

end type Dist2psData

contains

!########################################################################################################
!> Receive a message of unknown size from MPI communicator mMintercomm, get the message size, allocate the buffer 
!>(if the buffer is smaller than the received message the buffer is reallocated) and save the message
!> received in the buffer 
!########################################################################################################
!> \ingroup DCM
subroutine fesmpipackrecv
  implicit none
  integer :: packsize
  call MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mMintercomm, fesmpistatus, fesmpierr)
  call MPI_Get_count(fesmpistatus, MPI_PACKED, packsize, fesmpierr)
  if (packsize>fesmpirecvbufsize) then
    call fesmpireallocrecvbuf(packsize)
  endif
  call MPI_Recv(fesmpirecvbuf, packsize, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, mMintercomm, fesmpistatus, fesmpierr)
  fesmpirecvbufpos=0
end subroutine fesmpipackrecv

!########################################################################################################

!########################################################################################################
!> The buffer is reallocated before receiving the message through MPI communicator because the buffer is 
!> smaller than the message received fro MPI communicator
!########################################################################################################
!> \ingroup DCM
subroutine fesmpireallocrecvbuf(newbufsize)
  implicit none
  integer :: newbufsize, ierr
  character, pointer :: tmpbuf(:)

  newbufsize=2*newbufsize

  if (associated(fesmpirecvbuf)) then
    allocate(tmpbuf(newbufsize), stat=ierr)
    if(ierr/=0) stop 'FATAL ERROR: problem while allocating fesmpi receive buffer'
    tmpbuf(1:fesmpirecvbufsize)=fesmpirecvbuf(1:fesmpirecvbufsize)
    deallocate(fesmpirecvbuf)
    fesmpirecvbuf=>tmpbuf
  else
    allocate(fesmpirecvbuf(newbufsize))
  endif
  fesmpirecvbufsize=newbufsize
end subroutine fesmpireallocrecvbuf

!########################################################################################################


!########################################################################################################
!> Pack an array of doubles "buf" of size n to be send through MPI communicator mMintercomm (type PACKED)
!########################################################################################################
!> \ingroup DCM
subroutine fesmpipack_doublearray(buf,n)
  implicit none
  real(mytype) :: buf(:)
  integer :: n
  integer :: packsizeinc
  call MPI_Pack_size(n, MPI_DOUBLE_PRECISION, mMintercomm, packsizeinc, fesmpierr)
  if (fesmpisendbufpos+packsizeinc>fesmpisendbufsize) then
    call fesmpireallocsendbuf(fesmpisendbufpos+packsizeinc)
  endif
  call MPI_Pack(buf, n, MPI_DOUBLE_PRECISION, fesmpisendbuf, fesmpisendbufsize, fesmpisendbufpos, mMintercomm, fesmpierr)
end subroutine fesmpipack_doublearray

!########################################################################################################

!########################################################################################################
!> Send a pack (type PACKED) trough MPI communicator mMintercomm
!########################################################################################################
!> \ingroup DCM
subroutine fesmpipacksend
  implicit none
  call MPI_Send(fesmpisendbuf, fesmpisendbufpos, MPI_PACKED, 0, 0, mMintercomm, fesmpierr)
  fesmpisendbufpos=0
end subroutine fesmpipacksend

!########################################################################################################

!########################################################################################################
!> The buffer is reallocated before sending the message through MPI communicator because the buffer is 
!> smaller than the message to be sent through the  MPI communicator
!########################################################################################################
!> \ingroup DCM
subroutine fesmpireallocsendbuf(newbufsize)
  implicit none
  integer :: newbufsize, ierr
  character, pointer :: tmpbuf(:)

  newbufsize=2*newbufsize

  if (associated(fesmpisendbuf)) then
    allocate(tmpbuf(newbufsize), stat=ierr)
    if(ierr/=0) stop 'FATAL ERROR: problem while allocating fesmpi receive buffer'
    tmpbuf(1:fesmpisendbufsize)=fesmpisendbuf(1:fesmpisendbufsize)
    deallocate(fesmpisendbuf)
    fesmpisendbuf=>tmpbuf
  else
    allocate(fesmpisendbuf(newbufsize))
  endif
  fesmpisendbufsize=newbufsize
end subroutine fesmpireallocsendbuf

!########################################################################################################

!===================================================================================================
!                         SUBROUTINE OPEN_MM
!
!> Subroutine to open the communication with microMegas 
!!
!! 1) The number of processors that will be used by the DD code, and the path to the executable
!!    are passed as input in subroutine 
!! 3) The executable of the DD code is run by using MPI_SPAWN 
!! 4) The current working directory location is sent to the DD code
!! 5) A test number (dislo_ack) is received from amitex to test if the communication has 
!!    been successfull
!!
!! ICI:
!!      DD_num_procs    :      number of processors used by the DD code
!!      DD_name         :      location of the executable of the DD code
!!      cwd             :      current working directory path
!!      dislo_ack       :      number received from the DD code to test if the communication
!!                             between amitex and the DD code has been successfull  
!!
!!
!===================================================================================================
subroutine open_mM(DD_num_procs,DD_name)

    implicit none
    integer,intent(in)              :: DD_num_procs
    character(len=200),intent(in)   :: DD_name
    integer              :: dislo_ack=-1,pos=0,size_cwd,my_size
    integer              :: tag=100,testcomm=0
    character(len=1024)  :: cwd
    character(len=100)   :: buff
    integer, dimension(:), allocatable:: errCodes
    
    call getcwd(cwd)

    if (kkdebug>0) print *,"The current working directory is : ",cwd
    size_cwd = LEN_TRIM(cwd)
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD,my_size,ierror)
    allocate( errCodes(my_size))
    

    if (nrank==0) then
       !*** open_micromegas ***
       call MPI_COMM_SPAWN(trim(DD_name),MPI_ARGV_NULL,DD_num_procs,MPI_INFO_NULL,0, &
             MPI_COMM_SELF,mMintercomm,errCodes, ierror)
       if(ierror/=0) print * , "AmiteX could not spawn mM"

       if ( kkdebug > 0 )  print *,"spawn to MicroMegas. DD num procs",DD_num_procs," My rank=",nrank    

       call MPI_Pack(testcomm,1,MPI_INTEGER,buff,100,pos,mMintercomm,ierror)           ! pack 0 to test communications
       call MPI_Pack(size_cwd,1, MPI_INTEGER,buff,100,pos,mMintercomm,ierror)          ! Size of the string containing the current working directory
       call MPI_Pack(trim(cwd),size_cwd,MPI_CHARACTER,buff,100,pos,mMintercomm,ierror) ! Current working directory
    
       call MPI_Send(buff,100,MPI_PACKED,0,tag,mMintercomm,ierror) 
    
       call MPI_RECV(dislo_ack,1,MPI_INTEGER,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror)
       if( dislo_ack/=666)  print *,"MicroMegas returned the wrong magic number"    
    end if

end subroutine open_mM    
!===================================================================================================
  
  
!===================================================================================================
!                         SUBROUTINE SEND_MATRIX_OF_ELASTICITY_ISO
!
!> Subroutine to send elastic constants to the DD code microMegas 
!!
!! 1) Build matrix of elasticity from isotropic elastic constants lambda and mu
!! 2) Send matrix of elasticity to the DD code microMegas
!! 3) A test number (dislo_ack) is received from amitex to test if the communication has 
!!    been successfull
!!
!! ICI:
!!
!!      lambda,mu       :      isotropic elastic constants
!!      K(6,6)          :      full matrix of elastic constants built from the isotropic 
!!                             elastic constants lambda and mu
!!      dislo_ack       :      number received from the DD code to test if the communication
!!                             between amitex and the DD code has been successfull
!!
!===================================================================================================  
  subroutine send_matrix_of_elasticity(lambda,mu)

  implicit none

  real (kind=mytype),intent(in)     :: lambda,mu !lambda and mu are in Pa
  real (kind=mytype)                :: K(6,6)
  integer                           :: i,j
  character(len=100)                :: buff
  integer                           :: pos=0,dislo_ack=-1,tag=100
  print *,"lambda, mu: ", lambda,mu
  !! Assemble elasticity matrix from Lamé constants in MPa.
  !! This is done to maintain the interface with mM as it is right now
  do i=1,3
    K(i,i) = lambda*1d-6+2*mu*1.d-6
    K(i+3,i+3) = 2*mu*1.d-6
  end do
  do i=1,2
  K(1,1+i)=lambda*1.d-6
  K(i+1,1)=lambda*1.d-6   
  end do
  K(2,3)=lambda*1.d-6
  K(3,2)=lambda*1.d-6

  if (nrank==0) then
    do i=1,6
      do j=1,6
             call MPI_Pack(K(i,j),1, MPI_DOUBLE_PRECISION,buff,100,pos,mMintercomm,ierror)  
      end do
      call MPI_Send(buff,100,MPI_PACKED,0,tag,mMintercomm,ierror)
      pos = 0;
    end do
  
    call MPI_Recv(dislo_ack,1,MPI_INTEGER,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror)
    if( dislo_ack/=333)  print *,"MicroMegas returned the wrong magic number"
    
  end if
end subroutine send_matrix_of_elasticity
!===================================================================================================
 
!===================================================================================================
!                         SUBROUTINE TO SEND RESTART KEY 
!
!> Subroutine to send a restart key to the DD code microMegas 
!! N.B. : The restart_script should be used to restart the calculation with amitex
!!        Presently: (1) the geometry file should be rewritten assigning 1 zone per voxel,
!!        (2) the internal variables should be read as input from the res*.bin files 
!!        (3) the loading file should be adjusted to assure the continuity of the loading, and 
!!        (4) the restart_key should be set equal to 1 in the user algorithm parameter file  
!!        in order to restart the calculation 
!!
!! 1) Send key_restart to the DD code microMegas
!! 2) A test number (dislo_ack) is received from amitex to test if the
!!    communication has  been successfull
!!
!! ICI:
!!
!!      key_restart     :      integer value : 0=no restart, 1=restart
!!
!!      dislo_ack       :      number received from the DD code to test if the communication
!!                             between amitex and the DD code has been successfull
!!
!===================================================================================================  
  subroutine send_restart_key(restart_key)

  implicit none
  integer,optional,intent(in) :: restart_key
  integer                     :: tag=101
  integer                     :: dislo_ack,restart_key_loc

   if (.NOT. present (restart_key)) restart_key_loc=0
   if (      present (restart_key)) restart_key_loc=restart_key

   if (nrank==0) then 

     !-MPI- Begin Z ->DD No. 4
     call MPI_Send(restart_key_loc,1,MPI_INTEGER,0,tag,mMintercomm,ierror) 

     call MPI_Recv(dislo_ack,1,MPI_INTEGER,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror)
     if( dislo_ack/=444)  print *,"MicroMegas returned the wrong magic number"
    
   end if

  end subroutine send_restart_key
!===================================================================================================  


 
!===================================================================================================
!                         SUBROUTINE VERIFICATION_TIMESTEP
!
!> Subroutine to send time step and time step ratio to the DD code microMegas 
!!
!! 1) Send Amitex time step to the DD code microMegas 
!! 2) Send ratio between Amitex time step and DD time step to the DD code microMegas
!! 3) A test number (dislo_ack) is received from amitex to test if the communication has 
!!    been successfull
!!
!! ICI:

!!      dt              :      Amitex time step
!!      dt_ratio        :      ratio between microMegas et amitex time step
!!                             dt_ratio>=1 ;if dt_ratio>1 microMegas call amitex to update the 
!!                             stress field every n steps,where n=dt_ratio   
!!      dislo_ack       :      number received from the DD code to test if the communication
!!                             between amitex and the DD code has been successfull
!!
!===================================================================================================   
subroutine send_timestep(dt,dt_ratio)
  
  implicit none
  
  integer,intent(in)            ::  dt_ratio
  real(mytype),intent(in)       ::  dt
  integer                       ::  pos=0,dislo_ack,tag=100
  character(len=100)            ::  buff

  if (nrank==0) then 
    call MPI_Pack(dt_ratio,1, MPI_INTEGER,buff,100,pos,mMintercomm,ierror)  
    call MPI_Send(buff,pos,MPI_PACKED,0,tag,mMintercomm,ierror)
    if ( kkdebug > 0 ) print *,"AMITEX--------> send time step ratio",dt_ratio
    pos = 0
    call MPI_Pack(dt/real(dt_ratio,mytype),1,MPI_DOUBLE_PRECISION,buff,100,pos,mMintercomm,ierror)
    call MPI_Send(buff,pos,MPI_PACKED,0,100,mMintercomm,ierror)
    if ( kkdebug > 0 )  print *,"AMITEX--------> send time step ",dt/real(dt_ratio,mytype)
    
    dislo_ack=-1 
  
    call MPI_Recv(dislo_ack,1,MPI_INTEGER,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror)
    if(dislo_ack/=777) then
      print *,"Transferring deltatef to INCREMENT_DE_TEMPS"
    end if
    if ( kkdebug > 0 ) print *,"AMITEX--------> receive dislo_ack",dislo_ack
 
 end if
  

end subroutine send_timestep  
!===================================================================================================
  
!===================================================================================================
!                         SUBROUTINE GET_GLOBAL_DATA
!
!> Subroutine to receive some key parameters from the DD code microMegas 
!!
!! 1) The proc 0 receive the DD box size 
!! 2) The proc 0 test if the ratio between amitex box size and the DD code box size is the same
!!    along x, y and z direction. An error is sent and amitex is stopped if it's not the case
!! 3) The box size values and the ratio value are sent to all the procs
!! 4) The proc 0 receive the Burger's vector module, avalue and the parameter h that is used for
!!    stress regularization and eigenstrain distribution   
!! 5) These three parameters are sent to all the procs
!!
!! ICI:
!!
!!      box             :      DD box size
!!      b               :      Burger's vector in the units of the DD code micromegas 
!!      h               :      parameter that defines the thickness of the eigenstrain distribution 
!!                             volume and which is used in the stress regularization expression,
!!                             expressed in in the units of the DD code micromegas
!!      NOTE : b and h are definedi in the DD code units !!! 
!!             b is kept in DD code units because the coordinates of the amitex grid are
!!             rescaled by "ratio", i.e. the ratio between the amitex grid and the DD grid
!!             h is rescale by a factor ratio to spread correctly the eigenstrains on the 
!!             Amitex grid
!===================================================================================================   
subroutine get_global_data(ratio)
  
  implicit none
  
  real(mytype),intent(out)          :: ratio
  integer               :: pos=0
  character(len=1000)   :: buff
  character(len=100)    :: msg
  real(mytype)          :: r1,r2,r3,avalue

  if (nrank==0) then
  
    call MPI_Recv(buff,1000, MPI_PACKED, 0, MPI_ANY_TAG, mMintercomm, MPI_STATUS_IGNORE,ierror)
    call MPI_Unpack(buff,1000,pos,box,3,MPI_double_PRECISION,mMintercomm,ierror)

    r1=box(1)/grid%nx
    r2=box(2)/grid%ny
    r3=box(3)/grid%nz
    if ( abs(r1-r2)>1.e-3 .or. abs(r1-r3)>1.e-3 ) then
      msg="ratio between Amitex box size and DD box size is not the same along x,y and z"
      call amitex_abort( msg, 1, nrank)
    else
      ratio=box(1)/grid%nx
      print *," Ratio between mM grid and Amitex grid is : ",ratio
      print *," nx, ny, nz from Amitex: ",grid%nx,grid%ny,grid%nz
    end if 

  end if 
  call MPI_Bcast(ratio,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror);
  
  call MPI_Bcast(box,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror);
  
  if (nrank==0) then
    pos=0
    call MPI_Recv(buff,1000, MPI_PACKED, 0, MPI_ANY_TAG, mMintercomm, MPI_STATUS_IGNORE,ierror)
    call MPI_Unpack(buff,1000,pos,b,1,MPI_double_PRECISION,mMintercomm,ierror)
    call MPI_Unpack(buff,1000,pos,h,1,MPI_double_PRECISION,mMintercomm,ierror)
    call MPI_Unpack(buff,1000,pos,avalue,1,MPI_double_PRECISION,mMintercomm,ierror)
    b=b/avalue
    !call MPI_Recv(b,1,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror);
    !call MPI_Recv(h,1,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror);
    print *,">>>>>>>>>>>>>>> Thickness h : ",h
    !h=h*ratio
  endif
  
  call MPI_Bcast(b,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror);
  call MPI_Bcast(h,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror);

  if ( kkdebug > 0 ) print *," h ",h," b ",b," ratio ",ratio

end subroutine get_global_data  
!===================================================================================================

!===================================================================================================
!                         SUBROUTINE SEND_GLOBAL_DATA
!
!> Subroutine to send some key parameters to the DD code microMegas 
!> Here the volume box is send to the DD code mM
!!
!! 1) The proc 0 send the box volume  
!===================================================================================================
subroutine send_global_data(ratio)

  implicit none

  real(mytype),intent(in) :: ratio
  real(mytype)          :: mdcvol
  integer               :: pos=0,tag=100
  character(len=1000)   :: buff

  mdcvol=grid%nx*grid%ny*grid%nz*ratio*ratio*ratio
  if (nrank==0) then
      call MPI_Pack(mdcvol,1,MPI_DOUBLE_PRECISION,buff,1000,pos,mMintercomm,ierror)
      call MPI_Send(buff,pos,MPI_PACKED,0,tag,mMintercomm,ierror)
  endif

end subroutine send_global_data




!===================================================================================================
!                         SUBROUTINE EIGENSTRAIN_CALCULATION
!
!> Subroutine that receive the information about the areas swept by dislocation segments 
!> and calculate the corresponding eigenstrain by using a regularization procedure (Wei Cai et al. 2004)
!!
!! 1) The proc 0 receives the information about the areas swept by the dislocation segments 
!! 2) The proc 0 sends this information to all the procs
!! 3) Construct the table bxn (containing the cross product between the Burger's vector
!!    the vector normal to segment, for all the nbseg segments) and the vector dy (containing the norm of segment length 
!!    divided by the norm of the Burger's vector, for all the nbseg segments)
!! 4) The eigenstrain corresponding to the nbseg swept areas are calculated and save into ep   
!! 
!!
!! ICI:
!!
!!      nbseg           :      Number of dislocation segment that moved and induced a new 
!!                             deltaEigenstrain = Number of swept areas
!!      sweptsurfdata   :      Vector that contains the information about the swept areas 
!!                             - Dimansion = 16 * nbseg
!!                             - Each swept area is described by 16 parameters to reconstruct
!!                               its eigenstrain
!!                             - sweptsurfdata(1) = segment displacement amplitude perpendicular to the segment line
!!                             - sweptsurfdata(2:4) = segment origin coordinates 
!!                             - sweptsurfdata(5:7) = segment end coordinates
!!                             - sweptsurfdata(8:10) = displacement vector from the segment origin 
!!                             - sweptsurfdata(11:13) = displacement vector from the segment end
!!                             - sweptsurfdata(14:16) = Burger's vector 
!!      ep              :      Table to stock the eigeinstrain values produce by the nbseg 
!!      ratio           :      Ratio between amitex box size and the DD box size
!!      h               :      parameter that defines the thickness of the eigenstrain distribution 
!!                             volume and which is used in the stress regularization expression,
!!                             expressed in in the units of the DD code micromegas
!!      cutoff2         :
!!      bxn             :      vector containing the cross product between the Burger's vector and 
!!                             the vector normal to segment
!!      dy              :      norm of segment length divided by the norm of the Burger's vector
!! 
!===================================================================================================   

subroutine eigenstrain_calculation(ep,ratio,idt,istep)

    implicit none
    
    real(mytype),dimension(1:6,xsize(1)*xsize(2)*xsize(3)), &
    intent(out)                         :: ep 
    real(mytype),intent(in)             :: ratio
    integer,intent(in)                  :: idt,istep
    integer                             :: nbseg    
    integer,allocatable                 :: sweptsurfdata(:)
    real(mytype),allocatable            :: bxn(:),dy(:)
    real(mytype)                        :: h2,h4,cutoff2 
    real(mytype)                        :: S(6),buf(14),bxn_tmp(6),dy_tmp,bdS,x(3)
    type(Dist2psData)                   :: d2psdata
    integer                             :: absdep,ss(16),ix,iy,iz,xinf(3),xsup(3),iseg,i,j,k,indice_ep
    real(mytype)                        :: ttmp=0.
    

    if ( kkdebug ==istep ) print *,"Begin eigenstrain calculation on rank=",nrank

    if (key_cputime==1) then  
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      ttmp=MPI_WTIME()
    end if
    if(nrank==0) then 
      call fesmpipackrecv
      call MPI_Unpack(fesmpirecvbuf,fesmpirecvbufsize,fesmpirecvbufpos,nbseg,1,MPI_INTEGER,mMintercomm,fesmpierr)
    end if
    call MPI_Bcast(nbseg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror) 
      
    tcomm1=tcomm1+MPI_WTIME()-ttmp  
    
    
    if (nbseg>0) then
      allocate(sweptsurfdata(nbseg*16))
 
      if (key_cputime==1) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        ttmp=MPI_WTIME()
      end if

      if (nrank==0) call MPI_Unpack(fesmpirecvbuf,fesmpirecvbufsize,fesmpirecvbufpos,sweptsurfdata,16*nbseg, &
                         MPI_INTEGER,mMintercomm,fesmpierr)
      call MPI_Bcast(sweptsurfdata,16*nbseg,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

      tcomm1=tcomm1+MPI_WTIME()-ttmp
      
      if ( kkdebug == istep ) then
        if (nrank==0) then
          print *, "--------------> Number of areas received from mMegas = ", nbseg
          print *, "*************Received swept areas from mM"
          do iseg=0,nbseg-1
            print *,'amitex',iseg+1, sweptsurfdata(iseg*16+1:iseg*16+16)
          end do
        end if

      end if
 
      allocate(bxn(nbseg*6))
      allocate(dy(nbseg))
   
      bxn=0.
      dy=0.
    
      call setbxndy(sweptsurfdata,nbseg,bxn,dy)  !dy normalization (longueur du segment sur norme de b)

    end if 

    !! > eigenstrain computation
    h2=h*h;
    h4=h2*h2;
    cutoff2=(M_CUTOFF*M_CUTOFF)*h2;  

    if (idt==1) ep=0. 
    
    do iseg=0,nbseg-1
    
    ss(1:16)=sweptsurfdata(iseg*16+1:iseg*16+16)
    bxn_tmp(1:6)=bxn(iseg*6+1:iseg*6+6)
    dy_tmp=dy(iseg+1)
    absdep=ss(1) 
 
    ! boucle ix,iy,iz : gestion de la periodicite dans les 3 dimensions
    do ix=-1,1     
      do iy=-1,1   
        do iz=-1,1 
        
          ss(1)=sweptsurfdata(iseg*16+1)
          ss(2)=sweptsurfdata(iseg*16+2)+ix*int(box(1))
          ss(5)=sweptsurfdata(iseg*16+5)+ix*int(box(1))
          ss(3)=sweptsurfdata(iseg*16+3)+iy*int(box(2))
          ss(6)=sweptsurfdata(iseg*16+6)+iy*int(box(2))
          ss(4)=sweptsurfdata(iseg*16+4)+iz*int(box(3))
          ss(7)=sweptsurfdata(iseg*16+7)+iz*int(box(3))
    
            do while(ss(1)/=0) !(sweptsurfdata(iseg*16+1)/=0)
        
                call sweptsurfshift(ss(1:16),S,absdep) ! on discretise les "aires balayees"
        
                bdS = b * dy_tmp * sqrt( (S(4)-S(1))*(S(4)-S(1)) + (S(5)-S(2))*(S(5)-S(2)) + &
                                         (S(6)-S(3))*(S(6)-S(3)) )

                    if (ss(2)<=ss(5)) then
                      xinf(1)=nint((ss(2)-h)/ratio-0.5)
                      xsup(1)=nint((ss(5)+h)/ratio+0.5)
                    else
                      xinf(1)=nint((ss(5)-h)/ratio-0.5)
                      xsup(1)=nint((ss(2)+h)/ratio+0.5)
                    endif
                    if (xinf(1)<xstart(1)) xinf(1)=xstart(1)
                    if (xsup(1)>xend(1)) xsup(1)=xend(1)
    
                    if (ss(3)<=ss(6)) then
                      xinf(2)=nint((ss(3)-h)/ratio-0.5)
                      xsup(2)=nint((ss(6)+h)/ratio+0.5)
                    else
                      xinf(2)=nint((ss(6)-h)/ratio-0.5)
                      xsup(2)=nint((ss(3)+h)/ratio+0.5)
                    endif
                    if (xinf(2)<xstart(2)) xinf(2)=xstart(2)
                    if (xsup(2)>xend(2)) xsup(2)=xend(2)
    
                    if (ss(4)<=ss(7)) then
                      xinf(3)=nint((ss(4)-h)/ratio-0.5)
                      xsup(3)=(nint((ss(7)+h)/ratio+0.5))
                    else
                      xinf(3)=nint((ss(7)-h)/ratio-0.5)
                      xsup(3)=(nint((ss(4)+h)/ratio+0.5))
                    endif
                    if (xinf(3)<xstart(3)) xinf(3)=xstart(3)
                    if (xsup(3)>xend(3)) xsup(3)=xend(3)
 
                    do k=xinf(3),xsup(3)
                        
                      do j=xinf(2),xsup(2)
       
                        do i=xinf(1),xsup(1)
                        
                          indice_ep=(k-xstart(3))*xsize(1)*xsize(2)+(j-xstart(2))*xsize(1)+(i-xstart(1))+1
                          x(1)=(i-0.5)*ratio
                          x(2)=(j-0.5)*ratio
                          x(3)=(k-0.5)*ratio
                          
                          call dist2ps(x(1:3),S,buf,d2psdata)
                          call gammaplas(ep(1:6,indice_ep), &
                                             b, dy_tmp, bxn_tmp, d2psdata, h2, h4, cutoff2)
                  
                        end do

                      end do
 
                    end do
       
               
            end do
            
        end do
      end do 
    end do 
  
  end do
  
  if (nbseg>0) then
    deallocate(bxn)
    deallocate(dy)
    deallocate(sweptsurfdata)
  end if 

end subroutine eigenstrain_calculation
!===================================================================================================


!===================================================================================================
!                         SUBROUTINE SETBXNDY
!
!> Subroutine to calculate the cross product between the Burger's vector b and the vector normal to 
!> the dislocation segment n for all the dislocation segment that moved producing an eigenstrain : bxn
!> and the value of norm(n)/norm(b) for the same ensemble of segments : dy
!!
!! 1) calculation of bxn for the nbseg segments 
!! 2) calculation of dy for the nbseg segments
!!
!! ICI:
!!
!!      bxn       :      vector of  dimension 6*nbseg : 6 components of the product b x n
!!                       for the nbseg segments
!!      dy        :      vector of dimension nbseg : norm(n)/norm(b) for the nbseg segments
!!
!! Very Important:
!! Here Voigt notation is used (e1=exx, e2=eyy, e3=ezz, e4=2exy, e5=2exz, e6=2eyz).
!! Hence components 4, 5, and 6 of bxn are multiply by a factor 2 since they
!! contribute to off-diagonals components e4, e5, and e6!
!!
!===================================================================================================   
subroutine setbxndy(sweptsurfdata,nbseg,bxn,dy)
  implicit none
  
  integer,intent(in )                    :: sweptsurfdata(nbseg*16),nbseg
  real(mytype),intent(out)               :: bxn(nbseg*6),dy(nbseg)
  integer                                :: dep,a(3),b(3),c(3),iseg
  real(mytype)                           :: normn,normb,fac,l(3),n(3),absdep

  do iseg=0,nbseg-1
    ! > Calculation of the vector bxn

    dep = sweptsurfdata(16*iseg+1)
    absdep = abs(dep)

    a(1:3)=sweptsurfdata(16*iseg+2:16*iseg+4)
    b(1:3)=sweptsurfdata(16*iseg+5:16*iseg+7)
    c(1:3)=sweptsurfdata(16*iseg+8:16*iseg+10)

    l(1:3)=a(1:3)-b(1:3)

    n(1)=(c(2)*l(3)-c(3)*l(2))/absdep
    n(2)=(c(3)*l(1)-c(1)*l(3))/absdep
    n(3)=(c(1)*l(2)-c(2)*l(1))/absdep

    b(1:3)=sweptsurfdata(16*iseg+14:16*iseg+16) ! re-using vector b

    normn=sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
    normb=sqrt(1._mytype*(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)))
    fac=1./(normn*normb);

    bxn(iseg*6+1)=fac*b(1)*n(1)
    bxn(iseg*6+2)=fac*b(2)*n(2)
    bxn(iseg*6+3)=fac*b(3)*n(3)
    bxn(iseg*6+4)=fac*(b(1)*n(2)+b(2)*n(1))!*0.5
    bxn(iseg*6+5)=fac*(b(1)*n(3)+b(3)*n(1))!*0.5
    bxn(iseg*6+6)=fac*(b(2)*n(3)+b(3)*n(2))!*0.5
    ! < End calculation of the vectorbxn

    ! > Calculation of the vector dy
    normb=sqrt(l(1)*l(1)+l(2)*l(2)+l(3)*l(3)) ! re-using normb for the length of the segment

    dy(iseg+1)=normn/normb;
    ! < End of calculation of the vector dy
  end do


end subroutine setbxndy
!===================================================================================================

!===================================================================================================
!                         SUBROUTINE SWEPTSURFSHIFT
!
!> Subroutine to shift the position on the swept area. The Eigenstrain produced by a single 
!> segment is calculated as the sum of deltaeigenstrain each one corresponding to a discrete
!> steps on the DD grid of micromegas  
!> summation of eigenstrain
!!
!! 1) The proc 0 send the parameter endcalcul to the DD code 
!! 2) A test number (dislo_ack) is received from amitex to test if the communication has 
!!    been successfull
!!
!! ICI:
!!
!!      ss              :      portion of the vector sweptsurfdata corresponding to one 
!!                             segment i
!!      mseg            :      mseg(1:3) : middle position of the displacement vector starting at the 
!!                             origin of segment
!!                             mseg(4:6) : middle position of the displacement vector starting at the 
!!                             end of segment
!!
!===================================================================================================   

subroutine sweptsurfshift (ss,mseg,absdep)
    implicit none 
    integer,intent(inout)            :: ss(16)
    real(mytype),intent(out)         :: mseg(6)
    integer,intent(in)               :: absdep
    integer                          :: shift

    shift=sign(1,ss(1))     !shift is equal to 1 or -1 depending on the signof ss(1)    
    ss(1)=ss(1)-1*shift           !decrease by 1 remaining displacement amplitude
    
    !In the following abs(shift)=1 otherwise one's would need to multiply by abs(shift) the 
    !right term in the following 4 expression
    
    mseg(1:3)=ss(2:4)+ss(8:10)/absdep*0.5
    mseg(4:6)=ss(5:7)+ss(11:13)/absdep*0.5
    
    ss(2:4)=ss(2:4)+ss(8:10)/absdep
    ss(5:7)=ss(5:7)+ss(11:13)/absdep
    
end subroutine
!===================================================================================================

!===================================================================================================
!!                         SUBROUTINE DIST2PS
!!
!> Compute position of the point P with respect to the segment S!!
!!
!! ICI:
!!     P: input point (points to 3 double)
!!     S: input segment (points to 6 double [S0x, S0y, S0z, S1x, S1y, S1z])
!!     buf: buffer, for optimisation. MUST BE allocated at least 14 double
!!     d2: output data
!===================================================================================================
subroutine dist2ps(P,S,buf,d2)
  implicit none
  real(mytype),dimension(3),intent(in)   :: P
  real(mytype),dimension(6),intent(in)   :: S
  real(mytype),dimension(14),intent(out) :: buf
  type(Dist2psData),intent(out)          :: d2
  real(mytype)                           :: vtdtt
  
  buf(1:3)  = P(1:3)-(S(1:3)+S(4:6))*0.5 !buf(1:3) contains vector V=(center of S)->P
  buf(4:6) = S(4:6)-S(1:3) ! buf(4:6) contains T=S1-S0
  buf(7) =buf(4)*buf(4)+buf(5)*buf(5)+buf(6)*buf(6); ! T.T
  buf(8) =buf(1)*buf(4)+buf(2)*buf(5)+buf(3)*buf(6); ! V.T
  vtdtt=buf(8)/buf(7); ! (V.T)/(T.T)
  buf(9:11)=buf(4:6)*vtdtt ! buf(9:11) contains Vt=(V.T)T/(T.T)
  buf(12:14)=buf(1:3)-buf(9:11) ! buf(12:14) contains Vn=V-Vt
  d2%l2=buf(7);
  d2%dn2=buf(12)*buf(12)+buf(13)*buf(13)+buf(14)*buf(14);
  d2%dt2=buf(9)*buf(9)+buf(10)*buf(10)+buf(11)*buf(11);
  d2%sdt2=sign(1.D0, buf(8))
  
end subroutine dist2ps
!===================================================================================================

!===================================================================================================
!!                         SUBROUTINE GAMMAPLAS
!!
!> Compute plastic contribution to gamma of the segment iseg at the grid point x
!> The plastic strain at x ep_tmp is received as input, is incremented in this subroutine
!> and sent back as output
!!
!! ICI:
!!     ep_tmp       :   ep at the grid point x
!!     bxt_tmp      :   b x n vector relative to the segment iseg 
!!     dy_tmp       :   dy value relative to the segment iseg
!!     d2           :   information about the distance of the segment iseg from the
!!                      point x
!!     b            :   Burger's vector norm 
!!     h2           :   h * h, being h the regularization thickness
!!     h4           :   h2 * h2, being h the regularization thickness
!!     cutoff2      :   cutoff * cutoff, being cutoff ...
!===================================================================================================
subroutine gammaplas(ep_tmp,  b, dy_tmp, bxn_tmp, d2, h2, h4, cutoff2) 
  
  implicit none
  real (mytype),dimension(6),intent(inout) :: ep_tmp
  real (mytype),dimension(6),intent(in)     :: bxn_tmp
  real (mytype),intent(in)                  :: b,dy_tmp,h2,h4,cutoff2
  type(Dist2psData),intent(in)               :: d2 
  real (mytype)                             :: adt,dt,l,ll
  real (mytype)                             :: gplas

  
  adt=sqrt(d2%dt2);
  dt=d2%sdt2*adt;
  l=sqrt(d2%l2);
  ll=l/2.;

  ! check if point is inside a end disk : ((adt-ll)*(adt-ll)+d2%dn2<cutoff2);
  ! check if point is outside capsule   :
  if ( .not.( (d2%dn2<cutoff2 .and. adt<ll) .or.  ((adt-ll)*(adt-ll)+d2%dn2<cutoff2) ) ) then 
    return
  else
    gplas=b*dy_tmp*intsegwtilde(d2%dn2, dt, ll, h2, h4)
    ep_tmp(1:6)=ep_tmp(1:6)+gplas*bxn_tmp(1:6)
  end if

end subroutine gammaplas
!===================================================================================================


!===================================================================================================
!!                         FUNCTION POW2p5
!!
!> Compute x^(2.5) being x the input variable
!!   
!===================================================================================================
function pow2p5(x) !x^(5./2.)
  implicit none
  real(mytype) :: x,xx,pow2p5

  x=sqrt(x)
  xx=x*x
  pow2p5=xx*xx*x
end function pow2p5
!===================================================================================================

!===================================================================================================
!!                         FUNCTION POW3p5
!!
!> Compute x^(3.5) being x the input variable
!!   
!===================================================================================================
function pow3p5(x) !x^(7./2.)
  implicit none
  real(mytype) :: x,xx,pow3p5

  x=sqrt(x)
  xx=x*x
  pow3p5=xx*xx*xx*x
end function pow3p5
!===================================================================================================

!===================================================================================================
!!                         FUNCTION w,wtilde
!!
!> w, wtilde: fonctions to compute the Wei Cai distribution
!> cf: A non-singular continuum theory of dislocations, Cai, Arsenlis, Weinberger, Bulatov, Journal of the Mechanics and Physics of Solids, 2006   
!!
!===================================================================================================

function w(r2,a2,a3) 

  implicit none 
  real(mytype)      ::r2,a2,a3,w
  
  w=15./(8.*Pi*a3*pow3p5((r2/a2)+1))

end function w
!===================================================================================================


!===================================================================================================
function wtilde(r2, a2, a3) 

  implicit none
  real(mytype)      ::r2,a2,a3,wtilde
  
  wtilde=.3425*w(r2, .81685444*a2, .738273042872*a3)+ &
         .6575*w(r2, 0.29713401*a2, 0.161967748851*a3)

end function wtilde
!===================================================================================================

!===================================================================================================
!!                         FUNCTION intsegw, intsegwtilde
!!
!> intsegw, intsegwtilde: fonctions to compute the integral of w and wtilde over a segment
!> cf: A non-singular continuum theory of dislocations, Cai, Arsenlis, Weinberger, Bulatov, Journal of the Mechanics and Physics of Solids, 2006   
!!
!===================================================================================================
function intsegw(dn2, p, p2, p4, q, q2, q4, a2, a4) ! p=ll-dt, q=-ll-dt
  implicit none
  real (mytype) :: dn2,p,p2,p4,q,q2,q4,a2,a4
  real (mytype) :: s,s2,s3
  real (mytype) :: intsegw
  
  s=dn2+a2
  s2=s*s
  s3=s2*s

  intsegw = a4/(8.*Pi)*(                                                                                &
            p*(15.*s2 + 20.*s*p2 + 8.*p4) / (s3*pow2p5(s+p2))        &
            -q*(15.*s2 + 20.*s*q2 + 8.*q4) / (s3*pow2p5(s+q2))        &
                )
           
end function intsegw
!===================================================================================================

!===================================================================================================
function intsegwtilde(dn2, dt,ll, a2, a4) ! dt=sdt2*sqrt(dt2), ll=l/2
  implicit none
  real (mytype)  :: dn2,dt,ll,a2,a4
  real (mytype)  :: p,p2,p4,q,q2,q4
  real (mytype)  :: intsegwtilde
  
  p=ll-dt
  p2=p*p
  p4=p2*p2
  
  q=-ll-dt
  q2=q*q
  q4=q2*q2
  
  intsegwtilde=0.3425*intsegw(dn2, p, p2, p4, q, q2, q4, 0.81685444*a2, 0.6672511761477136*a4)+ &
               0.6575*intsegw(dn2, p, p2, p4, q, q2, q4, 0.29713401*a2, 0.0882886198986801*a4)
           
end function intsegwtilde
!===================================================================================================


subroutine calculation_stress_avg(s0)
  implicit none
  real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),1:6), &
      intent(in)                :: s0
  integer                       :: isig,nsig,kk
  real(mytype),dimension(6)     :: SigAvg,SigAvg_loc


  SigAvg=0.
  SigAvg_loc=0.
  nsig=xsize(1)*xsize(2)*xsize(3)

  ! do k=xstart(3),xend(3)
  !   do j=xstart(2),xend(2)
  !     do i=xstart(1),xend(1)
  !       isig=(k-xstart(3))*xsize(1)*xsize(2)+(j-xstart(2))*xsize(1)+(i-xstart(1))+1
   do isig=1,nsig
         do kk=1,6
           SigAvg_loc(kk)=SigAvg_loc(kk)+s0(isig,kk)
         end do
  end do
  !     end do  
  !  end do
  !end do    

  call MPI_Reduce(SigAvg_loc,SigAvg,6,MPI_DOUBLE_PRECISION,MPI_SUM, &
         0,MPI_Comm_world,ierror)
  SigAvg=SigAvg/real(grid%nx*grid%ny*grid%nz,mytype)
  
  if (nrank==0) then
  
    call fesmpipack_doublearray(SigAvg,6)
    call fesmpipacksend

    !write (40,*) 40, SigAvg
    
  end if


end subroutine calculation_stress_avg




!===================================================================================================
!                         SUBROUTINE CALCULATION_STRESS_SEG
!
!> Subroutine to interpolate the stress values at the dislocation segments position.
!!
!! 1) The proc 0 receive the positions of the Integration Points (IP), i.e. the points where the stress 
!!    is calculated (in mM the center of the dislocation segments)
!! 2) The IP positions are sent to all the procs
!! 2) The stress field is calculated at the IP by using a deconvolution on the voxel belonging to the 
!!    eigenstrain distribution volume
!! 3) The proc 0 send the stress tensor values at each IP to the dd code mM
!!
!! ICI:
!!
!!      ratio       :      ratio between the amitex grid and the dd grid
!!      s0          :      vector containing the stress tensor on the amitex grid (Pa)
!!      istep       :      current step in the main loop (use in the debug procedure) 
!!      nbseg       :      number of integration points (in micromegas nbseg=number of dislocation segments 
!!      seg_pos     :      vector containing the position of the IP (dislocation segment centers in micromegas)
!!      seg_str     :      vector containing the total stress tensor values at the IP (Pa converted into MPa
!!                         at the end of the subroutine)
!!      seg_str_loc :      vector containing the stress tensor values at the IP on the current proc:
!!                         only the segment positions inside the domain accessible to the current proc 
!!                         are considered 
!!      meanfactor  :      vector containing the weigth of the w function for all the IP 
!!      meanfactor_loc  :  vector containing the weigth of the w function at the IP on the current proc
!!                         only the segment positions inside the domain accessible to the current proc 
!!                         are considered 
!!
!===================================================================================================   
subroutine calculation_stress_seg(ratio,s0,istep)

    implicit none
    
    real (mytype),parameter     :: halfsqrt2 = 0.707106781185d0
    real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),1:6), &
      intent(in)                :: s0
    real(mytype),intent(in)     :: ratio
    integer,intent(in)          :: istep
    integer                     :: nbseg
    integer,allocatable         :: seg_pos(:)
    real(mytype),allocatable    :: seg_str(:),seg_str_loc(:),meanfactor(:),meanfactor_loc(:)
    integer                     :: ix,iy,iz,xinf(3),xsup(3),iseg,i,j,k,rp(3),kk,indice_sig
    real(mytype)                :: h2,h3,cutoff,cutoff2,w,x(3),r2,ttmp=0.
    
 
  

    !**********************************
    !receive coordinates of dislocation segment centers from mM 

    if (key_cputime==1) then    
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      ttmp=MPI_WTIME()
    end if
    
    if(nrank==0) then 
      call fesmpipackrecv
      call MPI_Unpack(fesmpirecvbuf,fesmpirecvbufsize,fesmpirecvbufpos,nbseg,1,MPI_INTEGER,mMintercomm,fesmpierr)
      if ( kkdebug == istep ) print *, "The number of segment centers received from mM is :",nbseg 
    end if
    call MPI_Bcast(nbseg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

    tcomm2=tcomm2+MPI_WTIME()-ttmp


    if (nbseg>0) then
      allocate(seg_pos(3*nbseg)) 

      if (key_cputime==1) then    
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        ttmp=MPI_WTIME()
      end if
      
      if (nrank==0) call MPI_Unpack(fesmpirecvbuf,fesmpirecvbufsize,fesmpirecvbufpos,seg_pos,3*nbseg, &
                         MPI_INTEGER,mMintercomm,fesmpierr)
      call MPI_Bcast(seg_pos,3*nbseg,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)  

      tcomm2=tcomm2+MPI_WTIME()-ttmp
      
      allocate(seg_str(6*nbseg))
      allocate(seg_str_loc(6*nbseg))
      allocate(meanfactor(nbseg))
      allocate(meanfactor_loc(nbseg))
      seg_str=0.
      seg_str_loc=0.
      meanfactor=0.
      meanfactor_loc=0.
    end if

    if ( kkdebug == istep ) then

      if (nrank==0) then
        print *,  "AMITEX---->received nbseg=",nbseg
        do i=0,nbseg-1
          print *,"AMITEX---->",i+1,seg_pos(i*3+1:i*3+3)
        end do
      end if

    end if

    !**********************************
    !Interpolation of the stress field at the IP, i.e. at dislocation segment centers in mM

    if ( kkdebug == istep )  print *, "Begin Stress interpolation at segment center position on proc with rank",nrank

    h2=h*h
    h3=h*h2
    cutoff=M_cutoff*h
    cutoff2=cutoff*cutoff
    
    do iseg=0,nbseg-1
      do ix=-1,1     
        do iy=-1,1   
          do iz=-1,1 
            
            rp(1)=seg_pos(iseg*3+1)+ix*int(box(1))
            rp(2)=seg_pos(iseg*3+2)+iy*int(box(2))
            rp(3)=seg_pos(iseg*3+3)+iz*int(box(3))
            !determination de la boite pour application du filtre
            do kk=1,3
              xinf(kk)=nint((rp(kk)-cutoff)/ratio-0.5)
              xsup(kk)=nint((rp(kk)+cutoff)/ratio-0.5)
              if (xinf(kk)<xstart(kk)) xinf(kk)=xstart(kk)
              if (xsup(kk)>xend(kk)) xsup(kk)=xend(kk)
            end do
            do k=xinf(3),xsup(3)
              do j=xinf(2),xsup(2)
                 do i=xinf(1),xsup(1)
                    indice_sig=(k-xstart(3))*xsize(1)*xsize(2)+(j-xstart(2))*xsize(1)+(i-xstart(1))+1
                    x(1)=(i-0.5)*ratio
                    x(2)=(j-0.5)*ratio
                    x(3)=(k-0.5)*ratio
                    r2=0.
                    do kk=1,3 
                     r2=r2+(x(kk)-rp(kk))*(x(kk)-rp(kk))
                    end do
                    if (r2<cutoff2) then
                      w=wtilde(r2,h2,h3)
                      do kk=1,4
                        seg_str_loc(iseg*6+kk)=seg_str_loc(iseg*6+kk)+s0(indice_sig,kk)*w     
                      end do
                      seg_str_loc(iseg*6+5)=seg_str_loc(iseg*6+5)+s0(indice_sig,6)*w     
                      seg_str_loc(iseg*6+6)=seg_str_loc(iseg*6+6)+s0(indice_sig,5)*w     
                      meanfactor_loc(iseg+1)=meanfactor_loc(iseg+1)+w             
                    end if 
      
                 end do
              end do
            end do
            
          end do
        end do
      end do
      
        
    end do
    
    !MPI_Reduce is used because only the proc 0  send the the stress values to mM    
    call MPI_Reduce(seg_str_loc,seg_str,nbseg*6,MPI_DOUBLE_PRECISION,MPI_SUM, &
         0,MPI_Comm_world,ierror)
    call MPI_Reduce(meanfactor_loc,meanfactor,nbseg,MPI_DOUBLE_PRECISION,MPI_SUM, &
         0,MPI_Comm_world,ierror)  

    !**********************************
    !Send the stress field at dislocation segment centers

    if ( kkdebug == istep ) then

      print *, "End of Stress calculation at segment center position on proc with rank",nrank
      if (nrank==0) then
        print *,  "AMITEX---->send 6*",nbseg,"=",6*nbseg,"stress values"
        do i=0,nbseg-1
          print *,"AMITEX---->",i+1,seg_str(i*6+1:i*6+6)
        end do
      end if

    end if
 
    if (key_cputime==1) then
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      ttmp=MPI_WTIME()
    end if
    
    ! N.B. the off-diagonals components of the stress tensor are divided by halfsqrt2 to 
    ! preserve the compatibility with the dd code micromegas 
    ! For the same reason all the stress tensor components are transformed from Pa (amitex)
    ! to MPa (micromegas) before sending the stress values to mM
    if (nrank==0) then
      do iseg=0,nbseg-1
        seg_str(iseg*6+1:iseg*6+3)=seg_str(iseg*6+1:iseg*6+3)*1e-6/meanfactor(iseg+1) 
        seg_str(iseg*6+4:iseg*6+6)=seg_str(iseg*6+4:iseg*6+6)*1e-6/halfsqrt2/meanfactor(iseg+1) 
      end do  
      call fesmpipack_doublearray(seg_str,6*nbseg)
      call fesmpipacksend

    end if
    
    tcomm2=tcomm2+MPI_WTIME()-ttmp


    if ( kkdebug == istep ) then
      if (nrank==0) then
        print *, "Sent stress at segment center postion to mM "
        write (60,*) nbseg
        do i=0,nbseg-1
          write (60,*)  seg_str(i*6+1:i*6+6)
        end do
      end if
    end if

    if (nbseg>0) then
      deallocate(seg_pos)
      deallocate(seg_str)
      deallocate(seg_str_loc)
      deallocate(meanfactor_loc)
      deallocate(meanfactor)
    end if 
end subroutine calculation_stress_seg
!===================================================================================================
  


!===================================================================================================
!                         SUBROUTINE TELL_DD_WHAT_TO_DO
!
!> Subroutine to tell the DD code microMegas to continue the time loop or to exit the time loop
!> The final time is set in the amitex code but not in the DD code, is amitex that tells the 
!> DD code to continue the simulation or to exit from the time loop
!!
!! 1) The proc 0 send the parameter endcalcul to the DD code 
!! 2) A test number (dislo_ack) is received from amitex to test if the communication has 
!!    been successfull
!!
!! ICI:
!!
!!      endcalcul       :      variable send to the DD code : 0 tells the DD code to 
!!                             continue, 1 tells the DD code to exit the time loop
!!      dislo_ack       :      number received from the DD code to test if the communication
!!                             between amitex and the DD code has been successfull
!!
!===================================================================================================   
subroutine tell_DD_whattodo(endcalcul) 

  implicit none
  
  integer,intent(in)  :: endcalcul 
  integer             :: tag=104,dislo_ack=-1

  !! Send simulation status to mM 
  if (nrank==0) then
  
    call MPI_Send(endcalcul,1,MPI_INTEGER,0,tag,mMintercomm,ierror)  

    call MPI_Recv(dislo_ack,1,MPI_INTEGER,0,MPI_ANY_TAG,mMintercomm,MPI_STATUS_IGNORE,ierror)
    if ( dislo_ack/=999) then 
      print *,"MicroMegas returned the wrong magic number to FIN_CALCUL"
   ! else
   !   print *, "Received the correct Magic number from whattodo" 
    endif

  endif
 
end subroutine tell_DD_whattodo
!===================================================================================================



!===================================================================================================
!                         SUBROUTINE EIGENSTRAIN_CALCULATION_FROM_FILE
!
!> Subroutine that read the information about the areas swept by dislocation segments from an input file
!> and calculate the corresponding eigenstrain by using a regularization procedure (Wei Cai et al. 2004)
!!
!! 1) The proc 0 reads the information about the areas swept by the dislocation segments from the file "mm.dat"
!! 2) The proc 0 sends this information to all the procs
!! 3) Construct the table bxn (containing the cross product between the Burger's vector
!!    the vector normal to segment, for all the nbseg segments) and the vector dy (containing the norm of segment length 
!!    divided by the norm of the Burger's vector, for all the nbseg segments)
!! 4) The eigenstrain corresponding to the nbseg swept areas are calculated and save into ep   
!! 
!!
!! ICI:
!!
!!      nbseg           :      Number of dislocation segment that moved and induced a new 
!!                             deltaEigenstrain = Number of swept areas
!!      sweptsurfdata   :      Vector that contains the information about the swept areas 
!!                             - Dimansion = 16 * nbseg
!!                             - Each swept area is described by 16 parameters to reconstruct
!!                               its eigenstrain
!!                             - sweptsurfdata(1) = segment displacement amplitude perpendicular to the segment line
!!                             - sweptsurfdata(2:4) = segment origin coordinates 
!!                             - sweptsurfdata(5:7) = segment end coordinates
!!                             - sweptsurfdata(8:10) = displacement vector from the segment origin 
!!                             - sweptsurfdata(11:13) = displacement vector from the segment end
!!                             - sweptsurfdata(14:16) = Burger's vector 
!!      ep1             :      Table to stock the eigeinstrain values produce by the nbseg 
!!      ratio           :      Ratio between amitex box size and the DD box size
!!      h               :      parameter that defines the thickness of the eigenstrain distribution 
!!                             volume and which is used in the stress regularization expression,
!!                             expressed in in the units of the DD code micromegas
!!      cutoff2         :
!!      bxn             :      vector containing the cross product between the Burger's vector and 
!!                             the vector normal to segment
!!      dy              :      norm of segment length divided by the norm of the Burger's vector
!! 
!!      N.B.  This subroutine has been written at the beginning of the process of coupling the dd
!!            code mM with the fft solver amitex. It is currently not used. 
!!            It should  work, but it is better to check it before using it !!!
!!
!!
!===================================================================================================   


subroutine EigS_fromfile(ep1,ratio)
!TO DO check maxshift, and positions : where have positions been rescaled? 

  implicit none
  integer                               :: indice_ep1,ix,iy,iz
  real(mytype)                          :: ratio!,h,b
  real(mytype),dimension(1:6,xsize(1)*xsize(2)*xsize(3))  :: ep1
  real(mytype),allocatable              :: bxn(:),dy(:)
  integer,allocatable                   :: sweptsurfdata(:)!,idx(:)
  integer                               :: nbseg
  real(mytype)                          :: h2,h4,cutoff2

  real(mytype)                          :: S(6),buf(14),bxn_tmp(6),dy_tmp,bdS
  type(Dist2psData)                     :: d2psdata
  integer                               :: absdep,ss(16),xinf(3),xsup(3),iseg,i,j,k!,box
  real(mytype)                          :: x(3),r1,r2,r3
  character(len=100)                    :: msg
  
  if (nrank==0) then
      open(51,file="../in/mm.dat",status='old',position='rewind')
  
      read(51,*) box(1),box(2),box(3),h,b,flagfreesurf
      read(51,*) nbseg
      ratio=100
      if ( kkdebug > 0 )  print *," h=",h," b=",b," ratio=",ratio,"nbseg",nbseg
      allocate(sweptsurfdata(nbseg*16))
      
      do i=0,nbseg-1
        read(51,*)sweptsurfdata(i*16+1:i*16+16)
      end do 
  
      close(51)
      
      r1=box(1)/grid%nx
      r2=box(2)/grid%ny
      r3=box(3)/grid%nz
      if ( abs(r1-r2)>1.e-3 .or. abs(r1-r3)>1.e-3 ) then
        msg="ratio between Amitex box size and DD box size is not the same along x,y and z"
      call amitex_abort( msg, 1, nrank)
      else
        ratio=box(1)/grid%nx
        print *," Ratio between mM grid and Amitex grid is : ",ratio
        print *," nx, ny, nz from Amitex: ",grid%nx,grid%ny,grid%nz
    end if 
      
  end if
  
  call MPI_Bcast(nbseg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror) 
  call MPI_Bcast(sweptsurfdata,16*nbseg,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  
  allocate(bxn(nbseg*6))
  allocate(dy(nbseg))
  
  bxn=0._mytype
  dy=0._mytype
  
  
  !calculate bx and dy for all the swept areas (in FEniCs_MDC/src/mio.Mdc function regularizer::sweptarea::setbxndy)
  call setbxndy(sweptsurfdata,nbseg,bxn,dy)
  !end calculation bx and dy 
  
  !eigenstrain computation
  h2=h*h
  h4=h2*h2;
  cutoff2=(M_CUTOFF*M_CUTOFF)*h2;
  ep1=0._mytype
  
  if ( kkdebug > 0 )  print *,"xsize(1)",xsize(1),"xstart(1)",xstart(1),"xend(1)",xend(1)
  do iseg=0,nbseg-1
           
              
    ss(1:16)=sweptsurfdata(iseg*16+1:iseg*16+16)
    bxn_tmp(1:6)=bxn(iseg*6+1:iseg*6+6)
    dy_tmp=dy(iseg+1)
    absdep=ss(1) !sweptsurfdata(iseg*16+1)
    
    do ix=-1,1     
      do iy=-1,1   
        do iz=-1,1 
          
          ss(1)=sweptsurfdata(iseg*16+1)
          ss(2)=sweptsurfdata(iseg*16+2)+ix*int(box(1))
          ss(5)=sweptsurfdata(iseg*16+5)+ix*int(box(1))
          ss(3)=sweptsurfdata(iseg*16+3)+iy*int(box(2))
          ss(6)=sweptsurfdata(iseg*16+6)+iy*int(box(2))
          ss(4)=sweptsurfdata(iseg*16+4)+iz*int(box(3))
          ss(7)=sweptsurfdata(iseg*16+7)+iz*int(box(3))
    
            do while(ss(1)/=0) !(sweptsurfdata(iseg*16+1)/=0)
        
                call sweptsurfshift(ss(1:16),S,absdep)
        
                bdS = b * dy_tmp * sqrt( (S(4)-S(1))*(S(4)-S(1)) + (S(5)-S(2))*(S(5)-S(2)) + &
                                         (S(6)-S(3))*(S(6)-S(3)) )
                    if (ss(2)<=ss(5)) then
                      xinf(1)=nint((ss(2)-h)/ratio-0.5)
                      xsup(1)=nint((ss(5)+h)/ratio+0.5)
                    else
                      xinf(1)=nint((ss(5)-h)/ratio-0.5)
                      xsup(1)=nint((ss(2)+h)/ratio+0.5)
                    endif
                    if (xinf(1)<xstart(1)) xinf(1)=xstart(1)
                    if (xsup(1)>xend(1)) xsup(1)=xend(1)
    
                    if (ss(3)<=ss(6)) then
                      xinf(2)=nint((ss(3)-h)/ratio-0.5)
                      xsup(2)=nint((ss(6)+h)/ratio+0.5)
                    else
                      xinf(2)=nint((ss(6)-h)/ratio-0.5)
                      xsup(2)=nint((ss(3)+h)/ratio+0.5)
                    endif
                    if (xinf(2)<xstart(2)) xinf(2)=xstart(2)
                    if (xsup(2)>xend(2)) xsup(2)=xend(2)
    
                    if (ss(4)<=ss(7)) then
                      xinf(3)=nint((ss(4)-h)/ratio-0.5)
                      xsup(3)=(nint((ss(7)+h)/ratio+0.5))
                    else
                      xinf(3)=nint((ss(7)-h)/ratio-0.5)
                      xsup(3)=(nint((ss(4)+h)/ratio+0.5))
                    endif
                    if (xinf(3)<xstart(3)) xinf(3)=xstart(3)
                    if (xsup(3)>xend(3)) xsup(3)=xend(3)
                    
                    do k=xinf(3),xsup(3)
                      
                      do j=xinf(2),xsup(2)
       
                        do i=xinf(1),xsup(1)
                          indice_ep1=(k-xstart(3))*xsize(1)*xsize(2)+(j-xstart(2))*xsize(1)+(i-xstart(1))+1
                          x(1)=(i-0.5)*ratio
                          x(2)=(j-0.5)*ratio
                          x(3)=(k-0.5)*ratio
                          call dist2ps(x(1:3),S,buf,d2psdata)
                          call gammaplas(ep1(1:6,indice_ep1), &
                                             b, dy_tmp, bxn_tmp, d2psdata, h2, h4, cutoff2)
                  
                        end do

                      end do
 
                    end do
       
               
            end do
            
        end do
      end do 
    end do 
                                

  end do
  
  deallocate(dy)
  deallocate(bxn)
  deallocate(sweptsurfdata)

end subroutine EigS_fromfile


end module dd_mod  
