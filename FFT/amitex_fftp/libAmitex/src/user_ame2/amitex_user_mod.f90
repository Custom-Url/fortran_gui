!===============================================================================
!
!       MODULE AMITEX_USER_MOD FOR AME (AMITEX MULTI-ECHELLES)
!>
!>   - Declaration and initialization of USER defined variables 
!>   - User-defined procedures uses in resolution_user_mod or standard_user_mod 
!>
!>   User-variables for AME
!!
!!         - SIMU2(1:Nsimuslave+1)
!!         - AME(1:Nsimuslave+1)
!!
!===============================================================================
module amitex_user_mod

  use ISO_FORTRAN_ENV

  use decomp_2d
  use mpi

  use amitex_mod
  use error_mod
  use simu_mod,  only : SIMU_AMITEX,SIMU,& ! type
                        init_simu_command_line,nullify_pointers_simu, associate_pointers_simu ! functions
  use param_algo_mod, only : user_param  ! pointer 
  use io2_amitex_mod, only : print_field_vtk ! function
  use field_mod, only : gene_cartesian_coord

  implicit none

  private

  public :: init_user_variables

!!------------------------------------------------------------------------------
!>                                                              PUBLIC VARIABLES 
 
  public :: SIMU2, AME  
  

!!------------------------------------------------------------------------------
!>                                                               TYPE MULTISCALE
!>
!>                                       Any 'current' simulation (I) (SIMU2(I)) 
!>              has a corresponding decriptor in the multiscale framework AME(I)
!>                 AME(I) provides informations on 'parent' and 'children' simus   
!>                                           SIMU2(1) is the 'master' simulation  
!>                                   SIMU2(2:Nsimuslave+1) are 'slave' simulations
!>  
!>
  type RANKS
     integer                            :: N   !> number of procs
     integer, allocatable, dimension(:) :: r   !> proc indexes (1:N)
  end type RANKS

  type PROLONGATION
     integer, dimension(3)              :: Imin,Imax !> indices (min,max) "Gros" des pinceaux "Fins"
                                                     !! 'appui large' pour permettre l'interpolation
     integer, allocatable, dimension(:) :: Scount,Sdisp,Rcount,Rdisp
                                                     !> Vecteurs pour le MPI_ALLTOALLV
     integer,allocatable,dimension(:,:) :: parent_data ! lignes [nrank, kmin(3),kmax(3)] utile pour recv
     integer,allocatable,dimension(:,:) :: child_data  ! lignes [nrank, kmin(3),kmax(3)] utile pour send
                                                       ! tableaux de taille (:,7)
  end type PROLONGATION
  
  type REDUCTION
     integer, dimension(3)              :: Imin,Imax !> indices (min,max) "Gros" des pinceaux "Fins"
     integer, allocatable, dimension(:) :: Scount,Sdisp,Rcount,Rdisp
                                                       !> Vecteurs pour le MPI_ALLTOALLV
     integer,allocatable,dimension(:,:) :: parent_data ! lignes [nrank, kmin(3),kmax(3)] utile pour recv
     integer,allocatable,dimension(:,:) :: child_data  ! lignes [nrank, kmin(3),kmax(3)] utile pour send
                                                       ! tableaux de taille (:,7)
  end type REDUCTION
  
  type MULTISCALE
  
     !> INPUT FILE for the simulation  
     character(len=1000)                :: input_file=""!> input file (simuxx.in in rep_ame)(except for AME(1))
     
     !> FAMILY DESCRIPTION (global)
     integer                            :: Nchildren    !> number of 'child' simus  ("0" if no simu 'child')
     integer, allocatable, dimension(:) :: children     !> index in SIMU2 of 'child' simus (1:Nchildren)
     integer                            :: parent       !> index in SIMU2 of the 'parent' simu (-1 for 'master' simu AME(1))
     integer                            :: child_number !> index in the 'children' list of 'parent' simu 
                                                        !!     for simu I, index in AME(AME(I)%parent)%children
                                                        !!     (-1 for 'master' simu AME(1))
     
     !> ALL INDEXES below (ie AME(I)%IJKxxx) are attached to the grid of the 'current' simu SIMU2(I)
     
     !> MIN,MAX INDEXES  (global)  
     integer, dimension(3)                  :: IJKpmin=0  !> min,max indexes of the 'parent' contour (indexes of voxel corners)
     integer, dimension(3)                  :: IJKpmax=0  !! remains 0 for "master" simu SIMU2(1), which has no parent  
     integer, allocatable, dimension(:,:)   :: IJKcmin    !> min,max indexes of the 'child' contours (indexes of voxel corners)
     integer, allocatable, dimension(:,:)   :: IJKcmax    !!      array(3,Nchildren) <- one min/max per child
     
     !> INDEXES OF THE BOX (surface boudary, and volume)(local per pencil)
     integer, allocatable, dimension(:,:)   :: IJKsp     !> indexes of the 'parent' contour surface (indexes of voxel corners)
                                                        !!      except for "master" which has no parent 
                                                        !!      (local per pencil), elements in IJKsp(:,1:3) 
     integer, allocatable, dimension(:,:,:) :: IJKsc     !> indexes of the 'child' contour surfaces (indexes of voxel corners)
                                                        !!      one IJKc(:,1:3,ichild) per 'child'
                                                        !!      (local per pencil), elements in IJKsc(:,1:3,1:Nchildren) 
     integer, allocatable, dimension(:,:)   :: IJKvp     !> indexes of the 'parent' volume (indexes of voxel corners)
                                                        !!      except for "master" which has no parent 
                                                        !!      (local per pencil), elements in IJKvp(:,1:3) 
     integer, allocatable, dimension(:,:,:) :: IJKvc     !> indexes of the 'child' volumes (indexes of voxel corners)
                                                        !!      one IJKc(:,1:3,ichild) per 'child'
                                                        !!      (local per pencil), elements in IJKvc(:,1:3,1:Nchildren) 

     !> REFINEMENT
     integer                                :: refine = 0      !> refinement coeff wrt the parent simu 
                                                               !! remains 0 for the "master" simu SIMU2(1)

     !> POSITIONS OF CONTOURS (global)
     real(mytype), dimension(3)                :: XYZpmin = 0 !> min coordinates of the 'parent' contour
     real(mytype), dimension(3)                :: XYZpmax = 0 !> max         -> remain 0 for the "master" simu SIMU2(1)
     real(mytype), allocatable, dimension(:,:) :: XYZcmin     !> min coordinates of the 'child' contours (1:3,1:Nchildren)
     real(mytype), allocatable, dimension(:,:) :: XYZcmax     !> max 
     
     !> FIELDS FOR INTERPOLATION  Current -> Parent (i.e. Fine -> Gros) (local)
     integer, dimension(3)                     :: deltafi     !> increment of fine voxel to make coincide the first gros voxel 
                                                              !! with a fine voxel in the pencil
                                                              !! deltafi=x_1er_voxel_gros - xstart (TO BE VEIFIED)
			                                                  !! (floor(xstart/refine)*refine + 1 + floor(refine/2)) - xstart
     integer, dimension(3)                     :: dimfi       !> dimensions of tmpfi 
                                                              !! dimfi = (floor((xend - xstart+deltafi)/refine)) (TO BE VERIFIED)
     !real(mytype), allocatable,dimension(3)    :: tmpfi       !> interpolated field  
                                                              !! tempfi=field(xstart+deltafi:xend:refine,idem,idem) (TO BE VERIFIED)
     
     !> MPI RANKS INFORMATIONS FOR SEND/RECEIVE (local per pencil)
     type(RANKS)                             :: ranks_parent !> ranks_parent%N : number of pencils of the 'parent' simu
                                                             !!                  overlying (at least partly) the pencil of the 'current' simu 
                                                             !!                  maximum 4, minimum 1, except 0 for "master" SIMU2(1)
                                                             !! ranks_parent%r : ranks in the 'parent' simu
     type(RANKS),allocatable, dimension(:)   :: ranks_child  !> ranks_child(ichild)%N : number of pencils of the 'child' simu
                                                             !!                  overlying (at least partly) the pencil of the 'current' simu
                                                             !!                  minimum 0, maximum ((refine of the child+1)^2) 

     !> MPI RANKS INFORMATIONS FOR MPI_ALLTOALLV (local per pencil)
     !integer,dimension(nproc)                :: Scount, Sdisp !! Send from current to child (PREVOIR plusieurs child)

     !> MPI SIZE INFORMATIONS FOR MPI_ALLTOALLV (global)
     !> N on every proc
     integer, allocatable, dimension(:)      :: sizeIJKsp_proc !> sizeIJKp_proc(nbproc) (size(IJKsp) then MPI_ALLGATHER)
     integer, allocatable, dimension(:,:)    :: sizeIJKsc_proc !> sizeIJKc_proc(nbproc,ichild) (size(IJKsc(:,:,ichild)) then MPI_ALLGATHER)
     integer, allocatable, dimension(:)      :: sizeIJKvp_proc !> sizeIJKp_proc(nbproc) (size(IJKvp) then MPI_ALLGATHER)
     integer, allocatable, dimension(:,:)    :: sizeIJKvc_proc !> sizeIJKc_proc(nbproc,ichild) (size(IJKvc(:,:,ichild)) then MPI_ALLGATHER)
                   
     !> MPI INFORMATION FOR PROLONGATION WITH MPI_ALLTOALLV
     type(PROLONGATION),allocatable,dimension(:):: prolong    
     type(REDUCTION)                            :: reduc                          
  end type MULTISCALE
  
  !!------------------------------------------------------------------------------
  !>                               IF ADDITIONNAL FIELDS REQUIRED BY NB (POST-DOC)
  type MULTISCALE_NB
  end type MULTISCALE_NB

!!------------------------------------------------------------------------------
!>                                                                 VARIABLES AME
  
!! SIMU2 vector of simulations: SIMU2(1) = simu master
!>                              SIMU2(2:Nsimuslave+1) : slave simulations
!>                              initialized in read_AME(SIMU2, AME, "ame_rep")
!>
!> AME vector of MULTISCALE descriptors : simulatin SIMU2(I) is decribed by AME(I)
  type(SIMU_AMITEX), dimension(:), allocatable, target :: SIMU2 
  type(MULTISCALE), dimension(:), allocatable          :: AME 
  type(MULTISCALE_NB), dimension(:), allocatable       :: AME_NB ! if required by NB (post-doc)
  integer, allocatable, dimension(:,:)                 :: Gmin, Gmax  ! MPI_ALLGATHER sur xstart,xend -> taille(3,nproc)
                                                                      ! peut etre supprime avec initialisation     
  real(mytype),allocatable,dimension(:)                :: Rbuf, Sbuf                                                        
  real(mytype),allocatable,dimension(:,:,:)            :: FPc ! Field Parent on Child pencil (size (Imin:Imax) on the 3 directions)  
                                                              ! intermediate field used in PROLONGATION                                                     
  real(mytype),allocatable,dimension(:,:,:)            :: FCp ! Field Child on parent pencil (size (Imin:Imax) on the 3 directions)   
                                                              ! intermediate field used in REDUCTION                                                   
contains

!===============================================================================
!>                  init_user_variables :  ALLOCATE AND INITIALIZE PUBLIC FIELDS 
!!
!!
!===============================================================================
subroutine init_user_variables()

  implicit none
  
  character(len=1000) :: rep_AME          !< from user_param%p_string(1)
  integer             :: Nsimu            !< to deduce from read_Nsimu_AME
  integer             :: Isimu            
  
  rep_AME = user_param%p_string(1)
if (nrank==0) print *,"rep_AME =",trim(rep_AME),nrank  
  !                                                              read_simus_AME
  !                                                   partial definition of AME
  !                                                      -> children and parent 
  !----------------------------------------------------------------------------
  call read_simus_AME0(rep_AME) ! simple example with two scale
  !call read_simus_AME(rep_AME) ! TODO : extension to any arbitrary rep_ame 
  
  !                                                                    read_AME
  !                                                       full definition SIMU2 
  !                                                 + partial definition of AME 
  !----------------------------------------------------------------------------
  call read_AME0(rep_AME) ! simple example with two scale
  !call read_AME(rep_AME) ! TODO : extension to any arbitrary rep_ame 
  
  !                                                                complete_AME
  !                      complete with MPI INFORMATIONS for SEND/RECEIVE in AME 
  !----------------------------------------------------------------------------
  call complete_AME()
  call print_AME()
  
  !                                Associate pointers to "master" simu SIMU2(1) 
  !----------------------------------------------------------------------------
  call nullify_pointers_simu()
  Isimu = 1
  call associate_pointers_simu(SIMU2,Isimu)

end subroutine init_user_variables


!===============================================================================
!>  read_simus_AME0 : read simulation names in rep_ame
!!                    "0" -> version for one parent / 1 children
!!  TODO : read_simus_ame, for any arbitraty "family tree"
!!
!!  Allocate SIMU2, assign SIMU2(1)=SIMU(1), deallocate SIMU  
!!  Allocate AME 
!!  Partial definitions of AME : 
!!        Nchildren, child, parent, child_number, input_file
!!   
!===============================================================================
subroutine read_simus_AME0(rep_ame)
  implicit none
  character(len=*), intent(in) :: rep_ame
  integer                      :: Nsimu,Isimu
    
  !            Allocate SIMU2/AME and deallocate SIMU
  !--------------------------------------------------  
  Nsimu = 2
  allocate(AME(1:Nsimu))
  allocate(SIMU2(1:Nsimu))
  SIMU2(1)=SIMU(1)
  deallocate(SIMU)

  !  define Nchildren, children, parent, child_number
  !--------------------------------------------------  

  !Cas SLAVE Isimu=2
  Isimu=2 
  AME(Isimu)%Nchildren     = 0      ! ici : slave n'a pas d'enfant
  allocate(AME(Isimu)%children(0))  ! 
  AME(Isimu)%parent        = 1      ! ici : parent = master = 1
  AME(Isimu)%child_number  = 1      !       slave est le 1er enfant de master
  AME(Isimu)%input_file    = "simu1.in"      


  !Cas MASTER Isimu=1
  Isimu=1
  AME(Isimu)%Nchildren    = 1       ! ici : 1 seul enfant
  allocate(AME(Isimu)%children(1))  !
  AME(Isimu)%children(1)  = 2       !       indice de la simulation 'child'
  AME(Isimu)%parent       = -1      ! ici : master n'a pas de parent
  AME(Isimu)%child_number = -1      !       master n'est donc le fils de personne!
  AME(Isimu)%input_file   = ""      !       master n'a pas de fichier input
  
end subroutine read_simus_AME0

!===============================================================================
!>  read_AME0 : partial definition of SIMU2 and AME from reading files in rep_ame
!!              "0" -> version for one parent / 1 children
!!  TODO : read_ame, for any arbitraty "family tree"
!!
!!  Complete definition of SIMU2 
!!  Partial definitions of AME : 
!!        IJKpmin, IJKpmax
!!        IJKcmin, IJKcmax
!!        XYZpmin, XYZpmax
!!        XYZcmin, XYZcmax
!!        refine
!!        IJKsp, IJKsc, IJKvp, IJKvc (local per pencil)
!!        
!===============================================================================
subroutine read_AME0(rep_ame)
  implicit none
  character(len=*), intent(in) :: rep_ame
  character(len=200)           :: err_msg, version = "user_ame2_0.0.0"  
  character(len=1000)          :: cmd_file         !< Command file name
  integer                      :: Fcmd
  
  integer, dimension(3)        :: IJKpmin, IJKpmax, Imin, Imax
  real(mytype),dimension(3)    :: refine,resol,resolp,corner,ratio
  integer,dimension(3)         :: Irefine
  
  integer                      :: Isimu, Ichild !indices de boucle
  integer                      :: Sparent, Schild
  integer                      :: ios
  logical                      :: PBratio,T1,T2,T3,T4,T5,T6,T7,T8
  integer                      :: i,j,k,I0,J0,K0, S0,S00,S1,S2,S3,S4,S5,S6,S7,S8,Stot,Stot0
  integer,dimension(:,:),allocatable :: tmp_alloc

  !!  3 variables inutiles dans init_simu_command_line (car presence cmd_file)
  integer            :: coeff2print =0   !< indice du coeff a sortir (ligne de commande) 
  integer            :: varint2print=0   !< indice de varint a sortir (ligne de commande) 
  logical            :: help=.false.     !< test pour avoir l'aide et sortir du code
  
  !! namelist pour lecture dans fichier simuxx.in
  namelist /AMEinfo/ IJKpmin, IJKpmax
  
  IJKpmin = 0
  IJKpmax = 0
  
!---------------------------------- TODO : extension from 2 scales to arbitary AME 
!NOW : two scales

!=======================================================================
do Isimu = 1, size(AME)                 !DEBUT 1ere boucle sur les simus

  !   LOAD 'slave' simulations in SIMU2 from cmd_file
  !--------------------------------------------------
  if (Isimu .ne. 1) then 
    cmd_file=trim(rep_ame)//"/"//trim(AME(Isimu)%input_file)
    call init_simu_command_line(SIMU2(Isimu),version, help,coeff2print,varint2print,&
                 trim(cmd_file))
    SIMU2(Isimu)%simu_name=trim(AME(Isimu)%input_file)  !rename simulation for output
  end if
  
  !                      READ parent contour boundary
  !                               AME%IJKpmin,IJKpmax
  !                               AME%XYZpmin,XYZpmax
  !-------------------------------------------------- 
  if (Isimu .ne. 1) then 
    ! pas de test d'existence de file_cmd : deja fait dans init_simu_command_line
    open(newunit=Fcmd, file=trim(cmd_file),form="formatted", status="old", action="read",iostat= ios)
    if ( ios /= 0 ) then
       write(err_msg,fmt="(3A,I0,A)") "Problem when opening file : ",trim(cmd_file),&
            " , ",ios," (amitex_user_mod / read_AME)"
       call amitex_abort(trim(err_msg),2,0) 
    end if
    read(Fcmd, nml=AMEinfo,iostat= ios)
    if ( ios /= 0 ) then
       write(err_msg,fmt="(3A,I0,A)") "Problem when reading file : ",trim(cmd_file),&
            " , ",ios," (amitex_user_mod / read_AME)"
       call amitex_abort(trim(err_msg),2,0) 
    end if
    close(Fcmd)
    
    AME(Isimu)%IJKpmin = IJKpmin
    AME(Isimu)%IJKpmax = IJKpmax
    
    resol =  (/ SIMU2(Isimu)%grid%dx, SIMU2(Isimu)%grid%dy, SIMU2(Isimu)%grid%dz /)
    corner = (/ SIMU2(Isimu)%grid%x0, SIMU2(Isimu)%grid%y0, SIMU2(Isimu)%grid%z0 /)
    
    AME(Isimu)%XYZpmin = corner + (IJKpmin-1) * resol
    AME(Isimu)%XYZpmax = corner + (IJKpmax-1) * resol   
    
  end if
  
end do ! FIN 1ere boucle sur les simus

!=======================================================================
do Isimu = 1, size(AME)                 !DEBUT 2eme boucle sur les simus

  !                                 DEDUCE refinement
  !                                        AME%refine
  !--------------------------------------------------   
  if (Isimu .ne. 1) then  
     Sparent = AME(Isimu)%parent
     resol = (/ SIMU2(Isimu)%grid%dx, SIMU2(Isimu)%grid%dy, SIMU2(Isimu)%grid%dz /)
     resolp = (/ SIMU2(Sparent)%grid%dx, SIMU2(Sparent)%grid%dy, SIMU2(Sparent)%grid%dz /)
     refine = resolp / resol

     ! arrondi entier
     Irefine = anint(refine)
     ! detections d'erreur : refinement must be integer
     !                       same refinements in x, y and z directions
     if (maxval(abs(refine - Irefine)) > 1e-4 * maxval(Irefine)) then
        !TODO call amitex_abort(a remplir)
     elseif ((Irefine(1) .ne. Irefine(2)) .or. (Irefine(2) .ne. Irefine(3))) then
        !TODO call amitex_abort(a remplir)     
     else
        AME(Isimu)%refine = Irefine(1)
     end if   
  end if

end do ! FIN 2eme boucle sur les simus

!=======================================================================
do Isimu = 1, size(AME)                 !DEBUT 3eme boucle sur les simus

  !                     DEDUCE MIN-MAX child CONTOURS
  !                              AME%XYZcmin, XYZcmax
  !                              AME%IJKcmin, IJKcmax
  !--------------------------------------------------   
  allocate(AME(Isimu)%IJKcmin(3,AME(Isimu)%Nchildren))
  allocate(AME(Isimu)%IJKcmax(3,AME(Isimu)%Nchildren))
  allocate(AME(Isimu)%XYZcmin(3,AME(Isimu)%Nchildren))
  allocate(AME(Isimu)%XYZcmax(3,AME(Isimu)%Nchildren))
    
  if (AME(Isimu)%Nchildren .ne. 0) then 
    do Ichild=1,AME(Isimu)%Nchildren
    
      Schild = AME(Isimu)%children(Ichild)
      resol = (/ SIMU2(Isimu)%grid%dx, SIMU2(Isimu)%grid%dy, SIMU2(Isimu)%grid%dz /)
      PBratio=.false.
      
      ! contour 'fils' de la simu courante = contour 'pere' de la simu fille!
      AME(Isimu)%XYZcmin(:,Ichild) =  AME(Schild)%XYZpmin
      AME(Isimu)%XYZcmax(:,Ichild) =  AME(Schild)%XYZpmax
      
      ! coordonnes du coin min
      ratio = AME(Isimu)%XYZcmin(:,Ichild) / resol
      Imin = anint(ratio) 
      if ((maxval(abs(ratio) - Imin)) > 1e-2) PBratio=.true.
      
      ! coordonnes du coin max
      ratio = AME(Isimu)%XYZcmax(:,Ichild) / resol
      Imax = anint(ratio) 
      if ((maxval(abs(ratio) - Imax)) > 1e-2) PBratio=.true.
      
      if (PBratio) then
         !TODO call amitex_abort
      end if
      
      AME(Isimu)%IJKcmin(:,Ichild) = Imin + 1
      AME(Isimu)%IJKcmax(:,Ichild) = Imax + 1
        
    end do
  end if
end do ! FIN 3eme boucle sur les simus

!=======================================================================
do Isimu = 1, size(AME)                 !DEBUT 4eme boucle sur les simus

  !            DEDUCE COMPLETE parent CONTOUR & VOLUMES
  !                                  -> AME%IJKsp,IJKvp
  !----------------------------------------------------
  if (Isimu .ne. 1) then
     call box_boundary(xstart,xend,AME(Isimu)%IJKpmin,AME(Isimu)%IJKpmax,AME(Isimu)%IJKsp)
     call box_volume(xstart,xend,AME(Isimu)%IJKpmin,AME(Isimu)%IJKpmax,AME(Isimu)%IJKvp)
  end if
  
  !            DEDUCE COMPLETE child CONTOURS & VOLUMES
  !                                  -> AME%IJKsc,IJKvc
  !----------------------------------------------------
  if (AME(Isimu)%Nchildren .ne. 0) then
  do Ichild=1,AME(Isimu)%Nchildren
     ! contour
     call box_boundary(xstart,xend,AME(Isimu)%IJKcmin(:,Ichild), &
                       AME(Isimu)%IJKcmax(:,Ichild),tmp_alloc)
     if (Ichild==1) allocate(AME(Isimu)%IJKsc(size(tmp_alloc,1),3,AME(Isimu)%Nchildren))
     AME(Isimu)%IJKsc(:,:,Ichild) = tmp_alloc
     deallocate(tmp_alloc)
     ! volume
     call box_volume(xstart,xend,AME(Isimu)%IJKcmin(:,Ichild), &
                       AME(Isimu)%IJKcmax(:,Ichild),tmp_alloc)
     if (Ichild==1) allocate(AME(Isimu)%IJKvc(size(tmp_alloc,1),3,AME(Isimu)%Nchildren))
     AME(Isimu)%IJKvc(:,:,Ichild) = tmp_alloc
     deallocate(tmp_alloc)
  end do    
  end if
  
end do
!---------------------------------------------------------------------------------------
!---- EXEMPLE D'UTILISATION DES INDICES IJK(s-v)(c-p)
!---- VERIF : AME(Isimu)%IJKsp, AME(Isimu)%IJKsc,AME(Isimu)%IJKvp, AME(Isimu)%IJKvc
!----         -> LES CHAMPS CONSTRUITS SEMBLENT CORRECTS 
!                (il faudrait modifier print_field_vtk pour prendre en compte
!                la position des grilles)
 
if (1 == 1) then
   block
   use green_mod
   real(kind=REAL32),pointer,dimension(:,:,:)     :: tmp32_3D
   
   Isimu = 2
   SIMU2(Isimu)%TEMPfield32 = 0.
   !pointer utile ici pour reprofilage
   tmp32_3D(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) => SIMU2(Isimu)%TEMPfield32
   do i=1,size(AME(Isimu)%IJKsp,1)
      tmp32_3D(AME(Isimu)%IJKsp(i,1), &
               AME(Isimu)%IJKsp(i,2), &
               AME(Isimu)%IJKsp(i,3)) = 1.
   end do
   do i=1,size(AME(Isimu)%IJKvp,1)
      tmp32_3D(AME(Isimu)%IJKvp(i,1), &
               AME(Isimu)%IJKvp(i,2), &
               AME(Isimu)%IJKvp(i,3)) = 2.
   end do
   grid => SIMU2(Isimu)%grid
   call print_field_vtk(SIMU2(Isimu)%TEMPfield32,"loc_boxboundary.vtk","ID")
   
   Isimu = 1; Ichild = 1
   SIMU2(Isimu)%TEMPfield32 = 0.
   !pointer utile ici pour reprofilage
   tmp32_3D(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) => SIMU2(Isimu)%TEMPfield32
   do i=1,size(AME(Isimu)%IJKsc,1)
      tmp32_3D(AME(Isimu)%IJKsc(i,1,Ichild), &
               AME(Isimu)%IJKsc(i,2,Ichild), &
               AME(Isimu)%IJKsc(i,3,Ichild)) = 1.
   end do
   do i=1,size(AME(Isimu)%IJKvc,1)
      tmp32_3D(AME(Isimu)%IJKvc(i,1,Ichild), &
               AME(Isimu)%IJKvc(i,2,Ichild), &
               AME(Isimu)%IJKvc(i,3,Ichild)) = 2.
   end do
   grid => SIMU2(Isimu)%grid
   call print_field_vtk(SIMU2(Isimu)%TEMPfield32,"mac_boxboundary.vtk","ID")
  
   end block
end if
!---- FIN VERIF -----------------------------------------------------------------------
!--------------------------------------------------------------------------------------

end subroutine read_AME0

!===============================================================================
!>  complete_AME : complete AME with informations for MPI_ALLTOALLV
!!                    
!===============================================================================
subroutine complete_AME()
  implicit none
  integer :: Isimu, Ichild, child
  integer :: tmp
  integer :: ierror
  
  !                                          PARTAGE DES DIMENSIONS DES PINCEAUX 
  !-----------------------------------------------------------------------------
  allocate(Gmin(3,nproc))
  allocate(Gmax(3,nproc))
  call MPI_ALLGATHER(xstart,3,MPI_INTEGER,Gmin,3,MPI_INTEGER,MPI_COMM_WORLD,ierror)
  call MPI_ALLGATHER(xend,3,MPI_INTEGER,Gmax,3,MPI_INTEGER,MPI_COMM_WORLD,ierror)
  
  !                                               ALLOCATION-AFFECTATION prolong  
  !-----------------------------------------------------------------------------
  do Isimu=1,size(AME)
     allocate(AME(Isimu)%prolong(AME(Isimu)%Nchildren))
  end do
  
  
  do Isimu=1,size(AME)
  if (AME(Isimu)%Nchildren .ne. 0) then
     do child = 1,AME(Isimu)%Nchildren
        call init_prolongIminImax(Isimu,child,1)
        call init_prolongRecv(Isimu,child)
        call init_prolongSend(Isimu,child)
     end do
  end if
  end do
  
  !                                                            AFFECTATION reduc  
  !-----------------------------------------------------------------------------
  do Isimu=2,size(AME) ! Isimu = 1 (master) ne realise pas de reduction
     call init_reducIminImax(Isimu)
     call init_reducRecv(Isimu)
     call init_reducSend(Isimu)
  end do
  
  !                                                                 VERIFICATION  
  !-----------------------------------------------------------------------------
  if (1 == 0) then
     Isimu=1
     child=1
     call print_prolong(Isimu,child)
  end if
  
  
  !                                                             TEST UTILISATION  
  !-----------------------------------------------------------------------------
  if (1 == 1) then
     Isimu = 1
     child = 1
     call test_prolong_reduc1(Isimu,child)
     call test_prolong_reduc2(Isimu,child)
  end if
  

end subroutine complete_AME



!===============================================================================
!>  print_AME : print the AME structure on stdout
!!
!!        
!===============================================================================
subroutine print_AME()
  implicit none
  
  integer :: Isimu,Ichild
  
  if (nrank==0) then
  
  write(output_unit,*) ACHAR(10)//"                               STRUCTURE AME"
  write(output_unit,*) "============================================"//ACHAR(10)

  do Isimu=1,size(AME)
     write(output_unit,*) "SIMULATION ",Isimu
     write(output_unit,*) "---------------------------------------------"
     write(output_unit,*) "   refine  ", AME(Isimu)%refine      
     write(output_unit,*) "FAMILY DESCRIPTION  "     
     write(output_unit,*) "   parent       ", AME(Isimu)%parent     
     write(output_unit,*) "   Nchildren    ", AME(Isimu)%Nchildren     
     write(output_unit,*) "   children     ", AME(Isimu)%children     
     write(output_unit,*) "   child_number ", AME(Isimu)%child_number   
     write(output_unit,*) "PARENT CONTOUR  "          
     write(output_unit,*) "   IJKpmin ", AME(Isimu)%IJKpmin
     write(output_unit,*) "   IJKpmax ", AME(Isimu)%IJKpmax
     write(output_unit,*) "   XYZpmin ", AME(Isimu)%XYZpmin
     write(output_unit,*) "   XYZpmax ", AME(Isimu)%XYZpmax
     do Ichild = 1, AME(Isimu)%Nchildren
     write(output_unit,*) "CHILD CONTOUR  ", Ichild          
     write(output_unit,*) "   IJKcmin ", AME(Isimu)%IJKcmin(:,Ichild)
     write(output_unit,*) "   IJKcmax ", AME(Isimu)%IJKcmax(:,Ichild)
     write(output_unit,*) "   XYZcmin ", AME(Isimu)%XYZcmin(:,Ichild)
     write(output_unit,*) "   XYZcmax ", AME(Isimu)%XYZcmax(:,Ichild)
     end do
     
     write(output_unit,*) ""

  end do
  
  end if !nrank==0
end subroutine print_AME



!===============================================================================
!>  box_boundary : defines points of a box boundary given by (min,max) indexes
!!
!!
!! /param[in]  xstart  starting coordinates of the pencil, integer(3)
!! /param[in]  xend    ending coordinates of the pencil, integer(3)
!! /param[in]  IJKmin  min corner indexes of the box, integer(3)
!! /param[in]  IJKmax  max corner indexes of the box, integer(3)
!! /param[out] IJK     coordinates the box boundary points, integer(:,3)
!!   
!===============================================================================
subroutine box_boundary(xstart,xend,IJKmin,IJKmax,IJK)
  implicit none
   
  integer, dimension(3),intent(in)                  :: xstart,xend
  integer, dimension(3),intent(in)                  :: IJKmin,IJKmax
  integer, allocatable, dimension(:,:),intent(out)  :: IJK
   
  integer :: I0,J0,K0, Stot,Stot0,S0,S00,S1,S2,S3,S4,S5,S6,S7,S8
  integer :: i,j,k
  logical :: T1,T2,T3,T4,T5,T6,T7,T8
  integer,dimension(3) :: Imin,Imax
  integer,dimension(:,:),allocatable :: tmp_alloc


  !                               ALLOCATION PAR DEFAUT
  !            (si le contour n'est pas sur le pinceau)
  !----------------------------------------------------
  allocate(IJK(0,3))

  
  ! SURFACES  min et max X, Y et Z (sans bord ni coin)
  Stot  = 0
  Stot0 = 0
  do I0 = 1,3
    if (I0 .eq. 1) then; J0=2 ; K0=3; end if
    if (I0 .eq. 2) then; J0=3 ; K0=1; end if
    if (I0 .eq. 3) then; J0=1 ; K0=2; end if
    
    T1    = xstart(I0) <= IJKmin(I0) .AND. &
            IJKmin(I0) <= xend(I0)
    T2    = xstart(I0) <= IJKmax(I0) .AND. &
            IJKmax(I0) <= xend(I0)

    Imin(J0) = maxval([IJKmin(J0)+1, xstart(J0)])
    Imax(J0) = minval([IJKmax(J0)-1, xend(J0)])
    Imin(K0) = maxval([IJKmin(K0)+1, xstart(K0)])
    Imax(K0) = minval([IJKmax(K0)-1, xend(K0)])
    
    S1 = 0
    S2 = 0
    if (T1) S1 = (Imax(J0) - Imin(J0) + 1) * (Imax(K0) - Imin(K0) + 1) 
    if (T2) S2 = (Imax(J0) - Imin(J0) + 1) * (Imax(K0) - Imin(K0) + 1) 
    Stot = Stot0 + S1 + S2

    if(allocated(IJK)) then
       tmp_alloc = IJK              ! reallocation "a la volee" (fortran 2003)
       deallocate(IJK)
       allocate(IJK(Stot,3))
       IJK(1:Stot0,:) = tmp_alloc(:,:)  
    else
       allocate(IJK(Stot,3))
    end if
    
    if (S1 >0) then
      IJK(Stot0+1:Stot0+S1,I0) = IJKmin(I0)
      IJK(Stot0+1:Stot0+S1,J0) = (/  ((j,j=Imin(J0),Imax(J0)),k=Imin(K0),Imax(K0))  /) 
      IJK(Stot0+1:Stot0+S1,K0) = (/  ((k,j=Imin(J0),Imax(J0)),k=Imin(K0),Imax(K0))  /) 
    end if
    if (S2 >0) then
      IJK(Stot0+S1+1:Stot,I0) = IJKmax(I0)
      IJK(Stot0+S1+1:Stot,J0) = (/  ((j,j=Imin(J0),Imax(J0)),k=Imin(K0),Imax(K0))  /) 
      IJK(Stot0+S1+1:Stot,K0) = (/  ((k,j=Imin(J0),Imax(J0)),k=Imin(K0),Imax(K0))  /) 
    end if
    
    Stot0=Stot
    
  end do

  ! LIGNES  // X, Y et Z (coin)
  do I0 = 1,3
    if (I0 .eq. 1) then; J0=2 ; K0=3; end if
    if (I0 .eq. 2) then; J0=3 ; K0=1; end if
    if (I0 .eq. 3) then; J0=1 ; K0=2; end if
    
    T1 = xstart(J0) <= IJKmin(J0) .AND. &
         IJKmin(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmin(K0) .AND. &
         IJKmin(K0) <= xend(K0) 
    T2 = xstart(J0) <= IJKmin(J0) .AND. &
         IJKmin(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmax(K0) .AND. &
         IJKmax(K0) <= xend(K0) 
    T3 = xstart(J0) <= IJKmax(J0) .AND. &
         IJKmax(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmax(K0) .AND. &
         IJKmax(K0) <= xend(K0) 
    T4 = xstart(J0) <= IJKmax(J0) .AND. &
         IJKmax(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmin(K0) .AND. &
         IJKmin(K0) <= xend(K0) 
      
    Imin(I0) = maxval([IJKmin(I0)+1, xstart(I0)])
    Imax(I0) = minval([IJKmax(I0)-1, xend(I0)])
    
    S1 = 0; S2 = 0; S3=0;S4=0
    if (T1) S1 = (Imax(I0) - Imin(I0) + 1) 
    if (T2) S2 = (Imax(I0) - Imin(I0) + 1)
    if (T3) S3 = (Imax(I0) - Imin(I0) + 1)
    if (T4) S4 = (Imax(I0) - Imin(I0) + 1)
    Stot = Stot0 + S1 + S2 + S3 + S4

    if(allocated(IJK)) then
       tmp_alloc = IJK              ! reallocation "a la volee" (fortran 2003)
       deallocate(IJK)
       allocate(IJK(Stot,3))
       IJK(1:Stot0,:) = tmp_alloc(:,:)  
    else
       allocate(IJK(Stot,3))
    end if
    
    if (S1 >0) then
      IJK(Stot0+1:Stot0+S1,I0) = (/  (i,i=Imin(I0),Imax(I0))  /)
      IJK(Stot0+1:Stot0+S1,J0) = IJKmin(J0) 
      IJK(Stot0+1:Stot0+S1,K0) = IJKmin(K0) 
    end if
    if (S2 >0) then
      IJK(Stot0+S1+1:Stot0+S1+S2,I0) = (/  (i,i=Imin(I0),Imax(I0))  /)
      IJK(Stot0+S1+1:Stot0+S1+S2,J0) = IJKmin(J0) 
      IJK(Stot0+S1+1:Stot0+S1+S2,K0) = IJKmax(K0) 
    end if
    if (S3 >0) then
      IJK(Stot0+S1+S2+1:Stot0+S1+S2+S3,I0) = (/  (i,i=Imin(I0),Imax(I0))  /)
      IJK(Stot0+S1+S2+1:Stot0+S1+S2+S3,J0) = IJKmax(J0) 
      IJK(Stot0+S1+S2+1:Stot0+S1+S2+S3,K0) = IJKmax(K0) 
    end if
    if (S4 >0) then
      IJK(Stot0+S1+S2+S3+1:Stot,I0) = (/  (i,i=Imin(I0),Imax(I0))  /)
      IJK(Stot0+S1+S2+S3+1:Stot,J0) = IJKmax(J0) 
      IJK(Stot0+S1+S2+S3+1:Stot,K0) = IJKmin(K0) 
    end if
    
    Stot0 = Stot
    
  end do
 
  ! COINS  
    I0=1; J0=2; K0=3
    T1 = xstart(I0) <= IJKmin(I0) .AND. &
         IJKmin(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmin(J0) .AND. &
         IJKmin(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmin(K0) .AND. &
         IJKmin(K0) <= xend(K0) 
    T2 = xstart(I0) <= IJKmin(I0) .AND. &
         IJKmin(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmin(J0) .AND. &
         IJKmin(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmax(K0) .AND. &
         IJKmax(K0) <= xend(K0) 
    T3 = xstart(I0) <= IJKmin(I0) .AND. &
         IJKmin(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmax(J0) .AND. &
         IJKmax(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmax(K0) .AND. &
         IJKmax(K0) <= xend(K0) 
    T4 = xstart(I0) <= IJKmin(I0) .AND. &
         IJKmin(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmax(J0) .AND. &
         IJKmax(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmin(K0) .AND. &
         IJKmin(K0) <= xend(K0) 
    T5 = xstart(I0) <= IJKmax(I0) .AND. &
         IJKmax(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmin(J0) .AND. &
         IJKmin(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmin(K0) .AND. &
         IJKmin(K0) <= xend(K0) 
    T6 = xstart(I0) <= IJKmax(I0) .AND. &
         IJKmax(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmin(J0) .AND. &
         IJKmin(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmax(K0) .AND. &
         IJKmax(K0) <= xend(K0) 
    T7 = xstart(I0) <= IJKmax(I0) .AND. &
         IJKmax(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmax(J0) .AND. &
         IJKmax(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmax(K0) .AND. &
         IJKmax(K0) <= xend(K0) 
    T8 = xstart(I0) <= IJKmax(I0) .AND. &
         IJKmax(I0) <= xend(I0)   .AND. &
         xstart(J0) <= IJKmax(J0) .AND. &
         IJKmax(J0) <= xend(J0)   .AND. &
         xstart(K0) <= IJKmin(K0) .AND. &
         IJKmin(K0) <= xend(K0) 


    S1 = 0; S2 = 0; S3=0;S4=0;S5=0;S6=0;S7=0;S8=0
    if (T1) S1 = 1
    if (T2) S2 = 1
    if (T3) S3 = 1
    if (T4) S4 = 1
    if (T5) S5 = 1
    if (T6) S6 = 1
    if (T7) S7 = 1
    if (T8) S8 = 1
    Stot = Stot0 + S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8

    if(allocated(IJK)) then
       tmp_alloc = IJK              ! reallocation "a la volee" (fortran 2003)
       deallocate(IJK)
       allocate(IJK(Stot,3))
       IJK(1:Stot0,:) = tmp_alloc(:,:)  
    else
       allocate(IJK(Stot,3))
    end if

    if (S1 >0) then
      S00=0; S0 = S1
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmin(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmin(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmin(K0) 
    end if
    if (S2 >0) then
      S00=S1; S0 = S1+S2
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmin(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmin(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmax(K0) 
    end if
    if (S3 >0) then
      S00=S1+S2; S0 = S1+S2+S3
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmin(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmax(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmax(K0) 
    end if
    if (S4 >0) then
      S00=S1+S2+S3; S0 = S1+S2+S3+S4
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmin(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmax(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmin(K0) 
    end if
    if (S5 >0) then
      S00=S1+S2+S3+S4; S0 = S1+S2+S3+S4+S5
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmax(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmin(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmin(K0) 
    end if
    if (S6 >0) then
      S00=S1+S2+S3+S4+S5; S0 = S1+S2+S3+S4+S5+S6
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmax(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmin(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmax(K0) 
    end if
    if (S7 >0) then
      S00=S1+S2+S3+S4+S5+S6; S0 = S1+S2+S3+S4+S5+S6+S7
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmax(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmax(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmax(K0) 
    end if
    if (S8 >0) then
      S00=S1+S2+S3+S4+S5+S6+S7; S0 = S1+S2+S3+S4+S5+S6+S7+S8
      IJK(Stot0+S00+1:Stot0+S0,I0) = IJKmax(J0)
      IJK(Stot0+S00+1:Stot0+S0,J0) = IJKmax(J0) 
      IJK(Stot0+S00+1:Stot0+S0,K0) = IJKmin(K0) 
    end if

end subroutine box_boundary

!===============================================================================
!>  box_volume : defines points strictly inside a box boundary 
!!               given by (min,max) indexes
!!
!! /param[in]  xstart  starting coordinates of the pencil, integer(3)
!! /param[in]  xend    ending coordinates of the pencil, integer(3)
!! /param[in]  IJKmin  min corner indexes of the box, integer(3)
!! /param[in]  IJKmax  max corner indexes of the box, integer(3)
!! /param[out] IJK     indexes of the points inside box, integer(:,3)
!!   
!===============================================================================
subroutine box_volume(xstart,xend,IJKmin,IJKmax,IJK)
  implicit none
  integer, dimension(3),intent(in)                  :: xstart,xend
  integer, dimension(3),intent(in)                  :: IJKmin,IJKmax
  integer, allocatable, dimension(:,:),intent(out)  :: IJK
  
  integer     :: i,j,k,K0 ! loop index
  
  !-------------------------------------------DECOMPTE DU NBRE DE PTS A ALLOUER
  K0 = 0
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)
     if (i > IJKmin(1) .AND. j > IJKmin(2) .AND. k > IJKmin(3) .AND. &
         i < IJKmax(1) .AND. j < IJKmax(2) .AND. k < IJKmax(3)) then
         K0 = K0 + 1
     end if 
  end do
  end do
  end do
  
  !-------------------------------------------ALLOCATION
  allocate(IJK(K0,3))

  !-------------------------------------------AFFECTATION
  K0 = 0
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)
     if (i > IJKmin(1) .AND. j > IJKmin(2) .AND. k > IJKmin(3) .AND. &
         i < IJKmax(1) .AND. j < IJKmax(2) .AND. k < IJKmax(3)) then
        k0 = K0 + 1
        IJK(K0,:) = (/ i, j , k /)
     end if 
  end do
  end do
  end do

end subroutine box_volume

!===============================================================================
!>  init_prolongIminImax : Imin,Imax indices "gros" des pinceaux "fins"
!!                         appui strict + halo-cells autour
!!                            -> necessaire pour interpolation par la suite
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!! /param[in]  halo    Nombre de couches de halo-cells (generalement 1)
!!   
!===============================================================================
subroutine init_prolongIminImax(Isimu,child,halo)
  implicit none
  integer, intent(in) :: Isimu, child, halo
  
  integer :: Ichild
  real(mytype), dimension(3) :: X0c,DXc,DX,Xmin,Xmax
  
  ! RECUPERATION DES COORD min,max du pinceau de la grille Fine (simu Child)
  Ichild = AME(Isimu)%children(child)
  X0c = (/ SIMU2(Ichild)%grid%x0,SIMU2(Ichild)%grid%y0,SIMU2(Ichild)%grid%z0 /) 
  DXc = (/ SIMU2(Ichild)%grid%dx,SIMU2(Ichild)%grid%dy,SIMU2(Ichild)%grid%dz /)

  Xmin = X0c + (xstart-1)*DXc
  Xmax = X0c + (xend-1)*DXc
  
  ! CONVERSION DANS LA GRILLE GROSSIERE (avec appui large)
  DX  = (/ SIMU2(Isimu)%grid%dx,SIMU2(Isimu)%grid%dy,SIMU2(Isimu)%grid%dz /)
  AME(Isimu)%prolong(child)%Imin = (ceiling((Xmin - 0.001*maxval(DXc)) / DX) + 1) - halo 
  AME(Isimu)%prolong(child)%Imax =  (floor((Xmax + 0.001*maxval(DXc)) / DX) + 1)  + halo

end subroutine init_prolongIminImax
!===============================================================================
!>  init_prolongIminImax : Imin,Imax indices "gros" des pinceaux "fins"
!!                         appuyee 'strict'
!!
!! /param[in]  Isimu   Index of the simulation 
!!   
!===============================================================================
subroutine init_reducIminImax(Isimu)
  implicit none
  integer, intent(in) :: Isimu
  
  integer :: Iparent
  real(mytype), dimension(3) :: X0c,DXc,DX,Xmin,Xmax
  
  Iparent = AME(Isimu)%parent
  
  ! RECUPERATION DES COORD min,max du pinceau de la grille Fine (simu courante)
  X0c = (/ SIMU2(Isimu)%grid%x0,SIMU2(Isimu)%grid%y0,SIMU2(Isimu)%grid%z0 /) 
  DXc = (/ SIMU2(Isimu)%grid%dx,SIMU2(Isimu)%grid%dy,SIMU2(Isimu)%grid%dz /)

  Xmin = X0c + (xstart-1)*DXc
  Xmax = X0c + (xend-1)*DXc
  
  ! CONVERSION DANS LA GRILLE GROSSIERE
  DX  = (/ SIMU2(Iparent)%grid%dx,SIMU2(Iparent)%grid%dy,SIMU2(Iparent)%grid%dz /)
  AME(Isimu)%reduc%Imin = ceiling((Xmin - 0.001*maxval(DXc)) / DX) + 1 
  AME(Isimu)%reduc%Imax =   floor((Xmax + 0.001*maxval(DXc)) / DX) + 1

end subroutine init_reducIminImax
!===============================================================================
!>  init_prolongRecv : determine les vecteurs Rcount,Rdisp et le tableau parent_data,
!!                   dans AME(Isimu)%prolong(child), pour la RECEPTION du ALLTOALLV 
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!!   
!===============================================================================
subroutine init_prolongRecv(Isimu,child)
  implicit none
  integer, intent(in) :: Isimu, child
  
  integer,dimension(3) :: Imin, Imax,kmin,kmax
  integer              :: i,j,K0,rank,count0,disp0
  logical              :: test
  
  ! RECUPERATION DES DONNEES
  Imin = AME(Isimu)%prolong(child)%Imin
  Imax = AME(Isimu)%prolong(child)%Imax
  
  ! AFFECTATION DE parent_data (proc d'ou vienne les data et portion de grille correspondante)
  call prolong_alloc_assign_Bdata(Gmin,Gmax,Imin,Imax,AME(Isimu)%prolong(child)%parent_data)

  ! ALLOCATION - AFFECTATION de Rcount, Rdisp (<-- parent_data)
  allocate(AME(Isimu)%prolong(child)%Rcount(nproc))
  allocate(AME(Isimu)%prolong(child)%Rdisp(nproc))
  call prolong_assign_count_disp(AME(Isimu)%prolong(child)%Rcount,&
                                 AME(Isimu)%prolong(child)%Rdisp,&
                                 AME(Isimu)%prolong(child)%parent_data)
                                   
end subroutine init_prolongRecv
!===============================================================================
!>  init_reducRecv : determine les vecteurs Rcount,Rdisp et le tableau parent_data,
!!                   dans AME(Isimu)%reduc, pour la RECEPTION du ALLTOALLV 
!!
!! /param[in]  Isimu   Index of the simulation 
!!   
!===============================================================================
subroutine init_reducRecv(Isimu)
  implicit none
  integer, intent(in) :: Isimu

  integer, dimension(3)       :: Imin,Imax
  integer, dimension(3,nproc) :: Fmin,Fmax  ! MPI_ALLGATHER sur Imin, Imax
  integer                     :: ierror
  
  ! ON DISTRIBUE LES "INDEX MIN-MAX DES PINCEAUX FINS DANS LA GRILLE GROSSE" 
  Imin = AME(Isimu)%reduc%Imin
  Imax = AME(Isimu)%reduc%Imax
  call MPI_ALLGATHER(Imin,3,MPI_INTEGER,Fmin,3,MPI_INTEGER,MPI_COMM_WORLD,ierror)
  call MPI_ALLGATHER(Imax,3,MPI_INTEGER,Fmax,3,MPI_INTEGER,MPI_COMM_WORLD,ierror)

  ! AFFECTATION DE parent_data (proc d'ou vienne les data et portion de grille correspondante)
  call prolong_alloc_assign_Bdata(Fmin,Fmax,xstart,xend,AME(Isimu)%reduc%parent_data)

  ! ALLOC/AFFECTATION de Rcount, Rdisp (<-- parent_data)
  allocate(AME(Isimu)%reduc%Rcount(nproc))
  allocate(AME(Isimu)%reduc%Rdisp(nproc))
  call prolong_assign_count_disp(AME(Isimu)%reduc%Rcount,&
                                 AME(Isimu)%reduc%Rdisp,&
                                 AME(Isimu)%reduc%parent_data)

end subroutine init_reducRecv
!===============================================================================
!>  init_reducSend : determine les vecteurs Scount,Sdisp et le tableau child_data,
!!                   dans AME(Isimu)%reduc, pour l'ENVOI du ALLTOALLV 
!!
!! /param[in]  Isimu   Index of the simulation 
!!   
!===============================================================================
subroutine init_reducSend(Isimu)
  implicit none
  integer, intent(in) :: Isimu
  
  integer, dimension(3)       :: Imin,Imax

  ! RECUPERATION DES DONNEES
  Imin = AME(Isimu)%reduc%Imin
  Imax = AME(Isimu)%reduc%Imax

  ! AFFECTATION DE child_data (proc vers qui envoyer et portion de grille a envoyer)
  call prolong_alloc_assign_Bdata(Gmin,Gmax,Imin,Imax,AME(Isimu)%reduc%child_data)
 
  ! ALLOC/AFFECTATION de Scount, Sdisp (<-- child_data)
  allocate(AME(Isimu)%reduc%Scount(nproc))
  allocate(AME(Isimu)%reduc%Sdisp(nproc))
  call prolong_assign_count_disp(AME(Isimu)%reduc%Scount,&
                                 AME(Isimu)%reduc%Sdisp,&
                                 AME(Isimu)%reduc%child_data)
end subroutine init_reducSend
!===============================================================================
!>  init_prolongSend : determine les vecteurs Scount,Sdisp et le tableau child_data,
!!                   dans AME(Isimu)%prolong(child), pour l'ENVOI du ALLTOALLV 
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!!   
!===============================================================================
subroutine init_prolongSend(Isimu,child)
  implicit none
  integer, intent(in)         :: Isimu, child
  
  integer, dimension(3)       :: Imin,Imax
  integer, dimension(3,nproc) :: Fmin,Fmax  ! MPI_ALLGATHER sur Imin, Imax
  integer                     :: ierror
  
  ! ON DISTRIBUE LES DIM DE LA GRILLE GROSSES SUPERPOSES AUX PINCEAUX "FINS" 
  Imin = AME(Isimu)%prolong(child)%Imin
  Imax = AME(Isimu)%prolong(child)%Imax
  call MPI_ALLGATHER(Imin,3,MPI_INTEGER,Fmin,3,MPI_INTEGER,MPI_COMM_WORLD,ierror)
  call MPI_ALLGATHER(Imax,3,MPI_INTEGER,Fmax,3,MPI_INTEGER,MPI_COMM_WORLD,ierror)

  ! AFFECTATION DE child_data (proc vers qui envoyer et portion de grille a envoyer)
  call prolong_alloc_assign_Bdata(Fmin,Fmax,xstart,xend,AME(Isimu)%prolong(child)%child_data)
  
  ! ALLOC/AFFECTATION de Scount, Sdisp (<-- child_data)
  allocate(AME(Isimu)%prolong(child)%Scount(nproc))
  allocate(AME(Isimu)%prolong(child)%Sdisp(nproc))
  call prolong_assign_count_disp(AME(Isimu)%prolong(child)%Scount,&
                                 AME(Isimu)%prolong(child)%Sdisp,&
                                 AME(Isimu)%prolong(child)%child_data)
end subroutine init_prolongSend
!===============================================================================
!>  prolong_alloc_assign_Bdata : 
!!          allocate and assign  prolong()%child_data or parent_data
!!
!!          identifies intersection between box I (current proc)
!!          and boxes H defined on all the procs.
!!
!!
!! /param[in]  Hmin    integer(3,nproc), min indexes of boxes, one per proc 
!! /param[in]  Hmax    integer(3,nproc), max indexes of boxes, one per proc 
!! /param[in]  Imin    integer(3), min indexes of a box (on the current proc) 
!! /param[in]  Imax    integer(3), max indexes of a box (on the current proc) 
!! /param[out] Bdata   lines of integer vectors [rk, kmin(3), kmax(3)]
!!                     one line if the box H of proc rk overlaps box I
!!                     kmin-kmax : min-max indexes of the overlap region
!!                     
!===============================================================================
subroutine prolong_alloc_assign_Bdata(Hmin,Hmax,Imin,Imax,Bdata)
  implicit none
  integer,dimension(3,nproc),intent(in)            :: Hmin,Hmax 
  integer,dimension(3),intent(in)                  :: Imin,Imax 
  integer, allocatable, dimension(:,:),intent(out) :: Bdata
  
  integer,dimension(3) :: kmin,kmax
  integer              :: K0,i,j,alloc_assign
  logical              :: test
  
  ! 1ere passage : on alloue
  ! 2eme passage : on affecte
  do alloc_assign=1,2
    K0 = 0
    test=.true.
    do i=1,nproc
       do j=1,3
          kmin(j)=maxval((/ Hmin(j,i), Imin(j) /))
          kmax(j)=minval((/ Hmax(j,i), Imax(j) /))
          if(kmin(j)>Hmax(j,i)) test=.false.
          if(kmax(j)<Hmin(j,i)) test=.false.
       end do
       if (test) then
         K0 = K0 + 1
         if (alloc_assign == 2) then ! assign
           Bdata(K0,1)   = i-1
           Bdata(K0,2:4) = kmin(:)
           Bdata(K0,5:7) = kmax(:)
         end if
       end if
       test=.true. 
    end do
    if (alloc_assign == 1) allocate(Bdata(K0,7)) !allocate
  end do
  
end subroutine prolong_alloc_assign_Bdata
!===============================================================================
!>  prolong_assign_count_disp : 
!!          read prolong()%child_data or parent_data
!!          and deduce the corresponding [Scount,Sdisp] or [Rcount,Rdisp]
!!
!! /param[in]  Bdata   lines of integer vectors [rk, kmin(3), kmax(3)]
!!                     one line if the box H of proc rk overlaps box I
!!                     kmin-kmax : min-max indexes of the overlap region
!! /param[out] Bcount   integer(nproc), count for mpi_alltoallv
!! /param[out] Bdisp    integer(nproc), disp for mpi_alltoallv
!!                     
!===============================================================================
subroutine prolong_assign_count_disp(Bcount,Bdisp,Bdata)
  implicit none
  integer, dimension(:,:),intent(in)  :: Bdata
  integer, dimension(:),intent(out)   :: Bcount, Bdisp
  
  integer,dimension(3) :: kmin,kmax
  integer              :: i,rank, count0, disp0
  
  Bcount = 0
  Bdisp  = 0
  disp0=0
  count0=0
  do i=1,size(Bdata,1)
     rank = Bdata(i,1)
     kmin = Bdata(i,2:4)
     kmax = Bdata(i,5:7)
     count0 = product(kmax-kmin+1)
     Bcount(rank+1) = count0
     Bdisp(rank+1)  = disp0
     disp0 = disp0 + count0
  end do
  
end subroutine prolong_assign_count_disp

!===============================================================================
!> test_prolong_reduc1 
!!
!!  TEST PROLONGATION
!!  On construit FP champ gros (Parent) constant par pinceau
!!                     sur la simu courante Isimu
!!                     On prolonge sur la grille fine de la simu (Child) pour 
!!                     construire FC.
!!                     Le champ est prolonge par "injection" (valeur identique 
!!                     aux noeuds communs entre grilles)
!!                     
!!  => Sorties VTK : comparaison FP FC, les résultats semblent OK (10/6/2021)
!!
!!  TEST REDUCTION
!!  On part du champ FC, -1 partout sauf aux noeuds "injectés"
!!  On initialise FP a -2
!!  On réalise la réduction pour construire un nouveau champ FP
!!  Apres reduction, on retrouve FP avec :
!!        FP identique au champ AVANT prolongation dans la zone "fine"
!!        FP = -2 ailleurs  
!!
!!  => Sorties VTK : comparaison FP FC, les résultats semblent OK (14/6/2021)
!!
!!
!! /param[in]  Isimu   Index of the current simulation 
!! /param[in]  child   Child number
!!   
!===============================================================================
subroutine test_prolong_reduc1(Isimu,child)
  implicit none
  integer, intent(in)         :: Isimu, child
  
  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) ::&
                                 FP,FC  ! Field Parent, Field Child
  
  integer :: Ichild,ierror,Isimu2
  
  !==============================================================================
  !                                                             TEST PROLONGATION
  !                                                                    FP -> FC
   
  !                             INITIALIZATION FP : constant par pinceau
  !---------------------------------------------------------------------
  FP = nrank
  FC = 0.
  
  !                                                          PACK_BUFFER
  !                                                           FP -> Sbuf
  !---------------------------------------------------------------------
  call pack_buffer_prolong(Isimu,child,FP)
  
  !                                                            ALLTOALLV   
  !                                                         Sbuf -> Rbuf    
  !---------------------------------------------------------------------
  call MPI_ALLTOALLV(Sbuf,AME(Isimu)%prolong(child)%Scount,AME(Isimu)%prolong(child)%Sdisp,real_type,&
                     Rbuf,AME(Isimu)%prolong(child)%Rcount,AME(Isimu)%prolong(child)%Rdisp,real_type,&
                     MPI_COMM_WORLD,ierror)

  !                                                        UNPACK_BUFFER
  !                                                          Rbuf -> FPc
  !---------------------------------------------------------------------
  call unpack_buffer_prolong(Isimu,child)
  
  !                                                          INTERPOLATE
  !                                                           FPc -> FC
  !---------------------------------------------------------------------
  call interp_prolong(Isimu,child,"injection",FC)
  
  !                                                          SORTIES VTK 
  !---------------------------------------------------------------------
  !   Rq : on doit basculer les pointeurs sur la bonne simu pour avoir 
  !        les bonnes valeurs de dx,dy,dz,x0,y0,z0
  Ichild = AME(Isimu)%children(child)
  call associate_pointers_simu(SIMU2,Isimu)
  call print_field_vtk(reshape(FP,(/xsize(1)*xsize(2)*xsize(3)/)),"parent.vtk","ID")
  
  call associate_pointers_simu(SIMU2,Ichild)
  call print_field_vtk(reshape(FC,(/xsize(1)*xsize(2)*xsize(3)/)),"child_prolong.vtk","ID")


  !==============================================================================
  !                                                                TEST REDUCTION
  !                                                                    FC -> FP

  ! Pour realiser une reduction, la simulation courante est ici "2"
  Isimu2 = 2
  
  !                                                          INTERPOLATE
  !                                                            FC -> FCp
  !---------------------------------------------------------------------
  call interp_reduc(Isimu2,"injection",FC)

  !                                                          PACK_BUFFER
  !                                                          FCp -> Sbuf
  !---------------------------------------------------------------------
  call pack_buffer_reduc(Isimu2)
  
  !                                                            ALLTOALLV  
  !                                                         Sbuf -> Rbuf     
  !---------------------------------------------------------------------
  call MPI_ALLTOALLV(Sbuf,AME(Isimu2)%reduc%Scount,AME(Isimu2)%reduc%Sdisp,real_type,&
                     Rbuf,AME(Isimu2)%reduc%Rcount,AME(Isimu2)%reduc%Rdisp,real_type,&
                     MPI_COMM_WORLD,ierror)

  !                                                        UNPACK_BUFFER
  !                                                           Rbuf -> FP
  !---------------------------------------------------------------------
  FP = -2
  call unpack_buffer_reduc(Isimu2,FP)

  !                                                          SORTIES VTK 
  !---------------------------------------------------------------------
  call associate_pointers_simu(SIMU2,Isimu) ! on place les pointers sur la simu "parente"
  call print_field_vtk(reshape(FP,(/xsize(1)*xsize(2)*xsize(3)/)),"parent_reduc.vtk","ID")

end subroutine test_prolong_reduc1

!===============================================================================
!> test_prolong_reduc2 
!!
!! IDEM test_prolon_reduc1 mais avec champs lineaire (pas ct/pinceau)
!!
!!
!! /param[in]  Isimu   Index of the current simulation 
!! /param[in]  child   Child number
!!   
!===============================================================================
subroutine test_prolong_reduc2(Isimu,child)
  implicit none
  integer, intent(in)         :: Isimu, child
  
  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) ::&
                                 FP,FC,X,Y,Z  ! Field Parent, Field Child
  
  integer :: Ichild,ierror,Isimu2
  
  !==============================================================================
  !                                                             TEST PROLONGATION
  !                                                                    FP -> FC
   
  !                             INITIALIZATION FP : lineaire
  !---------------------------------------------------------------------
  call gene_cartesian_coord(X,Y,Z,SIMU2(Isimu)%grid%dx,SIMU2(Isimu)%grid%dy,SIMU2(Isimu)%grid%dz)
  
  FP = X + 2._mytype * Y + 3._mytype * Z
  FC = 0.
  
  !                                                          PACK_BUFFER
  !                                                           FP -> Sbuf
  !---------------------------------------------------------------------
  call pack_buffer_prolong(Isimu,child,FP)
  
  !                                                            ALLTOALLV   
  !                                                         Sbuf -> Rbuf    
  !---------------------------------------------------------------------
  call MPI_ALLTOALLV(Sbuf,AME(Isimu)%prolong(child)%Scount,AME(Isimu)%prolong(child)%Sdisp,real_type,&
                     Rbuf,AME(Isimu)%prolong(child)%Rcount,AME(Isimu)%prolong(child)%Rdisp,real_type,&
                     MPI_COMM_WORLD,ierror)

  !                                                        UNPACK_BUFFER
  !                                                          Rbuf -> FPc
  !---------------------------------------------------------------------
  call unpack_buffer_prolong(Isimu,child)
  
  !                                                          INTERPOLATE
  !                                                           FPc -> FC
  !---------------------------------------------------------------------
  !call interp_prolong(Isimu,child,"injection",FC)
  call interp_prolong(Isimu,child,"linear",FC)
  
  !                                                          SORTIES VTK 
  !---------------------------------------------------------------------
  !   Rq : on doit basculer les pointeurs sur la bonne simu pour avoir 
  !        les bonnes valeurs de dx,dy,dz,x0,y0,z0
  Ichild = AME(Isimu)%children(child)
  call associate_pointers_simu(SIMU2,Isimu)
  call print_field_vtk(reshape(FP,(/xsize(1)*xsize(2)*xsize(3)/)),"parent2.vtk","ID")
  
  call associate_pointers_simu(SIMU2,Ichild)
  call print_field_vtk(reshape(FC,(/xsize(1)*xsize(2)*xsize(3)/)),"child_prolong2.vtk","ID")


  !==============================================================================
  !                                                                TEST REDUCTION
  !                                                                    FC -> FP

  ! Pour realiser une reduction, la simulation courante est ici "2"
  Isimu2 = 2
  
  !                                                          INTERPOLATE
  !                                                            FC -> FCp
  !---------------------------------------------------------------------
  call interp_reduc(Isimu2,"injection",FC)

  !                                                          PACK_BUFFER
  !                                                          FCp -> Sbuf
  !---------------------------------------------------------------------
  call pack_buffer_reduc(Isimu2)
  
  !                                                            ALLTOALLV  
  !                                                         Sbuf -> Rbuf     
  !---------------------------------------------------------------------
  call MPI_ALLTOALLV(Sbuf,AME(Isimu2)%reduc%Scount,AME(Isimu2)%reduc%Sdisp,real_type,&
                     Rbuf,AME(Isimu2)%reduc%Rcount,AME(Isimu2)%reduc%Rdisp,real_type,&
                     MPI_COMM_WORLD,ierror)

  !                                                        UNPACK_BUFFER
  !                                                           Rbuf -> FP
  !---------------------------------------------------------------------
  FP = -2
  call unpack_buffer_reduc(Isimu2,FP)

  !                                                          SORTIES VTK 
  !---------------------------------------------------------------------
  call associate_pointers_simu(SIMU2,Isimu) ! on place les pointers sur la simu "parente"
  call print_field_vtk(reshape(FP,(/xsize(1)*xsize(2)*xsize(3)/)),"parent_reduc2.vtk","ID")

end subroutine test_prolong_reduc2

!===============================================================================
!>  pack_buffer_prolong :  
!!        Allocate and assign Sbuf
!!        Allocate Rbuf
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!! /param[in]  FP      Parent Field (to be prolongated)
!!   
!===============================================================================
subroutine pack_buffer_prolong(Isimu,child,FP)
  implicit none
  integer, intent(in)                               :: Isimu, child
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                                         intent(in) :: FP
  
  integer, dimension(3)       :: kmin,kmax
  integer                     :: i,Npencil_child,rank,count0,disp0,ierror
  logical,dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) &
                              :: mask

  ! ALLOCATION DES BUFFERS
  if (allocated(Rbuf)) deallocate(Rbuf)
  if (allocated(Sbuf)) deallocate(Sbuf)

  allocate(Sbuf(sum(AME(Isimu)%prolong(child)%Scount)))
  allocate(Rbuf(sum(AME(Isimu)%prolong(child)%Rcount)))
  
  ! AFFECTATION BUFFER Sbuf (utilisation de pack)
  Npencil_child = size(AME(Isimu)%prolong(child)%child_data,1)
  mask=.false.
  do i = 1, Npencil_child
    rank = AME(Isimu)%prolong(child)%child_data(i,1)
    kmin = AME(Isimu)%prolong(child)%child_data(i,2:4)
    kmax = AME(Isimu)%prolong(child)%child_data(i,5:7)
    count0 = AME(Isimu)%prolong(child)%Scount(rank+1)
    disp0  = AME(Isimu)%prolong(child)%Sdisp(rank+1)
    !
    mask(kmin(1):kmax(1),kmin(2):kmax(2),kmin(3):kmax(3))=.true.
!print *, "rk ", nrank,"count0/count/disp ", count0,count(mask),disp0

    Sbuf(1+disp0:disp0+count0) = pack(FP,mask)
    mask=.false.
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
end subroutine pack_buffer_prolong  
!===============================================================================
!>  pack_buffer_reduc :  
!!        Allocate and assign Sbuf from FCp
!!        Allocate Rbuf
!!
!! /param[in]  Isimu   Index of the simulation 
!!   
!===============================================================================
subroutine pack_buffer_reduc(Isimu)
  implicit none
  integer, intent(in)                               :: Isimu
  
  integer, dimension(3)       :: kmin,kmax, Imin,Imax
  integer                     :: i,Npencil_child,rank,count0,disp0,ierror

  logical,allocatable, dimension(:,:,:) :: mask


  ! ALLOCATION DES BUFFERS
  if (allocated(Rbuf)) deallocate(Rbuf)
  if (allocated(Sbuf)) deallocate(Sbuf)

  allocate(Sbuf(sum(AME(Isimu)%reduc%Scount)))
  allocate(Rbuf(sum(AME(Isimu)%reduc%Rcount)))

  Imin    = AME(Isimu)%reduc%Imin
  Imax    = AME(Isimu)%reduc%Imax  
  allocate(mask(Imin(1):Imax(1),Imin(2):Imax(2),Imin(3):Imax(3)))
  
  ! AFFECTATION BUFFER Sbuf (utilisation de pack)
  Npencil_child = size(AME(Isimu)%reduc%child_data,1)
  mask=.false.
  do i = 1, Npencil_child
    rank = AME(Isimu)%reduc%child_data(i,1)
    kmin = AME(Isimu)%reduc%child_data(i,2:4)
    kmax = AME(Isimu)%reduc%child_data(i,5:7)
    count0 = AME(Isimu)%reduc%Scount(rank+1)
    disp0  = AME(Isimu)%reduc%Sdisp(rank+1)
    !TODO : definir correctement le mask ou passer par un champ FCp
    mask(kmin(1):kmax(1),kmin(2):kmax(2),kmin(3):kmax(3))=.true.
!print *, "rk ", nrank,"count0/count/disp ", count0,count(mask),disp0

    Sbuf(1+disp0:disp0+count0) = pack(FCp,mask)
    mask=.false.
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
end subroutine pack_buffer_reduc  

!===============================================================================
!>  unpack_buffer_prolong :  
!!        Allocate and assign FPc from Rbuf, just after MPI_ALLTOALLV
!!        FPc : Parent Field on Child pencil
!!        FPc can be interpolated to obtain FF the prolongation of FP
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!!   
!===============================================================================
subroutine unpack_buffer_prolong(Isimu,child)
  implicit none
  integer, intent(in)                               :: Isimu, child

  integer, dimension(3)       :: Imin,Imax,kmin,kmax
  integer                     :: i,Npencil_parent,rank,count0,disp0
  logical,allocatable, dimension(:,:,:) :: mask

  ! ALLOCATION FPc,mask
  if (allocated(FPc)) deallocate(FPc)
  Imin = AME(Isimu)%prolong(child)%Imin
  Imax = AME(Isimu)%prolong(child)%Imax
  allocate(FPc(Imin(1):Imax(1),Imin(2):Imax(2),Imin(3):Imax(3)))
  allocate(mask(Imin(1):Imax(1),Imin(2):Imax(2),Imin(3):Imax(3)))
  FPc = -1.
  mask = .false.
  
  ! AFFECTATION FPc
  Npencil_parent = size(AME(Isimu)%prolong(child)%parent_data,1)
  mask=.false.
  do i = 1, Npencil_parent
    rank = AME(Isimu)%prolong(child)%parent_data(i,1)
    kmin = AME(Isimu)%prolong(child)%parent_data(i,2:4)
    kmax = AME(Isimu)%prolong(child)%parent_data(i,5:7)
    count0 = AME(Isimu)%prolong(child)%Rcount(rank+1)
    disp0  = AME(Isimu)%prolong(child)%Rdisp(rank+1)
    mask(kmin(1):kmax(1),kmin(2):kmax(2),kmin(3):kmax(3))=.true.
    FPc = unpack(Rbuf(1+disp0:disp0+count0),mask,FPc)
    mask=.false. 
  end do

end subroutine unpack_buffer_prolong
!===============================================================================
!>  unpack_buffer_reduc :  
!!        Assign FP from Rbuf, just after MPI_ALLTOALLV
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[out] FP      Parent field with datas assigned from the 
!!                     reduction on the child field      
!!   
!===============================================================================
subroutine unpack_buffer_reduc(Isimu,FP)
  implicit none
  integer, intent(in)                               :: Isimu
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                                        intent(out) :: FP

  integer, dimension(3)       :: Imin,Imax,kmin,kmax
  integer                     :: i,Npencil_parent,rank,count0,disp0
  logical,dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: mask
  
  ! AFFECTATION FP
  Npencil_parent = size(AME(Isimu)%reduc%parent_data,1)
  mask=.false.
  do i = 1, Npencil_parent
    rank = AME(Isimu)%reduc%parent_data(i,1)
    kmin = AME(Isimu)%reduc%parent_data(i,2:4)
    kmax = AME(Isimu)%reduc%parent_data(i,5:7)
    count0 = AME(Isimu)%reduc%Rcount(rank+1)
    disp0  = AME(Isimu)%reduc%Rdisp(rank+1)
    mask(kmin(1):kmax(1),kmin(2):kmax(2),kmin(3):kmax(3))=.true.
    FP = unpack(Rbuf(1+disp0:disp0+count0),mask,FP)
    mask=.false. 
  end do

end subroutine unpack_buffer_reduc
!===============================================================================
!>  interp_prolong :  
!!        Allocate and assign FPc from Rbuf, just after MPI_ALLTOALLV
!!        FPc : Parent Field on Child pencil
!!        FPc can be interpolated to obtain FF the prolongation of FP
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!! /param[in]  mode    "injection" (-1 partout ailleurs)
!!                     "linear"  
!! /param[out] FC      Field interpolated from FPc (
!===============================================================================
subroutine interp_prolong(Isimu,child,mode,FC)
  implicit none
  integer, intent(in)                               :: Isimu, child
  character(len=*), intent(in)                      :: mode
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                                        intent(out) :: FC

  integer, dimension(3)       :: Imin,Imax,I0
  real(mytype), dimension(3)  :: X0c, DXc,X0, DX
  integer                     :: i,j,k,im,jm,km,Ichild,raf
  logical                     :: test


  Imin   = AME(Isimu)%prolong(child)%Imin
  Imax   = AME(Isimu)%prolong(child)%Imax
  Ichild = AME(Isimu)%children(child)
  raf    = AME(Ichild)%refine
  X0c = (/ SIMU2(Ichild)%grid%x0,SIMU2(Ichild)%grid%y0,SIMU2(Ichild)%grid%z0 /) 
  X0  = (/ SIMU2(Isimu)%grid%x0,SIMU2(Isimu)%grid%y0,SIMU2(Isimu)%grid%z0 /)
  DX  = (/ SIMU2(Isimu)%grid%dx,SIMU2(Isimu)%grid%dy,SIMU2(Isimu)%grid%dz /)
  
  ! translation entre les grilles
  I0 = nint((X0c - X0)/DX)
  I0 = I0*raf 
  
  !-----------------------------------TEST INJECTION
  if (mode=="injection") then
  FC = -1.
  do k=Imin(3),Imax(3)
  do j=Imin(2),Imax(2)
  do i=Imin(1),Imax(1)
     im = (i-1)*raf + 1 - I0(1)
     jm = (j-1)*raf + 1 - I0(2)
     km = (k-1)*raf + 1 - I0(3)
     ! verif that (im,jm,km) is on the pencil
     test = im >= xstart(1) .and. im <= xend(1) .and. &
            jm >= xstart(2) .and. jm <= xend(2) .and. &
            km >= xstart(3) .and. km <= xend(3)
     if (test) FC(im,jm,km) = FPc(i,j,k)
  end do
  end do
  end do
  end if
  
end subroutine interp_prolong

!===============================================================================
!>  interp_reduc :  
!!        Allocate and assign FCp before packing and use in MPI_ALLTOALLV
!!        FCp : Child Field on parent pencil
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!! /param[in]  mode    "injection" 
!! /param[in]  FC      Field to reduce in FCp
!===============================================================================
subroutine interp_reduc(Isimu,mode,FC)
  implicit none
  integer, intent(in)                               :: Isimu
  character(len=*) ,intent(in)                      :: mode
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                                        intent(in) :: FC

  integer, dimension(3)       :: Imin,Imax,I0
  real(mytype), dimension(3)  :: X0c, DXc,X0, DX
  integer                     :: i,j,k,im,jm,km,Iparent,raf
  logical                     :: test


  Imin    = AME(Isimu)%reduc%Imin
  Imax    = AME(Isimu)%reduc%Imax
  Iparent = AME(Isimu)%parent
  raf     = AME(Isimu)%refine
  X0c = (/ SIMU2(Isimu)%grid%x0,SIMU2(Isimu)%grid%y0,SIMU2(Isimu)%grid%z0 /) 
  X0  = (/ SIMU2(Iparent)%grid%x0,SIMU2(Iparent)%grid%y0,SIMU2(Iparent)%grid%z0 /)
  DX  = (/ SIMU2(Iparent)%grid%dx,SIMU2(Iparent)%grid%dy,SIMU2(Iparent)%grid%dz /)
  
  
  if (allocated(FCp)) deallocate(FCp)
  allocate(FCp(Imin(1):Imax(1),Imin(2):Imax(2),Imin(3):Imax(3)))

  ! translation entre les grilles
  I0 = nint((X0c - X0)/DX) ! decalage en indixes "gros"
  I0 = I0*raf              ! decalage en indices "fin"
  
  ! On affecte les sections regulières de tableaux
  if (mode=="injection") then
  FCp(Imin(1):Imax(1),Imin(2):Imax(2),Imin(3):Imax(3)) = FC( &
      (Imin(1)-1)*raf+1-I0(1):(Imax(1)-1)*raf+1-I0(1):raf,   &
      (Imin(2)-1)*raf+1-I0(2):(Imax(2)-1)*raf+1-I0(2):raf,   &
      (Imin(3)-1)*raf+1-I0(3):(Imax(3)-1)*raf+1-I0(3):raf)      
  else
  !message d'erreur
  end if
  
end subroutine interp_reduc

!===============================================================================
!>  print_prolong :  
!!        ecriture de la structure PROLONGATION pour (Isimu,child)
!!
!! /param[in]  Isimu   Index of the simulation 
!! /param[in]  child   Child number
!!   
!===============================================================================
subroutine print_prolong(Isimu,child)
  implicit none
  integer, intent(in)                               :: Isimu, child

  integer :: ierror  
  integer :: i,Npencil_child

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  Npencil_child = size(AME(Isimu)%prolong(child)%child_data,1)
  do i = 1,Npencil_child
     !print*, "rk ", nrank, "child_data", AME(Isimu)%prolong(child)%child_data
  end do
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

end subroutine print_prolong
end module amitex_user_mod
