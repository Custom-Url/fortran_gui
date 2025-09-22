!===============================================================================
!>                  PROGRAM OPTIMIZING 2DECOMP - DECOMPOSITION 
!! 

!!==============================================================================
!!
!! UTILISATION DE L'EXECUTABLE :
!!      mpirun -np nbproc optim_2decomp -NX nx -NY ny -NZ -nz
!!
!!==============================================================================
program optim_2decomp

!> CHARGEMENT DES DIFFERENTS MODULES

  use ISO_FORTRAN_ENV

  use MPI             
  use decomp_2d
  use decomp_2d_fft

  use error_mod
  use io_amitex_mod
  use field_mod

!!==============================================================================
!!                                                                  DECLARATIONS
!!==============================================================================
  implicit none

!!------------------------------------------------------------------------------
!>                                                       DIMENSIONS DE LA GRILLE

  integer               :: nx, ny, nz                   !< dimensions de la cellule

!!------------------------------------------------------------------------------
!>                                                                   CONTRAINTES 
!>                                                           tableaux paralleles

  !> dans l'espace reel : tableaux (ntot,Ntens) - pinceaux X
  !>                      notation 1D (1:xsize(1)*xsize(2)*xsize(3),6)
  real(mytype),allocatable,dimension(:,:)         :: Sig

  !> dans l'espace de Fourier : tableaux (nx/2+1,ny,nz,Ntens), ATTENTION AU SENS DE SigF - pinceaux Z
  !>                      notation 3D (fft_start(1):fft_end(1),...(2),...(3),6)
  complex(mytype),allocatable,dimension(:,:,:,:)  :: SigF      


!!------------------------------------------------------------------------------
!>                                                                  PARALLELISME 

  integer :: p_row, p_col,p_row_max,p_col_max, p_row_optim,p_col_optim
                                   !< nb. de lignes, de colonnes pour 2decomp
  integer :: np
  integer :: ierror                !< erreur relative au fonction MPI

  integer :: pmax
  logical :: col_max

!!------------------------------------------------------------------------------
!>                                                          SUIVI DE L'EXECUTION 

  double precision :: t0,t00,tfft,tfftmin,tfftmax,tifftmin,tifftmax, Tmax 
                                        !< double precision au lieu de real_mytype requis par MPI_WTIME
                                        !! t0 : date de depart de l'algorithme
                                        !! t1 : temps pour evaluer la duree de certaines fonctions

!!------------------------------------------------------------------------------
!>                                                                        DIVERS

  integer :: i,j,k                      !< indice de boucle
  integer :: FID                        !< file ID
  integer :: alloc_stat                 !< erreur lors d'une allocation memoire
  integer,dimension(4) :: test_arg      !< variable test
  character(len=200)   :: arg           !< argument de la ligne de commande
  character(len=200)   :: err_msg
  integer              :: Nsim = 1


!!==============================================================================
!!                                                ALLOCATIONS ET INITIALISATIONS
!!
!!              ATTENTION : L'ORDRE DES INITIALISATIONS CI-DESSOUS EST IMPORTANT
!!
!!==============================================================================
!trick to avoid gcc-warning
allocate(Sig(1,1));deallocate(Sig)
allocate(SigF(1,1,1,1));deallocate(SigF)

!!------------------------------------------------------------------------------
!>                                                            INITIALISATION MPI 
!>                                                        + INITIALISATION nrank
!>                                 -> utile tant que 2decomp n'est pas inialisee

  call MPI_INIT(ierror)
  t00 = MPI_WTIME()
  call MPI_COMM_RANK(mpi_comm_world,nrank,ierror)

!!------------------------------------------------------------------------------
!>                                               LECTURE DE LA LIGNE DE COMMANDE 
!>                                                                   -> nx,ny,nz

  test_arg=0
  nx = 0; ny = 0; nz = 0
  j=0
  do i=1,command_argument_count(),2
     j=j+1
     call getarg(i,arg)

     if(arg=="-NX")then
        call getarg(i+1,arg)
        read(arg,"(I8)") nx
        test_arg(1)=1
     else if(arg=="-NY")then
        call getarg(i+1,arg)
        read(arg,"(I8)") ny
        test_arg(2)=1
     else if(arg=="-NZ")then
        call getarg(i+1,arg)
        read(arg,"(I8)") nz
        test_arg(3)=1
     else if(arg=="-Nsim")then
        call getarg(i+1,arg)
        read(arg,"(I8)") Nsim
        test_arg(4)=1
     else
        call write_stdout0("invalid command line : unknown argument "//trim(arg))
        call write_stdout0("unknown argument "//trim(arg))
        call write_stdout0("valid command line : optim_2decomp -NX nx -NY ny -NZ nz (-Nsim nsim)")
        call amitex_abort("invalid command line",2,0)
     end if
  end do

  if ( .not. all(test_arg(1:3) == (/1,1,1/)) ) then
        call write_stdout0("invalid command line : missing arguments ")
        call write_stdout0("valid command line : optim_2decomp -NX nx -NY ny -NZ nz (-Nsim nsim)")
        call amitex_abort("invalid command line",2,0)
  end if 

  if (nrank==0) print *, "nx, ny, nz = ",nx,",", ny,",",nz 
  if (nrank==0) print *, "Nsim = ", Nsim
!!------------------------------------------------------------------------------
!>                                       CALCUL DE LA DECOMPOSITION p_row_,p_col 
!>                                                                   -

  p_row_max = min(nx/2+1,ny)
  p_col_max = min(ny,nz)
  if (nrank==0) print *, "p_row_max, p_col_max = ",p_row_max,",", p_col_max

  call MPI_Comm_size(MPI_COMM_WORLD, np,ierror)  !nombre de processus utilises

  if (np > p_row_max*p_col_max) then
     ! Nombre de processus demande trop eleve 
     write(err_msg,fmt="(A,I0)") "Number of process is too high : the maximum possible is = ", p_row_max*p_col_max
     call amitex_abort(err_msg,2,0)
  end if

  pmax= max(p_col_max,p_row_max)
  col_max=.false.
  if (pmax == p_col_max) col_max=.true.  

TMAX = huge(TMAX)
do k=1,pmax !###################################################################################

  if (col_max) then
     p_col = k
     p_row = np / k
  else
     p_row = k
     p_col = np / k
  end if

if ((p_row * p_col == np) .AND. (p_row .LE. p_row_max) .AND. (p_col .LE. p_col_max)) then;!===========

!!------------------------------------------------------------------------------
!>                                                        INITIALISATION 2DECOMP 
!>                                                                   -
  call decomp_2d_init(nx,ny,nz,p_row,p_col)               !initialisation decomp_2d
  call decomp_2d_fft_init                                 !initialisation decomp_2d_fft
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size) !fft_xxx declares dans le module field_mod.f90

  if (nrank==0) print *, "p_row, p_col = ",p_row,",", p_col

!!------------------------------------------------------------------------------
!>                                          INITIALISATION DES CHAMPS CONTRAINTE
!>                                       
  ! Contrainte
  allocate(Sig(xsize(1)*xsize(2)*xsize(3),1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex13)",2)
     
  ! Contrainte, espace spectral (pinceaux-Z, taille globale ~(nx/2+1)*ny*nz)
  allocate(SigF(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),&
       1:6),stat=alloc_stat)
  if(alloc_stat /=0) call amitex_abort("espace memoire disponible insuffisant (amitex20)",2)

  Sig=0.
  SigF=0.

!!------------------------------------------------------------------------------
!>                                                                           FFT
!>   
! on lance deux series et on ne conserve les temps que de la deuxieme serie
  do j = 1,2
  tfft = 0.
  do i=1,Nsim
    t0 = MPI_WTIME()
    call field_fft(Sig,SigF,6,1)
    tfft = tfft + MPI_WTIME() - t0
  end do
 
  call MPI_allreduce(tfft,tfftmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierror)
  call MPI_allreduce(tfft,tfftmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  end do

        
!!------------------------------------------------------------------------------
!>                                                                          iFFT
!>                                       

  tfft = 0.
  do i=1,Nsim
    t0 = MPI_WTIME()
    call field_ifft(Sig,SigF,6,1)
    tfft = tfft + MPI_WTIME() - t0
  end do
 
  call MPI_allreduce(tfft,tifftmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierror)
  call MPI_allreduce(tfft,tifftmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)

  if (nrank==0) print *, "time fft + ifft (max) = ",tifftmax + tifftmax


!!------------------------------------------------------------------------------
!>                                                                  OPTIMISATION
!>                                       
  
  if (tifftmax+tfftmax < TMAX) then
    TMAX = tifftmax+tfftmax
    p_row_optim = p_row
    p_col_optim = p_col
  end if


  deallocate(Sig)
  deallocate(SigF)
  call decomp_2d_fft_finalize
  call decomp_2d_finalize 

end if !==========================================================================================

end do !##########################################################################################

if (nrank==0) print *, "-------------------------------------------------------" 
if (nrank==0) print *, "" 
if (nrank==0) print *, "" 
if (nrank==0) print *, "         p_row_optim,   p_col_optim = ", p_row_optim, p_col_optim 

!!------------------------------------------------------------------------------
!>                
FID = 7489654                       
if(nrank==0)then
     open(unit=FID, file="user_2decomp.txt",form="formatted", status="replace",action="write")
     write(FID,*) p_row_optim
     write(FID,*) p_col_optim     
     close(FID)
end if

!OPTIMISATION 2decomp
!p_row=0
!p_col=0
!call decomp_2d_init(nx,ny,nz,p_row,p_col)               
!call decomp_2d_finalize 

if (nrank==0) print *, ""
if (nrank==0) print *, ""
if (nrank==0) print *, "elapse time (node 0) = ",MPI_WTIME() - t00

  call MPI_FINALIZE(ierror)

end program optim_2decomp

!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
