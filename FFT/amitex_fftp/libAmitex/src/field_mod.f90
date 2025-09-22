!===========================================================================================================
!
! MODULE FIELD_MOD : 
!
!> Regroupement de subroutines specifiques à la manipulation de champs (tableau 3D
!! a 1 ou plusieurs composantes)
!!
!! champ :  tableau reel (nxP,nyP,nzP,ncomp,nvar) ou (ntotP,ncomp,nvar) ou (ntotP,ncompxnvar) 
!!          tableau complexe (nxP/2,nyP,nzP,ncomp,nvar) ou (ntotPF,ncomp,nvar) ou (ntotPF,ncompxnvar) 
!! Dans Fourier (pinceaux-Z) :     indices (fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3))
!! Dans espace reel (pinceaux-X) : indices (xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))
!!
!! IMPORTANT : beaucoup de fonctions sont indifferentes au profil du tableau          
!! => utilisation du passage d'arguments par adresse et profil explicite pour "reformatter" les donnees
!! 
!! 
!!
!!  Subroutines
!! - add_gradgradU      :       ajoute un champ grad(grad(U)) a un champ de deformation
!! - field_Mean         :       moyenne d'un champ
!! - field_MeanF        :       moyenne d'un champ reel a partir de sa TF
!! - field_isNaN_real   :       verification des NaN sur un champ
!! - field_isNaN_complex:       verification des NaN sur un champ
!! - field_fft          :       fft sur un champ
!! - field_ifft         :       ifft sur un champ
!! - field_sum_squareF  :       somme des carres d'un champ reel partir de sa TF
!! - field_first_order_partial_centeredF  : derivee premiere d'un champ dans une direction (Fourier)
!! - field_second_order_partial_centeredF : derivee seconde
!! - field_filter       :       Filtre "moyenne" dasn l'espace reel
!! - field_filterF      :       Filtre "moyenne" dans l'espace de Fourier
!! - field_mask_outbox  :       Applique un mask a l'exterieur d'une "boite"
!! - laplace_scalarF    :       Laplacien d'un champ scalaire (Fourier)
!! - rotF               :       Rotationnel d'un champ (vecteur, tenseur2 symetrique ou non) (Fourier) 
!! - rot_vectorF        :       Rotationnel d'un champ (vecteur) (Fourier) 
!! - rot_sym_tensor2F   :       Rotationnel d'un champ (tenseur2 symetrique) (Fourier) 
!! - rot_tensor2F       :       Rotationnel d'un champ (tenseur2 non symetrique) (Fourier) 
!! - rotrotF            :       Rot. du rot. d'un champ (vecteur, tenseur2 symetrique ou non) (Fourier) 
!! - rotrot_vectorF     :       Rot. du rot. d'un champ (vecteur) (Fourier) 
!! - rotrot_sym_tensor2F:       Rot. du rot. d'un champ (tenseur2 symetrique) (Fourier) 
!! - rotrot_tensor2F    :       Rot. du rot. d'un champ (tenseur2 non symetrique) (Fourier) 
!! - 
!! - 
!===========================================================================================================
module field_mod

  use ISO_FORTRAN_ENV
  use mpi
  use decomp_2d,      only : mytype, real_type,DECOMP_INFO,&    !TYPES
                             nrank, &!              !variable
                             xstart, xend, xsize, & !pointers variable (since modif decomp_2d for multigrid)
                             update_halo            !functions
  use decomp_2d_fft,  only : decomp_2d_fft_3d
  use error_mod

  private

  !> Pointer "publiques" vers composantes de SIMU_AMITEX
  public :: sp_decomp,ph_decomp,fft_start,fft_end,fft_size,ntotP,ntotFP,&
            times_f


  !> Types publiques (pour definition de SIMU_AMITEX)
  public :: TIMES_FIELD

  !> Fonctions publiques
  public :: field_mean,field_meanF, field_sum_squareF, field_normcomp,  &
            field_isNaN_real,field_isNaN_complex, field_fft, field_ifft, &
            field_first_order_partial_centeredF,field_second_order_partial_centeredF,&
            rot_vectorF,rot_tensor2F,rot_sym_tensor2F,rotF,&
            rotrot_vectorF,rotrot_tensor2F,rotrot_sym_tensor2F,rotrotF,&
            laplace_scalarF, gene_cartesian_coord, add_gradient_scalar_field, add_gradgradU, &
            field_mask_outbox,field_Mean_box, field_filterF,field_filter,field_div,field_divF,&
            field_gradF, field_grad

  public :: var_to_std   ! attention, cette fonction ne travaille pas sur des champs 
                         ! utilisee egalement dans  sortie_std_mod


  type TIMES_FIELD
     double precision  :: fft  = 0  ! temps cumule dans FFT
     double precision  :: ifft = 0  ! temps cumule dans iFFT
  end type TIMES_FIELD
  
  type(DECOMP_INFO), pointer     :: ph_decomp   ! usefull for halo-cell
  type(DECOMP_INFO), pointer     :: sp_decomp        
  integer,dimension(:), pointer  :: fft_start,fft_end,fft_size   
  integer(kind=8),pointer        :: ntotP, ntotFP   
  type(TIMES_FIELD),pointer      :: times_f        
  
  ! besoin de ce parametre pour mettre a zero les variances trop proches du zero machine
  ! et qui soulevent des erreurs dans certains cas (homogene par exemple)
  double precision,parameter                     :: eps=1e-6_8 

  type HALO
    real(mytype),allocatable,dimension(:,:,:) :: Var 
  end type
  
  type(HALO), allocatable,dimension(:) :: halo_tmp
 
contains
!===========================================================================================================
!> field_Mean(field,ntot,ncomp,fieldMean,fieldstd) : MOYENNE D'UN CHAMP & ECART TYPE SI REQUIS
!!
!!   \param[in]  field          champ entree, tableau (nxP,nyP,nzP,ncomp,nvar) (profil indifferent)                           
!!   \param[in]  ntot           ntot, nombre de points du champ complet
!!   \param[in]  ncomp          ncomp, nombre de composantes du champ
!!   \param[out] fieldMean      moyenne du champ, tableau (ncomp,nvar) (profil indifferent)
!!                              
!!   entree optionnelle 
!!   \param[in] nvar            
!!
!!   sortie optionnelle 
!!   \param[out] fieldStd       Ecart-type du champ; tableau (ncomp,nvar) (profil indifferent)
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (field et fieldMoy)
!!
!===========================================================================================================
subroutine field_Mean(field,ntot,ncomp,fieldMean,fieldStd)
  implicit none
!!------------------------------------------------------------------------------
!>                                                                 ENTREE/SORTIE
  integer(KIND=8),intent(in)        :: ntot   
  integer,intent(in)                :: ncomp 
  real(mytype),dimension(xsize(1)*xsize(2)*xsize(3),ncomp),intent(in)           ::      field
  real(mytype),dimension(ncomp),intent(out)                                     ::      fieldMean
  real(mytype),dimension(ncomp),intent(out),optional                            ::      fieldStd

!!------------------------------------------------------------------------------
!>                                                             VARIABLES LOCALES
  integer                                   ::  ierror      !< erreur fonctions MPI
  real(mytype),dimension(ncomp)             ::  fieldSum
  real(mytype),dimension(ncomp)             ::  fieldSquareSum
  integer                                   ::  i

!!------------------------------------------------------------------------------
!!                                                                        CALCUL
  do i=1,ncomp
     fieldSum(i) = sum(field(:,i))
     if (present(fieldStd)) then
        fieldSquareSum(i) = dot_product(field(:,i),field(:,i))
     end if
  end do
  
  call MPI_Allreduce(fieldSum,fieldMean,ncomp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
  fieldMean = fieldMean / ntot

  if (present(fieldStd)) then
    call MPI_Allreduce(fieldSquareSum,fieldStd,ncomp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
    fieldStd = fieldStd / ntot
    do i=1,ncomp
       fieldStd(i) = fieldStd(i) - (fieldMean(i)*fieldMean(i))
    end do
    call var_to_std(fieldMean,fieldStd,ncomp,eps)
  end if

end subroutine field_Mean


!===========================================================================================================
subroutine field_MeanF(fieldF,fieldMean,ntot,ncomp)
  implicit none
!!------------------------------------------------------------------------------
!>                                                                 ENTREE/SORTIE
  integer(KIND=8),intent(in)    :: ntot   
  integer,intent(in)            :: ncomp    
  complex(mytype),dimension(fft_size(1)*fft_size(2)*fft_size(3),ncomp),intent(in)  :: fieldF
  real(mytype),dimension(ncomp),intent(out)                                        :: fieldMean

!!------------------------------------------------------------------------------
!>                                                             VARIABLES LOCALES
  integer                         :: ierror      !< erreur fonctions MPI
  real(mytype),dimension(ncomp)   :: fieldSum
  integer                         :: i

!!------------------------------------------------------------------------------
!!                                                                        CALCUL

  fieldMean=0
  fieldSum=0
  if (fft_start(1)==1 .and. fft_start(2) == 1 .and. fft_start(3)==1 ) then
  do i=1,ncomp
     fieldSum(i) = real(fieldF(1,i),mytype)
  end do
  end if
  
  call MPI_Allreduce(fieldSum(1:ncomp),fieldMean,ncomp,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  
  fieldMean = fieldMean / ntot

end subroutine field_MeanF



!===========================================================================================================
!> field_sum_squareF(fieldF,sum2,nx,ny,nz,ncomp,nvar) : 
!!   
!!   SOMME DES CARRES D'UN CHAMP A PARTIR DE SA TF, COMPOSANTE PAR COMPOSANTE
!!
!!   Definition et calcul dans l'espace de Fourier (Parseval + convention de normalisation fftw): 
!!
!!              sum2 = \sum_{n=0}^{N} |x[n]|^2  
!!
!!                      = (\sum_{k=0}^{N} |\hat{x}[k]|^2)/Ntot 
!!
!!
!!   \param[in]  fieldF         champ a normer (Fourier), 
!!                              tableau complexe (nxP/2,nyP,nzP,ncomp,nvar) (profil indifferent)
!!   \param[in]  nx,ny,nz       nombre de points dans chaque direction de l'espace REEL
!!   \param[in]  ncomp          nombre de composantes du champ
!!   \param[in]  nvar           nombre de champs
!!   \param[out] sum2           pour chaque composante, somme des carre du champ
!!                              tableau (ncomp,nvar) (profil indifferent)
!!
!! ATTENTION :  la transformee de Fourier est decrite sur une grille (nx/2,ny,nz)
!!              il faut doubler la somme sur toutes les frequences 
!!              et soustraire les contributions (1,:,:) et ((nx/2+1,:,:) dans le cas nx pair)
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF,sum2)
!!
!!
!===========================================================================================================
subroutine field_sum_squareF(fieldF,sum2,nx,ny,nz,ncomp,nvar)
  implicit none
!!------------------------------------------------------------------------------
!>                                                                 ENTREE/SORTIE
  complex(mytype),&
     dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
     intent(in)                                         :: fieldF
  real(mytype),dimension(ncomp,nvar),intent(out)        :: sum2
  integer,intent(in)                                    :: nx,ny,nz
  integer,intent(in)                                    :: ncomp,nvar   

!!------------------------------------------------------------------------------
!>                                                             VARIABLES LOCALES
  integer                       :: ierror      !< erreur fonctions MPI
  integer                       :: i,j
  real(mytype)                  :: tmp,tmp2
  integer(KIND=8)               :: ntot   

!!------------------------------------------------------------------------------
!>                                                               INITIALISATIONS

  ntot = (ny * nz)
  ntot = ntot * nx
  
!!------------------------------------------------------------------------------
!!                                                    BOUCLE SUR LES COMPOSANTES
  do j=1,nvar
  do i=1,ncomp

    ! 2 x somme des carres sur le 1/2 espace de Fourier
    tmp = 2._mytype*sum(real(fieldF(:,:,:,i,j)*conjg(fieldF(:,:,:,i,j)),mytype))
    call MPI_AllReduce(tmp,sum2(i,j),1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)

    ! correction des plans (1,:,:)
    tmp=0
    tmp2=0
    if(fft_start(1)==1) then
      tmp = sum(real(fieldF(1,:,:,i,j)*conjg(fieldF(1,:,:,i,j)),mytype))
    end if
    call MPI_AllReduce(tmp,tmp2,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    sum2(i,j) = sum2(i,j) - tmp2

    ! correction des plans (nx/2+1,:,:)
    tmp=0
    tmp2=0
    if(modulo(nx,2) == 0 .AND. fft_end(1)==nx/2+1) then
      tmp = sum(real(fieldF(nx/2+1,:,:,i,j)*conjg(fieldF(nx/2+1,:,:,i,j)),mytype))
    end if
    call MPI_AllReduce(tmp,tmp2,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    sum2(i,j) = sum2(i,j) - tmp2

  end do
  end do

  ! Normalisation fft
  sum2 = sum2 / ntot

end subroutine field_sum_squareF


!===========================================================================================================
!> function field_isNaN(field,nval,test) 
!!
!!   \param[in]  field         champ a moyenner, tableau (nxP,nyP,nzP,ncomp,nvar) (profil indifferent)
!!   \param[in]  nval          nombre de valeur dans le champ (ntotPxncompxnvar)
!!   \param[out] test          logique .true. si presence de NaN (.false. sinon)   
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (field)
!!
!===========================================================================================================
subroutine field_isNaN_real(field,nval,test)
  implicit none
  integer(kind=INT64),intent(in)                :: nval
  real(mytype),dimension(nval),intent(in)       :: field
  logical,intent(out)                           :: test
  integer(kind=INT64)                           :: i

  test=.false.
  do i=1,nval
     test = test .or. isnan(field(i)) 
  end do
end subroutine field_isNaN_real
!-----------------------------------------------------------------------------------------------------------
subroutine field_isNaN_complex(field,nval,test)
  implicit none
  integer(kind=INT64),intent(in)                :: nval
  complex(mytype),dimension(nval),intent(in)    :: field
  logical,intent(out)                           :: test
  integer(kind=INT64)                           :: i

  test=.false.
  do i=1,nval
     test = test .or. isnan(real(field(i))) 
     test = test .or. isnan(aimag(field(i))) 
  end do
end subroutine field_isNaN_complex

!!==========================================================================================================
!> fft(field,fieldF,ncomp,nvar) : fft d'un champ
!!
!!   \param[in]  field          champ entree (nxP,nyP,nzP,ncomp,nvar)  (profil indifferent)
!!   \param[out] fieldF         champ sortie (nxP/2,nyP,nzP,ncomp,nvar)(profil indifferent)
!!   \param[in]  ncomp          nombre de composantes du champ
!!   \param[in]  nvar           nombre de champs (1 en mecanique)
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (field,fieldF)
!!
!===========================================================================================================
subroutine field_fft(field,fieldF,ncomp,nvar)
  
  implicit none
  ! on utilise ici le passage par adresse pour reshaper field et fieldF
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncomp,nvar),&
               intent(inout)    :: field
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
               intent(out)      :: fieldF
  integer,intent(in)            :: ncomp,nvar
  integer                       :: i,j,ierror
  double precision              :: t1

  call mpi_barrier(mpi_comm_world,ierror)          
  t1 = MPI_WTIME()

  do j = 1,nvar
  do i = 1,ncomp
     call decomp_2d_fft_3d(field(:,:,:,i,j),fieldF(:,:,:,i,j))
  end do
  end do

  call mpi_barrier(mpi_comm_world,ierror)            
  times_f%fft =  times_f%fft + MPI_WTIME() - t1

end subroutine field_fft

!!==========================================================================================================
!> fft(field,fieldF,ncomp,nvar) : fft d'un champ
!!
!!   \param[in]  field          champ entree (nxP,nyP,nzP,ncomp,nvar)  (profil indifferent)
!!   \param[out] fieldF         champ sortie (nxP/2,nyP,nzP,ncomp,nvar)(profil indifferent)
!!   \param[in]  ncomp          nombre de composantes du champ
!!   \param[in]  nvar           nombre de champs (1 en mecanique)
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (field,fieldF)
!!
!===========================================================================================================
subroutine field_ifft(field,fieldF,ncomp,nvar)
  
  implicit none
  ! on utilise ici le passage par adresse pour reshaper field et fieldF
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncomp,nvar),&
               intent(out)      :: field
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
               intent(inout)    :: fieldF
  integer,intent(in)            :: ncomp,nvar
  integer                       :: i,j,ierror
  double precision              :: t1

  call mpi_barrier(mpi_comm_world,ierror)          
  t1 = MPI_WTIME()

  do j = 1,nvar
  do i = 1,ncomp
     call decomp_2d_fft_3d(fieldF(:,:,:,i,j),field(:,:,:,i,j))
  end do
  end do

  call mpi_barrier(mpi_comm_world,ierror)            
  times_f%ifft =  times_f%ifft + MPI_WTIME() - t1

end subroutine field_ifft

!!==========================================================================================================
!> field_normcomp(field,norm_field,ncomp,nvar) : norme euclidienne d'un champ a plusieurs composantes
!>                                     racine carree(somme des composantes au carre)
!!
!!   \param[in]  field          champ entree (nxP,nyP,nzP,ncomp,nvar) (profil indifferent)
!!   \param[out] norm_field     champ sortie (nxP,nyP,nzP,nvar)       (profil indifferent)
!!   \param[in]  ncomp          nombre de composantes du champ
!!   \param[in]  nvar           nombre de champs (1 en mecanique)
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (field,norm_field)
!!
!===========================================================================================================
subroutine field_normcomp(field,norm_field,ncomp,nvar)
  
  implicit none
  ! on utilise ici le passage par adresse pour reshaper field et norm_field
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncomp,nvar),&
               intent(in)      :: field
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),nvar),&
               intent(inout)    :: norm_field
  integer,intent(in)            :: ncomp,nvar
  integer :: i,j,ix,iy,iz
  
  norm_field(:,:,:,:) = 0._mytype

  do iz = xstart(3),xend(3)
  do iy = xstart(2),xend(2)
  do ix = xstart(1),xend(1)
      do j = 1,nvar
          do i = 1,ncomp
             norm_field(ix,iy,iz,j) = norm_field(ix,iy,iz,j) + field(ix,iy,iz,i,j)*field(ix,iy,iz,i,j)
          end do
             norm_field(ix,iy,iz,j) = sqrt(norm_field(ix,iy,iz,j))
      end do
  end do 
  end do
  end do  

end subroutine field_normcomp
!==============================================================================
!       SUBROUTINE VAR_TO_STD
!------------------------------------------------------------------------------
!> Evalue l'ecart-type a partir de la variance et de la moyenne 
!!           en annulant les éventuels arrondis négatifs
!!
!! \param[in]   ncomposantes   entier   
!!              Moyenne        tableau real (ncomposantes)
!!              EcartType      variance
!!                             tableau real (ncomposantes)
!!                 
!! \param[out]  EcartType      Ecart-Type
!!                             tableau real (ncomposantes) 

!! La variance est théoriquement toujours positive. Néanmoins lorsqu'elle est 
!! très petite par rapport à la moyenne, une erreur d'arrondi peut se produire
!! et résulter en une variance négative. 
!! Lorsque ces deux conditions sont réunies (variance très petite devant
!! la moyenne, et variance négative), la variance est mise à zéro.
!! 
!! Une variance négative et non petite devant la moyenne détectée entraine
!! un message d'erreur et un arrêt du code.
!==============================================================================
subroutine var_to_std(Moyenne,EcartType,ncomposantes,seuil)

  integer,intent(in)                                :: ncomposantes
  real(mytype),dimension(ncomposantes),intent(in)   :: Moyenne
  real(mytype),dimension(ncomposantes),intent(inout):: EcartType
  real(mytype),dimension(ncomposantes)              :: PosErreur
  real(mytype),intent(in)                           :: seuil
  character(96)                                     :: tmp_char0,tmp_char1 ! 96 = 6x16
  

  where(Moyenne .ne. 0) 
     PosErreur(:)=abs(EcartType/(Moyenne*Moyenne))
  elsewhere
     PosErreur(:)=seuil+1.
  end where
  where(PosErreur(:)<seuil )
     EcartType=0
  end where

  if (minval(EcartType) < 0) then
     write (tmp_char0,"(6E15.8)") EcartType
     write (tmp_char1,"(6E15.8)") Moyenne * Moyenne 
     call amitex_abort("Negative variance (field_mod, var_to_std) - variance and mean^2 : "&
                       //achar(10)//tmp_char0//achar(10)//tmp_char1,2,0)
  end if

  ! ecart type
  EcartType=sqrt(EcartType)

end subroutine var_to_std

!==============================================================================
!       SUBROUTINE FIELD_FIRST_ORDER_PARTIAL_CENTEREDF
!------------------------------------------------------------------------------
!> Evalue la derivee partielle de premier ordre dans l'espace de Fourier du 
!! champ en entree pour chacune de ses composantes
!!
!!              D_FieldF = D_fieldF + coeff x dfieldF/dxi
!!
!! schema de derivation : DIFFERENCE FINIE CENTREE
!!
!!   \param[inout] D_fieldF    champ entree-sortie : D_fieldF(entree) + coeff x derivee partielle du champ fieldF 
!!                             tableau complexe (nxP/2,nyP,nzP,ncomp,nvar) (profil indifferent)
!!
!!   \param[in] coeff          coefficient multiplicateur (reel)
!!
!!   \param[in] fieldF         champ entree, (nxP/2,nyP,nzP,ncomp,nvar)    (profil indifferent)
!!
!!   \param[in]  dir           dimension selon laquelle on derive 
!!                             (direction de derivation)
!!
!!   \param[in]  ncomp         nombre de composantes du champ 
!!   \param[in]  nvar          nombre de champs (1 en mecanique)
!!
!!   \param[in]  N             dimension de la cellule dans la direction dir
!!   \param[in]  d             taille des voxels dans la direction dir
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!!
!!   Derivation du premier ordre de chaque composante scalaire du champ correspondant a :
!!   D_fieldF(compi) = D_fieldF(compi) + coeff x dfieldF(compi)/dxdir dans l'espace de fourier
!!
!!   exemple : si dir = 1 (x1 = x)  D_fieldF(:) = D_fieldF(:) + coeff x dfieldF(:)/dx
!==============================================================================
subroutine field_first_order_partial_centeredF(D_fieldF,coeff,fieldF,dir,ncomp,nvar,N,d)

  implicit none
  
  !Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
                intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
               intent(inout)    :: D_fieldF 
  integer,intent(in)          :: dir,ncomp,nvar,N
  real(mytype),intent(in)     :: d, coeff
  
  ! Autre variables
  integer(kind=8)             :: i,imin=0,imax=0 ! variables de boucles
#ifdef DOUBLE_PREC
  real(mytype),parameter      :: PI = 4._mytype*DATAN(1._mytype)
#else
  real(mytype),parameter      :: PI = 4._mytype*ATAN(1._mytype)
#endif
  complex(mytype)             :: Frequence
  complex(mytype)             :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i


  ! Determination des parametres numeriques necessaires au calcul de 
  ! la frequence liee a l'approximation differences finies
  select case(dir)
    case(1)
        imin = fft_start(1)
        imax = fft_end(1)
    case(2)
        imin = fft_start(2)
        imax = fft_end(2)
    case(3)
        imin = fft_start(3)
        imax = fft_end(3)
    case default
       call amitex_abort("indices de la direction de derivation non compris entre 1 et 3&
            & (field_mod, field_first_order_partial_centeredF)",2,0) 
  end select

  ! calcul de la derivee dans l'espace de Fourier : 
  do i=imin,imax
     
     ! Calcul de la frequence liee a l'approximation
     ! differences finies
     Frequence = (1/d)*sin(2*PI*(i-1)/N)

     select case(dir)
       case(1)
           D_fieldF(i,:,:,:,:) = D_fieldF(i,:,:,:,:) + coeff*imP*Frequence*fieldF(i,:,:,:,:) 
       case(2)
           D_fieldF(:,i,:,:,:) = D_fieldF(:,i,:,:,:) + coeff*imP*Frequence*fieldF(:,i,:,:,:)
       case(3)
           D_fieldF(:,:,i,:,:) = D_fieldF(:,:,i,:,:) + coeff*imP*Frequence*fieldF(:,:,i,:,:)
       case default
          call amitex_abort("indices de la direction de derivation non compris entre 1 et 3&
               & (field_mod, field_first_order_partial_centeredF)",2,0) 
     end select
  end do

end subroutine field_first_order_partial_centeredF

!==============================================================================
!       SUBROUTINE FIELD_SECOND_ORDER_PARTIAL_CENTEREDF
!------------------------------------------------------------------------------
!> Evalue la derivee partielle de second ordre dans l'espace de Fourier du 
!! champ en entree pour chacune de ses composantes
!!
!!              D_FieldF = D_fieldF + coeff x d2fieldF/dxidxj
!!
!! schema de derivation : DIFFERENCE FINIE CENTREE "approximation 27 voxels"
!!
!!   \param[inout] D_fieldF     champ entree/sortie : D_fieldF(entree) + coeff x derivee partielle du champ fieldF
!!                              tableau complexe (nxP/2,nyP,nzP,ncomp,nvar) (profil indifferent)
!!
!!   \param[in] coeff           coefficient multiplicateur (reel)
!!
!!   \param[in] fieldF          champ entree 
!!                              tableau complexe (nxP/2,nyP,nzP,ncomp,nvar) (profil indifferent)
!!
!!
!!   \param[in]  direction      dimensions selon lesquelles on derive
!!                              vecteur de 2 entiers  (direction de derivation)
!!
!!   \param[in]  ncomp          nombre de composantes du champ 
!!   \param[in]  nvar           nombre de champs (1 en mecanique)
!!
!!   \param[in]  N1             dimension de la cellule dans la direction(1)
!!   \param[in]  N2             dimension de la cellule dans la direction(2)
!!
!!   \param[in]  d1             dimension des voxels dans la direction(1)
!!   \param[in]  d2             dimension des voxels dans la direction(2)
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!!
!!   Derivation du second ordre de chaque composante scalaire du champ correspondant a :
!!   D_fieldF(compi) = D_fieldF(compi) + coeff x d2fieldF(compi)/dxdir(1)dxdir(2) dans l'espace de fourier
!!
!!   exemple : si direction = (1 2)   D_fieldF(:) = D_fieldF + coeff x d2fieldF(:)/dx1dx2 
!========================================================================================
subroutine field_second_order_partial_centeredF(D_fieldF,coeff,fieldF,direction,ncomp,nvar,N1,N2,d1,d2)

  implicit none
  
  !Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
                            intent(inout)  :: D_fieldF 
  real(mytype),intent(in)                  :: coeff
  integer,intent(in),dimension(2)          :: direction
  integer,intent(in)                       :: ncomp,nvar,N1,N2
  real(mytype),intent(in)                  :: d1,d2
  
  ! Autre variables
  integer                     :: i,j,imin=0,imax=0,jmin=0,jmax=0    ! variables de boucles initialisees ici pour
                                                                     ! eviter warning gcc
#ifdef DOUBLE_PREC
  real(mytype),parameter      :: PI = 4._mytype*DATAN(1._mytype)
#else
  real(mytype),parameter      :: PI = 4._mytype*ATAN(1._mytype)
#endif
  complex(mytype)             :: Frequence
  integer,dimension(2)        :: dir
  integer                     :: NN1,NN2
  real(mytype)                :: dd1,dd2


  ! Mise dans l'ordre croissant des directions de derivation
  dir(1) = minval(direction)
  dir(2) = maxval(direction)

  ! Determination des parametres numeriques necessaires au calcul de 
  ! la frequence liee a l'approximation differences finies

  if (dir(1) == dir(2)) then
         ! Derivee seconde dans la meme dimension d²/dxdim(1)²
          
         ! dimension de la cellule et taille des voxels dans la direction de derivation (1=2)
          NN1  = N1
          dd1  = d1

    select case(dir(1))
      case(1)
          ! d²/dx² 
          imin = fft_start(1)
          imax = fft_end(1)
          jmin = 1
          jmax = 1
      case(2)
          ! d²/dy² 
          imin = fft_start(2)
          imax = fft_end(2)
          jmin = 1
          jmax = 1
      case(3)
          ! d²/dz² 
          imin = fft_start(3)
          imax = fft_end(3)
          jmin = 1
          jmax = 1
      case default
          call amitex_abort("indices des directions de derivation non compris entre 1 et 3&
               & (field_mod, field_second_order_partial_centeredF)",2,0)           
    end select
  else
         ! Derivee croisee par rapport a deux dimensions d²/dxdim(1)dxdim(2)
    ! dimension de la cellule et taille des voxels dans les dir. de derivation 1 et 2    
    if (direction(1) < direction(2)) then
      NN1  = N1
      dd1  = d1
      NN2  = N2
      dd2  = d2        
    else
      NN1  = N2
      dd1  = d2
      NN2  = N1
      dd2  = d1  
    end if
    
    select case(dir(1))
      case(1)
          imin = fft_start(1)
          imax = fft_end(1)

          select case(dir(2))
             case(2)
              ! d²/dxdy 
               jmin = fft_start(2)
               jmax = fft_end(2)
             case(3)
              ! d²/dxdz 
               jmin = fft_start(3)
               jmax = fft_end(3)
             case default
               call amitex_abort("indices des directions de derivation non compris entre 1 et 3&
                    &  (field_mod, field_second_order_partial_centeredF)",2,0)       
          end select
      case(2)
          imin = fft_start(2)
          imax = fft_end(2)

          select case(dir(2))
             case(3)
              ! d²/dydz
               jmin = fft_start(3)
               jmax = fft_end(3)
             case default
               call amitex_abort("indices des directions de derivation non compris entre 1 et 3&
                    &  (field_mod, field_second_order_partial_centeredF)",2,0)  
          end select
       ! 3 en premier n'est pas possible dans le cas ou dir(1) et dir(2) sont differents 
       ! car mise dans l'ordre croissant au depart
    end select

  end if


  ! calcul de la derivee dans l'espace de Fourier : 
  do i=imin,imax
    do j=jmin,jmax
      
       if (dir(1) == dir(2)) then
         ! Derivee seconde dans la meme dimension d²/dxdim(1)²

         Frequence = ( 2/(dd1*dd1) )*(cos( 2*PI*(i-1)/NN1 ) - 1)

         select case(dir(1))
           case(1)
               D_fieldF(i,:,:,:,:) = D_fieldF(i,:,:,:,:) + coeff * Frequence*fieldF(i,:,:,:,:) 
           case(2)
               D_fieldF(:,i,:,:,:) = D_fieldF(:,i,:,:,:) + coeff * Frequence*fieldF(:,i,:,:,:)
           case(3)
               D_fieldF(:,:,i,:,:) = D_fieldF(:,:,i,:,:) + coeff * Frequence*fieldF(:,:,i,:,:)
         end select
       else
         ! Derivee croisee par rapport a deux dimensions d²/dxdim(1)dxdim(2)

            ! oblige de mettre real(,mytype) car les termes entre parenthese dans les cos sont consideres comme des entiers 
            ! par l'operateur / qui renvoie donc 0 (quotient < 1)
           Frequence = ( 1/(2*dd1*dd2) )* &
                  ( cos( 2*PI*( real(i-1,mytype)/NN1 + real(j-1)/NN2 ) ) &
                  - cos( 2*PI*( real(i-1,mytype)/NN1 - real(j-1,mytype)/NN2 ) ) ) 

           if( (dir(1)==1 .and. dir(2)==2) ) then
                 ! derivee croisee 1/2 : d²/dxdy

                 D_fieldF(i,j,:,:,:) = D_fieldF(i,j,:,:,:) + coeff * Frequence*fieldF(i,j,:,:,:) 
            elseif ( (dir(1)==1 .and. dir(2)==3) ) then
                 ! derivee croisee 1/3 : d²/dxdz

                 D_fieldF(i,:,j,:,:) = D_fieldF(i,:,j,:,:) + coeff * Frequence*fieldF(i,:,j,:,:)
            elseif ( (dir(1)==2 .and. dir(2)==3)) then
                 ! derivee croisee 2/3 : d²/dydz

                 D_fieldF(:,i,j,:,:) = D_fieldF(:,i,j,:,:) + coeff * Frequence*fieldF(:,i,j,:,:)
           end if
       end if
    end do
  end do

end subroutine field_second_order_partial_centeredF

!==============================================================================
!       SUBROUTINE LAPLACE_scalarF
!------------------------------------------------------------------------------
!> Evalue le laplacien d'un champ scalaire en entree dans l'espace de Fourier
!!
!! Si les frequences FREQ sont passees en argument :
!!             On les utilise pour evaluer les derivees secondes
!! Sinon :
!!             On utilise l'operateur de derivee centree 
!!
!!
!!   \param[in] fieldF  Entree 
!!                      tableau complexe (nxP/2,nyP,nzP) (profil indifferent)
!!
!!   \param[out] Laplace_fieldF laplacien du champ entree 
!!                              tableau complexe (nxP/2,nyP,nzP) (profil indifferent)
!!
!!   \param[in] FREQ  tableau de fréquence, espace Fourier 
!!                    pinceaux-Z (nx/2+1,ny,nz,3) 
!!                    optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et Laplace_fieldF)
!!
!
!  formule :
!  --------
!
!  laplace(S) = S,11 + S,22 + S,33
!
!==============================================================================
subroutine laplace_scalarF(fieldF,Laplace_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1),&
                            intent(out)    :: Laplace_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: r,s,q        ! variables de boucle calcul derivee

  Laplace_fieldF(:,:,:,1) = 0._mytype

  if (present(FREQ)) then
     ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
     ! l'operateur d'ordre pour les derivees d'ordre 2


     do r=fft_start(3),fft_end(3)
        do s = fft_start(2),fft_end(2)
           do q = fft_start(1),fft_end(1)
              Laplace_fieldF(q,s,r,1) = Laplace_fieldF(q,s,r,1) -&
              (FREQ(q,s,r,1)*FREQ(q,s,r,1)+FREQ(q,s,r,2)*FREQ(q,s,r,2)+FREQ(q,s,r,3)*FREQ(q,s,r,3))*fieldF(q,s,r,1)
           end do
        end do
     end do

  else 
     ! Calcul des derivees a partir des operateurs differences finies sur cube de 27 voxels
     ! l'operateur field_second_order_partial_centeredF somme les contributions successives
     ! des derivees partielles (11, 22, 33) dans Laplace_fieldF
     call field_second_order_partial_centeredF(&
             Laplace_fieldF(:,:,:,1),1._mytype,fieldF(:,:,:,1),(/ 1,1 /),1,1,N(1),N(1),d(1),d(1))
     call field_second_order_partial_centeredF(&
             Laplace_fieldF(:,:,:,1),1._mytype,fieldF(:,:,:,1),(/ 2,2 /),1,1,N(2),N(2),d(2),d(2))
     call field_second_order_partial_centeredF(&
             Laplace_fieldF(:,:,:,1),1._mytype,fieldF(:,:,:,1),(/ 3,3 /),1,1,N(3),N(3),d(3),d(3))
  end if

end subroutine laplace_scalarF

!==============================================================================
!       SUBROUTINE ROTF
!------------------------------------------------------------------------------
!> 
!> Evalue ROT(xx) pour differents type de xx : 
!>                 vecteur, tenseur d'ordre 2 sym. et non sym. 
!>
!>  Verifie la taille du champ en entree et le renvoie vers la fonction rotationnel
!>  correspondante
!>
!!   \param[in] fieldF         	champ entree, tableau complexe (nxP/2,nyP,nzP,ncomp)
!!   \param[out] Rot_fieldF   	champ sortie, tableau complexe (nxP/2,nyP,nzP,ncomprot)
!!   \param[in] N               dimensions de la cellule [Nx Ny Nz] 
!!   \param[in] d               dimension des voxels [dx,dy,dz]
!!   \param[in] FREQ		frequences (nxP/2,nyP,nzP,3)
!!
!!
!! ATTENTION : si ncomp = 3,6,9 alors ncomprot = 3,9,9
!!           => Utilisation de tableaux a profil implicite
!!
!! TODO : A REVOIR pour passage tableau avec profil indifferent
!!        passage ncomp, nvar, ncomprot
!!
!==============================================================================
subroutine rotF(fieldF,Rot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(:,:,:,:),intent(in)     :: fieldF 
  complex(mytype),dimension(:,:,:,:),intent(out)    :: Rot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Variables locales
  integer                    :: size_tens, size_grad

  size_tens = size(fieldF,4)
  size_grad = size(Rot_fieldF,4)

  select case(size_tens)
      case(3) ! rot(vecteur)
         if (size_grad == 3) then
            call rot_vectorF(fieldF,Rot_fieldF,N,d,FREQ)
         else
            call amitex_abort("Taille du champ de sortie incorrecte, attendu=(:,:,:,3) (rotrot)" ,2,0)  
         end if
      case(6) ! rot(tenseur_sym)
         if (size_grad == 9) then
            call rot_sym_tensor2F(fieldF,Rot_fieldF,N,d,FREQ)
         else
            call amitex_abort("Taille du champ de sortie incorrecte, attendu=(:,:,:,9) (rotrot)" ,2,0)  
         end if
      case(9) ! rot(tenseur)
         if (size_grad == 9) then
            call rot_tensor2F(fieldF,Rot_fieldF,N,d,FREQ)
         else
            call amitex_abort("Taille du champ de sortie incorrecte, attendu=(:,:,:,9) (rotrot)" ,2,0)  
         end if   
      case default   
            call amitex_abort("Taille du champ en entree incorrecte, rot ne peut s'appliquer&
                 & que a des champs de taille 3 (vecteurs) ,6 (tenseurs sym)et 9 (tenseurs)",2,0)     
   end select

end subroutine rotF

!==============================================================================
!       SUBROUTINE ROT_VECTORF
!------------------------------------------------------------------------------
!> Evalue le rot() du champ de vecteurs en entree dans l'espace de Fourier
!!
!!   \param[in] fieldF         champ de vecteur entree 
!!                             tableau complexe (nxP/2,nyP,nzP,3) (profil indifferent)
!!
!!   \param[out] D_fieldF      rotationnel du champ entree : derivee partielle du champ entree 
!!                             tableau complexe(nxP/2,nyP,nzP,3) (profil indifferent)
!!
!!   \param[in] FREQ  tableau de fréquence, espace Fourier 
!!                   pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                    optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et Rot_fieldF)
!!
!
!  formule :
!  --------
!
!  rot(U)(k) = perm(i,j,k) * U(i),j
!!
!! Peut etre calcule de deux façons : 
!!      A) en utilisant les operateurs de derivation centree definis dans le field_mod
!          field_(first/second)_order_partial_centeredF
!          equivalent a une approximation differences finies sur un cube de 3 voxels
!          de cote
!
!       B) en utilisant les frequences initialisees dans le green_mod ce qui equivaut
!          a appliquer deux fois l'operateur de derivation du premier ordre 
!          correspondant au filtre choisi dans algo_param. 
!          APPLIQUE SI LE TABLEAU FREQ EST PASSE EN ARGUMENT
!
!==============================================================================
subroutine rot_vectorF(fieldF,Rot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3),&
                            intent(out)    :: Rot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: k,i,j        ! variables de boucles formule rotationnel
  integer                                  :: r,s,q        ! variables de boucle calcul derivee
  real(mytype)                             :: perm2
  complex(mytype)                          :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i

  Rot_fieldF(:,:,:,:) = 0.

  do k=1,3
  ! Calcul composante k  du rotationnel
     do j=1,3
        do i=1,3
             ! calcul de permutation_tensor(i,j,k)*fieldF(i),j


             perm2 = real(permutation_tensor(i,j,k),mytype)
            if (perm2 .ne. 0) then
              
               if (present(FREQ)) then
                 ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
                 ! l'operateur d'ordre pour les derivees d'ordre 2

                  
                 do r=fft_start(3),fft_end(3)
                   do s = fft_start(2),fft_end(2)
                      do q = fft_start(1),fft_end(1)
                          Rot_fieldF(q,s,r,k) = Rot_fieldF(q,s,r,k) + imP*perm2*FREQ(q,s,r,j)*fieldF(q,s,r,i)
                      end do
                   end do
                 end do

               else 
                 ! Calcul des derivees a partir des operateurs differences finies sur cube de 9 voxels
                 ! l'operateur field_second_order_partial_centeredF somme les contributions successives
                 ! des derivees partielles (,mj) dans Laplace_fieldF                   
                 call field_first_order_partial_centeredF(Rot_fieldF(:,:,:,k),&
                                   perm2,fieldF(:,:,:,i),j,1,1,N(j),d(j))

              end if
             end if

              ! fin calcul de la derivee seconde permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * dfieldF(l)/dXmdXj
        end do
     end do
  ! Fin calcul composante k
  end do  

end subroutine rot_vectorF

!==============================================================================
!       SUBROUTINE ROT_SYM_TENSOR2F
!------------------------------------------------------------------------------
!> Evalue le rot() du champ de tenseur symétrique d'ordre 2 en entree, dans l'espace de Fourier
!!
!!   \param[in] fieldF         champ de tenseur ordre 2 sym entree 
!!                             tableau complexe (nxP/2,nyP,nzP,6) (profil indifferent)
!!
!!   \param[out] D_fieldF      rotationnel du champ entree : derivee partielle du champ entree 
!!                             tableau complexe (nxP/2,nyP,nzP,9) (profil indifferent)
!!
!!   \param[in] FREQ           tableau de fréquence, espace Fourier 
!!                             pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                             optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et Rot_fieldF)
!!
!
!  formule :
!  --------
!
!  rot(T)(i,j) = perm(j,m,k) * T(i,m),k
!!
!! Peut etre calcule de deux façons : 
!!      A) en utilisant les operateurs de derivation centree definis dans le field_mod
!          field_(first/second)_order_partial_centeredF
!          equivalent a une approximation differences finies sur un cube de 3 voxels
!          de cote
!
!       B) en utilisant les frequences initialisees dans le green_mod ce qui equivaut
!          a appliquer deux fois l'operateur de derivation du premier ordre 
!          correspondant au filtre choisi dans algo_param. 
!          APPLIQUE SI LE TABLEAU FREQ EST PASSE EN ARGUMENT
!
!==============================================================================
subroutine rot_sym_tensor2F(fieldF,Rot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),6),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),&
                            intent(out)    :: Rot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: k,i,j,m  ! variables de boucles formule rotationnel
  integer                                  :: r,s,q        ! variables de boucle calcul derivee
  real(mytype)                             :: perm2
  complex(mytype)                          :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i

  Rot_fieldF(:,:,:,:) = 0._mytype

  do i=1,3
  do j=1,3
  ! Calcul composante (i,j) du rotationnel
     do k=1,3
       do m=1,3
         ! calcul de perm(j,m,k) * T(i,m),k
         ! on fait le lien entre la composante (i,j) / (i,m) des champs en ecriture conventionnelle
         ! et leur indice dans les Rot_field et fieldF a l'aide des fonctions indice_lin_st2 et indice_lin_t2

         perm2 = real(permutation_tensor(j,m,k),mytype)
         if (perm2 .ne. 0) then
          
           if (present(FREQ)) then
             ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
             ! l'operateur d'ordre pour les derivees d'ordre 2 

              
             do r=fft_start(3),fft_end(3)
               do s = fft_start(2),fft_end(2)
                  do q = fft_start(1),fft_end(1)

                      Rot_fieldF(q,s,r,indice_lin_t2(i,j)) = Rot_fieldF(q,s,r,indice_lin_t2(i,j)) &
                                     + imP*perm2*FREQ(q,s,r,k)*fieldF(q,s,r,indice_lin_st2(i,m))
                  end do
               end do
             end do

           else 
             ! Calcul des derivees a partir des operateurs differences finies sur cube de 27 voxels
             ! l'operateur field_first_order_partial_centeredF somme les contributions successives
             ! des derivees partielles (,im) dans Laplace_fieldF                                                   
                 call field_first_order_partial_centeredF(Rot_fieldF(:,:,:,indice_lin_t2(i,j)),&
                                   perm2,fieldF(:,:,:,indice_lin_st2(i,m)),k,1,1,N(k),d(k))
             
           end if
          end if

              ! fin calcul de la derivee seconde permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * dfieldF(l)/dXmdXj
        end do
     end do
  ! Fin calcul de la composante (i,j)
  end do
  end do

end subroutine rot_sym_tensor2F

!==============================================================================
!       SUBROUTINE ROT_TENSOR2F
!------------------------------------------------------------------------------
!> Evalue le rot(rot()) du champ de tenseur d'ordre 2 en entree, dans l'espace de Fourier
!!
!!   \param[in] fieldF         champ de tenseur ordre2 entree 
!!                             tableau complexe (nxP/2,nyP,nzP,9) (profil indifferent)
!!
!!   \param[out]rot_fieldF     rotationnel du champ entree : derivee partielle du champ entree 
!!                             tableau complexe (nxP/2,nyP,nzP,9) (profil indifferent)
!!
!!   \param[in] FREQ           tableau de fréquence, espace Fourier 
!!                             pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                             optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et Rot_fieldF)
!!
!
!  --------
!
!  rot(T)(i,j) = perm(j,m,k) * T(i,m),k
!!
!! Peut etre calcule de deux façons : 
!!      A) en utilisant les operateurs de derivations definis dans le field_mod
!          field_(first/second)_order_partial_centeredF
!          equivalent a une approximation differences finies sur un cube de 3 voxels
!          de cote
!
!       B) en utilisant les frequences initialisees dans le green_mod ce qui equivaut
!          a appliquer deux fois l'operateur de derivation du premier ordre 
!          correspondant au filtre choisi dans algo_param. 
!          Dans le cas difference finies cela correspond a un calcul sur une grille
!          de 5 voxels de cote
!          APPLIQUE SI LE TABLEAU FREQ EST PASSE EN ARGUMENT
!
!==============================================================================
subroutine rot_tensor2F(fieldF,Rot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),&
                            intent(out)    :: Rot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: k,i,j,m      ! variables de boucles formule rotationnel
  integer                                  :: r,s,q        ! variables de boucle calcul derivee
  real(mytype)                             :: perm2
  complex(mytype)                          :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i

 Rot_fieldF(:,:,:,:) = 0._mytype

  do i=1,3
  do j=1,3
  ! Calcul composante (i,j) du rotationnel
     do k=1,3
       do m=1,3
             ! calcul de perm(j,m,k) * T(i,m),k
             ! on fait le lien entre la composante (i,j) / (i,m) des champs en ecriture conventionnelle
             ! et leur indice dans les RoRot_field et fieldF a l'aide de la fonction indice_lin_t2

         perm2 = real(permutation_tensor(j,m,k),mytype)
            if (perm2 .ne. 0) then
              
               if (present(FREQ)) then
                 ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
                 ! l'operateur d'ordre pour les derivees d'ordre 2 

                  
                 do r=fft_start(3),fft_end(3)
                   do s = fft_start(2),fft_end(2)
                      do q = fft_start(1),fft_end(1)
                        Rot_fieldF(q,s,r,indice_lin_t2(i,j)) = Rot_fieldF(q,s,r,indice_lin_t2(i,j)) &
                                     + imP*perm2*FREQ(q,s,r,k)*fieldF(q,s,r,indice_lin_t2(i,m))
                      end do
                   end do
                 end do

               else 
                 ! Calcul des derivees a partir des operateurs differences finies sur cube de 9 voxels
                 ! l'operateur field_second_order_partial_centeredF somme les contributions successives
                 ! des derivees partielles (,im) dans Laplace_fieldF                                         
                 call field_first_order_partial_centeredF(Rot_fieldF(:,:,:,indice_lin_t2(i,j)),&
                                   perm2,fieldF(:,:,:,indice_lin_t2(i,m)),k,1,1,N(k),d(k))
                
              end if
             end if

              ! fin calcul de la derivee seconde permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * dfieldF(l)/dXmdXj
        end do
     end do
  ! Fin calcul de la composante (i,j)
  end do
  end do
end subroutine rot_tensor2F


!==============================================================================
!       SUBROUTINE ROTROTF
!------------------------------------------------------------------------------
!> 
!> Evalue ROT(ROT(xx)) pour differents type de xx : 
!>                 vecteur, tenseur d'ordre 2 sym. et non sym. 
!>
!>  Verifie la taille du champ en entree et le renvoie vers la fonction rotationnel
!>  correspondante
!!
!!   \param[in] fieldF         champ de tenseur entree 
!!                             tableau complexe (nxP/2,nyP,nzP,3 6 ou 9) (profil impose)
!!
!!   \param[out] RotRot_fieldFD rotationnel du rotationnel du champ entree 
!!                              tableau complexe (nxP/2,nyP,nzP,3 6 ou 9) (profil impose)
!!
!!   \param[in] FREQ           tableau de fréquence, espace Fourier 
!!                             pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                             optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!  
!==============================================================================
subroutine rotrotF(fieldF,RotRot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(:,:,:,:),intent(in)     :: fieldF 
  complex(mytype),dimension(:,:,:,:),intent(out)    :: RotRot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Variables locales
  integer                    :: size_tens, size_grad

  size_tens = size(fieldF,4)
  size_grad = size(RotRot_fieldF,4)

  select case(size_tens)
      case(3) ! rot(rot(vecteur))
         if (size_grad == 3) then
            call rotrot_vectorF(fieldF,RotRot_fieldF,N,d,FREQ)
         else
            call amitex_abort("Taille du champ de sortie incorrecte, attendu=(:,:,:,3) (rotrot)" ,2,0)  
         end if
      case(6) ! rot(rot(tenseur_sym))
         if (size_grad == 9) then
            call rotrot_sym_tensor2F(fieldF,RotRot_fieldF,N,d,FREQ)
         else
            call amitex_abort("Taille du champ de sortie incorrecte, attendu=(:,:,:,9) (rotrot)" ,2,0)  
         end if
      case(9) ! rot(rot(tenseur))
         if (size_grad == 9) then
            call rotrot_tensor2F(fieldF,RotRot_fieldF,N,d,FREQ)
         else
            call amitex_abort("Taille du champ de sortie incorrecte, attendu=(:,:,:,9) (rotrot)" ,2,0)  
         end if   
      case default   
            call amitex_abort("Taille du champ en entree incorrecte, rotrot ne peut s'appliquer&
                 & que a des champs de taille 3 (vecteurs) ,6 (tenseurs sym)et 9 (tenseurs)",2,0)     
   end select

end subroutine rotrotF

!==============================================================================
!       SUBROUTINE ROTROT_VECTORF
!------------------------------------------------------------------------------
!> Evalue le rot(rot()) du champ de vecteurs en entree dans l'espace de Fourier
!!
!!   \param[in] fieldF         champ de vecteur entree 
!!                             tableau complexe (nxP/2,nyP,nzP,3) (profil indifferent)
!!
!!   \param[out] RorRot_fieldF rotationnel du champ entree : derivee partielle du champ entree 
!!                             tableau complexe (nxP/2,nyP,nzP,3) (profil indifferent)
!!
!!   \param[in] FREQ           tableau de fréquence, espace Fourier 
!!                             pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                             optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et RotRot_fieldF)
!!
!
!  formule :
!  --------
!
!  rot(rot(U))(k) = perm(i,j,k) * perm(i,l,m) * U(l),mj
!!
!! Peut etre calcule de deux façons : 
!!      A) en utilisant les operateurs de derivation centree definis dans le field_mod
!          field_(first/second)_order_partial_centeredF
!          equivalent a une approximation differences finies sur un cube de 3 voxels
!          de cote
!
!       B) en utilisant les frequences initialisees dans le green_mod ce qui equivaut
!          a appliquer deux fois l'operateur de derivation du premier ordre 
!          correspondant au filtre choisi dans algo_param. 
!          APPLIQUE SI LE TABLEAU FREQ EST PASSE EN ARGUMENT
!
!       Test sur un champ de vecteur compose de sinus, sur une cellule de taille 1
!       pulsations : W2 = 4pi W3 = 6pi W4 = 8pi
!
!       u(x) = sin(w3y + w4z) u(y) = sin(w2x + w4z) u(z) = sin(w4x + w3y)
!
!       discretisation 100x100x100, erreur sur la derivee 
!       par rapport a la solution analytique :
!       methode A : 0.4% 
!       methode B : 1.5%
!       La methode B est 2x plus rapide que la A sur ce cas pour le calcul du rotrot
!
!==============================================================================
subroutine rotrot_vectorF(fieldF,RotRot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),3),&
                            intent(out)    :: RotRot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: k,i,j,l,m    ! variables de boucles formule rotationnel
  integer                                  :: r,s,q        ! variables de boucle calcul derivee
  real(mytype)                             :: perm2

  RotRot_fieldF(:,:,:,:) = 0._mytype

  do k=1,3
  ! Calcul composante k  du rotationnel
     do j=1,3
        do l=1,3
           do m=1,3
              do i=1,3
                 ! calcul de permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * fieldF(l),mj


                 perm2 = real(permutation_tensor(i,j,k)*permutation_tensor(i,l,m),mytype)
                if (perm2 .ne. 0) then
                  
                   if (present(FREQ)) then
                     ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
                     ! l'operateur d'ordre pour les derivees d'ordre 2

                      
                     do r=fft_start(3),fft_end(3)
                       do s = fft_start(2),fft_end(2)
                          do q = fft_start(1),fft_end(1)
                              RotRot_fieldF(q,s,r,k) = RotRot_fieldF(q,s,r,k) - perm2*FREQ(q,s,r,m)*FREQ(q,s,r,j)*fieldF(q,s,r,l)
                          end do
                       end do
                     end do

                   else 
                     ! Calcul des derivees a partir des operateurs differences finies sur cube de 9 voxels
                     ! l'operateur fieldF_second_order_partial_centered somme les contributions successives
                     ! des derivees partielles (,mj) dans Laplace_fieldF                   
                     call field_second_order_partial_centeredF(RotRot_fieldF(:,:,:,k),&
                                       perm2,fieldF(:,:,:,l),(/ m,j /),1,1,N(m),N(j),d(m),d(j))

                  end if
                 end if

                  ! fin calcul de la derivee seconde permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * dfieldF(l)/dXmdXj
              end do 
           end do
        end do
     end do
  ! Fin calcul composante k
  end do  

end subroutine rotrot_vectorF

!==============================================================================
!       SUBROUTINE ROTROT_SYM_TENSOR2F
!------------------------------------------------------------------------------
!> Evalue le rot(rot()) du champ de tenseur symmétrique d'ordre 2 en entree, dans l'espace de Fourier
!!
!!   \param[in] fieldF         champ de tenseur ordre 2 sym entree 
!!                             tableau complexe (nxP/2,nyP,nzP,6) (profil indifferent)
!!
!!   \param[out] RorRot_fieldF rotationnel du champ entree : derivee partielle du champ entree 
!!                             tableau complexe (nxP/2,nyP,nzP,9) (profil indifferent)
!!
!!   \param[in] FREQ           tableau de fréquence, espace Fourier 
!!                             pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                             optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et RotRot_fieldF)
!!
!
!  formule :
!  --------
!
!  rot(rot(T))(i,j) = perm(t,j,k) * perm(t,l,m) * T(i,m),lk
!!
!! Peut etre calcule de deux façons : 
!!      A) en utilisant les operateurs de derivation centree definis dans le field_mod
!          field_(first/second)_order_partial_centeredF
!          equivalent a une approximation differences finies sur un cube de 3 voxels
!          de cote
!
!       B) en utilisant les frequences initialisees dans le green_mod ce qui equivaut
!          a appliquer deux fois l'operateur de derivation du premier ordre 
!          correspondant au filtre choisi dans algo_param. 
!          APPLIQUE SI LE TABLEAU FREQ EST PASSE EN ARGUMENT
!
!==============================================================================
subroutine rotrot_sym_tensor2F(fieldF,RotRot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),6),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),&
                            intent(out)    :: RotRot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: k,i,j,l,m,t  ! variables de boucles formule rotationnel
  integer                                  :: r,s,q        ! variables de boucle calcul derivee
  real(mytype)                             :: perm2

  RotRot_fieldF(:,:,:,:) = 0._mytype

  do i=1,3
  do j=1,3
  ! Calcul composante (i,j) du rotationnel
     do k=1,3
        do l=1,3
           do m=1,3
              do t=1,3
                 ! calcul de permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * fieldF(i,m),lk
                 ! on fait le lien entre la composante (i,j) / (i,m) des champs en ecriture conventionnelle
                 ! et leur indice dans les RoRot_field et fieldF a l'aide des fonctions indice_lin_st2 et indice_lin_t2

                 perm2 = real(permutation_tensor(t,j,k)*permutation_tensor(t,l,m),mytype)
                if (perm2 .ne. 0) then
                  
                   if (present(FREQ)) then
                     ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
                     ! l'operateur d'ordre pour les derivees d'ordre 2 

                      
                     do r=fft_start(3),fft_end(3)
                       do s = fft_start(2),fft_end(2)
                          do q = fft_start(1),fft_end(1)

                              RotRot_fieldF(q,s,r,indice_lin_t2(i,j)) = RotRot_fieldF(q,s,r,indice_lin_t2(i,j)) &
                                             - perm2*FREQ(q,s,r,l)*FREQ(q,s,r,k)*fieldF(q,s,r,indice_lin_st2(i,m))
                          end do
                       end do
                     end do

                   else 
                     ! Calcul des derivees a partir des operateurs differences finies sur cube de 27 voxels
                     ! l'operateur field_second_order_partial_centeredF somme les contributions successives
                     ! des derivees partielles (,im) dans Laplace_fieldF                                        
                     call field_second_order_partial_centeredF(RotRot_fieldF(:,:,:,indice_lin_t2(i,j)),&
                                       perm2,fieldF(:,:,:,indice_lin_st2(i,m)),(/ l,k /),1,1,N(l),N(k),d(l),d(k))
                     
                  end if
                 end if

                  ! fin calcul de la derivee seconde permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * dfieldF(l)/dXmdXj
              end do 
           end do
        end do
     end do
  ! Fin calcul de la composante (i,j)
  end do
  end do

end subroutine rotrot_sym_tensor2F


!==============================================================================
!       SUBROUTINE ROTROT_TENSOR2F
!------------------------------------------------------------------------------
!> Evalue le rot(rot()) du champ de tenseur d'ordre 2 en entree, dans l'espace de Fourier
!!
!!   \param[in] fieldF         champ de tenseur ordre2 entree 
!!                             tableau complexe (nxP/2,nyP,nzP,9) (profil indifferent)
!!
!!   \param[in] D_fieldF       rotationnel du champ entree : derivee partielle du champ entree 
!!                             tableau complexe (nxP/2,nyP,nzP,9) (profil indifferent)
!!
!!   \param[in] FREQ           tableau de fréquence, espace Fourier 
!!                             pinceaux-Z (nxP/2,nyP,nzP,3) 
!!                             optionnel
!!
!!   \param[in]  N             dimensions de la cellule [Nx Ny Nz] 
!!   \param[in]  d             dimension des voxels [dx,dy,dz]
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et RotRot_fieldF)
!!
!
!  --------
!
!  rot(rot(T))(i,j) = perm(t,j,k) * perm(t,l,m) * T(i,m),lk
!!
!! Peut etre calcule de deux façons : 
!!      A) en utilisant les operateurs de derivations definis dans le field_mod
!          field_(first/second)_order_partial_centeredF
!          equivalent a une approximation differences finies sur un cube de 3 voxels
!          de cote
!
!       B) en utilisant les frequences initialisees dans le green_mod ce qui equivaut
!          a appliquer deux fois l'operateur de derivation du premier ordre 
!          correspondant au filtre choisi dans algo_param. 
!          Dans le cas difference finies cela correspond a un calcul sur une grille
!          de 5 voxels de cote
!          APPLIQUE SI LE TABLEAU FREQ EST PASSE EN ARGUMENT
!
!==============================================================================
subroutine rotrot_tensor2F(fieldF,RotRot_fieldF,N,d,FREQ)

  implicit none
  ! Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),&
                             intent(in)    :: fieldF 
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),9),&
                            intent(out)    :: RotRot_fieldF 
  integer,dimension(3),intent(in)          :: N
  real(mytype),dimension(3),intent(in)     :: d

  ! Entrees optionnelles
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
                      optional,intent(in)  :: FREQ

  ! Autre variables
  integer                                  :: k,i,j,l,m,t  ! variables de boucles formule rotationnel
  integer                                  :: r,s,q        ! variables de boucle calcul derivee
  real(mytype)                             :: perm2

  RotRot_fieldF(:,:,:,:) = 0._mytype

  do i=1,3
  do j=1,3
  ! Calcul composante (i,j) du rotationnel
     do k=1,3
        do l=1,3
           do m=1,3
              do t=1,3
                 ! calcul de permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * fieldF(i,m),lk
                 ! on fait le lien entre la composante (i,j) / (i,m) des champs en ecriture conventionnelle
                 ! et leur indice dans les RoRot_field et fieldF a l'aide de la fonction indice_lin_t2

                 perm2 = real(permutation_tensor(t,j,k)*permutation_tensor(t,l,m),mytype)
                if (perm2 .ne. 0) then
                  
                   if (present(FREQ)) then
                     ! Calcul des derivees a partir des frequences definies dans Green_mod = application successive de 
                     ! l'operateur d'ordre pour les derivees d'ordre 2 

                      
                     do r=fft_start(3),fft_end(3)
                       do s = fft_start(2),fft_end(2)
                          do q = fft_start(1),fft_end(1)

                              RotRot_fieldF(q,s,r,indice_lin_t2(i,j)) = RotRot_fieldF(q,s,r,indice_lin_t2(i,j)) &
                                            - perm2*FREQ(q,s,r,l)*FREQ(q,s,r,k)*fieldF(q,s,r,indice_lin_t2(i,m))
                          end do
                       end do
                     end do

                   else 
                     ! Calcul des derivees a partir des operateurs differences finies sur cube de 9 voxels
                     ! l'operateur field_second_order_partial_centeredF somme les contributions successives
                     ! des derivees partielles (,im) dans Laplace_fieldF                   
                     call field_second_order_partial_centeredF(RotRot_fieldF(:,:,:,indice_lin_t2(i,j)),&
                                      perm2,fieldF(:,:,:,indice_lin_t2(i,m)),(/ l,k /),1,1,N(l),N(k),d(l),d(k))
                    
                  end if
                 end if

                  ! fin calcul de la derivee seconde permutation_tensor(i,j,k)*permutation_tensor(i,l,m) * dfieldF(l)/dXmdXj
              end do 
           end do
        end do
     end do
  ! Fin calcul de la composante (i,j)
  end do
  end do
end subroutine rotrot_tensor2F

!==============================================================================
!       FUNCTION PERMUTATION_TENSOR
!------------------------------------------------------------------------------
!> Evalue la valeur (entier) du tenseur de permutation de Levi_Civita pour un triplet
!! d'entier (sa composante indicee par les trois entiers en entree)
!==============================================================================
function permutation_tensor(i,j,k)

  implicit none
  integer,intent(in)     :: i,j,k
  integer                :: permutation_tensor

  permutation_tensor = (i-j)*(j-k)*(k-i)/2

end function permutation_tensor

!==============================================================================
!       FUNCTION INDICE_LIN_T2
!------------------------------------------------------------------------------
!> Renvoie l'indice lineaire d'une composante d'un tenseur d'ordre 2
!  sous sa forme vectorielle a partir de ses deux indices de composantes
!  sous sa forme matricielle
!
!  forme matricielle T = T11 T12 T13  =  T1 T4 T5
!                        T21 T22 T23     T7 T2 T6
!                        T31 T32 T33     T8 T9 T3
!
!  forme vectorielle T = T1 T2 T3 T4 T5 T6 T7 T8 T9
!
!  Correspondance : 1=11 2=22 3=33 4=12 5=13 6=23 7=21 8=31 9=32
!
!  \param[in]  :  i,j  indices de la composante sous forme matricielle
!
!  \param[out]      indice de la composante sous forme vectorielle
!==============================================================================
function indice_lin_t2(i,j)

  implicit none
  integer,intent(in)     :: i,j
  integer                :: indice_lin_t2

  select case(i)
    case(1)
      select case(j)
         case(1)
            indice_lin_t2 = 1
         case(2)
            indice_lin_t2 = 4
         case(3)
            indice_lin_t2 = 5
      end select
    case(2)
      select case(j)
         case(1)
            indice_lin_t2 = 7
         case(2)
            indice_lin_t2 = 2
         case(3)
            indice_lin_t2 = 6
      end select
    case(3)
      select case(j)
         case(1)
            indice_lin_t2 = 8
         case(2)
            indice_lin_t2 = 9
         case(3)
            indice_lin_t2 = 3
      end select
  end select

end function indice_lin_t2

!==============================================================================
!       FUNCTION INDICE_LIN_ST2
!------------------------------------------------------------------------------
!> Renvoie l'indice lineaire d'une composante d'un tenseur d'ordre 2 symmetrique
!  sous sa forme vectorielle a partir de ses deux indices de composantes
!  sous sa forme matricielle
!
!  forme matricielle T = T11 T12 T13  =  T1 T4 T5
!                        T12 T22 T23     T4 T2 T6
!                        T13 T23 T33     T5 T6 T3
!
!  forme vectorielle T = T1 T2 T3 T4 T5 T6
!
!  Correspondance : 1=11 2=22 3=33 4=12 5=13 6=23
!
!  \param[in]  :  i,j  indices de la composante sous forme matricielle
!
!  \param[out] : indice_lin_st2  indice de la composante sous forme vectorielle
!==============================================================================
function indice_lin_st2(i,j)

  implicit none
  integer,intent(in)     :: i,j
  integer                :: indice_lin_st2

  select case(i)
    case(1)
      select case(j)
         case(1)
            indice_lin_st2 = 1
         case(2)
            indice_lin_st2 = 4
         case(3)
            indice_lin_st2 = 5
      end select
    case(2)
      select case(j)
         case(1)
            indice_lin_st2 = 4
         case(2)
            indice_lin_st2 = 2
         case(3)
            indice_lin_st2 = 6
      end select
    case(3)
      select case(j)
         case(1)
            indice_lin_st2 = 5
         case(2)
            indice_lin_st2 = 6
         case(3)
            indice_lin_st2 = 3
      end select
  end select


end function indice_lin_st2



!===========================================================================================================
!> gene_cartesian_coord(X,Y,Z,dx,dy,dz) generate X,Y,Z fields of coordinates of voxel's centers
!>                                      the lower coordinates corner being (0.,0.,0.)
!!
!!   \param[in]  dx,dy,dz       voxel sizes
!!   \param[out] X,Y,Z          coordinates of  grid centers or nodes, double prec array(nx,ny,nz)
!!
!===========================================================================================================
subroutine gene_cartesian_coord(X,Y,Z,dx,dy,dz)
  implicit none

!!------------------------------------------------------------------------------
!>                                                                        IN/OUT

  real(mytype),intent(out),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: X,Y,Z  
  real(mytype),intent(in)         :: dx,dy,dz

!!------------------------------------------------------------------------------
!>                                                                        OTHERS
  real(mytype)                    :: dx_2,dy_2,dz_2 
  integer                         :: i,j,k

!!------------------------------------------------------------------------------
!>                                                                        COORDS 

  dx_2 = dx / 2._mytype
  dy_2 = dy / 2._mytype
  dz_2 = dz / 2._mytype


  do k = xstart(3),xend(3)
  do j = xstart(2),xend(2)
  do i = xstart(1),xend(1)
     X(i,j,k) = i*dx 
     Y(i,j,k) = j*dy 
     Z(i,j,k) = k*dz 
  end do
  end do
  end do
  X = X - dx_2
  Y = Y - dy_2
  Z = Z - dz_2

end subroutine gene_cartesian_coord

!===========================================================================================================
!> add_gradient_scalar_field(field,Gx,Gy,Gz,dx,dy,dz) 
!>                      field = field + fieldG with fieldG = G.X
!>                                              and fieldG(1,1,1) = 0
!!
!!   \param[in]  Gx,Gy,Gz       gradient components
!!   \param[in]  dx,dy,dz       voxel sizes
!!   \param[inout] field          
!!
!===========================================================================================================
subroutine add_gradient_scalar_field(field, Gx,Gy,Gz,dx,dy,dz)
  implicit none

!!------------------------------------------------------------------------------
!>                                                                        IN/OUT
 
  real(mytype), intent(in) :: Gx, Gy, Gz, dx, dy, dz 
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: field  

!!------------------------------------------------------------------------------
!>                                                                        OTHERS

  integer                         :: i,j,k

!!------------------------------------------------------------------------------
!>                                                            GEN GRADIENT FIELD

  do k = xstart(3),xend(3)
  do j = xstart(2),xend(2)
  do i = xstart(1),xend(1)
     field(i,j,k) = field(i,j,k) + Gx*(i-1)*dx + Gy*(j-1)*dy + Gz*(k-1)*dz 
  end do
  end do
  end do

end subroutine add_gradient_scalar_field


!===========================================================================================================
!> add_gradgradU(Def,gradgradU,gradgradU_evol) : AJOUTE A DEF UN CHAMP DE GRADIENT DE U LINEAIRE
!!
!!   \param[inout]  Def         champ entree, tableau (nxP,nyP,nzP,6 ou 9) (profil indifferent)                           
!!   \param[in]  Ntens          nombre de variable pour la deformation (6 ou 9)
!!   \param[in]  gradgradU      composantes de gradgradU dimension (27 ou (3,3,3) dans la fonction)                             
!!   \param[in]  gradgradU_evol masque des composantes de gradgradU : valeur -1 = composante non imposee     
!!   \param[in]  dx,dy,dz       dimensions des voxels (reels)    
!!   \param[in]  nx,ny,nz       dimensions de la grille  (entiers)               
!!
!!  Dans le cas 6 composantes : on utlise la notation de Voigt pour les cisailments
!!                              et l'ordre 11,22,33,12,13,23
!!
!===========================================================================================================
subroutine add_gradgradU(Def,Ntens,gradgradU,gradgradU_evol,dx,dy,dz,nx,ny,nz)
  
  implicit none

  integer, intent(in)                       :: Ntens
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),Ntens), intent(inout) :: Def
  real(mytype), dimension(3,3,3),intent(in) :: gradgradU
  integer, dimension(:,:,:),intent(in)      :: gradgradU_evol
  real(mytype),intent(in)                   :: dx,dy,dz
  integer,intent(in)                        :: nx,ny,nz

  real(mytype),dimension(3)                 :: X
  integer                                   :: i,j,k     !< indices sur la grille
  integer                                   :: l,m,n,N0  !< indices sur les composantes de gradgradU
  integer, dimension(3,3)                   :: notation6, notation9

  
  do l = 1,3
  do m = 1,3
     if (l==m) then
       notation6(l,m) = l
       notation9(l,m) = l
     else
       notation6(l,m) = l + m + 1
       notation9(l,m) = notation6(l,m)
       if (l>m) then
          notation9(l,m) = notation9(l,m) + 3
       end if
     end if
  end do
  end do

  if (Ntens==6) then
  do k = xstart(3),xend(3)
  do j = xstart(2),xend(2)
  do i = xstart(1),xend(1)
     X(1) = i*dx - dx / 2._mytype - (nx*dx / 2._mytype)
     X(2) = j*dy - dy / 2._mytype - (ny*dy / 2._mytype)
     X(3) = k*dz - dz / 2._mytype - (nz*dz / 2._mytype)
     do n = 1,3
     do m = 1,3
     do l = 1,3
       if (gradgradU_evol(l,m,n) .ne. -1) then
          N0 = notation6(l,m)
          Def(i,j,k,N0) = Def(i,j,k,N0) + gradgradU(l,m,n)*X(n) ! valable pour toutes les composantes
       end if
     end do
     end do
     end do
  end do
  end do
  end do
  end if

  if (Ntens==9) then
  do k = xstart(3),xend(3)
  do j = xstart(2),xend(2)
  do i = xstart(1),xend(1)
     X(1) = i*dx - dx / 2._mytype - (nx*dx / 2._mytype)
     X(2) = j*dy - dy / 2._mytype - (ny*dy / 2._mytype)
     X(3) = k*dz - dz / 2._mytype - (nz*dz / 2._mytype)
     do n = 1,3
     do m = 1,3
     do l = 1,3
       if (gradgradU_evol(l,m,n) .ne. -1) then
          N0 = notation9(l,m)
          Def(i,j,k,N0) = Def(i,j,k,N0) + gradgradU(l,m,n)*X(n) ! valable pour toutes les composantes
       end if
     end do
     end do
     end do
  end do
  end do
  end do
  end if

end subroutine add_gradgradU
!===========================================================================================================
!> field_mask_outbox(F,ncomp,Imin,Imax,maskvalue) 
!>   apply a mask with maskvalue exclusively OUT of a box defined by its indices Imin,Imax
!>   
!!
!!   \param[inout]  F         field
!!   \param[in]     ncomp     number of components of F
!!   \param[in]     Imin,Imax integer(3), indices of the box
!!   \param[in]     maskvalue value assigned to the region of the mask 
!!
!===========================================================================================================
subroutine field_mask_outbox(F,ncomp,Imin,Imax,maskvalue) 
  implicit none

!!------------------------------------------------------------------------------
!>                                                                        IN/OUT
  integer,intent(in)              :: ncomp
  real(mytype),intent(inout),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncomp) :: F
  integer,dimension(3),intent(in) :: Imin,Imax
  real(mytype),intent(in)         :: maskvalue
  
!!------------------------------------------------------------------------------
!>                                                                        OTHERS
  integer  :: i

!!------------------------------------------------------------------------------
!>                                                                      MASK OUT

  
  do i = 1,3
    if (Imin(i) > xend(i))   F = maskvalue
    if (Imax(i) < xstart(i)) F = maskvalue
  end do   
  
  if (Imin(1) > xstart(1)  .AND. Imin(1) <= xend(1) ) F(xstart(1):Imin(1)-1,:,:,:) = maskvalue
  if (Imin(2) > xstart(2)  .AND. Imin(2) <= xend(2) ) F(:,xstart(2):Imin(2)-1,:,:) = maskvalue
  if (Imin(3) > xstart(3)  .AND. Imin(3) <= xend(3) ) F(:,:,xstart(3):Imin(3)-1,:) = maskvalue
    
  if (Imax(1) >= xstart(1) .AND. Imax(1) < xend(1) ) F(Imax(1)+1:xend(1),:,:,:)    = maskvalue
  if (Imax(2) >= xstart(2) .AND. Imax(2) < xend(2) ) F(:,Imax(2)+1:xend(2),:,:)    = maskvalue
  if (Imax(3) >= xstart(3) .AND. Imax(3) < xend(3) ) F(:,:,Imax(3)+1:xend(3),:)    = maskvalue

end subroutine field_mask_outbox

!==============================================================================
!       SUBROUTINE FIELD_FILTERF
!------------------------------------------------------------------------------
!> Applique un Filtre a un champ dans l'espace de Fourier 
!!
!!              field_filterF(fieldF,FILTER,ncomp,nvar,N,d)
!!
!!              FieldF = Filtre * fieldF
!!
!! FILTER : 1   : Filtre "Moyenne" 8 pts (avec changement de support noeuds->centres)
!!          -1  : Filtre "Moyenne" 8 pts (avec changement de support centres->noeuds)
!!          2   : Filtre "Moyenne" 8 pts applique 2 fois (retour support initial)
!!
!!          100 : Filtre de Gauss isotrope, param(1) = sigma 
!!                   g=Aexp(-X^2 / (2.sigma^2) (A = facteur de normalisation (integrale = 1)
!!
!!   \param[inout] fieldF      champ entree/sortie, (nxP/2,nyP,nzP,ncomp,nvar)    (profil indifferent)
!!   \param[in]  ncomp         nombre de composantes du champ 
!!   \param[in]  nvar          nombre de champs (1 en mecanique)
!!
!!   \param[in]  N(3)          dimension de la cellule dans la direction dir
!!   \param[in]  d(3)          taille des voxels dans la direction dir
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!! REMARQUE : les filtres 1, -1, 2 ne dependant pas de la dimension d des voxels  
!!
!==============================================================================
subroutine field_filterF(fieldF,FILTER,ncomp,nvar,N,d,param)

  implicit none
  
  !Entrees/Sorties :
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncomp,nvar),&
                intent(inout)          :: fieldF 
  integer,intent(in)                   :: ncomp,nvar,FILTER
  integer,dimension(3),intent(in)      :: N
  real(mytype),dimension(3),intent(in) :: d
  real(mytype),dimension(:),intent(in),optional :: param
  
  ! Autre variables
  integer                     :: i,j,k,l,m ! variables de boucles
#ifdef DOUBLE_PREC
  real(mytype),parameter      :: PI = 4._mytype*DATAN(1._mytype)
#else
  real(mytype),parameter      :: PI = 4._mytype*ATAN(1._mytype)
#endif
  real(mytype),dimension(3):: F          ! Frequence
  complex(mytype)             :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i
  
  integer,dimension(1:N(1)/2 + 1) :: idx    !> Utiles pour definition des frequences (no_filter)
  integer,dimension(1:N(2))       :: idy
  integer,dimension(1:N(3))       :: idz
  integer                         :: nx_2, ny_2, nz_2
  real(mytype)                    :: DF1,DF2,DF3
  real(mytype)                    :: fact, FiltreF


  !--- Definition de variables pour l'application de Filtre Gaussien et autres
  fact = 0.  
  DF1 = 2._mytype * PI / (float(N(1))*d(1))
  DF2 = 2._mytype * PI / (float(N(2))*d(2))
  DF3 = 2._mytype * PI / (float(N(3))*d(3))
  
  ! indices pour ramener les frequences nulles en (1,1,1)
  ! attention aux cas pairs et impairs
  nx_2 = N(1) / 2
  idx = (/(nx_2 + i,i=1,nx_2+1) /)

  ny_2 = N(2) / 2
  if (mod(N(2),2) .eq. 0) then
     idy = (/(ny_2 + i,i=1,ny_2),(i,i=1,ny_2)/)
  else
     idy = (/(ny_2 + i,i=1,ny_2+1),(i,i=1,ny_2)/)
  endif

  nz_2 = N(3) / 2
  if (mod(N(3),2) .eq. 0) then
     idz = (/(nz_2 + i,i=1,nz_2),(i,i=1,nz_2)/)
  else
     idz = (/(nz_2 + i,i=1,nz_2+1),(i,i=1,nz_2)/)
  endif

  !--- Filtre 8 voisins applique UNE fois (Centres => Noeuds)
  if (FILTER == -1) then 
  do m=1,nvar
  do l=1,ncomp
  do k=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do i=fft_start(1),fft_end(1)
    F(1) = PI*(i-1)/N(1)   ! attention ici F= F/2.
    F(2) = PI*(j-1)/N(2)
    F(3) = PI*(k-1)/N(3)
    fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * &
                        (cos(F(1)) * cos(F(2)) * cos(F(3))) * &! attention ici F= F/2.
                        exp(-imP*(F(1)+F(2)+F(3))) 
  end do
  end do
  end do
  end do
  end do
  
  !--- Filtre 8 voisins applique UNE fois (Noeuds => Centres)
  elseif (FILTER == 1) then 
  do m=1,nvar
  do l=1,ncomp
  do k=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do i=fft_start(1),fft_end(1)
    F(1) = PI*(i-1)/N(1)   ! attention ici F= F/2.
    F(2) = PI*(j-1)/N(2)
    F(3) = PI*(k-1)/N(3)
    fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * &
                        (cos(F(1)) * cos(F(2)) * cos(F(3))) * & ! attention ici F= F/2.
                        exp(imP*(F(1)+F(2)+F(3))) 
  end do
  end do
  end do
  end do
  end do
  
  !--- Filtre 8 voisins applique DEUX fois (support inchange : C => N => C ou N => C => N)
  elseif (FILTER == 2) then
  do m=1,nvar
  do l=1,ncomp
  do k=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do i=fft_start(1),fft_end(1)
    F(1) = PI*(i-1)/N(1)   ! attention ici F= F/2.
    F(2) = PI*(j-1)/N(2)
    F(3) = PI*(k-1)/N(3)
    fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * &
                         (cos(F(1)) * cos(F(2)) * cos(F(3))) ** 2 ! attention ici F= F/2.
  end do
  end do
  end do
  end do
  end do
  
  !--- Filtre de GAUSS
  !    
  !    g = Aexp(-alpha.X^2)   avec A = (alpha/pi)^3/2 (normalisation)
  !
  !      -> gF = exp(-freq^2/(4.alpha))
  ! 
  !      alpha = 1/(2.sigma^2)  -> 1/(4.alpha) = sigma^2 / 2
  !
  elseif (FILTER == 100) then
  
    fact = param(1)**2 / 2._mytype
    do m=1,nvar
    do l=1,ncomp
    do k=fft_start(3),fft_end(3)
    do j=fft_start(2),fft_end(2)
    do i=fft_start(1),fft_end(1)
      F(1) = (-nx_2 - 1 + idx(i)) * DF1
      F(2) = (-ny_2 - 1 + idy(j)) * DF2
      F(3) = (-nz_2 - 1 + idz(k)) * DF3
      FiltreF = exp(-sum(F*F) * fact)
      fieldF(i,j,k,l,m) = fieldF(i,j,k,l,m) * FiltreF
    end do
    end do
    end do
    end do
    end do
  
  else
     call amitex_abort("field_filterF : FILTER different from -1,1,2,100 - not yet implemented " ,2,0)  
  end if

end subroutine field_filterF

!==============================================================================
!       SUBROUTINE FIELD_FILTER
!------------------------------------------------------------------------------
!> Applique un Filtre "moyenn" a un champ dans l'espace reel
!!
!!              field_filter(field,FILTER,ncomp,nvar,N,d)
!!
!!              Field = Filtre * field
!!
!! FILTER : 1  : Filter "Mean" 8 pts (avec changement de support noeuds->centres)
!!          -1 : Filter "Mean" 8 pts (avec changement de support centres->noeuds)
!!          1000 : Filter "Sum" 6 pts - sum on 6 closest neighbors
!!                            WARNING : do not include the central voxel!!
!!
!! Pour FILTER = 2 analogue a field_filterF : 
!!               Appliquer le filtre successivement avec FILTER=1 puis FILTER=-1
!!
!!   \param[inout] field      champ entree/sortie, (nxP,nyP,nzP,ncomp,nvar)    (profil indifferent)
!!   \param[in]  ncomp         nombre de composantes du champ 
!!   \param[in]  nvar          nombre de champs (1 en mecanique)
!!
!!   optionel 
!!   \param[in]  N(3)          dimension de la cellule dans la direction dir
!!   \param[in]  d(3)          taille des voxels dans la direction dir
!!
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!! REMARQUE : les filtres 1, -1, ne dependant pas des entrees N et d
!!
!==============================================================================

subroutine field_filter(field,FILTER,ncomp,nvar,N,d)

  implicit none
  
  !Entrees/Sorties :
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncomp,nvar),&
                intent(inout)          :: field 
  integer,intent(in)                   :: ncomp,nvar,FILTER
  integer,dimension(3),intent(in),optional      :: N
  real(mytype),dimension(3),intent(in),optional :: d
  
  real(mytype),allocatable,dimension(:,:,:) :: Fh ! Hallo of Field   !TODO : utiliser halo_tmp (voir field_div)
  integer                                   :: i,j,k,l,m
  
  i=N(1)+int(d(1)) ! trick to avoid gcc-warning
  
  if (FILTER == 1) then
  !-----------------------------------Filtre MEAN 8 voisin Noeuds -> Centres
  do m=1,nvar
  do l=1,ncomp
  call update_halo(Field(:,:,:,l,m),Fh,level=1,opt_decomp=ph_decomp,opt_global=.true.) 

  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)-1
      field(i,j,k,l,m) = ( Fh(i  ,j  ,k  ) + &
                           Fh(i  ,j  ,k+1) + &
                           Fh(i  ,j+1,k  ) + &
                           Fh(i  ,j+1,k+1) + &
                           Fh(i+1,j  ,k  ) + &
                           Fh(i+1,j  ,k+1) + &
                           Fh(i+1,j+1,k  ) + &
                           Fh(i+1,j+1,k+1) )/8._mytype
  end do
  end do
  end do  
  !-- traitement specifique xend(1) car halo uniquement sur directions 2 et 3 
  i=xend(1)
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
        field(i,j,k,l,m) = ( Fh(i  ,j  ,k  ) + &
                           Fh(i  ,j  ,k+1) + &
                           Fh(i  ,j+1,k  ) + &
                           Fh(i  ,j+1,k+1) + &
                           Fh(1  ,j  ,k  ) + &
                           Fh(1  ,j  ,k+1) + &
                           Fh(1  ,j+1,k  ) + &
                           Fh(1  ,j+1,k+1) )/8._mytype
  end do
  end do  

  end do
  end do
  
  elseif (FILTER == -1) then
  !-----------------------------------Filtre MEAN 8 voisin Centres -> Noeuds
  do m=1,nvar
  do l=1,ncomp
  call update_halo(Field(:,:,:,l,m),Fh,level=1,opt_decomp=ph_decomp,opt_global=.true.) 

  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1)+1,xend(1)
       field(i,j,k,l,m) = ( Fh(i  ,j  ,k  ) + &
                            Fh(i  ,j  ,k-1) + &
                            Fh(i  ,j-1,k  ) + &
                            Fh(i  ,j-1,k-1) + &
                            Fh(i-1,j  ,k  ) + &
                            Fh(i-1,j  ,k-1) + &
                            Fh(i-1,j-1,k  ) + &
                            Fh(i-1,j-1,k-1) )/8._mytype
  end do
  end do
  end do  
  !-- traitement specifique xstart(1) car halo uniquement sur directions 2 et 3 
  i=xstart(1)
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
      field(i,j,k,l,m)=( Fh(i  ,j  ,k  ) + &
                            Fh(i  ,j  ,k-1) + &
                            Fh(i  ,j-1,k  ) + &
                            Fh(i  ,j-1,k-1) + &
                            Fh(xend(1),j  ,k  ) + &
                            Fh(xend(1),j  ,k-1) + &
                            Fh(xend(1),j-1,k  ) + &
                            Fh(xend(1),j-1,k-1) )/8._mytype
  end do
  end do  
  
  end do
  end do

  elseif (FILTER == 1000) then  
  !-----------------------------------Filtre SUM / 6 voisins (central voxel excluded)
  do m=1,nvar
  do l=1,ncomp
  call update_halo(Field(:,:,:,l,m),Fh,level=1,opt_decomp=ph_decomp,opt_global=.true.) 

  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1)+1,xend(1)-1
       field(i,j,k,l,m) =   Fh(i  ,j  ,k-1) + &
                            Fh(i  ,j  ,k+1) + &
                            Fh(i  ,j-1,k  ) + &
                            Fh(i  ,j+1,k  ) + &
                            Fh(i-1,j  ,k  ) + &
                            Fh(i+1,j  ,k  )
  end do
  end do
  end do  
  !-- traitement specifique xstart(1), xend(1) car halo uniquement sur directions 2 et 3 
  i=xstart(1)
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
       field(i,j,k,l,m) =   Fh(i  ,j  ,k-1) + &
                            Fh(i  ,j  ,k+1) + &
                            Fh(i  ,j-1,k  ) + &
                            Fh(i  ,j+1,k  ) + &
                            Fh(xend(1),j ,k) + &
                            Fh(i+1,j  ,k  )
  end do
  end do  
  i=xend(1)
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
       field(i,j,k,l,m) =   Fh(i  ,j  ,k-1) + &
                            Fh(i  ,j  ,k+1) + &
                            Fh(i  ,j-1,k  ) + &
                            Fh(i  ,j+1,k  ) + &
                            Fh(i-1,j  ,k  ) + &
                            Fh(1  ,j  ,k  )
  end do
  end do  
  
  end do
  end do
  
  else 
     call amitex_abort("field_filter : FILTER different from -1,1,1000 - not yet implemented " ,2,0)  
  end if
  

  if (allocated(Fh)) deallocate(Fh)

end subroutine field_filter


!==============================================================================
!       SUBROUTINE FIELD_DIV
!------------------------------------------------------------------------------
!> Evaluate div() in real space
!!
!!              field_div(fin,div,ncompin,ncompout,d,N,dtype,voigt_convention,BC)
!!
!!
!!   \param[in]  fin        field in (nxP,nyP,nzP,ncompin=3, 6 or 9)    (profil indifferent)
!!                          ncompin= 3 : vector field
!!                                   6 : symmetric 2nd order tensor
!!                                   9 : non-symmetric second order tensor
!!   \param[out] div        output field (nxP,nyP,nzP,ncompout=1 or 3)    (profil indifferent)
!!   \param[in]  ncompin    number of components of fin 
!!   \param[in]  ncompout   number of components of div 
!!   \param[in]  d          [dx,dy,dz] voxel sizes 
!!   \param[in]  N          [Nx,Ny,Nz] grid dimensions 
!!   \param[in]  dtype      derivation type 
!!   \param[in]  voigt_convention   optional "stress" or "strain" 
!!                                  required if ncompout=6 (symetric tensor)
!!                                  default = "stress"
!!   \param[in]  BC         optional [Ix, Iy, Iz] 
!!                          if I = 0 : periodic BC (default is periodic BC)
!!                             I = 1 : linear extrapolation
!!                                     (on plane 0, if dtype=-1, or plane N+1, if dtype=+1)
!!               For the moment only [0,0,0] or [1,1,1] are implemented
!!
!! dtype will be used to apply different derivative with FE-IR
!!         dtype = -1 => derivation from centers to nodes
!!         dtype = 1  => derivation from nodes to centers
!!         
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!==============================================================================
subroutine field_div(fin,div,ncompin,ncompout,d,N,dtype,voigt_convention,BC)

  implicit none
  
  ! Inputs/Outputs
  integer,intent(in)                   :: ncompin,ncompout
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncompin),&
                intent(in)             :: fin
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncompout),&
                intent(out)            :: div                                 
  real(mytype),dimension(3),intent(in) :: d
  integer,dimension(3),intent(in)      :: N
  integer,intent(in)                   :: dtype
  
  character(len=*),intent(in),optional :: voigt_convention
  integer, dimension(3), intent(in), optional :: BC
  
  ! Local Variables
  integer                              :: i,I1,J1,K1,i0,j0,k0,i00
  integer                              :: KSTART,KEND, JSTART,JEND,ISTART,IEND,DTRANS,IBORD,ITEST,I1mod,IP1mod,I0VAL
  real(mytype)                         :: dx,dy,dz, factor
  real(mytype),dimension(ncompout)     :: Floc 
  integer, dimension(3)                :: BC0
  real(mytype),allocatable,dimension(:,:,:) :: tmp_field ! to reallocate tmp_halo%var in X direction

  
  !--------------------------- CHECK INPUTS
  if (ncompin /=3 .AND. ncompin /= 6 .AND. ncompin /= 9) &
                 call amitex_abort("field_div : ncompin /= 3, 6 and 9",2,0)
  
  if (ncompin ==3 .AND. ncompout /= 1) &
                 call amitex_abort("field_div : ncompin=3 but ncompout /= 1",2,0)
                 
  if ((ncompin == 6 .OR. ncompin==9) .AND. ncompout /= 3) &
                 call amitex_abort("field_div : ncompin=6 or 9 but ncompout /= 3",2,0)
  
  if (dtype /= 1 .and. dtype /= -1) &
                 call amitex_abort("field_div : dtype different from 1 or -1",2,0)
                 
  if (ncompin == 6 .and. .not. present(voigt_convention)) &
                 call amitex_abort("field_div :  ncompin=6 but voigt_convention absent",2,0)
                                        
  factor = 1.                                    
  if (present(voigt_convention)) then
     if (voigt_convention == "stress") then
        factor = 1.
     elseif (voigt_convention == "strain") then
        factor = 2.
     else
        call amitex_abort("field_div:  voigt_convention different from 'stress' and 'strain'",2,0) 
     end if
  end if
 
  BC0 = 0 
  if (present(BC)) then
     BC0 = BC
     if (.not. all(BC==0) .and. .not. all(BC==1)) then
        call amitex_abort("field_grad :  BC different from [0,0,0] or [1,1,1], only available at the moment",2,0) 
     end if
  end if

  !--------------------------- INITIALIZATIONS
  dx = d(1)
  dy = d(2)
  dz = d(3)
  
  KSTART=0;KEND=0;JSTART=0;JEND=0;ISTART=0;IEND=0;DTRANS=0;IBORD=0;ITEST=0;I1mod=0;IP1mod=0;I0VAL=0
    
  !--------------------------- INITIALIZE halo_tmp(ncomp) with fin
  if (allocated(halo_tmp) .AND. size(halo_tmp) < ncompin) then
    deallocate(halo_tmp)
    allocate(halo_tmp(ncompin)) 
  end if
  if (.not. allocated(halo_tmp)) allocate(halo_tmp(ncompin)) 

  do i = 1,ncompin
     call update_halo(fin(:,:,:,i),halo_tmp(i)%var,level=1,opt_decomp=ph_decomp,opt_global=.true.)
  end do

  !---------------------------- TREATMENT FOR SYMETRIC TENSORS  
  if (ncompin==6 .and. voigt_convention=="strain") then
  do i = 4,6
     halo_tmp(i)%var=halo_tmp(i)%var / factor
  end do
  end if 

  !--------------------------- RE-ALLOCATE AND INITIALIZE tmp_halo, with X-BC
  allocate(tmp_field(xstart(1):xend(1),xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
  
  do i = 1,ncompin
     tmp_field = halo_tmp(i)%var
     deallocate(halo_tmp(i)%var)
     allocate(halo_tmp(i)%var(xstart(1)-1:xend(1)+1,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
     halo_tmp(i)%var(xstart(1):xend(1),xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1) = tmp_field
  end do
    
  if (all(BC0==0)) then
     do i = 1,ncompin
        halo_tmp(i)%var(xstart(1)-1,:,:) = halo_tmp(i)%var(xend(1),:,:)
        halo_tmp(i)%var(xend(1)+1,:,:)   = halo_tmp(i)%var(xstart(1),:,:)
     end do
  end if
  if (all(BC0==1)) then
     do i = 1,ncompin
        halo_tmp(i)%var(xstart(1)-1,:,:) = 2*halo_tmp(i)%var(xstart(1),:,:) - halo_tmp(i)%var(xstart(1)+1,:,:)
        halo_tmp(i)%var(xend(1)+1,:,:)   = 2*halo_tmp(i)%var(xend(1),:,:)   - halo_tmp(i)%var(xend(1)-1,:,:) 
     end do
  end if

  !--------------------------- DEFINE BC ON halo_tmp, for Y and Z-planes
  !if (all(BC0)==0) : nothing to do, case by default when initializing 2decomp
  if (all(BC0==1)) then
  do i = 1, ncompin
     if (dtype == 1) then
       if (xend(2) == N(2)) &
          halo_tmp(i)%var(:,N(2)+1,:) = 2 * halo_tmp(i)%var(:,N(2),:) - halo_tmp(i)%var(:,N(2)-1,:)
       if (xend(3) == N(3)) &       
          halo_tmp(i)%var(:,:,N(3)+1) = 2 * halo_tmp(i)%var(:,:,N(3)) - halo_tmp(i)%var(:,:,N(3)-1)
     else if (dtype == -1) then
       if (xstart(2) == 1) &
          halo_tmp(i)%var(:,0,:) = 2 * halo_tmp(i)%var(:,1,:) - halo_tmp(i)%var(:,2,:)
       if (xstart(3) == 1) &
          halo_tmp(i)%var(:,:,0) = 2 * halo_tmp(i)%var(:,:,1) - halo_tmp(i)%var(:,:,2)     
     end if
  end do
  end if
  
  !--------------------------- EVALUATE DIV(fin) - TODO : traitement sans 'boucle' I1,J1,K1
  !                                                       unification avec un vecteur comp=(/1 2 3 4 5 6 7 8 9/) en GD
  !                                                                                        (/1 2 3 4 5 6 4 5 6/)  
  ! 
  if (dtype == -1) then !default Centers -> nodes
     KSTART = xstart(3)-1
     KEND   = xend(3)-1
     JSTART = xstart(2)-1
     JEND   = xend(2)-1
     ISTART = xstart(1)-1
     IEND   = xend(1)-1
     DTRANS = 1
  elseif (dtype == 1) then ! nodes -> centers
     KSTART = xstart(3)
     KEND   = xend(3)
     JSTART = xstart(2)
     JEND   = xend(2)
     ISTART = xstart(1)
     IEND   = xend(1)
     DTRANS = 0     
  end if

  if (ncompin==3) then !---- Vector case
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    do I1=ISTART,IEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(I1+1,J1+j0,K1+k0) - halo_tmp(1)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(2)%var(I1+i0,J1+1,K1+k0) - halo_tmp(2)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         Floc(1) = Floc(1) + (halo_tmp(3)%var(I1+i0,J1+j0,K1+1) - halo_tmp(3)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      div(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
    end do
  elseif (ncompin==6) then !---- 2nd order symmetric tensor case
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    do I1=ISTART,IEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(I1+1,J1+j0,K1+k0) - halo_tmp(1)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(2) = Floc(2) + (halo_tmp(4)%var(I1+1,J1+j0,K1+k0) - halo_tmp(4)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(3) = Floc(3) + (halo_tmp(5)%var(I1+1,J1+j0,K1+k0) - halo_tmp(5)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(4)%var(I1+i0,J1+1,K1+k0) - halo_tmp(4)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i0,J1+1,K1+k0) - halo_tmp(2)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(3) = Floc(3) + (halo_tmp(6)%var(I1+i0,J1+1,K1+k0) - halo_tmp(6)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         Floc(1) = Floc(1) + (halo_tmp(5)%var(I1+i0,J1+j0,K1+1) - halo_tmp(5)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(2) = Floc(2) + (halo_tmp(6)%var(I1+i0,J1+j0,K1+1) - halo_tmp(6)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i0,J1+j0,K1+1) - halo_tmp(3)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      div(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
    end do
  elseif (ncompin==9) then !---- 2nd order NON-symmetric tensor case
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    do I1=ISTART,IEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(I1+1,J1+j0,K1+k0) - halo_tmp(1)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(2) = Floc(2) + (halo_tmp(7)%var(I1+1,J1+j0,K1+k0) - halo_tmp(7)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(3) = Floc(3) + (halo_tmp(8)%var(I1+1,J1+j0,K1+k0) - halo_tmp(8)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(4)%var(I1+i0,J1+1,K1+k0) - halo_tmp(4)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i0,J1+1,K1+k0) - halo_tmp(2)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(3) = Floc(3) + (halo_tmp(9)%var(I1+i0,J1+1,K1+k0) - halo_tmp(9)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         Floc(1) = Floc(1) + (halo_tmp(5)%var(I1+i0,J1+j0,K1+1) - halo_tmp(5)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(2) = Floc(2) + (halo_tmp(6)%var(I1+i0,J1+j0,K1+1) - halo_tmp(6)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i0,J1+j0,K1+1) - halo_tmp(3)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      div(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
    end do
  end if  
  
if (1 == 0) then ! shunt old treatment for X-direction  
  !--------------Traitement du cas particulier des bords IBORD car halo uniquement sur directions 2 et 3 
  if (dtype == -1) then !default Centers -> nodes
    IBORD = xstart(1)-1
    IP1mod = IBORD+1
    I1mod = xend(1)
    Itest = 0
    I0val = xend(1)
  elseif (dtype == 1) then ! nodes -> centers
    IBORD = xend(1)
    IP1mod = xstart(1)
    I1mod = IBORD
    Itest = 1
    I0val = xstart(1)
  end if

  if (ncompin==3) then !---- vector case
    I1 = IBORD
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(1)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(1) = Floc(1) + (halo_tmp(2)%var(I1+i00,J1+1,K1+k0) - halo_tmp(2)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(1) = Floc(1) + (halo_tmp(3)%var(I1+i00,J1+j0,K1+1) - halo_tmp(3)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      div(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
  elseif (ncompin==6) then !---- 2nd order symmetric tensor case
    I1 = IBORD
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(1)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(2) = Floc(2) + (halo_tmp(4)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(4)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(3) = Floc(3) + (halo_tmp(5)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(5)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(1) = Floc(1) + (halo_tmp(4)%var(I1+i00,J1+1,K1+k0) - halo_tmp(4)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i00,J1+1,K1+k0) - halo_tmp(2)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(3) = Floc(3) + (halo_tmp(6)%var(I1+i00,J1+1,K1+k0) - halo_tmp(6)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(1) = Floc(1) + (halo_tmp(5)%var(I1+i00,J1+j0,K1+1) - halo_tmp(5)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(2) = Floc(2) + (halo_tmp(6)%var(I1+i00,J1+j0,K1+1) - halo_tmp(6)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i00,J1+j0,K1+1) - halo_tmp(3)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      div(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
  elseif (ncompin==9) then !---- 2nd order NON-symmetric tensor case
    I1 = IBORD
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(1)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(2) = Floc(2) + (halo_tmp(7)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(7)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(3) = Floc(3) + (halo_tmp(8)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(8)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(1) = Floc(1) + (halo_tmp(4)%var(I1+i00,J1+1,K1+k0) - halo_tmp(4)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i00,J1+1,K1+k0) - halo_tmp(2)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(3) = Floc(3) + (halo_tmp(9)%var(I1+i00,J1+1,K1+k0) - halo_tmp(9)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(1) = Floc(1) + (halo_tmp(5)%var(I1+i00,J1+j0,K1+1) - halo_tmp(5)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(2) = Floc(2) + (halo_tmp(6)%var(I1+i00,J1+j0,K1+1) - halo_tmp(6)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i00,J1+j0,K1+1) - halo_tmp(3)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      div(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
  end if  
end if ! shunt old treatment in X-direction  
end subroutine field_div

!===========================================================================================================
!       SUBROUTINE FIELD_DIVF
!------------------------------------------------------------------------------
!> Evaluate div() in Fourier space
!!
!!              field_divF(finF,divF,ncompin,ncompout,FREQ,FREQ_2,dtype,voigt_convention)
!!
!!   \param[in]  finF       field in Fourier (.,.,.,ncompin=3, 6 or 9)    (profil indifferent)
!!                          ncompin= 3 : vector field
!!                                 6 : symmetric 2nd order tensor
!!                                 9 : non-symmetric second order tensor
!!   \param[out] divF       divergence field, size (.,.,.,ncompout=1 or 3) (profil indifferent)
!!   \param[in]  ncompin    number of component of finF
!!   \param[in]  ncompout   number of component of divF 
!!   \param[in]  FREQ       Frequencies real(.,.,.,3) 
!!   \param[in]  FREQ_2     Half-voxel translation in Fourier complex(.,.,.) 
!!   \param[in]  dtype      derivation type 
!!   \param[in]  voigt_convention (optional) required for ncompin=6
!!                          character "stress" or "strain"
!!   
!!
!! dtype will be used to apply different derivative with EF_IR
!!         dtype = -1 => derivation from centers to nodes
!!         dtype = 1  => derivation from nodes to centers
!!         
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!===========================================================================================================
subroutine field_divF(finF,divF,ncompin,ncompout,FREQ,FREQ_2,dtype,voigt_convention)
  implicit none

!!------------------------------------------------------------------------------
!>                                                                        IN/OUT
  integer, intent(in)      :: ncompin,ncompout
  complex(mytype), dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncompin),&
            intent(in)     :: finF
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
            intent(in)     :: FREQ
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)), &
            intent(in)     :: FREQ_2
  complex(mytype), dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncompout),&
               intent(out) :: divF                            
  integer,intent(in)       :: dtype
  character(len=*),intent(in),optional :: voigt_convention
 
  
!!------------------------------------------------------------------------------
!>                                                                        OTHERS
  integer               :: i,j,h,l
  real(mytype)          :: factor
  complex(mytype)       :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i

!!------------------------------------------------------------------------------
!>                                                               INITIALIZATIONS

  
!!------------------------------------------------------------------------------
!>                                                                   CHECK INPUT

  if (ncompin /=3 .AND. ncompin /= 6 .AND. ncompin /= 9) &
                 call amitex_abort("field_divF : ncompin /= 3, 6 and 9",2,0)
  
  if (ncompin ==3 .AND. ncompout /= 1) &
                 call amitex_abort("field_divF : ncompin=3 but ncompout /= 1",2,0)
                 
  if ((ncompin == 6 .OR. ncompin==9) .AND. ncompout /= 3) &
                 call amitex_abort("field_divF : ncompin=6 or 9 but ncompout /= 3",2,0)
  
  if (dtype /= 1 .and. dtype /= -1) &
                 call amitex_abort("field_divF : dtype different from 1 or -1",2,0)
                 
  if (ncompin == 6 .and. .not. present(voigt_convention)) &
                 call amitex_abort("field_divF :  ncompin=6 but voigt_convention absent",2,0)
                                        
  factor = 1.                                    
  if (present(voigt_convention)) then
     if (voigt_convention == "stress") then
        factor = 1.
     elseif (voigt_convention == "strain") then
        factor = 2.
     else
        call amitex_abort("field_divF:  voigt_convention different from 'stress' and 'strain'",2,0) 
     end if
  end if
  
!!------------------------------------------------------------------------------
!>                                                                    DIVERGENCE

  if (ncompin==3) then
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  divF(h,j,i,1)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,1) +FREQ(h,j,i,2)* finF(h,j,i,2) + FREQ(h,j,i,3)* finF(h,j,i,3) )

  end do
  end do
  end do

  elseif (ncompin==6) then
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  !divF(h,j,i,1)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,1) +FREQ(h,j,i,2)* finF(h,j,i,4) + FREQ(h,j,i,3)* finF(h,j,i,5) )
  !divF(h,j,i,2)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,4) +FREQ(h,j,i,2)* finF(h,j,i,2) + FREQ(h,j,i,3)* finF(h,j,i,6) )
  !divF(h,j,i,3)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,5) +FREQ(h,j,i,2)* finF(h,j,i,6) + FREQ(h,j,i,3)* finF(h,j,i,3) )
  divF(h,j,i,1)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,1) + (FREQ(h,j,i,2)* finF(h,j,i,4) + FREQ(h,j,i,3)* finF(h,j,i,5))/factor )
  divF(h,j,i,2)= imP * (  FREQ(h,j,i,2)* finF(h,j,i,2) + (FREQ(h,j,i,1)* finF(h,j,i,4) + FREQ(h,j,i,3)* finF(h,j,i,6))/factor )
  divF(h,j,i,3)= imP * (  FREQ(h,j,i,3)* finF(h,j,i,3) + (FREQ(h,j,i,2)* finF(h,j,i,6) + FREQ(h,j,i,1)* finF(h,j,i,5))/factor )

  end do
  end do
  end do
  
  elseif (ncompin==9) then
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  divF(h,j,i,1)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,1) +FREQ(h,j,i,2)* finF(h,j,i,4) + FREQ(h,j,i,3)* finF(h,j,i,5) )
  divF(h,j,i,2)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,7) +FREQ(h,j,i,2)* finF(h,j,i,2) + FREQ(h,j,i,3)* finF(h,j,i,6) )
  divF(h,j,i,3)= imP * (  FREQ(h,j,i,1)* finF(h,j,i,8) +FREQ(h,j,i,2)* finF(h,j,i,9) + FREQ(h,j,i,3)* finF(h,j,i,3) )

  end do
  end do
  end do
 
  end if
!!------------------------------------------------------------------------------
!>                                                        HALF-VOXEL TRANSLATION

  if (dtype == -1) then
  do l=1,ncompout
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  divF(h,j,i,l)= divF(h,j,i,l) / FREQ_2(h,j,i) 

  end do
  end do
  end do
  end do
  elseif (dtype == 1) then
  do l=1,ncompout
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  divF(h,j,i,l)= divF(h,j,i,l) * FREQ_2(h,j,i) 

  end do
  end do
  end do
  end do
  end if  
  
end subroutine field_divF

!===========================================================================================================
!       SUBROUTINE FIELD_GRADF
!------------------------------------------------------------------------------
!> Evaluate div() in Fourier space
!!
!!              field_gradF(finF,gradF,ncompin,ncompout,FREQ,FREQ_2,dtype,voigt_convention)
!!
!!   \param[in]  finF       field in Fourier (.,.,.,ncomp=3, 6 or 9)    (profil indifferent)
!!                          ncompin= 1 : scalar field
!!                          ncompin= 3 : vector field
!!   \param[out] gradF      gradient size (.,.,.,ncompout) (profil indifferent)
!!                          ncompout = 3, 6 or 9
!!   \param[in]  ncompin    Number of components for input field
!!   \param[in]  ncompout   Number of components for output field
!!   \param[in]  FREQ       Frequencies real(.,.,.,3) 
!!   \param[in]  FREQ_2     Half-voxel translation in Fourier complex(.,.,.) 
!!   \param[in]  dtype      derivation type 
!!   \patam[in]  voigt_convention    "stress" or "strain" 
!!                          required if ncompout=6 (symetric tensor)
!!                          default is "stress"
!! dtype will be used to apply different derivative with EF_IR
!!         dtype = 1   => derivation from nodes to centers
!!         dtype = -1  => derivation from centers to nodes
!!         
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!===========================================================================================================
subroutine field_gradF(finF,gradF,ncompin,ncompout,FREQ,FREQ_2,dtype,voigt_convention)
  implicit none

!!------------------------------------------------------------------------------
!>                                                                        IN/OUT
  integer, intent(in)           :: ncompin, ncompout
  complex(mytype), dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncompin),&
            intent(in)          :: finF
  real(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),1:3), &
            intent(in)          :: FREQ
  complex(mytype),dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3)), &
            intent(in)          :: FREQ_2
  complex(mytype), dimension(fft_start(1):fft_end(1),fft_start(2):fft_end(2),fft_start(3):fft_end(3),ncompout),&
               intent(out)      :: gradF                            
  integer,intent(in)            :: dtype
  character(len=*),intent(in),optional :: voigt_convention
 
  
!!------------------------------------------------------------------------------
!>                                                                        OTHERS
  real(mytype)          :: factor
  integer               :: i,j,h,l
  complex(mytype)       :: imP = cmplx(0,1,kind=mytype) ! le nombre complexe i

!!------------------------------------------------------------------------------
!>                                                               INITIALIZATIONS

  factor = 0.
    
!!------------------------------------------------------------------------------
!>                                                                   CHECK INPUT

  if (ncompin /= 3 .and. ncompin /= 1) call amitex_abort("field_gradF :  ncompin different from 1 and 3",2,0)

  if (ncompin == 1 .and. ncompout /= 3) call amitex_abort("field_gradF :  ncompin=1 but ncompout/=3",2,0)

  if (ncompin == 3 .and. (ncompout /= 6 .and. ncompout /=9)) &
                                        call amitex_abort("field_gradF :  ncompin=3 but ncompout/=6 and 9",2,0)

  if (dtype /= 1 .and. dtype /= -1) call amitex_abort("field_gradF : dtype different from 1 or -1",2,0)
   
  if (ncompout == 6 .and. .not. present(voigt_convention)) &
                                        call amitex_abort("field_gradF :  ncompout=6 but voigt_convention absent",2,0)
                                        
  factor = 1.
  if (present(voigt_convention)) then
     if (voigt_convention == "stress") then
        factor = 1.
     elseif (voigt_convention == "strain") then
        factor = 2.
     else
        call amitex_abort("field_gradF :  voigt_convention different fron 'stress' and 'strain'",2,0) 
     end if
  end if
 
!!------------------------------------------------------------------------------
!>                                                                      GRADIENT

  if (ncompin==1) then ! input = scalar, output = vector
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  gradF(h,j,i,1)= imP * finF(h,j,i,1) *  FREQ(h,j,i,1) 
  gradF(h,j,i,2)= imP * finF(h,j,i,1) *  FREQ(h,j,i,2) 
  gradF(h,j,i,3)= imP * finF(h,j,i,1) *  FREQ(h,j,i,3) 

  end do
  end do
  end do

  elseif (ncompin==3 .and. ncompout==6) then  ! input = vector, output = symmetrized tensor 2 (6 comp)
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  gradF(h,j,i,1)= imP * finF(h,j,i,1) *  FREQ(h,j,i,1) 
  gradF(h,j,i,2)= imP * finF(h,j,i,2) *  FREQ(h,j,i,2) 
  gradF(h,j,i,3)= imP * finF(h,j,i,3) *  FREQ(h,j,i,3) 

  gradF(h,j,i,4)= imP * ( finF(h,j,i,1) *  FREQ(h,j,i,2) + finF(h,j,i,2) *  FREQ(h,j,i,1)) * (factor/2._mytype)
  gradF(h,j,i,5)= imP * ( finF(h,j,i,1) *  FREQ(h,j,i,3) + finF(h,j,i,3) *  FREQ(h,j,i,1)) * (factor/2._mytype)
  gradF(h,j,i,6)= imP * ( finF(h,j,i,2) *  FREQ(h,j,i,3) + finF(h,j,i,3) *  FREQ(h,j,i,2)) * (factor/2._mytype)

  end do
  end do
  end do
  
  elseif (ncompin==3 .and. ncompout==9) then ! input = vector, output non symmetric tensor (9 comp)
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  gradF(h,j,i,1)= imP * finF(h,j,i,1) *  FREQ(h,j,i,1)
  gradF(h,j,i,2)= imP * finF(h,j,i,2) *  FREQ(h,j,i,2)
  gradF(h,j,i,3)= imP * finF(h,j,i,3) *  FREQ(h,j,i,3)

  gradF(h,j,i,4)= imP * finF(h,j,i,1) *  FREQ(h,j,i,2)
  gradF(h,j,i,5)= imP * finF(h,j,i,1) *  FREQ(h,j,i,3) 
  gradF(h,j,i,6)= imP * finF(h,j,i,2) *  FREQ(h,j,i,3)

  gradF(h,j,i,7)= imP * finF(h,j,i,2) *  FREQ(h,j,i,1)
  gradF(h,j,i,8)= imP * finF(h,j,i,3) *  FREQ(h,j,i,1)
  gradF(h,j,i,9)= imP * finF(h,j,i,3) *  FREQ(h,j,i,2)

  end do
  end do
  end do
 
  end if
!!------------------------------------------------------------------------------
!>                                                        HALF-VOXEL TRANSLATION

  if (dtype == -1) then
  do l=1,ncompout
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  gradF(h,j,i,l)= gradF(h,j,i,l) / FREQ_2(h,j,i) 

  end do
  end do
  end do
  end do
  elseif (dtype == 1) then
  do l=1,ncompout
  do i=fft_start(3),fft_end(3)
  do j=fft_start(2),fft_end(2)
  do h=fft_start(1),fft_end(1)

  gradF(h,j,i,l)= gradF(h,j,i,l) * FREQ_2(h,j,i) 

  end do
  end do
  end do
  end do
  end if  
  
end subroutine field_gradF

!==============================================================================
!       SUBROUTINE FIELD_GRAD
!------------------------------------------------------------------------------
!> Evaluate grad() in real space
!!
!!              field_grad(fin,grad,ncompin,ncompout,dtype,voigt_convention)
!!
!!   \param[in]  fin       field in real (.,.,.,ncomp=1 or 3 )    (profil indifferent)
!!                          ncompin= 1 : scalar field
!!                          ncompin= 3 : vector field
!!   \param[out] grad      gradient size (.,.,.,ncompout) (profil indifferent)
!!                          ncompout = 3, 6 or 9
!!   \param[in]  ncompin    Number of components for input field
!!   \param[in]  ncompout   Number of components for output field
!!   \param[in]  dtype      derivation type 
!!   \param[in]  d          [dx,dy,dz] voxel sizes 
!!   \param[in]  N          [Nx,Ny,Nz] grid dimensions 
!!   \param[in]  voigt_convention   optional "stress" or "strain" 
!!                                  required if ncompout=6 (symetric tensor)
!!                                  default = "stress"
!!   \param[in]  BC         optional [Ix, Iy, Iz] 
!!                          if I = 0 : periodic BC (default is periodic BC)
!!                             I = 1 : linear extrapolation
!!                                     (on plane 0, if dtype=-1, or plane N+1, if dtype=+1)
!!               For the moment only [0,0,0] or [1,1,1] are implemented
!!
!! dtype will be used to apply different derivative with EF_IR
!!         dtype = 1   => derivation from nodes to centers
!!         dtype = -1  => derivation from centers to nodes
!!         
!!
!! IMPORTANT : ON UTILISE ICI LE PASSAGE D'ARGUMENTS PAR ADRESSE SANS PROFIL IMPLICITE
!!             POUR TRAITER DES PROFILS DE TABLEAUX DIFFERENTS (fieldF et D_fieldF)
!!
!==============================================================================
subroutine field_grad(fin,grad,ncompin,ncompout,d,N,dtype,voigt_convention,BC)

  implicit none
  
  ! Inputs/Outputs
  integer,intent(in)                   :: ncompin,ncompout
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncompin),&
                intent(in)             :: fin
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncompout),&
                intent(out)            :: grad                                
  real(mytype),dimension(3),intent(in) :: d
  integer,dimension(3),intent(in)      :: N
  integer,intent(in)                   :: dtype
  
  character(len=*),intent(in),optional     :: voigt_convention
  integer,dimension(3),intent(in),optional :: BC
  
  ! Local Variables
  integer                              :: i,I1,J1,K1,i0,j0,k0,i00
  integer                              :: KSTART,KEND, JSTART,JEND,ISTART,IEND,DTRANS,IBORD,ITEST,I1mod,IP1mod,I0VAL
  real(mytype)                         :: dx,dy,dz
  real(mytype),dimension(ncompout)     :: Floc 
  real(mytype)                         :: factor
  real(mytype),allocatable,dimension(:,:,:) :: tmp_field ! to reallocate tmp_halo%var in X direction
  integer,dimension(3)                 :: BC0
  

    
!!------------------------------------------------------------------------------
!>                                                                   CHECK INPUT

  if (ncompin /= 3 .and. ncompin /= 1) call amitex_abort("field_grad :  ncompin different from 1 and 3",2,0)

  if (ncompin == 1 .and. ncompout /= 3) call amitex_abort("field_grad :  ncompin=1 but ncompout/=3",2,0)

  if (ncompin == 3 .and. (ncompout /= 6 .and. ncompout /=9)) &
                                        call amitex_abort("field_grad :  ncompin=3 but ncompout/=6 and 9",2,0)

  if (dtype /= 1 .and. dtype /= -1) call amitex_abort("field_grad : dtype different from 1 or -1",2,0)
   
  if (ncompout == 6 .and. .not. present(voigt_convention)) &
                                        call amitex_abort("field_grad :  ncompout=6 but voigt_convention absent",2,0)
  factor = 1.                                    
  if (present(voigt_convention)) then
     if (voigt_convention == "stress") then
        factor = 1.
     elseif (voigt_convention == "strain") then
        factor = 2.
     else
        call amitex_abort("field_grad :  voigt_convention different fron 'stress' and 'strain'",2,0) 
     end if
  end if

  BC0 = 0 
  if (present(BC)) then
     BC0 = BC
     if (.not. all(BC==0) .and. .not. all(BC==1)) then
        call amitex_abort("field_grad :  BC different from [0,0,0] or [1,1,1], only available at the moment",2,0) 
     end if
  end if

  !--------------------------- INITIALIZATIONS
  
  dx = d(1)
  dy = d(2)
  dz = d(3)
  
  KSTART=0;KEND=0;JSTART=0;JEND=0;ISTART=0;IEND=0;DTRANS=0;IBORD=0;ITEST=0;I1mod=0;IP1mod=0;I0VAL=0
  
  !--------------------------- INITIALIZE halo_tmp(ncompin) with fin
  if (allocated(halo_tmp) .AND. size(halo_tmp) < ncompin) then 
    deallocate(halo_tmp)
    allocate(halo_tmp(ncompin)) 
  end if
  if (.not. allocated(halo_tmp)) allocate(halo_tmp(ncompin)) 

  do i = 1,ncompin
     call update_halo(fin(:,:,:,i),halo_tmp(i)%var,level=1,opt_decomp=ph_decomp,opt_global=.true.) 
  end do
    
  !--------------------------- RE-ALLOCATE AND INITIALIZE tmp_halo, with X-BC
  allocate(tmp_field(xstart(1):xend(1),xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
  
  do i = 1,ncompin
     tmp_field = halo_tmp(i)%var
     deallocate(halo_tmp(i)%var)
     allocate(halo_tmp(i)%var(xstart(1)-1:xend(1)+1,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
     halo_tmp(i)%var(xstart(1):xend(1),xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1) = tmp_field
  end do
    
  if (all(BC0==0)) then
     do i = 1,ncompin
        halo_tmp(i)%var(xstart(1)-1,:,:) = halo_tmp(i)%var(xend(1),:,:)
        halo_tmp(i)%var(xend(1)+1,:,:)   = halo_tmp(i)%var(xstart(1),:,:)
     end do
  end if
  if (all(BC0==1)) then
     do i = 1,ncompin
        halo_tmp(i)%var(xstart(1)-1,:,:) = 2*halo_tmp(i)%var(xstart(1),:,:) - halo_tmp(i)%var(xstart(1)+1,:,:)
        halo_tmp(i)%var(xend(1)+1,:,:)   = 2*halo_tmp(i)%var(xend(1),:,:)   - halo_tmp(i)%var(xend(1)-1,:,:) 
     end do
  end if


  !--------------------------- DEFINE BC ON halo_tmp, for Y and Z-planes
  !if (all(BC0)==0) : nothing to do, case by default when initializing 2decomp
  if (all(BC0==1)) then
  do i = 1, ncompin
     if (dtype == 1) then
       if (xend(2) == N(2)) &
          halo_tmp(i)%var(:,N(2)+1,:) = 2 * halo_tmp(i)%var(:,N(2),:) - halo_tmp(i)%var(:,N(2)-1,:)
       if (xend(3) == N(3)) &       
          halo_tmp(i)%var(:,:,N(3)+1) = 2 * halo_tmp(i)%var(:,:,N(3)) - halo_tmp(i)%var(:,:,N(3)-1)
     else if (dtype == -1) then
       if (xstart(2) == 1) &
          halo_tmp(i)%var(:,0,:) = 2 * halo_tmp(i)%var(:,1,:) - halo_tmp(i)%var(:,2,:)
       if (xstart(3) == 1) &
          halo_tmp(i)%var(:,:,0) = 2 * halo_tmp(i)%var(:,:,1) - halo_tmp(i)%var(:,:,2)     
     end if
  end do
  end if
  
  !--------------------------- EVALUATE GRAD(fin) - TODO : traitement sans 'boucle' I1,J1,K1
  !                                                       unification avec un vecteur comp=(/1 2 3 4 5 6 7 8 9/) en GD
  !                                                                                        (/1 2 3 4 5 6 4 5 6/)  
  ! 
  if (dtype == -1) then !default Centers -> nodes
     KSTART = xstart(3)-1
     KEND   = xend(3)-1
     JSTART = xstart(2)-1
     JEND   = xend(2)-1
     ISTART = xstart(1)-1
     IEND   = xend(1)-1
     DTRANS = 1
  elseif (dtype == 1) then ! nodes -> centers
     KSTART = xstart(3)
     KEND   = xend(3)
     JSTART = xstart(2)
     JEND   = xend(2)
     ISTART = xstart(1)
     IEND   = xend(1)
     DTRANS = 0     
  end if

  if (ncompin==1) then !---- IN scalar OUT vector
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    do I1=ISTART,IEND
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(I1+1,J1+j0,K1+k0) - halo_tmp(1)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         Floc(2) = Floc(2) + (halo_tmp(1)%var(I1+i0,J1+1,K1+k0) - halo_tmp(1)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         Floc(3) = Floc(3) + (halo_tmp(1)%var(I1+i0,J1+j0,K1+1) - halo_tmp(1)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      grad(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    end do
    end do
    end do
  elseif (ncompin==3 .AND. ncompout==6) then !---- IN vector, OUT 2nd order symmetric tensor case
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    do I1=ISTART,IEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(I1+1,J1+j0,K1+k0) - halo_tmp(1)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(4) = Floc(4) + (halo_tmp(2)%var(I1+1,J1+j0,K1+k0) - halo_tmp(2)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(5) = Floc(5) + (halo_tmp(3)%var(I1+1,J1+j0,K1+k0) - halo_tmp(3)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         Floc(4) = Floc(4) + (halo_tmp(1)%var(I1+i0,J1+1,K1+k0) - halo_tmp(1)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i0,J1+1,K1+k0) - halo_tmp(2)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(6) = Floc(6) + (halo_tmp(3)%var(I1+i0,J1+1,K1+k0) - halo_tmp(3)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         Floc(5) = Floc(5) + (halo_tmp(1)%var(I1+i0,J1+j0,K1+1) - halo_tmp(1)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(6) = Floc(6) + (halo_tmp(2)%var(I1+i0,J1+j0,K1+1) - halo_tmp(2)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i0,J1+j0,K1+1) - halo_tmp(3)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      Floc(4:6) = Floc(4:6) * 0.5_mytype * factor
      grad(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
    end do
  elseif (ncompin==3 .AND. ncompout==9) then !---- IN Vector OUT : 2nd order NON-symmetric tensor case
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    do I1=ISTART,IEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(I1+1,J1+j0,K1+k0) - halo_tmp(1)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(7) = Floc(7) + (halo_tmp(2)%var(I1+1,J1+j0,K1+k0) - halo_tmp(2)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(8) = Floc(8) + (halo_tmp(3)%var(I1+1,J1+j0,K1+k0) - halo_tmp(3)%var(I1,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         Floc(4) = Floc(4) + (halo_tmp(1)%var(I1+i0,J1+1,K1+k0) - halo_tmp(1)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i0,J1+1,K1+k0) - halo_tmp(2)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
         Floc(9) = Floc(9) + (halo_tmp(3)%var(I1+i0,J1+1,K1+k0) - halo_tmp(3)%var(I1+i0,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         Floc(5) = Floc(5) + (halo_tmp(1)%var(I1+i0,J1+j0,K1+1) - halo_tmp(1)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(6) = Floc(6) + (halo_tmp(2)%var(I1+i0,J1+j0,K1+1) - halo_tmp(2)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i0,J1+j0,K1+1) - halo_tmp(3)%var(I1+i0,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      grad(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
    end do
  end if  

if (1==0) then  ! old version for PBC, TODO remove when everything OK 
  !--------------Traitement du cas particulier des bords IBORD car halo uniquement sur directions 2 et 3 
  if (dtype == -1) then !default Centers -> nodes
    IBORD = xstart(1)-1
    IP1mod = IBORD+1
    I1mod = xend(1)
    Itest = 0
    I0val = xend(1)
  elseif (dtype == 1) then ! nodes -> centers
    IBORD = xend(1)
    IP1mod = xstart(1)
    I1mod = IBORD
    Itest = 1
    I0val = xstart(1)
  end if

  if (ncompin==1) then !---- IN scalar OUT vector
    I1 = IBORD
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(1)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(2) = Floc(2) + (halo_tmp(1)%var(I1+i00,J1+1,K1+k0) - halo_tmp(1)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(3) = Floc(3) + (halo_tmp(1)%var(I1+i00,J1+j0,K1+1) - halo_tmp(1)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      grad(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
  elseif (ncompin==3 .AND. ncompout==6) then !---- IN vector OUT 2nd order symmetric tensor case
    I1 = IBORD
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(1)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(4) = Floc(4) + (halo_tmp(2)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(2)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(5) = Floc(5) + (halo_tmp(3)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(3)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(4) = Floc(4) + (halo_tmp(1)%var(I1+i00,J1+1,K1+k0) - halo_tmp(1)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i00,J1+1,K1+k0) - halo_tmp(2)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(6) = Floc(6) + (halo_tmp(3)%var(I1+i00,J1+1,K1+k0) - halo_tmp(3)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(5) = Floc(5) + (halo_tmp(1)%var(I1+i00,J1+j0,K1+1) - halo_tmp(1)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(6) = Floc(6) + (halo_tmp(2)%var(I1+i00,J1+j0,K1+1) - halo_tmp(2)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i00,J1+j0,K1+1) - halo_tmp(3)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      Floc(4:6) = Floc(4:6) * 0.5_mytype * factor
      grad(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
  elseif (ncompin==3 .AND. ncompout==9) then !---- IN vector OUT 2nd order NON-symmetric tensor case
    I1 = IBORD
    do K1=KSTART,KEND
    do J1=JSTART,JEND
    !do i = 1, npts
      Floc = 0._mytype
      !I1=IJK(i,1)-1;J1=IJK(i,2)-1;K1=IJK(i,3)-1;
      do j0=0,1
      do k0=0,1
         Floc(1) = Floc(1) + (halo_tmp(1)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(1)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(7) = Floc(7) + (halo_tmp(2)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(2)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
         Floc(8) = Floc(8) + (halo_tmp(3)%var(IP1mod,J1+j0,K1+k0) - halo_tmp(3)%var(I1mod,J1+j0,K1+k0)) / (4._mytype*dx)
      end do
      end do
      do i0=0,1
      do k0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(4) = Floc(4) + (halo_tmp(1)%var(I1+i00,J1+1,K1+k0) - halo_tmp(1)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(2) = Floc(2) + (halo_tmp(2)%var(I1+i00,J1+1,K1+k0) - halo_tmp(2)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
         Floc(9) = Floc(9) + (halo_tmp(3)%var(I1+i00,J1+1,K1+k0) - halo_tmp(3)%var(I1+i00,J1,K1+k0)) / (4._mytype*dy)
      end do
      end do      
      do i0=0,1
      do j0=0,1
         i00 = i0
         if (i00==Itest) i00 = I0val-I1
         Floc(5) = Floc(5) + (halo_tmp(1)%var(I1+i00,J1+j0,K1+1) - halo_tmp(1)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(6) = Floc(6) + (halo_tmp(2)%var(I1+i00,J1+j0,K1+1) - halo_tmp(2)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
         Floc(3) = Floc(3) + (halo_tmp(3)%var(I1+i00,J1+j0,K1+1) - halo_tmp(3)%var(I1+i00,J1+j0,K1)) / (4._mytype*dz)
      end do
      end do
      grad(I1+DTRANS,J1+DTRANS,K1+DTRANS,:) = Floc
    !end do
    end do
    end do
  end if  
end if !shunt traitement des bords  
end subroutine field_grad

 
!===============================================================================
!> field_Mean_box(field,ncomp,IJKmin,IJKmax,fieldMean)
!!  Compute the mean value of a field inside a region (box)
!!
!!   The box boundary is included in the computation
!!   
!!   \param[in]     field         input field (nx*ny*nz,ncomp)
!!   \param[in]     ncomp         number of components in the input field
!!   \param[in]     IJKmin        min corner indexes of the box, integer(3)
!!   \param[in]     IJKmax        max corner indexes of the box, integer(3)
!!   \param[out]    fieldMean     mean of the field inside the box 
!!
!================================================================================
subroutine field_Mean_box(field,ncomp,IJKmin,IJKmax,fieldMean)
  implicit none
  !!------------------------------------------------------------------------------
!>                                                                 ENTREE/SORTIE
  integer,intent(in)                                                     :: ncomp 
  integer, dimension(3),intent(in)                                       :: IJKmin,IJKmax
   ! reshape par le passage par adresse (norme 2003)
  real(mytype),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),ncomp),intent(in)    :: field
  real(mytype),dimension(ncomp),intent(out)                              :: fieldMean
  
  !!------------------------------------------------------------------------------
!>                                                             VARIABLES LOCALES
  integer(kind=INT64)                                                    :: Npts 
  integer                                                                :: i,j,k    ! loop index
  integer                                                                :: ierror   ! error for MPI subroutine
  
  !-------------------------------------------CHECK INPUTS
  if (minval(IJKmin) < 1) call amitex_abort("field_Mean_box : input minval(IJKmin) <0",2,0)
  if (minval(IJKmax) < 1) call amitex_abort("field_Mean_box : input minval(IJKmax) <0",2,0)

  !-------------------------------------------NUMBER OF POINTS IN THE BOX
  Npts = (IJKmax(1)-IJKmin(1)+1)*(IJKmax(2)-IJKmin(2)+1)*(IJKmax(3)-IJKmin(3)+1)
  
  !-------------------------------------------MEAN COMPUTATION

  fieldMean = 0._mytype
  
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)
     if (i >= IJKmin(1) .AND. j >= IJKmin(2) .AND. k >= IJKmin(3) .AND. &
         i <= IJKmax(1) .AND. j <= IJKmax(2) .AND. k <= IJKmax(3)) then
         fieldMean = fieldMean + field(i,j,k,:)
     end if 
  end do
  end do
  end do
  
  
  call MPI_ALLREDUCE(MPI_IN_PLACE,fieldMean, ncomp, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)

  fieldMean = fieldMean/Npts

end subroutine field_Mean_box



end module field_mod
