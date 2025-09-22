!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE STANDARD_USER_MOD 
!>
!>      Module for the standard algorithm : implement BEFORE and AFTER unpas user defined procedures 
!>
!!
!!
module standard_user_mod


  use ISO_FORTRAN_ENV


! MPI AND 2DECOMP MODULES
!------------------------
  use MPI             
  use decomp_2d
  use decomp_2d_fft

! ALL AMITEX MODULES (complete list of amitex modules)
!------------------------------------------------------
  use algo_functions_mod  
  use amitex_mod
  use amitex_user_mod     
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
  !!use resolution_mod       this module uses standard_user_mod
  !!use resolution_user_mod  no need
  use sortie_std_mod
  !!use standard_user_mod    this file
  
#ifdef OPENMP
  use omp_lib
#endif

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

  !trick to avoid gcc-warning
  integer :: bidon
  bidon = load_n+load_incr+ind_tps

end subroutine before_unpas_user


!===========================================================================================================
!!==========================================================================================================
!> after_unpas_user : write what must be done just AFTER unpas (using the standard algorithm)
!!                       
!!                       
!!==========================================================================================================
subroutine after_unpas_user(load_n,load_incr,ind_tps)
  
  implicit none

  integer,intent(in)     :: load_n,load_incr, ind_tps

  real(mytype)           :: dt, t1
  integer                :: ierror
  integer                :: i
  character(len=200)     :: tmp_char2,tmp_char1,tmp_char   !< variable caractere temporaire

  !-- Calcul du pas de temps
  if (load_incr==0) then
     dt = load(load_n)%time(load_incr)
  else
     dt = load(load_n)%time(load_incr) -  load(load_n)%time(load_incr-1)
  end if


  !================================================================================
  !                                                                      SORTIE VTK
  !             TODO AMITEX - mettre en place une structure pour analyser les temps 

  if(load_incr/=0) then ! Pas de sortie vtk pour les pas de temps fictifs
    if(extract%tpsVTK(ind_tps))then
      call mpi_barrier(mpi_comm_world,ierror)
      t1 = MPI_Wtime()

      ! Sortie vtk pour les variables de champ de phase
      !---------------------------------------
      ! a convertir si + de 1 PF
      !Sortie vtk des champs de phase
      write(tmp_char1,"(I4)") ind_tps
      do i=1,Neta
        write(tmp_char2,"(I4)") i
        tmp_char="_eta"//trim(adjustl(tmp_char2))//"_"//trim(adjustl(tmp_char1))
   
        call print_field_vtk(Eta0(:,i),trim(fic_vtk)//trim(tmp_char)//".vtk","eta"//trim(adjustl(tmp_char2)))
        !call write_header_vtk(trim(fic_vtk)//trim(tmp_char)//".vtk",1000, &
        !      grid%nx,grid%ny,grid%nz,grid%dx,grid%dy,grid%dz)
        !TEMPfield32=real(eta0(:,i))
        !call write_bin(1,TEMPfield32,"eta",trim(fic_vtk)//trim(tmp_char)//".vtk",1000)
      end do
          
      if (nrank == 0) write(OUTPUT_UNIT,*) "[PF VTK Write] dt, t = ", dt, load(load_n)%time(load_incr)

      call mpi_barrier(mpi_comm_world,ierror)
      if (nrank == 0) write(Flog,"(A,E15.8)")"Temps ecriture VTK PF : ",MPI_Wtime() - t1

      !t9_it = t9_it + MPI_Wtime() - t1 TODO
    end if
  end if


  !================================================================================
  !                                                  EVOLUTION DES CHAMPS DE PHASES

  ! get df/dEta(i=1->Neta) from UMAT
  !---------------------------------------
  call get_df_dEta(df_dEta,Neta)

  ! get alpha coefficient from UMAT
  !---------------------------------------
  call get_interface_energy_coef(alpha,Neta)

  ! get kinetic coefficient from UMAT
  !---------------------------------------
  call get_kinetic_coef(Kinetic,Neta)

  ! CALCUL (Eta0 - L.dt.df_dEta) -> Eta
  !---------------------------------------
  do i=1,Neta
    Eta(:,i) = Eta0(:,i) - Kinetic(i)*dt*df_dEta(:,i)
  end do

  ! TRANSFORMEE DE FOURIER -> EtaF
  !---------------------------------------
  do i=1,Neta
    call field_fft(Eta,EtaF,1,1) ! a convertir si + de 1 PF
  end do
  EtaF=EtaF/ real(grid%ntot,mytype)

  ! schema semi implicte, force motrice interface isotrope en Fourier = EtaF(:,:,:,1) * FREQ2 * alpha
  ! EtaF = EtaF  / (1+Freq^2 alpha.dt) ! ! a convertir si + de 1 PF
  !---------------------------------------------------------------
  do i=1,Neta
    EtaF(:,:,:,i) = EtaF(:,:,:,i)/(1._mytype + FREQ2(:,:,:) * Kinetic(i) * alpha(i) * dt)
  end do

  ! TRANSFORMEE DE FOURIER INVERSE -> Eta
  !---------------------------------------
  do i=1,Neta
    call field_ifft(Eta(:,i),EtaF(:,:,:,i),1,1)
  end do

  ! ACTUALISATION Eta -> Eta0
  !---------------------------------------
  eta0 = eta

  !================================================================================
  !                                  MISE A JOUR DES VARIABLES INTERNES DANS L UMAT
  call updateIntVar_JB(eta,Neta)

end subroutine after_unpas_user




end module standard_user_mod



