!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
!===================================================================================================
!
!       MODULE STANDARD_USER_MOD -> DDD_MM2 
!>
!>      Module for the standard algorithm : implement BEFORE and AFTER unpas user defined procedures 
!>
!!
!!
!! REMARQUE : lors de l'appel eu comportement Varint en sortie vaut varint en entree 
!!            => varint=varint0
!!            -> pas de souci avec la mise a jour des VI varint0=varint dans resolution_NL_base
!!            => Important : si on veut mettre à jour les variables internes avant l'appel à 
!!            -> resolution_un_pas l faut ajouter Deigs à VarInt0, au contraire, si on veut
!!            -> mettre à jour les variables internes après resolution_un_pas (cas present) 
!!            -> il faut ajouter DEigs à VarInt 

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
  
  use dd_mod

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
  bidon = load_n + load_incr + ind_tps
      
  
  if (key_cputime==1) then
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    tmdc0=MPI_WTIME()
    tmdc1=MPI_WTIME()
  end if  

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
  integer                :: i,nsize
  character(100)         :: tmp_tps, tmp_i, suffixe
  real(kind=mytype)      :: timesave1(1:6),timesave2(1:6)
  !real(kind=mytype)      :: local_loading(13),defMoy(6),sigMoy(6)
  !real(kind=mytype),dimension(-2:0) :: t_load

    if (key_cputime==1) then
      !tmdc2=time to solve one step amitex (resolution_un_pas)
      tmdc2 = tmdc2 + MPI_WTIME() - tmdc1

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      tmdc1=MPI_WTIME()
    end if  


        !!5****************************************************CALCUL CONTRAINTE ET ENVOI VERS MM

        !Receive segment positions, stress calculation via interpolation 
        !send stress values at seg position to mM 
        
        call calculation_stress_seg(ratio,Sig,ind_tps)
        
        if ( ind_tps==kkdebug  ) then
          call MPI_Barrier(MPI_COMM_WORLD,ierror)
          print *,"after sending stress to mM. My rank is ",nrank
        end if
        
        !call get_current_loading(load_n,load_incr,defMoy, sigMoy, t_load, local_loading)
        !if (nrank==0) write(*,*)"loading", local_loading(7:12)
        !if (nrank==0) write(*,*)"SigMoy", sigMoy(1:6)
        call calculation_stress_avg(Sig)
        
        if ( ind_tps==kkdebug  ) then
          call MPI_Barrier(MPI_COMM_WORLD,ierror)
          print *,"after sending average stress to mM. My rank is ",nrank
        end if

        
        !!5*********************************************FIN CALCUL CONTRAINTE ET ENVOI VERS MM
    if (key_cputime==1) then
      !tmdc3=time for stress interpolation at the Gauss points
      tmdc3 = tmdc3 + MPI_WTIME() - tmdc1

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      tmdc1=MPI_WTIME()
    end if
        !!6****************************************************CALCUL EIGENSTRAIN 
        ! Calcul eigenstrain
        !
        do idt=1,timestep_ratio
          call eigenstrain_calculation(DEigs,ratio,idt,ind_tps)
        end do
        
        if ( ind_tps==kkdebug  ) then
          call MPI_Barrier(MPI_COMM_WORLD,ierror)
          print *,"after getting new eigenstrains*timeste_ratio from mM. My rank is ",nrank
        end if

        !!6****************************************************FIN CALCUL EIGENSTRAIN 
  
        !!7**************************************************** TRANFERT EIGENSTRAIN VAR.INT MM
            
        ! transfer to MattoP%VarInt (after un_pas the eigenstrain transfer is done to VarInt  )
        !!>MattotP(1)%Varint(:,:) variables internes dimensions: (nVarInt,NPosition)
        do i=1,size(MattotP)
            MattotP(i)%VarInt(1:MattotP(i)%nVarInt,:)=MattotP(i)%VarInt(1:MattotP(i)%nVarInt,:)+DEigs(1:MattotP(i)%nVarInt,:)
        end do 

        !7*********************************************FIN  TRANFERT EIGENSTRAIN VAR.INT MM
        
    if (key_cputime==1) then
      !tmdc4= time eigeinstrain calculation from segment area list
      tmdc4 = tmdc4 + MPI_WTIME() - tmdc1

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      tmdc1=MPI_WTIME()
    endif
                  
        !8*********************************************SAVE INTERNAL VARIABLES FOR RESTART MM  (ATTENTION LG : VALABLE SI 1 SEUL MATERIAU)

        if(MOD(ind_tps,int(user_param%p_real(1))) .eq. 0 .AND.  ind_tps .NE. 0) then 
          do i=1,6
            write(tmp_i,"(I1)") i
            suffixe = trim(tmp_i)//".bin"
               if(nrank==0) call write_header_bin(trim(fic_vtk)//trim(suffixe),&
                                                  grid%nx,grid%ny,grid%nz)
               call convert_field_for_write_bin32(TEMPfield32,MattotP(1)%VarInt0(i,:)) ! conversion 32 bits + bigendian
               call write_bin32(1,TEMPfield32,"DEIGS_"//tmp_i,trim(fic_vtk)//trim(suffixe),key_ascii)
          end do
        end if  
        
        !8*********************************************FIN  SAVE INTERNAL VARIABLES FOR RESTART MM        
 
        !9*********************************************SAVE VTK FILES PERIODICALLY MM

        if(MOD(ind_tps,int(user_param%p_real(3))) .eq. 0 .AND.  ind_tps .NE. 0) then 
          write(tmp_tps,"(I0)") ind_tps
          do i=1,6
            write(tmp_i,"(I1)") i
            suffixe = "_sig_"//trim(tmp_i)//"_"//trim(tmp_tps)//".vtk"
            call print_field_vtk(Sig(:,i),trim(fic_vtk)//trim(suffixe),"Sig_"//trim(tmp_i))
          end do
        end if  

        !9*********************************************SAVE VTK FILES PERIODICALLY MM
        

        !!10***************************************TELL MM TO STAY OR TO LEAVE THE MAIN LOOP
        ! Send 0 to tell the DD code to continue the iterations
        ! Send 1 to tell the DD code to quit the time loop 

        if (load_n < size(load)) then  
          call tell_DD_whattodo(0)
        else 
          if (load_incr < load(load_n)%NIncr) then
            call tell_DD_whattodo(0)
          else
          ! last iteration of the double loop on the loading, amitex tell the DD code to exit
          ! from the time loop
            call tell_DD_whattodo(1)

            if ( ind_tps==kkdebug  ) then
               call MPI_Barrier(MPI_COMM_WORLD,ierror)
               print *,"time step=",ind_tps,"load_n=",load_n," load_incr=",load_incr
            end if
          end if  
        end if 
        !    
        !!10***********************************FIN TELL MM TO STAY OR TO LEAVE THE MAIN LOOP
        
    if (key_cputime==1) then
      !tmdc5= save internal variable for restart+ save vtk files periodically + tell wtd 
      tmdc5 = tmdc5 + MPI_WTIME() - tmdc1
      tmdc6 = tmdc6 + MPI_WTIME() - tmdc0
      !call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      

    

    
      !!9********************************************************start WRITE TIME SPENT BY AMITEX AND DD
      tcomm=tcomm1+tcomm2
        !if (nrank ==0) then     
      if ( ind_tps==0 .and. nrank==0) write(565,*) &
          "#(1)nstep, Time CPU: (2)tot (3)stress_int (4)eigen_definition (5)saving (6)communication (7)resolution " 

      if ( mod((ind_tps*timestep_ratio),100)==0 .and. nrank==0 ) write(565,"(I4,6E15.8)")ind_tps,tmdc6, &
             tmdc3-tcomm2,tmdc4-tcomm1,tmdc5,tcomm,tmdc2
      if ( (ind_tps*timestep_ratio)==100 ) timesave1=(/tmdc6,tmdc3-tcomm2,tmdc4-tcomm1,tmdc5,tcomm,tmdc2/)            
      if ( (ind_tps*timestep_ratio)==200 ) then 
            timesave2=(/tmdc6,tmdc3-tcomm2,tmdc4-tcomm1,tmdc5,tcomm,tmdc2/) 
            call MPI_Comm_Size(MPI_Comm_World,nsize,ierror)
            print *, "nproc DD=",nprocDD,"proc amitex=",nsize
            !compute max, min, and average timing statistics
            !call MPI_Reduce(tmdc6, maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD,ierror);
            !call MPI_Reduce(tmdc6, mintime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD,ierror);
            !call MPI_Reduce(tmdc6, avgtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierror);
            !avgtime =avgtime/dble(nsize);
            if (nrank==0) then
              !write(566,"(I4,I4,6E15.8,E15.8,E15.8,E15.8)") nsize,int(user_param%p_real(6)),(timesave2-timesave1)/100.d0, &
              !  mintime, maxtime, avgtime
              write(566,"(I4,I4,6E15.8)") nsize,int(user_param%p_real(6)),(timesave2-timesave1)/100.d0
            end if
      end if
        !!9**********************************************************end WRITE TIME USED BY AMITEX AND DD
    endif

end subroutine after_unpas_user




end module standard_user_mod



