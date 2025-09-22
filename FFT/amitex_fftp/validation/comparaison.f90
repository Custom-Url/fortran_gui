!===================================================================
!> Code pour comparer les resultats obtenus avec differentes versions
!! d'AmitexFFTP
!!
!! Deux types de verification :
!!  - Verifier que les resultats sont coherents avec ceux obtenus 
!!    precedemment (valeurs moyennes)
!!  - Verifier que le programme reste au moins aussi efficace
!! \todo Gerer la compatibilite avec les resultats Matlab
!! \todo Ne plus prendre comme reference la valeur moyenne sur toute la cellule mais une
!!       valeur maximale
!!       Ex : Si on a une valeur positive sur la moitie du domaine et une valeur negative
!!            equivalente sur l'autre moitie. Deux valeurs moyennes peuvent etre relativement
!!            eloignees mais les valeurs dans chaque voxel seront proches
!!==============================================================================
program compare

  use decomp_2d

  implicit none

  integer            :: int_ref, int_res, int_sum
  character(len=100) :: fic_ref, fic_res, fic_sum
  character(len=25)  :: test_name
  real(mytype),dimension(6) :: mean_sig_res, mean_def_res, &
                               mean_sig_ref, mean_def_ref
  integer          :: iterations_res, iterations_ref
  logical	   :: diffusion !test si pb diffusion ou mecanique

  call lire_commande(fic_ref, fic_res, fic_sum, test_name, int_ref, int_res, int_sum)  
  call getValues(fic_res, int_res, mean_def_res,mean_sig_res, iterations_res,diffusion)
  call getValues(fic_ref, int_ref, mean_def_ref,mean_sig_ref, iterations_ref,diffusion)
  call compareValues(mean_def_ref,mean_sig_ref, iterations_ref,mean_def_res,mean_sig_res, &
                     iterations_res,fic_sum,int_sum,test_name,diffusion)
end program compare

!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************




!==============================================================================
!
!>   Subroutine de lecture des arguments de la ligne de commande.
!!
!! Lit la ligne de commande qui doit etre de la forme\n
!! ./compare  <reference> <resultat> [-n <nom>] [-s <sortie>] 
!!
!!
!! - <reference>  nom des fichiers contenant les resultats de reference
!!   (suppose que les fichiers .std et .log ont le meme nom)
!! - <resultat>  nom des fichiers contenant les resultats avec la nouvelle version
!!   du code (suppose que les fichiers .std et .log ont le meme nom)
!!
!!optionnel:
!! - -n <nom>      nom du test
!! - -s <sortie>   nom de fichier de bilan des tests
!!
!! 
!! \param[out] fic_ref: (chaine de caracteres) nom des fichiers de reference
!! \param[out] fic_res: (chaine de caracteres) nom des fichiers de resultats
!!                      avec la nouvelle version du code
!! \param[out] int_ref, int_res: (entier) unite logique des fichiers ci-dessus
!!
!====================================================================
subroutine lire_commande(fic_ref, fic_res, fic_sum, test_name, int_ref, int_res, int_sum)


  implicit none

  integer,intent(out)           :: int_ref, int_res, int_sum
  character(len=*),intent(out)  :: fic_ref, fic_res, fic_sum, test_name
  integer                       :: FT,i
  logical                       :: existe
  character(len=100)            :: arg !< arguments du code fortran lance en ligne de commande

  fic_sum = "resume.log"


  !fichier d'entree: nom vides, unites logiques a -1
  ! en fin de subroutine:
  ! unite logique <0 => erreur dans la ligne de commande
  fic_res=""
  fic_ref=""
  int_ref = -1
  int_res = -1

  FT=7

  !! On verifie qu'il y a assez d'arguments
  if(command_argument_count() < 2) then
     write(*,"(A)") "La ligne de commande doit etre de la forme:"
     write(*,*) "./compare  <reference> <resultat> -s <sortie> -n <nom>"
     call abort
  end if


  !! Par defaut le nom du test est celui des
  !! fichiers de resultat
  test_name = fic_res
  !! unite du fichier ou sont resumes les resultats
  int_sum = FT
  FT = FT + 1
  if(command_argument_count()>2) then
     do i=3,command_argument_count(),2
        call getarg(i,arg)     
        if(arg=="-s")then
           !! on recupere le nom du fichier de sortie
           call getarg(i+1,arg)
           read(arg,"(A)") fic_sum
           fic_sum = trim(fic_sum) 
        elseif(arg=="-n")then
           !! on recupere le nom du test
           call getarg(i+1,arg)
           read(arg,*) test_name
           test_name = trim(test_name)
        else
           write(*,"(A)") "ligne de commande invalide."
           write(*,"(2A)") "argument inconnu: ",trim(arg)
           write(*,"(A)") "La ligne de commande doit etre de la forme:"
           write(*,*) "./compare  <reference> <resultat> -s <sortie> -n <nom>"
           call abort
        end if
     end do
  end if

  !! initialisation du fichier de sortie
  call initOutput(fic_sum,int_sum)

  !! On recupere les deux fichiers de resultats + ecriture messages d'erreur
  open(unit=int_sum, file=trim(fic_sum)&
       ,form="formatted",status="old",position="append",action="write")
  
  call get_file_output(1, fic_ref, int_ref, FT, existe)
  if((.NOT.existe)) then 
     write (*,"(A,A,A)")"Erreur (comparaison.f90) : Les fichiers de reference ",&
          trim(fic_ref)," n'existent pas"
     write (int_sum,"(A,A,A)")"ERREUR (comparaison.f90) : Les fichiers de reference ",&
          trim(fic_ref)," n'existent pas"
  end if 
  
  call get_file_output(2, fic_res, int_res, FT, existe)
  if((.NOT.existe)) then 
     write (*,"(A,A,A)")"Erreur (comparaison.f90) :Les fichiers de resultats ",&
          trim(fic_res)," n'existent pas"
     write (int_sum,"(A,A,A)")"ERREUR (comparaison.f90) :Les fichiers de resultats ",&
          trim(fic_res)," n'existent pas"
  end if
  
  close(int_sum)

  ! on verifie que les fichiers sont presents
  if(int_ref==-1) call abort
  if(int_res==-1) call abort
  if(int_sum==-1) call abort

end subroutine lire_commande
!-------------------------------------------------------------------


!===================================================================
!
!       SUBROUTINE GET_FILE_OUTPUT
!>  Lit les fichiers .std et .log associes au nom file_name.\n
!!  Retourne le nom du fichier
!!  et une unite de fichier s'il existe.
!!
!! \param[in]       ind_arg     (entier) numero de l'argument
!! \param[in]       cur_file_u  (entier) derniere unite de fichier utilisee
!! \param[out]      file_exist  (bool) existence du fichier
!! \param[out]      file_unit   (entier) unite du fichier .std (le fichier
!!                              .log est associe a l'unite de fichier suivante)
!! \param[out]      cur_file_u  (entier) unite de fichier disponible (file_unit+2)
!! \param[out]      file_name   (chaine de caracteres) nom du fichier en argument
!
!===================================================================
subroutine get_file_output( ind_arg, file_name, file_unit, cur_file_u, file_exist )

  implicit none

  integer, intent(in)           :: ind_arg
  integer, intent(inout)        :: cur_file_u
  logical, intent(out)          :: file_exist
  integer, intent(out)          :: file_unit
  character(len=*),intent(out)  :: file_name
  character(len=100)            :: arg !< arguments du code fortran lance en ligne de commande

  call getarg(ind_arg,arg)
  read(arg,"(A)") file_name
  INQUIRE(FILE=trim(file_name)//".std", EXIST=file_exist)
  if(file_exist) then
     INQUIRE(FILE=trim(file_name)//".log", EXIST=file_exist)
     if(file_exist) then
        file_unit=cur_file_u
        cur_file_u=cur_file_u +2
     else
        file_unit=-1
     end if
  else
     file_unit=-1
  end if
end subroutine get_file_output
!-------------------------------------------------------------------

!==================================================================================
!                           SUBROUTINE initOutput
!> Initialise le fichier de sortie regroupant les resultats.
!!
!! Si le fichier est vide ou ne respecte pas le format voulu, on efface le fichier
!! et on ecrit la premiere ligne donnant les noms des colonnes.\n
!! Sinon, on ajoute le resultat a la suite des precedents
!!
!! \param[in]  file_sum: (chaine de caracteres) fichier ou seront ecrits les resultats
!! \param[in]  int_sum: unite du fichier
!!
!==================================================================================
subroutine initOutput(file_sum,int_sum)

  implicit none

  character(len=*), intent(in)  :: file_sum
  integer, intent(in)           :: int_sum
  character(len=150), parameter :: str_begin = "Nom test                      Erreur moyenne&
       &                                                Nb iterations"
  character(len=150)            :: tmp_str_begin
  logical                       :: file_exist

  INQUIRE(FILE=trim(file_sum), EXIST=file_exist)
  if(file_exist) then 
     open(unit=int_sum, file=trim(file_sum)&
          ,form="formatted", action="read")
     read(int_sum,"(A)") tmp_str_begin
     if(trim(tmp_str_begin) == trim(str_begin)) then 
        close(int_sum)
     else
        close(int_sum)
        open(unit=int_sum, file=trim(file_sum)&
             ,form="formatted", action="write")
        write(int_sum,"(A)") str_begin
        close(int_sum)
     end if
  else
     open(unit=int_sum, file=trim(file_sum)&
          ,form="formatted",status="NEW", action="write")
     write(int_sum,"(A)") str_begin
     close(int_sum)
  end if
     
end subroutine initOutput

!==================================================================================
!                           SUBROUTINE GET VALUES
!> Recupere les valeurs interessantes pour une serie de resultats.
!!
!!  Lit dans les fichiers file_res.std et file_res.log puis
!!  recupere les valeurs qui vont etre comparees avec l'autre version du programme. \n
!!
!!
!! \param[in]  file_res: (chaine de caracteres) racine des fichiers de resultats
!! \param[in]  int_res: unite du premier fichier de resultats (.std)
!! \param[out]  mean_def: deformation moyenne 
!!		       ou gradQ moyen (mean_def(1:3))
!! \param[out]  mean_sig: contrainte moyenne 
!!		       ou Flux moyen (mean_sig(1:3))
!! \param[out]  iterations: nombre d'iterations
!! \param[out]  diffusion: test si on a un pb de diffusion ou de mecanique
!!		
!! REMARQUE : 
!!	ne prend en compte que la diffusion à une variable
!!
!==================================================================================
subroutine getValues(file_res,int_res, mean_def,mean_sig, iterations,diffusion)

  use decomp_2d

  implicit none

  !moyenne et ecart type des contraintes et deformation sur la cellule
  double precision,dimension(6) :: mean_sig_dp, mean_def_dp
  real(mytype),dimension(6),intent(out) :: mean_sig, mean_def
  character(len=*), intent(in)  :: file_res
  integer, intent(in)           :: int_res
  integer, intent(out)          :: iterations
  logical, intent(out)		:: DIFFUSION
  character(len=5),parameter    :: FMT_real="E15.8"
  !! Chaines de caracteres temporaires pour recuperer
  !! les valeurs importantes
  character(len=1)             :: tmp_1
  character(len=7)             :: tmp_7
  character(len=24)             :: tmp_str_ts
  character(len=100000)         :: tmp_str_long
  logical                       :: HPP

  integer                       :: nb_timestep, int_res_log, ios,compteur

  int_res_log = int_res +1
  open(unit=int_res, file=trim(file_res)//".std"&
       ,form="formatted", status="old", action="read")
  read(int_res,fmt="(/,/,A)") tmp_7
  if(tmp_7 == "#  8-13") then
     HPP = .TRUE.
     DIFFUSION=.FALSE.
  elseif(tmp_7 == "#  8-16") then
     HPP = .FALSE.
     DIFFUSION=.FALSE.
  elseif(tmp_7 == "#   2-4") then
     HPP = .FALSE.
     DIFFUSION=.TRUE.
  else
     print*, "Erreur inattendue "&
          //"(detection HPP, DIFFUSION, ou pas)."
     call abort
  end if
    
  tmp_1="#"
  compteur = 0
  !! On saute les lignes de commentaires
  do while (tmp_1 == "#")
     read(int_res,fmt="(A)") tmp_1
  end do
  !! On veut recuperer les valeurs moyennes de deformation et de contrainte
  !! On stocke donc la chaine de caracteres de la derniere ligne
  ios = 0
  nb_timestep =0
  do while (ios==0)
     nb_timestep = nb_timestep + 1
     read(int_res,fmt="(A)",iostat=ios) tmp_str_long
  end do
  
  !! On saute la colonne des temps (codee sur 15 caracteres) + 1 espace
  if(.not. DIFFUSION) then  
  read(tmp_str_long,fmt="(16X,6("//FMT_real//",1X))",iostat=ios) mean_sig_dp
  if(ios /= 0) then
     print*, "Erreur inattendue "&
          //"(recuperation de la derniere contrainte moyenne)."
     call abort
  end if
  end if

  if(DIFFUSION) then  
  read(tmp_str_long,fmt="(16X,6("//FMT_real//",1X))",iostat=ios) mean_sig_dp
  if(ios /= 0) then
     print*, "Erreur inattendue "&
          //"(recuperation de la derniere contrainte moyenne)."
     call abort
  end if
  end if
  
  !! Deformation
  if(.not. DIFFUSION) then    
  if(HPP) then
     !! Il faut sauter 7 colonnes (16*7 = 112)
     read(tmp_str_long,fmt="(112X,6("//FMT_real//",1X))",iostat=ios) mean_def_dp
  endif
  if(.not. HPP) then
     !! Il faut sauter 16 colonnes (16*16 = 256)
     read(tmp_str_long,fmt="(256X,6("//FMT_real//",1X))",iostat=ios) mean_def_dp
  end if
  end if
  if(DIFFUSION) then
     !! Il faut sauter 16 colonnes (16*4 = 64)
     read(tmp_str_long,fmt="(64X,6("//FMT_real//",1X))",iostat=ios) mean_def_dp
  end if

  if(ios /= 0) then
     print*, "Erreur inattendue "&
          //"(recuperation de la derniere deformation moyenne)."
     call abort
  end if
  close(int_res)
  !! On remet les valeurs dans le bon format
  mean_sig=real(mean_sig_dp,mytype)
  mean_def=real(mean_def_dp,mytype)
  
  !! On va recuperer le nombre d'iteration dans le .log
  open(unit=int_res_log, file=trim(file_res)//".log"&
       ,form="formatted", status="old", action="read")
  iterations = 0
  do while ( iterations == 0)
     read(int_res_log,"(A)", IOSTAT=ios) tmp_str_ts
     if(ios /= 0) then
        print*, "Fin de fichier "//trim(file_res)//".log inattendue "&
             //"(recherche du nombre d'iterations total)."
        call abort
     end if
     if(trim(tmp_str_ts) == trim("nombre total d'iteration")) then
        read(int_res_log,"(I12)", IOSTAT=ios) iterations
        if(ios /= 0) then
           print*, "Fin de fichier "//trim(file_res)//".log inattendue "&
                //"(recuperation du nombre d'iterations total)."
           call abort
        end if
     end if
  end do
  close(int_res_log)
end subroutine getValues

!==================================================================================
!                           SUBROUTINE COMPARE VALUES
!> Compare les valeurs trouvees avec la reference.
!!
!!
!! \param[in]  mean_def_ref, mean_def_res  deformation (ou gradQ) moyenne de la
!!                                         reference et du resultat obtenu
!! \param[in]  mean_sig_ref, mean_sig_res  contrainte (ou Flux) moyenne de la
!!                                         reference et du resultat obtenu
!! \param[in]  iterations_ref, iterations_res  combre d'iterations pour le cas
!!                                         de reference et pour le cas obtenu
!! \param[in]  file_sum     nom du fichier ou sont ecrits les resultats
!! \param[in]  int_sum      unite du fichier ou sont ecrits les resultats
!! \param[in]  test name    nom du test dans le fichier de resultats
!!
!==================================================================================
subroutine compareValues(mean_def_ref,mean_sig_ref, &
     iterations_ref,mean_def_res,&
     mean_sig_res, iterations_res,&
     file_sum,int_sum,test_name,diffusion)

  use decomp_2d

  implicit none
  real(mytype),dimension(6),intent(in) :: mean_sig_res, mean_def_res, &
                                          mean_sig_ref, mean_def_ref
  integer, intent(in)                  :: iterations_res, iterations_ref, int_sum
  character(len=*),intent(in)          :: file_sum, test_name
  logical,intent(in)		       :: diffusion
  !> Parametre a partir duquel on considere que l'ecart relatif
  !! des champs de deformation et de contrainte avec la reference
  !! est trop eleve
  real(mytype),parameter               :: tolerance = 1e-8_mytype
  !> indice
  integer                              :: i
  !> variables intermediaires pour stocker les erreurs relatives et les normes
  !! des champs de reference
  real(mytype)                         :: error_sig, error_def, norm
  !> Chaines de caracteres pour savoir si les resultats sont satisfaisants
  character(len=5)                     :: OK_sig, OK_def, OK_sig_def, OK_it
  !> Format d'un reel
  character(len=5),parameter    :: FMT_real="D15.8"

  OK_sig_def ="   OK"
  !! Calcul de la norme L2 de l'erreur entre les 2 champs de contrainte
  !! ainsi que la norme L2 de champs de contrainte de la reference
  error_sig = 0
  norm = 0
  do i=1,3
     error_sig = error_sig + (mean_sig_res(i) - mean_sig_ref(i))**2
     norm =  norm + mean_sig_ref(i)**2
  end do
  !! Prise en compte de la notation (et diffusion)
  if (.not. diffusion) then
  do i=4,6
     error_sig = error_sig + 2._mytype*(mean_sig_res(i) - mean_sig_ref(i))**2
     norm =  norm + 2._mytype*mean_sig_ref(i)**2
  end do
  end if
  norm = sqrt(norm)
  !! error_sig est l'erreur relative : on divise par
  !! la norme de la reference
  error_sig = sqrt(error_sig)/norm
  !! On regarde si l'erreur est trop elevee
  if(error_sig>tolerance) then 
     OK_sig = "ERROR"
     OK_sig_def = "ERROR"
  else
     OK_sig = "   OK"
  end if
  !! Calcul de la norme L2 de l'erreur entre les 2 champs de deformation
  !! ainsi que la norme L2 de champs de contrainte de la reference  
  error_def = 0
  norm = 0
  do i=1,3
     error_def = error_def + (mean_def_res(i) - mean_def_ref(i))**2
     norm =  norm + mean_def_ref(i)**2
  end do
  !! Prise en compte de la notation
  if (.not. diffusion) then  
  do i=4,6
     error_def = error_def + 0.5_mytype*(mean_def_res(i) - mean_def_ref(i))**2
     norm =  norm + 0.5_mytype*mean_def_ref(i)**2
  end do
  end if
  norm = sqrt(norm)
  !! error_def est l'erreur relative : on divise par
  !! la norme de la reference
  error_def = sqrt(error_def)/norm
  !! On regarde si l'erreur est trop elevee
  if(error_def>tolerance) then 
     OK_def = "ERROR"
     OK_sig_def = "ERROR"
  else
     OK_def = "   OK"
  end if

  !! Si le nombre d'iterations est superieur
  !! a celui trouver avec la reference : ERREUR
  if(iterations_res> iterations_ref) then
     OK_it = "ERROR"
  else
     OK_it = "   OK"
  end if
  !! On ecrit tous les resultats dans le fichier file_sum
  open(unit=int_sum, file=trim(file_sum)&
       ,form="formatted",status="old",position="append",action="write")
  write(int_sum, "(2(A,"//FMT_real//"),2(A,I4),A)")test_name//"    "//&
       OK_sig_def//" (err_sig: ",error_sig,", err_def: ",error_def,")    "//&
       OK_it//" (Ref: ",iterations_ref,", res: ",iterations_res,")"
  close(int_sum)
  write(6, "(A,"//FMT_real//",A)")"Sigma "//OK_sig//" (erreur relative : ",error_sig,")"
  write(6, "(A,"//FMT_real//",A)")"Deformation "//OK_def//" (erreur relative : ",error_def,")"
  write(6, "(2(A,I4),A)")"Iterations "//OK_it//" (Reference : ",iterations_ref,", resultat : ",iterations_res,")"
end subroutine compareValues
