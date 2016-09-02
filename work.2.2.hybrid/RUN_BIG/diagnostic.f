!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecriture d'une matrice ligne par ligne (fichier "UNFORMATTED")   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module bib_writemat
      interface
!
      subroutine write_mat(iu,mat)
!
      integer, intent(in) :: iu
      double precision, dimension(:,:), intent(in) :: mat
!
!
      end subroutine write_mat
!
      end interface
      end module bib_writemat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Lecture d'une matrice ligne par ligne (fichier "UNFORMATTED")    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module bib_readmat
      interface
!
      subroutine read_mat(iu,mat)
!
      integer, intent(in) :: iu
      double precision, dimension(:,:), intent(out) :: mat
      end subroutine read_mat
!
      end interface
      end module bib_readmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Numero logique de fichier disponible                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module bib_nofich
!
      contains
!
      function nofich(fich)
!
      implicit none
!
      character (len=*), optional, intent(in) :: fich
      integer :: nofich
!
      character (len=11) :: form
      character (len=26) :: forma
!
      nofich=99
      form="           "
      do while(nofich>0.and.form/="UNDEFINED  ")
         if (nofich/=5.and.nofich/=6) inquire(nofich,form=form)
         if (form/="UNDEFINED  ") nofich=nofich-1
      enddo
      if (form/="UNDEFINED  ") then
         write(6,*) "Aucun numero logique disponible"
         stop
      endif
!
      if (present(fich)) then
         write(forma,'("(""Fichier "",i2,"" : "",a",i3,")")') len_trim(fich)
         write(6,forma) nofich,fich
      endif
!
      return
      end function nofich
!
      end module bib_nofich
!###########################################################
!###########################################################
      program DIAGNOSTIC
!===================================================
!===================================================
	use bib_nofich
	use bib_writemat
	use bib_readmat
      IMPLICIT NONE 
      logical :: trouve,exis440,exis410
      character (len=100):: ligne
      integer :: iu,nr,npr,n,np,ierr,ip,itda,ii,jj
      double precision,dimension(:,:),allocatable::refA,refB,matA,matB
      double precision ::diff,diffb,sumdif
!==========================================================
      !lecture des matrice dans fort.440 et 410
!===================matrice de test=============================      
	inquire(file='fort.440',exist=exis440)
	IF(exis440)then
	  write(6,'("Je lis fort.440 matrice de test")')
	  iu=nofich()
	  open(iu,file='fort.440',form='unformatted',status='old',action='read')
	  read(iu)n,np
!	  read(iu)itda
          allocate(matA(n,n))
	  matA=0.D0
	  CALL READ_MAT(iu,matA)
!	  IF(ITDA.eq.2)then
!V0            allocate(matB(n,n))
!V0            matB=0.D0
!V0	   CALL READ_MAT(iu,matB)
!	  END IF !ITDA.eq.2
	  close(iu)
	 ELSE
	write(6,*)'PAS DE FICHIER  440'
	stop
	END IF !if (exis)
!===================matrice de reference=============================	
	inquire(file='fort.410',exist=exis410)
	IF(exis410)then
	  write(6,'("Je lis fort.410 matrice de reference")')
	  iu=nofich()
	  open(iu,file='fort.410',form='unformatted',status='old',action='read')
	  read(iu)nr,npr
!V0	  read(iu)itda
          allocate(refA(nr,nr))
	  refA=0.D0
	  CALL READ_MAT(iu,refA)
!V0	  IF(ITDA.eq.2)then
!V0            allocate(refB(nr,nr))
!V0            refB=0.D0
!V0	   CALL READ_MAT(iu,refB)
!V0	  END IF !ITDA.eq.2
	  close(iu)
	 ELSE
	write(6,*)'PAS DE FICHIER  410'
	stop
	END IF !if (exis)	
!===============================================================================    	
!V1        IF(nr.ne.n.or.ITDA.ne.2)THEN
        IF(nr.ne.n)THEN
	write(6,*)'les calculs ne sont pas conformes :',n,nr
	stop
	END IF	
!=================================================================================
!   comparaison
!
        sumdif=0.D0
        DIFFB=0.D0
        DO ii=1,N
	DO JJ=ii,N
	DIFF=refA(II,JJ)-matA(II,JJ)
	sumdif=sumdif+DIFF
!V0	DIFFB=refB(II,JJ)-matB(II,JJ)
!V0	if(diff.ge.1.D-14.or.DIFFB.ge.1.D-14)then
        if(diff.ge.1.D-14)then
        write(444,*)"A",II,JJ,DIFF,refA(ii,JJ),matA(II,JJ)
	DIFFB=DIFFB+DIFF
!if(diff.ge.1.D-14.or.DIFFB.ge.1.D-14)then1	write(444,*)"B",II,JJ,DIFFB,refB(II,JJ),matB(II,JJ)
        end if
	END DO
	END DO	
	write(6,*)"  DIFFERENCE = ",DIFFB
	write(6,*)"  somme de toutes les diff =",sumdif
       END ! fin du programme de test
!-----------------------------------------------------------------------------------
      subroutine write_mat(iu,mat)
!
      implicit none
!
      integer, intent(in) :: iu
      double precision, dimension(:,:), intent(in) :: mat
!
      integer :: i, n
!
      n=size(mat,2)
      do i=1,n
         write(iu) mat(:,i)
      enddo
!
      end subroutine write_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Lecture d'une matrice ligne par ligne (fichier "UNFORMATTED")    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
      subroutine read_mat(iu,mat)
!
      implicit none
!
      integer, intent(in) :: iu
      double precision, dimension(:,:), intent(out) :: mat
!
      integer :: i, n, ierr
!
      n=size(mat,2)
      do i=1,n
         read(iu,iostat=ierr) mat(:,i)
      enddo
      if (ierr/=0) then
         write(6,'("Erreur dans read_mat pour iu = ",i3,"   code erreur : ",i4)') iu,ierr
         stop
      endif
!
      end subroutine read_mat
!





