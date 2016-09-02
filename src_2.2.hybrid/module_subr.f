      module itf_INDICE
	interface
      SUBROUTINE INDICE
      IMPLICIT NONE
	END SUBROUTINE INDICE
	end interface
	end module itf_indice
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module itf_PL
	interface
	subroutine pl(n,rr,r2,s,w)
      IMPLICIT NONE
	INTEGER::N
	DOUBLE PRECISION::RR,R2,S,W
	end subroutine pl
	end interface
	end module itf_PL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module itf_LAGUERRE
	interface
	subroutine glag(ngr)
	implicit none
      INTEGER::NGR
	end subroutine glag
	end interface
	end module itf_LAGUERRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module itf_HERMITTE
	interface
	SUBROUTINE DGHERL(NGZ)
      IMPLICIT NONE
	INTEGER::NGZ
	end subroutine DGHERL
	end interface
	end module itf_HERMITTE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module itf_DFAC
	interface
      FUNCTION DFAC(N)
      IMPLICIT NONE
	INTEGER::N
	DOUBLE PRECISION::DFAC
	END FUNCTION DFAC
	end interface
	end module itf_DFAC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module itf_T2N
	interface
      FUNCTION T2N(NB,M,MB,MMU,MNU,NMU,NNU)
      IMPLICIT NONE
	INTEGER::NB,M
	INTEGER::MB,MMU,MNU,NMU,NNU
	DOUBLE PRECISION::T2N
	end FUNCTION T2N
	end interface
	end module itf_T2N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module itf_valeur
      interface
      function valeur(IJ1,ipn)
      implicit none
      integer,intent(in)::IJ1,ipn
      integer,dimension(11)::valeur
      end function valeur
      end interface
      end module itf_valeur
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        module itf_RESC_pn
	interface
       function RESC_pn(valeur12)
	use don_vecteur
	implicit none
	integer,dimension(22),intent(in)::valeur12
	double precision::RESC_pn
	end function RESC_pn
	end interface
	end module itf_RESC_pn
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        module itf_RESCB_pn
	interface
       function RESCB_pn(valeur12)
	use don_vecteur
	implicit none
	integer,dimension(22),intent(in)::valeur12
	double precision::RESCB_pn
	end function RESCB_pn
	end interface
	end module itf_RESCB_pn
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	module itf_DPAGL
	interface
c------------------------------------------------------------------------------
	SUBROUTINE DPAGL(X,W,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION,DIMENSION(:):: X,W
	end subroutine DPAGL
c-----------------------------------------------------------------------------
	end interface
	end module itf_DPAGL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        module itf_prepa_VABCD
	interface
      subroutine prepa_VABCD
        use don_parametres
	use don_indice
	use don_centrale
	use don_vabcd
	use don_mpi
      integer::I1,I2,NZNUM,NZNUX,NZNU,NZA,NZB,NZC,NZD
      double precision::S1,S2
      end  subroutine !prepa_VABCD
       end interface
       end module itf_prepa_VABCD

!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Ecriture d'une matrice ligne par ligne (fichier "UNFORMATTED")   !
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
