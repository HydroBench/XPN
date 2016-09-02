c-------------------------------------------------------------------------
      module don_mpi
      INTEGER,SAVE::IPRC,NPRC
      end module don_mpi
c-----------------------------------------------------------------------
      module don_matrice
      double precision,dimension(:,:),allocatable,save::RESA,RESB,RES,RESY  !  
      end module don_matrice
c-----------------------------------------------------------------------
      module don_parametres
	INTEGER,save::DIMMAX,NNTOP,MNTOP,MMTOP1,NPTOP,NTOP2,ITOP,DIMOH,Netat,N2TOP,NIOM,N6TOP,N4top
c-------------IOMMAX=2*NNTOP+1 pour base paire IOMMAX=2*NNTOP+1+mod(NNTOP,2) pour les bases impaires	 
	end module don_parametres
c-----------------------------------------------------------------------
      module don_OH
	DOUBLE PRECISION,SAVE::B,BZ,alpha,beta
        DOUBLE PRECISION,DIMENSION(:,:),allocatable,save::FOH1
	DOUBLE PRECISION,DIMENSION(:,:,:),allocatable,save::FOH2
	end module don_OH
c------------------------------------------------------------------------
      module don_constantes
	DOUBLE PRECISION,save::PI,STRUFI,HBARC,ehdeux
	end module don_constantes
c------------------------------------------------------------------------	
	module don_interaction
	DOUBLE PRECISION,SAVE::P1,P2,aw1,ab1,ah1,am1,am2,ah2,ab2,aw2,t3,wso
	end module don_interaction
c-----------------------------------------------------------------------	
	module don_dimensions
      INTEGER,SAVE::nqp,nbqpz,nbqpn,NGR,NGZ,NLD
!      INTEGER,SAVE::nbqpzold,nbqpnold,oldnqp  !!!pour la rerprise du calcul
	end module don_dimensions
c-----------------------------------------------------------------------
        module don_valeur
	INTEGER,DIMENSION(:),allocatable,save::ITABA,ITABB,KQP2
	integer,dimension(:,:),allocatable,save::JZ,IPAR
	DOUBLE PRECISION,dimension(:),allocatable,save::E2QP
	end module don_valeur
c-------------------------------------------------------------------------	
	module don_calcul
	INTEGER,SAVE::ITDA,Kimp,jPI
	end module don_calcul
c----------------------------------------------------------------------
      module don_ILN
	INTEGER,DIMENSION(:),allocatable,save::ILN1,ILN2
	INTEGER,DIMENSION(:,:),allocatable,save::ILN3,ILN4
CCC      COMMON/UNtab/ILN1,ILN2,ILN3,ILN4
	end module don_ILN	
c----------------------------------------------------------------------
      module don_indice
	INTEGER,DIMENSION(:),allocatable,save::IMS,INRS,INZS,IDD,IBL,ISPP,ISPM
      integer,dimension(:,:),allocatable,save::ID1,ID2,ID3,ID1N,ID3N,IDM   
	end module don_indice
c----------------------------------------------------------------------
      module don_centrale
	DOUBLE PRECISION,DIMENSION(:,:),allocatable,save::T111,XJ111,XJ112,T222,XJ221,XJ222  
	end module don_centrale  
c---------------------------------------------------------------------
!      module don_vabcd
!      double precision,dimension(:,:),allocatable,save::VABCD11,VABCD12
!      double precision,dimension(:,:),allocatable,save::VABCDF11,VABCDF12,VABCDF21,VABCDF22,VABCDF31,VABCDF32,VABCDF11c,
!     &VABCDF12c,VABCDF21c,VABCDF22c,VABCDF31c,VABCDF32c,VABCDF11e,VABCDF12e,VABCDF21e,VABCDF22e,VABCDF31e,VABCDF32e
!      double precision,dimension(:,:), allocatable,save::VABCDM1,VABCDH1,VABCDHM1,VABCDH2,VABCDM2,VABCDHM2,
!     &VABCDW1,VABCDB1,VABCDWB1,VABCDW2,VABCDB2,VABCDWB2	
!      end module don_vabcd
c---------------------------------------------------------------------
      module don_vabcd
!      double precision,dimension(:,:),allocatable,save::VABCD11,VABCD12,VABCD21,VABCD22
      double precision,dimension(:,:,:),allocatable,save::TABCD1,TABCD2
      double precision,dimension(:,:,:),allocatable,save::TABCD1HMBW,TABCD1HB,TABCD1MW
      end module don_vabcd      
c-----------------------------------------------------------------------	
      module don_colago
	DOUBLE PRECISION,DIMENSION(:),allocatable,save::X,W,Z,WZ
	end module don_colago
!c------------------------------------------------------------
!      module don_FALOG
!      double precision,dimension(:),allocatable,save::FLF
!	end module don_FALOG
c--------------------------------------------------------------------------	
      module don_vecteur
      double precision,dimension(:,:,:),allocatable,save::XPASU,XPASV,XPASUb,XPASVb
      end module don_vecteur
c------------------------------------------------------------------------------
      module don_phi
      DOUBLE PRECISION,dimension(:,:,:,:),allocatable,save:: PHIPU,PHIMU,PHIPV,PHIMV,APHIPU,APHIMU,APHIPV,APHIMV
      end module don_phi
