      program QRPA  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use don_matrice
	use don_parametres
	use don_dimensions
	use don_calcul
	use don_OH
	use don_interaction
	use don_colago
	use don_ILN
	use don_indice
	use don_centrale
	use don_constantes
	use don_vabcd
	use don_vecteur
	use don_valeur
	use don_phi
	use don_mpi
	use itf_INDICE
	use itf_LAGUERRE
	use itf_HERMITTE
	use itf_DFAC
	use itf_T2N
	use itf_prepa_VABCD
	use bib_nofich
	use bib_writemat
	use bib_readmat
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'mpif.h'
      CHARACTER(len=8)::date,fich,fich0
      CHARACTER(len=10)::time
      character (len=12) :: forma
      CHARACTER(len=50)::rephome
      DOUBLE PRECISION,DIMENSION(:,:),allocatable::XHERM
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable::XLAGUE
c------------------------------------------------------------------------------------
      DOUBLE PRECISION,dimension(:,:),allocatable:: XHERMTR,FOH1TR
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable::XLAGUETR,FOH2TR,RHOTR,aqabTR,aqabTRu,aqabTRv   !pour les densites de transition RHOTR p n et tot
c--------------------------------------------------------------------------------------
      integer::ok,sauv0
c-----------------------------------------------------------            
      INTEGER::ierr
c------------------------------------------------------------------------------------------------------
      DOUBLE PRECISION:: t1,t2
!
!  Initialisations MPI
!
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,nprc,ierr)
      call mpi_comm_rank(mpi_comm_world,iprc,ierr)
      t1 = MPI_WTIME()
c--------------------------GRANDEURS PHYSIQUES----------------------------------------------------------
      PI=DACOS(-1.D0)
c------------------------------------------------------------------------------------
      itda =1
      sauv0=440
      read(24)NNTOP,alpha,beta,DIMMAX,NETAT  !!!HNN->NNTOP
      if(IPRC==0)write(16,'(a,i5,a,2i8)')'base', nntop
      if(IPRC==0)write(6,*)'oscillator parameters alpha et beta ',alpha,beta,'DIMMAX NETAT',DIMMAX,NETAT
	read(24)NGR,NGZ,NLD
      if(IPRC==0)write(6,*)'NGR NGZ NLD',NGR,NGZ,NLD
      open (9, file='exmatin.dat', FORM='formatted')   !! fichier d'input pour mode batch
      if(IPRC==0)write(6,*)'Les entrees ecran seront dans le fichier exmatin pour batch'
      b=1/dsqrt(beta)
      bz=1/dsqrt(alpha)
C-------------------PARAMETRES-----------------------------------------------
       MNTOP=NNTOP/2
	MMTOP1=NNTOP-1
	NPTOP=NNTOP+1+mod(NNTOP,2)  !!+pour base impaire
	NTOP2=NPTOP*NPTOP
	ITOP=NPTOP*(MNTOP+1)+1   !!+pour base impaire
        N2TOP=2*NNTOP
	NIOM=NNTOP*2+1 !! NIOM=nombre de blocs
	N6TOP=6*NNTOP
	N4top=2*n2top
	DIMOH=NTOP2*dimmax+1
c------------------------------------------------------------------------------------
      CALL INDICE  !!ne calcule plus NETAT et DIMMAX ils sont lus
c------------------------------------------------------------------------------------
c      write(6,'(a,i3,a,i8)')'base ',nntop
      if(IPRC==0)write(16,'(a,i3)')'base_m',nntop
            if(IPRC==0)write(6,'(a,i3)')'sortie de INDICE base',nntop
      read(21)nqp,nqp,Kimp,jPI  !!double lecture de nqp pour compatibilite de fort.21 avec diago de Marc (premiere valeur nulle)
c-------------------------------------------------------------------------
      if(IPRC==0)write(6,*)'nombre de configurations protons/neutron',nqp
	nbqpz=nqp  !pour calculer rapidement la matrice pn sans changer tous les index
	nbqpn=0
c      if(IPRC==0)write(6,*)'nqp',nqp
      if(IPRC==0)write(16,*)nqp,'configurations a 2 qp proton/neutron'
      if(IPRC==0)write(16,'(i6,a,i4,a,i3)')nbqpz,'config protons/neutron  pour un Kspin',Kimp,' parite',jPI
c-----------------------------------------------------------------------------------
	allocate(ITABA(nqp),stat=ok)
	if(ok/=0)then
	write(6,*)"probablement pas assez de memoire pour ITABA"
	stop
	end if
	allocate(ITABB(nqp),stat=ok)
	if(ok/=0)then
	write(6,*)"probablement pas assez de memoire pour ITABB"
	stop
	end if
	read(21,iostat=ierrA)ITABA
	read(21,iostat=ierrB)ITABB
	if(ierrA.ne.0.or.ierrB.ne.0)then
	write(6,*)'vieux fichier pour ITABA et ITABB?'
	stop
	end if
c-----------------------------------------------------------------------
      allocate(KQP2(nqp),stat=ok)
	if(ok/=0)then
	write(6,*)"probablement pas assez de memoire pour KQP2"
	stop
	end if
c---------------------------------------------------------------------
       read(21)KQP2 !pour la diago de Marc il faut le tableau KQP2 dans fort.21 
c       mais ne sert a rien ici  (avril 2012)
c-----------------------------------------------------------------------
!	allocate(E2QP(nqp),stat=ok)
!	if(ok/=0)then
!	write(6,*)"probablement pas assez de memoire pour E2QP"
!	stop
!	end if
c-------------------------------------------------------------------
!      read(21)E2QP
      read(21)
      read(21)cut  
      if(IPRC==0)write(16,'(a,f12.8)')'pour un cut de ',cut
	allocate(JZ(Netat,2))
	allocate(IPAR(Netat,2))
	allocate(XPASU(DIMMAX,Netat,2))
	allocate(XPASV(DIMMAX,Netat,2))
	allocate(XPASUb(DIMMAX,Netat,2))
	allocate(XPASVb(DIMMAX,Netat,2))
!	allocate(XPAS(DIMMAX,Netat,2))
!	XPAS=0.D0  !initialisation
      read(23)XPASU,XPASV,JZ,IPAR
c-----------------------------------------------------------------------
       if(iprc==0)then
      call date_and_time(date,time)
      write(6,'("date:",a8,3x,"time:",a10)')date,time
      write(16,'("fin lecture exin date:",a8,3x,"time:",a10)')date,time
      end if
c--------------------------------------------------------------------------------
!      allocate(DTOT(NGR,NGZ),stat=ok)
!	if(ok/=0)then
!      write(6,*)"probablement pas assez de memoire pour DTOT"
!	stop
!	end if
!	allocate(DISO(NGR,NGZ,2),stat=ok)
!	if(ok/=0)then
!      write(6,*)"probablement pas assez de memoire pour DISO"
!	stop
!	end if
c-------------------------------------------------------------------
!      read(24)dtot
!      read(24)diso
       read(24)
       read(24)
      read(24,iostat=ierr)P1,P2,aw1,ab1,ah1,am1,am2,ah2,ab2,aw2,t3,wso
!	if(ierr.ne.0)write(6,*)'interaction D1S verouillee car vieux fichier'
      if(IPRC==0)write(6,*)'interaction',P1,P2,aw1,ab1,ah1,am1,am2,ah2,ab2,aw2,t3,wso
c-------------------------------------------------------------------------
      if(iprc==0)then
         allocate(RESA(nqp,nqp),stat=ok)
	 if(ok/=0)then
	    write(6,*)"probablement pas assez de memoire pour RESA"
	    stop
	 endif
	 RESA=0.D0
!V1	    allocate(RESB(nqp,nqp),stat=ok)				    
!V1	    if(ok/=0)then						    
!V1	       write(6,*)"probablement pas assez de memoire pour RESB"      
!V1	       stop							    
!V1	       end if							    
!V1	    RESB=0.D0							    
c==================================================================================
      END IF  !iprc==0
c--------------------------------------------------------------------------
      call glag(ngr)       !!!calcul des poids et des points d'integration
      call dgherl(ngz)
C==============================================================      
      do ipn=1,2   !!  boucle sur isospin
      jphase=0
c      write(16,*)'ipn',ipn
      DO I=1,NETAT/2                          !!!dedouble l'espace
         Ibar=I+NETAT/2                         !!!inf a IHF===T+
         JZ(Ibar,ipn)=-JZ(I,ipn)              !!!sup a IHF===T-
         IPAR(Ibar,ipn)=IPAR(I,ipn)           !!!T sym de renv du sens du temps
      END DO
      DO I=1,NETAT/2
         jphase=0
         Ibar=I+NETAT/2
         iphf=IPAR(I,ipn)
         IIO=JZ(I,ipn)-(IPAR(I,ipn)-1)/2      !!!IIO num du bloc omega pi
         NPI=ISPP(IIO)
         NPIV=ISPP(-IIO)
         NMI=ISPM(IIO)
         NMIV=ISPM(-IIO)
         IMP=IMS(IBL(IIO)+1)!!phase(-)^m des etats +
         IMM=IMP+1  !!phase des etats -
         iphasp=1-2*mod(imp,2)
         iphasm=1-2*mod(imm,2)
         iqsign=iphasp
c--------phase sur la signature des blos 3/2  7/2 11/2
         ijz=jz(i,ipn)/2.D0
         ijzbb=mod(ijz,2) !!!!ijzb=(-)**ijzbb
         ijzb=1-2*ijzbb
         jphase=(1-jz(i,ipn))/2.D0
         jphase=abs(jphase)
         jkphase=mod(jphase,2)
         jphase=1-2*jkphase
         Do j=1,NMI+NPI
         XPASV(J,I,ipn)=jphase*XPASV(J,I,ipn)
         end do
         DO J=1,NMI
           XPASU(J,Ibar,ipn)=jphase*iphasm*XPASU(J+NPI,I,ipn)
           XPASV(J+NPI,I,ipn)=iphasp*XPASV(J+NPI,I,ipn)
           XPASV(J,Ibar,ipn)=jphase*iphasp*XPASV(J+NPI,I,ipn) !!!iphasm
           XPASUb(J+NPI,I,ipn)=jphase*iphasp*XPASU(J+NPI,I,ipn)
           XPASUb(J,Ibar,ipn)=jphase*iphasp*XPASUb(J+NPI,I,ipn)
           XPASVb(J+NPI,I,ipn)=jphase*iphasm*XPASV(J+NPI,I,ipn)
           XPASVb(J,Ibar,ipn)=jphase*iphasm*XPASVb(J+NPI,I,ipn) !!!iphasp
         END DO
         DO J=1,NPI
           XPASU(J+NMI,Ibar,ipn)=jphase*iphasp*XPASU(J,I,ipn)
           XPASV(J,I,ipn)=iphasm*XPASV(J,I,ipn)
           XPASV(J+NMI,Ibar,ipn)=jphase*iphasm*XPASV(J,I,ipn) !! iphasp
           XPASUb(J,I,ipn)=jphase*iphasm*XPASU(J,I,ipn)
           XPASUb(J+NMI,Ibar,ipn)=jphase*iphasm*XPASUb(J,I,ipn)
           XPASVb(J,I,ipn)=jphase*iphasp*XPASV(J,I,ipn)
           XPASVb(J+NMI,Ibar,ipn)=jphase*iphasp*XPASVb(J,I,ipn) !!iphasm
         END DO
      END DO
c-------------------------------------------------------------       
      end do   !! boucle sur isospin
C==================================================================      
      IBAR=0
      I=0
C**************************************
C**************************************
C     CONSTRUCTION DES ETATS DE OH    *
C**************************************
C**************************************
      allocate(ID1(0:NNTOP,0:NNTOP))
	allocate(ID2(-NNTOP:NNTOP,0:NNTOP))
	allocate(ID3(0:ITOP,0:ITOP))
      IN=0
      DO I=0,NNTOP
         DO J=0,NNTOP
            IN=IN+1
            ID1(I,J)=IN
         END DO
      END DO
      if(in.gt.NTOP2)write(6,*)'ATTENTION! id1 max',IN,'NTOP2',NTOP2
      IN=0
      DO I=-NNTOP,NNTOP
         INRM=ILN2(IABS(I))-1
         DO J=0,INRM
            IN=IN+1
            ID2(I,J)=IN
         END DO
      END DO
      IDT2=IN
      IN=0
      ID3(0,0)=0
      DO I=1,IDT2
      ID3(i,0)=0
      ID3(0,i)=0
         DO J=1,IDT2
            IN=IN+1
            ID3(I,J)=IN
         END DO
      END DO
      if(in.gt.DIMOH)write(6,*)'ATTENTION!  DIMOH',IN,DIMOH
      if(in.gt.DIMOH)write(16,*)'ATTENTION!'
!      if (ichoixcoul.eq.1.or.ichoixCE.eq.1)then
C*************nouveau tableau d'indicage**********************
        allocate(ID1N(NETAT,NETAT))
	allocate(ID3N(NETAT,NETAT))
	allocate(IDM(NETAT,NETAT))
	Do IAA=1,NETAT
	   IMA=IMS(IAA)                  !!!IAA=num de l'etat parmis tous ceux de meme ipn
           INA=INRS(IAA)
           IZA=INZS(IAA)
           J1=ID2(IMA,INA)
	   Do ICC=1,NETAT
              IMC=IMS(ICC)
              INC=INRS(ICC)
              IZC=INZS(ICC)
              J3=ID2(IMC,INC)
              ID1N(IAA,ICC)=ID1(IZA,IZC)
              ID3N(IAA,ICC)=ID3(J1,J3)
	      IDM(IAA,ICC)=iabs(IMA-IMC)
	  end do !!ICC
	end do  !do IAA
!      END IF
C**************CALCUL des FONCTIONS d'ONDES***********************************************
        allocate(XHERM(NGZ,0:NPTOP))
	allocate(Xlague(NGR,0:NPTOP,0:MMTOP1))
	allocate(FOH1(NGZ,-1:NPTOP))
	allocate(FOH2(NGR,0:NPTOP,-1:MMTOP1))
c****************************************************************
      DO IH=1,NGZ
         XI=Z(IH)/DSQRT(2.0D0)
         DO INZ=0,NNTOP+1
            IF(INZ.EQ.0) XHERM(IH,INZ)=1.0D0
            IF(INZ.EQ.1) XHERM(IH,INZ)=2*XI
            IF(INZ.NE.0.AND.INZ.NE.1) THEN
               XHERM(IH,INZ)=2*XI*XHERM(IH,INZ-1)-
     &              2*DFLOAT(INZ-1)*XHERM(IH,INZ-2)
            END IF
         END DO
      END DO
      INRM=NNTOP/2+1
      DO IL=1,NGR
         XETA=X(IL)/2
         DO IM=0,NNTOP+1
            DO INR=0,INRM
               IF(INR.EQ.0) XLAGUE(IL,IM,INR)=1.0D0
               IF(INR.EQ.1) XLAGUE(IL,IM,INR)=1+DFLOAT(IABS(IM))-XETA
               IF(INR.NE.0.AND.INR.NE.1) THEN
                  XLAGUE(IL,IM,INR)=((DFLOAT(2*INR+IM)-1.0D0-XETA)
     &               *XLAGUE(IL,IM,INR-1)-DFLOAT(INR-1+IM)
     &               *XLAGUE(IL,IM,INR-2))/DFLOAT(INR)
               END IF
            END DO
         END DO
      END DO

      DO IH=1,NGZ
         FOH1(IH,-1)=0.0D0
      END DO
      DO IL=1,NGR
         DO IM=0,NNTOP+1
            FOH2(IL,IM,-1)=0.0D0
         END DO
      END DO
      DO INZ=0,NNTOP+1
         C2=1.D0/DSQRT(2**INZ*DFAC(INZ))
         DO IH=1,NGZ
            XI=Z(IH)/DSQRT(2.0D0)
            FOH1(IH,INZ)=C2*XHERM(IH,INZ)
         END DO
      END DO
      DO IL=1,NGR
         XETA=X(IL)/2
         DO IM=0,NNTOP+1
            C3=DSQRT(XETA)**IABS(IM)
            DO INR=0,INRM
               C5=DSQRT(DFAC(INR)/DFAC(INR+IABS(IM)))
               FOH2(IL,IM,INR)=C3*C5*XLAGUE(IL,IM,INR)
            END DO
         END DO
      END DO
      deallocate(XHERM)
      deallocate(Xlague)
!====FIN de CALCUL des FONCTIONS d'ONDE
!      if (ichoixcoul.eq.1.or.ichoixCE.eq.1)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CALCUL DES T ET DES J DE LA PARTIE CENTRALE   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      allocate(T111(0:N2TOP,1:NTOP2))
      allocate(T222(0:NNTOP,1:DIMOH))
C     CALCUL DES T commun au termes centrale et coulomb
      DO I=0,n2top
         DO J=1,ntop2
            T111(I,J)=0.D0
         END DO
      END DO
      DO I=0,NNTOP
         DO J=0,NNTOP
            IN=ID1(I,J)
            INZBM=IABS(I-J)
            INZBX=I+J
            DO K=INZBM,INZBX
               I1=K+I-J
               I2=K+J-I
               I3=I+J-K
               J1=MOD(I1,2)
               J2=MOD(I2,2)
               J3=MOD(I3,2)
               IF(J1.EQ.0.AND.J2.EQ.0.AND.J3.EQ.0) THEN
                  XA=DSQRT(DFAC(I)*DFAC(J)*DFAC(K))
                  XB=DFAC(I1/2)*DFAC(I2/2)*DFAC(I3/2)
                  T111(K,IN)=XA/XB
               END IF
            END DO
         END DO
      END DO
      IDT3=ID3(IDT2,IDT2)
      DO I=0,nntop
         DO J=1,DIMOH
            T222(I,J)=0.D0
         END DO
      END DO
      DO I1=-NNTOP,NNTOP
         INRM1=ILN2(IABS(I1))-1
         MMU=I1
         DO J1=0,INRM1
            NMU=J1
            K1=ID2(I1,J1)
            DO I2=-NNTOP,NNTOP
               MB=I2-I1
               INRM2=ILN2(IABS(I2))-1
               MNU=I2
               DO J2=0,INRM2
                  NNU=J2
                  K2=ID2(I2,J2)
                  K3=ID3(K1,K2)
                  INBX=(2*NMU+2*NNU+IABS(MMU)+IABS(MNU)-IABS(MB))/2
                  minbx=nmu-nnu+0.5D0*(iabs(mmu)-iabs(mnu))
                  minbx=iabs(minbx)-iabs(mb)*0.5D0
                  minbx=max0(0,minbx)
                  DO NB=minbx,INBX
                     X1=DFAC(NB)*DFAC(NB+IABS(MB))
                     X2=DFAC(NMU)*DFAC(NMU+IABS(MMU))
                     X3=DFAC(NNU)*DFAC(NNU+IABS(MNU))
                     CFN=DSQRT(X1*X2*X3)
                     XB=2D0*DFLOAT(NB)+DABS(DFLOAT(MB))
                     XMU=2D0*DFLOAT(NMU)+DABS(DFLOAT(MMU))
                     XNU=2D0*DFLOAT(NNU)+DABS(DFLOAT(MNU))
                     X1=.5D0*(XMU-XNU+XB)-DFLOAT(MMU)
                     X2=.5D0*(XB+XNU-XMU)-DFLOAT(MNU)
                     X3=.5D0*(XMU+XNU-XB)
                     IF(X1.NE.0D0) D1=X1/DABS(X1)
                     IF(X2.NE.0D0) D2=X2/DABS(X2)
                     IF(X3.NE.0D0) D3=X3/DABS(X3)
                     XF1=X1-DINT(X1)
                     XF2=X2-DINT(X2)
                     XF3=X3-DINT(X3)
                     MAX=MIN0(IDNINT(X1-D1*XF1)
     &                    ,IDNINT(X2-D2*XF2)
     &                    ,IDNINT(X3-D3*XF3))
                     X1=-.5D0*(XMU-XNU+XB)-DFLOAT(MMU)
                     X2=-.5D0*(XB+XNU-XMU)-DFLOAT(MNU)
                     X3=-X3
                     IF(X1.NE.0D0) D1=X1/DABS(X1)
                     IF(X2.NE.0D0) D2=X2/DABS(X2)
                     IF(X3.NE.0D0) D3=X3/DABS(X3)
                     XF1=X1-DINT(X1)
                     XF2=X2-DINT(X2)
                     XF3=X3-DINT(X3)
                     MIN=MAX0(IDNINT(X1+D1*XF1)
     &                    ,IDNINT(X2+D2*XF2)
     &                    ,IDNINT(X3+D3*XF3))
                     T=0D0
                     DO M=MIN,MAX
                        TI=T2N(NB,M,MB,MMU,MNU,NMU,NNU)
                        T=T+TI
                     END DO
                     I=NB+NMU+NNU
                     T222(NB,K3)=(-1)**I*T*CFN
                  END DO
               END DO
            END DO
         END DO
      END DO
C     CALCUL DES J  pour le central uniquement!!!!
        allocate(XJ111(0:N2TOP,1:NTOP2))
	allocate(XJ112(0:N2TOP,1:NTOP2))
      DO I=0,N2TOP
         DO J=1,NTOP2
            XJ111(I,J)=0.D0
            XJ112(I,J)=0.D0
         END DO
      END DO
      P=P1
C      PI=DACOS(-1.D0)
!      write(6,*)"N2TOP NTOP2",N2TOP,NTOP2
      BZP=DSQRT(BZ**2+.5D0*P**2)
      XK=BZ/BZP
      C=P*PI**.25D0/DSQRT(2.D0)
      CN=C/BZP
      DIVS2=1.D0/DSQRT(2.D0)
      DO I=0,NNTOP
         DO J=0,NNTOP
            IN1=ID1(I,J)
            NZNUMIN=IABS(I-J)
            NZNUMAX=I+J
            DO K=0,2*NNTOP
               NSS=NZNUMIN+K
               N=MOD(NSS,2)
               NZNUMIN1=NZNUMIN
               IF(N.NE.0) NZNUMIN1=NZNUMIN+1
               S=0.D0
               DO L=NZNUMIN1,NZNUMAX,2
                  NS=L+K
                  XA=DFAC(NS)
                  XB=DFAC(K)*DFAC(L)
                  C1=DSQRT(XA/XB)
                  RM10=(-1)**K*DIVS2**NS*C1
                  XC=2.D0**(DFLOAT(NS))*DSQRT(PI)
                  N2=NS/2
                  XX=DSQRT(DFAC(NS)/XC)
                  PHI10=(-1)**N2*XX/DFAC(N2)
                 S=S+T111(L,IN1)*XK**NS*RM10*PHI10
              END DO
               XJ111(K,IN1)=CN*S
            END DO
         END DO
      END DO
      P=P2
      BZP=DSQRT(BZ**2+.5D0*P**2)
      XK=BZ/BZP
      C=P*PI**.25D0/DSQRT(2.D0)
      CN=C/BZP
      DO I=0,NNTOP
         DO J=0,NNTOP
            IN1=ID1(I,J)
            NZNUMIN=IABS(I-J)
            NZNUMAX=I+J
            DO K=0,2*NNTOP
               NSS=NZNUMIN+K
               N=MOD(NSS,2)
               NZNUMIN1=NZNUMIN
               IF(N.NE.0) NZNUMIN1=NZNUMIN+1
               S=0.D0
               DO L=NZNUMIN1,NZNUMAX,2
                  NS=L+K
                  XA=DFAC(NS)
                  XB=DFAC(K)*DFAC(L)
                  C1=DSQRT(XA/XB)
                  RM10=(-1)**K*DIVS2**NS*C1
                  XC=2.D0**(DFLOAT(NS))*DSQRT(PI)
                  N2=NS/2
                  XX=DSQRT(DFAC(NS)/XC)
                  PHI10=(-1)**N2*XX/DFAC(N2)
                  S=S+T111(L,IN1)*XK**NS*RM10*PHI10
               END DO
               XJ112(K,IN1)=CN*S
            END DO
         END DO
      END DO
	allocate(XJ221(0:NNTOP,DIMOH))
	allocate(XJ222(0:NNTOP,DIMOH))
      DO I=0,nntop
         DO J=1,DIMOH
            XJ221(I,J)=0.0D0
            XJ222(I,J)=0.0D0
         END DO
      END DO
      MS=0
      P=P1
      A1=DSQRT(B**2+.5D0*P**2)
      A2=B/A1
      C=.5D0*P**2/A1**2
      DO I1=-NNTOP,NNTOP
         MA=I1
         INRM1=ILN2(IABS(I1))-1
          DO J1=0,INRM1
            NA=J1
            K1=ID2(I1,J1)
            DO I2=-NNTOP,NNTOP
               MC=I2
               MB=I1-I2
             MMU=MB
               MNU=-MMU
               INRM2=ILN2(IABS(I2))-1
               DO J2=0,INRM2
                  NC=J2
                  K2=ID2(I2,J2)
                  K4=ID3(K1,K2)
                  DO J3=0,NNTOP
                     NMU=J3
                     IXMU=2*NMU+IABS(MMU)
                     I2MU=NMU+IABS(MMU)
                     IXA=2*NA+IABS(MA)
                     IXC=2*NC+IABS(MC)
                     XA=DFLOAT(IXA)
                     XC=DFLOAT(IXC)
                     XI1=DABS(XA-XC)-DFLOAT(IABS(MNU))
                     XI2=XA+XC-DFLOAT(IABS(MNU))
                     IMIN=IDNINT(XI1*.5D0+.5D0*DMOD(XI1,2D0))
                     NNUMAX=IDNINT(XI2*.5D0-.5D0*DMOD(XI2,2D0))
                     NNUMIN=MAX0(0,IMIN)
                     S=0D0
                     DO NNU=NNUMIN,NNUMAX
                        NS=NMU+NNU+IABS(MMU)
                        IXNU=2*NNU+IABS(MNU)
                        IXS=2*NS
                        I=NMU+NNU-NS
                        I2NU=NNU+IABS(MNU)
                        XA=DFAC(NS)*DFAC(NS+IABS(MS))
                        XB=DFAC(NMU)*DFAC(I2MU)*DFAC(NNU)*DFAC(I2NU)
                        C2=(-1)**I*DSQRT(XA/XB)
                        RM20=(-1)**(IXNU)*DIVS2**IXS*C2
                        S=S+T222(NNU,K4)*A2**IXS*RM20
                     END DO
                     XJ221(J3,K4)=S*C
                  END DO
               END DO
            END DO
         END DO
      END DO
      MS=0
      P=P2
      A1=DSQRT(B**2+.5D0*P**2)
      A2=B/A1
      C=.5D0*P**2/A1**2
      DO I1=-NNTOP,NNTOP
         MA=I1
         INRM1=ILN2(IABS(I1))-1
         DO J1=0,INRM1
            NA=J1
            K1=ID2(I1,J1)
            DO I2=-NNTOP,NNTOP
               MC=I2
               MB=I1-I2
               MMU=MB
               MNU=-MMU
               INRM2=ILN2(IABS(I2))-1
               DO J2=0,INRM2
                  NC=J2
                  K2=ID2(I2,J2)
                  K4=ID3(K1,K2)
                  DO J3=0,NNTOP
                     NMU=J3
                     IXMU=2*NMU+IABS(MMU)
                     I2MU=NMU+IABS(MMU)
                     IXA=2*NA+IABS(MA)
                     IXC=2*NC+IABS(MC)
                     XA=DFLOAT(IXA)
                     XC=DFLOAT(IXC)
                     XI1=DABS(XA-XC)-DFLOAT(IABS(MNU))
                     XI2=XA+XC-DFLOAT(IABS(MNU))
                     IMIN=IDNINT(XI1*.5D0+.5D0*DMOD(XI1,2D0))
                     NNUMAX=IDNINT(XI2*.5D0-.5D0*DMOD(XI2,2D0))
                     NNUMIN=MAX0(0,IMIN)
                     S=0D0
                     DO NNU=NNUMIN,NNUMAX
                        NS=NMU+NNU+IABS(MMU)
                        IXNU=2*NNU+IABS(MNU)
                        IXS=2*NS
                        I=NMU+NNU-NS
                        I2NU=NNU+IABS(MNU)
                        XA=DFAC(NS)*DFAC(NS+IABS(MS))
                        XB=DFAC(NMU)*DFAC(I2MU)*DFAC(NNU)*DFAC(I2NU)
                        C2=(-1)**I*DSQRT(XA/XB)
                        RM20=(-1)**(IXNU)*DIVS2**IXS*C2
                        S=S+T222(NNU,K4)*A2**IXS*RM20
                     END DO
                     XJ222(J3,K4)=S*C
                  END DO
               END DO
            END DO
         END DO
      END DO
!      call prepa_centrale
      call prepa_VABCD
	call matrice 
!c-----------------------------------------------------------
	deallocate(XJ111)
	deallocate(XJ112)
	deallocate(XJ221)
	deallocate(XJ222)
        deallocate(ID1N)
        deallocate(ID3N)
        deallocate(IDM)
      IF(IPRC==0)THEN	!!!!!FIN  de multi proc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(6,*)iprc,'FIN DE CALCUL MATRICE'
      call date_and_time(date,time)
      write(6,'("date:",a8,3x,"time:",a10)')date,time
      write(16,'("FIN DE CALCUL MATRICE date:",a8,3x,"time:",a10)')date,time
         open(sauv0,file="fort.440",form='unformatted',status='unknown',action='write')
         write(sauv0)nbqpz,nbqpn
!	 write(sauv0)itda 
         call WRITE_MAT(sauv0,RESA)
!V1         call WRITE_MAT(sauv0,RESB)
	 close(sauv0)
	 open(sauv0,file="fort.440",form='unformatted',status='old',action='read')
	 read(sauv0)nbqpz,nbqpn	 
	 call READ_MAT(sauv0,RESA)
	 write(6,*)RESA(1,1), RESA(4,5)
         call date_and_time(date,time)
         write(6,'("date:",a8,3x,"time:",a10)')date,time
         write(16,'("FIN DE SAUVEGARDE DANS fort.440  date:",a8,3x,"time:",a10)')date,time
      do i=5,12
      write(6,'(8e12.4)')(resa(i,j),j=5,12)
      end do
      END IF  !!!if iprc==0  partie sequentielle
c-----------------------------------------------
  101 Format('vecteur ui',9f12.8)
  102 Format('vecteur vi',9f12.8)
  103 Format('vecteur uj',9f12.8)
  104 Format('vecteur vj',9f12.8)
  306 format('TERME DE PAIRING    ', 4f8.4,' pour (ij1 ij2)',2I4)
  319 format((f12.8),'   ', i4,'   ',i2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Finalisation MPI
!
      t2 = MPI_WTIME()
      write(6,340) t2-t1
  340 format('**** Duree totale : ',f12.3,' seconds ****')
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_finalize(ierr)
      stop
      END !!prog principal

      SUBROUTINE INDICE
	use don_parametres
	use don_ILN
	use don_indice
	use don_mpi
      IMPLICIT NONE
c------------VARIABLES LOCALES----------------
      INTEGER::IM,INMMAX,INRP,INR,IN1,INM,INZ,IL1,IL2,ITIN1,ITIN2,IOMMAXODD,
     &IOM,IM1,IM2,IN4,IN5,INRM,INZM,IP,IOMM,I,INRPtest,nper,nzed
      INTEGER::nbr
c------------------------------------------------------------------------------------
	allocate(IDD(NIOM+1))
	allocate(IBL(NIOM+1))
	allocate(ISPP(-NIOM-1:NIOM+1))
	allocate(ISPM(-NIOM-1:NIOM+1))
        allocate(ILN1(0:NNTOP))
	allocate(ILN2(0:NNTOP))
	allocate(ILN3(0:NNTOP,0:MNTOP))
	allocate(ILN4(0:NNTOP,0:MNTOP))
      iln1(0)=0
      IL2=0
      DO IM=0,NNTOP
         INMMAX=NNTOP-IABS(IM)
         IF(MOD(INMMAX,2).EQ.0) INRP=INMMAX/2+1
         IF(MOD(INMMAX,2).NE.0) INRP=(INMMAX-1)/2+1
         IL1=0
         DO INR=0,INRP-1
            IN1=0
            DO INM=0,INMMAX
               INZ=INM-2*INR
               IF(INZ.GE.0) IN1=IN1+1
            END DO
            ILN3(IM,INR)=IN1         !!pour n per  et m donnes nmbre de possibilites pour nz
            IL1=IL1+IN1
            IL2=IL2+IN1
         END DO
         IF(IM.LT.NNTOP) ILN1(IM+1)=IL2  !!indice=nmbre total d'etat (situes avant) pour m donne
         ILN2(IM)=INRP   !!valeur max de n per pour un m donne
      END DO
      DO IM=0,NNTOP
         ITIN1=ILN1(IM)
         ITIN2=0
         DO INR=0,ILN2(IM)-1
            ITIN2=ITIN2+ILN3(IM,INR)
            ILN4(IM,INR)=ITIN1+ITIN2 !!indice reperage dans classement total
         END DO
      END DO
      iommaxodd=NIOM+mod(NNTOP,2)
c==========================================================================
        allocate(IMS(NETAT))
	allocate(INRS(NETAT))
	allocate(INZS(NETAT))
c-----------------------------------------------------------------------------------------
      IN1=0
      DO IOM=1,NIOM,2
         IM1=(IOM-1)/2
         IM2=(IOM+1)/2
         IF(IOM.EQ.NIOM) IM2=NNTOP
         IN4=0
         IN5=0
         DO IM=IM1,IM2
            INRM=ILN2(IM)-1
            DO INR=0,INRM
               INZM=ILN3(IM,INR)
               DO INZ=0,INZM-1
                  IP=(-1)**(INZ+IABS(IM))
                  IF(IP.EQ.1) THEN
                     IN1=IN1+1
                     IMS(IN1)=IM
                     INRS(IN1)=INR
                     INZS(IN1)=INZ
                     IF(IM.EQ.IM1) IN4=IN4+1
                     IF(IM.EQ.IM2.AND.IM2.NE.IM1) IN5=IN5+1
                  END IF
               END DO
            END DO
         END DO
         ISPP(IOM)=IN4
         ISPM(IOM)=IN5
         IDD(IOM)=IN4+IN5
         IN4=0
         IN5=0
         IF(IOM.NE.IOMMAXodd)then    !!base impaires
         DO IM=IM1,IM2
            INRM=ILN2(IM)-1
            DO INR=0,INRM
               INZM=ILN3(IM,INR)
               DO INZ=0,INZM-1
                  IP=(-1)**(INZ+IABS(IM))
                  IF(IP.EQ.-1) THEN
                     IN1=IN1+1
                     IMS(IN1)=IM
                     INRS(IN1)=INR
                     INZS(IN1)=INZ
                     IF(IM.EQ.IM1) IN4=IN4+1
                     IF(IM.EQ.IM2.AND.IM2.NE.IM1) IN5=IN5+1
                  END IF
               END DO
            END DO
         END DO
         ISPP(IOM+1)=IN4
         ISPM(IOM+1)=IN5
         IDD(IOM+1)=IN4+IN5
         end if
      END DO
ccccccccccccccccccccccccccccccccccccccccccccc
      DO IOM=-NIOM,-1,2
         IOMM=-NIOM-IOM-1
         IM1=(IOMM-1)/2
         IM2=(IOMM+1)/2
         IF(IOM.EQ.-1) IM1=-NNTOP
         IN4=0
         IN5=0
         DO IM=IM1,IM2
            INRM=ILN2(IABS(IM))-1
            DO INR=0,INRM
               INZM=ILN3(IABS(IM),INR)
               DO INZ=0,INZM-1
                  IP=(-1)**(INZ+IABS(IM))
                  IF(IP.EQ.1) THEN
                     IN1=IN1+1
                     IMS(IN1)=IM
                     INRS(IN1)=INR
                     INZS(IN1)=INZ
                     IF(IM.EQ.IM1.AND.IM1.NE.IM2) IN4=IN4+1
                     IF(IM.EQ.IM2) IN5=IN5+1
                  END IF
               END DO
            END DO
         END DO
         ISPP(IOMM)=IN4
         ISPM(IOMM)=IN5
         IN4=0
         IN5=0
         IF(IOMM.NE.-iommaxodd)then
         DO IM=IM1,IM2
            INRM=ILN2(IABS(IM))-1
            DO INR=0,INRM
               INZM=ILN3(IABS(IM),INR)
               DO INZ=0,INZM-1
                  IP=(-1)**(INZ+IABS(IM))
                  IF(IP.EQ.-1) THEN
                     IN1=IN1+1
                     IMS(IN1)=IM
                     INRS(IN1)=INR
                     INZS(IN1)=INZ
                     IF(IM.EQ.IM1.AND.IM1.NE.IM2) IN4=IN4+1
                     IF(IM.EQ.IM2) IN5=IN5+1
                  END IF
               END DO
            END DO
         END DO
         ISPP(IOMM-1)=IN4
         ISPM(IOMM-1)=IN5
         END IF
      END DO
cccccccccccccccccccccccccccccccccccccccccccccccccc
      IBL(1)=0
      DO I=2,iommaxodd
         IBL(I)=IDD(I-1)+IBL(I-1)
      END DO
  766 format(a,i3,a,i3)
      END

      FUNCTION DFAC(N)
      IMPLICIT NONE
      INTEGER::N,I
	DOUBLE PRECISION::DFAC
      DFAC=1
      IF(N.NE.0) THEN
         DO I=1,N
            DFAC=DFAC*I
         END DO
      END IF
      END

      FUNCTION T2N(NB,M,MB,MMU,MNU,NMU,NNU)
	use itf_DFAC
      IMPLICIT NONE
	INTEGER::NB,M
	INTEGER::MB,MMU,MNU,NMU,NNU
	DOUBLE PRECISION::T2N
	INTEGER::I1,I2,I3,I4,I5,I6,J1,J2,J3,J4,J5,J6
	DOUBLE PRECISION::IXMU,IXB,IXNU,X1,X2,X3,X4,X5,X6
      IXMU=2*NMU+IABS(MMU)
      IXB=2*NB+IABS(MB)
      IXNU=2*NNU+IABS(MNU)
      I1=IXB+IXMU-IXNU-2*IABS(M+MMU)
      I2=IXB+IXMU-IXNU+2*IABS(M+MMU)
      I3=IXMU+IXNU-IXB-2*IABS(M)
      I4=IXMU+IXNU-IXB+2*IABS(M)
      I5=IXNU+IXB-IXMU-2*IABS(M+MNU)
      I6=IXNU+IXB-IXMU+2*IABS(M+MNU)
      J1=MOD(I1,4)
      J2=MOD(I2,4)
      J3=MOD(I3,4)
      J4=MOD(I4,4)
      J5=MOD(I5,4)
      J6=MOD(I6,4)
	IF(J1.EQ.0.AND.J2.EQ.0.AND.J3.EQ.0.AND.J4.EQ.0.AND.J5.EQ.0.AND.J6.EQ.0) THEN
      I1=I1/4
      I2=I2/4
      I3=I3/4
      I4=I4/4
      I5=I5/4
      I6=I6/4
      X1=DFAC(I1)
      X2=DFAC(I2)
      X3=DFAC(I3)
      X4=DFAC(I4)
      X5=DFAC(I5)
      X6=DFAC(I6)
      T2N=1D0/(X1*X2*X3*X4*X5*X6)
	ELSE
	T2N=0D0
	END IF
      END

cc subroutine pour les poids et les points de Laguerre puis Hermitte
      subroutine glag(ngr)
	use don_colago
	use itf_PL
	implicit none
      INTEGER::NGR
      DOUBLE PRECISION,dimension(:),allocatable:: abscisse,poids
	INTEGER::NMAX,NMAX2,M,N,I,J,NP,NM,N1
	DOUBLE PRECISION::X1,X2,Y1,Y2,Z1,Z2,t1,t2,X3,Y3,Z3,t3
      nmax=ngr
	nmax2=nmax+2
	allocate(abscisse(nmax2))
      abscisse(1)=2.q0
	allocate(poids(nmax2))
	allocate(x(ngr))
	allocate(w(ngr))
      do 8 n=2,nmax
      abscisse(n+1)=abscisse(n-1)*2
      do 1 i=2,n
      j=n-i+2
    1 abscisse(j)=abscisse(j-1)
      abscisse(1)=1.q0/dble(n)
      do 7 m=1,n
      np=0
      nm=0
      n1=2
      x1=abscisse(m)
      x2=abscisse(m+1)
      call pl(n,x1,y1,z1,t1)
      call pl(n,x2,y2,z2,t2)
      if (y1*y2.gt.0.q0) go to 9
    2 if (np+nm.le.4) x3=x1+(x2-x1)*y1/(y1-y2)
      if (np+nm.gt.4) x3=.5q0*(x1+x2)
      if ((n1.lt.9).and.(x3.eq.x1)) x3=.2q0*(4.q0*x1+x2)
      call pl(n,x3,y3,z3,t3)
      n1=n1+1
      if (y3.eq.0.q0) go to 5
      if (y1*y3.lt.0.q0) go to 3
      if (x1.eq.x3) go to 5
      nm=nm+1
      np=0
      x1=x3
      y1=y3
      go to 4
    3 if (x2.eq.x3) go to 5
      np=np+1
      nm=0
      x2=x3
      y2=y3
    4 if (x2-x1.gt.1.q-20) go to 2
    5 abscisse(m)=x3
      poids(m)=z3
      if(n.eq.nmax) then
c      write (6,1001) m,abscisse(m),poids(m)
      x(m)=abscisse(m)
      w(m)=poids(m)
      end if
 1001 format (i5,3d35.25)
    7 continue
    8 continue
      deallocate(abscisse)
	deallocate(poids)
      return
    9 write (16,*) '  Erreur'
      stop
      end
      subroutine pl(n,rr,r2,s,w)
c r2 est le  polynome de Laguerre Ln(rr),w saderivee
c et s le poids d'integration
      IMPLICIT NONE
	INTEGER::N
	DOUBLE PRECISION::RR,R2,S,W
c--------VARIABLES LOCALES-----------------
	INTEGER::I
	DOUBLE PRECISION::R1,R3
      r2=1.q0
      r1=1.q0-rr
      do 1 I=1,n
      r3=r2
      r2=r1
    1 r1=((2*i+1-rr)*r2-i*r3)/(i+1)
      s=rr/((n+1)*r1)**2
      w=(n+1)*(r1-r2)/rr+r2
      return
      end
      SUBROUTINE DGHERL(NGZ)
	use don_colago
      IMPLICIT NONE
	INTEGER::NGZ
      INTEGER::NMAX,I,J,L,IT
	DOUBLE PRECISION,dimension(:),allocatable:: absciss,poids
	DOUBLE PRECISION::C,XT,Q,P,PN,DP,DPN,B,D,DQ
      nmax=ngz+2
	allocate(absciss(nmax))
	allocate(poids(nmax))
	allocate(z(ngz))
	allocate(wz(ngz))
      C=1.7724538509055D0
      absciss(1)=0.D0
      poids(1)=C
      IF(Ngz.EQ.1) RETURN
      XT=1.27D0/(0.1D0+2.4D0*Ngz)
      DO 1 J=2,Ngz
   1  C=C*(J-1)*(J-1.5D0)
      DO 3 I=1,Ngz
      IF(I.GT.2) GO TO 5
      IF(I.EQ.2) XT=XT+11.88D0/(0.55D0+2.5D0*Ngz)
      GO TO 6
   5  L=I-2
      Q=absciss(L)
      XT=XT+(XT-Q*Q)*(1+6.05D0*L+7.728D0*L*L)/((1.625D0+5.653D0*L)*L)
   6  CONTINUE
      DO 2 IT=1,10
      P=1.D0
      PN=XT-0.5D0
      DP=0.D0
      DPN=1.D0
      DO 4 J=2,Ngz
      B=XT-2*J+1.5D0
      D=(J-1)*(J-1.5D0)
      Q=B*PN-D*P
      DQ=B*DPN-D*DP+PN
      P=PN
      PN=Q
      DP=DPN
   4  DPN=DQ
   2  XT=XT-PN/DPN
      absciss(I)=SQRT(XT)
      poids(I)=0.5D0*C/(DPN*P)
      z(i)=absciss(i)
      wz(i)=poids(i)
c      write (6,1001) i,absciss(i),poids(i)
    3 continue
      deallocate(absciss)
	deallocate(poids)
 1001 format (i5,2d35.25)
      RETURN
      END
c-----------------------------------------------------------------------------------
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


