      subroutine matrice
        use don_matrice
        use don_parametres
	use don_ILN
        use don_dimensions
	use don_calcul
	use don_valeur
	use don_mpi
	use itf_RESC_pn
	use itf_RESCB_pn
 	use itf_valeur
      IMPLICIT NONE
      include 'mpif.h'
      CHARACTER(len=8)::date
      CHARACTER(len=10)::time
!---------------VARIABLES LOCALES--------------------------------------------
      INTEGER,dimension(22)::valeur12
      INTEGER,dimension(mpi_status_size)::statu
      integer::ip,tag,tagend,ijob,njob,nrec,sender,tagjob,ierr,inout_par
      INTEGER,dimension(:),allocatable::tabij1,tabij2
      INTEGER,dimension(:),allocatable::tabij1G,tabij2G
      double precision::RESCA
      double precision,dimension(2)::RESCAB
      integer::IJ1,IJ2,idenspar
      inout_par=6
      if(iprc==0)THEN
          allocate(tabij1(nqp*(nqp+1)/2))
          allocate(tabij2(nqp*(nqp+1)/2))
	 DO IJ1=1,NBqpZ
	 DO IJ2=IJ1,nbqpz
	    tag=tag+1
	    tabij1(tag)=ij1
	    tabij2(tag)=ij2
	 END DO
         END DO
      END IF
!!      IF(ichoixCE.eq.1.or.idensc.eq.1.or.ichoixSO.eq.1)then
!C***************************************************
!C     CALCUL DU TERME CENTRAL DANS A !!!!!!et B ***
!C***************************************************
      if(iprc==0)then
      WRITE(16,*)'MATRICE A +B : new central+coul+dens+SO'
      write(16,'("calcul avec",i4," processeurs")')nprc
      call date_and_time(date,time)
      write(16,'("A B debut boucle ij1 ij2 date:",a8,3x,"time:",a10)')date,time
      END IF  !!iprc==0
      tagend=0
      tag=0
      njob=nbqpz*(nbqpz+1)/2  !!boucle sur les configurations : 2 QP(PROTON/NEUTRON)- 2QP(PROTON /NEUTRON)
      if (iprc==0)write(16,*)'NJOB proton proton',njob
      call mpi_bcast(njob,1,mpi_integer,0,mpi_comm_world,ierr)
      IF(NJOB/=0)THEN
      if (nprc > 1) then  ! GCdV
      if (iprc==0) then
         ip=0
         do while(ip<nprc-1.and.ip<njob)   !!!and ip<=njob
            ip=ip+1
	    tag=tag+1
            IJ1=tabij1(tag)
      	    valeur12(1:11)=valeur(IJ1,1)
	    IJ2=tabij2(tag)
	    valeur12(12:22)=valeur(IJ2,1)
            call mpi_send(valeur12,22,mpi_integer,ip,tag,mpi_comm_world,ierr)
         enddo
	 nrec=0
	 ijob=ip+1  !!ip=nprc-1 ou njob
	 Do while(ip<nprc-1)
	    ip=ip+1
	    call mpi_send(valeur12,22,mpi_integer,ip,tagend,mpi_comm_world,ierr)
	 end do
	 do while(ijob<=njob)
	    call mpi_recv(RESCAB,2,mpi_double_precision,mpi_any_source,mpi_any_tag,mpi_comm_world,statu,ierr)
	    sender=statu(mpi_source)
	    tagjob=statu(mpi_tag)
	    IJ1=tabij1(tagjob)
	    IJ2=tabij2(tagjob)
	    RESA(IJ1,IJ2)=RESCAB(1)  
	    nrec=nrec+1
	    tag=tag+1
            IJ1=tabij1(tag)
      	    valeur12(1:11)=valeur(IJ1,1)
	    IJ2=tabij2(tag)
	    valeur12(12:22)=valeur(IJ2,1)
            call mpi_send(valeur12,22,mpi_integer,sender,tag,mpi_comm_world,ierr)
	    ijob=ijob+1
	 enddo
	 do while(nrec<njob)
	    call mpi_recv(RESCAB,2,mpi_double_precision,mpi_any_source,mpi_any_tag,mpi_comm_world,statu,ierr)
	    sender=statu(mpi_source)
	    tagjob=statu(mpi_tag)
	    IJ1=tabij1(tagjob)
	    IJ2=tabij2(tagjob)
	    RESA(IJ1,IJ2)=RESCAB(1) 
	    nrec=nrec+1
	    call mpi_send(valeur12,22,mpi_integer,sender,tagend,mpi_comm_world,ierr)
	 enddo
      else
         tagjob=tagend+1  !!valeur arbitraire pour assurer l'entree dans la boucle
         do while(tagjob/=tagend)
            call mpi_recv(valeur12,22,mpi_integer,mpi_any_source,mpi_any_tag,mpi_comm_world,statu,ierr)
            sender=statu(mpi_source)
            tagjob=statu(mpi_tag)
            if (tagjob/=tagend) then
	    RESCAB(1)=0.D0
            RESCAB(2)=0.D0
	    RESCAB(1)=RESC_pn(valeur12)
               call mpi_send(RESCAB,2,mpi_double_precision,0,tagjob,mpi_comm_world,ierr)
            endif
         enddo
      endif
     else
              ! GCdV
              ! dans ce cas on travaille sur un seul processeur: pas la
              ! peine de
              ! faire une gestion fine
              do tagjob = 1, njob
                 IJ1 = tabij1(tagjob)
                 IJ2 = tabij2(tagjob)
                 valeur12(1:11) = valeur(IJ1, 1)
                 valeur12(12:22) = valeur(IJ2, 1)
                 RESCAB(1) = 0.D0
                 RESCAB(2) = 0.D0
                    RESCAB(1) = RESC_pn(valeur12) !central
                 RESA(IJ1, IJ2) = RESA(IJ1, IJ2) + RESCAB(1)
              end do
      endif ! nprc > 1
      END IF !njob/=0
      IF(IPRC==0)THEN
      call date_and_time(date,time)
      write(16,'(" A B fin boucle par ij1 ij2 date:",a8,3x,"time:",a10)')date,time
      END IF
      if(iprc==0)then
      deallocate(tabij1)
      deallocate(tabij2)
      END IF
       END
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
       function RESC_pn(valeur12)
        use don_calcul
	use don_vecteur
	use don_indice
	use don_centrale
	use don_vabcd
!	use don_interaction
	implicit none
	integer,dimension(22),intent(in)::valeur12
	integer::I,J,K,L,NPA,NMA,NPB,NMB,NPC,NMC,NPD,NMD,NNPB,NNMB,NNPC,NNMC,&
     itrIU,itrIV,itrJU,itrJV,itrKU,itrKV,itrLU,itrLV,ipn
	double precision::RESC_pn
	double precision::RES1A,RES1C,RES1B,RES1D,XTE1,XTE2,XTD1,XTD2,RESP
	integer::IA,IB,IC,ID,IAA,IBB,ICC,IDDD,IA1,IA1V,IB1,IB1V,IC1,IC1V,IDD1,IDD1V,I1,I2,I3,I4,K1,K2,K3,K4,&
     IAV,IBV,ICV,IDV
        integer:: NZB,NZD,NZNUM,NZNUX,NZNU
	double precision::S1,S2
       I=valeur12(1)
       J=valeur12(2)
       NPA=valeur12(3)
       NMA=valeur12(4)
       NPC=valeur12(5)
       NMC=valeur12(6)
       itrIU=valeur12(7)
       itrIV=valeur12(8)
       itrJU=valeur12(9)
       itrJV=valeur12(10)
       K=valeur12(12)
       L=valeur12(13)
       NPD=valeur12(14)
       NMD=valeur12(15)
       NPB=valeur12(16)
       NMB=valeur12(17)
       itrKU=valeur12(18)
       itrKV=valeur12(19)
       itrLU=valeur12(20)
       itrLV=valeur12(21)
       NNPB=NPC
       NNMB=NMC
       NNPC=NPB
       NNMC=NMB
       RESP=0.D0
!ccc!cccc BOUCLE I UV(ij) UV(kl)   +UV VU
       RES1A=0.D0
#ifdef OPENMP
   !$OMP PARALLEL DEFAULT( SHARED)                    &
   !$OMP             PRIVATE( IA,IA1,IAA,IA1V,        &          
   !$OMP                      RES1C,IC,IC1,ICC,ICV,I1,K1,IC1V,&
   !$OMP                      RES1B,IB,IB1,IBB,IBV,I4,K4,IB1V,&
   !$OMP                      RES1D,ID,ID1,IDD,IDV,IDD1,IDD1V,IDDD,I2, K2,I3,K3,&
   !$OMP                      XTD1,XTD2,XTE1,XTE2 )
   !$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RES1A)
#endif
       DO IA=1,NMA
          IA1=IA+NPA
          IAA=IA1+Itriu  !!!-(ISI-1)*NETAT/4=translation de NETAT/2 pour etats bar
          RES1C=0.D0
		  DO IC=1,NPC
		     ICC=IC+NMC+Itrjv
             I1=ID1N(IAA,ICC)
             K1=ID3N(IAA,ICC)
             RES1B=0.D0
		     DO IB=1,NPB
			    IBB=IB+NMB+Itrlv
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1=ID+NPD
                   IDDD=IDD1+Itrku
			       I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!-----------------------------B D---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,K1,K2)!*(ah1+am1)
			       XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,K1,K2)!*(ah2+am2)			    
!----------------------------B C-------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(ab1+aw1)
			       XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(ab2+aw2)								
!-----------------------------------------------------------------------------------
			       RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASU(IDD1,K,1)
                END DO  !d-  I
                RES1B=RES1B+RES1D*XPASV(IB,L,2)
             END DO  ! b-    I
             DO IB=1,NMB
                IB1V=IB+NPB
                IBB=IB+Itrlv
           	    RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+Itrku
                   I2=ID1N(IBB,IDDD)
				   K2=ID3N(IBB,IDDD)
!-----------------------------B D---------------------------------------------
	  		       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,K1,K2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,K2)!*ah2
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
			       XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2
!----------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASU(ID,K,1)
            	END DO  !d+ 	I
                RES1B=RES1B+RES1D*XPASV(IB1V,L,2)
             END DO !b+		I
             RES1C=RES1C+RES1B*XPASV(IC,J,2)
          END DO  !c-
          DO IC=1,NMC
             IC1V=IC+NPC
             ICC=IC+Itrjv
             I1=ID1N(IAA,ICC)
	    	 K1=ID3N(IAA,ICC)
             RES1B=0.D0
             DO IB=1,NMB
                IB1V=IB+NPB
                IBB=IB+Itrlv
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1=ID+NPD
                   IDDD=IDD1+Itrku
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD) 
!-----------------------------B D---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2
!---------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASU(IDD1,K,1)
                END DO !d-     I
                RES1B=RES1B+RES1D*XPASV(IB1V,L,2)
             END DO ! b+     I
             RES1C=RES1C+RES1B*XPASV(IC1V,J,2)
          END DO ! c+
          RES1A=RES1A+RES1C*XPASU(IA1,I,1)
       END DO !a-
#ifdef OPENMP
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RES1A)
#endif 
       DO IA=1,NPA
          IAA=IA+Itriu
          RES1C=0.D0
          DO IC=1,NPC
             ICC=IC+NMC+Itrjv
             I1=ID1N(IAA,ICC)
		     K1=ID3N(IAA,ICC)
             RES1B=0.D0
             DO IB=1,NPB
                IBB=IB+NMB+Itrlv
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+Itrku
                   I2=ID1N(IBB,IDDD)
		           K2=ID3N(IBB,IDDD)
!-----------------------------B D---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,K1,K2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,K1,K2)!*am2
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2
!---------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASU(ID,K,1)
                END DO  !d+   I
                RES1B=RES1B+RES1D*XPASV(IB,L,2)
             END DO     !b-   I
             RES1C=RES1C+RES1B*XPASV(IC,J,2)
          END DO   !c-
          DO IC=1,NMC
             IC1V=IC+NPC
             ICC=IC+Itrjv
             I1=ID1N(IAA,ICC)
             K1=ID3N(IAA,ICC)
             RES1B=0.D0
             DO IB=1,NPB
                IBB=IB+NMB+Itrlv
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1=ID+NPD
                   IDDD=IDD1+Itrku
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!-----------------------------B D---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,K1,K2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,K1,K2)!*ah2 
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
			       XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2
!--------------------------------------------------------- --------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASU(IDD1,K,1)
                END DO  !d- I
                RES1B=RES1B+RES1D*XPASV(IB,L,2)
             END DO     !b- I
             DO IB=1,NMB
                IB1V=IB+NPB
                IBB=IB+ITRLV
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!-----------------------------B D---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,K1,K2)!*(ah1+am1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,K1,K2)!*(ah2+am2)
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)  
!------------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASU(ID,K,1)
                END DO  !d+ I
                RES1B=RES1B+RES1D*XPASV(IB1V,L,2)
             END DO     !b+ I
             RES1C=RES1C+RES1B*XPASV(IC1V,J,2)
          END DO   !c+
          RES1A=RES1A+RES1C*XPASU(IA,I,1)
       END DO   !a+,
#ifdef OPENMP
!$OMP END DO
#endif 
!CCCC fin de la boucle UV UV sans les termes  UV VU
!CCC!CCC!CCCCc BOUCLE   VU UV + VU VU en enlevant les termes VU UV
#ifdef OPENMP
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RES1A)
#endif 
       DO IC=1,NMC         !  c-
          IC1=IC+NPC
          ICC=IC1+ITRJU
          RES1C=0.D0
          DO IA=1,NPA          !  a-
		     IAA=IA+NMA+ITRIV
             I1=ID1N(ICC,IAA)
		     K1=ID3N(ICC,IAA)
             RES1B=0.D0
             DO IB=1,NMB
                IB1=IB+NPB
                IBB=IB1+ITRLU
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+NMD+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!-----------------------------D  B ===  B D---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,k1,k2)!*(ah1+am1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,k1,k2)!*(ah2+am2)
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)
!---------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASV(ID,K,1)
                END DO    !d- IV
                RES1B=RES1B+RES1D*XPASU(IB1,L,2)
             END DO       !b- IV
             DO IB=1,NPB
                IBB=IB+ITRLU
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1V=ID+NPD
                   IDDD=ID+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!-----------------------------D  B ===  B D---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,k1,k2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,k2)!*ah2
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
                   XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2  
!------------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASV(IDD1V,K,1)
                END DO    !d+ IV
                RES1B=RES1B+RES1D*XPASU(IB,L,2)
             END DO       !b+ IV
             RES1C=RES1C+RES1B*XPASV(IA,I,1)
          END DO   !c-   donc a-
          DO IA=1,NMA         !a+
             IA1V=IA+NPA
             IAA=IA+ITRIV
             I1=ID1N(ICC,IAA)
		     K1=ID3N(ICC,IAA)
             RES1B=0.D0
             DO IB=1,NMB
                IB1=IB+NPB
                IBB=IB1+ITRLU
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1V=ID+NPD
                   IDDD=ID+ITRKV
                   I2=ID1N(IDDD,IBB)
		   	       K2=ID3N(IDDD,IBB)
!-----------------------------D  B ===  B D---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1  
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2  
!------------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASV(IDD1V,K,1)
                END DO    !d- IV
                RES1B=RES1B+RES1D*XPASU(IB1,L,2)
             END DO       !b+ IV
             RES1C=RES1C+RES1B*XPASV(IA1V,I,1)
          END DO   !c+    donc a+
          RES1A=RES1A+RES1C*XPASU(IC1,J,2)
       END DO   !a-    donc c-
#ifdef OPENMP
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RES1A)
#endif
       DO IC=1,NPC  !c+
          ICC=IC+ITRJU
          RES1C=0.D0
          DO IA=1,NPA          !  a-
             IA1=IA+NMA
             IAA=IA1+ITRIV
             I1=ID1N(ICC,IAA)
		     K1=ID3N(ICC,IAA)
             RES1B=0.D0
             DO IB=1,NPB
                IBB=IB+ITRLU
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+NMD+ITRKV
                   I2=ID1N(IDDD,IBB)
			       k2=ID3N(IDDD,IBB)
!-----------------------------D  B ===  B D---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2  
!--------------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASV(ID,K,1)
                END DO  !d+ IV
                RES1B=RES1B+RES1D*XPASU(IB,L,2)
             END DO     !b- IV
             RES1C=RES1C+RES1B*XPASV(IA,I,1)
          END DO !c-   donc a-
          DO IA=1,NMA
             IA1V=IA+NPA
             IAA=IA+ITRIV
             I1=ID1N(ICC,IAA)
	         K1=ID3N(ICC,IAA)
             RES1B=0.D0
             DO IB=1,NMB
                IB1=IB+NPB
                IBB=IB1+ITRLU
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+NMD+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!-----------------------------D  B ===  B D---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,k1,k2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,k2)!*ah2
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
                   XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2
!------------------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASV(ID,K,1)
                END DO    !d- IV
                RES1B=RES1B+RES1D*XPASU(IB1,L,2)
             END DO       !b- IV
             DO IB=1,NPB
                IBB=IB+ITRLU
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1V=ID+NPD
                   IDDD=ID+ITRKV
                   I2=ID1N(IDDD,IBB)
		           K2=ID3N(IDDD,IBB)
!-----------------------------D  B ===  B D---------------------------------------------
		           XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,k1,k2)!*(am1+ah1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,k1,k2)!*(am2+ah2)
!----------------------------D A--------------------------------------------------
		           XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)
!--------------------------------------------------------------------------------------------------
                   RES1D=RES1D-(XTD1+XTD2+XTE1+XTE2)*XPASV(IDD1V,K,1)
                END DO  !d+ IV
                RES1B=RES1B+RES1D*XPASU(IB,L,2)
             END DO     !b+ IV
             RES1C=RES1C+RES1B*XPASV(IA1V,I,1)
          END DO  !c+   donc a+
          RES1A=RES1A+RES1C*XPASU(IC,J,2)
       END DO   !a+ donc c+
#ifdef OPENMP
!$OMP END DO
#endif 
!CCC!CCC!CCCc fin boucle VU  UV +  VU VU   III + IV   
!CCC!CCC!CCC!CCCCC TERMES DE PAIRING!CCC!CCC!CCC!CCC!CCC!CCC!CCCC pairing si meme isospin
!ccc!ccc BOUCLE V UU(ij) UU(kl)
!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc
#ifdef OPENMP
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RESP)
#endif 
       DO IA=1,NMA          !!a-  b- c- d-
          IA1=IA+NPA
          IAA=IA1+ITRIU
          RES1B=0.D0
          DO IB=1,NNMB
             IB1=IB+NNPB
             IBB=IB1+ITRJU
             RES1C=0.D0
             DO IC=1,NNMC
                IC1=IC+NNPC
                ICC=IC1+ITRLU
                I1=ID1N(IAA,ICC)
		        K1=ID3N(IAA,ICC)
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1=ID+NPD
                   IDDD=IDD1+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!----------------------------- B D---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,k1,k2)!*(ah1+am1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,k1,k2)!*(ah2+am2)
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)
!--------------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASU(IDD1,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASU(IC1,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASU(IB1,J,2)
          END DO
          DO IB=1,NNPB    !a- b+  c- d+ (2)    ici on ne permute pas c et d cest bien la 2ieme bcle en sigma
             IBB=IB+ITRJU
             RES1C=0.D0
             DO IC=1,NNMC
                IC1=IC+NNPC
                ICC=IC1+ITRLU
                I1=ID1N(IAA,ICC)
			    K1=ID3N(IAA,ICC)
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!----------------------------- B D---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------B C--------------------------------------------------
	       	       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2  
!---------------------------------------------------------------------------------------------------
	    		   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASU(ID,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASU(IC1,L,2)
             END DO
             DO IC=1,NNPC    !a- b+ c+ d- (3)
                ICC=IC+ITRLU
                I1=ID1N(IAA,ICC)
			    K1=ID3N(IAA,ICC)
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1=ID+NPD
                   IDDD=IDD1+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!----------------------------- B D---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,k1,k2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,k2)!*ah2
!----------------------------B C--------------------------------------------------
		           XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
                   XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2
!-----------------------------------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASU(IDD1,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASU(IC,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASU(IB,J,2)
          END DO
	      RESP=RESP+RES1B*XPASU(IA1,I,1)   !!!TEST SIGNE
       END DO
#ifdef OPENMP
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RESP)
#endif
       DO IA=1,NPA   !a + b- c- d+ (4)
          IAA=IA+ITRIU
          RES1B=0.D0
          DO IB=1,NNMB
             IB1=IB+NNPB
             IBB=IB1+ITRJU
		     RES1C=0.D0
             DO IC=1,NNMC
                IC1=IC+NNPC
                ICC=IC1+ITRLU
                I1=ID1N(IAA,ICC)
			    K1=ID3N(IAA,ICC)
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!----------------------------- B D---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,k1,k2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,k2)!*ah2  
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
                   XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2
!----------------------------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASU(ID,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASU(IC1,L,2)
             END DO
             DO IC=1,NNPC        !a + b- c+ d- (5)
                ICC=IC+ITRLU
                I1=ID1N(IAA,ICC)
			    K1=ID3N(IAA,ICC)
                RES1D=0.D0
                DO ID=1,NMD
                   IDD1=ID+NPD
                   IDDD=IDD1+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!----------------------------- B D---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2
!------------------------------------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASU(IDD1,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASU(IC,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASU(IB1,J,2)
          END DO
          DO IB=1,NNPB
             IBB=IB+ITRJU
             RES1C=0.D0
             DO IC=1,NNPC
                ICC=IC+ITRLU
                I1=ID1N(IAA,ICC)
                K1=ID3N(IAA,ICC)
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+ITRKU
                   I2=ID1N(IBB,IDDD)
			       K2=ID3N(IBB,IDDD)
!----------------------------- B D---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,k1,k2)!*(ah1+am1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,k1,k2)!*(ah2+am2)
!----------------------------B C--------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)
!--------------------------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASU(ID,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASU(IC,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASU(IB,J,2)
          END DO
          RESP=RESP+RES1B*XPASU(IA,I,1)   	     
       END DO
#ifdef OPENMP
!$OMP END DO
#endif
!CCC!CCC!CCCC fin de la boucle UU UU  
!ccc!ccc BOUCLE VI Vabar Vbbar(ij) Vcbar Vdbar(kl)   <"kl"|v|IJ> =..= <lk||ij>::::<cd|v|ab>
!c                  res1V=0.D0
#ifdef OPENMP
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RESP)
#endif
       DO IA=1,NPA        ! a-b-c-d- boucle (1)
          IAA=IA+NMA+ITRIV
          RES1B=0.d0
          DO IB=1,NNPB
             IBB=IB+NNMB+ITRJV
             RES1C=0.d0
             DO IC=1,NNPC
                ICC=IC+NNMC+ITRLV
                I1=ID1N(ICC,IAA)
			    K1=ID3N(ICC,IAA)
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+NMD+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!----------------------------- D B ---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,k1,k2)!*(ah1+am1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,k1,k2)!*(ah2+am2)
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)
!-----------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASV(ID,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASV(IC,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASV(IB,J,2)
          END DO
          DO IB=1,NNMB      !a- b+ c- d+ == <c- d+ ||a- b+> (2)
             IBV=IB+NNPB
             IBB=IB+ITRJV
             RES1C=0.D0
             DO IC=1,NNPC
                ICC=IC+NNMC+ITRLV
                I1=ID1N(ICC,IAA)
			    K1=ID3N(ICC,IAA)
                RES1D=0.D0
                DO ID=1,NMD
                   IDV=ID+NPD
                   IDDD=ID+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!----------------------------- D B ---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2
!-----------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASV(IDV,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASV(IC,L,2)
             END DO
             DO IC=1,NNMC       !a- b+ c+ d- ===<c+ d- ||a- b+> (4)
                ICV=IC+NNPC
                ICC=IC+ITRLV
		  	    I1=ID1N(ICC,IAA)
			    K1=ID3N(ICC,IAA)
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+NMD+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!----------------------------- D B ---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,k1,k2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,k2)!*ah2
!----------------------------D A--------------------------------------------------
		           XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
                   XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2 
!-------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASV(ID,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASV(ICV,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASV(IBV,J,2)
          END DO  !!!B+
          RESP=RESP+RES1B*XPASV(IA,I,1)   		     
       END DO   !!!!A-
#ifdef OPENMP
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME) REDUCTION(+: RESP)
#endif
       DO IA=1,NMA     !!!!!!!!!!!!!!!!!!!!!!!   !a+ b- c- d+ ==<c- d+ ||a+ b-> (3)
          IAV=IA+NPA
          IAA=IA+ITRIV
          RES1B=0.D0
          DO IB=1,NNPB
             IBB=IB+NNMB+ITRJV         !!!!
             RES1C=0.D0
             DO IC=1,NNPC
                ICC=IC+NNMC+ITRLV
                I1=ID1N(ICC,IAA)
			    K1=ID3N(ICC,IAA)
                RES1D=0.D0
                DO ID=1,NMD
                   IDV=ID+NPD
                   IDDD=ID+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!----------------------------- D B ---------------------------------------------
			       XTD1=TABCD1HB(1,I1,I2)*TABCD2(1,k1,k2)!*ah1
                   XTD2=TABCD1HB(2,I1,I2)*TABCD2(2,k1,k2)!*ah2
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1HB(3,I1,I2)*TABCD2(3,K1,K2)!*ab1
                   XTE2=TABCD1HB(4,I1,I2)*TABCD2(4,K1,K2)!*ab2
!--------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASV(IDV,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASV(IC,L,2)
             END DO
             DO IC=1,NNMC         !a+ b- c+ d- ==< c+d-||a+ b->(5)
                ICV=IC+NNPC
                ICC=IC+ITRLV
                I1=ID1N(ICC,IAA)
			    K1=ID3N(ICC,IAA)
                RES1D=0.D0
                DO ID=1,NPD
                   IDDD=ID+NMD+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)   
!----------------------------- D B ---------------------------------------------
			       XTD1=TABCD1MW(1,I1,I2)*TABCD2(1,k1,k2)!*am1
                   XTD2=TABCD1MW(2,I1,I2)*TABCD2(2,k1,k2)!*am2
!----------------------------D A--------------------------------------------------
		           XTE1=TABCD1MW(3,I1,I2)*TABCD2(3,K1,K2)!*aw1
                   XTE2=TABCD1MW(4,I1,I2)*TABCD2(4,K1,K2)!*aw2 !
!---------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASV(ID,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASV(ICV,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASV(IB,J,2)
          END DO
          DO IB=1,NNMB           !a+ b+  c+ d+  (6)
             IBV=IB+NNPB
             IBB=IB+ITRJV
             RES1C=0.D0
             DO IC=1,NNMC
                ICV=IC+NNPC
                ICC=IC+ITRLV
                I1=ID1N(ICC,IAA)
		   	    K1=ID3N(ICC,IAA)
                RES1D=0.D0
                DO ID=1,NMD
                   IDV=ID+NPD
                   IDDD=ID+ITRKV
                   I2=ID1N(IDDD,IBB)
			       K2=ID3N(IDDD,IBB)
!----------------------------- D B ---------------------------------------------
			       XTD1=TABCD1HMBW(1,I1,I2)*TABCD2(1,k1,k2)!*(ah1+am1)
                   XTD2=TABCD1HMBW(2,I1,I2)*TABCD2(2,k1,k2)!*(ah2+am2)				
!----------------------------D A--------------------------------------------------
			       XTE1=TABCD1HMBW(3,I1,I2)*TABCD2(3,K1,K2)!*(aw1+ab1)
                   XTE2=TABCD1HMBW(4,I1,I2)*TABCD2(4,K1,K2)!*(aw2+ab2)
!-------------------------------------------------------------------------------------
                   RES1D=RES1D+(XTE1+XTE2+XTD1+XTD2)*XPASV(IDV,K,1)
                END DO
                RES1C=RES1C+RES1D*XPASV(ICV,L,2)
             END DO
             RES1B=RES1B+RES1C*XPASV(IBV,J,2)
          END DO
		  RESP=RESP+RES1B*XPASV(IAV,I,1)
       END DO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
!CCC!CCC!CCCC fin de la boucle VV VV
!CCC!CCC!CCCFIN TERME DE PAIRING!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCC!CCCC
       RES1A=RES1A-RESP  
       RESC_pn=RES1A
	END
!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc
!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc
      function valeur(IJ1,ipn)
      use don_valeur
      use don_parametres
      use don_indice
      implicit none
      integer,intent(in)::IJ1,ipn
      integer,dimension(11)::valeur
      integer::I,J,JZI,JZJ,iPRI,IPRJ,ISI,ISJ,II,JJ
      I=ITABA(IJ1)   !!indice du premier etat de "qp" de la configuration1
      J=ITABB(IJ1)   !!!indice du second etat de "qp" de la configuration1
         JZI=JZ(I,1)       !!!U est sur le bloc omega v sur omega bar
         iPRI=IPAR(I,1)
         ISI=SIGN(1,JZI)
         II=JZI-ISI*(iPRI-1)/2
            JZJ=JZ(J,1)
            iPRJ=IPAR(J,1)
            ISJ=SIGN(1,JZJ)
            JJ=JZJ-ISJ*(iPRJ-1)/2
	    valeur(1)=I
	    valeur(2)=J
	    valeur(3)=ISPP(II)!NPA
	    valeur(4)=ISPM(II)!NMA
	    valeur(5)=ISPP(JJ)!NPC
	    valeur(6)=ISPM(JJ)!NMC
	    valeur(7)=IBL(ABS(II))-((ISI-1)*NETAT)/4
	    valeur(8)=IBL(ABS(II))+((ISI+1)*NETAT)/4
	    valeur(9)=IBL(ABS(JJ))-((ISJ-1)*NETAT)/4
	    valeur(10)=IBL(ABS(JJ))+((ISJ+1)*NETAT)/4
	    valeur(11)=ipn
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine prepa_VABCD
        use don_parametres
	use don_indice
	use don_centrale
	use don_vabcd
	use don_mpi
	use don_ILN
	use don_interaction
      IMPLICIT NONE
      integer::I1,I2,NZNUM,NZNUX,NZNU,NZA,NZB,NZC,NZD,I3,I4
      integer::OK
      integer::MA,Na,mb,nb,mc,nc,md,nd,INRM4,INRM3,INRM2,INRM1,Ka,Kb,Kc,Kd,INRNU
      double precision::S1,S2
!C***************************************************
!C     CALCUL DU TERME CENTRAL DANS LA BASE DE OH ***
!C***************************************************
         allocate(TABCD1HMBW(4,NTOP2,NTOP2))
	 allocate(TABCD1HB(4,NTOP2,NTOP2))
	 allocate(TABCD1MW(4,NTOP2,NTOP2))
	 allocate(TABCD1(4,NTOP2,NTOP2))
!	 VABCD11=0.d0
!         VABCD12=0.d0
	 TABCD1=0.D0
      DO NZA=0,NNTOP
         DO NZC=0,NNTOP
            I1=ID1(NZA,NZC)
            DO NZB=0,NNTOP
	       I4=ID1(NZB,NZC)
               DO NZD=0,NNTOP
                  I2=ID1(NZB,NZD)
		  I3=ID1(NZA,NZD)
                  NZNUM=ABS(NZB-NZD)
                  NZNUX=NZB+NZD
                  S1=0.0D0
                  S2=0.0D0
                  DO NZNU=NZNUM,NZNUX  !!direct
                     S1=S1+T111(NZNU,I2)*XJ111(NZNU,I1)
                     S2=S2+T111(NZNU,I2)*XJ112(NZNU,I1)
                  END DO
!                  VABCD11(I2,I1)=S1
!                  VABCD12(I2,I1)=S2
		  TABCD1(1,I1,I2)=S1
		  TABCD1(2,I1,I2)=S2
		  TABCD1HMBW(1,I1,I2)=S1*(ah1+am1)
		  TABCD1HMBW(2,I1,I2)=S2*(ah2+am2)
		  TABCD1HB(1,I1,I2)=S1*ah1
		  TABCD1HB(2,I1,I2)=S2*ah2
		  TABCD1MW(1,I1,I2)=S1*am1
		  TABCD1MW(2,I1,I2)=S2*am2
		  NZNUM=ABS(NZB-NZC)
                  NZNUX=NZB+NZC
                  S1=0.0D0
                  S2=0.0D0
                  DO NZNU=NZNUM,NZNUX  !!echange
                     S1=S1+T111(NZNU,I4)*XJ111(NZNU,I3)
                     S2=S2+T111(NZNU,I4)*XJ112(NZNU,I3)
                  END DO
		  TABCD1(3,I1,I2)=S1
		  TABCD1(4,I1,I2)=S2
		  TABCD1HMBW(3,I1,I2)=S1*(ab1+aw1)
		  TABCD1HMBW(4,I1,I2)=S2*(ab2+aw2)
		  TABCD1HB(3,I1,I2)=S1*ab1
		  TABCD1HB(4,I1,I2)=S2*ab2
		  TABCD1MW(3,I1,I2)=S1*aw1
		  TABCD1MW(4,I1,I2)=S2*aw2
               END DO
            END DO
         END DO
      END DO
       allocate(TABCD2(4,DIMOH,DIMOH),stat=ok)
	if(ok/=0)then
	write(6,*)"probablement pas assez de memoire pour VABCD22"
	stop
	end if   
	TABCD2=0.D0
!T222(INRNU,k2Ia)*XJ221(INRNU,k1Ia)  :::(T222(0:NNTOP,1:DIMOH)) et (XJ222(0:NNTOP,DIMOH))   
      DO MA=-NNTOP,NNTOP
         INRM1=ILN2(IABS(MA))-1
         DO NA=0,INRM1
            KA=ID2(MA,NA)
            DO MB=-NNTOP,NNTOP
               INRM2=ILN2(IABS(MB))-1
               DO NB=0,INRM2
                  KB=ID2(MB,NB)
		  DO MC=-NNTOP,NNTOP
		     INRM3=ILN2(IABS(MC))-1
		     DO NC=0,INRM3
		        KC=ID2(MC,NC)
			I1=ID3(KA,KC)
			I4=ID3(KB,KC)
		        DO MD=-NNTOP,NNTOP
		           INRM4=ILN2(IABS(MD))-1
			   DO ND=0,INRM4
			      KD=ID2(MD,ND)
			      I2=ID3(KB,KD)
			      I3=ID3(KA,KD)
			      S1=0.D0  !!!DIRECT
			      S2=0.D0
			      DO INRNU=0,NNTOP
			      S1=S1+T222(INRNU,I2)*XJ221(INRNU,I1)
			      S2=S2+T222(INRNU,I2)*XJ222(INRNU,I1)
			      END DO
			      TABCD2(1,I1,I2)=S1
			      TABCD2(2,I1,I2)=S2
			      S1=0.D0  !!!ECHANGE
			      S2=0.D0
			      DO INRNU=0,NNTOP
			      S1=S1+T222(INRNU,I4)*XJ221(INRNU,I3)
			      S2=S2+T222(INRNU,I4)*XJ222(INRNU,I3)
			      END DO
			      TABCD2(3,I1,I2)=S1
			      TABCD2(4,I1,I2)=S2
			   END DO
		        END DO
		     END DO
		  END DO
	       END DO
	    END DO   
	 END DO 
      END DO 
	IF(IPRC==0)write(6,*)'fin de prepa central'
      end  !!!! FIN de prepa_VABCD
