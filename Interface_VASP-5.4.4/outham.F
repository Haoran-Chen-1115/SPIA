#include "symbol.inc"

  SUBROUTINE OUT_HAM(GRID, NONL_S_ORI, PP, WDES, LATT_CUR, T_INFO, INFO, SV, &
                LMDIM, CDIJ, CQIJ, HAMILTONIAN, EFERMI, TOTEN, NSTEP_ori,DYN)
    USE prec
    USE mgrid
    USE nonl
    USE pseudo
    USE wave
    USE lattice
    USE base
    USE poscar
    USE mpimy
    USE fock
    USE mymath
    USE hamil_high

    IMPLICIT NONE
    ! input
    TYPE (grid_3d) :: GRID
    TYPE (nonl_struct) :: NONL_S_ORI
    TYPE (potcar) :: PP(NONL_S_ORI%NTYP)
    TYPE (wavedes) :: WDES
    TYPE (latt) :: LATT_CUR
    TYPE (type_info) :: T_INFO
    TYPE (info_struct) :: INFO
    RGRID :: SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
    INTEGER :: LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q) :: EFERMI
    REAL(q) :: TOTEN

    ! extra check
    LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION
    TYPE (ham_handle) :: HAMILTONIAN

    ! local
    TYPE (nonl_struct) :: NONL_S
    INTEGER :: ISP, NK, NT, LMMAXC, NIS,CH
    INTEGER :: NPLWV, NPLMAX, NPL
    INTEGER, ALLOCATABLE :: G_INDEX(:,:),G_INDEX_LOC(:,:)
    INTEGER :: NP, NNP, ISPINOR
    COMPLEX(q), ALLOCATABLE :: FHAM(:,:),CWORK(:)
    REAL(q),ALLOCATABLE :: WAE_PS(:,:),WAE(:,:),WPS(:,:)

    ! io
    INTEGER :: IUHAM, i, j,CH0,CH1,CH2,CH3, i_in1, i_in2
    CHARACTER(5) :: NKINDEX,NISP,NTINDEX,NB_C,INB_INDEX
    CHARACTER(40) :: FORM

    ! time
    INTEGER :: t1,t2
    
    ! molecular dynamics
    TYPE (dynamics)    DYN
    INTEGER :: NSTEP_ori, NSTEP
    CHARACTER(10) :: NSTEP_C
    CHARACTER(40) :: MK_DIR
    CHARACTER(200) :: MV
    INTEGER :: NKPTS
    
    LOGICAL :: LINDEX    

    LINDEX = INFO%LINDEX

    ! Hybrid functional not supported
    ! meta-GGA not supported
    IF (USEFOCK_CONTRIBUTION()) THEN
        WRITE(*,*) 'Hybrid functional not supported, quit'
        CALL M_exit(); stop
    ENDIF

    IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
        WRITE(*,*) 'meta-GGA not supported, quit'
        CALL M_exit(); stop
    ENDIF

    IF (WDES%LSPIRAL) THEN
        WRITE(*,*) 'spinor not supported, quit'
        CALL M_exit(); stop
    ENDIF

    ! make a copy for non_local structure
    !NONL_S=NONL_S_ORI
    !IF (INFO%LREAL) THEN
        ! if projector function was calculated in real space
        ! turn to reciprocal space now
        CALL NONL_ALLOC(NONL_S, T_INFO, PP, WDES, .FALSE.)
        CALL SPHER(GRID, NONL_S, PP, WDES, LATT_CUR,  1)
        CALL PHASE(WDES,NONL_S,0)
    !ENDIF

    WRITE(NB_C,'(g0)') INFO%NB

    NPLWV=GRID%NGX*GRID%NGY*GRID%NGZ
    NPLMAX=maxval(WDES%NGVECTOR)

    NSTEP=NSTEP_ori+INFO%NCONT
    IF (MOD(NSTEP_ori+INFO%NCONT-1,INFO%NJUMP).EQ.0) THEN

    IF (DYN%IBRION==0 .OR. DYN%IBRION==23) THEN
        WRITE(NSTEP_C,'(g0)') NSTEP
        MK_DIR='mkdir ./BEAD_'//TRIM(NSTEP_C)//' 1>/dev/null 2>&1'
        CALL system( MK_DIR )
    ENDIF

    IF (MOD(0,WDES%COMM_KIN%NCPU) .NE. (WDES%COMM_KIN%NODE_ME-1)) THEN
        GOTO 202
    ENDIF

    CALL system_clock(t1)
    IF (INFO%NB>1 .OR. NSTEP>1) LINDEX=.FALSE.

    IF (INFO%LKPOINTS) THEN
        NKPTS = WDES%NKPTS
    ELSE
        NKPTS = 1
    ENDIF


      IF (MOD(0,WDES%COMM_KINTER%NCPU) .EQ. (WDES%COMM_KINTER%NODE_ME-1)) THEN
          IF (LINDEX) THEN
          IUHAM=540
          OPEN(IUHAM,FILE='./VKPT',&
                  STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8)
          DO NK=1,WDES%NKPTS
            DO i=1,3
              WRITE(IUHAM,REC=3*(NK-1)+i) WDES%VKPT(i,NK)
            ENDDO
          ENDDO
          CLOSE(IUHAM)
          ENDIF

          IUHAM=540
          OPEN(IUHAM,FILE='HEAD_'//TRIM(NB_C),&
                  STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8)
          DO i=1,3
          DO j=1,3
            WRITE(IUHAM,REC=j+3*(i-1)) LATT_CUR%A(j,i)
          ENDDO
          ENDDO
          WRITE(IUHAM,REC=10) EFERMI
          WRITE(IUHAM,REC=11) TOTEN
          CLOSE(IUHAM)
          OPEN(IUHAM,FILE='HEAD_'//TRIM(NB_C),&
                  STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4)
          WRITE(IUHAM,REC=23) WDES%NKPTS
          WRITE(IUHAM,REC=24) WDES%ISPIN
          WRITE(IUHAM,REC=25) NONL_S%NTYP
          DO i=1,NONL_S%NTYP
            WRITE(IUHAM,REC=25+i) NONL_S%NITYP(i)
          ENDDO
          WRITE(IUHAM,REC=26+NONL_S%NTYP) WDES%COMM_INB%NCPU
          CLOSE(IUHAM)
 
          IF (DYN%IBRION==0 .OR. DYN%IBRION==23) THEN
              MV='mv ./HEAD_'//TRIM(NB_C)//' ./BEAD_'//TRIM(NSTEP_C)//' 1>/dev/null 2>&1'
              CALL system(MV)
          ENDIF
      ENDIF
 
      ! index for calculated G-grids in the whole grid
      IF (LINDEX) THEN
        ALLOCATE(G_INDEX(NPLMAX,WDES%NKPTS),G_INDEX_LOC(NPLWV,WDES%NKPTS))
        CALL GEN_G_INDEX(G_INDEX,G_INDEX_LOC,NPLWV,NPLMAX,WDES,GRID)
      ENDIF
 
      IF (MOD(0,WDES%COMM_KINTER%NCPU) .EQ. (WDES%COMM_KINTER%NODE_ME-1)) THEN
          NIS=1
          typ2:DO NT=1,NONL_S%NTYP
            LMMAXC=NONL_S%LMMAX(NT)
  
            IF (LMMAXC.NE.0) THEN
              WRITE(NTINDEX,'(g0)') NT
              ! take a lot of space, write position instead.
              ! phase
              IUHAM=540
              OPEN(IUHAM,FILE='POS_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                      STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                      RECL=24*NONL_S%NITYP(NT))
  
              WRITE(IUHAM,REC=1) NONL_S%POSION(1:3,NIS:NIS+NONL_S%NITYP(NT)-1)
  
              CLOSE(IUHAM)
  
              ! information needed for transformation opertor T
              ! radial grids and simpson coefficient
              IF (LINDEX) THEN
                  ! INDEX for lm 
                  OPEN(IUHAM,FILE='LM_INDEX_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=4*PP(NT)%LMAX)
                  WRITE(IUHAM,Rec=1) PP(NT)%LPS(1:PP(NT)%LMAX)
         
                  CLOSE(IUHAM)

                  CALL SET_SIMP(PP(NT)%R)
                  OPEN(IUHAM,FILE='SIMPI_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8)
                  WRITE(IUHAM,REC=1) PP(NT)%R%NMAX,NT
               
                  DO i=1,PP(NT)%R%NMAX
                    WRITE(IUHAM,REC=i+1) PP(NT)%R%R(i)
                  ENDDO
               
                  DO i=1,PP(NT)%R%NMAX
                    WRITE(IUHAM,REC=i+PP(NT)%R%NMAX+1) PP(NT)%R%SI(i)
                  ENDDO
               
                  CLOSE(IUHAM)
               
                  ALLOCATE(WAE_PS(PP(NT)%R%NMAX,PP(NT)%LMAX))
                  DO CH=1,PP(NT)%LMAX
                    WAE_PS(1:PP(NT)%R%NMAX,CH)=PP(NT)%R%R(1:PP(NT)%R%NMAX)*&
                            (PP(NT)%WAE(1:PP(NT)%R%NMAX,CH)-&
                            PP(NT)%WPS(1:PP(NT)%R%NMAX,CH))
                  ENDDO
               
                  OPEN(IUHAM,FILE='PWAV_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*PP(NT)%LMAX*PP(NT)%R%NMAX)
                  WRITE(IUHAM,REC=1) WAE_PS
                  CLOSE(IUHAM)
                  DEALLOCATE(WAE_PS)
               
                  ALLOCATE(WAE(PP(NT)%R%NMAX,PP(NT)%LMAX))
                  DO CH=1,PP(NT)%LMAX
                    WAE(1:PP(NT)%R%NMAX,CH)=&
                            PP(NT)%WAE(1:PP(NT)%R%NMAX,CH)
                  ENDDO
               
                  OPEN(IUHAM+1,FILE='WAE_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*PP(NT)%LMAX*PP(NT)%R%NMAX)
                  WRITE(IUHAM+1,REC=1) WAE
                  CLOSE(IUHAM+1)
                  DEALLOCATE(WAE)
  
                  ALLOCATE(WPS(PP(NT)%R%NMAX,PP(NT)%LMAX))
                  DO CH=1,PP(NT)%LMAX
                    WPS(1:PP(NT)%R%NMAX,CH)=&
                            PP(NT)%WPS(1:PP(NT)%R%NMAX,CH)
                  ENDDO
               
                  OPEN(IUHAM+2,FILE='WPS_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*PP(NT)%LMAX*PP(NT)%R%NMAX)
                  WRITE(IUHAM+2,REC=1) WPS
                  CLOSE(IUHAM+2)
                  DEALLOCATE(WPS)
               
                  ! Projector for interpolation 
                  OPEN(IUHAM+3,FILE='PSPNL_'//TRIM(NB_C)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8)
                  WRITE(IUHAM+3,REC=1) PP(NT)%PSMAXN
                  DO CH=1,PP(NT)%LDIM
                    DO i=1,NPSNL+1
                      WRITE(IUHAM+3,REC=i+(CH-1)*(NPSNL+1)+1) PP(NT)%PSPNL(i-1,CH)
                    ENDDO
                  ENDDO
                  CLOSE(IUHAM+3)
              ENDIF
  
            ENDIF
            NIS=NIS+NONL_S%NITYP(NT)
          ENDDO typ2
 
          IF (DYN%IBRION==0 .OR. DYN%IBRION==23) THEN
              MV='mv ./POS_'//TRIM(NB_C)//'_* ./SIMPI_'//TRIM(NB_C)//'_* ./PWAV_'//TRIM(NB_C)//'_* ./WAE_'//TRIM(NB_C)//'_* ./WPS_'//TRIM(NB_C)//'_* ./BEAD_'//TRIM(NSTEP_C)//' 1>/dev/null 2>&1'
              CALL system(MV)
          ENDIF
      ENDIF
 
 
      WRITE(INB_INDEX,'(g0)') WDES%COMM_INB%NODE_ME
      spin: DO ISP=1,WDES%NCDIJ
 
        IUHAM=540
        WRITE(NISP,'(g0)') ISP
        ! information needed for local potential and kinetic energy
        IF (MOD(0,WDES%COMM_KINTER%NCPU) .EQ. (WDES%COMM_KINTER%NODE_ME-1)) THEN
          IF (WDES%COMM_INB%NCPU.NE.1) THEN
#ifdef realmode
              OPEN(IUHAM,FILE='SV_'//TRIM(NB_C)//'_'//TRIM(NISP)//'_'//TRIM(INB_INDEX),&
                      STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*GRID%RL%NP)
#else
              OPEN(IUHAM,FILE='SV_'//TRIM(NB_C)//'_'//TRIM(NISP)//'_'//TRIM(INB_INDEX),&
                      STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=16*GRID%RL%NP)
#endif
          ELSE
#ifdef realmode
              OPEN(IUHAM,FILE='SV_'//TRIM(NB_C)//'_'//TRIM(NISP),&
                      STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*GRID%RL%NP)
#else
              OPEN(IUHAM,FILE='SV_'//TRIM(NB_C)//'_'//TRIM(NISP),&
                      STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=16*GRID%RL%NP)
#endif
          ENDIF
          WRITE(IUHAM,REC=1) SV(1:GRID%RL%NP,ISP)
          CLOSE(IUHAM)

          NIS=1;i_in1=1;i_in2=1
          typ1: DO NT=1,NONL_S%NTYP
          LMMAXC=NONL_S%LMMAX(NT)
     
            IF (LMMAXC.NE.0) THEN
              WRITE(NTINDEX,'(g0)') NT
              ! nonlocal pseudopotential strength
              IF (WDES%COMM_INB%NCPU.NE.1) THEN
#ifdef realmode
                  OPEN(IUHAM,FILE='./CDIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*LMMAXC)
#else
                  OPEN(IUHAM,FILE='./CDIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=16*LMMAXC)
#endif
              ELSE
#ifdef realmode
                  OPEN(IUHAM,FILE='./CDIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*LMMAXC)
#else
                  OPEN(IUHAM,FILE='./CDIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=16*LMMAXC)
#endif
              ENDIF
              ! DO i=1,NONL_S%NITYP(NT)
              DO i=NIS,NIS+NONL_S%NITYP(NT)-1
                IF (NI_LOCAL(i,WDES%COMM_INB) /=0 ) THEN
                  DO j=1,LMMAXC
                    WRITE(IUHAM,Rec=(i-NIS)*LMMAXC+j) CDIJ(1:LMMAXC,j,i_in1,ISP)
                  ENDDO
                  i_in1=i_in1+1
                ENDIF
              ENDDO
     
              CLOSE(IUHAM)
     
              ! overlap strength
              IF (WDES%COMM_INB%NCPU.NE.1) THEN
#ifdef realmode
                  OPEN(IUHAM,FILE='./CQIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*LMMAXC)
#else
                  OPEN(IUHAM,FILE='./CQIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=16*LMMAXC)
#endif
              ELSE
#ifdef realmode
                  OPEN(IUHAM,FILE='./CQIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8*LMMAXC)
#else
                  OPEN(IUHAM,FILE='./CQIJ_'//TRIM(NB_C)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=16*LMMAXC)
#endif
              ENDIF
              ! DO i=1,NONL_S%NITYP(NT)
              DO i=NIS,NIS+NONL_S%NITYP(NT)-1
                IF (NI_LOCAL(i,WDES%COMM_INB) /=0 ) THEN
                  DO j=1,LMMAXC
                    WRITE(IUHAM,Rec=(i-NIS)*LMMAXC+j) CQIJ(1:LMMAXC,j,i_in2,ISP)
                  ENDDO
                  i_in2=i_in2+1
                ENDIF
              ENDDO
     
              CLOSE(IUHAM)
              NIS=NIS+NONL_S%NITYP(NT)
            ENDIF
          ENDDO typ1
        ENDIF
 
        kpoint: DO NK=1,NKPTS
          ! used to allocate kpoints to different nodes
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU) .NE. (WDES%COMM_KINTER%NODE_ME-1)) &
                  CYCLE kpoint
        ! kpoint: DO NK=1,1

          CALL PHASE(WDES,NONL_S,NK)
          NPL=WDES%NGVECTOR(NK)

          ! io information
          IUHAM=540
          WRITE(NKINDEX,'(g0)') NK

#ifdef debug
          ! used for checking
          ALLOCATE(FHAM(NPL,NPL),CWORK(NPL*WDES%NRSPINORS))
          DO NP=1, NPL
              CALL CAC_FHAM_LOCAL(SV(1,ISP),CWORK,NP,GRID,NPL,NK,&
                      WDES%NINDPW(1,NK),WDES,G_INDEX_LOC(1,NK))
              DO NNP=1, NPL
              DO ISPINOR=0,WDES%NRSPINORS-1
                  FHAM(NNP,NP)=FHAM(NNP,NP)+CWORK(NNP+NPL*ISPINOR)
              ENDDO
              ENDDO
          ENDDO
 
          OPEN(IUHAM,FILE='LOC_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_'//TRIM(NISP),&
                  STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=16*NPL*NPL)
          WRITE(IUHAM,REC=1) FHAM
          CLOSE(IUHAM)
          DEALLOCATE(FHAM,CWORK)
#endif

          ! Since it's the same for all configuartions, only output once
          IF (LINDEX) THEN
              ! index for calculated G-grids in the whole grid
              OPEN(IUHAM,FILE='INDEX_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_'//TRIM(NISP),&
                      STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=16)
              WRITE(IUHAM,REC=1) GRID%NGX,GRID%NGY,GRID%NGZ,NPL
            
              DO i=1,NPL
                WRITE(IUHAM,REC=i+1) WDES%IGX(i,NK),WDES%IGY(i,NK),WDES%IGZ(i,NK),G_INDEX(i,NK)
              ENDDO
            
              WRITE(IUHAM,REC=NPL+2) WDES%VKPT(1,NK),WDES%VKPT(2,NK)
              WRITE(IUHAM,REC=NPL+3) WDES%VKPT(3,NK),0._q
            
              CLOSE(IUHAM)
         
              ! information needed for non local potential
              NIS=1
              typ:DO NT=1,NONL_S%NTYP
                LMMAXC=NONL_S%LMMAX(NT)
           
                IF (LMMAXC.NE.0) THEN
                  WRITE(NTINDEX,'(g0)') NT
#ifdef wrtphase
                  ! take a lot of space, write position instead.
                  ! phase
                  OPEN(IUHAM,FILE='PHASE_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=16*NPL*NONL_S%NITYP(NT))
           
                  WRITE(IUHAM,REC=1) NONL_S%CREXP(1:NPL,&
                          NIS:NIS+NONL_S%NITYP(NT)-1)
           
                  CLOSE(IUHAM)
#endif
                  ! projection function
                  OPEN(IUHAM,FILE='PROJ_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_'//&
                          TRIM(NISP)//'_'//TRIM(NTINDEX),&
                          STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                          RECL=8)
                  WRITE(IUHAM,Rec=1) PP(NT)%LMAX,LMMAXC
                  CH0=0;CH1=0;CH2=0;CH3=0;
                  DO i=1,PP(NT)%LMAX
                    IF (PP(NT)%LPS(i)==0) THEN
                      CH0=CH0+1
                    ELSEIF (PP(NT)%LPS(i)==1) THEN
                      CH1=CH1+1
                    ELSEIF (PP(NT)%LPS(i)==2) THEN
                      CH2=CH2+1
                    ELSEIF (PP(NT)%LPS(i)==3) THEN
                      CH3=CH3+1
                    ELSE
                      WRITE(*,*) 'warning:L>3 not supported now'
                    ENDIF
                  ENDDO
                  WRITE(IUHAM,Rec=2) 0,CH0
                  WRITE(IUHAM,Rec=3) 1,CH1
                  WRITE(IUHAM,Rec=4) 2,CH2
                  WRITE(IUHAM,Rec=5) 3,CH3
                 
                  DO i=1,LMMAXC
                    DO j=1,NPL
                      WRITE(IUHAM,Rec=j+5+(i-1)*NPL) NONL_S%QPROJ(j,i,NT,NK,1)
                    ENDDO
                  ENDDO
              
                  ! not used
                  IF (WDES%LSPIRAL) THEN
                    DO i=LMMAXC+1,2*LMMAXC
                      DO j=1,NPL
                        WRITE(IUHAM,Rec=j+5+(i-1)*NPL) NONL_S%QPROJ(1:NPL,i,NT,NK,2)
                      ENDDO
                    ENDDO
                  ENDIF
         
                  CLOSE(IUHAM)
         
         
                ENDIF
                NIS=NIS+NONL_S%NITYP(NT)
              ENDDO typ
          ENDIF
 
          IF (DYN%IBRION==0 .OR. DYN%IBRION==23) THEN
              MV='mv ./INDEX_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_* ./PROJ_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_* ./CDIJ_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_* ./CQIJ_'//TRIM(NB_C)//'_'//TRIM(NKINDEX)//'_* ./BEAD_'//TRIM(NSTEP_C)//' 1>/dev/null 2>&1'
              CALL system(MV)
          ENDIF
 
        ENDDO kpoint
 
      ENDDO spin

      IF (DYN%IBRION==0 .OR. DYN%IBRION==23) THEN
          MV='mv ./SV_'//TRIM(NB_C)//'_* ./BEAD_'//TRIM(NSTEP_C)//' 1>/dev/null 2>&1'
          CALL system(MV)
      ENDIF
     
      CALL system_clock(t2)
      IF (MOD(0,WDES%COMM_KINTER%NCPU) .EQ. (WDES%COMM_KINTER%NODE_ME-1)) THEN
        WRITE(*,*) ''
        WRITE(*,*) 'OUTPUT OVER:',(t2-t1)/10._q,'ms'
      ENDIF
     
      IF (LINDEX) THEN
        DEALLOCATE(G_INDEX,G_INDEX_LOC)
      ENDIF

  202 IF (MOD(0,WDES%COMM_INTER%NCPU) .EQ. (WDES%COMM_INTER%NODE_ME-1)) THEN
          IF (MOD(0,WDES%COMM_INB%NCPU) .NE. (WDES%COMM_INB%NODE_ME-1)) THEN
            spin2: DO ISP=1,WDES%NCDIJ
  
              IUHAM=540+WDES%NKPTS*ISP
              WRITE(NISP,'(g0)') ISP
              WRITE(INB_INDEX,'(g0)') WDES%COMM_INB%NODE_ME
              ! information needed for local potential and kinetic energy
              IF (MOD(0,WDES%COMM_KINTER%NCPU) .EQ. (WDES%COMM_KINTER%NODE_ME-1)) THEN
#ifdef realmode
                OPEN(IUHAM,FILE='./SV_'//TRIM(NB_C)//'_'//TRIM(NISP)//'_'//TRIM(INB_INDEX),&
                        STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*GRID%RL%NP)
#else
                OPEN(IUHAM,FILE='./SV_'//TRIM(NB_C)//'_'//TRIM(NISP)//'_'//TRIM(INB_INDEX),&
                        STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=16*GRID%RL%NP)
#endif
                WRITE(IUHAM,REC=1) SV(1:GRID%RL%NP,ISP)
                CLOSE(IUHAM)

                NIS=1;i_in1=1;i_in2=1
                typ3: DO NT=1,NONL_S%NTYP
                LMMAXC=NONL_S%LMMAX(NT)
           
                  IF (LMMAXC.NE.0) THEN
                    WRITE(NTINDEX,'(g0)') NT
                    ! nonlocal pseudopotential strength
#ifdef realmode
                    OPEN(IUHAM,FILE='./CDIJ_'//TRIM(NB_C)//'_'//&
                            TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                            STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                            RECL=8*LMMAXC)
#else
                    OPEN(IUHAM,FILE='./CDIJ_'//TRIM(NB_C)//'_'//&
                            TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                            STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                            RECL=16*LMMAXC)
#endif
                    ! DO i=1,NONL_S%NITYP(NT)
                    DO i=NIS,NIS+NONL_S%NITYP(NT)-1
                      IF (NI_LOCAL(i,WDES%COMM_INB) /=0 ) THEN
                        DO j=1,LMMAXC
                          WRITE(IUHAM,Rec=(i-NIS)*LMMAXC+j) CDIJ(1:LMMAXC,j,i_in1,ISP)
                        ENDDO
                        i_in1=i_in1+1
                      ENDIF
                    ENDDO
           
                    CLOSE(IUHAM)
           
                    ! overlap strength
#ifdef realmode
                    OPEN(IUHAM,FILE='./CQIJ_'//TRIM(NB_C)//'_'//&
                            TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                            STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                            RECL=8*LMMAXC)
#else
                    OPEN(IUHAM,FILE='./CQIJ_'//TRIM(NB_C)//'_'//&
                            TRIM(NISP)//'_'//TRIM(NTINDEX)//'_'//TRIM(INB_INDEX),&
                            STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='DIRECT',&
                            RECL=16*LMMAXC)
#endif
                    ! DO i=1,NONL_S%NITYP(NT)
                    DO i=NIS,NIS+NONL_S%NITYP(NT)-1
                      IF (NI_LOCAL(i,WDES%COMM_INB) /=0 ) THEN
                        DO j=1,LMMAXC
                          WRITE(IUHAM,Rec=(i-NIS)*LMMAXC+j) CQIJ(1:LMMAXC,j,i_in2,ISP)
                        ENDDO
                        i_in2=i_in2+1
                      ENDIF
                    ENDDO
           
                    CLOSE(IUHAM)
                    NIS=NIS+NONL_S%NITYP(NT)
                  ENDIF
                ENDDO typ3
              ENDIF
            ENDDO spin2


          ENDIF

      ENDIF

      ENDIF
      CALL NONL_DEALLOC(NONL_S)
    
      IF (INFO%LCLOSE) THEN
        STOP
      ENDIF
  END SUBROUTINE

  SUBROUTINE GEN_G_INDEX(G_INDEX, G_INDEX_LOC, NPLWV, NPLMAX, WDES, GRID)
    USE prec
    USE mgrid
    USE wave

    IMPLICIT NONE
    TYPE (wavedes) :: WDES
    TYPE (grid_3d) :: GRID
    INTEGER :: NPLWV, NPLMAX
    INTEGER :: G_INDEX(NPLMAX,WDES%NKPTS)
    INTEGER :: G_INDEX_LOC(NPLWV,WDES%NKPTS)

    INTEGER :: NK, IND, NG

    G_INDEX=0
    G_INDEX_LOC=0
    
    kpoint:DO NK=1,WDES%NKPTS
        ! used to allocate kpoints to different nodes
        IF (MOD(NK-1,WDES%COMM_KINTER%NCPU) .NE. (WDES%COMM_KINTER%NODE_ME-1)) &
                CYCLE kpoint

        ! give the relationship in position between G points and grid points
        !OPEN(UNIT=2000,FILE='G_INDEX',STATUS='UNKNOWN')
        !REWIND 2000
        DO IND=1, WDES%NGVECTOR(NK)
            NG= MOD(WDES%IGX(IND,NK)+GRID%NGX,GRID%NGX)+1 +&
                    MOD(WDES%IGY(IND,NK)+GRID%NGY,GRID%NGY)*GRID%NGX_rd + &
                    MOD(WDES%IGZ(IND,NK)+GRID%NGZ,GRID%NGZ)*GRID%NGY*GRID%NGX_rd
            G_INDEX(IND,NK)=NG
            G_INDEX_LOC(NG,NK)=IND
            !WRITE(2000,*) IND,WDES%IGX(IND,NK),WDES%IGY(IND,NK),WDES%IGZ(IND,NK),NG
        ENDDO
        !CLOSE(2000)
    ENDDO kpoint

  END SUBROUTINE
