!***********************************************************************
      subroutine ReadCardFile_EXP
!***********************************************************************
      use EXP_Global_def,  only: drydep, dt, iadv, imix, a0, thetac, iripple, advcoef, rainfall, isedform
      use EXP_Global_def,  only: dtsed, dtmorph, dtsalt, dtheat, rf_filename, maxunit, rf_unit, rf_frac_ro
      USE EXP_bndcond_def, only: modify_h, hmod
      USE EXP_transport_def,  only: rate_avalanche, chparms, saltsimd, cohes_flow_bc_file, cohes_flow_bc
      use EXP_Structures_def, only: srm, structures, srm_on, srm, cul, cul_on
      use prec_def, only: ikind
      use flow_def, only: hmin
      use comvarbl, only: ctlfile,dtime  
                 
      implicit none
      !local vars
      integer i,j,itimes,ios, ierr
      real(ikind) Adummy
      character :: astring*50   
      CHARACTER*30  CARDNAME    
      CHARACTER*12  STRINGNAME   
      CHARACTER*3   CSWITCH   
      CHARACTER*200  ADUM      
      LOGICAL   FOUND, bool 
      LOGICAL :: UNRECOGNIZED = .FALSE.
 
      !rhow = densit   !rhow now used in all code - MEB 11/13
      drydep = hmin   !hmin used in implicit code
      dt = dtime      !dtime used in explicit code

!********************
!    Set Defaults
!********************

      !hydrodyanmics
      IADV = 1
      IMIX = 1  
      !sediment transport
      A0 = -1.D0
      THETAC = -1.D0
      IRIPPLE = 1   
      rate_avalanche = 0.01     
      !BC extrabploation
      ADVCOEF = 0.0  !0.90d0 !1.0d0 
      !cohesive
      CHparms%Dfac = 1.0
      !COHD_PATH = "COHES_DIFF"    
      !RAINFALL
      rainfall = .false.   

      OPEN(1, FILE=CTLFILE, STATUS='OLD') 
    
!********************
!*  Read Cards
!********************
      DO   !UNTIL END OF FILE
        READ (1,*,IOSTAT=IOS) CARDNAME
        IF (IOS.NE.0) EXIT
        IF (CARDNAME(1:1).EQ."!" .OR. CARDNAME(1:1).EQ." ") CYCLE
        BACKSPACE(1)  !ALWAYS BACK UP TO READ THE REST OF THE LINE FOR OTHER OPTIONS

        SELECT CASE (CARDNAME)
             
!--START OF HYDRO CARDS
          CASE ('USE_ADVECTION_TERMS')
            call card_boolean(1,bool,ierr)
            !READ(1,*) CARDNAME, CSWITCH
            IF (.not.BOOL) IADV = 0            !DEFAULT = 'ON'

          CASE ('USE_MIXING_TERMS')
            call card_boolean(1,bool,ierr)              
            !READ(1,*) CARDNAME, CSWITCH
            IF (.not.BOOL) IMIX = 0            !DEFAULT = 'ON'

!--START OF SED TRANS CARDS   
          CASE ('SED_TRAN_FORMULATION')
            READ(1,*) CARDNAME, STRINGNAME
            IF     (STRINGNAME.EQ.'WATANABE')  THEN 
              ISEDFORM=1
            ELSEIF (STRINGNAME.EQ.'LUND_CIRP') THEN 
              ISEDFORM=2
            ELSEIF (STRINGNAME.EQ.'EXNER')     THEN 
              ISEDFORM=2              
            ELSEIF (STRINGNAME.EQ.'A-D')       THEN 
              ISEDFORM=3
            ELSEIF (STRINGNAME.EQ.'COHES')     THEN  
              ISEDFORM=5
            ENDIF
            STRINGNAME = '        '
          
          CASE ('SED_TRAN_CALC_INTERVAL')
            READ(1,*) CARDNAME, DTSED                                   
!
          CASE ('MORPH_UPDATE_INTERVAL')
            READ(1,*) CARDNAME, DTMORPH                                        
          
          CASE ('A_COEFFICIENT_WATANABE')
            READ(1,*) CARDNAME, A0

          CASE ('THETA_CRITICAL')
            READ(1,*) CARDNAME, THETAC
            
          CASE ('AVALANCHE_RATE')                          
            READ(1,*) CARDNAME, rate_avalanche
            
!--START OF SALINITY CARDS - MEB 11/11/2008          
        
          CASE ('SALINITY_DENSITY')                                               !Added 11/20/2009 - Chris Reed
            READ (1,*) CARDNAME, CSWITCH
            if(CSWITCH.eq.'OFF') then
              saltsimD = .false. !Turn off density changes due to salinity
            else
              saltsimD = .true.  !Default is TRUE, so if anything other than OFF, leave it on.
            endif
            CSWITCH='   '           

          CASE ('SALINITY_CALC_INTERVAL')
            READ (1,*) CARDNAME, DTSALT
          
          CASE ('TEMPERATURE_CALC_INTERVAL')
            READ (1,*) CARDNAME, DTHEAT
  
!--START OF BC CARDS

          CASE ('BOUND_ADVECTION_EXTRAP_COEFFIC','BOUND_ADVECTION_EXTRAP_COEFFICIENT')
            READ(1,*) CARDNAME, ADVCOEF           

          CASE ('MODIFY_BC_CONSTANT')
            READ(1,*) CARDNAME, Adummy
            MODIFY_H = .TRUE.
            ALLOCATE(HMOD%DTIME(2),HMOD%DVALUE(2))
            HMOD%DTIME(1)  = 0.0    ; HMOD%DTIME(2)  = 99999999.0
            HMOD%DVALUE(1) = Adummy ; HMOD%DVALUE(2) = Adummy
            HMOD%INC = 1
            
          CASE ('MODIFY_BC_VARIABLE')
            READ(1,*) CARDNAME, ASTRING
            INQUIRE(FILE=TRIM(ASTRING),EXIST=FOUND)
            IF (FOUND) THEN
              MODIFY_H = .TRUE.
              OPEN(2,FILE=TRIM(ASTRING),STATUS='OLD')
              ITIMES = 0
              DO 
                READ(2,*,IOSTAT=IOS) ADUM
                IF (IOS.NE.0) EXIT
                ITIMES = ITIMES + 1
              ENDDO
              CLOSE(2)
              OPEN(2,FILE=TRIM(ASTRING),STATUS='OLD')
              ALLOCATE(HMOD%DTIME(ITIMES),HMOD%DVALUE(ITIMES))
              DO I=1,ITIMES
                READ(2,*) HMOD%DTIME(I),HMOD%DVALUE(I)
              ENDDO
              CLOSE(2)
              HMOD%INC=1
            ENDIF
            
!--START OF COHESIVE CARDS
          CASE ('COHESIVE_E')
            READ (1,*) CARDNAME,CHparms%E
            ISEDFORM = 5                           !added 05/15/09 reed
          
          CASE ('COHESIVE_TCRIT_EROSION')
            READ (1,*) CARDNAME,CHparms%Tcrit_E
          
          CASE ('COHESIVE_TCRIT_DEPOSITION')
            READ (1,*) CARDNAME,CHparms%Tcrit_D
          
          CASE ('COHESIVE_WS_MAX')
            READ (1,*) CARDNAME,CHparms%ws_max   

          CASE ('COHESIVE_CONC_MAX')
            READ (1,*) CARDNAME,CHparms%c_max

          CASE ('COHESIVE_CONC_PEAK')
            READ (1,*) CARDNAME,CHparms%c_peak
          
          CASE ('COHESIVE_WSE_BC')
            READ (1,*) CARDNAME,CHparms%wse_bc  
               
          CASE ('COHESIVE_DIFF_FAC')
            READ (1,*) CARDNAME,CHparms%Dfac 
 
          CASE ('COHES_FLOW_BC_INPUT')
            READ(1,*) CARDNAME,COHES_FLOW_BC_FILE
            cohes_flow_bc = .true. 
          
          CASE ('COHES_TCRIT_VARIABLE')
            READ(1,*) CARDNAME, CHparms%numDEPTHS
            BACKSPACE(1)
            CHparms%Tcrit_variable = .true.
            IF (.NOT.ALLOCATED(CHparms%DEPTH)) ALLOCATE(CHparms%DEPTH(CHparms%numDEPTHS))
            IF (.NOT.ALLOCATED(CHparms%Tcrit_Ev)) ALLOCATE(CHparms%Tcrit_Ev(CHparms%numDEPTHS))
            IF (.NOT.ALLOCATED(CHparms%Tcrit_Dv)) ALLOCATE(CHparms%Tcrit_Dv(CHparms%numDEPTHS))        
            READ(1,*) CARDNAME, CHparms%numDEPTHS, (CHparms%DEPTH(J),CHparms%Tcrit_Ev(J),CHparms%Tcrit_Dv(J),J=1,CHparms%numDEPTHS)
            write(*,*)"COHESIVE: Variable E and D Invoked"
            DO J=1,CHparms%numDEPTHS
              write(*,*)CHparms%DEPTH(J),CHparms%Tcrit_Ev(J),CHparms%Tcrit_Dv(J)
            ENDDO

!--------START OF STRUCTURE CARDS - RUBBLE MOUND
          CASE ('STRUCT_RM')
            READ(1,*) CARDNAME, SRM%NCELLS
            BACKSPACE(1)
            STRUCTURES = .true.
            SRM_ON = .true.
            IF (.NOT.ALLOCATED(SRM%CELLS)) ALLOCATE(SRM%CELLS(SRM%NCELLS))
            READ(1,*)CARDNAME, SRM%NCELLS, (SRM%CELLS(J),J=1,SRM%NCELLS)
            IF (.NOT.ALLOCATED(SRM%POR)) ALLOCATE(SRM%POR(SRM%NCELLS))
            SRM%POR = 0.2             
            IF (.NOT.ALLOCATED(SRM%HGT)) ALLOCATE(SRM%HGT(SRM%NCELLS))
            SRM%HGT = -1.0           
            IF (.NOT.ALLOCATED(SRM%HC))ALLOCATE(SRM%HC(SRM%NCELLS))
            SRM%HC = 2000   
            IF (.NOT.ALLOCATED(SRM%A)) ALLOCATE(SRM%A(SRM%NCELLS))
            SRM%a = 0.02  
            IF (.NOT.ALLOCATED(SRM%B)) ALLOCATE(SRM%B(SRM%NCELLS))
            SRM%B = 3.0                                    
            IF (.NOT.ALLOCATED(SRM%ival)) ALLOCATE(SRM%ival(SRM%NCELLS))
            write(*,*)'Rubble Mound Structures (RM) Invoked'
            write(*,*)"number of RM cells = ",SRM%ncells
                          
          CASE ('STRUCT_RM_POR')
            READ(1,*) CARDNAME, Adummy
            do J=1,SRM%NCELLS
              SRM%POR(J)= Adummy
            ENDDO
            
          CASE ('STRUCT_RM_BASE')
            READ(1,*) CARDNAME, SRM%BASE
            
          CASE ('STRUCT_RM_HC')
            READ(1,*) CARDNAME, Adummy
            do J=1,SRM%NCELLS
              SRM%HC(J)= Adummy
            ENDDO   
          
          CASE ('STRUCT_RM_A')
            READ(1,*) CARDNAME, Adummy
            do J=1,SRM%NCELLS
              SRM%A(J)= Adummy
            ENDDO         
          
          CASE ('STRUCT_RM_B')
            READ(1,*) CARDNAME, Adummy
            do J=1,SRM%NCELLS
              SRM%B(J)= Adummy
            ENDDO          

!--------START OF STRUCTURE CARDS - CULVERT   
          CASE ('STRUCT_CUL')
            READ(1,*) CARDNAME, cul%NUM
            BACKSPACE(1)
            STRUCTURES = .true.
            CUL_ON = .true.
            IF (.NOT.ALLOCATED(cul%CELLS1)) ALLOCATE(cul%CELLS1(cul%NUM))
            IF (.NOT.ALLOCATED(cul%CELLS2)) ALLOCATE(cul%CELLS2(cul%NUM)) 
            READ(1,*)CARDNAME,ADUMMY, (cul%CELLS1(J),cul%CELLS2(J),J=1,cul%NUM)
            write(*,*)cul%num,cul_on
            do i=1,cul%num
              write(*,*)'cul pair = ',cul%cells1(i),cul%cells2(i)
            enddo
            IF (.NOT.ALLOCATED(cul%INVERT)) ALLOCATE(cul%INVERT(cul%NUM))
            IF (.NOT.ALLOCATED(cul%DIA))    ALLOCATE(cul%DIA(cul%NUM))
            IF (.NOT.ALLOCATED(cul%FRIC))   ALLOCATE(cul%FRIC(cul%NUM))
            IF (.NOT.ALLOCATED(cul%LENGTH)) ALLOCATE(cul%LENGTH(cul%NUM))
            cul%INVERT = 1.0            
            cul%DIA = 0.67           
            cul%FRIC = 0.025           
            cul%LENGTH = 5
                                       
          CASE ('STRUCT_CUL_INVERT')
            READ(1,*) CARDNAME, Adummy
            do J=1,cul%NUM
              cul%INVERT(J)= Adummy
            ENDDO
            
          CASE ('STRUCT_CUL_DIA')
            READ(1,*) CARDNAME, Adummy
            do J=1,cul%NUM
              cul%DIA(J)= Adummy
            ENDDO
            
          CASE ('STRUCT_CUL_FRIC')
            READ(1,*) CARDNAME, Adummy
            do J=1,cul%NUM
              cul%FRIC(J)= Adummy
            ENDDO  
            
          CASE ('STRUCT_CUL_LEN')
            READ(1,*) CARDNAME, Adummy
            do J=1,cul%NUM
              cul%LENGTH(J)= Adummy
            ENDDO             
                     
! START OF RAINFALL CARDS
          CASE ('RF_FILENAME')
            READ (1,*) CARDNAME, RF_FILENAME
            maxunit=maxunit+1
            RF_unit=maxunit
            RAINFALL=.true.
            open(unit=RF_unit,FILE=trim(RF_FILENAME),STATUS='old')
            write(*,*)'rainfall enabled using file =',trim(RF_FILENAME)
            write(*,*)'thus rainfall set to ',rainfall
          
          CASE ('RF_FRAC_RO')
            READ (1,*) CARDNAME, RF_FRAC_RO
            if(RF_FRAC_RO.gt.1.0) then
              write(*,*)'RF_FRAC_RO 0 <= and <= 1.0'
              stop
              write(*,*)'dry cell rainfall ro set to ',rf_frac_ro
            endif
 
          CASE ('END_PARAMETERS      ')
            READ(1,*) CARDNAME

          CASE DEFAULT
            READ(1,*) CARDNAME
            !WRITE(*,*)' WARNING: UNRECOGNIZED CARD: '//TRIM(CARDNAME)
            UNRECOGNIZED = .TRUE.
         
        END SELECT
      ENDDO

      CLOSE(1)
        
      return
      end subroutine
    
      