!***********************************************************************
      subroutine ST_lundcirp()
!***********************************************************************
      use EXP_Global_def,    only: iripple, waves, ncn, nce, etan, cdx, cdy
      USE EXP_transport_def, only: qsx, qsy
      use wave_flowgrid_def, only: wper, wavediss, whgt, wang
      use flow_def, only: rhow, iwet
      use sed_def,  only: rhosed,d50,scalesus,scalebed 
      use size_def, only: ncells
      use geo_def,  only: zb,cell2cell     
            
      implicit none     
      real*8 uvel,vvel,arg,c_angle,w_angle,taucwtb,taucwmtb,epscw,fcf,fwf
      real*8  H_lc,period_lc,db_lc,phi_lc,DEP_lc,vnu_lc,qtot_lc
      real*8  qbs_lc,qss_lc,U0_lc,BDpart,CRCW,USTC,USTW,single_d50  
      real*8 xRHOW,xRHOSED,xWSFALL,xTAUCR 
      integer i,iiripple
                   
!lund-cirp formula
! I = COUNTER (?; NOT USED)
! H = WAVE HEIGHT
! PERIOD = WAVE PERIOD
! DB = WAVE ENERGY DISSIPATION DUE TO BREAKING
! U0 = CURRENT SPEED (MAGNITUDE)
! DEP = WATER DEPTH
! D50 = MEDIAN GRAIN SIZE   ( use D50() )
! PHI_M = ANGLE BETWEEN WAVE AND CURRENT DIRECTION
! VNU = KINEMATIC VISCOSITY  (set to 0.000001 in old to HS map file)
! RHOW = DENSITY OF WATER
! RHOS = DENSITY OF SEDIMENT  (use RHOSED)
! IRIPPLE = SWITCH FOR RIPPLE CALCULATIONS (=1 IMPLIES RIPPLES)
! SCALEBED = CALIBRATION FACTOR FOR BED LOAD (=1.0 IMPLIES ORIGINAL FORMULA)
! SCALESUS = CALIBRATION FACTOR FOR SUSPENDED LOAD (=1.0 IMPLIES ORIGINAL FORMULA)

      !need local variables for subroutine calls
      xRHOW = RHOW
      xRHOSED = RHOSED
      iiripple = iripple

      call update_Q_forSedTrans()   !this fixes WSE BC cells lower vlocity from reducing Sed Trans when
                                    !extrap < 1.0 is used                                     
      IF(waves) then
!!$omp parallel do     
!!$omp+ private(H_lc,PERIOD_lc,DB_lc,DEP_lc,uvel,vvel,U0_lc,hrms,gammatt,
!     . arg,c_angle,w_angle,PHI_lc,ncn,nce,single_d50,WLtt,qtot_lc,
!     . qba,qsm,qsr,qbm,qall,alfa,um,ur,vnu_lc,qbs_lc,Qss_lc,
!     . bdpart,crcw,xWSFALL,ustc,ustw,xtaucr,taucwtb,
!     . taucwmtb,epscw,fcf,fwf)
        DO I=1,NCELLS
          qsx(i) = 0.
          qsy(i) = 0.
          ncn = cell2cell(1,i)     !cwr- 072909
          nce = cell2cell(2,i)     !cwr- 072909
          IF(iwet(i) .eq. 1) then
            PERIOD_lc = Wper(i)
            DB_lc = abs(wavediss(i))          
            DEP_lc = -zb(i) + etan(i)
            H_lc = Whgt(i) !min(Whgt(i),0.78*Dep_lc)             
            uvel = cdx(i) !((qxn(i)+qxn(nce))/2. + 1.e-10)/dep_lc
            vvel = cdy(i) !((qyn(i)+qyn(ncn))/2.)/dep_lc
            U0_lc = sqrt(uvel**2+vvel**2)
            H_lc = min(0.7*Dep_lc,H_lc)
            arg = vvel/uvel
            c_angle = atan2(vvel,uvel)
            w_angle = Wang(i) !Wdir(i)*pi/180.0
            PHI_lc = abs(w_angle - c_angle)
            SINGLE_D50 = D50(I) !CHANGE SINGLED50 FOR EACH VARIABLE D50(I)
            VNU_lc = .000001d0
            CALL LUNDCIRP_(H_lc,PERIOD_lc,DB_lc,U0_lc,DEP_lc,     &
                 single_D50,PHI_lc,VNU_lc,xRHOW,xRHOSED,QBS_lc,   & !Variable D50 here
                 QSS_lc,iIRIPPLE,BDpart,CRCW,xWSFALL,             &
                 USTC,USTW,xTAUCR,TAUCWTB,TAUCWMTB,EPSCW,FCF,FWF)
            QSS_lc = SCALESUS*QSS_lc
            QBS_lc = SCALEBED*QBS_lc
            QTOT_lc = QSS_lc + QBS_lc   
                  
            !regular LC total load - in direction of depth-averaged current
            qsx(i) = (uvel/U0_lc)*qtot_lc
            qsy(i) = (vvel/U0_lc)*qtot_lc
          endif !wet         
        ENDDO
!!$omp end parallel do        
      ELSE
!!$omp parallel do     
!!$omp+ private(DEP_lc,uvel,vvel,U0_lc,ncn,nce,QTOT_lc)
!!$omp+ private(vnu_lc,h_lc,period_lc,db_lc,phi_lc,single_D50)
!!$omp+ private(bdpart,crcw,xWSFALL,ustc,ustw,xtaucr,taucwtb)
!!$omp+ private(taucwmtb,epscw,fcf,fwf,qbs_lc,Qss_lc)
        DO I=1,NCELLS
          qsx(i) = 0.
          qsy(i) = 0.
          ncn = cell2cell(1,i)    !cwr- 072909
          nce = cell2cell(2,i)    !cwr- 072909
          IF(iwet(i) .eq. 1) then
            VNU_lc = .000001d0
            H_lc = 0.0                                  
            PERIOD_lc = 0.0
            DB_lc = 0.0
            PHI_lc = 0.0          
            DEP_lc = -zb(i) + etan(i)
            uvel = cdx(i) !((qxn(i)+qxn(nce))/2. + 1.e-10)/dep_lc
            vvel = cdy(i) !((qyn(i)+qyn(ncn))/2.)/dep_lc
            U0_lc = sqrt(uvel**2+vvel**2)
            SINGLE_D50 = D50(I)!CHANGE SINGLED50 FOR EACH VARIABLE D50(I)
            CALL LUNDCIRP_(H_lc,PERIOD_lc,DB_lc,U0_lc,DEP_lc,    &
                 single_D50,PHI_lc,VNU_lc,xRHOW,xRHOSED,QBS_lc,  &  !Variable D50 here
                 QSS_lc,iIRIPPLE,BDpart,CRCW,xWSFALL,            &
                 USTC,USTW,xTAUCR,TAUCWTB,TAUCWMTB,EPSCW,FCF,FWF)
            QSS_lc = SCALESUS*QSS_lc
            QBS_lc = SCALEBED*QBS_lc
            QTOT_lc = QSS_lc + QBS_lc
         
            qsx(i) = (uvel/U0_lc)*qtot_lc
            qsy(i) = (vvel/U0_lc)*qtot_lc              
          endif
        ENDDO    
        
!!$omp end parallel do        
      ENDIF

      return
      end subroutine

