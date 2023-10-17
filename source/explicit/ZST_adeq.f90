!***********************************************************************
    subroutine ST_adeq()
!***********************************************************************
      use EXP_Global_def,    only: iripple, waves, etan, cdx, cdy, linktodummies
      USE EXP_transport_def, only: qsx, qsy, adss, tsed_elapse
      use wave_flowgrid_def, only: whgt, wper, wavediss, wang
      use flow_def, only: rhow, iwet, vis
      use size_def, only: ncells, ncellsd
      use sed_def, only: rhosed,d50,scalesus,scalebed 
      use geo_def, only: zb
                
      implicit none
      real*8 uvel,vvel,arg,c_angle,w_angle,taucwtb,taucwmtb,epscw,fcf,fwf,term
      integer i,ii,jj    
      real*8  H_lc,period_lc,db_lc,phi_lc,DEP_lc,vnu_lc,qtot_lc
      real*8  qbs_lc,qss_lc,U0_lc,BDpart,CRCW,USTC,USTW,single_D50       
      real*8 xRHOW,xRHOSED,xWSFALL,xTAUCR 
      integer iiripple
            
      !need local variables for subroutine calls
      xRHOW = RHOW
      xRHOSED = RHOSED            
      iiripple = iripple
      
      call update_Q_forSedTrans()   !this fixes WSE BC cells lower vlocity from reducing 
                                    !Sed Trans when extrap < 1.0 is used                               
                                    
      !advection-diffusion formula for sediment tranport
!$omp parallel 
      !IF(waves .or. wavescwr_sim) then
      IF(waves) then   

!$omp do                                                           &
!$omp& private(H_lc,PERIOD_lc,DB_lc,DEP_lc,uvel,vvel,U0_lc)        &
!$omp& private(arg,c_angle,w_angle,PHI_lc,term,single_d50,qtot_lc) &
!$omp& private(vnu_lc,qbs_lc,Qss_lc)                               &
!$omp& private(bdpart,crcw,xWSFALL,ustc,ustw,xtaucr,taucwtb)       &
!$omp& private(taucwmtb,epscw,fcf,fwf)
       DO I=1,NCELLS
         qsx(i) = 0                                                     !cr 8/25/08
         qsy(i) = 0                                                     !cr 8/25/08
         adss(i)%eros=0                                                 !cr 8/25/08
         adss(i)%depo =0                                                !cr 8/25/08
         adss(i)%diffC=0                                                !cr 8/25/08
         if(iwet(i) .eq. 1) then                                        !cr 8/25/08
           H_lc = Whgt(i)
           PERIOD_lc = Wper(i)
           DB_lc = abs(wavediss(i))
           DEP_lc = -zb(i) + etan(i)
           uvel = cdx(i) !((qxn(i)+qxn(cell2cell(2,i)))/2. + 1.e-10)/dep_lc
           vvel = cdy(i) !((qyn(i)+qyn(cell2cell(1,i)))/2.)/dep_lc
           U0_lc = sqrt(uvel**2+vvel**2)
           arg = vvel/uvel
           c_angle = atan2(vvel,uvel)
           w_angle = Wang(i) !Wdir(i)*pi/180.0
           PHI_lc = abs(w_angle - c_angle)
           SINGLE_D50 = D50(I)
           VNU_lc = .000001d0
           
           H_lc = min(0.7*Dep_lc,H_lc) 
          
           CALL LUNDCIRP_(H_lc,PERIOD_lc,DB_lc,U0_lc,DEP_lc,         &    !01/08/09 - Variable D50
             SINGLE_D50,PHI_lc,VNU_lc,xRHOW,xRHOSED,QBS_lc,QSS_lc,   &
             iIRIPPLE,BDpart,CRCW,xWSFALL,USTC,USTW,                 &
             xTAUCR,TAUCWTB,TAUCWMTB,EPSCW,FCF,FWF)
     
           QSS_lc = SCALESUS*QSS_lc
           QBS_lc = SCALEBED*QBS_lc
           QTOT_lc = QSS_lc + QBS_lc
           
           !regular LC total load - in direction of depth-averaged current    
           qsx(i) = (uvel/U0_lc)*qbs_lc
           qsy(i) = (vvel/U0_lc)*qbs_lc
                
           BDpart = BDpart/DEP_lc
           term = xWSFALL/BDPART
           IF(TSED_ELAPSE*xWSFALL/(BDPART*Dep_lc).GE.1.0)   term = 0.99*Dep_lc/tsed_elapse
           ADSS(i)%depo = ADSS(i)%conc*TERM
           ADSS(i)%eros = CRCW*xWSFALL
           ADSS(i)%DiffC = vis(i)
         endif  ! if wet                                                !cr 8/25/08
       ENDDO
      
!$omp end do

      ELSE
!$omp do                                                       &
!$omp& private(DEP_lc,uvel,vvel,U0_lc,term,QTOT_lc)            &
!$omp& private(vnu_lc,h_lc,period_lc,db_lc,phi_lc,single_D50)  &
!$omp& private(bdpart,crcw,xWSFALL,ustc,ustw,xtaucr,taucwtb)   &
!$omp& private(taucwmtb,epscw,fcf,fwf,qbs_lc,Qss_lc)                      
       DO I=1,NCELLS
         qsx(i) = 0                                                     !cr 8/25/08
         qsy(i) = 0                                                     !cr 8/25/08
         adss(i)%eros = 0                                               !cr 8/25/08
         adss(i)%depo =0                                                !cr 8/25/08
         adss(i)%diffC=0                                                !cr 8/25/08
         if(iwet(i) .eq. 1) then                                                !cr 8/25/08
           VNU_lc = .000001d0
           H_lc = 0.0                                   
           PERIOD_lc = 0.0
           DB_lc = 0.0
           PHI_lc = 0.0        
           DEP_lc = -zb(i) + etan(i)
           uvel = cdx(i) !((qxn(i)+qxn(cell2cell(2,i)))/2. + 1.e-10)/dep_lc
           vvel = cdy(i) !((qyn(i)+qyn(cell2cell(1,i)))/2.)/dep_lc
           U0_lc = sqrt(uvel**2+vvel**2)
           SINGLE_D50 = D50(I)!CHANGE SINGLED50 FOR EACH VARIABLE D50(I)
            
           CALL LUNDCIRP_(H_lc,PERIOD_lc,DB_lc,U0_lc,DEP_lc,    &    
                single_D50,PHI_lc,VNU_lc,xRHOW,xRHOSED,QBS_lc,  &  !Variable D50 here
                QSS_lc,iIRIPPLE,BDpart,CRCW,xWSFALL,            & 
                USTC,USTW,xTAUCR,TAUCWTB,TAUCWMTB,EPSCW,FCF,FWF)
      
           QSS_lc = SCALESUS*QSS_lc
           QBS_lc = SCALEBED*QBS_lc
           QTOT_lc = QSS_lc + QBS_lc
               
           qsx(i) = (uvel/U0_lc)*qbs_lc
           qsy(i) = (vvel/U0_lc)*qbs_lc  
           BDpart = BDpart/DEP_lc
           term = xWSFALL/(BDPART)
           IF(TSED_ELAPSE*xWSFALL/(BDPART*Dep_lc).GE.1.0)    term =0.99*Dep_lc/tsed_elapse
           ADSS(i)%depo = ADSS(i)%conc*TERM
           ADSS(i)%eros = CRCW*xWSFALL    
           ADSS(i)%diffC = vis(i)/0.7 !5.93*dep_lc*USTC
         endif !if wet                                                  !cr 8/25/08
       ENDDO  
!$omp end do       
      ENDIF
!$omp end parallel      
      
      !copy diffC values to dummy cells
      ii=0
      do i=ncells+1,ncellsD
        ii=ii+1
        jj = linktodummies(ii)
        ADSS(i)%diffC = ADSS(jj)%diffC
      enddo

      return
      end subroutine