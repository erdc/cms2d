!***********************************************************************
    subroutine Init_BalanceCheck_tel()
!***********************************************************************
    use BalanceCheck_def
    use sed_def,  only: sedtrans
    use sal_def,  only: saltrans,sal  
    use comvarbl, only: timehrs,dtime    
    
    implicit none
    character*80 filename
    real*8  Error,Error_Percent
     
    Error = 0
    Error_Percent = 0
    
    write(*,*)'Check mass/mom balance = ',Balcheck
    write(*,*)'Check interval = ',BalCheck_int,' seconds'
    BalanceCheck_time = 0
    MassFluxin = 0
    MassFluxout = 0
    MassTotal = 0
    MassTotal_old = 0
          
    MomFluxin = 0
    MomFluxout = 0
    MomTotal = 0
    MomTotal_old = 0         
    MomFricLoss = 0
    MomWindAdd = 0
    MomPressAdd = 0
          
    SMassFluxin = 0
    SMassFluxout = 0
    SMassTotal = 0
    SMassTotal_old = 0
    SbedAdd = 0
    SbedLoss = 0    
          
    Salfluxin = 0
    SalFluxout = 0
    SalmassTotal = 0
    SalmassTotal_old = 0
          
    filename = 'Check_MassBalance.csv'
    call balancecheck_tel_mass()  !get intial water volume
    open(unit=961,file=trim(filename))
    write(961,*)'Time(hrs),FluxIn(m^3),Fluxout(m^3),TotMass(m^3),Erorr(m^3),%Error(%)'
    write(961,"(6(e19.9,','))")timehrs,MassFluxin,MassFluxout,MassTotal,Error,Error_percent
    MassTotal_old = MassTotal
          
    filename = 'Check_MomentumBalance.csv'
    call balancecheck_tel_mom()  !get intial momentum          
    open(unit=962,file=trim(filename))  
    write(962,"(a130)")'Time(hrs),Fluxin(m^4/s),Fluxout(m^4/s),TotMom(m^4/s),FricLoss(m^4/s),WindAdd(m^4/s),PressAdd(m^4/s),Erorr(m^4/s),%Error(%)'
    write(962,"(9(e19.9,','))")timehrs,MomFluxin,MomFluxout,MomTotal,MomFricLoss,MomWindAdd,MomPressAdd,Error,Error_percent 
    MomTotal_old = MomTotal          
          
    if(Saltrans) then
      filename = 'Check_SalinityBalance.csv'
      call balancecheck_tel_sal_mass()  !get initial salinty mass
      open(unit=963,file=trim(filename)) 
      write(963,*)'Time(hrs),SalIn(mg),Salout(mg),TotSal(mg),Erorr(mg),%Error(m^3)'
      write(963,"(6(e19.9,','))")timehrs,Salfluxin,SalFluxout,SalmassTotal,Error,Error_percent 
      SalMassTotal_old = SalMassTotal
    endif
          
    if(Sedtrans)then
      filename = 'Check_SedimentBalance.csv'
      !call balancecheck_tel_sed()  !get initial sediment mass          
      open(unit=964,file=trim(filename)) 
      write(964,*)'Time(hrs),SedIn(*),Sedout(*),TotSed(*),SedDep(*),SedEros(*),Erorr(*),%Error(%)'
      write(964,"(8(e19.9,','))")timehrs,SMassFluxin,SMassFluxout,SMassTotal,SbedAdd,SbedLoss,Error,Error_percent 
      SMassTotal_old = SMassTotal
    endif 
    
    return
    end subroutine !Init_BalanceCheck

!***********************************************************************
    subroutine balancecheck_tel()
!***********************************************************************
    use EXP_Global_def,  only: num_ext_n, num_ext_s, num_ext_w, num_ext_e, fac_dw, fac_uw
    USE EXP_bndcond_def, only: qstringexp, ext_n, ext_s, ext_w, ext_e
    use EXP_TELESCOPING, only: xface_qn, yface_qn, xface_length, yface_length, xface_q, yface_q, xface_vel, yface_vel
    use comvarbl, only: timehrs,dtime
    use bnd_def,  only: nqstr, q_str
    use geo_def,  only: dx,dy
    use size_def, only: ncells,ncellsD
    use geo_def,  only: zb 
    use sed_def,  only: sedtrans
    use sal_def,  only: saltrans,sal       
    use BalanceCheck_def
      
    implicit none
    integer i,j,ii,iid,IDO
    real*8 Error,Error_percent
      
    BalanceCheck_time = BalanceCheck_time + dtime
      
!________________________Water Volume____________________________________
    call balancecheck_tel_mass()     
!_________________water fluxes at Q BC strings___________________________    
    if(nQstr .gt. 0) then    
      do i = 1,nQstr  !for each cell string           
        IF(QstringEXP(i)%vface) THEN          
          if(QstringEXP(I)%sgn .eq. -1) then !flow on north side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_qn(IID) .gt. 0) then
                MassFluxout = MassFluxout + yface_qn(IID)*dtime*yface_Length(IID)
              else
                MassFluxin = MassFluxin - yface_qn(IID)*dtime*yface_Length(IID)
              endif
            enddo
          else !flow on south side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_qn(IID) .gt. 0) then
                MassFluxin = MassFluxin + yface_qn(IID)*dtime*yface_Length(IID)
              else
                MassFluxout = MassFluxout - yface_qn(IID)*dtime*yface_Length(IID)
              endif
            enddo                  
          endif  
        ELSE            
          if(QstringEXP(I)%sgn .eq. -1) then !flow on east side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_qn(IID) .gt. 0) then
                MassFluxout = MassFluxout + xface_qn(IID)*dtime*xface_Length(IID)
              else
                MassFluxin = MassFluxin - xface_qn(IID)*dtime*xface_Length(IID)
              endif     
            enddo
          else  !flow on west side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_qn(IID) .gt. 0) then
                MassFluxin = MassFluxin + xface_qn(IID)*dtime*xface_Length(IID)
              else
                MassFluxout = MassFluxout - xface_qn(IID)*dtime*xface_Length(IID)
              endif     
            enddo                
          endif
        ENDIF    
      enddo ! end of NQdriver
    endif  !Q_single    

!___________________water fluxes at tidal BCs_______________________________         
    do i=1,num_ext_N
      if(yface_q(ext_N(i,2)) .gt. 0) then
         MassFluxout = MassFluxout + yface_q(ext_N(i,2))*dtime*yface_Length(ext_N(i,2))
      else
         MassFluxin = MassFluxin - yface_q(ext_N(i,2))*dtime*yface_Length(ext_N(i,2))
      endif
    enddo
    do i=1,num_ext_S
      if(yface_q(ext_S(i,2)) .gt. 0) then
        MassFluxin = MassFluxin + yface_q(ext_S(i,2))*dtime*yface_Length(ext_S(i,2))
      else
        MassFluxout = MassFluxout - yface_q(ext_S(i,2))*dtime*yface_Length(ext_S(i,2))
      endif
    enddo
    do i=1,num_ext_W
      if(xface_q(ext_W(i,2)) .gt. 0) then
        MassFluxin = MassFluxin + xface_q(ext_W(i,2))*dtime*xface_Length(ext_W(i,2))
      else
        MassFluxout = MassFluxout - xface_q(ext_W(i,2))*dtime*xface_Length(ext_W(i,2))
      endif
    enddo
    do i=1,num_ext_E
      if(xface_q(ext_E(i,2)) .gt. 0) then
        MassFluxout = MassFluxout + xface_q(ext_E(i,2))*dtime*xface_Length(ext_E(i,2))
      else
        MassFluxin = MassFluxin - xface_q(ext_E(i,2))*dtime*xface_Length(ext_E(i,2))
      endif
    enddo  
    
!________________________Total Momentum____________________________________
    call balancecheck_tel_mom()     
!_________________momentum fluxes at Q BC strings______NOT DONE_____________ 
    if(nQstr .gt. 0) then
      do i = 1,nQstr  !for each cell string           
        IF(QstringEXP(i)%vface) THEN          
          if(QstringEXP(I)%sgn .eq. -1) then !flow on north side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_qn(IID) .gt. 0) then
                MassFluxout = MassFluxout + yface_qn(IID)*dtime*yface_Length(IID)*yface_vel(IID)
              else
                MassFluxin = MassFluxin + yface_qn(IID)*dtime*yface_Length(IID)*yface_vel(IID)
              endif
            enddo
          else !flow on south side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_qn(IID) .gt. 0) then
                MassFluxin = MassFluxin + yface_qn(IID)*dtime*yface_Length(IID)*yface_vel(IID)
              else
                MassFluxout = MassFluxout + yface_qn(IID)*dtime*yface_Length(IID)*yface_vel(IID)
              endif
            enddo                  
          endif  
        ELSE            
          if(QstringEXP(I)%sgn .eq. -1) then !flow on east side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_qn(IID) .gt. 0) then
                MassFluxout = MassFluxout + xface_qn(IID)*dtime*xface_Length(IID)*xface_vel(IID)
              else
                MassFluxin = MassFluxin + xface_qn(IID)*dtime*xface_Length(IID)*xface_vel(IID)
              endif     
            enddo
          else  !flow on west side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_qn(IID) .gt. 0) then
                MassFluxin = MassFluxin + xface_qn(IID)*dtime*xface_Length(IID)*xface_vel(IID)
              else
                MassFluxout = MassFluxout + xface_qn(IID)*dtime*xface_Length(IID)*xface_vel(IID)
              endif     
            enddo                
          endif
        ENDIF    
      enddo ! end of NQdriver
    endif  !Q_single    

!___________________momentum fluxes at tidal BCs_______________________________  
    do i=1,num_ext_N
      if(yface_q(ext_N(i,2)) .gt. 0) then
        MomFluxout = MomFluxout + yface_q(ext_N(i,2))*(1.0+fac_DW)*yface_vel(ext_N(i,2))*(1.0+fac_DW)*dtime*yface_Length(ext_N(i,2))/4.0
      else
        MomFluxout = MomFluxout + yface_q(ext_N(i,2))*(1.0+fac_UW)*yface_vel(ext_N(i,2))*(1.0+fac_UW)*dtime*yface_Length(ext_N(i,2))/4.0
      endif
    enddo
    do i=1,num_ext_s
      if(yface_q(ext_s(i,2)) .gt. 0) then
        MomFluxin = MomFluxin + yface_q(ext_s(i,2))*(1.0+fac_DW)*yface_vel(ext_s(i,2))*(1.0+fac_DW)*dtime*yface_Length(ext_s(i,2))/4.0
      else
        MomFluxin = MomFluxin + yface_q(ext_s(i,2))*(1.0+fac_UW)*yface_vel(ext_s(i,2))*(1.0+fac_UW)*dtime*yface_Length(ext_s(i,2))/4.0
      endif
    enddo
    do i=1,num_ext_W
      if(xface_q(ext_W(i,2)) .gt. 0) then
        MomFluxin = MomFluxin + xface_q(ext_w(i,2))*(1.0+fac_DW)*xface_vel(ext_w(i,2))*(1.0+fac_DW)*dtime*xface_Length(ext_w(i,2))/4.0
      else        
        MomFluxin = MomFluxin + xface_q(ext_w(i,2))*(1.0+fac_UW)*xface_vel(ext_w(i,2))*(1.0+fac_UW)*dtime*xface_Length(ext_w(i,2))/4.0
      endif
    enddo
    do i=1,num_ext_E
      if(xface_q(ext_E(i,2)) .gt. 0) then
        MomFluxout = MomFluxout + xface_q(ext_e(i,2))*(1.0+fac_DW)*xface_vel(ext_e(i,2))*(1.0+fac_DW)*dtime*xface_Length(ext_e(i,2))/4.0
      else
        MomFluxout = MomFluxout + xface_q(ext_e(i,2))*(1.0+fac_UW)*xface_vel(ext_e(i,2))*(1.0+fac_UW)*dtime*xface_Length(ext_e(i,2))/4.0
      endif
    enddo
    
!___________________momentum fluxes due to friction_______________________________
                     !computed in momentum solvers
!___________________momentum fluxes due to Wind_______________________________
                     !computed in momentum solvers
  
!________________output at prescribe interval_____________________________________
    if(balancecheck_time .ge. BalCheck_int) then
      balancecheck_time = 0
      Error = MassTotal - MassTotal_old + MassFluxout - MassFluxin
      Error_percent = Error*100/(MassTotal + MassFluxin + MassFluxout)
      write(961,"(6(e19.9,','))")timehrs,MassFluxin,MassFluxout,MassTotal,Error,Error_percent
      MassFluxin = 0
      MassFluxout = 0
      MassTotal_old = MassTotal
      MassTotal = 0 
       
      Error = MomTotal - MomTotal_old + MomFluxout - MomFluxin + MomFricLoss - MomWindAdd
      Error_percent = Error*100/(abs(MomTotal) + abs(MomFluxin) + abs(MomFluxout) + abs(MomFricLoss) + abs(MomWindAdd))        
      write(962,"(9(e19.9,','))")timehrs,MomFluxin,MomFluxout,MomTotal,MomFricLoss,MomWindAdd,MomPressAdd,Error,Error_percent         
      MomFluxin = 0
      MomFluxout = 0
      MomTotal_old = MomTotal        
      MomTotal = 0
      MomFricLoss = 0
      MomWindAdd = 0
      MomPressAdd = 0
        
      if(saltrans) then
        Error = SalMassTotal - SalMassTotal_old + SalFluxout - SalFluxin
        Error_percent = Error*100/(SalMassTotal + SalFluxin + SalFluxout)            
        write(963,"(6(e19.9,','))")timehrs,Salfluxin,SalFluxout,SalmassTotal,Error,Error_percent             
        Salfluxin=0
        SalFluxout=0
        SalmassTotal_old = SalmassTotal      
        SalmassTotal=0
      endif
    endif   
    
    return
    end subroutine !balancecheck

!***********************************************************************
    subroutine balancecheck_tel_mass()
!***********************************************************************
    use EXP_Global_def,  only: active
    use EXP_TELESCOPING, only: numregcells, regcells, numtbcells, tbcells
    use flow_def, only: eta
    use geo_def,  only: dx,dy
    use geo_def,  only: zb 
    use BalanceCheck_def, only: masstotal
      
    implicit none
    integer i,j,ii,iid 
      
!________________Water Volume__________________________________________
    MassTotal = 0
    do ii=1,numREGCells
      i=REGCells(ii) 
      if(active(i,3))  MassTotal = MassTotal + (eta(i)-zb(i))*dx(i)*dy(i)  
    enddo       
    do ii=1,numTBCells
      i=TBCells(ii) 
      if(active(i,3))  MassTotal = MassTotal + (eta(i)-zb(i))*dx(i)*dy(i) 
    enddo  
      
    return
    end subroutine  !balancecheck_tel_mass

!***********************************************************************
    subroutine balancecheck_tel_mom()
!***********************************************************************
    use EXP_TELESCOPING, only: numtbxfaces, tbxfaces, xface_wall, xface_cells, xface_qn, xface_length
    use EXP_TELESCOPING, only: numtbyfaces, tbyfaces, yface_wall, yface_cells, yface_qn, yface_length
    use EXP_TELESCOPING, only: numregxfaces, numregyfaces, regxfaces, regyfaces
    use flow_def, only: ikind
    use geo_def,  only: dx,dy
    use geo_def,  only: zb 
    use sed_def,  only: sedtrans
    use sal_def,  only: saltrans,sal       
    use BalanceCheck_def, only: momtotal
      
    implicit none
    integer i,j,ii,iid 
    real*8 DELTAY,DELTAX
      
!________________total momentum__________________________________________
    MomTotal = 0
    do ii=1,numTBXfaces
      i=TBXfaces(ii)
      if(.not. xface_wall(i)) then
        DELTAX = (DX(xface_cells(1,i))+DX(xface_cells(2,i)))/2.0_ikind                
        MomTotal = MomTotal + xface_Qn(i)*xface_length(i)*DELTAX   
      endif
    enddo
    do ii=1,numTBYfaces
      i=TBYfaces(ii)
      if(.not. yface_wall(i)) then
        DELTAY = (DY(yface_cells(1,i))+DY(yface_cells(2,i)))/2.0_ikind               
        MomTotal = MomTotal + yface_Qn(i)*yface_length(i)*DELTAY      
      endif  
    enddo
    do ii=1,numREGXfaces
      i=REGXfaces(ii) 
      if(.not. xface_wall(i)) then         
        DELTAX = (DX(xface_cells(1,i))+DX(xface_cells(2,i)))/2.0_ikind        
        MomTotal = MomTotal + xface_Qn(i)*xface_length(i)*DELTAX
      endif
    enddo
    do ii=1,numREGYfaces
      i=REGYfaces(ii)  
      if(.not. yface_wall(i)) then
        DELTAY = (DY(yface_cells(1,i))+DY(yface_cells(2,i)))/2.0_ikind       
        MomTotal = MomTotal + yface_Qn(i)*yface_length(i)*DELTAY 
      endif
    enddo

    return
    end subroutine  !balancecheck_tel_mom    

!***********************************************************************
    subroutine balancecheck_tel_sal_mass()
!***********************************************************************
    use EXP_Global_def,  only: active 
    use EXP_TELESCOPING, only: numregcells, numtbcells, regcells, tbcells
    use flow_def, only: eta
    use geo_def,  only: dx,dy
    use geo_def,  only: zb 
    use sal_def,  only: saltrans,sal       
    use BalanceCheck_def, only: salmasstotal 
      
    implicit none
    integer i,j,ii,iid 
    
!______________Salinity Mass ___________________________________________________
    if(saltrans) then
      SalMassTotal = 0
      do ii=1,numREGCells
        i=REGCells(ii) 
        if(active(i,3))  SalMassTotal = SalMassTotal + (eta(i)-zb(i))*dx(i)*dy(i)*sal(i) 
      enddo       
      do ii=1,numTBCells
        i=TBCells(ii) 
        if(active(i,3))  SalMassTotal = SalMassTotal + (eta(i)-zb(i))*dx(i)*dy(i)*sal(i)
      enddo    
    endif   

    return
    end subroutine !balancecheck_tel_mass_sal   

!***********************************************************************
    subroutine balancecheck_tel_sal_flux()   
!***********************************************************************
    use EXP_Global_def,    only: num_ext_n, num_ext_s, num_ext_w, num_ext_e
    USE EXP_bndcond_def,   only: qstringexp, ext_n, ext_s, ext_w, ext_e
    USE EXP_transport_def, only: tsalt_elapse
    use EXP_TELESCOPING,   only: xface_flux, yface_flux
    use bnd_def, only: nqstr, q_str
    use geo_def, only: dx,dy
    use geo_def, only: zb 
    use sal_def, only: saltrans,sal       
    use BalanceCheck_def, only: salfluxout, salfluxin
      
    implicit none
    integer i,j,ii,iid 
    
!________________salinity fluxes at Q BC strings___________________________    
    if(nQstr .gt. 0) then    
      do i = 1,nQstr  !for each cell string           
        IF(QstringEXP(i)%vface) THEN          
          if(QstringEXP(I)%sgn .eq. -1) then !flow on north side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_flux(IID) .gt. 0) then
                SalFluxout = SalFluxout + yface_flux(IID)*tsalt_elapse
              else
                SalFluxin = SalFluxin + yface_flux(IID)*tsalt_elapse
              endif
            enddo
          else !flow on south side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_flux(IID) .gt. 0) then
                SalFluxin = SalFluxin + yface_flux(IID)*tsalt_elapse
              else
                SalFluxout = SalFluxout + yface_flux(IID)*tsalt_elapse
              endif
            enddo                  
          endif  
        ELSE            
          if(QstringEXP(I)%sgn .eq. -1) then !flow on east side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_flux(IID) .gt. 0) then
                SalFluxout = SalFluxout + xface_flux(IID)*tsalt_elapse
              else
                SalFluxin = SalFluxin + xface_flux(IID)*tsalt_elapse
              endif     
            enddo
          else  !flow on west side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_flux(IID) .gt. 0) then
                SalFluxin = SalFluxin + xface_flux(IID)*tsalt_elapse
              else
                SalFluxout = SalFluxout + xface_flux(IID)*tsalt_elapse
              endif     
            enddo                
          endif
        ENDIF    
      enddo ! end of NQdriver
    endif  !Q_single    
    
!__________________salinity fluxes at tidal BCs____________________________             
    do i=1,num_ext_N
      if(yface_flux(ext_N(i,2)) .gt. 0) then
        SalFluxout = SalFluxout + yface_flux(ext_N(i,2))*tsalt_elapse
      else
        SalFluxin = SalFluxin - yface_flux(ext_N(i,2))*tsalt_elapse
      endif
    enddo
    do i=1,num_ext_S
      if(yface_flux(ext_S(i,2)) .gt. 0) then
        SalFluxin = SalFluxin + yface_flux(ext_S(i,2))*tsalt_elapse
      else
        SalFluxout = SalFluxout - yface_flux(ext_S(i,2))*tsalt_elapse
      endif
    enddo
    do i=1,num_ext_W
      if(xface_flux(ext_W(i,2)) .gt. 0) then
        SalFluxin = SalFluxin + xface_flux(ext_W(i,2))*tsalt_elapse
      else
        SalFluxout = SalFluxout - xface_flux(ext_W(i,2))*tsalt_elapse
      endif
    enddo
    do i=1,num_ext_E
      if(xface_flux(ext_E(i,2)) .gt. 0) then
        SalFluxout = SalFluxout + xface_flux(ext_E(i,2))*tsalt_elapse
      else
        SalFluxin = SalFluxin - xface_flux(ext_E(i,2))*tsalt_elapse
      endif
    enddo  
    
    return
    end subroutine !balancecheck_tel_sal_flux()      
    
!***********************************************************************
    subroutine balancecheck_tel_sed_mass()
!***********************************************************************
    use EXP_TELESCOPING, only: numregcells, numtbcells, regcells, tbcells
    use geo_def, only: dx,dy
    use geo_def, only: zb 
    use sed_def, only: sedtrans
    use sal_def, only: saltrans
    use BalanceCheck_def, only: salmasstotal
      
    implicit none
    integer i,j,ii,iid 
    
!______________Sediment Mass ___________________________________________________
    if(saltrans) then
      SalMassTotal = 0
      do ii=1,numREGCells
        i=REGCells(ii) 
        ! if(active(i,3))  SMassTotal = SMassTotal + (eta(i)-zb(i))*dx(i)*dy(i)*sed(i) 
      enddo       
      do ii=1,numTBCells
        i=TBCells(ii) 
        ! if(active(i,3))  SMassTotal = SMassTotal + (eta(i)-zb(i))*dx(i)*dy(i)*sed(i)
      enddo    
    endif   

    return
    end subroutine !balancecheck_tel_sed_mass   

!***********************************************************************
   subroutine balancecheck_tel_sed_flux()   
!***********************************************************************
    use EXP_Global_def,    only: num_ext_n, num_ext_s, num_ext_w, num_ext_e
    USE EXP_bndcond_def,   only: qstringexp, ext_n, ext_s, ext_w, ext_e
    USE EXP_transport_def, only: tsalt_elapse
    use EXP_TELESCOPING,   only: xface_flux, yface_flux
    use bnd_def, only: nqstr, q_str
    use geo_def, only: dx,dy
    use geo_def, only: zb 
    use BalanceCheck_def, only: salfluxout, salfluxin
      
    implicit none
    integer i,j,ii,iid 
    
!________________sediment fluxes at Q BC strings___________________________    
    if(nQstr .gt. 0) then    
      do i = 1,nQstr  !for each cell string           
        IF(QstringEXP(i)%vface) THEN          
          if(QstringEXP(I)%sgn .eq. -1) then !flow on north side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_flux(IID) .gt. 0) then
                SalFluxout = SalFluxout + yface_flux(IID)*tsalt_elapse
              else
                SalFluxin = SalFluxin + yface_flux(IID)*tsalt_elapse
              endif
            enddo
          else !flow on south side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j) 
              if(yface_flux(IID) .gt. 0) then
                SalFluxin = SalFluxin + yface_flux(IID)*tsalt_elapse
              else
                SalFluxout = SalFluxout + yface_flux(IID)*tsalt_elapse
              endif
            enddo                  
          endif  
        ELSE            
          if(QstringEXP(I)%sgn .eq. -1) then !flow on east side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_flux(IID) .gt. 0) then
                SalFluxout = SalFluxout + xface_flux(IID)*tsalt_elapse
              else
                SalFluxin = SalFluxin + xface_flux(IID)*tsalt_elapse
              endif     
            enddo
         else  !flow on west side
            do j=1,Q_str(i)%NCells    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              if(xface_flux(IID) .gt. 0) then
                SalFluxin = SalFluxin + xface_flux(IID)*tsalt_elapse
              else
                SalFluxout = SalFluxout + xface_flux(IID)*tsalt_elapse
              endif     
            enddo                
          endif
        ENDIF    
      enddo ! end of NQdriver
    endif  !Q_single    
    
!__________________sediment fluxes at tidal BCs____________________________             
    do i=1,num_ext_N
      if(yface_flux(ext_N(i,2)) .gt. 0) then
        SalFluxout = SalFluxout + yface_flux(ext_N(i,2))*tsalt_elapse
      else
        SalFluxin = SalFluxin - yface_flux(ext_N(i,2))*tsalt_elapse
      endif
    enddo
    do i=1,num_ext_S
      if(yface_flux(ext_S(i,2)) .gt. 0) then
        SalFluxin = SalFluxin + yface_flux(ext_S(i,2))*tsalt_elapse
      else
        SalFluxout = SalFluxout - yface_flux(ext_S(i,2))*tsalt_elapse
      endif
    enddo
    do i=1,num_ext_W
      if(xface_flux(ext_W(i,2)) .gt. 0) then
        SalFluxin = SalFluxin + xface_flux(ext_W(i,2))*tsalt_elapse
      else
        SalFluxout = SalFluxout - xface_flux(ext_W(i,2))*tsalt_elapse
      endif
    enddo
    do i=1,num_ext_E
      if(xface_flux(ext_E(i,2)) .gt. 0) then
        SalFluxout = SalFluxout + xface_flux(ext_E(i,2))*tsalt_elapse
      else
        SalFluxin = SalFluxin - xface_flux(ext_E(i,2))*tsalt_elapse
      endif
    enddo  
    
    return
    end subroutine !balancecheck_tel_sed_flux()          