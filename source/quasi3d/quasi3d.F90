!================================================================================
! Quasi-3D CMS routines (BETA)
!
! Description:
!   The Quasi-3D module is divided into a variable definitions fortran module 
!   (q3d_def.f90), a library module (q3d_lib.F90), and CMS driver routines
!   (quasi3d.F90). In the quasi-3D method, the vertical structure of the current
!   velocity and sediment concentrations are approximated by assuming analytical
!   expressions, or by solving simple 1DV equations, and then either analytically
!   or numerically integrating the 1DV dispersion expersions. 
!
! Contains:
!   q3d_default - Sets the default values for the Quasi3D module 
!   q3d_cards - Reads the Quasi3D module cards 
!   q3d_init  - Initializes Quasi3D module variables 
!   q3d_print - Prints Quasi3D module variables to the screen and the diagnositic file
!   q3d_flow  - Solves the vertical variation of the current velocity    
!   q3d_velpro - Calculates the vertical current velocity profiles 
!   q3d_sed - Solves the vertical variation of the sediment concentration
!   surface_windwave_stress - Calculates the surface stress including wind and waves
!
! written by Alex Sanchez, USACE-CHL
!================================================================================

!*********************************************************************************
    subroutine q3d_default
!Sets the default values for the Quasi3D module    
!*********************************************************************************    
    use q3d_def
    implicit none
    
    q3d = .false. !On/Off
    q3d_to_flow = .true.
    q3d_to_sed = .true.
    
    !Wave-current interaction term
    wavcurint = .false.
    facwci = 0.5
    
    !Velocity for wave model
    ivelwav = 1 !1-mean, 2-surface, 3-weighted
        
    !Streamwise profile due to bottom friction and wind (outside of surf zone)
    istreampro = 1 !0-None, 1-Log, 2-Power
    
    !Transverse profile (normal) due to helical flow and Coriolis
    itranspro = 0 !0-None, 1-Linear    
    
    !Profile in the surf zone
    isurfpro = 0 !0-None, 1-Quad, 2-Log
    
    !Vertical eddy viscosity
    visvcon = 1.0e-4 !Base value
    cvisvv = 0.067 !Shear
    cvisvhc = 0.0 !Celerity x depth
    verteff = 0.5
    visvslp0 = 0.05  !Mean value
    cvisvslpmu = 8.5e-3 !From Steezel
        
    !Output
    nzpsta = 50 !Stations
    q3d_lay = .false.
    nzplay = 3  !Layers (normally usefull to output bottom, middle, and surface layers)
    
    return
    end subroutine q3d_default

!*********************************************************************************
    subroutine q3d_cards(cardname,foundcard)
! Reads the Quasi3D module cards  
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************    
    use q3d_def
    implicit none
    integer :: i,ierr
    character(len=37) :: cardname,cdum
    logical :: foundcard
    
    foundcard = .true.
    selectcase(cardname)          
    case('QUASI3D_MODE','Q3D_MODE')
      call card_boolean(77,q3d,ierr)
      
    case('QUASI3D_TO_FLOW')  
      call card_boolean(77,q3d_to_flow,ierr)
      
    !--- Wave-current interactio ---------------------------------------
    case('QUASI3D_WAVE_CURRENT','QUASI3D_WAVE_CURRENT_INTERACTION',&
      'WAVE_CUR_DEV_INTER')
      call card_boolean(77,wavcurint,ierr)  
      
    case('WAVE_CUR_DEV_INTER_FACTOR')
      backspace(77)
      read(77,*) cardname,facwci
      if(facwci>1.e-6)then
        wavcurint = .true.  
      else
        wavcurint = .false.  
      endif
      
    !--- Velocity for Wave Model ------------------------------------
    case('WAVE_EFFECTIVE_VELOCITY')
      backspace(77)
      read(77,*) cardname, cdum 
      do i=0,3
        if(cdum==avelwav(i))then
          ivelwav = i
          exit
        endif
      enddo
      !if(ivelwav==0) noptvel = 0 !None
      
    !--- Profiles -----------------------
    case('STREAMWISE_PROFILE')
      backspace(77)
      read(77,*) cardname, cdum 
      do i=0,2
        if(cdum==astreampro(i))then
          istreampro = i
          exit
        endif
      enddo
        
    case('TRANSVERSE_PROFILE')
      backspace(77)
      read(77,*) cardname, cdum 
      do i=0,1
        if(cdum==atranspro(i))then
          itranspro = i
          exit
        endif
      enddo
      
    case('SURFZONE_PROFILE')
      backspace(77)
      read(77,*) cardname, cdum 
      do i=0,2
        if(cdum==asurfpro(i))then
          isurfpro = i
          exit
        endif
      enddo
    
    !--- Vertical Eddy Viscosity -----------------------------------------------  
    case('EDDY_VISCOSITY_VERTICAL_SLOPE','EDDY_VISCOSITY_VERTICAL_SLOPE_BASE')
      backspace(77)
      read(77,*) cardname, visvslp0
      visvslp0 = max(visvslp0,1.0e-6)
      
    case('EDDY_VISCOSITY_VERTICAL_SLOPE_COEF','EDDY_VISCOSITY_VERTICAL_SLOPE_COEFF')
      backspace(77)
      read(77,*) cardname, cvisvslpmu
      cvisvslpmu = max(cvisvslpmu,0.0)
    
    case('EDDY_VISCOSITY_VERTICAL_SHEAR')
      backspace(77)
      read(77,*) cardname, cvisvv
      cvisvv = max(cvisvv,0.0)
        
    case('EDDY_VISCOSITY_VERTICAL_CELERITY')
      backspace(77)
      read(77,*) cardname, cvisvhc
      cvisvhc = max(cvisvhc,0.0)
      
    case('EDDY_VISCOSITY_VERTICAL_BASE','EDDY_VISCOSITY_VERTICAL_CONSTANT','EDDY_VISCOSITY_VERTICAL_BACKGROUND')
      backspace(77)
      read(77,*) cardname, visvcon
      visvcon = max(visvcon,1.0e-6) 
      
    case('EDDY_VISCOSITY_VERTICAL_EFFICIENCY','EDDY_VISCOSITY_VERTICAL_EFF')
      backspace(77)
      read(77,*) cardname, verteff    
      
    !--- Output ----------------------------------------------------
    case('VERTICAL_GRID_NUMBER','VERTICAL_GRID_NUMBER_STATION',&
      'STATION_VERTICAL_GRID_NUMBER','SAVEPOINT_VERTICAL_GRID_NUMBER')
      backspace(77)
      read(77,*) cardname, nzpsta
      nzpsta = min(max(nzpsta,5),100)
      
    case('VERTICAL_GRID_NUMBER_LAYERS','GLOBAL_VERTICAL_GRID_NUMBER',&
      'VERTICAL_LAYERS_NUMBER')
      backspace(77)
      read(77,*) cardname, nzplay
      nzplay = min(max(nzplay,2),10)  
      q3d_lay = .true.
    
    case('QUASI3D_LAYERS')  
      call card_boolean(77,q3d_lay,ierr) 
      
    case default
      foundcard = .false.  
               
    endselect
    
    return
    end subroutine q3d_cards

!*********************************************************************************
    subroutine q3d_init
! Initializes Quasi3D module variables 
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************    
    use size_def, only: ncellsD
    use flow_def, only: viscos
    use q3d_def
    use sed_def, only: sedtrans,nsed
    implicit none
    integer :: i
    
    !--- Flow diserpsion ---------------------------------
    allocate(f3dxx(ncellsD),f3dxy(ncellsD),f3dyy(ncellsD))
    f3dxx = 0.0; f3dxy = 0.0; f3dyy = 0.0
    allocate(f3du(ncellsD),f3dv(ncellsD))
    f3du = 0.0; f3dv = 0.0
    
    if(wavcurint)then
      allocate(udsx(ncellsD),udsy(ncellsD))
      udsx = 0.0; udsy = 0.0
    endif
    
    !--- Sediment dispersion -----------------------------------
    if(sedtrans)then
      allocate(s3dxk(ncellsD,nsed),s3dyk(ncellsD,nsed))
      s3dxk=0.0; s3dyk = 0.0
    endif
            
    !--- Vertical Eddy Viscosity ------------
    !Base (constant) value
    visvcon = max(visvcon,viscos)
    
    !Celerity model
    if(cvisvhc>1.0e-6)then
      cvisvv = 0.0
    endif
    
    !Shear model
    if(cvisvv>1.0e-6)then
      cvisvhc = 0.0 
    endif
    
    !--- Current velocity profiles (temporary variables for output) ---
    !Stations (Save Points)
    allocate(zpsta(nzpsta),uzpsta(nzpsta),vzpsta(nzpsta))
    uzpsta = 0.0; vzpsta = 0.0
    do i=1,nzpsta
      zpsta(i)=real(i,ikind)/real(nzpsta+1,ikind)
    enddo
    
    !Layers
    if(.not.allocated(zplay))then        
      allocate(zplay(nzplay))  
      do i=1,nzplay
        zplay(i)=real(i-1,ikind)/real(nzplay-1,ikind)
        zplay(i)=max(zplay(i),0.02)
        zplay(i)=min(zplay(i),0.98)
        !zplay(i)=real(i,ikind)/real(nzplay+1,ikind)
      enddo
    endif
    allocate(uzplay(nzplay,ncellsD),vzplay(nzplay,ncellsD))
    uzplay = 0.0; vzplay = 0.0
    
    return
    end subroutine q3d_init
    
!*********************************************************************************
    subroutine q3d_print()
! Prints Quasi3D module variables to the screen and the diagnositic file
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************    
    use diag_def, only: dgunit,dgfile
    use q3d_def
    implicit none
    integer :: iunit(2),i
    
    iunit=(/6,dgunit/)
        
111 format(' ',A,T40,A)    
778 format(' ',A,T40,F0.5)
    
    open(dgunit,file=dgfile,access='append') 
    write(*,*)
    
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),111)   'Quasi-3D Mode:','ON (BETA)'
      if(q3d_to_flow)then
        write(iunit(i),111) '  Q3D to 2DH Interaction:','ON'
      else
        write(iunit(i),111) '  Q3D to 2DH Interaction:','OFF'
      endif
      if(wavcurint)then
        write(iunit(i),111) '  Wave-current Interaction:','ON'
        write(iunit(i),778) '    Factor:',facwci
      else
        write(iunit(i),111) '  Wave-current Interaction:','OFF'
      endif
      write(iunit(i),111)   '  Profiles:'
      write(iunit(i),111)   '    Streamwise:',astreampro(istreampro)
      write(iunit(i),111)   '    Transverse:',atranspro(itranspro)
      write(iunit(i),111)   '    Surf zone:',asurfpro(isurfpro)    
      write(iunit(i),111)   '  Vertical Viscosity: '
      write(iunit(i),778)   '    Base Value:',visvcon
      write(iunit(i),778)   '    Surface Stress Factor:',verteff
      write(iunit(i),778)   '    Shear Coefficient:',cvisvv
      write(iunit(i),778)   '    Celerity Coefficient:',cvisvhc
      if(isurfpro==2)then
        write(iunit(i),778) '    Slope Base Value:',visvslp0
        write(iunit(i),778) '    Slope Coefficient:',cvisvslpmu
      endif
    enddo
    
    close(dgunit)
    
    return
    end subroutine q3d_print  
    
!*********************************************************************************
    subroutine q3d_flow
! Solves the vertical variation of the current velocity    
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************    
    use size_def, only: ncells
    use flow_def, only: h,vis,iwet,u,v,uv,us,vs,dux,duy,dvx,dvy,fc,rhow
    use flow_lib, only: streamwise_curvature
    use fric_def, only: z0,bsxy,uelwc,cfrict,z0,bbl_stream
    use fric_lib, only: fric_normapprough,bed_current_wave_stress,&
      fric_streaming_stress,fric_conv_nlength2pow
    use met_def, only: windconst,windvar,wndx,wndy,uwind,vwind
    use met_lib, only: wind_normrough
    use wave_flowgrid_def, only: wavediss,wunitx,wunity,wlen,wper,worbrep,whgt
    use cms_def, only: noptset
    use geo_def, only: mapid
    use q3d_def
    use q3d_lib
    implicit none
    integer :: i
    real(ikind):: ucx,ucy,wx,wy,unx,uny !,udsx,udsy
    real(ikind):: tausx,tausy,taubx,tauby,zwp,zap,ustars,taustr,Kc
    real(ikind):: za,uvc,visvm,visv0,visvslp,ax,ay,bx,by,dx,dy,c,cm,cn
    logical :: isnankind
    
    !Initialize variables
    ucx=0.0; ucy=0.0    !Eulerian mean current velocities [m/s]
    wx=0.0;  wy=0.0     !Lagrangian wind velocities [m/s]
    unx=0.0; uny=0.0    !Secondary flow (normal) current velocities [m/s]
    !udsx=0.0; udsy=0.0  !Current velocities deviations at the surface [m/s]
    
    do i=1,ncells
      if(iwet(i)==0)then
        f3dxx(i)=0.0; f3dxy(i)=0.0; f3dyy(i)=0.0
        cycle
      endif
      call bed_current_wave_stress(bsxy(i),u(i),v(i),us(i),vs(i),ucx,ucy,uvc,taubx,tauby)  
      if(noptset>=3)then !Waves        
        if(bbl_stream)then
          !zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
          !za = zap*h(i)
          !za = 0.00533*max(worbrep(i),1.0e-5)**2.25
          !za = 1.6667e-05 !~d50/12.0
          za = z0(i)
          taustr = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i))
          taubx = taubx - taustr*wunitx(i) !*ramp
          tauby = tauby - taustr*wunity(i) !*ramp
        endif  
        if(wavediss(i)>0.001)then !Surf zone with wind
        !if(.true.)then
          call surface_windwave_stress(i,tausx,tausy,ustars)  
          visvm = q3d_eddyvert_mean(h(i),bsxy(i),ustars)          
          if(isurfpro==1)then !Quadratic
            call q3d_flow_coef_surfquad(h(i),tausx,tausy,&
              taubx,tauby,visvm,ax,ay,bx,by)
            call q3d_flow_disp_surfquad(h(i),ax,ay,bx,by,&
              f3dxx(i),f3dxy(i),f3dyy(i))
            if(wavcurint)then
              call q3d_flow_uds_surfquad(ax,ay,bx,by,udsx(i),udsy(i))
              call q3d_flow_wavcur(h(i),udsx(i),udsy(i),us(i),vs(i),&
                  f3dxx(i),f3dxy(i),f3dyy(i))
              !if(mapid(i)==2390)then
              !  continue
              !endif
            endif
          else !Logarithmic
            call q3d_eddyvert_bottom_slope(h(i),whgt(i),visvm,visv0,visvslp)  
            call q3d_flow_coef_surflog(h(i),tausx,tausy,taubx,tauby,&
                 visv0,visvslp,ax,ay,bx,by,dx,dy,c)
            call q3d_flow_disp_surflog(h(i),ax,ay,bx,by,dx,dy,c,&
                 f3dxx(i),f3dxy(i),f3dyy(i))
            if(wavcurint)then
              call q3d_flow_uds_surflog(ax,ay,bx,by,dx,dy,c,udsx(i),udsy(i))
              call q3d_flow_wavcur(h(i),udsx(i),udsy(i),us(i),vs(i),&
                  f3dxx(i),f3dxy(i),f3dyy(i))
            endif           
          endif
          cycle !Go to next cell
        endif !Surf zone
      endif !waves 
      
      !Not surf zone
      if(windvar)then !Lagrangian wind velocities
        wx = uwind(i)-ucx
        wy = vwind(i)-ucy
      else
        wx = wndx-ucx
        wy = wndy-ucy
      endif  
      zwp = wind_normrough(wx,wy)
      zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
      if(itranspro>0)then
        Kc = streamwise_curvature(ucx,ucy,uv(i),duy(i),duy(i),dvx(i),dvy(i))
        call q3d_flow_normal(h(i),ucx,ucy,fc(i),zap,Kc,unx,uny)            
      endif
      if(istreampro==1)then !Log-log outside of surf zone with wind and waves                    
        call q3d_flow_disp_logloglin(h(i),ucx,ucy,unx,uny,zap,wx,wy,zwp, &
             f3dxx(i),f3dxy(i),f3dyy(i))
        if(wavcurint)then
          call q3d_flow_uds_logloglin(ucx,ucy,unx,uny,zap,wx,wy,zwp,udsx(i),udsy(i))
          call q3d_flow_wavcur(h(i),udsx(i),udsy(i),us(i),vs(i),&
            f3dxx(i),f3dxy(i),f3dyy(i))
        endif    
      elseif(istreampro==2)then !Power-Power outside of surf zone with wind and waves  
        cm = fric_conv_nlength2pow(zap)
        cn = fric_conv_nlength2pow(zwp)  
        call q3d_flow_disp_powpowlin(h(i),ucx,ucy,unx,uny,cm,wx,wy,cn, &
             f3dxx(i),f3dxy(i),f3dyy(i))
        if(wavcurint)then
          call q3d_flow_uds_powpowlin(ucx,ucy,unx,uny,cm,wx,wy,udsx(i),udsy(i))
          call q3d_flow_wavcur(h(i),udsx(i),udsy(i),us(i),vs(i),&
            f3dxx(i),f3dxy(i),f3dyy(i))
        endif
      endif
      if(isnankind(f3dxx(i)) .or. isnankind(f3dxy(i)) .or. isnankind(f3dyy(i)))then
        continue
      endif
    enddo    
    
    return
    end subroutine q3d_flow
    
!****************************************************************************************************
    subroutine q3d_velpro(i,nzp,zp,iprofile,ucx,ucy,unx,uny,wx,wy,zap,zwp,Kc,visvm,visv0,visvslp,uzp,vzp)
! Calculates the vertical current velocity profiles    
!
! Input:
!   i - Cell id [-]
!   nzp - # of vertical coordinates [-]
!   zp - Normalized vertical coordinates [-]
!
! Output:
!   ucx,ucy - Stream-wise current velocity components [m/s]
!   wx,wy - Lagrangian wind velocity components [m/s]
!   unx,uny - Normal-to-stream current velocity components [m/s]
!   zap   - Normalized bed roughness [-]
!   zwp   - Normalized surface roughness [-]
!   Kc    - Stream-wise curvature [1/m]
!   visvm - Mean vertical eddy viscosity [m^2/s] 
!   visv0 - Bottom vertical eddy viscosity [m^2/s]
!   visvslp - Slope of vertical eddy viscosity [m/s]
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************************************
    use size_def, only: ncells
    use flow_def, only: h,vis,iwet,u,v,uv,us,vs,dux,duy,dvx,dvy,fc,rhow
    use flow_lib, only: streamwise_curvature
    use fric_def, only: z0,bsxy,uelwc,cfrict,bbl_stream
    use fric_lib, only: fric_normapprough,bed_current_wave_stress,&
      fric_streaming_stress,fric_conv_nlength2pow
    use met_def, only: windconst,windvar,wndx,wndy,uwind,vwind
    use met_lib, only: wind_normrough
    use wave_flowgrid_def, only: wavediss,wunitx,wunity,worbrep,wlen,wper,whgt
    use cms_def, only: noptset
    use q3d_def, only: isurfpro,itranspro,istreampro
    use q3d_lib
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in):: i,nzp
    real(ikind),intent(in):: zp(nzp)
    !Output
    integer,    intent(out):: iprofile     !0-none, 1-Logloglin, 2-powpowlin, 3-quad, 4-log    
    real(ikind),intent(out):: ucx,ucy,wx,wy,unx,uny,zap,zwp,Kc,visvm,visv0,visvslp
    real(ikind),intent(out):: uzp(nzp),vzp(nzp)
    !Internal    
    real(ikind):: tausx,tausy,taubx,tauby
    real(ikind):: ustars,uvc,taustr,za,ax,ay,bx,by,dx,dy,c,cm,cn
    
    !Early exit for dry cells
    if(iwet(i)==0)then
      uzp=0.0; vzp=0.0
      return
    endif
    
    !Initialize variables
    ucx=0.0; ucy=0.0 !Stream-wise current velocities
    wx=0.0; wy=0.0   !Lagrangian wind velocity components
    unx=0.0; uny=0.0 !Normal-to-stream current velocity components   
    iprofile = 0
    zap = 0.0
    zwp = 0.0
    visvm = 0.0
    visv0 = 0.0
    visvslp = 0.0
    
    call bed_current_wave_stress(bsxy(i),u(i),v(i),us(i),vs(i),ucx,ucy,uvc,taubx,tauby)
    if(noptset>=3)then
      if(bbl_stream)then
        !zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness
        !za = zap*h(i)  
        !za = 0.00533*max(worbrep(i),1.0e-5)**2.25
        !za = 1.6667e-05 !=0.2/1000/12
        za = z0(i)
        taustr = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i))
        taubx = taubx - taustr*wunitx(i)
        tauby = tauby - taustr*wunity(i)
      endif
      if(wavediss(i)>0.001)then !Surf zone with wind
        call surface_windwave_stress(i,tausx,tausy,ustars)
        visvm = q3d_eddyvert_mean(h(i),bsxy(i),ustars)
        !!val = rhow*6.0*visvm/h(i)
        !!tausx = tausx + (u(i)-us(i))*val-2.0*taubx
        !!tausy = tausy + (v(i)-vs(i))*val-2.0*tauby
        if(isurfpro==1)then !Quadratic
          call q3d_flow_coef_surfquad(h(i),tausx,tausy,taubx,tauby,visvm,ax,ay,bx,by)
          call q3d_velpro_surfquad(nzp,zp,ucx,ucy,ax,ay,bx,by,uzp,vzp)
          visv0 = visvm
          visvslp = 0.0
          iprofile = 3
        else !Logarithminc
          call q3d_eddyvert_bottom_slope(h(i),whgt(i),visvm,visv0,visvslp)
          call q3d_flow_coef_surflog(h(i),tausx,tausy,taubx,tauby,&
            visv0,visvslp,ax,ay,bx,by,dx,dy,c)
          call q3d_velpro_surflog(nzp,zp,ucx,ucy,ax,ay,bx,by,dx,dy,c,uzp,vzp)
          iprofile = 4
        endif
        return
      endif !In surf zone
    endif !Waves
    
    !Not surf zone
    visv0 = 0.0
    visvslp = 0.0
    if(windvar)then !Lagrangian wind velocities
      wx = uwind(i)-ucx
      wy = vwind(i)-ucy
    else
      wx = wndx-ucx
      wy = wndy-ucy
    endif  
    zwp = wind_normrough(wx,wy)  
    zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness
    if(itranspro>0)then
      Kc = streamwise_curvature(ucx,ucy,uv(i),duy(i),duy(i),dvx(i),dvy(i))
      call q3d_flow_normal(h(i),ucx,ucy,fc(i),zap,Kc,unx,uny)
    endif
    if(istreampro==1)then !Log-log outside of surf zone with wind and waves        
      call q3d_velpro_logloglin(nzp,zp,ucx,ucy,unx,uny,zap,wx,wy,zwp,uzp,vzp)
      iprofile = 1
    elseif(istreampro==2)then !Power-Power outside of surf zone with wind and waves   
      cm = fric_conv_nlength2pow(zap)
      cn = fric_conv_nlength2pow(zwp)
      call q3d_velpro_powpowlin(nzp,zp,ucx,ucy,unx,uny,cm,wx,wy,cn,uzp,vzp)
      iprofile = 2
    endif

    return
    end subroutine q3d_velpro

!*********************************************************************************
    subroutine q3d_sed
! Solves the vertical variation of the sediment concentration
!*********************************************************************************    
    use size_def, only: ncells
    use flow_def, only: h,vis,iwet,u,v,uv,us,vs,dux,duy,dvx,dvy
    use fric_def, only: z0,bsxy,uelwc,cfrict
    use fric_lib, only: bed_current_wave_stress
    use met_def, only: windconst,windvar,wndx,wndy,uwind,vwind
    use met_lib, only: wind_normrough
    use wave_flowgrid_def, only: wavediss,wunitx,wunity
    use cms_def, only: noptset
    use sed_def, only: nsed,wsfall,cak,epsvk
    use q3d_def
    use q3d_lib
    implicit none
    integer :: i,ks
    real(ikind):: ucx,ucy        !Stream-wise current velocity components, m/s
    real(ikind):: wx,wy        !Lagrangian wind velocity components, m/s
    real(ikind):: unx,uny        !Normal-to-stream current velocity components, m/s
    real(ikind):: tausx,tausy  !Surface stresses, N/m^2
    real(ikind):: taubx,tauby  !Bed stresses, N/m^2
    real(ikind):: ustars       !Surface friction velocity [m/s]
    real(ikind):: uvc,visvm
    real(ikind):: ax,ay,bx,by
    
    !Initialize variables
    ucx=0.0; ucy=0.0
    wx=0.0; wy=0.0
    unx=0.0; uny=0.0
    
    do i=1,ncells
      if(iwet(i)==0)then
        s3dxk(i,:)=0.0; s3dyk(i,:)=0.0
        cycle
      endif
      call bed_current_wave_stress(bsxy(i),u(i),v(i),us(i),vs(i),ucx,ucy,uvc,taubx,tauby)
!      if(wavediss(i)>0.001)then !Surf zone with wind        
        call surface_windwave_stress(i,tausx,tausy,ustars)
        visvm = q3d_eddyvert_mean(h(i),bsxy(i),ustars)
!       if(isurfpro==1)then
        call q3d_flow_coef_surfquad(h(i),tausx,tausy,taubx,tauby,visvm,ax,ay,bx,by)
        do ks=1,nsed
          call q3d_sed_surfquad(h(i),ax,ay,bx,by,&
             cak(i,ks),wsfall(ks),epsvk(i,ks),s3dxk(i,ks),s3dyk(i,ks))
        enddo !ks
!       else
!         call s3d_surflog(h(i),tausx,tausy,taubx,tauby,visv,us(i),vs(i),s3dxk(i,:),s3dyk(i,:))
!       endif
!      else !Not surf zone
!        if(windvar)then
!          wx = u(i)-uwind(i); wy = v(i)-vwind(i) !Lagrangian velocities
!        else
!          wx = u(i)-wndx; wy = v(i)-wndy !Lagrangian velocities  
!        endif   
!        zwp = wind_normrough(wx,wy)
!        zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
!        if(itranspro>0)then
!          Kc = streamwise_curvature(ucx,ucy,uv(i),duy(i),duy(i),dvx(i),dvy(i))
!          call q3d_flow_normal(h(i),ucx,ucy,fc(i),zap,Kc,unx,uny)            
!        endif
!        if(istreampro==1)then !Log-log outside of surf zone with wind and waves                    
!          call s3d_fricwindlog_helcorlin_wave(h(i),u(i),v(i),unx,uny,zap,wx,wy,zwp, &
!               us(i),vs(i),s3dxk(i,:),s3dyk(i,:))
!        elseif(istreampro==2)then !Power-Power outside of surf zone with wind and waves              
!          call s3d_fricwindpow_helcorlin_wave(h(i),u(i),v(i),unx,uny,zap,wx,wy,zwp, &
!               us(i),vs(i),s3dxk(i,:),s3dyk(i,:))
!        endif  
!      endif
    enddo
    
    return
    end subroutine q3d_sed

!******************************************************************
    subroutine surface_windwave_stress(i,tausx,tausy,ustars)
! Calculates the surface stress including wind and waves
!******************************************************************  
    use geo_def, only: areap
    use flow_def, only: h,u,v,us,vs,rhow,dpx,dpy,grav
    use q3d_def, only: verteff
    use const_def, only: twopi
    use wave_flowgrid_def, only: wavediss,wunitx,wunity,wlen,Ssr,&
       wavestrx,wavestry,Whgt
    use met_def, only: iwndlagr,windconst,windvar,cdWndareap,&
        wndx,wndy,uwind,vwind,tauwindx,tauwindy,tauwx,tauwy
    use rol_def, only: br,ceff
    use cms_def, only: noptset
    use prec_def
    implicit none
    !Input
    integer,intent(in) :: i
    !Output
    real(ikind):: tausx,tausy,ustars
    !Internal
    real(ikind):: val,fac
    
    !Wave stresses
    if(noptset>=3)then
      !fac = 1.0
      fac = 0.5*Whgt(i)/sqrt(2.0)/h(i)
      !val = verteff*wavediss(i)/sqrt(grav*h(i))
      !val = verteff*grav*br*Ssr(i)/sqrt(grav*h(i)) !verteff*Dr/c    
      val = verteff*(ceff*grav*br*Ssr(i)+(1.0-ceff)*wavediss(i))/sqrt(grav*h(i)) !verteff*(fe*Dr+(1-fe)*Db)/c, more general works with/without roller
      !val = verteff*wavediss(i)*twopi/wlen(i)   
      tausx = val*wunitx(i)*fac
      tausy = val*wunity(i)*fac
      !tausx = verteff*rhow*wavestrx(i)*fac
      !tausy = verteff*rhow*wavestry(i)*fac    
    endif
    !tausx = 0.0
    !tausy = 0.0
    !tausx = tausx - rhow*h(i)*dpx(i)*fac
    !tausy = tausy - rhow*h(i)*dpy(i)*fac
    
    !Wind surface stresses
    if(iwndlagr==0)then
      if(windconst)then
        tausx = tausx + tauwx
        tausy = tausy + tauwy
      elseif(windvar)then  
        tausx = tausx + tauwindx(i)
        tausy = tausy + tauwindy(i)
      endif
    elseif(iwndlagr==1)then
      if(windconst)then 
        !Note: better to use near-surface current velocity in the future
        tausx = tausx + cdWndareap(i)*(wndx-u(i)+us(i))/areap(i)
        tausy = tausy + cdWndareap(i)*(wndy-v(i)+vs(i))/areap(i)
      elseif(windvar)then
        tausx = tausx + cdWndareap(i)*(uwind(i)-u(i)+us(i))/areap(i)
        tausy = tausy + cdWndareap(i)*(vwind(i)-v(i)+vs(i))/areap(i)
      endif
    endif
    
    val = sqrt(tausx*tausx+tausy*tausy) !Magnitude of surface stress
    ustars = sqrt(val/rhow) !Surface friction velocity
    
    return
    end subroutine surface_windwave_stress
    
!******************************************************************    
    subroutine q3d_flow_vel_surface(ueff,veff)
!******************************************************************           
    use size_def, only: ncells,ncellsD
    use flow_def, only: h,vis,iwet,u,v,uv,us,vs,dux,duy,dvx,dvy,fc,rhow
    use flow_lib, only: streamwise_curvature
    use fric_def, only: z0,bsxy,uelwc,cfrict,z0,bbl_stream
    use fric_lib, only: fric_normapprough,bed_current_wave_stress,&
      fric_streaming_stress,fric_conv_nlength2pow
    use met_def, only: windconst,windvar,wndx,wndy,uwind,vwind
    use met_lib, only: wind_normrough
    use wave_flowgrid_def, only: wavediss,wunitx,wunity,wlen,wper,worbrep,whgt
    use cms_def, only: noptset
    use geo_def, only: mapid
    use q3d_def
    use q3d_lib
    implicit none
    !Output
    real(ikind),intent(out):: ueff(ncellsD),veff(ncellsD)
    !Internal Variables
    integer :: i
    real(ikind):: ucx,ucy,wx,wy,unx,uny !,udsx,udsy
    real(ikind):: tausx,tausy,taubx,tauby,zwp,zap,ustars,taustr,Kc
    real(ikind):: za,uvc,visvm,visv0,visvslp,ax,ay,bx,by,dx,dy,c,cm
    
    !Initialize variables
    ucx=0.0; ucy=0.0    !Eulerian mean current velocities [m/s]
    wx=0.0;  wy=0.0     !Lagrangian wind velocities [m/s]
    unx=0.0; uny=0.0    !Secondary flow (normal) current velocities [m/s]
    !udsx=0.0; udsy=0.0  !Current velocities deviations at the surface [m/s]
    
    do i=1,ncells
      if(iwet(i)==0)then
        ueff(i)=0.0; veff(i)=0.0
        cycle
      endif
      call bed_current_wave_stress(bsxy(i),u(i),v(i),us(i),vs(i),ucx,ucy,uvc,taubx,tauby)
      if(noptset>=3)then
        if(bbl_stream)then
          !zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
          !za = zap*h(i)
          !za = 0.00533*max(worbrep(i),1.0e-5)**2.25
          !za = 1.6667e-05 !~d50/12.0
          za = z0(i)
          taustr = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i))
          taubx = taubx - taustr*wunitx(i) !*ramp
          tauby = tauby - taustr*wunity(i) !*ramp
        endif
        !if(wavediss(i)>0.001)then !Surf zone with wind
        if(.true.)then
          call surface_windwave_stress(i,tausx,tausy,ustars)
          visvm = q3d_eddyvert_mean(h(i),bsxy(i),ustars)          
          if(isurfpro==1)then !Quadratic profile
            call q3d_flow_coef_surfquad(h(i),tausx,tausy,&
              taubx,tauby,visvm,ax,ay,bx,by)
            call q3d_flow_uds_surfquad(ax,ay,bx,by,udsx(i),udsy(i))          
          else !Logarithmic profile
           call q3d_eddyvert_bottom_slope(h(i),whgt(i),visvm,visv0,visvslp)  
           call q3d_flow_coef_surflog(h(i),tausx,tausy,taubx,tauby,&
             visv0,visvslp,ax,ay,bx,by,dx,dy,c)
           call q3d_flow_uds_surflog(ax,ay,bx,by,dx,dy,c,udsx(i),udsy(i))         
          endif
        endif !Surf zone
        ueff(i) = ucx + udsx(i)
        veff(i) = ucy + udsy(i)
        cycle
      endif !Waves 
      
      !Not surf zone
      if(windvar)then !Lagrangian wind velocities
        wx = uwind(i)-ucx
        wy = vwind(i)-ucy
      else
        wx = wndx-ucx
        wy = wndy-ucy
      endif  
      zwp = wind_normrough(wx,wy)
      zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
      if(itranspro>0)then
        Kc = streamwise_curvature(ucx,ucy,uv(i),duy(i),duy(i),dvx(i),dvy(i))
        call q3d_flow_normal(h(i),ucx,ucy,fc(i),zap,Kc,unx,uny)            
      endif
      if(istreampro==1)then !Log-log outside of surf zone with wind and waves                    
        call q3d_flow_uds_logloglin(ucx,ucy,unx,uny,zap,wx,wy,zwp,udsx(i),udsy(i))
      elseif(istreampro==2)then !Power-Power outside of surf zone with wind and waves  
        cm = fric_conv_nlength2pow(zap)
        call q3d_flow_uds_powpowlin(ucx,ucy,unx,uny,cm,wx,wy,udsx(i),udsy(i))
      endif
      ueff(i) = ucx + udsx(i)
      veff(i) = ucy + udsy(i)
    enddo
    
    return
    end subroutine q3d_flow_vel_surface
    
!******************************************************************    
    subroutine q3d_flow_vel_weighted(ueff,veff)
!******************************************************************           
    use size_def, only: ncells,ncellsD
    use flow_def, only: h,vis,iwet,u,v,uv,us,vs,dux,duy,dvx,dvy,fc,rhow
    use flow_lib, only: streamwise_curvature
    use fric_def, only: z0,bsxy,uelwc,cfrict,z0,bbl_stream
    use fric_lib, only: fric_normapprough,bed_current_wave_stress,&
      fric_streaming_stress,fric_conv_nlength2pow
    use met_def, only: windconst,windvar,wndx,wndy,uwind,vwind
    use met_lib, only: wind_normrough
    use wave_flowgrid_def, only: wavediss,wunitx,wunity,wlen,wper,worbrep,whgt
    use cms_def, only: noptset
    use geo_def, only: mapid
    use const_def, only: pi
    use q3d_def
    use q3d_lib
    implicit none
    !Output
    real(ikind),intent(out):: ueff(ncellsD),veff(ncellsD)
    !Internal Variables
    integer :: i
    real(ikind):: ucx,ucy,wx,wy,unx,uny !,udsx,udsy
    real(ikind):: tausx,tausy,taubx,tauby,zwp,zap,ustars,taustr,Kc,kw
    real(ikind):: za,uvc,visvm,visv0,visvslp,ax,ay,bx,by,dx,dy,c,cm,udex,udey
    
    !Initialize variables
    ucx=0.0; ucy=0.0    !Eulerian mean current velocities [m/s]
    wx=0.0;  wy=0.0     !Lagrangian wind velocities [m/s]
    unx=0.0; uny=0.0    !Secondary flow (normal) current velocities [m/s]
    udex=0.0; udsy=0.0  !Current velocities deviations [m/s]
    
    do i=1,ncells
      if(iwet(i)==0)then
        ueff(i)=0.0; veff(i)=0.0
        cycle
      endif
      kw = 2.0*pi/wlen(i)
      call bed_current_wave_stress(bsxy(i),u(i),v(i),us(i),vs(i),ucx,ucy,uvc,taubx,tauby)
      if(bbl_stream)then
        !zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
        !za = zap*h(i)
        !za = 0.00533*max(worbrep(i),1.0e-5)**2.25
        !za = 1.6667e-05 !~d50/12.0
        za = z0(i)
        taustr = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i))
        taubx = taubx - taustr*wunitx(i) !*ramp
        tauby = tauby - taustr*wunity(i) !*ramp
      endif
      !if(wavediss(i)>0.001)then !Surf zone with wind
      if(.true.)then
        visvm = q3d_eddyvert_mean(h(i),bsxy(i),ustars)
        call surface_windwave_stress(i,tausx,tausy,ustars)
        if(isurfpro==1)then
          call q3d_flow_coef_surfquad(h(i),tausx,tausy,&
            taubx,tauby,visvm,ax,ay,bx,by)          
          call q3d_flow_ude_surfquad(h(i),ax,ay,bx,by,kw,udex,udey)          
        else
         call q3d_eddyvert_bottom_slope(h(i),whgt(i),visvm,visv0,visvslp)  
         call q3d_flow_coef_surflog(h(i),tausx,tausy,taubx,tauby,&
           visv0,visvslp,ax,ay,bx,by,dx,dy,c)
         call q3d_flow_ude_surflog(h(i),ax,ay,bx,by,dx,dy,c,kw,udex,udey)
        endif
      else !Not surf zone
        if(windvar)then !Lagrangian wind velocities
          wx = uwind(i)-ucx
          wy = vwind(i)-ucy
        else
          wx = wndx-ucx
          wy = wndy-ucy
        endif  
        zwp = wind_normrough(wx,wy)
        zap = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness  
        if(itranspro>0)then
          Kc = streamwise_curvature(ucx,ucy,uv(i),duy(i),duy(i),dvx(i),dvy(i))
          call q3d_flow_normal(h(i),ucx,ucy,fc(i),zap,Kc,unx,uny)            
        endif
        if(istreampro==1)then !Log-log outside of surf zone with wind and waves                    
          !call q3d_flow_ude_logloglin(ucx,ucy,unx,uny,zap,wx,wy,zwp,udex,udsy)
        elseif(istreampro==2)then !Power-Power outside of surf zone with wind and waves  
          cm = fric_conv_nlength2pow(zap)
          !call q3d_flow_ude_powpowlin(ucx,ucy,unx,uny,cm,wx,wy,udex,udey)
        endif  
      endif !Outside of surfzone
      ueff(i) = ucx + udex
      veff(i) = ucy + udey
    enddo
    
    return
    end subroutine q3d_flow_vel_weighted 
    
!********************************************************************************   
    subroutine q3d_vel_lay
! Writes an observayion cell profile to its corresponding file
! written by Alex Sanchez, USACE-ERDC-CHL
!********************************************************************************   
    use size_def, only: ncells
    use q3d_def, only: nzplay,zplay,uzplay,vzplay
    use prec_def
    implicit none
    !Internal variables
    integer :: i,iprofile
    real(ikind) :: uc,vc,vx,vy,wx,wy,zap,zwp
    real(ikind) :: Kc,visvm,visv0,visvslp
    
    do i=1,ncells
      call q3d_velpro(i,nzplay,zplay,iprofile,uc,vc,vx,vy,wx,wy,zap,zwp,&
        Kc,visvm,visv0,visvslp,uzplay(:,i),vzplay(:,i)) !Compute velocity vertical profile
    enddo
    
    return
    end subroutine q3d_vel_lay  