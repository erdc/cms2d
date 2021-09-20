!================================================================
! CMS Surface Roller routines
!
! Description:
!   The surface roller model solves the 2D steady-state surface
!   roller energy balance equation. The current implimentation
!   uses the same grid as the wave model for convenience. 
!   The governing equation is solves using an iterative
!   finite difference method. The model runs after each
!   spectral wave model run, adds a roller contribution to 
!   the wave radiation stresses which are then passes to the
!   flow model. 
!
! Contains:
!   rol_default - Sets roller model default parameters
!   rol_cards   - Reads the roller cards from the control file
!   rol_init    - Initializes the roller model
!   rol_print   - Prints the roller model settings to the screen
!                 and diagnositic file
!   rol_solve   - Solves the roller model
!   rol_write   - Writes the output for the roller model
!
! written by Alex Sanchez, USACE-CHL
!================================================================
    
!*************************************************************   
    subroutine rol_default
! Sets roller model default parameters
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use rol_def
    implicit none   
    
    roller = .false.  !Turns on the roller model
    rolflux = .false. !Turns on the roller mass flux
    br = 0.1          !Breaking dissipation coefficient
    ceff = 1.0        !Efficiency coefficient
    roltol = 0.001    !Norm tolerance
    irolscheme = 1    !Transport/advection scheme
    rol_courant = 0.5
   
    return
    endsubroutine rol_default
   
!*************************************************************   
    subroutine rol_cards(cardname,foundcard)
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use rol_def
    implicit none
    integer :: ierr
    character(len=37) :: cardname,cdum
    logical :: foundcard
    
    foundcard = .true.
    selectcase (cardname)      
    case('CALC_ROLLER','SURFACE_ROLLER','ROLLER')
      call card_boolean(77,roller,ierr)  
        
    case('ROLLER_DISSIPATION_COEFFICIENT','ROLLER_DISSIPATION_COEF')
      backspace(77)
      read(77,*) cardname,br
          
    case('ROLLER_EFFICIENCY_COEFFICIENT','ROLLER_EFFICIENCY_COEF')
      backspace(77)
      read(77,*) cardname,ceff  
        
    case('ROLLER_NORM_TOLERANCE')
      backspace(77)
      read(77,*) cardname,roltol
       
    case('ROLLER_SCHEME','ROLLER_TRANSPORT_SCHEME','ROLLER_ADVECTION_SCHEME')
      backspace(77)
      read(77,*) cardname,cdum     
      if(cdum(1:3)=='LAX')then
        irolscheme = 2 !Lax-Wendroff scheme
      elseif(cdum(1:7)=='UPWIND2' .or. cdum(1:3)=='SOU')then
        irolscheme = 3 !Second order upwind scheme
      else !if(cdum(1:7)=='UPWIND1')then
        irolscheme = 1 !First order upwind scheme
      endif
      
    case('ROLLER_COURANT','ROLLER_COURANT_NUMBER')  
      backspace(77)
      read(77,*) cardname,rol_courant   
      rol_courant = max(min(rol_courant,1.0),0.01)
      
    case default
      foundcard = .false.
          
    endselect

    return
    endsubroutine rol_cards

!****************************************************************
    subroutine rol_init
! Initialize surface roller model    
! written by Alex Sanchez, USACE-CHL
!****************************************************************
    use flow_def, only: waveflux
    use rol_def
    use diag_lib
    use wave_wavegrid_def, only: nwavei,nwavej,xwave,ywave
    implicit none
    
    !if(roller .or. waveflux)then
      !x and y double distances      
      allocate(d2x(nwavei),d2y(nwavej))    
      d2x(1)=xwave(2)-xwave(1)
      d2x(nwavei)=xwave(nwavei)-xwave(nwavei-1)
      d2x(2:nwavei-1)=xwave(3:nwavei)-xwave(1:nwavei-2)
      d2y(1)=ywave(2)-ywave(1)
      d2y(nwavej)=ywave(nwavej)-ywave(nwavej-1)
      d2y(2:nwavej-1)=ywave(3:nwavej)-ywave(1:nwavej-2)      
    !endif
    
    !Cell center distances
    allocate(delx(nwavei),dely(nwavej))
    delx(1:nwavei-1)=xwave(2:nwavei)-xwave(1:nwavei-1)
    delx(nwavei)=delx(nwavei-1)
    dely(1:nwavej-1)=ywave(2:nwavej)-ywave(1:nwavej-1) 
    dely(nwavej)=dely(nwavej-1)   
    
    if(rolflux .and. .not.roller)then
      call diag_print_warning('Roller mass flux turned on but roller model turned off',&
                              '  Ignoring roller mass flux.')
    endif
    
    if(.not.roller)then
      rolflux = .false.
      return
    endif
    
    allocate(Sr(nwavei,nwavej)) !2x the roller energy per unit area (Sr=2*Er) *********
    Sr=0.0
    
    allocate(rxrs(nwavei,nwavej),ryrs(nwavei,nwavej))
    rxrs=0.0; ryrs=0.0
    
    allocate(roldiss(nwavei,nwavej))
    roldiss=0.0
    
    return
    endsubroutine rol_init
    
!****************************************************************
    subroutine rol_print()
! Print roller settings to screen and diagnostic file   
! written by Alex Sanchez, USACE-CHL
!****************************************************************       
    use rol_def
    use diag_def, only: dgunit,dgfile
    use tool_def, only: vstrlz
    implicit none
    integer :: iunit(2),i
    
341 format(' ',A,T40,A,A)     !Added for vstrlz function results
888 format(' ',A,T40,A)    
222 format(' ',A,T40,F0.2,A)
    
    open(dgunit,file=dgfile,access='append')
    iunit = (/6, dgunit/)
    do i=1,2
        write(iunit(i),*)      
      if(roller)then
         write(iunit(i),888)     'Roller Model:','ON'
        if(rolflux)then
          write(iunit(i),888)    '  Mass Flux:','ON'
        else
          write(iunit(i),888)    '  Mass Flux:','OFF'
        endif
        write(iunit(i),341)      '  Dissipation Coefficient:',trim(vstrlz(br,'(F0.2)'))
        write(iunit(i),341)      '  Wave Breaking Efficiency:',trim(vstrlz(ceff,'(F0.2)'))
        write(iunit(i),341)      '  Courant Number:',trim(vstrlz(rol_courant,'(F0.2)'))
        if(irolscheme==1)then
          write(iunit(i),888)    '  Transport Scheme:','UPWIND1'
        elseif(irolscheme==2)then        
          write(iunit(i),888)    '  Transport Scheme:','LAX'
        else
          write(iunit(i),888)    '  Transport Scheme:','UPWIND2'
        endif
      else
        write(iunit(i),888)      'Roller Model:','OFF'
        endif    
    enddo
    close(dgunit)
    
    return
    endsubroutine rol_print        

!****************************************************************
    subroutine rol_solve
! Simplified roller model    
! written by Alex Sanchez, USACE-CHL
!****************************************************************          
    use rol_def
    use flow_def, only: rhow,hdry,grav
    use diag_lib
    use flow_wavegrid_def, only: depwave,uwave,vwave,etawave,hwave
    use wave_wavegrid_def, only: nwaveij,nwavei,nwavej,xwave,ywave,&
        dxwav,dywav,wheight,wperiod,wdiss,wcos,wsin,wxrs1,wyrs1,wibr
    !use wavestress3D     !Wu, 8/5/2011
    !use fl3d, only: dsigma   !Wu
    use const_def, only: twopi,pi,deg2rad,small
    use prec_def  
    implicit none
    integer :: i,j,k,i1,i2,j1,j2
    integer :: irxy(nwavei,nwavej)
    real(ikind) :: Sxx(nwavei,nwavej),Sxy(nwavei,nwavej),Syy(nwavei,nwavej)
    real(ikind) :: c(nwavei,nwavej),cx(nwavei,nwavej),cy(nwavei,nwavej) !Roller velocities
    real(ikind) :: cySrdy,cxSrdx,fac,hminlim,val    
    real(ikind) :: rn,rn1,dtrol,densitinv
    character(len=100) :: msg
    !real(ikind) :: ratiozroller(numlaywavegrid)   !Wu
    !real(ikind) :: ratioz    !Wu,
    !integer :: klay   !Wu
    
    call diag_print_message(' ','*** Starting Roller ****')

!    hminlim=5.0*hdry   !Depth to limit stresses ************************
    hminlim=max(5.0*hdry,0.05)   !Depth to limit stresses ************************
    
    fac=0.0625*rhow*grav
!$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,nwavej
      do i=1,nwavei
        !Roller speeds  
!        wlen=wavelength(wperiod(i,j),hwave(i,j)) !No wave-current interaction
!        kh=twopi/wlen*hwave(i,j)
!!        wdiss(i,j)=wave_diss(hwave(i,j),wheight(i,j),wperiod(i,j),wlen) 
!        c(i,j)=wlen/wperiod(i,j)
        c(i,j)=sqrt(grav*max(hwave(i,j),hdry)) !long wave approximation
!!        ECg(i,j)=(0.5+kh/sinh(2.0*kh))*c(i,j)*fac*wheight(i,j)*wheight(i,j) !Energy flux
        cx(i,j)=c(i,j)*wcos(i,j)
        cy(i,j)=c(i,j)*wsin(i,j)
        cx(i,j)=cx(i,j)+uwave(i,j)
        cy(i,j)=cy(i,j)+vwave(i,j)
        !Determine wet and dry areas
        !if(hwave>hdry .and. abs(wxrs1)>small .and. abs(wyrs1)>small .and. wheight>1.0e-6)then
        if(hwave(i,j)>hdry .and. wheight(i,j)>1.0e-6)then
          irxy(i,j) = 1
        else
          irxy(i,j) = 0  
        endif
      enddo
    enddo    
!$OMP END PARALLEL DO

    !Time step    
    dtrol = 1.e20
    do i=1,nwavei-1
      do j=1,nwavej   
        val=min(delx(i),dely(j))/max(c(i,j),0.1)
        dtrol=min(dtrol,val)
      enddo
    enddo
    dtrol = rol_courant*dtrol
    
    call diag_print_message('Iteration  Error')
    
661 format(3x,i4,2x,1pe12.3)
    
!    Sr=0.0 !Initialize  
    do k=1,1000      
      rn=0.0
      do i=2,nwavei-1
      do j=1,nwavej    
!!        if(i==37 .and. j==1)then
!!          continue
!!        endif    
        if(irxy(i,j)==0) cycle
        selectcase(irolscheme)
          case(1) !First Order Upwind scheme
            if(j>1 .and. cy(i,j)>0.0)then
              cySrdy=(cy(i,j)*Sr(i,j)-cy(i,j-1)*Sr(i,j-1))/dely(j)      
            elseif(j<nwavej .and. cy(i,j)<=0.0)then
              cySrdy=(cy(i,j+1)*Sr(i,j+1)-cy(i,j)*Sr(i,j))/dely(j+1)
            else
              cySrdy=0.0  
            endif  
            if(cx(i,j)>0.0)then
              cxSrdx=(cx(i,j)*Sr(i,j)-cx(i-1,j)*Sr(i-1,j))/delx(i)  
            elseif(cx(i,j)<0.0)then
              cxSrdx=(cx(i+1,j)*Sr(i+1,j)-cx(i,j)*Sr(i,j))/delx(i+1)
            else
              cxSrdx=0.0
            endif
          case(2) !Lax-Wendroff scheme
            cxSrdx=(cx(i+1,j)*Sr(i+1,j)-cx(i-1,j)*Sr(i-1,j))/d2x(i) &
                  -0.5*dtrol*(cx(i+1,j)*cx(i+1,j)*Sr(i+1,j) &
                         -2.0*cx(i  ,j)*cx(i  ,j)*Sr(i  ,j) &
                             +cx(i-1,j)*cx(i-1,j)*Sr(i-1,j))/delx(i+1)/delx(i) 
            if(j>1 .and. j<nwavej)then                                                       
              cySrdy=(cy(i,j+1)*Sr(i,j+1)-cy(i,j-1)*Sr(i,j-1))/d2y(j) &
                    -0.5*dtrol*(cy(i,j+1)*cy(i,j+1)*Sr(i,j+1) &
                           -2.0*cy(i,j  )*cy(i,j  )*Sr(i,j  ) &
                               +cy(i,j-1)*cy(i,j-1)*Sr(i,j-1))/dely(j+1)/dely(j)
            else
              if(cy(i,j)>0.0 .and. j>1)then
                cySrdy=(cy(i,j)*Sr(i,j)-cy(i,j-1)*Sr(i,j-1))/dely(j)      
              elseif(cy(i,j)<=0.0 .and. j<nwavej)then
                cySrdy=(cy(i,j+1)*Sr(i,j+1)-cy(i,j)*Sr(i,j))/dely(j+1)
              else
                cySrdy=0.0  
              endif 
            endif                            
        case default !Second Order Upwind Scheme          
          if(cy(i,j)>0.0)then
            if(j>2)then
              cySrdy=(3.0*cy(i,j)*Sr(i,j)-4.0*cy(i,j-1)*Sr(i,j-1) &
                      +cy(i,j-2)*Sr(i,j-2))/(2.0*dely(j))
            elseif(j>1)then
              cySrdy=(cy(i,j)*Sr(i,j)-cy(i,j-1)*Sr(i,j-1))/dely(j)     
            else
              cySrdy=0.0
            endif    
          elseif(cy(i,j)<=0.0)then
            if(j<nwavej-1)then
              cySrdy=(-3.0*cy(i,j)*Sr(i,j)+4.0*cy(i,j+1)*Sr(i,j+1) &
                      -cy(i,j+2)*Sr(i,j+2))/(2.0*dely(j))                
            elseif(j<nwavej)then
              cySrdy=(cy(i,j+1)*Sr(i,j+1)-cy(i,j)*Sr(i,j))/dely(j+1)
            else
              cySrdy=0.0  
            endif    
          else
            cySrdy=0.0  
          endif 
          if(cx(i,j)>0.0)then
            if(i>2)then
              cxSrdx=(3.0*cx(i,j)*Sr(i,j)-4.0*cx(i-1,j)*Sr(i-1,j) &
                      +cx(i-2,j)*Sr(i-2,j))/(2.0*delx(i))
            elseif(j>1)then  
              cxSrdx=(cx(i,j)*Sr(i,j)-cx(i-1,j)*Sr(i-1,j))/delx(i)  
            else
              cxSrdx=0.0
            endif    
          elseif(cx(i,j)<=0.0)then
            if(i<nwavei-1)then
              cxSrdx=(-3.0*cx(i,j)*Sr(i,j)+4.0*cx(i+1,j)*Sr(i+1,j) &
                      -cx(i+2,j)*Sr(i+2,j))/(2.0*delx(i))                
            elseif(i<nwavei)then
              cxSrdx=(cx(i+1,j)*Sr(i+1,j)-cx(i,j)*Sr(i,j))/delx(i+1)
            else
              cxSrdx=0.0  
            endif    
          else
            cxSrdx=0.0  
          endif                         
        endselect      
        roldiss(i,j)=grav*br*Sr(i,j)/c(i,j)                     !Stive and De Vriend, Note: Sr=2*Er
        val=dtrol*(-cxSrdx-cySrdy-roldiss(i,j)+ceff*wdiss(i,j)) !Wave dissipation units should be N/m/s=kg/s^3
        Sr(i,j)=Sr(i,j)+val  
        Sr(i,j)=max(Sr(i,j),0.0)        
        rn=rn+val*val
      enddo !j
      enddo !i
      rn=sqrt(rn/real(nwaveij)) !Normalized L-2 Norm            
      if(k==1)then
        rn1=rn
        write(msg,661) k,rn
        call diag_print_message(msg)
        cycle
      endif
      if(k<20) cycle
      if(rn<roltol) exit !Error
      val=abs(rn-rn1)
      if(val<1.e-6) exit !Change    
      val=val/max(rn,roltol)
      if(val<0.01) exit !Relative change  
      if(mod(k,50)==0)then
        write(msg,661) k,rn
        call diag_print_message(msg)
        if(rn>rn1)then
          call diag_print_warning('Roller Model Divergent',&
                                  '  Roller Stresses Ignored',&
                                  '*** Roller Finished ****')
          return
        endif
        rn1=rn      
      endif       
    enddo !k
    write(msg,661) k,rn
    call diag_print_message(msg)
        
!$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,nwavej
      do i=1,nwavei
        !Stresses due to roller  
        Sxx(i,j) = Sr(i,j)*wcos(i,j)*wcos(i,j)
        Sxy(i,j) = Sr(i,j)*wcos(i,j)*wsin(i,j)
        Syy(i,j) = Sr(i,j)*wsin(i,j)*wsin(i,j)
        !Total dissipation with roller dissipation
        !wdiss(i,j) = wdiss(i,j) + roldiss(i,j) !********** IMPORTANT ****************** commented by bdj 2021-01-13 
        !wdiss(i,j) = roldiss(i,j)                        ! 1st fix, replaced by next. by bdj 2021-01-13
        wdiss(i,j) = (1.0-ceff)*wdiss(i,j) + roldiss(i,j) ! new statement by bdj 2021-02-03 as suggested by A Sanchez
        !write(*,*) 'bdj ceff = ',ceff
      enddo
    enddo
!$OMP END PARALLEL DO
    
    !Compute gradients of stresses and update wave radiation stress gradients           
    do i=1,nwavei
      i1=min(i+1,nwavei); i2=max(i-1,1)
      do j=1,nwavej
        j1=min(j+1,nwavej); j2=max(j-1,1)
!        if(.true.)then !Central difference
          rxrs(i,j)=-(Sxx(i1,j)-Sxx(i2,j))/d2x(i)-(Sxy(i,j1)-Sxy(i,j2))/d2y(j) 
          ryrs(i,j)=-(Syy(i,j1)-Syy(i,j2))/d2y(j)-(Sxy(i1,j)-Sxy(i2,j))/d2x(i)         
!        else        !Upwind
!          rxrs(i,j)=0.0
!          ryrs(i,j)=0.0   
!          if(cx(i,j)>0.0)then
!            rxrs(i,j)=rxrs(i,j)-(Sxx(i,j)-Sxx(i2,j))/delx(i)
!            ryrs(i,j)=ryrs(i,j)-(Sxy(i,j)-Sxy(i2,j))/delx(i) 
!          elseif(cx(i,j)<0.0)then
!            rxrs(i,j)=rxrs(i,j)-(Sxx(i1,j)-Sxx(i,j))/delx(i+1)
!            ryrs(i,j)=ryrs(i,j)-(Sxy(i1,j)-Sxy(i,j))/delx(i+1) 
!          endif
!          if(j>1 .and. cy(i,j)>0.0)then
!            rxrs(i,j)=rxrs(i,j)-(Sxy(i,j)-Sxy(i,j2))/dely(j)     
!            ryrs(i,j)=ryrs(i,j)-(Syy(i,j)-Syy(i,j2))/dely(j)
!          elseif(j<nwavej .and. cy(i,j)<=0.0)then
!            rxrs(i,j)=rxrs(i,j)-(Sxy(i,j1)-Sxy(i,j))/dely(j+1)     
!            ryrs(i,j)=ryrs(i,j)-(Syy(i,j1)-Syy(i,j))/dely(j+1)
!          endif 
!        endif  
!        fac=(0.5-0.5*cos(pi*min(hwave(i,j),hminlim)/hminlim))**5.0 !To reduce in very shallow water
!        fac=0.5-0.5*cos(pi*min(hwave(i,j),hminlim)/hminlim) !To reduce in very shallow water
        fac=0.5-0.5*cos(pi*min(max(hwave(i,j)-hminlim,0.0),hminlim)/hminlim);
        rxrs(i,j)=rxrs(i,j)*fac !*irxy(i,j)
        ryrs(i,j)=ryrs(i,j)*fac !*irxy(i,j)
      enddo
    enddo
    
    !Boundaries
!!    rxrs(:,1)=rxrs(:,1)
!!    rxrs(:,nwavej)=rxrs(:,nwavej-1)
!!    ryrs(:,1)=ryrs(:,1)
!!    ryrs(:,nwavej)=ryrs(:,nwavej-1)    
    
!    do i=1,nwavei
!     i1=min(i+1,nwavei); i2=max(i,1)
!      do j=1,nwavej
!        j1=min(j+1,nwavej); j2=max(j-1,1)
!        rxrs(i,j)=-(Sxx(i1,j)-Sxx(i2,j))/delx(i)-(Sxy(i,j1)-Sxy(i,j2))/d2y(j) 
!        ryrs(i,j)=-(Syy(i,j1)-Syy(i,j2))/d2y(j)-(Sxy(i1,j)-Sxy(i2,j))/delx(i)
!      enddo
!    enddo    
    
    densitinv = 1.0/rhow
!$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,nwavej
      do i=1,nwavei    
        rxrs(i,j) = rxrs(i,j)*densitinv
        ryrs(i,j) = ryrs(i,j)*densitinv
        wxrs1(i,j) = wxrs1(i,j) + rxrs(i,j)
        wyrs1(i,j) = wyrs1(i,j) + ryrs(i,j)
      enddo
    enddo
!$OMP END PARALLEL DO

    call rol_write
    
    !do i=1,nwavei    !Wu, 8/5/2011
    !   do j=1,nwavej
    !       ratioz=0.0
    !   do klay=1,numlaywavegrid
    !      ratiozroller(klay)=1.0-tanh((2.0*sigmawave(klay)*hwave(i,j)   &
    !                                   /(wheight(i,j)+0.00001))**4)+0.00000000001
    !      ratioz=ratioz+ratiozroller(klay)*dsigma(klay)
    !   enddo  
    !   ratiozroller(:)=ratiozroller(:)/ratioz
    !   do klay=1,numlaywavegrid
    !      wxrs3D1(i,j,klay)=wxrs3D1(i,j,klay)+rxrs(i,j)/hwave(i,j)*ratiozroller(klay)        
    !      wyrs3D1(i,j,klay)=wyrs3D1(i,j,klay)+ryrs(i,j)/hwave(i,j)*ratiozroller(klay)        
    !   enddo    
    !enddo    
    !enddo    !Wu

    call diag_print_message('*** Roller Finished ****',' ')
    
    return
    endsubroutine rol_solve
       
!*****************************************************************************
    subroutine rol_write
!*****************************************************************************
#include "CMS_cpp.h"
    use wave_wavegrid_def, only: nwavei,nwavej,wxrs1,wyrs1,wheight,wdiss
    use comvarbl, only: reftime
    use cms_def, only: nsteer,dtsteer
    use rol_def, only: Sr,rxrs,ryrs
#ifdef XMDF_IO
    use xmdf
#endif
    use prec_def
    implicit none
    integer :: i,j,k,l,IGMX,JGMX
    real :: rol(nwavei*nwavej),radstr(2*nwavei*nwavej)
    real :: rrs(nwavei*nwavej),dis(nwavei*nwavej) !Output variables, must be single precision
    common /FileNames/ OptsFile, DepFile, CurrFile, EngInFile,    &
                         WaveFile, ObsFile, EngOutFile, NestFile,   &
                         BreakFile, RadsFile, StrucFile, SurgeFile, &
                         MudFile, FricFile, FrflFile, BrflFile,     &
                         SpecFile, WindFile, XMDFFile,SetupFile
    character(len=180) :: WaveFile,ObsFile,EngOutFile,BreakFile,RadsFile
    character(len=180) :: OptsFile,DepFile,CurrFile,EngInFile
    character(len=180) :: NestFile,StrucFile,SurgeFile
    character(len=180) :: MudFile,FricFile,FrflFile,BrflFile
    character(len=180) :: SpecFile,WindFile,XMDFFile,SetupFile
    character(len=20) :: PREFIX
    real(8) :: TIME2
    integer :: ierr,PID,DID,DGID,NIJ
    common /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
    integer :: imod,iprp,island,imd,iprpp,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd
    real :: PAI2,PAI,HPAI,RAD,akap,depmin0
    
    IGMX=nwavei
    JGMX=nwavej
    NIJ=IGMX*JGMX 
      
    K=1; L=1      
    do j=1,jgmx
      do i=1,igmx                            
        rol(K)=Sr(i,j)/2.0 !Bug fix, added 2.0
        dis(K)=wdiss(i,j)
        rrs(K)=sqrt(rxrs(i,j)*rxrs(i,j)+ryrs(i,j)*ryrs(i,j))                 
        radstr(L)=rxrs(i,j)
        radstr(L+1)=ryrs(i,j)        
        K=K+1
        L=L+2
      enddo              
    enddo
    TIME2=float(nsteer-1)*dtsteer/3600.0 
    
#ifdef XMDF_IO
    if(ixmdf==1)then
      call XF_OPEN_FILE(trim(XMDFFile),READWRITE,PID,ierr)
      if(ierr<0)then
        call XF_CREATE_FILE (trim(XMDFFile),READWRITE,PID,ierr)
        write(*,*) ' '
        if(ierr<0)then
          write(*,*) 'ERROR CREATING XMDF FILE ',trim(XMDFFile)
          write(*,*) 'Press any key to continue.'
          read(*,*)
          stop
         else
          write(*,*) 'BINARY FILE CREATED: ',trim(XMDFFile)
        endif
      endif
      
      call XF_OPEN_GROUP(PID,'Dataset',DGID,ierr)
      if(ierr<0)then 
        call XF_CREATE_GENERIC_GROUP(PID,'Dataset',DGID,ierr)
        if(ierr < 0)then
          write(*,*) 'COULD NOT CREATE DATASET - '//'Dataset'
          write(*,*) 'Press any key to continue.'
          read(*,*)
          stop
        endif        
      endif
      
      prefix='Roller_Energy'
      call XF_OPEN_GROUP(DGID,trim(prefix),DID,ierr)
      if(ierr<0)then
        call XF_CREATE_SCALAR_DATASET(DGID,trim(prefix),'none', &
               TS_HOURS,0,DID,ierr)
        if(ierr>0)then
          write(*,*) 'CREATED DATASET - '//trim(prefix)
          else
            write(*,*) 'COULD NOT CREATE DATASET - '//trim(prefix)
            write(*,*) 'Press any key to continue.'
            read(*,*)
          stop
        endif
        call XF_DATASET_REFTIME(DID,reftime,ierr)
        call XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ierr)
      endif
      call XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,rol,ierr)
      call XF_CLOSE_GROUP(DID, ierr)
      
      prefix='Roller_Stress_Mag'
      call XF_OPEN_GROUP (DGID,trim(prefix),DID,ierr)
      if(ierr<0)then
        call XF_CREATE_SCALAR_DATASET(DGID,trim(prefix),'none',  &
               TS_HOURS,0, DID,ierr)
        if(ierr>0)then
          write(*,*) 'CREATED DATASET - '//trim(prefix)
          else
            write(*,*) 'COULD NOT CREATE DATASET - '//trim(prefix)
            write(*,*) 'Press any key to continue.'
            read(*,*)
          stop
        endif
        call XF_DATASET_REFTIME(DID,reftime,ierr)
        call XF_SCALAR_DATA_LOCATION(DID,GRID_LOC_CENTER,ierr)
      endif
      call XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,rrs,ierr)
      call XF_CLOSE_GROUP(DID,ierr)
      
      PREFIX='Dissipation'
      call XF_OPEN_GROUP(DGID,trim(prefix),DID,ierr)
      if(ierr<0)then
        call XF_CREATE_SCALAR_DATASET(DGID,trim(prefix),'none',  &
               TS_HOURS,0, DID,ierr)
        if(ierr>0)then
          write(*,*) 'CREATED DATASET - '//trim(prefix)
          else
            write(*,*) 'COULD NOT CREATE DATASET - '//trim(prefix)
          stop
        endif
        call XF_DATASET_REFTIME(DID,reftime,ierr)
        call XF_SCALAR_DATA_LOCATION(DID,GRID_LOC_CENTER,ierr)
      endif
      call XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,dis,ierr)
      call XF_CLOSE_GROUP(DID,ierr)
      
      prefix='Roller_Stress'
      call XF_OPEN_GROUP(DGID,trim(prefix),DID,ierr)
      if(ierr<0)then
        call XF_CREATE_VECTOR_DATASET(DGID,trim(prefix),'none',  &
               TS_HOURS,0, DID,ierr)
        if(ierr>0)then
          write(*,*) 'CREATED DATASET - '//trim(prefix)
        else
          write(*,*) 'COULD NOT CREATE DATASET - '//trim(prefix)
          write(*,*) 'Press any key to continue.'
          read(*,*)
          stop
        endif
        call XF_DATASET_REFTIME(DID,reftime,ierr)
! removed XF_VECTORS_IN_LOCAL_COORDS (DID,ierr) for proper vector direction
        call XF_VECTOR_2D_DATA_LOCS(DID,GRID_LOC_FACE_I,GRID_LOC_FACE_J,ierr)
      endif
      call XF_WRITE_VECTOR_TIMESTEP(DID,TIME2,NIJ,2,radstr,ierr)
      call XF_CLOSE_GROUP(DID,ierr)
      
      call XF_CLOSE_GROUP(DGID,ierr)
      call XF_CLOSE_FILE(PID,ierr)  
    endif
#endif
      
      return
      endsubroutine rol_write
