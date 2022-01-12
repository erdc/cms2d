!****************************************************
    subroutine steering_default
! Sets the default values for the steering variables
! Author: Alex Sanchez, USACE-CHL
!****************************************************
    use cms_def
    use flow_def
    use comvarbl, only: flowpath
    implicit none
    
    noptset = 3          !Option identifying the mode of simulation
    noptwse = 3          !0-none, 1-eta2=eta1, 2-eta2=tide2, 3-eta2=tide2+eta1-tide1
    nopttime = 0         !0-Interval-based, 1-Time-based
    noptvel = 1          !0-none, 1-vel2=vel1
    noptzb  = 1          !0-none, 1-zb2=zb,2-bed change,3-interpolate from datasets
    nbrksm = 3           !Breaking smoothing iterations
    ndissm = 1           !Dissipation smoothing iterations   !Lihwa changed from 3 to 1    12/18/2020
    nradsm = 1           !Radiation stresses smoothing iterations
    npersm = 4           !Wave period smoothing iterations
    nsteer = 0           !Wave steering interval number
    noptxtrpfl = 2       !0-No extrapolation, 1-User specified, 2-Automatic    
    xtrpdistfl = 0.0     !Flow extrapolation distance
    xtrpdistwav = 0.0    !Wave extrapolation distance    
    noptxtrpwav = 2      !0-No extrapolation, 1-User specified, 2-Automatic   
    dtsteer = -1.0       !hours, steering interval
    cmsflow = .false.    !toggle for running CMS-Flow
    cmswave = .false.    !toggle for running CMS-Wave
    wave_interp = .true. !toggle for temporal wave interpolation   MEB  01/22/2014
    flowpath = ''
    wavepath = ''
    wavedisstol = 0.001 !Tolerance for wave breaking
    waveradfac = 0.0
    
    return
    endsubroutine steering_default
    
!**************************************************************
    subroutine steering_cards(cardname)    
! Reads the CMS Cards related to steering/coupling
! between CMS-Flow and CMS-Wave
!
! History:
!  Version 2, 11/21/2012, Alex Sanchez, USACE-CHL
!    Change - Moved wave smoothing variables from flow cards
!**************************************************************
    use cms_def
    use geo_def, only: grdfile,wgrdfile
    use flow_def
    implicit none
    integer :: ierr
    character(len=10) :: aext
    character(len=37) :: cardname
    character(len=100) :: cdum
    character(len=200) :: astr
    
    selectcase(cardname)                
      case('CMS-WAVE_SIM_FILE','CMS_WAVE_SIM_FILE',&
           'CMSWAVE_SIM_FILE','WAVE_SIM_FILE')
        backspace(77)
        read(77,*) cardname,astr
        call fileparts(astr,Wavepath,wavename,aext)           
        WavSimFile = trim(wavename) // '.sim'     
        cmswave = .true.   
        noptset = 3

      case('STEERING_INTERVAL','CMS-STEERING_INTERVAL')
        call card_scalar(77,'hrs','sec',dtsteer,ierr)
        cmswave = .true.
        
      !Temporal Wave Interpolation toggle - MEB 01/22/2014
      case('TEMPORAL_WAVE_INTERPOLATION')   
        call card_boolean(77,wave_interp,ierr)
      
      case('WAVE_WATER_LEVEL','FLOW-TO-WAVE_WATER_LEVEL','FLOW_TO_WAVE_WATER_LEVEL',&
        'FLOW-TO-WAVE_WATER_ELEVATION','FLOW_TO_WAVE_WATER_ELEVATION')
        backspace(77)   
        read(77,*) cardname, cdum
        if(cdum(1:1)=='1' .or. cdum(1:1)=='N' .or. cdum(1:1)=='n')then
          noptwse = 0   !NONE
        elseif(cdum(1:1)=='1' .or. cdum(1:1)=='L')then
          noptwse = 1   !LAST
        elseif(cdum(1:1)=='3' .or. cdum(6:6)=='_')then
          noptwse = 3   !TIDAL_PLUS_VARIATION   
        else !if(cdum(1:1)=='T')then
          noptwse = 2   !TIDAL 
        endif
      
      case('WAVE_CURRENT_VELOCITY','FLOW-TO-WAVE_CURRENT_VELOCITY','FLOW_TO_WAVE_CURRENT_VELOCITY')
        backspace(77)   
        read(77,*) cardname, cdum
        if(cdum(1:1)=='1' .or. cdum(1:1)=='N' .or. cdum(1:1)=='n')then
          noptvel = 0   !NONE
        else !if(cdum(1:1)=='1' .or. cdum(1:1)=='L' .or. cdum(1:1)=='l')then
          noptvel = 1   !LAST
        endif
        
      case('WAVE_BED_ELEVATION','FLOW-TO-WAVE_BED_ELEVATION','FLOW_TO_WAVE_BED_ELEVATION',&
          'WAVE_WATER_DEPTH','FLOW-TO-WAVE_WATER_DEPTH','FLOW_TO_WAVE_WATER_DEPTH','FLOW_TO_WAVE_DEPTH')
        backspace(77)   
        read(77,*) cardname, cdum
        if(cdum(1:1)=='0' .or. cdum(1:1)=='N' .or. cdum(1:1)=='n')then
          noptzb = 0   !NONE
        elseif(cdum(1:1)=='1' .or. cdum(1:1)=='L' .or. cdum(1:1)=='l')then            
          noptzb = 1   !LAST
        elseif(cdum(1:1)=='2' .or. cdum(1:1)=='C' .or. cdum(1:1)=='c')then            
          noptzb = 2   !CHANGE
        elseif(cdum(1:1)=='3' .or. cdum(1:1)=='I' .or. cdum(1:1)=='i')then   
          noptzb = 3   !INTERP  
        endif

      case('FLOW_EXTRAPOLATION_DISTANCE','FLOW_TO_WAVE_EXTRAPOLATION_DISTANCE','FLOW-TO-WAVE_EXTRAPOLATION_DISTANCE')
        call card_scalar(77,'m','m',xtrpdistfl,ierr)
        if (xtrpdistfl < 0) then
          noptxtrpfl = 2 !0-No extrapolation, 1-User specified, 2-Automatic  
        elseif (xtrpdistfl == 0) then
          noptxtrpfl = 0 
        else
          noptxtrpfl = 1
        endif
        
      case('WAVE_EXTRAPOLATION_DISTANCE','WAVE_TO_FLOW_EXTRAPOLATION_DISTANCE','WAVE-TO-FLOW_EXTRAPOLATION_DISTANCE')
        call card_scalar(77,'m','m',xtrpdistwav,ierr) 
        if (xtrpdistwav < 0) then
          noptxtrpwav = 2 !0-No extrapolation, 1-User specified, 2-Automatic  
        elseif (xtrpdistwav == 0) then
          noptxtrpwav = 0 
        else
          noptxtrpwav = 1
        endif
        
      case('WAVE_RSTRESS_DATASET')
        backspace(77)   
        read(77,*) cardname,wgrdfile,radpath    
        noptset = 4
        cmswave = .true.  
        
      case('WAVE_HEIGHT_DATASET')
        backspace(77)   
        read(77,*) cardname,wgrdfile,wavpath    
        noptset = 4  
        cmswave = .true.  
        
      case('WAVE_PERIOD_DATASET')
        backspace(77)   
        read(77,*) cardname,wgrdfile,perpath    
        noptset = 4    
        cmswave = .true.  
        
      case('WAVE_DIRECTION_DATASET')
        backspace(77)   
        read(77,*) cardname,wgrdfile,dirpath    
        noptset = 4    
        cmswave = .true.    
        
      case('WAVE_DISS_DATASET')
        backspace(77)   
        read(77,*) cardname,wgrdfile,disspath    
        noptset = 4    
        cmswave = .true.   
        
      !==== Wave Smoothing =========================
      !case('WAVE_BREAKING_SMOOTHING_ITERATIONS',&
      !  'WAVE_BREAKING_SMOOTHING_ITER','WAVE_BREAKING_SMOOTH_ITER')
      !  backspace(77)  
      !  read(77,*) cardname, nbrksm
      !  nbrksm = min(max(nbrksm,0),10)
          
      case('WAVE_DISSIPATION_SMOOTHING_ITERATIONS',&
        'WAVE_DISSIPATION_SMOOTHING_ITER','WAVE_DISSIPATION_SMOOTH_ITER')
        backspace(77)  
        read(77,*) cardname, ndissm
        ndissm = min(max(ndissm,0),10)  
          
      !case('WAVE_STRESSES_SMOOTHING_ITERATIONS','WAVE_STRESS_SMOOTHING_ITERATIONS',&
      !  'WAVE_STRESSES_SMOOTHING_ITER','WAVE_STRESSES_SMOOTH_ITER','WAVE_STRESS_SMOOTH_ITER')
      !  backspace(77)  
      !  read(77,*) cardname, nradsm
      !  nradsm = min(max(nradsm,0),10)
          
      case('WAVE_PERIOD_SMOOTHING_ITERATIONS',&
        'WAVE_PERIOD_SMOOTHING_ITER','WAVE_PERIOD_SMOOTH_ITER')
        backspace(77)  
        read(77,*) cardname, npersm
        npersm = min(max(npersm,0),10)    
        
    endselect    

    return
    endsubroutine steering_cards

!*******************************************************************
    subroutine steer_init()
! Initializes the steering variables
!*******************************************************************    
    use size_def
    use geo_def, only: xOrigin,yOrigin,azimuth_fl
    use wave_flowgrid_def, only: wavfl_intpcoef_file
    use flow_wavegrid_def, only: flwav_intpcoef_file,depwave,depwave0,hwave
    use wave_wavegrid_def, only: nwavei,nwavej,nwaveij,&
       dxwav,dywav,xwave,ywave,azimuth_wav,xwav0,ywav0
    !use global !CMS-Wave module
    use global_inline    !Lihwa recommends these instead      12/18/2020
    use sed_def, only: sedtrans
    use comvarbl, only: flowpath
    !use cms_def, only: wavepath,noptzb,nbrksm,ndissm,nradsm,npersm
    use cms_def, only: wavepath,noptzb,ndissm,npersm
    use const_def, only: deg2rad
    use geo_def, only: bathydata
    use prec_def
    implicit none
    integer :: i,j,ni,nj
    integer :: isteer,iidate,imod,iprp,island,imd,iprpp
    integer :: nonln,igrav,isolv,ixmdf,iproc,imud,iwnd
    real*4 :: depin,etain,uin,vin,x0,y0,azimuth,sinaz,cosaz
    real*4 :: PAI2,PAI,HPAI,RAD,akap,dvarxx,dvaryy,depmin0
    common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz
    common/wavegrid/ni,nj
    common /dxxdyy/dvarxx(ipmx),dvaryy(ipmx) 
    common /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
    common /fl2wav/depin(ipmx,jpmx),etain(ipmx,jpmx), &        !Alex
             uin(ipmx,jpmx),vin(ipmx,jpmx)                       !Alex
    
    !Initialize interpolation coefficient files     
    wavfl_intpcoef_file = trim(wavepath) // 'Intpcoef_wavfl.bin'          
    flwav_intpcoef_file = trim(flowpath) // 'Intpcoef_flwav.bin'    

    !Wave variables on the flow grid
    call wave_flgrid_init

    !nwavei and nwavej are the mesh dismension before rotation in CMS-Wave
    nwavei = ni
    nwavej = nj
    nwaveij = nwavei*nwavej

    !Wave variables on the wave grid (used in getwave and roller)
    call wave_wavegrid_init           
    
    !Flow variables on the wave grid
    call flow_wavegrid_init   
    
    do i=1,nwavei
      do j=1,nwavej
        depwave(i,j)=depin(i,j)
        hwave(i,j)=depwave(i,j)
      enddo
    enddo
    
    !Flow to wave still water depths
    if(noptzb<0) noptzb = 1 !LAST
    !if(noptzb<0)then
      !if(bathydata%ison)then
      !  noptzb = 3 !INTERP
      !else
      !  noptzb = 1 !LAST
      !endif
    !endif
    if(.not.sedtrans .and. .not.bathydata%ison) noptzb = 0 !NONE; If sed or bathy not on, depths are not passed
    
    if(noptzb>0)then
      allocate(depwave0(nwavei,nwavej))
      do i=1,nwavei
        do j=1,nwavej
          depwave0(i,j)=depin(i,j)
        enddo
      enddo
    endif

    do i=1,nwavei
      dxwav(i)=dvarxx(i)
    enddo
    do j=1,nwavej
      dywav(j)=dvaryy(j)
    enddo
    
    xwave(1)=0.5*dxwav(1)
    do i=2,nwavei
      xwave(i)=xwave(i-1)+(dxwav(i-1)+dxwav(i))/2.0
    enddo        
       ywave(1)=0.5*dywav(1)
    do j=2,nwavej
      ywave(j)=ywave(j-1)+(dywav(j-1)+dywav(j))/2.0    
    enddo

    azimuth_wav = azimuth
    xwav0 = x0
    ywav0 = y0

!    cosflow=cos(azimuth_fl*deg2rad)
!    sinflow=sin(azimuth_fl*deg2rad)
!    coswave=cos(azimuth_wav*deg2rad)
!    sinwave=sin(azimuth_wav*deg2rad)
!    
!    !Get flow coordinates with respect to wave cooridinate system
!    do i=1,ncellsD
!      x_global=xOrigin+x(i)*cosflow-y(i)*sinflow
!      y_global=yOrigin+x(i)*sinflow+y(i)*cosflow
!      xflwav(i)= (x_global-xwav0)*coswave+(y_global-ywav0)*sinwave
!      yflwav(i)=-(x_global-xwav0)*sinwave+(y_global-ywav0)*coswave      
!    enddo
!    
!    !Get wave coordinates with respect to flow coordinate system
!    do i=1,nwavei
!      do j=1,nwavej
!        x_global=xwav0+xwave(i)*coswave-ywave(j)*sinwave
!        y_global=ywav0+xwave(i)*sinwave+ywave(j)*coswave
!        xwavfl(i,j)= (x_global-xOrigin)*cosflow+(y_global-yOrigin)*sinflow
!        ywavfl(i,j)=-(x_global-xOrigin)*sinflow+(y_global-yOrigin)*cosflow     
!      enddo
!    enddo
    
    return
    endsubroutine steer_init

!*******************************************************************
    subroutine wave_flgrid_init
! Allocates and initializes the wave variables on the flow grid
! written by Weiming Wu, NCCHE
! modified by Alex Sanchez, USACE-CHL
!*******************************************************************    
    use size_def
    use geo_def, only: zb
    use flow_def, only: hdry
    use const_def, only: deg2rad
    use wave_flowgrid_def
    use rol_def, only: roller
    use prec_def
    implicit none
    integer :: ii
        
    !test for allocation, do not reallocate if already allocated
    if(.not.allocated(Whgt)) then
      !allocate(xflwav(ncellsD),yflwav(ncellsD))                        !Flow grid on wave coordinate system
      allocate(Whgt(ncellsD),Whgt1(ncellsD),Whgt2(ncellsD))             !Significant wave height [m]
      allocate(Wper(ncellsD),Wper1(ncellsD),Wper2(ncellsD))             !Peak wave period [s]
      allocate(Wang(ncellsD))                                           !Wave angle [rad]
      !allocate(waveibr(ncellsD),waveibr1(ncellsD),waveibr2(ncellsD))    !Wave breaking index [-]
      allocate(wavediss(ncellsD),wavediss1(ncellsD),wavediss2(ncellsD)) !Wave breaking dissipation                 
      allocate(Wlen(ncellsD))                                           !Wave length [m]
      allocate(Worb(ncellsD))                                           !Wave bottom orbital velocity based on Whgt and Wper [m/s]
      allocate(Worbrep(ncellsD))                                        !Representative bottom orbital velocity [m/s]
      allocate(wavestrx(ncellsD),wavestry(ncellsD))                     !Wave forcing
      allocate(wavestrx1(ncellsD),wavestry1(ncellsD))
      allocate(wavestrx2(ncellsD),wavestry2(ncellsD))
      allocate(ijwavcell(2,ncellsD),coefintp_wavfl(4,ncellsD))          !Interpolation indeces and coefficients
      allocate(Wunitx(ncellsD),Wunitx1(ncellsD),Wunitx2(ncellsD))       !Wave unit vectors [-]
      allocate(Wunity(ncellsD),Wunity1(ncellsD),Wunity2(ncellsD))       !Wave unit vectors [-]
      allocate(ueff(ncellsD),veff(ncellsD))

      !Moved from below  MEB
      allocate(wetsteer(ncellsD))  
      allocate(Ssr(ncellsD))       !Surface roller energy
    endif

    !Initialize variables, Alex, Sep. 1, 2009, needed for plotting purposes
    Wang  = 0.0; !Wang1 = 0.0; Wang2 = 0.0
    Wlen = 1.0; Worb = 0.0; Worbrep = 0.0
    Whgt  = 0.0; Whgt1 = 0.0; Whgt2 = 0.0; 
    Wper  = 1.0; Wper1 = 1.0; Wper2 = 1.0;
    wavestrx = 0.0; wavestrx1 = 0.0; wavestrx2 = 0.0;
    wavestry = 0.0; wavestry1 = 0.0; wavestry2 = 0.0;
    wavediss = 0.0; wavediss1 = 0.0; wavediss2 = 0.0; 
    !waveibr = 0.0; waveibr1 = 0.0; waveibr2 = 0.0; 
    ijwavcell = 0; coefintp_wavfl = 0.0
    Wunitx = 0.0; Wunitx1 = 0.0; Wunitx2 = 0.0
    Wunity = 0.0; Wunity1 = 0.0; Wunity2 = 0.0         
    ueff = 0.0; veff = 0.0
    
    !Moving allocate(wetsteer) into group above  MEB  01/12/2022
    wetsteer = 0.0
    do ii=1,ncells
      if(-zb(ii)>hdry)then !wet
        wetsteer = 1.0
      endif
    enddo  
          
    !Moving allocate(Ssr) into group above  MEB  01/12/2022
    Ssr = 0.0
    
    return
    endsubroutine wave_flgrid_init
    
!***********************************************************************
    subroutine flow_wavegrid_init
!***********************************************************************    
    use size_def, only: nmaxcells
    use flow_wavegrid_def
    use wave_wavegrid_def, only: nwavei,nwavej
    implicit none
    
    !Grid
    !allocate(xwavfl(nwavei,nwavej),ywavfl(nwavei,nwavej)) !Wave grid on flow coordinate system          
    allocate(depwave(nwavei,nwavej)) !Alex, added depwave
    
    !Flow    
    allocate(hwave(nwavei,nwavej),etawave(nwavei,nwavej)) !Alex, added hwave June 22, 2012
    allocate(uwave(nwavei,nwavej),vwave(nwavei,nwavej))
    etawave=0.0; uwave=0.0; vwave=0.0     
            
    !Interpolation
    allocate(iiflcell(0:nmaxcells,nwavei,nwavej))
    allocate(coefintp_flwav(nmaxcells,nwavei,nwavej))  
    iiflcell = 0; coefintp_flwav = 0.0
    
    return
    endsubroutine flow_wavegrid_init

!***********************************************************************
    subroutine wave_wavegrid_init
!***********************************************************************    
    use wave_wavegrid_def
    use prec_def
    implicit none
    
    allocate(xwave(nwavei),ywave(nwavej)) !Wave grid on the wave grid coordinate system    
    allocate(dxwav(nwavei),dywav(nwavej)) !Alex
    allocate(wxrs1(nwavei,nwavej),wyrs1(nwavei,nwavej))
    allocate(wheight(nwavei,nwavej),wperiod(nwavei,nwavej))
    !allocate(wibr(nwavei,nwavej),wdiss(nwavei,nwavej))
    allocate(wdiss(nwavei,nwavej))							 !Lihwa 12/18/2020
    allocate(wcos(nwavei,nwavej),wsin(nwavei,nwavej))  
    xwave=0.0; ywave=0.0
    wxrs1=0.0; wyrs1=0.0
    wheight=0.0; wperiod=0.0
    !wibr=0.0; wdiss=0.0
    wdiss=0.0												 !Lihwa 12/18/2020
    wcos=0.0; wsin=0.0
    
    return
    endsubroutine wave_wavegrid_init
    
!***********************************************************************
    subroutine interp_coef_flwav()
! Calculates interpolation coefficients from flow to waves
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces,nmaxcells,&
        nnodes,ncellsimple,ncellpoly
    use geo_def, only: cell2cell,idirface,ncface,x,y,dx,dy,&
        nflowgrd,xOrigin,yOrigin,azimuth_fl,areaavg,nncell,node2cell
    use comvarbl, only: flowpath
    use cms_def, only: noptxtrpfl,xtrpdistfl,wavepath
    use flow_wavegrid_def
    use wave_flowgrid_def
    use wave_wavegrid_def, only: nwavei,nwavej,xwav0,ywav0,&
          azimuth_wav,xwave,ywave
    use bnd_def, only: nbndstr,bnd_str
    use interp_lib, only: interp_coef_tel2cart,interp_coef_poly2cart
    use const_def, only: deg2rad
    use diag_lib
    use prec_def
    implicit none
    integer :: i,j,k,nk,iversionchk
    integer, parameter :: iversion = 7 !flow-to-wave interpolation file version number
    integer :: nflowgrdchk,nchk,nDchk,nichk,njchk,ierr,nmchk
    real(ikind) :: xtrpdistflchk
    character(len=200) :: msg
    logical :: foundfile
    !real(ikind) :: coswave,sinwave,xglobal,yglobal
    !character(len=200) :: aname
    
    !Calculate extrapolation distance for flow
    if(noptxtrpfl==2)then
        xtrpdistfl = 2.0*sqrt(areaavg) !Extrapolation distance
!        xtrpdistfl=0.0
!        do i=1,nwavei
!          xtrpdistfl=max(xtrpdistfl,dxwav(i))
!        enddo  
!        do j=1,nwavej
!          xtrpdistfl=max(xtrpdistfl,dywav(j))
!        enddo
!        xtrpdistfl=2.0*xtrpdistfl
    endif 
    
    inquire(file=flwav_intpcoef_file,exist=foundfile)
    if(foundfile)then
      open(92,file=flwav_intpcoef_file,form='unformatted')
      read(92,iostat=ierr) iversionchk
      if(iversionchk==iversion .and. ierr==0)then
        read(92,iostat=ierr) nflowgrdchk
        if(nflowgrdchk==nflowgrd .and. ierr==0)then
          read(92,iostat=ierr) nchk,nDchk,nichk,njchk
          if(nchk==ncells .and. nDchk==ncellsD .and. &
            nichk==nwavei .and. njchk==nwavej .and. ierr==0)then
            read(92,iostat=ierr) xtrpdistflchk
            if(abs(xtrpdistflchk-xtrpdistfl)<1.0e-3 .and. ierr==0)then
              read(92,iostat=ierr) nmchk  
              if(abs(nmchk-nmaxcells)==0 .and. ierr==0)then
                call diag_print_message(' ','Reading flow-to-wave interpolation coefficients',' ')
         loop1: do i=1,nwavei
                  do j=1,nwavej        
                    read(92,iostat=ierr) iiflcell(0,i,j)
                    if(ierr/=0) exit loop1
                    nk = max(iiflcell(0,i,j),1)
                    read(92,iostat=ierr) (iiflcell(k,i,j),k=1,nk)
                    if(ierr/=0) exit loop1
                    read(92,iostat=ierr) (coefintp_flwav(k,i,j),k=1,nk)
                    if(ierr/=0) exit loop1
                  enddo
                enddo loop1
                if(ierr==0)then
                  close(92) 
                  return !Read file successfully
                endif !ierr
              endif !nmchk
            endif !xtrpdistflchk
          endif !nchk,nDchk
        endif !nflowgrdchk
      endif !iversionchk
      close(92,status='delete')
    endif !foundfile

    write(msg,'(A,F9.3,A)') '   Extrapolation distance: ',xtrpdistfl,' m'
    call diag_print_message(' ','Calculating flow-to-wave interpolation coefficients',msg)
    
    if(ncellsimple>0)then
      call interp_coef_tel2cart(ncells,ncellsD,nmaxfaces,&
        xOrigin,yOrigin,azimuth_fl,x,y,dx,dy,cell2cell,idirface,ncface, &            
        nwavei,nwavej,xwav0,ywav0,azimuth_wav,xwave,ywave,&
        xtrpdistfl,nbndstr,bnd_str,iiflcell,coefintp_flwav)
    elseif(ncellpoly>0)then
      call interp_coef_poly2cart(ncells,ncellsD,nmaxfaces,x,y,&
        ncface,cell2cell,nnodes,nmaxcells,nncell,node2cell, &
        nwavei,nwavej,xwav0,ywav0,azimuth_wav,xwave,ywave,&
        xtrpdistfl,iiflcell,coefintp_flwav)
    endif
    
    call diag_print_message(' ','Saving flow-to-wave interpolation coefficients',' ')
    open(92,file=flwav_intpcoef_file,form='unformatted')  
    write(92) iversion
    write(92) nflowgrd
    write(92) ncells,ncellsD,nwavei,nwavej
    write(92) xtrpdistfl
    write(92) nmaxcells
    do i=1,nwavei
      do j=1,nwavej 
        write(92) iiflcell(0,i,j)
        nk = max(iiflcell(0,i,j),1)
        write(92) (iiflcell(k,i,j),k=1,nk)
        write(92) (coefintp_flwav(k,i,j),k=1,nk)
      enddo
    enddo
    close(92)
    
    !!ASCII File
    !flwav_intpcoef_file = trim(flowpath) // 'Intpcoef_flwav.dat'  
    !open(92,file=flwav_intpcoef_file)  
    !write(92,*) iversion
    !write(92,*) nflowgrd
    !write(92,*) ncells,ncellsD,nwavei,nwavej
    !write(92,*) xtrpdistfl
    !write(92,*) nmaxcells
    !do i=1,nwavei
    !  do j=1,nwavej        
    !    write(92,*) (iiflcell(k,i,j),k=0,nmaxcells)
    !    write(92,*) (coefintp_flwav(k,i,j),k=1,nmaxcells)
    !  enddo
    !enddo
    !close(92)
    
    !return
    !coswave=cos(azimuth_wav*deg2rad)
    !sinwave=sin(azimuth_wav*deg2rad)
    !!Output wave grid
    !aname = trim(wavepath) // 'Coord_Wave.xyz'
    !open(4234,file=aname)
    !do i=1,nwavei
    !  do j=1,nwavej
    !    xglobal=xwav0+xwave(i)*coswave-ywave(j)*sinwave
    !    yglobal=ywav0+xwave(i)*sinwave+ywave(j)*coswave  
    !    write(4234,*) xglobal,yglobal,0.0 !,depwave(i,j)
    !  enddo
    !enddo
    !close(4234)
    
    !aname = trim(flowpath) // 'Coord_Flow.xyz'
    !open(4234,file=aname)
    !do i=1,ncells
    !  val=-zb(i)
    !  write(4234,*) x(i),y(i),val
    !enddo
    !close(4234)

    return 
    endsubroutine interp_coef_flwav
    
!***********************************************************************
    subroutine interp_coef_wavfl()
!   Calculates interpolation coefficients from wave to flow grids
!***********************************************************************
    use size_def, only: ncells,ncellsD,ncellsimple,ncellpoly
    use geo_def, only: nflowgrd,xOrigin,yOrigin,azimuth_fl,x,y
    use cms_def, only: noptxtrpwav,xtrpdistwav
    use wave_flowgrid_def, only: ijwavcell,coefintp_wavfl,&
                           wavfl_intpcoef_file
    use wave_wavegrid_def, only: nwavei,nwavej,dxwav,dywav,xwave,ywave,&
        xwav0,ywav0,azimuth_wav
    use comvarbl, only: flowpath
    use interp_lib, only: interp_coef_cart2tel,interp_coef_cart2pts
    use prec_def
    use diag_lib
    use bnd_def
    implicit none
    integer :: i,j,ierr
    integer :: nchk,nDchk,nichk,njchk,nflowgrdchk,iversionchk
    integer, parameter :: iversion = 6
    real(ikind) :: xtrpdistwavchk
    character(len=200) :: msg
    logical :: foundfile
    
    !Calculate extrapolation distance for waves
    if(noptxtrpwav==2)then
!      xtrpdistwav = 10.0*sqrt(areaavg) !Extrapolation distance      
      xtrpdistwav=0.0
!      do ii=1,ncells
!        xtrpdistwav=max(xtrapdist,dx(ii),dy(ii))
!      enddo
!      xtrpdistwav=10.0*xtrapdist    
      do i=1,nwavei      
        xtrpdistwav = max(xtrpdistwav,dxwav(i))
      enddo  
      do j=1,nwavej
        xtrpdistwav = max(xtrpdistwav,dywav(j))
      enddo 
      xtrpdistwav = 2.0*xtrpdistwav
    endif  
    
    inquire(file=wavfl_intpcoef_file,exist=foundfile)
    if(foundfile)then
      open(91,file=wavfl_intpcoef_file,form='unformatted')
      read(91,iostat=ierr) iversionchk
      if(iversionchk==iversion .and. ierr==0)then
        read(91,iostat=ierr) nflowgrdchk
        if(nflowgrdchk==nflowgrd .and. ierr==0)then
          read(91,iostat=ierr) nchk,nDchk,nichk,njchk
          if(nchk==ncells .and. nDchk==ncellsD .and. &
            nichk==nwavei .and. njchk==nwavej .and. ierr==0)then
            read(91,iostat=ierr) xtrpdistwavchk
            if(abs(xtrpdistwavchk-xtrpdistwav)<1.0e-4 .and. ierr==0)then   
              call diag_print_message(' ','Reading wave-to-flow interpolation coefficients',' ')
              do i=1,ncellsD
                read(91,iostat=ierr) (ijwavcell(j,i),j=1,2)
                if(ierr/=0) exit
                read(91,iostat=ierr) (coefintp_wavfl(j,i),j=1,4)
                if(ierr/=0) exit
              enddo
              close(91)
              return !File read successfully
            endif
          endif  
        endif     
      endif           
      close(91,status='delete')
    endif
    
    write(msg,'(A,F9.3,A)')'   Extrapolation distance: ',xtrpdistwav,' m'
    call diag_print_message('Calculating wave-to-flow interpolation coefficients',msg)
    
    !Note: All points are interpolated to Cartesian grid
    if(ncellsimple>0)then
      call interp_coef_cart2tel(nwavei,nwavej,xwav0,ywav0,&
          azimuth_wav,xwave,ywave,dxwav,dywav,&
          ncellsD,ncellsD,xOrigin,yOrigin,azimuth_fl,x,y,&
          xtrpdistwav,ijwavcell,coefintp_wavfl)
    elseif(ncellpoly>0)then
      call interp_coef_cart2pts(nwavei,nwavej,xwav0,ywav0,&
           azimuth_wav,xwave,ywave,dxwav,dywav,&
           ncellsD,ncellsD,x,y,xtrpdistwav,ijwavcell,coefintp_wavfl)
    endif
    
    call diag_print_message(' ','Saving wave-to-flow interpolation coefficients',' ')
    open(91,file=wavfl_intpcoef_file,form='unformatted')
    write(91) iversion
    write(91) nflowgrd !Grid modification number
    write(91) ncells,ncellsD,nwavei,nwavej !Grid dimensions
    write(91) xtrpdistwav !Extrapolation distance
    do i=1,ncellsD
      write(91) (ijwavcell(j,i),j=1,2)
      write(91) (coefintp_wavfl(j,i),j=1,4)
    enddo
    close(91)
    
    !!!ASCII file
    !!wavfl_intpcoef_file = trim(flowpath) // 'Intpcoef_wavfl.dat'
    !!open(91,file=wavfl_intpcoef_file)
    !!write(91,*) iversion
    !!write(91,*) nflowgrd !Grid modification number
    !!write(91,*) ncells,ncellsD,nwavei,nwavej !Grid dimensions
    !!write(91,*) xtrpdistwav !Extrapolation distance
    !!do i=1,ncellsD
    !!  write(91,*) (ijwavcell(i,j),j=1,2)
    !!  write(91,*) (coefintp_wavfl(i,j),j=1,4)
    !!enddo
    !!close(91) 
              
    return 
    endsubroutine interp_coef_wavfl    

!***********************************************************************
    subroutine interp_scal_flwav(varflow,varwave,iextrap)
! Interpolate from flow mesh to wave mesh 
!
!
! Authors: Mingliang, Z., Alex Sanchez, Weiming, Wu
!***********************************************************************
    use size_def, only: ncellsD
    use flow_def, only: iwet
    use wave_wavegrid_def, only: nwavei,nwavej
    use flow_wavegrid_def, only: iiflcell,coefintp_flwav
    use interp_lib, only: interp_scal_tel2cart
    use prec_def
    implicit none
    !Input
    integer,    intent(in) :: iextrap
    real(ikind),intent(in) :: varflow(ncellsD)
    !Output
    real(ikind),intent(out) :: varwave(nwavei,nwavej)
    !Internal
    real(ikind) :: valdry
            
    valdry = -999.0
    
    call interp_scal_tel2cart(ncellsD,varflow,iwet, &
           nwavei,nwavej,iiflcell,coefintp_flwav,&
           varwave,valdry,iextrap)

    return 
    endsubroutine interp_scal_flwav   
    
!***********************************************************************
    subroutine interp_vec_flwav(vecxflow,vecyflow,&
           vecxwave,vecywave,iextrap)
! Interpolate from flow mesh to wave mesh 
!
!
! Authors: Mingliang, Z., Alex Sanchez, Weiming, Wu
!***********************************************************************
    use size_def, only: ncellsD
    use flow_def, only: iwet
    use wave_wavegrid_def, only: nwavei,nwavej
    use flow_wavegrid_def, only: iiflcell,coefintp_flwav
    use interp_lib, only: interp_vec_tel2cart
    use prec_def
    implicit none
    !Input
    integer,    intent(in) :: iextrap
    real(ikind),intent(in),dimension(ncellsD) :: vecxflow,vecyflow
    !Output
    real(ikind),intent(out),dimension(nwavei,nwavej) :: vecxwave,vecywave
    !Internal
    real(ikind) :: valdry
            
    valdry = -999.0
    
    call interp_vec_tel2cart(ncellsD,vecxflow,vecyflow,iwet, &
           nwavei,nwavej,iiflcell,coefintp_flwav,&
           vecxwave,vecywave,valdry,iextrap)

    return 
    endsubroutine interp_vec_flwav    

!***********************************************************************
    subroutine setwave()
! Calculates the hydro variables on the wave grid
!
!Input:
! u - u-velocity on flow grid and time step
! v - v-velocity on flow grid and time step
! zb - bed elevation (negative is wet)
! wetsteer - 1 if cell is wet during whole steering internal and 0 otherwise
!
!Output:
! uwave - u-velocity on wave grid and time step
! vwave - v-velocity on wave grid and time step
! depwave - still water depth on wave grid and time step
! etawave - mean water level on wave grid and time step
! hwave = depwave + etawave
!
! written by Alex Sanchez, USACE-CHL    
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: zb,azimuth_fl,zb0
    use flow_def, only: u,v,p,hdry,hdry,us,vs,grav,gravinv
    use comvarbl, only: ctime
    use cms_def, only: noptwse,noptvel,noptzb,nsteer,dtsteer
    use flow_wavegrid_def
    use wave_flowgrid_def
    use wave_wavegrid_def, only: nwavei,nwavej
    use sed_def, only: sedtrans
    use global_inline, only: IGPX,JGPX,IPMX,JPMX,KOMX,NPF   !Lihwa recommends  12/18/2020
    use const_def, only: deg2rad
#ifdef DEV_MODE
    use q3d_def, only: ivelwav
#endif    
    use diag_lib
    use prec_def
    implicit none
    integer :: ii,i,j,idate
    real(ikind) :: theta,coswave,sinwave,uwaveij,vwaveij
    real(ikind) :: wetwave(nwavei,nwavej)
    real(ikind) :: varflow(ncellsD)
    real*4 :: x0,y0,azimuth,isteer,iidate,sinaz,cosaz
    real*4 :: depin,etain,uin,vin
    common /fl2wav/depin(ipmx,jpmx),etain(ipmx,jpmx), &        !Alex
           uin(ipmx,jpmx),vin(ipmx,jpmx)                       !Alex             
    common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz
    common /FileNames/ OptsFile, DepFile, CurrFile, EngInFile,    &
                       WaveFile, ObsFile, EngOutFile, NestFile,   &
                       BreakFile, RadsFile, StrucFile, SurgeFile, &
                       MudFile, FricFile, FrflFile, BrflFile,     &
                       WindFile, XMDFFile, SetupFile
    character(len=180)  OptsFile, DepFile, CurrFile, EngInFile,        &
                   WaveFile,ObsFile,EngOutFile,NestFile,          &
                   BreakFile,RadsFile, StrucFile, SurgeFile,      & 
                   MudFile, FricFile, FrflFile, BrflFile,         &
                   WindFile, XMDFFile, SetupFile
    character(len=200) :: msg
    
    write(*,*)
    write(*,'(A)') ' Starting flow-to-wave interpolation'                     

    IDate = nsteer
    
    !--- Wetting and drying --------------
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)    
      do j=1,nwavej
        do i=1,nwavei
          wetwave(i,j) = 1.0 !By default cells are wet
        enddo           
      enddo  
!$OMP END PARALLEL DO 
    call interp_scal_flwav(wetsteer,wetwave,0)  !NOT extrapolated. Values outside domain are left unchanged
    
    !--- Still Water Depths ----------------------------
    !(bathymetry only updated in wave model for wet cells)
    !If bathymetry if being passed from flow to wave model
    selectcase(noptzb)
    case(1) !Last bed elevation
!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j) COLLAPSE(2)    
      do j=1,nwavej
        do i=1,nwavei
          depwave(i,j) = depwave0(i,j) !Important to use initial depths for extrapolation
        enddo           
      enddo  
!$OMP END DO 
!$OMP DO PRIVATE(ii)       
      do ii=1,ncellsD
        varflow(ii) = -zb(ii)
      enddo 
!$OMP END DO
!$OMP END PARALLEL
      call interp_scal_flwav(varflow,depwave,2) !Extrapolate to original values in depwave0
      
    case(2) !Bed elevation change
!$OMP PARALLEL
!$OMP DO PRIVATE(i,j) COLLAPSE(2)    
      do j=1,nwavej
        do i=1,nwavei
          depwave(i,j) = 0.0 !!Bed change. Important to initialize
        enddo           
      enddo  
!$OMP END DO 
!$OMP DO PRIVATE(ii)       
      do ii=1,ncellsD
        varflow(ii) = zb(ii)-zb0(ii) !Bed change
      enddo 
!$OMP END DO
!$OMP END PARALLEL
      call interp_scal_flwav(varflow,depwave,1) !Extrapolate to zero. depwave contains bed change
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)    
      do j=1,nwavej
        do i=1,nwavei
          depwave(i,j) = depwave0(i,j) + depwave(i,j) !Important to use initial depths for extrapolation
        enddo           
      enddo  
!$OMP END PARALLEL DO

    case(3) !Temporal interpolation in time using bathymetry datasets
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)    
      do j=1,nwavej
        do i=1,nwavei
          depwave(i,j) = depwave0(i,j) !Important to use initial depths for extrapolation
        enddo           
      enddo  
!$OMP END PARALLEL DO
      call geo_bathy_update(varflow,1)
      call interp_scal_flwav(varflow,depwave,2) !Extrapolate to original values in depwave0
      
    endselect
    
    !---Water levels ------------
    selectcase(noptwse)
    case(0) !None
      call diag_print_message('   wse(wave_time,wave_grid) = 0.0 m')
      !Convert p to eta and check dry cells        
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)
      do j=1,nwavej  
        do i=1,nwavei
          etawave(i,j) = 0.0
        enddo           
      enddo  
!$OMP END PARALLEL DO
        
    case(1) !Last time step
      call diag_print_message('   wse(wave_time,wave_grid) = wse(flow_time,flow_grid)')
      
      call interp_scal_flwav(p,etawave,3) !Notes: ***** etawave is pressure here ****, Extrapolated out
      !Convert p to eta and check dry cells
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)        
      do j=1,nwavej
        do i=1,nwavei
          etawave(i,j) = etawave(i,j)*gravinv  !Convert from p to eta
        enddo           
      enddo
!$OMP END PARALLEL DO

    case(2) !Tidal
      write(msg,*) '   wse(wave_time,wave_grid) = tide(flow_time) = ',tide2,' m'
      call diag_print_message(msg)
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)        
      do i=1,nwavei
        do j=1,nwavej
          etawave(i,j) = tide2
        enddo           
      enddo  
!$OMP END PARALLEL DO

    case default !(3) !Tidal Plus Variation
      write(msg,*)      '   tide(flow_time) = ',tide1,'tide(wave_time) = ', tide2,' m'
      call diag_print_message('   wse(wave_time,wave_grid) = wse(flow_time,flow_grid)',&
                              '            + tide(wave_time) - tide(flow_time)',msg)
      call interp_scal_flwav(p,etawave,3) !Notes: ***** etawave is pressure here ****, Extrapolated
      !Convert p to eta and check dry cells
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)         
      do j=1,nwavej
        do i=1,nwavei
          etawave(i,j) = tide2 + (etawave(i,j)*gravinv-tide1) !Convert from p to eta
        enddo           
      enddo
!$OMP END PARALLEL DO      
    endselect  
    
    !---Current Velocities ------------
    selectcase(noptvel)  
    case(0) !None
      call diag_print_message('   vel(wave_time,wave_grid)=0.0')
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)        
      do j=1,nwavej
        do i=1,nwavei
          uwave(i,j)=0.0
          vwave(i,j)=0.0  
        enddo
      enddo  
!$OMP END PARALLEL DO        

    case(1) !Last  
      call diag_print_message('   vel(wave_time,wave_grid)=vel(flow_time,flow_grid)')
#ifdef DEV_MODE
      selectcase(ivelwav)   
      case(1) !Mean
#endif
!$OMP PARALLEL DO PRIVATE(i)  
        do i=1,ncellsD
          ueff(i)=u(i)-us(i)
          veff(i)=v(i)-vs(i)
        enddo
!$OMP END PARALLEL DO
#ifdef DEV_MODE
      case(2) !Surfaces
        call q3d_flow_vel_surface(ueff,veff)
      case(3) !Weighted
        !call q3d_flow_vel_weighted!(ueff,veff)
      endselect
#endif
      call interp_vec_flwav(ueff,veff,uwave,vwave,1) !Extrapolated to zero
      theta=(azimuth_fl-azimuth)*deg2rad
      coswave=cos(theta)
      sinwave=sin(theta)      
!$OMP PARALLEL DO PRIVATE(i,j,uwaveij,vwaveij) COLLAPSE(2)      
      do j=1,nwavej
        do i=1,nwavei
          !Rotate velocities to wave grid orientation
          uwaveij=uwave(i,j)*coswave-vwave(i,j)*sinwave
          vwaveij=uwave(i,j)*sinwave+vwave(i,j)*coswave
          uwave(i,j)=uwaveij
          vwave(i,j)=vwaveij
        enddo
      enddo
!$OMP END PARALLEL DO        
    endselect
      
    !--- Water depths and wetting and drying ------
!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)     
    do j=1,nwavej
      do i=1,nwavei
        hwave(i,j) = depwave(i,j) + etawave(i,j)
        !Note: wetwave indicates whether a cell was wet during the flow simulation
        if(hwave(i,j)<hdry .or. wetwave(i,j)<0.001)then
          wetwave(i,j) = 0.0  
          uwave(i,j) = 0.0
          vwave(i,j) = 0.0
          etawave(i,j) = 0.0
          hwave(i,j) = 1.0e-5 !To avoid divide by zero
        endif
        !Save communication variables
        uin(i,j) = uwave(i,j) 
        vin(i,j) = vwave(i,j) !Save communication variable
        etain(i,j) = etawave(i,j) !Save communication variable
        if(noptzb>0 .and. wetwave(i,j)>0.95)then !only update cells that were wet during steering period
          depin(i,j)=depwave(i,j)
        endif    
      enddo
    enddo
!$OMP END PARALLEL DO
    
!    ierr=0
!    call check_wave_var(nwavei,nwavej,uwave,'uwave',ierr)
!    call check_wave_var(nwavei,nwavej,vwave,'vwave',ierr)
!    call check_wave_var(nwavei,nwavej,etawave,'etawave',ierr)
!    call check_wave_var(nwavei,nwavej,depwave,'depwave',ierr)
!    if(ierr/=0) stop    

    write(*,'(A)') ' Flow-to-wave interpolation complete'   
    write(*,*)   
    
    return
    
!    !--- Write communication files --------
!    open (unit=216,file= CurrFile, status='replace')   
!    write(216,*) nwavei,nwavej,0
!    write(216,*) IDate
!    do j=nwavej,1,-1
!      write (216,*) (uwave(i,j),vwave(i,j),i=1,nwavei) 
!    enddo
!    close(216)    
!
!    open (unit=221,file= SurgeFile, status='replace')       
!    write(221,*) nwavei,nwavej,0
!    write(221,*) IDate
!    do j=nwavej,1,-1
!      write(221,*) (etawave(i,j),i=1,nwavei)
!    enddo
!    close(221)
!    return
!    
!    if(sedtrans)then
!      nn=len_trim(DepFile)
!      MorphFile = DepFile(1:nn-4) // '_dep.dep'
!      open(unit=254,file=MorphFile, status='replace')       
!      write(254,*) nwavei,nwavej,0
!      do j=nwavej,1,-1
!        write(254,*) (depwave(i,j),i=1,nwavei)
!      enddo
!      close(254)
!    endif
!    return
!    
!    nn=len_trim(DepFile)
!    text=DepFile(1:nn-4) // '_dry.dep'
!    open(unit=254,file=text, status='replace')       
!    write(254,*) nwavei,nwavej,0
!    do j=nwavej,1,-1
!      write(254,*) (wetwave(i,j),i=1,nwavei)
!    enddo
!    close(254)
             
    return
    endsubroutine setwave

!***********************************************************************
    subroutine interp_scal_wavfl(varwave,varflow)
! Interpolates a variable from the wave mesh to the flow mesh
! written by Weiming Wu
! modified by Alex Sanchez
!***********************************************************************
    use size_def, only: ncellsD
    use wave_wavegrid_def, only: nwavei,nwavej
    use wave_flowgrid_def, only: ijwavcell,coefintp_wavfl
    use interp_lib
    use prec_def
    implicit none    
    !Input
    real(ikind),intent(in) :: varwave(nwavei,nwavej)   
    !Output
    real(ikind),intent(out) :: varflow(ncellsD)
    !Internal
    integer :: iextrap  
    integer :: iwetwave(nwavei,nwavej)
    real(ikind) :: valdry    
    
    iextrap = 1  !All wave variables are extrapolated to zero
    valdry = -999.0 !Value for dry cells (Not used because all wave cells are considered wet/active)
    iwetwave = 1 !All cells are considered active for interolation   
    call interp_scal_cart2tel(nwavei,nwavej,varwave,iwetwave,ijwavcell,coefintp_wavfl,&
           ncellsD,ncellsD,varflow,valdry,iextrap) !Note: Ghost cells are also interpolated

    return 
    endsubroutine interp_scal_wavfl
    
!***********************************************************************
    subroutine interp_vec_wavfl(vecxwave,vecywave,vecxflow,vecyflow)
! Interpolates a variable from the wave mesh to the flow mesh
! written by Weiming Wu
! modified by Alex Sanchez
!***********************************************************************
    use size_def, only: ncellsD
    use wave_wavegrid_def, only: nwavei,nwavej
    use wave_flowgrid_def, only: ijwavcell,coefintp_wavfl
    use interp_lib
    use prec_def
    implicit none    
    !Input
    real(ikind),intent(in),dimension(nwavei,nwavej) :: vecxwave,vecywave
    !Output
    real(ikind),intent(out),dimension(ncellsD) :: vecxflow,vecyflow
    !Internal
    integer :: iextrap  
    integer :: iwetwave(nwavei,nwavej)
    real(ikind) :: valdry    
    
    iextrap = 1  !All wave variables are extrapolated to zero
    valdry = -999.0 !Value for dry cells (Not used because all wave cells are considered wet/active)
    iwetwave = 1 !All cells are considered active for interolation   
    call interp_vec_cart2tel(nwavei,nwavej,vecxwave,vecywave,iwetwave,ijwavcell,coefintp_wavfl,&
           ncellsD,ncellsD,vecxflow,vecyflow,valdry,iextrap) !Note: Ghost cells are also interpolated

    return 
    endsubroutine interp_vec_wavfl    
            
!***********************************************************************
    subroutine getwave()
! Interface between CMS-Wave and CMS-Flow
! Gets wave variables from CMS-Wave, passes them to CMS-Flow
! variables on the wave grid and then interpolates them to the flow grid
! 
! written by Weiming Wu, and Alex Sanchez
!***********************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD
    use flow_def, only: rhow,waveflux,hdry,grav
    use comvarbl, only: nfsch
    use geo_def, only: azimuth_fl
    use flow_wavegrid_def, only: hwave
    use wave_flowgrid_def, only: wavestrx2,wavestry2,whgt2,wper2,&
       wunitx2,wunity2,wavediss2
    use wave_wavegrid_def, only: nwavei,nwavej,azimuth_wav,&
       wxrs1,wyrs1,wheight,wperiod,wcos,wsin,wdiss,idatewave
    use cms_def, only: ndissm,npersm,wavedisstol
    use global_inline
    use rol_def, only: roller
    use sed_def, only: sedtrans,wavesedtrans
    use const_def, only: deg2rad
    use time_lib, only: index2calendar,calendar2julian
    use diag_lib
    use prec_def
    implicit none
    common /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
    common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz
    common /rsb/ wxrs(ipmx,jpmx),wyrs(ipmx,jpmx)
    common /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
    common /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
    common /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
    integer :: i,j,iwave,jwave !,ierr
    real(ikind) :: thetadgr, g
    !real :: sfac(nwavei,nwavej),svar(nwavei,nwavej)
    real(4) :: DX,DY,DXX,DMESH,DTH,depmax,depmax0
    integer :: kdate,idate
    integer :: imod,iprp,island,imd,iprpp,nonln,igrav
    integer :: isolv,ixmdf,iproc,imud,iwnd
    integer :: isteer,iidate,IBR
    real(4) :: x0,y0,azimuth,sinaz,cosaz
    real(4) :: wxrs,wyrs,H13,T13,DMN,H13S,DISS
    real(4) :: PAI2,PAI,HPAI,RAD,akap,depmin0
    integer :: iyrwav,imowav,idaywav,ihrwav,iminwav,isecwav,ierr
    real(ikind) :: tjulwav
    real(ikind), parameter :: wpermin = 1.0
    real(ikind), parameter :: valdry = -999.0
#ifdef DIAG_MODE
    character(len=100) :: msg2
    logical :: isnankind
#endif
    
    write(*,'(A)') ' Starting wave-to-flow interpolation'
    write(*,'(A,I3)') '  imod=', imod
    
    idatewave = idate
    
    call index2calendar(idatewave,&
      iyrwav,imowav,idaywav,ihrwav,iminwav,isecwav,ierr)
    call calendar2julian(iyrwav,imowav,idaywav,ihrwav,iminwav,isecwav,tjulwav) 
    
    if(iyrwav .ge. 2000 .or. idatewave .gt. 10000000)then !Otherwise it is not a date
      write(*,'(A,I4,5(1x,I2))') '  Wave date = ',&
        iyrwav,imowav,idaywav,ihrwav,iminwav,isecwav
    else
      write(*,'(A,I10)') '  Wave index = ',idatewave
    endif
    
    if(imod==0)then
!$OMP PARALLEL DO PRIVATE(i,j,iwave,jwave)             
      do iwave=1,nwavei
        do jwave=1,nwavej
          i=iwave
          j=jwave
          wxrs1(iwave,jwave)=wxrs(i,j) !positive in x-direction of wave grid
          wyrs1(iwave,jwave)=wyrs(i,j) !positive in y-direction of wave grid
          wheight(iwave,jwave)=H13S(i,j)    
          wperiod(iwave,jwave)=max(T13(i,j),wpermin)
          wcos(iwave,jwave)=cos(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          wsin(iwave,jwave)=sin(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          !wibr(iwave,jwave)=max(real(ibr(i,j)),0.0)
          !wdiss(iwave,jwave)=-rhow*diss(i,j)*ibr(i,j) !Alex, wdiss=[N/m/s]
          wdiss(iwave,jwave) = diss(i,j)   !Lihwa  12/18/2020
        enddo
      enddo  
!$OMP END PARALLEL DO
    elseif(imod==2)then
!$OMP PARALLEL DO PRIVATE(i,j,iwave,jwave)
      do iwave=1,nwavei
        do jwave=1,nwavej
          i=nwavei-iwave+1
          j=nwavej-jwave+1
          wxrs1(iwave,jwave)=wxrs(i,j) !positive in x-direction of wave grid
          wyrs1(iwave,jwave)=wyrs(i,j) !positive in y-direction of wave grid
          wheight(iwave,jwave)=H13S(i,j)    
          wperiod(iwave,jwave)=max(T13(i,j),wpermin)
          wcos(iwave,jwave)=cos(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          wsin(iwave,jwave)=sin(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          !wibr(iwave,jwave)=max(real(ibr(i,j)),0.0)
          !wdiss(iwave,jwave)=-rhow*diss(i,j)*ibr(i,j) !Alex, wdiss=[N/m/s]
          wdiss(iwave,jwave)=diss(i,j)   !Lihwa  12/18/2020
        enddo
      enddo  
!$OMP END PARALLEL DO
    elseif(imod==1)then
!$OMP PARALLEL DO PRIVATE(i,j,iwave,jwave)          
      do iwave=1,nwavei
        do jwave=1,nwavej
          i=jwave
          j=nwavei-iwave+1
          wxrs1(iwave,jwave)=wxrs(i,j) !positive in x-direction of wave grid
          wyrs1(iwave,jwave)=wyrs(i,j) !positive in y-direction of wave grid
          wheight(iwave,jwave)=H13S(i,j)    
          wperiod(iwave,jwave)=max(T13(i,j),wpermin)
          wcos(iwave,jwave)=cos(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          wsin(iwave,jwave)=sin(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          !wibr(iwave,jwave)=max(real(ibr(i,j)),0.0)
          !wdiss(iwave,jwave)=-rhow*diss(i,j)*ibr(i,j) !Alex, wdiss=[N/m/s]
          wdiss(iwave,jwave)=diss(i,j)   !Lihwa  12/18/2020
        enddo
      enddo  
!$OMP END PARALLEL DO
    else !if(imod==3)then
!$OMP PARALLEL DO PRIVATE(i,j,iwave,jwave)
      do iwave=1,nwavei
        do jwave=1,nwavej
          i=nwavej-jwave+1
          j=iwave
          wxrs1(iwave,jwave)=wxrs(i,j) !positive in x-direction of wave grid
          wyrs1(iwave,jwave)=wyrs(i,j) !positive in y-direction of wave grid
          wheight(iwave,jwave)=H13S(i,j)    
          wperiod(iwave,jwave)=max(T13(i,j),wpermin)
          wcos(iwave,jwave)=cos(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          wsin(iwave,jwave)=sin(dmn(i,j)*deg2rad) !going, Cartesian, 0-east,90-north,180-west,270-south
          !wibr(iwave,jwave)=max(real(ibr(i,j)),0.0)
          !wdiss(iwave,jwave)=-rhow*diss(i,j)*ibr(i,j) !Alex, wdiss=[N/m/s]
          wdiss(iwave,jwave)=diss(i,j)   !Lihwa  12/18/2020
        enddo
      enddo 
!$OMP END PARALLEL DO
    endif    

    wdiss=-rhow*wdiss*grav		!Lihwa  12/18/2020
    
#ifdef DIAG_MODE
    ierr=0
    call check_wave_var(nwavei,nwavej,wxrs1,'wxrs1',ierr)
    call check_wave_var(nwavei,nwavej,wyrs1,'wyrs1',ierr)
    call check_wave_var(nwavei,nwavej,wheight,'wheight',ierr)
    call check_wave_var(nwavei,nwavej,wperiod,'wperiod',ierr)
    call check_wave_var(nwavei,nwavej,wcos,'wcos',ierr)
    call check_wave_var(nwavei,nwavej,wsin,'wsin',ierr)    
    if(ierr/=0)then
      call diag_print_error('Instability in wave model',&
        '  Check wave model solution')
    endif
#endif
    
    !call wave_rad_stress !Overrides the wave radiation stresses from the wave model
    !if(waveflux) call stokes_stress !Adds the stokes stress gradients to the wave forcing
    !call wave_dissipation !Over-rides the wave dissipation from the wave model

    !where(wdiss<=wavedisstol)
    !  wibr=0.0 !If no dissipation than there must be no breaking
    !endwhere  

    !Smoothing factor for transition from breaking to nonbreaking
    !!call smooth_wavegrid_scal(wibr,nbrksm,0) !Replaced by below procedure
    !do i=1,nbrksm
    !  where(hwave<=hdry)
    !    wibr=1.0 !Set all dry cells as breaking
    !  endwhere  
    !  call smooth_wavegrid_scal(wibr,1,0)
    !enddo
    if(ndissm>0) call smooth_wavegrid_scal(wdiss,ndissm,0)
    !if(nradsm>0) call smooth_wavegrid_vec(wxrs1,wyrs1,nradsm)
    !wxrs1=wxrs1*wibr !Limit radiation stresses to breaking zone
    !wyrs1=wyrs1*wibr !Limit radiation stresses to breaking zone
    if(npersm>0) call smooth_wavegrid_scal(wperiod,npersm,0) !Wave period smoothing

    !--- Add roller contribution to stresses ----    
    if(roller) call rol_solve !Note: roller model is on wave grid
    
    !--- Spatial Interpolation from Wave to Flow ------
    call interp_vec_wavfl(wxrs1,wyrs1,wavestrx2,wavestry2) !Extrapolate to zero
    call interp_scal_wavfl(wheight,Whgt2)             !Extrapolate to zero
    call interp_scal_wavfl(wperiod,Wper2)             !Extrapolate to zero
    call interp_vec_wavfl(wcos,wsin,Wunitx2,Wunity2)  !Extrapolate to zero
    !call interp_scal_wavfl(wibr,waveibr2)             !Extrapolate to zero
    call interp_scal_wavfl(wdiss,wavediss2)           !Extrapolate to zero

    !Check for negative wave heights from bad interpolations
!$OMP PARALLEL DO PRIVATE(i)    
    do i=1,ncells
      if(Whgt2(i)<0.0)    Whgt2(i) = 0.0
      !if(waveibr2(i)<0.0) waveibr2(i) = 0.0
      if(Wper2(i)<1.0)    Wper2(i) = 1.0
#ifdef DIAG_MODE
      if(isnankind(wavestrx2(i)) .or. isnankind(wavestry2(i)))then
        write(msg2,*) 'i = ',i,wavestrx2(i),wavestry2(i)
        call diag_print_error('Found NaN in Radiation stress',msg2,&
           '  Setting value to zero.')
      endif
#endif
    enddo  
!$OMP END PARALLEL DO  

    !--- Wave stresses ---
    thetadgr=azimuth_fl-azimuth_wav
    if(imod==1)then
      thetadgr=thetadgr-90.0      !Lihwa
    elseif(imod==2)then
      thetadgr=thetadgr-180.0     !Lihwa
    elseif(imod==3)then
      thetadgr=thetadgr+90.0      !Lihwa        
    endif
    call rotate_vector(ncellsD,ncellsD,thetadgr,wavestrx2,wavestry2) 

    !--- Wave unit vectors -----
    !!thetadgr=azimuth_fl-azimuth !Note: should use the same angle as for the wave stresses
    call rotate_vector(ncellsD,ncellsD,thetadgr,Wunitx2,Wunity2)    
    
    write(*,'(A)') ' Wave-to-flow interpolation complete'

    return 
    endsubroutine getwave

!**********************************************************************
    subroutine check_wave_var(ni,nj,val,name,ierr)
!***********************************************************************       
    use diag_def, only: dgfile,dgunit
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: ni,nj
    integer, intent(out) :: ierr
    real(ikind), intent(inout) :: val(ni,nj)
    character(len=*), intent(in) :: name
    !Internal variables
    integer :: i,j    
    logical :: isnankind

    ierr = 0
    open(dgunit,file=dgfile,access='append')
431 format(A,'(',I5,I5,')=',F12.5)
    do i=1,ni
      do j=1,nj
        if(isnankind(val(i,j)))then
          ierr = -1
          write(*,431)      name,i,j,val(i,j)
          write(dgunit,431) name,i,j,val(i,j)
          val(i,j) = 0.0
        endif
      enddo      
    enddo
    close(dgunit)
    
    return 
    endsubroutine check_wave_var

!**********************************************************************
    subroutine smooth_wavegrid_scal(val,niter,ibc)
! Smooths a wave variables on the wave grid
! by using the following stencil
!        1/6
!   1/6  1/3  1/6
!        1/6
! written by Alex Sanchez, USACE-CHL    
!***********************************************************************        
    use wave_wavegrid_def, only: nwavei,nwavej
    use flow_wavegrid_def, only: hwave
    use flow_def, only: hdry
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: niter,ibc
    real(ikind), intent(inout) :: val(nwavei,nwavej)
    !Internal variables
    integer :: i,i1,i2,j,j1,j2,k
    real(ikind) :: val2(nwavei,nwavej)

    selectcase(ibc)
    case(0)
    do k=1,niter
!$OMP PARALLEL        
!$OMP DO PRIVATE(i,i1,i2,j,j1,j2)    
      do i=1,nwavei
        i1=min(i+1,nwavei)
        i2=max(i-1,1)
        do j=1,nwavej
          j1=min(j+1,nwavej)
          j2=max(j-1,1)
          val2(i,j)=(2.0_ikind*val(i,j)+val(i1,j)+val(i2,j)+val(i,j1)+val(i,j2))/6.0_ikind
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j) COLLAPSE(2) 
      do i=1,nwavei
        do j=1,nwavej
          val(i,j) = val2(i,j)
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
    enddo
    
    case(1)    
    do k=1,niter
!$OMP PARALLEL        
!$OMP DO PRIVATE(i,i1,i2,j,j1,j2)    
      do i=1,nwavei
        i1=min(i+1,nwavei)
        i2=max(i-1,1)
        do j=1,nwavej
          j1=min(j+1,nwavej)
          j2=max(j-1,1)
          if(hwave(i,j)>hdry)then
            val2(i,j)=val(i,j)/3.0_ikind  
            if(hwave(i1,j)>hdry)then
              val2(i,j)=val2(i,j)+val(i1,j)/6.0_ikind
            else
              val2(i,j)=val2(i,j)+val(i,j)/6.0_ikind
            endif
            if(hwave(i2,j)>hdry)then
              val2(i,j)=val2(i,j)+val(i2,j)/6.0_ikind
            else
              val2(i,j)=val2(i,j)+val(i,j)/6.0_ikind
            endif
            if(hwave(i,j1)>hdry)then
              val2(i,j)=val2(i,j)+val(i,j1)/6.0_ikind
            else
              val2(i,j)=val2(i,j)+val(i,j)/6.0_ikind
            endif
            if(hwave(i,j2)>hdry)then
              val2(i,j)=val2(i,j)+val(i,j2)/6.0_ikind
            else
              val2(i,j)=val2(i,j)+val(i,j)/6.0_ikind  
            endif
          endif
          !val2(i,j)=(2.0_ikind*val(i,j)+val(i1,j)+val(i2,j)+val(i,j1)+val(i,j2))/6.0_ikind
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j) COLLAPSE(2) 
      do i=1,nwavei
        do j=1,nwavej
          val(i,j) = val2(i,j)
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
    enddo  
    
    endselect
    
    return
    endsubroutine smooth_wavegrid_scal
    
!**********************************************************************
    subroutine smooth_wavegrid_vec(vecx,vecy,niter)
! Smooths a wave variables on the wave grid
! by using the following stencil
!        1/6
!   1/6  1/3  1/6
!        1/6
! written by Alex Sanchez, USACE-CHL    
!***********************************************************************        
    use wave_wavegrid_def, only: nwavei,nwavej
    use prec_def
    implicit none
    !Input/Output    
    integer, intent(in) :: niter
    real(ikind), intent(inout) :: vecx(nwavei,nwavej),vecy(nwavei,nwavej)
    !Internal
    integer :: i,i1,i2,j,j1,j2,k
    real(ikind) :: vecx2(nwavei,nwavej),vecy2(nwavei,nwavej)

    do k=1,niter
!$OMP PARALLEL        
!$OMP DO PRIVATE(i,i1,i2,j,j1,j2)    
      do i=1,nwavei
        i1=min(i+1,nwavei)
        i2=max(i-1,1)
        do j=1,nwavej
          j1=min(j+1,nwavej)
          j2=max(j-1,1)    
          vecx2(i,j)=(2.0_ikind*vecx(i,j)+vecx(i1,j)+vecx(i2,j)+vecx(i,j1)+vecx(i,j2))/6.0_ikind
          vecy2(i,j)=(2.0_ikind*vecy(i,j)+vecy(i1,j)+vecy(i2,j)+vecy(i,j1)+vecy(i,j2))/6.0_ikind
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j) COLLAPSE(2) 
      do i=1,nwavei
        do j=1,nwavej
          vecx(i,j) = vecx2(i,j)
          vecy(i,j) = vecy2(i,j)
        enddo
      enddo
!$OMP END DO
!!$OMP DO PRIVATE(i,i1,i2,j,j1,j2)    
!      do i=2,nwavei-1
!        i1=i+1
!        i2=i-1
!        do j=2,nwavej-1
!          j1=j+1
!          j2=j-1   
!          vecx2(i,j)=(2.0_ikind*vecx(i,j)+vecx(i1,j)+vecx(i2,j)+vecx(i,j1)+vecx(i,j2))/6.0_ikind
!          vecy2(i,j)=(2.0_ikind*vecy(i,j)+vecy(i1,j)+vecy(i2,j)+vecy(i,j1)+vecy(i,j2))/6.0_ikind
!        enddo
!      enddo
!!$OMP END DO
!!$OMP DO PRIVATE(i,j) COLLAPSE(2) 
!      do i=2,nwavei-1
!        do j=2,nwavej-1
!          vecx(i,j) = vecx2(i,j)
!          vecy(i,j) = vecy2(i,j)
!        enddo
!      enddo
!!$OMP END DO
!$OMP END PARALLEL
    enddo    
    
    return
    endsubroutine smooth_wavegrid_vec    
    
!**********************************************************************
    subroutine smooth_flowgrid_scal(val,niter)
! Smooths a wave variables on the flow grid
! written by Alex Sanchez, USACE-CHL       
!***********************************************************************    
    use size_def
    use geo_def, only: cell2cell,ncface
    use flow_def, only: iwet
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: niter
    real(ikind), intent(inout) :: val(ncellsD)
    !Internal variables
    integer :: i,j,k,nck
    real(ikind) :: valout(ncellsD),cn

    do j=1,niter
!$OMP PARALLEL        
!$OMP DO PRIVATE(i,k,cn,nck)    
      do i=1,ncells
        if(iwet(i)==0) cycle  
        cn = 1.0
        valout(i) = val(i)
        do k=1,ncface(i)
          nck = cell2cell(k,i)
          if(nck<=ncells .and. iwet(nck)==1)then
            cn = cn + 1.0
            valout(i) = valout(i) + val(nck)
          endif
        enddo
        valout(i)=valout(i)/cn 
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i)    
      do i=1,ncells
        if(iwet(i)==0) cycle   
        val(i) = valout(i)
      enddo
!$OMP END DO
!$OMP END PARALLEL      
    enddo
    
    return 
    endsubroutine smooth_flowgrid_scal
    
!**********************************************************************
    subroutine smooth_flowgrid_vec(vecx,vecy,niter)
! Smooths a wave variables on the flow grid
! written by Alex Sanchez, USACE-CHL       
!***********************************************************************    
    use size_def
    use geo_def, only: cell2cell,ncface
    use flow_def, only: iwet
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: niter
    real(ikind),intent(inout) :: vecx(ncellsD),vecy(ncellsD)
    !Internal variables
    integer :: i,j,k,nck    
    real(ikind) :: vecxout(ncellsD),vecyout(ncellsD),cn

    do j=1,niter
!$OMP PARALLEL        
!$OMP DO PRIVATE(i,k,cn,nck)    
      do i=1,ncells
        if(iwet(i)==0) cycle  
        cn = 1.0
        vecxout(i) = vecx(i)
        vecyout(i) = vecy(i)
        do k=1,ncface(i)
          nck = cell2cell(k,i)
          if(nck<=ncells .and. iwet(nck)==1)then
            cn = cn + 1.0
            vecxout(i) = vecxout(i) + vecx(nck)
            vecyout(i) = vecyout(i) + vecy(nck)    
          endif
        enddo
        vecxout(i) = vecxout(i)/cn
        vecyout(i) = vecyout(i)/cn 
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i)    
      do i=1,ncells
        if(iwet(i)==0) cycle   
        vecx(i) = vecxout(i)
        vecy(i) = vecyout(i)
      enddo
!$OMP END DO
!$OMP END PARALLEL      
    enddo
    
    return 
    endsubroutine smooth_flowgrid_vec    

!***********************************************************************
    subroutine wave_eval()
! Temporal interpolation of wave variables
! written by Weiming Wu, NCCHE (now at Clarkson University)
! and Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h" 
    use size_def
    use geo_def
    use flow_def
    use comvarbl
    use wave_flowgrid_def
    use cms_def
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalsteph5,readvecsteph5
#endif   
    use rol_def, only: roller,Sr,rolflux
    use sed_def, only: sedtrans,wavesedtrans,dzb
    use const_def, only: deg2rad,small
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    integer :: i,ii,ierr !,k
    real(ikind) :: facintep,val !,wavestrnorm
 
    if(ctime>tswave2 .and. noptset==3)then
      call diag_print_message(' ','*** CMS-Flow run successful ***')
            
      nsteer=nsteer+1
      tswave1=tswave2
!!        tide1=tide2
      call tidevalue(ctime,tide1)   !Tide1 is tidal level at ctime, not at previous steering interval, 
                                          !because etawave is at ctime
!$OMP PARALLEL DO PRIVATE(ii)         
      do ii=1,ncellsD
        wavestrx1(ii)=wavestrx2(ii)
        wavestry1(ii)=wavestry2(ii)
        Whgt1(ii)=Whgt2(ii)
        Wper1(ii)=Wper2(ii)     
        Wunitx1(ii)=Wunitx2(ii)
        Wunity1(ii)=Wunity2(ii)
        !waveibr1(ii)=waveibr2(ii)
        wavediss1(ii)=wavediss2(ii)  !Alex, Aug 3, 2009
      enddo
!$OMP END PARALLEL DO
!!      if(roller .and. sedtrans .and. wavesedtrans) Esr1=Esr2
            
      tswave2=tswave2+dtsteer
      call tidevalue(tswave2,tide2)
      call setwave
      
      call diag_print_message('*** Starting CMS-Wave run ***',' ')      
      open(dgunit,file=dgfile,access='append') !Open before CMS-Wave run
      
      call cms_wave_inline                         !call cms_wave !(noptset,nsteer)  !modified 10/15/2018
      
      close(dgunit) !Close after CMS-Wave run
      call diag_print_message(' ','*** CMS-Wave run successful ***')
        
!Lihwa added  12/18/2020
!$  if(nthr>=1) call omp_set_num_threads(nthr)

      write(*,*) 'NTHR =', nthr

      call getwave
      
      if(sedtrans) wetsteer = 0.0 !Restart wet count for wet cells, Alex

      call diag_print_message(' ','*** Starting CMS-Flow Run ***')
      if(tswave2<tswave1) then
        write(msg,*) IDate0 
        call diag_print_error('Wave boundary input time series is wrong at: ',msg)
      endif
    endif
    
    if(ctime>(tswave2-1.0e-6) .and. noptset==4)then
      nsteer=nsteer+1
#ifdef XMDF_IO
      call readscalsteph5(wgrdfile,wavpath,nsteer,tswave2,Whgt2,ierr)                 !Updated with 'wgrdfile' to use since there is no _grid.h5 file anymore.  MEB  06/10/2021 
#else
      call diag_print_error('ERROR: Cannot read initial wave condition without XMDF')
#endif     
      if(ierr>=0)then
#ifdef XMDF_IO
        call readscalsteph5(wgrdfile,perpath,nsteer,tswave2,Wper2,ierr)               !Updated with 'wgrdfile' to use since there is no _grid.h5 file anymore.  MEB  06/10/2021
        call readscalsteph5(wgrdfile,dirpath,nsteer,tswave2,Wang,ierr)
        call readscalsteph5(wgrdfile,wavpath,nsteer,tswave2,wavediss2,ierr)
        call readvecsteph5 (wgrdfile,wavpath,nsteer,tswave2,wavestrx2,wavestry2,ierr)
#endif       
!$OMP PARALLEL DO PRIVATE(i)
          do i=1,ncellsD
            Wunitx2(i) = cos((Wang(i)-azimuth_fl)*deg2rad)
          Wunity2(i) = sin((Wang(i)-azimuth_fl)*deg2rad)
          wavediss2(i) = -rhow*wavediss2(i)  !Flip sign and change units [N/m/s]
          !if(wavediss2(i)>wavedisstol)then
          !  waveibr2(i) = 1.0
          !endif
        enddo
!$OMP END PARALLEL DO    
!!        call rotate_vector(ncellsD,ncellsD,azimuth_fl,wavestrx2,wavestry2)                   
!!        call rotate_vector(ncellsD,ncellsD,azimuth_fl,Wunitx2,Wunity2)   

        !Smoothing
        !do i=1,nbrksm
        !  where(iwet==0)
        !    waveibr2 = 1.0 !Set all dry cells as breaking
        !  endwhere  
        !  call smooth_flowgrid_scal(waveibr2,1)
        !enddo
        if(ndissm>0) call smooth_flowgrid_scal(wavediss2,ndissm)
        !if(nradsm>0) call smooth_flowgrid_vec (wavestrx2,wavestry2,nradsm)
        if(npersm>0) call smooth_flowgrid_scal(Wper2,npersm) !Wave period smoothing
      
!$OMP PARALLEL DO PRIVATE(i)
        do i=1,ncellsD
          if(wavediss2(i)<wavedisstol)then
            wavediss2(i) = 0.0
          endif
          !wavestrx2(i) = wavestrx2(i)*waveibr2(i) !Limit radiation stresses to breaking zone
          !wavestry2(i) = wavestry2(i)*waveibr2(i) !Limit radiation stresses to breaking zone
        enddo
!$OMP END PARALLEL DO 
      endif  
    endif
    
    if(ctime<=tswave1) then
!$OMP PARALLEL DO PRIVATE(ii,val)      
      do ii=1,ncellsD
        wavestrx(ii)=wavestrx1(ii)
        wavestry(ii)=wavestry1(ii)
        Whgt(ii)=Whgt1(ii)
        Wper(ii)=Wper1(ii)
        Wunitx(ii)=Wunitx1(ii)
        Wunity(ii)=Wunity1(ii)
        val=sqrt(Wunitx(ii)**2+Wunity(ii)**2)
        
        val=max(val,1.0e-5)               !To fix a divide by zero issue, unknown date of implementation.
        Wunitx(ii)=Wunitx(ii)/val
        Wunity(ii)=Wunity(ii)/val
        
        Wang(ii)=atan2(Wunity(ii),Wunitx(ii))
        !waveibr(ii)=waveibr1(ii)
        wavediss(ii)=wavediss1(ii)  !Alex, Aug 3, 2009
      enddo
!$OMP END PARALLEL DO
    elseif(ctime>=tswave2) then
!$OMP PARALLEL DO PRIVATE(ii,val)        
      do ii=1,ncellsD
        wavestrx(ii)=wavestrx2(ii)
        wavestry(ii)=wavestry2(ii)
        Whgt(ii)=Whgt2(ii)
        Wper(ii)=Wper2(ii)
        Wunitx(ii)=Wunitx2(ii) !Bug fix, replaced Wunitx1 with Wunitx
        Wunity(ii)=Wunity2(ii) !Bug fix, replaced Wunity1 with Wunity
        val=sqrt(Wunitx(ii)**2+Wunity(ii)**2)
        val=max(val,1.0e-5)
        Wunitx(ii)=Wunitx(ii)/val
        Wunity(ii)=Wunity(ii)/val
        Wang(ii)=atan2(Wunity(ii),Wunitx(ii))
        !waveibr(ii)=waveibr2(ii)
        wavediss(ii)=wavediss2(ii)  !Alex, Aug 3, 2009
      enddo
!$OMP END PARALLEL DO 
!!      if(roller .and. sedtrans .and. wavesedtrans) Esr=Esr2
    else                                                                   !Flow time is between Wave intervals  MEB  01/21/2014
      facintep = (ctime-tswave1)/(tswave2-tswave1+small)
      if (.not. wave_interp) facintep = 0.0                                !Added to make weighting of second wave condition always equal to zero in this section.  MEB  01/22/2014
!$OMP PARALLEL DO PRIVATE(ii,val)         
      do ii=1,ncellsD
        wavestrx(ii)=wavestrx1(ii)+facintep*(wavestrx2(ii)-wavestrx1(ii))
        wavestry(ii)=wavestry1(ii)+facintep*(wavestry2(ii)-wavestry1(ii))
        Whgt(ii)=Whgt1(ii)+facintep*(Whgt2(ii)-Whgt1(ii))
        Wper(ii)=Wper1(ii)+facintep*(Wper2(ii)-Wper1(ii))
        Wunitx(ii)=Wunitx1(ii)+facintep*(Wunitx2(ii)-Wunitx1(ii))
        Wunity(ii)=Wunity1(ii)+facintep*(Wunity2(ii)-Wunity1(ii))
        val=sqrt(Wunitx(ii)**2+Wunity(ii)**2)
        val=max(val,1.0e-5)
        Wunitx(ii)=Wunitx(ii)/val
        Wunity(ii)=Wunity(ii)/val
        Wang(ii)=atan2(Wunity(ii),Wunitx(ii))
        !waveibr(ii)=waveibr1(ii)+facintep*(waveibr2(ii)-waveibr1(ii))
        wavediss(ii)=wavediss1(ii)+facintep*(wavediss2(ii)-wavediss1(ii))  !Alex, Aug 3, 2009
      enddo
!$OMP END PARALLEL DO 
!!      if(roller .and. sedtrans .and. wavesedtrans) Esr=Esr1+facintep*(Esr2-Esr1)
    endif

    !!do ii=1,ncellsD
    !!  if(iwet(ii)==0) then
    !!    wavestrx(ii)=0.0
    !!    wavestry(ii)=0.0  
    !!  endif
    !!enddo

    if(roller) call interp_scal_wavfl(Sr,Ssr) !Extrapolate to zero
    
    if(ramp<1.0e-10)then
      val=0.0
    else  
      val=sqrt(ramp)
    endif  
    
    !Apply ramp and check values
!$OMP PARALLEL    
!$OMP DO PRIVATE(i)
    do i=1,ncells
       Wper(i) = max(Wper(i),1.0)
       wavestrx(i) = wavestrx(i)*ramp !Alex, Mar 27, 2010 
       wavestry(i) = wavestry(i)*ramp
       Whgt(i) = Whgt(i)*val
    enddo
!$OMP END DO

    if(waveflux .and. roller)then
!$OMP DO PRIVATE(i)
      do i=1,ncells
        Ssr(i)=Ssr(i)*ramp
      enddo
    endif  
!$OMP END PARALLEL

    return 
    endsubroutine wave_eval
    
!***********************************************************************
    subroutine wave_wetdry()
! Impose wetting and drying on wave variables
! and limits variables
! written by Alex Sanchez, USAC-CHL
! last modified 10/01/2013
!***********************************************************************
    use size_def
    use flow_def, only: h,u,v,iwet,waveflux,us,vs,rhow,grav
    use wave_flowgrid_def
    use rol_def, only: roller,Sr,rolflux
    use sed_def, only: sedtrans
    use const_def, only: twopi
    use wave_lib, only: waveorb_linear,waveorbrep_ss87
    use prec_def
    implicit none
    integer :: i !,k
    !real(ikind) :: wavestrnorm
    real*4 :: uu,vv,om,d,q,cw,cg,sg,akk !used for wccg, must be single     

!$OMP PARALLEL    
!$OMP DO PRIVATE(i,q,d,uu,vv,om,cw,cg,sg,akk)
    do i=1,ncells       
       !Wave length and Orbital velocity 
       Whgt(i)=min(Whgt(i),0.88*h(i)) !limit wave height ******************
!      Wlen(i)=wavelength(Wper(i),h(i)) !No wave-current interaction
       om=twopi/Wper(i)       
       q=Wang(i); d=h(i); uu=u(i); vv=v(i);
       !call wccg(q,d,uu,vv,om,cw,cg,sg,akk) !Wave-current interaction
       call wccg_inline(q,d,uu,vv,om,cw,cg,sg,akk) !Wave-current interaction   !Lihwa 12/18/2020
       !call wccg3(q,d,uu,vv,om,cw,cg,sg,akk) !Wave-current interaction
       Wlen(i)=twopi/akk
!       f = 1.0/Wper(i); d=h(i)
!       uu1 = u(i)*Wunitx(i)+v(i)*Wunity(i)
!       call wkcgen(f,grav,pi,uu1,d,i,akk)
!       Wlen(i)=twopi/akk
       
       !Apply wet/dry
       wavestrx(i)=iwet(i)*wavestrx(i) !Alex, Mar 27, 2010 
       wavestry(i)=iwet(i)*wavestry(i)
       Whgt(i)=iwet(i)*Whgt(i)
       
       !Orbital velocities (includes ramp and wet/dry through Whgt)
       Worb(i)=waveorb_linear(h(i),Whgt(i),Wper(i),Wlen(i))
!!       Worbrep(i)=Worb(i)/sqrttwo !Old approach
       Worbrep(i)=waveorbrep_ss87(grav,h(i),Whgt(i),Wper(i)) !More accurate   
       !if(Worb(i)>Whgt(i) .or. Worbrep(i)>Worb(i))then
       !  write(*,*)
       !  write(*,*) 'Problem calculating orbital velocities'
       !  write(*,*) 'h(i) ',h(i)
       !  write(*,*) 'u(i) ', u(i)
       !  write(*,*) 'v(i) ', v(i)
       !  write(*,*) 'Worb(i) ', Worb(i)
       !  write(*,*) 'Worbrep(i) ',Worbrep(i)
       !  write(*,*) 'Whgt(i) ', Whgt(i)
       !  write(*,*) 'Wper(i) ', Wper(i)
       !  write(*,*) 'Wlen(i) ', Wlen(i)
       !  write(*,*) 'Wang(i) ', Wang(i)
       !endif
    enddo
!$OMP END DO

    !Wave flux velocity
    if(waveflux)then
      if(roller .and. rolflux)then
!$OMP DO PRIVATE(i,q,cw)
        do i=1,ncells
          cw=Wlen(i)/Wper(i)  
          Ssr(i)=iwet(i)*Ssr(i)   
          q=iwet(i)*(0.0625*grav*Whgt(i)*Whgt(i)+Ssr(i)/rhow)/h(i)/cw
          us(i)=q*Wunitx(i)
          vs(i)=q*Wunity(i)
        enddo
!$OMP END DO
      else
!$OMP DO PRIVATE(i,q,cw)
        do i=1,ncells
          cw=Wlen(i)/Wper(i)  
          q=iwet(i)*0.0625*grav*Whgt(i)*Whgt(i)/h(i)/cw !Stokes velocity 
          us(i)=q*Wunitx(i)
          vs(i)=q*Wunity(i)
        enddo
!$OMP END DO
      endif
    endif
!$OMP END PARALLEL        
        

    if(sedtrans)then
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        wetsteer(i) = max(wetsteer(i),real(iwet(i),kind=ikind))
      enddo
!$OMP END PARALLEL DO   
    endif
    
    return
    endsubroutine wave_wetdry

!***********************************************************************
    subroutine tidevalue(tswave2,tide2)
!   evaluation of Tide at boundaris for steering 
!   Made by Weiming Wu, NCCHE, April 2009
!***********************************************************************
    use size_def
    use flow_def
    use bnd_def
    use comvarbl, only: ramp
    use const_def, only: small
    use plagr_lib
    use prec_def
    implicit none
    integer :: j,k,iwse,ibnd,np,inc
    real(ikind) :: tswave2,tide2,tideiwse,timehrswave
    integer,parameter :: nb = 3 !Must be larger than maximum order of interpolation
    real(ikind) :: lb(nb)

    tide2=0.0
    ibnd=0
    timehrswave=tswave2/3600.0

!--- Tidal/Harmonic BC --------------------------------------------------------    
    do iwse=1,nTHstr
      tideiwse = TH_str(iwse)%wseoffset
      do k=1,TH_str(iwse)%ntc  
        tideiwse = tideiwse + TH_str(iwse)%amp(k) &
                 *cos(TH_str(iwse)%speed(k)*timehrswave + TH_str(iwse)%phase(k))  
      enddo
      j = TH_str(iwse)%ncells/2 !Use center cell as approximate value
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*TH_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo

!--- Single Water Level BC ---------------------------------------------------------        
    do iwse=1,nHstr  !for each cell string
      if(H_STR(iwse)%ntimes>0)then !Time series  
        inc = H_str(iwse)%inc  
        call plagr_fit(H_str(iwse)%ntimes,H_str(iwse)%times,timehrswave,nb,lb,H_str(iwse)%nti,np,inc)  
        tideiwse = sum(lb(1:np+1)*H_str(iwse)%wsecurv(H_str(iwse)%inc:H_str(iwse)%inc+np))
      else
        tideiwse = H_str(iwse)%wseconst  
      endif
      j = H_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = tideiwse + H_str(iwse)%wseoffset
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*H_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo ! end of each cell string

!--- Multiple Water Level BC ------------------------------------------------------
    do iwse=1,nMHstr  !for each cell string
      inc = MH_str(iwse)%inc !To protect MH_str(iwse)%inc from changing below
      call plagr_fit(MH_str(iwse)%ntimes,MH_str(iwse)%times,&
         timehrswave,nb,lb,MH_str(iwse)%nti,np,inc)
      j = MH_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = sum(lb(1:np+1)*MH_str(iwse)%wsedata(inc:inc+np,j))
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*MH_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo ! end of each cell string    
    
!--- Multiple Water Level and Velocity BC -----------------------------------------------
    do iwse=1,nMHVstr  !for each cell string
      inc = MHV_str(iwse)%incwse !To protect MHV_str(iwse)%incwse from changing below
      call plagr_fit(MHV_str(iwse)%ntimeswse,MHV_str(iwse)%timeswse,&
         timehrswave,nb,lb,MHV_str(iwse)%ntiwse,np,inc)
      j = MHV_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = sum(lb(1:np+1)*MHV_str(iwse)%wsedata(inc:inc+np,j))  
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*MHV_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo ! end of each cell string

!*********** TEMPORARY *******************************    
!--- Nested Water Level BC -----------------------------------------------
    do iwse=1,nNHstr
      j = NH_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = NHV_str(iwse)%wsebnd(j)
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*NH_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo
    
!--- Nested Water Level and Velocity BC -----------------------------------------------
    do iwse=1,nNHVstr
      j = NHV_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = NHV_str(iwse)%wsebnd(j)
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*NHV_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo    
!*********** TEMPORARY *******************************

!--- Nested Tidal Database Water Level BC (Type 9-NTH) --------------------------------    
    do iwse=1,nNTHstr  !for each cell string
      j=NTH_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = NTH_str(iwse)%wseoffset
      do k=1,NTH_str(iwse)%ntc  
       tideiwse = tideiwse + NTH_str(iwse)%amp(j,k) &
          *cos(NTH_str(iwse)%speed(k)*timehrswave + NTH_str(iwse)%phase(j,k))
      enddo
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*NTH_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo

!--- Nested Tidal Database WSE and Velocity BC (Type 10-NTHV) --------------------------------    
    do iwse=1,nNTHVstr  !for each cell string
      j=NTH_str(iwse)%ncells/2 !Use center cell as approximate value
      tideiwse = NTHV_str(iwse)%wseoffset
      do k=1,NTHV_str(iwse)%ntc            
        tideiwse = tideiwse + NTHV_str(iwse)%amp(j,k) &
           *cos(NTHV_str(iwse)%speed(k)*timehrswave + NTHV_str(iwse)%phase(j,k))
      enddo
      tide2 = tide2 + ramp*tideiwse + (1.0-ramp)*NTHV_str(iwse)%wsebnd0(j)
      ibnd = ibnd + 1
    enddo        
    
    tide2 = tide2/max(real(ibnd,kind=ikind),1.0)
    
    return
    endsubroutine tidevalue
    
!****************************************    
    subroutine freememory_fl_wav
!****************************************
    use flow_wavegrid_def
    use wave_wavegrid_def
    use wave_flowgrid_def
    implicit none
    
    !deallocate(xflwav,yflwav)
    deallocate(xwave,ywave)
    !deallocate(xwavfl,ywavfl)
    deallocate(dxwav,dywav)    
    
    return
    endsubroutine freememory_fl_wav
