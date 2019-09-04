!====================================================================================
module met_lib
! Meteorological Library
!
! Contains the following subroutines and functions:
! ~ Parameterizations ~
!     wind_drag_hsu - Wind drag coefficient based on Hsu (1988)
!     wind_drag_teeter - Wind drag coefficient based on Teeter (2001)
!     wind_drag_coeff_smith_banke - Wind drag coefficient based Smith and Banke (1975)
!     wind_drag_garratt - Wind drag coefficient based on Garratt (1977)
!     wind_drag_large_pond - Wind drag coefficient based on Large and Pond (1981)
!     wind_normrough - Calculates the normalized surface roughness = zwp z_w/h  
!     wind_heightcorr - Correction for an anemometer height other than 10 m
! ~ Input file formats ~ 
!     read_stdmetfile - Reads a NOAA Standard Meteorological file
!     read_wind_owi - Reads an Oceanweather wind file
!     read_pres_owi - Reads an Oceanweather pressure file
!     read_blended_winds - Reads a NOAA Blended Sea Winds file 
!     read_met_uvp - Reads the wind components in the x and y directions and the
!                   atmospheric pressure from meteorologic single file.
!                   The file is equivalent to the ADCIRC NWS=6 fort.22 format.
!     read_met_spddir - Reads wind speed and directions from the fleet wind file
!                   The file is equivalent to the ADCIRC NWS=3 fort.22 format.
!
! written by Alex Sanchez, USACE-CHL
!====================================================================================
    use prec_def
    implicit none

contains

!*********************************************************    
    function wind_drag_hsu(W) result(Cd)
! Wind drag coefficient based on Hsu (1988)
! and modified for high wind speeds based on data from 
! Powell et al. (2003) 
!
! Input:
!  W - 10-m wind speed [m/s]
!
! Output:
!  Cd - Wind drag coefficient
!
! References:
!   Hsu, S.A. 1988. Coastal meteorology. Academic Press, 
!      San Diego, CA.
!   Powell, M.D., Vickery, P.J., and Reinhold, T.A. (2003). 
!      Reduced drag coefficient for high wind speeds in 
!      tropical cyclones, Nature, 422, 279-283.
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************    
    use met_def, only: wkappa
    use prec_def
    implicit none
    real(ikind),intent(in) :: W
    real(ikind) :: Cd
    
    if(W<=30.0)then
      Cd = (wkappa/(14.56 - 2.0*log(max(W,0.01))))**2
    else !Reduction of drag at high speeds
      Cd = max(3.86e-3-4.0e-5*W,0.0015) !Empirical based on data by Powell
    endif    
    
    return
    endfunction wind_drag_hsu

!*********************************************************    
    function wind_drag_teeter(W,dep) result(Cd)
! Wind Drag Coefficient based on that of Teeter (2001)
!
! Input:
!  W - 10-m wind speed [m/s]
!  dep - Effective depth minimum depth of water and the
!        depth to cannopy height
!
! Output:
!  Cd - Wind drag coefficient
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: W,dep
    real(ikind) :: Cd
    
    Cd = (0.4/(16.11 - 0.5*log(max(dep,0.001)) &
       - 2.48*log(max(W,0.01))))**2
    
    return
    endfunction wind_drag_teeter    
    
!********************************************************************    
    function wind_drag_coeff_smith_banke(W,WA,WB,CdA,CdB) result(Cd)
! Wind Drag Coefficient based on Smith and Banke (1975)
! The formulation is used by Delft3D is convenient in that
! it is highly adjustable.
!
! Input:
!  W - 10-m wind speed [m/s]
!  WA - 10-m wind speed corresponding to drag coefficient CdA [m/s]
!       Range=0.0-0.1, Default=0.00063
!  WB - 10-m wind speed corresponding to drag coefficient CdB [m/s]
!       Range=0.0-0.1, Default=0.00723
!  CdA - Drag coefficient A corresponding to WA. 
!       Range=0.0-100.0 m/s, Default=0.0 m/s
!  CdB - Drag Coefficient B corresponding to WB. 
!       Range=0.0-100.0 m/s, Default=100.0 m/s
!
! Output:
!  Cd - Wind drag coefficient corresponding to wind speed W.
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************** 
    use prec_def
    implicit none
    real(ikind),intent(in) :: W,WA,WB,CdA,CdB
    real(ikind) :: Cd
    
    if(W<=WA)then
      Cd = CdA
    elseif(W<=WB)then
      Cd = CdA + (CdB-CdA)*(W-WA)/(WB-WA)
    else
      Cd = CdB  
    endif 
    
    return
    endfunction wind_drag_coeff_smith_banke    
    
!******************************************************************** 
    function wind_drag_garratt(W) result(Cd)
! Wind drag coefficient based on Garratt (1977)
!
! Description:
!   Calculates a wind drag coefficient based on Garratt (1977). 
!   A lower limit of Cd = 0.75e-3 is based on
!   Donelan et al. (2004). The upper limit of Cd = 2.64e-3 
!   is based on data from Powell et al. (2003).
!
! Usage:
!   Cd = wind_drag_large_pond(W)
!
! Input:
!   W - 10-m wind speed [m/s]
!
! Output:
!   Cd - Wind drag coefficient [-]
!
! References:
!   Donelan, M.A., Haus, B.K., Reul, N., Plant, W.J., Stiassnie, M., 
!     Graber, H.C., Brown, O.B., and Saltzman, E.S. 2004. On the limiting 
!     aerodynamic roughness of the ocean in very strong winds, Geophysical 
!     Research Letters, 31, L18306.
!   Garratt, J. R., 1977: Review of drag coefficients over oceans and
!     continents. Mon. Wea. Rev., 105, 915–929
!   Powell, M.D., Vickery, P.J., and Reinhold, T.A. 2003. Reduced drag 
!     coefficient for high wind speeds in tropical cyclones. Nature, 
!     422, 279-283.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************** 
    use prec_def
    implicit none
    real(ikind),intent(in) :: W
    real(ikind) :: Cd

    Cd = 0.001*(0.75 + 0.067*W)
    Cd = min(max(Cd,0.75e-3),2.64e-3)

    return
    endfunction wind_drag_garratt 
    
!********************************************************************     
    function wind_drag_large_pond(W) result(Cd)
! Wind drag coefficient based on Large and Pond (1981)
!
! Description:
!   Calculates a wind drag coefficient based on Large and Pond (1981). 
!   A lower limit of Cd = 0.75e-3 below W = 4 m/s is based on
!   Donelan et al. (2004). The upper limit of Cd = 2.64e-3 above W = 33 m/s 
!   is based on data from Powell et al. (2003).
!
! Usage:
!   Cd = wind_drag_large_pond(W);
!
! Input:
!   W - 10-m wind speed [m/s]
!
! Output:
!   Cd - Wind drag coefficient [-]
!
! References:
!   Donelan, M.A., Haus, B.K., Reul, N., Plant, W.J., Stiassnie, M., 
!     Graber, H.C., Brown, O.B., and Saltzman, E.S. 2004. On the limiting 
!     aerodynamic roughness of the ocean in very strong winds, Geophysical 
!     Research Letters, 31, L18306.
!   Large, W.G., and Pond, S., 1981. Open ocean flux measurements in 
!     moderate to storng winds. Journal of Physical Oceanography, 
!     11, 324-336.
!   Powell, M.D., Vickery, P.J., and Reinhold, T.A. 2003. Reduced drag 
!     coefficient for high wind speeds in tropical cyclones. Nature, 
!     422, 279-283.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************** 
    use prec_def
    implicit none
    real(ikind),intent(in) :: W
    real(ikind) :: Cd

    Cd = 0.001*(0.49 + 0.065*W)
    Cd = min(max(Cd,0.75e-3),2.64e-3)

    return
    endfunction wind_drag_large_pond
    
!*****************************************************************
    function wind_normrough(wx,wy) result(zwp)
! Calculates the normalized surface roughness = zwp z_w/h  
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: wx,wy
    real(ikind) :: zwp
    !Internal Variables
    real(ikind) :: Cd,wndspd
    
    wndspd = sqrt(wx*wx+wy*wy)
    Cd = wind_drag_hsu(wndspd)
    zwp = exp(1.0_ikind - 0.4_ikind/sqrt(Cd))
    
    return
    endfunction wind_normrough    
    
!*****************************************************************************   
    function wind_heightcorr(wndhgt,wndhgtrefin,alphain) result(fac)
! Calculates a wind velocity height correction 
!
! Description:
!  Calculates a correction factor which accounts for different wind
!  velocity measurement heights. The correction factor can be used 
!  to convert measurements between different elevations or heights. 
!  The correction formula is given by
!    w = w0*fac
!  where 
!    w0 = reference wind velocity measured at ht0
!    w = wind velocity corrected to height ht
!    fac = correction factor = (ht0/ht)**alpha
!    alpha = wind shear coefficient [-]
!
!  When a wind shear coefficient of 1/7 is used the formula is known
!  as the 1/7th rule.
!
! Input:
!  wndhgt - Wind output height [m]
!  wndhgtrefin - Reference wind height (optional, default = 10.0 m)  [m]
!  alphain - Wind shear coefficient (optional, default = 1/7) [-]
!
! Output:
!  fac - Correction factor to convert wind speed from wndhgt to wndhgtref [-]
!
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: wndhgt
    real(ikind),intent(in),optional :: wndhgtrefin,alphain
    !Internal
    real(ikind) :: wndhgtref,alpha,fac
    real(ikind), parameter :: wndhgtmin = 0.01_ikind
    
    if(present(wndhgtrefin))then
      wndhgtref = wndhgtrefin
    else
      wndhgtref = 10.0 ![m]
    endif
    if(present(alphain))then
      alpha = alphain ![-]
    else
      alpha = 1.0_ikind/7.0_ikind
    endif
    
    fac = (wndhgtref/max(wndhgt,wndhgtmin))**alpha
    
    return
    endfunction wind_heightcorr
    
!***********************************************************************
    subroutine read_stdmetfile(aname,apath,tjs,tje,nt,t,wspd,wdir,pres)
! Reads a NOAA NDBC Standard Meteorological Data File *.txt    
! The files can be downloaded from http://www.ndbc.noaa.gov/
!
!#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS PTDY  TIDE
!#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC  nmi  hPa    ft
!2007 04 15 13 50 120  4.0  6.0   0.4     3    MM  MM 1023.4  20.6  22.5  10.8   MM +1.7    MM
!
! WDIR - Wind direction (the direction the wind is coming from in degrees clockwise from true N) 
!        during the same period used for WSPD. See Wind Averaging Methods
! WSPD - Wind speed (m/s) averaged over an eight-minute period for buoys and a two-minute period 
!        for land stations. Reported Hourly. See Wind Averaging Methods.
! GST - Peak 5 or 8 second gust speed (m/s) measured during the eight-minute or two-minute period. 
!       The 5 or 8 second period can be determined by payload, See the Sensor Reporting, 
!       Sampling, and Accuracy section.
! WVHT - Significant wave height (meters) is calculated as the average of the highest one-third of 
!        all of the wave heights during the 20-minute sampling period. See the Wave Measurements section.
! DPD - Dominant wave period (seconds) is the period with the maximum wave energy. See the Wave Measurements section.
! APD - Average wave period (seconds) of all waves during the 20-minute period. See the Wave Measurements section.
! MWD - The direction from which the waves at the dominant period (DPD) are coming. The units are degrees 
!       from true North, increasing clockwise, with North as 0 (zero) degrees and East as 90 degrees. 
!       See the Wave Measurements section.
! PRES - Sea level pressure (hPa). For C-MAN sites and Great Lakes buoys, the recorded pressure is 
!        reduced to sea level using the method described in NWS Technical Procedures Bulletin 291 (11/14/80). 
!        (labeled BAR in Historical files)
! ATMP - Air temperature (Celsius). For sensor heights on buoys, see Hull Descriptions. 
!        For sensor heights at C-MAN stations, see C-MAN Sensor Locations
! WTMP - Sea surface temperature (Celsius). For sensor depth, see Hull Description.
! DEWP - Dewpoint temperature taken at the same height as the air temperature measurement.
! VIS - Station visibility (nautica miles). Note that buoy stations are limited to reports from 0 to 1.6 nmi.
! PTDY - Pressure Tendency is the direction (plus or minus) and the amount of pressure change (hPa) 
!        for a three hour period ending at the time of observation. (not in Historical files).
! TIDE - The water level in feet above or below Mean Lower Low Water (MLLW). 
!************************************************************************************
    use time_lib, only: calendar2julian,julian2calendar
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    character(len=*),intent(in) :: aname,apath !NOAA NDBC Std Met Data File (*.txt) including path
    real(ikind),intent(in) :: tjs,tje !Starting and ending dates in Julian days
    integer,intent(out) :: nt !Number of times in series
    real(ikind),intent(out),pointer :: t(:),wspd(:),wdir(:) !Reference Julian time in days
    real(ikind),intent(out),pointer,optional :: pres(:)
    !Internal
    integer :: i,jyr,ierr,nttemp,iyre,iyrs,nyr
    integer :: iyrm,imom,idaym,ihrm,iminm,isecm
    real(ikind), allocatable :: ttemp(:),wspdtemp(:),wdirtemp(:),prestemp(:)
    real(ikind) :: tji,wspdi,wdiri,wgsti,whgti,dpdi,apdi,wdmi,presi
    character(len=200) :: afile
    character(len=4) :: ayear
    logical :: foundfile,outpres
    
    ierr=0
    !Determine name(s) of files that need to be read in
    call julian2calendar(tje,iyrm,imom,idaym,ihrm,iminm,isecm)
    iyre = iyrm !End year
    call julian2calendar(tjs,iyrm,imom,idaym,ihrm,iminm,isecm)
    iyrs = iyrm !Start year
    nyr = iyre - iyrs + 1 !number of years to read
    
    !Over-allocate temporary time series arrays
    nttemp = (tje-tjs)*24 + 24
    allocate(ttemp(nttemp),wspdtemp(nttemp),wdirtemp(nttemp))
    
    outpres = .false.
    if(present(pres))then
      allocate(prestemp(nttemp))
      outpres = .true.
    endif
    
    nt = 0
    do jyr=iyrs,iyre
      write(ayear,'(I4)') jyr
      afile = trim(apath) // trim(aname) // 'h' // trim(ayear) // '.txt'  
      inquire(file=afile,exist=foundfile)
      if(.not.foundfile)then
        call diag_print_error('Could not find Met Station File: ',afile)
      endif
      !write(*,*) 'Reading Met Station File: ',trim(afile)
      open(121,file=afile) 
      read(121,*,iostat=ierr) !Skip first two lines
      read(121,*,iostat=ierr)
      do i=1,9000 !hrly records, therefore cannot be more thatn 24*365
        read(121,*,iostat=ierr) iyrm,imom,idaym,ihrm,iminm,wdiri,wspdi,&
           wgsti,whgti,dpdi,apdi,wdmi,presi
        if(ierr/=0) exit
        if(abs(wdiri-999.0)<1.0e-4 .or. abs(wspdi-99.0)<1.0e-4) cycle
        call calendar2julian(iyrm,imom,idaym,ihrm,iminm,0,tji)
        if(tji>tje+0.041667)then
          exit  
        elseif(tji>=tjs-0.041667)then
          nt = nt + 1
          ttemp(nt) = (tji - tjs)*24.0 !Hours relative to start of simulation
          ttemp(nt) = max(ttemp(nt),0.0)
          wspdtemp(nt) = wspdi
          wdirtemp(nt) = wdiri
          if(outpres) prestemp(nt) = presi
        endif          
      enddo !i-record 
      close(121)
    enddo !jyr
    
    if(nt==0)then
      call diag_print_error('No valid data found for Met Station: ',aname,&        
        ' Remove station or replace with another')
    endif
    
    !Save to output pointers
    allocate(t(nt),wspd(nt),wdir(nt))
    do i=1,nt
      t(i) = ttemp(i)
      wspd(i) = wspdtemp(i)
      wdir(i) = wdirtemp(i)
    enddo    
    if(outpres)then
      do i=1,nt  
        pres(i) = prestemp(i)
      enddo
    endif  
    
    return
    endsubroutine read_stdmetfile

!*********************************************************************************
    subroutine read_wind_owi(wunit,nwindi,nwindj,windhr2,wndspdx2,wndspdy2,ierr)
! Reads an Oceanweather wind file
!
! written by Alex Sanchez, USACE-CHL    
!*********************************************************************************
    use geo_def, only: azimuth_fl
    use comvarbl, only: tjulhr0
    use time_lib, only: julday
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: wunit,nwindi,nwindj
    integer,intent(out) :: ierr
    real(ikind),intent(out) :: windhr2
    real(ikind),intent(out) :: wndspdx2(nwindi,nwindj),wndspdy2(nwindi,nwindj)
    !Internal Variables
    integer :: i,j
    integer :: iyear,imonth,iday,ihour,imin
    logical :: fileopen
    
13  format(68x,I4,4(I2))    
    !write(*,*) ' Reading wind data'
    
    inquire(unit=wunit,opened=fileopen)
    if(.not.fileopen)then !Reached end of file  
      ierr = -4
      return 
    endif
    
    read(wunit,13,iostat=ierr) iyear,imonth,iday,ihour,imin !Read time stamp
    if(ierr/=0)then
      close(wunit)  
      call diag_print_warning('End of Oceanweather Wind File Reached')
      ierr = -1  
      return
    endif
    windhr2 = julday(iyear,imonth,iday)*24.0 + ihour + imin/60.0 - tjulhr0
    read(wunit,*,iostat=ierr) ((wndspdx2(i,j),j=1,nwindj),i=1,nwindi)
    if(ierr/=0)then
      close(wunit)  
      call diag_print_warning('Missing data in Oceanweather Wind File')
      ierr = -2
      return
    endif
    read(wunit,*,iostat=ierr) ((wndspdy2(i,j),j=1,nwindj),i=1,nwindi)
    if(ierr/=0)then
      close(wunit)  
      call diag_print_warning('Missing data in Oceanweather Wind File')
      ierr = -3
      return
    endif
    ierr = 0 !Finished successfully
    
    return
    endsubroutine read_wind_owi
    
!************************************************************************
    subroutine read_pres_owi(punit,nwindi,nwindj,preshr2,atmpres2,ierr)
! Reads an Oceanweather pressure file
!
! written by Alex Sanchez, USACE-CHL    
!************************************************************************
    use comvarbl, only: tjulhr0
    use time_lib, only: julday
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: punit,nwindi,nwindj
    integer, intent(out) :: ierr
    real(ikind), intent(out) :: preshr2
    real(ikind), intent(out) :: atmpres2(nwindi,nwindj)
    !Internal Variables
    integer :: i,j
    integer :: iyear,imonth,iday,ihour,imin
    logical :: fileopen
    
13  format(68x,I4,4(I2))
    !write(*,*) ' Reading pressure field'
    
    inquire(unit=punit,opened=fileopen)
    if(.not.fileopen)then !Reached end of file  
      ierr = -4
      return 
    endif
    
    read(punit,13,iostat=ierr) iyear,imonth,iday,ihour,imin
    if(ierr/=0)then
      close(punit)  
      call diag_print_warning('End of Oceanweather Pressure File Reached')
      ierr = -1
      return
    endif
    preshr2 = julday(iyear,imonth,iday)*24.0 + ihour + imin/60.0 - tjulhr0
    read(punit,*,iostat=ierr) ((atmpres2(i,j),j=1,nwindj),i=1,nwindi)
    if(ierr/=0)then
      close(punit)  
      call diag_print_warning('Missing data in Oceanweather Pressure File')
      ierr = -2
      return
    endif
    atmpres2 = atmpres2*100.0 !Convert from mbar to Pascal (N/m^2)      
    ierr = 0 !Finished successfully
    
    return
    endsubroutine read_pres_owi
    
!**************************************************************************
	subroutine read_blended_winds(wunit,thrs,uw,vw,ierr)
! Reads a NOAA blended winds file
!**************************************************************************
    use comvarbl, only: tjulday0
    use time_lib, only: julian2calendar
    use prec_def
    implicit none
    !Parameters
    integer, parameter :: im=1440, jm=719
    !Input/Output
    integer, intent(in) :: wunit
    integer, intent(out) :: ierr
    real(ikind),intent(in) :: thrs
    real(ikind),intent(out) :: uw(im,jm),vw(im,jm)
    !Internal Variables
    integer :: k,iyrw,imow,idayw,ihrw,iminw,isecw
    character(len=100) :: fname
    character(len=200) :: fpath,fin
    real(ikind) :: tjuldayw
        
    tjuldayw = tjulday0 + thrs/24.0
    call julian2calendar(tjuldayw,iyrw,imow,idayw,ihrw,iminw,isecw)
    
456 format(I4.4,I2.2,I2.2)
    write(fname(3:10),456) iyrw,imow,idayw	!for daily files, each contains 4 (6-hrly) records
    fin = trim(fpath) // trim(fname)
    write(*,*) ' Reading Blended Sea Winds File: ',fin
    open(wunit,file=fin,form='unformatted',status='old',convert='big_endian') !daily files,each contains 4 (6-hrly) records 
    do k=1,4
      read(wunit) uw(:,:)
	  read(wunit) vw(:,:)
      if(ihrw<=(k-1)*6) exit 
    enddo
    close(wunit)
    ierr = 0    
    
    return
    endsubroutine read_blended_winds
    
!*******************************************************************
    subroutine read_met_uvp(wunit,nwindi,nwindj,windhr2,&
                wndspdx2,wndspdy2,atmpres2,wtiminc,ierr)
! Reads the wind components in the x and y directions and the
! atmospheric pressure from meteorologic single file.
! The file format is equivalent to the ADCIRC NWS=6 format.
! written by Alex Sanchez, USACE-CHL  
!*******************************************************************
    use comvarbl, only: tjulhr0
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: wunit,nwindi,nwindj
    integer, intent(out) :: ierr
    real(ikind), intent(in) :: wtiminc
    real(ikind), intent(out) :: windhr2
    real(ikind), intent(out) :: wndspdx2(nwindi,nwindj),wndspdy2(nwindi,nwindj)
    real(ikind), intent(out) :: atmpres2(nwindi,nwindj)
    !Internal Variables
    integer :: i,j
    logical :: isopen
        
    ierr=0
    windhr2 = windhr2 + wtiminc
    
    inquire(unit=wunit,opened=isopen)
    if(.not.isopen)then !Reached end of file  
      return 
    endif
    
di: do i=1,nwindi
      do j=1,nwindj
        read(wunit,*,iostat=ierr) wndspdx2(i,j),wndspdy2(i,j),atmpres2(i,j) ![m/s], [m/s], [Pa]
        if(ierr/=0) exit di
      enddo !j 
    enddo di !i
    
    if(ierr<0)then
      close(wunit) !In case end of file is reached 
      call diag_print_warning('End of Meteorologic Single File Reached')
    elseif(ierr>0)then
      close(wunit)
      call diag_print_error('Invalid format for Meteorologic Single File')
    endif
    
    return
    endsubroutine read_met_uvp
                
!*******************************************************************
    subroutine read_met_spddir(wunit,nwindi,nwindj,windhr2,&
                wndspdx2,wndspdy2,error)
! Reads wind speed and directions from the fleet wind file
! The file is equivalent to the ADCIRC NWS=3 fort.22 format.
! The wind speed is in m/s and the wind direction is in degrees 
! coming from clockwise from true North.
! The speed and direction are converted to the wind components.
! written by Alex Sanchez, USACE-CHL  
!*******************************************************************
    use comvarbl, only: tjulhr0
    use const_def, only: deg2rad
    use time_lib, only: julday
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: wunit,nwindi,nwindj
    integer, intent(out) :: error
    real(ikind), intent(out) :: windhr2
    real(ikind), intent(out) :: wndspdx2(nwindi,nwindj),wndspdy2(nwindi,nwindj)
    !Internal Variables
    integer :: i,j
    integer :: iyear,imonth,iday,ihour
    real(ikind) :: wspd,wdir

    !Read time stamp
13  format(I4,4(I2))    
    read(wunit,13,iostat=error) iyear,imonth,iday,ihour ! = year*1000000 + month*10000 + day*100 + hr
    windhr2 = julday(iyear,imonth,iday)*24.0 + ihour - tjulhr0
        
d1: do i=1,nwindi
      do j=1,nwindj
        read(wunit,*,iostat=error) wndspdx2(i,j) !Temp
        if(error/=0) exit d1     
      enddo !j 
    enddo d1 !i   

    if(error==0)then
d2:   do i=1,nwindi
        do j=1,nwindj
          read(wunit,*,iostat=error) wdir
          if(error/=0) exit d2
          wspd=wndspdx2(i,j)
          wdir=wdir*deg2rad
          wndspdx2(i,j)=-wspd*sin(wdir)
          wndspdy2(i,j)=-wspd*cos(wdir)        
        enddo !j 
      enddo d2 !i
    endif

    if(error<0)then
      close(wunit) !In case end of file is reached 
      call diag_print_warning('End of Meteorologic Fleet File Reached')
    elseif(error>0)then
      close(wunit)
      call diag_print_error('Invalid format for Meteorologic Fleet File')
    endif
    
    return
    endsubroutine read_met_spddir    

endmodule met_lib
