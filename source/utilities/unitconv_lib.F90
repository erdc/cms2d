!=====================================================================
module unitconv_lib
! Unit Conversion library
!
! Calculates the variables for a general linear conversion of the form: 
!    Y = X*fac + con
!  where
!    X = Input variable
!    Y = Output variable
!    fac = conversion factor (coefficient)
!    con = conversion constant
!
!   Supports converting to and from the following units:
!     Angles - rad, deg
!     Time - sec, min, hrs, days, weeks, years    
!     Length - um, mm, cm, m, km, in, ft, yd
!     Velocity - mm/s, cm/s, m/s, km/hr, ft/s, knot, mph,
!                m/min, m/hr, ft/min, ft/hr
!     Temperature - celcius, fahrenheit, kelvin, rankine
!     Mass - mg, cm, g, kg, oz, lb, tonUS
!     Density - g/cm^3, g/m^3, kg/m^3, oz/gal, oz/in^3, 
!               lb/in^3, lb/ft^3, lb/yd^3, ton/yd^3
!     Pressure - Pa, hPa, kPa, kg/cm^2, kg/m^2, mb, 
!                bar, mmHg, inHg, psf, psi
!     Salinity - pph, ppt, ppm, w 
!     Volume Flux - m^3/s, m^3/hr, m^3/day, 
!                  ft^3/s, ft^3/hr, ft^3/day,
!                  yd^3/s, yd^2/hr, yd^3/day
!     Volume - mm^3, cm^3, m^3, km^3, in^3, ft^3, yd^3, 
!             mi^3, nmi^3, gal, ml, l
!
!   Units may be input and a variety of formats. 
!   For a list of supported formats see tables below.
!
! Contains:
!   unitconv_scal - Scalar unit conversion
!   unitconv_vec  - Vector array unit conversion
!   unitconv_var  - Calculates unit conversion factor and constant
!   unitconv_fac  - Calculates unit conversion factor
!
! Author: Alex Sanchez, USACE-CHL
!=====================================================================
    use prec_def
    implicit none
        
contains

!*********************************************************************    
    subroutine unitconv_scal(fromunits,tounits,scalar)
! Scalar unit conversion
!
! Description:
!   Computes the the general linear conversion of the form: 
!      v = v*fac + con
!   where
!      v = Input vector array
!      fac = conversion factor (coefficient)
!      con = conversion constant
!  
! Input:
!   tounits = Input unit string
!   fromunits = Output unit strings
!   n = Size of vector array fromunits
!   v = Vector array with units tounits
!
! Output:
!   v = Vector array with units tounits
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
    implicit none
    !Input/Output
    character(len=*),intent(in) :: fromunits,tounits
    real(ikind),     intent(inout) :: scalar
    !Internal Variables
    real(ikind) :: fac,con
    
    call unitconv_var(fromunits,tounits,fac,con)
    scalar = scalar*fac + con !General 
    
    return
    end subroutine unitconv_scal
    
!*********************************************************************    
    subroutine unitconv_vec(fromunits,tounits,n,v)
! Vector array unit conversion
!
! Description:
!   Computes the the general linear conversion of the form: 
!      v = v*fac + con
!   where
!      v = Input vector array
!      fac = conversion factor (coefficient)
!      con = conversion constant
!  
! Input:
!   tounits = Input unit string
!   fromunits = Output unit strings
!   n = Size of vector array fromunits
!   v = Vector array with units tounits
!
! Output:
!   v = Vector array with units tounits
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
    implicit none
    !Input/Output
    character(len=*),intent(in) :: fromunits,tounits
    integer,         intent(in) :: n
    real(ikind),     intent(inout) :: v(n)
    !Internal Variables
    integer :: i
    real(ikind) :: fac,con
    
    call unitconv_var(fromunits,tounits,fac,con)

!$OMP PARALLEL DO PRIVATE(i) IF(n>5000)
    do i=1,n
      v(i) = v(i)*fac + con
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine unitconv_vec
    
!*********************************************************************    
    function unitconv_fac(fromunits,tounits) result(fac)
! Calculates the unit conversion factor for the general conversion
!   Y = X*fac
!
! Input:
!   tounits = Input unit string
!   fromunits = Output unit strings
!
! Output:
!   fac = factor
!
! Author: Alex Sanchez, USACE-CHL
!*********************************************************************    
    use prec_def
    implicit none
    !Input/Output
    character(len=*),intent(in) :: fromunits,tounits
    real(ikind) :: fac
    
    call unitconv_var(fromunits,tounits,fac)
    
    return
    end function unitconv_fac
    
!*********************************************************************    
    subroutine unitconv_var(fromunits,tounits,fac,con)
! Calculates the unit conversion variables
!
! Description:
!   Computes the variables for a general linear conversion of the form: 
!      Y = X*fac + con
!   where
!      X = Input variable
!      Y = Output variable
!      fac = conversion factor (coefficient)
!      con = conversion constant
!  
! Input:
!   tounits = Input unit string
!   fromunits = Output unit strings
!
! Output:
!   fac = factor
!   con = constant (optional)
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
    !use diag_lib
    implicit none
    !Input/Output
    character(len=*),intent(in) :: fromunits,tounits
    real(ikind),intent(out) :: fac
    real(ikind),intent(out),optional :: con
    !Internal variables
    integer :: ifromtype,itotype,itounits,ifromunits,ierr
    integer, parameter :: nunits = 12
    real(ikind) :: factable(nunits,nunits),contable(nunits,nunits)
    
    call unitcodes(fromunits,ifromunits,ifromtype,ierr)
    if(ierr==-1)then
      write(*,*) 'Unit conversion error: Unknown "from" unit type: ',trim(fromunits)  
      read(*,*)
      stop
    endif
    call unitcodes(tounits,itounits,itotype,ierr)
    if(ierr==-1)then
      write(*,*) 'Unit conversion error: Unknown "to" unit type: ',trim(tounits)
      read(*,*)
      stop
    endif
    
    !Check unit types
    if(ifromtype/=itotype)then
      write(*,*) 'Unit conversion error: "To" and "from" units do not match '
      write(*,*)  '   "To" Units: ',trim(tounits)
      write(*,*) '   "From" Units: ',trim(fromunits)
      read(*,*)
      stop
    endif
    
    !Early exit for dimensionless variables
    if(itounits==0 .or. ifromunits==0)then
      fac = 1.0  
      if(present(con)) con = 0.0
      return !Dimensionless
    endif  
    
    !Creat Unit Conversion Table
    factable = 0.0 !Initialize
    contable = 0.0 !Initialize
    select case(ifromtype)    
    case(1) !Angles
    !from               rad                deg
    factable(1,1:2) = (/1.0               ,0.017453292519943/) !to rad
    factable(2,1:2) = (/57.295779513082323,1.0/)               !to deg
    
    case(2) !Time    
    !from              sec           , min       , hrs       , days      , weeks     , years    
    factable(1,1:6) = (/1.0          ,60.0       ,3600.0     ,86400.0    ,604800.0   ,31536000.0 /) !to sec
    factable(2,1:6) = (/0.01666666667,1.0        ,60.0       ,1440.0     ,10080.0    ,525600.0   /) !to min
    factable(3,1:6) = (/0.00027777778,0.016666667,1.0        ,24.0       ,168.0      ,8760.0     /) !to hrs
    factable(4,1:6) = (/1.15741E-05  ,0.000694444,0.041666667,1.0        ,7.0        ,365.0      /) !to days
    factable(5,1:6) = (/1.65344E-06  ,9.92063E-05,0.005952381,0.142857143,1.0        ,52.14285714/) !to weeks
    factable(6,1:6) = (/3.17098E-08  ,1.90259E-06,0.000114155,0.002739726,0.019178082,1.0        /) !to years
    
    case(3) !Length
    !from               um           ,  mm         , cm          , m         , km          , in       , ft     ,yd      , mile  , nautical mile
    factable(1 ,1:10) = (/       1.0,        1000.0,       10000.0,       1.0e+06,       1.0e+09,       25400.0,      304800.0,      914400.0,    1.60934e+09,  1.852e+09    /)  !to um
    factable(2 ,1:10) = (/     0.001,          1.0,       10.0,        1000.0,       1.0e+06,          25.4,         304.8,         914.4,    1.60934e+06,  1.852e+06    /)  !to mm
    factable(3 ,1:10) = (/    0.0001,           0.1,           1.0,         100.0,       1.0E+05,          2.54,         30.48,         91.44,      160934.0,      185200.0    /)  !to cm
    factable(4 ,1:10) = (/   1.0e-06,         0.001,          0.01,           1.0,        1000.0,        0.0254,        0.3048,        0.9144,       1609.34,        1852.0    /)  !to m
    factable(5 ,1:10) = (/   1.0e-09,       1.0e-06,       1.0e-05,         0.001,           1.0,      2.54e-05,     0.0003048,     0.0009144,       1.60934,         1.852    /)  !km
    factable(6 ,1:10) = (/3.93701e-05,     0.0393701,      0.393701,       39.3701,       39370.1,           1.0,          12.0,          36.0,       63360.0,       72913.4    /)  !in
    factable(7 ,1:10) = (/3.28084e-06,    0.00328084,     0.0328084,       3.28084,       3280.84,     0.0833333,           1.0,           3.0,        5280.0,       6076.12    /)  !ft
    factable(8 ,1:10) = (/1.09361e-06,    0.00109361,     0.0109361,       1.09361,       1093.61,     0.0277778,      0.333333,           1.0,        1760.0,       2025.37    /)  !yd
    factable(9 ,1:10) = (/6.21371e-10,  6.21371e-07, 6.21371e-06, 0.000621371, 0.621371, 1.57828e-05, 0.000189394, 0.000568182,    1.0,       1.15078    /)  !mi
    factable(10,1:10) = (/5.39957e-10,  5.39957e-07, 5.39957e-06, 0.000539957, 0.539957, 1.37149e-05, 0.000164579, 0.000493737, 0.868976,      1.0    /)  !nmi
    
    case(4) !Velocity
    !from               mm/s       , cm/s      , m/s        , km/hr     , ft/s      ,knot       ,mph        , m/min    , m/hr      , ft/min    , ft/hr
    factable(1 ,1:11) = (/1.0        ,10.0       ,1.0E+03    ,277.7777778,304.80     ,514.4444   ,447.04     ,16.6666667,0.277777778,5.08       ,0.084666667/) !to mm/s
    factable(2 ,1:11) = (/0.1        ,1.0        ,1.0E+02    ,27.77778   ,30.480     ,51.44444   ,4.4704     ,1.6666667 ,0.027777778,0.508      ,0.008466667/) !to cm/s
    factable(3 ,1:11) = (/1.0E-03    ,0.01       ,1.0        ,0.277777778,0.3048     ,0.514444444,0.44704    ,0.0166667 ,0.002777778,0.0508     ,0.000846667/) !to m/s
    factable(4 ,1:11) = (/3.60E-03   ,3.60E-02   ,3.60       ,1.0        ,1.097279968,1.852      ,1.609344   ,0.06      ,0.001      ,0.018288   ,0.0003048  /) !to km/hr
    factable(5 ,1:11) = (/3.28084E-03,0.0328084  ,3.28084    ,0.911344442,1.0        ,1.68781    ,1.466666666,0.0546807 ,0.000911344,0.016666667,0.000277778/) !to ft/s
    factable(6 ,1:11) = (/1.94384E-03,1.94384E-02,1.94384    ,0.539956803,0.592483801,1.0        ,0.868976242,0.0323974 ,0.000539957,0.00987473 ,0.000164579/) !to knot
    factable(7 ,1:11) = (/2.23694E-03,2.23694E-02,2.23694    ,0.621371192,0.681818182,1.150779448,1.0        ,0.0372823 ,0.000621371,0.011363636,0.000189394/) !to mph
    factable(8 ,1:11) = (/0.06        ,0.6        ,60.0       ,16.66666667,18.288     ,30.86666664,26.8224    ,1.0       ,0.016666667,3.280839895,0.054680665/) !to m/min
    factable(9 ,1:11) = (/3.6        ,36.0       ,3600.00    ,1000.0     ,1097.28    ,1852.00    ,1609.344   ,60.0      ,1.0        ,196.8503937,3.280839895/) !to m/hr
    factable(10,1:11) = (/0.196850394,1.9685039  ,196.8503937,54.68066492,60.0       ,101.2685914,88.0       ,0.3048    ,0.00508    ,1.0        ,0.0166667  /) !to ft/min
    factable(11,1:11) = (/11.81102362,118.1102362,11811.02362,3280.839895,3600.0     ,6076.115486,5280.0     ,18.288    ,0.3048     ,60.0       ,1.0        /) !to ft/hr
    
    case(5) !Temperature
    !from               celcius,fahrenheit, kelvin  , rankine
    factable(1,1:4) = (/1.0    ,0.5555555 ,1.0      ,0.5555556/) !to celcius
    factable(2,1:4) = (/1.8    ,1.0       ,0.5555556,1.0      /) !to fahrenheit
    factable(3,1:4) = (/1.0    ,1.8       ,1.0      ,1.8      /) !to kelvin
    factable(4,1:4) = (/1.8    ,1.0       ,1.0      ,1.0      /) !to rankine  
    
    contable(1,1:4) = (/0.0    ,-17.777778,0.0      ,0.0      /) !to celcius
    contable(2,1:4) = (/57.6   ,0.0       ,255.37222,459.67   /) !to fahrenheit
    contable(3,1:4) = (/-273.15,-459.67   ,0.0      ,0.0      /) !to kelvin
    contable(4,1:4) = (/-273.15,-459.67   ,0.0      ,0.0      /) !to rankine    
    
    case(6) !Mass
    !from               mg           , cm           , g            , kg          , oz           , lb           , tonUK        , tonUS         , tonne
    factable(1,1:9) = (/1.0          , 10.0         , 1000.0       , 1000000.0   , 2.8349523E+04, 4.5359237E+05, 1.0160469E+09, 907184740.0003, 1.0E+09/) !to mg
    factable(2,1:9) = (/0.1          , 1.0          , 100.0        , 10000.0     , 2.8349523E+02, 4.5359237E+03, 1.0160469E+07, 90718474.0    ,    1.0E+08/) !to cg
    factable(3,1:9) = (/0.001        , 0.010        , 1.0          , 1000.0      , 2.8349523E+01, 4.5359237E+02, 1.0160469E+06, 907184.74     ,    1.0E+06/) !to g
    factable(4,1:9) = (/0.000001     , 0.00001      , 0.001        , 1.0         , 0.028349523  , 0.45359237   , 1016.046909  , 907.18474     ,    1.0E+03/) !to kg
    factable(5,1:9) = (/3.5273962E-05, 3.5273962E-04, 3.5273962E-02, 35.27396195 , 1.0          , 16.0         , 35840.0      , 3.2E+04         , 3.5273962E+04/) !to oz
    factable(6,1:9) = (/2.2046226E-06, 2.2046226E-05, 2.2046226E-03, 2.20462262  , 0.0625       , 1.0          , 2240.0       , 2.0E+03          , 2204.622622/) !to lb
    factable(7,1:9) = (/9.8420653E-10, 9.8420653E-09, 9.8420653E-07, 9.842065E-04, 2.7901786E-05, 4.4642857E-04, 1.0          , 0.89286          , 0.984206528/) !to ton UK (long)
    factable(8,1:9) = (/1.1023113E-09, 1.1023113E-08,    1.1023113E-06, 1.102311E-03, 3.1250000E-05, 5.0000000E-04, 1.12         , 1.0           ,    1.102311311/) !to ton US (short)
    factable(9,1:9) = (/1.0000000E-09, 1.0000000E-08, 1.0000000E-06, 1.000000E-03, 2.8349523E-05, 4.5359237E-04, 1.0160469088, 0.90718474     , 1.0/)         !to tonne (SI units)

    case(7) !Density
    !from               g/cm^3      , g/m^3      , kg/m^3      , oz/gal      , oz/in^3    , lb/in^3  , lb/ft^3, lb/yd^3      , ton/yd^3
    factable(1,1:9) = (/1.0         ,1.0E-06     , 0.001       , 0.00748915  , 1.72999    , 27.68    , 1.602E-02, 5.93E-04   ,0.693592385 /) !to g/cm^3
    factable(2,1:9) = (/1.0E+06     ,1.0         , 1000.0      , 7489.15     , 1.72999E+06, 2.768E+07, 1.602E+04, 5.93276E+02, 6.9359E+05 /) !to g/m^3
    factable(3,1:9) = (/1.0E+03     ,1.00E-03    , 1.0         , 7.48915     , 1.72999E+03, 2.768E+04, 1.602E+01, 5.93276E-01, 6.9359E+11 /) !to kg/m^3
    factable(4,1:9) = (/133.5265    ,1.335265E-04, 0.1335265   , 1.0         , 231.0      , 3696.0   , 2.14     , 7.92181E-02, 6.9359E+17 /) !to oz/gal
    factable(5,1:9) = (/0.5780368   ,5.780368E-07, 0.0005780   , 0.0043290   , 1.0        , 16.0     , 0.0093   , 0.0003429  , 0.7680492  /) !to oz/in^3
    factable(6,1:9) = (/3.612729E-02,3.612729E-08, 3.612729E-05, 2.705628E-04, 0.0625     , 1.0      , 5.787E-04, 2.14335E-05, 4.80307E-02/) !to lb/in^3
    factable(7,1:9) = (/62.42796    ,6.242796E-05, 0.06242796  , 0.4675325   , 108.0      , 1728.0   , 1.0      , 0.037      , 82.96      /) !to lb/ft^3
    factable(8,1:9) = (/1.685555E+03,1.685555E-03, 1.685555E+00, 1.262338E+01, 2.916E+03  , 46656.0  , 27.0     , 1.0        , 2240.0     /) !to lb/yd^3
    factable(9,1:9) = (/1.441769    ,1.441769E-06, 1.441769E-12, 1.441769E-18, 1.302      , 20.82    , 0.0121   , 4.46429E-04, 1.0        /) !to ton/yd^3
    
    case(8) !Pressure
    !from                 Pa         , hPa         , kPa         , kg/cm^2    , kg/m^2     , mb         , bar        , mmHg       , inHg       , psf        , psi
    factable(1 ,1:11) = (/1.0        , 100.0       , 1000.0      , 9.80665E+04, 9.80665    , 100.0      , 1.0E+05    , 1.33322E+02, 3.38639E+03, 6.89476E+03, 4.78803E+01/) !to Pa
    factable(2 ,1:11) = (/0.01       , 1.0         , 10.0        , 9.80665E+03, 9.80665E-01, 10.0       , 1.0E+04    , 1.33322E+01, 3.38639E+02, 6.89476E+02, 4.78803    /) !to hPa
    factable(3 ,1:11) = (/0.001      , 0.1         , 1.0         , 9.80665E+02, 9.80665E-02, 1.0        , 1.0E+03    , 1.33322    , 3.38639E+01, 6.89476E+01, 4.78803E-01/) !to kPa
    factable(4 ,1:11) = (/1.01972E-05, 1.01972E-03 , 1.01972E-02 , 1.0        , 1.0E-04    , 1.01972E-03, 1.01972    , 1.35951E-03, 3.45316E-02, 7.0307E-02 , 4.88242E-04/) !to kg/cm^2
    factable(5 ,1:11) = (/0.10197    , 1.01972E+01 , 1.019716E+02, 1.0E+04    , 1.0        , 1.01972E+01, 1.01972E+04, 13.5950981 , 345.3155   , 7.0307E+02 , 4.88242    /) !to kg/m^2
    factable(6 ,1:11) = (/0.01       , 1.0         , 10.0        , 9.80665E+02, 9.80665E-02, 1.0        , 1.00000E+03, 1.33322    , 33.86388   , 6.89476E+01, 4.78802E-01/) !to mb
    factable(7 ,1:11) = (/1.0E-05    , 1.0E-03     , 0.01        , 9.80665E-01, 9.80665E-05, 0.001      , 1.0        , 1.33322E-03, 0.03386388 , 6.89476E-02, 4.78802E-04/) !to bar
    factable(8 ,1:11) = (/7.50062E-03, 7.50062E-01 , 7.50062     , 7.35559E+02, 7.35559E-02, 7.50062E-01, 7.50062E+02, 1.0        , 25.4       , 51.71493   , 0.359131096/) !to mmHg
    factable(9 ,1:11) = (/2.95300E-04, 2.952999E-02, 2.952999E-01, 28.95903   , 2.89590E-03, 2.953E-02  , 2.95300E+01, 3.93701E-02, 1.0        , 2.036021   , 0.014139034/) !to inHg
    factable(10,1:11) = (/1.45038E-04, 1.450377E-02, 1.450377E-01, 14.22336   , 1.42234E-03, 1.45038E-02, 1.45038E+01, 1.93368E-02, 0.4911541  , 1.0        , 0.006944444/) !to psf
    factable(11,1:11) = (/2.08854E-02, 2.088543E+00, 2.088543E+01, 2.04816E+03, 2.04816E-01, 2.08855    , 2.08855E+03, 2.7845     , 70.72619   , 144.0      , 1.0        /) !to psi

    case(9) !Salinity
    !from               pph      ppt      ppm      kg/kg
    factable(1,1:4) = (/1.0       , 0.1    , 1.0E-04, 100.0  /) !to pph 
    factable(2,1:4) = (/10.0   , 1.0    , 1.0E-03, 1.0E+03/) !to ppt
    factable(3,1:4) = (/1.0E+04, 1.0E+03, 1.0    , 1.0E+06/) !to ppm
    factable(4,1:4) = (/1.0E-02, 1.0E-03, 1.0E-06, 1.0    /) !to kg/kg
    
    case(10) !Volume flux (L^3/T)
    !From                   m^3/s           m^3/hr        m^3/day          ft^3/s        ft^3/hr          ft^3/day          L/s           L/hr             L/day           yd^3/s        yd^3/hr          yd^3/day
    factable(1 ,1:12) = (/1.000000E+00, 2.777780E-04, 1.157400E-05,    2.831685E-02, 7.865791E-06,    3.277413E-07, 1.000000E-03,    3.600000E+00, 8.640000E+01,    7.645549E-01,2.752397E+03,    6.605754E+04/) !to m^3/s
    factable(2 ,1:12) = (/3.600000E+03,    1.000000E+00, 4.166667E-02,    1.019406E+02, 2.831685E-02,    1.179869E-03, 3.600000E+00,    1.296000E+04, 3.110400E+05,    2.752397E+03, 9.908631E+06,    2.378071E+08/) !to m^3/hr
    factable(3 ,1:12) = (/8.640000E+04,    2.400000E+01, 1.000000E+00,    2.446576E+03, 6.796043E-01,    2.831685E-02, 8.640000E+01,    3.110400E+05, 7.464960E+06,    6.605754E+04, 2.378071E+08,    5.707371E+09/) !to m^3/day
    factable(4 ,1:12) = (/3.531467E+01,    9.809630E-03, 4.087350E-04,    1.000000E+00, 2.777778E-04,    1.157407E-05, 3.531467E-02,    1.271328E+02, 3.051187E+03,    2.700000E+01, 9.720000E+04,    2.332800E+06/) !to ft^3/s
    factable(5 ,1:12) = (/1.271328E+05,    3.531467E+01, 1.471444E+00,    3.600000E+04, 1.000000E+00,    4.166667E-02, 1.271328E+02,    4.576781E+05, 1.098427E+07,    9.720000E+04, 3.499200E+08,    8.398080E+09/) !to ft^3/hr
    factable(6 ,1:12) = (/3.051187E+06,    8.475520E+02, 3.531467E+01,    8.640000E+05, 2.400000E+01,    1.000000E+00, 3.051187E+03,    1.098427E+07, 2.636226E+08,    2.332800E+06, 8.398080E+09,    2.015539E+11/) !to ft^3/day
    factable(7 ,1:12) = (/1.000000E+03,    2.777778E-01, 1.157407E-02,    2.831685E+01, 7.865791E-03,    3.277413E-04, 1.000000E+00,    3.600000E+03, 8.640000E+04,    7.645549E+02, 2.752397E+06,    6.605754E+07/) !to L/s
    factable(8 ,1:12) = (/3.600000E+06,    1.000000E+03, 4.166667E+01,    1.019406E+05, 2.831685E+01,    1.179869E+00, 3.600000E+03,    1.000000E+00, 3.110400E+08,    2.752397E+06, 9.908631E+09,    2.378071E+11/) !to L/hr
    factable(9 ,1:12) = (/8.640000E+07,    2.400000E+04, 1.000000E+03,    2.446576E+06, 6.796043E+02,    2.831685E+01, 8.640000E+04,    2.400000E+01, 1.000000E+00,    6.605754E+07, 2.378071E+11,    5.707371E+12/) !to L/day
    factable(10,1:12) = (/1.307923E+00,    3.633120E-04, 1.513800E-05,    3.703704E-02, 1.028800E-05,    4.286667E-07, 1.307951E-03,    4.708624E+00, 1.130070E+02,    1.000000E+00, 3.600000E+03,    8.640000E+04/) !to yd^3/s
    factable(11,1:12) = (/4.708622E+03,    1.307951E+00, 5.449794E-02,    1.333333E+02, 3.703680E-02,    1.543200E-03, 4.708624E+00,    1.695104E+04, 4.068251E+05,    3.600000E+03, 1.000000E+00,    3.110400E+08/) !to yd^3/hr
    factable(12,1:12) = (/1.130069E+05,    3.139081E+01, 1.307951E+00,    3.200000E+03, 8.888832E-01,    3.703680E-02, 1.130070E+02,    4.068251E+05, 9.763802E+06,    8.640000E+04, 2.400000E+01,    1.000000E+00/) !to yd^3/day
        
    case(11) !Volume
    !From                   mm^3           cm^3              m^3            km^3          in^3            ft^3          yd^3            mi^3          nmi^3            gal               ml              l
    factable(1 ,1:12) = (/1.000000E+00, 1.000000E+03, 1.000000E+09, 1.000000E+18, 1.638710E+04, 2.831680E+07, 7.645550E+08, 4.168180E+18, 6.352182E+18, 3.785410E+06, 1.000000E+03, 1.000000E+06/) !to mm^3
    factable(2 ,1:12) = (/1.000000E-03, 1.000000E+00, 1.000000E+06, 1.000000E+15, 1.638710E+01, 2.831680E+04, 7.645550E+05, 4.168180E+15, 6.352182E+15, 3.785410E+03, 1.000000E+00, 1.000000E+03/) !to cm^3
    factable(3 ,1:12) = (/1.000000E-09,    1.000000E-06, 1.000000E+00, 1.000000E+09, 1.638710E-05, 2.831680E-02, 7.645550E-01, 4.168180E+09, 6.352182E+09, 3.785410E-03, 1.000000E-06, 1.000000E-03/) !to m^3
    factable(4 ,1:12) = (/1.000000E-18,    1.000000E-15, 1.000000E-09, 1.000000E+00, 1.638710E-14, 2.831680E-11, 7.645550E-10, 4.168180E+00, 6.352182E+00, 3.785410E-12, 1.000000E-15, 1.000000E-12/) !to km^3
    factable(5 ,1:12) = (/6.102361E-05,    6.102361E-02, 6.102361E+04, 6.102361E+13, 1.000000E+00, 1.727993E+03, 4.665591E+04, 2.543574E+14, 3.876331E+14, 2.309994E+02, 6.102361E-05, 6.102361E-02/) !to in^3
    factable(6 ,1:12) = (/3.531472E-08,    3.531472E-05, 3.531472E+01, 3.531472E+10, 5.787059E-04, 1.000000E+00, 2.700005E+01, 1.471981E+11, 2.243256E+11, 1.336807E-01, 3.531472E-08, 3.531472E-05/) !to ft^3
    factable(7 ,1:12) = (/1.307950E-09,    1.307950E-06, 1.307950E+00, 1.307950E+09, 2.143351E-05, 3.703697E-02, 1.000000E+00, 5.451773E+09, 8.308339E+09, 4.951128E-03, 1.307950E-09, 1.307950E-06/) !to yd^3
    factable(8 ,1:12) = (/2.399129E-19,    2.399129E-16, 2.399129E-10, 2.399129E-01, 3.931476E-15, 6.793565E-12, 1.834266E-10, 1.000000E+00, 1.523970E+00, 9.081686E-13, 2.399129E-19, 2.399129E-16/) !to mi^3
    factable(9 ,1:12) = (/1.574262E-19,    1.574262E-16, 1.574262E-10, 1.574262E-01, 2.579759E-15, 4.457807E-12, 1.203610E-10, 6.561808E-01, 1.000000E+00, 5.959228E-13, 1.574262E-19, 1.574262E-16/) !to nmi^3
    factable(10,1:12) = (/2.641722E-07,    2.641722E-04, 2.641722E+02, 2.641722E+11, 4.329016E-03, 7.480511E+00, 2.019742E+02, 1.101117E+12, 1.678070E+12, 1.000000E+00, 2.641722E-07, 2.641722E-04/) !to gal
    factable(11,1:12) = (/1.000000E-03,    1.000000E+00, 1.000000E+06, 1.000000E+15, 1.638710E+01, 2.831680E+04, 7.645550E+05, 4.168180E+15, 6.352181E+15, 3.785410E+06, 1.000000E+00, 1.000000E+03/) !to ml
    factable(12,1:12) = (/1.000000E-06,    1.000000E-03, 1.000000E+03, 1.000000E+12, 1.638710E-02, 2.831680E+01, 7.645550E+02, 4.168180E+12, 6.352181E+12, 3.785410E+03, 1.000000E-03, 1.000000E+00/) !to l
    
    case(12) !Acceleration
    !From               m/s^2    ft/s^2
    factable(1,1:2) = (/1.0,    0.3048/) !m/s^2
    factable(2,1:2) = (/3.28084, 1.0/) !ft/s^2
    
    case(13) !Diffusivity
    !From               m^2/s    ft^2/s
    factable(1,1:2) = (/1.0       , 0.09290304/) !m^2/s
    factable(2,1:2) = (/10.7639103, 1.0       /) !ft^2/s
    
    case(14) !Specific weight
    !From              kg/m^2/s^2      lb/ft^2/s^2
    factable(1,1:2) = (/1.0        , 4.882427636/) !kg/m^2/s^2
    factable(2,1:2) = (/0.204816144, 1.0        /) !lb/ft^2/s^2
      
    end select
    
    fac = factable(itounits,ifromunits)
    if(present(con)) con = contable(itounits,ifromunits)
    
    return
    contains
!=============================================================
    subroutine unitcodes(aunits,iunitid,iunittype,ierr)
! For a unit string outputs the unit ID and type codes.
! written by Alex Sanchez, USACE-CHL    
!=============================================================
    implicit none
    !Input/Output    
    character(len=*),intent(in) :: aunits
    integer,         intent(out) :: iunitid,iunittype,ierr
    !Internal Variables
    character(len=len(aunits)) :: aunitstemp
    
    ierr = 0
    aunitstemp = aunits
    call lowercase(aunitstemp)
    select case(aunitstemp)
    !--- None ----------------------
    case(' ','!','#','-')
      iunittype = 0; iunitid = 0
      
    !--- Angles --------------------
    case('rad','radian','radians')
      iunittype = 1; iunitid = 1
    case('deg','º','degree','degrees')
      iunittype = 1; iunitid = 2
      
    !--- Time Units --------------------  
    case('s','sec','second','seconds')
      iunittype = 2; iunitid = 1
    case('min','minute','minutes')
      iunittype = 2; iunitid = 2
    case('h','hr','hrs','hour','hours')
      iunittype = 2; iunitid = 3
    case('d','dy','day','days')
      iunittype = 2; iunitid = 4
    case('w','wk','week','weeks','wks')
      iunittype = 2; iunitid = 5
    case('y','yr','yrs','year','years')
      iunittype = 2; iunitid = 6
      
    !---Length Units ---------------------
    case('um','micron','microns')
      iunittype = 3; iunitid = 1
    case('mm','millimeter','millimeters')
      iunittype = 3; iunitid = 2
    case('cm','centimeter','centimeters')
      iunittype = 3; iunitid = 3
    case('m','meter','meters')
      iunittype = 3; iunitid = 4
    case('km','kilometer','kilometers')
      iunittype = 3; iunitid = 5
    case('in','inch','inches')
      iunittype = 3; iunitid = 6
    case('ft','feet','foot')
      iunittype = 3; iunitid = 7
    case('yd','yds','yard','yards')
      iunittype = 3; iunitid = 8
    case('sm','smi','mile','mi','statutemile','statmi')
      iunittype = 3; iunitid = 9
    case('nm','nmi','nauticalmile','nautmi')
      iunittype = 3; iunitid = 10
      
    !--- Velocity Units ------------------------------------------
    case('mmps','mm/s','mm/sec','millimeter/second','millimeters/second')
      iunittype = 4; iunitid = 1
    case('cmps','cm/s','cm/sec','centimeter/second','centimeters/second')
      iunittype = 4; iunitid = 2
    case('mps','m/s','m/sec','meter/second','meters/second')
      iunittype = 4; iunitid = 3
    case('kph','km/hr','km/h','kilometer/hour','kilometers/hour')
      iunittype = 4; iunitid = 4
    case('fps','ftps','ft/s','ft/sec','feet/sec','feet/second','foot/second')
      iunittype = 4; iunitid = 5
    case('knot','knots')
      iunittype = 4; iunitid = 6
    case('mph','mile/hr','mi/hr','miles/hour','mile/hour')
      iunittype = 4; iunitid = 7
    case('m/min','meter/minute','meters/minute')
      iunittype = 4; iunitid = 8  
    case('m/hr','m/hrs','meter/hour','meters/hour','meters/hours')
      iunittype = 4; iunitid = 9    
    case('ft/min','fpm','feet/min','foot/min')
      iunittype = 4; iunitid = 10
    case('ft/hr','ft/hour','feet/hr','feet/hour')
      iunittype = 4; iunitid = 11
      
    !--- Temperature Units ---------------------------
    case('ºc','c','celcius','cel','centigrade','cen','deg c')
      iunittype = 5; iunitid = 1
    case('ºf','f','fahrenheit','fah','deg f')
      iunittype = 5; iunitid = 2
    case('ºk','k','kelvin','kel')
      iunittype = 5; iunitid = 3
    case('ºr','r','rankine','ran')
      iunittype = 5; iunitid = 4
    
    !--- Mass -----------------------
    case('mg','milligram')
      iunittype = 6; iunitid = 1      
    case('cg','centigram')
      iunittype = 6; iunitid = 2
    case('g','gram')
      iunittype = 6; iunitid = 3
    case('kg','kilogram')
      iunittype = 6; iunitid = 4
    case('oz','once')
      iunittype = 6; iunitid = 5
    case('lb','pound')
      iunittype = 6; iunitid = 6
     case('ton','tonUS')
      iunittype = 6; iunitid = 7
    
    !--- Density ------------------
    case('g/cm^3','g/cm3','gram/cm^3','gram/centimeter^3')
      iunittype = 7; iunitid = 1
    case('g/m^3','g/m3','gram/meter^3')
      iunittype = 7; iunitid = 2
    case('kg/m^3','kg/m3','kilogram/meter^3')
      iunittype = 7; iunitid = 3
    case('oz/gal')
      iunittype = 7; iunitid = 4
    case('oz/in^3','oz/in3','oz/cin')
      iunittype = 7; iunitid = 5
    case('lb/in^3','lbs/in^3','lb/in3','lbs/in3','lb/cin','lbs/cin')
      iunittype = 7; iunitid = 6
    case('lb/ft^3','lb/ft3','lb/cft','lb/cf','lbs/ft^3','lbs/ft3','lbs/cft','lbs/cf')
      iunittype = 7; iunitid = 7
    case('lb/yd^3','lb/yd3','lb/cyd','lbs/yd^3','lbs/yd3','lbs/cyd')
      iunittype = 8; iunitid = 6
    case('ton/yd^3','ton/yd3','ton/cyd','tons/yd^3','tons/yd3','tons/cyd')
      iunittype = 7; iunitid = 9
    
    !--- Pressure -------------------------
    case('pa','pascal','pascals')
      iunittype = 8; iunitid = 1      
    case('hpa','heptopascal','heptopascals')
      iunittype = 8; iunitid = 2      
    case('kpa','kilopascal','kilopascals')
      iunittype = 8; iunitid = 3
    case('kg/cm^2')
      iunittype = 8; iunitid = 4
    case('kg/m^2')
      iunittype = 8; iunitid = 5
    case('mb','milibar')
      iunittype = 8; iunitid = 6
    case('b','bar')
      iunittype = 8; iunitid = 7
    case('mmhg')
      iunittype = 8; iunitid = 8
    case('inhg')
      iunittype = 8; iunitid = 9
    case('psf','lb/ft^2','lb/ft2','lb/foot^2','lb/foot2')
      iunittype = 8; iunitid = 10
    case('psi','lb/in^2','lb/inch^2''lb/in2','lb/inch2')
      iunittype = 8; iunitid = 11
      
    !--- Salinity ---------------------------------
    case('pph','%') !parts per hundred
      iunittype = 9; iunitid = 1
    case('ppt','g/kg','psu') !parts per thousand or practical salinity units
      iunittype = 9; iunitid = 2
    case('ppm','mg/kg','mgr/kg') !parts per million
      iunittype = 9; iunitid = 3
    case('kg/kg','g/g') !mass fraction 
      iunittype = 9; iunitid = 4
      
    !--- Volume Flux [L^3/T] ---------------------------------------------------------------
    case('cms','m^3/s','cu m/s','m^3/sec','meter^3/second','meter^3/seconds')
      iunittype = 10; iunitid = 1
    case('cmh','m^3/hr','cu m/hr','m^3/hour','meter^3/hour','meter^3/hours')
      iunittype = 10; iunitid = 2
    case('cmd','m^3/d','cu m/d','m^3/day','meter^3/day','meter^3/days')
      iunittype = 10; iunitid = 3
    case('cfs','ft^3/s','cu ft/s','ft^3/sec','feet^3/second','feet^3/seconds')
      iunittype = 10; iunitid = 4
    case('cfh','ft^3/hr','cu ft/hr','ft^3/hour','feet^3/hour','feet^3/hours')
      iunittype = 10; iunitid = 5
    case('cfd','ft^3/d','cu ft/d','ft^3/day','feet^3/day','feet^3/days')
      iunittype = 10; iunitid = 6
    case('ls','l/s','l/sec','l/second','l/seconds')
      iunittype = 10; iunitid = 7
    case('lh','l/h','l/hr','l/hour','l/hours')
      iunittype = 10; iunitid = 8
    case('ld','l/d','l/day','l/days')
      iunittype = 10; iunitid = 9
    case('cys','yd^3/s','cu yd/s','yd^3/sec','yard^3/second','yard^3/seconds')
      iunittype = 10; iunitid = 10  
    case('cyh','yd^3/hr','cu yd/hr','yard^3/hour','yard^3/hours')
      iunittype = 10; iunitid = 11
    case('cyd','yd^3/d','cu yd/d','yd^3/day','yard^3/day','yard^3/days')
      iunittype = 10; iunitid = 12
      
    !--- Volume Units [L^3] ------------------------------------------------------
    case('mm^3','millimeter^3','cu millimeter','cu millimeters')
      iunittype = 11; iunitid = 1
    case('cm^3','centimeter^3','cu centimeter','cu centimeters')
      iunittype = 11; iunitid = 2
    case('m^3','cu m','cu meter','cu meters')
      iunittype = 11; iunitid = 3
    case('km^3','cu kilometer','cu kilometers')
      iunittype = 11; iunitid = 4
    case('in^3','inch^3','cu inch','cu inches')
      iunittype = 11; iunitid = 5
    case('ft^3','feet^3','foot^3')
      iunittype = 11; iunitid = 6
    case('yd^3','yds^3','yard^3','cu yards')
      iunittype = 11; iunitid = 7
    case('sm^3','smi^3','mile^3','cu mile','mi^3','cu mi','statutemile^3','statmi^3')
      iunittype = 11; iunitid = 8
    case('nm^3','nmi^3','nauticalmile^3','cu nauticalmile','nautmi^3')
      iunittype = 11; iunitid = 9
    case('gal','gallon','gallons')
      iunittype = 11; iunitid = 10  
    case('ml','mililiter','mililitre')
      iunittype = 11; iunitid = 11  
    case('l','liter','litre')
      iunittype = 11; iunitid = 12
      
    !--- Acceleration [L/T^2] ----------------------------------------------------------
    case('m/s^2','meter/second^2','m/sec^2')
      iunittype = 12; iunitid = 1
    case('ft/s^2','feet/second^2','ft/sec^2')
      iunittype = 12; iunitid = 2
      
    !--- Diffusivity [L^2/T] ----------------------------------------------------------
    case('m^2/s','meter^2/second','m^2/sec')
      iunittype = 13; iunitid = 1
    case('ft^2/s','feet^2/second','ft^2/sec')
      iunittype = 13; iunitid = 2
      
    !--- Specific Weight [M/L^2/T^2] -----------------------------------------------------  
    case('kg/m^2/s^2','kg/m^2/sec^2','kg/(m^2/s)','kg/(m^2/sec)')
      iunittype = 14; iunitid = 1
    case('lb/ft^2/s^2','lb/ft^2/sec^2','lb/(ft^2/s)','lb/(ft^2/sec)')
      iunittype = 14; iunitid = 2
      
    !--- Erosion Rate [M/L^2/T] --------------------------------------------------------    
    
    !--- Units not found --------------------
    case default
      iunittype = 0; iunitid = 0; ierr = -1
      
    end select 
    
    return
    end subroutine unitcodes    
    end subroutine unitconv_var

end module unitconv_lib    