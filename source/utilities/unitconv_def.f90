!=====================================================================
module unitconv_def
! UNIT CONVersion variable DEFinitions
!
! Author: Alex Sanchez, USACE-CHL
!=====================================================================
    use prec_def
    implicit none

    !--- Distances ----
    !per meter
    real(ikind),parameter :: ftpermm = 0.00328084_ikind   !foot per milimeter
    real(ikind),parameter :: ftperm = 3.28084_ikind       !foot per meter
    real(ikind),parameter :: inperm = 39.3701_ikind       !Inch per meter
    real(ikind),parameter :: ydperm = 1.09361_ikind       !Yard per meter
    real(ikind),parameter :: miperm = 0.000621371_ikind   !Mile per meter
    real(ikind),parameter :: nmiperm = 0.000539957_ikind  !Nautical mile per meter
    !per foot
    real(ikind),parameter :: inperft = 12.0_ikind         !Inch per foot
    real(ikind),parameter :: ydperft = 0.333333_ikind     !Yard per foot
    real(ikind),parameter :: mperft = 0.3048_ikind        !Meter per foot
    real(ikind),parameter :: mmperft = 304.8_ikind        !Millimeter per foot
    real(ikind),parameter :: miperft = 0.000189394_ikind  !Mile per foot
    real(ikind),parameter :: nmiperft = 0.000164579_ikind !Nautical mile per foot
    
    !--- Time ---
    real(ikind),parameter :: secpermin = 60.0_ikind       !Seconds per minute
    real(ikind),parameter :: secperhr = 3600.0_ikind      !Seconds per hour
    real(ikind),parameter :: secperday = 86400.0_ikind    !Seconds per day
    real(ikind),parameter :: secperyr = 31536000.0_ikind  !Seconds per year
    real(ikind),parameter :: minperhr = 60.0_ikind        !Minutes per hour
    real(ikind),parameter :: minperday = 1440.0_ikind     !Minutes per day
    
    !!--- Angle ---
    !real(ikind),parameter :: pikind = acos(-1.0_ikind)
    !real(ikind),parameter :: deg2rad = pikind/180.0_ikind
    !real(ikind),parameter :: rad2deg = 180.0_ikind/pikind
    
    !--- Mass ---
    real(ikind),parameter :: lbperkg = 2.20462262_ikind      !lb per kg
    real(ikind),parameter :: kgperlb = 0.45359237_ikind      !kg per lb
    real(ikind),parameter :: kgpertonUK = 1016.046909_ikind  !kg per tonUK
    real(ikind),parameter :: kgpertonLong = kgpertonUK       !kg per tonlong
    real(ikind),parameter :: lbpertonUK = 2240.0_ikind       !lb per tonUK
    real(ikind),parameter :: lbpertonLong = lbpertonUK       !lb per tonLong
    real(ikind),parameter :: kgpertonUS = 907.18474_ikind    !kg per tonUS
    real(ikind),parameter :: kgpertonShort = kgpertonUS      !kg per tonShort
    real(ikind),parameter :: lbpertonUS = 2000.0_ikind       !lb per tonUS
    real(ikind),parameter :: lbpertonShort = lbpertonUS      !lb per tonShort
    real(ikind),parameter :: kgpertonne = 1000.0_ikind       !kg per tonne
    real(ikind),parameter :: lbpertonne = 2204.622622_ikind  !lb per tonne
    
    !--- Area ---
    real(ikind),parameter :: ft2perm2 = 10.76391042_ikind    !ft^2 per m^2
    real(ikind),parameter :: m2perft2 = 0.09290304_ikind     !m^2 per ft^2
    
    !--- Volume ---
    real(ikind),parameter :: ft3perm3 = 35.31472_ikind       !ft^3 per m^3
    real(ikind),parameter :: m3perft3 = 0.02831680_ikind     !m^3 per ft^3
    
    !--- Density ----
    real(ikind),parameter :: kgft3perlbm3 = 16.02_ikind      !kg/m^3 to lb/ft^3
    real(ikind),parameter :: lbm3perkgft3 = 0.06242796_ikind !lb/ft^3 to kg/m^3
    
end module unitconv_def