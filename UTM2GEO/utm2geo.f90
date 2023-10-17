      MODULE UTM2GEO

! Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39180
! June 7, 2011
!
! Added a module called UTM2GEO that contains a set of 
! subroutines for converting to and from UTM to Geographic
! horizontal coordinate systems.  The original routines were
! converted to fixed-format, fixed bugs where mixed mode
! arthimetic was being used, changed all real values to .D0
! type instead of default, and commented out the call to update
! the grid_zone in the subroutine ll2utm.
!
      CONTAINS

!*************************************************************************
!/*
! * Peter Daly
! * MIT Ocean Acoustics
! * pmd@mit.edu
! * 25-MAY-1998
! * 
! Revisions:
!   Jan. 25, 1999 DHG  Port to Fortran 90
!   Mar. 23, 1999 DHG  To add Lewis Dozier's fix to "rr1" calculation 
! * 
! Description:
! * 
! * These routines convert UTM to Lat/Longitude and vice-versa,
! * using the WGS-84 (GPS standard) or Clarke 1866 Datums.
! * 
! * The formulae for these routines were originally taken from
! * Chapter 10 of "GPS: Theory and Practice," by B. Hofmann-Wellenhof,
! * H. Lictenegger, and J. Collins. (3rd ed) ISBN: 3-211-82591-6,
! * however, several errors were present in the text which
! * made their formulae incorrect.
! *
! * Instead, the formulae for these routines was taken from
! * "Map Projections: A Working Manual," by John P. Snyder
! * (US Geological Survey Professional Paper 1395)
! *
! * Copyright (C) 1998 Massachusetts Institute of Technology
! *               All Rights Reserved
! *
! * RCS ID: $Id: convert_datum.c,v 1.2 1998/06/04 20:50:47 pmd Exp pmd $
! */
!
!*************************************************************************
!
      subroutine get_grid_zone (longitude, latitude, grid_zone, lambda0)

      IMPLICIT NONE

      real (kind=8) longitude, latitude
      integer       grid_zone(2)
      real (kind=8) lambda0

      integer  zone_long, zone_lat

      real (kind=8) M_PI
!!!   parameter (M_PI = 3.141592654)

!-------------------------------------------------------------------------

      m_pi = ACOS (-1.0d0)

!  /* Solve for the grid zone, returns the central meridian */

      zone_long = INT ((longitude + 180.0d0) / 6.0d0) + 1
      zone_lat = NINT ((latitude + 80.0d0) / 8.0d0)
      grid_zone(1) = zone_long
      grid_zone(2) = zone_lat

!  /* First, let's take care of the polar regions */

      if ((latitude < -80.0d0) .OR. (latitude > 84.0d0)) then
         lambda0 = 0.0d0 * M_PI / 180.0d0
         return
      endif

!  /* Now the special "X" grid */

      if (latitude .GT. 72.0d0 .AND. longitude .GT. 0.0d0 .AND. longitude .LT. 42.0d0) then
         if (longitude .LT. 9.0d0) then
            lambda0 = 4.5d0 * M_PI / 180.0d0
         elseif (longitude .LT. 21.0d0) then
            lambda0 = 15.0d0 * M_PI / 180.0d0
         elseif (longitude .LT. 33.0d0) then
            lambda0 = 27.0d0 * M_PI / 180.0d0
         elseif (longitude .LT. 42.0d0) then
            lambda0 = 37.5d0 * M_PI / 180.0d0
         endif
         return
      endif

!  /* Handle the special "V" grid */

      if (latitude .GT. 56.0d0 .AND. latitude .LT. 64.0d0 .AND. longitude .GT. 0.0d0 .AND. longitude .LT. 12.0d0) then
         if (longitude .LT. 3.0d0) then
            lambda0 = 1.5d0 * M_PI / 180.0d0
         elseif (longitude .LT. 12.0d0) then
            lambda0 = 7.5d0 * M_PI / 180.0d0
         endif
         return
      endif

!  /* The remainder of the grids follow the standard rule */

      lambda0 = (FLOAT (zone_long - 1) * 6.0d0 + (-180.0d0) + 3.0d0) * M_PI / 180.0d0
  
      return
      end subroutine get_grid_zone

!*************************************************************************

      subroutine get_lambda0 (grid_zone, lambda0, ierr)

      IMPLICIT NONE

      integer       grid_zone(2)
      real (kind=8) lambda0
      integer       ierr

      integer zone_long
      integer zone_lat
      real (kind=8) latitude, longitude

      real (kind=8) M_PI
!!!   parameter (M_PI = 3.141592654)

!---------------------------------------------------------------------------


      m_pi = ACOS (-1.0d0)

  !/* Given the grid zone, then set the central meridian, lambda0 */

  !/* Check the grid zone format */

      zone_long = grid_zone(1)
      zone_lat = grid_zone(2)
      if ((zone_long .LT. 1) .OR. (zone_long .GT. 61)) then
         write (*,*) 'Invalid grid zone format: ', zone_long, zone_lat
         ierr = -1
         return 
      endif

      longitude = (FLOAT (zone_long - 1) * 6.0d0) - 180.0d0
      latitude = (FLOAT (zone_lat) * 8.0d0) - 80.0d0

  !/* Take care of special cases */

      if ((latitude .LT. -80.0d0) .OR. (latitude .GT. 84.0d0)) then
         lambda0 = 0.0d0
         ierr = 0
         return 
      endif

      if (latitude .GT. 56.0d0 .AND. latitude .LT. 64.0d0 .AND. longitude .GT. 0.0d0 .AND. longitude .LT. 12.0d0) then
         if (longitude .LT. 3.0d0) then
            lambda0 = 1.5 * M_PI / 180.0d0
         elseif (longitude .LT. 12.0d0) then
            lambda0 = 7.5d0 * M_PI / 180.0d0
         endif
         ierr = 0
         return
      endif
  
      if (latitude .GT. 72.0d0 .AND. longitude .GT. 0.0d0 .AND. longitude < 42.0d0) then
         if (longitude .LT. 9.0d0) then
            lambda0 = 4.5d0 * M_PI / 180.0d0
         elseif (longitude .LT. 21.0d0) then
            lambda0 = 15.0d0 * M_PI / 180.0d0
         elseif (longitude .LT. 33.0d0) then
            lambda0 = 27.0d0 * M_PI / 180.0d0
         elseif (longitude .LT. 42.0d0) then
            lambda0 = 37.5d0 * M_PI / 180.0d0
         endif
         ierr = 0
         return
      endif

  !/* Now handle standard cases */

      lambda0 = (FLOAT (zone_long - 1) * 6.0d0 + (-180.0d0) + 3.0d0) * M_PI / 180.0d0

  !/* All done */

      ierr = 0
      return
      end subroutine get_lambda0

!*************************************************************************
      subroutine ll2utm (longitude, latitude, utm_x, utm_y, grid_zone, lambda0, datum)

! tcm 20110608 -- commented out the call to get_grid_zone and added lambda0 to input data

      IMPLICIT NONE

      real (kind=8) latitude, longitude
      real (kind=8) utm_x, utm_y
      integer       grid_zone(2)
      integer       datum

      real (kind=8)  a, b, f, e, e2, e4, e6
      real (kind=8)  phi, lambda, lambda0, phi0, k0
      real (kind=8)  t, rho, x, y, mm, mm0                !m, k   Never used, so commented out from declaration  MEB  01/26/2022
      real (kind=8)  aa, aa2, aa3, aa4, aa5, aa6
      real (kind=8)  ep2, nn, tt, cc

      real (kind=8) M_PI
!!!   parameter (M_PI = 3.141592654)

      integer CLARKE_1866_DATUM
      parameter (CLARKE_1866_DATUM = 1)
      integer GRS_80_DATUM
      parameter (GRS_80_DATUM = 2)
      integer WGS_84_DATUM
      parameter (WGS_84_DATUM = 3)

!---------------------------------------------------------------------------


      m_pi = ACOS (-1.0d0)

  !/* Converts lat/long to UTM, using the specified datum */

      if (datum == CLARKE_1866_DATUM) then      ! CLARKE_1866_DATUM:
        a = 6378206.4d0
        b = 6356583.8d0
      elseif (datum == GRS_80_DATUM) then      ! GRS_80_DATUM:
        a = 6378137.0d0
        b = 6356752.3d0
      elseif (datum == WGS_84_DATUM) then      ! WGS_84_DATUM:
        a = 6378137.0d0           !/* semimajor axis of ellipsoid (meters) */
        b = 6356752.31425d0       !/* semiminor axis of ellipsoid (meters) */
      else
        write (*,*) 'Unknown datum: ', datum
        return
      endif

  !/* Calculate flatness and eccentricity */

      f = 1.d0 - (b / a)
      e2 = 2.d0 * f - f * f
      e = sqrt (e2)
      e4 = e2 * e2
      e6 = e4 * e2

  !/* Convert latitude/longitude to radians */
  
      phi = latitude * M_PI / 180.0d0
      lambda = longitude * M_PI / 180.0d0

  !/* Figure out the UTM zone, as well as lambda0 */

! tcm 20110608 -- comment out the call to get_grid_zone
!      call get_grid_zone (longitude, latitude, grid_zone, lambda0)

      phi0 = 0.0d0

  !/* See if this will use UTM or UPS */

      if (latitude .GT. 84.0d0) then

    !/* use Universal Polar Stereographic Projection (north polar aspect) */

         k0 = 0.994d0

         t = sqrt ( ((1.0d0 - sin (phi)) / (1.0d0 + sin (phi))) * (((1.0d0 + e * sin (phi)) / (1. - e * sin (phi))) ** e) )
         rho = 2.0 * a * k0 * t / sqrt ( ((1.0d0 + e) ** (1.0d0 + e)) * ((1.0d0 - e) ** (1.0d0 - e)) )
    !!! Not needed (dhg) m = cos (phi) / sqrt (1.0 - e2 * sin (phi) * sin (phi))

         x = rho * sin (lambda - lambda0)
         y = -rho * cos (lambda - lambda0)
    !!! Not needed (dhg) k = rho * a * m

    !/* Apply false easting/northing */

         x = x + 2000000.0d0
         y = y + 2000000.0d0

      elseif (latitude .LT. -80.0d0) then

    !/* use Universal Polar Stereographic Projection (south polar aspect) */

         phi = -phi
         lambda = -lambda
         lambda0 = -lambda0

         k0 = 0.994d0

         t = sqrt (((1.0d0 - sin (phi)) / (1.0d0 + sin (phi))) * ( ( (1.0d0 + e * sin (phi)) / (1.0d0 - e * sin (phi)) ** e) ) )
         rho = 2.0d0 * a * k0 * t / sqrt ( ((1.0d0+e) ** (1.0d0+e)) * ((1.0d0-e) ** (1.0d0-e)) )
    !!! Not needed (dhg) m = cos (phi) / sqrt (1.0 - e2 * sin (phi) * sin (phi))

         x = rho * sin (lambda - lambda0)
         y = -rho * cos (lambda - lambda0)
    !!! Not needed (dhg) k = rho * a * m

         x = -x
         y = -y

    !/* Apply false easting/northing */

         x = x + 2000000.0d0
         y = y + 2000000.0d0

      else

    !/* Use UTM */

    !/* set scale on central median (0.9996 for UTM) */
    
         k0 = 0.9996d0

         mm = a * ((1.0d0-e2/4.0d0 - 3.0d0*e4/64.0d0 - 5.0d0*e6/256.0d0) * phi -              &
              (3.0d0*e2/8.0d0 + 3.0d0*e4/32.0d0 + 45.0d0*e6/1024.0d0) * sin (2.0d0*phi) +     &
              (15.0d0*e4/256.0d0 + 45.0d0*e6/1024.0d0) * sin (4.0d0*phi) -                    &
              (35.0d0*e6/3072.0d0) * sin (6.0d0*phi))

         mm0 = a * ((1.0d0-e2/4.0d0 - 3.0d0*e4/64.0d0 - 5.0d0*e6/256.0d0) * phi0 -            &
               (3.0d0*e2/8.0d0 + 3.0d0*e4/32.0d0 + 45.0d0*e6/1024.0d0) * sin (2.0d0*phi0) +   &
               (15.0d0*e4/256.0d0 + 45.0d0*e6/1024.0d0) * sin (4.0d0*phi0) -                  &
               (35.0d0*e6/3072.0d0) * sin (6.0d0*phi0))

         aa = (lambda - lambda0) * cos(phi)
         aa2 = aa * aa
         aa3 = aa2 * aa
         aa4 = aa2 * aa2
         aa5 = aa4 * aa
         aa6 = aa3 * aa3

         ep2 = e2 / (1.0d0 - e2)
         nn = a / sqrt (1.0 - e2 * sin (phi) * sin (phi))
         tt = tan (phi) * tan (phi)
         cc = ep2 * cos (phi) * cos (phi)

    !!! Not needed (dhg) k = k0 * (1 + (1+cc)*aa2/2 + (5-4*tt+42*cc+13*cc*cc-28*ep2) * aa4 / 24.0 + &
    !!! Not needed (dhg)          (61-148*tt+16*tt*tt) * aa6 / 720.0)
         x = k0 * nn * (aa + (1.0d0-tt+cc) * aa3 / 6.0d0 + (5.0d0-18.0d0*tt+tt*tt+72.0d0*cc-58.0d0*ep2) * aa5 / 120.0d0)
         y = k0 * (mm - mm0 + nn * tan (phi) *                                               &
                  (aa2 / 2.0d0 + (5.0d0-tt+9.0d0*cc+4.0d0*cc*cc) * aa4 / 24.0d0 +            &
                  (61.0d0 - 58.0d0*tt + tt*tt + 600.0d0*cc - 330.0d0*ep2) * aa6 / 720.0d0))

    !/* Apply false easting and northing */

         x = x + 500000.0d0
         if (y .LT. 0.0d0) then
            y = y + 10000000.0d0
         endif
      endif

  !/* Set entries in UTM structure */

      utm_x = x
      utm_y = y

  !/* done */

      return
      end subroutine ll2utm

!*************************************************************************
      subroutine utm2ll (utm_x, utm_y, longitude, latitude, grid_zone, datum)

      IMPLICIT NONE
 
      real (kind=8) utm_x, utm_y
      real (kind=8) latitude, longitude
      integer       grid_zone(2)
      integer       datum

      integer ierr
      real (kind=8)  a, b, f, e, e2, e4, e6, e8
      real (kind=8)  lambda0, x, y, k0, rho, t, chi, phi, phi1, phit
      real (kind=8)  lambda, phi0, e1, e12, e13, e14
      real (kind=8)  mm, mm0, mu, ep2, cc1, tt1, nn1, rr1
      real (kind=8)  dd, dd2, dd3, dd4, dd5, dd6

      real (kind=8) M_PI
!!!   parameter (M_PI = 3.141592654)
      real (kind=8) LOWER_EPS_LIMIT
      parameter (LOWER_EPS_LIMIT = 1.0d-14)
      real (kind=8) M_PI_2

      integer CLARKE_1866_DATUM
      parameter (CLARKE_1866_DATUM = 1)
      integer GRS_80_DATUM
      parameter (GRS_80_DATUM = 2)
      integer WGS_84_DATUM
      parameter (WGS_84_DATUM = 3)

!---------------------------------------------------------------------------


      m_pi = ACOS (-1.0d0)

      M_PI_2 = M_PI * 2.0d0

  !/* Converts UTM to lat/long, using the specified datum */

      if (datum == CLARKE_1866_DATUM) then      ! CLARKE_1866_DATUM:
         a = 6378206.4d0
         b = 6356583.8d0
      elseif (datum == GRS_80_DATUM) then       ! GRS_80_DATUM:
         a = 6378137.0d0
         b = 6356752.3d0
      elseif (datum == WGS_84_DATUM) then       ! WGS_84_DATUM:
         a = 6378137.0d0             !/* semimajor axis of ellipsoid (meters) */
         b = 6356752.31425d0         !/* semiminor axis of ellipsoid (meters) */
      else
         write (*,*) 'Unknown datum: ', datum
         return
      endif

  !/* Calculate flatness and eccentricity */

      f = 1.0d0 - (b / a)
      e2 = (2.0d0 * f) - (f * f)
      e = sqrt (e2)
      e4 = e2 * e2
      e6 = e4 * e2
      e8 = e4 * e4

  !/* Given the UTM grid zone, generate a baseline lambda0 */

      call get_lambda0 (grid_zone, lambda0, ierr)
      if (ierr .NE. 0) then
         write (*,*) 'Unable to translate UTM to LL'
         return
      endif

      latitude = ( FLOAT(grid_zone(2)) * 8.0d0 ) - 80.0d0

  !/* Take care of the polar regions first. */

      if (latitude .GT. 84.0d0) then !/* north polar aspect */

    !/* Subtract the false easting/northing */

         x = utm_x - 2000000.0d0
         y = utm_y - 2000000.0d0

    !/* Solve for inverse equations */

         k0 = 0.994d0
         rho = sqrt (x*x + y*y)
         t = rho * sqrt ( ((1.0d0+e) ** (1.0d0+e)) * ((1.0d0-e) ** (1.0d0-e)) ) / (2.0d0*a*k0)

    !/* Solve for latitude and longitude */

         chi = M_PI_2 - 2.0d0 * atan (t)
         phit = chi + (e2/2.0d0 + 5.0d0*e4/24.0d0 + e6/12.0d0 + 13.0d0*e8/360.0d0) * sin(2.0d0*chi) + &
                  (7.0d0*e4/48.0d0 + 29.0d0*e6/240.0d0 + 811.0d0*e8/11520.0d0) * sin(4.0d0*chi) +     &
                  (7.0d0*e6/120 + 81.0d0*e8/1120.0d0) * sin(6.0d0*chi) +                              &
                  (4279.0d0*e8/161280.0d0) * sin(8.0d0*chi)

         do while (ABS (phi-phit) .GT. LOWER_EPS_LIMIT)
            phi = phit
            phit = M_PI_2 - 2.0d0 * atan ( t * (((1.0d0 - e * sin (phi)) / (1.0d0 + e * sin (phi))) ** (e / 2.0d0)) )
         enddo

         lambda = lambda0 + atan2 (x, -y)

      elseif( latitude .LT. -80.0d0 ) then !/* south polar aspect */

    !/* Subtract the false easting/northing */

         x = -(utm_x - 2000000.0d0)
         y = -(utm_y - 2000000.0d0)

    !/* Solve for inverse equations */

         k0 = 0.994d0
         rho = sqrt (x*x + y*y)
         t = rho * sqrt ( ((1.0d0+e) ** (1.0d0+e)) * ((1.0d0-e) ** (1.0d0-e)) ) / (2.0d0*a*k0)

    !/* Solve for latitude and longitude */

         chi = M_PI_2 - 2.0d0 * atan (t)
         phit = chi + (e2/2.0d0 + 5.0d0*e4/24.0d0 + e6/12.0d0 + 13.0d0*e8/360.0d0) * sin (2.0d0*chi) + &
            (7.0d0*e4/48.0d0 + 29.0d0*e6/240.0d0 + 811.0d0*e8/11520.0d0) * sin (4.0d0*chi) +           &
            (7.0d0*e6/120.0d0 + 81.0d0*e8/1120.0d0) * sin (6.0d0*chi) +                                &
            (4279.0d0*e8/161280.0d0) * sin (8.0d0*chi)

         do while (ABS (phi-phit) .GT. LOWER_EPS_LIMIT)
            phi = phit;
            phit = M_PI_2 - 2.0d0 * atan (t * ( ((1.0d0-e*sin(phi)) / (1.0d0+e*sin(phi)) ) ** (e/2.0d0)))
         enddo

         phi = -phi
         lambda = -(-lambda0 + atan2 (x , -y))

      else

    !/* Now take care of the UTM locations */

         k0 = 0.9996d0

    !/* Remove false eastings/northings */

         x = utm_x - 500000.0d0
         y = utm_y

         if (latitude .LT. 0.0d0) then  !/* southern hemisphere */
            y = y - 10000000.0d0
         endif

    !/* Calculate the footpoint latitude */

         phi0 = 0.0d0
         e1 = (1.0d0 - sqrt (1.0d0-e2)) / (1.0d0 + sqrt (1.0d0-e2))
         e12 = e1 * e1
         e13 = e1 * e12
         e14 = e12 * e12

         mm0 = a * ((1.0d0-e2/4.0d0 - 3.0d0*e4/64.0d0 - 5.0d0*e6/256.0d0) * phi0 -         &
             (3.0d0*e2/8.0d0 + 3.0d0*e4/32.0d0 + 45.0d0*e6/1024.0d0) * sin (2.0d0*phi0) +  &
             (15.0d0*e4/256.0d0 + 45.0d0*e6/1024.0d0) * sin (4.0d0*phi0) -                 &
             (35.0d0*e6/3072.0d0) * sin (6.0d0*phi0))
         mm = mm0 + y/k0;
         mu = mm / (a * (1.0d0-e2/4.0d0-3.0d0*e4/64.0d0-5.0d0*e6/256.0d0))

         phi1 = mu + (3.0d0*e1/2.0d0 - 27.0d0*e13/32.0d0) * sin (2.0d0*mu) +  &
             (21.0d0*e12/16.0d0 - 55.0d0*e14/32.0d0) * sin (4.0d0*mu) +       &
             (151.0d0*e13/96.0d0) * sin (6.0d0*mu) +                          &
             (1097.0d0*e14/512.0d0) * sin (8.0d0*mu)

    !/* Now calculate lambda and phi */

         ep2 = e2 / (1.0d0 - e2)
         cc1 = ep2 * cos (phi1) * cos (phi1)
         tt1 = tan (phi1) * tan (phi1)
         nn1 = a / sqrt (1.0d0 - e2 * sin (phi1) * sin (phi1))
    !!!DHG Old Code rr1 = a * (1.0 - e2) / ((1.0 - e2 * sin (phi) * sin (phi)) ** 1.5)
    !!!DHG L.Dozier's fix is next
         rr1 = a * (1.0d0 - e2) / ((1.0d0 - e2 * sin (phi1) * sin (phi1)) ** 1.5d0)
         dd = x / (nn1 * k0)

         dd2 = dd * dd
         dd3 = dd * dd2
         dd4 = dd2 * dd2
         dd5 = dd3 * dd2
         dd6 = dd4 * dd2

         phi = phi1 - (nn1 * tan (phi1) / rr1) *                                                             &
           (dd2/2.0d0 - (5.0d0+3.0d0*tt1+10.0d0*cc1-4.0d0*cc1*cc1-9.0d0*ep2) * dd4 / 24.0d0 +                &
           (61.0d0+90.0d0*tt1+298.0d0*cc1+45.0d0*tt1*tt1-252.0d0*ep2-3.0d0*cc1*cc1) * dd6 / 720.0d0)
         lambda = lambda0 + (dd - (1.0d0+2.0d0*tt1+cc1) * dd3 / 6.0d0 +                                      &
           (5.0d0-2.0d0*cc1+28.0d0*tt1-3.0d0*cc1*cc1+8.0d0*ep2+24.0d0*tt1*tt1) * dd5 / 120.0d0) / cos (phi1)
      endif

  !/* Convert phi/lambda to degrees */
  
      latitude = phi * 180.0d0 / M_PI
      longitude = lambda * 180.0d0 / M_PI
  
  !/* All done */

      return
      end subroutine utm2ll


      END MODULE UTM2GEO