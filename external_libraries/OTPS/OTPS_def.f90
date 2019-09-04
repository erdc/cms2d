!================================================================
module otps_var
!================================================================    
    implicit none
    character*1 :: zuv
    character*80 :: modname !Model name
    character*80 :: lltname !Latitude and longitude
    character*80 :: outname !Output name
    logical :: APRI !true = output amplitude and phase (GMT), false = real/imaginary parts.
    logical :: Geo  !extract ocean (false) or geocentric (true) HC for elevations only
    logical :: interp_micon !true = correct for minor constituents    
    
endmodule otps_var

!================================================================
module otps_constit
!   This file contains the standard parameters which define the
!   amplitudes, frequencies, etc. for the primary tidal constituents
!   Currently knows about 20 tidal constituents:
!      NOW KNOWS about 21: M4 also
!================================================================
      integer, parameter :: ncmx = 29
      character*4 constid(ncmx)
      data constid/'m2  ','s2  ','k1  ','o1  ',&
                   'n2  ','p1  ','k2  ','q1  ',&
                   '2n2 ','mu2 ','nu2 ','l2  ',&
                   't2  ','j1  ','m1  ','oo1 ',&
                   'rho1','mf  ','mm  ','ssa ',&
                   'm4  ','ms4 ','mn4 ','m6  ',&
                   'm8  ','mk3 ','s6  ','2sm2',&
                   '2mk3'/

!    FOR EACH POSSIBLE CONSTIUENT, these parameters are given:
!    alpha = correction factor for first order load tides
!    amp = amplitude of equilibrium tide in m
!    ph = Currently set to zero ...   phases for
!             each constituent are referred to the time
!             when the phase of the forcing for that
!             constituent is zero on the Greenich meridian.)

!    omega = angular frequency of constituent, in radians
      real :: alpha_d(ncmx),ph_d(ncmx),amp_d(ncmx),omega_d(ncmx),&
           phase_mkB(ncmx),beta_SE(ncmx)
      integer :: ispec_d(ncmx)

!     Tidal parameters taken from Rodney's constituent.h, 2/23/96:
!     (except for ispec).
      data ispec_d/ &
          2,2,1,1,&
          2,1,2,1,&
          2,2,2,2,&
          2,1,1,1,&
          1,0,0,0,&
          0,0,0,0,&
          0,0,0,0,&
          0/
!     note: for now I am just leaving ispec for M4 set to 0 (ispec
!     is only used to define forcing in atgf, and this is always  0
!     for M4)

      data alpha_d/ &
          0.693,0.693,0.736,0.695,&
          0.693,0.706,0.693,0.695,&
          0.693,0.693,0.693,0.693,&
          0.693,0.695,0.695,0.695,&
          0.695,0.693,0.693,0.693,&
          0.693,0.693,0.693,0.693,&
          0.693,0.693,0.693,0.693,&
          0.693/

      data omega_d/ &
          1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,&
          1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,&
          1.352405e-04,1.355937e-04,1.382329e-04,1.431581e-04,&
          1.452450e-04,7.556036e-05,7.028195e-05,7.824458e-05,&
          6.531174e-05,0.053234e-04,0.026392e-04,0.003982e-04,&
          2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04,&
          5.620755e-04,2.134402e-04,4.363323e-04,1.503693e-04,&
          2.081166e-04/

      data ph_d/29*0.0/

      data amp_d/ &
         0.242334,0.112743,0.141565,0.100661,&
         0.046397,0.046848,0.030684,0.019273,&
         0.006141,0.007408,0.008811,0.006931,&
         0.006608,0.007915,0.007915,0.004338,&
         0.003661,0.042041,0.022191,0.019567,&
!       amplitude for M4 etc. is zero
         0.,0.,0.,0.,&
         0.,0.,0.,0.,&
         0./
 
! Astronomical arguments, obtained with Richard Ray's
! "arguments" and "astrol", for Jan 1, 1992, 00:00 Greenwich time
! Corrected July 12, 2000  
       data phase_mkB/ &
         1.731557546,0.000000000,0.173003674,1.558553872,&
         6.050721243,6.110181633,3.487600001,5.877717569,&
         4.086699633,3.463115091,5.427136701,0.553986502,&
         0.052841931,2.137025284,2.436575100,1.929046130,&
         5.254133027,1.756042456,1.964021610,3.487600001,&
         3.463115091,1.731557546,1.499093481,5.194672637,&
         6.926230184,1.904561220,0.000000000,4.551627762,&
         3.809122439/
! I am putting 0 for ms2,mn4 etc. for now: correct later
! Now this correction is done using the SAL file (h_TPXO3_90-90.load)
! I replace beta_SE with units for now (on case we decide to switch back
! to old version) and comment the old numbers - this way I do NOT change
! anything in subroutines  
! This was in weights.h before - placed here not to mix with w!
! to remove solid Earth tide multily by beta:
       data beta_SE/ &
         0.9540,0.9540,0.9400,0.9400,&
         0.9540,0.9400,0.9540,0.9400,&
         0.9540,0.9540,0.9540,0.9540,&
         0.9540,0.9400,0.9400,0.9400,&
         0.9400,0.9400,0.9400,0.9400,&
!      for M4 just using value for semi-diurnals (no good reason!)
         0.9540,0.9540,0.9540,0.954,&
         0.9540,0.9540,0.9540,0.954,&
         0.9540/
!       data beta_SE/29*1./


endmodule otps_constit

!================================================================
module otps_weights
! Gary's approach to minor constituents interpolation
! DISABLED in this package
! This is comments from old weights.h
!/* 
! * weights.h:: original file from Gary for NLP=4 
! *
! * $Id: weights.h,v 1.4 1996/02/13 22:54:26 rodney Exp $
! */
! Redone by Lana 1998.02.23 for a new FORTRAN version
! This is good only for 8 constituents in order:
! 'm2  ','s2  ','k1  ','o1  ','n2  ','p1  ','k2  ','q1  '
! 
! The same order is supported in constit.h - so
! do not need to care about the corresponding indices
!================================================================
    real w(17,8),beta(17)
    data w(1,:)/1.0,  .00,  .00,  .00,  .00,  .00,  .00,  .00/
    data w(2,:)/0.0, 1.00,  .00,  .00,  .00,  .00,  .00,  .00/
    data w(3,:)/0.0,  .00,  1.0,  .00,  .00,  .00,  .00,  .00/
    data w(4,:)/0.0,  .00,  .00, 1.00,  .00,  .00,  .00,  .00/
    data w(5,:)/0.0,  .00,  .00,  .00, 1.00,  .00,  .00,  .00/
    data w(6,:)/0.0,  .00,  .00,  .00,  .00, 1.00,  .00,  .00/
    data w(7,:)/0.0,  .00,  .00,  .00,  .00,  .00, 1.00,  .00/
    data w(8,:)/0.0,  .00,  .00,  .00,  .00,  .00,  .00, 1.00/
    data w(9,:)/-0.0379, .0,.00,  .00,  .30859 ,0.0, .03289,.0/
    data w(10,:)/-0.03961,.0,.00,  .00,  .34380, 0.0, .03436,.0/
    data w(11,:)/.00696,  .0,.00,  .00,  .15719, 0.0, -.00547,.0/
    data w(12,:)/.02884,  .0,.00,  .00, -.05036, 0.0,  .07424,.0/
    data w(13,:)/.00854,  .0,.00,  .00, -.01913, 0.0,  .17685,.0/
    data w(14,:)/.0,  .0, -.00571, .11234, .0, .05285, .0, -.26257/
    data w(15,:)/.0,  .0,  .00749, .07474, .0, .03904, .0, -.12959/
    data w(16,:)/.0,  .0, -.03748, .12419, .0, .05843, .0, -.29027/
    data w(17,:)/.0,  .0,  .00842, .01002, .0,-.03064, .0,  .15028/

    
endmodule otps_weights
