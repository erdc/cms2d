!*************************************************
module CSHORE_def
! CSHORE Variable Definitions
!
! Author: Brad Johnson, USACE-CHL
! Modifications:
!  Chris Reed, URS
!  Alex Sanchez, USACE-CHL
!*************************************************
#include "CMS_cpp.h"
#ifdef CSHORE
    use prec_def,ONLY: ikind
    implicit none
    save
    
    !Note: Brad, please add comments below to describe each variable including the units
    logical CSHORE_ON
    real(ikind), allocatable :: QCSx(:),QCSy(:)  !Sediment transport rates [kg/m/s]
    real(ikind), allocatable :: Calpha(:),Salpha(:),DTLoc(:)
    real(ikind), allocatable :: BSLOPE(:),GSLOPERAW(:),GSLOPE(:),ASLOPERAW(:)
    real(ikind), allocatable :: ASLOPE(:),ASLOPEX(:),ASLOPEY(:)  !Smoothed bedslope
    real(ikind), allocatable :: CP(:),Re(:),RBETA(:),DBSTA(:),RQ(:),RQnew(:),RQf(:,:)
    real(ikind):: CTHETA, STHETA, SQR2, RBZERO, SQRG1,SQRG2, EFFF, EFFB
    real(ikind):: sgm1,sporo1,gsgm1,blp,shield, bld
    real(ikind):: GSLMAX,TANPHI,SLP,GAMMA_CS
    real(ikind):: SG, SPORO, SQR8, FB2
    
#endif
endmodule CSHORE_def