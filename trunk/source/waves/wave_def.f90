!======================================================================
module global_inline
! CMS-Wave 
!======================================================================
    implicit none
    save
    
    real, parameter :: A=0.17  !A=0.15 suggested by Li's experiments
    integer, parameter :: NPF=50,MPD=50,IPMX=2500,JPMX=2500  !Alex: NPF should be maybe increased
    integer, parameter :: IGPX=IPMX-1,JGPX=JPMX-1,MPMX=MPD*JGPX
    integer, parameter :: IJPMX=IPMX*JPMX,IJGPX=IGPX*JGPX
    integer, parameter :: KOMX=4000,NOMX=80000
    
endmodule global_inline
    