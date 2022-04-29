!======================================================================
    module watch_def
! CMS Watch Library
! 
! Description:
!   Used to time sections of the code and profile the code efficiency
!
! Authors:
!   Alex Sanchez, UASCE-CHL
!======================================================================    
    implicit none   
    
    integer :: nwatches      !# of active watches
    integer, parameter :: nmaxwatches = 50 !Maximum # of watches    
    integer, parameter :: watchcharlen = 30 !character length for watches

    type chronotype
      real(8) :: lap
      real(8) :: split
    endtype chronotype
    
    type watchtype
      logical :: active
      logical :: running   
      character(len=watchcharlen) :: name
      type(chronotype) :: cpu
      type(chronotype) :: clock
    endtype watchtype
    
    type(watchtype), allocatable, dimension(:) :: watch
    type(watchtype) :: watchdef
   
end module watch_def