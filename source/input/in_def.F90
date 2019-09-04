!=======================================================
module in_def    
! Variable definitions for input files
!written by Alex Sanchez, USACE-CHL
!=======================================================    
    use prec_def
    implicit none
    
!--- SMS ASCII --------------------------------------------------    
    !Dataset Type
    type dat_type  
      character(len=50) :: name   
      integer :: id
      integer :: nd   !# of data values
      integer :: nc   !# of cells or elements
      integer :: nt   !# of time steps
      logical :: usereftime
      real(ikind) :: reftime  !Reference time as a Julian number
      integer :: istatus
      real(ikind), allocatable :: time(:)
      character(len=10) :: time_units
      integer :: itype   !0-applied to nodes, 1-applied to elements/cells    
      integer :: ndim    !1-Scalar, 2-2D Vector, 3-2D Vector
      integer,     allocatable :: stat(:,:)  !(cell,time)
      real(ikind), allocatable :: val(:,:,:) !(cell,comp,time)
      !real(ikind), allocatable :: scal(:,:)  !(cell,time)
      !real(ikind), allocatable :: vec(:,:,:) !(cell,comp,time)
      character(len=10) :: units
    endtype dat_type  

    !Scalar Type
    type scaldattype 
      character(len=50) :: name            
      integer :: id
      integer :: nd   !# of data values
      integer :: nc   !# of cells or elements
      integer :: nt   !# of time stepss
      logical :: usereftime
      real(ikind) :: reftime  !Reference time as a Julian number
      integer :: istatus
      real(ikind), allocatable :: time(:)
      character(len=10) :: time_units
      integer,    allocatable :: stat(:,:) !(cell,time)
      real(ikind), allocatable :: val(:,:)  !(cell,time)
      character(len=10) :: units
    endtype scaldattype    
    
    !Vector Type
    type vecdattype  
      character(len=50) :: name   
      integer :: id
      integer :: nd   !# of data values
      integer :: nc   !# of cells or elements
      integer :: nt   !# of time steps
      logical :: usereftime
      real(ikind) :: reftime  !Reference time as a Julian number
      integer :: istatus
      real(ikind), allocatable :: time(:)
      character(len=10) :: time_units
      integer :: vectype  !0-applied to nodes, 1-applied to elements/cells      
      integer,     allocatable :: stat(:,:)  !(cell,time)
      real(ikind), allocatable :: val(:,:,:) !(cell,comp,time)
      character(len=10) :: units
    endtype vecdattype    

endmodule in_def