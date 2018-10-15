!===================================================================
module sal_def
! Salinity transport variables
! written by Alex Sanchez, USACE-CHL
!===================================================================
    use prec_def
    implicit none
    save
    
    logical :: saltrans, salinitconst
    integer :: nsalobs,numsalbnd,itermaxsal,icsal
    real(ikind) :: salic,tolsal,schmidtsal
    !real(ikind), allocatable :: salbnd(:)  !Salinity for each cellstring, constant over each cellstring
    real(ikind), allocatable :: sal(:),sal1(:),sal2(:) 
    real(ikind), allocatable :: dsalx(:),dsaly(:) !Gradients
    real(ikind), allocatable :: susal0(:)
    real(ikind), allocatable :: rssal(:)
    real(ikind), allocatable :: salobspts(:,:) !Contains x,y,z,id    
    character(len=200) :: salobsfile
    character(len=200), allocatable :: saldrvpath(:)    
    character(len=200) :: salfile,salpath

    integer :: nsalstr
    type sal_driver
      integer              :: ibndstr !Boundary cell string
      integer              :: ncells
      integer, allocatable :: cells(:)
      integer, allocatable :: faces(:) 
      character(len=200)        :: bidfile,bidpath
      integer :: ntimes
      integer :: inc
      real(ikind) :: salbnd              !Salinity for each cellstring. Note: Constant over each cellsstring    
      real(4), allocatable :: timesal(:) !Time [hours]
      real(4), allocatable :: val(:)     !Salinity(time), Time-series of salinity values
      real(ikind), allocatable :: salbnd0(:) !Salinity(cell), initial boundary salinity values
      character(len=200) :: salfile,salpath
    endtype sal_driver
    type(sal_driver), allocatable :: sal_str(:)
    
endmodule sal_def