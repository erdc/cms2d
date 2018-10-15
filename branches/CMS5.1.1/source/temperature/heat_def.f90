!===================================================================
module heat_def
! Heat transfer variables
!  USACE-CHL
!===================================================================
    use prec_def
    implicit none
    save
    
    logical :: heattrans, heatinitconst
    integer :: nheatobs,numheatbnd,itermaxheat,icheat
    real(ikind) :: heatic,tolheat,schmidtheat
    !real(ikind), allocatable :: heatbnd(:)  !Temperature for each cellstring, constant over each cellstring
    real(ikind), allocatable :: heat(:),heat1(:),heat2(:) 
    real(ikind), allocatable :: dheatx(:),dheaty(:) !Gradients
    real(ikind), allocatable :: suheat0(:)
    real(ikind), allocatable :: rsheat(:)
    real(ikind), allocatable :: heatobspts(:,:) !Contains x,y,z,id    
    character(len=200) :: heatobsfile
    character(len=200), allocatable :: heatdrvpath(:)    
    character(len=200) :: heatfile,heatpath

    integer :: nheatstr
    type heat_driver
      integer              :: ibndstr !Boundary cell string
      integer              :: ncells
      integer, allocatable :: cells(:)
      integer, allocatable :: faces(:) 
      character(len=200)        :: bidfile,bidpath
      integer :: ntimes
      integer :: inc
      real(ikind) :: heatbnd              !Temperature for each cellstring. Note: Constant over each cellsstring    
      real(4), allocatable :: timeheat(:) !Time [hours]
      real(4), allocatable :: val(:)     !Temperature(time), Time-series of temperature values
      real(ikind), allocatable :: heatbnd0(:) !Temperature(cell), initial boundary temperature values
      character(len=200) :: heatfile,heatpath
    endtype heat_driver
    type(heat_driver), allocatable :: heat_str(:)

    !Time Series of Input Parameters
    integer numhtflux
    real cloud, dewpt
    real airtmp, velwind
    real solarms
    real,allocatable :: clouddat(:),dewptdat(:)
    real,allocatable :: airtmpdat(:),dathtflux(:)
    real,allocatable :: solarmsdat(:)
    
 endmodule heat_def
