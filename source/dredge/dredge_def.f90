!===================================================================
module dredge_def
! dredging simulation variables
! written by Chris Reed
!===================================================================
    use prec_def, only: ikind
    
    implicit none
    save
    
    logical :: dredging, write_dredge_diag
    integer :: ndredge_operations,num_place_areas  
    integer, allocatable  :: DredgeUnit(:)
    integer :: num_trigger_TS
    real(ikind) :: dredge_interval,dredge_time_lapse
    real(ikind), allocatable :: DredgeTS_vars(:,:)
    real(ikind), allocatable :: trigger_start(:),trigger_finish(:)
    real(ikind), allocatable :: bed_change(:)
    real(ikind), allocatable :: dredge_mat_gradation(:)
    character(len=32)  :: METHOD
    character(len=100) :: dredge_diag_file = 'dredge_module_diagnostics.txt'
    
    type dredge_operations_type
      logical :: active
      character(len=200) :: name
      character(len=200) :: DredgeSourceAreaFile,DredgeSourceAreaPath
      character(len=32) :: DredgeSourceArea
      character(len=200), allocatable :: DredgePlaceAreaFile(:),DredgePlaceAreaPath(:)   
      integer :: ndredge_placement_areas
      integer :: NumDredgeAreaCells    
      integer, allocatable :: DredgeAreaCells(:)
      real(ikind), allocatable:: dredge_depth(:)

      !used for trigger methods
      integer :: Trigger_Approach,num_trigger_intervals   
      real(ikind) :: Trigger_Depth,Trigger_vol,Trigger_Percentage
      real(ikind), allocatable :: trigger_start(:),trigger_finish(:)
      
      !used for controling dredge operations (all approaches)
      integer :: Dredge_Approach
      real(ikind) :: Rate
      
      !used for dredge approach number 2
      integer :: Dredge_start_cell, cell_inc
      integer, allocatable :: cellorder(:)
      
      !used for placment of dredged material
      character*32, allocatable :: DredgePlacementArea(:)    
      integer :: NumPlacementAreas,PAmethod
      integer :: num_dredge_intervals
      integer, allocatable :: Placement_Approach(:),placement_cell_inc(:),Placement_Start_Cell(:)
      integer, allocatable :: PlacementAreaCells(:,:)
      integer, allocatable :: NumPlacementAreaCells(:)
      real(ikind) :: total_dredged_vol,total_placed_vol
      real(ikind) :: Placement_vol,Placement_excess
      real(ikind), allocatable :: Placement_Limit(:,:),Placement_area(:),Placement_thickness(:),Placement_Depth(:) 
      real(ikind), allocatable :: DredgePlacementAllocation(:) 
      real(ikind), allocatable :: total_placed_by_area(:),placed_vol_by_area(:)
      real(4), allocatable :: start(:),finish(:)

    end type dredge_operations_type
    type(dredge_operations_type), allocatable :: dredge_operations(:)      
    
    end module dredge_def
