!===================================================================
module struct_def
! Structures Module
! written by Weiming Wu, NCCHE
! modified by Alex Sanchez
! updated by Mitchell Brown, 10/27/2023
!===================================================================
    use prec_def
    save
        
!==================================================================    
!Tidal gates
    integer numtidegate
    integer,allocatable :: ntidegate(:),idtidegate(:),opentidegate(:),ntidegateopt(:),mtidegateopt(:), &
                           orienttgbay(:),orienttgsea(:),orienttidegate(:),methtidegate(:)
    real(ikind),allocatable :: coeftidegate(:,:),coeftglateral(:),elevtidegate(:),openhgttidegate(:),   &
                               qtidegate(:),dqtgdzdown(:),dqtgdzup(:),tidegateopt(:),Qtottidegate(:) 
    !Added meb 10/23/2024
    character(len=12), allocatable :: aschedtype(:)
    character(len=5), allocatable  :: aorienttgsea(:)
    integer :: itidegate = 0
    type TG_type
    !This type facilitates reading in of the information from multiple blocks.
    !Information from this type is assigned other defined arrays.
      integer :: ncells
      integer, allocatable :: cells(:)
      real    :: coefTGlateral 
      integer :: orientTGsea       !1= North, 2= East, 3= South, 4= West
      real    :: coefTG_b2s        !Bay to Sea
      real    :: coefTG_s2b        !Sea to Bay
      real    :: open_height       !Opening height
      real    :: bottom_elev       !Bottom elevation - Positive upward
      integer :: method            !1= Approach 1, 2= Approach 2      
      integer :: mtidegateopt      !1= REG, 2= DES, 3= EBB, 4= UCG
      integer :: reg_start_time    !Start time for REGular gates
      integer :: reg_open_freq     !Open frequency for REGular gates
      integer :: reg_open_dur      !Open duration for REGular gates
      integer :: ntidegateopt      !Number of times for all gates
      integer, allocatable :: des_start_time(:) !Start time for DESignated gates
      integer, allocatable :: des_open_dur(:)   !Open duration for DESignated gates
    endtype TG_type
    type(TG_type),allocatable :: TG_struct(:)

!==================================================================    
!Weirs
    integer :: numweir
    integer,allocatable :: nweir(:),idweir(:),iweirtype(:),methweir(:),  &
                           orientweirbay(:),orientweirsea(:),orientweir(:) 
    real(ikind),allocatable :: coefweir(:,:),coefweirlateral(:),elevweir(:),qweir(:),   &
                               dqweirdzdown(:),dqweirdzup(:),Qtotweir(:)  
    character(len=5),allocatable :: aweirtype(:),aorientweirsea(:)
!Added meb 01/28/2019
    integer :: iweir = 0
    type WR_type
    !This type facilitates reading in of the information from multiple blocks.
    !Information from this type is assigned other defined arrays.
      integer :: ncells
      integer, allocatable :: cells(:)
      real    :: coefWeirLateral
      integer :: orientWeir       !1= North, 2= East, 3= South, 4= West
      integer :: weirType         !1= Sharp-crested, 2=Broad-crested
      real    :: coefWeir_b2s     !Bay to Sea
      real    :: coefWeir_s2b     !Sea to Bay
      real    :: elevWeir         !Positive is upward
      integer :: methWeir         !1= Approach 1, 2= Approach 2
    endtype WR_type
    type(WR_type),allocatable :: WR_struct(:)

!==================================================================    
!Culverts
    integer :: numculvert
    integer,allocatable :: nculvert(:),idculvert(:,:),iculverttype(:),iculvertflap(:) 
    real(ikind),allocatable :: culvertrad(:),culvertwidth(:),culvertheight(:),culvertelevbay(:),&
                               culvertelevsea(:),culvertlength(:),cvheadlossbayentr(:),  &
                               cvheadlossbayexit(:),cvheadlossseaentr(:),cvheadlossseaexit(:), &
                               culvertfrict(:),culvertmann(:),qculvert(:),dqcvdzdown(:),dqcvdzup(:), &
                               uvculvert(:),angleculvertbay(:),angleculvertsea(:)
    character(len=3), allocatable :: aculverttype(:)
    
    !Added meb 10/22/2024
    integer :: iculvert = 0
    type CV_type
    !This type facilitates reading in of the information from multiple blocks.
    !Information from this type is assigned other defined arrays.
      integer :: cell_bay          !Bayside cell id
      integer :: cell_sea          !Seaside cell id
      integer :: culvertType       !1= Circle, 2=Box
      logical :: culvertFlap       !True if ON, False if OFF
      real    :: radius            !Radius if Type is Circular
      real    :: height            !Height if Type is Box
      real    :: width             !Width if Type is Box
      real    :: length            !Length of culvert
      real    :: darcyCoef         !
      real    :: manningsCoef
      real    :: invertElevBay
      real    :: invertElevSea
      real    :: entryCoefWeir_b2s !Bay to Sea
      real    :: entryCoefWeir_s2b !Sea to Bay
      real    :: exitCoefWeir_b2s  !Bay to Sea
      real    :: exitCoefWeir_s2b  !Sea to Bay
      real    :: outflowAngle_bay
      real    :: outflowAngle_sea
    endtype CV_type
    type(CV_type),allocatable :: CV_struct(:)

!==================================================================    
!Rubble mounds
    integer :: nrubmoundcells      !refactored 'numrubmound' to 'nrubmoundcells' for clarity.  All of these look so similar.  10/27/2023
    integer,allocatable :: nrubmound(:),idrubmound(:)
    real(ikind),allocatable :: rubmounddia(:),rubmoundporo(:),rubmounda(:),rubmoundb(:),rubmoundbotm(:)
    real(ikind),allocatable :: permeability(:),rockdiam(:),structporo(:),structbaseD(:)  !hli(12/11/12)
    integer,allocatable     :: methrubmoundab(:)
    real(ikind),allocatable :: structmeth(:)  !hli(12/11/12)
    character(len=200),allocatable :: rubmoundname(:)

 !Added meb 1/18/2019
    !This type facilitates reading in of the information from multiple blocks.
    !Information from this type is assigned other defined arrays.
    integer :: irubmound = 0  !present index of rubble mound inputs from block
    type RM_type
      integer              :: ncells
      character (len=200)  :: name
      integer, allocatable :: cells(:)
      real                 :: rockdia_const
      real                 :: structporo_const
      real                 :: structbaseD_const
      character(len=200)   :: rockdia_dset(2)      !Added meb 10/27/2023
      character(len=200)   :: structporo_dset(2)   !Added meb 10/27/2023
      character(len=200)   :: structbaseD_dset(2)  !Added meb 10/27/2023 
      character(len=200)   :: structmeth_dset(2)   !Added meb 10/27/2023
      integer              :: rubmoundmeth
      integer              :: inc                  !Added meb 10/27/2023
    endtype RM_type
    type(RM_type), allocatable :: RM_struct(:)
    logical :: rm_dset
    
    integer :: rmblock     !hli(12/11/12)
    character(len=200) :: arubmoundfile,    arubmoundpath='Empty'   !Points to file that contains non-zero values in cells where there is permeability.  Each number denotes a different structure. 10/27/2023
    character(len=200) :: arockdiamfile,    arockdiampath='Empty'           
    character(len=200) :: astructporofile,  astructporopath='Empty'
    character(len=200) :: astructbaseDfile, astructbaseDpath='Empty'
    character(len=200) :: astructmethfile,  astructmethpath='Empty'
    character(len=34)  :: methrm           !hli(12/11/12)

end module struct_def
