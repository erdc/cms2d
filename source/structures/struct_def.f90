!===================================================================
module struct_def
! Structures Module
! written by Weiming Wu, NCCHE
! modified by Alex Sanchez
!===================================================================
    use prec_def
    save
        
    !Tidal gate
    integer numtidegate
    integer,allocatable :: ntidegate(:),idtidegate(:),opentidegate(:),ntidegateopt(:),mtidegateopt(:), &
                           orienttgbay(:),orienttgsea(:),orienttidegate(:),methtidegate(:)
    real(ikind),allocatable :: coeftidegate(:,:),coeftglateral(:),elevtidegate(:),openhgttidegate(:),   &
                               qtidegate(:),dqtgdzdown(:),dqtgdzup(:),tidegateopt(:),Qtottidegate(:) 

    !Weir
    integer :: numweir
    integer,allocatable :: nweir(:),idweir(:),iweirtype(:),methweir(:),  &
                           orientweirbay(:),orientweirsea(:),orientweir(:) 
    real(ikind),allocatable :: coefweir(:,:),coefweirlateral(:),elevweir(:),qweir(:),   &
                               dqweirdzdown(:),dqweirdzup(:),Qtotweir(:)  
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
!-------------

    !Culvert
    integer :: numculvert
    integer,allocatable :: nculvert(:),idculvert(:,:),iculverttype(:),iculvertflap(:) 
    real(ikind),allocatable :: culvertrad(:),culvertwidth(:),culvertheight(:),culvertelevbay(:),&
                               culvertelevsea(:),culvertlength(:),cvheadlossbayentr(:),  &
                               cvheadlossbayexit(:),cvheadlossseaentr(:),cvheadlossseaexit(:), &
                               culvertfrict(:),culvertmann(:),qculvert(:),dqcvdzdown(:),dqcvdzup(:), &
                               uvculvert(:),angleculvertbay(:),angleculvertsea(:)

    !Rubble mound
    integer :: numrubmound
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
      integer              :: rubmoundmeth
    endtype RM_type
    type(RM_type), allocatable :: RM_struct(:)
    logical :: rm_dset
    !real :: rockdia_const,structporo_const,structbaseD_const,rubmoundmeth  !hli(12/11/12)
 !-------
    
    !
    integer :: rmblock  !hli(12/11/12)
    character(len=200) :: arubmoundfile,arubmoundpath
    character(len=200) :: arockdiamfile,arockdiampath
    character(len=200) :: astructporofile,astructporopath
    character(len=200) :: astructbaseDfile,astructbaseDpath
    character(len=200) :: astructmethfile,astructmethpath
    character(len=200) :: astructIDfile, astructIDpath          !meb incomplete - finish or remove  09/25/23
    character(len=34)  :: methrm           !hli(12/11/12)

end module struct_def
