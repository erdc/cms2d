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
    real(ikind),allocatable :: methrubmoundab(:),structmeth(:)  !hli(12/11/12)
    real :: rockdia_const,structporo_const,structbaseD_const,rubmoundmeth  !hli(12/11/12)

    integer :: rmblock  !hli(12/11/12)
    character(len=200) :: arubmoundfile,arubmoundpath
    character(len=200) :: arockdiamfile,arockdiampath
    character(len=200) :: astructporofile,astructporopath
    character(len=200) :: astructbaseDfile,astructbaseDpath
    character(len=200) :: astructmethfile,astructmethpath
    character(len=34)  :: methrm           !hli(12/11/12)

endmodule struct_def
