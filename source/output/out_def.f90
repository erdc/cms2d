!===================================================================    
module out_def
!Output Module
!===================================================================
    use prec_def
    implicit none
    save
    
!--- Simulation label -------------
   character(len=100) :: simlabel  !Used for XMDF output path
    
!--- File Names -------------------------------------------------
    character(len=200) :: hisfile    !Tecplot history file
    character(len=200) :: datfile    !Tecplot snapshot file
    character(len=200) :: goutfile   !Global XMDF Solution File _sol.h5
    character(len=200) :: dgoutfile  !Diagnostic global output file _diag.h5 (for model crashes)
    character(len=200) :: outpath    !Path directory for output files
    character(len=200) :: outprefix  !Prefix of output files including path
    
!--- Output Times List ---------------------------------   
   type output_time_list
     character(len=11) :: name
     logical :: write_dat              !Write output list group
     logical :: use_default=.true.     !To use the default name or not
     integer :: ilist                  !List id number
     integer :: inc                    !Counter for current output time
     integer :: ntimes                 !Number of output times
     real(8), allocatable :: times(:)  !Output times
     character(len=200) :: afile       !Output file name and path
   endtype output_time_list
   integer, parameter :: noutlist = 14 !One for each group   
   type(output_time_list) :: outlist(noutlist)
   logical :: use_common = .false.     !Use common solution file named '*_sol.h5'
   
!--- Output Time Series -------------------------
   type input_time_series
     integer :: ntimes
     real(8), allocatable :: times(:)
   endtype input_time_series       
   integer, parameter :: nmaxinputlist = 4    
   type(input_time_series) :: outseries(nmaxinputlist)         
   
!--- Toggles -------------------------------------------------  
    integer, parameter :: incXMDF_IO = 1 !Toggle for using compression in XMDF output files 
    logical :: write_xmdf         !Turns on or off XMDF output files
    integer :: ixmdfcomp          !Toggle for using compression in XMDF output files 
    logical :: write_sup          !Turns on or off super ascii output files
    logical :: write_tecplot      !Turns on or off the tecplot output files
    
    !Geometry
    logical :: write_areap        !Turns on or off the current velocity magnitude = uv
    logical :: write_cell2cell    !Turns on or off the forward cell to cell connectivity
    logical :: write_llec2llec    !Turns on or off the backwards cell to cell connectivity
    logical :: write_cell2node    !Turns on or off the forward cell to node connectivity    
    logical :: write_node2cell    !Turns on or off the forward node to cell connectivity    
    logical :: write_geocells     !Turns on or off the geometry file for all cells
    
    !Hydro
    logical :: write_velmag       !Turns on or off the current velocity magnitude = uv
    logical :: write_veltotal     !Turns on or off the total flux velocities = u,v
    logical :: write_veltotalmag  !Turns on or off the total flux velocity magnitude
    logical :: write_pres         !Turns on or off the wse pressure = p
    logical :: write_totdep       !Turns on or off the total depth = h
    logical :: write_hpred        !Turns on or off the total depth = hpred
    logical :: write_wsepred      !Turns on or off the water level = hpred+zb
    logical :: write_velpred      !Turns on or off the current velocities = upred and vpred
    logical :: write_velcurv      !Turns on or off the velocity curvature = (d2u/dx2,d2u/dy2,d2v/dx2,d2v/dy2)
    logical :: write_wsegrad      !Turns on or off the wse gradient = (deta/dx,deta/dy)
    logical :: write_velgrad      !Turns on or off the wse gradient = (du/dx,du/dy,dv/dx,dv/dy)
    logical :: write_prescor      !Turns on or off the pressure correction = pp
    logical :: write_prescorgrad  !Turns on or off the pressure correction gradients = (dppx,dppy)
    logical :: write_presgrad     !Turns on or off the pressure gradients = (dpx,dpy)
    logical :: write_presres      !Turns on or off the normalized residuals for p
    logical :: write_velres       !Turns on or off the normalized residuals for u,v
    logical :: write_streamcurv   !Turns on or off the streamwise curvature
    logical :: write_volres       !Turns on or off the volume residual
    
    !Sediment transport
    logical :: write_conc         !Turns on or off the total-load sediment concentration output
    logical :: write_capac        !Turns on or off the total-load sediment capacity output
    logical :: write_concsusp     !Turns on or off the suspended-load concentration output
    !logical :: write_concbed      !Turns on or off the suspended-load concentration output
    logical :: write_fracsusp     !Turns on or off the fraction of suspended transport output
    logical :: write_concfrac     !Turns on or off the fracational sediment concenctrations
    logical :: write_bedload      !Turns on or off the bed-load transport rate
    logical :: write_suspload     !Turns on or off the suspended-load transport rate
    logical :: write_sizefrac     !Turns on or off the fractional bed composition
    logical :: write_thickness    !Turns on or off the bed layer thickness
    logical :: write_morph        !Turns on or off the morphology change output
    logical :: write_alphatot     !Turns on or off the adaptation coefficient for total-load
    logical :: write_lengthtot    !Turns on or off the adaptation length for total-load
    logical :: write_lengthsusp   !Turns on or off the adaptation length for suspended-load
    logical :: write_lengthbed    !Turns on or off the adaptation length for bed-load
    logical :: write_betatot      !Turns on or off the correction factor for suspended-load
    logical :: write_betasusp     !Turns on or off the correction factor for bed-load
    logical :: write_bedvel       !Turns on or off the bed-load velocity
    
    !Waves
    logical :: write_wave_detail  !Turns on or off the wave detailed variables    
    logical :: write_wavbrkind    !Turns on or off the wave breaking index
    logical :: write_wavbrkdiss   !Turns on or off the wave dissipation output
    logical :: write_wavstress    !Turns on or off the wave radiation stress gradient output
    logical :: write_fracbreak    !Turns on or off the fraction of breaking waves
    logical :: write_wavdisscoef  !Turns on or off the fraction of breaking waves
    logical :: write_wavorbvel    !Turns on or off the representative wave bottom orbital velocity
    logical :: write_wavlen       !Turns on or off the wave length
    logical :: write_wavvel       !Turns on or off the wave mass transport velocity vector
    logical :: write_wavvelmag    !Turns on or off the wave mass transport velocity magnitude
    
    !Roller
    logical :: write_rollenergy   !Turns on or off the roller energy
    logical :: write_rollstress   !Turns on or off the roller stress
    logical :: write_rolldiss     !Turns on or off the roller dissipation
    
    !Wind and Atmospheric Pressure
    logical :: write_wndvel       !Turns on or off the wind velocity magnitude output
    logical :: write_wndmag       !Turns on or off the wind velocity magnitude output
    logical :: write_wndstress    !Turns on or off the wind stress vector
    logical :: write_wndstressmag !Turns on or off the wind stress magnitude
    logical :: write_atmpres      !Turns on or off the atmospheric pressure
    logical :: write_atmpresgrad  !Turns on or off the atmospheric pressure gradient
    
    !Bed roughness
    logical :: write_meanbedshear !Turns on or off the mean bed shear stress
    logical :: write_normrough    !Turns on or off the normalized bed roughness length
    logical :: write_normapprough !Turns on or off the normalized apparent roughness length
    logical :: write_rough        !Turns on or off the total bed roughness
    logical :: write_roughrip     !Turns on or off the roughness due to rippples
    logical :: write_roughmegarip !Turns on or off the roughness due to megarippples
    logical :: write_roughdune    !Turns on or off the roughness due to dunes
    logical :: write_roughgrain   !Turns on or off the roughness due to grains
    logical :: write_roughtrans   !Turns on or off the roughness due to sediment transport
    logical :: write_mannings     !Turns on or off the Manning's roughness coefficient
    
    !ASCII/XMDF Controls
    logical :: write_ascii_input  !Turns on or off outputting ASCII input files
    logical :: write_xmdf_output  !Turns on or off outputting XMDF input files
    
!--- Observation cells (save points) --------------------------------
    logical :: obs_cell = .false.   !Any observation cell
    type obs_cell_driver     
      character(len=20) :: group
      logical :: active     
      real(ikind) :: time_inc             !Output time increment for each group [sec]
      integer :: nvar                     !Number of output variables in each group
      integer,allocatable :: units(:)     !Output file unit for each variable     
      character(len=10),allocatable:: names(:) !Name of each variable in group
      character(len=200),allocatable :: files(:)     
      integer :: ncells !Number of cells
      integer,allocatable :: cells(:)
      character(len=37),allocatable:: identifiers(:) !Identifier for each cell     
    endtype obs_cell_driver   
    type(obs_cell_driver) :: obs(4) !Groups
    
!--- New save point vars ------------   5/7/12 - meb
    logical :: save_point = .false.
    character(len=100) :: splabel = ' ' !Used for Save Points output path
    
    type save_pt_group
      character(len=10) :: group
      logical :: active
      real(ikind) :: time_inc             !Output time increment for each group
      integer :: ncells                        !Number of cells
      integer :: nvar   !Number of output variables in each group
      integer,           allocatable :: funits(:)   !Output file unit for each variable     
      character(len=200),allocatable :: files(:)    !File names for each variable 
      character(len=10), allocatable :: ounit(:)    !Output units for each variable, ie. "m/sec"
      integer,           allocatable :: vals(:)     !Scalar (1) or Vector (2)
      character(len=10), allocatable :: names(:)    !Name of each variable in group
      integer,           allocatable :: cell(:)     !Cell number for each cell. 
      character(len=37), allocatable :: id(:)       !Identifier for each cell.  Doesn't have to be the Cell #     
      real(ikind),       allocatable :: x(:), y(:)  !Store X and Y for each cell 
    endtype save_pt_group
    integer, parameter :: ngroups = 5
    type(save_pt_group) :: savept(ngroups) !Groups  1=Hydro, 2=Sediment, 3=Salinity, 4=Waves, 5=Morphology

endmodule out_def
