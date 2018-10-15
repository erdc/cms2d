cmake_minimum_required(VERSION 3.2)

# CMake control file for building the CMS model

# Gather library requirements-- these are a Fortran compiler, a recent
# version of the HDF5 library and the XMDF library.
ENABLE_LANGUAGE( Fortran )

# Source
SET( CMS_SOURCE_FILES
  ./main/CMS_cpp.h
  ./main/definitions.F90
  ./boundaries/bnd_def.f90
  ./boundaries/bnd_block.f90
  ./boundaries/bnd_lib.f90
  ./boundaries/bnd_vel.f90
  ./boundaries/boundaries.F90
  ./boundaries/nesting.F90
  ./boundaries/nest_lib.f90
  ./boundaries/nest_par.f90
  ./boundaries/tide_lib.F90
  ./boundaries/tide_fac.f
  ./derivatives/der_def.f90
  ./derivatives/der_lib.f90
  ./derivatives/derivatives.f90
  ./dredge/dredge_def.f90
  ./dredge/dredge_op.F90
  ./dredge/dredge_setup.F90
  ./explicit/initialize_flow_arrays.f90
  ./explicit/initialize_flow_arrays_Tel_v02.f90
  ./explicit/ZCMS_FLOW_EXP_C.f90
  ./explicit/ZCMS_FLOW_EXP_Main.f90
  ./explicit/ZCMS_FLOW_EXP_T.f90
  ./explicit/ZderivativesEXP.f90
  ./explicit/ZdiswallEXP.f90
  ./explicit/Zimplicit_deallocate.f90
  ./explicit/ZInit_Activity.f90
  ./explicit/ZInit_Activity_tel.f90
  ./explicit/Zinitialize_extrapolations.f90
  ./explicit/Zinitialize_extrapolations_Tel.f90
  ./explicit/Zinitialize_structures.f90
  ./explicit/ZModules.F90
  ./explicit/Zprestart_EXP.f90
  ./explicit/Zread_card_file.f90
  ./explicit/Zsed_exp.F90
  ./explicit/Zsed_exp_tel.F90
  ./explicit/ZST_adeq.f90
  ./explicit/ZST_avalanche.f90
  ./explicit/ZST_bed_changes.f90
  ./explicit/ZST_cohesive.f90
  ./explicit/ZST_hardbottom.f90
  ./explicit/ZST_hardbottom_AD.f90
  ./explicit/ZST_hardbottom_COHES.f90
  ./explicit/ZST_lundcirp.f90
  ./explicit/ZST_morphology.f90
  ./explicit/ZST_watanabe.f90
  ./explicit/ZST_wet-dry_check.f90
  ./explicit/ZSub_update_rainf_evap.f90
  ./explicit/ZSubs-Init.f90
  ./explicit/ZSubs-Init_tel.f90
  ./explicit/ZSubs-Other.f90
  ./explicit/Zupdate_ADSS_conc.f90
  ./explicit/Zupdate_advection_Tel_Reg_v00.f90
  ./explicit/Zupdate_advection_Tel_v00.f90
  ./explicit/Zupdate_COHES_conc.f90
  ./explicit/Zupdate_cohes_flow_bc.f90
  ./explicit/Zupdate_culverts.f90
  ./explicit/Zupdate_momentum.f90
  ./explicit/Zupdate_momentum_Tel_reg_v00.f90
  ./explicit/Zupdate_momentum_Tel_v00.f90
  ./explicit/Zupdate_q_and_u_extrapolations.f90
  ./explicit/Zupdate_q_and_u_extrapolations_Tel.f90
  ./explicit/Zupdate_Q_bc.f90
  ./explicit/Zupdate_Q_bc_Tel.f90
  ./explicit/Zupdate_Q_cellcenter_forSedTrans.f90
  ./explicit/Zupdate_salinity.f90
  ./explicit/Zupdate_salinity_bc.f90
  ./explicit/Zupdate_salinity_bc_tel.f90
  ./explicit/Zupdate_salinity_tel.f90
  ./explicit/Zupdate_salt_conc.f90
  ./explicit/Zupdate_salt_conc_tel.f90
  ./explicit/Zupdate_sedtrans.f90
  ./explicit/Zupdate_uv_from_qxqy.f90
  ./explicit/Zupdate_uv_from_qxqy_Tel.f90
  ./explicit/Zupdate_variables.f90
  ./explicit/Zupdate_variables_tel.f90
  ./explicit/Zupdate_WABC.f90
  ./explicit/Zupdate_wetdry.f90
  ./explicit/Zupdate_wetdry_tel.f90
  ./explicit/Zupdate_wse.f90
  ./explicit/Zupdate_wse_bc.f90
  ./explicit/Zupdate_wse_Tel.f90
  ./flow/flow.F90
  ./flow/flow_def.f90
  ./flow/flow_eddyvis.f90
  ./flow/flow_imp.F90
  ./flow/flow_imp_flux.f90
  ./flow/flow_lib.F90
  ./friction/friction.F90
  ./friction/fric_def.f90
  ./friction/fric_lib.f90
  ./geospatial/geo_def.f90
  ./geospatial/geo_lib.f90
  ./geospatial/geospatial.F90
  ./geospatial/geo_test.F90
  ./hot_start/hot_def.f90
  ./hot_start/hot_start.F90
  ./input/in_def.F90
  ./input/in_lib.f90
  ./input/input.F90
  ./interpolation/interp_def.f90
  ./interpolation/interp_lib.f90
  ./interpolation/interpolation.F90
  ./interpolation/plagr_lib.f90
  ./main/CMS2D_main.f90
  ./main/CMS_Flow.f90
  ./main/CMS_test.F90
  ./main/prestart.f90
  ./meteorological/met_def.f90
  ./meteorological/met_lib.f90
  ./meteorological/meteorological.F90
  ./output/out_def.f90
  ./output/out_lib.f90
  ./output/output.F90
  ./quasi3d/q3d_def.F90
  ./quasi3d/q3d_lib.F90
  ./quasi3d/quasi3d.F90
  ./roller/rol_def.f90
  ./roller/roller.F90
  ./salinity/sal_bnd.f90
  ./salinity/sal_def.f90
  ./salinity/salinity.F90
  ./schemes/comp_lib.f90
  ./schemes/computational.F90
  ./sediment/sed_def.f90
  ./sediment/sed_lib.f90
  ./sediment/lund-cirp.F90
  ./sediment/sediment.F90
  ./sediment/sed_adapt.f90
  ./sediment/sed_aval.f90
  ./sediment/sed_bed.f90
  ./sediment/sed_beta.f90
  ./sediment/sed_beta_def.F90
  ./sediment/sed_beta_lib.F90
  ./sediment/sed_bnd.f90
  ./sediment/sed_capac.f90
  ./sediment/sed_hidexp.f90
  ./sediment/sed_imp.f90
  ./sediment/sed_slope.f90
  ./sediment/wuetal.F90
  ./sediment/sed_dry.f90
  ./solvers/solv_def.f90
  ./solvers/gauss-seidel.f90
  ./solvers/iccg.f90
  ./solvers/iccgstab.f90
  ./solvers/ilu0.f90
  ./solvers/ilutp.f90
  ./solvers/LU_lib.F90
  ./solvers/pbicgstab.f90
  ./solvers/pgmres.f90
  ./solvers/secant_lib.f90
  ./solvers/sip.f90
  ./solvers/solvers.f90
  ./solvers/sparse.f90
  ./stacking/solution_scheme_option.f90
  ./temperature/heat_def.f90
  ./temperature/heat_bnd.f90
  ./temperature/heat.f90
  ./temperature/heatflux.f90
  ./statistics/stat_def.f90
  ./statistics/statistics.F90
  ./steering/steer_def.f90
  ./steering/steering.f90
  ./structures/struct_def.f90
  ./structures/struct_lib.f90
  ./structures/structures.f90
  ./utilities/diag_def.f90
  ./utilities/diag_lib.f90
  ./utilities/diagnostics.f90
  ./utilities/math_lib.F90
  ./utilities/math_test.F90
  ./utilities/time_lib.F90
  ./utilities/time_test.F90
  ./utilities/tools.F90
  ./utilities/unitconv_lib.F90
  ./utilities/watch_def.F90
  ./utilities/watch_lib.F90
  ./utilities/watch_test.F90
  ./vegetation/veg_def.f90
  ./vegetation/vegetation.F90
  ./waves/wave_def.f90
  ./waves/wave_lib.F90
  ./waves/wave.F90
  ./waves/CMS-Wave_v3-2W_01Dec2014.f90
  ./waves/wave_rad_stress_fl.f90
  ../external_libraries/proj_conv/utm2geo.F
  ../external_libraries/proj_conv/spcs83_combined.F
)

IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")

  SET( CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-cpp -ffree-line-length-none -ffixed-line-length-none -O3" )

ELSEIF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")

  # Setting -static so that ifort runtime libraries are linked into the executable.
  # Setting -assume no underscore so C code compiled with GCC can be linked by ifort.
  SET( CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-fpp -extend_source 132 -static -assume nounderscore" )

ELSE()

  MESSAGE("Warning: Fortran compiler ${CMAKE_Fortran_COMPILER_ID} untested!")

ENDIF()

INCLUDE_DIRECTORIES( main )

# Does not work
# SET(Fortran_MODULE_DIRECTORY 
#	${PROJECT_BINARY_DIR}/build
# 	CACHE PATH
# 	"Single directory for all modules.")

ADD_EXECUTABLE( cms ${CMS_SOURCE_FILES} )

INSTALL( 
  TARGETS cms
  RUNTIME DESTINATION build
)
