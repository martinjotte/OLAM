#!/bin/bash
sed 's@Module consts_coms@Module ed_consts_coms@' <../ED2/memory/ed_consts_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_consts_coms.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_max_dims.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_max_dims.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/grid_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/grid_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/pft_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/pft_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/rk4_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/rk4_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/soil_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/soil_coms.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_misc_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_misc_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/allometry.f90 >dum.f90
mv dum.f90 ../ED2/utils/allometry.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/ed_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_therm_lib.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/dp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/dp_therm_lib.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/sp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/sp_therm_lib.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/c34constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/c34constants.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/canopy_air_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_air_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/canopy_layer_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_layer_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/canopy_radiation_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_radiation_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/decomp_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/decomp_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/disturb_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/disturb_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_mem_alloc.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_alloc.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_mem_grid_dim_defs.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_grid_dim_defs.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_state_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_state_vars.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_var_tables.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_var_tables.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_work_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_work_vars.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ename_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ename_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/fusion_fission_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/fusion_fission_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/ed_hdf5_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_hdf5_coms.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/hydrology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/hydrology_constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_constants.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/mem_polygons.f90 >dum.f90
mv dum.f90 ../ED2/memory/mem_polygons.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/met_driver_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/met_driver_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/optimiz_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/optimiz_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/phenology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/phenology_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/memory/physiology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/physiology_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/mpi/ed_mpass_init.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_mpass_init.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/mpi/ed_node_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_node_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/mpi/ed_para_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_coms.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/mpi/ed_para_init.F90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_init.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/driver/ed_1st.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_1st.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/driver/ed_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_driver.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/driver/ed_met_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_met_driver.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/driver/ed_model.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_model.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/budget_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/budget_utils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/charutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/charutils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/ed_dateutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_dateutils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/ed_filelist.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_filelist.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/ed_grid.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_grid.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/fatal_error.f90 >dum.f90
mv dum.f90 ../ED2/utils/fatal_error.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/fuse_fiss_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/fuse_fiss_utils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/great_circle.f90 >dum.f90
mv dum.f90 ../ED2/utils/great_circle.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/ed_hdf5_utils.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_hdf5_utils.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/invmondays.f90 >dum.f90
mv dum.f90 ../ED2/utils/invmondays.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/lapse.f90 >dum.f90
mv dum.f90 ../ED2/utils/lapse.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/libxml2f90.f90_pp.f90 >dum.f90
mv dum.f90 ../ED2/utils/libxml2f90.f90_pp.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/numutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/numutils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/radiate_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/radiate_utils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/rsys.F90 >dum.f90
mv dum.f90 ../ED2/utils/rsys.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/stable_cohorts.f90 >dum.f90
mv dum.f90 ../ED2/utils/stable_cohorts.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/update_derived_props.f90 >dum.f90
mv dum.f90 ../ED2/utils/update_derived_props.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/utils_c.c >dum.f90
mv dum.f90 ../ED2/utils/utils_c.c
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/utils/utils_f.f90 >dum.f90
mv dum.f90 ../ED2/utils/utils_f.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/canopy_struct_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/canopy_struct_dynamics.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/disturbance.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/disturbance.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/euler_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/euler_driver.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/events.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/events.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/farq_leuning.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/farq_leuning.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/fire.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/fire.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/forestry.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/forestry.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/growth_balive.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/growth_balive.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/heun_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/heun_driver.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/lsm_hyd.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/lsm_hyd.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/mortality.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/mortality.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/phenology_aux.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_aux.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/phenology_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_driv.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/photosyn_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/photosyn_driv.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/radiate_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/radiate_driver.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/reproduction.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/reproduction.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/rk4_derivs.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_derivs.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/rk4_driver.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_driver.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/rk4_integ_utils.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_integ_utils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/rk4_misc.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_misc.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/rk4_stepper.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_stepper.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/soil_respiration.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/soil_respiration.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/structural_growth.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/structural_growth.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/twostream_rad.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/twostream_rad.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/dynamics/vegetation_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/vegetation_dynamics.f90

sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/ed_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_init.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/ed_init_atm.F90 >dum.f90
mv dum.f90 ../ED2/init/ed_init_atm.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/ed_nbg_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_nbg_init.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/ed_params.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_params.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/ed_type_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_type_init.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/init_hydro_sites.f90 >dum.f90
mv dum.f90 ../ED2/init/init_hydro_sites.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/landuse_init.f90 >dum.f90
mv dum.f90 ../ED2/init/landuse_init.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/init/phenology_startup.f90 >dum.f90
mv dum.f90 ../ED2/init/phenology_startup.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/an_header.f90 >dum.f90
mv dum.f90 ../ED2/io/an_header.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/average_utils.f90 >dum.f90
mv dum.f90 ../ED2/io/average_utils.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_init_full_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_init_full_history.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_load_namelist.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_load_namelist.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_opspec.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_opspec.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_print.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_print.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_read_ed10_20_history.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed10_20_history.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_read_ed21_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed21_history.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/ed_xml_config.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_xml_config.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/edio.f90 >dum.f90
mv dum.f90 ../ED2/io/edio.f90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/h5_output.F90 >dum.f90
mv dum.f90 ../ED2/io/h5_output.F90
sed 's@use consts_coms@use ed_consts_coms@' <../ED2/io/leaf_database.f90 >dum.f90
mv dum.f90 ../ED2/io/leaf_database.f90

##
## hdf5_coms
##

sed 's@module hdf5_coms@module ed_hdf5_coms@' <../ED2/memory/ed_hdf5_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_hdf5_coms.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_max_dims.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_max_dims.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/grid_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/grid_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/pft_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/pft_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/rk4_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/rk4_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/soil_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/soil_coms.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_misc_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_misc_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/allometry.f90 >dum.f90
mv dum.f90 ../ED2/utils/allometry.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/ed_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_therm_lib.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/dp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/dp_therm_lib.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/sp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/sp_therm_lib.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/c34constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/c34constants.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/canopy_air_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_air_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/canopy_layer_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_layer_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/canopy_radiation_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_radiation_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/decomp_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/decomp_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/disturb_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/disturb_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_mem_alloc.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_alloc.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_mem_grid_dim_defs.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_grid_dim_defs.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_state_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_state_vars.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_var_tables.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_var_tables.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ed_work_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_work_vars.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/ename_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ename_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/fusion_fission_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/fusion_fission_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/hydrology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/hydrology_constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_constants.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/mem_polygons.f90 >dum.f90
mv dum.f90 ../ED2/memory/mem_polygons.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/met_driver_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/met_driver_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/optimiz_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/optimiz_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/phenology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/phenology_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/memory/physiology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/physiology_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/mpi/ed_mpass_init.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_mpass_init.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/mpi/ed_node_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_node_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/mpi/ed_para_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_coms.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/mpi/ed_para_init.F90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_init.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/driver/ed_1st.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_1st.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/driver/ed_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_driver.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/driver/ed_met_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_met_driver.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/driver/ed_model.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_model.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/budget_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/budget_utils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/charutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/charutils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/ed_dateutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_dateutils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/ed_filelist.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_filelist.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/ed_grid.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_grid.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/fatal_error.f90 >dum.f90
mv dum.f90 ../ED2/utils/fatal_error.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/fuse_fiss_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/fuse_fiss_utils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/great_circle.f90 >dum.f90
mv dum.f90 ../ED2/utils/great_circle.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/ed_hdf5_utils.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_hdf5_utils.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/invmondays.f90 >dum.f90
mv dum.f90 ../ED2/utils/invmondays.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/lapse.f90 >dum.f90
mv dum.f90 ../ED2/utils/lapse.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/libxml2f90.f90_pp.f90 >dum.f90
mv dum.f90 ../ED2/utils/libxml2f90.f90_pp.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/numutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/numutils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/radiate_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/radiate_utils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/rsys.F90 >dum.f90
mv dum.f90 ../ED2/utils/rsys.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/stable_cohorts.f90 >dum.f90
mv dum.f90 ../ED2/utils/stable_cohorts.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/update_derived_props.f90 >dum.f90
mv dum.f90 ../ED2/utils/update_derived_props.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/utils_c.c >dum.f90
mv dum.f90 ../ED2/utils/utils_c.c
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/utils/utils_f.f90 >dum.f90
mv dum.f90 ../ED2/utils/utils_f.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/canopy_struct_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/canopy_struct_dynamics.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/disturbance.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/disturbance.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/euler_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/euler_driver.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/events.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/events.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/farq_leuning.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/farq_leuning.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/fire.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/fire.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/forestry.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/forestry.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/growth_balive.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/growth_balive.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/heun_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/heun_driver.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/lsm_hyd.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/lsm_hyd.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/mortality.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/mortality.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/phenology_aux.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_aux.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/phenology_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_driv.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/photosyn_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/photosyn_driv.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/radiate_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/radiate_driver.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/reproduction.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/reproduction.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/rk4_derivs.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_derivs.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/rk4_driver.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_driver.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/rk4_integ_utils.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_integ_utils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/rk4_misc.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_misc.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/rk4_stepper.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_stepper.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/soil_respiration.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/soil_respiration.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/structural_growth.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/structural_growth.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/twostream_rad.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/twostream_rad.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/dynamics/vegetation_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/vegetation_dynamics.f90

sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/ed_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_init.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/ed_init_atm.F90 >dum.f90
mv dum.f90 ../ED2/init/ed_init_atm.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/ed_nbg_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_nbg_init.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/ed_params.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_params.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/ed_type_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_type_init.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/init_hydro_sites.f90 >dum.f90
mv dum.f90 ../ED2/init/init_hydro_sites.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/landuse_init.f90 >dum.f90
mv dum.f90 ../ED2/init/landuse_init.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/init/phenology_startup.f90 >dum.f90
mv dum.f90 ../ED2/init/phenology_startup.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/an_header.f90 >dum.f90
mv dum.f90 ../ED2/io/an_header.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/average_utils.f90 >dum.f90
mv dum.f90 ../ED2/io/average_utils.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_init_full_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_init_full_history.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_load_namelist.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_load_namelist.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_opspec.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_opspec.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_print.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_print.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_read_ed10_20_history.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed10_20_history.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_read_ed21_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed21_history.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/ed_xml_config.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_xml_config.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/edio.f90 >dum.f90
mv dum.f90 ../ED2/io/edio.f90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/h5_output.F90 >dum.f90
mv dum.f90 ../ED2/io/h5_output.F90
sed 's@use hdf5_coms@use ed_hdf5_coms@' <../ED2/io/leaf_database.f90 >dum.f90
mv dum.f90 ../ED2/io/leaf_database.f90
##
## hdf5_utils
##

sed 's@module hdf5_utils@module ed_hdf5_utils@' <../ED2/utils/ed_hdf5_utils.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_hdf5_utils.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_max_dims.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_max_dims.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/grid_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/grid_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/pft_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/pft_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/rk4_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/rk4_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/soil_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/soil_coms.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_misc_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_misc_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/allometry.f90 >dum.f90
mv dum.f90 ../ED2/utils/allometry.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/ed_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_therm_lib.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/dp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/dp_therm_lib.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/sp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/sp_therm_lib.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/c34constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/c34constants.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/canopy_air_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_air_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/canopy_layer_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_layer_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/canopy_radiation_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_radiation_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/decomp_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/decomp_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/disturb_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/disturb_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_mem_alloc.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_alloc.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_mem_grid_dim_defs.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_grid_dim_defs.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_state_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_state_vars.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_var_tables.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_var_tables.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_work_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_work_vars.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ename_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ename_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/fusion_fission_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/fusion_fission_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/ed_hdf5_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_hdf5_coms.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/hydrology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/hydrology_constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_constants.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/mem_polygons.f90 >dum.f90
mv dum.f90 ../ED2/memory/mem_polygons.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/met_driver_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/met_driver_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/optimiz_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/optimiz_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/phenology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/phenology_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/memory/physiology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/physiology_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/mpi/ed_mpass_init.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_mpass_init.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/mpi/ed_node_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_node_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/mpi/ed_para_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_coms.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/mpi/ed_para_init.F90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_init.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/driver/ed_1st.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_1st.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/driver/ed_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_driver.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/driver/ed_met_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_met_driver.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/driver/ed_model.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_model.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/budget_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/budget_utils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/charutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/charutils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/ed_dateutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_dateutils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/ed_filelist.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_filelist.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/ed_grid.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_grid.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/fatal_error.f90 >dum.f90
mv dum.f90 ../ED2/utils/fatal_error.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/fuse_fiss_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/fuse_fiss_utils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/great_circle.f90 >dum.f90
mv dum.f90 ../ED2/utils/great_circle.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/invmondays.f90 >dum.f90
mv dum.f90 ../ED2/utils/invmondays.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/lapse.f90 >dum.f90
mv dum.f90 ../ED2/utils/lapse.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/libxml2f90.f90_pp.f90 >dum.f90
mv dum.f90 ../ED2/utils/libxml2f90.f90_pp.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/numutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/numutils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/radiate_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/radiate_utils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/rsys.F90 >dum.f90
mv dum.f90 ../ED2/utils/rsys.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/stable_cohorts.f90 >dum.f90
mv dum.f90 ../ED2/utils/stable_cohorts.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/update_derived_props.f90 >dum.f90
mv dum.f90 ../ED2/utils/update_derived_props.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/utils_c.c >dum.f90
mv dum.f90 ../ED2/utils/utils_c.c
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/utils/utils_f.f90 >dum.f90
mv dum.f90 ../ED2/utils/utils_f.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/canopy_struct_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/canopy_struct_dynamics.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/disturbance.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/disturbance.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/euler_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/euler_driver.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/events.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/events.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/farq_leuning.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/farq_leuning.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/fire.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/fire.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/forestry.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/forestry.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/growth_balive.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/growth_balive.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/heun_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/heun_driver.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/lsm_hyd.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/lsm_hyd.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/mortality.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/mortality.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/phenology_aux.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_aux.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/phenology_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_driv.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/photosyn_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/photosyn_driv.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/radiate_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/radiate_driver.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/reproduction.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/reproduction.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/rk4_derivs.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_derivs.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/rk4_driver.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_driver.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/rk4_integ_utils.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_integ_utils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/rk4_misc.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_misc.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/rk4_stepper.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_stepper.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/soil_respiration.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/soil_respiration.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/structural_growth.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/structural_growth.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/twostream_rad.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/twostream_rad.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/dynamics/vegetation_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/vegetation_dynamics.f90

sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/ed_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_init.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/ed_init_atm.F90 >dum.f90
mv dum.f90 ../ED2/init/ed_init_atm.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/ed_nbg_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_nbg_init.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/ed_params.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_params.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/ed_type_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_type_init.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/init_hydro_sites.f90 >dum.f90
mv dum.f90 ../ED2/init/init_hydro_sites.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/landuse_init.f90 >dum.f90
mv dum.f90 ../ED2/init/landuse_init.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/init/phenology_startup.f90 >dum.f90
mv dum.f90 ../ED2/init/phenology_startup.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/an_header.f90 >dum.f90
mv dum.f90 ../ED2/io/an_header.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/average_utils.f90 >dum.f90
mv dum.f90 ../ED2/io/average_utils.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_init_full_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_init_full_history.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_load_namelist.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_load_namelist.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_opspec.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_opspec.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_print.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_print.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_read_ed10_20_history.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed10_20_history.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_read_ed21_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed21_history.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/ed_xml_config.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_xml_config.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/edio.f90 >dum.f90
mv dum.f90 ../ED2/io/edio.f90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/h5_output.F90 >dum.f90
mv dum.f90 ../ED2/io/h5_output.F90
sed 's@use hdf5_utils@use ed_hdf5_utils@' <../ED2/io/leaf_database.f90 >dum.f90
mv dum.f90 ../ED2/io/leaf_database.f90

##
## Convert therm_lib8 to something safe.
##

sed 's@module therm_lib8@module dp_therm_lib@' <../ED2/utils/dp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/dp_therm_lib.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/ed_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_therm_lib.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/ed_hdf5_utils.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_hdf5_utils.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/driver/ed_met_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_met_driver.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_max_dims.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_max_dims.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/grid_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/grid_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/pft_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/pft_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/rk4_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/rk4_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/soil_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/soil_coms.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_misc_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_misc_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/allometry.f90 >dum.f90
mv dum.f90 ../ED2/utils/allometry.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/sp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/sp_therm_lib.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/c34constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/c34constants.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/canopy_air_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_air_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/canopy_layer_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_layer_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/canopy_radiation_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_radiation_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/decomp_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/decomp_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/disturb_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/disturb_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_mem_alloc.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_alloc.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_mem_grid_dim_defs.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_grid_dim_defs.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_state_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_state_vars.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_var_tables.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_var_tables.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_work_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_work_vars.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ename_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ename_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/fusion_fission_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/fusion_fission_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/ed_hdf5_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_hdf5_coms.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/hydrology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/hydrology_constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_constants.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/mem_polygons.f90 >dum.f90
mv dum.f90 ../ED2/memory/mem_polygons.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/met_driver_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/met_driver_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/optimiz_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/optimiz_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/phenology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/phenology_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/memory/physiology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/physiology_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/mpi/ed_mpass_init.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_mpass_init.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/mpi/ed_node_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_node_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/mpi/ed_para_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_coms.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/mpi/ed_para_init.F90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_init.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/driver/ed_1st.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_1st.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/driver/ed_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_driver.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/driver/ed_model.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_model.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/budget_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/budget_utils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/charutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/charutils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/ed_dateutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_dateutils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/ed_filelist.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_filelist.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/ed_grid.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_grid.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/fatal_error.f90 >dum.f90
mv dum.f90 ../ED2/utils/fatal_error.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/fuse_fiss_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/fuse_fiss_utils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/great_circle.f90 >dum.f90
mv dum.f90 ../ED2/utils/great_circle.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/invmondays.f90 >dum.f90
mv dum.f90 ../ED2/utils/invmondays.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/lapse.f90 >dum.f90
mv dum.f90 ../ED2/utils/lapse.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/libxml2f90.f90_pp.f90 >dum.f90
mv dum.f90 ../ED2/utils/libxml2f90.f90_pp.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/numutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/numutils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/radiate_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/radiate_utils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/rsys.F90 >dum.f90
mv dum.f90 ../ED2/utils/rsys.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/stable_cohorts.f90 >dum.f90
mv dum.f90 ../ED2/utils/stable_cohorts.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/update_derived_props.f90 >dum.f90
mv dum.f90 ../ED2/utils/update_derived_props.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/utils_c.c >dum.f90
mv dum.f90 ../ED2/utils/utils_c.c
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/utils/utils_f.f90 >dum.f90
mv dum.f90 ../ED2/utils/utils_f.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/canopy_struct_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/canopy_struct_dynamics.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/disturbance.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/disturbance.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/euler_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/euler_driver.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/events.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/events.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/farq_leuning.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/farq_leuning.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/fire.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/fire.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/forestry.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/forestry.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/growth_balive.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/growth_balive.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/heun_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/heun_driver.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/lsm_hyd.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/lsm_hyd.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/mortality.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/mortality.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/phenology_aux.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_aux.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/phenology_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_driv.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/photosyn_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/photosyn_driv.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/radiate_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/radiate_driver.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/reproduction.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/reproduction.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/rk4_derivs.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_derivs.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/rk4_driver.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_driver.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/rk4_integ_utils.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_integ_utils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/rk4_misc.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_misc.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/rk4_stepper.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_stepper.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/soil_respiration.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/soil_respiration.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/structural_growth.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/structural_growth.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/twostream_rad.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/twostream_rad.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/dynamics/vegetation_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/vegetation_dynamics.f90

sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/ed_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_init.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/ed_init_atm.F90 >dum.f90
mv dum.f90 ../ED2/init/ed_init_atm.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/ed_nbg_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_nbg_init.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/ed_params.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_params.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/ed_type_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_type_init.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/init_hydro_sites.f90 >dum.f90
mv dum.f90 ../ED2/init/init_hydro_sites.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/landuse_init.f90 >dum.f90
mv dum.f90 ../ED2/init/landuse_init.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/init/phenology_startup.f90 >dum.f90
mv dum.f90 ../ED2/init/phenology_startup.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/an_header.f90 >dum.f90
mv dum.f90 ../ED2/io/an_header.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/average_utils.f90 >dum.f90
mv dum.f90 ../ED2/io/average_utils.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_init_full_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_init_full_history.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_load_namelist.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_load_namelist.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_opspec.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_opspec.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_print.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_print.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_read_ed10_20_history.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed10_20_history.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_read_ed21_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed21_history.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/ed_xml_config.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_xml_config.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/edio.f90 >dum.f90
mv dum.f90 ../ED2/io/edio.f90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/h5_output.F90 >dum.f90
mv dum.f90 ../ED2/io/h5_output.F90
sed 's@use therm_lib8@use dp_therm_lib@' <../ED2/io/leaf_database.f90 >dum.f90
mv dum.f90 ../ED2/io/leaf_database.f90

##
## sp_therm_lib
##

sed 's@module therm_lib@module sp_therm_lib@' <../ED2/utils/sp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/sp_therm_lib.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/ed_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_therm_lib.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/ed_hdf5_utils.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_hdf5_utils.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_max_dims.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_max_dims.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/grid_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/grid_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/pft_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/pft_coms.f90
 sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/rk4_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/rk4_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/soil_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/soil_coms.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_misc_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_misc_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/allometry.f90 >dum.f90
mv dum.f90 ../ED2/utils/allometry.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/dp_therm_lib.f90 >dum.f90
mv dum.f90 ../ED2/utils/dp_therm_lib.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/c34constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/c34constants.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/canopy_air_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_air_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/canopy_layer_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_layer_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/canopy_radiation_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/canopy_radiation_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/decomp_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/decomp_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/disturb_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/disturb_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_mem_alloc.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_alloc.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_mem_grid_dim_defs.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_mem_grid_dim_defs.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_state_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_state_vars.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_var_tables.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_var_tables.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_work_vars.f90 >dum.f90
mv dum.f90 ../ED2/memory/ed_work_vars.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ename_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/ename_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/fusion_fission_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/fusion_fission_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/ed_hdf5_coms.F90 >dum.f90
mv dum.f90 ../ED2/memory/ed_hdf5_coms.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/hydrology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/hydrology_constants.f90 >dum.f90
mv dum.f90 ../ED2/memory/hydrology_constants.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/mem_polygons.f90 >dum.f90
mv dum.f90 ../ED2/memory/mem_polygons.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/met_driver_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/met_driver_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/optimiz_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/optimiz_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/phenology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/phenology_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/memory/physiology_coms.f90 >dum.f90
mv dum.f90 ../ED2/memory/physiology_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/mpi/ed_mpass_init.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_mpass_init.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/mpi/ed_node_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_node_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/mpi/ed_para_coms.f90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_coms.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/mpi/ed_para_init.F90 >dum.f90
mv dum.f90 ../ED2/mpi/ed_para_init.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/driver/ed_1st.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_1st.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/driver/ed_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_driver.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/driver/ed_met_driver.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_met_driver.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/driver/ed_model.f90 >dum.f90
mv dum.f90 ../ED2/driver/ed_model.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/budget_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/budget_utils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/charutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/charutils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/ed_dateutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_dateutils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/ed_filelist.F90 >dum.f90
mv dum.f90 ../ED2/utils/ed_filelist.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/ed_grid.f90 >dum.f90
mv dum.f90 ../ED2/utils/ed_grid.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/fatal_error.f90 >dum.f90
mv dum.f90 ../ED2/utils/fatal_error.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/fuse_fiss_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/fuse_fiss_utils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/great_circle.f90 >dum.f90
mv dum.f90 ../ED2/utils/great_circle.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/invmondays.f90 >dum.f90
mv dum.f90 ../ED2/utils/invmondays.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/lapse.f90 >dum.f90
mv dum.f90 ../ED2/utils/lapse.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/libxml2f90.f90_pp.f90 >dum.f90
mv dum.f90 ../ED2/utils/libxml2f90.f90_pp.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/numutils.f90 >dum.f90
mv dum.f90 ../ED2/utils/numutils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/radiate_utils.f90 >dum.f90
mv dum.f90 ../ED2/utils/radiate_utils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/rsys.F90 >dum.f90
mv dum.f90 ../ED2/utils/rsys.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/stable_cohorts.f90 >dum.f90
mv dum.f90 ../ED2/utils/stable_cohorts.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/update_derived_props.f90 >dum.f90
mv dum.f90 ../ED2/utils/update_derived_props.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/utils_c.c >dum.f90
mv dum.f90 ../ED2/utils/utils_c.c
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/utils/utils_f.f90 >dum.f90
mv dum.f90 ../ED2/utils/utils_f.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/canopy_struct_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/canopy_struct_dynamics.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/disturbance.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/disturbance.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/euler_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/euler_driver.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/events.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/events.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/farq_leuning.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/farq_leuning.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/fire.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/fire.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/forestry.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/forestry.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/growth_balive.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/growth_balive.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/heun_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/heun_driver.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/lsm_hyd.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/lsm_hyd.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/mortality.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/mortality.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/phenology_aux.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_aux.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/phenology_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/phenology_driv.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/photosyn_driv.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/photosyn_driv.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/radiate_driver.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/radiate_driver.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/reproduction.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/reproduction.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/rk4_derivs.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_derivs.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/rk4_driver.F90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_driver.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/rk4_integ_utils.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_integ_utils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/rk4_misc.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_misc.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/rk4_stepper.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/rk4_stepper.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/soil_respiration.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/soil_respiration.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/structural_growth.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/structural_growth.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/twostream_rad.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/twostream_rad.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/dynamics/vegetation_dynamics.f90 >dum.f90
mv dum.f90 ../ED2/dynamics/vegetation_dynamics.f90

sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/ed_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_init.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/ed_init_atm.F90 >dum.f90
mv dum.f90 ../ED2/init/ed_init_atm.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/ed_nbg_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_nbg_init.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/ed_params.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_params.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/ed_type_init.f90 >dum.f90
mv dum.f90 ../ED2/init/ed_type_init.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/init_hydro_sites.f90 >dum.f90
mv dum.f90 ../ED2/init/init_hydro_sites.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/landuse_init.f90 >dum.f90
mv dum.f90 ../ED2/init/landuse_init.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/init/phenology_startup.f90 >dum.f90
mv dum.f90 ../ED2/init/phenology_startup.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/an_header.f90 >dum.f90
mv dum.f90 ../ED2/io/an_header.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/average_utils.f90 >dum.f90
mv dum.f90 ../ED2/io/average_utils.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_init_full_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_init_full_history.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_load_namelist.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_load_namelist.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_opspec.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_opspec.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_print.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_print.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_read_ed10_20_history.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed10_20_history.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_read_ed21_history.F90 >dum.f90
mv dum.f90 ../ED2/io/ed_read_ed21_history.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/ed_xml_config.f90 >dum.f90
mv dum.f90 ../ED2/io/ed_xml_config.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/edio.f90 >dum.f90
mv dum.f90 ../ED2/io/edio.f90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/h5_output.F90 >dum.f90
mv dum.f90 ../ED2/io/h5_output.F90
sed 's@use therm_lib@use sp_therm_lib@' <../ED2/io/leaf_database.f90 >dum.f90
mv dum.f90 ../ED2/io/leaf_database.f90
