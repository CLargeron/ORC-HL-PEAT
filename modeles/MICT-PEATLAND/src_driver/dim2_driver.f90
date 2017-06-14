PROGRAM driver
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_driver/dim2_driver.f90 $ 
!< $Date: 2012-10-26 16:32:10 +0200 (Fri, 26 Oct 2012) $
!< $Author: josefine.ghattas $
!< $Revision: 1042 $
!- IPSL (2006)
!-  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!---------------------------------------------------------------------
!- This PROGRAM is the driver of the dim2 version of SECHIBA. It's
!- main use is for PILPS type experiments. To run it you need to have
!- the following software :
!- - SECHIBA
!- - ioipsl
!- - netcdf
!- - F90 compiler
!- - tk (optional but very useful for configuring the model)
!- Obviously also the associated Makefiles.
!-
!- The forcing data needs to be in netCDF format and should
!- contain the following variables :
!- - Incoming SW radiation
!- - Incoming LW radiation
!- - Precipitation
!- - Air temperature at a reference level
!- - Air humidity at the same level
!- - wind at the same level
!- - surface pressure
!-
!- Once you have all this and compiled the model you can configure it.
!- Type make config and set all the values presented in the menu.
!- This tool will create a run.def which will be used by the model.
!- The run.def can also be edited by hand but it is more tedious.
!- Once this run.def is created the model is ready to run.
!---------------------------------------------------------------------
  USE netcdf
  USE ioipsl
  USE grid
  USE intersurf,  ONLY : intersurf_main
  USE constantes
  USE readdim2
  USE parallel
  USE timer
!-
  IMPLICIT NONE
!-
  INTEGER :: iim, jjm, llm
  INTEGER :: im, jm, lm, tm, is, force_id, itest, jtest
  REAL :: dt, dt_force, dt_rest, date0, date0_rest
  REAL :: zlflu
  REAL :: alpha
!-
  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
 & swdown, sinang, precip_rain, precip_snow, tair_obs, u, v, &
 & qair_obs, pb, lwdown, &
 & eair_obs, zlev_vec, zlevuv_vec, relax
!- Variables which are forcings for SECHIBA
  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
 & petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, &
 & for_u, for_v, for_swnet, for_swdown, for_sinang, for_lwdown, &
 & for_psurf, for_qair, for_tair, for_eair, &
 & for_ccanopy, for_rau
!!$, tmp_eair, tmp_tair, &
!!$ & tmp_qair, tmp_pb
!-
  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
 & for_contfrac, old_zlev, old_qair, old_eair, tsol_rad, vevapp, &
 & temp_sol_NEW, temp_sol_old, qsurf, dtdt, coastalflow, riverflow, &
 & fluxsens, fluxlat, emis, z0
!!$, epot_sol
!-
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: for_neighbours
!-
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: for_resolution
!-
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: albedo
  REAL, ALLOCATABLE, DIMENSION(:,:) :: albedo_vis
  REAL, ALLOCATABLE, DIMENSION(:,:) :: albedo_nir
!-
  INTEGER, ALLOCATABLE, DIMENSION(:) :: kindex
  REAL, ALLOCATABLE, DIMENSION(:,:) :: lon, lat
!-
  REAL :: old_tair
  REAL :: atmco2
  INTEGER :: nbindex
  REAL,DIMENSION(:), ALLOCATABLE :: lev, levuv
  REAL :: julian, ss
  INTEGER :: yy, mm, dd
!-
  LOGICAL :: relaxation, lower_wind
!-
  CHARACTER(LEN=80) :: filename, restname_in, restname_out
  CHARACTER(LEN=30) :: time_str, var_name
!-
  INTEGER :: it, istp, istp_old, rest_id, it_force
!-
  INTEGER :: split, split_start, nb_spread, for_offset
  INTEGER :: itau_dep, itau_dep_rest, itau_fin, itau_skip, itau_len

  INTEGER,DIMENSION(2) :: ml
!-
  LOGICAL :: last_CALL, first_CALL, debug
  LOGICAL :: no_inter, inter_lin, netrad_cons
!-
  ! to check variables passed to intersurf
  INTEGER :: ik
  INTEGER :: i,j
  LOGICAL, PARAMETER :: longprint = .FALSE.
  INTEGER :: ierr

  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
  & fluxsens_g,vevapp_g,old_zlev_g,old_qair_g,old_eair_g,for_rau_g, &
  & petAcoef_g, petBcoef_g,peqAcoef_g,peqBcoef_g,albedo_g,z0_g
  LOGICAL :: Flag

  REAL :: fill_init

  fill_init=REAL(nf90_fill_real,r_std)
  CALL ioconf_expval(val_exp)
!-
! Init parallelism

  CALL init_para(.FALSE.)
  CALL init_timer
  
! driver only for process root
  
!---------------------------------------------------------------------
!-
! set debug to have more information
!-
  !Config Key   = DEBUG_INFO
  !Config Desc  = Flag for debug information
  !Config If    = [-]
  !Config Def   = n
  !Config Help  = This option allows to switch on the output of debug
  !Config         information without recompiling the code.
  !Config Units = [FLAG] 
!-
  debug = .FALSE.
  CALL getin_p('DEBUG_INFO',debug)
!=====================================================================
!- 1.0 This section defines the general set-up of the experiment :
!-   - Forcing data to be used with its time step
!-   - Restart file to be used
!-   - The time step that will be used
!-   - The starting date of the simulation
!-   - Length of integration
!-   - Basic flags for SSIPSL
!=====================================================================
!- 1.1 Initialize the driving variables. It essentialy gets the mode
!-     used and the size of the driving variables.
!=====================================================================
  IF (debug) WRITE(numout,*) 'Reading name of the forcing file'
!- 
 !Config Key   = FORCING_FILE
 !Config Desc  = Name of file containing the forcing data
 !Config If    = [-]
 !Config Def   = forcing_file.nc
 !Config Help  = This is the name of the file which should be opened
 !Config         for reading the forcing data of the dim0 model.
 !Config         The format of the file has to be netCDF and COADS
 !Config         compliant.
 !Config Units = [FILE] 
!- 
  filename='forcing_file.nc'
  CALL getin_p('FORCING_FILE',filename)
!-
  IF (debug) WRITE(numout,*) 'Opening forcing file'
!-
! We call flininfo to obtain the dimensions
! of iim, jjm and others.
! This information will allow us to allocate all the space needed.
!-
  CALL forcing_info &
 &  (filename, iim, jjm, llm, tm, date0, dt_force, force_id) 
!-
  WRITE(numout,*) 'Out of flininfo : date0 ', date0, &
       'iim, jjm, llm, tm',iim,jjm,llm,tm,' dt_force ',dt_force
!-
  IF (debug) THEN
    WRITE(numout,*) 'Allocate memory for the driver :', iim, jjm, llm
  ENDIF
!-
  ALLOCATE(lev(llm),levuv(llm))
  ALLOCATE &
 & (swdown(iim,jjm), sinang(iim,jjm), precip_rain(iim,jjm), precip_snow(iim,jjm), &
 &  tair_obs(iim,jjm), u(iim,jjm), v(iim,jjm), qair_obs(iim,jjm), &
 &  pb(iim,jjm), lwdown(iim,jjm), &
 &  eair_obs(iim,jjm), zlev_vec(iim,jjm), zlevuv_vec(iim,jjm), relax(iim,jjm))
!- 
  ALLOCATE &
 & (petAcoef(iim,jjm), peqAcoef(iim,jjm), &
 &  petBcoef(iim,jjm), peqBcoef(iim,jjm), &
 &  cdrag(iim,jjm), for_u(iim,jjm), for_v(iim,jjm), &
 &  for_swnet(iim,jjm), for_swdown(iim,jjm), for_sinang(iim,jjm), for_lwdown(iim,jjm), &
 &  for_psurf(iim,jjm), for_qair(iim,jjm), for_tair(iim,jjm), &
 &  for_eair(iim,jjm), for_ccanopy(iim,jjm), for_rau(iim,jjm))
!!$, &
!!$ &  tmp_eair(iim,jjm), tmp_tair(iim,jjm), &
!!$ &  tmp_qair(iim,jjm), tmp_pb(iim,jjm))
!-
  ALLOCATE &
 & (for_contfrac(iim,jjm), for_neighbours(iim,jjm,8), for_resolution(iim,jjm,2), &
 &  old_zlev(iim,jjm), old_qair(iim,jjm), old_eair(iim,jjm), &
 &  tsol_rad(iim,jjm), vevapp(iim,jjm), &
 &  temp_sol_NEW(iim,jjm), temp_sol_old(iim,jjm), &
 &  dtdt(iim,jjm), coastalflow(iim,jjm), riverflow(iim,jjm), &
 &  fluxsens(iim,jjm), fluxlat(iim,jjm), emis(iim,jjm), &
 &  z0(iim,jjm), qsurf(iim,jjm))
!!$, epot_sol(iim,jjm)
!-
  ALLOCATE(albedo(iim,jjm,2))
  ALLOCATE(albedo_vis(iim,jjm),albedo_nir(iim,jjm))
!-
  ALLOCATE(kindex(iim*jjm))
  ALLOCATE(lon(iim,jjm), lat(iim,jjm))
!-- 
  lev(:) = 0.0
  levuv(:) = 0.0
  swdown(:,:) = fill_init
  precip_rain(:,:) = 0.0
  precip_snow(:,:) = 0.0
  tair_obs(:,:) = 0.0
  u(:,:) = fill_init
  v(:,:) = fill_init
  qair_obs(:,:) = fill_init
  pb(:,:) = fill_init
  lwdown(:,:) = fill_init
  eair_obs(:,:) = fill_init
  zlev_vec(:,:) = 0.0
  zlevuv_vec(:,:) = 0.0
  relax(:,:) = 0.0
  petAcoef(:,:) = 0.0
  peqAcoef(:,:) = 0.0
  petBcoef(:,:) = 0.0
  peqBcoef(:,:) = 0.0
  cdrag(:,:) = 0.0
  for_u(:,:) = fill_init
  for_v(:,:) = fill_init
  for_swnet(:,:) = fill_init
  for_swdown(:,:) = fill_init
  for_lwdown(:,:) = fill_init
  for_psurf(:,:) = fill_init
  for_qair(:,:) = fill_init
  for_tair(:,:) = fill_init
  for_eair(:,:) = fill_init
  for_ccanopy(:,:) = 0.0
  for_rau(:,:) = fill_init
  for_contfrac(:,:) = 0.0
  for_neighbours(:,:,:) = 0
  for_resolution(:,:,:) = 0.0
  old_zlev(:,:) = 0.0
  old_qair(:,:) = 0.0
  old_eair(:,:) = 0.0
  tsol_rad(:,:) = 0.0
  vevapp(:,:) = 0.0
  temp_sol_NEW(:,:) = fill_init
  temp_sol_old(:,:) = fill_init
  dtdt(:,:) = 0.0
  coastalflow(:,:) = 0.0
  riverflow(:,:) = 0.0
  fluxsens(:,:) = 0.0
  fluxlat(:,:) = 0.0
  emis(:,:) = 0.0
  z0(:,:) = fill_init
  qsurf(:,:) = 0.0
  albedo(:,:,:) = fill_init
  albedo_vis(:,:) = fill_init
  albedo_nir(:,:) = fill_init
  kindex(:) = 0
  lon(:,:) = 0.0
  lat(:,:) = 0.0
!-
! We need to know the grid.
! Then we can initialize the restart files, and then we
! can give read the restart files in the forcing subroutine.
!-
  CALL forcing_grid (iim,jjm,llm,lon,lat,lev,levuv,init_f=.FALSE.)
!=====================================================================
!- 1.2 Time step to be used.
!-     This is relative to the time step of the  forcing data
!=====================================================================
  IF ( .NOT. weathergen ) THEN
     !Config Key   = SPLIT_DT
     !Config Desc  = splits the timestep imposed by the forcing
     !Config If    = NOT(WEATHERGEN)
     !Config Def   = 12
     !Config Help  = With this value the time step of the forcing
     !Config         will be devided. In principle this can be run
     !Config         in explicit mode but it is strongly suggested
     !Config         to use the implicit method so that the
     !Config         atmospheric forcing has a smooth evolution.
     !Config Units = [-]
!-
     split = 12
     CALL getin_p('SPLIT_DT', split)
  ELSE
     split = 1
  ENDIF
!-
! The time step which is going to be used is computed
! from the forcing time step and the user input.
!-
  IF (split > 0) THEN
     dt = dt_force/split
  ELSE
     split = 1
     dt = dt_force
  ENDIF
!=====================================================================
!- 1.3 Initialize the restart file for the driver
!=====================================================================
  !Config Key   = RESTART_FILEIN
  !Config Desc  = Name of restart to READ for initial conditions
  !Config If    = [-]
  !Config Def   = NONE
  !Config Help  = This is the name of the file which will be opened
  !Config         to extract the initial values of all prognostic
  !Config         values of the model. This has to be a netCDF file.
  !Config         Not truly COADS compliant. NONE will mean that
  !Config         no restart file is to be expected.
  !Config Units = [FILE]
!-
  restname_in = 'NONE'
  CALL getin_p('RESTART_FILEIN',restname_in)
  if (debug) WRITE(numout,*) 'INPUT RESTART_FILE : ',TRIM(restname_in)
!-
  !Config Key   = RESTART_FILEOUT
  !Config Desc  = Name of restart files to be created by the driver
  !Config If    = [-]
  !Config Def   = driver_rest_out.nc
  !Config Help  = This variable give the  name for
  !Config         the restart files. The restart software within
  !Config         IOIPSL will add .nc if needed
  !Config Units = [FILE]
!-
  restname_out = 'driver_rest_out.nc'
  CALL getin_p('RESTART_FILEOUT', restname_out)
  if (debug) WRITE(numout,*) 'OUTPUT RESTART_FILE : ',TRIM(restname_out)
!-
! We choose some default values for the chronology
! of the restart file.
!-
  dt_rest = dt
  date0_rest = date0
!-
! Set default values for the start and end of the simulation
! in the forcing chronology.
!-
  itau_dep = 0
  itau_dep_rest = 0
  itau_fin = tm
!-
  CALL gather2D(lon,lon_g)
  CALL gather2D(lat,lat_g)
!-
  if (debug) WRITE(numout,*) &
 &   'Before Restart file initialization : itau_dep, itau_fin, date0, dt', &
 &   itau_dep, itau_fin, date0, dt
  
  IF (is_root_prc) THEN
     CALL restini &
 &       (restname_in, iim_g, jjm_g, lon_g, lat_g, llm, lev, &
 &        restname_out, itau_dep_rest, date0_rest, dt_rest, rest_id)
!-
     IF (debug) WRITE(numout,*) &
 &       'For Restart file initialization : itau_dep_rest, date0_rest, dt_rest', &
 &       itau_dep_rest, date0_rest, dt_rest
!-
!    MM When itau_dep_rest is not equal to itau_dep
     IF (itau_dep /= itau_dep_rest) THEN
        itau_dep = itau_dep_rest
        itau_fin = itau_dep+tm
     ENDIF
  ENDIF
  CALL bcast (itau_dep_rest)
  CALL bcast (itau_dep)
  CALL bcast (itau_fin)
  CALL bcast (date0_rest)
  CALL bcast (dt_rest)
!-
  WRITE(numout,*) &
 & 'Restart file initialized : itau_dep, itau_fin, date0_rest, dt_rest', &
 & itau_dep, itau_fin, date0_rest, dt_rest
!=====================================================================
!- 1.4 Here we do the first real reading of the driving. It only
!-     gets a few variables.
!=====================================================================
!-
  zlev_vec(:,:) = lev(1)
!-
!- In most forcing datas obtained from observations, the wind will be available
!- another height than temperature and humidity. We will thus need to lower the
!- wind to take this effect into account.
!-
  lower_wind=.FALSE.
  IF ( ABS(lev(1) - levuv(1)) > EPSILON(lev) ) THEN
     lower_wind=.TRUE.
     zlevuv_vec(:,:) = levuv(1)
  ELSE
     lower_wind=.FALSE.
     zlevuv_vec(:,:) = lev(1)
  ENDIF
!-
! prepares kindex table from the information obtained
! from the forcing data and reads some initial values for
! temperature, etc.
!- 
  kindex(1) = 1
!-
  CALL forcing_READ &
 &  (filename, rest_id, .TRUE., .FALSE., &
 &   0, itau_dep, 0, split, nb_spread, netrad_cons, &
 &   date0, dt_force, iim, jjm, lon, lat, zlev_vec, zlevuv_vec, tm, &
 &   swdown, sinang, precip_rain, precip_snow, tair_obs, &
 &   u, v, qair_obs, pb, lwdown, for_contfrac, for_neighbours, for_resolution, &
 &   for_swnet, eair_obs, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, for_ccanopy, &
 &   kindex, nbindex, force_id)
!-

  WRITE (numout,*) ">> Number of land points =",nbindex
  IF (nbindex == 0) THEN
     WRITE(numout,*) "Limits : (W,E / N,S)", limit_west, limit_east, &
          &                             limit_north, limit_south
     CALL ipslerr ( 3, 'dim2_driver','number of land points error.', &
          &         ' is zero !','stop driver')
  ENDIF

  jlandindex = (((kindex(1:nbindex)-1)/iim) + 1)
  ilandindex = (kindex(1:nbindex) - (jlandindex(1:nbindex)-1)*iim)
  IF (longprint) THEN
     WRITE(numout,*) "kindex of land points : ", kindex(1:nbindex)
     WRITE(numout,*) "index i of land points : ", ilandindex
     WRITE(numout,*) "index j of land points : ", jlandindex 
  ENDIF

  IF (is_watchout) THEN
     lower_wind = .FALSE.
  ENDIF
  IF (debug) WRITE(numout,*) "lower wind : ",lower_wind

!-
  im = iim; jm = jjm; lm = llm;
  IF ( (iim > 1).AND.(jjm > 1) ) THEN
    jtest = INT((kindex(INT(nbindex/2))-1)/iim)+1
    itest = MAX( 1, kindex(INT(nbindex/2))-(jtest-1)*iim )
  ELSE
    jtest = 1
    itest = 1
  ENDIF
  IF (debug) WRITE(numout,*) "test point in dim2_driver : ",itest,jtest
!-
  IF ((im /= iim) .AND. (jm /= jjm) .AND. (lm /= llm))  THEN
    WRITE(numout,*) ' dimensions are not good. Verify FILE :'
    WRITE(numout,*) ' filename = ',filename
    WRITE(numout,*) ' im, jm, lm lus         = ', im, jm, lm
    WRITE(numout,*) ' iim, jjm, llm demandes = ', iim, jjm, llm
    STOP 'dim2_driver'
  ENDIF  
!=====================================================================
!- 1.5  Configures the time-steps and other parameters 
!-      of the run.
!=====================================================================
!-
! If the time steping of the restart is different from the one
! of the forcing we need to convert the itau_dep into the
! chronology of the forcing. This ensures that the forcing will
! start at the date of the restart file. Obviously the itau_fin
! needs to be shifted as well !
!-
  IF ( (dt_rest /= dt_force).AND.(itau_dep > 1) ) THEN
    itau_dep = NINT((itau_dep*dt_rest )/dt_force)
    itau_fin = itau_dep+tm
    if (debug) WRITE(numout,*) &
 & 'The time steping of the restart is different from the one ',&
 & 'of the forcing we need to convert the itau_dep into the ',&
 & 'chronology of the forcing. This ensures that the forcing will ',&
 & 'start at the date of the restart file. Obviously the itau_fin ',&
 & 'needs to be shifted as well : itau_dep, itau_fin ', &
 & itau_dep, itau_fin
  ENDIF
!-
! Same things if the starting dates are not the same.
! Everything should look as if we had only one forcing file !
!-
  IF (date0 /= date0_rest) THEN
    WRITE(numout,*) 'date0_rest , date0 : ',date0_rest , date0
    for_offset = NINT((date0_rest-date0)*one_day/dt_force)
  ELSE
    for_offset = 0
  ENDIF
  WRITE(numout,*) 'OFFSET FOR THE data read :', for_offset

  CALL ioconf_startdate(date0_rest)
!-
  !Config Key   = TIME_SKIP
  !Config Desc  = Time in the forcing file at which the model is started.
  !Config If    = [-]
  !Config Def   = 0
  !Config Help  = This time give the point in time at which the model
  !Config         should be started. If exists, the date of the restart file is use.
  !Config         The FORMAT of this date can be either of the following :
  !Config         n   : time step n within the forcing file
  !Config         nS  : n seconds after the first time-step in the file
  !Config         nD  : n days after the first time-step
  !Config         nM  : n month after the first time-step (year of 365 days)
  !Config         nY  : n years after the first time-step (year of 365 days)
  !Config         Or combinations :
  !Config         nYmM: n years and m month
  !Config Units = [seconds, days, months, years]
!-
  itau_skip = 0
  WRITE(time_str,'(I10)') itau_skip
  CALL getin_p('TIME_SKIP', time_str)
!-
! Transform into itau
!-
  CALL tlen2itau (time_str, dt_force, date0, itau_skip)
!-
  itau_dep = itau_dep+itau_skip
!-
! We need to select the right position of the splited time steps.
!-
  istp = itau_dep*split+1
  IF (MOD(istp-1,split) /= 0) THEN
    split_start = MOD(istp-1,split)+1
  ELSE
    split_start = 1
  ENDIF
  istp_old = itau_dep_rest
!!$  it_force = MAX(MOD(itau_dep,INT(one_day*one_year/dt_force)),1)+itau_skip
!!$  WRITE(numout,*) 'itau_dep, it_force :', itau_dep, it_force
!-
  itau_len = itau_fin-itau_dep
!-
  !Config Key   = TIME_LENGTH
  !Config Desc  = Length of the integration in time.
  !Config If    = [-]
  !Config Def   = DEF
  !Config Help  = Length of integration. By default the entire length
  !Config         of the forcing is used. The FORMAT of this date can
  !Config         be either of the following :
  !Config         n   : time step n within the forcing file
  !Config         nS  : n seconds after the first time-step in the file
  !Config         nD  : n days after the first time-step
  !Config         nM  : n month after the first time-step (year of 365 days)
  !Config         nY  : n years after the first time-step (year of 365 days)
  !Config         Or combinations :
  !Config         nYmM: n years and m month
  !Config Units = [seconds, days, months, years]
!-
  WRITE(time_str,'(I10)') itau_len
  CALL getin_p('TIME_LENGTH', time_str)
!-
! Transform into itau
!-
  CALL tlen2itau (time_str, dt_force, date0, itau_len)
!-
  itau_fin = itau_dep+itau_len
!-
  WRITE(numout,*) '>> Time origine in the forcing file :', date0
  WRITE(numout,*) '>> Time origine in the restart file :', date0_rest
  WRITE(numout,*) '>> Simulate starts at forcing time-step : ', itau_dep
  WRITE(numout,*) '>> The splited time-steps start at (Sets the '
  WRITE(numout,*) '>>  chronology for the history and restart files):',istp
  WRITE(numout,*) '>> The time spliting starts at :', split_start
  WRITE(numout,*) '>> Simulation ends at forcing time-step: ', itau_fin
  WRITE(numout,*) '>> Length of the simulation is thus :', itau_len
  WRITE(numout,*) '>> Length of the forcing data is in time-steps : ', tm
  IF (tm < itau_len) THEN
     CALL ipslerr ( 2, 'dim2_driver','Length of the simulation is greater than.', &
 &         ' Length of the forcing data is in time-steps','verify TIME_LENGTH parameter.')
  ENDIF
  WRITE(numout,*) '>> Time steps : true, forcing and restart : ', dt,dt_force,dt_rest

!=====================================================================
!- 2.0 This section is going to define the details by which
!-     the input data is going to be used to force the
!-     land-surface scheme. The tasks are the following :
!-   - Is the coupling going to be explicit or implicit
!-   - Type of interpolation to be used.
!-   - At which height are the atmospheric forcings going to be used ?
!-   - Is a relaxation method going to be used on the forcing
!-   - Does net radiation in the interpolated data need to be conserved
!-   - How do we distribute precipitation.
!=====================================================================
  !Config Key   = RELAXATION
  !Config Desc  = method of forcing
  !Config If    = [-]
  !Config Def   = n
  !Config Help  = A method is proposed by which the first atmospheric
  !Config         level is not directly forced by observations but
  !Config         relaxed with a time constant towards observations.
  !Config         For the moment the methods tends to smooth too much
  !Config         the diurnal cycle and introduces a time shift.
  !Config         A more sophisticated method is needed.
  !Config Units = [FLAG]
!-
  relaxation = .FALSE.
  CALL getin_p('RELAXATION', relaxation)  
  IF ( relaxation ) THEN
     WRITE(numout,*) 'dim2_driver : The relaxation option is temporarily disabled as it does not'
     WRITE(numout,*) '              produce energy conservation in ORCHIDEE. If you intend to use it'
     WRITE(numout,*) '              you should validate it.'
     STOP 'dim2_driver'
!- 
     !Config Key   = RELAX_A
     !Config Desc  = Time constant of the relaxation layer
     !Config If    = RELAXATION
     !Config Def   = 1.0
     !Config Help  = The time constant associated to the atmospheric
     !Config         conditions which are going to be computed
     !Config         in the relaxed layer. To avoid too much
     !Config         damping the value should be larger than 1000.
     !Config Units = [days?]
!-
     alpha = 1000.0
     CALL getin_p('RELAX_A', alpha)
  ENDIF
!-
  no_inter = .TRUE.
  inter_lin = .FALSE.

  !Config Key   = NO_INTER 
  !Config Desc  = No interpolation IF split is larger than 1
  !Config If    = [-]
  !Config Def   = y
  !Config Help  = Choose IF you do not wish to interpolate linearly.
  !Config Units = [FLAG]
  CALL getin_p('NO_INTER', no_inter)

  !Config Key   = INTER_LIN
  !Config Desc  = Interpolation IF split is larger than 1
  !Config If    = [-]
  !Config Def   = n
  !Config Help  = Choose IF you wish to interpolate linearly.
  !Config Units = [FLAG]
  CALL getin_p('INTER_LIN', inter_lin)
!-
  IF (inter_lin) THEN
     no_inter = .FALSE.
  ELSE
     no_inter = .TRUE.
  ENDIF
!-
  IF (split == 1) THEN
    no_inter = .TRUE.
    inter_lin = .FALSE.
  ENDIF
!-
  !Config Key   = SPRED_PREC
  !Config Desc  = Spread the precipitation.
  !Config If    = [-]
  !Config Def   = 1
  !Config Help  = Spread the precipitation over n steps of the splited forcing 
  !Config         time step. This ONLY applied if the forcing time step has been splited.
  !Config         If the value indicated is greater than SPLIT_DT, SPLIT_DT is used for it.
  !Config Units = [-]
!-
  nb_spread = 1
  CALL getin_p('SPRED_PREC', nb_spread)  
  IF (nb_spread > split) THEN
    WRITE(numout,*) 'WARNING : nb_spread is too large it will be '
    WRITE(numout,*) '          set to the value of split'
    nb_spread = split
  ELSE IF (split == 1) THEN
    nb_spread = 1
  ENDIF
!-
  IF (inter_lin) THEN
     !Config Key   = NETRAD_CONS
     !Config Desc  = Conserve net radiation in the forcing
     !Config Def   = y
     !Config If    = INTER_LIN
     !Config Help  = When the interpolation is used the net radiation
     !Config         provided by the forcing is not conserved anymore.
     !Config         This should be avoided and thus this option should
     !Config         be TRUE (y).
     !Config         This option is not used for short-wave if the
     !Config         time-step of the forcing is longer than an hour.
     !Config         It does not make sense to try and reconstruct
     !Config         a diurnal cycle and at the same time conserve the 
     !Config         incoming solar radiation.
     !Config Units = [FLAG]
     !-
     netrad_cons = .TRUE.
     CALL getin_p('NETRAD_CONS', netrad_cons)
  ELSE
     netrad_cons = .FALSE.
  ENDIF
!=====================================================================
!- 3.0 Finaly we can prepare all the variables before launching
!-     the simulation !
!=====================================================================
! Initialize LOGICAL and the length of the integration
!-
  last_CALL = .FALSE.
  first_CALL = .TRUE.
!-
  temp_sol_NEW(:,:) = tp_00
!-  
  !Config Key   = ATM_CO2
  !Config Desc  = Value for atm CO2
  !Config If    = [-]
  !Config Def   = 350.
  !Config Help  = Value to prescribe the atm CO2.
  !Config         For pre-industrial simulations, the value is 286.2 .
  !Config         348. for 1990 year.
  !Config Units = [ppm]
!-
  atmco2=350.
  CALL getin_p('ATM_CO2',atmco2)
  for_ccanopy(:,:)=atmco2
!-
! Preparing for the implicit scheme.
! This means loading the prognostic variables from the restart file.
!-
  Flag=.FALSE.
  IF (is_root_prc) THEN
     ALLOCATE(fluxsens_g(iim_g,jjm_g))
     var_name= 'fluxsens'
     CALL restget &
 &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., fluxsens_g)
     IF (ALL(fluxsens_g(:,:) == val_exp)) THEN
        Flag=.TRUE.
     ELSE
        Flag=.FALSE.
     ENDIF
  ELSE
     ALLOCATE(fluxsens_g(0,1))
  ENDIF
  CALL bcast(Flag)
  IF (.NOT. Flag) THEN
     CALL scatter2D(fluxsens_g,fluxsens)
  ELSE
     fluxsens(:,:) = zero
  ENDIF
  DEALLOCATE(fluxsens_g)
!-
  IF (is_root_prc) THEN
     ALLOCATE(vevapp_g(iim_g,jjm_g))
     var_name= 'vevapp'
     CALL restget &
 &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., vevapp_g)
     IF (ALL(vevapp_g(:,:) == val_exp)) THEN
        Flag=.TRUE.
     ELSE
        Flag=.FALSE.
     ENDIF
  ELSE
     ALLOCATE(vevapp_g(0,1))
  ENDIF
  CALL bcast(Flag)
  IF (.NOT. Flag) THEN
     CALL scatter2D(vevapp_g,vevapp)
  ELSE
     vevapp(:,:) = zero
  ENDIF
  DEALLOCATE(vevapp_g)
!-
  IF (is_root_prc) THEN
     ALLOCATE(old_zlev_g(iim_g,jjm_g))
     var_name= 'zlev_old'
     CALL restget &
 &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., old_zlev_g)
     IF (ALL(old_zlev_g(:,:) == val_exp)) THEN
        Flag=.TRUE.
     ELSE
        Flag=.FALSE.
     ENDIF
  ELSE
     ALLOCATE(old_zlev_g(0,1))
  ENDIF
  CALL bcast(Flag)
  IF ( .NOT. Flag ) THEN
     CALL scatter2D(old_zlev_g,old_zlev)
  ELSE
     old_zlev(:,:)=zlev_vec(:,:)
  ENDIF
  DEALLOCATE(old_zlev_g)
!-
  IF (is_root_prc) THEN
     ALLOCATE(old_qair_g(iim_g,jjm_g))
     var_name= 'qair_old'
     CALL restget &
 &       (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., old_qair_g)
    IF (ALL(old_qair_g(:,:) == val_exp)) THEN
      Flag=.TRUE.
    ELSE
      Flag=.FALSE.
    ENDIF
  ELSE
     ALLOCATE(old_qair_g(0,1))
  ENDIF
  CALL bcast(Flag)
  IF ( .NOT. Flag ) THEN
     CALL scatter2D(old_qair_g,old_qair)
  ELSE
     old_qair(:,:) = qair_obs(:,:)
  ENDIF
  DEALLOCATE(old_qair_g)
!-
  IF (is_root_prc) THEN
     ALLOCATE(old_eair_g(iim_g,jjm_g))
     var_name= 'eair_old'
     CALL restget &
 &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., old_eair_g)
    IF (ALL(old_eair_g(:,:) == val_exp)) THEN
      Flag=.TRUE.
    ELSE
      Flag=.FALSE.
    ENDIF
  ELSE
     ALLOCATE(old_eair_g(0,1))
  ENDIF
  CALL bcast(Flag)
  IF ( .NOT. Flag ) THEN
     CALL scatter2D(old_eair_g,old_eair)
  ELSE
     DO ik=1,nbindex
        i=ilandindex(ik)
        j=jlandindex(ik)
        old_eair(i,j) = cp_air * tair_obs(i,j) + cte_grav*zlev_vec(i,j)
     ENDDO
  ENDIF
  DEALLOCATE(old_eair_g)
!-
! old density is also needed because we do not yet have the right pb
!-
!=> obsol�te ??!! (tjrs calcul� apr�s forcing_read) 
  IF (is_root_prc) THEN
     ALLOCATE(for_rau_g(iim_g,jjm_g))
     var_name= 'rau_old'
     CALL restget &
 &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., for_rau_g)
    IF (ALL(for_rau_g(:,:) == val_exp)) THEN
      Flag=.TRUE.
    ELSE
      Flag=.FALSE.
    ENDIF
  ELSE
     ALLOCATE(for_rau_g(0,1))
  ENDIF
  CALL bcast(Flag)
  IF ( .NOT. Flag ) THEN
     CALL scatter2D(for_rau_g,for_rau)
  ELSE
     DO ik=1,nbindex
        i=ilandindex(ik)
        j=jlandindex(ik)
        for_rau(i,j) = pb(i,j)/(cte_molr*(tair_obs(i,j)))
     ENDDO
  ENDIF
  DEALLOCATE(for_rau_g)
!-
! For this variable the restart is extracted by SECHIBA
!-
  temp_sol_NEW(:,:) = tair_obs(:,:)
!- 
  if (.NOT. is_watchout) THEN
!-
!   This does not yield a correct restart in the case of relaxation
!-
     IF (is_root_prc) THEN
        ALLOCATE(petAcoef_g(iim_g,jjm_g))
        var_name= 'petAcoef'
        CALL restget &
 &           (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., petAcoef_g)
        IF (ALL(petAcoef_g(:,:) == val_exp)) THEN
           Flag=.TRUE.
        ELSE
           Flag=.FALSE.
        ENDIF
     ELSE
        ALLOCATE(petAcoef_g(0,1))
     ENDIF
     CALL bcast(Flag)
     IF ( .NOT. Flag ) THEN
        CALL scatter2D(petAcoef_g,petAcoef)
     ELSE
        petAcoef(:,:) = zero
     ENDIF
     DEALLOCATE(petAcoef_g)
!--
     IF (is_root_prc) THEN
        ALLOCATE(petBcoef_g(iim_g,jjm_g))
        var_name= 'petBcoef'
        CALL restget &
 &           (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., petBcoef_g)
        IF (ALL(petBcoef_g(:,:) == val_exp)) THEN
           Flag=.TRUE.
        ELSE
           Flag=.FALSE.
        ENDIF
     ELSE
        ALLOCATE(petBcoef_g(0,1))
     ENDIF
     CALL bcast(Flag)
     IF ( .NOT. Flag ) THEN
        CALL scatter2D(petBcoef_g,petBcoef)
     ELSE
        petBcoef(:,:) = old_eair(:,:)
     ENDIF
     DEALLOCATE(petBcoef_g)
!--
     IF (is_root_prc) THEN
        ALLOCATE(peqAcoef_g(iim_g,jjm_g))
        var_name= 'peqAcoef'
        CALL restget &
 &           (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., peqAcoef_g)
        IF (ALL(peqAcoef_g(:,:) == val_exp)) THEN
           Flag=.TRUE.
        ELSE
           Flag=.FALSE.
        ENDIF
     ELSE
        ALLOCATE(peqAcoef_g(0,1))
     ENDIF
     CALL bcast(Flag)
     IF ( .NOT. Flag ) THEN
        CALL scatter2D(peqAcoef_g,peqAcoef)
     ELSE
        peqAcoef(:,:) = zero
     ENDIF
     DEALLOCATE(peqAcoef_g)
!--
     IF (is_root_prc) THEN
        ALLOCATE(peqBcoef_g(iim_g,jjm_g))
        var_name= 'peqBcoef'
        CALL restget &
 &           (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., peqBcoef_g)
        IF (ALL(peqBcoef_g(:,:) == val_exp)) THEN
           Flag=.TRUE.
        ELSE
           Flag=.FALSE.
        ENDIF
     ELSE
        ALLOCATE(peqBcoef_g(0,1))
     ENDIF
     CALL bcast(Flag)
     IF ( .NOT. Flag ) THEN
        CALL scatter2D(peqBcoef_g,peqBcoef)
     ELSE
        peqBcoef(:,:) = old_qair(:,:)
     ENDIF
     DEALLOCATE(peqBcoef_g)
  ENDIF
!-
! And other variables which need initial variables. These variables
! will get properly initialized by ORCHIDEE when it is called for 
! the first time.
!-
  albedo(:,:,:) = 0.13
  emis(:,:) = 1.0
  z0(:,:) = 0.1
!--
!=====================================================================
!- 4.0 START THE TIME LOOP
!=====================================================================

  CALL init_timer   

  it = itau_dep+1
  DO WHILE ( it <= itau_fin )
!----
    it_force = it+for_offset
    IF (it_force < 0) THEN
      WRITE(numout,*) 'The day is not in the forcing file :', &
 &               it_force, it, for_offset
      STOP 'dim2_driver'
    ENDIF
!!$    IF (it_force > itau_dep+tm) THEN
!!$       WRITE(numout,*) 'ERROR : more time-steps than data'
!!$       WRITE(numout,*) 'it_force : ', it_force, ' itau_dep+tm : ', itau_dep+tm
!!$      STOP 'dim2_driver'
!!$    ENDIF 
!---
    is=split_start
    DO WHILE ( is <= split )
!-----
      julian = itau2date(istp, date0_rest, dt)
      CALL ju2ymds(julian, yy, mm, dd, ss)
      IF (debug) THEN
         WRITE(numout,*) "=============================================================="
         WRITE(numout,"('New iteration at date : ',I4,'-',I2.2,'-',I2.2,':',F8.4)") &
              &               yy,mm,dd,ss/3600.
#ifdef CPP_PARA
         IF (is_root_prc) THEN
            WRITE(*,*) "=============================================================="
            WRITE(*,"('New iteration at date : ',I4,'-',I2.2,'-',I2.2,':',F8.4)") &
              &               yy,mm,dd,ss/3600.
         ENDIF
#endif
      ENDIF
!----- 
      IF ( (it == itau_fin).AND.(is == split) ) THEN
        last_CALL = .TRUE.
      ENDIF
!-----
      IF (debug) WRITE(numout,*) 'Into forcing_read'
!-----
      CALL forcing_READ &
 &      (filename, rest_id, .FALSE., last_call, &
 &       it_force, istp, is, split, nb_spread, netrad_cons, &
 &       date0_rest, dt_force, iim, jjm, lon, lat, zlev_vec, zlevuv_vec, tm, &
 &       swdown, sinang, precip_rain, precip_snow, tair_obs, &
 &       u, v, qair_obs, pb, for_lwdown, for_contfrac, for_neighbours, for_resolution, &
 &       for_swnet, eair_obs, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, for_ccanopy, &
 &       kindex, nbindex, force_id)

!-----
!---- SECHIBA expects surface pressure in hPa
!-----
      for_psurf(:,:)  = pb(:,:)/100.

      IF (longprint) THEN
         WRITE(numout,*) "dim2_driver 0 ",it_force 
         WRITE(numout,*) ">> Index of land points =",kindex(1:nbindex)
         WRITE(numout,*) "Lowest level wind speed North = ", &
              & (/ ( u(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Lowest level wind speed East = ", &
              & (/ ( v(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "z0            ; Surface roughness = ", &
              & (/ ( z0(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Height of first layer = ", &
              & (/ ( zlev_vec(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Lowest level specific humidity = ", &
              & (/ ( qair_obs(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Rain precipitation = ", &
              & (/ ( precip_rain(ilandindex(ik), jlandindex(ik))*dt,ik=1,nbindex ) /)
         WRITE(numout,*) "Snow precipitation = ", &
              & (/ ( precip_snow(ilandindex(ik), jlandindex(ik))*dt,ik=1,nbindex ) /)
         WRITE(numout,*) "Down-welling long-wave flux = ", &
              & (/ ( for_lwdown(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Net surface short-wave flux = ", &
              & (/ ( for_swnet(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Downwelling surface short-wave flux = ", &
              & (/ ( swdown(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Air temperature in Kelvin = ", &
              & (/ ( tair_obs(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Air potential energy = ", &
              & (/ ( eair_obs(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "CO2 concentration in the canopy = ", &
              & (/ ( for_ccanopy(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Coeficients A from the PBL resolution = ", &
              & (/ ( petAcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "One for T and another for q = ", &
              & (/ ( peqAcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Coeficients B from the PBL resolution = ", &
              & (/ ( petBcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "One for T and another for q = ", &
              & (/ ( peqBcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Cdrag = ", &
              & (/ ( cdrag(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /) 
         WRITE(numout,*) "CO2 concentration in the canopy = ", &
              & (/ ( for_ccanopy(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /) 
         WRITE(numout,*) "Lowest level pressure = ", &
              & (/ ( for_psurf(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Geographical coordinates lon = ", &
              & (/ (  lon(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Geographical coordinates lat = ", &
              & (/ (  lat(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Fraction of continent in the grid = ", &
              & (/ ( for_contfrac(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
      ENDIF
!-----
!---- Prepare : tmp_qair, tmp_eair, tmp_tair, tmp_pb
!---- and     : for_u, for_v, for_lwdown, for_swnet, for_swdown
!---- All the work is done in forcing_read
!---- We do not need the no_inter options as we will
!---- allways interpolate
!----- 
      IF (debug) WRITE(numout,*) 'Prepare the atmospheric forcing'
!----- 
      IF (.NOT. is_watchout) THEN
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            eair_obs(i,j) = cp_air*tair_obs(i,j)+cte_grav*zlev_vec(i,j)
            for_swnet(i,j) = (1.-(albedo(i,j,1)+albedo(i,j,2))/2.)*swdown(i,j)
         ENDDO
      ENDIF
      DO ik=1,nbindex
         i=ilandindex(ik)
         j=jlandindex(ik)
         for_swdown(i,j) = swdown(i,j)
         for_sinang(i,j) = sinang(i,j)
      ENDDO
!----- 
!---- Computing the buffer zone !
!----- 
      IF (relaxation) THEN
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_qair(i,j) = peqAcoef(i,j)*(-1.) * vevapp(i,j)*dt+peqBcoef(i,j)
!-------
            for_eair(i,j) = petAcoef(i,j)*(-1.) * fluxsens(i,j)+petBcoef(i,j)
!-------
         ENDDO
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_tair(i,j) = (for_eair(i,j) - cte_grav*zlev_vec(i,j))/cp_air
!-------
!!$        if (.NOT. is_watchout) &
!!$             epot_sol(:,:) = cp_air*temp_sol_NEW(:,:)
!-------
         ENDDO
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_rau(i,j) = pb(i,j) / (cte_molr*for_tair(i,j))
!-------
            relax(i,j) = for_rau(i,j)*alpha 
         ENDDO
         zlflu = lev(1)/2.0*dt
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            peqAcoef(i,j) = 1.0/(zlflu+relax(i,j))
            peqBcoef(i,j) = (relax(i,j) * qair_obs(i,j)/(zlflu+relax(i,j))) + & 
                 & for_qair(i,j)/(1.0+relax(i,j)/zlflu)
         ENDDO
!-------
!        relax(:,:) = for_rau(:,:)*alpha 
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            petAcoef(i,j) = 1.0/(zlflu+relax(i,j))
            petBcoef(i,j) = ( relax(i,j) * eair_obs(i,j) / (zlflu+relax(i,j)) ) &
                 & + for_eair(i,j)/(1.0+relax(i,j)/zlflu)
         ENDDO
      ELSE
         for_qair(:,:) = fill_init
         for_eair(:,:) = fill_init
         for_tair(:,:) = fill_init
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_qair(i,j) = qair_obs(i,j)
            for_eair(i,j) = eair_obs(i,j)
            for_tair(i,j) = tair_obs(i,j)
         ENDDO
!-------
!!$        if (.NOT. is_watchout) &
!!$             epot_sol(:,:) =  cp_air*temp_sol_NEW(:,:)
!-------
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_rau(i,j) = pb(i,j) / (cte_molr*for_tair(i,j))
         ENDDO
!-------
         IF (.NOT. is_watchout) THEN
           petAcoef(:,:) = 0.0
           peqAcoef(:,:) = 0.0
           DO ik=1,nbindex
              i=ilandindex(ik)
              j=jlandindex(ik)
              petBcoef(i,j) = eair_obs(i,j)
              peqBcoef(i,j) = qair_obs(i,j)
           ENDDO
        ENDIF
      ENDIF
!-----
      IF (.NOT. is_watchout) &
           cdrag(:,:)  = 0.0
      for_ccanopy(:,:)=atmco2
!-----
!---- SECHIBA expects wind, temperature and humidity at the same height.
!---- If this is not the case then we need to correct for that.
!-----
      IF ( lower_wind .AND. .NOT. is_watchout ) THEN
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_u(i,j) = u(i,j)*LOG(zlev_vec(i,j)/z0(i,j)) / &
                 &              LOG(zlevuv_vec(i,j)/z0(i,j))
            for_v(i,j) = v(i,j)*LOG(zlev_vec(i,j)/z0(i,j)) / &
                 &              LOG(zlevuv_vec(i,j)/z0(i,j))
         ENDDO
      ELSE
         DO ik=1,nbindex
            i=ilandindex(ik)
            j=jlandindex(ik)
            for_u(i,j) = u(i,j)
            for_v(i,j) = v(i,j)
         ENDDO
      ENDIF
!-----
!---- Prepare the other variables WITH the special CASE
!---- of splited time steps
!----
!---- PRINT input value for debug
!-----
      IF (debug) THEN
        WRITE(numout,*) ' >>>>>> time step it_force = ',it_force
        WRITE(numout,*) &
 &       ' tair, qair, eair = ', &
 &       for_tair(itest,jtest),for_qair(itest,jtest), &
 &       for_eair(itest,jtest)
        WRITE(numout,*) &
 &       ' OBS : tair, qair, eair = ', &
 &       tair_obs(itest,jtest),qair_obs(itest,jtest), &
 &       eair_obs(itest,jtest)
        WRITE(numout,*) ' u et v = ',for_u(itest,jtest),for_v(itest,jtest)
        WRITE(numout,*) ' precip rain et snow = ', &
        & precip_rain(itest,jtest),precip_snow(itest,jtest)
        WRITE(numout,*) ' lwdown et swnet = ', &
        & for_lwdown(itest,jtest),for_swnet(itest,jtest)
        WRITE(numout,*) ' petAcoef et peqAcoef = ', &
        & petAcoef(itest,jtest), peqAcoef(itest,jtest)
        WRITE(numout,*) ' petBcoef et peqAcoef = ', &
        & petBcoef(itest,jtest),peqBcoef(itest,jtest)
        WRITE(numout,*) ' zlev = ',zlev_vec(itest,jtest)
      ENDIF
!-----
      IF (first_CALL) THEN

        DO ik=1,nbindex
           i=ilandindex(ik)
           j=jlandindex(ik)
           for_swdown(i,j) = swdown(i,j)
           for_sinang(i,j) = sinang(i,j)
        ENDDO
        IF (longprint) THEN
           WRITE(numout,*) "dim2_driver first_CALL ",it_force 
           WRITE(numout,*) ">> Index of land points =",kindex(1:nbindex)
           WRITE(numout,*) "Lowest level wind speed North = ", &
             &     (/ ( for_u(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Lowest level wind speed East = ", &
             &     (/ ( for_v(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "z0            ; Surface roughness = ", &
             &     (/ ( z0(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Height of first layer = ", &
             &     (/ ( zlev_vec(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Lowest level specific humidity = ", &
             &     (/ ( for_qair(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Rain precipitation = ", &
             &     (/ ( precip_rain(ilandindex(ik), jlandindex(ik))*dt,ik=1,nbindex ) /)
           WRITE(numout,*) "Snow precipitation = ", &
             &     (/ ( precip_snow(ilandindex(ik), jlandindex(ik))*dt,ik=1,nbindex ) /)
           WRITE(numout,*) "Down-welling long-wave flux = ", &
             &     (/ ( for_lwdown(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Net surface short-wave flux = ", &
             &     (/ ( for_swnet(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Downwelling surface short-wave flux = ", &
             &     (/ ( for_swdown(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Air temperature in Kelvin = ", &
             &     (/ ( for_tair(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Air potential energy = ", &
             &     (/ ( for_eair(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "CO2 concentration in the canopy = ", &
             &     (/ ( for_ccanopy(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Coeficients A from the PBL resolution = ", &
             &     (/ ( petAcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "One for T and another for q = ", &
             &     (/ ( peqAcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Coeficients B from the PBL resolution = ", &
             &     (/ ( petBcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "One for T and another for q = ", &
             &     (/ ( peqBcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Cdrag = ", &
             &     (/ ( cdrag(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Lowest level pressure = ", &
             &     (/ ( for_psurf(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Geographical coordinates lon = ", &
             &     (/ (  lon(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Geographical coordinates lat = ", &
             &     (/ (  lat(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Fraction of continent in the grid = ", &
             &     (/ ( for_contfrac(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
        ENDIF
!-------
!------ CALL sechiba to initialize fields
!------ and have some initial results: emis, albedo, z0
!-------
        CALL intersurf_main &
 &        (istp_old, iim, jjm, nbindex, kindex, dt, &
 &         first_CALL, .FALSE., lon, lat, for_contfrac, for_neighbours, for_resolution, date0_rest, &
!       first level conditions
 &         zlev_vec, for_u, for_v, &
 &         for_qair, for_tair, for_eair, for_ccanopy, &
!       Variables for the implicit coupling
 &         cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
!       Rain, snow, radiation and surface pressure
 &         precip_rain, precip_snow, &
 &         for_lwdown, for_swnet, for_swdown, for_psurf, &
!       Output : Fluxes
 &         vevapp, fluxsens, fluxlat, coastalflow, riverflow,  &
!       Surface temperatures and surface properties
 &         tsol_rad, temp_sol_NEW, qsurf, albedo, emis, z0, & 
!       VOC : radiation
 &         for_sinang)

        !
        first_CALL = .FALSE.
        !
        ! Get Restart values for albedo and z0, 
        ! as they modify forcing variables swnet and wind.
!-------
        ! albedo 
        IF (is_root_prc) THEN
           ALLOCATE(albedo_g(iim_g,jjm_g))
        ELSE
           ALLOCATE(albedo_g(0,1))
        ENDIF
        !
        IF (is_root_prc) THEN
           var_name= 'albedo_vis'
           CALL restget &
                &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., albedo_g)
           IF (ALL(albedo_g(:,:) == val_exp)) THEN
              Flag=.TRUE.
           ELSE
              Flag=.FALSE.
           ENDIF
        ENDIF
        CALL bcast(Flag)
        IF ( .NOT. Flag ) THEN
           CALL scatter2D(albedo_g,albedo_vis)
           albedo(:,:,1)=albedo_vis(:,:)
        ELSE
           albedo_vis(:,:)=albedo(:,:,1)
        ENDIF
        !
        IF (is_root_prc) THEN
           var_name= 'albedo_nir'
           CALL restget &
                &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., albedo_g)
           IF (ALL(albedo_g(:,:) == val_exp)) THEN
              Flag=.TRUE.
           ELSE
              Flag=.FALSE.
           ENDIF
        ENDIF
        CALL bcast(Flag)
        IF ( .NOT. Flag ) THEN
           CALL scatter2D(albedo_g,albedo_nir)
           albedo(:,:,2)=albedo_nir(:,:)
        ELSE
           albedo_nir(:,:)=albedo(:,:,2)
        ENDIF
        !
        DEALLOCATE(albedo_g)
        !--
        ! z0 
        IF (is_root_prc) THEN
           ALLOCATE(z0_g(iim_g,jjm_g))
           var_name= 'z0'
           CALL restget &
                &        (rest_id, var_name, iim_g, jjm_g, 1, istp_old, .TRUE., z0_g)
           IF (ALL(z0_g(:,:) == val_exp)) THEN
              Flag=.TRUE.
           ELSE
              Flag=.FALSE.
           ENDIF
        ELSE
           ALLOCATE(z0_g(0,1))
        ENDIF
        CALL bcast(Flag)
        IF (.NOT. Flag) &
             CALL scatter2D(z0_g,z0)
        DEALLOCATE(z0_g)
!-------
        DO ik=1,nbindex
           i=ilandindex(ik)
           j=jlandindex(ik)
           temp_sol_old(i,j) = temp_sol_NEW(i,j)
           for_swnet(i,j) = (1.- (albedo(i,j,1)+albedo(i,j,2))/2.)*swdown(i,j)
           for_swdown(i,j) = swdown(i,j)
           for_sinang(i,j) = sinang(i,j)
        ENDDO
!
!     MM : z0 have been modified then we must lower the wind again
!-----
!---- SECHIBA expects wind, temperature and humidity at the same height.
!---- If this is not the case then we need to correct for that.
!-----
        IF ( lower_wind .AND. .NOT. is_watchout ) THEN
           DO ik=1,nbindex
              i=ilandindex(ik)
              j=jlandindex(ik)
              for_u(i,j) = u(i,j) * LOG(zlev_vec(i,j)/z0(i,j)) / &
                   &                LOG(zlevuv_vec(i,j)/z0(i,j))
              for_v(i,j) = v(i,j) * LOG(zlev_vec(i,j)/z0(i,j)) / &
                   &                LOG(zlevuv_vec(i,j)/z0(i,j))
           ENDDO
        ELSE
           DO ik=1,nbindex
              i=ilandindex(ik)
              j=jlandindex(ik)
              for_u(i,j) = u(i,j)
              for_v(i,j) = v(i,j)
           ENDDO
        ENDIF
!-----
!---- PRINT input value after first_CALL for debug
!-----
        IF (debug) THEN
           WRITE(numout,*) ' >>>>>> after first_CALL = ',first_CALL
           WRITE(numout,*) ' u et v = ',for_u(itest,jtest),for_v(itest,jtest)
           WRITE(numout,*) ' swnet = ', for_swnet(itest,jtest)
        ENDIF
!-------
        IF (longprint) THEN
           WRITE(numout,*) "dim2_driver first_CALL outputs"
           !       Output : Fluxes
           WRITE(numout,*) "vevapp        ; Total of evaporation = ", &
             &     (/ ( vevapp(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Sensible heat flux = ", &
             &     (/ ( fluxsens(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "Latent heat flux = ", &
             &     (/ ( fluxlat(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "coastalflow   ; Diffuse flow of water into the ocean (m^3/dt) = ", &
             &     (/ ( coastalflow(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "riverflow     ; Largest rivers flowing into the ocean (m^3/dt) = ", &
             &     (/ ( riverflow(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           !       Surface temperatures and surface properties
           WRITE(numout,*) "tsol_rad      ; Radiative surface temperature = ", &
             &     (/ ( tsol_rad(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "temp_sol_new  ; New soil temperature = ", &
             &     (/ ( temp_sol_NEW(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "qsurf         ; Surface specific humidity = ", &
             &     (/ ( qsurf(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "albedoVIS = ", &
             &     (/ ( albedo(ilandindex(ik), jlandindex(ik), 1),ik=1,nbindex ) /)
           WRITE(numout,*) "albedoNIR = ", &
             &     (/ ( albedo(ilandindex(ik), jlandindex(ik), 2),ik=1,nbindex ) /)
           WRITE(numout,*) "emis          ; Emissivity = ", &
             &     (/ ( emis(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
           WRITE(numout,*) "z0            ; Surface roughness = ", &
             &     (/ ( z0(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
        ENDIF
!-------
        IF (debug) THEN
          WRITE(numout,*) &
 &         ' OUT rest : z0, albedoVIS, albedoNIR, emis = ', &
 &         z0(itest,jtest),albedo(itest,jtest,1), &
 &                         albedo(itest,jtest,2),emis(itest,jtest)
          WRITE(numout,*) ' OUT rest : coastal and river flow = ', &
 &         coastalflow(itest,jtest), riverflow(itest,jtest)
          WRITE(numout,*) ' OUT rest : tsol_rad, vevapp = ', &
 &         tsol_rad(itest,jtest), vevapp(itest,jtest)
          WRITE(numout,*) ' OUT rest : temp_sol_new =', &
 &         temp_sol_NEW(itest,jtest)
        ENDIF
    
        CALL barrier_para
        CALL start_timer(timer_global)
        CALL start_timer(timer_mpi)
   
      ENDIF
!-----
!---- Calling SECHIBA and doing the number crunching.
!---- Note that for the first time step SECHIBA is called twice.
!----
!---- All H_2O fluxes are now in Kg/m^2s
!-----
      IF (longprint) THEN
         WRITE(numout,*) "dim2_driver ",it_force 
         WRITE(numout,*) ">> Index of land points =",kindex(1:nbindex)
         WRITE(numout,*) "Lowest level wind speed North = ", &
           &     (/ ( for_u(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Lowest level wind speed East = ", &
           &     (/ ( for_v(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "z0            ; Surface roughness = ", &
           &     (/ ( z0(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Height of first layer = ", &
           &     (/ ( zlev_vec(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Lowest level specific humidity = ", &
           &     (/ ( for_qair(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Rain precipitation = ", &
           &     (/ ( precip_rain(ilandindex(ik), jlandindex(ik))*dt,ik=1,nbindex ) /)
         WRITE(numout,*) "Snow precipitation = ", &
           &     (/ ( precip_snow(ilandindex(ik), jlandindex(ik))*dt,ik=1,nbindex ) /)
         WRITE(numout,*) "Down-welling long-wave flux = ", &
           &     (/ ( for_lwdown(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Net surface short-wave flux = ", &
           &     (/ ( for_swnet(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Downwelling surface short-wave flux = ", &
           &     (/ ( for_swdown(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Air temperature in Kelvin = ", &
           &     (/ ( for_tair(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Air potential energy = ", &
           &     (/ ( for_eair(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "CO2 concentration in the canopy = ", &
           &     (/ ( for_ccanopy(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Coeficients A from the PBL resolution = ", &
           &     (/ ( petAcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "One for T and another for q = ", &
           &     (/ ( peqAcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Coeficients B from the PBL resolution = ", &
           &     (/ ( petBcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "One for T and another for q = ", &
           &     (/ ( peqBcoef(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Cdrag = ", &
           &     (/ ( cdrag(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Lowest level pressure = ", &
           &     (/ ( for_psurf(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Geographical coordinates lon = ", &
           &     (/ (  lon(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Geographical coordinates lat = ", &
           &     (/ (  lat(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Fraction of continent in the grid = ", &
           &     (/ ( for_contfrac(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
      ENDIF
      
      CALL intersurf_main &
 &      (istp, iim, jjm, nbindex, kindex, dt, &
 &       first_CALL, last_CALL, lon, lat, for_contfrac, for_neighbours, for_resolution, date0_rest, &
!     first level conditions
 &       zlev_vec, for_u, for_v, &
 &       for_qair, for_tair, for_eair, for_ccanopy, &
!     Variables for the implicit coupling
 &       cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
!     Rain, snow, radiation and surface pressure
 &       precip_rain, precip_snow, &
 &       for_lwdown, for_swnet, for_swdown, for_psurf, &
!     Output : Fluxes
 &       vevapp, fluxsens, fluxlat, coastalflow, riverflow,  &
!     Surface temperatures and surface properties
 &       tsol_rad, temp_sol_NEW, qsurf, albedo, emis, z0, &
!       VOC : radiation
 &       for_sinang)

!-------
      IF (longprint) THEN
         WRITE(numout,*) "dim2_driver outputs"
         !       Output : Fluxes
         WRITE(numout,*) "vevapp        ; Total of evaporation = ", &
           &     (/ ( vevapp(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Sensible heat flux = ", &
           &     (/ ( fluxsens(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "Latent heat flux = ", &
           &     (/ ( fluxlat(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "coastalflow   ; Diffuse flow of water into the ocean (m^3/dt) = ", &
           &     (/ ( coastalflow(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "riverflow     ; Largest rivers flowing into the ocean (m^3/dt) = ", &
           &     (/ ( riverflow(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         !       Surface temperatures and surface properties
         WRITE(numout,*) "tsol_rad      ; Radiative surface temperature = ", &
           &     (/ ( tsol_rad(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "temp_sol_new  ; New soil temperature = ", &
           &     (/ ( temp_sol_NEW(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "qsurf         ; Surface specific humidity = ", &
           &     (/ ( qsurf(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "albedoVIS = ", &
           &     (/ ( albedo(ilandindex(ik), jlandindex(ik), 1),ik=1,nbindex ) /)
         WRITE(numout,*) "albedoNIR = ", &
           &     (/ ( albedo(ilandindex(ik), jlandindex(ik), 2),ik=1,nbindex ) /)
         WRITE(numout,*) "emis          ; Emissivity = ", &
           &     (/ ( emis(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
         WRITE(numout,*) "z0            ; Surface roughness = ", &
           &     (/ ( z0(ilandindex(ik), jlandindex(ik)),ik=1,nbindex ) /)
      ENDIF
!-----
      dtdt(:,:) = zero
      DO ik=1,nbindex
         i=ilandindex(ik)
         j=jlandindex(ik)
         dtdt(i,j) = ABS(temp_sol_NEW(i,j)-temp_sol_old(i,j))/dt
      ENDDO
!-----
!---- Test if the point with the largest change has more than 5K per dt
!-----
      IF (MAXVAL(dtdt(:,:)) > 5./dt .AND. .FALSE.) THEN ! Chloe 05082014 ajout .AND. .false.
        ml = MAXLOC(dtdt)
        CALL ju2ymds(julian, yy, mm, dd, ss)
        WRITE(numout,"('ATT :',I5,' big temperature jumps on ', &
 &       I4,'-',I2.2,'-',I2.2,':',F8.4)") &
 &       COUNT(dtdt(:,:) > 5./dt),yy,mm,dd,ss/3600.
        WRITE(numout,*) &
 &       'Maximum change of surface temperature located at :', &
 &       lon(ml(1),ml(2)),lat(ml(1),ml(2))
        WRITE(numout,*) 'Coordinates in grid space: ',ml(1),ml(2)
        WRITE(numout,*) 'Change from ',temp_sol_old(ml(1),ml(2)), &
 &       ' to ',temp_sol_new(ml(1),ml(2)),&
 &       'with sw_in = ',for_swnet(ml(1),ml(2))
        old_tair = &
 &        (old_eair(ml(1),ml(2))-cte_grav*old_zlev(ml(1),ml(2)))/cp_air
        WRITE(numout,*) 'Air temperature change from ',old_tair, &
 &       ' to ',for_tair(ml(1),ml(2))
        WRITE(numout,*) 'Max of dtdt : ',dtdt(ml(1),ml(2)),' with dt = ',dt
      ENDIF
!-----
      temp_sol_old(:,:) = temp_sol_NEW(:,:)
!-----
!---- PRINT output value for debug
!-----
      IF (debug) THEN
        WRITE(numout,*) ' OUT : z0, albedoVIS, albedoNIR, emis = ', &
 &       z0(itest,jtest),albedo(itest,jtest,1), &
 &                       albedo(itest,jtest,2),emis(itest,jtest)
        WRITE(numout,*) ' OUT : coastal and river flow = ',&
 &       coastalflow(itest,jtest), riverflow(itest,jtest)
        WRITE(numout,*) ' OUT : tsol_rad, vevapp = ', &
 &       tsol_rad(itest,jtest), vevapp(itest,jtest)
        WRITE(numout,*) ' OUT : temp_sol_new =', temp_sol_NEW(itest,jtest)
      ENDIF
!-----
!---- Give some variables to the output package
!---- for saving on the history tape
!-----
      IF (debug) WRITE(numout,*) 'history written for ', istp
!-----
      istp_old = istp
      istp = istp+1
!-----
      old_zlev(:,:) = zlev_vec(:,:)
      old_qair(:,:) = for_qair(:,:)
      old_eair(:,:) = for_eair(:,:)
!-----
      is = is + 1
    ENDDO
    split_start = 1
!!$    it_force =it_force+1
    IF (it==itau_fin-1) THEN
     CALL stop_timer(timer_global)
     CALL stop_timer(timer_mpi)
    ENDIF
    it = it + 1
!!! CALL histsync
  ENDDO
!-
! Archive in restart file the prognostic variables
!-
  IF (debug) WRITE(numout,*) 'Write the restart for the driver', istp_old
!-
  var_name = 'fluxsens'
  IF (is_root_prc) THEN
     ALLOCATE(fluxsens_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(fluxsens_g(0,1))
  ENDIF
  CALL gather2D(fluxsens , fluxsens_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, fluxsens_g)
  DEALLOCATE(fluxsens_g)
  
  var_name = 'vevapp'
  IF (is_root_prc) THEN
     ALLOCATE(vevapp_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(vevapp_g(0,1))
  ENDIF
  CALL gather2D( vevapp, vevapp_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, vevapp_g)
  DEALLOCATE(vevapp_g)
  
  var_name = 'zlev_old'
  IF (is_root_prc) THEN
     ALLOCATE(old_zlev_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(old_zlev_g(0,1))
  ENDIF
  CALL gather2D( old_zlev, old_zlev_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, old_zlev_g)
  DEALLOCATE(old_zlev_g)
  
  var_name = 'qair_old'
  IF (is_root_prc) THEN
     ALLOCATE(old_qair_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(old_qair_g(0,1))
  ENDIF
  CALL gather2D( old_qair, old_qair_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, old_qair_g)
  DEALLOCATE(old_qair_g)
  
  var_name = 'eair_old'
  IF (is_root_prc) THEN
     ALLOCATE(old_eair_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(old_eair_g(0,1))
  ENDIF
  CALL gather2D( old_eair, old_eair_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, old_eair_g)
  DEALLOCATE(old_eair_g)
  
  var_name = 'rau_old'
  IF (is_root_prc) THEN
     ALLOCATE(for_rau_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(for_rau_g(0,1))
  ENDIF
  CALL gather2D( for_rau, for_rau_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, for_rau_g)
  DEALLOCATE(for_rau_g)
  
  IF (is_root_prc) THEN
     ALLOCATE(albedo_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(albedo_g(0,1))
  ENDIF
  var_name= 'albedo_vis'
  albedo_vis(:,:)=albedo(:,:,1)
  CALL gather2D(albedo_vis,albedo_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, albedo_g)  
  !
  var_name= 'albedo_nir'
  albedo_nir(:,:)=albedo(:,:,2)
  CALL gather2D(albedo_nir,albedo_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, albedo_g)  
  DEALLOCATE(albedo_g)

  IF (is_root_prc) THEN
     ALLOCATE(z0_g(iim_g,jjm_g))
  ELSE
     ALLOCATE(z0_g(0,1))
  ENDIF
  var_name= 'z0'
  CALL gather2D(z0,z0_g)
  IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, z0_g)  
  DEALLOCATE(z0_g)

  if (.NOT. is_watchout) THEN
     var_name = 'petAcoef'
     IF (is_root_prc) THEN
        ALLOCATE(petAcoef_g(iim_g,jjm_g))
     ELSE
        ALLOCATE(petAcoef_g(0,1))
     ENDIF
     CALL gather2D( petAcoef, petAcoef_g)
     IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, petAcoef_g)
     DEALLOCATE(petAcoef_g)
  
     var_name = 'petBcoef'
     IF (is_root_prc) THEN
        ALLOCATE(petBcoef_g(iim_g,jjm_g))
     ELSE
        ALLOCATE(petBcoef_g(0,1))
     ENDIF
     CALL gather2D( petBcoef, petBcoef_g)
     IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, petBcoef_g)
     DEALLOCATE(petBcoef_g)
  
     var_name = 'peqAcoef'
     IF (is_root_prc) THEN
        ALLOCATE(peqAcoef_g(iim_g,jjm_g))
     ELSE
        ALLOCATE(peqAcoef_g(0,1))
     ENDIF
     CALL gather2D( peqAcoef, peqAcoef_g)
     IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, peqAcoef_g)
     DEALLOCATE(peqAcoef_g)
  
     var_name = 'peqBcoef'
     IF (is_root_prc) THEN
        ALLOCATE(peqBcoef_g(iim_g,jjm_g))
     ELSE
        ALLOCATE(peqBcoef_g(0,1))
     ENDIF
     CALL gather2D( peqBcoef, peqBcoef_g)
     IF(is_root_prc) CALL restput (rest_id, var_name, iim_g, jjm_g, 1, istp_old, peqBcoef_g)
     DEALLOCATE(peqBcoef_g)
  ENDIF
!-
  IF (debug) WRITE(numout,*) 'Restart for the driver written'
!=====================================================================
!- 5.0 Closing all files
!=====================================================================
  CALL flinclo(force_id)
  IF ( debug )  WRITE(numout,*) 'FLIN CLOSED'
  CALL histclo
  IF ( debug )   WRITE(numout,*) 'HIST CLOSED'      
    
  IF(is_root_prc) THEN
     CALL restclo
     IF ( debug )  WRITE(numout,*) 'REST CLOSED'
     CALL getin_dump
     IF ( debug )  WRITE(numout,*) 'GETIN CLOSED'
  ENDIF
!---------------
  CALL finalize_para(timer_global,timer_mpi)
!-
  STOP 'END of dim2_driver'
!---------------
END PROGRAM driver
