! =================================================================================================================================
! MODULE       : constantes
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        "constantes" module contains some public technical constants like
!! pi, Earth radius, etc...and non-PFTs externalized parameters.
!!
!!\n DESCRIPTION: In this module, you can set the flag diag_qsat in order to detect the pixel where the
!!                temperature is out of range (look qsatcalc and dev_qsatcalc in qsat_moisture.f90).\n
!!                The Earth radius is approximated by the Equatorial radius.The Earth's equatorial radius a,
!!                or semi-major axis, is the distance from its center to the equator and equals 6,378.1370 km.
!!                The equatorial radius is often used to compare Earth with other planets.\n
!!                The meridional mean is well approximated by the semicubic mean of the two axe yielding 
!!                6367.4491 km or less accurately by the quadratic mean of the two axes about 6,367.454 km
!!                or even just the mean of the two axes about 6,367.445 km.\n
!!
!! RECENT CHANGE(S): Didier Solyga : This module contains now all the externalized parameters of ORCHIDEE 
!!                   listed by modules which are not pft-dependent  
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2013-03-19 18:15:15 +0100 (Tue, 19 Mar 2013) $
!! $Revision: 1223 $
!! \n
!_ ================================================================================================================================

MODULE constantes

  USE defprec
  USE parallel
!-
  IMPLICIT NONE
!-

                         !-----------------------!
                         !  ORCHIDEE CONSTANTS   !
                         !-----------------------!


  !
  ! FLAGS 
  !
  TYPE control_type
    LOGICAL :: river_routing      !! activate river routing (true/false)
    LOGICAL :: hydrol_cwrr        !! activate 11 layers hydrolgy model (true/false)
    LOGICAL :: do_floodplains
    LOGICAL :: do_irrigation
    LOGICAL :: ok_sechiba         !! activate physic of the model (true/false)
    LOGICAL :: ok_co2             !! activate photosynthesis (true/false)
    LOGICAL :: ok_stomate         !! activate carbon cycle (true/false)
    LOGICAL :: ok_dgvm            !! activate dynamic vegetation (true/false)
    LOGICAL :: stomate_watchout   !! activate the creation of restart files for STOMATE even if STOMATE is not activated (true/false)
    LOGICAL :: ok_pheno           !! activate the calculation of lai using stomate rather than a prescription (true/false)
!Isa
!    LOGICAL :: ok_freeze_choisnel ! ce flag est inutile en fait
     LOGICAL :: ok_freeze ! si oui, ce flag active toutes les modifs d'Isabelle en meme temps
     LOGICAL :: ok_converge_isaorig ! si oui, quelques details pour converger avec isaorig
    ! end isa
    LOGICAL :: do_land_use
    LOGICAL :: ok_inca            !! activate biogenic volatile organic coumpounds ? (true/false)
    LOGICAL :: ok_leafage         !! activate leafage? (true/false)
    LOGICAL :: ok_radcanopy       !! use canopy radiative transfer model (true/false)
    LOGICAL :: ok_multilayer      !! use canopy radiative transfer model with multi-layers (true/false)
    LOGICAL :: ok_pulse_NOx       !! calculate NOx emissions with pulse (true/false)
    LOGICAL :: ok_bbgfertil_NOx   !! calculate NOx emissions with bbg fertilizing effect (true/false)
    LOGICAL :: ok_cropsfertil_NOx !! calculate NOx emissions with fertilizers use (true/false)
  END TYPE control_type

  !-
  TYPE(control_type), SAVE :: control  !! Flags that (de)activate parts of the model

  LOGICAL, SAVE :: OFF_LINE_MODE = .FALSE.  !! ORCHIDEE detects if it is coupled with a GCM or 
                                            !! just use with one driver in OFF-LINE. (true/false)

  !
  ! TIME
  !
  REAL(r_std), SAVE :: one_day  !! One day in seconds (s)
  REAL(r_std), SAVE :: one_year !! One year in seconds (s)
  REAL(r_std), PARAMETER :: one_hour = 3600.0  !! One hour in seconds (s)


  !
  ! SPECIAL VALUES 
  !
  INTEGER(i_std), PARAMETER :: undef_int = 999999999     !! undef integer for integer arrays (unitless)
  !-
  REAL(r_std), SAVE :: val_exp = 999999.                 !! Specific value if no restart value  (unitless)
  REAL(r_std), PARAMETER :: undef = -9999.               !! Special value for stomate (unitless)
  !-
  REAL(r_std), PARAMETER :: min_sechiba = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: undef_sechiba = 1.E+20_r_std !! The undef value used in SECHIBA (unitless)
  !-
  REAL(r_std), PARAMETER :: min_stomate = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: large_value = 1.E33_r_std    !! some large value (for stomate) (unitless)


  !
  !  DIMENSIONING AND INDICES PARAMETERS  
  !
  INTEGER(i_std), PARAMETER :: ivis = 1          !! index for albedo in visible range (unitless)
  INTEGER(i_std), PARAMETER :: inir = 2          !! index for albeod i near-infrared range (unitless) 
  INTEGER(i_std), PARAMETER :: nnobio = 1        !! Number of other surface types: land ice (lakes,cities, ...) (unitless)
  INTEGER(i_std), PARAMETER :: iice = 1          !! Index for land ice (see nnobio) (unitless)
  !-
  !! Soil
!!$  INTEGER(i_std), PARAMETER :: ngrnd = 7         !! Number of soil level (unitless)
!!$  INTEGER(i_std), PARAMETER :: nbdl = 11         !! Number of diagnostic levels in the soil (unitless)
!!$!MM : if you want to compare hydrology variables with old TAG 1.6 and lower, 
!!$!     you must set the Number of diagnostic levels in the soil to 6 :
!!$!    INTEGER(i_std),PARAMETER :: nbdl=6
!!$  INTEGER(i_std), PARAMETER :: nslm = 11         !! Number of levels in CWRR (unitless)
!!$  INTEGER(i_std), PARAMETER :: nstm = 3          !! Number of soil types (unitless)
  INTEGER(i_std), PARAMETER :: classnb = 9       !! Levels of soil colour classification (unitless)
  !-
  INTEGER(i_std), PARAMETER :: nleafages = 4     !! leaf age discretisation ( 1 = no discretisation )(unitless)
  !-
  !! litter fractions: indices (unitless)
  INTEGER(i_std), PARAMETER :: ileaf = 1         !! Index for leaf compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapabove = 2     !! Index for sapwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapbelow = 3     !! Index for sapwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartabove = 4   !! Index for heartwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartbelow = 5   !! Index for heartwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iroot = 6         !! Index for roots compartment (unitless)
  INTEGER(i_std), PARAMETER :: ifruit = 7        !! Index for fruits compartment (unitless)
  INTEGER(i_std), PARAMETER :: icarbres = 8      !! Index for reserve compartment (unitless)
  INTEGER(i_std), PARAMETER :: nparts = 8        !! Number of biomass compartments (unitless)
  !-
  !! indices for assimilation parameters 
  INTEGER(i_std), PARAMETER :: itmin = 1         !! Index for minimum photosynthesis temperature (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: itopt = 2         !! Index for optimal photosynthesis temperature (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: itmax = 3         !! Index for maxmimum photosynthesis temperature (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: ivcmax = 4        !! Index for vcmax (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: ivjmax = 5        !! Index for vjmax (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: npco2 = 5         !! Number of assimilation parameters (unitless)
  !-
  !! trees and litter: indices for the parts of heart-
  !! and sapwood above and below the ground 
  INTEGER(i_std), PARAMETER :: iabove = 1       !! Index for above part (unitless)
  INTEGER(i_std), PARAMETER :: ibelow = 2       !! Index for below part (unitless)
  INTEGER(i_std), PARAMETER :: nlevs = 2        !! Number of levels for trees and litter (unitless)
  !-
  !! litter: indices for metabolic and structural part
  INTEGER(i_std), PARAMETER :: imetabolic = 1   !! Index for metabolic litter (unitless)
  INTEGER(i_std), PARAMETER :: istructural = 2  !! Index for structural litter (unitless)
  INTEGER(i_std), PARAMETER :: nlitt = 2        !! Number of levels for litter compartments (unitless)
  !-
  !! carbon pools: indices
  INTEGER(i_std), PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ncarb = 3        !! Number of soil carbon pools (unitless)


  !
  ! NUMERICAL AND PHYSICS CONSTANTS
  !
  !

  !-
  ! 1. Mathematical and numerical constants
  !-
  REAL(r_std), PARAMETER :: pi = 3.141592653589793238   !! pi souce : http://mathworld.wolfram.com/Pi.html (unitless)
  REAL(r_std), PARAMETER :: euler = 2.71828182845904523 !! e source : http://mathworld.wolfram.com/e.html (unitless)
  REAL(r_std), PARAMETER :: zero = 0._r_std             !! Numerical constant set to 0 (unitless)
  REAL(r_std), PARAMETER :: undemi = 0.5_r_std          !! Numerical constant set to 1/2 (unitless)
  REAL(r_std), PARAMETER :: un = 1._r_std               !! Numerical constant set to 1 (unitless)
  REAL(r_std), PARAMETER :: moins_un = -1._r_std        !! Numerical constant set to -1 (unitless)
  REAL(r_std), PARAMETER :: deux = 2._r_std             !! Numerical constant set to 2 (unitless)
  REAL(r_std), PARAMETER :: trois = 3._r_std            !! Numerical constant set to 3 (unitless)
  REAL(r_std), PARAMETER :: quatre = 4._r_std           !! Numerical constant set to 4 (unitless)
  REAL(r_std), PARAMETER :: cinq = 5._r_std             !![DISPENSABLE] Numerical constant set to 5 (unitless)
  REAL(r_std), PARAMETER :: six = 6._r_std              !![DISPENSABLE] Numerical constant set to 6 (unitless)
  REAL(r_std), PARAMETER :: huit = 8._r_std             !! Numerical constant set to 8 (unitless)
  REAL(r_std), PARAMETER :: mille = 1000._r_std         !! Numerical constant set to 1000 (unitless)

  !-
  ! 2 . Physics
  !-
  REAL(r_std), PARAMETER :: R_Earth = 6378000.              !! radius of the Earth : Earth radius ~= Equatorial radius (m)
  REAL(r_std), PARAMETER :: mincos  = 0.0001                !! Minimum cosine value used for interpolation (unitless) 
  REAL(r_std), PARAMETER :: pb_std = 1013.                  !! standard pressure (hPa)
  REAL(r_std), PARAMETER :: ZeroCelsius = 273.15            !! Freezing point (K)
  REAL(r_std), PARAMETER :: tp_00 = 273.15                  !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: chalsu0 = 2.8345E06             !! Latent heat of sublimation (J.kg^{-1})
  REAL(r_std), PARAMETER :: chalev0 = 2.5008E06             !! Latent heat of evaporation (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: chalfu0 = chalsu0-chalev0       !! Latent heat of fusion (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: c_stefan = 5.6697E-8            !! Stefan-Boltzman constant (W.m^{-2}.K^{-4})
  REAL(r_std), PARAMETER :: cp_air = 1004.675               !! Specific heat of dry air (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: cte_molr = 287.05               !! Specific constant of dry air (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: kappa = cte_molr/cp_air         !! Kappa : ratio between specific constant and specific heat 
                                                            !! of dry air (unitless)
  REAL(r_std), PARAMETER :: msmlr_air = 28.964E-03          !! Molecular weight of dry air (kg.mol^{-1})
  REAL(r_std), PARAMETER :: msmlr_h2o = 18.02E-03           !! Molecular weight of water vapor (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: cp_h2o = &                      !! Specific heat of water vapor (J.kg^{-1}.K^{-1}) 
       & cp_air*(quatre*msmlr_air)/( 3.5_r_std*msmlr_h2o) 
  REAL(r_std), PARAMETER :: cte_molr_h2o = cte_molr/quatre  !! Specific constant of water vapor (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: retv = msmlr_air/msmlr_h2o-un   !! Ratio between molecular weight of dry air and water 
                                                            !! vapor minus 1(unitless)  
  REAL(r_std), PARAMETER :: rvtmp2 = cp_h2o/cp_air-un       !! Ratio between specific heat of water vapor and dry air
                                                            !! minus 1 (unitless)
  REAL(r_std), PARAMETER :: cepdu2 = (0.1_r_std)**2         !! Squared wind shear (m^2.s^{-2}) 
  REAL(r_std), PARAMETER :: ct_karman = 0.35_r_std          !! Van Karmann Constant (unitless)
  REAL(r_std), PARAMETER :: cte_grav = 9.80665_r_std        !! Acceleration of the gravity (m.s^{-2})
  REAL(r_std), PARAMETER :: pa_par_hpa = 100._r_std         !! Transform pascal into hectopascal (unitless)
  REAL(r_std), PARAMETER :: R = 8.314                       !! Ideal gas constant (J.mol^{-1}.K^{-1})
  REAL(r_std), PARAMETER :: Sct = 1370.                     !! Solar constant (W.m^{-2}) 


  !-
  ! 3. Climatic constants
  !-
  !! Constantes of the Louis scheme 
  REAL(r_std), PARAMETER :: cb = 5._r_std         !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
  REAL(r_std), PARAMETER :: cc = 5._r_std         !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
  REAL(r_std), PARAMETER :: cd = 5._r_std         !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
  !-
  REAL(r_std), PARAMETER :: rayt_cste = 125.      !! Constant in the computation of surface resistance (W.m^{-2})
  REAL(r_std), PARAMETER :: defc_plus = 23.E-3    !! Constant in the computation of surface resistance (K.W^{-1})
  REAL(r_std), PARAMETER :: defc_mult = 1.5       !! Constant in the computation of surface resistance (K.W^{-1})

  !-
  ! 4. Soil thermodynamics constants
  !-
!!$  REAL(r_std), PARAMETER :: so_cond = 1.5396       !! Average Thermal Conductivity of soils (W.m^{-2}.K^{-1})
!!$  REAL(r_std), PARAMETER :: so_capa = 2.0514e+6    !! Average Heat capacity of soils (J.m^{-3}.K^{-1}) 
!!$  !-
!!$  !! Values taken from : PIELKE,'MESOSCALE METEOROLOGICAL MODELING',P.384
!!$  !  Dry soil heat capacity was decreased and conductivity increased.
!!$  REAL(r_std), PARAMETER :: so_capa_dry = 1.80e+6  !! Dry soil Heat capacity of soils (J.m^{-3}.K^{-1}) 
!!$  REAL(r_std), PARAMETER :: so_cond_dry = 0.40     !! Dry soil Thermal Conductivity of soils (W.m^{-2}.K^{-1})
!!$  REAL(r_std), PARAMETER :: so_capa_wet = 3.03e+6  !! Wet soil Heat capacity of soils (J.m^{-3}.K^{-1})
!!$  REAL(r_std), PARAMETER :: so_cond_wet = 1.89     !! Wet soil Thermal Conductivity of soils (W.m^{-2}.K^{-1})
!!$  REAL(r_std), PARAMETER :: sn_cond = 0.3          !! Thermal Conductivity of snow (W.m^{-2}.K^{-1})
!!$  REAL(r_std), PARAMETER :: sn_dens = 330.0        !! Snow density for the soil thermodynamics (unitless)
!!$  REAL(r_std), PARAMETER :: sn_capa = 2100.0_r_std*sn_dens !! Heat capacity for snow (J.m^{-3}.K^{-1}) 

  !
  ! OPTIONAL PARTS OF THE MODEL
  !
  LOGICAL, SAVE     :: long_print = .FALSE.       !! To set for more printing
  LOGICAL,PARAMETER :: diag_qsat = .TRUE.         !! One of the most frequent problems is a temperature out of range
                                                  !! we provide here a way to catch that in the calling procedure. 
                                                  !! (from Jan Polcher)(true/false) 
  LOGICAL            :: almaoutput                !! Selects the type of output for the model.(true/false)
                                                  !! Value is read from run.def in intersurf_history

!!$  !
!!$  ! DIAGNOSTIC VARIABLES
!!$  !
!!$  REAL(r_std),DIMENSION(nbdl),SAVE :: diaglev     !! The lower limit of the layer on which soil moisture (relative)
!!$                                                  !! and temperature are going to be diagnosed.
!!$                                                  !! These variables are made for transfering the information
!!$                                                  !! to the biogeophyical processes modelled in STOMATE. (unitless)
  !
  ! DIVERSE
  !
  CHARACTER(LEN=100) :: stomate_forcing_name='NONE'  !! NV080800 Name of STOMATE forcing file (unitless)
                                                     ! Compatibility with Nicolas Viovy driver.
  CHARACTER(LEN=100) :: stomate_Cforcing_name='NONE' !! NV080800 Name of soil forcing file (unitless)
                                                     ! Compatibility with Nicolas Viovy driver.
  INTEGER(i_std), SAVE :: forcing_id                 !! Index of the forcing file (unitless)




                         !------------------------!
                         !  SECHIBA PARAMETERS    !
                         !------------------------!
 

  !
  ! GLOBAL PARAMETERS   
  !
  REAL(r_std), SAVE :: min_wind = 0.1      !! The minimum wind (m.s^{-1})
  REAL(r_std), SAVE :: snowcri = 1.5       !! Sets the amount above which only sublimation occures (kg.m^{-2})
!!$  REAL(r_std), SAVE :: qsintcst = 0.1      !! Transforms leaf area index into size of interception reservoir (unitless)

  !
  ! FLAGS ACTIVATING SUB-MODELS
  !
  LOGICAL, SAVE :: treat_expansion = .FALSE.   !! Do we treat PFT expansion across a grid point after introduction? (true/false)
  LOGICAL, SAVE :: ok_herbivores = .FALSE.     !! flag to activate herbivores (true/false)
  LOGICAL, SAVE :: harvest_agri = .TRUE.       !! flag to harvest aboveground biomass from agricultural PFTs)(true/false)
  LOGICAL, SAVE :: lpj_gap_const_mort = .TRUE. !! constant moratlity (true/false)
  LOGICAL, SAVE :: disable_fire = .FALSE.      !! flag that disable fire (true/false)

  !
  ! CONFIGURATION VEGETATION
  !
  LOGICAL, SAVE :: agriculture = .TRUE.    !! allow agricultural PFTs (true/false)
  LOGICAL, SAVE :: impveg = .FALSE.        !! Impose vegetation ? (true/false)
  LOGICAL, SAVE :: impsoilt = .FALSE.      !! Impose soil ? (true/false)
  LOGICAL, SAVE :: lcchange = .FALSE.      !! Land cover change flag (true/false)
  LOGICAL, SAVE :: read_lai = .FALSE.      !! Flag to read a map of LAI if STOMATE is not activated (true/false)
  LOGICAL, SAVE :: old_lai = .FALSE.       !! Flag for the old LAI map interpolation (SHOULD BE DROPED ??)(true/false)
  LOGICAL, SAVE :: old_veget = .FALSE.     !! Flag to use the old vegetation Map interpolation (SHOULD BE DROPED ?)(true/false)
  LOGICAL, SAVE :: land_use = .TRUE.       !! flag to account or not for Land Use  (true/false)
  LOGICAL, SAVE :: veget_reinit = .TRUE.   !! To change LAND USE file in a run. (true/false)

  !
  ! PARAMETERS USED BY BOTH HYDROLOGY MODELS
  !
  REAL(r_std), SAVE :: max_snow_age = 50._r_std !! Maximum period of snow aging (days)
  REAL(r_std), SAVE :: snow_trans = 0.3_r_std   !! Transformation time constant for snow (m)
  REAL(r_std), SAVE :: sneige                   !! Lower limit of snow amount (kg.m^{-2})
  REAL(r_std), SAVE :: maxmass_glacier = 3000.  !! The maximum mass of a glacier (kg.m^{-2})
!!$  REAL(r_std), SAVE :: mx_eau_eau = 150.        !! Maximum quantity of water (kg.m^{-3})
 
  !
  ! BVOC : Biogenic activity  for each age class
  !
  REAL(r_std), SAVE, DIMENSION(nleafages) :: iso_activity = (/0.5, 1.5, 1.5, 0.5/)     !! Biogenic activity for each 
                                                                                       !! age class : isoprene (unitless)
  REAL(r_std), SAVE, DIMENSION(nleafages) :: methanol_activity = (/1., 1., 0.5, 0.5/)  !! Biogenic activity for each
                                                                                       !! age class : methanol (unnitless)

  !
  ! condveg.f90
  !

  ! 1. Scalar

  ! 1.1 Flags used inside the module

  LOGICAL, SAVE :: alb_bare_model = .FALSE. !! Switch for choosing values of bare soil 
                                            !! albedo (see header of subroutine)
                                            !! (true/false)
  LOGICAL, SAVE :: impaze = .FALSE.         !! Switch for choosing surface parameters
                                            !! (see header of subroutine).  
                                            !! (true/false)
  LOGICAL, SAVE :: z0cdrag_ave = .TRUE.     !! Chooses between two methods to calculate the 
                                            !! grid average of the roughness (see header of subroutine)   
                                            !! (true/false)

  ! 1.2 Others 

  REAL(r_std), SAVE :: z0_over_height = un/16.           !! Factor to calculate roughness height from 
                                                         !! vegetation height (unitless)   
  REAL(r_std), SAVE :: height_displacement = 0.75        !! Factor to calculate the zero-plane displacement
                                                         !! height from vegetation height (m)
  REAL(r_std), SAVE :: z0_bare = 0.01                    !! bare soil roughness length (m)
  REAL(r_std), SAVE :: z0_ice = 0.001                    !! ice roughness length (m)
  REAL(r_std), SAVE :: tcst_snowa = 5.0                  !! Time constant of the albedo decay of snow (days)
  REAL(r_std), SAVE :: snowcri_alb = 10.                 !! Critical value for computation of snow albedo (kg.m^{-2})
  REAL(r_std), SAVE :: fixed_snow_albedo = undef_sechiba !! To choose a fixed snow albedo value (unitless)
  REAL(r_std), SAVE :: z0_scal = 0.15                    !! Surface roughness height imposed (m)
  REAL(r_std), SAVE :: roughheight_scal = zero           !! Effective roughness Height depending on zero-plane 
                                                         !! displacement height (m) (imposed)

  REAL(r_std), SAVE :: emis_scal = 1.0                   !! Surface emissivity imposed (unitless)

  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(2) :: alb_deadleaf = (/ .12, .35/)    !! albedo of dead leaves, VIS+NIR (unitless)
  REAL(r_std), SAVE, DIMENSION(2) :: alb_ice = (/ .60, .20/)         !! albedo of ice, VIS+NIR (unitless)
  REAL(r_std), SAVE, DIMENSION(2) :: albedo_scal = (/ 0.25, 0.25 /)  !! Albedo values for visible and near-infrared 
                                                                     !! used imposed (unitless) 
  REAL(r_std), DIMENSION(classnb) :: vis_dry = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in visible range
  REAL(r_std), DIMENSION(classnb) :: nir_dry = (/0.48,&
       &0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in near-infrared range 
  REAL(r_std), DIMENSION(classnb) :: vis_wet = (/0.12,&
       &0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in visible range 
  REAL(r_std), DIMENSION(classnb) :: nir_wet = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in near-infrared range
  REAL(r_std), DIMENSION(classnb) :: albsoil_vis = (/ &
       &0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25/)   !! Soil albedo values to soil colour classification:
                                                                   !! Averaged of wet and dry soil albedo values
                                                                   !! in visible and near-infrared range
  REAL(r_std), DIMENSION(classnb) :: albsoil_nir = (/ &
       &0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45/)  !! Soil albedo values to soil colour classification:
                                                                !! Averaged of wet and dry soil albedo values
                                                                !! in visible and near-infrared range


  !
  ! diffuco.f90
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Tetens_1 = 0.622         !! Ratio between molecular weight of water vapor and molecular weight  
                                                     !! of dry air (unitless)
  REAL(r_std), PARAMETER :: Tetens_2 = 0.378         !!
  REAL(r_std), PARAMETER :: std_ci_frac = 0.667      !!
  REAL(r_std), PARAMETER :: alpha_j = 0.8855         !! Quantum yield of RuBP regeneration 
  REAL(r_std), PARAMETER :: curve_assim = 0.7        !! Curvature of the quantum response (unitless)
  REAL(r_std), PARAMETER :: WJ_coeff1 = 4.5          !! First coefficient for calculating the generation-limited rate RuBP (unitless)
  REAL(r_std), PARAMETER :: WJ_coeff2 = 10.5         !! Second coefficient for calculating the generation-limited rate RuBP (unitless)
  REAL(r_std), PARAMETER :: Vc_to_Rd_ratio = 0.011   !!
  REAL(r_std), PARAMETER :: O2toCO2_stoechio = 1.6   !! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)
  REAL(r_std), PARAMETER :: mmol_to_m_1 = 0.0244     !!
  REAL(r_std), PARAMETER :: RG_to_PAR = 0.5          !!
  REAL(r_std), PARAMETER :: W_to_mmol = 4.6          !! W_to_mmol * RG_to_PAR = 2.3

  ! 1. Scalar

  INTEGER(i_std), SAVE :: nlai = 20             !! Number of LAI levels (unitless)
  LOGICAL, SAVE :: ldq_cdrag_from_gcm = .FALSE. !! Set to .TRUE. if you want q_cdrag coming from GCM
  REAL(r_std), SAVE :: laimax = 12.             !! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
  REAL(r_std), SAVE :: xc4_1 = 0.83             !! Factor in the first Collatz equation for C4 plants (unitless)
  REAL(r_std), SAVE :: xc4_2 = 0.93             !! Factor in the second Collatz equation for C4 plants (unitless)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: lai_level_depth = 0.15  !!
  REAL(r_std), SAVE :: x1_coef =  0.177        !! Multiplicative factor for calculating the pseudo first order rate constant 
                                               !! of assimilation response to co2 kt (unitless)
  REAL(r_std), SAVE :: x1_Q10 =  0.069         !! Exponential factor in the equation defining kt (unitless)
  REAL(r_std), SAVE :: quantum_yield =  0.092  !!
  REAL(r_std), SAVE :: kt_coef = 0.7           !! Multiplicative factor in the equation defining kt (unitless)
  REAL(r_std), SAVE :: kc_coef = 39.09         !! Multiplicative factor for calculating the Michaelis-Menten 
                                               !! coefficient Kc (unitless)
  REAL(r_std), SAVE :: Ko_Q10 = 0.085          !! Exponential factor for calculating the Michaelis-Menten coefficients 
                                               !! Kc and Ko (unitless)
  REAL(r_std), SAVE :: Oa = 210000.            !! Intercellular concentration of O2 (ppm)
  REAL(r_std), SAVE :: Ko_coef =  2.412        !! Multiplicative factor for calculating the Michaelis-Menten coefficient Ko (unitless)
  REAL(r_std), SAVE :: CP_0 = 42.              !! Multiplicative factor for calculating the CO2 compensation point CP (unitless)
  REAL(r_std), SAVE :: CP_temp_coef = 9.46     !! Exponential factor for calculating the CO2 compensation point CP (unitless)
  REAL(r_std), SAVE :: CP_temp_ref = 25.       !! Reference temperature for the CO2 compensation point CP (C)
  !
  REAL(r_std), SAVE, DIMENSION(2) :: rt_coef = (/ 0.8, 1.3 /)    !! 
  REAL(r_std), SAVE, DIMENSION(2) :: vc_coef = (/ 0.39, 0.3 /)   !!
  !
  REAL(r_std), SAVE, DIMENSION(6) :: dew_veg_poly_coeff = &            !! coefficients of the 5 degree polynomomial used
  & (/ 0.887773, 0.205673, 0.110112, 0.014843,  0.000824,  0.000017 /) !! in the equation of coeff_dew_veg



  !
  ! hydrolc.f90
  !

  ! 1. Scalar

!!$  LOGICAL, SAVE     :: ok_hdiff  = .FALSE.        !! do horizontal diffusion? (true/false)
!!$  REAL(r_std), SAVE :: qwilt = 5.0                !! Wilting point (Has a numerical role for the moment) (unitless)
!!$  REAL(r_std), SAVE :: min_resdis = 2.e-5         !! The minimal size we allow for the upper reservoir (m)
!!$  REAL(r_std), SAVE :: min_drain = 0.001          !! Diffusion constant for the slow regime (kg.m^{-2}.dt^{-1}) 
!!$                                                  !! (This is for the diffusion between reservoirs)
!!$  REAL(r_std), SAVE :: max_drain = 0.1            !! Diffusion constant for the fast regime (kg.m^{-2}.dt^{-1}) 
!!$  REAL(r_std), SAVE :: exp_drain = 1.5            !! The exponential in the diffusion law (unitless)
!!$  REAL(r_std), SAVE :: rsol_cste = 33.E3          !! Constant in the computation of resistance for bare soil evaporation (s.m^{-2})
!!$  REAL(r_std), SAVE :: hcrit_litter = 0.08_r_std  !! Scaling depth for litter humidity (m)



  !
  ! hydrol.f90
  !

  ! 0. Constants
!Chloe : défini dans constantes_soil
  !!$INTEGER(i_std),PARAMETER :: imin = 1       !! CWRR linearisation (unitless)
  !!$INTEGER(i_std),PARAMETER :: nbint = 100    !! number of interval for CWRR (unitless)
  !!$INTEGER(i_std),PARAMETER :: imax = nbint+1 !! number of points for CWRR (unitless)

  ! 1. Scalar

!!$  REAL(r_std), SAVE :: w_time = un  !! Time weighting for discretisation (unitless)

  !
  ! routing.f90
  !

  ! 1. Scalar
!!$
!!$  REAL(r_std), SAVE :: crop_coef = 1.5   !! Empirical crop coefficient dependent on vegetation characteristics
!!$                                         !! according to Kassel irrigation parametrization.
!!$                                         !! When potential transpiration is used this coefficient has another interpretation (unitless)

  !
  ! slowproc.f90 
  !

  ! 1. Scalar

  INTEGER(i_std), SAVE :: veget_year_orig = 0        !!  first year for landuse (number)
  REAL(r_std), SAVE :: clayfraction_default = 0.2    !! Default value for clay fraction (0-1, unitless)
  REAL(r_std), SAVE :: min_vegfrac = 0.001           !! Minimal fraction of mesh a vegetation type can occupy (0-1, unitless)
  REAL(r_std), SAVE :: frac_nobio_fixed_test_1 = 0.0 !! Value for frac_nobio for tests in 0-dim simulations (0-1, unitless)
                                                     
  REAL(r_std), SAVE :: stempdiag_bid = 280.          !! only needed for an initial LAI if there is no restart file



                           !-----------------------------!
                           !  STOMATE AND LPJ PARAMETERS !
                           !-----------------------------!


  !
  ! lpj_constraints.f90
  !
  
  ! 1. Scalar

  REAL(r_std), SAVE  :: too_long = 5.      !! longest sustainable time without 
                                           !! regeneration (vernalization) (years)



  !
  ! lpj_establish.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: estab_max_tree = 0.12   !! Maximum tree establishment rate (0-1, unitless)
  REAL(r_std), SAVE :: estab_max_grass = 0.12  !! Maximum grass establishment rate (0-1, unitless)
  
  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: establish_scal_fact = 5.  !!
  REAL(r_std), SAVE :: max_tree_coverage = 0.98  !! (0-1, unitless)
  REAL(r_std), SAVE :: ind_0_estab = 0.2         !! = ind_0 * 10.



  !
  ! lpj_fire.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tau_fire = 30.           !! Time scale for memory of the fire index (days).
  REAL(r_std), SAVE :: litter_crit = 200.       !! Critical litter quantity for fire
                                                !! below which iginitions extinguish 
                                                !! @tex $(gC m^{-2})$ @endtex
  REAL(r_std), SAVE :: fire_resist_struct = 0.5 !!

  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(nparts) :: co2frac = &    !! The fraction of the different biomass 
       & (/ .95, .95, 0., 0.3, 0., 0., .95, .95 /)       !! compartments emitted to the atmosphere 
                                                         !! when burned (unitless, 0-1)  

  ! 3. Coefficients of equations

  REAL(r_std), SAVE, DIMENSION(3) :: bcfrac_coeff = (/ .3,  1.3,  88.2 /)         !! (unitless)
  REAL(r_std), SAVE, DIMENSION(4) :: firefrac_coeff = (/ 0.45, 0.8, 0.6, 0.13 /)  !! (unitless)


  !
  ! lpj_gap.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: ref_greff = 0.035         !! Asymptotic maximum mortality rate
                                                 !! @tex $(year^{-1})$ @endtex

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: availability_fact = 0.1   !!


  !               
  ! lpj_light.f90 
  !              

  ! 1. Scalar
  
  LOGICAL, SAVE :: annual_increase = .TRUE. !! for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or
                                            !! to fpc of last time step (F)? (true/false)
  REAL(r_std), SAVE :: min_cover = 0.05     !! For trees, minimum fraction of crown area occupied
                                            !! (due to its branches etc.) (0-1, unitless)
                                            !! This means that only a small fraction of its crown area
                                            !! can be invaded by other trees.
  !
  ! lpj_pftinout.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: min_avail = 0.01         !! minimum availability
  REAL(r_std), SAVE :: ind_0 = 0.02             !! initial density of individuals

  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: RIP_time_min = 1.25      !! test whether the PFT has been eliminated lately (years)
  REAL(r_std), SAVE :: npp_longterm_init = 10.  !! Initialisation value for npp_longterm (gC.m^{-2}.year^{-1})
  REAL(r_std), SAVE :: everywhere_init = 0.05   !!



  !
  ! stomate_alloc.f90
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: max_possible_lai = 10. !! (m^2.m^{-2})
  REAL(r_std), PARAMETER :: Nlim_Q10 = 10.         !!

  ! 1. Scalar

  LOGICAL, SAVE :: ok_minres = .TRUE.              !! [DISPENSABLE] Do we try to reach a minimum reservoir even if
                                                   !! we are severely stressed? (true/false)

  REAL(r_std), SAVE :: tau_leafinit = 10.          !! Time required to develop a minimal LAI
                                                   !! using the carbohydrate reserve (days)
  REAL(r_std), SAVE :: reserve_time_tree = 30.     !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! trees (days)
  REAL(r_std), SAVE :: reserve_time_grass = 20.    !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! grasses (days)


  REAL(r_std), SAVE :: R0 = 0.3                    !! Default root allocation (0-1, unitless)
  REAL(r_std), SAVE :: S0 = 0.3                    !! Default sapwood allocation (0-1, unitless)
  REAL(r_std), SAVE :: L0                          !! Default leaf allocation (0-1, unitless)
  REAL(r_std), SAVE :: f_fruit = 0.1               !! Default fruit allocation (0-1, unitless)
  REAL(r_std), SAVE :: alloc_sap_above_grass = 1.0 !! fraction of sapwood allocation above ground
                                                   !! for grass (0-1, unitless)
  REAL(r_std), SAVE :: min_LtoLSR = 0.2            !! Prescribed lower bounds for leaf 
                                                   !! allocation (0-1, unitless)
  REAL(r_std), SAVE :: max_LtoLSR = 0.5            !! Prescribed upper bounds for leaf 
                                                   !! allocation (0-1, unitless)
  REAL(r_std), SAVE :: z_nitrogen = 0.2            !! Curvature of the root profile (m)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: lai_max_to_happy = 0.5      !!
  REAL(r_std), SAVE :: Nlim_tref = 25.             !! (C)



  !
  ! stomate_data.f90 
  !

  ! 1. Scalar 

  ! 1.1 Parameters for the pipe model

  REAL(r_std), SAVE :: pipe_tune1 = 100.0        !! crown area = pipe_tune1. stem diameter**(1.6) (Reinicke's theory) (unitless)
  REAL(r_std), SAVE :: pipe_tune2 = 40.0         !! height=pipe_tune2 * diameter**pipe_tune3 (unitless)
  REAL(r_std), SAVE :: pipe_tune3 = 0.5          !! height=pipe_tune2 * diameter**pipe_tune3 (unitless)
  REAL(r_std), SAVE :: pipe_tune4 = 0.3          !! needed for stem diameter (unitless)
  REAL(r_std), SAVE :: pipe_density = 2.e5       !! Density
  REAL(r_std), SAVE :: pipe_k1 = 8.e3            !! one more SAVE
  REAL(r_std), SAVE :: pipe_tune_exp_coeff = 1.6 !! pipe tune exponential coeff (unitless)

  ! 1.2 climatic parameters 

  REAL(r_std), SAVE :: precip_crit = 100.        !! minimum precip, in (mm/year)
  REAL(r_std), SAVE :: gdd_crit_estab = 150.     !! minimum gdd for establishment of saplings
  REAL(r_std), SAVE :: fpc_crit = 0.95           !! critical fpc, needed for light competition and establishment (0-1, unitless)

  ! 1.3 sapling characteristics

  REAL(r_std), SAVE :: alpha_grass = 0.5         !! alpha coefficient for grasses (unitless)
  REAL(r_std), SAVE :: alpha_tree = 1.           !! alpha coefficient for trees (unitless)
  REAL(r_std), SAVE :: mass_ratio_heart_sap = 3. !! mass ratio (heartwood+sapwood)/sapwood (unitless)
  REAL(r_std), SAVE :: frac_growthresp = 0.28    !! fraction of GPP which is lost as growth respiration (0-1, unitless)

  ! 1.4  time scales for phenology and other processes (in days)

  REAL(r_std), SAVE :: tau_hum_month = 20.        !! (days)       
  REAL(r_std), SAVE :: tau_hum_week = 7.          !! (days)  
  REAL(r_std), SAVE :: tau_t2m_month = 20.        !! (days)      
  REAL(r_std), SAVE :: tau_t2m_week = 7.          !! (days)  
  REAL(r_std), SAVE :: tau_tsoil_month = 20.      !! (days)     
  REAL(r_std), SAVE :: tau_soilhum_month = 20.    !! (days)     
  REAL(r_std), SAVE :: tau_gpp_week = 7.          !! (days)  
  REAL(r_std), SAVE :: tau_gdd = 40.              !! (days)  
  REAL(r_std), SAVE :: tau_ngd = 50.              !! (days)  
  REAL(r_std), SAVE :: coeff_tau_longterm = 3.    !! (unitless)
  REAL(r_std), SAVE :: tau_longterm               !! (days)  

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_carbres = 5.             !!
  REAL(r_std), SAVE :: bm_sapl_sapabove = 0.5           !!
  REAL(r_std), SAVE :: bm_sapl_heartabove = 2.          !!
  REAL(r_std), SAVE :: bm_sapl_heartbelow = 2.          !!
  REAL(r_std), SAVE :: init_sapl_mass_leaf_nat = 0.1    !!
  REAL(r_std), SAVE :: init_sapl_mass_leaf_agri = 1.    !!
  REAL(r_std), SAVE :: init_sapl_mass_carbres = 5.      !!
  REAL(r_std), SAVE :: init_sapl_mass_root = 0.1        !!
  REAL(r_std), SAVE :: init_sapl_mass_fruit = 0.3       !!  
  REAL(r_std), SAVE :: cn_sapl_init = 0.5               !!
  REAL(r_std), SAVE :: migrate_tree = 10.*1.E3          !!
  REAL(r_std), SAVE :: migrate_grass = 10.*1.E3         !!
  REAL(r_std), SAVE :: lai_initmin_tree = 0.3           !!
  REAL(r_std), SAVE :: lai_initmin_grass = 0.1          !!
  REAL(r_std), SAVE, DIMENSION(2) :: dia_coeff = (/ 4., 0.5 /)            !!
  REAL(r_std), SAVE, DIMENSION(2) :: maxdia_coeff =(/ 100., 0.01/)        !!
  REAL(r_std), SAVE, DIMENSION(4) :: bm_sapl_leaf = (/ 4., 4., 0.8, 5./)  !!



  !
  ! stomate_litter.f90 
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Q10 = 10.               !!

  ! 1. Scalar

  REAL(r_std), SAVE :: z_decomp = 0.2               !!  Maximum depth for soil decomposer's activity (m)

  ! 2. Arrays

  REAL(r_std), SAVE :: frac_soil_struct_aa = 0.55   !! corresponding to frac_soil(istructural,iactive,iabove) 
  REAL(r_std), SAVE :: frac_soil_struct_ab = 0.45   !! corresponding to frac_soil(istructural,iactive,ibelow)
  REAL(r_std), SAVE :: frac_soil_struct_sa = 0.7    !! corresponding to frac_soil(istructural,islow,iabove)
  REAL(r_std), SAVE :: frac_soil_struct_sb = 0.7    !! corresponding to frac_soil(istructural,islow,ibelow)
  REAL(r_std), SAVE :: frac_soil_metab_aa = 0.45    !! corresponding to frac_soil(imetabolic,iactive,iabove)
  REAL(r_std), SAVE :: frac_soil_metab_ab = 0.45    !! corresponding to frac_soil(imetabolic,iactive,ibelow)
  REAL(r_std), SAVE, DIMENSION(nparts) :: CN = &    !! C/N ratio of each plant pool (0-100, unitless)
       & (/ 40., 40., 40., 40., 40., 40., 40., 40. /) 
  REAL(r_std), SAVE, DIMENSION(nparts) :: LC = &    !! Lignin/C ratio of different plant parts (0,22-0,35, unitless)
       & (/ 0.22, 0.35, 0.35, 0.35, 0.35, 0.22, 0.22, 0.22 /)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: metabolic_ref_frac = 0.85    !! used by litter and soilcarbon (0-1, unitless)
  REAL(r_std), SAVE :: metabolic_LN_ratio = 0.018   !! (0-1, unitless)   
  REAL(r_std), SAVE :: tau_metabolic = 0.066        !!
  REAL(r_std), SAVE :: tau_struct = 0.245           !!
  REAL(r_std), SAVE :: soil_Q10 = 0.69              !!= ln 2
  REAL(r_std), SAVE :: tsoil_ref = 30.              !!
  REAL(r_std), SAVE :: litter_struct_coef = 3.      !! 
  REAL(r_std), SAVE, DIMENSION(3) :: moist_coeff = (/ 1.1,  2.4,  0.29 /) !!



  !
  ! stomate_lpj.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: frac_turnover_daily = 0.55  !! (0-1, unitless)



  !
  ! stomate_npp.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tax_max = 0.8 !! Maximum fraction of allocatable biomass used 
                                     !! for maintenance respiration (0-1, unitless)



  !
  ! stomate_phenology.f90
  !

  ! 1. Scalar

  LOGICAL, SAVE :: always_init = .FALSE.           !! take carbon from atmosphere if carbohydrate reserve too small? (true/false)
  REAL(r_std), SAVE :: min_growthinit_time = 300.  !! minimum time since last beginning of a growing season (days)
  REAL(r_std), SAVE :: moiavail_always_tree = 1.0  !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !!  - for trees (0-1, unitless)
  REAL(r_std), SAVE :: moiavail_always_grass = 0.6 !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !! - for grass (0-1, unitless)
  REAL(r_std), SAVE :: t_always                    !! monthly temp. above which temp. tendency doesn't matter
  REAL(r_std), SAVE :: t_always_add = 10.          !! monthly temp. above which temp. tendency doesn't matter (C)

  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: gddncd_ref = 603.           !!
  REAL(r_std), SAVE :: gddncd_curve = 0.0091       !!
  REAL(r_std), SAVE :: gddncd_offset = 64.         !!



  !
  ! stomate_prescribe.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_rescale = 40.       !!



  !
  ! stomate_resp.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: maint_resp_min_vmax = 0.3   !!
  REAL(r_std), SAVE :: maint_resp_coeff = 1.4      !!



  !
  ! stomate_soilcarbon.f90 
  !

  ! 2. Arrays 

  ! 2.1 frac_carb_coefficients

  REAL(r_std), SAVE :: frac_carb_ap = 0.004  !! from active pool: depends on clay content  (0-1, unitless)
                                             !! corresponding to frac_carb(:,iactive,ipassive)
  REAL(r_std), SAVE :: frac_carb_sa = 0.42   !! from slow pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,islow,iactive)
  REAL(r_std), SAVE :: frac_carb_sp = 0.03   !! from slow pool (0-1, unitless) 
                                             !! corresponding to frac_carb(:,islow,ipassive)
  REAL(r_std), SAVE :: frac_carb_pa = 0.45   !! from passive pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,ipassive,iactive)
  REAL(r_std), SAVE :: frac_carb_ps = 0.0    !! from passive pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,ipassive,islow)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: active_to_pass_clay_frac = 0.68  !! (0-1, unitless)
  !! residence times in carbon pools (days)
  REAL(r_std), SAVE :: carbon_tau_iactive = 0.149   !! residence times in active pool (days)
  REAL(r_std), SAVE :: carbon_tau_islow = 5.48      !! residence times in slow pool (days)
  REAL(r_std), SAVE :: carbon_tau_ipassive = 241.   !! residence times in passive pool (days)
  REAL(r_std), SAVE, DIMENSION(3) :: flux_tot_coeff = (/ 1.2, 1.4, .75/)


  !
  ! stomate_turnover.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: new_turnover_time_ref = 20. !!(days)
  REAL(r_std), SAVE :: dt_turnover_time = 10.      !!(days)
  REAL(r_std), SAVE :: leaf_age_crit_tref = 20.    !! (C)
  REAL(r_std), SAVE, DIMENSION(3) :: leaf_age_crit_coeff = (/ 1.5, 0.75, 10./) !! (unitless)



  !
  ! stomate_vmax.f90
  !
 
  ! 1. Scalar

  REAL(r_std), SAVE :: vmax_offset = 0.3        !! minimum leaf efficiency (unitless)
  REAL(r_std), SAVE :: leafage_firstmax = 0.03  !! relative leaf age at which efficiency
                                                !! reaches 1 (unitless)
  REAL(r_std), SAVE :: leafage_lastmax = 0.5    !! relative leaf age at which efficiency
                                                !! falls below 1 (unitless)
  REAL(r_std), SAVE :: leafage_old = 1.         !! relative leaf age at which efficiency
                                                !! reaches its minimum (vmax_offset) 
                                                !! (unitless)

  !
  ! stomate_season.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: gppfrac_dormance = 0.2  !! report maximal GPP/GGP_max for dormance (0-1, unitless)
  REAL(r_std), SAVE :: min_gpp_allowed = 0.3   !! minimum gpp considered as not "lowgpp" (gC.m^{-2}.year^{-1})
  REAL(r_std), SAVE :: tau_climatology = 20.   !! tau for "climatologic variables (years)
  REAL(r_std), SAVE :: hvc1 = 0.019            !! parameters for herbivore activity (unitless)
  REAL(r_std), SAVE :: hvc2 = 1.38             !! parameters for herbivore activity (unitless)
  REAL(r_std), SAVE :: leaf_frac_hvc = 0.33    !! leaf fraction (0-1, unitless)
  REAL(r_std), SAVE :: tlong_ref_max = 303.1   !! maximum reference long term temperature (K)
  REAL(r_std), SAVE :: tlong_ref_min = 253.1   !! minimum reference long term temperature (K)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: ncd_max_year = 3.
  REAL(r_std), SAVE :: gdd_threshold = 5.
  REAL(r_std), SAVE :: green_age_ever = 2.
  REAL(r_std), SAVE :: green_age_dec = 0.5


!!!!Chloe Peatland !!!
!!!! Ajout de la routine stomate_cste_wetland (version HIGH_LAT)
!!!! qui est dans constantes_var dans la nouvelle structure d'orchidee
!!!! MICT3 = ancienne structure d'orchidee, donc a mettre ici dans le module constantes.f90

!From modele Bruno Ringeval (based on Walter model)
!ancien para.h sans ngrid et ntime
  INTEGER(i_std),SAVE  :: nvert = 171
  INTEGER(i_std),SAVE  :: ns = 151
  INTEGER(i_std),SAVE  :: nday = 24 !Chloé change le 24 en 240
  REAL(r_std),SAVE  :: h = 0.1 !delta z (en dm)
  REAL(r_std),SAVE  :: rk = 1 !1 !delta t Attention changer rk (et donc rkh aussi) avec
!  nday si nday=240, rk=0.1, si nday=2400 rk=0.01
  REAL(r_std),SAVE  :: rkh = 100 !rk/h**2 !dépendance delta t/(delta z)^2 !Chloe change 100 en 10
  REAL(r_std),SAVE  :: diffair = 7.2 
  REAL(r_std),SAVE  :: pox = 0.5
  REAL(r_std),SAVE  :: dveg = 0.001
  REAL(r_std),SAVE  :: rkm = 5.0
  REAL(r_std),SAVE  :: xvmax = 20.0
  REAL(r_std),SAVE  :: oxq10 = 2.0
!  REAL(r_std),SAVE  :: catm = 0.0033
  REAL(r_std),SAVE  :: funit = 3.84 !3.84/rk ! Chloe change 3.84 en 38.4
  REAL(r_std),SAVE  :: scmax = 500.
  REAL(r_std),SAVE  :: sr0pl = 600.

!valeur de WTD pour les routines de calcul de densite de flux de CH4
  REAL(r_std),SAVE :: pwater_peat=-3
  !REAL(r_std),PARAMETER :: pwater_wet3=-15
  !REAL(r_std),PARAMETER :: pwater_wet4=-21
!

  REAL(r_std),SAVE :: rpv = 0.5
  REAL(r_std),SAVE :: iother = -1.0
  !!!!plus necessaire maintenant que je mets 2 subroutines differents pour des wt differentes
  !REAL(r_std),PARAMETER :: pwater = 0.0

!!definit vecteurs tve

!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
!!!ATTENTION!!!!
!!!PFT differents selon que l on soit dans Stomate ou Sechiba!!!
!!!a verifier!!!
!!!!!!!!!!!!!!!
! Chloe : Intégré dans pft_parameters :
!  INTEGER(i_std),SAVE,DIMENSION(nvm) :: sdepth_v  = (/0,129,129,129,129,129,129,129,129,79,79,162,162,79/)
!  INTEGER(i_std),SAVE,DIMENSION(nvm) :: rdepth_v  = (/0,64,64,64,64,64,64,64,64,39,39,81,81,39/)
!  INTEGER(i_std),SAVE,DIMENSION(nvm) :: tveg_v  = (/0,1,1,1,1,1,1,1,1,10,10,15,15,10/)

!!rq10 et alpha
!!pour l instant je les mets constantes pour toutes les latitudes
  REAL(r_std),SAVE :: rq10 = 3.0  !Chloe :rq10=3 si depend de carbon et vaut 6 si depend de la NPP
!  REAL(r_std),PARAMETER :: alpha = 0.010
!  REAL(r_std),PARAMETER,DIMENSION(3) :: alpha = (/0.009,0.004,0.021/)
  REAL(r_std),SAVE,DIMENSION(3) :: alpha_CH4 = (/0.004,0.003,0.018/)

!!!!Chloe 
!!!end of bruno's add.



 CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : activate_sub_models
!!
!>\BRIEF         This subroutine reads the flags in the configuration file to
!! activate some sub-models like routing, irrigation, fire, herbivory, ...  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

   SUBROUTINE activate_sub_models(active_flags)

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.1 Input variables

     TYPE(control_type), INTENT(in) :: active_flags     !! What parts of the code are activated ?

     !! 0.4 Local variables

     LOGICAL, SAVE ::  first_call = .TRUE.             !! To keep first call trace (true/false)

!_ ================================================================================================================================

     IF (first_call) THEN 

!!$        IF(active_flags%ok_sechiba .AND. active_flags%river_routing) THEN
!!$           
!!$           !Config Key   = DO_IRRIGATION
!!$           !Config Desc  = Should we compute an irrigation flux 
!!$           !Config If    = RIVER_ROUTING 
!!$           !Config Def   = n
!!$           !Config Help  = This parameters allows the user to ask the model
!!$           !Config         to compute an irigation flux. This performed for the
!!$           !Config         on very simple hypothesis. The idea is to have a good
!!$           !Config         map of irrigated areas and a simple function which estimates
!!$           !Config         the need to irrigate.
!!$           !Config Units = [FLAG]
!!$           CALL getin_p('DO_IRRIGATION', doirrigation)
!!$           !
!!$           !Config Key   = DO_FLOODPLAINS
!!$           !Config Desc  = Should we include floodplains 
!!$           !Config If    = RIVER_ROUTING 
!!$           !Config Def   = n
!!$           !Config Help  = This parameters allows the user to ask the model
!!$           !Config         to take into account the flood plains and return 
!!$           !Config         the water into the soil moisture. It then can go 
!!$           !Config         back to the atmopshere. This tried to simulate 
!!$           !Config         internal deltas of rivers.
!!$           !Config Units = [FLAG]
!!$           CALL getin_p('DO_FLOODPLAINS', dofloodplains)
!!$        
!!$        ENDIF
!!$
           
        IF(active_flags%ok_stomate) THEN

           !Config Key   = HERBIVORES
           !Config Desc  = herbivores allowed?
           !Config If    = OK_STOMATE 
           !Config Def   = n
           !Config Help  = With this variable, you can determine
           !Config         if herbivores are activated
           !Config Units = [FLAG]
           CALL getin_p('HERBIVORES', ok_herbivores)
           !
           !Config Key   = TREAT_EXPANSION
           !Config Desc  = treat expansion of PFTs across a grid cell?
           !Config If    = OK_STOMATE 
           !Config Def   = n
           !Config Help  = With this variable, you can determine
           !Config         whether we treat expansion of PFTs across a
           !Config         grid cell.
           !Config Units = [FLAG]
           CALL getin_p('TREAT_EXPANSION', treat_expansion)
           !
           !Config Key   = LPJ_GAP_CONST_MORT
           !Config Desc  = prescribe mortality if not using DGVM?
           !Config If    = OK_STOMATE
           !Config Def   = y
           !Config Help  = set to TRUE if constant mortality is to be activated
           !Config         ignored if DGVM=true!
           !Config Units = [FLAG]
           CALL getin_p('LPJ_GAP_CONST_MORT', lpj_gap_const_mort)
           !
           !Config Key   = HARVEST_AGRI
           !Config Desc  = Harvest model for agricultural PFTs.
           !Config If    = OK_STOMATE 
           !Config Def   = y
           !Config Help  = Compute harvest above ground biomass for agriculture.
           !Config         Change daily turnover.
           !Config Units = [FLAG]
           CALL getin_p('HARVEST_AGRI', harvest_agri)
           !
           !Config Key   = FIRE_DISABLE
           !Config Desc  = no fire allowed
           !Config If    = OK_STOMATE 
           !Config Def   = n
           !Config Help  = With this variable, you can allow or not
           !Config         the estimation of CO2 lost by fire
           !Config Units = [FLAG]
           CALL getin_p('FIRE_DISABLE', disable_fire)

        ENDIF

        !
        ! Check consistency (see later)
        !
!!$        IF(.NOT.(ok_routing) .AND. (doirrigation .OR. dofloodplains)) THEN
!!$           CALL ipslerr (2,'activate_sub_models', &
!!$               &     'Problem :you tried to activate the irrigation and floodplains without activating the routing',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF
       
!!$        IF(.NOT.(ok_stomate) .AND. (ok_herbivores .OR. treat_expansion .OR. lpj_gap_const_mort &
!!$            & .OR. harvest_agri .OR. disable_fire)) THEN
!!$          CALL ipslerr (2,'activate_sub_models', &
!!$               &     'Problem : try to activate the following options : herbivory, treat_expansion, fire,',&
!!$               &     'harvest_agri and constant mortality without stomate activated.',&
!!$               &     '(check your parameters).')
!!$        ENDIF
            
        first_call =.FALSE.

     ENDIF

   END SUBROUTINE activate_sub_models
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : veget_config
!!
!>\BRIEF         This subroutine reads the flags controlling the configuration for
!! the vegetation : impose_veg, veget_mpa, lai_map, etc...       
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE veget_config

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.4 Local variables  

     LOGICAL, SAVE ::  first_call = .TRUE.        !! To keep first call trace (true/false)  

!_ ================================================================================================================================
     
     IF (first_call) THEN 

        !Config Key   = AGRICULTURE
        !Config Desc  = agriculture allowed?
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = y
        !Config Help  = With this variable, you can determine
        !Config         whether agriculture is allowed
        !Config Units = [FLAG]
        CALL getin_p('AGRICULTURE', agriculture)
        !
        !Config Key   = IMPOSE_VEG
        !Config Desc  = Should the vegetation be prescribed ?
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = n
        !Config Help  = This flag allows the user to impose a vegetation distribution
        !Config         and its characteristics. It is espacially interesting for 0D
        !Config         simulations. On the globe it does not make too much sense as
        !Config         it imposes the same vegetation everywhere
        !Config Units = [FLAG]
        CALL getin_p('IMPOSE_VEG', impveg)

        IF(impveg) THEN
           !Config Key   = IMPOSE_SOILT
           !Config Desc  = Should the soil type be prescribed ?
           !Config Def   = n
           !Config If    = IMPOSE_VEG
           !Config Help  = This flag allows the user to impose a soil type distribution.
           !Config         It is espacially interesting for 0D
           !Config         simulations. On the globe it does not make too much sense as
           !Config         it imposes the same soil everywhere
           !Config Units = [FLAG]
           CALL getin_p('IMPOSE_SOILT', impsoilt)     
        ENDIF

        !Config Key   = LAI_MAP
        !Config Desc  = Read the LAI map
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = n
        !Config Help  = It is possible to read a 12 month LAI map which will
        !Config         then be interpolated to daily values as needed.
        !Config Units = [FLAG]
        CALL getin_p('LAI_MAP',read_lai)

        IF(read_lai) THEN
           !Config Key   = SLOWPROC_LAI_OLD_INTERPOL
           !Config Desc  = Flag to use old "interpolation" of LAI
           !Config If    = LAI_MAP
           !Config Def   = n
           !Config Help  = If you want to recover the old (ie orchidee_1_2 branch) 
           !Config         "interpolation" of LAI map.
           !Config Units = [FLAG]
           CALL getin_p('SLOWPROC_LAI_OLD_INTERPOL',old_lai)
        ENDIF
 
        !Config Key   = LAND_USE
        !Config Desc  = Read a land_use vegetation map
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = y
        !Config Help  = pft values are needed, max time axis is 293
        !Config Units = [FLAG]
        CALL getin_p('LAND_USE',land_use)

        IF(land_use) THEN
           !Config Key   = VEGET_REINIT
           !Config Desc  = booleen to indicate that a new LAND USE file will be used.
           !Config If    = LAND_USE
           !Config Def   = y
           !Config Help  = The parameter is used to bypass veget_year count 
           !Config         and reinitialize it with VEGET_YEAR parameter.
           !Config         Then it is possible to change LAND USE file.
           !Config Units = [FLAG] 
           CALL getin_p('VEGET_REINIT', veget_reinit)
           !
           !Config Key   = LAND_COVER_CHANGE
           !Config Desc  = treat land use modifications
           !Config If    = LAND_USE
           !Config Def   = n
           !Config Help  = With this variable, you can use a Land Use map
           !Config         to simulate anthropic modifications such as 
           !Config         deforestation.
           !Config Units = [FLAG] 
           CALL getin_p('LAND_COVER_CHANGE', lcchange)
           !
           !Config Key   = VEGET_YEAR
           !Config Desc  = Year of the land_use vegetation map to be read
           !Config If    = LAND_USE
           !Config Def   = 1
           !Config Help  = First year for landuse vegetation (2D map by pft).
           !Config         If VEGET_YEAR is set to 0, this means there is no time axis.
           !Config Units = [FLAG] 
           CALL getin_p('VEGET_YEAR', veget_year_orig)
        ENDIF

        IF(.NOT. impveg .AND. .NOT. land_use) THEN
           !Config Key   = SLOWPROC_VEGET_OLD_INTERPOL
           !Config Desc  = Flag to use old "interpolation" of vegetation map.
           !Config If    = NOT(IMPOSE_VEG) and NOT(LAND_USE)
           !Config Def   = n
           !Config Help  = If you want to recover the old (ie orchidee_1_2 branch) 
           !Config         "interpolation" of vegetation map.
           !Config Units = [FLAG] 
           CALL getin_p('SLOWPROC_VEGET_OLD_INTERPOL',old_veget)
         ENDIF  

         !
         ! Check consistency
         !
         ! 1. You have to activate agriculture and land_use
         IF ( .NOT. agriculture .AND. land_use ) THEN 
            CALL ipslerr (2,'veget_config', &
                 &     'Problem with agriculture desactivated and Land Use activated.',&
                 &     'Are you sure ?', &
                 &     '(check your parameters).')
         ENDIF


        first_call = .FALSE.

     ENDIF

!!$        ! DS : Add warning in case of a wrong configuration (need to be discussed)
!!$        ! 2. 
!!$        IF (.NOT.(read_lai) .AND. old_lai) THEN
!!$           CALL ipslerr (2,'veget_config', &
!!$               &     'Problem with lai_map desactivated and old_lai activated.',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF
!!$    
!!$        ! 3.
!!$        IF ((impveg .OR. land_use) .AND. old_veget) THEN
!!$           CALL ipslerr (2,'veget_config', &
!!$                &     'Problem : try to use the old interpolation with a land use map or in impose_veg.',&
!!$                &     'Are you sure ?', &
!!$                &     '(check your parameters).')
!!$        ENDIF
!!$
!!$        ! 4.
!!$        IF ( .NOT.(impveg) .AND. impsoilt) THEN
!!$           CALL ipslerr (2,'veget_config', &
!!$               &     'Problem : try to activate impose_soilt without activating impose_veg.',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF
!!$
!!$        ! 5.
!!$        IF (.NOT.(land_use) .AND. (veget_reinit)) THEN
!!$           CALL ipslerr (2,'veget_config', &
!!$                &     'Problem : try to use a land_use map without activating land_use.',&
!!$                &     'Are you sure ?', &
!!$                &     '(check your parameters).')        
!!$        ENDIF
!!$
!!$        ! 6.
!!$        IF (.NOT.(land_use) .AND. lcchange) THEN
!!$           CALL ipslerr (2,'veget_config', &
!!$                &     'Problem : lcchange is activated without activating land_use.',&
!!$                &     'Are you sure ?', &
!!$                &     '(check your parameters).')        
!!$        ENDIF
           
   END SUBROUTINE veget_config
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : veget_config
!!
!>\BRIEF         This subroutine reads in the configuration file the imposed values of the parameters for all SECHIBA modules.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE config_sechiba_parameters

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.4 Local variables
     
     LOGICAL, SAVE ::  first_call = .TRUE.    !! To keep first call trace (true/false)

!_ ================================================================================================================================
     
     IF(first_call) THEN 
        
        ! Global : parameters used by many modules
        !
        !Config Key   = MAXMASS_GLACIER
        !Config Desc  = The maximum mass of a glacier
        !Config If    = OK_SECHIBA or HYDROL_CWRR
        !Config Def   = 3000.
        !Config Help  = 
        !Config Units = [kg/m^2]  
        CALL getin_p('MAXMASS_GLACIER',maxmass_glacier)
        !
        !Config Key   = SNOWCRI
        !Config Desc  = Sets the amount above which only sublimation occures 
        !Config If    = OK_SECHIBA or HYDROL_CWRR
        !Config Def   = 1.5
        !Config Help  = 
        !Config Units = [kg/m^2]  
        CALL getin_p('SNOWCRI',snowcri)
        !
        !! Initialization of sneige
        sneige = snowcri/mille
!!$        !
!!$        !Config Key   = SECHIBA_QSINT 
!!$        !Config Desc  = Interception reservoir coefficient
!!$        !Config If    = OK_SECHIBA 
!!$        !Config Def   = 0.1
!!$        !Config Help  = Transforms leaf area index into size of interception reservoir
!!$        !Config         for slowproc_derivvar or stomate
!!$        !Config Units = [m]
!!$        CALL getin_p('SECHIBA_QSINT',qsintcst)
!!$        !
!!$        !Config Key   = HYDROL_SOIL_DEPTH
!!$        !Config Desc  = Total depth of soil reservoir
!!$        !Config If    = OK_SECHIBA 
!!$        !Config Def   = 4.
!!$        !Config Help  =
!!$        !Config Units = [m]
!!$        CALL getin_p("HYDROL_SOIL_DEPTH",dpu_max)
        !
        !Config Key   = MIN_WIND
        !Config Desc  = Minimum wind speed
        !Config If    = OK_SECHIBA
        !Config Def   = 0.1
        !Config Help  = 
        !Config Units = [m/s]
        CALL getin_p('MIN_WIND',min_wind)
        !
        !Config Key   = MAX_SNOW_AGE
        !Config Desc  = Maximum period of snow aging 
        !Config If    = OK_SECHIBA
        !Config Def   = 50.
        !Config Help  = 
        !Config Units = [days?]
        CALL getin_p('MAX_SNOW_AGE',max_snow_age)
        !
        !Config Key   = SNOW_TRANS
        !Config Desc  = Transformation time constant for snow
        !Config If    = OK_SECHIBA
        !Config Def   = 0.3
        !Config Help  = 
        !Config Units = [m]   
        CALL getin_p('SNOW_TRANS',snow_trans)
        !
!!$        !Config Key   = MX_EAU_EAU
!!$        !Config Desc  = Maximum quantity of water 
!!$        !Config If    = OK_SECHIBA 
!!$        !Config Def   = 150.
!!$        !Config Help  = 
!!$        !Config Units = [kg/m^3]  
!!$        CALL getin_p('MX_EAU_EAU',mx_eau_eau)
        !-
        ! condveg
        !-
        !
        !Config Key   = Z0_OVER_HEIGHT
        !Config Desc  = to get z0 from height 
        !Config If    = OK_SECHIBA 
        !Config Def   = 1/16.
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('Z0_OVER_HEIGHT',z0_over_height)
        !
        !Config Key   = HEIGHT_DISPLACEMENT
        !Config Desc  = Magic number which relates the height to the displacement height.
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.75
        !Config Help  = 
        !Config Units = [m]  
        CALL getin_p('HEIGHT_DISPLACEMENT',height_displacement)
        !
        !Config Key   = Z0_BARE
        !Config Desc  = bare soil roughness length
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.01 
        !Config Help  = 
        !Config Units = [m]   
        CALL getin_p('Z0_BARE',z0_bare)
        !
        !Config Key   = Z0_ICE
        !Config Desc  = ice roughness length
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.001
        !Config Help  = 
        !Config Units = [m]   
        CALL getin_p('Z0_ICE',z0_ice)
        !
        !Config Key   = TCST_SNOWA
        !Config Desc  = Time constant of the albedo decay of snow
        !Config If    = OK_SECHIBA 
        !Config Def   = 5.0 
        !Config Help  = 
        !Config Units = [days]
        CALL getin_p('TCST_SNOWA',tcst_snowa)
        !
        !Config Key   = SNOWCRI_ALB
        !Config Desc  = Critical value for computation of snow albedo
        !Config If    = OK_SECHIBA
        !Config Def   = 10. 
        !Config Help  = 
        !Config Units = [kg/m^2]  
        CALL getin_p('SNOWCRI_ALB',snowcri_alb)
        !
        !
        !Config Key   = VIS_DRY
        !Config Desc  = The correspondance table for the soil color numbers and their albedo 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('VIS_DRY',vis_dry)
        !
        !Config Key   = NIR_DRY
        !Config Desc  = The correspondance table for the soil color numbers and their albedo 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.48, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('NIR_DRY',nir_dry)
        !
        !Config Key   = VIS_WET 
        !Config Desc  = The correspondance table for the soil color numbers and their albedo
        !Config If    = OK_SECHIBA  
        !Config Def   = 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('VIS_WET',vis_wet)
        !
        !Config Key   = NIR_WET
        !Config Desc  = The correspondance table for the soil color numbers and their albedo 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31
        !Config Help  = 
        !Config Units = [-]    
        CALL getin_p('NIR_WET',nir_wet)
        !
        !Config Key   = ALBSOIL_VIS
        !Config Desc  = 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('ALBSOIL_VIS',albsoil_vis)
        !
        !Config Key   = ALBSOIL_NIR 
        !Config Desc  = 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('ALBSOIL_NIR',albsoil_nir)
        !-
        !
        !Config Key   = ALB_DEADLEAF 
        !Config Desc  = albedo of dead leaves, VIS+NIR 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.12, 0.35
        !Config Help  = 
        !Config Units = [-]     
        CALL getin_p('ALB_DEADLEAF',alb_deadleaf)
        !
        !Config Key   = ALB_ICE
        !Config Desc  = albedo of ice, VIS+NIR
        !Config If    = OK_SECHIBA
        !Config Def   = 0.60, 0.20
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('ALB_ICE',alb_ice)
        !
        ! Get the fixed snow albedo if needed
        !
        !Config Key   = CONDVEG_SNOWA
        !Config Desc  = The snow albedo used by SECHIBA
        !Config Def   = 1.E+20
        !Config if    = OK_SECHIBA
        !Config Help  = This option allows the user to impose a snow albedo.
        !Config         Default behaviour is to use the model of snow albedo
        !Config         developed by Chalita (1993).
        !Config Units = [-]
        CALL getin_p('CONDVEG_SNOWA',fixed_snow_albedo)
        !
        !Config Key   = ALB_BARE_MODEL
        !Config Desc  = Switch bare soil albedo dependent (if TRUE) on soil wetness
        !Config Def   = n
        !Config if    = OK_SECHIBA
        !Config Help  = If TRUE, the model for bare soil albedo is the old formulation.
        !Config         Then it depend on the soil dry or wetness. If FALSE, it is the 
        !Config         new computation that is taken, it is the mean of soil albedo.
        !Config Units = [FLAG]
        CALL getin_p('ALB_BARE_MODEL',alb_bare_model)
        !
        !Config Key   = Z0CDRAG_AVE
        !Config Desc  = Average method for z0
        !Config Def   = y
        !Config if    = OK_SECHIBA
        !Config Help  = If this flag is set to true (y) then the neutral Cdrag
        !Config         is averaged instead of the log(z0). This should be
        !Config         the prefered option. We still wish to keep the other
        !Config         option so we can come back if needed. If this is
        !Config         desired then one should set Z0CDRAG_AVE=n
        !Config Units = [FLAG]
        CALL getin_p('Z0CDRAG_AVE',z0cdrag_ave)
        !
        !Config Key   = IMPOSE_AZE
        !Config Desc  = Should the surface parameters be prescribed
        !Config Def   = n
        !Config if    = OK_SECHIBA
        !Config Help  = This flag allows the user to impose the surface parameters
        !Config         (Albedo Roughness and Emissivity). It is espacially interesting for 0D
        !Config         simulations. On the globe it does not make too much sense as
        !Config         it imposes the same vegetation everywhere
        !Config Units = [FLAG]
        CALL getin_p('IMPOSE_AZE',impaze)
        !
        IF(impaze) THEN
           !
           !Config Key   = CONDVEG_Z0
           !Config Desc  = Surface roughness
           !Config Def   = 0.15
           !Config If    = IMPOSE_AZE
           !Config Help  = Surface rougness to be used on the point if a 0-dim version
           !Config         of SECHIBA is used. Look at the description of the forcing  
           !Config         data for the correct value.
           !Config Units = [m]
           CALL getin_p('CONDVEG_Z0', z0_scal)  
           !
           !Config Key   = ROUGHHEIGHT
           !Config Desc  = Height to be added to the height of the first level
           !Config Def   = 0.0
           !Config If    = IMPOSE_AZE
           !Config Help  = ORCHIDEE assumes that the atmospheric level height is counted
           !Config         from the zero wind level. Thus to take into account the roughness
           !Config         of tall vegetation we need to correct this by a certain fraction
           !Config         of the vegetation height. This is called the roughness height in
           !Config         ORCHIDEE talk.
           !Config Units = [m] 
           CALL getin_p('ROUGHHEIGHT', roughheight_scal)
           ! 
           !Config Key   = CONDVEG_ALBVIS
           !Config Desc  = SW visible albedo for the surface
           !Config Def   = 0.25
           !Config If    = IMPOSE_AZE
           !Config Help  = Surface albedo in visible wavelengths to be used 
           !Config         on the point if a 0-dim version of SECHIBA is used. 
           !Config         Look at the description of the forcing data for 
           !Config         the correct value.
           !Config Units = [-]
           CALL getin_p('CONDVEG_ALBVIS', albedo_scal(ivis))
           !
           !Config Key   = CONDVEG_ALBNIR
           !Config Desc  = SW near infrared albedo for the surface
           !Config Def   = 0.25
           !Config If    = IMPOSE_AZE
           !Config Help  = Surface albedo in near infrared wavelengths to be used 
           !Config         on the point if a 0-dim version of SECHIBA is used. 
           !Config         Look at the description of the forcing data for 
           !Config         the correct value.
           !Config Units = [-]  
           CALL getin_p('CONDVEG_ALBNIR', albedo_scal(inir))
           !
           !Config Key   = CONDVEG_EMIS
           !Config Desc  = Emissivity of the surface for LW radiation
           !Config Def   = 1.0
           !Config If    = IMPOSE_AZE
           !Config Help  = The surface emissivity used for compution the LE emission
           !Config         of the surface in a 0-dim version. Values range between 
           !Config         0.97 and 1.. The GCM uses 0.98.
           !Config Units = [-] 
           CALL getin_p('CONDVEG_EMIS', emis_scal)
        ENDIF
        !
        !-
        ! diffuco 
        !-
        !
        !Config Key   = NLAI
        !Config Desc  = Number of LAI levels
        !Config If    = OK_SECHIBA
        !Config Def   = 20
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('NLAI',nlai)
        !
        !Config Key   = LAIMAX
        !Config Desc  = Maximum LAI
        !Config If    = OK_SECHIBA
        !Config Def   = 
        !Config Help  = 
        !Config Units = [m^2/m^2]   
        CALL getin_p('LAIMAX',laimax)
        !
        !Config Key   = XC4_1 
        !Config Desc  = Factor in the first Collatz equation for C4 plants 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.83
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('XC4_1',xc4_1)
        !
        !Config Key   = XC4_2
        !Config Desc  = Factor in the second Collatz equation for C4 plants
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.93
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('XC4_2',xc4_2)
        !
        !Config Key   = DEW_VEG_POLY_COEFF
        !Config Desc  = coefficients of the polynome of degree 5 for the dew
        !Config If    = OK_SECHIBA
        !Config Def   = 0.887773, 0.205673, 0.110112, 0.014843, 0.000824, 0.000017 
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('DEW_VEG_POLY_COEFF',dew_veg_poly_coeff)
        !-
        ! slowproc
        !-
        !
        !Config Key   = CLAYFRACTION_DEFAULT
        !Config Desc  = default fraction of clay
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.2 
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('CLAYFRACTION_DEFAULT',clayfraction_default)
        !
        !Config Key   = MIN_VEGFRAC 
        !Config Desc  = Minimal fraction of mesh a vegetation type can occupy 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.001 
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('MIN_VEGFRAC',min_vegfrac)
        !
        !Config Key   = STEMPDIAG_BID 
        !Config Desc  = only needed for an initial LAI if there is no restart file
        !Config If    = OK_SECHIBA 
        !Config Def   = 280.
        !Config Help  = 
        !Config Units = [K]
        CALL getin_p('STEMPDIAG_BID',stempdiag_bid)
        !
        first_call =.FALSE.
        
     ENDIF
     
   END SUBROUTINE config_sechiba_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_co2_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when OK_CO2 is set to true. (ie : when the photosynthesis is activated) 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE config_co2_parameters
     
     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.4 Local variables
     
     LOGICAL, SAVE ::  first_call = .TRUE.      !! To keep first call trace (true/false)

!_ ================================================================================================================================
     
     IF(first_call) THEN
        
        !
        !Config Key   = LAI_LEVEL_DEPTH
        !Config Desc  = 
        !Config If    = OK_CO2
        !Config Def   = 0.15
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('LAI_LEVEL_DEPTH',lai_level_depth)
        !
        !Config Key   = X1_COEF
        !Config Desc  = Multiplicative factor in the equation defining kt 
        !Config If    = OK_CO2
        !Config Def   = 0.177
        !Config Help  = Multiplicative factor for calculating the pseudo first order rate constant 
        !Config         of assimilation response to co2 kt
        !Config Units = [-]  
        CALL getin_p('X1_COEF',x1_coef)
        !
        !Config Key   = X1_Q10
        !Config Desc  = Exponential factor in the equation defining kt
        !Config If    = OK_CO2
        !Config Def   = 0.069
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('X1_Q10',x1_Q10)
        !
        !Config Key   = QUANTUM_YIELD
        !Config Desc  = 
        !Config If    = OK_CO2
        !Config Def   = 0.092 
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('QUANTUM_YIELD',quantum_yield)
        !
        !Config Key   = KT_COEF
        !Config Desc  = Multiplicative factor in the equation defining kt
        !Config If    = OK_CO2
        !Config Def   = 0.7 
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('KT_COEF',kt_coef)
        !
        !Config Key   = KC_COEF
        !Config Desc  = Multiplicative factor for calculating Kc
        !Config If    = OK_CO2
        !Config Def   = 39.09 
        !Config Help  = Multiplicative factor for calculating the Michaelis-Menten
        !Config         coefficient Kc  
        !Config Units = [-]  
        CALL getin_p('KC_COEF',kc_coef)
        !
        !Config Key   = KO_Q10 
        !Config Desc  = Exponential factor for calculating Kc and Ko 
        !Config If    = OK_CO2 
        !Config Def   = 0.085 
        !Config Help  = Exponential factor for calculating the Michaelis-Menten coefficients
        !Config         Kc and Ko 
        !Config Units = [-]  
        CALL getin_p('KO_Q10',Ko_Q10)
        !
        !Config Key   = OA
        !Config Desc  = Intercellular concentration of O2
        !Config If    = OK_CO2
        !Config Def   = 210000. 
        !Config Help  = 
        !Config Units = [ppm]  
        CALL getin_p('OA',Oa)
        !
        !Config Key   = KO_COEF
        !Config Desc  = Multiplicative factor for calculating Ko
        !Config If    = OK_CO2
        !Config Def   = 2.412 
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('KO_COEF',Ko_coef)
        !
        !Config Key   = CP_0
        !Config Desc  = Multiplicative factor for calculating the CO2 compensation point
        !Config If    = OK_CO2
        !Config Def   = 42. 
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('CP_0',CP_0)
        !
        !Config Key   = CP_TEMP_COEF
        !Config Desc  = Exponential factor for calculating the CO2 compensation point
        !Config If    = OK_CO2
        !Config Def   = 9.46
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('CP_TEMP_COEF',cp_temp_coef)
        !
        !Config Key   = CP_TEMP_REF
        !Config Desc  = Reference temperature for the CO2 compensation point CP
        !Config If    = OK_CO2
        !Config Def   = 25.
        !Config Help  = 
        !Config Units = [C]  
        CALL getin_p('CP_TEMP_REF',cp_temp_ref)
        !
        !Config Key   = RT_COEF
        !Config Desc  = 
        !Config If    = OK_CO2
        !Config Def   = 0.8, 1.3
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('RT_COEF',rt_coef)
        !
        !Config Key   = VC_COEF
        !Config Desc  = 
        !Config If    = OK_CO2
        !Config Def   = 0.39, 0.3 
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('VC_COEF',vc_coef)
        
        first_call =.FALSE.
        
     ENDIF
     
   END SUBROUTINE config_co2_parameters
!
!=
!

!!$!! ================================================================================================================================
!!$!! SUBROUTINE   : config_hydrolc_parameters
!!$!!
!!$!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!!$!! needed when the Choisnel hydrology model is activated.
!!$!!
!!$!! DESCRIPTION  : None
!!$!!
!!$!! RECENT CHANGE(S): None
!!$!!
!!$!! MAIN OUTPUT VARIABLE(S): 
!!$!!
!!$!! REFERENCE(S) :
!!$!!
!!$!! FLOWCHART    :
!!$!! \n
!!$!_ ================================================================================================================================
!!$
!!$   SUBROUTINE config_hydrolc_parameters
!!$     
!!$     IMPLICIT NONE
!!$
!!$    !! 0. Variables and parameters declaration
!!$
!!$    !! 0.4 Local variables 
!!$
!!$     LOGICAL, SAVE ::  first_call = .TRUE.     !! To keep first call trace (true/false)
!!$
!!$!_ ================================================================================================================================
!!$     
!!$     IF(first_call) THEN 
!!$        !
!!$        !Config Key   = QWILT
!!$        !Config Desc  = Wilting point
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
!!$        !Config Def   = 5.0
!!$        !Config Help  = Has a numerical role for the moment
!!$        !Config Units = [-]
!!$        CALL getin_p('QWILT',qwilt)
!!$        !
!!$        !Config Key   = MIN_RESDIS
!!$        !Config Desc  = The minimal size we allow for the upper reservoir
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
!!$        !Config Def   = 2.e-5
!!$        !Config Help  = 
!!$        !Config Units = [m]
!!$        CALL getin_p('MIN_RESDIS',min_resdis)
!!$        !
!!$        !Config Key   = MIN_DRAIN
!!$        !Config Desc  = Diffusion constant for the slow regime
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
!!$        !Config Def   = 0.001
!!$        !Config Help  = 
!!$        !Config Units = [kg/m^2/dt]
!!$        CALL getin_p('MIN_DRAIN',min_drain)
!!$        !
!!$        !Config Key   = MAX_DRAIN
!!$        !Config Desc  = Diffusion constant for the fast regime
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
!!$        !Config Def   = 0.1
!!$        !Config Help  = 
!!$        !Config Units = [kg/m^2/dt]
!!$        CALL getin_p('MAX_DRAIN',max_drain)
!!$        !
!!$        !Config Key   = EXP_DRAIN
!!$        !Config Desc  = The exponential in the diffusion law
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
!!$        !Config Def   = 1.5
!!$        !Config Help  = 
!!$        !Config Units = [-]
!!$        CALL getin_p('EXP_DRAIN',exp_drain)
!!$        !
!!$        !Config Key   = RSOL_CSTE
!!$        !Config Desc  = Constant in the computation of resistance for bare  soil evaporation 
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
!!$        !Config Def   = 33.E3
!!$        !Config Help  = 
!!$        !Config Units = [s/m^2]
!!$        CALL getin_p('RSOL_CSTE',rsol_cste)
!!$        !
!!$        !Config Key   = HCRIT_LITTER
!!$        !Config Desc  = Scaling depth for litter humidity
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR) 
!!$        !Config Def   = 0.08 
!!$        !Config Help  = 
!!$        !Config Units = [m]
!!$        CALL getin_p('HCRIT_LITTER',hcrit_litter)
!!$        !
!!$        !Config Key   = HYDROL_OK_HDIFF
!!$        !Config Desc  = do horizontal diffusion?
!!$        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)  
!!$        !Config Def   = n
!!$        !Config Help  = If TRUE, then water can diffuse horizontally between
!!$        !Config         the PFTs' water reservoirs.
!!$        !Config Units = [FLAG]
!!$        CALL getin_p('HYDROL_OK_HDIFF',ok_hdiff)         
!!$
!!$        first_call =.FALSE.
!!$        
!!$     ENDIF
!!$     
!!$   END SUBROUTINE config_hydrolc_parameters
!!$   
!!$!
!!$!=
!!$!
!!$
!!$!! ================================================================================================================================
!!$!! SUBROUTINE   : config_hydrol_cwrr_parameters
!!$!!
!!$!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!!$!! needed when the 11-layers hydrology model is activated.
!!$!!
!!$!! DESCRIPTION  : None
!!$!!
!!$!! RECENT CHANGE(S): None
!!$!!
!!$!! MAIN OUTPUT VARIABLE(S): 
!!$!!
!!$!! REFERENCE(S) :
!!$!!
!!$!! FLOWCHART    :
!!$!! \n
!!$!_ ================================================================================================================================
!!$
!!$   SUBROUTINE config_hydrol_cwrr_parameters
!!$     
!!$     IMPLICIT NONE
!!$
!!$     !! 0. Variables and parameters declaration
!!$
!!$     !! 0.4 Local variables    
!!$
!!$     LOGICAL, SAVE ::  first_call = .TRUE.       !! To keep first call trace (true/false)
!!$
!!$!_ ================================================================================================================================
!!$     
!!$     IF (first_call) THEN
!!$
!!$        !
!!$        !Config Key   = W_TIME 
!!$        !Config Desc  = Time weighting for discretisation
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 1.
!!$        !Config Help  = 
!!$        !Config Units = [-]
!!$        CALL getin_p('W_TIME',w_time)
!!$        !
!!$        !Config Key   = NVAN
!!$        !Config Desc  = Van genuchten coefficient n
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 1.89, 1.56, 1.31
!!$        !Config Help  = 
!!$        !Config Units = [-]
!!$        CALL getin_p('NVAN',nvan)
!!$        !
!!$        !Config Key   = AVAN
!!$        !Config Desc  = Van genuchten coefficient a
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.0075, 0.0036, 0.0019
!!$        !Config Help  = 
!!$        !Config Units = [1/mm]  
!!$        CALL getin_p('AVAN',avan)
!!$        !
!!$        !Config Key   = MCR 
!!$        !Config Desc  = Residual soil water content
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.065, 0.078, 0.095
!!$        !Config Help  = 
!!$        !Config Units = [mm]  
!!$        CALL getin_p('MCR',mcr)
!!$        !
!!$        !Config Key   = MCS
!!$        !Config Desc  = Saturated soil water content
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.41, 0.43, 0.41
!!$        !Config Help  = 
!!$        !Config Units = [-]  
!!$        CALL getin_p('MCS',mcs)     
!!$        !
!!$        !Config Key   = KS 
!!$        !Config Desc  = Hydraulic conductivity Saturation
!!$        !Config If    = HYDROL_CWRR 
!!$        !Config Def   = 1060.8, 249.6, 62.4
!!$        !Config Help  = 
!!$        !Config Units = [mm/d]   
!!$        CALL getin_p('KS',ks)
!!$        !
!!$        !Config Key   = PCENT
!!$        !Config Desc  = Soil moisture above which transpir is max
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.5, 0.5, 0.5
!!$        !Config Help  = 
!!$        !Config Units = [-]    
!!$        CALL getin_p('PCENT',pcent)
!!$        !
!!$        !Config Key   = FREE_DRAIN_MAX 
!!$        !Config Desc  = Max value of the permeability coeff at the bottom of the soil 
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 1.0, 1.0, 1.0 
!!$        !Config Help  = 
!!$        !Config Units = [-]   
!!$        CALL getin_p('FREE_DRAIN_MAX',free_drain_max)
!!$        !
!!$        !Config Key   = MCF 
!!$        !Config Desc  = Volumetric water content field capacity
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.32, 0.32, 0.32
!!$        !Config Help  = 
!!$        !Config Units = [-]   
!!$        CALL getin_p('MCF',mcf)
!!$        !
!!$        !Config Key   = MCW
!!$        !Config Desc  = Volumetric water content Wilting pt
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.10, 0.10, 0.10 
!!$        !Config Help  = 
!!$        !Config Units = [-]   
!!$        CALL getin_p('MCW',mcw)
!!$        !
!!$        !Config Key   = MC_AWET
!!$        !Config Desc  = Vol. wat. cont. above which albedo is cst
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.25, 0.25, 0.25
!!$        !Config Help  = 
!!$        !Config Units = [-]   
!!$        CALL getin_p('MC_AWET',mc_awet)
!!$        !
!!$        !Config Key   = MC_ADRY
!!$        !Config Desc  = Vol. wat. cont. below which albedo is cst
!!$        !Config If    = HYDROL_CWRR
!!$        !Config Def   = 0.1, 0.1, 0.1
!!$        !Config Help  = 
!!$        !Config Units = [-]   
!!$        CALL getin_p('MC_ADRY',mc_adry)
!!$         
!!$        first_call =.FALSE.
!!$        
!!$     ENDIF
!!$
!!$   END SUBROUTINE config_hydrol_cwrr_parameters
!!$!
!!$!=
!!$!
!!$
!!$!! ================================================================================================================================
!!$!! SUBROUTINE   : config_routing_parameters 
!!$!!
!!$!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!!$!! needed when the routing is activated.
!!$!!
!!$!! DESCRIPTION  : None
!!$!!
!!$!! RECENT CHANGE(S): None
!!$!!
!!$!! MAIN OUTPUT VARIABLE(S): 
!!$!!
!!$!! REFERENCE(S) :
!!$!!
!!$!! FLOWCHART    :
!!$!! \n
!!$!_ ================================================================================================================================
!!$
!!$   SUBROUTINE config_routing_parameters
!!$     
!!$     IMPLICIT NONE
!!$     
!!$     !! 0. Variables and parameters declaration
!!$     
!!$     !! 0.4 Local variables
!!$
!!$     LOGICAL, SAVE ::  first_call = .TRUE.    !! To keep first call trace (true/false)
!!$
!!$!_ ================================================================================================================================
!!$     
!!$     IF(first_call) THEN
!!$        !
!!$        !Config Key   = CROP_COEF 
!!$        !Config Desc  = Parameter for the Kassel irrigation parametrization linked to the crops
!!$        !Config If    = OK_ROUTING
!!$        !Config Def   = 1.5
!!$        !Config Help  = Empirical crop coefficient dependent on vegetation characteristics
!!$        !Config         according to Kassel irrigation parametrization.
!!$        !Config         When potential transpiration is used this coefficient has another interpretation
!!$        !Config Units = [-]  
!!$        CALL getin_p('CROP_COEF',crop_coef)
!!$        
!!$        first_call =.FALSE.
!!$        
!!$     ENDIF
!!$     
!!$   END SUBROUTINE config_routing_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when stomate is activated (ie : when OK_STOMATE is set to true).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE config_stomate_parameters
     
    IMPLICIT NONE
    
    !! 0. Variables and parameters declaration

    !! 0.4 Local variables   

    LOGICAL, SAVE ::  first_call = .TRUE.  !! To keep first call trace (true/false)

!_ ================================================================================================================================
    
    IF(first_call) THEN
       !-
       ! constraints_parameters
       !-
       !
       !Config Key   = TOO_LONG 
       !Config Desc  = longest sustainable time without regeneration (vernalization)
       !Config If    = OK_STOMATE
       !Config Def   = 5.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TOO_LONG',too_long)

       !-
       ! fire parameters
       !-
       !
       !Config Key   = TAU_FIRE 
       !Config Desc  = Time scale for memory of the fire index (days). Validated for one year in the DGVM. 
       !Config If    = OK_STOMATE 
       !Config Def   = 30.
       !Config Help  = 
       !Config Units = [days]    
       CALL getin_p('TAU_FIRE',tau_fire)
       !
       !Config Key   = LITTER_CRIT
       !Config Desc  = Critical litter quantity for fire
       !Config If    = OK_STOMATE 
       !Config Def   = 200.
       !Config Help  = 
       !Config Units = [gC/m^2]  
       CALL getin_p('LITTER_CRIT',litter_crit)
       !
       !Config Key   = FIRE_RESIST_STRUCT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('FIRE_RESIST_STRUCT',fire_resist_struct)
       !
       !
       !Config Key   = CO2FRAC
       !Config Desc  = What fraction of a burned plant compartment goes into the atmosphere
       !Config If    = OK_STOMATE 
       !Config Def   = 0.95, 0.95, 0., 0.3, 0., 0., 0.95, 0.95
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('CO2FRAC',co2frac)
       !
       !Config Key   = BCFRAC_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3, 1.3, 88.2 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('BCFRAC_COEFF',bcfrac_coeff)
       !
       !Config Key   = FIREFRAC_COEFF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45, 0.8, 0.6, 0.13
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FIREFRAC_COEFF',firefrac_coeff)

       !-
       ! gap parameters (+ lpj_const_mort)
       !-
       !
       !Config Key   = AVAILABILITY_FACT 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('AVAILABILITY_FACT', availability_fact)  
       !
       !Config Key   = REF_GREFF
       !Config Desc  = Asymptotic maximum mortality rate
       !Config If    = OK_STOMATE 
       !Config Def   = 0.035
       !Config Help  = Set asymptotic maximum mortality rate from Sitch 2003
       !Config         (they use 0.01) (year^{-1})
       !Config Units = [1/year]  
       CALL getin_p('REF_GREFF',ref_greff)
       !-
       ! allocation parameters
       !-
       !
       !Config Key   = OK_MINRES
       !Config Desc  = Do we try to reach a minimum reservoir even if we are severely stressed?
       !Config If    = OK_STOMATE 
       !Config Def   = y
       !Config Help  = 
       !Config Units = [FLAG]
       CALL getin_p('OK_MINRES',ok_minres)
       !
       !Config Key   = TAU_LEAFINIT
       !Config Desc  = time to attain the initial foliage using the carbohydrate reserve 
       !Config If    = OK_STOMATE 
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('TAU_LEAFINIT', tau_leafinit)
       !
       !Config Key   = RESERVE_TIME_TREE 
       !Config Desc  = maximum time during which reserve is used (trees) 
       !Config If    = OK_STOMATE 
       !Config Def   = 30.
       !Config Help  = 
       !Config Units = [days]    
       CALL getin_p('RESERVE_TIME_TREE',reserve_time_tree)
       !
       !Config Key   = RESERVE_TIME_GRASS 
       !Config Desc  = maximum time during which reserve is used (grasses) 
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('RESERVE_TIME_GRASS',reserve_time_grass)
       !
       !Config Key   = R0
       !Config Desc  = Standard root allocation 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('R0',R0)
       !
       !Config Key   = S0 
       !Config Desc  = Standard sapwood allocation 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('S0',S0)
       !
       !Config Key   = F_FRUIT
       !Config Desc  = Standard fruit allocation
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('F_FRUIT',f_fruit)
       !
       !Config Key   = ALLOC_SAP_ABOVE_GRASS 
       !Config Desc  = fraction of sapwood allocation above ground 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALLOC_SAP_ABOVE_GRASS',alloc_sap_above_grass)
       !
       !Config Key   = MIN_LTOLSR 
       !Config Desc  = extrema of leaf allocation fraction 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MIN_LTOLSR',min_LtoLSR)
       !
       !Config Key   = MAX_LTOLSR
       !Config Desc  = extrema of leaf allocation fraction
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MAX_LTOLSR',max_LtoLSR)
       !
       !Config Key   = Z_NITROGEN
       !Config Desc  = scaling depth for nitrogen limitation 
       !Config If    = OK_STOMATE
       !Config Def   = 0.2 
       !Config Help  =
       !Config Units = [m]  
       CALL getin_p('Z_NITROGEN',z_nitrogen)
       !
       !Config Key   = LAI_MAX_TO_HAPPY
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('LAI_MAX_TO_HAPPY',lai_max_to_happy)
       !
       !Config Key   = NLIM_TREF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 25. 
       !Config Help  = 
       !Config Units = [C]  
       CALL getin_p('NLIM_TREF',Nlim_tref) 
  
       !-
       ! data parameters
       !-
       !
       !Config Key   = PIPE_TUNE1
       !Config Desc  = crown area = pipe_tune1. stem diameter**(1.6) (Reinicke's theory)
       !Config If    = OK_STOMATE 
       !Config Def   = 100.0
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE1',pipe_tune1)
       !
       !Config Key   = PIPE_TUNE2 
       !Config Desc  = height=pipe_tune2 * diameter**pipe_tune3
       !Config If    = OK_STOMATE 
       !Config Def   = 40.0 
       !Config Help  = 
       !Config Units = [-]      
       CALL getin_p('PIPE_TUNE2',pipe_tune2) 
        !
       !Config Key   = PIPE_TUNE3
       !Config Desc  = height=pipe_tune2 * diameter**pipe_tune3
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE3',pipe_tune3)
       !
       !Config Key   = PIPE_TUNE4
       !Config Desc  = needed for stem diameter
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('PIPE_TUNE4',pipe_tune4)
       !
       !Config Key   = PIPE_DENSITY 
       !Config Desc  = Density
       !Config If    = OK_STOMATE 
       !Config Def   = 2.e5 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('PIPE_DENSITY',pipe_density)
       !
       !Config Key   = PIPE_K1 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 8.e3 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_K1',pipe_k1)
       !
       !Config Key   = PIPE_TUNE_EXP_COEFF 
       !Config Desc  = pipe tune exponential coeff 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.6 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_TUNE_EXP_COEFF',pipe_tune_exp_coeff)
       !
       !
       !Config Key   = PRECIP_CRIT 
       !Config Desc  = minimum precip
       !Config If    = OK_STOMATE 
       !Config Def   = 100.
       !Config Help  = 
       !Config Units = [mm/year]  
       CALL getin_p('PRECIP_CRIT',precip_crit)
       !
       !Config Key   = GDD_CRIT_ESTAB
       !Config Desc  = minimum gdd for establishment of saplings
       !Config If    = OK_STOMATE 
       !Config Def   = 150. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GDD_CRIT_ESTAB',gdd_crit_estab)
        !
       !Config Key   = FPC_CRIT
       !Config Desc  = critical fpc, needed for light competition and establishment
       !Config If    = OK_STOMATE 
       !Config Def   = 0.95
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('FPC_CRIT',fpc_crit)
       !
       !Config Key   = ALPHA_GRASS
       !Config Desc  = sapling characteristics : alpha's
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_GRASS',alpha_grass)
       !
       !Config Key   = ALPHA_TREE
       !Config Desc  = sapling characteristics : alpha's 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_TREE',alpha_tree)
       !-
       !
       !Config Key   = MASS_RATIO_HEART_SAP
       !Config Desc  = mass ratio (heartwood+sapwood)/sapwood
       !Config If    = OK_STOMATE 
       !Config Def   = 3.
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MASS_RATIO_HEART_SAP',mass_ratio_heart_sap)
       !
       !Config Key   = FRAC_GROWTHRESP
       !Config Desc  = fraction of GPP which is lost as growth respiration
       !Config If    = OK_STOMATE 
       !Config Def   = 0.28
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_GROWTHRESP',frac_growthresp)
       !
       !Config Key   = TAU_HUM_MONTH
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('TAU_HUM_MONTH',tau_hum_month)
       !
       !Config Key   = TAU_HUM_WEEK
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 7.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_HUM_WEEK',tau_hum_week)
       !
       !Config Key   = TAU_T2M_MONTH
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20.
       !Config Help  = 
       !Config Units = [days]     
       CALL getin_p('TAU_T2M_MONTH',tau_t2m_month)
       !
       !Config Key   = TAU_T2M_WEEK
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 7.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_T2M_WEEK',tau_t2m_week)
       !
       !Config Key   = TAU_TSOIL_MONTH 
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]     
       CALL getin_p('TAU_TSOIL_MONTH',tau_tsoil_month)
       !
       !Config Key   = TAU_SOILHUM_MONTH
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_SOILHUM_MONTH',tau_soilhum_month)
       !
       !Config Key   = TAU_GPP_WEEK 
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 7. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_GPP_WEEK',tau_gpp_week)
       !
       !Config Key   = TAU_GDD
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 40. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_GDD',tau_gdd)
       !
       !Config Key   = TAU_NGD
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 50.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_NGD',tau_ngd)
       !
       !Config Key   = COEFF_TAU_LONGTERM
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 3. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('COEFF_TAU_LONGTERM',coeff_tau_longterm)
       !-
       !
       !Config Key   = BM_SAPL_CARBRES 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 5. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('BM_SAPL_CARBRES',bm_sapl_carbres)
       !
       !Config Key   = BM_SAPL_SAPABOVE
       !Config Desc  = 
       !Config If    = OK_STOMATE
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('BM_SAPL_SAPABOVE',bm_sapl_sapabove)
       !
       !Config Key   = BM_SAPL_HEARTABOVE 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 2.
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('BM_SAPL_HEARTABOVE',bm_sapl_heartabove)
       !
       !Config Key   = BM_SAPL_HEARTBELOW 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 2. 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('BM_SAPL_HEARTBELOW',bm_sapl_heartbelow)
       !
       !Config Key   = INIT_SAPL_MASS_LEAF_NAT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_LEAF_NAT',init_sapl_mass_leaf_nat)
       !
       !Config Key   = INIT_SAPL_MASS_LEAF_AGRI
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1. 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_LEAF_AGRI',init_sapl_mass_leaf_agri)
       !
       !Config Key   = INIT_SAPL_MASS_CARBRES
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 5. 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_CARBRES',init_sapl_mass_carbres)
       !
       !Config Key   = INIT_SAPL_MASS_ROOT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('INIT_SAPL_MASS_ROOT',init_sapl_mass_root)
       !
       !Config Key   = INIT_SAPL_MASS_FRUIT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_FRUIT',init_sapl_mass_fruit)
       !
       !Config Key   = CN_SAPL_INIT 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('CN_SAPL_INIT',cn_sapl_init)
       !
       !Config Key   = MIGRATE_TREE 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10000.
       !Config Help  = 
       !Config Units = [m/year]   
       CALL getin_p('MIGRATE_TREE',migrate_tree)
       !
       !Config Key   = MIGRATE_GRASS
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10000.
       !Config Help  = 
       !Config Units = [m/year]   
       CALL getin_p('MIGRATE_GRASS',migrate_grass)
       !
       !Config Key   = LAI_INITMIN_TREE
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [m^2/m^2]  
       CALL getin_p('LAI_INITMIN_TREE',lai_initmin_tree)
       !
       !Config Key   = LAI_INITMIN_GRASS 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [m^2/m^2]    
       CALL getin_p('LAI_INITMIN_GRASS',lai_initmin_grass)
       !
       !Config Key   = DIA_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 4., 0.5
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DIA_COEFF',dia_coeff)
       !
       !Config Key   = MAXDIA_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 100., 0.01 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MAXDIA_COEFF',maxdia_coeff)
       !
       !Config Key   = BM_SAPL_LEAF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 4., 4., 0.8, 5. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('BM_SAPL_LEAF',bm_sapl_leaf)

       !-
       ! litter parameters
       !-
       !
       !Config Key   = METABOLIC_REF_FRAC
       !Config Desc  =
       !Config If    = OK_STOMATE 
       !Config Def   = 0.85  
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('METABOLIC_REF_FRAC',metabolic_ref_frac)
       !
       !Config Key   = Z_DECOMP
       !Config Desc  = scaling depth for soil activity
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [m]   
       CALL getin_p('Z_DECOMP',z_decomp)
       !
       !Config Key   = CN
       !Config Desc  = C/N ratio
       !Config If    = OK_STOMATE 
       !Config Def   = 40., 40., 40., 40., 40., 40., 40., 40.
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('CN',CN)
       !
       !Config Key   = LC 
       !Config Desc  = Lignine/C ratio of the different plant parts
       !Config If    = OK_STOMATE 
       !Config Def   = 0.22, 0.35, 0.35, 0.35, 0.35, 0.22, 0.22, 0.22
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('LC',LC)
       !
       !Config Key   = FRAC_SOIL_STRUCT_AA
       !Config Desc  = frac_soil(istructural,iactive,iabove)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.55
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_SOIL_STRUCT_AA',frac_soil_struct_aa)
       !
       !Config Key   = FRAC_SOIL_STRUCT_A 
       !Config Desc  = frac_soil(istructural,iactive,ibelow)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_SOIL_STRUCT_AB',frac_soil_struct_ab)
       !
       !Config Key   = FRAC_SOIL_STRUCT_SA
       !Config Desc  = frac_soil(istructural,islow,iabove)
       !Config If    = OK_STOMATE
       !Config Def   = 0.7  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_STRUCT_SA',frac_soil_struct_sa)
       !
       !Config Key   = FRAC_SOIL_STRUCT_SB
       !Config Desc  = frac_soil(istructural,islow,ibelow) 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.7  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_STRUCT_SB',frac_soil_struct_sb)
       !
       !Config Key   = FRAC_SOIL_METAB_AA 
       !Config Desc  = frac_soil(imetabolic,iactive,iabove) 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_METAB_AA',frac_soil_metab_aa)
       !
       !Config Key   = FRAC_SOIL_METAB_AB 
       !Config Desc  = frac_soil(imetabolic,iactive,ibelow)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_METAB_AB',frac_soil_metab_ab)
       !
       !
       !Config Key   = METABOLIC_LN_RATIO
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.018  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('METABOLIC_LN_RATIO',metabolic_LN_ratio) 
       !
       !Config Key   = TAU_METABOLIC
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.066
       !Config Help  = 
       !Config Units = [days] 
       CALL getin_p('TAU_METABOLIC',tau_metabolic)
       !
       !Config Key   = TAU_STRUCT 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.245 
       !Config Help  = 
       !Config Units = [days]
       CALL getin_p('TAU_STRUCT',tau_struct)
       !
       !Config Key   = SOIL_Q10
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.69 (=ln2)
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('SOIL_Q10',soil_Q10)
       !
       !Config Key   = TSOIL_REF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 30. 
       !Config Help  = 
       !Config Units = [C]   
       CALL getin_p('TSOIL_REF',tsoil_ref)
       !
       !Config Key   = LITTER_STRUCT_COEF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 3. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('LITTER_STRUCT_COEF',litter_struct_coef)
       !
       !Config Key   = MOIST_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.1, 2.4, 0.29
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIST_COEFF',moist_coeff)

       !-
       ! lpj parameters
       !-
       !
       !Config Key   = FRAC_TURNOVER_DAILY 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.55
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_TURNOVER_DAILY',frac_turnover_daily)   

       !-
       ! npp parameters
       !-
       !
       !Config Key   = TAX_MAX
       !Config Desc  = maximum fraction of allocatable biomass used for maintenance respiration
       !Config If    = OK_STOMATE 
       !Config Def   = 0.8
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('TAX_MAX',tax_max) 

       !-
       ! phenology parameters
       !-
       !
       !Config Key   = ALWAYS_INIT
       !Config Desc  = take carbon from atmosphere if carbohydrate reserve too small? 
       !Config If    = OK_STOMATE 
       !Config Def   = n 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALWAYS_INIT',always_init)
       !
       !Config Key   = MIN_GROWTHINIT_TIME 
       !Config Desc  = minimum time since last beginning of a growing season
       !Config If    = OK_STOMATE 
       !Config Def   = 300. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('MIN_GROWTHINIT_TIME',min_growthinit_time)
       !
       !Config Key   = MOIAVAIL_ALWAYS_TREE
       !Config Desc  = moisture availability above which moisture tendency doesn't matter 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIAVAIL_ALWAYS_TREE',moiavail_always_tree)
       !
       !Config Key   = MOIAVAIL_ALWAYS_GRASS 
       !Config Desc  = moisture availability above which moisture tendency doesn't matter
       !Config If    = OK_STOMATE 
       !Config Def   = 0.6 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIAVAIL_ALWAYS_GRASS',moiavail_always_grass)
       !
       !Config Key   = T_ALWAYS_ADD
       !Config Desc  = monthly temp. above which temp. tendency doesn't matter 
       !Config If    = OK_STOMATE 
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [C]    
       CALL getin_p('T_ALWAYS_ADD',t_always_add)
       !
       !
       !Config Key   = GDDNCD_REF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 603. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('GDDNCD_REF',gddncd_ref)
       !
       !Config Key   = GDDNCD_CURVE
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.0091 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GDDNCD_CURVE',gddncd_curve)
       !
       !Config Key   = GDDNCD_OFFSET
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 64. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GDDNCD_OFFSET',gddncd_offset)
       !-
       ! prescribe parameters
       !-
       !
       !Config Key   = BM_SAPL_RESCALE 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 40. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('BM_SAPL_RESCALE',bm_sapl_rescale)

       !-
       ! respiration parameters
       !-
       !
       !Config Key   = MAINT_RESP_MIN_VMAX
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('MAINT_RESP_MIN_VMAX',maint_resp_min_vmax)  
       !
       !Config Key   = MAINT_RESP_COEFF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.4 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('MAINT_RESP_COEFF',maint_resp_coeff)

       !-
       ! soilcarbon parameters 
       !-
       !
       !Config Key   = FRAC_CARB_AP
       !Config Desc  = frac carb coefficients from active pool: depends on clay content
       !Config if    = OK_STOMATE 
       !Config Def   = 0.004
       !Config Help  = fraction of the active pool going into the passive pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_AP',frac_carb_ap)  
       !
       !Config Key   = FRAC_CARB_SA
       !Config Desc  = frac_carb_coefficients from slow pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.42
       !Config Help  = fraction of the slow pool going into the active pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_SA',frac_carb_sa)
       !
       !Config Key   = FRAC_CARB_SP
       !Config Desc  = frac_carb_coefficients from slow pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.03
       !Config Help  = fraction of the slow pool going into the passive pool
       !Config Units = [-] 
       CALL getin_p('FRAC_CARB_SP',frac_carb_sp)
       !
       !Config Key   = FRAC_CARB_PA
       !Config Desc  = frac_carb_coefficients from passive pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.45
       !Config Help  = fraction of the passive pool going into the active pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_PA',frac_carb_pa)
       !
       !Config Key   = FRAC_CARB_PS
       !Config Desc  = frac_carb_coefficients from passive pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.0
       !Config Help  = fraction of the passive pool going into the slow pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_PS',frac_carb_ps)
       !
       !Config Key   = ACTIVE_TO_PASS_CLAY_FRAC
       !Config Desc  = 
       !Config if    = OK_STOMATE 
       !Config Def   = 0.68  
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('ACTIVE_TO_PASS_CLAY_FRAC',active_to_pass_clay_frac)
       !
       !Config Key   = CARBON_TAU_IACTIVE
       !Config Desc  = residence times in carbon pools
       !Config if    = OK_STOMATE 
       !Config Def   = 0.149
       !Config Help  =
       !Config Units =  [days] 
       CALL getin_p('CARBON_TAU_IACTIVE',carbon_tau_iactive)
       !
       !Config Key   = CARBON_TAU_ISLOW
       !Config Desc  = residence times in carbon pools
       !Config if    = OK_STOMATE 
       !Config Def   = 5.48
       !Config Help  =
       !Config Units = [days]
       CALL getin_p('CARBON_TAU_ISLOW',carbon_tau_islow)
       !
       !Config Key   = CARBON_TAU_IPASSIVE
       !Config Desc  = residence times in carbon pools
       !Config if    = OK_STOMATE 
       !Config Def   = 241.
       !Config Help  = residence time in the passive pool
       !Config Units = [days] 
       CALL getin_p('CARBON_TAU_IPASSIVE',carbon_tau_ipassive)
       !
       !Config Key   = FLUX_TOT_COEFF
       !Config Desc  =
       !Config if    = OK_STOMATE 
       !Config Def   = 1.2, 1.4,.75
       !Config Help  =
       !Config Units = [days] 
       CALL getin_p('FLUX_TOT_COEFF',flux_tot_coeff)

       !-
       ! turnover parameters
       !-
       !
       !Config Key   = NEW_TURNOVER_TIME_REF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('NEW_TURNOVER_TIME_REF',new_turnover_time_ref)
       !
       !Config Key   = DT_TURNOVER_TIME 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('DT_TURNOVER_TIME',dt_turnover_time)
       !
       !Config Key   = LEAF_AGE_CRIT_TREF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('LEAF_AGE_CRIT_TREF',leaf_age_crit_tref)
       !
       !Config Key   = LEAF_AGE_CRIT_COEFF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.5, 0.75, 10. 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('LEAF_AGE_CRIT_COEFF',leaf_age_crit_coeff)

       !-
       ! vmax parameters
       !-
       !
       !Config Key   = VMAX_OFFSET 
       !Config Desc  = offset (minimum relative vcmax)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = offset (minimum vcmax/vmax_opt)
       !Config Units = [-]  
       CALL getin_p('VMAX_OFFSET',vmax_offset)
       !
       !Config Key   = LEAFAGE_FIRSTMAX
       !Config Desc  = leaf age at which vmax attains vcmax_opt (in fraction of critical leaf age)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.03 
       !Config Help  = relative leaf age at which vmax attains vcmax_opt
       !Config Units = [-] 
       CALL getin_p('LEAFAGE_FIRSTMAX',leafage_firstmax)
       !
       !Config Key   = LEAFAGE_LASTMAX 
       !Config Desc  = leaf age at which vmax falls below vcmax_opt (in fraction of critical leaf age) 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = relative leaf age at which vmax falls below vcmax_opt
       !Config Units = [-]  
       CALL getin_p('LEAFAGE_LASTMAX',leafage_lastmax)
       !
       !Config Key   = LEAFAGE_OLD 
       !Config Desc  = leaf age at which vmax attains its minimum (in fraction of critical leaf age)
       !Config If    = OK_STOMATE 
       !Config Def   = 1.
       !Config Help  = relative leaf age at which vmax attains its minimum
       !Config Units = [-]  
       CALL getin_p('LEAFAGE_OLD',leafage_old)

       !-
       ! season parameters
       !-
       !
       !Config Key   = GPPFRAC_DORMANCE 
       !Config Desc  = rapport maximal GPP/GGP_max pour dormance
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2 
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('GPPFRAC_DORMANCE',gppfrac_dormance)
       !
       !Config Key   = MIN_GPP_ALLOWED
       !Config Desc  = minimum gpp considered as not "lowgpp"
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3 
       !Config Help  = 
       !Config Units = [gC/m^2/year] 
       CALL getin_p('MIN_GPP_ALLOWED',min_gpp_allowed)
       !
       !Config Key   = TAU_CLIMATOLOGY
       !Config Desc  = tau for "climatologic variables 
       !Config If    = OK_STOMATE 
       !Config Def   = 20
       !Config Help  = 
       !Config Units = [days]
       CALL getin_p('TAU_CLIMATOLOGY',tau_climatology)
       !
       !Config Key   = HVC1 
       !Config Desc  = parameters for herbivore activity
       !Config If    = OK_STOMATE 
       !Config Def   = 0.019
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('HVC1',hvc1)
       !
       !Config Key   = HVC2 
       !Config Desc  = parameters for herbivore activity 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.38
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('HVC2',hvc2)
       !
       !Config Key   = LEAF_FRAC_HVC
       !Config Desc  = parameters for herbivore activity 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.33
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('LEAF_FRAC_HVC',leaf_frac_hvc)
       !
       !Config Key   = TLONG_REF_MAX
       !Config Desc  = maximum reference long term temperature 
       !Config If    = OK_STOMATE 
       !Config Def   = 303.1
       !Config Help  = 
       !Config Units = [K]  
       CALL getin_p('TLONG_REF_MAX',tlong_ref_max)
       !
       !Config Key   = TLONG_REF_MIN 
       !Config Desc  = minimum reference long term temperature 
       !Config If    = OK_STOMATE 
       !Config Def   = 253.1
       !Config Help  = 
       !Config Units = [K]  
       CALL getin_p('TLONG_REF_MIN',tlong_ref_min)
       !
       !Config Key   = NCD_MAX_YEAR
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 3. 
       !Config Help  = NCD : Number of Chilling Days
       !Config Units = [days]
       CALL getin_p('NCD_MAX_YEAR',ncd_max_year)
       !
       !Config Key   = GDD_THRESHOLD 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 5. 
       !Config Help  = GDD : Growing-Degree-Day
       !Config Units = [days] 
       CALL getin_p('GDD_THRESHOLD',gdd_threshold)
       !
       !Config Key   = GREEN_AGE_EVER 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 2. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GREEN_AGE_EVER',green_age_ever)
       !
       !Config Key   = GREEN_AGE_DEC
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('GREEN_AGE_DEC',green_age_dec)
   
!Chloe Methane model : 
   
       !-
       ! establish WETLAND CH4 methane parameters
       !

       !Config Key   = nvert
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 171 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NVERT',nvert)

       !Config Key   = ns
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 151 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NS',ns)

       !Config Key   = nday
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 24
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NDAY',nday)

       !Config Key   = h
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('H',h)

       !Config Key   = rk
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RK',rk)

       !Config Key   = diffair
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 7.2 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DIFFAIR',diffair)
       !Config Key   = pox
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('POX',pox)

       !Config Key   = dveg
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.001 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DVEG',dveg)

       !Config Key   = rkm
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 5.0
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RKM',rkm)

       !Config Key   = xvmax
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 20.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('XVMAX',xvmax)

       !Config Key   = oxq10
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 2.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('OXQ10',oxq10)

       !Config Key   = scmax
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 500. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('SCMAX',scmax)

       !Config Key   = sr0pl
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 600. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('SR0PL',sr0pl)

       !Config Key   = pwater_peat
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -3 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_PEAT',pwater_peat)

       !Config Key   = rpv
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RPV',rpv)

       !Config Key   = iother
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = -1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('IOTHER',iother)

  !Config Key   = rq10
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 3.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RQ10',rq10)

       !Config Key   = alpha_CH4
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = /0.009,0.004,0.021/
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_CH4',alpha_CH4)

!Chloe -- End Methane def parameters


       first_call = .FALSE.
       
    ENDIF
    
  END SUBROUTINE config_stomate_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_dgvm_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when the DGVM model is activated (ie : when ok_dgvm is set to true).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_dgvm_parameters   
    
    IMPLICIT NONE
    
    !! 0. Variables and parameters declaration

    !! 0.4 Local variables

    LOGICAL, SAVE ::  first_call = .TRUE.         !! To keep first call trace (true/false)

!_ ================================================================================================================================    

    IF(first_call) THEN
  
       !-
       ! establish parameters
       !-
       !
       !Config Key   = ESTAB_MAX_TREE
       !Config Desc  = Maximum tree establishment rate 
       !Config If    = OK_DGVM
       !Config Def   = 0.12 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ESTAB_MAX_TREE',estab_max_tree)
       !
       !Config Key   = ESTAB_MAX_GRASS
       !Config Desc  = Maximum grass establishment rate
       !Config If    = OK_DGVM
       !Config Def   = 0.12 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('ESTAB_MAX_GRASS',estab_max_grass)
       !
       !Config Key   = ESTABLISH_SCAL_FACT
       !Config Desc  = 
       !Config If    = OK_DGVM 
       !Config Def   = 5.
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('ESTABLISH_SCAL_FACT',establish_scal_fact)
       !
       !Config Key   = MAX_TREE_COVERAGE 
       !Config Desc  = 
       !Config If    = OK_DGVM 
       !Config Def   = 0.98
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('MAX_TREE_COVERAGE',max_tree_coverage)
       !
       !Config Key   = IND_0_ESTAB
       !Config Desc  = 
       !Config If    = OK_DGVM 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('IND_0_ESTAB',ind_0_estab)

       !-
       ! light parameters
       !-
       !
       !Config Key   = ANNUAL_INCREASE
       !Config Desc  = for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or to fpc of last time step (F)?
       !Config If    = OK_DGVM
       !Config Def   = y
       !Config Help  = 
       !Config Units = [FLAG]
       CALL getin_p('ANNUAL_INCREASE',annual_increase)
       !
       !Config Key   = MIN_COVER 
       !Config Desc  = For trees, minimum fraction of crown area occupied 
       !Config If    = OK_DGVM
       !Config Def   = 0.05 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('MIN_COVER',min_cover)

       !-
       ! pftinout parameters
       !
       !Config Key   = IND_0 
       !Config Desc  = initial density of individuals
       !Config If    = OK_DGVM
       !Config Def   = 0.02 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('IND_0',ind_0)
       !
       !Config Key   = MIN_AVAIL
       !Config Desc  = minimum availability
       !Config If    = OK_DGVM
       !Config Def   = 0.01
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('MIN_AVAIL',min_avail)
       !
       !Config Key   = RIP_TIME_MIN
       !Config Desc  = 
       !Config If    = OK_DGVM
       !Config Def   = 1.25 
       !Config Help  = 
       !Config Units = [year]  
       CALL getin_p('RIP_TIME_MIN',RIP_time_min)
       !
       !Config Key   = NPP_LONGTERM_INIT
       !Config Desc  = 
       !Config If    = OK_DGVM
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [gC/m^2/year]
       CALL getin_p('NPP_LONGTERM_INIT',npp_longterm_init)
       !
       !Config Key   = EVERYWHERE_INIT
       !Config Desc  = 
       !Config If    = OK_DGVM
       !Config Def   = 0.05 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('EVERYWHERE_INIT',everywhere_init)
       
       first_call = .FALSE.
       
    ENDIF
    
    
  END SUBROUTINE config_dgvm_parameters


END MODULE constantes
