! ================================================================================================================================
!  MODULE       : diffuco
!
!  CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module calculates the limiting coefficients, both aerodynamic
!! and hydrological, for the turbulent heat fluxes.
!!
!!\n DESCRIPTION: The aerodynamic resistance R_a is used to limit
!! the transport of fluxes from the surface layer of vegetation to the point in the atmosphere at which
!! interaction with the LMDZ atmospheric circulation model takes place. The aerodynamic resistance is
!! calculated either within the module r_aerod (if the surface drag coefficient is provided by the LMDZ, and 
!! if the flag 'ldq_cdrag_from_gcm' is set to TRUE) or r_aero (if the surface drag coefficient must be calculated).\n
!!
!! Within ORCHIDEE, evapotranspiration is a function of the Evaporation Potential, but is modulated by a
!! series of resistances (canopy and aerodynamic) of the surface layer, here represented by beta.\n
!!
!! DESCRIPTION	:
!! \latexonly 
!!     \input{diffuco_intro.tex}
!! \endlatexonly
!! \n
!!
!! This module calculates the beta for several different scenarios: 
!! - diffuco_snow calculates the beta coefficient for sublimation by snow, 
!! - diffuco_inter calculates the beta coefficient for interception loss by each type of vegetation, 
!! - diffuco_bare calculates the beta coefficient for bare soil, 
!! - diffuco_trans or diffuco_trans_co2 both calculate the beta coefficient for transpiration for each type
!!   of vegetation. Those routines differ by the formulation used to calculate the canopy resistance (Jarvis in 
!!   diffuco_trans, Farqhar in diffuco_trans_co2)
!! - diffuco_inca calculates the beta coefficient for emissions of biogenic compounds \n
!!
!! Finally, the module diffuco_comb computes the combined $\alpha$ and $\beta$ coefficients for use 
!! elsewhere in the module. \n

!! RECENT CHANGE(S): Nathalie le 28 mars 2006 - sur proposition de Fred Hourdin, ajout
!! d'un potentiometre pour regler la resistance de la vegetation (rveg is now in pft_parameters)
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/diffuco.f90 $
!! $Date: 2013-06-20 15:26:25 +0200 (Thu, 20 Jun 2013) $
!! $Revision: 1322 $
!! \n
!_ ================================================================================================================================

MODULE diffuco

  ! modules used :
  USE constantes
  USE qsat_moisture
  USE sechiba_io
  USE ioipsl
  USE pft_parameters
  USE grid
  USE parallel
  USE slowproc
!  USE WRITE_FIELD_p

  IMPLICIT NONE

  ! public routines :
  ! diffuco_main only
  PRIVATE
  PUBLIC :: diffuco_main,diffuco_clear

  !
  ! variables used inside diffuco module : declaration and initialisation
  !
  LOGICAL, SAVE                                      :: l_first_diffuco = .TRUE.  !! Initialisation has to be done one time
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: leaf_ci                   !! intercellular CO2 concentration (ppm)
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: rstruct                   !! architectural resistance (s m^{-1})
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: raero                     !! Aerodynamic resistance (s m^{-1})
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: qsatt 
!Isa
  LOGICAL, SAVE      ::  ok_snowfact ! facteur divisant le z0 en cas de neige _ pour limiter ma subli
                    !! Surface saturated humidity (kg kg^{-1})
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: wind                      !! Wind module (m s^{-1})

  ! variables used inside diffuco_inca module 
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)     :: ok_siesta        !! Flag for controlling post-pulse period (true/false)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)     :: allow_pulse      !! Flag for controlling pulse period (true/false)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pulse            !! Pulse fonction 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pulseday         !! Counter for pulse period
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: siestaday        !! Counter for post-pulse period
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pulselim         !! Pulse period length
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: siestalim        !! Post-pulse period length
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: area2            !! Grid cell area (m^2)
  REAL(r_std), SAVE                            :: nbre_precip 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: flx_co2_bbg_year !! CO2 emissions from bbg, 
                                                                   !! read in a file (kgC.m^{-2}.year^{-1})
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: N_qt_WRICE_year  !! N fertilizers on wetland rice,
                                                                   !! read in a file expressed in kgN/year/grid cell
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: N_qt_OTHER_year  !! N fertilizers on other crops and grasses,
                                                                   !! read in a file expressed in kgN/year/grid cell

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_main
!!
!>\BRIEF	 The root subroutine for the module, which calls all other required
!! subroutines.
!! 
!! DESCRIPTION   : 

!! This is the root subroutine for the module. Following
!! initialisation (if required) and preparation of the restart file, the module first of all calculates
!! the surface drag coefficient (via a call to diffuco_aero), using available parameters to determine
!! the stability of air in the surface layer by calculating the Richardson Nubmber. If a value for the 
!! surface drag coefficient is passed down from the atmospheric model and and if the flag 'ldq_cdrag_from_gcm' 
!! is set to TRUE, then the subroutine diffuco_aerod is called instead. This calculates the aerodynamic coefficient. \n
!!
!! Following this, an estimation of the saturated humidity at the surface is made (via a call
!! to qsatcalc in the module qsat_moisture). Following this the beta coefficients for sublimation (via 
!! diffuco_snow), interception (diffuco_inter), bare soil (diffuco_bare), and transpiration (via 
!! diffuco_trans_co2 if co2 is considered, diffuco_trans otherwise) are calculated in sequence. Finally 
!! the alpha and beta coefficients are combined (diffuco_comb). \n
!!
!! The surface drag coefficient is calculated for use within the module enerbil. It is required to to
!! calculate the aerodynamic coefficient for every flux. \n
!!
!! The various beta coefficients are used within the module enerbil for modifying the rate of evaporation, 
!! as appropriate for the surface. As explained in Chapter 2 of Guimberteau (2010), that module (enerbil) 
!! calculates the rate of evaporation essentially according to the expression $E = /beta E_{pot}$, where
!! E is the total evaporation and $E_{pot}$ the evaporation potential. If $\beta = 1$, there would be
!! essentially no resistance to evaporation, whereas should $\beta = 0$, there would be no evaporation and
!! the surface layer would be subject to some very stong hydrological stress. \n
!!
!! ATTENTION ::valpha [DISPOSABLE] is not used anymore ! Its value is always 1.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): humrel, q_cdrag, vbeta, valpha, vbeta1, vbeta4,
!! vbetaco2, vbeta2, vbeta3, rveget, cimean   
!!
!! REFERENCE(S) :				        
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273.
!! - de Rosnay, P, 1999. Représentation des interactions sol-plante-atmosphère dans le modèle de circulation générale
!! du LMD, 1999. PhD Thesis, Université Paris 6, available (25/01/12): 
!! http://www.ecmwf.int/staff/patricia_de_rosnay/publications.html#8
!! - Ducharne, A, 1997. Le cycle de l'eau: modélisation de l'hydrologie continentale, étude de ses interactions avec 
!! le climat, PhD Thesis, Université Paris 6
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available (25/01/12):
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Lathière, J, 2005. Evolution des émissions de composés organiques et azotés par la biosphère continentale dans le 
!! modèle LMDz-INCA-ORCHIDEE, Université Paris 6
!!
!! FLOWCHART	:
!! \latexonly 
!!     \includegraphics[scale=0.5]{diffuco_main_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================


! Main routine for *diffuco* module.
! - called only one time for initialisation
! - called every time step for calculations (also the first time step)
! - called one more time at last time step for writing _restart_ file
!
! The following processes are calculated:
! - call diffuco_aero for aerodynamic transfer coeficient
! - call diffuco_snow for partial beta coefficient: sublimation
! - call diffuco_inter for partial beta coefficient: interception for each type of vegetation
! - call diffuco_bare for partial beta coefficient: bare soil
! - call diffuco_trans for partial beta coefficient: transpiration for each type of vegetation, using Jarvis formula
! - call diffuco_trans_co2 for partial beta coefficient: transpiration for each type of vegetation, using Farqhar's formula
! - call diffuco_comb for alpha and beta coefficient
! - call diffuco_inca for alpha and beta coefficients for biogenic emissions
  SUBROUTINE diffuco_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, indexveg, indexlai, u, v, &
! Ajout Nathalie - Juin 2006 - passage q2m/t2m pour calcul Rveget
!     & zlev, z0, roughheight, temp_sol, temp_air, rau, q_cdrag, qsurf, qair, pb, &
     & zlev, z0, roughheight, temp_sol, temp_air, rau, q_cdrag, qsurf, qair, q2m, t2m, pb, &
     & rsol, evap_bare_lim, evapot, evapot_corr, snow, flood_frac, flood_res, frac_nobio, snow_nobio, totfrac_nobio, &
     & swnet, swdown, sinang, ccanopy, humrel, veget, veget_max, lai, qsintveg, qsintmax, assim_param, &
     & vbeta , valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, vbetaco2, rveget, rstruct, cimean, &
     & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, &
     & rest_id, hist_id, hist2_id)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! _Restart_ file identifier (-)
    INTEGER(i_std),INTENT (in)                         :: hist_id          !! _History_ file identifier (-)
    INTEGER(i_std),INTENT (in)                         :: hist2_id         !! _History_ file 2 identifier (-)
    REAL(r_std), INTENT (in)                           :: dtradia          !! Time step (s)
    LOGICAL, INTENT(in)                                :: ldrestart_read   !! Logical for restart file to read (-)
    LOGICAL, INTENT(in)                                :: ldrestart_write  !! Logical for restart file to write (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: index          !! Indeces of the points on the map (-)
    INTEGER(i_std),DIMENSION (kjpindex*(nlai+1)), INTENT (in) :: indexlai  !! Indeces of the points on the 3D map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in) :: indexveg       !! Indeces of the points on the 3D map (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u                !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v                !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: zlev             !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: z0               !! Surface roughness Length (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: roughheight      !! Effective height for roughness (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol         !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_air         !! Lowest level temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rau              !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qsurf            !! Near surface air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qair             !! Lowest level air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q2m              !! 2m air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: t2m              !! 2m air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow             !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac       !! Fraction of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_res        !! Reservoir in floodplains (estimation to avoid over-evaporation)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb               !! Surface level pressure (Pa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rsol             !! Bare soil evaporation resistance (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evap_bare_lim    !! Limit to the bare soil evaporation when the 
                                                                           !! 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot           !! Soil Potential Evaporation (mm day^{-1}) 
                                                                           !! NdN This variable does not seem to be used at 
                                                                           !! all in diffuco
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_corr      !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice,lakes,cities,... (-)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Snow on ice,lakes,cities,... (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio    !! Total fraction of ice+lakes+cities+... (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet            !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swdown           !! Down-welling surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: sinang           !! Sinus of Solar Angle (as computed in read_dim2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ccanopy          !! CO2 concentration inside the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI->infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: lai              !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg         !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax         !! Maximum water on vegetation for interception 
                                                                           !! (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2), INTENT (in) :: assim_param !! min+max+opt temps, vcmax, vjmax
                                                                           !! for photosynthesis (K ??)
    REAL(r_std),DIMENSION (kjpindex,2),   INTENT (in)  :: lalo               !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,8),INTENT (in)  :: neighbours         !! Vector of neighbours for each 
                                                                             !! grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)     :: resolution         !! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ptnlev1            !! 1st level of soil temperature (Kelvin)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain        !! Rain precipitation expressed in mm/tstep
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT (in)  :: frac_age !! Age efficiency for isoprene emissions (from STOMATE)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta            !! Total beta coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: valpha           !! Total alpha coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta1           !! Beta for sublimation (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta4           !! Beta for bare soil evaporation (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta5           !! Beta for floodplains
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbetaco2         !! STOMATE: Beta for CO2 assimilation (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta2           !! Beta for interception loss (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3           !! Beta for transpiration (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3pot        !! Beta for potential transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rveget           !! Stomatal resistance for the whole canopy (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rstruct          !! Structural resistance for the vegetation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: cimean           !! STOMATE: mean intercellular ci (\mumole m^{-2} s^{-1}) 

    !! 0.3 Modified variables
 
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (inout) :: humrel        !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: q_cdrag       !! Product of drag coefficient and wind speed (m s^{-1})

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex,nvm)               :: vbeta23          !! Beta for fraction of wetted foliage that will
                                                                           !! transpire once intercepted water has evaporated (-)
    INTEGER(i_std)                                    :: ilai
    CHARACTER(LEN=4)                                  :: laistring
    CHARACTER(LEN=80)                                 :: var_name                  !! To store variables names for I/O

    ! Biogenic emissions
    REAL(r_std),DIMENSION(kjpindex)          :: PAR                !! Photosynthetic active radiation, half of swdown
                                                                   !! @tex ($\mu mol photons. m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: PARsun             !! PAR received by sun leaves
                                                                   !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: PARsh              !! PAR received by shaded leaves 
                                                                   !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: laisun             !! Leaf area index of Sun leaves (m^2.m^{-2})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: laish              !! Leaf area index of Shaded leaves (m^2.m^{-2}) 
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_iso            !! Biogenic isoprene emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_mono           !! Biogenic monoterpene emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_ORVOC          !! Biogenic ORVOC emission - (kgC.m^{-2}.s^{-1}) 
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_MBO            !! Biogenic MBO emission -
                                                                   !! MethylButanOl (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_methanol       !! Biogenic methanol emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_acetone        !! Biogenic acetone emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_acetal         !! Biogenic Acetaldehyde emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_formal         !! Biogenic Formaldehyde emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_acetic         !! Biogenic Acetic Acid emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_formic         !! Biogenic Formic Acid emission (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_no_soil        !! Biogenic NO emission by soil (kgC.m^{-2}.s^{-1})               
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_no             !! Biogenic net NO emission (kgC.m^{-2}.s^{-1})                 
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: CRF                !! Canopy reduction factor for net NO flux calculation    
    REAL(r_std),DIMENSION(kjpindex,nvm)      :: flx_fertil_no      !! Biogenic NO emission due to N-fertilisation 
                                                                   !! (kgC.m^{-2}.s^{-1})
    REAL(r_std),DIMENSION(kjpindex)          :: Fdf                !! Diffuse Fraction of the radiation (0-1, unitless)
    REAL(r_std),DIMENSION(kjpindex,nlai+1)   :: PARsuntab          !! PAR received by sun leaves 
                                                                   !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex,nlai+1)   :: PARshtab           !! PAR received by shaded leaves 
                                                                   !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex)          :: PARdf              !! Diffuse PAR
                                                                   !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex)          :: PARdr              !! Direct PAR 
                                                                   !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex)          :: Trans              !! Atmospheric Transmissivity (unitless)

    REAL(r_std),DIMENSION(kjpindex)          :: julian_diff_2d     !!
    REAL(r_std),DIMENSION(kjpindex)          :: year_length_2d     !!

!_ ================================================================================================================================
    
  !! 1. Perform initialisation, if required
    
    IF (l_first_diffuco) THEN

        !Config Key   = CDRAG_FROM_GCM
        !Config Desc  = Keep cdrag coefficient from gcm.
        !Config If    = OK_SECHIBA
        !Config Def   = y
        !Config Help  = Set to .TRUE. if you want q_cdrag coming from GCM (if q_cdrag on initialization is non zero).
        !Config         Keep cdrag coefficient from gcm for latent and sensible heat fluxes.
        !Config Units = [FLAG]
        IF ( ABS(MAXVAL(q_cdrag)) .LE. EPSILON(q_cdrag)) THEN
           ldq_cdrag_from_gcm = .FALSE.
        ELSE
           ldq_cdrag_from_gcm = .TRUE.
        ENDIF
        
        !?? q_cdrag is always 0 on initialization ??
        CALL getin_p('CDRAG_from_GCM', ldq_cdrag_from_gcm)
        
        WRITE(numout,*) "ldq_cdrag_from_gcm = ",ldq_cdrag_from_gcm

        IF (long_print) WRITE (numout,*) ' call diffuco_init '

        ! If cdrag is 
        CALL diffuco_init(kjit, ldrestart_read, kjpindex, index, rest_id, q_cdrag, rstruct)

        WRITE(numout,*) "control%ok_inca:",control%ok_inca
        IF ( control%ok_inca ) CALL diffuco_inca_init(kjpindex, lalo, neighbours, resolution, dtradia)

        RETURN

    ENDIF
    
  !! 2. Prepare the restart file for the next simulation

    IF (ldrestart_write) THEN

        IF (long_print) WRITE (numout,*) ' we have to complete restart file with DIFFUCO variables '

        var_name= 'rstruct'
	CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, rstruct, 'scatter',  nbp_glo, index_g)
  
        var_name= 'raero'
	CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, raero, 'scatter',  nbp_glo, index_g)

        var_name= 'qsatt'
	CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, qsatt, 'scatter',  nbp_glo, index_g)
  
        ! the following variable is written only if CO2 was calculated
        IF ( control%ok_co2 ) THEN

          DO ilai = 1, nlai
  
            ! variable name is somewhat complicated as ioipsl does not allow 3d variables for the moment...
            write(laistring,'(i4)') ilai
            laistring=ADJUSTL(laistring)
            var_name='leaf_ci_'//laistring(1:LEN_TRIM(laistring))
  
            CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, leaf_ci(:,:,ilai), 'scatter',  nbp_glo, index_g)
  
          ENDDO

        ENDIF
  
        IF (.NOT.ldq_cdrag_from_gcm) THEN
           var_name= 'cdrag'
           CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, q_cdrag, 'scatter',  nbp_glo, index_g)
        ENDIF
  
      RETURN
  
    END IF
    
    wind(:) = SQRT (u(:)*u(:) + v(:)*v(:))

    
  !! 3. Calculate the different coefficients

    IF (.NOT.ldq_cdrag_from_gcm) THEN
        ! Case 3a)
        !Isa: ajout snow
       CALL diffuco_aero (kjpindex, kjit, u, v, zlev, z0, roughheight, temp_sol, temp_air, &
            &             qsurf, qair, q_cdrag, snow)
    ENDIF

    ! Case 3b)
    CALL diffuco_raerod (kjpindex, u, v, q_cdrag, raero)

  !! 4. Make an estimation of the saturated humidity at the surface

    CALL qsatcalc (kjpindex, temp_sol, pb, qsatt)

  !! 5. Calculate the beta coefficient for sublimation
  
    CALL diffuco_snow (kjpindex, dtradia, qair, qsatt, rau, u, v, q_cdrag, &
         & snow, frac_nobio, totfrac_nobio, snow_nobio, vbeta1)


    CALL diffuco_flood (kjpindex, dtradia, qair, qsatt, rau, u, v, q_cdrag, evapot, evapot_corr, &
         & flood_frac, flood_res, vbeta5)

  !! 6. Calculate the beta coefficient for interception

    ! Correction Nathalie - Juin 2006 - introduction d'un terme vbeta23
    !CALL diffuco_inter (kjpindex, dtradia, qair, qsatt, rau, u, v, q_cdrag, veget, &
    !   & qsintveg, qsintmax, rstruct, vbeta2) 
    CALL diffuco_inter (kjpindex, dtradia, qair, qsatt, rau, u, v, q_cdrag, veget, &
       & qsintveg, qsintmax, rstruct, vbeta2, vbeta23) 


  !! 8. Calculate the beta coefficient for transpiration

    IF ( control%ok_co2 ) THEN

      ! Ajout Nathalie - Juin 2006 - passage q2m/t2m pour calcul Rveget
      ! Correction Nathalie - Juin 2006 - introduction d'un terme vbeta23
      !CALL diffuco_trans_co2 (kjpindex, dtradia, swdown, temp_air, pb, qair, rau, u, v, q_cdrag, humrel, &
      !                        assim_param, ccanopy, &
      !                        veget, veget_max, lai, qsintveg, qsintmax, vbeta3, rveget, rstruct, cimean, vbetaco2)
 
      ! case 8a)
      CALL diffuco_trans_co2 (kjpindex, dtradia, swdown, temp_air, pb, qair, q2m, t2m, rau, u, v, q_cdrag, humrel, &
                              assim_param, ccanopy, &
                              veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, &
                              rveget, rstruct, cimean, vbetaco2, vbeta23)

    ELSE

      ! Correction Nathalie - Juin 2006 - introduction d'un terme vbeta23
      !CALL diffuco_trans (kjpindex, dtradia, swnet, temp_air, pb, qair, rau, u, v, q_cdrag, humrel, &
      !                    veget, veget_max, lai, qsintveg, qsintmax, vbeta3, rveget, rstruct, cimean, vbetaco2)
      
      ! case 8b) 
      CALL diffuco_trans (kjpindex, dtradia, swnet, temp_air, pb, qair, rau, u, v, q_cdrag, humrel, &
           & veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, rveget, rstruct, cimean, &
           & vbetaco2, vbeta23)

    ENDIF


    !
    !biogenic emissions
    !
    IF ( control%ok_inca ) THEN
       CALL diffuco_inca (kjpindex, dtradia, swdown, sinang, temp_air, temp_sol, ptnlev1, precip_rain, humrel, &
                        veget_max, lai, frac_age, &
                        lalo, &
                        PAR, PARsun, PARsh, laisun, laish, &
                        flx_iso, flx_mono, flx_ORVOC, flx_MBO, flx_methanol, flx_acetone, flx_acetal, &
                        flx_formal, flx_acetic, flx_formic, &
                        flx_no_soil, flx_no,CRF, flx_fertil_no, Fdf, PARsuntab, PARshtab, PARdf, PARdr, Trans)
    ENDIF
    !
    ! combination of coefficient : alpha and beta coefficient
    ! beta coefficient for bare soil
    !

    CALL diffuco_bare (kjpindex, dtradia, u, v, q_cdrag, rsol, evap_bare_lim, evapot, humrel, &
         veget, veget_max, vbeta2, vbeta3, vbeta4) 

  !! 9. Combine the alpha and beta coefficients

    ! Ajout qsintmax dans les arguments de la routine.... Nathalie / le 13-03-2006
    CALL diffuco_comb (kjpindex, dtradia, humrel, rau, u, v, q_cdrag, pb, qair, temp_sol, temp_air, snow, &
       & veget, lai, vbeta1, vbeta2, vbeta3, vbeta4, valpha, vbeta, qsintmax)    

    julian_diff_2d(:) = julian_diff
    year_length_2d(:) = year_length

    IF ( .NOT. almaoutput ) THEN
       CALL histwrite(hist_id, 'raero', kjit, raero, kjpindex, index)
       ! Ajouts Nathalie - novembre 2006
       CALL histwrite(hist_id, 'cdrag', kjit, q_cdrag, kjpindex, index)
       CALL histwrite(hist_id, 'Wind', kjit, wind, kjpindex, index)
       ! Fin ajouts Nathalie
       CALL histwrite(hist_id, 'qsatt', kjit, qsatt, kjpindex, index)

       IF ( control%ok_inca ) THEN
          CALL histwrite(hist_id, 'PAR', kjit, PAR, kjpindex, index)
          IF ( control%ok_radcanopy ) THEN
             CALL histwrite(hist_id, 'PARsun', kjit, PARsun, kjpindex*nvm, indexveg)
             CALL histwrite(hist_id, 'PARsh', kjit, PARsh, kjpindex*nvm, indexveg)
             CALL histwrite(hist_id, 'laisun', kjit, laisun, kjpindex*nvm, indexveg)
             CALL histwrite(hist_id, 'laish', kjit, laish, kjpindex*nvm, indexveg)
             CALL histwrite(hist_id, 'Fdf', kjit, Fdf, kjpindex, index)
             CALL histwrite(hist_id, 'PARsuntab', kjit, PARsuntab, kjpindex*(nlai+1), indexlai)
             CALL histwrite(hist_id, 'PARshtab', kjit, PARshtab, kjpindex*(nlai+1), indexlai)
             CALL histwrite(hist_id, 'Sinang', kjit, Sinang, kjpindex, index)
             CALL histwrite(hist_id, 'PARdf', kjit, PARdf, kjpindex, index)
             CALL histwrite(hist_id, 'PARdr', kjit, PARdr, kjpindex, index)
             CALL histwrite(hist_id, 'Trans', kjit, Trans, kjpindex, index)
             CALL histwrite(hist_id, 'Day', kjit, julian_diff_2d,kjpindex,index)
             CALL histwrite(hist_id, 'Year_length', kjit, year_length_2d,kjpindex,index)

          END IF
          CALL histwrite(hist_id, 'flx_fertil_no', kjit, flx_fertil_no, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'CRF', kjit, CRF, kjpindex*nvm, indexveg)
          IF ( control%ok_bbgfertil_Nox ) THEN
             CALL histwrite(hist_id, 'flx_co2_bbg_year', 1, flx_co2_bbg_year, kjpindex, index)
          ENDIF
          IF ( control%ok_cropsfertil_Nox ) THEN
             CALL histwrite(hist_id, 'N_qt_WRICE_year', 1, N_qt_WRICE_year, kjpindex, index)
             CALL histwrite(hist_id, 'N_qt_OTHER_year', 1, N_qt_OTHER_year, kjpindex, index)
          ENDIF
          CALL histwrite(hist_id, 'ptnlev1', kjit, ptnlev1, kjpindex, index)
          CALL histwrite(hist_id, 'flx_iso', kjit, flx_iso, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_mono', kjit, flx_mono, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_ORVOC', kjit, flx_ORVOC, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_MBO', kjit, flx_MBO, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_methanol', kjit, flx_methanol, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_acetone', kjit, flx_acetone, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_acetal', kjit, flx_acetal, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_formal', kjit, flx_formal, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_acetic', kjit, flx_acetic, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_formic', kjit, flx_formic, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_no_soil', kjit, flx_no_soil, kjpindex*nvm, indexveg)
          CALL histwrite(hist_id, 'flx_no', kjit, flx_no, kjpindex*nvm, indexveg)
       END IF

       IF ( hist2_id > 0 ) THEN
          CALL histwrite(hist2_id, 'raero', kjit, raero, kjpindex, index)
          CALL histwrite(hist2_id, 'cdrag', kjit, q_cdrag, kjpindex, index)
          CALL histwrite(hist2_id, 'Wind', kjit, wind, kjpindex, index)
          CALL histwrite(hist2_id, 'qsatt', kjit, qsatt, kjpindex, index)

          IF ( control%ok_inca ) THEN
             CALL histwrite(hist2_id, 'PAR', kjit, PAR, kjpindex, index)
             IF ( control%ok_radcanopy ) THEN
                CALL histwrite(hist2_id, 'PARsun', kjit, PARsun, kjpindex*nvm, indexveg)
                CALL histwrite(hist2_id, 'PARsh', kjit, PARsh, kjpindex*nvm, indexveg)
                CALL histwrite(hist2_id, 'laisun', kjit, laisun, kjpindex*nvm, indexveg)
                CALL histwrite(hist2_id, 'laish', kjit, laish, kjpindex*nvm, indexveg)
             ENDIF
             CALL histwrite(hist2_id, 'flx_fertil_no', kjit, flx_fertil_no, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'CRF', kjit, CRF, kjpindex*nvm, indexveg)
             IF ( control%ok_bbgfertil_Nox ) THEN
                CALL histwrite(hist2_id, 'flx_co2_bbg_year', 1, flx_co2_bbg_year, kjpindex, index)
             ENDIF
             IF ( control%ok_cropsfertil_Nox ) THEN
                CALL histwrite(hist2_id, 'N_qt_WRICE_year', 1, N_qt_WRICE_year, kjpindex, index)
                CALL histwrite(hist2_id, 'N_qt_OTHER_year', 1, N_qt_OTHER_year, kjpindex, index)
             ENDIF
             CALL histwrite(hist2_id, 'ptnlev1', kjit, ptnlev1, kjpindex, index)
             CALL histwrite(hist2_id, 'flx_iso', kjit, flx_iso, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_mono', kjit, flx_mono, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_ORVOC', kjit, flx_ORVOC, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_MBO', kjit, flx_MBO, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_methanol', kjit, flx_methanol, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_acetone', kjit, flx_acetone, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_acetal', kjit, flx_acetal, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_formal', kjit, flx_formal, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_acetic', kjit, flx_acetic, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_formic', kjit, flx_formic, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_no_soil', kjit, flx_no_soil, kjpindex*nvm, indexveg)
             CALL histwrite(hist2_id, 'flx_no', kjit, flx_no, kjpindex*nvm, indexveg)
          ENDIF
       ENDIF
    ELSE

    ENDIF
   
    IF (long_print) WRITE (numout,*) ' diffuco_main done '

  END SUBROUTINE diffuco_main


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_init
!!
!>\BRIEF					Dynamic allocation algorithm for local arrays 
!!
!! DESCRIPTION				        : Housekeeping module to allocate local arrays
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S)	                : q_cdrag
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_init(kjit, ldrestart_read, kjpindex, index, rest_id, q_cdrag, rstruct)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                       :: kjit               !! Time step number  (-)
    LOGICAL,INTENT (in)                               :: ldrestart_read     !! Logical for restart file to read (-)
    INTEGER(i_std), INTENT (in)                       :: kjpindex           !! Domain size (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: index              !! Indeces of the points on the map (-)
    INTEGER(i_std), INTENT (in)                       :: rest_id            !! _Restart_ file identifier (-)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rstruct           !! STOMATE: architectural resistance

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)  :: q_cdrag            !! Product of Surface drag and wind speed (m s^{-1})

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: ier, jv
    INTEGER(i_std)                                    :: ilai
    CHARACTER(LEN=4)                                  :: laistring
    REAL(r_std),DIMENSION (kjpindex)                  :: temp
    CHARACTER(LEN=80)                                 :: var_name            !! To store variables names for I/O
!_ ================================================================================================================================
    
  !! 1. Initialisation
    
    IF (l_first_diffuco) THEN 
        l_first_diffuco=.FALSE.
    ELSE 
        WRITE (numout,*) ' l_first_diffuco false . we stop '
        STOP 'diffuco_init'
    ENDIF

    ! allocate only if CO2 is calculated
    IF ( control%ok_co2 ) THEN

      ALLOCATE (leaf_ci(kjpindex,nvm,nlai),stat=ier)
      IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in leaf_ci allocation. We stop. We need kjpindex*nvm*nlai words = ',&
            kjpindex*nvm*nlai
          STOP 'diffuco_init'
      END IF

    ENDIF

    ALLOCATE (raero(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in raero allocation. We stop. We need kjpindex x nvm words = ', kjpindex
        STOP 'diffuco_init'
    END IF

    ALLOCATE (qsatt(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in qsatt allocation. We stop. We need kjpindex x nvm words = ', kjpindex
        STOP 'diffuco_init'
    END IF

    ALLOCATE (wind(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in wind allocation. We stop. We need kjpindex x nvm words = ', kjpindex
        STOP 'diffuco_init'
    END IF

    IF (ldrestart_read) THEN

        IF (long_print) WRITE (numout,*) ' we have to read a restart file for DIFFUCO variables'

        var_name='rstruct'
        CALL ioconf_setatt('UNITS', 's/m')
        CALL ioconf_setatt('LONG_NAME','Structural resistance')
        CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., rstruct, "gather", nbp_glo, index_g)
        
        IF ( MINVAL(rstruct) .EQ. MAXVAL(rstruct) .AND.  MAXVAL(rstruct) .EQ. val_exp ) THEN
        
           DO jv = 1, nvm
              rstruct(:,jv) = rstruct_const(jv)
           ENDDO

        ENDIF

        var_name='raero' ;
        CALL ioconf_setatt('UNITS', 's/m')
        CALL ioconf_setatt('LONG_NAME','Aerodynamic resistance')
        IF ( ok_var(var_name) ) THEN
           CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
           IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
              raero(:) = temp(:)
           ENDIF
        ENDIF
        
        var_name='qsatt' ;
        CALL ioconf_setatt('UNITS', 'g/g')
        CALL ioconf_setatt('LONG_NAME','Surface saturated humidity')
        IF ( ok_var(var_name) ) THEN
           CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
           IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
              qsatt(:) = temp(:)
           ENDIF
        ENDIF

        ! the following variable is read only if CO2 is calculated
        IF ( control%ok_co2 ) THEN

           CALL ioconf_setatt('UNITS', 'ppm')
           CALL ioconf_setatt('LONG_NAME','Leaf CO2')

           DO ilai = 1, nlai

              ! variable name is somewhat complicated as ioipsl does not allow 3d variables for the moment...
              write(laistring,'(i4)') ilai
              laistring=ADJUSTL(laistring)
              var_name='leaf_ci_'//laistring(1:LEN_TRIM(laistring))
              
              CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE.,leaf_ci(:,:,ilai), "gather", nbp_glo, index_g)
              
           ENDDO
           
           !Config Key   = DIFFUCO_LEAFCI
           !Config Desc  = Initial leaf CO2 level if not found in restart
           !Config If    = OK_SECHIBA
           !Config Def   = 233.
           !Config Help  = The initial value of leaf_ci if its value is not found
           !Config         in the restart file. This should only be used if the model is
           !Config         started without a restart file.
           !Config Units = [ppm]
           CALL setvar_p (leaf_ci, val_exp,'DIFFUCO_LEAFCI', 233._r_std)
        ENDIF
        
        IF (.NOT.ldq_cdrag_from_gcm) THEN
           var_name= 'cdrag'
           CALL ioconf_setatt('LONG_NAME','Drag coefficient for LE and SH')
           CALL ioconf_setatt('UNITS', '-')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 q_cdrag(:) = temp(:)
              ENDIF
           ENDIF
           
        ENDIF
        
    ENDIF
    !
    !Config Key  = RVEG_PFT
    !Config Desc = Artificial parameter to increase or decrease canopy resistance.
    !Config Def  = 1.
    !Config Help = This parameter is set by PFT.
    rveg_pft(:) = 1.
    CALL getin_p('RVEG_PFT', rveg_pft)

    WRITE(numout,*) 'DANS DIFFUCO_INIT , RVEG_PFT=',rveg_pft

    !Isa
    if (control%ok_freeze) then
        ! on active par défaut les options d'Isabelle
        write(*,*) 'diffuco 761: ok_freeze active, on active ok_snowfact par defaut'
        ok_snowfact = .TRUE.
    else !if (ok_freeze) then
        write(*,*) 'diffuco 764: ok_freeze desactive, on lit ok_snowfact'
        ok_snowfact = .FALSE.
        CALL getin_p('OK_SNOWFACT', ok_snowfact)
    endif !if (ok_freeze) then
! end isa

    IF (long_print) WRITE (numout,*) ' diffuco_init done '

  END SUBROUTINE diffuco_init


!! ================================================================================================================================
!! SUBROUTINE   : diffuco_inca_init
!!
!>\BRIEF         This subroutine initializes diffuco_inca-related state variables
!!
!! DESCRIPTION  : Some of the variables and flags used in diffuco_inca are allocated and initialised here.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

  SUBROUTINE diffuco_inca_init(kjpindex, lalo, neighbours, resolution, dtradia)
    
    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (unitless) 
    REAL(r_std), DIMENSION(kjpindex,2), INTENT (in)    :: lalo             !! Geographical coordinates
    INTEGER(i_std), DIMENSION(kjpindex,8), INTENT (in) :: neighbours       !! Vector of neighbours for each 
                                                                           !! grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)     :: resolution       !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)                            :: dtradia          !! Time step in seconds

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    LOGICAL                             :: allow_weathergen, density
    CHARACTER(LEN=80)                   :: filename, filename2, fieldname
    INTEGER(i_std)                      :: iml, jml, lml, tml, force_id
    INTEGER(i_std)                      :: ier

!_ ================================================================================================================================

    ALLOCATE (pulse(kjpindex),stat=ier)
    IF (ier /= 0) THEN
       WRITE (numout,*) ' error in pulse allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'diffuco_inca_init'
    END IF
    pulse(:) = un

    ! If we acount for NOx pulse emissions
    IF (control%ok_pulse_NOx) THEN

       ALLOCATE (ok_siesta(kjpindex),stat=ier)
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in ok_siesta allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       ok_siesta(:) = .FALSE.

       ALLOCATE (allow_pulse(kjpindex),stat=ier)
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in allow_pulse allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       allow_pulse(:) = .FALSE.

       ALLOCATE (pulseday(kjpindex),stat=ier)
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in pulseday allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       pulseday(:) = zero

       ALLOCATE (siestaday(kjpindex),stat=ier)
       IF (ier /=0 ) THEN
          WRITE (numout,*) ' error in siestaday allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       siestaday(:) = zero

       ALLOCATE (pulselim(kjpindex),stat=ier)
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in pulselim allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       pulselim(:) = zero

       ALLOCATE (siestalim(kjpindex),stat=ier)
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in siestalim allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       siestalim(:) = zero

    END IF ! (control%ok_pulse_NOx) 

    ! If we acount for NOx emissions by N-fertilizers
    IF (control%ok_cropsfertil_NOx) THEN

       ALLOCATE (area2(kjpindex),stat=ier)
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in area2 allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       area2(:) = resolution(:,1)*resolution(:,2) 

       ALLOCATE (N_qt_WRICE_year(kjpindex),stat=ier)  !! N fertilizers on wetland rice, read in file 
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in N_qt_WRICE_year allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       N_qt_WRICE_year(:) = zero
    
       ALLOCATE (N_qt_OTHER_year(kjpindex),stat=ier)  !! N fertilizers on other crops and grasses, read in file 
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in N_qt_OTHER_year allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       N_qt_OTHER_year(:) = zero

       WRITE (numout,*) ' *********************** Interpolating N fertilizers files for NOx emissions... '
       filename = 'orchidee_fertilizer_1995.nc'
       filename2= 'N_FERTIL_FILE'
       density  = .FALSE. 
       fieldname= 'N_qt_WRICE_year'
       CALL diffuco_inca_read (kjpindex, lalo, neighbours, resolution, N_qt_WRICE_year, filename, filename2, fieldname, density)
       fieldname= 'N_qt_OTHER_year'
       CALL diffuco_inca_read (kjpindex, lalo, neighbours, resolution, N_qt_OTHER_year, filename, filename2, fieldname, density)
    END IF

    ! If we acount for NOx emissions due to Biomass Burning
    IF (control%ok_bbgfertil_NOx) THEN

       ALLOCATE (flx_co2_bbg_year(kjpindex),stat=ier) !! CO2 emissions from bbg, read in file 
       IF (ier /= 0) THEN
          WRITE (numout,*) ' error in flx_co2_bbg_year allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'diffuco_inca_init'
       END IF
       flx_co2_bbg_year(:) = zero    

       WRITE (numout,*) ' *********************** Interpolating CO2 bbg files for NOx emissions... '
       filename = 'orchidee_bbg_clim.nc'
       filename2= 'CO2_BBG_FILE'
       fieldname= 'flx_co2_bbg_year'
       density  = .TRUE.
       CALL diffuco_inca_read (kjpindex, lalo, neighbours, resolution, flx_co2_bbg_year,filename,filename2,fieldname,density)
    END IF

    IF ( OFF_LINE_MODE ) THEN

       !-
       !- What are the alowed options for the temportal interpolation
       !-
       !Config Key   = ALLOW_WEATHERGEN
       !Config Desc  = Allow weather generator to create data
       !Config If    = 
       !Config Def   = n
       !Config Help  = This flag allows the forcing-reader to generate
       !Config         synthetic data if the data in the file is too sparse
       !Config         and the temporal resolution would not be enough to
       !Config         run the model.
       !Config Units = [FLAG]
       !-
       allow_weathergen = .FALSE.
       CALL getin_p('ALLOW_WEATHERGEN',allow_weathergen)
       
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
       ! Filename should be parameterized
       ! To Be Updated
       CALL getin_p('FORCING_FILE',filename)
       CALL flininfo(filename,iml, jml, lml, tml, force_id)   
       WRITE(*,*) 'Number of data per year in forcing file :', tml 
       CALL flinclo(force_id)
       WRITE(*,*) 'Forcing file closed in INCA'
       
       
       IF ( allow_weathergen ) THEN
          WRITE(*,*) '**INCA: Using weather generator, careful to precip division for NOx '
          nbre_precip = un
          WRITE(*,*) 'Division pour les precip, NOx:', nbre_precip
       ELSE
          WRITE(*,*) 'DTRADIA :', dtradia 
          nbre_precip = (one_day/dtradia)/(tml/365.)
          WRITE(*,*) 'Division pour les precip, NOx:', nbre_precip
       END IF

    ELSE ! (in coupled mode)

       nbre_precip = un
       
    END IF  ! (OFF_LINE_MODE)

       
  END SUBROUTINE diffuco_inca_init


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_clear
!!
!>\BRIEF					Housekeeping module to deallocate the variables
!! leaf_ci, rstruct and raero
!!
!! DESCRIPTION				        : Housekeeping module to deallocate the variables
!! leaf_ci, rstruct and raero
!!
!! RECENT CHANGE(S)                             : None
!!
!! MAIN OUTPUT VARIABLE(S)	                : None
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_clear()

    l_first_diffuco=.TRUE.

    IF (ALLOCATED (leaf_ci)) DEALLOCATE (leaf_ci)
    IF (ALLOCATED (raero)) DEALLOCATE (raero)

  END SUBROUTINE diffuco_clear


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_aero
!!
!>\BRIEF	This module first calculates the surface drag 
!! coefficient, for cases in which the surface drag coefficient is NOT provided by the coupled 
!! atmospheric model LMDZ or when the flag ldq_cdrag_from_gcm is set to FALSE 
!!
!! DESCRIPTION	: Computes the surface drag coefficient, for cases 
!! in which it is NOT provided by the coupled atmospheric model LMDZ. The module first uses the 
!! meteorolgical input to calculate the Richardson Number, which is an indicator of atmospheric 
!! stability in the surface layer. The formulation used to find this surface drag coefficient is 
!! dependent on the stability determined. \n
!!
!! Designation of wind speed
!! \latexonly 
!!     \input{diffucoaero1.tex}
!! \endlatexonly
!!
!! Calculation of geopotential. This is the definition of Geopotential height (e.g. Jacobson 
!! eqn.4.47, 2005). (required for calculation of the Richardson Number)
!! \latexonly 
!!     \input{diffucoaero2.tex}
!! \endlatexonly
!! 
!! \latexonly 
!!     \input{diffucoaero3.tex}
!! \endlatexonly
!!
!! Calculation of the virtual air temperature at the surface (required for calculation
!! of the Richardson Number)
!! \latexonly 
!!     \input{diffucoaero4.tex}
!! \endlatexonly
!!
!! Calculation of the virtual surface temperature (required for calculation of th
!! Richardson Number)
!! \latexonly 
!!     \input{diffucoaero5.tex}
!! \endlatexonly
!!
!! Calculation of the squared wind shear (required for calculation of the Richardson
!! Number)
!! \latexonly 
!!     \input{diffucoaero6.tex}
!! \endlatexonly
!! 
!! Calculation of the Richardson Number. The Richardson Number is defined as the ratio 
!! of potential to kinetic energy, or, in the context of atmospheric science, of the
!! generation of energy by wind shear against consumption
!! by static stability and is an indicator of flow stability (i.e. for when laminar flow 
!! becomes turbulent and vise versa). It is approximated using the expression below:
!! \latexonly 
!!     \input{diffucoaero7.tex}
!! \endlatexonly
!!
!! The Richardson Number hence calculated is subject to a minimum value:
!! \latexonly 
!!     \input{diffucoaero8.tex}
!! \endlatexonly
!! 
!! Computing the drag coefficient. We add the add the height of the vegetation to the 
!! level height to take into account that the level 'seen' by the vegetation is actually 
!! the top of the vegetation. Then we we can subtract the displacement height.
!! \latexonly 
!!     \input{diffucoaero9.tex}
!! \endlatexonly
!! 
!! For the stable case (i.e $R_i$ $\geq$ 0)
!! \latexonly 
!!     \input{diffucoaero10.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucoaero11.tex}
!! \endlatexonly
!!          
!! For the unstable case (i.e. $R_i$ < 0)
!! \latexonly 
!!     \input{diffucoaero12.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucoaero13.tex}
!! \endlatexonly
!!               
!! If the Drag Coefficient becomes too small than the surface may uncouple from the atmosphere.
!! To prevent this, a minimum limit to the drag coefficient is defined as:
!!
!! \latexonly 
!!     \input{diffucoaero14.tex}
!! \endlatexonly
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): q_cdrag
!!
!! REFERENCE(S)	: 
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Jacobson M.Z., Fundamentals of Atmospheric Modeling (2nd Edition), published Cambridge 
!! University Press, ISBN 0-521-54865-9
!!
!! FLOWCHART	:
!! \latexonly 
!!     \includegraphics[scale=0.5]{diffuco_aero_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_aero (kjpindex, kjit, u, v, zlev, z0, roughheight, temp_sol, temp_air, &
       &                   qsurf, qair, q_cdrag, snow)
! Isa: ajout snow
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex, kjit   !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u                !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: v                !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev             !! Height of first atmospheric layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: z0               !! Surface roughness Length (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: roughheight      !! Effective roughness height (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol         !! Ground temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_air         !! Lowest level temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qsurf            !! near surface specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qair             !! Lowest level specific air humidity (kg kg^{-1})
!Isa
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: snow             !! snow

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)    :: q_cdrag          !! Product of Surface drag coefficient and wind speed 
                                                                            !! (m s^{-1})

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv
    REAL(r_std)                                         :: speed, zg, zdphi, ztvd, ztvs, zdu2
    REAL(r_std)                                         :: zri, cd_neut, zscf, cd_tmp
!_ ================================================================================================================================
!Isa
    REAL(r_std)                                    :: snowfact

  !! 1. Initialisation

    ! test if we have to work with q_cdrag or to calcul it
    DO ji=1,kjpindex
       
       !! 1a).1 Designation of wind speed
       !! \latexonly 
       !!     \input{diffucoaero1.tex}
       !! \endlatexonly
       speed = wind(ji)
    
       !! 1a).2 Calculation of geopotentiel
       !! This is the definition of Geopotential height (e.g. Jacobson eqn.4.47, 2005). (required
       !! for calculation of the Richardson Number)
       !! \latexonly 
       !!     \input{diffucoaero2.tex}
       !! \endlatexonly
       zg = zlev(ji) * cte_grav
      
       !! \latexonly 
       !!     \input{diffucoaero3.tex}
       !! \endlatexonly
       zdphi = zg/cp_air
       
       !! 1a).3 Calculation of the virtual air temperature at the surface 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero4.tex}
       !! \endlatexonly
       ztvd = (temp_air(ji) + zdphi / (un + rvtmp2 * qair(ji))) * (un + retv * qair(ji)) 
       
       !! 1a).4 Calculation of the virtual surface temperature 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero5.tex}
       !! \endlatexonly
       ztvs = temp_sol(ji) * (un + retv * qsurf(ji))
     
       !! 1a).5 Calculation of the squared wind shear 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero6.tex}
       !! \endlatexonly
       zdu2 = MAX(cepdu2,speed**2)
       
       !! 1a).6 Calculation of the Richardson Number
       !!  The Richardson Number is defined as the ratio of potential to kinetic energy, or, in the 
       !!  context of atmospheric science, of the generation of energy by wind shear against consumption
       !!  by static stability and is an indicator of flow stability (i.e. for when laminar flow 
       !!  becomes turbulent and vise versa).\n
       !!  It is approximated using the expression below:
       !!  \latexonly 
       !!     \input{diffucoaero7.tex}
       !! \endlatexonly
       zri = zg * (ztvd - ztvs) / (zdu2 * ztvd)
      
       !! The Richardson Number hence calculated is subject to a minimum value:
       !! \latexonly 
       !!     \input{diffucoaero8.tex}
       !! \endlatexonly       
       zri = MAX(MIN(zri,5.),-5.)
       
       !! 1a).7 Computing the drag coefficient
       !!  We add the add the height of the vegetation to the level height to take into account
       !!  that the level 'seen' by the vegetation is actually the top of the vegetation. Then we 
       !!  we can subtract the displacement height.
       !! \latexonly 
       !!     \input{diffucoaero9.tex}
       !! \endlatexonly
       !Isa
	! 7.0 snow smoothering
	! Snow induces low levels of turbulence.
	! Sensible heat fluxes can therefore be reduced of ~1/3. Pomeroy et al., 1998
	snowfact=1.
	if (snow(ji).GT.snowcri.AND.ok_snowfact)  snowfact=10.

       cd_neut = (ct_karman / LOG( (zlev(ji) + roughheight(ji)) / z0(ji) )) ** 2 
       
       !! 1a).7.1 - for the stable case (i.e $R_i$ $\geq$ 0)
       IF (zri .GE. zero) THEN
          
          !! \latexonly 
          !!     \input{diffucoaero10.tex}
          !! \endlatexonly
          zscf = SQRT(un + cd * ABS(zri))
         
          !! \latexonly 
          !!     \input{diffucoaero11.tex}
          !! \endlatexonly          
          cd_tmp=cd_neut/(un + trois * cb * zri * zscf)
       ELSE
          
          !! 1a).7.2 - for the unstable case (i.e. $R_i$ < 0)
          !! \latexonly 
          !!     \input{diffucoaero12.tex}
          !! \endlatexonly
          zscf = un / (un + trois * cb * cc * cd_neut * SQRT(ABS(zri) * &
               & ((zlev(ji) + roughheight(ji)) / z0(ji)/snowfact)))
          
          !! \latexonly 
          !!     \input{diffucoaero13.tex}
          !! \endlatexonly               
          cd_tmp=cd_neut * (un - trois * cb * zri * zscf)
       ENDIF
       
       !! If the Drag Coefficient becomes too small than the surface may uncouple from the atmosphere.
       !! To prevent this, a minimum limit to the drag coefficient is defined as:
       
       !! \latexonly 
       !!     \input{diffucoaero14.tex}
       !! \endlatexonly
       !!
       q_cdrag(ji) = MAX(cd_tmp, 1.e-4/MAX(speed,min_wind))

       ! In some situations it might be useful to give an upper limit on the cdrag as well. 
       ! The line here should then be uncommented.
      !q_cdrag(ji) = MIN(q_cdrag(ji), 0.5/MAX(speed,min_wind))

    END DO

    IF (long_print) WRITE (numout,*) ' not ldqcdrag_from_gcm : diffuco_aero done '

  END SUBROUTINE diffuco_aero


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_snow
!!
!>\BRIEF         This subroutine computes the beta coefficient for snow sublimation.
!!
!! DESCRIPTION   : This routine computes beta coefficient for snow sublimation, which
!! integrates the snow on both vegetation and other surface types (e.g. ice, lakes,
!! cities etc.) \n
!!
!! A critical depth of snow (snowcri) is defined to calculate the fraction of each grid-cell
!! that is covered with snow (snow/snowcri) while the remaining part is snow-free.
!! We also carry out a first calculation of sublimation (subtest) to lower down the beta
!! coefficient if necessary (if subtest > snow). This is a predictor-corrector test. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta1 
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp. 248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
SUBROUTINE diffuco_snow (kjpindex, dtradia, qair, qsatt, rau, u, v,q_cdrag, &
       & snow, frac_nobio, totfrac_nobio, snow_nobio, vbeta1)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size (-)
    REAL(r_std), INTENT (in)                             :: dtradia        !! Time step in seconds (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair           !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qsatt          !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau            !! Air density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u              !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v              !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag        !! Product of surface drag coefficient and wind speed (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow           !! Snow mass (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Snow on ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: totfrac_nobio  !! Total fraction of ice, lakes, cities etc. (-)
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vbeta1         !! Beta for sublimation (dimensionless ratio) 
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std)                                          :: subtest        !! Sublimation for test (kg m^{-2})
    REAL(r_std)                                          :: zrapp          !! Modified factor (ratio)
    REAL(r_std)                                          :: speed          !! Wind speed (m s^{-1})
    REAL(r_std)                                          :: vbeta1_add     !! Beta for sublimation (ratio)
    INTEGER(i_std)                                       :: ji, jv         !! Indices (-)
!_ ================================================================================================================================

  !! 1. Calculate beta coefficient for snow sublimation on the vegetation\n

    DO ji=1,kjpindex  ! Loop over # pixels - domain size

       ! Fraction of mesh that can sublimate snow
       vbeta1(ji) = (un - totfrac_nobio(ji)) * MAX( MIN(snow(ji)/snowcri,un),zero)

       ! Limitation of sublimation in case of snow amounts smaller than the atmospheric demand. 
       speed = MAX(min_wind, wind(ji))

       subtest = dtradia * vbeta1(ji) * speed * q_cdrag(ji) * rau(ji) * &
               & ( qsatt(ji) - qair(ji) )

       IF ( subtest .GT. zero ) THEN
          zrapp = snow(ji) / subtest
          IF ( zrapp .LT. un ) THEN
             vbeta1(ji) = vbeta1(ji) * zrapp
          ENDIF
       ENDIF

    END DO ! Loop over # pixels - domain size

  !! 2. Add the beta coefficients calculated from other surfaces types (snow on ice,lakes, cities...)

    DO jv = 1, nnobio ! Loop over # other surface types
!!$      !
!!$      IF ( jv .EQ. iice ) THEN
!!$        !
!!$        !  Land ice is of course a particular case
!!$        !
!!$        DO ji=1,kjpindex
!!$          vbeta1(ji) = vbeta1(ji) + frac_nobio(ji,jv)
!!$        ENDDO
!!$        !
!!$      ELSE
        !
        DO ji=1,kjpindex ! Loop over # pixels - domain size

           vbeta1_add = frac_nobio(ji,jv) * MAX(MIN(snow_nobio(ji,jv)/snowcri,un), zero)

           ! Limitation of sublimation in case of snow amounts smaller than
           ! the atmospheric demand. 
           speed = MAX(min_wind, wind(ji))
            
            !!     Limitation of sublimation by the snow accumulated on the ground 
            !!     A first approximation is obtained with the old values of
            !!     qair and qsol_sat: function of temp-sol and pb. (see call of qsatcalc)
           subtest = dtradia * vbeta1_add * speed * q_cdrag(ji) * rau(ji) * &
                & ( qsatt(ji) - qair(ji) )

           IF ( subtest .GT. zero ) THEN
              zrapp = snow_nobio(ji,jv) / subtest
              IF ( zrapp .LT. un ) THEN
                 vbeta1_add = vbeta1_add * zrapp
              ENDIF
           ENDIF

           vbeta1(ji) = vbeta1(ji) + vbeta1_add

        ENDDO ! Loop over # pixels - domain size

!!$      ENDIF
      
    ENDDO ! Loop over # other surface types

    IF (long_print) WRITE (numout,*) ' diffuco_snow done '

  END SUBROUTINE diffuco_snow


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_flood 
!!
!>\BRIEF				       	This routine computes partial beta coefficient : floodplains
!!
!! DESCRIPTION				        : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S)	                : vbeta5
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_flood (kjpindex, dtradia, qair, qsatt, rau, u, v, q_cdrag, evapot, evapot_corr, &
       & flood_frac, flood_res, vbeta5)

    ! interface description
    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex   !! Domain size
    REAL(r_std), INTENT (in)                                 :: dtradia    !! Time step in seconds
    ! input fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair       !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qsatt      !! Surface saturated humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau        !! Density
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u          !! Lowest level wind speed 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v          !! Lowest level wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag    !! Surface drag
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: flood_res  !! water mass in flood reservoir
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: flood_frac !! fraction of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: evapot     !! Potential evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: evapot_corr!! Potential evaporation2
    ! output fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vbeta5     !! Beta for floodplains

    ! local declaration
    REAL(r_std)                                              :: subtest, zrapp, speed
    INTEGER(i_std)                                           :: ji, jv

!_ ================================================================================================================================
    !
    ! beta coefficient for sublimation for floodplains
    !
    DO ji=1,kjpindex
       !
       IF (evapot(ji) .GT. min_sechiba) THEN
          vbeta5(ji) = flood_frac(ji) *evapot_corr(ji)/evapot(ji)
       ELSE
          vbeta5(ji) = flood_frac(ji)
       ENDIF
       !
       ! -- Limitation of evaporation in case of water amounts smaller than
       !    the atmospheric demand. 
       
       !
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
       !
       subtest = dtradia * vbeta5(ji) * speed * q_cdrag(ji) * rau(ji) * &
               & ( qsatt(ji) - qair(ji) )
       !  
       IF ( subtest .GT. zero ) THEN
          zrapp = flood_res(ji) / subtest
          IF ( zrapp .LT. un ) THEN
             vbeta5(ji) = vbeta5(ji) * zrapp
          ENDIF
       ENDIF
       !
    END DO

    IF (long_print) WRITE (numout,*) ' diffuco_flood done '

  END SUBROUTINE diffuco_flood


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_inter
!!
!>\BRIEF	 This routine computes the partial beta coefficient
!! for the interception for each type of vegetation
!!
!! DESCRIPTION   : We first calculate the dry and wet parts of each PFT (wet part = qsintveg/qsintmax).
!! The former is submitted to transpiration only (vbeta3 coefficient, calculated in 
!! diffuco_trans or diffuco_trans_co2), while the latter is first submitted to interception loss 
!! (vbeta2 coefficient) and then to transpiration once all the intercepted water has been evaporated 
!! (vbeta23 coefficient). Interception loss is also submitted to a predictor-corrector test, 
!! as for snow sublimation. \n
!!
!! \latexonly 
!!     \input{diffucointer1.tex}
!! \endlatexonly
!! Calculate the wet fraction of vegetation as  the ration between the intercepted water and the maximum water 
!! on the vegetation. This ratio defines the wet portion of vegetation that will be submitted to interception loss.
!!
!! \latexonly 
!!     \input{diffucointer2.tex}
!! \endlatexonly
!!
!! Calculation of $\beta_3$, the canopy transpiration resistance
!! \latexonly 
!!     \input{diffucointer3.tex}
!! \endlatexonly            
!! 
!! We here determine the limitation of interception loss by the water stored on the leaf. 
!! A first approximation of interception loss is obtained using the old values of
!! qair and qsol_sat, which are functions of temp-sol and pb. (see call of 'qsatcalc')
!! \latexonly 
!!     \input{diffucointer4.tex}
!! \endlatexonly
!!
!! \latexonly
!!     \input{diffucointer5.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucointer6.tex}
!! \endlatexonly
!!
!! Once the whole water stored on foliage has evaporated, transpiration can take place on the fraction
!! 'zqsvegrap'.
!! \latexonly 
!!     \input{diffucointer7.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta2, ::vbeta23
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp. 248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Perrier, A, 1975. Etude physique de l'évaporation dans les conditions naturelles. Annales 
!! Agronomiques, 26(1-18): pp. 105-123, pp. 229-243
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_inter (kjpindex, dtradia, qair, qsatt, rau, u, v, q_cdrag, veget, &
     & qsintveg, qsintmax, rstruct, vbeta2, vbeta23)
   
  !! 0 Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjpindex   !! Domain size (-)
    REAL(r_std), INTENT (in)                             :: dtradia    !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qsatt      !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u          !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v          !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag    !! Product of Surface drag coefficient and wind 
                                                                       !! speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget      !! vegetation fraction for each type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintveg   !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: rstruct    !! architectural resistance (s m^{-1})
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vbeta2     !! Beta for interception loss (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vbeta23    !! Beta for fraction of wetted foliage that will 
                                                                       !! transpire (-)

    !! 0.4 Local variables

    INTEGER(i_std)                                       :: ji, jv                               !! (-), (-)
    REAL(r_std)                                          :: zqsvegrap, ziltest, zrapp, speed     !!
!_ ================================================================================================================================

  !! 1. Initialize

    vbeta2(:,:) = zero
    vbeta23(:,:) = zero
   
  !! 2. The beta coefficient for interception by vegetation. 
    
    DO jv = 2,nvm

      DO ji=1,kjpindex

         IF (veget(ji,jv) .GT. min_sechiba .AND. qsintveg(ji,jv) .GT. zero ) THEN

            zqsvegrap = zero
            IF (qsintmax(ji,jv) .GT. min_sechiba ) THEN

            !! \latexonly 
            !!     \input{diffucointer1.tex}
            !! \endlatexonly
            !!
            !! We calculate the wet fraction of vegetation as  the ration between the intercepted water and the maximum water 
            !! on the vegetation. This ratio defines the wet portion of vegetation that will be submitted to interception loss.
            !!
                zqsvegrap = MAX(zero, qsintveg(ji,jv) / qsintmax(ji,jv))
            END IF

            !! \latexonly 
            !!     \input{diffucointer2.tex}
            !! \endlatexonly
            speed = MAX(min_wind, wind(ji))

            !! Calculation of $\beta_3$, the canopy transpiration resistance
            !! \latexonly 
            !!     \input{diffucointer3.tex}
            !! \endlatexonly
            vbeta2(ji,jv) = veget(ji,jv) * zqsvegrap * (un / (un + speed * q_cdrag(ji) * rstruct(ji,jv)))
            
            !! We here determine the limitation of interception loss by the water stored on the leaf. 
            !! A first approximation of interception loss is obtained using the old values of
            !! qair and qsol_sat, which are functions of temp-sol and pb. (see call of 'qsatcalc')
            !! \latexonly 
            !!     \input{diffucointer4.tex}
            !! \endlatexonly
            ziltest = dtradia * vbeta2(ji,jv) * speed * q_cdrag(ji) * rau(ji) * &
               & ( qsatt(ji) - qair(ji) )

            IF ( ziltest .GT. zero ) THEN

                !! \latexonly 
                !!     \input{diffucointer5.tex}
                !! \endlatexonly
                zrapp = qsintveg(ji,jv) / ziltest
                IF ( zrapp .LT. un ) THEN
                   
                    !! \latexonly 
                    !!     \input{diffucointer6.tex}
                    !! \endlatexonly
                    !!
		    !! Once the whole water stored on foliage has evaporated, transpiration can take place on the fraction
                    !! 'zqsvegrap'.
                    vbeta23(ji,jv) = MAX(vbeta2(ji,jv) - vbeta2(ji,jv) * zrapp, zero)
                    
                    !! \latexonly 
                    !!     \input{diffucointer7.tex}
                    !! \endlatexonly
                    vbeta2(ji,jv) = vbeta2(ji,jv) * zrapp
                ENDIF
            ENDIF
        END IF
!        ! Autre formulation possible pour l'evaporation permettant une transpiration sur tout le feuillage
!        !commenter si formulation Nathalie sinon Tristan
!        speed = MAX(min_wind, wind(ji))
!        
!        vbeta23(ji,jv) = MAX(zero, veget(ji,jv) * (un / (un + speed * q_cdrag(ji) * rstruct(ji,jv))) - vbeta2(ji,jv))

      END DO

    END DO

    IF (long_print) WRITE (numout,*) ' diffuco_inter done '

  END SUBROUTINE diffuco_inter


!! ==============================================================================================================================
!! SUBROUTINE      : diffuco_bare
!!
!>\BRIEF	   This routine computes the partial beta coefficient corresponding to
!! bare soil
!!
!! DESCRIPTION	   : Bare soil evaporation is either limited by a soil resistance 
!! (rsol) that is calculated in hydrolc.f90, when Choisnel hydrology is used or 
!! submitted to a maximum possible flow (evap_bare_lim) if the 11-layer hydrology is used.\n
!! 
!! Calculation of wind speed
!! \latexonly 
!!     \input{diffucobare1.tex}
!! \endlatexonly
!!             
!! The calculation of $\beta_4$
!! \latexonly 
!!     \input{diffucobare2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta4
!!
!! REFERENCE(S)	 :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_bare (kjpindex, dtradia, u, v, q_cdrag, rsol, evap_bare_lim, evapot, humrel, &
       & veget, veget_max, vbeta2, vbeta3, vbeta4)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjpindex       !! Domain size (-)
    REAL(r_std), INTENT (in)                           :: dtradia        !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u              !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v              !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q_cdrag        !! Product of Surface drag coefficient and wind speed 
                                                                         !! (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rsol           !! resistance for bare soil evaporation  (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evap_bare_lim  !! limiting factor for bare soil evaporation when the 
                                                                         !! 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot         !! Soil Potential Evaporation (??mm day^{-1}??) 
                                                                         !! NdN Variable that does not seem to be used at 
                                                                         !! all in diffuco!!
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel         !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget          !! Type of vegetation fraction (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max      !! Type of vegetation max fraction
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vbeta2         !! Beta for Interception 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vbeta3         !! Beta for Transpiration 

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: vbeta4         !! Beta for bare soil evaporation (-)
    
    !! 0.3 Modified variables
        
    !! 0.4 Local variables
    REAL(r_std)                                    :: humveg_prod

    INTEGER(i_std)                                     :: ji, jv
    REAL(r_std)                                        :: speed          !! Surface wind speed (m s^{-1})
!_ ================================================================================================================================

  !! 1. Calculation of the soil resistance and the beta (beta_4) for bare soil

    IF ( .NOT. control%hydrol_cwrr ) THEN
       DO ji = 1, kjpindex
          
          vbeta4(ji) = zero
          !     
          ! 1.   Soil resistance and beta for bare soil
          !      note: tot_bare_soil contains the fraction of bare soil
          !            see slowproc module
          !
          speed = MAX(min_wind, wind(ji))
!          write(*,*) 'diffuco 1736: ji,speed=',ji,speed
          !
          humveg_prod = tot_bare_soil(ji) * humrel(ji,1)
          !
          DO jv = 2, nvm
             humveg_prod = humveg_prod + veget(ji,jv) * humrel(ji,jv)
          ENDDO
            
             !! \latexonly 
             !!     \input{diffucobare1.tex}
             !! \endlatexonly
          IF (tot_bare_soil(ji) .GE. min_sechiba) THEN
             
             ! Correction Nathalie de Noblet - le 27 Mars 2006
             ! Selon recommandation de Frederic Hourdin: supprimer humrel dans formulation vbeta4
             !vbeta4(ji) = tot_bare_soil(ji) *humrel(ji,1)* (un / (un + speed * q_cdrag(ji) * rsol(ji)))
             ! Nathalie - le 28 mars 2006 - vbeta4 n'etait pas calcule en fonction de
             ! rsol mais de rsol_cste * hdry! Dans ce cas inutile de calculer rsol(ji)!!
             vbeta4(ji) = tot_bare_soil(ji) * (un / (un + speed * q_cdrag(ji) * rsol(ji)))
             
          ENDIF
          !Commenter la ligne ci-dessous si calcul Nathalie sinon Tristan
!          vbeta4(ji) = MIN(humveg_prod * (un / (un + speed * q_cdrag(ji) * rsol(ji))), &
!               & un - SUM(vbeta2(ji,:)+vbeta3(ji,:)))
          
       END DO
    ELSE
       DO ji = 1, kjpindex

          ! The limitation by 1-beta2-beta3 is due to the fact that evaporation under vegetation is possible
          !! \latexonly 
          !!     \input{diffucobare3.tex}
          !! \endlatexonly
          vbeta4(ji) = MIN(evap_bare_lim(ji), un - SUM(vbeta2(ji,:)+vbeta3(ji,:)))
       END DO
    ENDIF
    
    
    IF (long_print) WRITE (numout,*) ' diffuco_bare done '
    
  END SUBROUTINE diffuco_bare


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_trans 
!!
!>\BRIEF        This routine computes the partial beta coefficient 
!! corresponding to transpiration for each vegetation type.
!!
!! DESCRIPTION  : Beta coefficients for transpiration are calculated 
!! here using Jarvis formulation for stomatal resistance and
!! the structural resistance to represent the vertical gradient of 
!! transpiration within the canopy. \n
!!
!! The Jarvis formulation as used here is derived by Lohanner et al. (1980) from Jarvis (1976). This formulation is
!! semi-empirical: \n
!!
!! \latexonly 
!!     \input{diffucotrans4.tex}
!! \endlatexonly
!! \n
!!            
!! where in this expression LAI is the single sided Leaf Area Index, R_{new}^{SW} the net shortwave radiation,
!! R_{SO} the half light saturation factor, \delta c the water vapour concentration deficit, and a, k_0 and \lambda
!! are all parameters that are derived from extensive measurement of surface layer vegetation. \n
!!
!! Structural resistance (or architectural resistance) is a function of vegetation type and is assigned based on the
!! particular Plant Functional Type (PFT) in question. The range of values for the structural resistance are listed
!! in the module 'pft_parameters', and are described in de Noblet-Ducoudré et al (1993). \n
!!
!! vbetaco2 is here set to zero as this way to compute canopy resistances is only used 
!! without STOMATE, and there is therefore no photosynthesis. \n
!!
!! Moisture concentration at the leaf level.
!! \latexonly 
!!     \input{diffucotrans1.tex}
!! \endlatexonly
!!   
!! Calulation of the beta coefficient for vegetation transpiration, beta_3.
!! \latexonly 
!!     \input{diffucotrans2.tex}
!! \endlatexonly
!! \latexonly             
!!     \input{diffucotrans3.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucotrans4.tex}
!! \endlatexonly            
!!
!! This is the formulation for beta_3.
!! \latexonly 
!!     \input{diffucotrans5.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucotrans6.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vbeta3, ::rveget, ::cimean and ::vbetaco2
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Jarvis, PG, 1976. The interpretation of the variations in leaf water potential and stomatal
!! conductance found in canopies in the fields. Philosophical Transactions of the Royal Society of
!! London, Series B, 273, pp. 593-610
!! - Lohammer T, Larsson S, Linder S & Falk SO, 1980. Simulation models of gaseous exchange in Scotch
!! pine. Structure and function of Northern Coniferous Forest, Ecological Bulletin, 32, pp. 505-523
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_trans (kjpindex, dtradia, swnet, temp_air, pb, qair, rau, u, v, q_cdrag, humrel, &
                            veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, rveget, rstruct, &
                            cimean, vbetaco2, vbeta23)  

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                         :: kjpindex   !! Domain size (-)
    REAL(r_std), INTENT (in)                           :: dtradia    !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet      !! Short wave net flux at surface (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_air   !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb         !! Lowest level pressure (Pa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u          !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v          !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q_cdrag    !! Product of Surface drag coefficient and wind speed 
                                                                     !! (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel     !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget      !! Type of vegetation (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max  !! Max. vegetation fraction (LAI->infty) (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: lai        !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg   !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: rstruct    !! Structural resistance (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vbeta23    !! Beta for wetted foliage fraction that will transpire (-)
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3     !! Beta for Transpiration (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3pot  !! Beta for Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rveget     !! Stomatal resistance of the whole canopy (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: cimean     !! STOMATE: mean intercellular ci (see enerbil) 
                                                                     !! (\mumole m^{-2} s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbetaco2   !! STOMATE: Beta for CO2 (-)

    !! 0.3 Modified variables
  
    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji, jv
    REAL(r_std)                                        :: speed
    REAL(r_std), DIMENSION(kjpindex)                   :: zdefconc, zqsvegrap
    REAL(r_std), DIMENSION(kjpindex)                   :: qsatt
    REAL(r_std),DIMENSION (kjpindex,nvm)               :: rveget_min           !! Minimal Surface resistance of vegetation
!_ ================================================================================================================================


  !! 1.  Moisture concentration at the leaf level.
    
    CALL qsatcalc (kjpindex, temp_air, pb, qsatt)
    
    !! \latexonly 
    !!     \input{diffucotrans1.tex}
    !! \endlatexonly
    zdefconc(:) = rau(:) * MAX( qsatt(:) - qair(:), zero )

   
  !! 2. Calulation of the beta coefficient for vegetation transpiration, beta_3.

    rveget(:,:) = undef_sechiba
    rveget_min(:,:) = undef_sechiba
    vbeta3(:,:) = zero
    vbeta3pot(:,:) = zero

    DO jv = 2,nvm

      zqsvegrap(:) = zero

      DO ji = 1, kjpindex

         !! \latexonly 
         !!     \input{diffucotrans2.tex}
         !! \endlatexonly
         speed = MAX(min_wind, wind(ji))

         IF (qsintmax(ji,jv) .GT. min_sechiba) THEN
        
            !! \latexonly 
            !!     \input{diffucotrans3.tex}
            !! \endlatexonly
            zqsvegrap(ji) = MAX(zero, qsintveg(ji,jv) / qsintmax(ji,jv))
         ENDIF
         
         IF ( ( veget(ji,jv)*lai(ji,jv) .GT. min_sechiba ) .AND. &
              ( kzero(jv) .GT. min_sechiba ) .AND. &
              ( swnet(ji) .GT. min_sechiba ) ) THEN

            !! \latexonly 
            !!     \input{diffucotrans4.tex}
            !! \endlatexonly            
            rveget(ji,jv) = (( swnet(ji) + rayt_cste ) / swnet(ji) ) &
                 * ((defc_plus + (defc_mult * zdefconc(ji) )) / kzero(jv)) * (un / lai(ji,jv))

            rveget_min(ji,jv) = (defc_plus / kzero(jv)) * (un / lai(ji,jv))

            ! Corrections Nathalie - le 28 mars 2006 - sur conseils Fred Hourdin
            ! Introduction d'un potentiometre (rveg_pft) pour regler la somme rveg+rstruct
            ! vbeta3(ji,jv) = veget(ji,jv) * (un - zqsvegrap(ji)) * humrel(ji,jv) * &
            !     (un / (un + speed * q_cdrag(ji) * (rveget(ji,jv) + rstruct(ji,jv))))
            
            !! This is the formulation for $beta_3$.
            !! \latexonly 
            !!     \input{diffucotrans5.tex}
            !! \endlatexonly
            vbeta3(ji,jv) = veget(ji,jv) * (un - zqsvegrap(ji)) * humrel(ji,jv) * &
                 (un / (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget(ji,jv) + rstruct(ji,jv)))))
           
            ! Fin ajout Nathalie
            ! Ajout Nathalie - Juin 2006

            !! \latexonly 
            !!     \input{diffucotrans6.tex}
            !! \endlatexonly
            vbeta3(ji,jv) = vbeta3(ji,jv) + MIN( vbeta23(ji,jv), &
                 veget(ji,jv) * zqsvegrap(ji) * humrel(ji,jv) * &
                 (un / (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget(ji,jv) + rstruct(ji,jv))))))
            ! Fin ajout Nathalie
            ! Autre possibilite permettant la transpiration sur toute la canopee
            !Commenter si formulation Nathalie sinon Tristan
!            vbeta3(ji,jv) = MAX(zero, MIN(vbeta23(ji,jv), &
!                 & veget_max(ji,jv) * humrel(ji,jv) / &
!                 & (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget(ji,jv) + rstruct(ji,jv))))))

           ! vbeta3pot for computation of potential transpiration (needed for irrigation)
            vbeta3pot(ji,jv) = &
                 &  MAX(zero, veget_max(ji,jv) / &
                 & (un + speed * q_cdrag(ji) * (rveg_pft(jv)*(rveget_min(ji,jv) + rstruct(ji,jv)))))
         ENDIF

      ENDDO

    ENDDO

    ! STOMATE
    cimean(:,:) = zero
    vbetaco2(:,:) = zero

    IF (long_print) WRITE (numout,*) ' diffuco_trans done '

  END SUBROUTINE diffuco_trans


!! ==============================================================================================================================
!! SUBROUTINE   : diffuco_trans_co2
!!
!>\BRIEF        This subroutine computes carbon assimilation and stomatal 
!! conductance, following respectively Farqhuar et al. (1980) and Ball et al. (1987).
!!
!! DESCRIPTION  :\n
!! *** General:\n 
!! The equations are different depending on the photosynthesis mode (C3 versus C4).\n 
!! Assimilation and conductance are computed over 20 levels of LAI and then 
!! integrated at the canopy level.\n 
!! This routine also computes partial beta coefficient: transpiration for each 
!! type of vegetation.\n
!! There is a main loop on the PFTs, then inner loops on the points where 
!! assimilation has to be calculated.\n
!! This subroutine is called by diffuco_main only if photosynthesis is activated
!! for sechiba (flag STOMATE_OK_CO2=TRUE), otherwise diffuco_trans is called.\n
!! This subroutine is called at each sechiba time step by sechiba_main.\n
!! *** Details:
!! - Integration at the canopy level\n
!! \latexonly
!! \input{diffuco_trans_co2_1.1.tex}
!! \endlatexonly
!! - Light''s extinction \n
!! The available light follows a simple Beer extinction law. 
!! The extinction coefficients (ext_coef) are PFT-dependant constants and are defined in constant_co2.f90.\n
!! \latexonly
!! \input{diffuco_trans_co2_1.2.tex}
!! \endlatexonly
!! - Estimation of relative humidity of air (for calculation of the stomatal conductance)\n
!! \latexonly
!! \input{diffuco_trans_co2_1.3.tex}
!! \endlatexonly
!! - Calculation of the water limitation factor\n
!! \latexonly
!! \input{diffuco_trans_co2_2.1.tex}
!! \endlatexonly
!! - Calculation of temperature dependent parameters for C4 plants\n
!! \latexonly
!! \input{diffuco_trans_co2_2.2.tex}
!! \endlatexonly
!! - Calculation of temperature dependent parameters for C3 plants\n
!! \latexonly
!! \input{diffuco_trans_co2_2.3.tex}
!! \endlatexonly
!! - Vmax scaling\n 
!! Vmax is scaled into the canopy due to reduction of nitrogen 
!! (Johnson and Thornley,1984).\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.1.tex}
!! \endlatexonly
!! - Assimilation for C4 plants (Collatz et al., 1992)\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.2.tex}
!! \endlatexonly         
!! - Assimilation for C3 plants (Farqhuar et al., 1980)\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.3.tex}
!! \endlatexonly
!! - Estimation of the stomatal conductance (Ball et al., 1987)\n
!! \latexonly
!! \input{diffuco_trans_co2_2.4.4.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): N. de Noblet          2006/06
!!                - addition of q2m and t2m as input parameters for the 
!!                calculation of Rveget
!!                - introduction of vbeta23
!!
!! MAIN OUTPUT VARIABLE(S): beta coefficients, resistances, CO2 intercellular 
!! concentration
!!
!! REFERENCE(S) :
!! - Ball, J., T. Woodrow, and J. Berry (1987), A model predicting stomatal 
!! conductance and its contribution to the control of photosynthesis under 
!! different environmental conditions, Prog. Photosynthesis, 4, 221– 224.
!! - Collatz, G., M. Ribas-Carbo, and J. Berry (1992), Coupled photosynthesis 
!! stomatal conductance model for leaves of C4 plants, Aust. J. Plant Physiol.,
!! 19, 519–538.
!! - Farquhar, G., S. von Caemmener, and J. Berry (1980), A biochemical model of 
!! photosynthesis CO2 fixation in leaves of C3 species, Planta, 149, 78–90.
!! - Johnson, I. R., and J. Thornley (1984), A model of instantaneous and daily
!! canopy photosynthesis, J Theor. Biol., 107, 531 545
!! - McMurtrie, R.E., Rook, D.A. and Kelliher, F.M., 1990. Modelling the yield of Pinus radiata on a
!! site limited by water and nitrogen. For. Ecol. Manage., 30: 381-413
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

SUBROUTINE diffuco_trans_co2 (kjpindex, dtradia, swdown, temp_air, pb, qair, q2m, t2m, rau, u, v, q_cdrag, humrel, &
                                assim_param, ccanopy, &
                                veget, veget_max, lai, qsintveg, qsintmax, vbeta3, vbeta3pot, rveget, rstruct, &
                                cimean, vbetaco2, vbeta23)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !! Domain size (unitless)
    REAL(r_std), INTENT (in)                                 :: dtradia          !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown           !! Downwelling short wave flux 
                                                                                 !! @tex ($W m^{-2}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air         !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb               !! Lowest level pressure (Pa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair             !! Lowest level specific humidity 
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
! N. de Noblet - 2006/06 - addition of q2m and t2m
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q2m              !! 2m specific humidity 
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
! In off-line mode q2m and qair are the same.
! In off-line mode t2m and temp_air are the same.
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: t2m              !! 2m air temperature (K)
! N. de Noblet - 2006/06 - addition of q2m and t2m
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau              !! air density @tex ($kg m^{-3}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                !! Lowest level wind speed 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                !! Lowest level wind speed 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag          !! Surface drag 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2), INTENT (in)  :: assim_param      !! min+max+opt temps (K), vcmax, vjmax for 
                                                                                 !! photosynthesis 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: ccanopy          !! CO2 concentration inside the canopy (ppm) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: humrel           !! Soil moisture stress (0-1,unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget            !! Coverage fraction of vegetation for each PFT 
                                                                                 !! depending on LAI (0-1, unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max        !! Maximum vegetation fraction of each PFT inside 
                                                                                 !! the grid box (0-1, unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: lai              !! Leaf area index @tex ($m^2 m^{-2}$) @endtex
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: qsintveg         !! Water on vegetation due to interception 
                                                                                 !! @tex ($kg m^{-2}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: qsintmax         !! Maximum water on vegetation
                                                                                 !! @tex ($kg m^{-2}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: vbeta23          !! Beta for fraction of wetted foliage that will 
                                                                                 !! transpire (unitless) 
    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: vbeta3           !! Beta for Transpiration (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: vbeta3pot        !! Beta for Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: rveget           !! stomatal resistance of vegetation 
                                                                                 !! @tex ($s m^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: rstruct          !! structural resistance @tex ($s m^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: cimean           !! mean intercellular CO2 concentration 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: vbetaco2         !! beta for CO2 transfer (unitless)

    !
    !! 0.3 Modified variables
    !
 
    !
    !! 0.4 Local variables
    !
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: vcmax                               !! maximum rate of carboxylation 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: vjmax                               !! maximum rate of Rubisco regeneration 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: tmin                                !! minimum temperature for photosynthesis (K)
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: topt                                !! optimal temperature for photosynthesis (K)
    REAL(r_std),DIMENSION (kjpindex,nvm)  :: tmax                                !! maximum temperature for photosynthesis (K)
    INTEGER(i_std)                        :: ji, jv, jl                          !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: leaf_ci_lowest                      !! intercellular CO2 concentration at the lowest 
                                                                                 !! LAI level (ppm) 
    INTEGER(i_std), DIMENSION(kjpindex)   :: ilai                                !! counter for loops on LAI levels (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: zqsvegrap                           !! relative water quantity in the water 
                                                                                 !! interception reservoir (0-1,unitless) 
    REAL(r_std)                           :: speed                               !! wind speed @tex ($m s^{-1}$) @endtex
    ! Assimilation
    LOGICAL, DIMENSION(kjpindex)          :: assimilate                          !! where assimilation is to be calculated 
                                                                                 !! (unitless) 
    LOGICAL, DIMENSION(kjpindex)          :: calculate                           !! where assimilation is to be calculated for 
                                                                                 !! in the PFTs loop (unitless) 
    INTEGER(i_std)                        :: nic,inic,icinic                     !! counter/indices (unitless)
    INTEGER(i_std), DIMENSION(kjpindex)   :: index_calc                          !! index (unitless)
    INTEGER(i_std)                        :: nia,inia,nina,inina,iainia          !! counter/indices (unitless)
    INTEGER(i_std), DIMENSION(kjpindex)   :: index_assi,index_non_assi           !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: vc2                                 !! rate of carboxylation (at a specific LAI level) 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)      :: vj2                                 !! rate of Rubisco regeneration (at a specific LAI 
                                                                                 !! level) @tex ($\mu mol m^{-2} s^{-1}$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)      :: assimi                              !! assimilation (at a specific LAI level) 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std)                           :: x_1,x_2,x_3,x_4,x_5,x_6             !! terms of second order polynomials
                                                                                 !! (temporary variables)
    REAL(r_std), DIMENSION(kjpindex)      :: gstop                               !! stomatal conductance at topmost level 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: gs                                  !! stomatal conductance 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: Kc                                  !! Michaelis-Menten constant for the rubisco 
                                                                                 !! enzyme catalytic activity for CO2 (ppm) 
    REAL(r_std), DIMENSION(kjpindex)      :: Ko                                  !! Michaelis-Menten constant for the rubisco 
                                                                                 !! enzyme catalytic activity for O2 (ppm) 
    REAL(r_std), DIMENSION(kjpindex)      :: CP                                  !! CO2 compensation point (ppm)
    REAL(r_std), DIMENSION(kjpindex)      :: vc                                  !! rate of carboxylation 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: vj                                  !! rate of Rubisco regeneration 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: kt                                  !! temperature-dependent pseudo-first-order rate 
                                                                                 !! constant of assimilation response to CO2 (C4) 
                                                                                 !! (??) 
    REAL(r_std), DIMENSION(kjpindex)      :: rt                                  !! temperature-dependent non-photosynthetic 
                                                                                 !! respiration (C4) @tex ($\mu mol $) @endtex ?? 
    REAL(r_std), DIMENSION(kjpindex)      :: air_relhum                          !! air relative humidity at 2m 
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: water_lim                           !! water limitation factor (0-1,unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: temp_lim                            !! temperature limitation factor (0-1,unitles)
    REAL(r_std), DIMENSION(kjpindex)      :: gstot                               !! total stomatal conductance 
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: assimtot                            !! total assimilation 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: leaf_gs_top                         !! leaf stomatal conductance at topmost level 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(nlai+1)        :: laitab                              !! tabulated LAI steps @tex ($m^2 m^{-2}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)      :: qsatt                               !! surface saturated humidity at 2m (??) 
                                                                                 !! @tex ($g g^{-1}$) @endtex
    REAL(r_std), DIMENSION(nvm,nlai)      :: light                               !! fraction of light that gets through upper LAI 
                                                                                 !! levels (0-1,unitless) 
    REAL(r_std), DIMENSION(kjpindex)      :: ci_gs                               !! intermediate calculus variable 
                                                                                 !! for vectorization
    REAL(r_std)                           :: cresist                             !! coefficient for resistances (??)

! @defgroup Photosynthesis Photosynthesis
! @{   
    ! 1. Preliminary calculations\n
    REAL(r_std),DIMENSION (kjpindex,nvm)           :: rveget_min                 !! Minimal Surface resistance of vegetation

!_ ================================================================================================================================

    !
    ! 1.1 Calculate LAI steps\n
    ! The integration at the canopy level is done over nlai fixed levels.
    !! \latexonly
    !! \input{diffuco_trans_co2_1.1.tex}
    !! \endlatexonly
! @}
! @codeinc
    DO jl = 1, nlai+1
      laitab(jl) = laimax*(EXP(lai_level_depth*REAL(jl-1,r_std))-1.)/(EXP(lai_level_depth*REAL(nlai,r_std))-un)
    ENDDO
! @endcodeinc

! @addtogroup Photosynthesis
! @{   
    !
    ! 1.2 Calculate light fraction for each LAI step\n
    ! The available light follows a simple Beer extinction law. 
    ! The extinction coefficients (ext_coef) are PFT-dependant constants and are defined in constant_co2.f90.
    !! \latexonly
    !! \input{diffuco_trans_co2_1.2.tex}
    !! \endlatexonly
! @}
! @codeinc
    DO jl = 1, nlai
      DO jv = 1, nvm
        light(jv,jl) = exp( -ext_coeff(jv)*laitab(jl) )
      ENDDO
    ENDDO 
! @endcodeinc
    !
    ! Photosynthesis parameters
    !
    ! temperatures in K
    !
    tmin(:,:) = assim_param(:,:,itmin)
    tmax(:,:) = assim_param(:,:,itmax)
    topt(:,:) = assim_param(:,:,itopt)
    !
    vcmax(:,:) = assim_param(:,:,ivcmax)
    vjmax(:,:) = assim_param(:,:,ivjmax)

! @addtogroup Photosynthesis
! @{   
    !
    ! 1.3 Estimate relative humidity of air (for calculation of the stomatal conductance).\n
    !! \latexonly
    !! \input{diffuco_trans_co2_1.3.tex}
    !! \endlatexonly
! @}
    !
! N. de Noblet - 2006/06 - We use q2m/t2m instead of qair.
!    CALL qsatcalc (kjpindex, temp_air, pb, qsatt)
!    air_relhum(:) = &
!      ( qair(:) * pb(:) / (0.622+qair(:)*0.378) ) / &
!      ( qsatt(:)*pb(:) / (0.622+qsatt(:)*0.378 ) )
! @codeinc
    CALL qsatcalc (kjpindex, t2m, pb, qsatt)
    air_relhum(:) = &
      ( q2m(:) * pb(:) / (Tetens_1+q2m(:)* Tetens_2) ) / &
      ( qsatt(:)*pb(:) / (Tetens_1+qsatt(:)*Tetens_2 ) )
! @endcodeinc
! N. de Noblet - 2006/06 
    !
    !
    ! 2. beta coefficient for vegetation transpiration
    !
    rstruct(:,1) = rstruct_const(1)
    rveget(:,:) = undef_sechiba
    rveget_min(:,:) = undef_sechiba
    !
    vbeta3(:,:) = zero
    vbeta3pot(:,:) = zero
    vbetaco2(:,:) = zero
    !
    cimean(:,1) = ccanopy(:)
    !
    ! @addtogroup Photosynthesis
    ! @{   
    ! 2. Loop over vegetation types\n
    ! @} 
    !
    DO jv = 2,nvm
      !
      ! @addtogroup Photosynthesis
      ! @{   
      !
      ! 2.1 Initializations\n
      !! \latexonly
      !! \input{diffuco_trans_co2_2.1.tex}
      !! \endlatexonly
      ! @}      
      !
      ! beta coefficient for vegetation transpiration
      !
      rstruct(:,jv) = rstruct_const(jv)
      cimean(:,jv) = ccanopy(:)
      !
      !! mask that contains points where there is photosynthesis
      !! For the sake of vectorisation [DISPENSABLE], computations are done only for convenient points.
      !! nia is the number of points where the assimilation is calculated and nina the number of points where photosynthesis is not
      !! calculated (based on criteria on minimum or maximum values on LAI, vegetation fraction, shortwave incoming radiation, 
      !! temperature and relative humidity).
      !! For the points where assimilation is not calculated, variables are initialized to specific values. 
      !! The assimilate(kjpindex) array contains the logical value (TRUE/FALSE) relative to this photosynthesis calculation.
      !! The index_assi(kjpindex) array indexes the nia points with assimilation, whereas the index_no_assi(kjpindex) array indexes
      !! the nina points with no assimilation.
      nia=0
      nina=0
      !
      DO ji=1,kjpindex
        !
        
        IF ( ( lai(ji,jv) .GT. 0.01 ) .AND. &
              ( veget_max(ji,jv) .GT. min_sechiba ) ) THEN

          IF ( ( veget(ji,jv) .GT. min_sechiba ) .AND. &
               ( swdown(ji) .GT. min_sechiba ) .AND. &
               ( temp_air(ji) .GT. tmin(ji,jv) ) .AND. &
               ( temp_air(ji) .LT. tmax(ji,jv) ) .AND. &
               ( humrel(ji,jv) .GT. min_sechiba )       ) THEN
             !
             assimilate(ji) = .TRUE.
             nia=nia+1
             index_assi(nia)=ji
             !
          ELSE
             !
             assimilate(ji) = .FALSE.
             nina=nina+1
             index_non_assi(nina)=ji
             !
          ENDIF
        ELSE
          !
          assimilate(ji) = .FALSE.
          nina=nina+1
          index_non_assi(nina)=ji
          !
        ENDIF
        !
      ENDDO
      !
      gstot(:) = zero
      assimtot(:) = zero
      !
      zqsvegrap(:) = zero
      WHERE (qsintmax(:,jv) .GT. min_sechiba)
      !! relative water quantity in the water interception reservoir
          zqsvegrap(:) = MAX(zero, qsintveg(:,jv) / qsintmax(:,jv))
      ENDWHERE
      !
      !! Calculates the water limitation factor.
      WHERE ( assimilate(:) )
        water_lim(:) = MIN( 2.*humrel(:,jv), un )
      ENDWHERE
      ! give a default value of ci for all pixel that do not assimilate
      DO jl=1,nlai
         DO inina=1,nina
            leaf_ci(index_non_assi(inina),jv,jl) = ccanopy(index_non_assi(inina)) * std_ci_frac
         ENDDO
      ENDDO
      !
      ilai(:) = 1
      !
      ! Here is calculated photosynthesis (Farqhuar et al., 1980)
      ! and stomatal conductance (Ball et al., 1987)
      !
      ! Calculates temperature dependent parameters.
      !
      IF ( is_c4(jv) )  THEN
        !
        ! @addtogroup Photosynthesis
        ! @{   
        !
        ! 2.2 Calculates temperature dependent parameters for C4 plants.\n
        !! \latexonly
        !! \input{diffuco_trans_co2_2.2.tex}
        !! \endlatexonly
        ! @}           
        !
        IF (nia .GT. 0) then
          DO inia=1,nia
            !
            ! @codeinc
            x_1 = x1_coef * EXP( x1_Q10*(temp_air(index_assi(inia))-tp_00) ) 

            ! = 2.0**(((temp_air(index_assi(inia))-tp_00)-25.0)/10.0)
            !
            kt(index_assi(inia)) = kt_coef * x_1 * 1.e6
            rt(index_assi(inia)) = rt_coef(1)* x_1 / &
              ( un + EXP(rt_coef(2)*(temp_air(index_assi(inia))-tmax(index_assi(inia),jv))) )
            !
            vc(index_assi(inia)) = vcmax(index_assi(inia),jv) & 
                * vc_coef(1) * x_1 * water_lim(index_assi(inia)) / &
!                 * 0.39 * x_1  / &
                ( (un+EXP(vc_coef(2)*(tmin(index_assi(inia),jv)-temp_air(index_assi(inia))))) &
                * (un+EXP(vc_coef(2)*(temp_air(index_assi(inia))-topt(index_assi(inia),jv)))) )
            ! @endcodeinc
            !
          ENDDO
        ENDIF
        !
        IF (nina .GT. 0) then
          ! For the points where assimilation is not calculated, variables are initialized to specific values. 
          DO inina=1,nina
            kt(index_non_assi(inina)) = zero
            rt(index_non_assi(inina)) = zero
            vc(index_non_assi(inina)) = zero
           ENDDO
        ENDIF
        !
      ELSE
        !
        ! @addtogroup Photosynthesis
        ! @{   
        !
        ! 2.3 Calculates temperature dependent parameters for C3 plants.\n
        !! \latexonly
        !! \input{diffuco_trans_co2_2.3.tex}
        !! \endlatexonly
        ! @}           
        !
        IF (nia .GT. 0) then
          !
          DO inia=1,nia
            !
            ! @codeinc
            temp_lim(index_assi(inia)) = &
              (temp_air(index_assi(inia))-tmin(index_assi(inia),jv)) * &
              (temp_air(index_assi(inia))-tmax(index_assi(inia),jv))
            temp_lim(index_assi(inia)) = temp_lim(index_assi(inia)) / &
              (temp_lim(index_assi(inia))-(temp_air(index_assi(inia))-&
               topt(index_assi(inia),jv))**2)
            !
            Kc(index_assi(inia)) = kc_coef * EXP(Ko_Q10*(temp_air(index_assi(inia))-tp_00))
            !
            Ko(index_assi(inia)) = Ko_coef * Oa &
                * EXP(Ko_Q10*(temp_air(index_assi(inia))-tmin(index_assi(inia),jv))) / &
                  (temp_air(index_assi(inia))-tmin(index_assi(inia),jv))
            !
            CP(index_assi(inia)) = CP_0 * EXP( CP_temp_coef *(temp_air(index_assi(inia))-tp_00 - CP_temp_ref)/&
                                                           temp_air(index_assi(inia)) )
            !
            vc(index_assi(inia)) = vcmax(index_assi(inia),jv) * &
               temp_lim(index_assi(inia)) * water_lim(index_assi(inia))
!               temp_lim(index_assi(inia)) 
            vj(index_assi(inia)) = vjmax(index_assi(inia),jv) * &
               temp_lim(index_assi(inia)) * water_lim(index_assi(inia))
!                temp_lim(index_assi(inia))
            ! @endcodeinc
            !
          ENDDO
          !
        ENDIF
        !
        IF (nina .GT. 0) then
          ! For the points where assimilation is not calculated, variables are initialized to specific values. 
          DO inina=1,nina
            temp_lim(index_non_assi(inina)) = zero
            Kc(index_non_assi(inina)) = zero
            Ko(index_non_assi(inina)) = zero
            CP(index_non_assi(inina)) = zero
            !
            vc(index_non_assi(inina)) = zero
            vj(index_non_assi(inina)) = zero
          ENDDO
        ENDIF
      ENDIF      ! C3/C4
      !
      ! @addtogroup Photosynthesis
      ! @{   
      !
      ! 2.4 Loop over LAI discretized levels to estimate assimilation and conductance\n
      ! @}           
      !
      !! The calculate(kjpindex) array is of type logical to indicate wether we have to sum over this LAI fixed level (the LAI of
      !! the point for the PFT is lower or equal to the LAI level value). The number of such points is incremented in nic and the 
      !! corresponding point is indexed in the index_calc array.
      DO jl = 1, nlai
        !
        nic=0
        calculate(:) = .FALSE.
        !
        IF (nia .GT. 0) then
          DO inia=1,nia
            calculate(index_assi(inia)) = (laitab(jl) .LE. lai(index_assi(inia),jv) )
            IF ( calculate(index_assi(inia)) ) THEN
              nic=nic+1
              index_calc(nic)=index_assi(inia)
            ENDIF
          ENDDO
        ENDIF
        !
        ! @addtogroup Photosynthesis
        ! @{   
        !
        ! 2.4.1 Vmax is scaled into the canopy due to reduction of nitrogen 
        !! (Johnson and Thornley,1984).\n
        !! \latexonly
        !! \input{diffuco_trans_co2_2.4.1.tex}
        !! \endlatexonly
        ! @}           
        !
        x_1 = ( un - .7_r_std * ( un - light(jv,jl) ) )
        !
        IF ( nic .GT. 0 ) THEN
          DO inic=1,nic
            ! @codeinc
            vc2(index_calc(inic)) = vc(index_calc(inic)) * x_1
            vj2(index_calc(inic)) = vj(index_calc(inic)) * x_1
            ! @endcodeinc
          ENDDO
        ENDIF
        !
        IF ( is_c4(jv) )  THEN
          !
          ! @addtogroup Photosynthesis
          ! @{   
          !
          ! 2.4.2 Assimilation for C4 plants (Collatz et al., 1992)\n
          !! \latexonly
          !! \input{diffuco_trans_co2_2.4.2.tex}
          !! \endlatexonly
          ! @}           
          !
          DO ji = 1, kjpindex
            assimi(ji) = zero
          ENDDO
          !
          IF (nic .GT. 0) THEN
             DO inic=1,nic
              ! @codeinc
              ! W_to_mmol * RG_to_PAR = 2.3
              !
              x_1 = - ( vc2(index_calc(inic)) + quantum_yield * W_to_mmol * RG_to_PAR* swdown(index_calc(inic)) * &
                    ext_coeff(jv) * light(jv,jl) )
              x_2 = vc2(index_calc(inic)) * quantum_yield *W_to_mmol * RG_to_PAR * swdown(index_calc(inic)) * &
                    ext_coeff(jv) * light(jv,jl) 
              x_3 = ( -x_1 - sqrt( x_1*x_1 - 4.0 * xc4_1 * x_2 ) ) / (2.0*xc4_1)
              x_4 = - ( x_3 + kt(index_calc(inic)) * leaf_ci(index_calc(inic),jv,jl) * &
                    1.0e-6 )
              x_5 = x_3 * kt(index_calc(inic)) * leaf_ci(index_calc(inic),jv,jl) * 1.0e-6
              assimi(index_calc(inic)) = ( -x_4 - sqrt( x_4*x_4 - 4. * xc4_2 * x_5 ) ) / (2.*xc4_2)
              assimi(index_calc(inic)) = assimi(index_calc(inic)) - &
                                                      rt(index_calc(inic))
              ! @endcodeinc
            ENDDO
          ENDIF
        ELSE
          !
          ! @addtogroup Photosynthesis
          ! @{   
          !
          ! 2.4.3 Assimilation for C3 plants (Farqhuar et al., 1980)\n
          !! \latexonly
          !! \input{diffuco_trans_co2_2.4.3.tex}
          !! \endlatexonly
          ! @}           
          !
          DO ji = 1, kjpindex
            assimi(ji) = zero
          ENDDO
          !
          IF (nic .GT. 0) THEN
            DO inic=1,nic
              ! @codeinc
              x_1 = vc2(index_calc(inic)) * leaf_ci(index_calc(inic),jv,jl) / &
                    ( leaf_ci(index_calc(inic),jv,jl) + Kc(index_calc(inic)) * &
                    ( un + Oa / Ko(index_calc(inic)) ) )

              ! The value alpha_j = 0.8855 is equal to 0.385*4.6/2.
              ! The value 0.385 comes from MacMutrie et al. (1990), to calculate J based on the equation:
	      ! theta*J^2 - (0.385*theta+Jmax)*J+O.385*theta*Jmax=0 equation that comes
              ! from Farqhuar and Wong (1984).
              ! For C4 plants PAR is calculated from SWdown by dividing by 2 and multiplying by
              ! by 4.6 to convert from W/m^2 to micromol/s/m^2.
              x_2 = alpha_j*swdown(index_calc(inic))*ext_coeff(jv)*light(jv,jl)

              x_3 = x_2 + vj2(index_calc(inic))
              x_4 = ( x_3 - sqrt( x_3*x_3 - (quatre*curve_assim*x_2*vj2(index_calc(inic))) ) ) / &
                    (deux*curve_assim)
              x_5 = x_4 * leaf_ci(index_calc(inic),jv,jl) / &
                    ( WJ_coeff1 *  leaf_ci(index_calc(inic),jv,jl) + &
                     WJ_coeff2 *CP(index_calc(inic)) ) 
              x_6 = MIN( x_1, x_5 )
              assimi(index_calc(inic)) = x_6 * ( un  - CP(index_calc(inic))/&
                 leaf_ci(index_calc(inic),jv,jl) ) - Vc_to_Rd_ratio * vc2(index_calc(inic))
              ! @endcodeinc
            ENDDO
          ENDIF
        ENDIF
        !
        IF (nic .GT. 0) THEN
          !
          DO inic=1,nic
            !
            ! @addtogroup Photosynthesis
            ! @{   
            !
            !! 2.4.4 Estimatation of the stomatal conductance (Ball et al., 1987).\n
            !! \latexonly
            !! \input{diffuco_trans_co2_2.4.4.tex}
            !! \endlatexonly
            ! @}           
            !
            icinic=index_calc(inic)
!            gs(icinic) = water_lim(icinic) * &
            gs(icinic) = &
                 ( gsslope(jv) * assimi(icinic) * &
                 air_relhum(icinic) / ccanopy(icinic) ) &
                  + gsoffset(jv)
            gs(icinic) = MAX( gs(icinic), gsoffset(jv) )
          ENDDO
          !
          DO inic=1,nic
            icinic=index_calc(inic)
            !
            ! the new ci is calculated with
            ! dci/dt=(ccanopy-ci)*gs/1.6-A
            ! ci=ci+((ccanopy(icinic)-ci)*gs/1.6-&
            !    assimi(icinic))*dtradia
            ! we verify that ci is not out of possible values
            !
            ci_gs(icinic) = MIN( ccanopy(icinic), MAX( CP(icinic), &
                   ( ccanopy(icinic) - O2toCO2_stoechio * assimi(icinic) / &
                     gs(icinic) ) ) ) - leaf_ci(icinic,jv,jl)
          ENDDO
          DO inic=1,nic
            icinic=index_calc(inic)
            !to avoid some problem of numerical stability, the leaf_ci is bufferized
            !! [DISPENSABLE] This trick makes the code uncomprehensible:
            !! ci_gs = (Ca - 1.6*A)/gs -leaf_ci
            !! leaf_ci = leaf_ci + ci_gs/6. = leaf_ci + ((Ca - 1.6*A)/gs - leaf_ci)/6. = leaf_ci*5/6. + (Ca - 1.6*A)/gs*1/6.
            leaf_ci(icinic,jv,jl) = leaf_ci(icinic,jv,jl) + ci_gs(icinic)/6.
          ENDDO
          !
          DO inic=1,nic
            icinic=index_calc(inic)
            !
            ! this might be the last level for which Ci is calculated. Store it for 
            ! initialization of the remaining levels of the Ci array.
            !
            leaf_ci_lowest(icinic) = leaf_ci(icinic,jv,jl)
          ENDDO
          !
          DO inic=1,nic
            icinic=index_calc(inic)
            !
            ! @addtogroup Photosynthesis
            ! @{   
            !
            !! 2.4.5 Integration at the canopy level\n
            !! \latexonly
            !! \input{diffuco_trans_co2_2.4.5.tex}
            !! \endlatexonly
            ! @}           
            ! total assimilation and conductance
            assimtot(icinic) = assimtot(icinic) + &
              assimi(icinic) * (laitab(jl+1)-laitab(jl))
            gstot(icinic) = gstot(icinic) + &
              gs(icinic) * (laitab(jl+1)-laitab(jl))
            !
            ilai(icinic) = jl
            !
          ENDDO
          !
        ENDIF
        !
        ! keep stomatal conductance of topmost level
        !
        IF ( jl .EQ. 1 ) THEN
          !
          leaf_gs_top(:) = zero
          !
          IF ( nic .GT. 0 ) then
            DO inic=1,nic
              leaf_gs_top(index_calc(inic)) = gs(index_calc(inic))
             ENDDO
          ENDIF
          !
        ENDIF
        !
        IF (nia .GT. 0) THEN
          !
          DO inia=1,nia
            !
            IF ( .NOT. calculate(index_assi(inia)) ) THEN
              !
              !   a) for plants that are doing photosynthesis, but whose LAI is lower
              !      than the present LAI step, initialize it to the Ci of the lowest
              !      canopy level
              !
              leaf_ci(index_assi(inia),jv,jl) = leaf_ci_lowest(index_assi(inia))
              !
            ENDIF
            !
          ENDDO
          !
        ENDIF
        !
      ENDDO  ! loop over LAI steps
      !
      !! 2.5 Calculate resistances
      !
      IF (nia .GT. 0) THEN
        !
        DO inia=1,nia
          !
          iainia=index_assi(inia)
          !
          ! conversion from mmol/m^2/s to m/s
          !
          gstot(iainia) =  mmol_to_m_1 *(temp_air(iainia)/tp_00)*&
             (pb_std/pb(iainia))*gstot(iainia)
          gstop(iainia) =  mmol_to_m_1 * (temp_air(iainia)/tp_00)*&
             (pb_std/pb(iainia))*leaf_gs_top(iainia)*&
              laitab(ilai(iainia)+1)
          !
          rveget(iainia,jv) = un/gstop(iainia)
          rveget_min(iainia,jv) = (defc_plus / kzero(jv)) * (un / lai(iainia,jv))

          !
        ENDDO
        !
        DO inia=1,nia
          !
          iainia=index_assi(inia)
          !
          ! rstruct is the difference between rtot (=1./gstot) and rveget
          !
          ! Correction Nathalie - le 27 Mars 2006 - Interdire a rstruct d'etre negatif
          !rstruct(iainia,jv) = un/gstot(iainia) - &
          !     rveget(iainia,jv)
          rstruct(iainia,jv) = MAX( un/gstot(iainia) - &
               rveget(iainia,jv), min_sechiba)
          !
        ENDDO
        !
        DO inia=1,nia
          !
          iainia=index_assi(inia)
          !
          !! wind is a global variable of the diffuco module.
          speed = MAX(min_wind, wind(index_assi(inia)))
          !
          ! beta for transpiration
          !
          ! Corrections Nathalie - 28 March 2006 - on advices of Fred Hourdin
          !! Introduction of a potentiometer rveg_pft to settle the rveg+rstruct sum problem in the coupled mode.
          !! rveg_pft=1 in the offline mode. rveg_pft is a global variable declared in the diffuco module.
          !vbeta3(iainia,jv) = veget_max(iainia,jv) * &
          !  (un - zqsvegrap(iainia)) * &
          !  (un / (un + speed * q_cdrag(iainia) * (rveget(iainia,jv) + &
          !   rstruct(iainia,jv))))
          !! Global resistance of the canopy to evaporation
          cresist=(un / (un + speed * q_cdrag(iainia) * &
               (rveg_pft(jv)*(rveget(iainia,jv) + rstruct(iainia,jv)))))


          vbeta3(iainia,jv) = veget_max(iainia,jv) * &
            (un - zqsvegrap(iainia)) * cresist + &
!!$          ! Addition Nathalie - June 2006
!!$          vbeta3(iainia,jv) = vbeta3(iainia,jv) + &
             ! Corrections Nathalie - 09 November 2009 : veget => veget_max
!               MIN( vbeta23(iainia,jv), veget(iainia,jv) * &
               MIN( vbeta23(iainia,jv), veget_max(iainia,jv) * &
!               zqsvegrap(iainia) * humrel(iainia,jv) * &
               zqsvegrap(iainia) * cresist )
          ! Fin ajout Nathalie


          ! vbeta3pot for computation of potential transpiration (needed for irrigation)
          vbeta3pot(iainia,jv) = MAX(zero, veget_max(iainia,jv) * cresist)
          !
          ! beta for assimilation. The difference is that surface 
          ! covered by rain (un - zqsvegrap(iainia)) is not taken into account
          ! The factor 1.6 is conversion for H2O to CO2 conductance.
          ! vbetaco2(iainia,jv) = veget_max(iainia,jv) * &
          !   (un / (un + q_cdrag(iainia) * &
          !   (rveget(iainia,jv))))/1.6
          !

          !! Global resistance of the canopy to CO2 transfer
          vbetaco2(iainia,jv) = veget_max(iainia,jv) * &
             (un / (un + speed * q_cdrag(iainia) * &
             (rveget(iainia,jv) + rstruct(iainia,jv)))) / O2toCO2_stoechio
          !
          ! cimean is the "mean ci" calculated in such a way that assimilation 
          ! calculated in enerbil is equivalent to assimtot
          !
          cimean(iainia,jv) = ccanopy(iainia) - &
             assimtot(iainia) / &
             ( vbetaco2(iainia,jv)/veget_max(iainia,jv) * &
             rau(iainia) * speed * q_cdrag(iainia))
          !
        ENDDO
        !
      ENDIF
      !
    END DO         ! loop over vegetation types
    !

    IF (long_print) WRITE (numout,*) ' diffuco_trans_co2 done '


END SUBROUTINE diffuco_trans_co2


!! ================================================================================================================================
!! SUBROUTINE	   : diffuco_comb
!!
!>\BRIEF           This routine combines the previous partial beta 
!! coefficients and calculates the total alpha and complete beta coefficients.
!!
!! DESCRIPTION	   : Those integrated coefficients are used to calculate (in enerbil.f90) the total evapotranspiration 
!! from the grid-cell. \n
!!
!! In the case that air is more humid than surface, dew deposition can occur (negative latent heat flux). 
!! In this instance, for temperature above zero, all of the beta coefficients are set to 0, except for 
!! interception (vbeta2) and bare soil (vbeta4 with zero soil resistance). The amount of water that is 
!! intercepted by leaves is calculated based on the value of LAI of the surface. In the case of freezing 
!! temperatures, water is added to the snow reservoir, and so vbeta4 and vbeta2 are set to 0, and the 
!! total vbeta and valpha is set to 1.\n
!!
!! \latexonly 
!!     \input{diffucocomb1.tex}
!! \endlatexonly
!!
!! The beta and alpha coefficients are initially set to 1.
!! \latexonly 
!!     \input{diffucocomb2.tex}
!! \endlatexonly
!!
!! If snow is lower than the critical value:
!! \latexonly 
!!     \input{diffucocomb3.tex}
!! \endlatexonly
!! If in the presence of dew:
!! \latexonly 
!!     \input{diffucocomb4.tex}
!! \endlatexonly
!!
!! Determine where the water goes (soil, vegetation, or snow)
!! when air moisture exceeds saturation.
!! \latexonly 
!!     \input{diffucocomb5.tex}
!! \endlatexonly
!!
!! If it is not freezing dew is put into the interception reservoir and onto the bare soil. If it is freezing, 
!! water is put into the snow reservoir. 
!! Now modify valpha and vbetas where necessary: for soil and snow
!! \latexonly 
!!     \input{diffucocomb6.tex}
!! \endlatexonly
!!
!! and for vegetation
!! \latexonly 
!!     \input{diffucocomb7.tex}
!! \endlatexonly
!!
!! Then compute part of dew that can be intercepted by leafs.
!!
!! There will be no transpiration when air moisture is too high, under any circumstance
!! \latexonly 
!!     \input{diffucocomb8.tex}
!! \endlatexonly
!!
!! There will also be no interception loss on bare soil, under any circumstance.
!! \latexonly 
!!     \input{diffucocomb9.tex}
!! \endlatexonly
!!
!! The flowchart details the 'decision tree' which underlies the module. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): vbeta1, vbeta4, humrel, vbeta2, vbeta3, valpha, vbeta
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    :
!! \latexonly 
!!     \includegraphics[scale=0.25]{diffuco_comb_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_comb (kjpindex, dtradia, humrel, rau, u, v, q_cdrag, pb, qair, temp_sol, temp_air, &
       & snow, veget, lai, vbeta1, vbeta2, vbeta3 , vbeta4, valpha, vbeta, qsintmax)    
    
    ! Ajout qsintmax dans les arguments de la routine Nathalie / le 13-03-2006

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                           :: kjpindex   !! Domain size (-)
    REAL(r_std), INTENT (in)                             :: dtradia    !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u          !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v          !! Nortward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag    !! Product of Surface drag coefficient and wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: pb         !! Lowest level pressure (Pa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol   !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_air   !! Lower air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow       !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget      !! Fraction of vegetation type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: lai        !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintmax   !! Maximum water on vegetation (kg m^{-2})

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: valpha     !! Total Alpha coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vbeta      !! Total beta coefficient (-)

    !! 0.3 Modified variables 
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: vbeta1     !! Beta for sublimation (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: vbeta4     !! Beta for Bare soil evaporation (-) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: humrel     !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: vbeta2     !! Beta for interception loss (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: vbeta3     !! Beta for Transpiration (-)
    
    !! 0.4 Local variables
    
    INTEGER(i_std)                                       :: ji, jv
    REAL(r_std)                                          :: zevtest, zsoil_moist, zrapp
    REAL(r_std), DIMENSION(kjpindex)                     :: vbeta2sum, vbeta3sum
!!    REAL(r_std), DIMENSION(kjpindex)                     :: vegetsum, vegetsum2
    REAL(r_std), DIMENSION(kjpindex)                     :: qsatt
    LOGICAL, DIMENSION(kjpindex)                         :: toveg, tosnow
    REAL(r_std)                                          :: coeff_dew_veg
!_ ================================================================================================================================
    
    !! \latexonly 
    !!     \input{diffucocomb1.tex}
    !! \endlatexonly
    vbeta2sum(:) = zero
    vbeta3sum(:) = zero
    DO jv = 1, nvm
      vbeta2sum(:) = vbeta2sum(:) + vbeta2(:,jv)
      vbeta3sum(:) = vbeta3sum(:) + vbeta3(:,jv)
    ENDDO 

  !! 1. The beta and alpha coefficients are initially set to 1.
     
    !! \latexonly 
    !!     \input{diffucocomb2.tex}
    !! \endlatexonly
    vbeta(:) = un
    valpha(:) = un

    
  !! 2. if snow is lower than the critical value:
    
    !! \latexonly 
    !!     \input{diffucocomb3.tex}
    !! \endlatexonly
    DO ji = 1, kjpindex

      IF  (snow(ji) .LT. snowcri) THEN

          vbeta(ji) = vbeta4(ji) + vbeta2sum(ji) + vbeta3sum(ji)

          IF (vbeta(ji) .LT. min_sechiba) THEN
             vbeta(ji) = zero
          END IF

      END IF

    ENDDO

    
  !! 3 If we are in presence of dew:
     
    CALL qsatcalc (kjpindex, temp_sol, pb, qsatt)

    
    !! 9.3.1 Determine where the water goes 
    !! Determine where the water goes (soil, vegetation, or snow)
    !! when air moisture exceeds saturation.
    !! \latexonly 
    !!     \input{diffucocomb5.tex}
    !! \endlatexonly
    toveg(:) = .FALSE.
    tosnow(:) = .FALSE.
    DO ji = 1, kjpindex
      IF ( qsatt(ji) .LT. qair(ji) ) THEN
          IF (temp_air(ji) .GT. tp_00) THEN

              !! 9.3.1.1  If it is not freezing, 
              !! If it is not freezing dew is put into the 
              !! interception reservoir and onto the bare soil.
              toveg(ji) = .TRUE.
          ELSE

              !! 9.3.1.2  If it is freezing, 
              !! If it is freezing water is put into the 
              !! snow reservoir.
              tosnow(ji) = .TRUE.
          ENDIF
      ENDIF
    END DO


    !! 9.3.1.3 Now modify valpha and vbetas where necessary.
    
    !! 9.3.1.3.1 Soil and snow (2d)
    !! \latexonly 
    !!     \input{diffucocomb6.tex}
    !! \endlatexonly
    DO ji = 1, kjpindex
      IF ( toveg(ji) ) THEN
        vbeta1(ji) = zero
        vbeta4(ji) = tot_bare_soil(ji)
        ! Correction Nathalie - le 13-03-2006: le vbeta ne sera calcule qu'une fois tous les vbeta2 redefinis
        !vbeta(ji) = vegetsum(ji)
        vbeta(ji) = vbeta4(ji)
        valpha(ji) = un
      ENDIF
      IF ( tosnow(ji) ) THEN
        vbeta1(ji) = un
        vbeta4(ji) = zero
        vbeta(ji) = un
        valpha(ji) = un
      ENDIF
    ENDDO

    !! 9.3.1.3.2 Vegetation (3d)
    !! \latexonly 
    !!     \input{diffucocomb7.tex}
    !! \endlatexonly
    DO jv = 1, nvm
      
      DO ji = 1, kjpindex
        
        ! Correction Nathalie - 13-03-2006 / si qsintmax=0, vbeta2=0
        IF ( toveg(ji) ) THEN
           IF (qsintmax(ji,jv) .GT. min_sechiba) THEN
              
              ! Compute part of dew that can be intercepted by leafs.
              IF ( lai(ji,jv) .GT. min_sechiba) THEN
                IF (lai(ji,jv) .GT. 1.5) THEN
                   coeff_dew_veg= &
                         &   dew_veg_poly_coeff(6)*lai(ji,jv)**5 &
                         & - dew_veg_poly_coeff(5)*lai(ji,jv)**4 &
                         & + dew_veg_poly_coeff(4)*lai(ji,jv)**3 &
                         & - dew_veg_poly_coeff(3)*lai(ji,jv)**2 &
                         & + dew_veg_poly_coeff(2)*lai(ji,jv) &
                         & + dew_veg_poly_coeff(1)
                 ELSE
                    coeff_dew_veg=un
                 ENDIF
              ELSE
                 coeff_dew_veg=zero
              ENDIF
              IF (jv .EQ. 1) THEN
                 vbeta2(ji,jv) = coeff_dew_veg*tot_bare_soil(ji)
              ELSE
                 vbeta2(ji,jv) = coeff_dew_veg*veget(ji,jv)
             !vbeta2(ji,jv) = veget(ji,jv)
              ENDIF
           ELSE
              vbeta2(ji,jv) = zero
           ENDIF
           vbeta(ji) = vbeta(ji) + vbeta2(ji,jv)
        ENDIF
        IF ( tosnow(ji) ) vbeta2(ji,jv) = zero
        
      ENDDO
      
    ENDDO

    !! 9.3.2a  There will be no transpiration when air moisture is too high, under any circumstance
    !! \latexonly 
    !!     \input{diffucocomb8.tex}
    !! \endlatexonly
    DO jv = 1, nvm
      DO ji = 1, kjpindex
        IF ( qsatt(ji) .LT. qair(ji) ) THEN
          vbeta3(ji,jv) = zero
          humrel(ji,jv) = zero
        ENDIF
      ENDDO
    ENDDO

    
    !! 9.3.2b  There will also be no interception loss on bare soil, under any circumstance.
    !! \latexonly 
    !!     \input{diffucocomb9.tex}
    !! \endlatexonly
    DO ji = 1, kjpindex
       IF ( qsatt(ji) .LT. qair(ji) ) THEN
          vbeta2(ji,1) = zero
       ENDIF
    ENDDO

    IF (long_print) WRITE (numout,*) ' diffuco_comb done '

  END SUBROUTINE diffuco_comb


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_raerod
!!
!>\BRIEF	Computes the aerodynamic resistance, for cases in which the
!! surface drag coefficient is provided by the coupled atmospheric model LMDZ and  when the flag
!! 'ldq_cdrag_from_gcm' is set to TRUE
!!
!! DESCRIPTION	: Simply computes the aerodynamic resistance, for cases in which the
!! surface drag coefficient is provided by the coupled atmospheric model LMDZ. If the surface drag coefficient
!! is not provided by the LMDZ or signalled by the flag 'ldq_cdrag_from_gcm' set to FALSE, then the subroutine
!! diffuco_aero is called instead of this one.
!!
!! Calculation of the aerodynamic resistance, for diganostic purposes. First calculate wind speed:
!! \latexonly 
!!     \input{diffucoaerod1.tex}
!! \endlatexonly       
!!
!! next calculate ::raero
!! \latexonly 
!!     \input{diffucoaerod2.tex}
!! \endlatexonly
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::raero
!!
!! REFERENCE(S)	:
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influence de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    :  None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_raerod (kjpindex, u, v, q_cdrag, raero)
    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: u            !! Eastward Lowest level wind velocity (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: v            !! Northward Lowest level wind velocity (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: q_cdrag      !! Product of Surface drag coefficient and wind speed (m s^{-1})
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out) :: raero        !! Aerodynamic resistance (s m^{-1})
     
    !! 0.3 Modified variables

    !! 0.4 Local variables
    
    INTEGER(i_std)                                 :: ji           !! (-)
    REAL(r_std)                                    :: speed        !! (m s^{-1})
!_ ================================================================================================================================
   
  !! 1. Simple calculation of the aerodynamic resistance, for diganostic purposes.

    DO ji=1,kjpindex

       !! \latexonly 
       !!     \input{diffucoaerod1.tex}
       !! \endlatexonly       
       speed = MAX(min_wind, wind(ji))

       !! \latexonly 
       !!     \input{diffucoaerod2.tex}
       !! \endlatexonly
       raero(ji) = un / (q_cdrag(ji)*speed)
       
    ENDDO
  
  END SUBROUTINE diffuco_raerod

!! ================================================================================================================================
!! SUBROUTINE   : diffuco_ok_inca
!!
!>\BRIEF         This subroutine computes biogenic emissions of reactive compounds, that is of
!!               VOCs (volatile organic compounds) from vegetation and NOx (nitrogen oxides) from soils.
!!               Calculation are mostly based on the works by Guenther et al. (1995) and Yienger and Levy (1995).\n 
!!
!! DESCRIPTION  : Biogenic VOC emissions from vegetation are based on the parameterisations developped by
!!                Guenther et al. (1995). Biogenic VOCs considered here are: isoprene, monoterpenes, OVOC and ORVOC 
!!                as bulked emissions, methanol, acetone, acetaldehyde, formaldehyde, acetic acid, formic acid
!!                as single emissions.\n
!!                For every biogenic VOCs an emission factor (EF), depending on the PFT considered, is used.\n
!!                Isoprene emissions depend on temperature and radiation. A partition between sunlit and shaded
!!                leaves is taken into account and either one (if ok_multilayer = FALSE) or several layers
!!                (if ok_multilayer = TRUE) in the canopy can be used.\n
!!                When radiation extinction is considered, the canopy radiative transfer model takes into 
!!                account light extinction through canopy, calculating first need diffuse and direct radiation
!!                based on Andrew Friend 2001 radiative model and Spitters et al. 1986. The calculation of lai, 
!!                parscat, parsh and parsun, laisun and laishabsed based on Guenther et al.(JGR, 1995) and Norman (1982).\n
!!                Emissions for other BVOCs (monoterpenes, OVOC, ORVOC and other single compounds such as
!!                methanol, acetone...) depend only on temperature.\n   
!!                The impact of leaf age, using an emission activity prescribed for each of the 4 leaf age
!!                classes can also be considered for isoprene and methanol emissions when ok_leafage = TRUE.\n
!!                NOx emissions from soils are based on Yienger and Levy (1995) and depend on soil moisture
!!                and temperature and PFT. The pulse effect, related to significant rain occuring after severe
!!                drought can also be considered (ok_pulse_NOx = TRUE), as well as the increase in emissions related to
!!                biomass buring (ok_bbgfertil_NOx = TRUE) or use of fertilizers (ok_cropsfertil_NOx = TRUE). 
!!                A net NO flux is eventually calculated taking into account loss by deposition on the surface, using
!!                a Canopy Reduction Factor (CRF) based on stomatal and leaf area.\n 
!!                This subroutine is called by diffuco_main only if biogenic emissions are activated
!!                for sechiba (flag DIFFUCO_OK_INCA=TRUE).\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: PAR, :: PARsun, :: PARsh, :: laisun, :: laish,
!!                          :: flx_iso, :: flx_mono, :: flx_ORVOC, :: flx_MBO,
!!                          :: flx_methanol, :: flx_acetone, :: flx_acetal, :: flx_formal,
!!                          :: flx_acetic, :: flx_formic, :: flx_no_soil, :: flx_no,
!!                          :: CRF, :: flx_fertil_no, :: Trans, :: Fdf,
!!                          :: PARdf, :: PARdr, :: PARsuntab, :: PARshtab
!!
!! REFERENCE(S) :
!! - Andrew Friend (2001), Modelling canopy CO2 fluxes: are 'big-leaf' simplifications justified? 
!! Global Ecology and Biogeography, 10, 6, 603-619, doi: 10.1046/j.1466-822x.2001.00268.x 
!! - Spitters, C.J.T, Toussaint, H.A.J.M, Groudriaan, J. (1986), Separating the diffuse and direct
!! component of global radiation and its implications for modeling canopy photosynthesis, Agricultural
!! and Forest Meteorology, 38, 1-3, 217-229, doi:10.1016/0168-1923(86)90060-2
!! - Norman JM (1982) Simulation of microclimates. In: Hatfield JL, Thomason IJ (eds)
!!  Biometeorology in integrated pest management. Academic, New York, pp 65–99
!! - Guenther, A., Hewitt, C. N., Erickson, D., Fall, R., Geron, C., Graedel, T., Harley, P.,
!! Klinger, L., Lerdau, M., McKay, W. A., Pierce, T., Scholes, B., Steinbrecher, R., Tallamraju,
!! R., Taylor, J. et Zimmerman, P. (1995), A global model of natural volatile organic compound
!! emissions, J. Geophys. Res., 100, 8873-8892.
!! - MacDonald, R. et Fall, R. (1993), Detection of substantial emissions of methanol from
!! plants to the atmosphere, Atmos. Environ., 27A, 1709-1713.
!! - Guenther, A., Geron, C., Pierce, T., Lamb, B., Harley, P. et Fall, R. (2000), Natural emissions
!! of non-methane volatile organic compounds, carbon monoxide, and oxides of nitrogen from
!! North America, Atmos. Environ., 34, 2205-2230.
!! - Yienger, J. J. et Levy II, H. (1995), Empirical model of global soil-biogenic NOx emissions,
!! J. Geophys. Res., 100, 11,447-11,464.
!! - Lathiere, J., D.A. Hauglustaine, A. Friend, N. De Noblet-Ducoudre, N. Viovy, and
!!  G. Folberth (2006), Impact of climate variability and land use changes on global biogenic volatile 
!! organic compound emissions, Atmospheric Chemistry and Physics, 6, 2129-2146. 
!! - Lathiere, J., D.A. Hauglustaine, N. De Noblet-Ducoudre, G. Krinner et G.A. Folberth (2005),
!! Past and future changes in biogenic volatile organic compound emissions simulated with a global
!! dynamic vegetation model, Geophysical Research Letters, 32, doi: 10.1029/2005GL024164.
!! - Lathiere, J. (2005), Evolution des emissions de composes organiques et azotes par la biosphere
!!  continentale dans le modele LMDz-INCA-ORCHIDEE, These de doctorat, Universite Paris VI.
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

  SUBROUTINE diffuco_inca (kjpindex, dtradia, swdown, sinang, temp_air, temp_sol, ptnlev1, precip_rain, humrel, &
                        veget_max, lai, frac_age, &
                        lalo, &
                        PAR, PARsun, PARsh, laisun, laish, &
                        flx_iso, flx_mono, flx_ORVOC, flx_MBO, flx_methanol, flx_acetone, flx_acetal, &
                        flx_formal, flx_acetic, flx_formic, &
                        flx_no_soil, flx_no,CRF, flx_fertil_no, Fdf, PARsuntab, PARshtab, PARdf, PARdr, Trans)


    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                 :: kjpindex         !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std), INTENT(in)                                    :: dtradia          !! Time step (seconds)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)               :: swdown           !! Down-welling surface short-wave flux 
                                                                                   !! (W.m^{-2})
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)               :: sinang           !! Sinus of Solar Angle 
                                                                                   !! (as computed in read_dim2) (unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)               :: temp_air         !! Air temperature (K)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)               :: temp_sol         !! Skin temperature (K)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)               :: ptnlev1          !! 1st level of soil temperature (K)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)               :: precip_rain      !! Rain precipitation !!?? init
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)           :: humrel           !! Soil moisture stress (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)           :: veget_max        !! Max. vegetation fraction (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)           :: lai              !! Leaf area index (m^2.m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm,nleafages), INTENT(in) :: frac_age         !! Age efficacity from STOMATE for iso 
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)             :: lalo             !! Geographical coordinates for pixels (degrees)

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(kjpindex), INTENT(out)        :: PAR              !! Photosynthetic active radiation, half of swdown
                                                                             !! @tex ($\mu mol photons. m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: PARsun           !! PAR received by sun leaves 
                                                                             !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: PARsh            !! PAR received by shaded leaves 
                                                                             !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: laisun           !! Leaf area index of Sun leaves 
                                                                             !! (m^2.m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: laish            !! Leaf area index of Shaded leaves
                                                                             !! (m^2.m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_iso          !! Biogenic isoprene emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_mono         !! Biogenic monoterpene emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_ORVOC        !! Biogenic ORVOC emission - 
                                                                             !! Other Reactive Volatil Organic Components 
                                                                             !! (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_MBO          !! Biogenic MBO emission 
                                                                             !! MethylButanOl (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_methanol     !! Biogenic methanol emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_acetone      !! Biogenic acetone emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_acetal       !! Biogenic Acetaldehyde emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_formal       !! Biogenic Formaldehyde emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_acetic       !! Biogenic Acetic Acid emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_formic       !! Biogenic Formic Acid emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_no_soil      !! Biogenic NO emission by soil (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: flx_no           !! Biogenic total NO emission (kgC.m^{-2}.s^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)    :: CRF              !! Canopy reduction factor 
                                                                             !! for net NO flux calculation 
   REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)     :: flx_fertil_no    !! Biogenic NO emission due to N-fertilisation 
                                                                             !! (kgC.m^{-2}.s^{-1})
    !!Canopy radiative transfer model
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)         :: Trans           !! Atmospheric Transmissivity (unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)         :: Fdf             !! Diffuse Fraction of the radiation (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)         :: PARdf           !! Diffuse PAR
                                                                             !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)         :: PARdr           !! Direct PAR 
                                                                             !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    !! for multilayer canopy model for iso flux
    REAL(r_std), DIMENSION(kjpindex,nlai+1), INTENT(out)  :: PARsuntab       !! PAR received by sun leaves 
                                                                             !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nlai+1), INTENT(out)  :: PARshtab        !! PAR received by shaded leaves 
                                                                             !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    !! 0.3 Modified variables

    !! 0.4 Local variables 

    INTEGER(i_std)                             :: ji, jv, jf, jl    !! Indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: fol_dens          !! foliar density (gDM.m^{-2})
    REAL(r_std), DIMENSION(kjpindex)           :: tleaf             !! Foliar temperature (K)
    REAL(r_std), DIMENSION(kjpindex)           :: t_no              !! Temperature used for soil NO emissions (C)
    REAL(r_std), DIMENSION(kjpindex)           :: exp_1             !! First exponential used in the calculation of 
                                                                    !! isoprene dependancy to Temperature 
    REAL(r_std), DIMENSION(kjpindex)           :: exp_2             !! Second exponential used in the calculation of 
                                                                    !! Isoprene dependancy to Temperature
    REAL(r_std), DIMENSION(kjpindex)           :: Ct_iso            !! Isoprene dependancy to Temperature
    REAL(r_std), DIMENSION(kjpindex)           :: Cl_iso            !! Isoprene dependancy to Light 
    REAL(r_std), DIMENSION(kjpindex)           :: Ct_mono           !! Monoterpene dependancy to Temperature
    REAL(r_std), DIMENSION(kjpindex)           :: Ct_MBO            !! MBO dependance to Temperature
    REAL(r_std), DIMENSION(kjpindex)           :: Cl_MBO            !! MBO dependance to Light
    REAL(r_std), DIMENSION(kjpindex)           :: Xvar              !! Parameter used in the calculation 
                                                                    !! of MBO dependance to Temperature
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: flx_OVOC          !! Biogenic OVOC emission - 
                                                                    !! Other Volatil Organic Components (kgC.m^{-2}.s^{-1})
    !!Canopy radiative transfer model
    REAL(r_std)                                :: day               !! Day of The Year
    REAL(r_std), DIMENSION(kjpindex)           :: So                !! Maximum radiation at the Earth surface (W.m^{-2})
    REAL(r_std), DIMENSION(kjpindex)           :: Rfrac             !! Parameter in the regression of diffuse 
                                                                    !! share on transmission
    REAL(r_std), DIMENSION(kjpindex)           :: Kfrac             !! Parameter in the regression of diffuse 
                                                                    !! share on transmission
    REAL(r_std), DIMENSION(kjpindex)           :: swdf              !! Sw diffuse radiation (W.m^{-2}) 
    REAL(r_std), DIMENSION(kjpindex)           :: swdr              !! Sw direct radiation (W.m^{-2}) 
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: PARscat           !! Scatter PAR @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: Clsun_iso         !! Isoprene dependance to light for sun leaves 
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: Clsh_iso          !! Isoprene dependance to light for shaded leaves
    !! for multilayer canopy model for iso flux
    REAL(r_std), DIMENSION(kjpindex,nlai+1)    :: PARscattab        !! Scatter PAR @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(nlai+1)             :: laitab            !! LAI per layer (m^2.m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nlai+1)    :: laisuntab         !! LAI of sun leaves per layer (m^2.m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nlai+1)    :: laishtab          !! LAI of shaded leaves per layer 
                                                                    !! (m^2.m^{-2})
    REAL(r_std)                                :: Clsun_iso_tab     !! Isoprene dependance to light 
                                                                    !! for sun leaves and per layer 
    REAL(r_std)                                :: Clsh_iso_tab      !! Isoprene dependance to light 
                                                                    !! for shaded leaves and per layer
    !!Leaf age
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: Eff_age_iso       !! Isoprene emission dependance to Leaf Age 
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: Eff_age_meth      !! Methanol emission dependance to Leaf Age 
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: Eff_age_VOC       !! Other VOC emission dependance to Leaf Age 
    !!BBG and Fertilizers for NOx soil emission
    REAL(r_std), DIMENSION(kjpindex)           :: veget_max_nowoody !! sum of veget_max for non-woody PFT
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: N_qt_WRICE_pft    !! N fertiliser on RICE 
                                                                    !! (kgN per year per grid cell)
    REAL(r_std), DIMENSION(kjpindex,nvm)       :: N_qt_OTHER_pft    !! N fertiliser on other veg 
                                                                    !! (kgN per year per grid cell)


    !! 0.5 Parameters values

    REAL(r_std), PARAMETER :: CT1 = 95000.0       !! Empirical coeffcient (see Guenther .et. al, 1995, eq(10)) (J.mol^{-1})
    REAL(r_std), PARAMETER :: CT2 = 230000.0      !! Empirical coefficient (see Guenther .et. al, 1995, eq(10)) (J.mol^{-1})
    REAL(r_std), PARAMETER :: TS = 303.0          !! Leaf temperature at standard condition
                                                  !! (see Guenther .et. al, 1995, eq(10)) (K)
    REAL(r_std), PARAMETER :: TM = 314.0          !! Leaf temperature (see Guenther .et. al, 1995, eq(10)) (K)

    REAL(r_std), PARAMETER :: alpha_ = 0.0027     !! Empirical coeffcient (see Guenther .et. al, 1995, eq(9)) (unitless)
    REAL(r_std), PARAMETER :: CL1 = 1.066         !! Empirical coeffcient (see Guenther .et. al, 1995, eq(9)) (unitless)
    REAL(r_std), PARAMETER :: beta = 0.09         !! Empirical coeffcient (see Guenther .et. al, 1995, eq(11)) (K^{-1})
    REAL(r_std), PARAMETER :: lai_threshold = 11. !! Lai threshold for the calculation of scattered radiation
                                                  !! based on Guenther .et. al (1995) (m^2.m^{-2})

!_ ================================================================================================================================

    !! 1. Canopy radiative transfer model

    !! Canopy radiative transfer model: takes into account light extinction through canopy
    !! First need to calculate diffuse and direct radiation
    !! Based on Andrew Friend radiative model (Global Ecology & Biogeography, 2001)
    !! And Spitters et al. (Agricultural and Forest Meteorology, 1986)

    IF ( control%ok_radcanopy ) THEN

       DO ji = 1, kjpindex
          IF (sinang(ji) .GT. zero) THEN
             day = julian_diff
             !! 1.1 Extra-terrestrial solar irradiance at a plan parallel to Earh's surface
             So(ji) = Sct*( un + 0.033*COS(360.*pi/180.*day/365.))*sinang(ji)
             !! 1.2 Atmospheric transmissivity
             Trans(ji) = swdown(ji)/So(ji)
             !! 1.3 Empirical calculation of fraction diffuse from transmission based on Spitters et al. (1986)
             Rfrac(ji) = 0.847 - 1.61*sinang(ji) + 1.04*(sinang(ji)**2.)
             Kfrac(ji) = (1.47 - Rfrac(ji)*Rfrac(ji))/1.66      
             IF (Trans(ji) .LE. 0.22) THEN
                Fdf(ji) = un
             ELSE IF (Trans(ji) .LE. 0.35) THEN
                Fdf(ji) = un - 6.4*((Trans(ji) - 0.22)**2.) 
             ELSE IF (Trans(ji) .LE. Kfrac(ji)) THEN
                Fdf(ji) = 1.47 - 1.66*Trans(ji)
             ELSE
                Fdf(ji) = Rfrac(ji)
             END IF
             !! 1.4 Direct and diffuse sw radiation in W.m^{-2}
             swdf(ji) = swdown(ji)*Fdf(ji)
             swdr(ji) = swdown(ji)*(un-Fdf(ji))
          ELSE
             swdf(ji) = zero
             swdr(ji) = zero
          END IF

          !! 1.5 PAR diffuse and direct in umol/m^2/s
          PARdf(ji) = swdf(ji) * W_to_mmol * RG_to_PAR
          PARdr(ji) = swdr(ji) * W_to_mmol * RG_to_PAR 
       END DO

       !! 1.6 Calculation of lai, parscat, parsh and parsun, laisun and laish !!?? define the terms
       !! Based on Guenther et al. (JGR, 1995) and Norman (1982)
       ! One-layer canopy model or multi-layer canopy model
       IF (control%ok_multilayer) THEN 


          ! Calculation PER LAYER
          DO jl = 1, nlai+1
             laitab(jl) = laimax*(EXP(lai_level_depth*REAL(jl-1,r_std)) - un)/(EXP(lai_level_depth*REAL(nlai,r_std)) - un)

             DO ji = 1, kjpindex
                IF (laitab(jl) .LE. lai_threshold) THEN
                   PARscattab(ji,jl) = 0.07*PARdr(ji)*(1.1 - 0.1*laitab(jl))*exp(-sinang(ji))
                ELSE
                   PARscattab(ji,jl) = zero
                END IF

                IF (sinang(ji) .NE. zero ) THEN
                   PARshtab(ji,jl) = PARdf(ji)*exp(-0.5*((laitab(jl))**0.7)) + PARscattab(ji,jl)
                   PARsuntab(ji,jl) = PARdr(ji)*COS(60.*pi/180.)/sinang(ji) + PARshtab(ji,jl)
                ELSE
                   PARshtab(ji,jl) = PARdf(ji)*exp(-0.5*(laitab(jl)**0.7)) + PARscattab(ji,jl)
                   PARsuntab(ji,jl) = zero 
                END IF
                IF (ABS(ACOS(sinang(ji))) .LT. pi/2. .AND. sinang(ji) .NE. zero) THEN 
                   ! calculation corrected for multi-layer canopy model from Friend (2001)
                   laisuntab(ji,jl) = laitab(jl)*exp(-0.5*laitab(jl)) 
                   laishtab(ji,jl) = laitab(jl) - laisuntab(ji,jl)
                ELSE
                   laisuntab(ji,jl) = zero
                   laishtab(ji,jl) = laitab(jl)
                END IF
             END DO
          END DO
       ELSE
          ! Calculation FOR one layer
          DO jv = 1, nvm
             DO ji = 1, kjpindex
                IF (lai(ji,jv) .LE. lai_threshold) THEN
                   PARscat(ji,jv) = 0.07*PARdr(ji)*(1.1 - 0.1*lai(ji,jv))*exp(-sinang(ji))
                ELSE
                   PARscat(ji,jv) = zero
                END IF

                IF (sinang(ji) .NE. zero ) THEN
                   PARsh(ji,jv) = PARdf(ji)*exp(-0.5*((lai(ji,jv))**0.7)) + PARscat(ji,jv)
                   PARsun(ji,jv) = PARdr(ji)*COS(60.*pi/180.)/sinang(ji) + PARsh(ji,jv)
                ELSE
                   PARsh(ji,jv) = PARdf(ji)*exp(-0.5*(lai(ji,jv)**0.7)) + PARscat(ji,jv)
                   PARsun(ji,jv) = zero 
                END IF
                IF (ABS(ACOS(sinang(ji))) .LT. pi/2. .AND. sinang(ji) .NE. zero) THEN 
                   ! calculation as in Lathiere (2005) = with correction removing lai in Guenther (1995)
                   laisun(ji,jv) = (un - exp(-0.5*lai(ji,jv)/(sinang(ji))))*sinang(ji)/COS(60.*pi/180.)
                   laish(ji,jv) = lai(ji,jv) - laisun(ji,jv)
                ELSE
                   laisun(ji,jv) = zero
                   laish(ji,jv) = lai(ji,jv)
                END IF
             END DO
          END DO
       ENDIF
    END IF


    !! 2. Calculation of non-PFT dependant parameters used for VOC emissions
    DO ji = 1, kjpindex ! (loop over # pixels)
       !! 2.1 Calculation of Tleaf (based on Lathiere, 2005)
       IF ((temp_sol(ji) - temp_air(ji)) .GT. 2) THEN
          tleaf(ji) = temp_air(ji) + 2
       ELSE IF ((temp_sol(ji) - temp_air(ji)) .LT. -2) THEN
          tleaf(ji) = temp_air(ji) - 2
       ELSE
          tleaf(ji) = temp_sol(ji)
       END IF

       !! 2.2 Isoprene emission dependency - with no PARsun/PARshaded partitioning - Guenther et al. (1995) and Lathiere (2005)
       !> @codeinc $$?? ecrire les equation en latex ? 
       exp_1(ji) = exp( (CT1 * ( tleaf(ji) - TS )) / (R*TS*tleaf(ji)) )
       exp_2(ji) = exp( (CT2 *( tleaf(ji) - TM )) / (R*TS*tleaf(ji)) )
       PAR(ji)   = swdown(ji) * W_to_mmol * RG_to_PAR        ! from W/m^2 to umol photons/m^2/s and half of sw for PAR
       Ct_iso(ji)    = exp_1(ji)/(un + exp_2(ji))            ! temperature dependance  
       Cl_iso(ji)    = alpha_*CL1*PAR(ji)/sqrt(un + (alpha_**2) * (PAR(ji)**2) ) ! light dependance
       !> @endcodeinc
       !! 2.3 Monoterpene emission dependency to Temperature
       !> @codeinc
       Ct_mono(ji) = exp(beta*(tleaf(ji) - TS))
       !> @endcodeinc
       !! 2.4 MBO biogenic emissions dependency, only from PFT7 and PFT4 for location of vegetation emitter
       ! but in fact MBO fluxes only in America (ponderosa and lodgepole pines only found in these areas)
       !> @codeinc
       Xvar(ji) = ((un/312.3) - (un/tleaf(ji)))/R
       !> @endcodeinc
       !! 2.4.1 temperature dependency
       !> @codeinc
       Ct_MBO(ji)    = (1.52*209000.0*exp(67000.0*Xvar(ji)))/(209000.0 - 67000.0*(un - exp(209000.0*Xvar(ji))))
       !> @endcodeinc
       !! 2.4.2 light dependency
       Cl_MBO(ji)    = (0.0011*1.44*PAR(ji))/(sqrt(un + (0.0011**2)*(PAR(ji)**2)))
       !! 2.5 NO biogenic emissions given in ngN/m^2/s, emission factor in ngN/m^2/s too
       !! calculation of temperature used for NO soil emissions
       t_no(ji) = ptnlev1(ji) - ZeroCelsius  !!temp must be in celsius to calculate no emissions
       !! 2.6 calculation of non-woody veget_max fraction
       IF (control%ok_cropsfertil_NOx) THEN
          veget_max_nowoody(ji) = zero
          DO jv = 1,nvm
             IF ( is_c3(jv) .OR. is_c4(jv) ) THEN
                veget_max_nowoody(ji) = veget_max_nowoody(ji) + veget_max(ji,jv)
             ENDIF
          ENDDO
       END IF
    END DO ! (loop over # pixels)

    !! 3. Calculation of PFT dependant parameters and
    ! Calculation of VOC emissions flux

    Eff_age_iso(:,:) = zero
    Eff_age_meth(:,:) = zero

    DO jv = 1, nvm
       DO ji = 1, kjpindex
          ! 6-Calculation of Leaf Age Function (Lathiere 2005)
          IF ( control%ok_leafage ) THEN
             DO jf = 1, nleafages
                !> @codeinc
                Eff_age_iso(ji,jv) = Eff_age_iso(ji,jv) + frac_age(ji,jv,jf)*iso_activity(jf)
                Eff_age_meth(ji,jv) = Eff_age_meth(ji,jv) + frac_age(ji,jv,jf)*methanol_activity(jf)
                !> @endcodeinc 
             END DO
             !> @codeinc
             Eff_age_VOC(ji,jv) = un
             !> @endcodeinc
          ELSE
             Eff_age_iso(ji,jv) = un
             Eff_age_meth(ji,jv) = un
             Eff_age_VOC(ji,jv) = un
          END IF
          !! 5. Calculation of foliar density
          IF ( sla(jv) .eq. zero ) THEN
             fol_dens(ji,jv) = zero
          ELSE
             ! 2 factor for conversion from gC to gDM
             fol_dens(ji,jv) = 2 * lai(ji,jv)/sla(jv)
          ENDIF
          !! 6. Calculation of VOC emissions from vegetation
          IF ( control%ok_radcanopy ) THEN
             ! if multi-layer canopy model
             IF (control%ok_multilayer) THEN 
                flx_iso(ji,jv) = zero
                laisun(ji,jv) = zero
                laish(ji,jv) = zero
                ! loop over the NLAI canopy layers
                DO jl = 1, nlai
                   IF ((laitab(jl) .LE. lai(ji,jv)).AND.(lai(ji,jv).NE.zero)) THEN
                      !sunlit vegetation 
                      Clsun_iso_tab   = alpha_*CL1*PARsuntab(ji,jl)/sqrt(un + (alpha_**2) * (PARsuntab(ji,jl)**2) )
                      ! shaded vegetation
                      Clsh_iso_tab    = alpha_*CL1*PARshtab(ji,jl)/sqrt(un + (alpha_**2) * (PARshtab(ji,jl)**2) ) 

                      !                  flx_iso(ji,jv) = flx_iso(ji,jv) + ((laisuntab(ji,jl+1)-laisuntab(ji,jl))*Clsun_iso_tab+ &
                      !                       & (laishtab(ji,jl+1)-laishtab(ji,jl))*Clsh_iso_tab)* &
                      !                       & fol_dens(ji,jv)*Ct_iso(ji)*em_factor_isoprene(jv)* &
                      !                       & Eff_age_iso(ji,jv)*1e-9/one_hour
                      !
                      flx_iso(ji,jv) = flx_iso(ji,jv) + ((laisuntab(ji,jl+1) - laisuntab(ji,jl))*Clsun_iso_tab+ &
                           & (laishtab(ji,jl+1) - laishtab(ji,jl))*Clsh_iso_tab)* &
                           & fol_dens(ji,jv)/lai(ji,jv)*Ct_iso(ji)*em_factor_isoprene(jv)* &
                           & Eff_age_iso(ji,jv)*1e-9/one_hour

                      laisun(ji,jv) = laisuntab(ji,jl)
                      laish(ji,jv)  = laishtab(ji,jl)
                   END IF
                END DO
                ! if mono-layer canopy model
             ELSE
                !sunlit vegetation 
                Clsun_iso(ji,jv)   = alpha_*CL1*PARsun(ji,jv)/sqrt(un + (alpha_**2) * (PARsun(ji,jv)**2) )
                ! shaded vegetation      
                Clsh_iso(ji,jv)    = alpha_*CL1*PARsh(ji,jv)/sqrt(un + (alpha_**2) * (PARsh(ji,jv)**2) )       
                IF (lai(ji,jv) .NE. zero) THEN
                   flx_iso(ji,jv) = (laisun(ji,jv)*Clsun_iso(ji,jv) + laish(ji,jv)*Clsh_iso(ji,jv))* &
                        & fol_dens(ji,jv)/lai(ji,jv)*Ct_iso(ji)*em_factor_isoprene(jv)* &
                        & Eff_age_iso(ji,jv)*1e-9/one_hour
                ELSE
                   ! 
                   flx_iso(ji,jv) = zero
                END IF
             END IF
             ! if no light extinction due to vegetation  
          ELSE
             !! Isoprene emissions - general equation
             !> @codeinc
             flx_iso(ji,jv) = fol_dens(ji,jv)*Ct_iso(ji)*Cl_iso(ji)*Eff_age_iso(ji,jv)*em_factor_isoprene(jv)*1e-9/one_hour
             !> @endcodeinc
          END IF
          !! 6.2 Calculation of monoterpene biogenic emissions
          !> @codeinc
          flx_mono(ji,jv) = fol_dens(ji,jv)*em_factor_monoterpene(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !> @endcodeinc
          !! 6.3 Calculation of ORVOC biogenic emissions
          !! Other Reactive Volatile Organic Compounds
          !> @codeinc
          flx_ORVOC(ji,jv) = fol_dens(ji,jv)*em_factor_ORVOC(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !> @endcodeinc
          !! 6.4 Calculation of OVOC biogenic emissions
          !! Other Volatile Organic Compounds
          flx_OVOC(ji,jv) = fol_dens(ji,jv)*em_factor_OVOC(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !! 6.5 Calculation of MBO biogenic emissions
          !! 2-Methyl-3-Buten-2-ol 
          IF(lalo(ji,1) .GE. 20. .AND. lalo(ji,2) .LE. -100) THEN
             flx_MBO(ji,jv) = fol_dens(ji,jv)*em_factor_MBO(jv)*Ct_MBO(ji)*Cl_MBO(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          ELSE
             flx_MBO(ji,jv) = zero
          END IF
          !! 6.6 Calculation of methanol biogenic emissions
          flx_methanol(ji,jv) = fol_dens(ji,jv)*em_factor_methanol(jv)*Ct_mono(ji)*Eff_age_meth(ji,jv)*1e-9/one_hour
          !! 6.7 Calculation of acetone biogenic emissions
          flx_acetone(ji,jv) = fol_dens(ji,jv)*em_factor_acetone(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !! 6.8 Calculation of acetaldehyde biogenic emissions
          flx_acetal(ji,jv) = fol_dens(ji,jv)*em_factor_acetal(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !! 6.9 Calculation of formaldehyde biogenic emissions
          flx_formal(ji,jv) = fol_dens(ji,jv)*em_factor_formal(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !! 6.10 Calculation of acetic acid biogenic emissions
          flx_acetic(ji,jv) = fol_dens(ji,jv)*em_factor_acetic(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
          !! 6.11 Calculation of formic acid biogenic emissions 
          flx_formic(ji,jv) = fol_dens(ji,jv)*em_factor_formic(jv)*Ct_mono(ji)*Eff_age_VOC(ji,jv)*1e-9/one_hour
       END DO

    END DO


    !! 7. Calculation of NOx emissions from soils
    ! Based on Yienger & Levy (1995) and Lathiere (2005, chapter 3)
    DO ji = 1, kjpindex
       !! 7.1 Precipitation-related pulse function
       IF (control%ok_pulse_NOx) THEN
          ! if we are during a period where pulses are not allowed
          IF (ok_siesta(ji)) THEN
             ! if this period is not over 
             IF (FLOOR(siestaday(ji)) .LE. siestalim(ji)) THEN
                siestaday(ji) = siestaday(ji) + (dtradia/one_day)
                ! if this period is over
             ELSE
                ok_siesta(ji) = .FALSE.
                siestaday(ji) = zero
             END IF
          END IF
          ! if we are during a period where pulses are allowed
          IF ((.NOT. ok_siesta(ji)) .AND. (.NOT. allow_pulse(ji))) THEN
             IF (humrel(ji,1) .LT. 0.15) THEN
                ! if precip exceeds 1 mm/day over one time step => a pulse occurs
                IF(precip_rain(ji)/nbre_precip .GE. un/(one_day/dtradia)) THEN
                   ! if precip is up to 5 mm/day => pulse length is 3 days
                   IF (precip_rain(ji)/nbre_precip .LT. 5./(one_day/dtradia)) THEN
                      pulselim(ji) = 3.
                      ! if precip is up to 15 mm/day => pulse length is 7 days
                   ELSE IF (precip_rain(ji)/nbre_precip .LT. 15./(one_day/dtradia)) THEN
                      pulselim(ji) = 7.
                      ! if precip is upper than 15 mm/day => pulse length is 14 days
                   ELSE IF (precip_rain(ji)/nbre_precip .GE. 15./(one_day/dtradia)) THEN
                      pulselim(ji) = 14.
                   END IF
                   allow_pulse(ji)=.TRUE.
                   pulseday(ji) = un
                END IF
             END IF
          END IF
          ! if we were during a pulse period
          IF (allow_pulse(ji)) THEN
             ! if we are still during the pulse period
             ! 16/06/2010 We assume a (pulselim-1) days for the pulse length (NVui+Jlath)
             IF(FLOOR(pulseday(ji)) .LT. pulselim(ji)) THEN
                ! calculation of the pulse function
                IF (pulselim(ji).EQ.3) THEN
                   pulse(ji) = 11.19*exp(-0.805*pulseday(ji))
                ELSE IF (pulselim(ji).EQ.7) THEN
                   pulse(ji) = 14.68*exp(-0.384*pulseday(ji))
                ELSE IF (pulselim(ji).EQ.14) THEN 
                   pulse(ji) = 18.46*exp(-0.208*pulseday(ji))
                END IF
                pulseday(ji) = pulseday(ji) + (dtradia/one_day)
                ! if the pulse period is over
             ELSE
                ! pulse function is set to 1 
                pulse(ji) = un
                allow_pulse(ji) = .FALSE.
                siestaday(ji) = un
                siestalim(ji) = pulselim(ji)
                ok_siesta(ji) = .TRUE. 
             END IF
          END IF
          ! no precipitation-related pulse function
       ELSE
          pulse(ji) = un
       END IF
    END DO

    !! 7.2 Calculation of NO basal emissions including pulse effect
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          !Tropical forests
          IF ( is_tropical(jv) .AND. is_evergreen(jv) ) THEN
             ! Wet soils
             IF (humrel(ji,1) .GT. 0.3) THEN
                flx_no_soil(ji,jv) = 2.6*pulse(ji)
                ! Dry soils
             ELSE
                flx_no_soil(ji,jv) = 8.6*pulse(ji)
             END IF
             !Else If agricultural lands OR Wet soils
          ELSE IF ( ( .NOT.(natural(jv)) ) .OR. ( humrel(ji,1) .GT. 0.3 ) ) THEN
             ! Calculation of NO emissions depending of Temperature
             IF (t_no(ji) .LT. zero) THEN
                flx_no_soil(ji,jv) = zero
             ELSE IF (t_no(ji) .LE. 10.) THEN
                flx_no_soil(ji,jv) = 0.28*em_factor_no_wet(jv)*t_no(ji)*pulse(ji)
             ELSE IF (t_no(ji) .LE. 30.) THEN
                flx_no_soil(ji,jv) = em_factor_no_wet(jv)*exp(0.103*t_no(ji))*pulse(ji)
             ELSE
                flx_no_soil(ji,jv) = 21.97*em_factor_no_wet(jv)*pulse(ji)
             END IF
             !Else if Temp negative
          ELSE IF (t_no(ji) .LT. zero) THEN
             flx_no_soil(ji,jv) = zero
             !Else if Temp <= 30
          ELSE IF (t_no(ji) .LE. 30.) THEN
             flx_no_soil(ji,jv) = (em_factor_no_dry(jv)*t_no(ji))/30.*pulse(ji)
          ELSE
             flx_no_soil(ji,jv) = em_factor_no_dry(jv)*pulse(ji)
          END IF

          !! 7.3 IF ACTIVATED (ok_bbgfertil_NOx = TRUE) calculation of NOx soil emission increase due to biomass burning
          ! Calculation of Biomass-Burning-induced NOx emissions (Lathiere, 2005)
          ! => NOx emissions 3-fold increase
          IF (control%ok_bbgfertil_NOx) THEN
             IF ( natural(jv) ) THEN
                ! North Tropical zone from May to June
                IF ((lalo(ji,1) .LE. 30. .AND. lalo(ji,1) .GE. zero).AND. &
                     (day .GE. 121. .AND. day .LE. 181).AND.(flx_co2_bbg_year(ji) .GT. 0.1)) THEN
                   flx_no_soil(ji,jv) = flx_no_soil(ji,jv)*3.
                   ! South Tropical zone from November to December
                ELSE IF ((lalo(ji,1) .GE. -30. .AND. lalo(ji,1) .LT. zero).AND.(day .GE. 305.).AND. & 
                        (flx_co2_bbg_year(ji) .GT. 0.1)) THEN
                   flx_no_soil(ji,jv) = flx_no_soil(ji,jv)*3.
                END IF
             END IF
          END IF

          !! 7.4 IF ACTIVATED (ok_cropsfertil_NOx = TRUE) calculation of NOx soil emission increase due to fertilizer use 
          ! Calculation of N-fertiliser-induced NOx emissions
          flx_fertil_no(ji,jv) = zero
          IF (control%ok_cropsfertil_NOx) THEN
             IF (veget_max_nowoody(ji) .NE. zero) THEN
                ! Non-agricultural lands
                IF ( natural(jv) .AND. .NOT.(is_c3(jv)) .AND. .NOT.(is_c4(jv)) ) THEN
                   N_qt_WRICE_pft(ji,jv) = zero
                   N_qt_OTHER_pft(ji,jv) = zero
                END IF
                ! Grasslands or Croplands
                IF ( is_c3(jv) .OR. is_c4 (jv) ) THEN
                   N_qt_WRICE_pft(ji,jv) = N_qt_WRICE_year(ji)*veget_max(ji,jv)/veget_max_nowoody(ji)
                   N_qt_OTHER_pft(ji,jv) = N_qt_OTHER_year(ji)*veget_max(ji,jv)/veget_max_nowoody(ji)
                ENDIF
             ELSE
                N_qt_WRICE_pft(ji,jv) = zero
                N_qt_OTHER_pft(ji,jv) = zero
             END IF

             ! North temperate regions from May to August
             ! OR South Temperate regions from November to February
             IF (((lalo(ji,1) .GT. 30.) .AND. (day .GE. 121. .AND. day .LE. 243.).AND.(veget_max(ji,jv) .NE. zero)) .OR. & 
             &  ((lalo(ji,1) .LT. -30.) .AND. (day .GE. 305. .OR. day .LE. 59.) .AND.(veget_max(ji,jv) .NE. zero))) THEN
                ! 1e12 for conversion from kg to Ng
                ! 1/(365/12*24*60*60*4) for conversion from year to seconds corrected for 4 months of emissions
                flx_fertil_no(ji,jv) = (N_qt_WRICE_pft(ji,jv)*(1/30.)+N_qt_OTHER_pft(ji,jv))*(2.5/100.) &
                     & *1e12/(365*24*60*60*4/12)/(area2(ji)*veget_max(ji,jv))
                ! OR Tropical regions all the year
             ELSE IF ((lalo(ji,1) .GE. -30.).AND.(lalo(ji,1) .LE. 30.).AND.(veget_max(ji,jv) .NE. zero)) THEN
                flx_fertil_no(ji,jv) = (N_qt_WRICE_pft(ji,jv)*(1/30.)+N_qt_OTHER_pft(ji,jv))*(2.5/100.) &
                     & *1e12/(365*24*60*60)/(area2(ji)*veget_max(ji,jv))
             END IF
             flx_no_soil(ji,jv) = flx_no_soil(ji,jv) + flx_fertil_no(ji,jv)
          END IF

          !! 7.5 Calculation of net NO flux above soil accounting for surface deposition, 
          !! based on the Canopy Reduction Factor (CRF), calculated using Leaf Area and Stomatal Area
          !kc=cuticle absorptivity = 0.24m^2/m^2
          !ks=stomatal absorptivity = 8.75m^2/m^2
          !Larch=Larcher SAI/LAI ratio given in Larcher 1991
          !> @codeinc
          CRF(ji,jv) = (exp(-8.75*Larch(jv)*lai(ji,jv)) + exp(-0.24*lai(ji,jv)))/2.
          flx_no(ji,jv) = flx_no_soil(ji,jv)*CRF(ji,jv)
          !> @endcodeinc
       END DO
    END DO
    IF (long_print) WRITE(numout,*) 'OK diffuco inca'


  END SUBROUTINE diffuco_inca

!! ================================================================================================================================
!! SUBROUTINE   : diffuco_inca_read
!!
!>\BRIEF         This subroutine will read and interpolate the data set of CO2 emissions from biomass burning.\n
!!
!! DESCRIPTION  :  This subroutine will interpolate the 1 x 1 deg based data set of CO2 emissions from biomass burning
!!                 expresssed in kgC/m^2 per year (=> density = TRUE), N fertiliser amount expressed in kgN per year and
!!                 per grid cell (=> density = FALSE) to the resolution of the model.\ n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: data_year
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

  SUBROUTINE diffuco_inca_read (nbpt, lalo, neighbours, resolution, data_year, filename, filename2, fieldname, density)

    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)       :: nbpt                !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)          :: lalo(nbpt,2)        !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)       :: neighbours(nbpt,8)  !! Vector of neighbours for each grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std), INTENT(in)          :: resolution(nbpt,2)  !! The size in m of each grid-box in X and Y 
    CHARACTER(LEN=80), INTENT(in)    :: filename2
    CHARACTER(LEN=80), INTENT(in)    :: fieldname
    LOGICAL, INTENT(in)              :: density

    !! 0.2 Output variables

    REAL(r_std), INTENT(out)          :: data_year(:)       !! data per year (climatology) per unit area 

    !! 0.3 Modified variables

    CHARACTER(LEN=80), INTENT(inout)  :: filename


    !! 0.4 Local variables 

    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, ilf, lastjp, nbexp
    REAL(r_std) :: lev(1), date, dt, coslat
    REAL(r_std) :: lon_up, lon_low, lat_up, lat_low
    INTEGER(i_std) :: itau(1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_rel, lon_rel, lat_ful, lon_ful
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: data_year_file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: loup_rel, lolow_rel, laup_rel, lalow_rel
    REAL(r_std) :: area2                                     !! total area of the final grid box 
    REAL(r_std) :: data_grid                                 !! emissions for the final grid point before being 
                                                             !! ponderated by area2
    REAL(r_std) :: data_global                               !! global emissions to check
    REAL(r_std) :: ax, ay, sgn, surp, ax_file, ay_file
    REAL(r_std) :: lonrel, louprel, lolowrel

!_ ================================================================================================================================

    !
    CALL getin_p(filename2,filename)
    !
    CALL flininfo(filename,iml, jml, lml, tml, fid)
    !
    !
    ALLOCATE (lat_rel(iml,jml))
    ALLOCATE (lon_rel(iml,jml))
    ALLOCATE (laup_rel(iml,jml))
    ALLOCATE (loup_rel(iml,jml))
    ALLOCATE (lalow_rel(iml,jml))
    ALLOCATE (lolow_rel(iml,jml))
    ALLOCATE (lat_ful(iml+2,jml+2))
    ALLOCATE (lon_ful(iml+2,jml+2))
    ALLOCATE (data_year_file(iml,jml))
    !
    CALL flinopen(filename, .FALSE., iml, jml, lml, lon_rel, lat_rel, lev, tml, itau, date, dt, fid)
    !
    CALL flinget(fid, fieldname, iml, jml, lml, tml, 1, 1, data_year_file)
    !
    CALL flinclo(fid)
    !
    !
    WRITE(*,*) 'lon_rel : ', MAXVAL(lon_rel), MINVAL(lon_rel)
    WRITE(*,*) 'lat_rel : ', MAXVAL(lat_rel), MINVAL(lat_rel)
    WRITE(*,*) 'data_year_file min,max: ', MINVAL(data_year_file, MASK=data_year_file .GT. 0), &
         &                      MAXVAL(data_year_file, MASK=data_year_file .LT. undef_sechiba)
    !
    nbexp = 0
    !
    !    Duplicate the border assuming we have a global grid going from west to east
    !
    lon_ful(2:iml+1,2:jml+1) = lon_rel(1:iml,1:jml)
    lat_ful(2:iml+1,2:jml+1) = lat_rel(1:iml,1:jml)
    !
    IF ( lon_rel(iml,1) .LT. lon_ful(2,2)) THEN
       lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)
       lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
    ELSE
       lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)-360
       lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
    ENDIF

    IF ( lon_rel(1,1) .GT. lon_ful(iml+1,2)) THEN
       lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)
       lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
    ELSE
       lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)+360
       lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
    ENDIF
    !
    sgn = lat_rel(1,1)/ABS(lat_rel(1,1))
    lat_ful(2:iml+1,1) = sgn*180 - lat_rel(1:iml,1)
    sgn = lat_rel(1,jml)/ABS(lat_rel(1,jml))
    lat_ful(2:iml+1,jml+2) = sgn*180 - lat_rel(1:iml,jml)
    lat_ful(1,1) = lat_ful(iml+1,1)
    lat_ful(iml+2,1) = lat_ful(2,1)
    lat_ful(1,jml+2) = lat_ful(iml+1,jml+2)
    lat_ful(iml+2,jml+2) = lat_ful(2,jml+2)
    !
    ! Add the longitude lines to the top and bottom
    !
    lon_ful(:,1) = lon_ful(:,2)
    lon_ful(:,jml+2) = lon_ful(:,jml+1)
    !
    !  Get the upper and lower limits of each grid box
    !
    DO ip = 1,iml
       DO jp = 1,jml
          loup_rel(ip,jp) =MAX(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), 0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
          lolow_rel(ip,jp) =MIN(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), 0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
          laup_rel(ip,jp) =MAX(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), 0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
          lalow_rel(ip,jp) =MIN(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), 0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
       ENDDO
    ENDDO
    !
    !   Now we take each grid point and find out which values from the forcing we need to average
    !
    DO ib =1, nbpt
       !
       !  We find the 4 limits of the grid-box. As we transform the resolution of the model
       !  into longitudes and latitudes we do not have the problem of periodicity.
       ! coslat is a help variable here !
       !
       coslat = MAX(COS(lalo(ib,1) * pi/180. ), mincos )*pi/180. * R_Earth
       !
       lon_up = lalo(ib,2)+ resolution(ib,1)/(2.0*coslat)
       lon_low =lalo(ib,2) - resolution(ib,1)/(2.0*coslat)
       !
       coslat = pi/180. * R_Earth
       !
       lat_up =lalo(ib,1)+resolution(ib,2)/(2.0*coslat)
       lat_low =lalo(ib,1)-resolution(ib,2)/(2.0*coslat)
       !
       !
       !  Find the grid boxes from the data that go into the model's boxes
       !  We still work as if we had a regular grid ! Well it needs to be localy regular so
       !  so that the longitude at the latitude of the last found point is close to the one of the next point.
       !
       lastjp = 1
       data_grid = zero
       area2 = zero

       DO ip = 1,iml
          !
          !  Either the center of the data grid point is in the interval of the model grid or
          !  the East and West limits of the data grid point are on either sides of the border of
          !  the data grid.
          !
          !
          !  We find the 4 limits of the grid-box. As we transform the resolution of the model
          !  into longitudes and latitudes we do not have the problem of periodicity.
          ! coslat is a help variable here !
          !
          !
          !  To do that correctly we have to check if the grid box sits on the date-line.
          !
          IF ( lon_low < -180.0 ) THEN
             lonrel = MOD( lon_rel(ip,lastjp) - 360.0, 360.0)
             lolowrel = MOD( lolow_rel(ip,lastjp) - 360.0, 360.0)
             louprel = MOD( loup_rel(ip,lastjp) - 360.0, 360.0)
             !
          ELSE IF ( lon_up > 180.0 ) THEN
             lonrel = MOD( 360. - lon_rel(ip,lastjp), 360.0)
             lolowrel = MOD( 360. - lolow_rel(ip,lastjp), 360.0)
             louprel = MOD( 360. - loup_rel(ip,lastjp), 360.0)
          ELSE
             lonrel = lon_rel(ip,lastjp)
             lolowrel = lolow_rel(ip,lastjp)
             louprel = loup_rel(ip,lastjp)
          ENDIF
          !
          !
          !
          IF ( lonrel > lon_low .AND. lonrel < lon_up .OR. &
               & lolowrel < lon_low .AND.  louprel > lon_low .OR. &
               & lolowrel < lon_up  .AND.  louprel > lon_up ) THEN
             ! 
             DO jp = 1, jml
                !
                ! Now that we have the longitude let us find the latitude
                !
                IF ( lat_rel(ip,jp) > lat_low .AND. lat_rel(ip,jp) < lat_up .OR. &
                     & lalow_rel(ip,jp) < lat_low .AND. laup_rel(ip,jp) > lat_low .OR.&
                     & lalow_rel(ip,jp) < lat_up .AND. laup_rel(ip,jp) > lat_up) THEN
                   !
                   lastjp = jp
                   !
                   ! Mising values in the file are assumed to be 1e20
                   !
                   IF ( lon_low < -180.0 ) THEN
                      lolowrel = MOD( lolow_rel(ip,jp) - 360.0, 360.0)
                      louprel = MOD( loup_rel(ip,jp) - 360.0, 360.0)
                      !
                   ELSE IF ( lon_up > 180.0 ) THEN
                      lolowrel = MOD( 360. - lolow_rel(ip,jp), 360.0)
                      louprel = MOD( 360. - loup_rel(ip,jp), 360.0)
                   ELSE
                      lolowrel = lolow_rel(ip,jp)
                      louprel = loup_rel(ip,jp)
                   ENDIF
                   !
                   ! Get the area of the fine grid in the model grid
                   !
                   coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )
                   ax = (MIN(lon_up,louprel)-MAX(lon_low, lolowrel))*pi/180. * R_Earth * coslat
                   ay = (MIN(lat_up, laup_rel(ip,jp))-MAX(lat_low,lalow_rel(ip,jp)))*pi/180. * R_Earth
                   !!to calculate the surface of the initial grid box, ie data one
                   ax_file = (louprel-lolowrel)*pi/180. * R_Earth * coslat
                   ay_file = (laup_rel(ip,jp)-lalow_rel(ip,jp))*pi/180. * R_Earth
                   !
                   IF (data_year_file(ip,jp) .LT. undef_sechiba-un) THEN
                      IF(density) THEN
                         data_grid = data_grid + ax*ay*data_year_file(ip,jp)
                      ELSE
                         data_grid = data_grid + ax*ay*(data_year_file(ip,jp)/(ax_file*ay_file))
                      ENDIF
                      area2 = area2 + ax*ay
                   ENDIF
                   !
                ENDIF
                !
                !
             ENDDO
             !
          ENDIF
          !
       ENDDO
       !
       ! Put the total data_year areas in the output variables
       !
       IF(density) THEN
          data_year(ib) = data_grid/area2 
       ELSE
          data_year(ib) = data_grid
       ENDIF


       IF ( data_year(ib) < 0 ) THEN
          WRITE(*,*) 'We have a problem here : ', data_year(ib)
          WRITE(*,*) 'resolution :', resolution(ib,1), resolution(ib,2)
          WRITE(*,*) area2
          STOP
       ENDIF
       !
       !!Accumulates emissions (in Kg/year) per grid box to calculate global emissions
       data_global = data_global + data_grid
    ENDDO
    !
    WRITE(*,*) 'RESULT data_year : ', MINVAL(data_year), MAXVAL(data_year)
    WRITE(*,*) 'RESULT Global data emissions : ', data_global 
    !

  END SUBROUTINE diffuco_inca_read


END MODULE diffuco
