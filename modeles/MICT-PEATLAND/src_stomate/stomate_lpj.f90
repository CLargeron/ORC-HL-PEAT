! ================================================================================================================================
! MODULE       : stomate_lpj
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Main entry point for daily processes in STOMATE and LPJ (phenology, 
!! allocation, npp_calc, kill, turn, light, establish, crown, cover, lcchange)
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_lpj.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_lpj

  ! modules used:

  USE ioipsl
  USE grid
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE lpj_constraints
  USE lpj_pftinout
  USE lpj_kill
  USE lpj_crown
  USE lpj_fire
  USE lpj_gap
  USE lpj_light
  USE lpj_establish
  USE lpj_cover
  USE stomate_prescribe
  USE stomate_phenology
  USE stomate_alloc
  USE stomate_npp
  USE stomate_turnover
  USE stomate_litter
  USE stomate_soilcarbon
  USE stomate_vmax
  USE stomate_assimtemp
  USE stomate_lcchange
  !  USE Write_Field_p
  
  !Chloe peat :
  USE stomate_wet_ch4_pt_ter_peat
  USE stomate_wet_ch4_pt_ter_0

   IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC StomateLpj,StomateLpj_clear

  LOGICAL, SAVE                         :: firstcall = .TRUE.             !! first call

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : StomateLpj_clear
!!
!>\BRIEF        Re-initialisation of variable
!!
!! DESCRIPTION  : This subroutine reinitializes variables. To be used if we want to relaunch 
!! ORCHIDEE but the routine is not used in current version.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE StomateLpj_clear

    CALL prescribe_clear
    CALL phenology_clear
    CALL npp_calc_clear
    CALL turn_clear
    CALL soilcarbon_clear
    CALL constraints_clear
    CALL establish_clear
    CALL fire_clear
    CALL gap_clear
    CALL light_clear
    CALL pftinout_clear
    CALL alloc_clear
    !Chloe modele peat
    CALL ch4_wet_flux_density_clear_peat
 CALL ch4_wet_flux_density_clear_0
  END SUBROUTINE StomateLpj_clear


!! ================================================================================================================================
!! SUBROUTINE   : StomateLPJ
!!
!>\BRIEF        Main entry point for daily processes in STOMATE and LPJ, structures the call sequence 
!!              to the different processes such as dispersion, establishment, competition and mortality of PFT's.
!! 
!! DESCRIPTION  : This routine is the main entry point to all processes calculated on a 
!! daily time step. Is mainly devoted to call the different STOMATE and LPJ routines 
!! depending of the ok_dgvm (is dynamic veg used) and lpj_constant_mortality (is background mortality used).
!! It also prepares the cumulative 
!! fluxes or pools (e.g TOTAL_M TOTAL_BM_LITTER etc...)
!!
!! This routine makes frequent use of "weekly", "monthly" and "long term" variables. Quotion is used because
!! by default "weekly" denotes 7 days, by default "monthly" denotes 20 days and by default "Long term" denotes
!! 3 years. dtslow refers to 24 hours (1 day).
!!
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): All variables related to stomate and required for LPJ dynamic vegetation mode.
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudré, J. Ogeé, J. Polcher, P. Friedlingstein, P. Ciais, S. Sitch, 
!! and I. C. Prentice. 2005. A dynamic global vegetation model for studies of the coupled atmosphere-biosphere 
!! system. Global Biogeochemical Cycles 19:GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, I. C. Prentice, A. Arneth, A. Bondeau, W. Cramer, J. O. Kaplan, S. Levis, W. Lucht, 
!! M. T. Sykes, K. Thonicke, and S. Venevsky. 2003. Evaluation of ecosystem dynamics, plant geography and 
!! terrestrial carbon cycling in the LPJ dynamic global vegetation model. Global Change Biology 9:161-185.
!!
!! FLOWCHART    : Update with existing flowchart from N Viovy (Jan 19, 2012)
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE StomateLpj (npts,dt_days, EndOfYear, EndOfMonth, &
       neighbours, resolution, &
       clay, herbivores, &
       tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
       litterhum_daily, soilhum_daily, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       gdd0_lastyear, precip_lastyear, &
       moiavail_month, moiavail_week, tlong_ref, t2m_month, t2m_week, &
       tsoil_month, soilhum_month, &
       gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
       turnover_longterm, gpp_daily, time_lowgpp, &
       time_hum_min, maxfpc_lastyear, resp_maint_part, &
       PFTpresent, age, fireindex, firelitter, &
       leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
       senescence, when_growthinit, &
       litterpart, litter, dead_leaves, carbon, black_carbon, lignin_struc, &
       veget_max, npp_longterm, lm_lastyearmax, veget_lastlight, &
       everywhere, need_adjacent, RIP_time, &
       lai, rprof,npp_daily, turnover_daily, turnover_time,&
       control_moist, control_temp, soilcarbon_input, &
       co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, &
       height, deadleaf_cover, vcmax, vjmax, &
       t_photo_min, t_photo_opt, t_photo_max,bm_to_litter, &
       prod10,prod100,flux10, flux100, veget_max_new, &
       convflux,cflux_prod10,cflux_prod100, harvest_above, carb_mass_total, lcchange, &
       fpc_max, & !) 
       !Chloe :
       ch4_flux_density_tot_peat,ch4_flux_density_dif_peat,ch4_flux_density_bub_peat,ch4_flux_density_pla_peat,&
       ch4_flux_density_goxid_peat,wpro_conso_carb,uold2_peat,tsurf_year,&
       !Chloe provisoire : 
        ch4_flux_density_tot_0, ch4_flux_density_dif_0, ch4_flux_density_bub_0,& 
        ch4_flux_density_pla_0)
    
  !! 0. Variable and parameter declaration

    !! 0.1 input

    INTEGER(i_std),PARAMETER                         :: n = 171
!    INTEGER(i_std), INTENT(in)                                 :: nvert                 !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL(r_std), INTENT(in)                                    :: dt_days              !! Time step of Stomate (days)
    INTEGER(i_std), DIMENSION(npts,8), INTENT(in)              :: neighbours           !! Indices of the 8 neighbours of each grid 
                                                                                       !! point [1=N, 2=NE, 3=E, 4=SE,
                                                                                       !!  5=S, 6=SW, 7=W, 8=NW] 
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                 :: resolution           !! Resolution at each grid point (m)  
                                                                                       !! [1=E-W, 2=N-S] 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: clay                 !! Clay fraction (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! Time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_daily          !! Daily surface temperatures (K)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: tsoil_daily          !! Daily soil temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_daily            !! Daily 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_min_daily        !! Daily minimum 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: litterhum_daily      !! Daily litter humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: soilhum_daily        !! Daily soil humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! Last year's maximum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! Last year's minimum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: gdd0_lastyear        !! Last year's GDD0 (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: precip_lastyear      !! Lastyear's precipitation 
                                                                                       !! @tex $(mm year^{-1})$ @endtex
										       !! to determine if establishment possible
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
                                                                                       !! unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "Weekly" moisture availability 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tlong_ref            !! "Long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "Monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "Weekly" 2-meter temperatures (K)
!Chloe peat :
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_year           !! annual surface temperatures (K)
!Chloe-
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: tsoil_month          !! "Monthly" soil temperatures (K)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: soilhum_month        !! "Monthly" soil humidity
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gdd_m5_dormance      !! Growing degree days (K), threshold -5 deg 
                                                                                       !! C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: gdd_midwinter        !! Growing degree days (K), since midwinter 
                                                                                       !! (for phenology) - this is written to the history files 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: ncd_dormance         !! Number of chilling days (days), since 
                                                                                       !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: ngd_minus5           !! Number of growing days (days), threshold 
                                                                                       !! -5 deg C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: turnover_longterm    !! "Long term" turnover rate  
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gpp_daily            !! Daily gross primary productivity  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
!    REAL(r_std), DIMENSION(npts,nstm), INTENT(in)               :: wtsoil_daily         !! daily mean of wt_soil2
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: time_lowgpp          !! Duration of dormance (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: time_hum_min         !! Time elapsed since strongest moisture 
                                                                                       !! availability (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxfpc_lastyear      !! Last year's maximum foliage projected
                                                                                       !! coverage for each natural PFT,
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: resp_maint_part      !! Maintenance respiration of different 
                                                                                       !! plant parts  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: fpc_max              !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    LOGICAL, INTENT(in)                                        :: lcchange             !! Land cover change flag
    LOGICAL, INTENT(in)                                        :: EndOfYear            !! Flag set on the last day of the year used 
                                                                                       !! to update "yearly variables". This 
                                                                                       !! variable must be .TRUE. once a year
    LOGICAL, INTENT(in)                                        :: EndOfMonth           !! Flag set at end of each month to update 
                                                                                       !! monthly variable 

  !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: npp_daily            !! Net primary productivity 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out)       :: turnover_daily       !! Turnover rates 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_to_bm            !! CO2 taken up from atmosphere when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                                       !! fire (living and dead biomass)  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero          !! Heterotrophic respiration
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_maint           !! Maintenance respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_growth          !! Growth respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: deadleaf_cover       !! Fraction of soil covered by dead leaves 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax                !! Maximum rate of carboxylation 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vjmax                !! Maximum rate of RUbp regeneration  
                                                                                       !! @tex $(\mumol CO_2 m^{-2}s^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: t_photo_min          !! Minimum temperature for photosynthesis 
                                                                                       !! (C) This is not K because it is a scalar.
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: t_photo_opt          !! Optimum temperature for photosynthesis 
                                                                                       !! (C) This is not K because it is a scalar.
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: t_photo_max          !! Maximum temperature for photosynthesis 
                                                                                       !! (C) This is not K because it is a scalar.
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out)       :: bm_to_litter         !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
!Chloe provisoire
  REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_0
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_0
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_0
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_0


!Chloe peat :
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_peat
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_peat
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_peat
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_peat
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_goxid_peat
     REAL(r_std), DIMENSION(npts), INTENT(in)            :: wpro_conso_carb
     REAL(r_std), DIMENSION(npts,n), INTENT(in)            :: uold2_peat
!Chloe--


    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: height               !! Height of vegetation (m) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_moist        !! Moisture control of heterotrophic 
                                                                                       !! respiration (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_temp         !! Temperature control of heterotrophic 
                                                                                       !! respiration, above and below 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: soilcarbon_input     !! Quantity of carbon going into carbon 
                                                                                       !! pools from litter decomposition  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lai                  !! Leaf area index OF AN INDIVIDUAL PLANT,
										       !! where a PFT contains n indentical plants
										       !! i.e., using the mean individual approach 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: rprof                !! Prescribed root depth (m) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: PFTpresent           !! Tab indicating which PFTs are present in 
                                                                                       !! each pixel 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! Age (years)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: fireindex            !! Probability of fire (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: firelitter           !! Longer term litter above the ground that 
                                                                                       !! can be burned, @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age             !! Leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac            !! Fraction of leaves in leaf age class, 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)     :: biomass              !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: ind                  !! Density of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: adapted              !! Adaptation of PFT (killed if too cold) 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: regenerate           !! "Fitness": Winter sufficiently cold for 
                                                                                       !! PFT regeneration ? (0 to 1, unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: senescence           !! Flag for setting senescence stage (only 
                                                                                       !! for deciduous trees) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: when_growthinit      !! How many days ago was the beginning of 
                                                                                       !! the growing season (days) 
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: litterpart           !! Fraction of litter above the ground 
                                                                                       !! belonging to different PFTs
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs), INTENT(inout):: litter               !! Metabolic and structural litter, above 
                                                                                       !! and below ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: dead_leaves          !! Dead leaves on ground, per PFT, metabolic 
                                                                                       !! and structural,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon               !! Carbon pool: active, slow, or passive, 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: black_carbon         !! Black carbon on the ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)      :: lignin_struc         !! Ratio of Lignin/Carbon in structural 
                                                                                       !! litter, above and below ground,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_max            !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: npp_longterm         !! "Long term" mean yearly primary 
                                                                                       !! productivity 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lm_lastyearmax       !! Last year's maximum leaf mass, for each 
                                                                                       !! PFT @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_lastlight      !! Vegetation fractions (on ground) after 
                                                                                       !! last light competition  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: everywhere           !! Is the PFT everywhere in the grid box or 
                                                                                       !! very localized (after its introduction) 
                                                                                       !! (unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: need_adjacent        !! In order for this PFT to be introduced, 
                                                                                       !! does it have to be present in an 
                                                                                       !! adjacent grid box? 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: RIP_time             !! How much time ago was the PFT eliminated 
                                                                                       !! for the last time (y) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! Turnover_time of leaves for grasses 
                                                                                       !! (days)
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)             :: veget_max_new        !! New "maximal" coverage fraction of a PFT 
                                                                                       !! (LAI -> infinity) (unitless) 
    REAL(r_std),DIMENSION(npts,0:10), INTENT(inout)            :: prod10               !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,0:100), INTENT(inout)           :: prod100              !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,10), INTENT(inout)              :: flux10               !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,100), INTENT(inout)             :: flux100              !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: convflux             !! Release during first year following land 
                                                                                       !! cover change @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod10         !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod100        !! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: harvest_above        !! Harvest above ground biomass for 
                                                                                       !! agriculture @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: carb_mass_total      !! Carbon Mass total (soil, litter, veg) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_bm_to_litter    !! Total conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_live_biomass    !! Total living biomass  
                                                                                       !! @tex $(gC m{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts)                     :: bm_alloc            !! Biomass increase, i.e. NPP per plant part 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_turnover        !! Total turnover rate  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_soil_carb!! Total soil and litter carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_carb     !! Total litter carbon 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_soil_carb       !! Total soil carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: carb_mass_variation !! Carbon Mass variation  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: cn_ind              !! Crown area of individuals 
                                                                                       !! @tex $(m^{2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: woodmass_ind        !! Woodmass of individuals (gC) 
    REAL(r_std), DIMENSION(npts,nvm,nparts)                     :: f_alloc             !! Fraction that goes into plant part 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: avail_tree          !! Space availability for trees 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: avail_grass         !! Space availability for grasses 
                                                                                       !! (0 to 1, unitless) 
    INTEGER                                                     :: j
    REAL(r_std),DIMENSION(npts)                                 :: prod10_total        !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod100_total       !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: cflux_prod_total    !! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm)                             :: veget_max_old       !! "Maximal" coverage fraction of a PFT  
                                                                                       !! (LAI-> infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,nvm)                            :: mortality           !! Fraction of individual dying this time 
                                                                                       !! step (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: vartmp              !! Temporary variable used to add history
    REAL(r_std), DIMENSION(npts,nvm)                            :: histvar             !! History variables
!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering stomate_lpj'

  
  !! 1. Initializations
    
    !! 1.1 Initialize variables to zero
    co2_to_bm(:,:) = zero
    co2_fire(:,:) = zero
    npp_daily(:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    harvest_above(:) = zero
    bm_to_litter(:,:,:) = zero
    cn_ind(:,:) = zero
    woodmass_ind(:,:) = zero
    turnover_daily(:,:,:) = zero
    
    !! 1.2  Initialize variables to veget_max
    veget_max_old(:,:) = veget_max(:,:)

    !! 1.3 Calculate some vegetation characteristics
    
    !! 1.3.1 Calculate some vegetation characteristics 
    !        Calculate cn_ind (individual crown mass) and individual height from
    !        state variables if running DGVM or dynamic mortality in static cover mode
    !??        Explain (maybe in the header once) why you mulitply with veget_max in the DGVM
    !??        and why you don't multiply with veget_max in stomate.
    IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       IF(control%ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove)+biomass(:,:,isapbelow) &
                  +biomass(:,:,iheartabove)+biomass(:,:,iheartbelow)) & 
                  *veget_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove)+biomass(:,:,isapbelow) &
                  +biomass(:,:,iheartabove)+biomass(:,:,iheartbelow))/ind(:,:)
          ENDWHERE
       ENDIF

       CALL crown (npts,  PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height)
    ENDIF

    !! 1.3.2 Prescribe characteristics if the vegetation is not dynamic
    !        IF the DGVM is not activated, the density of individuals and their crown
    !        areas don't matter, but they should be defined for the case we switch on
    !        the DGVM afterwards. At the first call, if the DGVM is not activated, 
    !        impose a minimum biomass for prescribed PFTs and declare them present.
    CALL prescribe (npts, &
         veget_max, PFTpresent, everywhere, when_growthinit, &
         biomass, leaf_frac, ind, cn_ind)


  !! 2. Climatic constraints for PFT presence and regenerativeness

    !   Call this even when DGVM is not activated so that "adapted" and "regenerate"
    !   are kept up to date for the moment when the DGVM is activated.
    CALL constraints (npts, dt_days, &
         t2m_month, t2m_min_daily,when_growthinit, &
         adapted, regenerate)

    
  !! 3. Determine introduction and elimination of PTS based on climate criteria
 
    IF ( control%ok_dgvm ) THEN
      
       !! 3.1 Calculate introduction and elimination
       CALL pftinout (npts, dt_days, adapted, regenerate, &
            neighbours, veget_max, &
            biomass, ind, cn_ind, age, leaf_frac, npp_longterm, lm_lastyearmax, senescence, &
            PFTpresent, everywhere, when_growthinit, need_adjacent, RIP_time, &
            co2_to_bm, &
            avail_tree, avail_grass)

       !! 3.2 Reset attributes for eliminated PFTs.
       !     This also kills PFTs that had 0 leafmass during the last year. The message
       !     "... after pftinout" is misleading in this case.
       CALL kill (npts, 'pftinout  ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter)

       
       !! 3.3 Calculate new crown area and diameter 
       !      Calculate new crown area, diameter and maximum vegetation cover**[No longer used in the subroutine]
       !      unsure whether this is really required
       !      - in theory this could ONLY be done at the END of stomate_lpj
       !      calculate woodmass of individual tree
       WHERE ((ind(:,:).GT.min_stomate))
          WHERE  ( veget_max(:,:) .GT. min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove)+biomass(:,:,isapbelow) &
                  +biomass(:,:,iheartabove)+biomass(:,:,iheartbelow))*veget_max(:,:))/ind(:,:)
          ELSEWHERE
             woodmass_ind(:,:) =(biomass(:,:,isapabove)+biomass(:,:,isapbelow) &
                  +biomass(:,:,iheartabove)+biomass(:,:,iheartbelow))/ind(:,:)
          ENDWHERE

       ENDWHERE
       
       ! Calculate crown area and diameter for all PFTs (including the newly established)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height)

    ENDIF
    
  !! 4. Phenology

    !! 4.1 Write values to history file
    !      Current values for ::when_growthinit and time_lowGPP
    CALL histwrite (hist_id_stomate, 'WHEN_GROWTHINIT', itime, when_growthinit, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'TIME_LOWGPP', itime, time_lowgpp, npts*nvm, horipft_index)

    ! Set and write values for ::PFTpresent
    WHERE(PFTpresent)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE
    CALL histwrite (hist_id_stomate, 'PFTPRESENT', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_midwinter
    WHERE(gdd_midwinter.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_midwinter
    ENDWHERE
    CALL histwrite (hist_id_stomate, 'GDD_MIDWINTER', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for ncd_dormance
    WHERE(ncd_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=ncd_dormance
    ENDWHERE
    CALL histwrite (hist_id_stomate, 'NCD_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    !! 4.2 Calculate phenology
    CALL phenology (npts, dt_days, PFTpresent, &
         veget_max, &
         tlong_ref, t2m_month, t2m_week, gpp_daily, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_month, moiavail_week, &
         gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
         senescence, time_lowgpp, time_hum_min, &
         biomass, leaf_frac, leaf_age, &
         when_growthinit, co2_to_bm)
    
  !! 5. Allocate C to different plant parts
    
    CALL alloc (npts, dt_days, &
         lai, veget_max, senescence, when_growthinit, &
         moiavail_week, tsoil_month, soilhum_month, &
         biomass, age, leaf_age, leaf_frac, rprof, f_alloc)

  !! 6. NPP, maintenance and growth respiration

    !! 6.1 Calculate NPP and respiration terms
    CALL npp_calc (npts, dt_days, &
         PFTpresent, &
         tlong_ref, t2m_daily, tsoil_daily, lai, rprof, &
         gpp_daily, f_alloc, bm_alloc, resp_maint_part,&
         biomass, leaf_age, leaf_frac, age, &
         resp_maint, resp_growth, npp_daily)

    !! 6.2 Kill slow growing PFTs in DGVM or STOMATE with constant mortality
    IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       CALL kill (npts, 'npp       ', lm_lastyearmax,  &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter)

       !! 6.2.1 Update wood biomass      
       !        For the DGVM
       IF(control%ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove)+biomass(:,:,isapbelow) &
                  +biomass(:,:,iheartabove)+biomass(:,:,iheartbelow)) & 
                  *veget_max(:,:))/ind(:,:)
          ENDWHERE

       ! For all pixels with individuals
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove)+biomass(:,:,isapbelow) &
                  +biomass(:,:,iheartabove)+biomass(:,:,iheartbelow))/ind(:,:)
          ENDWHERE
       ENDIF ! control%ok_dgvm

       !! 6.2.2 New crown area and maximum vegetation cover after growth
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind,&
            veget_max, cn_ind, height)

    ENDIF ! control%ok_dgvm
    
  !! 7. fire

    !! 7.1. Burn PFTs
    CALL fire (npts, dt_days, litterpart, &
         litterhum_daily, t2m_daily, lignin_struc, veget_max, &
         fireindex, firelitter, biomass, ind, &
         litter, dead_leaves, bm_to_litter, black_carbon, &
         co2_fire)

    !! 7.2 Kill PFTs in DGVM
    IF ( control%ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'fire      ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter)

    ENDIF ! control%ok_dgvm
 
  !! 8. Tree mortality

    ! Does not depend on age, therefore does not change crown area.
    CALL gap (npts, dt_days, &
         npp_longterm, turnover_longterm, lm_lastyearmax, &
         PFTpresent, biomass, ind, bm_to_litter, mortality)

    IF ( control%ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'gap       ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter)

    ENDIF

  !! 9. Calculate vcmax, vjmax and photosynthesis temperatures

    CALL vmax (npts, dt_days, &
         leaf_age, leaf_frac, &
         vcmax, vjmax)

    CALL assim_temp (npts, tlong_ref, t2m_month, &
         t_photo_min, t_photo_opt, t_photo_max)

  !! 10. Leaf senescence, new lai and other turnover processes

    CALL turn (npts, dt_days, PFTpresent, &
         herbivores, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_week,  moiavail_month,tlong_ref, t2m_month, t2m_week, veget_max, &
         leaf_age, leaf_frac, age, lai, biomass, &
         turnover_daily, senescence,turnover_time)

    !! 11. Light competition
    
    !! If not using constant mortality then kill with light competition
!    IF ( control%ok_dgvm .OR. .NOT.(lpj_gap_const_mort) ) THEN
    IF ( control%ok_dgvm ) THEN
 
       !! 11.1 Light competition
       CALL light (npts, dt_days, &
            veget_max, fpc_max, PFTpresent, cn_ind, lai, maxfpc_lastyear, &
            lm_lastyearmax, ind, biomass, veget_lastlight, bm_to_litter, mortality)
       
       !! 11.2 Reset attributes for eliminated PFTs
       CALL kill (npts, 'light     ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter)

    ENDIF

    
  !! 12. Establishment of saplings
    
    IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort ) THEN

       !! 12.1 Establish new plants
       CALL establish (npts, dt_days, PFTpresent, regenerate, &
            neighbours, resolution, need_adjacent, herbivores, &
            precip_lastyear, gdd0_lastyear, lm_lastyearmax, &
            cn_ind, lai, avail_tree, avail_grass, npp_longterm, &
            leaf_age, leaf_frac, &
            ind, biomass, age, everywhere, co2_to_bm, veget_max, woodmass_ind)

       !! 12.2 Calculate new crown area (and maximum vegetation cover)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height)

    ENDIF

  !! 13. Calculate final LAI and vegetation cover
    
    CALL cover (npts, cn_ind, ind, biomass, &
         veget_max, veget_max_old, lai, &
         litter, carbon, turnover_daily, bm_to_litter)

  !! 14. Update litter pools to account for harvest
 
    ! the whole litter stuff:
    !    litter update, lignin content, PFT parts, litter decay, 
    !    litter heterotrophic respiration, dead leaf soil cover.
    !    No vertical discretisation in the soil for litter decay.\n
    ! added by shilong for harvest
    IF(harvest_agri) THEN
       CALL harvest(npts, dt_days, veget_max, &
            bm_to_litter, turnover_daily, &
            harvest_above)
    ENDIF

  !! 15. Land cover change

    !shilong adde turnover_daily
    IF(EndOfYear) THEN
       IF (lcchange) THEN
          CALL lcchange_main (npts, dt_days, veget_max, veget_max_new, &
               biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
               co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, tree, cn_ind,flux10,flux100, &
!$               prod10,prod100,prod10_total,prod100_total,&
!$               convflux,cflux_prod_total,cflux_prod10,cflux_prod100,leaf_frac,&
               prod10,prod100,convflux,cflux_prod10,cflux_prod100,leaf_frac,&
               npp_longterm, lm_lastyearmax, litter, carbon)
       ENDIF
    ENDIF
    !MM déplacement pour initialisation correcte des grandeurs cumulées :
    cflux_prod_total(:) = convflux(:) + cflux_prod10(:) + cflux_prod100(:)
    prod10_total(:)=SUM(prod10,dim=2)
    prod100_total(:)=SUM(prod100,dim=2)
    
  !! 16. Total heterotrophic respiration

    tot_soil_carb=zero
    tot_litter_carb=zero
    DO j=2,nvm

    !!! Chloe test again :
!    IF (j .EQ. 14) THEN 
! Stomate LPJ est en journalier 
!        carbon(:,iactive,14)=carbon(:,iactive,14) - wpro_conso_carb(:)/48
!            carbon(:,iactive,14)=carbon(:,iactive,14) - wpro_conso_carb(:)
!           write(*,*) 'Are you daily? Chloe',carbon(:,iactive,14)
!    ENDIF
       tot_litter_carb(:,j) = tot_litter_carb(:,j) + (litter(:,istructural,j,iabove) + &
            &          litter(:,imetabolic,j,iabove) + &
            &          litter(:,istructural,j,ibelow) + litter(:,imetabolic,j,ibelow))

       tot_soil_carb(:,j) = tot_soil_carb(:,j) + (carbon(:,iactive,j) + &
            &          carbon(:,islow,j)+  carbon(:,ipassive,j))

    ENDDO
    tot_litter_soil_carb = tot_litter_carb + tot_soil_carb

    tot_live_biomass = biomass(:,:,ileaf) + biomass(:,:,isapabove) + biomass(:,:,isapbelow) +&
         &             biomass(:,:,iheartabove) + biomass(:,:,iheartbelow) + &
         &             biomass(:,:,iroot)+ biomass(:,:,ifruit)+ biomass(:,:,icarbres)

    tot_turnover = turnover_daily(:,:,ileaf) + turnover_daily(:,:,isapabove) + &
         &         turnover_daily(:,:,isapbelow) + turnover_daily(:,:,iheartabove) + &
         &         turnover_daily(:,:,iheartbelow) + turnover_daily(:,:,iroot) + &
         &         turnover_daily(:,:,ifruit) + turnover_daily(:,:,icarbres)

    tot_bm_to_litter = bm_to_litter(:,:,ileaf) + bm_to_litter(:,:,isapabove) +&
         &             bm_to_litter(:,:,isapbelow) + bm_to_litter(:,:,iheartbelow) +&
         &             bm_to_litter(:,:,iheartabove) + bm_to_litter(:,:,iroot) + &
         &             bm_to_litter(:,:,ifruit) + bm_to_litter(:,:,icarbres)

    carb_mass_variation(:)=-carb_mass_total(:)
    carb_mass_total(:)=SUM((tot_live_biomass+tot_litter_carb+tot_soil_carb)*veget_max,dim=2) + &
         &                 (prod10_total + prod100_total)
    carb_mass_variation(:)=carb_mass_total(:)+carb_mass_variation(:)
    
  !! 17. Write history

    CALL histwrite (hist_id_stomate, 'RESOLUTION_X', itime, &
         resolution(:,1), npts, hori_index)
    CALL histwrite (hist_id_stomate, 'RESOLUTION_Y', itime, &
         resolution(:,2), npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CONTFRAC', itime, &
         contfrac(:), npts, hori_index)

    CALL histwrite (hist_id_stomate, 'LITTER_STR_AB', itime, &
         litter(:,istructural,:,iabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LITTER_MET_AB', itime, &
         litter(:,imetabolic,:,iabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LITTER_STR_BE', itime, &
         litter(:,istructural,:,ibelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LITTER_MET_BE', itime, &
         litter(:,imetabolic,:,ibelow), npts*nvm, horipft_index)

    CALL histwrite (hist_id_stomate, 'DEADLEAF_COVER', itime, &
         deadleaf_cover, npts, hori_index)

    CALL histwrite (hist_id_stomate, 'TOTAL_SOIL_CARB', itime, &
         tot_litter_soil_carb, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'CARBON_ACTIVE', itime, &
         carbon(:,iactive,:), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'CARBON_SLOW', itime, &
         carbon(:,islow,:), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'CARBON_PASSIVE', itime, &
         carbon(:,ipassive,:), npts*nvm, horipft_index)
!Chloe provisoire
   CALL histwrite (hist_id_stomate, 'CH4_FLUX_TOT_0', itime, &
                    ch4_flux_density_tot_0, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_DIF_0', itime, &
                    ch4_flux_density_dif_0, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_BUB_0', itime, &
                    ch4_flux_density_bub_0, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_PLA_0', itime, &
                    ch4_flux_density_pla_0, npts, hori_index)

!
!Chloe peat :
   CALL histwrite (hist_id_stomate, 'CH4_FLUX_TOT_peat', itime, &
                    ch4_flux_density_tot_peat, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_DIF_peat', itime, &
                    ch4_flux_density_dif_peat, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_BUB_peat', itime, &
                    ch4_flux_density_bub_peat, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_PLA_peat', itime, &
                    ch4_flux_density_pla_peat, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CH4_FLUX_GOXID_peat', itime, &
                    ch4_flux_density_goxid_peat, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'WPRO_CH4_peat', itime, &
                        wpro_conso_carb, npts, hori_index)
 
    CALL histwrite (hist_id_stomate, 'TSURF_YEAR', itime, &
                    tsurf_year, npts, hori_index)

!Chloe--

    CALL histwrite (hist_id_stomate, 'T2M_MONTH', itime, &
         t2m_month, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'T2M_WEEK', itime, &
         t2m_week, npts, hori_index)

    CALL histwrite (hist_id_stomate, 'HET_RESP', itime, &
         resp_hetero(:,:), npts*nvm, horipft_index)

    CALL histwrite (hist_id_stomate, 'BLACK_CARBON', itime, &
         black_carbon, npts, hori_index)

    CALL histwrite (hist_id_stomate, 'FIREINDEX', itime, &
         fireindex(:,:), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LITTERHUM', itime, &
         litterhum_daily, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CO2_FIRE', itime, &
         co2_fire, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'CO2_TAKEN', itime, &
         co2_to_bm, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite (hist_id_stomate, 'CONVFLUX', itime, &
         convflux, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CFLUX_PROD10', itime, &
         cflux_prod10, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'CFLUX_PROD100', itime, &
         cflux_prod100, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'HARVEST_ABOVE', itime, &
         harvest_above, npts, hori_index)

    CALL histwrite (hist_id_stomate, 'LAI', itime, &
         lai, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'VEGET_MAX', itime, &
         veget_max, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'NPP', itime, &
         npp_daily, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'GPP', itime, &
         gpp_daily, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'IND', itime, &
         ind, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'CN_IND', itime, &
         cn_ind, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'WOODMASS_IND', itime, &
         woodmass_ind, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'TOTAL_M', itime, &
         tot_live_biomass, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LEAF_M', itime, &
         biomass(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SAP_M_AB', itime, &
         biomass(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SAP_M_BE', itime, &
         biomass(:,:,isapbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'HEART_M_AB', itime, &
         biomass(:,:,iheartabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'HEART_M_BE', itime, &
         biomass(:,:,iheartbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'ROOT_M', itime, &
         biomass(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'FRUIT_M', itime, &
         biomass(:,:,ifruit), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'RESERVE_M', itime, &
         biomass(:,:,icarbres), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'TOTAL_TURN', itime, &
         tot_turnover, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LEAF_TURN', itime, &
         turnover_daily(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SAP_AB_TURN', itime, &
         turnover_daily(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'ROOT_TURN', itime, &
         turnover_daily(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'FRUIT_TURN', itime, &
         turnover_daily(:,:,ifruit), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'TOTAL_BM_LITTER', itime, &
         tot_bm_to_litter, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'LEAF_BM_LITTER', itime, &
         bm_to_litter(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SAP_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'SAP_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,isapbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'HEART_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'HEART_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'ROOT_BM_LITTER', itime, &
         bm_to_litter(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'FRUIT_BM_LITTER', itime, &
         bm_to_litter(:,:,ifruit), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'RESERVE_BM_LITTER', itime, &
         bm_to_litter(:,:,icarbres), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'MAINT_RESP', itime, &
         resp_maint, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'GROWTH_RESP', itime, &
         resp_growth, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'AGE', itime, &
         age, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'HEIGHT', itime, &
         height, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'MOISTRESS', itime, &
         moiavail_week, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'VCMAX', itime, &
         vcmax, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'TURNOVER_TIME', itime, &
         turnover_time, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite (hist_id_stomate, 'PROD10', itime, &
         prod10, npts*11, horip11_index)
    CALL histwrite (hist_id_stomate, 'PROD100', itime, &
         prod100, npts*101, horip101_index)
    CALL histwrite (hist_id_stomate, 'FLUX10', itime, &
         flux10, npts*10, horip10_index)
    CALL histwrite (hist_id_stomate, 'FLUX100', itime, &
         flux100, npts*100, horip100_index)

    IF ( hist_id_stomate_IPCC > 0 ) THEN
       vartmp(:)=SUM(tot_live_biomass*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cVeg", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_litter_carb*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_soil_carb*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(prod10_total + prod100_total)/1e3
       CALL histwrite (hist_id_stomate_IPCC, "cProduct", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=carb_mass_variation/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cMassVariation", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(lai*veget_max,dim=2)*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "lai", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(gpp_daily*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "gpp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((resp_maint+resp_growth)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "ra", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(npp_daily*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "npp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_hetero*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "rh", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_fire*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "fFire", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=harvest_above/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "fHarvest", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_total/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "fLuc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
            &        *veget_max,dim=2)-cflux_prod_total-harvest_above)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "nbp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((tot_bm_to_litter + tot_turnover)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "fVegLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(soilcarbon_input,dim=2)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "fLitterSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(biomass(:,:,ileaf)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((biomass(:,:,isapabove)+biomass(:,:,iheartabove))*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,iroot) + biomass(:,:,isapbelow) + biomass(:,:,iheartbelow) ) &
            &        *veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cRoot", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,icarbres) + biomass(:,:,ifruit))*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cMisc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,iabove)+litter(:,imetabolic,:,iabove))*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cLitterAbove", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,ibelow)+litter(:,imetabolic,:,ibelow))*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cLitterBelow", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,iactive,:)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cSoilFast", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,islow,:)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cSoilMedium", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,ipassive,:)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "cSoilSlow", itime, &
            vartmp, npts, hori_index)
       DO j=1,nvm
          histvar(:,j)=veget_max(:,j)*contfrac(:)*100
       ENDDO
       CALL histwrite (hist_id_stomate_IPCC, "landCoverFrac", itime, &
            histvar, npts*nvm, horipft_index)
       !-
       vartmp(:)=zero
       DO j=2,nvm
          IF(is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite (hist_id_stomate_IPCC, "treeFracPrimDec", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j=2,nvm
          IF(is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite (hist_id_stomate_IPCC, "treeFracPrimEver", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j=2,nvm
          IF(is_c3(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite (hist_id_stomate_IPCC, "c3PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j=2,nvm
          IF(is_c4(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite (hist_id_stomate_IPCC, "c4PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=SUM(resp_growth*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "rGrowth", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_maint*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "rMaint", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,ileaf)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "nppLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,isapabove)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "nppWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( bm_alloc(:,:,isapbelow) + bm_alloc(:,:,iroot) )*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite (hist_id_stomate_IPCC, "nppRoot", itime, &
            vartmp, npts, hori_index)

       CALL histwrite (hist_id_stomate_IPCC, 'RESOLUTION_X', itime, &
            resolution(:,1), npts, hori_index)
       CALL histwrite (hist_id_stomate_IPCC, 'RESOLUTION_Y', itime, &
            resolution(:,2), npts, hori_index)
       CALL histwrite (hist_id_stomate_IPCC, 'CONTFRAC', itime, &
            contfrac(:), npts, hori_index)

    ENDIF

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving stomate_lpj'

  END SUBROUTINE StomateLpj


!! ================================================================================================================================
!! SUBROUTINE   : harvest
!!
!>\BRIEF        Harvest of croplands
!!
!! DESCRIPTION  : To take into account biomass harvest from crop (mainly to take 
!! into account for the reduced litter input and then decreased soil carbon. it is a 
!! constant (40\%) fraction of above ground biomass.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::harvest_above the harvested biomass
!!
!! REFERENCE(S) :
!! - Piao, S., P. Ciais, P. Friedlingstein, N. de Noblet-Ducoudre, P. Cadule, N. Viovy, and T. Wang. 2009. 
!!   Spatiotemporal patterns of terrestrial carbon cycle during the 20th century. Global Biogeochemical 
!!   Cycles 23:doi:10.1029/2008GB003339.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE harvest(npts, dt_days, veget_max, &
       bm_to_litter, turnover_daily, &
       harvest_above)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts            !! Domain size (unitless) 
    REAL(r_std), INTENT(in)                                :: dt_days         !! Time step (days)                               
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_max       !! new "maximal" coverage fraction of a PFT (LAI -> 
                                                                              !! infinity) on ground @tex $(m^2 m^{-2})$ @endtex 
    
   !! 0.2 Output variables
   
   !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout) :: bm_to_litter     !! [DISPENSABLE] conversion of biomass to litter 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout) :: turnover_daily   !! Turnover rates 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: harvest_above    !! harvest above ground biomass for agriculture 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    !! 0.4 Local variables

    INTEGER(i_std)                                         :: i, j, k, l, m    !! indices                       
    REAL(r_std)                                            :: above_old        !! biomass of previous time step 
                                                                               !! @tex $(gC m^{-2})$ @endtex 
!_ ================================================================================================================================

  !! 1. Yearly initialisation

    above_old             = zero
    harvest_above         = zero

    DO i = 1, npts
       DO j = 1,nvm
          IF (.NOT. natural(j)) THEN
             above_old = turnover_daily(i,j,ileaf) + turnover_daily(i,j,isapabove) + &
                  &       turnover_daily(i,j,iheartabove) + turnover_daily(i,j,ifruit) + &
                  &       turnover_daily(i,j,icarbres) + turnover_daily(i,j,isapbelow) + &
                  &       turnover_daily(i,j,iheartbelow) + turnover_daily(i,j,iroot)

             turnover_daily(i,j,ileaf) = turnover_daily(i,j,ileaf)*frac_turnover_daily
             turnover_daily(i,j,isapabove) = turnover_daily(i,j,isapabove)*frac_turnover_daily
             turnover_daily(i,j,isapbelow) = turnover_daily(i,j,isapbelow)*frac_turnover_daily
             turnover_daily(i,j,iheartabove) = turnover_daily(i,j,iheartabove)*frac_turnover_daily
             turnover_daily(i,j,iheartbelow) = turnover_daily(i,j,iheartbelow)*frac_turnover_daily
             turnover_daily(i,j,iroot) = turnover_daily(i,j,iroot)*frac_turnover_daily
             turnover_daily(i,j,ifruit) = turnover_daily(i,j,ifruit)*frac_turnover_daily
             turnover_daily(i,j,icarbres) = turnover_daily(i,j,icarbres)*frac_turnover_daily
             harvest_above(i)  = harvest_above(i) + veget_max(i,j) * above_old *(un - frac_turnover_daily)
          ENDIF
       ENDDO
    ENDDO

!$    harvest_above = harvest_above
  END SUBROUTINE harvest
END MODULE stomate_lpj
