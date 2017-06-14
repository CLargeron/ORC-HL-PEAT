!x =================================================================================================================================
! MODULE       : stomate
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines that: (1) initialize all variables in 
!! stomate, (2) read and write forcing files of stomate and the soil component,
!! (3) aggregates and convert variables to handle the different time steps 
!! between sechiba and stomate, (4) call subroutines that govern major stomate
!! processes (litter, soil, and vegetation dynamics) and (5) structures these tasks 
!! in stomate_main
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate.f90 $
!! $Date: 2012-08-02 15:40:27 +0200 (Thu, 02 Aug 2012) $
!! $Revision: 965 $
!! \n
!_ ================================================================================================================================

MODULE stomate

  ! Modules used:
  USE netcdf
  USE ioipsl
  USE defprec
  USE grid
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE stomate_io
  USE stomate_data
  USE stomate_season
  USE stomate_lpj
  USE stomate_assimtemp
  USE stomate_litter
  USE stomate_vmax
  USE stomate_soilcarbon
  USE stomate_resp
  USE parallel
  ! USE Write_field_p
  !Chloe peat :
  USE stomate_wet_ch4_pt_ter_peat
  !Test Chloe : 
   USE stomate_wet_ch4_pt_ter_0
  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC stomate_main,stomate_clear,init_forcing,forcing_read


  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: biomass              !! Biomass per ground area @tex $(gC m^{-2})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: veget_cov_max        !! Maximal fractional coverage: maximum share of a pixel
                                                                         !! taken by a PFT 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: veget_cov_max_new    !! Updated maximal fractional coverage (unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ind                  !! Vegetation density, number of individuals per unit 
                                                                         !! ground area @tex $(m^{-2})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: age                  !! Age of PFT it normalized by biomass - can increase and
                                                                         !! decrease - (years)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: adapted              !! Winter too cold for PFT to survive (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: regenerate           !! Winter sufficiently cold to produce viable seeds 
                                                                         !! (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: everywhere           !! Is the PFT everywhere in the grid box or very localized 
                                                                         !! (after its intoduction)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: fireindex            !! Probability of fire (unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: veget_lastlight      !! Vegetation fractions (on ground) after last light 
                                                                         !! competition (unitless) 
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)   :: fpc_max              !! "maximal" coverage fraction of a grid box (LAI -> 
                                                                         !! infinity) on ground. [??CHECK??] It's set to zero here, 
                                                                         !! and then is used once in lpj_light.f90 to test if 
                                                                         !! fpc_nat is greater than it. Something seems missing
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: PFTpresent           !! PFT exists (equivalent to veget > 0 for natural PFTs)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: senescence           !! The PFT is senescent
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: need_adjacent        !! This PFT needs to be in present in an adjacent gridbox 
                                                                         !! if it is to be introduced in a new gridbox
!--
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: humrel_daily         !! Daily plant available water -root profile weighted 
                                                                         !! (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: humrel_week          !! "Weekly" plant available water -root profile weighted
                                                                         !! (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: humrel_month         !! "Monthly" plant available water -root profile weighted
                                                                         !! (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxhumrel_lastyear   !! Last year's max plant available water -root profile 
                                                                         !! weighted (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxhumrel_thisyear   !! This year's max plant available water -root profile 
                                                                         !! weighted (0-1, unitless) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: minhumrel_lastyear   !! Last year's min plant available water -root profile 
                                                                         !! weighted (0-1, unitless)  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: minhumrel_thisyear   !! This year's minimum plant available water -root profile
                                                                         !! weighted (0-1, unitless)
!---  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_daily            !! Daily air temperature at 2 meter (K)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_week             !! Mean "weekly" (default 7 days) air temperature at 2 
                                                                         !! meter (K)  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_month            !! Mean "monthly" (default 20 days) air temperature at 2 
                                                                         !! meter (K)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_longterm         !! Mean "Long term" (default 3 years) air temperature at 
!---
!Chloe peat
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: tsurf_year           !! "annual" surface temperatures (K)                                                    
                                                                         !! 2 meter (K) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_min_daily        !! Daily minimum air temperature at 2 meter (K)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: tlong_ref            !! "Long term" (default 3 years) reference temperature at
                                                                         !! 2 meter (K) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: tsurf_daily          !! Daily surface temperatures (K)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_daily         !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_lastyear      !! Last year's annual precipitation sum 
                                                                         !! @tex $??(mm year^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_thisyear      !! This year's annual precipitation sum 
                                                                         !! @tex $??(mm year^{-1})$ @endtex 
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: soilhum_daily        !! Daily soil humidity (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: soilhum_month        !! Soil humidity - integrated over a month (0-1, unitless) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: tsoil_daily          !! Daily soil temperatures (K)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: tsoil_month          !! Soil temperatures at each soil layer integrated over a
                                                                         !! month (K) 
!--- 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: litterhum_daily      !! Daily litter humidity (0-1, unitless)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: control_moist        !! Moisture control of heterotrophic respiration 
                                                                         !! (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: control_temp         !! Temperature control of heterotrophic respiration at the
                                                                         !! different soil levels (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: control_moist_daily  !! Moisture control of heterotrophic respiration daily 
                                                                         !! (0-1, unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: control_temp_daily   !! Temperature control of heterotrophic respiration, above
                                                                         !! and below daily (0-1, unitless)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: gdd0_lastyear        !! Last year's annual Growing Degree Days,
                                                                         !! threshold 0 deg C (K) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: gdd0_thisyear        !! This year's annual Growing Degree Days,
                                                                         !! threshold 0 deg C (K)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_m5_dormance      !! Growing degree days for onset of growing season, 
                                                                         !! threshold -5 deg C (K)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_midwinter        !! Growing degree days for onset of growing season, 
                                                                         !! since midwinter (K)

!Chloe++ : 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: wtsoil_daily         !! Daily mean wt_soil used for Walter model.     
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: wtold                !! Water table in previous day (used in bruno's code)  
!Chloe--                        

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ncd_dormance         !! Number of chilling days since leaves were lost (days) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ngd_minus5           !! Number of growing days, threshold -5 deg C (days)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: hum_min_dormance     !! Minimum moisture during dormance (0-1, unitless) 
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_daily            !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_week             !! Mean "weekly" (default 7 days) GPP  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_week             !! Mean "weekly" (default 7 days) GPP  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxgppweek_lastyear  !! Last year's maximum "weekly" GPP  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxgppweek_thisyear  !! This year's maximum "weekly" GPP  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex  
!Chloe :
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxnppweek_thisyear  !! This year's maximum "weekly" NPP  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex  
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_daily            !! Daily net primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_longterm         !! "Long term" (default 3 years) net primary productivity 
                                                                         !! per ground area  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex   
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_equil            !! Equilibrium NPP written to forcesoil 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: npp_tot              !! Total NPP written to forcesoil 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: resp_maint_part_radia!! Maintenance respiration of different plant parts per 
                                                                         !! total ground area at Sechiba time step  
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: resp_maint_part      !! Maintenance respiration of different plant parts per
                                                                         !! total ground area at Stomate time step 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_radia     !! Maintenance respiration per ground area at Sechiba time
                                                                         !! step   
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_d         !! Maintenance respiration per ground area at Stomate time 
                                                                         !! step  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_growth_d        !! Growth respiration per ground area 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_d        !! Heterotrophic respiration per ground area 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_radia    !! Heterothrophic respiration per ground area at Sechiba
                                                                         !! time step 
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex 
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: turnover_time        !! Turnover time of grasses 
                                                                         !! @tex $(dt_slow^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: turnover_daily       !! Senescence-driven turnover (better: mortality) of 
                                                                         !! leaves and roots  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: turnover_littercalc  !! Senescence-driven turnover (better: mortality) of 
                                                                         !! leaves and roots at Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: turnover_longterm    !! "Long term" (default 3 years) senescence-driven 
                                                                         !! turnover (better: mortality) of leaves and roots 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: bm_to_litter         !! Background (not senescence-driven) mortality of biomass
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: bm_to_littercalc     !! conversion of biomass to litter per ground area at 
                                                                         !! Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: dead_leaves          !! Metabolic and structural pools of dead leaves on ground
                                                                         !! per PFT @tex $(gC m^{-2})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: litter               !! Above and below ground metabolic and structural litter 
                                                                         !! per ground area 
                                                                         !! @tex $(gC m^{-2})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: litterpart           !! Fraction of litter above the ground belonging to 
                                                                         !! different litter pools (unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: firelitter           !! Total litter above the ground that could potentially 
                                                                         !! burn @tex $(gC m^{-2})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: soilcarbon_input     !! Quantity of carbon going into carbon pools from litter
                                                                         !! decomposition per ground area  at Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: soilcarbon_input_daily !! Daily quantity of carbon going into carbon pools from
                                                                           !! litter decomposition per ground area 
                                                                           !! @tex $(gC m^{-2} day^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: carbon               !! Soil carbon pools per ground area: active, slow, or 
                                                                         !! passive, @tex $(gC m^{-2})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: lignin_struc         !! Ratio Lignine/Carbon in structural litter for above and
                                                                         !! below ground compartments (unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: black_carbon         !! Black carbon on the ground 
                                                                         !! @tex $(gC m^{-2})$ @endtex   
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lm_lastyearmax       !! Last year's maximum leaf mass per ground area for each
                                                                         !! PFT @tex $(gC m^{-2})$ @endtex  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lm_thisyearmax       !! This year's maximum leaf mass per ground area for each
                                                                         !! PFT @tex $(gC m^{-2})$ @endtex  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxfpc_lastyear      !! Last year's maximum fpc for each natural PFT, on ground
                                                                         !! [??CHECK] fpc but this ones look ok (computed in 
                                                                         !! season, used in light)?? 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxfpc_thisyear      !! This year's maximum fpc for each PFT, on ground (see 
                                                                         !! stomate_season), [??CHECK] fpc but this ones look ok 
                                                                         !! (computed in season, used in light)??
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: leaf_age             !! Age of different leaf classes (days)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: leaf_frac            !! PFT fraction of leaf mass in leaf age class (0-1, 
                                                                         !! unitless) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: when_growthinit      !! Days since beginning of growing season (days)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: herbivores           !! Time constant of probability of a leaf to be eaten by a
                                                                         !! herbivore (days)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: RIP_time             !! How much time ago was the PFT eliminated for the last 
                                                                         !! time (year)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: time_lowgpp          !! Duration of dormancy (days)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: time_hum_min         !! Time elapsed since strongest moisture limitation (days) 
!---
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: clay_fm              !! Soil clay content (0-1, unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: clay_fm_g            !! Soil clay content (0-1, unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: precip_fm            !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex,
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: precip_fm_g          !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex,
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: litterhum_daily_fm   !! Daily relative humidity of litter (0-1, unitless), 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: litterhum_daily_fm_g !! Daily relative humidity of litter (0-1, unitless), 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: t2m_daily_fm         !! Daily air temperature at 2 meter (K), parallel 
                                                                         !! computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: t2m_daily_fm_g       !! Daily air temperature at 2 meter (K), parallel 
                                                                         !! computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: t2m_min_daily_fm     !! Daily minimum air temperature at 2 meter (K), 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: t2m_min_daily_fm_g   !! Daily minimum air temperature at 2 meter (K), 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: tsurf_daily_fm       !! Daily surface temperatures (K), parallel 
                                                                         !! computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)         :: tsurf_daily_fm_g     !! Daily surface temperatures (K), parallel 
                                                                         !! computing 
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: tsoil_daily_fm       !! Daily soil temperatures (K), parallel computing 
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: tsoil_daily_fm_g     !! Daily soil temperatures (K), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: soilhum_daily_fm     !! Daily soil humidity (0-1, unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: soilhum_daily_fm_g   !! Daily soil humidity (0-1, unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: humrel_daily_fm      !! Daily relative humidity of atmosphere (0-1, unitless), 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: humrel_daily_fm_g    !! Daily relative humidity of atmosphere (0-1, unitless), 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: gpp_daily_fm         !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex, 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: gpp_daily_fm_g       !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} day^{-1})$ @endtex, 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: veget_fm             !! Vegetation coverage taking into account non-biological
                                                                         !! coverage (unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: veget_fm_g           !! Vegetation coverage taking into account non-biological
                                                                         !! coverage (unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: veget_max_fm         !! Maximum vegetation coverage taking into account 
                                                                         !! non-biological coverage (unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: veget_max_fm_g       !! Maximum vegetation coverage taking into account none 
                                                                         !! biological coverage (unitless), parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: lai_fm               !! Leaf area index @tex $@tex $(m^2 m^{-2})$ @endtex$ @endtex, 
                                                                         !! parallel computing
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)       :: lai_fm_g             !! Leaf area index @tex $@tex $(m^2 m^{-2})$ @endtex$ @endtex, 
                                                                         !! parallel computing
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_fire             !! Carbon emitted to the atmosphere by burning living 
                                                                         !! and dead biomass 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_to_bm_dgvm       !! Psuedo-photosynthesis,C used to provide seedlings with
                                                                         !! an initial biomass, arbitrarily removed from the 
                                                                         !! atmosphere  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_flux_daily       !! Daily net CO2 flux between atmosphere and biosphere 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
                                                                         !! [??CHECK] sign convention?
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_flux_monthly     !! Monthly net CO2 flux between atmosphere and biosphere 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex
                                                                         !! [??CHECK] sign convention? 

!Chloe test provisoire
!pour wetland avc Water Table Depth (WTD) = 0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_0
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_0 !concentration dim = (npts,n)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_0 !concentration au pas de temps precedent
!Chloe fin du test provisoire

!---
!Chloe
!pour wetland avc WTD = -x1
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_tot_peat
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_dif_peat
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_bub_peat
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_pla_peat
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: ch4_flux_density_goxid_peat
!  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: wpro_conso_carb 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uo_peat
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: uold2_peat

!FLAG for CH4 from wetland
  LOGICAL,SAVE        :: CH4_calcul, CH4_WTD1 !atmoshpere methane concentration/ near surface methane concentration
  REAL(r_std), SAVE   :: CH4atmo_CONC

!Chloe--

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: prod10               !! Wood products remaining in the 10 year-turnover pool 
                                                                         !! after the annual release for each compartment 
                                                                         !! @tex $(gC m^{-2})$ @endtex    
                                                                         !! (0:10 input from year of land cover change),
                                                                         !! dimension(#pixels,0:10 years)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: prod100              !! Wood products remaining in the 100 year-turnover pool
                                                                         !! after the annual release for each compartment
                                                                         !! @tex $(gC m^{-2})$ @endtex  
                                                                         !! (0:100 input from year of land cover change), 
                                                                         !! dimension(#pixels,0:100 years)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: flux10               !! Wood decomposition from the 10 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
                                                                         !! dimension(#pixels,0:10)  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: flux100              !! Wood decomposition from the 100 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
                                                                         !! dimension(#pixels,0:100)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: convflux             !! Release during first year following land cover change 
                                                                         !! (paper, burned, etc...) 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod10         !! Total annual release from the 10 year-turnover pool
                                                                         !! sum of flux10  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod100        !! Total annual release from the 100 year-turnover pool 
                                                                         !! sum of flux100 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: harvest_above        !! Harvest of above ground biomass for agriculture -not 
                                                                         !! just from land use change 
                                                                         !! @tex $(??gC m^{-2} dt_slow^{-1})$ @endtex
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: carb_mass_total      !! Total on-site and off-site C pool 
                                                                         !! @tex $(??gC m^{-2})$ @endtex                        
!---
  REAL(r_std),SAVE                               :: dt_days=zero         !! Time step of STOMATE (days) 
  REAL(r_std),SAVE                               :: day_counter=zero     !! Time step each day sechiba (dtradia)   
  INTEGER(i_std),SAVE                            :: date=0               !! Date (days) 
  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:)   :: nforce               !! Number of states calculated for the soil forcing 
                                                                         !! variables (unitless), dimension(::nparan*::nbyear) both 
                                                                         !! given in the run definition file    
  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:)   :: isf                  !! Index for number of time steps that can be stored in 
                                                                         !! memory (unitless), dimension (#nsfm)
  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:)   :: nf_cumul             !! Number of years over which the average is calculated in
                                                                         !! forcesoil when cumul flag is set, dimension (#nsft)
                                                                         !! [??CHECK] definition the dimension is number of 
                                                                         !! timesteps in a year?
  INTEGER,PARAMETER                              :: r_typ = nf90_real4   !! Specify data format (server dependent)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:)          :: nf_written           !! Flag indicating whether the forcing data have been 
                                                                         !! written
!---
  LOGICAL, SAVE                                  :: do_slow=.FALSE.      !! Flag that determines whether stomate_accu calculates
                                                                         !! the sum(do_slow=.FALSE.) or the mean 
                                                                         !! (do_slow=.TRUE.)
  LOGICAL, SAVE                                  :: EndOfYear=.FALSE.    !! Update annual variables? This variable must be 
                                                                         !! .TRUE. once a year
  LOGICAL, SAVE                                  :: EndOfMonth=.FALSE.   !! Update monthly variables? This variable must be 
                                                                         !!.TRUE. once a month 
  LOGICAL, SAVE                                  :: l_first_stomate = .TRUE.!! Is this the first call of stomate?
  LOGICAL, SAVE                                  :: cumul_forcing=.FALSE.!! flag for cumul of forcing if teststomate
  LOGICAL, SAVE                                  :: cumul_Cforcing=.FALSE.  !! Flag, if internal parameter cumul_Cforcing is 
                                                                            !! TRUE then ::nbyear (defined in run definition 
                                                                            !! file will be forced to 1 later in this module. If 
                                                                            !! FALSE the mean over ::nbyear is written in forcesoil
!---   
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: harvest_above_monthly   !! [??CHECK] post-processing - should be removed?
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: cflux_prod_monthly      !! [??CHECK] post-processing - should be removed?
!---
 
PUBLIC clay_fm, humrel_daily_fm, litterhum_daily_fm, t2m_daily_fm, &
   & t2m_min_daily_fm, tsurf_daily_fm, tsoil_daily_fm, soilhum_daily_fm, &
   & precip_fm, gpp_daily_fm, veget_fm, veget_max_fm, lai_fm
PUBLIC  dt_days, day_counter, date, do_slow, EndOfYear
PUBLIC isf, nf_written

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_main
!!
!>\BRIEF        Manages variable initialisation, reading and writing forcing 
!! files, aggregating data at stomate's time step (dt_slow), aggregating data
!! at longer time scale (i.e. for phenology) and uses these forcing to calculate
!! CO2 fluxes (NPP and respirations) and C-pools (litter, soil, biomass, ...)
!!
!! DESCRIPTION  : The subroutine manages 
!! divers tasks:
!! (1) Initializing all variables of stomate (first call)
!! (2) Reading and writing forcing data (last call)
!! (3) Adding CO2 fluxes to the IPCC history files
!! (4) Converting the time steps of variables to maintain consistency between
!! sechiba and stomate
!! (5) Use these variables to call stomate_lpj, maint_respiration, littercalc,
!! soilcarbon. The called subroutines handle: climate constraints 
!! for PFTs, PFT dynamics, Phenology, Allocation, NPP (based on GPP and
!! authothropic respiration), fire, mortality, vmax, assimilation temperatures,
!! all turnover processes, light competition, sapling establishment, lai,  
!! land cover change and litter and soil dynamics.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): deadleaf_cover, assim_param, lai, height, veget, 
!! veget_max, veget_max_new, totfrac_nobio_new, resp_maint, 
!! resp_hetero,resp_growth, co2_flux, fco2_lu.
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale=0.5]{stomatemainflow.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
  
SUBROUTINE stomate_main &
       & (kjit, kjpij, kjpindex, dtradia, dt_slow, &
       &  ldrestart_read, ldrestart_write, ldforcing_write, ldcarbon_write, &
       &  index, lalo, neighbours, resolution, contfrac, totfrac_nobio, clay, &
       &  t2m, t2m_min, temp_sol, stempdiag, &
       &  humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, &
       &  gpp, deadleaf_cover, assim_param, &
       &  lai, frac_age, height, veget, veget_max, &
       &  veget_max_new, totfrac_nobio_new, &
       &  hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       &  co2_flux, fco2_lu,resp_maint,resp_hetero,resp_growth,nstm, &
       !Chloe
       & wt_soil,wt_soil2)
         
    
    IMPLICIT NONE

    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

!Chloe copie de shushi :
!pss:+ nb de niveau vertical pour le calcul de la concentration de CH4
    !INTEGER(i_std),PARAMETER                         :: n = 371
!Chloe change le 371 en 171 :
!    INTEGER(i_std),PARAMETER                     :: n = 171
!Chloe--
!INTEGER(i_std),INTENT(in)                       :: nvert   

!Chloe add nstm et wt_soil : 
    INTEGER(i_std),INTENT(in)                       :: nstm              !! Time step number (unitless)
!Chloe test debug a changer : inout en in 
    REAL(r_std),DIMENSION(kjpindex,nstm),INTENT(in) :: wt_soil          !! Position nappe phréatique utile pour em methane     
   REAL(r_std),DIMENSION(kjpindex,nstm),INTENT(in) :: wt_soil2          !! Position nappe phréatique utile pour em methane      
!    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: wpro_conso_carb  !! Amount of Carbon consumed by methanogenese
!Chloe--    
    INTEGER(i_std),INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpij             !! Total size of the un-compressed grid (unitless)
    REAL(r_std),INTENT(in)                          :: dtradia           !! Time step of SECHIBA (seconds)
    REAL(r_std),INTENT(in)                          :: dt_slow           !! Time step of STOMATE (days)
    LOGICAL,INTENT(in)                              :: ldrestart_read    !! Logical for _restart_ file to read
    LOGICAL,INTENT(in)                              :: ldrestart_write   !! Flag to write restart file for stomate
    LOGICAL,INTENT(in)                              :: ldforcing_write   !! Flag to write forcing file for main processes in 
                                                                         !! stomate 
    LOGICAL,INTENT(in)                              :: ldcarbon_write    !! Flag to write forcing file for soil processes in 
                                                                         !! stomate 
    INTEGER(i_std),INTENT(in)                       :: rest_id_stom      !! STOMATE's _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT(in)                       :: hist_id_stom      !! STOMATE's _history_ file identifier (unitless)
    INTEGER(i_std),INTENT(in)                       :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier 
                                                                         !! (unitless) 
    INTEGER(i_std),DIMENSION(kjpindex),INTENT(in)   :: index             !! Indices of the pixels on the map. Stomate uses a 
                                                                         !! reduced grid excluding oceans. ::index contains 
                                                                         !! the indices of the terrestrial pixels only 
                                                                         !! (unitless) 
    INTEGER(i_std),DIMENSION(kjpindex,8),INTENT(in) :: neighbours        !! Neighoring grid points if land for the DGVM 
                                                                         !! (unitless) 
    REAL(r_std),DIMENSION(kjpindex,2),INTENT(in)    :: lalo              !! Geographical coordinates (latitude,longitude) 
                                                                         !! for pixels (degrees) 
    REAL(r_std),DIMENSION(kjpindex,2),INTENT(in)    :: resolution        !! Size in x an y of the grid (m) - surface area of 
                                                                         !! the gridbox 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)   :: contfrac          !! Fraction of continent in the grid cell (unitless)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: totfrac_nobio     !! Fraction of grid cell covered by lakes, land 
                                                                         !! ice, cities, ... (unitless) 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)  :: humrel            !! Relative humidity ("moisture availability") 
                                                                         !! (0-1, unitless) 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: t2m               !! 2 m air temperature (K)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: t2m_min           !! Minimum 2 m air temp. during forcing time step 
                                                                         !! (K) 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: temp_sol          !! Surface temperature (K)
    REAL(r_std),DIMENSION(kjpindex,nbdl),INTENT(in) :: stempdiag         !! Soil humidity (0-1, unitless)
    REAL(r_std),DIMENSION(kjpindex,nbdl),INTENT(in) :: shumdiag          !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: litterhumdiag     !! Litter humidity (0-1, unitless)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: precip_rain       !! Rain precipitation  
                                                                         !! @tex $(mm dt_slow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: precip_snow       !! Snow precipitation  
                                                                         !! @tex $(mm dt_slow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)  :: gpp               !! GPP of total ground area  
                                                                         !! @tex $(gC m^{-2} time step^{-1})$ @endtex 
                                                                         !! Calculated in sechiba, account for vegetation 
                                                                         !! cover and effective time step to obtain ::gpp_d 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)  :: veget_max_new     !! New "maximal" coverage fraction of a PFT (LAI -> 
                                                                         !! infinity) on ground only if EndOfYear is 
                                                                         !! activated (unitless) 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)      :: totfrac_nobio_new !! New fraction of grid cell covered by lakes, land 
                                                                         !! ice, cities, ... (unitless) 
    INTEGER(i_std),INTENT(in)                       :: hist_id           !! ?? [DISPENSABLE] SECHIBA's _history_ file 
                                                                         !! identifier 
    INTEGER(i_std),INTENT(in)                       :: hist2_id          !! ?? [DISPENSABLE] SECHIBA's _history_ file 2 
                                                                         !! identifier 

    !! 0.2 Output variables

    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(out) :: co2_flux          !! CO2 flux between atmosphere and biosphere per 
                                                                         !! average ground area 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex  
                                                                         !! [??CHECK] sign convention? 
    REAL(r_std),DIMENSION(kjpindex),INTENT(out)     :: fco2_lu           !! CO2 flux between atmosphere and biosphere from 
                                                                         !! land-use (without forest management)  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex  
                                                                         !! [??CHECK] sign convention? 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(out) :: resp_maint        !! Maitenance component of autotrophic respiration in 
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(out) :: resp_growth       !! Growth component of autotrophic respiration in 
                                                                         !! @tex ($gC m^{-2} dt_slow^{-1}$) @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(out) :: resp_hetero       !! Heterotrophic respiration in  
                                                                         !! @tex $(gC m^{-2} dt_slow^{-1})$ @endtex  

    !! 0.3 Modified
   
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(inout)       :: lai            !! Leaf area inex @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)          :: veget          !! Fraction of vegetation type including 
                                                                              !! non-biological fraction (unitless) 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(inout)       :: veget_max      !! Maximum fraction of vegetation type including 
                                                                              !! non-biological fraction (unitless) 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(inout)       :: height         !! Height of vegetation (m)
    REAL(r_std),DIMENSION(kjpindex,nvm,npco2),INTENT(inout) :: assim_param    !! min+max+opt temperatures (K) & vmax for 
                                                                              !! photosynthesis  
                                                                              !! @tex $(\mu mol m^{-2}s^{-1})$ @endtex  
    REAL(r_std),DIMENSION(kjpindex),INTENT(inout)           :: deadleaf_cover !! Fraction of soil covered by dead leaves 
                                                                              !! (unitless) 
    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages),INTENT(inout):: frac_age    !! Age efficacity from STOMATE 

    !! 0.4 local variables
    
    REAL(r_std)                                   :: day_counter_read         !! Day read in restart file (days)
    REAL(r_std)                                   :: dt_days_read             !! STOMATE time step read in restart file (days)
    INTEGER(i_std)                                :: l,k,ji, jv, i, j         !! indices    
    INTEGER(i_std)                                :: date_read                !! STOMATE date read in restart file (days)
    REAL(r_std),PARAMETER                         :: max_dt_days = 5.         !! Maximum STOMATE time step (days)
    REAL(r_std)                                   :: hist_days                !! Writing frequency for history file (days)
    REAL(r_std),DIMENSION(0:nbdl)                 :: z_soil                   !! Variable to store depth of the different soil 
                                                                              !! layers (m) 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: rprof                    !! Coefficient of the exponential functions that 
                                                                              !! relates root density to soil depth (unitless) 
    REAL(r_std),DIMENSION(kjpindex)               :: cvegtot                  !! Total "vegetation" cover (unitless)
    REAL(r_std),DIMENSION(kjpindex)               :: precip                   !! Total liquid and solid precipitation  
                                                                              !! @tex $(??mm dt_slow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_d                    !! Gross primary productivity per ground area 
                                                                              !! @tex $(??gC m^{-2} dt_slow^{-1})$ @endtex  
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_daily_x              !! "Daily" gpp for teststomate  
                                                                              !! @tex $(??gC m^{-2} dt_slow^{-1})$ @endtex 
   
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: resp_hetero_litter       !! Litter heterotrophic respiration per ground area 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex  
                                                                              !! ??Same variable is also used to 
                                                                              !! store heterotrophic respiration per ground area 
                                                                              !! over ::dtradia?? 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: resp_hetero_soil         !! soil heterotrophic respiration  
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex
!Chloe : 
    REAL(r_std),DIMENSION(kjpindex)           :: wpro_conso_carb   

    REAL(r_std),DIMENSION(kjpindex,nvm)           :: veget_cov                !! Fractional coverage: actually share of the pixel 
                                                                              !! covered by a PFT (fraction of ground area), 
                                                                              !! taking into account LAI ??(= grid scale fpc)?? 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: vcmax                    !! Maximum rate of carboxylation
                                                                              !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: vjmax                    !! Maximum rate of RUbp regeneration
                                                                              !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: t_photo_min              !! Min temperature for photosynthesis (K)
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: t_photo_opt              !! Opt temperature for photosynthesis (K)
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: t_photo_max              !! Max temperature for photosynthesis (K)
    REAL(r_std),DIMENSION(kjpindex,nlevs)         :: control_moist_inst       !! Moisture control of heterotrophic respiration 
                                                                              !! (0-1, unitless) 
    REAL(r_std),DIMENSION(kjpindex,nlevs)         :: control_temp_inst        !! Temperature control of heterotrophic 
                                                                              !! respiration, above and below (0-1, unitless) 
    REAL(r_std),DIMENSION(kjpindex,ncarb,nvm)     :: soilcarbon_input_inst    !! Quantity of carbon going into carbon pools from 
                                                                              !! litter decomposition 
                                                                              !! @tex $(gC m^{-2} day^{-1})$ @endtex 
    
    INTEGER(i_std)                                :: ier                      !! Check errors in netcdf call (unitless)
    REAL(r_std),SAVE                              :: dt_forcesoil             !! Time step of soil forcing file (days)
                                                                              !!?? Unit?
    INTEGER(i_std),PARAMETER                      :: nparanmax=366            !! Maximum number of time steps per year for 
                                                                              !! forcesoil 
    INTEGER(i_std),SAVE                           :: nparan                   !! Number of time steps per year for forcesoil - 
                                                                              !! read from run definition (unitless) 
    INTEGER(i_std),SAVE                           :: nbyear                   !! Number of years saved for forcesoil - read from 
                                                                              !! run definition (unitless) 
    INTEGER(i_std),SAVE                           :: iatt                     !! Time step of forcing of soil processes (iatt = 1 
                                                                              !! to ::nparan*::nbyear) 
    INTEGER(i_std),SAVE                           :: iatt_old=1               !! Previous ::iatt
    REAL(r_std)                                   :: sf_time                  !! Intermediate variable to calculate current time 
                                                                              !! step 
    INTEGER(i_std)                                :: max_totsize              !! Memory management - maximum memory size (Mb)
    INTEGER(i_std)                                :: totsize_1step            !! Memory management - memory required to store one 
                                                                              !! time step on one processor (Mb) 
    INTEGER(i_std)                                :: totsize_tmp              !! Memory management - memory required to store one 
                                                                              !! time step on all processors(Mb) 
    REAL(r_std)                                   :: xn                       !! How many times have we treated in this forcing 
                                                                              !! state 
    INTEGER(i_std),SAVE                           :: nsfm                     !! Number of time steps that can be stored in 
                                                                              !! memory (unitless) 
    INTEGER(i_std),SAVE                           :: nsft                     !! Number of time steps in a year (unitless)
 
   INTEGER(i_std),SAVE                            :: iisf                     !! Current pointer for teststomate (unitless)
    REAL(r_std), DIMENSION(kjpindex)              :: vartmp                   !! Temporary variable
    CHARACTER(LEN=100), SAVE                      :: forcing_name             !! Name of forcing file 1
    CHARACTER(LEN=100), SAVE                      :: Cforcing_name            !! Name of forcing file 2
    INTEGER(i_std),SAVE                           :: Cforcing_id              !! File identifer of file 2
    
    INTEGER(i_std),PARAMETER                      :: ndm = 10                 !! Maximum number of dimensions (unitless)
    INTEGER(i_std)                                :: vid                      !! Variable identifer of netCDF (unitless)
    INTEGER(i_std)                                :: nneigh                   !! Number of neighbouring pixels
    INTEGER(i_std)                                :: direct                   !! ??
    INTEGER(i_std),DIMENSION(ndm)                 :: d_id                     !! ??

    REAL(r_std),DIMENSION(kjpindex,nvm)           :: t_root                   !! ??[DISPENSABLE] Root temperature (K) convolution 
                                                                              !! of root and soil temperature profiles 
    REAL(r_std),DIMENSION(kjpindex,nvm,nparts)    :: coeff_maint              !! ??[DISPENSABLE]
    REAL(r_std),DIMENSION(kjpindex,nparts)        :: t_maint_radia            !! ??[DISPENSABLE] Temperature which is pertinent 
                                                                              !! for maintenance respiration (K) 
    REAL(r_std),DIMENSION(kjpindex)               :: rpc                      !! ??[DISPENSABLE] Integration constant for root 
                                                                              !! profile (unitless) 
    REAL(r_std),DIMENSION(kjpindex)               :: tl                       !! ??[DISPENSABLE] long term annual mean 
                                                                              !! temperature, C 
    REAL(r_std),DIMENSION(kjpindex)               :: slope                    !! ??[DISPENSABLE] slope of maintenance respiration 
                                                                              !! coefficient (1/K) 
    INTEGER(i_std)                                :: iyear                    !! ??[DISPENSABLE]
    REAL(r_std)                                   :: net_co2_flux_monthly     !! ??[DISPENSABLE]
    REAL(r_std)                                   :: net_co2_flux_monthly_sum !! ??[DISPENSABLE]
    INTEGER                                       :: ios                      !! ??[DISPENSABLE]
    REAL(r_std)                                   :: trans_veg                !! ??[DISPENSABLE]
    REAL(r_std)                                   :: tmp_day(1)               !! ??[DISPENSABLE]
    INTEGER(i_std),SAVE                           :: lcanop                   !! ??[DISPENSABLE] soil level used for LAI
 
    REAL(r_std),DIMENSION(nbp_glo)                :: clay_g                   !! Clay fraction of soil (0-1, unitless), parallel 
                                                                              !! computing 
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:,:)    :: soilcarbon_input_g       !! Quantity of carbon going into carbon pools from 
                                                                              !! litter decomposition  
                                                                              !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex, parallel 
                                                                              !! computing 
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)      :: control_moist_g          !! Moisture control of heterotrophic respiration 
                                                                              !! (0-1, unitless), parallel computing 
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:,:)      :: control_temp_g           !! Temperature control of heterotrophic respiration 
                                                                              !! (0-1, unitless), parallel computing 
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:)        :: npp_equil_g              !! Equilibrium NPP written to forcesoil 
                                                                              !! @tex $(gC m^{-2} year^{-1})$ @endtex, parallel 
                                                                              !! computing 

    REAL(r_std)                                   :: net_cflux_prod_monthly_sum    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL(r_std)                                   :: net_cflux_prod_monthly_tot    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL(r_std)                                   :: net_harvest_above_monthly_sum !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL(r_std)                                   :: net_harvest_above_monthly_tot !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL(r_std)                                   :: net_biosp_prod_monthly_sum    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
    REAL(r_std)                                   :: net_biosp_prod_monthly_tot    !! AR5 output?? gC m2 month-1 (one variable for 
                                                                                   !! reduce_sum and one for bcast??), parallel 
                                                                                   !! computing 
!_ ================================================================================================================================
    
  !! 1. Initialize variables

    !! 1.1 Store current time step in a common variable
    itime = kjit
    
    !![DISPENSABLE] 1.2 Copy the depth of the different soil layers from diaglev specified in slow_proc
															
    !! 1.3 PFT rooting depth across pixels, humescte is pre-defined 
    ! (constantes_veg.f90). It is defined as the coefficient of an exponential 
    ! function relating root density to depth 
    DO j=1,nvm
       rprof(:,j) = 1./humcste(j)
    ENDDO
    
    !! 1.4 Initialize first call
    ! Set growth respiration to zero
    resp_growth=zero
    IF (l_first_stomate) THEN
       IF (long_print) THEN
          WRITE (numout,*) ' l_first_stomate : call stomate_init'
       ENDIF
       
       !! 1.4.0 Initialization of PFT specific parameters 
       ! Initialization of PFT specific parameters that have no value
       ! for the bare soil PFT i.e. fire resistance, flamability, maximum lai,
       ! settings for growing degree days (GDD), settings for senescence, 
       ! respiration coefficients, photosynthesis, etc.
       ! [DISPENSABLE]
       
       !! 1.4.1 Allocate memory for all variables in stomate
       ! Allocate memory for all variables in stomate, build new index
       ! tables accounting for the PFTs, read and check flags and set file
       ! identifier for restart and history files.
       CALL stomate_init (kjpij, kjpindex, index, ldforcing_write, lalo, &
            rest_id_stom, hist_id_stom, hist_id_stom_IPCC)

       !! 1.4.2 Initialization of PFT specific parameters
       ! Initialization of PFT specific parameters i.e. sla from leaf life, 
       ! sapling characteristics (biomass), migration speed, critical diameter,
       ! coldest tolerable temperature, critical values for phenology, maximum
       ! life time of leaves, respiration coefficients and photosynthesis.
       ! The subroutine also communicates settings read by stomate_constant_init.
       CALL data (kjpindex, lalo)

       !! 1.4.3 Initial conditions
       
       !! 1.4.3.1 Read initial values for STOMATE's variables from the _restart_ file
       ! ??Shouldn't this be included in stomate_init?? Looks like an initialization!
       co2_flux(:,:) = zero
       fco2_lu(:) = zero
       
       ! Get values from _restart_ file. Note that only ::kjpindex, ::index, ::lalo 
       ! and ::resolution are input variables, all others are output variables.
       CALL readstart &
            &        (kjpindex, index, lalo, resolution, &
            &         day_counter_read, dt_days_read, date_read, &
            &         ind, adapted, regenerate, &
            &         humrel_daily, litterhum_daily, &
            &         t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
            &         soilhum_daily, precip_daily, &
            &         gpp_daily, npp_daily, turnover_daily, &
            &         humrel_month, humrel_week, &
            &         t2m_longterm, tlong_ref, t2m_month, t2m_week, &
            &         tsoil_month, soilhum_month, fireindex, firelitter, &
            &         maxhumrel_lastyear, maxhumrel_thisyear, &
            &         minhumrel_lastyear, minhumrel_thisyear, &
            &         maxgppweek_lastyear, maxgppweek_thisyear, maxnppweek_thisyear, &
            &         gdd0_lastyear, gdd0_thisyear, &
            &         precip_lastyear, precip_thisyear, &
            &         gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &         PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
            &         maxfpc_lastyear, maxfpc_thisyear, &
            &         turnover_longterm, gpp_week,npp_week, biomass, resp_maint_part, &
            &         leaf_age, leaf_frac, &
            &         senescence, when_growthinit, age, &
            &         resp_hetero_d, resp_maint_d, resp_growth_d, co2_fire, co2_to_bm_dgvm, &
            &         veget_lastlight, everywhere, need_adjacent, RIP_time, time_lowgpp, &
            &         time_hum_min, hum_min_dormance, &
            &         litterpart, litter, dead_leaves, &
            &         carbon, black_carbon, lignin_struc,turnover_time,&
            &         prod10,prod100,flux10, flux100, &
            &         convflux, cflux_prod10, cflux_prod100, bm_to_litter, carb_mass_total, & !)
!Chloe
            &         uo_peat, uold2_peat, tsurf_year, uo_0, uold2_0, wtsoil_daily)

       !! 1.4.4 If the vegetation is dynamic
       ! If the vegetation is dynamic, long-term reference temperature was 
       ! read by the call above. If vegetation is static then the long-term 
       ! reference temperature is a boundary condition and read here.  
       CALL readbc (kjpindex, lalo, resolution, tlong_ref)
       
       !! 1.4.5 Check time step
       
       !! 1.4.5.1 Allow STOMATE's time step to change although this is dangerous
       IF (dt_days /= dt_days_read) THEN
          WRITE(numout,*) 'slow_processes: STOMATE time step changes:', &
           & dt_days_read,' -> ',dt_days
       ENDIF

       !! 1.4.5.2 Time step has to be a multiple of a full day
       IF ( ( dt_days-REAL(NINT(dt_days),r_std) ) > min_stomate ) THEN
          WRITE(numout,*) 'slow_processes: STOMATE time step is not a mutiple of a full day:', &
           & dt_days,' days.'
          STOP
       ENDIF

       !! 1.4.5.3 upper limit to STOMATE's time step
       IF ( dt_days > max_dt_days ) THEN
          WRITE(numout,*) 'slow_processes: STOMATE time step exceeds the maximum value:', &
           & dt_days,' days > ', max_dt_days, ' days.'  
          STOP
       ENDIF

       !! 1.4.5.4 STOMATE time step must not be less than the forcing time step
       IF ( dtradia > dt_days*one_day ) THEN
          WRITE(numout,*) &
            & 'slow_processes: STOMATE time step ::dt_days smaller than forcing time step ::dtradia'
          STOP
       ENDIF

       !! 1.4.5.5 For teststomate : test day_counter
       IF ( abs(day_counter - day_counter_read) > min_stomate ) THEN
          WRITE(numout,*) 'slow_processes: STOMATE day counter changes:', &
               day_counter_read,' -> ',day_counter
       ENDIF

       !! 1.4.5.6 Date check
       IF (date /= date_read) THEN
          WRITE(numout,*) 'Slow_processes: STOMATE date changes:', &
               date_read,' -> ',date
       ENDIF

       !! 1.4.5.6 Final message on time step
       WRITE(numout,*) 'Slow_processes, STOMATE time step (days): ', dt_days

       !! 1.4.6 Write forcing file for stomate?
       IF (ldforcing_write) THEN

          !Config Key   = STOMATE_FORCING_NAME
          !Config Desc  = Name of STOMATE's forcing file
          !Config If    = OK_STOMATE
          !Config Def   = NONE
          !Config Help  = Name that will be given
          !Config         to STOMATE's offline forcing file
          !Config         Compatible with Nicolas Viovy's driver
          !Config Units = [FILE]
          forcing_name = stomate_forcing_name
          CALL getin_p('STOMATE_FORCING_NAME',forcing_name)

          IF ( TRIM(forcing_name) /= 'NONE' ) THEN
             
             !! 1.4.6.1 Calculate steps that can be stored in memory
             ! Action for the root processor only (parallel computing)  
             IF (is_root_prc) CALL SYSTEM ('rm -f '//TRIM(forcing_name))
             WRITE(numout,*) 'writing a forcing file for STOMATE.'

             !Config Key   = STOMATE_FORCING_MEMSIZE
             !Config Desc  = Size of STOMATE forcing data in memory 
             !Config If    = OK_STOMATE
             !Config Def   = 50
             !Config Help  = This variable determines how many
             !Config         forcing states will be kept in memory.
             !Config         Must be a compromise between memory
             !Config         use and frequeny of disk access.
             !Config Units = [MegaBytes]
             max_totsize = 50
             CALL getin_p('STOMATE_FORCING_MEMSIZE', max_totsize)      
             max_totsize = max_totsize*1000000

             totsize_1step = &
                  &      SIZE(clay)*KIND(clay) &
                  &           +SIZE(humrel_daily)*KIND(humrel_daily) &
                  &     +SIZE(litterhum_daily)*KIND(litterhum_daily) &
                  &     +SIZE(t2m_daily)*KIND(t2m_daily) &
                  &     +SIZE(t2m_min_daily)*KIND(t2m_min_daily) &
                  &     +SIZE(tsurf_daily)*KIND(tsurf_daily) &
                  &     +SIZE(tsoil_daily)*KIND(tsoil_daily) &
                  &     +SIZE(soilhum_daily)*KIND(soilhum_daily) &
                  &     +SIZE(precip_daily)*KIND(precip_daily) &
                  &     +SIZE(gpp_daily_x)*KIND(gpp_daily_x) &
                  &     +SIZE(veget)*KIND(veget) &
                  &     +SIZE(veget_max)*KIND(veget_max) &
                  &     +SIZE(lai)*KIND(lai)
             
             ! Totsize_1step is the size on a single processor, sum
             ! all processors and send to all processors
             CALL reduce_sum(totsize_1step,totsize_tmp)
             CALL bcast(totsize_tmp)
             totsize_1step=totsize_tmp
 
             ! Total number of forcing steps
             nsft = INT(one_year/(dt_slow/one_day))

             ! Number of forcing steps in memory
             nsfm = MIN(nsft, &
                  &       MAX(1,NINT( REAL(max_totsize,r_std) &
                  &                  /REAL(totsize_1step,r_std))))
            
             
	     !! 1.6.4.2 Allocate memory for variables containing forcing data  
             ! and initialize variables (set to zero).
             CALL init_forcing (kjpindex,nsfm,nsft)
             
             ! Indexing for writing forcing file
             isf(:) = (/ (i,i=1,nsfm) /)
             nf_written(:) = .FALSE.
             nf_cumul(:) = 0
             iisf = 0

             !! 1.6.4.3 Create netcdf file
             ! Create, define and populate a netcdf file containing the forcing data.
             ! For the root processor only (parallel computing). NF90_ are functions
             ! from and external library.  
             IF (is_root_prc) THEN

                ! Create new netCDF dataset
                ier = NF90_CREATE (TRIM(forcing_name),NF90_SHARE,forcing_id)

                ! Add variable attribute
                ! Note ::iim_g and ::jjm_g are dimensions of the global field and 
                ! ::nbp_glo is the number of global continental points
                ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL,'dtradia',dtradia)
                ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL,'dt_slow',dt_slow)
                ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL, &
                     & 'nsft',REAL(nsft,r_std))
                ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL, &
                     & 'kjpij',REAL(iim_g*jjm_g,r_std))
                ier = NF90_PUT_ATT (forcing_id,NF90_GLOBAL, &
                     & 'kjpindex',REAL(nbp_glo,r_std))
 
                ! Add new dimension
                ier = NF90_DEF_DIM (forcing_id,'points',nbp_glo,d_id(1))
                ier = NF90_DEF_DIM (forcing_id,'layers',nbdl,d_id(2))
                ier = NF90_DEF_DIM (forcing_id,'pft',nvm,d_id(3))
                direct=2
                ier = NF90_DEF_DIM (forcing_id,'direction',direct,d_id(4))
                nneigh=8
                ier = NF90_DEF_DIM (forcing_id,'nneigh',nneigh,d_id(5))
                ier = NF90_DEF_DIM (forcing_id,'time',nsft,d_id(6))
                ier = NF90_DEF_DIM (forcing_id,'nbparts',nparts,d_id(7))

                ! Add new variable
                ier = NF90_DEF_VAR (forcing_id,'points',    r_typ,d_id(1),vid)
                ier = NF90_DEF_VAR (forcing_id,'layers',    r_typ,d_id(2),vid)
                ier = NF90_DEF_VAR (forcing_id,'pft',       r_typ,d_id(3),vid)
                ier = NF90_DEF_VAR (forcing_id,'direction', r_typ,d_id(4),vid)
                ier = NF90_DEF_VAR (forcing_id,'nneigh',    r_typ,d_id(5),vid)
                ier = NF90_DEF_VAR (forcing_id,'time',      r_typ,d_id(6),vid)
                ier = NF90_DEF_VAR (forcing_id,'nbparts',   r_typ,d_id(7),vid)
                ier = NF90_DEF_VAR (forcing_id,'index',     r_typ,d_id(1),vid)
                ier = NF90_DEF_VAR (forcing_id,'contfrac',  r_typ,d_id(1),vid) 
                ier = NF90_DEF_VAR (forcing_id,'lalo', &
                     & r_typ,(/ d_id(1),d_id(4) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'neighbours', &
                     & r_typ,(/ d_id(1),d_id(5) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'resolution', &
                     & r_typ,(/ d_id(1),d_id(4) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'clay', &
                     & r_typ,(/ d_id(1),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'humrel', &
                     & r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'litterhum', &
                     & r_typ,(/ d_id(1),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'t2m', &
                     & r_typ,(/ d_id(1),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'t2m_min', &
                     & r_typ,(/ d_id(1),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'tsurf', &
                     & r_typ,(/ d_id(1),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'tsoil', &
                     & r_typ,(/ d_id(1),d_id(2),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'soilhum', &
                     & r_typ,(/ d_id(1),d_id(2),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'precip', &
                     & r_typ,(/ d_id(1),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'gpp', &
                     & r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'veget', &
                     & r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'veget_max', &
                     & r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
                ier = NF90_DEF_VAR (forcing_id,'lai', &
                     & r_typ,(/ d_id(1),d_id(3),d_id(6) /),vid)
                ier = NF90_ENDDEF (forcing_id)
                
		! Given the name of a varaible, nf90_inq_varid finds the variable 
                ! ID (::vid). Put data value(s) into variable ::vid
                ier = NF90_INQ_VARID (forcing_id,'points',vid)
                ier = NF90_PUT_VAR (forcing_id,vid, &
                     & (/(REAL(i,r_std),i=1,nbp_glo) /))
                ier = NF90_INQ_VARID (forcing_id,'layers',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nbdl)/))
                ier = NF90_INQ_VARID (forcing_id,'pft',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nvm)/))
                ier = NF90_INQ_VARID (forcing_id,'direction',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,2)/))
                ier = NF90_INQ_VARID (forcing_id,'nneigh',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,8)/))
                ier = NF90_INQ_VARID (forcing_id,'time',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nsft)/))
                ier = NF90_INQ_VARID (forcing_id,'nbparts',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,(/(REAL(i,r_std),i=1,nparts)/))
                ier = NF90_INQ_VARID (forcing_id,'index',vid)  
                ier = NF90_PUT_VAR (forcing_id,vid,REAL(index_g,r_std))
                ier = NF90_INQ_VARID (forcing_id,'contfrac',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,REAL(contfrac_g,r_std))
                ier = NF90_INQ_VARID (forcing_id,'lalo',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,lalo_g)
                !ym attention a neighbours, a modifier plus tard      
                ier = NF90_INQ_VARID (forcing_id,'neighbours',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,REAL(neighbours_g,r_std))
                ier = NF90_INQ_VARID (forcing_id,'resolution',vid)
                ier = NF90_PUT_VAR (forcing_id,vid,resolution_g)
             ENDIF ! is_root_prc
          ENDIF ! (forcing_name) /= 'NONE'
       ENDIF ! ldforcing_write =.TRUE.
       
       !! 1.4.7 write forcing file for the soil?
       IF (ldcarbon_write) THEN
          
          !! 1.4.7.1 Initialize
          !Config Key   = STOMATE_CFORCING_NAME
          !Config Desc  = Name of STOMATE's carbon forcing file
          !Config If    = OK_STOMATE
          !Config Def   = NONE
          !Config Help  = Name that will be given to STOMATE's carbon
          !Config         offline forcing file
          !Config         Compatible with Nicolas Viovy's driver
          !Config Units = [FILE]
          Cforcing_name = stomate_Cforcing_name
          CALL getin_p('STOMATE_CFORCING_NAME',Cforcing_name)

          IF ( TRIM(Cforcing_name) /= 'NONE' ) THEN
  
             ! For root processor only (parallel computing)
             IF (is_root_prc) CALL SYSTEM ('rm -f '//TRIM(Cforcing_name))
             
             ! Time step of forcesoil
             !Config Key   = FORCESOIL_STEP_PER_YEAR
             !Config Desc  = Number of time steps per year for carbon spinup.
             !Config If    = OK_STOMATE
             !Config Def   = 365
             !Config Help  = Number of time steps per year for carbon spinup.
             !Config Units = [days, months, year]
             nparan = 365
             CALL getin_p('FORCESOIL_STEP_PER_YEAR', nparan)
             
             ! Correct if setting is out of bounds 
             IF ( nparan < 1 ) nparan = 1

             !Config Key   = FORCESOIL_NB_YEAR
             !Config Desc  = Number of years saved for carbon spinup.
             !Config If    = OK_STOMATE
             !Config Def   = 1
             !Config Help  = Number of years saved for carbon spinup. If internal parameter cumul_Cforcing is TRUE in stomate.f90
             !Config         Then this parameter is forced to one.
             !Config Units = [years]
             nbyear=1
             CALL getin_p('FORCESOIL_NB_YEAR', nbyear)
             
             ! Set ::nbyear to 1. if ::cumul_Cforcing=.TRUE.
             IF ( cumul_Cforcing ) THEN
                CALL ipslerr (1,'stomate', &
                     &          'Internal parameter cumul_Cforcing is TRUE in stomate.f90', &
                     &          'Parameter FORCESOIL_NB_YEAR is therefore forced to 1.', &
                     &          '::nbyear is thus set to 1.')
                nbyear=1
             ENDIF

             ! Make use of ::nparan to calculate ::dt_forcesoil
             dt_forcesoil = zero
             nparan = nparan+1
             DO WHILE ( dt_forcesoil < dt_slow/one_day )
                nparan = nparan-1
                IF ( nparan < 1 ) THEN
                   STOP 'Problem with number of soil forcing time steps ::nparan < 1.'
                ENDIF
                dt_forcesoil = one_year/REAL(nparan,r_std)
             ENDDO
             IF ( nparan > nparanmax ) THEN
                STOP 'Problem with number of soil forcing time steps ::nparan > ::nparanmax'
             ENDIF
             WRITE(numout,*) 'Time step of soil forcing (d): ',dt_forcesoil

             ! Allocate memory for the forcing variables of soil dynamics
             ALLOCATE( nforce(nparan*nbyear))
             nforce(:) = 0
             ALLOCATE(control_moist(kjpindex,nlevs,nparan*nbyear))
             ALLOCATE(npp_equil(kjpindex,nparan*nbyear))
             ALLOCATE(npp_tot(kjpindex))
             ALLOCATE(control_temp(kjpindex,nlevs,nparan*nbyear))
             ALLOCATE(soilcarbon_input(kjpindex,ncarb,nvm,nparan*nbyear)) 
             
             ! Initialize variables, set to zero
             control_moist(:,:,:) = zero
             npp_equil(:,:) = zero
             npp_tot(:) = zero
             control_temp(:,:,:) = zero
             soilcarbon_input(:,:,:,:) = zero

          ENDIF ! Cforcing_name) /= 'NONE'
       ENDIF ! ::ldcarbon_write!.TRUE.
       
       !! 1.4.8 Calculate STOMATE's vegetation fractions from veget, veget_max
       DO j=1,nvm
          WHERE ((1.-totfrac_nobio(:)) > min_sechiba)       
             ! Pixels with vegetation
             veget_cov(:,j) = veget(:,j)/( 1.-totfrac_nobio(:) )
             veget_cov_max(:,j) = veget_max(:,j)/( 1.-totfrac_nobio(:) )
          ELSEWHERE
             ! Pixels without vegetation
             veget_cov(:,j) = zero
             veget_cov_max(:,j) = zero
          ENDWHERE
       ENDDO ! Loop over PFTs
       IF ( EndOfYear ) THEN
          DO j=1,nvm
             WHERE ((1.-totfrac_nobio_new(:)) > min_sechiba)
                ! Update pixels with vegetation
                veget_cov_max_new(:,j) = veget_max_new(:,j)/( 1.-totfrac_nobio_new(:) )
             ELSEWHERE
                ! Update pixels without vegetation
                veget_cov_max_new(:,j) = zero
             ENDWHERE
          ENDDO
       ENDIF
       
       !! 1.4.9 Initialize non-zero variables
       IF (control%ok_stomate) THEN
          CALL stomate_var_init &
               &         (kjpindex, veget_cov_max, leaf_age, leaf_frac, &
               &          tlong_ref, t2m_month, dead_leaves, &
               &          veget, lai, deadleaf_cover, assim_param)
          
          ! Initialize land cover change variable
          ! ??Should be integrated in the subroutine?? 
          harvest_above = zero
       ENDIF
      
       !! 1.4.10 Update flag
       l_first_stomate = .FALSE.
      
       RETURN ! Out of stomate main
    ENDIF  ! First call


    IF (bavard >= 4) THEN
       WRITE(numout,*) 'DATE ',date,' ymds', year, month, day, sec, '-- stp --', itime, do_slow
    ENDIF


  !! 2. Prepares restart file for the next simulation

    !! 2.1 Write restart file for stomate
    IF (ldrestart_write) THEN
       IF (long_print) THEN
          WRITE (numout,*) &
               &      ' Write a restart file for STOMATE'
       ENDIF
       CALL writerestart &
            &         (kjpindex, index, &
            &          day_counter, dt_days, date, &
            &          ind, adapted, regenerate, &
            &          humrel_daily, litterhum_daily, &
            &          t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
            &          soilhum_daily, precip_daily, &
            &          gpp_daily, npp_daily, turnover_daily, &
            &          humrel_month, humrel_week, &
            &          t2m_longterm, tlong_ref, t2m_month, t2m_week, &
            &          tsoil_month, soilhum_month, fireindex, firelitter, &
            &          maxhumrel_lastyear, maxhumrel_thisyear, &
            &          minhumrel_lastyear, minhumrel_thisyear, &
            &          maxgppweek_lastyear, maxgppweek_thisyear,maxnppweek_thisyear, &
            &          gdd0_lastyear, gdd0_thisyear, &
            &          precip_lastyear, precip_thisyear, &
            &          gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &          PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
            &          maxfpc_lastyear, maxfpc_thisyear, &
            &          turnover_longterm, gpp_week,npp_week, biomass, resp_maint_part, &
            &          leaf_age, leaf_frac, &
            &          senescence, when_growthinit, age, &
            &          resp_hetero_d, resp_maint_d, resp_growth_d, co2_fire, co2_to_bm_dgvm, &
            &          veget_lastlight, everywhere, need_adjacent, &
            &          RIP_time, time_lowgpp, &
            &          time_hum_min, hum_min_dormance, &
            &          litterpart, litter, dead_leaves, &
            &          carbon, black_carbon, lignin_struc,turnover_time,&
            &          prod10,prod100,flux10, flux100, &
            &          convflux, cflux_prod10, cflux_prod100, bm_to_litter, carb_mass_total, & !)
! Chloe
            &          uo_peat, uold2_peat, tsurf_year, uo_0, uold2_0, wtsoil_daily) !Chloe--
       
       !! 2.2 Write file with variables that force general processes in stomate
       IF (ldforcing_write .AND. TRIM(forcing_name) /= 'NONE' ) THEN  
          CALL forcing_write(forcing_id,1,iisf)
          ! Close forcing file
          IF (is_root_prc) ier = NF90_CLOSE (forcing_id)
          forcing_id=-1
       ENDIF

       !! 2.3 Collect variables that force the soil processes in stomate
       IF (ldcarbon_write .AND. TRIM(Cforcing_name) /= 'NONE' ) THEN 
          
          !! 2.3.1 Collet variables 
          WRITE(numout,*) &
               &      'stomate: writing the forcing file for carbon spinup'
          DO iatt = 1, nparan*nbyear
             IF ( nforce(iatt) > 0 ) THEN
                soilcarbon_input(:,:,:,iatt) = &
                     & soilcarbon_input(:,:,:,iatt)/REAL(nforce(iatt),r_std)
                control_moist(:,:,iatt) = &
                     & control_moist(:,:,iatt)/REAL(nforce(iatt),r_std)
                control_temp(:,:,iatt) = &
                     & control_temp(:,:,iatt)/REAL(nforce(iatt),r_std)
                npp_equil(:,iatt) = &
                     & npp_equil(:,iatt)/REAL(nforce(iatt),r_std)
             ELSE
                WRITE(numout,*) &
                     &         'We have no soil carbon forcing data for this time step:', &
                     &         iatt
                WRITE(numout,*) ' -> we set them to zero'
                soilcarbon_input(:,:,:,iatt) = zero
                control_moist(:,:,iatt) = zero
                control_temp(:,:,iatt) = zero
                npp_equil(:,iatt) = zero
             ENDIF
          ENDDO

          ! Allocate memory for parallel computing
          IF (is_root_prc) THEN
             ALLOCATE(soilcarbon_input_g(nbp_glo,ncarb,nvm,nparan*nbyear))
             ALLOCATE(control_moist_g(nbp_glo,nlevs,nparan*nbyear))
             ALLOCATE(control_temp_g(nbp_glo,nlevs,nparan*nbyear))
             ALLOCATE(npp_equil_g(nbp_glo,nparan*nbyear))
          ENDIF
          
          ! Gather distributed variables
          CALL gather(clay,clay_g)
          CALL gather(soilcarbon_input,soilcarbon_input_g)
          CALL gather(control_moist,control_moist_g)
          CALL gather(control_temp,control_temp_g)
          CALL gather(npp_equil,npp_equil_g)
          
          !! 2.3.2 Create netcdf
          ! Create, define and populate a netcdf file containing the forcing data.
          ! For the root processor only (parallel computing). NF90_ are functions
          ! from and external library.  
          IF (is_root_prc) THEN

             ! Create new netCDF dataset
             ier = NF90_CREATE (TRIM(Cforcing_name),NF90_WRITE,Cforcing_id)

             ! Add variable attribute
             ! Note ::nbp_glo is the number of global continental points
             ier = NF90_PUT_ATT (Cforcing_id,NF90_GLOBAL, &
                  &                        'kjpindex',REAL(nbp_glo,r_std))
             ier = NF90_PUT_ATT (Cforcing_id,NF90_GLOBAL, &
                  &                        'nparan',REAL(nparan,r_std))
             ier = NF90_PUT_ATT (Cforcing_id,NF90_GLOBAL, &
                  &                        'nbyear',REAL(nbyear,r_std))
             
             ! Add new dimension
             ier = NF90_DEF_DIM (Cforcing_id,'points',nbp_glo,d_id(1))
             ier = NF90_DEF_DIM (Cforcing_id,'carbtype',ncarb,d_id(2))
             ier = NF90_DEF_DIM (Cforcing_id,'vegtype',nvm,d_id(3))
             ier = NF90_DEF_DIM (Cforcing_id,'level',nlevs,d_id(4))
             ier = NF90_DEF_DIM (Cforcing_id,'time_step',nparan*nbyear,d_id(5))
             
             ! Add new variable
             ier = NF90_DEF_VAR (Cforcing_id,'points',    r_typ,d_id(1),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'carbtype',  r_typ,d_id(2),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'vegtype',   r_typ,d_id(3),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'level',     r_typ,d_id(4),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'time_step', r_typ,d_id(5),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'index',     r_typ,d_id(1),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'clay',      r_typ,d_id(1),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'soilcarbon_input',r_typ, &
                  &                        (/ d_id(1),d_id(2),d_id(3),d_id(5) /),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'control_moist',r_typ, &
                  &                        (/ d_id(1),d_id(4),d_id(5) /),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'control_temp',r_typ, &
                  &                        (/ d_id(1),d_id(4),d_id(5) /),vid)
             ier = NF90_DEF_VAR (Cforcing_id,'npp_equil',r_typ, &
                  &                        (/ d_id(1),d_id(5) /),vid)
             ier = NF90_ENDDEF (Cforcing_id)
             
             ! Given the name of a varaible, nf90_inq_varid finds the variable 
             ! ID (::vid). Put data value(s) into variable ::vid 
             ier = NF90_INQ_VARID (Cforcing_id,'points',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, &
                  &                          (/(REAL(i,r_std),i=1,nbp_glo)/))
             ier = NF90_INQ_VARID (Cforcing_id,'carbtype',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, &
                  &                        (/(REAL(i,r_std),i=1,ncarb)/))
             ier = NF90_INQ_VARID (Cforcing_id,'vegtype',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, &
                  &                            (/(REAL(i,r_std),i=1,nvm)/))
             ier = NF90_INQ_VARID (Cforcing_id,'level',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, &
                  &                          (/(REAL(i,r_std),i=1,nlevs)/))
             ier = NF90_INQ_VARID (Cforcing_id,'time_step',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, &
                  &                          (/(REAL(i,r_std),i=1,nparan*nbyear)/))
             ier = NF90_INQ_VARID (Cforcing_id,'index',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, REAL(index_g,r_std) )
             ier = NF90_INQ_VARID (Cforcing_id,'clay',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, clay_g )
             ier = NF90_INQ_VARID (Cforcing_id,'soilcarbon_input',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, soilcarbon_input_g )
             ier = NF90_INQ_VARID (Cforcing_id,'control_moist',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, control_moist_g )
             ier = NF90_INQ_VARID (Cforcing_id,'control_temp',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, control_temp_g )
             ier = NF90_INQ_VARID (Cforcing_id,'npp_equil',vid)
             ier = NF90_PUT_VAR (Cforcing_id,vid, npp_equil_g )
             
             ! Close netCDF
             ier = NF90_CLOSE (Cforcing_id)
             Cforcing_id = -1
          ENDIF

          ! Clear memory
          IF (is_root_prc) THEN
             DEALLOCATE(soilcarbon_input_g)
             DEALLOCATE(control_moist_g)
             DEALLOCATE(control_temp_g)
             DEALLOCATE(npp_equil_g)
          ENDIF
       
       ENDIF ! ldcarbon_write=.TRUE.
       RETURN ! This routine is only called at the last time step to write a netcdf for forcesoil, 
              ! the rest of the code does not need to be executed.
    ENDIF  ! write restart-files
    
!! 3. Special treatment for some input arrays.
    
    !! 3.1 Sum of liquid and solid precipitation
    precip = ( precip_rain+precip_snow )*one_day/dtradia
    
    !! 3.2 Calculate STOMATE's vegetation fractions from veget and veget_max
    DO j=1,nvm 
       WHERE ((1.-totfrac_nobio(:)) > min_sechiba)
          ! Pixels with vegetation
          veget_cov(:,j) = veget(:,j)/( 1.-totfrac_nobio(:) )
          veget_cov_max(:,j) = veget_max(:,j)/( 1.-totfrac_nobio(:) )
       ELSEWHERE
          ! Pixels without vegetation
          veget_cov(:,j) = zero
          veget_cov_max(:,j) = zero
       ENDWHERE
    ENDDO
    IF ( EndOfYear ) THEN
       DO j=1,nvm
          WHERE ((1.-totfrac_nobio(:)) > min_sechiba)
             ! Pixels with vegetation
             veget_cov_max_new(:,j)=veget_max_new(:,j)/( 1.-totfrac_nobio(:) )
          ELSEWHERE
             ! Pixels without vegetation
             veget_cov_max_new(:,j) = zero
          ENDWHERE
       ENDDO
    ENDIF

    !! 3.3 Adjust time step of GPP 
    ! No GPP for bare soil
    gpp_d(:,1) = zero
    ! GPP per PFT
    DO j = 2,nvm   
       WHERE (veget_cov_max(:,j) > min_stomate)
          ! The PFT is available on the pixel
          gpp_d(:,j) =  gpp(:,j)/ veget_cov_max(:,j)* one_day/dtradia  
       ELSEWHERE
          ! The PFT is absent on the pixel
          gpp_d(:,j) = zero
       ENDWHERE
    ENDDO

    !! 3.4 The first time step of the first day of the month  
    ! implies that the month is over
    IF ( day == 1 .AND. sec .LT. dtradia ) THEN
       EndOfMonth=.TRUE.
    ELSE
       EndOfMonth=.FALSE.
    ENDIF
    

  !! 4. Calculate variables for dt_slow (i.e. "daily")

    ! Note: If dt_days /= 1, then variables 'xx_daily' (eg. half-daily or bi-daily) are by definition
    ! not expressed on a daily basis. This is not a problem but could be
    ! confusing

    !! 4.1 Accumulate instantaneous variables (do_slow=.FALSE.) 
    ! Accumulate instantaneous variables (do_slow=.FALSE.) and eventually 
    ! calculate daily mean value (do_slow=.TRUE.) 
    CALL stomate_accu (kjpindex, nvm, dt_slow, dtradia, &
         & do_slow, humrel, humrel_daily)
    CALL stomate_accu (kjpindex,    1, dt_slow, dtradia, &
         & do_slow, litterhumdiag, litterhum_daily)
    CALL stomate_accu (kjpindex,    1, dt_slow, dtradia, &
         & do_slow, t2m, t2m_daily)
    CALL stomate_accu (kjpindex,    1, dt_slow, dtradia, &
         & do_slow, temp_sol, tsurf_daily)
    CALL stomate_accu (kjpindex, nbdl, dt_slow, dtradia, &
         & do_slow, stempdiag, tsoil_daily)
    CALL stomate_accu (kjpindex, nbdl, dt_slow, dtradia, &
         & do_slow, shumdiag, soilhum_daily)
    CALL stomate_accu (kjpindex,    1, dt_slow, dtradia, &
         & do_slow, precip, precip_daily)
    CALL stomate_accu (kjpindex, nvm, dt_slow, dtradia, &
         & do_slow, gpp_d, gpp_daily)
!! Chloe wt_soil mean daily :
!! On utilise wt_soil2 pour faire la moyenne journaliere
    CALL stomate_accu (kjpindex, nstm, dt_slow, dtradia, &
         & do_slow, wt_soil2, wtsoil_daily) 

!    CALL getin_p('WTSOIL_DAILY', wtsoil_daily)

    !! 4.2 Daily minimum temperature
    t2m_min_daily(:) = MIN( t2m_min(:), t2m_min_daily(:) )

    !! 4.3 Calculate maintenance respiration
    ! Note: lai is passed as input argument to overcome previous problems with 
    ! natural and agricultural vegetation types. 
    CALL maint_respiration &
         & (kjpindex,dtradia,lai,t2m,tlong_ref,stempdiag,height,veget_cov_max, &
         & rprof,biomass,resp_maint_part_radia)
    
    ! Aggregate maintenance respiration across the different plant parts 
    resp_maint_radia(:,:) = zero
    DO j=2,nvm
       DO k= 1, nparts
          resp_maint_radia(:,j) = resp_maint_radia(:,j) &
               & + resp_maint_part_radia(:,j,k)
       ENDDO
    ENDDO
    
    ! Maintenance respiration separated by plant parts
    resp_maint_part(:,:,:) = resp_maint_part(:,:,:) &
         & + resp_maint_part_radia(:,:,:)
    
    !! 4.4 Litter dynamics and litter heterothropic respiration 
    ! Including: litter update, lignin content, PFT parts, litter decay,
    ! litter heterotrophic respiration, dead leaf soil cover.
    ! Note: there is no vertical discretisation in the soil for litter decay.
    turnover_littercalc = turnover_daily * dtradia/one_day
    bm_to_littercalc    = bm_to_litter*dtradia/one_day       
    CALL littercalc (kjpindex, dtradia/one_day, &
         turnover_littercalc, bm_to_littercalc, &
         veget_cov_max, temp_sol, stempdiag, shumdiag, litterhumdiag, &
         litterpart, litter, dead_leaves, lignin_struc, &
         deadleaf_cover, resp_hetero_litter, &
         soilcarbon_input_inst, control_temp_inst,&
         control_moist_inst,wpro_conso_carb)
         
    
    ! Heterothropic litter respiration during time step ::dtradia @tex $(gC m^{-2})$ @endtex
    resp_hetero_litter=resp_hetero_litter*dtradia/one_day
    
    !! 4.5 Soil carbon dynamics and soil heterotrophic respiration
    ! Note: there is no vertical discretisation in the soil for litter decay.
    CALL soilcarbon (kjpindex, dtradia/one_day, clay, &
         soilcarbon_input_inst, control_temp_inst, control_moist_inst, &
         carbon, resp_hetero_soil)
    
    ! Heterothropic soil respiration during time step ::dtradia @tex $(gC m^{-2})$ @endtex 
    resp_hetero_soil=resp_hetero_soil*dtradia/one_day

    ! Total heterothrophic respiration during time step ::dtradia @tex $(gC m^{-2})$ @endtex
    resp_hetero_radia = resp_hetero_litter+resp_hetero_soil
    resp_hetero_d = resp_hetero_d + resp_hetero_radia
    
    !! 4.6 Accumulate instantaneous variables (do_slow=.FALSE.) 
    ! Accumulate instantaneous variables (do_slow=.FALSE.) and eventually 
    ! calculate daily mean value (do_slow=.TRUE.) 
    CALL stomate_accu (kjpindex, nlevs, dt_slow, dtradia, &
     & do_slow, control_moist_inst, control_moist_daily)
    CALL stomate_accu (kjpindex, nlevs, dt_slow, dtradia, &
     & do_slow, control_temp_inst, control_temp_daily)
    DO i=1,ncarb
       CALL stomate_accu (kjpindex, nvm, dt_slow, dtradia, &
            & do_slow, soilcarbon_input_inst(:,i,:), soilcarbon_input_daily(:,i,:))
    ENDDO

!! 5. Daily processes - performed at the end of the day
    
    IF (do_slow) THEN
       ! WRITE(*,*) "Number of days: ", date !   Chloe : nombre de jour
!5.0 Chloe peat :
       !appel routines pour calcul des densites de flux de CH4

       CH4_calcul  = .FALSE.
       CALL getin_p('CH4_CALCUL', CH4_calcul)
       IF ( .NOT. control%ok_stomate ) CH4_calcul = .FALSE. ! GK08032012

       CH4atmo_CONC=0.0033
       CALL getin_p('CH4atmo_CONC', CH4atmo_CONC)

       CH4_WTD1  = .TRUE.
       CALL getin_p('CH4_WTD1', CH4_WTD1)

       IF(CH4_calcul) THEN
         !Chloe test ter0 : a supprimer ensuite :  
           CALL ch4_wet_flux_density_0 (kjpindex,dtradia,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
                    & carbon,lai,uo_0,uold2_0, ch4_flux_density_tot_0, ch4_flux_density_dif_0,&
                    & ch4_flux_density_bub_0,ch4_flux_density_pla_0, CH4atmo_CONC)
          !Chloe test fin


        IF (CH4_WTD1) THEN
!routine calcule densite de flux d un wetland ou WTD = pwater_peat (cf.stomate_cste_wetlands.f90) 
        CALL ch4_wet_flux_density_peat (kjpindex,dtradia,stempdiag,tsurf_daily,tsurf_year,veget_cov_max,veget,&
                        & carbon,lai,uo_peat,uold2_peat,ch4_flux_density_tot_peat, ch4_flux_density_dif_peat, &
                        & ch4_flux_density_bub_peat,ch4_flux_density_pla_peat,ch4_flux_density_goxid_peat, CH4atmo_CONC, & 
                        & wpro_conso_carb, nstm, wtsoil_daily ,npp_longterm, &
                        & npp_week,maxnppweek_thisyear, wtold)
          ELSE
             ch4_flux_density_tot_peat=0.0
             ch4_flux_density_dif_peat=0.0
             ch4_flux_density_bub_peat=0.0
             ch4_flux_density_pla_peat=0.0
             ch4_flux_density_goxid_peat=0.0
             wpro_conso_carb=0.0
        
         ENDIF

       ELSE

          ch4_flux_density_tot_peat=0.0
          ch4_flux_density_dif_peat=0.0
          ch4_flux_density_bub_peat=0.0
          ch4_flux_density_pla_peat=0.0
          ch4_flux_density_goxid_peat=0.0
          wpro_conso_carb=0.0
           
          ch4_flux_density_tot_0=0.0
          ch4_flux_density_dif_0=0.0
          ch4_flux_density_bub_0=0.0
          ch4_flux_density_pla_0=0.0


       ENDIF
!!!Chloe--

 !! 5.1 Update lai
       ! Use lai from stomate
       ! ?? check if this is the only time control%ok_pheno is used??
       ! ?? Looks like it is the only time. But this variables probably is defined 
       ! in stomate_constants or something, in which case, it is difficult to track.
       IF (control%ok_pheno) THEN
          !! 5.1.1 Update LAI 
          ! Set lai of bare soil to zero
          lai(:,ibare_sechiba) = zero
          ! lai for all PFTs
          DO j = 2, nvm
             lai(:,j) = biomass(:,j,ileaf)*sla(j)
          ENDDO
          frac_age(:,:,:) = leaf_frac(:,:,:)
       ELSE 
          ! 5.1.2 Use a prescribed lai
          ! WARNING: code in setlai is identical to the lines above
          ! Update subroutine if LAI has to be forced 
          CALL  setlai(kjpindex,lai) 
          frac_age(:,:,:) = zero
       ENDIF

       !! 5.2 Calculate long-term "meteorological" and biological parameters
       ! mainly in support of calculating phenology. If ::EndOfYear=.TRUE.
       ! annual values are update (i.e. xx_lastyear).
       CALL season &
            &          (kjpindex, dt_days, EndOfYear, &
            &           veget_cov, veget_cov_max, &
            &           humrel_daily, t2m_daily, tsoil_daily, soilhum_daily, &
            &           precip_daily, npp_daily, biomass, &
            &           turnover_daily, gpp_daily, when_growthinit, &
            &           maxhumrel_lastyear, maxhumrel_thisyear, &
            &           minhumrel_lastyear, minhumrel_thisyear, &
            &           maxgppweek_lastyear, maxgppweek_thisyear,maxnppweek_thisyear, &
            &           gdd0_lastyear, gdd0_thisyear, &
            &           precip_lastyear, precip_thisyear, &
            &           lm_lastyearmax, lm_thisyearmax, &
            &           maxfpc_lastyear, maxfpc_thisyear, &
            &           humrel_month, humrel_week, t2m_longterm, &
            &           tlong_ref, t2m_month, t2m_week, tsoil_month, soilhum_month, &
            &           npp_longterm, turnover_longterm, gpp_week,npp_week, &
            &           gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &           time_lowgpp, time_hum_min, hum_min_dormance, herbivores, & !)
!Chloe
            &           tsurf_daily, tsurf_year)
       
       !! 5.3 Use all processes included in stomate
       IF (control%ok_stomate) THEN

          !! 5.3.1  Activate stomate processes 
          ! Activate stomate processes (the complete list of processes depends 
          ! on whether the DGVM is used or not). Processes include: climate constraints 
          ! for PFTs, PFT dynamics, Phenology, Allocation, NPP (based on GPP and
          ! authothropic respiration), fire, mortality, vmax, assimilation temperatures,
          ! all turnover processes, light competition, sapling establishment, lai and 
          ! land cover change.
          CALL StomateLpj &
               &            (kjpindex, dt_days, EndOfYear, EndOfMonth, &
               &             neighbours, resolution, &
               &             clay, herbivores, &
               &             tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
               &             litterhum_daily, soilhum_daily, &
               &             maxhumrel_lastyear, minhumrel_lastyear, &
               &             gdd0_lastyear, precip_lastyear, &
               &             humrel_month, humrel_week, tlong_ref, t2m_month, t2m_week, &
               &             tsoil_month, soilhum_month, &
               &             gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
               &             turnover_longterm, gpp_daily, time_lowgpp, &
               &             time_hum_min, maxfpc_lastyear, resp_maint_part,&
               &             PFTpresent, age, fireindex, firelitter, &
               &             leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
               &             senescence, when_growthinit, litterpart, litter, &
               &             dead_leaves, carbon, black_carbon, lignin_struc, &
               &             veget_cov_max, npp_longterm, lm_lastyearmax, &
               &             veget_lastlight, everywhere, need_adjacent, RIP_time, &
               &             lai, rprof,npp_daily, turnover_daily, turnover_time,&
               &             control_moist_inst, control_temp_inst, soilcarbon_input_inst, &
               &             co2_to_bm_dgvm, co2_fire, &
               &             resp_hetero_d, resp_maint_d, resp_growth_d, &
               &             height, deadleaf_cover, vcmax, vjmax, &
               &             t_photo_min, t_photo_opt, t_photo_max,bm_to_litter,&
               &             prod10, prod100, flux10, flux100, veget_cov_max_new,&
               &             convflux, cflux_prod10, cflux_prod100, harvest_above, carb_mass_total, lcchange,&
               &             fpc_max, & !)
!Chloe
               &             ch4_flux_density_tot_peat,ch4_flux_density_dif_peat, ch4_flux_density_bub_peat, &
               &             ch4_flux_density_pla_peat,ch4_flux_density_goxid_peat,wpro_conso_carb,uold2_peat, tsurf_year, & !)
!Chloe
!Chloe Test provisoire :
               &  ch4_flux_density_tot_0, ch4_flux_density_dif_0, ch4_flux_density_bub_0, &
               &             ch4_flux_density_pla_0)

          !! 5.3.2 Calculate the total CO2 flux from land use change
          fco2_lu(:) = convflux(:) &
               &             + cflux_prod10(:)  &
               &             + cflux_prod100(:) &
               &             + harvest_above(:)

          !! 5.4 Calculate veget and veget_max
          veget_max(:,:) = zero 
          DO j = 1, nvm
             veget_max(:,j) = veget_max(:,j) + &
                  & veget_cov_max(:,j) * ( 1.-totfrac_nobio(:) )
          ENDDO
          
          !! 5.5 Photosynthesis parameters
          assim_param(:,:,ivcmax) = zero
          assim_param(:,:,ivjmax) = zero
          assim_param(:,:,itmin) = zero
          assim_param(:,:,itopt) = zero
          assim_param(:,:,itmax) = zero
          DO j = 2,nvm
             assim_param(:,j,ivcmax) = vcmax(:,j)
          ENDDO
          DO j = 2, nvm
             assim_param(:,j,ivjmax) = vjmax(:,j)
          ENDDO
          DO j = 2, nvm
             assim_param(:,j,itmin) = t_photo_min(:,j)
          ENDDO
          DO j = 2, nvm
             assim_param(:,j,itopt) = t_photo_opt(:,j)
          ENDDO
          DO j = 2, nvm
             assim_param(:,j,itmax) = t_photo_max(:,j)
          ENDDO

          !! 5.6 Update forcing variables for soil carbon
          IF (ldcarbon_write  .AND. TRIM(Cforcing_name) /= 'NONE') THEN
             npp_tot(:) = 0
             DO j=2,nvm
                npp_tot(:) = npp_tot(:) + npp_daily(:,j)
             ENDDO
             ! ::nbyear Number of years saved for carbon spinup
             sf_time = MODULO(REAL(date,r_std)-1,one_year*REAL(nbyear,r_std))
             iatt=FLOOR(sf_time/dt_forcesoil) + 1
             IF (iatt == 0) iatt = iatt_old + 1
             IF ((iatt<iatt_old) .and. (.not. cumul_Cforcing)) THEN
                nforce(:)=0
                soilcarbon_input(:,:,:,:) = zero
                control_moist(:,:,:) = zero
                control_temp(:,:,:) = zero
                npp_equil(:,:) = zero
             ENDIF
             iatt_old = iatt
             ! Update forcing
             nforce(iatt) = nforce(iatt)+1
             soilcarbon_input(:,:,:,iatt) = soilcarbon_input(:,:,:,iatt) + soilcarbon_input_daily(:,:,:)
             control_moist(:,:,iatt) = control_moist(:,:,iatt) + control_moist_daily(:,:)
             control_temp(:,:,iatt) = control_temp(:,:,iatt) + control_temp_daily(:,:)
             npp_equil(:,iatt) = npp_equil(:,iatt) + npp_tot(:)
          ENDIF

       ENDIF ! control%ok_stomate
       
       !! 5.8 Write forcing file if ::ldforcing_write=.TRUE.
       ! Note: if STOMATE is run in coupled mode the forcing file is written
       ! If run in stand-alone mode, the forcing file is read!
       IF ( ldforcing_write .AND. TRIM(forcing_name) /= 'NONE' ) THEN
          
          !! 5.8.1 Convert GPP to sechiba time steps
          ! GPP is multiplied by coverage to obtain forcing @tex $(gC m^{-2} dt_slow^{-1})$\f \end@tex $(m^2 m^{-2})$ @endtexonly
          ! @tex$ m^{-2}$ @endtex remains in the units because ::veget_cov_max is a fraction, not a 
          ! surface area. In sechiba values are ponderated by surface and frac_no_bio. 
          ! At the beginning of stomate, the units are converted. 
          ! When we use forcesoil we call sechiba_main and so we need the have the same units as in sechiba.
          gpp_daily_x(:,:) = zero
          DO j = 2, nvm             
             gpp_daily_x(:,j) = gpp_daily_x(:,j) + &
              & gpp_daily(:,j) * dt_slow / one_day * veget_cov_max(:,j)
          ENDDO
          
          ! Bare soil moisture availability has not been treated
          ! in STOMATE, update it here
          humrel_daily(:,ibare_sechiba) = humrel(:,ibare_sechiba)   

          ! Update index to store the next forcing step in memory
          iisf = iisf+1

          ! How many times have we treated this forcing state
          xn = REAL(nf_cumul(isf(iisf)),r_std)
          
          !! 5.8.2 Cumulate forcing variables
          ! Cumulate forcing variables (calculate average)
          ! Note: precipitation is multiplied by dt_slow/one_day to be consistent with 
          ! the units in sechiba
          IF (cumul_forcing) THEN
             clay_fm(:,iisf) = (xn*clay_fm(:,iisf)+clay(:))/(xn+1.)
             humrel_daily_fm(:,:,iisf) = &
                  & (xn*humrel_daily_fm(:,:,iisf) + humrel_daily(:,:))/(xn+1.)
             litterhum_daily_fm(:,iisf) = &
                  & (xn*litterhum_daily_fm(:,iisf)+litterhum_daily(:))/(xn+1.)
             t2m_daily_fm(:,iisf) = &
                  & (xn*t2m_daily_fm(:,iisf)+t2m_daily(:))/(xn+1.)
             t2m_min_daily_fm(:,iisf) = &
                  & (xn*t2m_min_daily_fm(:,iisf)+t2m_min_daily(:))/(xn+1.)
             tsurf_daily_fm(:,iisf) = &
                  & (xn*tsurf_daily_fm(:,iisf)+tsurf_daily(:))/(xn+1.)
             tsoil_daily_fm(:,:,iisf) = &
                  & (xn*tsoil_daily_fm(:,:,iisf)+tsoil_daily(:,:))/(xn+1.)
             soilhum_daily_fm(:,:,iisf) = &
                  & (xn*soilhum_daily_fm(:,:,iisf)+soilhum_daily(:,:))/(xn+1.)
             precip_fm(:,iisf) = &
                  & (xn*precip_fm(:,iisf)+precip_daily(:)*dt_slow/one_day)/(xn+1.)
             gpp_daily_fm(:,:,iisf) = &
                  & (xn*gpp_daily_fm(:,:,iisf) + gpp_daily_x(:,:))/(xn+1.)
             veget_fm(:,:,iisf) = &
                  & (xn*veget_fm(:,:,iisf) + veget(:,:) )/(xn+1.)
             veget_max_fm(:,:,iisf) = &
                  & (xn*veget_max_fm(:,:,iisf) + veget_max(:,:) )/(xn+1.)
             lai_fm(:,:,iisf) = &
                  & (xn*lai_fm(:,:,iisf) + lai(:,:) )/(xn+1.)
          ELSE
             ! Here we just calculate the values
             clay_fm(:,iisf) = clay(:)
             humrel_daily_fm(:,:,iisf) = humrel_daily(:,:)
             litterhum_daily_fm(:,iisf) = litterhum_daily(:)
             t2m_daily_fm(:,iisf) = t2m_daily(:)
             t2m_min_daily_fm(:,iisf) =t2m_min_daily(:)
             tsurf_daily_fm(:,iisf) = tsurf_daily(:)
             tsoil_daily_fm(:,:,iisf) =tsoil_daily(:,:)
             soilhum_daily_fm(:,:,iisf) =soilhum_daily(:,:)
             precip_fm(:,iisf) = precip_daily(:)
             gpp_daily_fm(:,:,iisf) =gpp_daily_x(:,:)
             veget_fm(:,:,iisf) = veget(:,:)
             veget_max_fm(:,:,iisf) =veget_max(:,:)
             lai_fm(:,:,iisf) =lai(:,:)
          ENDIF
          nf_cumul(isf(iisf)) = nf_cumul(isf(iisf))+1

          ! 5.8.3 Do we have to write the forcing states?
          IF (iisf == nsfm) THEN

             !! 5.8.3.1 Write these forcing states
             CALL forcing_write(forcing_id,1,nsfm)
             ! determine which forcing states must be read
             isf(1) = isf(nsfm)+1
             IF ( isf(1) > nsft ) isf(1) = 1
             DO iisf = 2, nsfm
                isf(iisf) = isf(iisf-1)+1
                IF (isf(iisf) > nsft)  isf(iisf) = 1
             ENDDO

             ! Read forcing variables - for debug use only
             ! CALL forcing_read(forcing_id,nsfm)
             iisf = 0

          ENDIF

       ENDIF


       !! 5.9 Compute daily CO2 flux (AR5 output - not essential)
       ! CO2 flux in @tex $gC m^{-2} s^{-1}$ @endtex (positive towards the atmosphere) is sum of:
       ! (1) heterotrophic respiration from ground + (2) maintenance respiration 
       ! from the plants + (3) growth respiration from the plants + (4) co2 
       ! emissions from fire - (5) co2 taken up in the DGVM to establish 
       ! saplings - (6) co2 taken up by photosyntyhesis
       co2_flux_daily(:,:)=   &
            & resp_maint_d(:,:) + resp_growth_d(:,:) + resp_hetero_d(:,:) + &
            & co2_fire(:,:) - co2_to_bm_dgvm(:,:) - gpp_daily(:,:)
       IF ( hist_id_stom_IPCC > 0 ) THEN
          vartmp(:) = SUM(co2_flux_daily*veget_cov_max,dim=2)/1e3/one_day*contfrac
          CALL histwrite (hist_id_stom_IPCC, "nep", itime, &
               vartmp, kjpindex, hori_index)
       ENDIF

       ! See 5.9 for details on NEP + fire. At the monthly time step also 
       ! harvest and land use change are calculated
       co2_flux_monthly(:,:) = co2_flux_monthly(:,:) + co2_flux_daily(:,:)
       harvest_above_monthly(:) = harvest_above_monthly(:) + harvest_above(:)
       cflux_prod_monthly(:) = cflux_prod_monthly(:) + convflux(:) + & 
        & cflux_prod10(:) + cflux_prod100(:)
      
       !! 5.10 Compute monthly CO2 fluxes 
       IF ( EndOfMonth ) THEN
          IF ( control%ok_stomate ) THEN
             
             !! 5.10.1 Write history file for monthly fluxes
             CALL histwrite (hist_id_stomate, 'CO2FLUX', itime, &
                  co2_flux_monthly, kjpindex*nvm, horipft_index)
          ENDIF
          !?? I (=VB) translated the French, but the whole stuff does not make sense to me.
          ! If one deletes the montly cumulation,
          ! one should not forget this change in resolution(:,1)*resolution(:,2)*contfrac(:)
          ! Si on supprimer le cumul par mois, 
          ! il ne faut pas oublier cette modif resolution(:,1)*resolution(:,2)*contfrac(:)
          ! Should be supressed, this is post-processing 
          DO j=2, nvm
             co2_flux_monthly(:,j) = co2_flux_monthly(:,j)* &
                  resolution(:,1)*resolution(:,2)*contfrac(:)
          ENDDO

          ! Should be supressed, this is post-processing
          ! ?? How does it differ from co2_flux_monthly??
          net_co2_flux_monthly = zero
          DO ji=1,kjpindex
             DO j=2,nvm
                net_co2_flux_monthly = net_co2_flux_monthly + &
                     &  co2_flux_monthly(ji,j)*veget_cov_max(ji,j)
             ENDDO
          ENDDO

     
          !! 5.10.2 Cumulative fluxes of land use cover change, harvest and net biosphere production
          ! Parallel processing, gather the information from different processors. first argument is the lo
          ! local variable, the second argument is the global variable. bcast send it to all processors.
          net_cflux_prod_monthly_sum = &
              &  SUM(cflux_prod_monthly(:)*resolution(:,1)*resolution(:,2)*contfrac(:))*1e-15
          CALL reduce_sum(net_cflux_prod_monthly_sum,net_cflux_prod_monthly_tot)
          CALL bcast(net_cflux_prod_monthly_tot)
          net_harvest_above_monthly_sum = &
             &   SUM(harvest_above_monthly(:)*resolution(:,1)*resolution(:,2)*contfrac(:))*1e-15
          CALL reduce_sum(net_harvest_above_monthly_sum,net_harvest_above_monthly_tot)
          CALL bcast(net_harvest_above_monthly_tot)
          net_co2_flux_monthly = net_co2_flux_monthly*1e-15
          CALL reduce_sum(net_co2_flux_monthly,net_co2_flux_monthly_sum)
          CALL bcast(net_co2_flux_monthly_sum)
          net_biosp_prod_monthly_tot =  &
             & ( net_co2_flux_monthly_sum + net_cflux_prod_monthly_tot + &
             & net_harvest_above_monthly_tot )
          
          WRITE(numout,9010) 'GLOBAL net_cflux_prod_monthly    (Peta gC/month)  = ',net_cflux_prod_monthly_tot
          WRITE(numout,9010) 'GLOBAL net_harvest_above_monthly (Peta gC/month)  = ',net_harvest_above_monthly_tot
          WRITE(numout,9010) 'GLOBAL net_co2_flux_monthly      (Peta gC/month)  = ',net_co2_flux_monthly_sum
          WRITE(numout,9010) 'GLOBAL net_biosp_prod_monthly    (Peta gC/month)  = ',net_biosp_prod_monthly_tot

9010  FORMAT(A52,F17.14)

!! DELETE
!!$          IF ( control%ok_stomate ) THEN
!!$             vartmp(:)=net_co2_flux_monthly_sum
!!$             CALL histwrite (hist_id_stomate, 'CO2FLUX_MONTHLY_SUM', itime, &
!!$                  vartmp, kjpindex, hori_index )
!!$          ENDIF
!!$          IF (is_root_prc) THEN
!!$             OPEN( unit=39,              &
!!$                  file="stomate_co2flux.data", &
!!$                  action="write",              &
!!$                  position="append",           &
!!$                  iostat=ios  )
!!$             IF ( ios /= 0 ) THEN
!!$                STOP "Erreur lors de la lecture/ecriture du fichier stomate_co2flux.data"
!!$             ELSE
!!$                WRITE(numout,*)
!!$                WRITE(numout,*) "Ecriture du fichier stomate_co2flux.data"
!!$                WRITE(numout,*)
!!$             END IF
!!$             WRITE(39,*) net_co2_flux_monthly_sum
!!$             CLOSE( unit=39 )
!!$          ENDIF
          
          ! Reset Monthly values
          co2_flux_monthly(:,:) = zero
          harvest_above_monthly(:) = zero
          cflux_prod_monthly(:)    = zero

       ENDIF ! Monthly processes - at the end of the month
       
       !! 5.11 Reset daily variables
       humrel_daily(:,:) = zero
       litterhum_daily(:) = zero
       t2m_daily(:) = zero
       t2m_min_daily(:) = large_value
       tsurf_daily(:) = zero
       tsoil_daily(:,:) = zero
       soilhum_daily(:,:) = zero
       precip_daily(:) = zero
       gpp_daily(:,:) = zero
       !Chloe:
       wtsoil_daily(:,:)=zero
       resp_maint_part(:,:,:)=zero
       resp_hetero_d=zero
       IF (bavard >= 3) THEN
          WRITE(numout,*) 'stomate_main: daily processes done'
       ENDIF

    ENDIF  ! Daily processes - at the end of the day
    
  !! 6. Outputs from Stomate

    ! co2_flux receives a value from STOMATE only if STOMATE is activated.
    ! Otherwise, the calling hydrological module must do this itself.
    IF ( control%ok_stomate ) THEN

       !! 6.1 Respiration and fluxes
       resp_maint(:,:) = resp_maint_radia(:,:)*veget_cov_max(:,:)
       resp_maint(:,ibare_sechiba) = zero
       resp_growth(:,:)= resp_growth_d(:,:)*veget_cov_max(:,:)*dtradia/one_day
       resp_hetero(:,:) = resp_hetero_radia(:,:)*veget_cov_max(:,:)

       !! 6.2 Derived CO2 fluxes
       ! CO2 flux in gC m^{-2} s^{-1} (positive towards the atmosphere) is sum of:
       ! (1) heterotrophic respiration from ground + (2) maintenance respiration 
       ! from the plants + (3) growth respiration from the plants + (4) co2 
       ! emissions from fire - (5) co2 taken up in the DGVM to establish 
       ! saplings - (6) co2 taken up by photosyntyhesis
       co2_flux(:,:) = resp_hetero(:,:) + resp_maint(:,:) + resp_growth(:,:) &
            & + (co2_fire(:,:)-co2_to_bm_dgvm(:,:))*veget_cov_max(:,:)/one_day &
            & - gpp(:,:)
    ENDIF
   
  !! 7. Messages
    
    IF ( (bavard >= 2).AND.EndOfYear.AND.do_slow) THEN
       WRITE(numout,*) 'stomate: EndOfYear'
    ENDIF
    IF (bavard >= 4) WRITE(numout,*) 'Leaving stomate_main'
    IF (long_print) WRITE (numout,*) ' stomate_main done '


  END SUBROUTINE stomate_main


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_init
!!
!>\BRIEF        The routine is called only at the first simulation. At that 
!! time settings and flags are read and checked for internal consistency and 
!! memory is allocated for the variables in stomate.
!!
!! DESCRIPTION  : The routine reads the 
!! following flags from the run definition file:
!! -bavard (level of online diagnostic in stomate)\n
!! -ipd (index of grid point for online diagnostics)\n
!! -ok_herbivores (flag to activate herbivores)\n
!! -treat_expansion (flag to activate PFT expansion across a pixel\n
!! -harvest_agri (flag to harvest aboveground biomass from agricultural PFTs)\n
!! \n
!! Check for inconsistent setting between the following flags:
!! -control%ok_stomate\n
!! -control%ok_dgvm\n
!! -control%ok_co2\n
!! -ldforcing_write\n
!! \n
!! Memory is allocated for all the variables of stomate and new indexing tables 
!! are build. New indexing tables are needed because a single pixel can conatin 
!! several PFTs. The new indexing tables have separate indices for the different 
!! PFTs. Similar index tables are build for land use cover change.\n
!! \n
!! Several global variables and land cover change variables are initialized to 
!! zero.\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory and builds new indexing 
!! variables for later use.\n 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_init &
       &  (kjpij, kjpindex, index, ldforcing_write, lalo, &
       &   rest_id_stom, hist_id_stom, hist_id_stom_IPCC)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
!Chloe
   !, nb of vertical layers for CH4 diffussion
!Chloe : et oh, les constantes doivent etre defini dans constantes.f90 !
!Chloe change le 371 en 171 :
    INTEGER(i_std),PARAMETER                         :: n = 171


    INTEGER(i_std),INTENT(in)                    :: kjpij             !! Total size of the un-compressed grid, including 
                                                                      !! oceans (unitless) 
    INTEGER(i_std),INTENT(in)                    :: kjpindex          !! Domain size - number of terrestrial pixels 
                                                                      !! (unitless) 
    INTEGER(i_std),INTENT(in)                    :: rest_id_stom      !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT(in)                    :: hist_id_stom      !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                    :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier 
    INTEGER(i_std),DIMENSION(kjpindex),INTENT(in):: index             !! Indices of the terrestrial pixels on the global 
                                                                      !! map 
    REAL(r_std),DIMENSION(kjpindex,2),INTENT(in) :: lalo              !! Geogr. coordinates (latitude,longitude) (degrees)
    LOGICAL,INTENT(in)                           :: ldforcing_write   !! Flag activates writing of _forcing_ file
   
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    LOGICAL                                      :: l_error           !! Check errors in netcdf call
    INTEGER(i_std)                               :: ier               !! Check errors in netcdf call
    INTEGER(i_std)                               :: ji,j,ipd,l        !! Indices

    REAL(r_std)                                  :: tmp_day(1)        !! ??[DISPENSABLE]
    REAL(r_std)                                  :: zcanop            !! ??[DISPENSABLE] soil depth taken for canopy
    REAL(r_std),DIMENSION(nbdl)                  :: zsoil             !! ??[DISPENSABLE] soil depths at diagnostic levels 
!_ ================================================================================================================================
    
  !! 1. Online diagnostics

    !Config Key   = BAVARD
    !Config Desc  = level of online diagnostics in STOMATE (0-4)
    !Config If    = OK_STOMATE
    !Config Def   = 1
    !Config Help  = With this variable, you can determine
    !Config         how much online information STOMATE
    !Config         gives during the run. 0 means
    !Config         virtually no info.
    !Config Units = [-]
    !
    ! By default, ::bavard is set to 1 in stomate_constants
    bavard = 1
    ! Get ::bavard from run definition file
    CALL getin_p('BAVARD', bavard)

    IF ( kjpindex > 0 ) THEN
       !Config  Key  = STOMATE_DIAGPT
       !Config  Desc = Index of grid point for online diagnostics
       !Config If    = OK_STOMATE
       !Config  Def  = 1
       !Config  Help = This is the index of the grid point which
       !               will be used for online diagnostics.
       !Config Units = [-]
       ! By default ::ipd is set to 1
       ipd = 1
       ! Get ::ipd from run definition file
       CALL getin_p('STOMATE_DIAGPT',ipd)
       ipd = MIN( ipd, kjpindex )
       WRITE(numout,*) 'Stomate: '
       WRITE(numout,*) '  Index of grid point for online diagnostics: ',ipd
       WRITE(numout,*) '  Lon, lat:',lalo(ipd,2),lalo(ipd,1)
       WRITE(numout,*) '  Index of this point on GCM grid: ',index(ipd)
       !
    ENDIF
    
  !! 2. Check consistency of flags

    IF ( ( .NOT. control%ok_stomate ) .AND. control%ok_dgvm ) THEN
       WRITE(numout,*) 'Cannot do dynamical vegetation without STOMATE.'
       WRITE(numout,*) 'Inconsistency between ::control%ok_stomate and ::control%ok_dgvm'
       WRITE(numout,*) 'Stop: fatal error'
       STOP
    ENDIF

    IF ((.NOT.control%ok_co2).AND.control%ok_stomate) THEN
       WRITE(numout,*) 'Cannot call STOMATE without GPP.'
       WRITE(numout,*) 'Inconsistency between ::control%ok_stomate and ::control%ok_co2'
       WRITE(numout,*) 'Stop: fatal error'
       STOP
    ENDIF

    IF ( ( .NOT. control%ok_co2 ) .AND. ldforcing_write ) THEN
       WRITE(numout,*)'Cannot write forcing file if photosynthesis is not activated'
       WRITE(numout,*)'Inconsistency between ::ldforcing_write and ::control%ok_co2'
       WRITE(numout,*) 'Stop: fatal error'
       STOP
    ENDIF
    
  !! 3. Communicate settings
    
    WRITE(numout,*) 'stomate first call - overview of the activated flags:'
    WRITE(numout,*) '  Photosynthesis: ', control%ok_co2
    WRITE(numout,*) '  STOMATE: ', control%ok_stomate
    WRITE(numout,*) '  LPJ: ', control%ok_dgvm
    
  !! 4. Allocate memory for STOMATE's variables

    l_error = .FALSE.

    ALLOCATE(veget_cov_max(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for veget_cov_max. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(veget_cov_max_new(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for veget_cov_max_new. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ind(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ind. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(adapted(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for adapted. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(regenerate(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for regenerate. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(humrel_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for humrel_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(litterhum_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litterhum_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_min_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_min_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsurf_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsurf_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsoil_daily(kjpindex,nbdl),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsoil_daily. We stop. We need kjpindex*nbdl words',kjpindex,nbdl
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(soilhum_daily(kjpindex,nbdl),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soilhum_daily. We stop. We need kjpindex*nbdl words',kjpindex,nbdl
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_daily. We stop. We need kjpindex words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
!Chloe++
    ALLOCATE(wtsoil_daily(kjpindex,nstm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for wtsoil_daily. We stop. We need kjpindex*nvm words',kjpindex,nstm
       STOP 'stomate_init'
    ENDIF
!Chloe--

    ALLOCATE(npp_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_daily(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_daily. We stop. We need kjpindex*nvm*nparts words', &
       &   kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_littercalc(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_littercalc. We stop. We need kjpindex*nvm*nparts words', & 
        &  kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(humrel_month(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for humrel_month. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(humrel_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for humrel_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_longterm(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tlong_ref(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tlong_ref. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_month(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_month. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_week(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_week. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
!Chloe
    ALLOCATE(tsurf_year(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turf_year. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
!Chloe

    ALLOCATE(tsoil_month(kjpindex,nbdl),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsoil_month. We stop. We need kjpindex*nbdl words',kjpindex,nbdl
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(soilhum_month(kjpindex,nbdl),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soilhum_month. We stop. We need kjpindex*nbdl words',kjpindex,nbdl
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(fireindex(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fireindex. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(firelitter(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for firelitter. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxhumrel_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxhumrel_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxhumrel_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxhumrel_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(minhumrel_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for minhumrel_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(minhumrel_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for minhumrel_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxgppweek_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxgppweek_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxgppweek_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxgppweek_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

!Chloe :
    ALLOCATE(maxnppweek_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxnppweek_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF


!!!Chloe++:
    ALLOCATE(wtold(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for wtold. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
!!Chloe--
    ALLOCATE(gdd0_lastyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd0_lastyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd0_thisyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd0_thisyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_lastyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_lastyear. We stop. We need kjpindex*nvm words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_thisyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_thisyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_m5_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_m5_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_midwinter(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_midwinter. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ncd_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ncd_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ngd_minus5(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ngd_minus5. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(PFTpresent(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for PFTpresent. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(npp_longterm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lm_lastyearmax(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lm_lastyearmax. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lm_thisyearmax(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lm_thisyearmax. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxfpc_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxfpc_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxfpc_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxfpc_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_longterm(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_longterm. We stop. We need kjpindex*nvm*nparts words', & 
       &    kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
!Chloe
    ALLOCATE(npp_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(biomass(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for biomass. We stop. We need kjpindex*nvm*nparts words',kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(senescence(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for senescence. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(when_growthinit(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for when_growthinit. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(age(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_radia(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_radia. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_growth_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_growth_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

!Chloe Test provisoire : 
    ALLOCATE(ch4_flux_density_tot_0(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)

    ALLOCATE(ch4_flux_density_dif_0(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(ch4_flux_density_bub_0(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(ch4_flux_density_pla_0(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)


!Chloe
    ALLOCATE(ch4_flux_density_tot_peat(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(ch4_flux_density_dif_peat(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(ch4_flux_density_bub_peat(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(ch4_flux_density_pla_peat(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
!Chloe add goxid methanotrophy
    ALLOCATE(ch4_flux_density_goxid_peat(kjpindex),stat=ier)
        l_error = l_error .OR. (ier /= 0)
    !ALLOCATE(wpro_conso_carb(kjpindex),stat=ier)
    !        l_error = l_error .OR. (ier /= 0)
! Chloe --

     ALLOCATE(uo_0(kjpindex,n),stat=ier)
     l_error = l_error .OR. (ier /= 0)

     ALLOCATE(uold2_0(kjpindex,n),stat=ier)
     l_error = l_error .OR. (ier /= 0)

    ALLOCATE(uo_peat(kjpindex,n),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(uold2_peat(kjpindex,n),stat=ier)
    l_error = l_error .OR. (ier /= 0)
!Chloe--

    ALLOCATE(co2_fire(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_fire. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(co2_to_bm_dgvm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_to_bm_dgvm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(veget_lastlight(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for veget_lastlight. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(everywhere(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for everywhere. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(need_adjacent(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for need_adjacent. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(leaf_age(kjpindex,nvm,nleafages),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_age. We stop. We need kjpindex*nvm*nleafages words', & 
       &      kjpindex,nvm,nleafages
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(leaf_frac(kjpindex,nvm,nleafages),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_frac. We stop. We need kjpindex*nvm*nleafages words', & 
       &      kjpindex,nvm,nleafages
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(RIP_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for RIP_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(time_lowgpp(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for time_lowgpp. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(time_hum_min(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for time_hum_min. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(hum_min_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for hum_min_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(litterpart(kjpindex,nvm,nlitt),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litterpart. We stop. We need kjpindex*nvm*nlitt words',  &
       &  kjpindex,nvm,nlitt
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(litter(kjpindex,nlitt,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litter. We stop. We need kjpindex*nlitt*nvm*nlevs words', & 
       &    kjpindex,nlitt,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(dead_leaves(kjpindex,nvm,nlitt),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for dead_leaves. We stop. We need kjpindex*nvm*nlitt words', & 
       &   kjpindex,nvm,nlitt
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(carbon(kjpindex,ncarb,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for carbon. We stop. We need kjpindex*ncarb*nvm words',kjpindex,ncarb,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(black_carbon(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for black_carbon. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lignin_struc(kjpindex,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lignin_struc. We stop. We need kjpindex*nvm*nlevs words',kjpindex,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(co2_flux_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_flux_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(co2_flux_monthly(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_flux_monthly. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod_monthly(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod_monthly. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
 
    ALLOCATE (harvest_above_monthly(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for harvest_above_monthly. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_litter(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_litter. We stop. We need kjpindex*nvm*nparts words', & 
       &    kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_littercalc(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_littercalc. We stop. We need kjpindex*nvm*nparts words', &
       &   kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(herbivores(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for herbivores. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(hori_index(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for hori_index. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(horipft_index(kjpindex*nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horipft_index. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_part_radia(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_part_radia. We stop. We need kjpindex*nvm*nparts words', &
       &  kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_radia(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_radia. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_part(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_part. We stop. We need kjpindex*nvm*nparts words', &
       &    kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF
    resp_maint_part(:,:,:) = zero

    ALLOCATE (horip10_index(kjpindex*10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip10_index. We stop. We need kjpindex*10 words',kjpindex,10
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (horip100_index(kjpindex*100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip100_index. We stop. We need kjpindex*100 words',kjpindex,100
       STOP 'stomate_init'
    ENDIF

!Chloe++ : 
    ALLOCATE (horip171_index(kjpindex*171), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip171_index. We stop. We need kjpindex*100 words',kjpindex,100
       STOP 'stomate_init'
    ENDIF
!Chloe--

    ALLOCATE (horip11_index(kjpindex*11), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip11_index. We stop. We need kjpindex*11 words',kjpindex,11
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (horip101_index(kjpindex*101), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip101_index. We stop. We need kjpindex*101 words',kjpindex,101
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (prod10(kjpindex,0:10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod10. We stop. We need kjpindex*11 words',kjpindex,11
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (prod100(kjpindex,0:100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod100. We stop. We need kjpindex*101 words',kjpindex,101
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (flux10(kjpindex,10), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux10. We stop. We need kjpindex*10 words',kjpindex,10
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (flux100(kjpindex,100), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux100. We stop. We need kjpindex*100 words',kjpindex,100
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (convflux(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for convflux. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod10(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod10. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (cflux_prod100(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for cflux_prod100. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (harvest_above(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for harvest_above. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (carb_mass_total(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for carb_mass_total. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (soilcarbon_input_daily(kjpindex,ncarb,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soilcarbon_input_daily. We stop. We need kjpindex*ncarb*nvm words', & 
       &    kjpindex,ncarb,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (control_temp_daily(kjpindex,nlevs), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for control_temp_daily. We stop. We need kjpindex*nlevs words',kjpindex,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (control_moist_daily(kjpindex,nlevs), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for control_moist_daily. We stop. We need kjpindex*nlevs words',kjpindex,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fpc_max(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fpc_max. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF
    
  !! 5. File definitions

    ! Store history and restart files in common variables
    hist_id_stomate = hist_id_stom
    hist_id_stomate_IPCC = hist_id_stom_IPCC
    rest_id_stomate = rest_id_stom
    
    ! In STOMATE reduced grids are used containing only terrestrial pixels.
    ! Build a new indexing table for the vegetation fields separating 
    ! between the different PFTs. Note that ::index has dimension (kjpindex) 
    ! wheras ::indexpft has dimension (kjpindex*nvm). 

    hori_index(:) = index(:)

    DO j = 1, nvm
       DO ji = 1, kjpindex
          horipft_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij
       ENDDO
    ENDDO

    ! Similar index tables are build for the land cover change variables
!Chloe build for nvert index 
    DO j=1,n
        DO ji=1,kjpindex
            horip171_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij
       ENDDO
    ENDDO

     DO j = 1, 10
       DO ji = 1, kjpindex
          horip10_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij
       ENDDO
    ENDDO

    DO j = 1, 100
       DO ji = 1, kjpindex
          horip100_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij
       ENDDO
    ENDDO

    DO j = 1, 11
       DO ji = 1, kjpindex
          horip11_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij
       ENDDO
    ENDDO

    DO j = 1, 101
       DO ji = 1, kjpindex
          horip101_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij
       ENDDO
    ENDDO

  !! 6. Initialization of global and land cover change variables. 

    ! All variables are cumulative variables. bm_to_litter is not and is therefore
    ! excluded
    !   bm_to_litter(:,:,:) = zero
    turnover_daily(:,:,:) = zero
    resp_hetero_d(:,:) = zero
    co2_flux_daily(:,:) = zero
    co2_flux_monthly(:,:) = zero
    cflux_prod_monthly(:) = zero
    harvest_above_monthly(:) = zero
    control_moist_daily(:,:) = zero
    control_temp_daily(:,:) = zero
    soilcarbon_input_daily(:,:,:) = zero
    ! Land cover change variables
    prod10(:,:)  = zero
    prod100(:,:) = zero
    flux10(:,:)  = zero
    flux100(:,:) = zero
    convflux(:)  = zero
    cflux_prod10(:) = zero
    cflux_prod100(:) = zero
    fpc_max(:,:)=zero
    
  END SUBROUTINE stomate_init


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_clear
!!
!>\BRIEF        Deallocate memory of the stomate variables.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_clear

  !! 1. Deallocate all dynamics variables

    IF (ALLOCATED(veget_cov_max)) DEALLOCATE(veget_cov_max)
    IF (ALLOCATED(veget_cov_max_new)) DEALLOCATE(veget_cov_max_new)
    IF (ALLOCATED(ind)) DEALLOCATE(ind)
    IF (ALLOCATED(adapted)) DEALLOCATE(adapted)
    IF (ALLOCATED(regenerate)) DEALLOCATE(regenerate)
    IF (ALLOCATED(humrel_daily)) DEALLOCATE(humrel_daily)
    IF (ALLOCATED(litterhum_daily)) DEALLOCATE(litterhum_daily)
    IF (ALLOCATED(t2m_daily))  DEALLOCATE(t2m_daily)
    IF (ALLOCATED(t2m_min_daily))  DEALLOCATE(t2m_min_daily)
    IF (ALLOCATED(tsurf_daily))  DEALLOCATE(tsurf_daily)
    IF (ALLOCATED(tsoil_daily)) DEALLOCATE(tsoil_daily)
    IF (ALLOCATED(soilhum_daily)) DEALLOCATE(soilhum_daily)
    IF (ALLOCATED(precip_daily)) DEALLOCATE(precip_daily)
    IF (ALLOCATED(gpp_daily)) DEALLOCATE(gpp_daily)
    IF (ALLOCATED(npp_daily)) DEALLOCATE(npp_daily)
    IF (ALLOCATED(turnover_daily)) DEALLOCATE(turnover_daily)
    IF (ALLOCATED(turnover_littercalc)) DEALLOCATE(turnover_littercalc)
    IF (ALLOCATED(humrel_month)) DEALLOCATE(humrel_month)
    IF (ALLOCATED(humrel_week)) DEALLOCATE(humrel_week)
    IF (ALLOCATED(t2m_longterm)) DEALLOCATE(t2m_longterm)
    IF (ALLOCATED(tlong_ref)) DEALLOCATE(tlong_ref)
    IF (ALLOCATED(t2m_month)) DEALLOCATE(t2m_month)
    IF (ALLOCATED(t2m_week)) DEALLOCATE(t2m_week)
    IF (ALLOCATED(tsoil_month)) DEALLOCATE(tsoil_month)
    IF (ALLOCATED(soilhum_month)) DEALLOCATE(soilhum_month)
    IF (ALLOCATED(fireindex)) DEALLOCATE(fireindex)
    IF (ALLOCATED(firelitter)) DEALLOCATE(firelitter)
    IF (ALLOCATED(maxhumrel_lastyear)) DEALLOCATE(maxhumrel_lastyear)
    IF (ALLOCATED(maxhumrel_thisyear)) DEALLOCATE(maxhumrel_thisyear)
    IF (ALLOCATED(minhumrel_lastyear)) DEALLOCATE(minhumrel_lastyear)
    IF (ALLOCATED(minhumrel_thisyear)) DEALLOCATE(minhumrel_thisyear)
    IF (ALLOCATED(maxgppweek_lastyear)) DEALLOCATE(maxgppweek_lastyear)
    IF (ALLOCATED(maxgppweek_thisyear)) DEALLOCATE(maxgppweek_thisyear)
!Chloe:
    IF (ALLOCATED(maxnppweek_thisyear)) DEALLOCATE(maxnppweek_thisyear)
    IF (ALLOCATED(wtold)) DEALLOCATE(wtold)
    IF (ALLOCATED(wtsoil_daily)) DEALLOCATE(wtsoil_daily)
!Chloe--
    IF (ALLOCATED(gdd0_lastyear)) DEALLOCATE(gdd0_lastyear)
    IF (ALLOCATED(gdd0_thisyear)) DEALLOCATE(gdd0_thisyear)
    IF (ALLOCATED(precip_lastyear)) DEALLOCATE(precip_lastyear)
    IF (ALLOCATED(precip_thisyear)) DEALLOCATE(precip_thisyear)
    IF (ALLOCATED(gdd_m5_dormance)) DEALLOCATE(gdd_m5_dormance)
    IF (ALLOCATED(gdd_midwinter)) DEALLOCATE(gdd_midwinter)
    IF (ALLOCATED(ncd_dormance)) DEALLOCATE(ncd_dormance)
    IF (ALLOCATED(ngd_minus5))  DEALLOCATE(ngd_minus5)
    IF (ALLOCATED(PFTpresent)) DEALLOCATE(PFTpresent)
    IF (ALLOCATED(npp_longterm)) DEALLOCATE(npp_longterm)
    IF (ALLOCATED(lm_lastyearmax)) DEALLOCATE(lm_lastyearmax)
    IF (ALLOCATED(lm_thisyearmax)) DEALLOCATE(lm_thisyearmax)
    IF (ALLOCATED(maxfpc_lastyear)) DEALLOCATE(maxfpc_lastyear)
    IF (ALLOCATED(maxfpc_thisyear)) DEALLOCATE(maxfpc_thisyear)
    IF (ALLOCATED(turnover_longterm)) DEALLOCATE(turnover_longterm)
    IF (ALLOCATED(gpp_week)) DEALLOCATE(gpp_week)
!Chloe :
    IF (ALLOCATED(npp_week)) DEALLOCATE(npp_week)
    IF (ALLOCATED(biomass)) DEALLOCATE(biomass)
    IF (ALLOCATED(senescence)) DEALLOCATE(senescence)
    IF (ALLOCATED(when_growthinit)) DEALLOCATE(when_growthinit)
    IF (ALLOCATED(age))  DEALLOCATE(age)
    IF (ALLOCATED(resp_hetero_d)) DEALLOCATE(resp_hetero_d)
    IF (ALLOCATED(resp_hetero_radia)) DEALLOCATE(resp_hetero_radia)
    IF (ALLOCATED(resp_maint_d)) DEALLOCATE(resp_maint_d)
    IF (ALLOCATED(resp_growth_d)) DEALLOCATE(resp_growth_d)
    IF (ALLOCATED(co2_fire)) DEALLOCATE(co2_fire)
    IF (ALLOCATED(co2_to_bm_dgvm)) DEALLOCATE(co2_to_bm_dgvm)
    IF (ALLOCATED(veget_lastlight)) DEALLOCATE(veget_lastlight)
    IF (ALLOCATED(everywhere)) DEALLOCATE(everywhere)
    IF (ALLOCATED(need_adjacent)) DEALLOCATE(need_adjacent)
    IF (ALLOCATED(leaf_age)) DEALLOCATE(leaf_age)
    IF (ALLOCATED(leaf_frac)) DEALLOCATE(leaf_frac)
    IF (ALLOCATED(RIP_time)) DEALLOCATE(RIP_time)
    IF (ALLOCATED(time_lowgpp)) DEALLOCATE(time_lowgpp)
    IF (ALLOCATED(time_hum_min)) DEALLOCATE(time_hum_min)
    IF (ALLOCATED(hum_min_dormance)) DEALLOCATE(hum_min_dormance)
    IF (ALLOCATED(litterpart)) DEALLOCATE(litterpart)
    IF (ALLOCATED(litter)) DEALLOCATE(litter)
    IF (ALLOCATED(dead_leaves)) DEALLOCATE(dead_leaves)
    IF (ALLOCATED(carbon)) DEALLOCATE(carbon)
    IF (ALLOCATED(black_carbon)) DEALLOCATE(black_carbon)
    IF (ALLOCATED(lignin_struc)) DEALLOCATE(lignin_struc)
    IF (ALLOCATED(turnover_time)) DEALLOCATE(turnover_time)
    IF (ALLOCATED(co2_flux_daily)) DEALLOCATE(co2_flux_daily)
    IF (ALLOCATED(co2_flux_monthly)) DEALLOCATE(co2_flux_monthly)
    IF (ALLOCATED(harvest_above_monthly)) DEALLOCATE (harvest_above_monthly)
    IF (ALLOCATED(cflux_prod_monthly)) DEALLOCATE (cflux_prod_monthly)
    IF (ALLOCATED(bm_to_litter)) DEALLOCATE(bm_to_litter)
    IF (ALLOCATED(bm_to_littercalc)) DEALLOCATE(bm_to_littercalc)
    IF (ALLOCATED(herbivores)) DEALLOCATE(herbivores)
    IF (ALLOCATED(resp_maint_part_radia)) DEALLOCATE(resp_maint_part_radia)
    IF (ALLOCATED(resp_maint_radia)) DEALLOCATE(resp_maint_radia)
    IF (ALLOCATED(resp_maint_part)) DEALLOCATE(resp_maint_part)


!Chloe test provisoire
    IF (ALLOCATED(uo_0)) DEALLOCATE(uo_0)
    IF (ALLOCATED(uold2_0)) DEALLOCATE(uold2_0)
   
    IF (ALLOCATED(ch4_flux_density_tot_0)) DEALLOCATE(ch4_flux_density_tot_0)
    IF (ALLOCATED(ch4_flux_density_dif_0)) DEALLOCATE(ch4_flux_density_dif_0)
    IF (ALLOCATED(ch4_flux_density_bub_0)) DEALLOCATE(ch4_flux_density_bub_0)
    IF (ALLOCATED(ch4_flux_density_pla_0)) DEALLOCATE(ch4_flux_density_pla_0)
!Chloe fin du test provisoire

!Chloe
    IF (ALLOCATED(uo_peat)) DEALLOCATE(uo_peat)
    IF (ALLOCATED(uold2_peat)) DEALLOCATE(uold2_peat)

    IF (ALLOCATED(ch4_flux_density_tot_peat)) DEALLOCATE(ch4_flux_density_tot_peat)
    IF (ALLOCATED(ch4_flux_density_dif_peat)) DEALLOCATE(ch4_flux_density_dif_peat)
    IF (ALLOCATED(ch4_flux_density_bub_peat)) DEALLOCATE(ch4_flux_density_bub_peat)
    IF (ALLOCATED(ch4_flux_density_pla_peat)) DEALLOCATE(ch4_flux_density_pla_peat)
    IF (ALLOCATED(ch4_flux_density_goxid_peat)) DEALLOCATE(ch4_flux_density_goxid_peat)
   ! IF (ALLOCATED(wpro_conso_carb)) DEALLOCATE(wpro_conso_carb)
!Chloe--  

    IF (ALLOCATED(hori_index)) DEALLOCATE(hori_index)
    IF (ALLOCATED(horipft_index)) DEALLOCATE(horipft_index)
    IF (ALLOCATED(clay_fm)) DEALLOCATE(clay_fm)
    IF (ALLOCATED(humrel_daily_fm)) DEALLOCATE(humrel_daily_fm)
    IF (ALLOCATED(litterhum_daily_fm))  DEALLOCATE(litterhum_daily_fm)
    IF (ALLOCATED(t2m_daily_fm))  DEALLOCATE(t2m_daily_fm)
    IF (ALLOCATED(t2m_min_daily_fm))  DEALLOCATE(t2m_min_daily_fm)
    IF (ALLOCATED(tsurf_daily_fm)) DEALLOCATE(tsurf_daily_fm)
    IF (ALLOCATED(tsoil_daily_fm)) DEALLOCATE(tsoil_daily_fm)
    IF (ALLOCATED(soilhum_daily_fm))  DEALLOCATE(soilhum_daily_fm)
    IF (ALLOCATED(precip_fm)) DEALLOCATE(precip_fm)
    IF (ALLOCATED(gpp_daily_fm))  DEALLOCATE(gpp_daily_fm)
    IF (ALLOCATED(veget_fm)) DEALLOCATE(veget_fm)
    IF (ALLOCATED(veget_max_fm)) DEALLOCATE(veget_max_fm)
    IF (ALLOCATED(lai_fm))  DEALLOCATE(lai_fm)

    IF (is_root_prc) THEN
       IF (ALLOCATED(clay_fm_g)) DEALLOCATE(clay_fm_g)
       IF (ALLOCATED(humrel_daily_fm_g)) DEALLOCATE(humrel_daily_fm_g)
       IF (ALLOCATED(litterhum_daily_fm_g))  DEALLOCATE(litterhum_daily_fm_g)
       IF (ALLOCATED(t2m_daily_fm_g))  DEALLOCATE(t2m_daily_fm_g)
       IF (ALLOCATED(t2m_min_daily_fm_g))  DEALLOCATE(t2m_min_daily_fm_g)
       IF (ALLOCATED(tsurf_daily_fm_g)) DEALLOCATE(tsurf_daily_fm_g)
       IF (ALLOCATED(tsoil_daily_fm_g)) DEALLOCATE(tsoil_daily_fm_g)
       IF (ALLOCATED(soilhum_daily_fm_g))  DEALLOCATE(soilhum_daily_fm_g)
       IF (ALLOCATED(precip_fm_g)) DEALLOCATE(precip_fm_g)
       IF (ALLOCATED(gpp_daily_fm_g))  DEALLOCATE(gpp_daily_fm_g)
       IF (ALLOCATED(veget_fm_g)) DEALLOCATE(veget_fm_g)
       IF (ALLOCATED(veget_max_fm_g)) DEALLOCATE(veget_max_fm_g)
       IF (ALLOCATED(lai_fm_g))  DEALLOCATE(lai_fm_g)
    ENDIF

    IF (ALLOCATED(isf)) DEALLOCATE(isf)
    IF (ALLOCATED(nf_written)) DEALLOCATE(nf_written)
    IF (ALLOCATED(nf_cumul)) DEALLOCATE(nf_cumul)
    IF (ALLOCATED(nforce)) DEALLOCATE(nforce)
    IF (ALLOCATED(control_moist)) DEALLOCATE(control_moist)
    IF (ALLOCATED(control_temp)) DEALLOCATE(control_temp)
    IF (ALLOCATED(soilcarbon_input)) DEALLOCATE(soilcarbon_input)
    IF ( ALLOCATED (horip10_index)) DEALLOCATE (horip10_index)
    IF ( ALLOCATED (horip100_index)) DEALLOCATE (horip100_index)
!Chloe:
    IF ( ALLOCATED (horip171_index)) DEALLOCATE (horip171_index)
    IF ( ALLOCATED (horip11_index)) DEALLOCATE (horip11_index)
    IF ( ALLOCATED (horip101_index)) DEALLOCATE (horip101_index)
    IF ( ALLOCATED (prod10)) DEALLOCATE (prod10)
    IF ( ALLOCATED (prod100)) DEALLOCATE (prod100)
    IF ( ALLOCATED (flux10)) DEALLOCATE (flux10)
    IF ( ALLOCATED (flux100)) DEALLOCATE (flux100)
    IF ( ALLOCATED (convflux)) DEALLOCATE (convflux)
    IF ( ALLOCATED (cflux_prod10)) DEALLOCATE (cflux_prod10)
    IF ( ALLOCATED (cflux_prod100)) DEALLOCATE (cflux_prod100)
    IF ( ALLOCATED (harvest_above)) DEALLOCATE (harvest_above)
    IF ( ALLOCATED (soilcarbon_input_daily)) DEALLOCATE (soilcarbon_input_daily)
    IF ( ALLOCATED (control_temp_daily)) DEALLOCATE (control_temp_daily)
    IF ( ALLOCATED (control_moist_daily)) DEALLOCATE (control_moist_daily)

    IF ( ALLOCATED (fpc_max)) DEALLOCATE (fpc_max)

 !! 2. reset l_first

    l_first_stomate=.TRUE.

 !! 3. call to clear functions

    CALL get_reftemp_clear
    CALL season_clear
    CALL stomatelpj_clear
    CALL littercalc_clear
    CALL vmax_clear
 
  END SUBROUTINE stomate_clear


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_var_init
!!
!>\BRIEF        Initialize variables of stomate with a none-zero initial value.
!! Subroutine is called only if ::ok_stomate = .TRUE. STOMATE diagnoses some 
!! variables for SECHIBA : assim_param, deadleaf_cover, etc. These variables can 
!! be recalculated from STOMATE's prognostic variables. Note that height is
!! saved in SECHIBA.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): leaf age (::leaf_age) and fraction of leaves in leaf 
!! age class (::leaf_frac). The maximum water on vegetation available for 
!! interception, fraction of soil covered by dead leaves
!! (::deadleaf_cover) and assimilation parameters (:: assim_param).
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_var_init &
       &  (kjpindex, veget_cov_max, leaf_age, leaf_frac, &
       &   tlong_ref, t2m_month, dead_leaves, &
       &   veget, lai, deadleaf_cover, assim_param)


  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                             :: kjpindex        !! Domain size - terrestrial pixels only
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)        :: veget           !! Fraction of pixel covered by PFT. Fraction 
                                                                             !! accounts for none-biological land covers 
                                                                             !! (unitless) 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)        :: veget_cov_max   !! Fractional coverage: maximum share of the pixel 
                                                                             !! covered by a PFT (unitless) 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)            :: tlong_ref       !! "long term" 2 meter reference temperatures (K)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)            :: t2m_month       !! Air temperature at 2 meter integrated over a 
                                                                             !! month (K) 
    REAL(r_std),DIMENSION(kjpindex,nvm,nlitt),INTENT(in)  :: dead_leaves     !! Metabolic and structural fraction of dead leaves 
                                                                             !! per ground area 
                                                                             !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)        :: lai             !! Leaf area index 
                                                                             !! @tex $(m^2 m{-2})$ @endtex 
    
    !! 0.2 Modified variables

    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages),INTENT(inout) :: leaf_age  !! Age of different leaf classes per PFT (days)
    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages),INTENT(inout) :: leaf_frac !! Fraction of leaves in leaf age class per PFT 
                                                                             !! (unitless; 1) 
    
    !! 0.3 Output variables

    REAL(r_std),DIMENSION(kjpindex), INTENT (out)         :: deadleaf_cover  !! Fraction of soil covered by dead leaves 
                                                                             !! (unitless) 
    REAL(r_std),DIMENSION(kjpindex,nvm,npco2),INTENT(out) :: assim_param    !! min+max+opt temperatures (K) & vmax for 
                                                                            !! photosynthesis  
                                                                            !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
    ! 0.4 Local variables
   
    REAL(r_std),PARAMETER                                 :: dt_0 = zero     !! Dummy time step, must be zero
    REAL(r_std),DIMENSION(kjpindex,nvm)                   :: vcmax           !! Dummy vcmax 
                                                                             !! @tex $(\mu mol m^{-2} s^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)                   :: vjmax           !! Dummy vjmax 
                                                                             !! @tex $(\mu mol m^{-2} s^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm)                   :: t_photo_min     !! Min temperature for photosynthesis (K)
    REAL(r_std),DIMENSION(kjpindex,nvm)                   :: t_photo_opt     !! Opt temperature for photosynthesis (K)
    REAL(r_std),DIMENSION(kjpindex,nvm)                   :: t_photo_max     !! Max temperature for photosynthesis (K)  
    INTEGER(i_std)                                        :: j               !! Index (untiless)
    
!_ ================================================================================================================================   

    ! Only if stomate is activated
    IF (control%ok_stomate) THEN
     
  !! 1. photosynthesis parameters

       !! 1.1 vcmax (stomate_vmax.f90)
       CALL vmax (kjpindex, dt_0, leaf_age, leaf_frac, vcmax, vjmax)
       
       !! 1.2 assimilation temperatures (stomate_assimtemp.f90)
       CALL assim_temp(kjpindex, tlong_ref, t2m_month, &
            &                    t_photo_min, t_photo_opt, t_photo_max)
       
       !! 1.3 transform into nvm vegetation types
       assim_param(:,:,ivcmax) = zero
       DO j = 2, nvm
          assim_param(:,j,ivcmax)=vcmax(:,j)
       ENDDO
       assim_param(:,:,ivjmax) = zero
       DO j = 2, nvm
          assim_param(:,j,ivjmax)=vjmax(:,j)
       ENDDO
       assim_param(:,:,itmin) = zero
       DO j = 2, nvm
          assim_param(:,j,itmin)=t_photo_min(:,j)
       ENDDO
       assim_param(:,:,itopt) = zero
       DO j = 2, nvm
          assim_param(:,j,itopt)=t_photo_opt(:,j)
       ENDDO
       assim_param(:,:,itmax) = zero
       DO j = 2, nvm
          assim_param(:,j,itmax)=t_photo_max(:,j)
       ENDDO
       
  !! 2. Dead leaf cover (stomate_litter.f90)
       
       CALL deadleaf (kjpindex, veget_cov_max, dead_leaves, deadleaf_cover)     
       
    ENDIF ! ok_stomate = .TRUE.
    
  END SUBROUTINE stomate_var_init


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_accu
!!
!>\BRIEF        Accumulate a variable for the time period specified by 
!! ::dt_tot or calculate the mean value over ::dt_tot.
!! 
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): accumulated or mean variable ::field_out:: 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_accu &
       &  (npts, n_dim2, dt_tot, dt, ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std),INTENT(in)                        :: npts      !! Domain size (unitless)
    INTEGER(i_std),INTENT(in)                        :: n_dim2    !! 2nd dimension (1 or nvm)
    REAL(r_std),INTENT(in)                           :: dt_tot    !! Current time step (days)
    REAL(r_std),INTENT(in)                           :: dt        !! Accumulated time step (days)
    LOGICAL,INTENT(in)                               :: ldmean    !! Flag to calculate the mean over
    REAL(r_std),DIMENSION(npts,n_dim2),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Output variables
    
    !! 0.3 Modified variables

    REAL(r_std),DIMENSION(npts,n_dim2),INTENT(inout) :: field_out !! Accumulated or mean field

    !! 0.4 Local variables

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:) = field_out(:,:)+field_in(:,:)*dt
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:) = field_out(:,:)/dt_tot
    ENDIF

  END SUBROUTINE stomate_accu


!! ================================================================================================================================
!! SUBROUTINE 	: init_forcing
!!
!>\BRIEF        Allocate memory for the variables containing the forcing data.
!! The maximum size of the allocated memory is specified in run definition file
!! (::max_totsize) and needs to be a compromise between charging the memory and 
!! accessing disks to get the forcing data.
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory for later use. 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE init_forcing (kjpindex,nsfm,nsft)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    
    INTEGER(i_std),INTENT(in) :: kjpindex !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std),INTENT(in) :: nsfm     !! Number of time steps that can be stored in memory (unitless)
    INTEGER(i_std),INTENT(in) :: nsft     !! Number of time steps in a year (unitless)

   !! 0.2 Output variables

   !! 0.3 Modified variables

   !! 0.4 Local variables

    LOGICAL                   :: l_error  !! Check errors in netcdf call
    INTEGER(i_std)            :: ier      !! Check errors in netcdf call
!_ ================================================================================================================================
    
  !! 1. Allocate memory

    ! Note ::nvm is number of PFTs and ::nbdl is number of soil layers
    l_error = .FALSE.
    ALLOCATE(clay_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables clay_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(humrel_daily_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables humrel_daily_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(litterhum_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables litterhum_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(t2m_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(t2m_min_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_min_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(tsurf_daily_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables tsurf_daily_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(tsoil_daily_fm(kjpindex,nbdl,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables tsoil_daily_fm ',kjpindex,nbdl,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(soilhum_daily_fm(kjpindex,nbdl,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables soilhum_daily_fm ',kjpindex,nbdl,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(precip_fm(kjpindex,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables precip_fm ',kjpindex,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(gpp_daily_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables gpp_daily_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(veget_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(veget_max_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_max_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(lai_fm(kjpindex,nvm,nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables lai_fm ',kjpindex,nvm,nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(isf(nsfm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables isf ',nsfm
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(nf_written(nsft),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables nf_written ',nsft
       STOP 'init_forcing'
    ENDIF
    ALLOCATE(nf_cumul(nsft),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables nf_cumul ',nsft
       STOP 'init_forcing'
    ENDIF
    
  !! 2. Allocate memory for the root processor only (parallel computing)

    ! Where, ::nbp_glo is the number of global continental points
    IF (is_root_prc) THEN
       ALLOCATE(clay_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables clay_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(humrel_daily_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables humrel_daily_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(litterhum_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables litterhum_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(t2m_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(t2m_min_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables t2m_min_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(tsurf_daily_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables tsurf_daily_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(tsoil_daily_fm_g(nbp_glo,nbdl,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables tsoil_daily_fm_g ',nbp_glo,nbdl,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(soilhum_daily_fm_g(nbp_glo,nbdl,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables soilhum_daily_fm_g ',nbp_glo,nbdl,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(precip_fm_g(nbp_glo,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables precip_fm_g ',nbp_glo,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(gpp_daily_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables gpp_daily_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(veget_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(veget_max_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables veget_max_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
       ALLOCATE(lai_fm_g(nbp_glo,nvm,nsfm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) 'Problem with memory allocation: forcing variables lai_fm_g ',nbp_glo,nvm,nsfm
          STOP 'init_forcing'
       ENDIF
    ELSE
       ! Allocate memory for co-processors
       ALLOCATE(clay_fm_g(0,nsfm),stat=ier)
       ALLOCATE(humrel_daily_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(litterhum_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(t2m_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(t2m_min_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(tsurf_daily_fm_g(0,nsfm),stat=ier)
       ALLOCATE(tsoil_daily_fm_g(0,nbdl,nsfm),stat=ier)
       ALLOCATE(soilhum_daily_fm_g(0,nbdl,nsfm),stat=ier)
       ALLOCATE(precip_fm_g(0,nsfm),stat=ier)
       ALLOCATE(gpp_daily_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(veget_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(veget_max_fm_g(0,nvm,nsfm),stat=ier)
       ALLOCATE(lai_fm_g(0,nvm,nsfm),stat=ier)
    ENDIF ! is_root_proc
    
    IF (l_error) THEN
       WRITE(numout,*) 'Problem with memory allocation: forcing variables'
       STOP 'init_forcing'
    ENDIF

  !! 3. Initilaize variables

    CALL forcing_zero
    
  END SUBROUTINE init_forcing


!! ================================================================================================================================
!! SUBROUTINE 	: forcing_zero
!!
!>\BRIEF        Initialize variables containing the forcing data; variables are 
!! set to zero.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE forcing_zero
    
    clay_fm(:,:) = zero
    humrel_daily_fm(:,:,:) = zero
    litterhum_daily_fm(:,:) = zero
    t2m_daily_fm(:,:) = zero
    t2m_min_daily_fm(:,:) = zero
    tsurf_daily_fm(:,:) = zero
    tsoil_daily_fm(:,:,:) = zero
    soilhum_daily_fm(:,:,:) = zero
    precip_fm(:,:) = zero
    gpp_daily_fm(:,:,:) = zero
    veget_fm(:,:,:) = zero
    veget_max_fm(:,:,:) = zero
    lai_fm(:,:,:) = zero
    
  END SUBROUTINE forcing_zero


!! ================================================================================================================================
!! SUBROUTINE 	: forcing_write
!!
!>\BRIEF        Appends data values to a netCDF file containing the forcing 
!! variables of the general processes in stomate.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): netCDF file
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE forcing_write(forcing_id,ibeg,iend)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)      :: forcing_id  !! File identifer of forcing file, assigned when netcdf is created
    INTEGER(i_std),INTENT(in)      :: ibeg, iend  !! First and last time step to be written

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                 :: iisf        !! Index of isf where isf is the number of time steps that can be 
                                                  !! stored in memory 
    INTEGER(i_std)                 :: iblocks     !! Index of block that is written
    INTEGER(i_std)                 :: nblocks     !! Number of blocks that needs to be written
    INTEGER(i_std)                 :: ier         !! Check errors in netcdf call
    INTEGER(i_std),DIMENSION(0:2)  :: ifirst      !! First block in memory - changes with iblocks
    INTEGER(i_std),DIMENSION(0:2)  :: ilast       !! Last block in memory - changes with iblocks
    INTEGER(i_std),PARAMETER       :: ndm = 10    !! Maximum number of dimensions
    INTEGER(i_std),DIMENSION(ndm)  :: start       !! First block to write
    INTEGER(i_std)                 :: ndim        !! Dimensions of forcing to be added to the netCDF
    INTEGER(i_std),DIMENSION(ndm)  :: count_force !! Number of elements in each dimension  
    INTEGER(i_std)                 :: vid         !! Variable identifer of netCDF
!_ ================================================================================================================================
    
  !! 1. Determine number of blocks of forcing variables that are stored in memory

    nblocks = 0
    ifirst(:) = 1
    ilast(:) = 1
    DO iisf = ibeg, iend
       IF (     (nblocks /= 0) &
            &      .AND.(isf(iisf) == isf(ilast(nblocks))+1)) THEN
          ! Last block found
          ilast(nblocks) = iisf
       ELSE
          ! First block found
          nblocks = nblocks+1
          IF (nblocks > 2)  STOP 'Problem in forcing_write'
          ifirst(nblocks) = iisf
          ilast(nblocks) = iisf
       ENDIF
    ENDDO

  !! 2. Gather distributed variables (parallel computing)

    CALL gather(clay_fm,clay_fm_g)
    CALL gather(humrel_daily_fm,humrel_daily_fm_g)
    CALL gather(litterhum_daily_fm,litterhum_daily_fm_g)
    CALL gather(t2m_daily_fm,t2m_daily_fm_g)
    CALL gather(t2m_min_daily_fm,t2m_min_daily_fm_g)
    CALL gather(tsurf_daily_fm,tsurf_daily_fm_g)
    CALL gather(tsoil_daily_fm,tsoil_daily_fm_g)
    CALL gather(soilhum_daily_fm,soilhum_daily_fm_g)
    CALL gather(precip_fm,precip_fm_g)
    CALL gather(gpp_daily_fm,gpp_daily_fm_g)
    CALL gather(veget_fm,veget_fm_g)
    CALL gather(veget_max_fm,veget_max_fm_g)
    CALL gather(lai_fm,lai_fm_g)
 
 !! 3. Append data to netCDF file
   
    IF (is_root_prc) THEN
       ! The netCDF file has been created earlier in this module, a file ID is available 
       ! and variables and dimensions have already been defined
       DO iblocks = 1, nblocks
          IF (ifirst(iblocks) /= ilast(iblocks)) THEN
             ndim = 2
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(clay_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'clay',vid)
             ier = NF90_PUT_VAR (forcing_id,vid, &
                  &              clay_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(humrel_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'humrel',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            humrel_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(litterhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'litterhum',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            litterhum_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            t2m_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_min_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m_min',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            t2m_min_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsurf_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsurf',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            tsurf_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsoil_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsoil',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            tsoil_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(soilhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'soilhum',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            soilhum_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(precip_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'precip',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            precip_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  & start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(gpp_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'gpp',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            gpp_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            veget_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_max_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget_max',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            veget_max_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(lai_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'lai',vid)
             ier = NF90_PUT_VAR (forcing_id, vid, &
                  &            lai_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
          ENDIF
       ENDDO
    ENDIF
    
  !! 4. Adjust flag of forcing file
    nf_written(isf(:)) = .TRUE.

  END SUBROUTINE forcing_write

  
!! ================================================================================================================================
!! SUBROUTINE 	: forcing_read
!!
!>\BRIEF        Read forcing file.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE forcing_read(forcing_id,nsfm)
   
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)  :: forcing_id           !! File identifer of forcing file, assigned when netcdf is created
    INTEGER(i_std),INTENT(in)  :: nsfm                 !! Number of time steps stored in memory        
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                 :: iisf              !! Index of isf where isf is the number of time steps that can be stored in 
                                                        !! memory 
    INTEGER(i_std)                 :: iblocks           !! Index of block that is written
    INTEGER(i_std)                 :: nblocks           !! Number of blocks that needs to be written
    INTEGER(i_std)                 :: ier               !! Check error of netcdf call
    INTEGER(i_std),DIMENSION(0:2)  :: ifirst            !! First block in memory - changes with iblocks
    INTEGER(i_std),DIMENSION(0:2)  :: ilast             !! Last block in memory - changes with iblocks
    INTEGER(i_std),PARAMETER       :: ndm = 10          !! Maximum number of dimensions
    INTEGER(i_std),DIMENSION(ndm)  :: start             !! First block to write
    INTEGER(i_std)                 :: ndim              !! Dimensions of forcing to be added to the netCDF
    INTEGER(i_std),DIMENSION(ndm)  :: count_force       !! Number of elements in each dimension
    INTEGER(i_std)                 :: vid               !! Variable identifer of netCDF
    LOGICAL, PARAMETER             :: check=.FALSE.     !! Flag for debugging 
    LOGICAL                        :: a_er=.FALSE.      !! Error catching from netcdf file
!_ ================================================================================================================================

    IF (check) WRITE(numout,*) "forcing_read "
    
  !! 1. Set to zero if the corresponding forcing state

    ! has not yet been written into the file  
    DO iisf = 1, nsfm
       IF (.NOT.nf_written(isf(iisf))) THEN
          clay_fm(:,iisf) = zero
          humrel_daily_fm(:,:,iisf) = zero
          litterhum_daily_fm(:,iisf) = zero
          t2m_daily_fm(:,iisf) = zero
          t2m_min_daily_fm(:,iisf) = zero
          tsurf_daily_fm(:,iisf) = zero
          tsoil_daily_fm(:,:,iisf) = zero
          soilhum_daily_fm(:,:,iisf) = zero
          precip_fm(:,iisf) = zero
          gpp_daily_fm(:,:,iisf) = zero
          veget_fm(:,:,iisf) = zero
          veget_max_fm(:,:,iisf) = zero
          lai_fm(:,:,iisf) = zero
       ENDIF
    ENDDO
    
  !! 2. determine blocks of forcing states that are contiguous in memory

    nblocks = 0
    ifirst(:) = 1
    ilast(:) = 1
    
    DO iisf = 1, nsfm
       IF (nf_written(isf(iisf))) THEN
          IF (     (nblocks /= 0) &
               &        .AND.(isf(iisf) == isf(ilast(nblocks))+1)) THEN

             ! element is contiguous with last element found
             ilast(nblocks) = iisf
          ELSE

             ! found first element of new block
             nblocks = nblocks+1
             IF (nblocks > 2)  STOP 'Problem in forcing_read'
             
             ifirst(nblocks) = iisf
             ilast(nblocks) = iisf
          ENDIF
       ENDIF
    ENDDO
    IF (check) WRITE(numout,*) "forcing_read nblocks, ifirst, ilast",nblocks, ifirst, ilast
    
  !! 3. Read variable values

    IF (is_root_prc) THEN
       DO iblocks = 1, nblocks
          IF (check) WRITE(numout,*) "forcing_read iblocks, ifirst(iblocks), ilast(iblocks)",iblocks, &
               ifirst(iblocks), ilast(iblocks)
          IF (ifirst(iblocks) /= ilast(iblocks)) THEN
             a_er=.FALSE.
             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(clay_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'clay',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            clay_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(humrel_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'humrel',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            humrel_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(litterhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'litterhum',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              litterhum_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              t2m_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(t2m_min_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'t2m_min',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              t2m_min_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsurf_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsurf',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              tsurf_daily_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(tsoil_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'tsoil',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              tsoil_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(soilhum_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'soilhum',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              soilhum_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 2;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(precip_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'precip',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &              precip_fm_g(:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(gpp_daily_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'gpp',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            gpp_daily_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            veget_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(veget_max_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'veget_max',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            veget_max_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)

             ndim = 3;
             start(1:ndim) = 1; start(ndim) = isf(ifirst(iblocks));
             count_force(1:ndim) = SHAPE(lai_fm_g)
             count_force(ndim) = isf(ilast(iblocks))-isf(ifirst(iblocks))+1
             ier = NF90_INQ_VARID (forcing_id,'lai',vid)
             a_er = a_er.OR.(ier /= 0)
             ier = NF90_GET_VAR (forcing_id, vid, &
                  &            lai_fm_g(:,:,ifirst(iblocks):ilast(iblocks)), &
                  &            start=start(1:ndim), count=count_force(1:ndim))
             a_er = a_er.OR.(ier /= 0)
             IF (a_er) THEN
                CALL ipslerr (3,'forcing_read', &
                     &        'PROBLEM when read forcing file', &
                     &        '','')
             ENDIF

          ENDIF ! (ifirst(iblocks) /= ilast(iblocks))
       ENDDO ! iblocks
    ENDIF ! is_root_prc

  !! 4. Distribute the variable over several processors

    CALL scatter(clay_fm_g,clay_fm)
    CALL scatter(humrel_daily_fm_g,humrel_daily_fm)
    CALL scatter(litterhum_daily_fm_g,litterhum_daily_fm)
    CALL scatter(t2m_daily_fm_g,t2m_daily_fm)
    CALL scatter(t2m_min_daily_fm_g,t2m_min_daily_fm)
    CALL scatter(tsurf_daily_fm_g,tsurf_daily_fm)
    CALL scatter(tsoil_daily_fm_g,tsoil_daily_fm)
    CALL scatter(soilhum_daily_fm_g,soilhum_daily_fm)
    CALL scatter(precip_fm_g,precip_fm)
    CALL scatter(gpp_daily_fm_g,gpp_daily_fm)
    CALL scatter(veget_fm_g,veget_fm)
    CALL scatter(veget_max_fm_g,veget_max_fm)
    CALL scatter(lai_fm_g,lai_fm)
  
  END SUBROUTINE forcing_read


!! ================================================================================================================================
!! SUBROUTINE 	: setlai
!!
!>\BRIEF        Routine to force the lai in STOMATE. The code in this routine
!! simply CALCULATES lai and is therefore not functional. The routine should be 
!! rewritten if one wants to force lai.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE setlai(npts,lai)

  !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                    :: npts !! Domain size - number of pixels (unitless)
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)  :: lai  !! PFT leaf area index @tex $(m^{2} m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                               :: j    !! index (unitless)
!_ ================================================================================================================================
    
    !! 1. Set lai for bare soil to zero

    lai(:,ibare_sechiba) = zero

    !! 2. Multiply foliage biomass by sla to calculate lai for all PFTs and pixels

    DO j=2,nvm
       lai(:,j) = biomass(:,j,ileaf)*sla(j)
    ENDDO
    
  END SUBROUTINE setlai

END MODULE stomate
