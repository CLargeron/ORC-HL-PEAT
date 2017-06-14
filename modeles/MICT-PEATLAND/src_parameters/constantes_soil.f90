
! MODULE 	: constantes_soil
!
! CONTACT       : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil" module contains public data for the soils.
!!
!!\n DESCRIPTION : The non saturated hydraulic properties are defined from the  
!!                 formulations of van Genuchten (1980) and Mualem (1976), combined as  
!!                 explained in d'Orgeval (2006). \n
!!                 The related parameters for three main  
!!                 soil textures (coarse, medium and fine) come from Carsel and Parrish  
!!                 (1988).
!!
!! RECENT CHANGE(S): Sonke Zaehle changed hcrit_litter value according to Shilong Piao from 0.03 to 0.08, 080806
!!                   Chloe Largeron ks_peat value for soiltile 4 corresponding to peat soil 
!!
!! REFERENCE(S)	:
!!- Roger A.Pielke, (2002), Mesoscale meteorological modeling, Academic Press Inc. 
!!- Polcher, J., Laval, K., Dümenil, L., Lean, J., et Rowntree, P. R. (1996).
!! Comparing three land surface schemes used in general circulation models. Journal of Hydrology, 180(1-4), 373--394.
!!- Ducharne, A., Laval, K., et Polcher, J. (1998). Sensitivity of the hydrological cycle
!! to the parametrization of soil hydrology in a GCM. Climate Dynamics, 14, 307--327. 
!!- Rosnay, P. de et Polcher, J. (1999). Modelling root water uptake in a complex land surface
!! scheme coupled to a GCM. Hydrol. Earth Syst. Sci., 2(2/3), 239--255.
!!- d'Orgeval, T. et Polcher, J. (2008). Impacts of precipitation events and land-use changes
!! on West African river discharges during the years 1951--2000. Climate Dynamics, 31(2), 249--262. 
!!- Carsel, R. and Parrish, R.: Developing joint probability distributions of soil water
!! retention characteristics, Water Resour. Res.,24, 755–769, 1988.
!!- Mualem Y (1976) A new model for predicting the hydraulic conductivity  
!! of unsaturated porous media. Water Resources Research 12(3):513-522
!!- Van Genuchten M (1980) A closed-form equation for predicting the  
!! hydraulic conductivity of unsaturated soils. Soil Sci Soc Am J, 44(5):892-898
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil

  USE constantes
!-
  IMPLICIT NONE
!
  LOGICAL, SAVE                             :: check_waterbal=.TRUE. !! The check the water balance
  LOGICAL, SAVE                             :: waterbal_error=.FALSE. !! If true the water balance is bad
!
!-
! Dimensioning parameters
!-
  INTEGER(i_std),PARAMETER :: ngrnd=7           !! Number of soil level (unitless)
!Chloe
!   INTEGER(i_std),PARAMETER :: ngrnd=11

    INTEGER(i_std),PARAMETER :: nbdl=11           !! Number of diagnostic levels in the soil (unitless)
!MM : if you want to compare hydrology variables with old TAG 1.6 and lower, 
!     you must set the Number of diagnostic levels in the soil to 6 :
!  INTEGER(i_std),PARAMETER :: nbdl=6

  INTEGER(i_std),PARAMETER :: nslm=11              !! Number of levels in CWRR (unitless)

!!$  REAL(r_std),SAVE :: dpu_max=2.0_r_std            !! Maximum depth of soil reservoir. (m)
!!$                                                   !! If a depth map is given, nbdl, and nslm will be put
!!$                                                   !! to nslm_max and dpu_max will be changed  
!!$                                                   !! in intersurf

  REAL(r_std),SAVE :: dpu_max = 4.0_r_std          !! Maximum depth of soil reservoir. (m)
                                                   !! If a depth map is given, nbdl, and nslm will be put
                                                   !! to nslm_max and dpu_max will be changed  
                                                   !! in intersurf (value for Choisnel/AR5)

  INTEGER(i_std),PARAMETER :: nslm_min=10          !! Number of levels min in CWRR (unitless)
                                                   !! (used only when a depth map is given)
  INTEGER(i_std),PARAMETER :: nslm_max=13          !! Number of levels max in CWRR (unitless)
                                                   !! (used only when a depth map is given)
!
!Number of soil classes
!Chloe : 4ème soiltile : peat for PFT14
  INTEGER(i_std),PARAMETER :: ntext=3           !! Number of soil textures (Silt, Sand, Clay)
  INTEGER(i_std),PARAMETER :: nstm=4            !! Number of soil types (unitless)
!
  CHARACTER(LEN=30) :: soil_classif           !! Type of classification used for the map of soil types.
                                              !! It must be consistent with soil file given by 
                                              !! SOILCLASS_FILE parameter.
  INTEGER(i_std),PARAMETER :: nscm_fao=3     !! For FAO Classification (unitless)
  INTEGER(i_std),PARAMETER :: nscm_usda=12    !! For USDA Classification (unitless)
  INTEGER(i_std), SAVE     :: nscm=nscm_fao   !! Default value for nscm

!-
!- Parameters for soil thermodynamics
!-
  REAL(r_std),PARAMETER :: so_cond = 1.5396     !! Average Thermal Conductivity of soils  @tex $(W.m^{-2}.K^{-1})$ @endtex
  REAL(r_std),PARAMETER :: so_capa = 2.0514e+6  !! Average Heat capacity of soils @tex $(J.m^{-3}.K^{-1})$ @endtex
  REAL(r_std), PARAMETER                :: so_capa_ice = 2.6e6
  ! water and ice density (kg m-3)
  REAL(r_std), PARAMETER                :: rho_water = 1000.
  REAL(r_std), PARAMETER                :: rho_ice = 920.  
  ! Heat capacities (J kg**-1 K**-1)
  ! Isa : useless here, just for info
!  REAL(r_std), PARAMETER		 :: capa_water = 4.188*1.E3
!  REAL(r_std), PARAMETER		 :: capa_ice = 2.228*1.E3
  ! thermal conductivities (W m**-1 K**-1)
  REAL(r_std), PARAMETER		 :: cond_water = 0.6
  REAL(r_std), PARAMETER		 :: cond_ice = 2.2
  REAL(r_std), PARAMETER		 :: cond_solid = 2.32 ! adjusted to obtain csat of 1.89 specified for so_cond_wet at poros 0.15
  ! time constant of long-term soil humidity (s)
  REAL(r_std), PARAMETER	   	 :: tau_freezesoil = 30.*86400.  !1M
  REAL(r_std), PARAMETER                :: tzero = 273.15
  REAL(r_std), PARAMETER                :: lhf = 0.3336*1.E6

  !Chloe+
  REAL(r_std),PARAMETER :: ks_peat = 2120.        !! Saturated Hydraulic conductivity of reference (mm/d) in the case of peat (soiltile 4)
  REAL(r_std),PARAMETER :: mcs_peat = 0.41    !mcs_peat = 0.9002  !! Saturated Moisture content in case of peat 
  REAL(r_std),PARAMETER :: mcr_peat = 0.15        !! Ref PhD Thesis of Q. Dawson p118 - MF Semi Fibrous Peat
  REAL(r_std),PARAMETER :: nvan_peat = 1.38       !! mcr=0.15 : vient de Lett et al pour hemic peat (middle layer)
  REAL(r_std),PARAMETER :: avan_peat = 0.00507       !! 5.07 m-1 = 5.07x10^(-3) mm-1 (unité mm-1 dans orchidee) 
  !Chloe-

  !CHLOE TEST le faux peat : 
  !REAL(r_std),PARAMETER :: ks_peat = 249.6        !! Saturated Hydraulic conductivity of reference (mm/d) in the case of peat (soiltile 4)
  !REAL(r_std),PARAMETER :: mcs_peat = 0.43      !! Saturated Moisture content in case of peat 
  !REAL(r_std),PARAMETER :: mcr_peat = 0.078        !! Ref PhD Thesis of Q. Dawson p118 - MF Semi Fibrous Peat
  !REAL(r_std),PARAMETER :: nvan_peat = 1.56       !! mcr=0.15 : vient de Lett et al pour hemic peat (middle layer)
  !REAL(r_std),PARAMETER :: avan_peat = 0.0036       !! 5.07 m-1 = 5.07x10^(-3) mm-1 (unité mm-1 dans orchidee) 
  !Chloe-



!-
! Values taken from : PIELKE,'MESOSCALE METEOROLOGICAL MODELING',P.384
! Dry soil heat capacity was decreased and conductivity increased.
!-
!
!*REAL(r_std),PARAMETER :: so_capa_dry = 1.35e+6
  REAL(r_std),PARAMETER :: so_capa_dry = 1.80e+6            !! Dry soil Heat capacity of soils @tex $(J.m^{-3}.K^{-1})$ @endtex 
!*REAL(r_std),PARAMETER :: so_cond_dry = 0.28
  REAL(r_std),PARAMETER :: so_cond_dry = 0.40               !! Dry soil Thermal Conductivity of soils @tex $(W.m^{-2}.K^{-1})$ @endtex
!-
  REAL(r_std),PARAMETER :: so_capa_wet = 3.03e+6            !! Wet soil Heat capacity of soils  @tex $(J.m^{-3}.K^{-1})$ @endtex
  REAL(r_std),PARAMETER :: so_cond_wet = 1.89               !! Wet soil Thermal Conductivity of soils @tex $(W.m^{-2}.K^{-1})$ @endtex 
!-
  REAL(r_std),PARAMETER :: sn_cond = 0.3                    !! Thermal Conductivity of snow @tex $(W.m^{-2}.K^{-1})$ @endtex  
  REAL(r_std),PARAMETER :: sn_dens = 330.0                  !! Snow density for the soil thermodynamics (unitless)
  REAL(r_std),PARAMETER :: sn_capa = 2100.0_r_std*sn_dens   !! Heat capacity for snow @tex $(J.m^{-3}.K^{-1})$ @endtex

!-
! Constantes from the Choisnel hydrology
!-
  REAL(r_std),PARAMETER :: qwilt = 5.0              !! Wilting point (Has a numerical role for the moment) (??units??)
  REAL(r_std),PARAMETER :: min_resdis = 2.e-5       !! The minimal size we allow for the upper reservoir (m)
  REAL(r_std),PARAMETER :: min_drain = 0.001        !! Diffusion constant for the slow regime @tex $(kg.m^{-2}.dt^{-1})$ @endtex
                                                    ! (This is for the diffusion between reservoirs)

  REAL(r_std),PARAMETER :: max_drain = 0.1          !! Diffusion constant for the fast regime @tex $(kg.m^{-2}.dt^{-1})$ @endtex
  REAL(r_std),PARAMETER :: exp_drain = 1.5          !! The exponential in the diffusion law (unitless)
  REAL(r_std),SAVE      :: qsintcst = 0.1           !! Transforms leaf area index into size of interception reservoir (unitless)
  REAL(r_std),PARAMETER :: mx_eau_eau = 150.        !! Maximum quantity of water @tex $(kg.m^{-3} of soil)$ @endtex
!-

  REAL(r_std),PARAMETER :: rsol_cste = 33.E3        !! Constant in the computation of resistance for bare soil evaporation
                                                    !! @tex $(s.m^{-2})$ @endtex
  REAL(r_std),PARAMETER :: hcrit_litter=0.08_r_std  !! Scaling depth for litter humidity (m)
!-
! Parameters for soil type distribution
!-
!!$ [DISPENSABLE] REAL(r_std),DIMENSION(nstm),SAVE :: soiltype_default = &    !! Default soil texture distribution in the following order 
!!$ & (/ 0.0, 1.0, 0.0 /)                                        !! sand, loam and clay (0-1, unitless)

!-
! Parameters specific for the CWRR hydrology.
!-
!-
!- 1. Parameters for FAO Classification
!-

!-
! Parameters for soil type distribution
!-
! Default soil texture distribution in the following order :
!   COARSE, MEDIUM, FINE
  REAL(r_std),DIMENSION(nscm_fao),SAVE :: soilclass_default_fao = &
 & (/ 0.28, 0.52, 0.20 /)

 
  REAL(r_std),DIMENSION(nscm_fao*2-1),SAVE :: soilclass_default_fao2 = &
 & (/ 0.28, 0.0, 0.52, 0.0, 0.20 /)

!!$[DISPENSABLE]  INTEGER, SAVE :: jsc_default = 2



 !GROS FAKE TEST CHLOE POUR AGNES : 
 ! REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcr_fao = &            !! Residual soil water content
 !& (/ 0.15_r_std, 0.15_r_std, 0.15_r_std /)

  !REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcs_fao = &            !! Saturated soil water content
 !& (/ 0.9002_r_std, 0.9002_r_std, 0.9002_r_std /)

  ! REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: ks_fao = &             !! Hydraulic conductivity Saturation (mm/d)
 !& (/ 2120.0_r_std, 2120.0_r_std, 2120.0_r_std /)

 !REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: nvan_fao = &          !! Van genuchten coefficient n
 !& (/ 1.38_r_std, 1.38_r_std, 1.38_r_std /)
! Van genuchten coefficient a (cm^{-1}) BIG BUG -> mm^{-1}
  !REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: avan_fao = &
  !& (/ 0.00507_r_std,  0.00507_r_std,  0.00507_r_std /) 



  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: nvan_fao = &          !! Van genuchten coefficient n
 & (/ 1.89_r_std, 1.56_r_std, 1.31_r_std /)
! Van genuchten coefficient a (cm^{-1}) BIG BUG -> mm^{-1}
 REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: avan_fao = &
  & (/ 0.0075_r_std, 0.0036_r_std, 0.0019_r_std /) 
!  & (/ 0.075_r_std, 0.036_r_std, 0.019_r_std /) 
! & (/ 0.036_r_std, 0.036_r_std, 0.036_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcr_fao = &            !! Residual soil water content
 & (/ 0.065_r_std, 0.078_r_std, 0.095_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcs_fao = &            !! Saturated soil water content
 & (/ 0.41_r_std, 0.43_r_std, 0.41_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: ks_fao = &             !! Hydraulic conductivity Saturation (mm/d)
 & (/ 1060.8_r_std, 249.6_r_std, 62.4_r_std /)





! Fraction of saturated volumetric soil moisture above which transpir is max
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: pcent_fao = &
 & (/ 0.5_r_std, 0.5_r_std, 0.5_r_std /)

! Max value of the permeability coeff at the bottom of the soil
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: free_drain_max_fao = &
 & (/ 1.0_r_std, 1.0_r_std, 1.0_r_std /)

! Volumetric water content field capacity
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcf_fao = &
 & (/ 0.32_r_std, 0.32_r_std, 0.32_r_std /)

! Volumetric water content Wilting pt
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcw_fao = &
 & (/ 0.10_r_std, 0.10_r_std, 0.10_r_std /)
! & (/ 0.07_r_std, 0.085_r_std, 0.10_r_std /)

! Vol. wat. cont. above which albedo is cst
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mc_awet_fao = &
 & (/ 0.25_r_std, 0.25_r_std, 0.25_r_std /)

! Vol. wat. cont. below which albedo is cst
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mc_adry_fao = &
 & (/ 0.1_r_std, 0.1_r_std, 0.1_r_std /)

! Matrix potential at saturation (mm)
  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: psis_fao = &
 & (/ -300.0_r_std, -300.0_r_std, -300.0_r_std /)
!-
!- 2. Parameters for USDA Classification
!-

!-
! Parameters for soil type distribution
!-
! Default soil texture distribution in the following order :
!    sand, loam and clay ??? OR COARSE, MEDIUM, FINE???
  REAL(r_std),DIMENSION(nscm_usda),SAVE :: soilclass_default_usda = &
 & (/ 0.28, 0.52, 0.20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

! Van genuchten coefficient n
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: nvan_usda = &
 & (/ 2.68_r_std, 2.28_r_std, 1.89_r_std, 1.41_r_std, &
 &    1.37_r_std, 1.56_r_std, 1.48_r_std, 1.23_r_std, &
 &    1.31_r_std, 1.23_r_std, 1.09_r_std, 1.09_r_std /)
! Van genuchten coefficient a (cm^{-1}) BIG BUG!!! -> mm^{-1}
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: avan_usda = &
 & (/ 0.0145_r_std, 0.0124_r_std, 0.0075_r_std, 0.0020_r_std, &
 &    0.0016_r_std, 0.0036_r_std, 0.0059_r_std, 0.0010_r_std, &
 &    0.0019_r_std, 0.0027_r_std, 0.0005_r_std, 0.0008_r_std /)
! & (/ 0.145_r_std, 0.124_r_std, 0.075_r_std, 0.020_r_std, &
! &    0.016_r_std, 0.036_r_std, 0.059_r_std, 0.010_r_std, &
! &    0.019_r_std, 0.027_r_std, 0.005_r_std, 0.008_r_std /)
! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay
! Residual soil water content
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcr_usda = &
 & (/ 0.045_r_std, 0.057_r_std, 0.065_r_std, 0.067_r_std, &
 &    0.034_r_std, 0.078_r_std, 0.100_r_std, 0.089_r_std, &
 &    0.095_r_std, 0.100_r_std, 0.070_r_std, 0.068_r_std /)
! Saturated soil water content
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcs_usda = &
 & (/ 0.43_r_std, 0.41_r_std, 0.41_r_std, 0.45_r_std, &
 &    0.46_r_std, 0.43_r_std, 0.39_r_std, 0.43_r_std, &
 &    0.41_r_std, 0.38_r_std, 0.36_r_std, 0.38_r_std /)
! Hydraulic conductivity Saturation (mm/d)
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: ks_usda = &
 & (/ 7128.0_r_std, 3501.6_r_std, 1060.8_r_std, 108.0_r_std, &
 &    60.0_r_std, 249.6_r_std, 314.4_r_std, 16.8_r_std, &
 &    62.4_r_std, 28.8_r_std, 4.8_r_std, 48.0_r_std /)
! Soil moisture above which transpir is max
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: pcent_usda = &
 & (/ 0.5_r_std, 0.5_r_std, 0.5_r_std, 0.5_r_std, &
 &    0.5_r_std, 0.5_r_std, 0.5_r_std, 0.5_r_std, &
 &    0.5_r_std, 0.5_r_std, 0.5_r_std, 0.5_r_std /)
! Max value of the permeability coeff at the bottom of the soil
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: free_drain_max_usda = &
 & (/ 1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std /)
! Volumetric water content field capacity
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcf_usda = &
 & (/ 0.32_r_std, 0.32_r_std, 0.32_r_std, 0.32_r_std, &
 &    0.32_r_std, 0.32_r_std, 0.32_r_std, 0.32_r_std, &
 &    0.32_r_std, 0.32_r_std, 0.32_r_std, 0.32_r_std /)
! Volumetric water content Wilting pt
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcw_usda = &
 & (/ 0.10_r_std, 0.10_r_std, 0.10_r_std, 0.10_r_std, &
 &    0.10_r_std, 0.10_r_std, 0.10_r_std, 0.10_r_std, &
 &    0.10_r_std, 0.10_r_std, 0.10_r_std, 0.10_r_std /)
! Vol. wat. cont. above which albedo is cst
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_awet_usda = &
 & (/ 0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std /)
! Vol. wat. cont. below which albedo is cst
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_adry_usda = &
 & (/ 0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std /)
! Matrix potential at saturation (mm)
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: psis_usda = &
 & (/ -300.0_r_std, -300.0_r_std, -300.0_r_std, -300.0_r_std, &
 &    -300.0_r_std, -300.0_r_std, -300.0_r_std, -300.0_r_std, &
 &    -300.0_r_std, -300.0_r_std, -300.0_r_std, -300.0_r_std /)

  INTEGER(i_std),PARAMETER :: imin = 1                        !! CWRR linearisation (unitless)
  INTEGER(i_std),PARAMETER :: nbint = 50                      !! Number of interval for CWRR (unitless)
  INTEGER(i_std),PARAMETER :: imax = nbint+1                  !! Number of points for CWRR (unitless)
  REAL(r_std),PARAMETER :: w_time = 1.0_r_std                 !! Time weighting for discretisation (unitless)
!-
! Diagnostic variables
!-
  REAL(r_std),DIMENSION(nbdl),SAVE :: diaglev     !! The lower limit of the layer on which soil moisture (relative)
                                                  !! and temperature are going to be diagnosed.
                                                  !! These variables are made for transfering the information
                                                  !! to the biogeophyical processes modelled in STOMATE. (unitless)
!-------------------------
END MODULE constantes_soil
