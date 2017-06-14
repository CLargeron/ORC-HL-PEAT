! =================================================================================================================================
! MODULE       : pft_parameters
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module initializes all the pft parameters in function of the
!! number of vegetation types and of the values chosen by the user.
!!
!!\n DESCRIPTION:  This module allocates and initializes the pft parameters in function of the number of pfts
!!                 and the values of the parameters. \n
!!                 The number of PFTs is read in intersurf.f90 (subroutine intsurf_config). \n
!!                 Then we can initialize the parameters. \n
!!                 This module is the result of the merge of constantes_co2, constantes_veg, stomate_constants.\n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2012-07-26 16:16:59 +0200 (Thu, 26 Jul 2012) $
!! $Revision: 956 $
!! \n
!_ ================================================================================================================================

MODULE pft_parameters

  USE constantes_mtc
  USE constantes
  USE ioipsl
  USE parallel
  USE defprec
  
  IMPLICIT NONE


  !
  ! PFT GLOBAL
  !
  INTEGER(i_std), SAVE :: nvm = 13                               !! Number of vegetation types (2-N, unitless)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pft_to_mtc  !! Table of conversion : we associate one pft to one metaclass 
                                                                 !! (1-13, unitless)

  CHARACTER(LEN=34), ALLOCATABLE, SAVE, DIMENSION(:) :: PFT_name !! Description of the PFT (unitless)

  LOGICAL, SAVE   :: l_first_pft_parameters = .TRUE.             !! To keep first call trace of the module (true/false)
  LOGICAL, SAVE   :: ok_throughfall_by_pft = .FALSE.             !! Flag to use the parameter PERCENT_THROUGHFALL_PFT (true/false) 


  !
  ! VEGETATION STRUCTURE
  !
  !-
  !  1. Sechiba
  !
  !   1.1 Labels - Characteristics
  !-
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_tree               !! Is the vegetation type a tree ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_deciduous          !! Is PFT deciduous ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_evergreen          !! Is PFT evegreen ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_summergreen        !! Is PFT summergreen ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_needleleaf         !! Is PFT needleleaf ? (true/false)
 
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_broadleaf          !! Is PFT broadleaf ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_tropical           !! Is PFT tropical ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_temperate          !! Is PFT temperate ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_boreal             !! Is PFT boreal ? (true/false)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_c3                 !! Is PFT c3 ? (true/false)

  CHARACTER(len=5), ALLOCATABLE, SAVE, DIMENSION(:) :: type_of_lai  !! Type of behaviour of the LAI evolution algorithm 
                                                                    !! for each vegetation type.
                                                                    !! Value of type_of_lai, one for each vegetation type :
                                                                    !! mean or interp
  !-
  !   1.2 Prescribed Values
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: veget_ori_fixed_test_1 !! Value for veget_ori for tests in 0-dim simulations 
                                                                         !! (0-1, unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: llaimax                !! laimax for maximum lai see also type of lai 
                                                                         !! interpolation (m^2.m^{-2}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: llaimin                !! laimin for minimum lai see also type of lai 
                                                                         !! interpolation (m^2.m^{-2})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: height_presc           !! prescribed height of vegetation.(m)
                                                                         !! Value for height_presc : one for each vegetation type

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) ::  rveg_pft              !! Potentiometer to set vegetation resistance (unitless)
                                                                         !! Nathalie on March 28th, 2006 - from Fred Hourdin,
  !-
  !  2. Stomate 
  !
  !   2.1 Labels - Characteristics
  !- 
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: natural                    !! natural? (true/false)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_tab            !! leaf type (1-4, unitless)
                                                                         !! 1=broad leaved tree, 2=needle leaved tree, 
                                                                         !! 3=grass 4=bare ground 
  !-
  !   2.2 Prescribed Values
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: sla                    !! specif leaf area (m^2.gC^{-1}) 


  !
  ! EVAPOTRANSPIRATION (sechiba)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rstruct_const          !! Structural resistance. (s.m^{-1})
                                                                         !! Value for rstruct_const : one for each vegetation type

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: kzero                  !! A vegetation dependent constant used in the calculation
                                                                         !! of the surface resistance. (kg.m^2.s^{-1})
                                                                         !! Value for kzero one for each vegetation type


  !
  ! WATER (sechiba)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: wmax_veg  !! Maximum field capacity for each of the vegetations (Temporary).
                                                            !! Value of wmax_veg : max quantity of water :
                                                            !! one for each vegetation type (kg.m^{-3})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: humcste   !! Root profile description for the different vegetation types.
                                                            !! These are the factor in the exponential which gets
                                                            !! the root density as a function of depth (m^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: throughfall_by_pft !! Fraction of rain intercepted by the canopy (0-100, unitless)


  !
  ! ALBEDO (sechiba)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_ini     !! Initial snow albedo value for each vegetation type
                                                                !! as it will be used in condveg_snow (unitless)
                                                                !! Source : Values are from the Thesis of S. Chalita (1992)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_dec     !! Decay rate of snow albedo value for each vegetation type
                                                                !! as it will be used in condveg_snow (unitless)
                                                                !! Source : Values are from the Thesis of S. Chalita (1992)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alb_leaf_vis  !! leaf albedo of vegetation type, visible albedo (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alb_leaf_nir  !! leaf albedo of vegetation type, near infrared albedo (unitless)


  !
  ! SOIL - VEGETATION
  !
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pref_soil_veg      !! Table which contains the correlation between the soil
                                                                        !! types and vegetation type. Two modes exist :
                                                                        !! 1) pref_soil_veg = 0 then we have an equidistribution
                                                                        !!    of vegetation on soil types
                                                                        !! 2) Else for each pft the prefered soil type is given :
                                                                        !!    1=sand, 2=loan, 3=clay
                                                                        !! This variable is initialized in slowproc.(1-3, unitless)

  !
  ! PHOTOSYNTHESIS
  !
  !-
  ! 1. CO2
  !-
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_c4             !! flag for C4 vegetation types (true/false)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gsslope       !! Slope of the gs/A relation (Ball & al.) (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gsoffset      !! intercept of the gs/A relation (Ball & al.) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: vcmax_fix     !! values used for vcmax when STOMATE is not activated
                                                                !! (�mol.m^{-2}.s^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: vjmax_fix     !! values used for vjmax when STOMATE is not activated 
                                                                !! (�mol.m^{-2}.s^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: co2_tmin_fix  !! values used for photosynthesis tmin when STOMATE 
                                                                !! is not activated (C)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: co2_topt_fix  !! values used for photosynthesis topt when STOMATE 
                                                                !! is not activated (C)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: co2_tmax_fix  !! values used for photosynthesis tmax when STOMATE 
                                                                !! is not activated (C)
  !-
  ! 2. Stomate
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff     !! extinction coefficient of the Monsi&Saeki relationship (1953)
                                                                !! (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: vcmax_opt     !! Maximum rate of carboxylation (�mol.m^{-2}.s^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: vjmax_opt     !! Maximum rate of RUbp regeneration (�mol.m^{-2}.s^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_min_a  !! minimum photosynthesis temperature, 
                                                                !! constant a of ax^2+bx+c (deg C),tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_min_b  !! minimum photosynthesis temperature, 
                                                                !! constant b of ax^2+bx+c (deg C),tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_min_c  !! minimum photosynthesis temperature, 
                                                                !! constant c of ax^2+bx+c (deg C),tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_opt_a  !! optimum photosynthesis temperature, 
                                                                !! constant a of ax^2+bx+c (deg C),tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_opt_b  !! optimum photosynthesis temperature,
                                                                !! constant b of ax^2+bx+c (deg C),tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_opt_c  !! optimum photosynthesis temperature, 
                                                                !! constant c of ax^2+bx+c (deg C),tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_max_a  !! maximum photosynthesis temperature, 
                                                                !! constant a of ax^2+bx+c (deg C), tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_max_b  !! maximum photosynthesis temperature, 
                                                                !! constant b of ax^2+bx+c (deg C), tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_max_c  !! maximum photosynthesis temperature, 
                                                                !! constant c of ax^2+bx+c (deg C), tabulated (unitless)


  !
  ! RESPIRATION (stomate)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: maint_resp_slope  !! slope of maintenance respiration coefficient 
                                                                      !! (1/K, 1/K^2, 1/K^3), used in the code

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_c  !! slope of maintenance respiration coefficient (1/K),
                                                                      !! constant c of aT^2+bT+c , tabulated

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_b  !! slope of maintenance respiration coefficient (1/K), 
                                                                      !! constant b of aT^2+bT+c , tabulated

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_a  !! slope of maintenance respiration coefficient (1/K), 
                                                                      !! constant a of aT^2+bT+c , tabulated

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: coeff_maint_zero  !! maintenance respiration coefficient at 0 deg C, 
                                                                      !! used in the code (gC.gC^{-1}.day^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_leaf        !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for leaves, tabulated (gC.gC^{-1}.day^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_sapabove    !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for sapwood above, tabulated (gC.gC^{-1}.day^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_sapbelow    !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for sapwood below, tabulated (gC.gC^{-1}.day^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_heartabove  !! maintenance respiration coefficient at 0 deg C
                                                                      !! for heartwood above, tabulated (gC.gC^{-1}.day^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_heartbelow  !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for heartwood below, tabulated (gC.gC^{-1}.day^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_root        !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for roots, tabulated (gC.gC^{-1}.day^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_fruit       !! maintenance respiration coefficient  at 0 deg C,
                                                                      !! for fruits, tabulated (gC.gC^{-1}.day^{-1})  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cm_zero_carbres     !! maintenance respiration coefficient at 0 deg C,
                                                                      !! for carbohydrate reserve, tabulated (gC.gC^{-1}.day^{-1}) 

 
  !
  ! FIRE (stomate)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: flam              !! flamability : critical fraction of water holding 
                                                                    !! capacity (0-1, unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: resist            !! fire resistance (0-1, unitless)


  !
  ! FLUX - LUC (Land Use Change)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_1   !! Coeff of biomass export for the year (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_10  !! Coeff of biomass export for the decade (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_100 !! Coeff of biomass export for the century (unitless)
 
 
  !
  ! PHENOLOGY
  !
  !-
  ! 1. Stomate
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_max           !! maximum LAI, PFT-specific (m^2.m^{-2})

  CHARACTER(len=6), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_model  !! which phenology model is used? (tabulated) (unitless)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_type     !! type of phenology (0-4, unitless)
                                                                    !! 0=bare ground 1=evergreen,  2=summergreen, 
                                                                    !! 3=raingreen,  4=perennial
                                                                    !! For the moment, the bare ground phenotype is not managed, 
                                                                    !! so it is considered as "evergreen"
  !-
  ! 2. Leaf Onset
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: pheno_gdd_crit   !! critical gdd,tabulated (C), used in the code

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_c   !! critical gdd,tabulated (C), 
                                                                     !! constant c of aT^2+bT+c (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_b   !! critical gdd,tabulated (C), 
                                                                     !! constant b of aT^2+bT+c (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_a   !! critical gdd,tabulated (C), 
                                                                     !! constant a of aT^2+bT+c (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ngd_crit           !! critical ngd,tabulated. Threshold -5 degrees (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ncdgdd_temp        !! critical temperature for the ncd vs. gdd function
                                                                     !! in phenology (C)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: hum_frac           !! critical humidity (relative to min/max) for phenology
                                                                     !! (0-1, unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lowgpp_time        !! minimum duration of dormance (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: hum_min_time       !! minimum time elapsed since moisture minimum (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tau_sap            !! sapwood -> heartwood conversion time (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tau_fruit          !! fruit lifetime (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ecureuil           !! fraction of primary leaf and root allocation put
                                                                     !! into reserve (0-1, unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alloc_min          !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alloc_max          !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: demi_alloc         !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaflife_tab       !! leaf longevity, tabulated (??units??)
  !-
  ! 3. Senescence
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaffall              !! length of death of leaves,tabulated (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leafagecrit           !! critical leaf age,tabulated (days)

  CHARACTER(len=6), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_type  !! type of senescence,tabulated (unitless)
                                                                        !! List of avaible types of senescence :
                                                                        !! 'cold  ', 'dry   ', 'mixed ', 'none  '

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_hum        !! critical relative moisture availability for senescence
                                                                        !! (0-1, unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: nosenescence_hum      !! relative moisture availability above which there is
                                                                        !! no humidity-related senescence (0-1, unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_turnover_time     !! maximum turnover time for grasses (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_turnover_time     !! minimum turnover time for grasses (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_leaf_age_for_senescence  !! minimum leaf age to allow senescence g (days)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: senescence_temp     !! critical temperature for senescence (C),
                                                                        !! used in the code

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_c     !! critical temperature for senescence (C), 
                                                                        !! constant c of aT^2+bT+c , tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_b     !! critical temperature for senescence (C), 
                                                                        !! constant b of aT^2+bT+c , tabulated (unitless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_a     !! critical temperature for senescence (C),
                                                                        !! constant a of aT^2+bT+c , tabulated (unitless)


  !
  ! DGVM
  !

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: residence_time        !! residence time of trees (y) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tmin_crit             !! critical tmin, tabulated (C)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tcm_crit              !! critical tcm, tabulated (C)

  !
  ! Biogenic Volatile Organic Compounds
  !

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_isoprene       !! Isoprene emission factor
                                                                           !! (\mu gC.g^{-1}.h^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_monoterpene    !! Monoterpene emission factor
                                                                           !! (\mu gC.g^{-1}.h^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_ORVOC          !! ORVOC emissions factor
                                                                           !! (\mu gC.g^{-1}.h^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_OVOC           !! OVOC emissions factor
                                                                           !! (\mu gC.g^{-1}.h^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_MBO            !! MBO emissions factor 
                                                                           !! (\mu gC.g^{-1}.h^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_methanol       !! Methanol emissions factor
                                                                           !! (\mu gC.g^{-1}.h^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetone        !! Acetone emissions factor 
                                                                           !! (\mu gC.g^{-1}.h^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetal         !! Acetaldehyde emissions factor
                                                                           !! (\mu gC.g^{-1}.h^{-1}) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_formal         !! Formaldehyde emissions factor
                                                                           !! (\mu gC.g^{-1}.h^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetic         !! Acetic Acid emissions factor 
                                                                           !! (\mu gC.g^{-1}.h^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_formic         !! Formic Acid emissions factor 
                                                                           !! (\mu gC.g^{-1}.h^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_no_wet         !! NOx emissions factor soil emissions and 
                                                                           !! exponential dependancy factor for wet soils
                                                                           !! (ngN.m^{-2}.s^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_no_dry         !! NOx emissions factor soil emissions and
                                                                           !! exponential dependancy factor for dry soils
                                                                           !! (ngN.m^{-2}.s^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Larch                    !! Larcher 1991 SAI/LAI ratio (unitless)

  !
  ! INTERNAL PARAMETERS USED IN STOMATE_DATA
  !

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_initmin   !! Initial lai for trees/grass  (m^2.m^{-2})

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: tree              !! is pft a tree? (used for consistency with is_tree)
                                                                !! (true/false)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: bm_sapl     !! sapling biomass (gC.ind^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: migrate       !! migration speed (m.year^{-1})

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maxdia        !! maximum stem diameter from which on crown area no longer 
                                                                !! increases (m)m

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cn_sapl       !! crown of tree when sapling (m^2)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_timecst  !! time constant for leaf age discretisation (days)


CONTAINS
 !

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_main
!!
!>\BRIEF          This subroutine initializes all the pft parameters in function of the
!! number of vegetation types chosen by the user.
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the options 
!!                activated by the user in the configuration files. (structure active_flags) \n
!!                The allocation is done just before reading the correspondence table  between PFTs and MTCs
!!                defined by the user in the configuration file.\n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, ok_co2, routing, new_hydrol...) in order to
!!                optimize the memory allocation. \n
!!                If the number of PFTs and pft_to_mtc are not found, the standard configuration will be used
!!                (13 PFTs, PFT = MTC). \n 
!!                Some restrictions : the pft 1 can only be the bare soil and it is unique. \n
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

 SUBROUTINE pft_parameters_main(active_flags)

   IMPLICIT NONE

   !! 0. Variables and parameters declaration

   !! 0.1 Input variables 

   TYPE(control_type),INTENT(in) :: active_flags   !! What parts of the code are activated ? (true/false)
   
   !! 0.4 Local variables  

   INTEGER(i_std) :: j                             !! Index (unitless)

!_ ================================================================================================================================ 
   
   !
   ! PFT global
   !

   IF(l_first_pft_parameters) THEN

      !! 1. First time step
      IF(long_print) THEN
         WRITE(numout,*) 'l_first_pft_parameters :we read the parameters from the def files'
      ENDIF

      IF ( active_flags%hydrol_cwrr ) THEN
         
         !! 2.1 Read the flag ok_throughfall_by_pft to know if 
         !!      we have to use the parameter throughfall_by_pft

         !Config Key   = OK_THROUGHFALL_PFT
         !Config Desc  = Activate use of PERCENT_THROUGHFALL_PFT
         !Config If    = HYDROL_CWRR
         !Config Def   = FALSE
         !Config Help  = If NOT OFF_LINE_MODE it is always TRUE (coupled with a GCM)
         !Config Units = [FLAG]
         IF ( .NOT. OFF_LINE_MODE ) ok_throughfall_by_pft = .TRUE.
         CALL getin_p('OK_THROUGHFALL_PFT',ok_throughfall_by_pft)   

      END IF
   
      !! 2.2 Memory allocation for the pfts-parameters
      CALL pft_parameters_alloc(active_flags)

      !! 3. Correspondance table 
      
      !! 3.1 Initialisation of the correspondance table
      !! Initialisation of the correspondance table
      IF (nvm == nvmc) THEN
         WRITE(numout,*) 'Message to the user : we will use ORCHIDEE to its standard configuration' 
         pft_to_mtc = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 /)
      ELSE
         pft_to_mtc(:) = undef_int
      ENDIF !(nvm  == nvmc)
      
      !! 3.2 Reading of the conrrespondance table in the .def file
      !
      !Config Key   = PFT_TO_MTC
      !Config Desc  = correspondance array linking a PFT to MTC
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('PFT_TO_MTC',pft_to_mtc)
      
      !! 3.3 If the user want to use the standard configuration, he needn't to fill the correspondance array
      !!     If the configuration is wrong, send a error message to the user.
      IF(nvm /= nvmc ) THEN
         !
         IF(pft_to_mtc(1) == undef_int) THEN
            STOP ' The array PFT_TO_MTC is empty : we stop'
         ENDIF !(pft_to_mtc(1) == undef_int)
         !
      ENDIF !(nvm /= nvmc )

      !! 3.4 Some error messages

      !! 3.4.1 What happened if pft_to_mtc(j) > nvmc or pft_to_mtc(j) <=0 (if the mtc doesn't exist)?
       DO j = 1, nvm ! Loop over # PFTs  
          !
          IF( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) ) THEN
             WRITE(numout,*) 'the metaclass chosen does not exist'
             STOP 'we stop reading pft_to_mtc'
          ENDIF !( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) )
          !
       ENDDO  ! Loop over # PFTs  


       !! 3.4.2 Check if pft_to_mtc(1) = 1 
       IF(pft_to_mtc(1) /= 1) THEN
          !
          WRITE(numout,*) 'the first pft has to be the bare soil'
          STOP 'we stop reading next values of pft_to_mtc'
          !
       ELSE
          !
          DO j = 2,nvm ! Loop over # PFTs different from bare soil
             !
             IF(pft_to_mtc(j) == 1) THEN
                WRITE(numout,*) 'only pft_to_mtc(1) has to be the bare soil'
                STOP 'we stop reading pft_to_mtc'
             ENDIF ! (pft_to_mtc(j) == 1)
             !
          ENDDO ! Loop over # PFTs different from bare soil
          !
       ENDIF !(pft_to_mtc(1) /= 1)
      

      !! 4.Initialisation of the pfts-parameters
      CALL pft_parameters_init(active_flags)

      !! 5. Useful data

      !! 5.1 Read the name of the PFTs given by the user
      !
      !Config Key   = PFT_NAME
      !Config Desc  = Name of a PFT
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = bare ground, tropical broad-leaved evergreen, tropical broad-leaved raingreen, 
      !Config         temperate needleleaf evergreen, temperate broad-leaved evergreen temperate broad-leaved summergreen,
      !Config         boreal needleleaf evergreen, boreal broad-leaved summergreen, boreal needleleaf summergreen,
      !Config         C3 grass, C4 grass, C3 agriculture, C4 agriculture    
      !Config Help  = the user can name the new PFTs he/she introducing for new species
      !Config Units = [-]
      CALL getin_p('PFT_NAME',pft_name)

      !! 5.2 A useful message to the user: correspondance between the number of the pft
      !! and the name of the associated mtc 
      DO j = 1,nvm ! Loop over # PFTs
         !
         WRITE(numout,*) 'the PFT',j, 'called  ', PFT_name(j),'corresponds to the MTC : ',MTC_name(pft_to_mtc(j))
         !
      ENDDO ! Loop over # PFTs


      !! 6. End message
      IF(long_print) THEN
         WRITE(numout,*) 'pft_parameters_done'
      ENDIF

      !! 8. Reset flag
      l_first_pft_parameters = .FALSE.

   ELSE 

      RETURN

   ENDIF !(l_first_pft_parameters)

 END SUBROUTINE pft_parameters_main
 !
 !=
 !

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_init 
!!
!>\BRIEF          This subroutine initializes all the pft parameters by the default values
!! of the corresponding metaclasse. 
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the correspondence
!!                table defined by the user in the configuration files. \n
!!                With the correspondence table, the subroutine can search the default values for the parameter
!!                even if the PFTs are classified in a random order (except bare soil). \n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, ok_co2, routing, new_hydrol...).\n
!!
!! RECENT CHANGE(S): Didier Solyga : Simplified PFT loops : use vector notation. 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE pft_parameters_init(active_flags)
  
   IMPLICIT NONE
   
   !! 0. Variables and parameters declaration

   !! 0.1 Input variables
   
   TYPE(control_type),INTENT(in) :: active_flags  !! What parts of the code are activated ? (true/false)

   !! 0.4 Local variables

!_ ================================================================================================================================ 

   !
   ! 1. Correspondance between the PFTs values and thes MTCs values 
   !
 

   ! 1.1 For parameters used anytime
   
   PFT_name(:) = MTC_name(pft_to_mtc(:))
   !
   ! Vegetation structure 
   !
   veget_ori_fixed_test_1(:) = veget_ori_fixed_mtc(pft_to_mtc(:))
   llaimax(:) = llaimax_mtc(pft_to_mtc(:))
   llaimin(:) = llaimin_mtc(pft_to_mtc(:))
   height_presc(:) = height_presc_mtc(pft_to_mtc(:))
   type_of_lai(:) = type_of_lai_mtc(pft_to_mtc(:))
   is_tree(:) = is_tree_mtc(pft_to_mtc(:))
   natural(:) = natural_mtc(pft_to_mtc(:))
   is_deciduous(:) = is_deciduous_mtc(pft_to_mtc(:))
   is_evergreen(:) = is_evergreen_mtc(pft_to_mtc(:))
   is_c3(:) = is_c3_mtc(pft_to_mtc(:))
   !
   ! Water - sechiba
   !
   If (active_flags%hydrol_cwrr ) THEN
      humcste(:) = humcste_cwrr(pft_to_mtc(:)) ! values for 2m soil depth
   ELSE
      humcste(:) = humcste_mtc(pft_to_mtc(:))  ! values for 4m soil depth 
   END IF
   !
   ! Soil - vegetation
   !
   pref_soil_veg(:) = pref_soil_veg_mtc(pft_to_mtc(:))
   !
   ! Photosynthesis
   !
   is_c4(:) = is_c4_mtc(pft_to_mtc(:))
   gsslope(:) = gsslope_mtc(pft_to_mtc(:))
   gsoffset(:) = gsoffset_mtc(pft_to_mtc(:))
   vcmax_fix(:) = vcmax_fix_mtc(pft_to_mtc(:))
   vjmax_fix(:) = vjmax_fix_mtc(pft_to_mtc(:))
   co2_tmin_fix(:) = co2_tmin_fix_mtc(pft_to_mtc(:))
   co2_topt_fix(:) = co2_topt_fix_mtc(pft_to_mtc(:))
   co2_tmax_fix(:) = co2_tmax_fix_mtc(pft_to_mtc(:))
   ext_coeff(:) = ext_coeff_mtc(pft_to_mtc(:))

   ! 1.2 For sechiba parameters

   IF (active_flags%ok_sechiba) THEN
      !
      ! Vegetation structure - sechiba
      !
      rveg_pft(:) = rveg_mtc(pft_to_mtc(:))
      !
      ! Evapotranspiration -  sechiba
      !
      rstruct_const(:) = rstruct_const_mtc(pft_to_mtc(:))
      kzero(:) = kzero_mtc(pft_to_mtc(:))
      !
      ! Water - sechiba
      !
      wmax_veg(:) = wmax_veg_mtc(pft_to_mtc(:))
      IF ( .NOT.(active_flags%hydrol_cwrr) .OR.  (active_flags%hydrol_cwrr .AND. ok_throughfall_by_pft) ) THEN
         throughfall_by_pft(:) = throughfall_by_mtc(pft_to_mtc(:))
      ENDIF
      !
      ! Albedo - sechiba
      !
      snowa_ini(:) = snowa_ini_mtc(pft_to_mtc(:))
      snowa_dec(:) = snowa_dec_mtc(pft_to_mtc(:)) 
      alb_leaf_vis(:) = alb_leaf_vis_mtc(pft_to_mtc(:))  
      alb_leaf_nir(:) = alb_leaf_nir_mtc(pft_to_mtc(:))
      !-
   ENDIF !(active_flags%ok_sechiba)

   ! 1.3 For BVOC parameters
   
   IF (active_flags%ok_inca) THEN
      !
      ! Biogenic Volatile Organic Compounds
      !
      em_factor_isoprene(:) = em_factor_isoprene_mtc(pft_to_mtc(:))
      em_factor_monoterpene(:) = em_factor_monoterpene_mtc(pft_to_mtc(:))
      em_factor_ORVOC(:) = em_factor_ORVOC_mtc(pft_to_mtc(:)) 
      em_factor_OVOC(:) = em_factor_OVOC_mtc(pft_to_mtc(:))
      em_factor_MBO(:) = em_factor_MBO_mtc(pft_to_mtc(:))
      em_factor_methanol(:) = em_factor_methanol_mtc(pft_to_mtc(:))
      em_factor_acetone(:) = em_factor_acetone_mtc(pft_to_mtc(:)) 
      em_factor_acetal(:) = em_factor_acetal_mtc(pft_to_mtc(:))
      em_factor_formal(:) = em_factor_formal_mtc(pft_to_mtc(:))
      em_factor_acetic(:) = em_factor_acetic_mtc(pft_to_mtc(:))
      em_factor_formic(:) = em_factor_formic_mtc(pft_to_mtc(:))
      em_factor_no_wet(:) = em_factor_no_wet_mtc(pft_to_mtc(:))
      em_factor_no_dry(:) = em_factor_no_dry_mtc(pft_to_mtc(:))
      Larch(:) = Larch_mtc(pft_to_mtc(:)) 
      !-
   ENDIF !(active_flags%ok_inca)

   ! 1.4 For stomate parameters

   IF (active_flags%ok_stomate) THEN
      !
      ! Vegetation structure - stomate
      !
      leaf_tab(:) = leaf_tab_mtc(pft_to_mtc(:))
      sla(:) = sla_mtc(pft_to_mtc(:))
      !
      ! Photosynthesis
      !
      vcmax_opt(:) = vcmax_opt_mtc(pft_to_mtc(:))
      vjmax_opt(:) = vjmax_opt_mtc(pft_to_mtc(:)) 
      tphoto_min_a(:) = tphoto_min_a_mtc(pft_to_mtc(:)) 
      tphoto_min_b(:) = tphoto_min_b_mtc(pft_to_mtc(:))
      tphoto_min_c(:) = tphoto_min_c_mtc(pft_to_mtc(:))
      tphoto_opt_a(:) = tphoto_opt_a_mtc(pft_to_mtc(:))
      tphoto_opt_b(:) = tphoto_opt_b_mtc(pft_to_mtc(:))
      tphoto_opt_c(:) = tphoto_opt_c_mtc(pft_to_mtc(:))
      tphoto_max_a(:) = tphoto_max_a_mtc(pft_to_mtc(:))
      tphoto_max_b(:) = tphoto_max_b_mtc(pft_to_mtc(:))
      tphoto_max_c(:) = tphoto_max_c_mtc(pft_to_mtc(:))
      !
      ! Respiration - stomate
      !
      maint_resp_slope_c(:) = maint_resp_slope_c_mtc(pft_to_mtc(:))               
      maint_resp_slope_b(:) = maint_resp_slope_b_mtc(pft_to_mtc(:))
      maint_resp_slope_a(:) = maint_resp_slope_a_mtc(pft_to_mtc(:))
      cm_zero_leaf(:) = cm_zero_leaf_mtc(pft_to_mtc(:))
      cm_zero_sapabove(:) = cm_zero_sapabove_mtc(pft_to_mtc(:))
      cm_zero_sapbelow(:) = cm_zero_sapbelow_mtc(pft_to_mtc(:)) 
      cm_zero_heartabove(:) = cm_zero_heartabove_mtc(pft_to_mtc(:)) 
      cm_zero_heartbelow(:) = cm_zero_heartbelow_mtc(pft_to_mtc(:))
      cm_zero_root(:) = cm_zero_root_mtc(pft_to_mtc(:))
      cm_zero_fruit(:) = cm_zero_fruit_mtc(pft_to_mtc(:))
      cm_zero_carbres(:) = cm_zero_carbres_mtc(pft_to_mtc(:))
      !
      ! Fire - stomate
      !
      flam(:) = flam_mtc(pft_to_mtc(:))
      resist(:) = resist_mtc(pft_to_mtc(:))
      !
      ! Flux - LUC
      !
      coeff_lcchange_1(:) = coeff_lcchange_1_mtc(pft_to_mtc(:))
      coeff_lcchange_10(:) = coeff_lcchange_10_mtc(pft_to_mtc(:))
      coeff_lcchange_100(:) = coeff_lcchange_100_mtc(pft_to_mtc(:))
      !
      ! Phenology
      !
      !
      ! 1. Stomate
      !
      lai_max(:) = lai_max_mtc(pft_to_mtc(:))
      pheno_model(:) = pheno_model_mtc(pft_to_mtc(:))
      pheno_type(:) = pheno_type_mtc(pft_to_mtc(:))
      !
      ! 2. Leaf Onset
      !
      pheno_gdd_crit_c(:) = pheno_gdd_crit_c_mtc(pft_to_mtc(:))
      pheno_gdd_crit_b(:) = pheno_gdd_crit_b_mtc(pft_to_mtc(:))         
      pheno_gdd_crit_a(:) = pheno_gdd_crit_a_mtc(pft_to_mtc(:))
      ngd_crit(:) =  ngd_crit_mtc(pft_to_mtc(:))
      ncdgdd_temp(:) = ncdgdd_temp_mtc(pft_to_mtc(:)) 
      hum_frac(:) = hum_frac_mtc(pft_to_mtc(:))
      lowgpp_time(:) = lowgpp_time_mtc(pft_to_mtc(:))
      hum_min_time(:) = hum_min_time_mtc(pft_to_mtc(:))
      tau_sap(:) = tau_sap_mtc(pft_to_mtc(:))
      tau_fruit(:) = tau_fruit_mtc(pft_to_mtc(:))
      ecureuil(:) = ecureuil_mtc(pft_to_mtc(:))
      alloc_min(:) = alloc_min_mtc(pft_to_mtc(:))
      alloc_max(:) = alloc_max_mtc(pft_to_mtc(:))
      demi_alloc(:) = demi_alloc_mtc(pft_to_mtc(:))
      leaflife_tab(:) = leaflife_mtc(pft_to_mtc(:))
      !
      ! 3. Senescence
      !
      leaffall(:) = leaffall_mtc(pft_to_mtc(:))
      leafagecrit(:) = leafagecrit_mtc(pft_to_mtc(:))
      senescence_type(:) = senescence_type_mtc(pft_to_mtc(:)) 
      senescence_hum(:) = senescence_hum_mtc(pft_to_mtc(:)) 
      nosenescence_hum(:) = nosenescence_hum_mtc(pft_to_mtc(:)) 
      max_turnover_time(:) = max_turnover_time_mtc(pft_to_mtc(:))
      min_turnover_time(:) = min_turnover_time_mtc(pft_to_mtc(:))
      min_leaf_age_for_senescence(:) = min_leaf_age_for_senescence_mtc(pft_to_mtc(:))
      senescence_temp_c(:) = senescence_temp_c_mtc(pft_to_mtc(:))
      senescence_temp_b(:) = senescence_temp_b_mtc(pft_to_mtc(:))
      senescence_temp_a(:) = senescence_temp_a_mtc(pft_to_mtc(:))
      !
      ! DGVM
      !
      residence_time(:) = residence_time_mtc(pft_to_mtc(:))
      tmin_crit(:) = tmin_crit_mtc(pft_to_mtc(:))
      tcm_crit(:) = tcm_crit_mtc(pft_to_mtc(:))
      !-
   ENDIF !(active_flags%ok_stomate)

 END SUBROUTINE pft_parameters_init
 !
 !=
 !

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_alloc
!!
!>\BRIEF         This subroutine allocates memory needed for the PFT parameters 
!! in function  of the flags activated.  
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

 SUBROUTINE pft_parameters_alloc(active_flags)

   IMPLICIT NONE

   !! 0. Variables and parameters declaration

   !! 0.1 Input variables 
   
   TYPE(control_type),INTENT(in) :: active_flags  !! What parts of the code are activated ? (true/false)

   !! 0.4 Local variables
   
   LOGICAL :: l_error                             !! Diagnostic boolean for error allocation (true/false) 
   INTEGER :: ier                                 !! Return value for memory allocation (0-N, unitless)

!_ ================================================================================================================================


   !
   ! 1. Parameters used anytime
   !

   l_error = .FALSE.

   ALLOCATE(pft_to_mtc(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for pft_to_mtc. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(PFT_name(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for PFT_name. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(height_presc(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for height_presc. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_tree(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_tree. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(natural(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for natural. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_c4(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_c4. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(gsslope(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for gsslope. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(gsoffset(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for gsoffset. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(humcste(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for humcste. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(ext_coeff(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for ext_coeff. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(veget_ori_fixed_test_1(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for veget_ori_fixed_test_1. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(llaimax(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for llaimax. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(llaimin(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for llaimin. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(type_of_lai(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for type_of_lai. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(vcmax_fix(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for vcmax_fix. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(vjmax_fix(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for vjmax_fix. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(co2_tmin_fix(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for co2_tmin_fix. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(co2_topt_fix(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for co2_topt_fix. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(co2_tmax_fix(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for co2_tmax_fix. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(pref_soil_veg(nvm),stat=ier)
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for pref_soil_veg. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_deciduous(nvm),stat=ier) 
   l_error = l_error .OR. (ier /= 0) 
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_deciduous. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_evergreen(nvm),stat=ier) 
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_evergreen. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_c3(nvm),stat=ier) 
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_c3. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_summergreen(nvm),stat=ier)   
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_summergreen. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_needleleaf(nvm),stat=ier)  
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_needleleaf. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_broadleaf(nvm),stat=ier)  
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_broadleaf. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_tropical(nvm),stat=ier)   
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_tropical. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_temperate(nvm),stat=ier)  
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_temperate. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   ALLOCATE(is_boreal(nvm),stat=ier)  
   l_error = l_error .OR. (ier /= 0)
   IF (l_error) THEN
      WRITE(numout,*) ' Memory allocation error for is_boreal. We stop. We need nvm words = ',nvm
      STOP 'pft_parameters_alloc'
   END IF

   !
   ! 2. Parameters used if ok_sechiba only
   !
   IF ( active_flags%ok_sechiba ) THEN

      l_error = .FALSE.

      ALLOCATE(rstruct_const(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for rstruct_const. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(kzero(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for kzero. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(rveg_pft(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for rveg_pft. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(wmax_veg(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for wmax_veg. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      IF ( .NOT.(active_flags%hydrol_cwrr) .OR. (active_flags%hydrol_cwrr .AND. ok_throughfall_by_pft) ) THEN
         ALLOCATE(throughfall_by_pft(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0)
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for throughfall_by_pft. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF
      END IF

      ALLOCATE(snowa_ini(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for snowa_ini. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(snowa_dec(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for snowa_dec. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(alb_leaf_vis(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for alb_leaf_vis. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(alb_leaf_nir(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for alb_leaf_nir. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      IF( active_flags%ok_inca ) THEN
         
         l_error = .FALSE.
         
         ALLOCATE(em_factor_isoprene(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_isoprene. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_monoterpene(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_monoterpene. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_ORVOC(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_ORVOC. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_OVOC(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0)       
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_OVOC. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_MBO(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_MBO. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_methanol(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_methanol. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_acetone(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_acetone. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_acetal(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_acetal. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_formal(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_formal. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_acetic(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0)       
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_acetic. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_formic(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_formic. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_no_wet(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0)
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_no_wet. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(em_factor_no_dry(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0)       
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for em_factor_no_dry. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

         ALLOCATE(Larch(nvm),stat=ier)
         l_error = l_error .OR. (ier /= 0) 
         IF (l_error) THEN
            WRITE(numout,*) ' Memory allocation error for Larch. We stop. We need nvm words = ',nvm
            STOP 'pft_parameters_alloc'
         END IF

      ENDIF ! (active_flags%ok_inca) 

   ENDIF !(active_flags%ok_sechiba)

   !
   ! 3. Parameters used if ok_stomate only
   !
   IF ( active_flags%ok_stomate ) THEN

      l_error = .FALSE.

      ALLOCATE(leaf_tab(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for leaf_tab. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(sla(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for sla. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(vcmax_opt(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for vcmax_opt. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(vjmax_opt(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for vjmax_opt. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_min_a(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_min_a. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_min_b(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_min_b. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_min_c(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_min_c. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_opt_a(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_opt_a. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_opt_b(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_opt_b. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_opt_c(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_opt_c. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_max_a(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_max_a. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_max_b(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_max_b. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tphoto_max_c(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tphoto_max_c. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(pheno_gdd_crit_c(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_c. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(pheno_gdd_crit_b(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_b. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(pheno_gdd_crit_a(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_a. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(pheno_gdd_crit(nvm,3),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit. We stop. We need nvm words = ',nvm*3
         STOP 'pft_parameters_alloc'
      END IF
      pheno_gdd_crit(:,:) = zero

      ALLOCATE(ngd_crit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ngd_crit. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ncdgdd_temp(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ncdgdd_temp. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(hum_frac(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for hum_frac. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(lowgpp_time(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for lowgpp_time. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(hum_min_time(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for hum_min_time. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tau_sap(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tau_sap. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tau_fruit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tau_fruit. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(ecureuil(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for ecureuil. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(alloc_min(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for alloc_min. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(alloc_max(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for alloc_max. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(demi_alloc(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(maint_resp_slope(nvm,3),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for maint_resp_slope. We stop. We need nvm*3 words = ',nvm*3
         STOP 'pft_parameters_alloc'
      END IF
      maint_resp_slope(:,:) = zero

      ALLOCATE(maint_resp_slope_c(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for maint_resp_slope_c. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(maint_resp_slope_b(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for maint_resp_slope_b. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(maint_resp_slope_a(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for maint_resp_slope_a. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(coeff_maint_zero(nvm,nparts),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for coeff_maint_zero. We stop. We need nvm*nparts words = ',nvm*nparts
         STOP 'pft_parameters_alloc'
      END IF
      coeff_maint_zero(:,:) = zero

      ALLOCATE(cm_zero_leaf(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_leaf. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_sapabove(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_sapabove. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_sapbelow(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_sapbelow. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_heartabove(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_heartabove. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_heartbelow(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_heartbelow. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_root(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_root. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_fruit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_fruit. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cm_zero_carbres(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cm_zero_carbres. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(flam(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(resist(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for resist. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(coeff_lcchange_1(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for coeff_lcchange_1. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(coeff_lcchange_10(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for coeff_lcchange_10. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(coeff_lcchange_100(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for coeff_lcchange_100. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(lai_max(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for lai_max. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(pheno_model(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for pheno_model. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(pheno_type(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for pheno_type. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(leaffall(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for leaffall. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(leafagecrit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for leafagecrit. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(senescence_type(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(senescence_hum(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for senescence_hum. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(nosenescence_hum(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for nosenescence_hum. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(max_turnover_time(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for max_turnover_time. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(min_turnover_time(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for min_turnover_time. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(min_leaf_age_for_senescence(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for min_leaf_age_for_senescence. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(senescence_temp_c(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for senescence_temp_c. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(senescence_temp_b(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for senescence_temp_b. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(senescence_temp_a(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for senescence_temp_a. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(senescence_temp(nvm,3),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for senescence_temp. We stop. We need nvm*3 words = ',nvm*3
         STOP 'pft_parameters_alloc'
      END IF
      senescence_temp(:,:) = zero

      ALLOCATE(residence_time(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for residence_time. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tmin_crit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tmin_crit. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tcm_crit(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tcm_crit. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(lai_initmin(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(tree(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for tree. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(bm_sapl(nvm,nparts),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for bm_sapl. We stop. We need nvm*nparts words = ',nvm*nparts
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(migrate(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for migrate. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(maxdia(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for maxdia. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(cn_sapl(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for cn_sapl. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(leaf_timecst(nvm),stat=ier)
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for leaf_timecst. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

      ALLOCATE(leaflife_tab(nvm),stat=ier)   
      l_error = l_error .OR. (ier /= 0)
      IF (l_error) THEN
         WRITE(numout,*) ' Memory allocation error for leaflife_tab. We stop. We need nvm words = ',nvm
         STOP 'pft_parameters_alloc'
      END IF

   ENDIF ! (active_flags%ok_stomate)

 END SUBROUTINE pft_parameters_alloc
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_pft_parameters 
!!
!>\BRIEF          This subroutine will read the imposed values for the global pft
!! parameters (sechiba + stomate). It is not called if IMPOSE_PARAM is set to NO.
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

 SUBROUTINE config_pft_parameters
   
   IMPLICIT NONE

   !! 0. Variables and parameters declaration
  
   !! 0.4 Local variable

   LOGICAL, SAVE ::  first_call = .TRUE.  !! To keep first call trace (true/false)

!_ ================================================================================================================================ 

   IF (first_call) THEN

      !
      ! Vegetation structure
      !

      !Config Key   = SECHIBA_LAI
      !Config Desc  = laimax for maximum lai(see also type of lai interpolation)
      !Config if    = OK_SECHIBA or IMPOSE_VEG
      !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
      !Config Help  = Maximum values of lai used for interpolation of the lai map
      !Config Units = [m^2/m^2]
      CALL getin_p('SECHIBA_LAI',llaimax)

      !Config Key   = LLAIMIN
      !Config Desc  = laimin for minimum lai(see also type of lai interpolation)
      !Config if    = OK_SECHIBA or IMPOSE_VEG
      !Config Def   = 0., 8., 0., 4., 4.5, 0., 4., 0., 0., 0., 0., 0., 0.
      !Config Help  = Minimum values of lai used for interpolation of the lai map
      !Config Units = [m^2/m^2]
      CALL getin_p('LLAIMIN',llaimin)

      !Config Key   = SLOWPROC_HEIGHT
      !Config Desc  = prescribed height of vegetation 
      !Config if    = OK_SECHIBA
      !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1., 1.
      !Config Help  =
      !Config Units = [m] 
      CALL getin_p('SLOWPROC_HEIGHT',height_presc)

      !Config Key   = TYPE_OF_LAI
      !Config Desc  = Type of behaviour of the LAI evolution algorithm 
      !Config if    = OK_SECHIBA
      !Config Def   = inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TYPE_OF_LAI',type_of_lai)

      !Config Key   = IS_TREE
      !Config Desc  = Is the vegetation type a tree ?
      !Config if    = OK_SECHIBA
      !Config Def   = n, y, y, y, y, y, y, y, y, n, n, n, n
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('IS_TREE',is_tree)

      !Config Key   = NATURAL
      !Config Desc  = natural? 
      !Config if    = OK_SECHIBA, OK_STOMATE
      !Config Def   = y, y, y, y, y, y, y, y, y, y, y, n, n 
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('NATURAL',natural)

      !Config Key   = IS_DECIDUOUS
      !Config Desc  = is PFT deciduous ?
      !Config if    = OK_STOMATE
      !Config Def   = n, n, y, n, n, y, n, y, y, n, n, n, n
      !Config Help  =
      !Config Units = [BOOLEAN] 
      CALL getin_p('IS_DECIDUOUS',is_deciduous)

      !Config Key   = IS_EVERGREEN
      !Config Desc  = is PFT evergreen ?
      !Config if    = OK_STOMATE
      !Config Def   = n, y, n, y, y, n, y, n, n, n, n, n, n
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('IS_EVERGREEN',is_evergreen)

      !Config Key   = IS_C3
      !Config Desc  = is PFT C3 ?
      !Config if    = OK_SECHIBA, OK_STOMATE
      !Config Def   = n, n, n, n, n, n, n, n, n, y, n, y, n
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('IS_C3',is_c3)   

!!$      !Config Key   = IS_SUMMERGREEN
!!$      !Config Desc  = Is PFT summergreen ?
!!$      !Config if    = OK_SECHIBA
!!$      !Config Def   = n, n, n, n, n, y, n, y, y, n, n, n, n 
!!$      !Config Help  =
!!$      !Config Units = [BOOLEAN]
!!$      CALL getin_p('IS_SUMMERGREEN',is_summergreen)
!!$
!!$      !Config Key   = IS_NEEDLELEAF
!!$      !Config Desc  = Is PFT needleleaf ?
!!$      !Config if    = OK_SECHIBA
!!$      !Config Def   = n, n, n, y, n, n, y, n, y, n, n, n, n 
!!$      !Config Help  =
!!$      !Config Units = [BOOLEAN]
!!$      CALL getin_p('IS_NEEDLELEAF',is_needleleaf)
!!$
!!$      !Config Key   = IS_BROADLEAF
!!$      !Config Desc  = Is PFT broadleaf ?
!!$      !Config if    = OK_SECHIBA
!!$      !Config Def   = n, y, y, n, y, y, n, y, n, n, n, n, n
!!$      !Config Help  =
!!$      !Config Units = [BOOLEAN]
!!$      CALL getin_p('IS_BROADLEAF',is_broadleaf)
      
      !Config Key   = IS_TROPICAL
      !Config Desc  = Is PFT tropical ?
      !Config if    = OK_SECHIBA
      !Config Def   = n, y, y, n, n, n, n, n, n, n, n, n, n
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('IS_TROPICAL',is_tropical)

!!$      !Config Key   = IS_TEMPERATE 
!!$      !Config Desc  = Is PFT temperate ?
!!$      !Config if    = OK_SECHIBA
!!$      !Config Def   = n, n, n, y, y, y, n, n, n, n, n, n, n
!!$      !Config Help  =
!!$      !Config Units = [BOOLEAN]
!!$      CALL getin_p('IS_TEMPERATE',is_temperate)
!!$
!!$      !Config Key   = IS_BOREAL 
!!$      !Config Desc  = Is PFT boreal ?
!!$      !Config if    = OK_SECHIBA
!!$      !Config Def   = n, n, n, n, n, n, y, y, y, n, n, n, n
!!$      !Config Help  =
!!$      !Config Units = [BOOLEAN]
!!$      CALL getin_p('IS_BOREAL',is_boreal)
      
      !
      ! Photosynthesis
      !

      !Config Key   = IS_C4
      !Config Desc  = flag for C4 vegetation types
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = n, n, n, n, n, n, n, n, n, n, n, y, n, y
      !Config Help  =
      !Config Units = [BOOLEAN]
      CALL getin_p('IS_C4',is_c4)

      !Config Key   = GSSLOPE
      !Config Desc  = Slope of the gs/A relation (Ball & al.)
      !Config if    = OK_CO2
      !Config Def   = 0., 9., 9., 9., 9., 9., 9., 9., 9., 9., 3., 9., 3.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('GSSLOPE',gsslope)

      !Config Key   = GSOFFSET
      !Config Desc  = intercept of the gs/A relation (Ball & al.)
      !Config if    = OK_CO2 or OK_STOMATE
      !Config Def   = 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.01, 0.03
      !Config Help  =
      !Config Units = [-] 
      CALL getin_p('GSOFFSET',gsoffset)

      !Config Key   = VCMAX_FIX
      !Config Desc  = values used for vcmax when STOMATE is not activated
      !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
      !Config Def   = 0., 40., 50., 30., 35., 40.,30., 40., 35., 60., 60., 70., 70.
      !Config Help  =
      !Config Units = [micromol/m^2/s] 
      CALL getin_p('VCMAX_FIX',vcmax_fix)

      !Config Key   = VJMAX_FIX
      !Config Desc  = values used for vjmax when STOMATE is not activated
      !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
      !Config Def   = 0., 80., 100., 60., 70., 80.,  60., 80., 70., 120., 120., 140., 140.
      !Config Help  =
      !Config Units = [micromol/m^2/s]
      CALL getin_p('VJMAX_FIX',vjmax_fix)

      !Config Key   = CO2_TMIN_FIX
      !Config Desc  = values used for photosynthesis tmin when STOMATE is not activated
      !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
      !Config Def   = 0.,  2.,  2., -4., -3., -2., -4., -4., -4., -5.,  6., -5.,  6.
      !Config Help  =
      !Config Units = [C] 
      CALL getin_p('CO2_TMIN_FIX',co2_tmin_fix)

      !Config Key   = CO2_TOPT_FIX 
      !Config Desc  = values used for photosynthesis topt when STOMATE is not activated 
      !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
      !Config Def   = 0., 27.5, 27.5, 17.5, 25., 20.,17.5, 17.5, 17.5, 20.,  32.5, 20.,  32.5
      !Config Help  =
      !Config Units = [C]
      CALL getin_p('CO2_TOPT_FIX',co2_topt_fix)

      !Config Key   = CO2_TMAX_FIX
      !Config Desc  = values used for photosynthesis tmax when STOMATE is not activated 
      !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
      !Config Def   = 0., 55., 55., 38., 48., 38.,38., 38., 38., 45., 55., 45., 55.
      !Config Help  =
      !Config Units = [C]
      CALL getin_p('CO2_TMAX_FIX',co2_tmax_fix)

      !Config Key   = EXT_COEFF
      !Config Desc  = extinction coefficient of the Monsi&Seaki relationship (1953)
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('EXT_COEFF',ext_coeff)
     
      !
      ! Water-hydrology - sechiba
      !

      !Config Key   = HYDROL_HUMCSTE
      !Config Desc  = Root profile
      !Config Def   = 5., .4, .4, 1., .8, .8, 1., 1., .8, 4., 1., 4., 1. 
      !Config if    = OK_SECHIBA
      !Config Help  = Default values were defined for 4 meters soil depth.
      !Config         For 2 meters soil depth, you may use those ones :
      !Config         5., .8, .8, 1., .8, .8, 1., 1., .8, 4., 4., 4., 4.
      !Config Units = [m]
      CALL getin_p('HYDROL_HUMCSTE',humcste)

      !
      ! Soil - vegetation
      !

      !Config Key   = PREF_SOIL_VEG
      !Config Desc  = The soil tile number for each vegetation
      !Config if    = OK_SECHIBA or OK_STOMATE
      !Config Def   = 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3
      !Config Help  = Gives the number of the soil tile on which we will
      !Config         put each vegetation. This allows to divide the hydrological column
      !Config Units = [-]        
      CALL getin_p('PREF_SOIL_VEG',pref_soil_veg)
      
      first_call = .FALSE.

   ENDIF !(first_call)

 END SUBROUTINE config_pft_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_sechiba_pft_parameters
!!
!>\BRIEF        This subroutine will read the imposed values for the sechiba pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO. 
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

 SUBROUTINE config_sechiba_pft_parameters(active_flags)

   IMPLICIT NONE
  
   !! 0. Variables and parameters declaration

   !! 0.1 Input variables

   TYPE(control_type), INTENT(in) :: active_flags     !! What parts of the code are activated ?

   !! 0.4 Local variable

   LOGICAL, SAVE ::  first_call = .TRUE.   !! To keep first call trace (true/false)

!_ ================================================================================================================================ 

   IF (first_call) THEN

      !
      ! Evapotranspiration -  sechiba
      !
      
      !Config Key   = RSTRUCT_CONST
      !Config Desc  = Structural resistance 
      !Config if    = OK_SECHIBA
      !Config Def   = 0.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0,  2.5,  2.0,  2.0,  2.0
      !Config Help  =
      !Config Units = [s/m]
      CALL getin_p('RSTRUCT_CONST',rstruct_const)
      
      !Config Key   = KZERO
      !Config Desc  = A vegetation dependent constant used in the calculation of the surface resistance.
      !Config if    = OK_SECHIBA
      !Config Def   = 0.0, 12.E-5, 12.E-5, 12.e-5, 12.e-5, 25.e-5, 12.e-5,25.e-5, 25.e-5, 30.e-5, 30.e-5, 30.e-5, 30.e-5 
      !Config Help  =
      !Config Units = [kg/m^2/s]
      CALL getin_p('KZERO',kzero)
      
      !Config Key   = RVEG_PFT
      !Config Desc  = Artificial parameter to increase or decrease canopy resistance.
      !Config if    = OK_SECHIBA
      !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
      !Config Help  = This parameter is set by PFT.
      !Config Units = [-]
      CALL getin_p('RVEG_PFT',rveg_pft)    
      
      !
      ! Water-hydrology - sechiba
      !

      !Config Key   = WMAX_VEG
      !Config Desc  = Maximum field capacity for each of the vegetations (Temporary): max quantity of water
      !Config if    = OK_SECHIBA
      !Config Def   = 150., 150., 150., 150., 150., 150., 150.,150., 150., 150., 150., 150., 150.
      !Config Help  =
      !Config Units = [kg/m^3]
      CALL getin_p('WMAX_VEG',wmax_veg)
      !
      IF ( .NOT.(active_flags%hydrol_cwrr) .OR. (active_flags%hydrol_cwrr .AND. ok_throughfall_by_pft) ) THEN
         !Config Key   = PERCENT_THROUGHFALL_PFT
         !Config Desc  = Percent by PFT of precip that is not intercepted by the canopy
         !Config if    = OK_SECHIBA OR HYDROL_CWRR
         !Config Def   = 30. 30. 30. 30. 30. 30. 30. 30. 30. 30. 30. 30. 30.
         !Config Help  = During one rainfall event, PERCENT_THROUGHFALL_PFT% of the incident rainfall
         !Config         will get directly to the ground without being intercepted, for each PFT.
         !Config Units = [%]
         CALL getin_p('PERCENT_THROUGHFALL_PFT',throughfall_by_pft)
         throughfall_by_pft(:) = throughfall_by_pft(:) / 100. 
      END IF
      
      !
      ! Albedo - sechiba
      !

      !Config Key   = SNOWA_INI
      !Config Desc  = Initial snow albedo value for each vegetation type as it will be used in condveg_snow
      !Config if    = OK_SECHIBA
      !Config Def   = 0.35, 0., 0., 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.18, 0.18, 0.18, 0.18
      !Config Help  = Values are from the Thesis of S. Chalita (1992)
      !Config Units = [-]
      CALL getin_p('SNOWA_INI',snowa_ini)

      !Config Key   = SNOWA_DEC
      !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow
      !Config if    = OK_SECHIBA
      !Config Def   = 0.45, 0.,  0., 0.06, 0.06, 0.11, 0.06, 0.11, 0.11, 0.52,0.52, 0.52, 0.52
      !Config Help  = Values are from the Thesis of S. Chalita (1992)
      !Config Units = [-]
      CALL getin_p('SNOWA_DEC',snowa_dec)

      !Config Key   = ALB_LEAF_VIS
      !Config Desc  = leaf albedo of vegetation type, visible albedo
      !Config if    = OK_SECHIBA
      !Config Def   = .00, .04, .06, .06, .06,.06, .06, .06, .06, .10, .10, .10, .10
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('ALB_LEAF_VIS',alb_leaf_vis)

      !Config Key   = ALB_LEAF_NIR
      !Config Desc  = leaf albedo of vegetation type, near infrared albedo
      !Config if    = OK_SECHIBA
      !Config Def   = .00, .20, .22, .22, .22,.22, .22, .22, .22, .30, .30, .30, .30 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('ALB_LEAF_NIR',alb_leaf_nir)
      
      IF ( active_flags%ok_inca ) THEN
         !
         ! BVOC
         !

         !Config Key   = ISO_ACTIVITY
         !Config Desc  = Biogenic activity for each age class : isoprene
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0.5, 1.5, 1.5, 0.5
         !Config Help  =
         !Config Units = [-]
         CALL getin_p('ISO_ACTIVITY',iso_activity)

         !Config Key   = METHANOL_ACTIVITY
         !Config Desc  = Isoprene emission factor for each age class : methanol
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 1., 1., 0.5, 0.5
         !Config Help  =
         !Config Units = [-]
         CALL getin_p('METHANOL_ACTIVITY',methanol_activity)

         !Config Key   = EM_FACTOR_ISOPRENE
         !Config Desc  = Isoprene emission factor
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0., 24., 24., 8., 16., 45., 8., 8., 8., 16., 24., 5., 5.
         !Config Help  =
         !Config Units = [ugC/g/h] 
         CALL getin_p('EM_FACTOR_ISOPRENE',em_factor_isoprene)

         !Config Key   = EM_FACTOR_MONOTERPENE
         !Config Desc  = Monoterpene emission factor 
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0.8, 0.8, 2.4, 1.2, 0.8, 2.4, 2.4, 2.4, 0.8, 1.2, 0.2, 0.2
         !Config Help  =
         !Config Units = [ugC/g/h] 
         CALL getin_p('EM_FACTOR_MONOTERPENE',em_factor_monoterpene)

         !Config Key   = EM_FACTOR_ORVOC
         !Config Desc  = ORVOC emissions factor 
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
         !Config Help  =
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_ORVOC',em_factor_ORVOC)

         !Config Key   = EM_FACTOR_OVOC
         !Config Desc  = OVOC emissions factor
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
         !Config Help  =
         !Config Units = [ugC/g/h]        
         CALL getin_p('EM_FACTOR_OVOC',em_factor_OVOC)

         !Config Key   = EM_FACTOR_MBO
         !Config Desc  = MBO emissions factor 
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0., 0., 20.0, 0., 0., 0., 0., 0., 0., 0., 0., 0.
         !Config Help  =
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_MBO',em_factor_MBO)

         !Config Key   = EM_FACTOR_METHANOL
         !Config Desc  = Methanol emissions factor 
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0.6, 0.6, 1.8, 0.9, 0.6, 1.8, 1.8, 1.8, 0.6, 0.9, 2., 2.
         !Config Help  =
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_METHANOL',em_factor_methanol)

         !Config Key   = EM_FACTOR_ACETONE
         !Config Desc  = Acetone emissions factor
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0.29, 0.29, 0.87, 0.43, 0.29, 0.87, 0.87, 0.87, 0.29, 0.43, 0.07, 0.07 
         !Config Help  =
         !Config Units = [ugC/g/h]     
         CALL getin_p('EM_FACTOR_ACETONE',em_factor_acetone)

         !Config Key   = EM_FACTOR_ACETAL
         !Config Desc  = Acetaldehyde emissions factor 
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0., 0.1, 0.1, 0.3, 0.15, 0.1, 0.3, 0.3, 0.3, 0.1, 0.15, 0.025, 0.025
         !Config Help  =
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_ACETAL',em_factor_acetal)

         !Config Key   = EM_FACTOR_FORMAL
         !Config Desc  = Formaldehyde emissions factor
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0.07, 0.07, 0.2, 0.1, 0.07, 0.2, 0.2, 0.2, 0.07, 0.1, 0.017, 0.017
         !Config Help  = 
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_FORMAL',em_factor_formal)

         !Config Key   = EM_FACTOR_ACETIC
         !Config Desc  = Acetic Acid emissions factor
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0.002, 0.002, 0.006, 0.003, 0.002, 0.006, 0.006, 0.006, 0.002, 0.003, 0.0005, 0.0005
         !Config Help  =
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_ACETIC',em_factor_acetic)

         !Config Key   = EM_FACTOR_FORMIC
         !Config Desc  = Formic Acid emissions factor
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0., 0.01, 0.01, 0.03, 0.015, 0.01, 0.03, 0.03, 0.03, 0.01, 0.015, 0.0025, 0.0025 
         !Config Help  =
         !Config Units = [ugC/g/h]  
         CALL getin_p('EM_FACTOR_FORMIC',em_factor_formic)

         !Config Key   = EM_FACTOR_NO_WET
         !Config Desc  = NOx emissions factor wet soil emissions and exponential dependancy factor 
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0., 2.6, 0.06, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.36, 0.36, 0.36, 0.36
         !Config Help  =
         !Config Units = [ngN/m^2/s]
         CALL getin_p('EM_FACTOR_NO_WET',em_factor_no_wet)

         !Config Key   = EM_FACTOR_NO_DRY
         !Config Desc  = NOx emissions factor dry soil emissions and exponential dependancy factor 
         !Config if    = DIFFUCO_OK_INCA
         !Config Def   = 0., 8.60, 0.40, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 2.65, 2.65, 2.65, 2.65
         !Config Help  =
         !Config Units = [ngN/m^2/s] 
         CALL getin_p('EM_FACTOR_NO_DRY',em_factor_no_dry)

         !Config Key   = LARCH
         !Config Desc  = Larcher 1991 SAI/LAI ratio
         !Config if    = DIFFUCO_OK_INCA 
         !Config Def   = 0., 0.015, 0.015, 0.003, 0.005, 0.005, 0.003, 0.005, 0.003, 0.005, 0.005, 0.008, 0.008
         !Config Help  =
         !Config Units = [-]  
         CALL getin_p('LARCH',Larch)
         
      ENDIF ! (active_flags%ok_inca)

      first_call = .FALSE.

   ENDIF !(first_call)

 END SUBROUTINE config_sechiba_pft_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_pft_parameters 
!!
!>\BRIEF         This subroutine will read the imposed values for the stomate pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO.
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

 SUBROUTINE config_stomate_pft_parameters

   IMPLICIT NONE
   
   !! 0. Variables and parameters declaration

   !! 0.4 Local variable

   LOGICAL, SAVE ::  first_call = .TRUE.   !! To keep first call trace (true/false)

!_ ================================================================================================================================

   IF (first_call) THEN
      
      !
      ! Vegetation structure
      !

      !Config Key   = LEAF_TAB
      !Config Desc  = leaf type : 1=broad leaved tree, 2=needle leaved tree, 3=grass 4=bare ground
      !Config if    = OK_STOMATE
      !Config Def   = 4, 1, 1, 2, 1, 1, 2, 1, 2, 3, 3, 3, 3 
      !Config Help  = 
      !Config Units = [-] 
      CALL getin_p('LEAF_TAB',leaf_tab)

      !Config Key   = SLA
      !Config Desc  = specif leaf area 
      !Config if    = OK_STOMATE
      !Config Def   = 1.5E-2, 1.53E-2, 2.6E-2, 9.26E-3, 2E-2, 2.6E-2, 9.26E-3, 2.6E-2, 1.9E-2, 2.6E-2, 2.6E-2, 2.6E-2, 2.6E-2
      !Config Help  =
      !Config Units = [m^2/gC]
      CALL getin_p('SLA',sla)

      !
      ! Photosynthesis
      !

      !Config Key   = VCMAX_OPT
      !Config Desc  = Maximum rate of carboxylation
      !Config if    = OK_STOMATE
      !Config Def   = undef, 65., 65., 35., 45., 55., 35., 45., 35., 70., 70., 70., 70.
      !Config Help  =
      !Config Units = [micromol/m^2/s]
      CALL getin_p('VCMAX_OPT',vcmax_opt)

      !Config Key   = VJMAX_OPT
      !Config Desc  = Maximum rate of RUbp regeneration
      !Config if    = OK_STOMATE
      !Config Def   = undef, 130., 130., 70., 80., 110., 70., 90., 70., 160., 160., 200., 200.
      !Config Help  =
      !Config Units = [micromol/m^2/s]
      CALL getin_p('VJMAX_OPT',vjmax_opt)

      !Config Key   = TPHOTO_MIN_A
      !Config Desc  = minimum photosynthesis temperature, constant a of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef,  0., 0., 0., 0., 0., 0.,  0., 0.,  0.0025, 0., 0., 0. 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_MIN_A',tphoto_min_a)

      !Config Key   = TPHOTO_MIN_B
      !Config Desc  = minimum photosynthesis temperature, constant b of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef,  0.,  0., 0., 0., 0., 0., 0., 0., 0.1, 0.,0.,0.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_MIN_B',tphoto_min_b)

      !Config Key   = TPHOTO_MIN_C
      !Config Desc  = minimum photosynthesis temperature, constant c of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef,  2., 2., -4., -3.,-2.,-4., -4., -4., -3.25, 13.,-5.,13.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_MIN_C',tphoto_min_c)

      !Config Key   = TPHOTO_OPT_A
      !Config Desc  = optimum photosynthesis temperature, constant a of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0.0025, 0., 0., 0. 
      !Config Help  = 
      !Config Units = [-]
      CALL getin_p('TPHOTO_OPT_A',tphoto_opt_a)

      !Config Key   = TPHOTO_OPT_B
      !Config Desc  = optimum photosynthesis temperature, constant b of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0., 0., 0.  
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_OPT_B',tphoto_opt_b)

      !Config Key   = TPHOTO_OPT_C
      !Config Desc  = optimum photosynthesis temperature, constant c of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 37., 37., 25., 32., 26., 25., 25., 25., 27.25, 36., 30., 36.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_OPT_C',tphoto_opt_c)

      !Config Key   = TPHOTO_MAX_A
      !Config Desc  = maximum photosynthesis temperature, constant a of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef,  0., 0., 0., 0., 0., 0., 0., 0., 0.00375, 0., 0., 0.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_MAX_A',tphoto_max_a)

      !Config Key   = TPHOTO_MAX_B
      !Config Desc  = maximum photosynthesis temperature, constant b of ax^2+bx+c (deg C), tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0.,0.35, 0., 0., 0.   
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_MAX_B',tphoto_max_b)

      !Config Key   = TPHOTO_MAX_C
      !Config Desc  = maximum photosynthesis temperature, constant c of ax^2+bx+c (deg C), tabulated 
      !Config if    = OK_STOMATE
      !Config Def   = undef, 55., 55.,38., 48.,38.,38., 38., 38., 41.125, 55., 45., 55.  
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('TPHOTO_MAX_C',tphoto_max_c)

      !
      ! Respiration - stomate
      !

      !Config Key   = MAINT_RESP_SLOPE_C
      !Config Desc  = slope of maintenance respiration coefficient (1/K), constant c of aT^2+bT+c , tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, .12, .12, .16, .16, .16, .16, .16, .16, .16, .12, .16, .12 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('MAINT_RESP_SLOPE_C',maint_resp_slope_c) 

      !Config Key   = MAINT_RESP_SLOPE_B
      !Config Desc  = slope of maintenance respiration coefficient (1/K), constant b of aT^2+bT+c , tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, -.00133, .0, -.00133, .0 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('MAINT_RESP_SLOPE_B',maint_resp_slope_b)

      !Config Key   = MAINT_RESP_SLOPE_A
      !Config Desc  = slope of maintenance respiration coefficient (1/K), constant a of aT^2+bT+c , tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0    
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('MAINT_RESP_SLOPE_A',maint_resp_slope_a)

      !Config Key   = CM_ZERO_LEAF
      !Config Desc  = maintenance respiration coefficient at 0 deg C, for leaves, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 2.35E-3, 2.62E-3, 1.01E-3, 2.35E-3, 2.62E-3, 1.01E-3,2.62E-3, 2.05E-3, 2.62E-3, 2.62E-3, 2.62E-3, 2.62E-3
      !Config Help  =
      !Config Units = [g/g/day]
      CALL getin_p('CM_ZERO_LEAF',cm_zero_leaf)

      !Config Key   = CM_ZERO_SAPABOVE
      !Config Desc  = maintenance respiration coefficient at 0 deg C,for sapwood above, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4
      !Config Help  =
      !Config Units = [g/g/day]
      CALL getin_p('CM_ZERO_SAPABOVE',cm_zero_sapabove)

      !Config Key   = CM_ZERO_SAPBELOW
      !Config Desc  = maintenance respiration coefficient at 0 deg C, for sapwood below, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4 
      !Config Help  =
      !Config Units = [g/g/day]
      CALL getin_p('CM_ZERO_SAPBELOW',cm_zero_sapbelow)

      !Config Key   = CM_ZERO_HEARTABOVE
      !Config Desc  = maintenance respiration coefficient at 0 deg C, for heartwood above, tabulated
      !Config if    = OK_STOMATE 
      !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 
      !Config Help  =
      !Config Units = [g/g/day]
      CALL getin_p('CM_ZERO_HEARTABOVE',cm_zero_heartabove)

      !Config Key   = CM_ZERO_HEARTBELOW
      !Config Desc  = maintenance respiration coefficient at 0 deg C,for heartwood below, tabulated
      !Config if    = OK_STOMATE 
      !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 
      !Config Help  =
      !Config Units = [g/g/day] 
      CALL getin_p('CM_ZERO_HEARTBELOW',cm_zero_heartbelow)

      !Config Key   = CM_ZERO_ROOT
      !Config Desc  = maintenance respiration coefficient at 0 deg C, for roots, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef,1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3,1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3
      !Config Help  =
      !Config Units = [g/g/day] 
      CALL getin_p('CM_ZERO_ROOT',cm_zero_root)

      !Config Key   = CM_ZERO_FRUIT
      !Config Desc  = maintenance respiration coefficient at 0 deg C, for fruits, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4,1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4    
      !Config Help  =
      !Config Units = [g/g/day] 
      CALL getin_p('CM_ZERO_FRUIT',cm_zero_fruit)

      !Config Key   = CM_ZERO_CARBRES
      !Config Desc  = maintenance respiration coefficient at 0 deg C, for carbohydrate reserve, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4,1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4
      !Config Help  =
      !Config Units = [g/g/day] 
      CALL getin_p('CM_ZERO_CARBRES',cm_zero_carbres)
      
      !
      ! Fire - stomate
      !

      !Config Key   = FLAM
      !Config Desc  = flamability: critical fraction of water holding capacity
      !Config if    = OK_STOMATE
      !Config Def   = undef, .15, .25, .25, .25, .25, .25, .25, .25, .25, .25, .35, .35
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('FLAM',flam)

      !Config Key   = RESIST
      !Config Desc  = fire resistance
      !Config if    = OK_STOMATE
      !Config Def   = undef, .95, .90, .12, .50, .12, .12, .12, .12, .0, .0, .0, .0 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('RESIST',resist)
     
      !
      ! Flux - LUC
      !

      !Config Key   = COEFF_LCCHANGE_1
      !Config Desc  = Coeff of biomass export for the year
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('COEFF_LCCHANGE_1',coeff_lcchange_1)

      !Config Key   = COEFF_LCCHANGE_10
      !Config Desc  = Coeff of biomass export for the decade
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0.403, 0.403, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.403, 0.299, 0.403
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('COEFF_LCCHANGE_10',coeff_lcchange_10)

      !Config Key   = COEFF_LCCHANGE_100
      !Config Desc  = Coeff of biomass export for the century
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0., 0., 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0., 0.104, 0.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('COEFF_LCCHANGE_100',coeff_lcchange_100)
      
      !
      ! Phenology
      !

      !Config Key   = LAI_MAX
      !Config Desc  = maximum LAI, PFT-specific
      !Config if    = OK_STOMATE
      !Config Def   = undef, 7., 7., 5., 5., 5., 4.5, 4.5, 3.0, 2.5, 2.5, 5.,5. 
      !Config Help  =
      !Config Units = [m^2/m^2]
      CALL getin_p('LAI_MAX',lai_max)

      !Config Key   = PHENO_MODEL
      !Config Desc  = which phenology model is used? (tabulated) 
      !Config if    = OK_STOMATE
      !Config Def   = none, none, moi, none, none, ncdgdd, none, ncdgdd, ngd, moigdd, moigdd, moigdd, moigdd
      !Config Help  =
      !Config Units = [-] 
      CALL getin_p('PHENO_MODEL',pheno_model)

      !Config Key   = PHENO_TYPE
      !Config Desc  = type of phenology, 0=bare ground 1=evergreen,  2=summergreen,  3=raingreen,  4=perennial
      !Config if    = OK_STOMATE
      !Config Def   = 0, 1, 3, 1, 1, 2, 1, 2, 2, 4, 4, 2, 3
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('PHENO_TYPE',pheno_type)

      !
      ! Phenology : Leaf Onset
      !

      !Config Key   = PHENO_GDD_CRIT_C
      !Config Desc  = critical gdd, tabulated (C), constant c of aT^2+bT+c
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 270., 400., 125., 400.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('PHENO_GDD_CRIT_C',pheno_gdd_crit_c)

      !Config Key   = PHENO_GDD_CRIT_B
      !Config Desc  = critical gdd, tabulated (C), constant b of aT^2+bT+c
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef,undef, undef, 6.25, 0., 0., 0.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('PHENO_GDD_CRIT_B',pheno_gdd_crit_b)

      !Config Key   = PHENO_GDD_CRIT_A
      !Config Desc  = critical gdd, tabulated (C), constant a of aT^2+bT+c
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.03125,  0., 0., 0.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('PHENO_GDD_CRIT_A',pheno_gdd_crit_a)

      !Config Key   = NGD_CRIT
      !Config Desc  = critical ngd, tabulated. Threshold -5 degrees
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef, 0., undef, undef, undef, undef, undef
      !Config Help  = NGD : Number of Growing Days.
      !Config Units = [days]
      CALL getin_p('NGD_CRIT',ngd_crit)

      !Config Key   = NCDGDD_TEMP
      !Config Desc  = critical temperature for the ncd vs. gdd function in phenology
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, 5., undef, 0., undef, undef, undef, undef, undef
      !Config Help  =
      !Config Units = [C] 
      CALL getin_p('NCDGDD_TEMP',ncdgdd_temp)

      !Config Key   = HUM_FRAC
      !Config Desc  = critical humidity (relative to min/max) for phenology
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, .5, undef, undef, undef, undef, undef,  undef, .5, .5, .5,.5     
      !Config Help  =
      !Config Units = [%]
      CALL getin_p('HUM_FRAC',hum_frac)

      !Config Key   = LOWGPP_TIME
      !Config Desc  = minimum duration of dormance for phenology 
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, 30., undef, undef, 30., undef, 30., 30., 30., 30., 30., 30.  
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('LOWGPP_TIME',lowgpp_time)

      !Config Key   = HUM_MIN_TIME
      !Config Desc  = minimum time elapsed since moisture minimum
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, 50., undef, undef, undef, undef, undef, undef, 35., 35., 75., 75.
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('HUM_MIN_TIME',hum_min_time)

      !Config Key   = TAU_SAP
      !Config Desc  = sapwood -> heartwood conversion time
      !Config if    = OK_STOMATE
      !Config Def   = undef, 730., 730., 730., 730., 730., 730., 730., 730., undef, undef, undef, undef
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('TAU_SAP',tau_sap)

      !Config Key   = TAU_FRUIT
      !Config Desc  = fruit lifetime
      !Config if    = OK_STOMATE
      !Config Def   = undef, 90., 90., 90., 90., 90., 90., 90., 90., undef, undef, undef, undef
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('TAU_FRUIT',tau_fruit)

      !Config Key   = ECUREUIL
      !Config Desc  = fraction of primary leaf and root allocation put into reserve
      !Config if    = OK_STOMATE
      !Config Def   = undef, .0, 1., .0, .0, 1., .0, 1., 1., 1., 1., 1., 1.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('ECUREUIL',ecureuil)

      !Config Key   = ALLOC_MIN
      !Config Desc  = minimum allocation above/below = f(age) - 30/01/04 NV/JO/PF
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, undef, undef, undef, undef 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('ALLOC_MIN',alloc_min)

      !Config Key   = ALLOC_MAX
      !Config Desc  = maximum allocation above/below = f(age) - 30/01/04 NV/JO/PF
      !Config if    = OK_STOMATE
      !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, undef, undef, undef, undef
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('ALLOC_MAX',alloc_max)

      !Config Key   = DEMI_ALLOC 
      !Config Desc  = mean allocation above/below = f(age) - 30/01/04 NV/JO/PF
      !Config if    = OK_STOMATE
      !Config Def   = undef, 5., 5., 5., 5., 5., 5., 5., 5., undef, undef, undef, undef
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('DEMI_ALLOC',demi_alloc)

      !Config Key   = LEAFLIFE_TAB
      !Config Desc  = leaf longevity
      !Config if    = OK_STOMATE
      !Config Def   = undef, .5, 2., .33, 1., 2., .33, 2., 2., 2., 2., 2., 2. 
      !Config Help  =
      !Config Units = [years]
      CALL getin_p('LEAFLIFE_TAB',leaflife_tab)

      !
      ! Phenology : Senescence
      !
      !
      !Config Key   = LEAFFALL
      !Config Desc  = length of death of leaves, tabulated 
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, 10., undef, undef, 10., undef, 10., 10., 10., 10., 10., 10. 
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('LEAFFALL',leaffall)

      !Config Key   = LEAFAGECRIT
      !Config Desc  = critical leaf age, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, 730., 180., 910., 730., 180., 910., 180., 180., 120., 120., 90., 90.  
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('LEAFAGECRIT',leafagecrit) 

      !Config Key   = SENESCENCE_TYPE
      !Config Desc  = type of senescence, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = none, none, dry, none, none, cold, none, cold, cold, mixed, mixed, mixed, mixed 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('SENESCENCE_TYPE',senescence_type) 

      !Config Key   = SENESCENCE_HUM
      !Config Desc  = critical relative moisture availability for senescence
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, .3, undef, undef, undef, undef, undef, undef, .2, .2, .3, .2 
      !Config Help  =
      !Config Units = [-] 
      CALL getin_p('SENESCENCE_HUM',senescence_hum)

      !Config Key   = NOSENESCENCE_HUM
      !Config Desc  = relative moisture availability above which there is no humidity-related senescence
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, .8, undef, undef, undef, undef, undef, undef, .3, .3, .3, .3 
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('NOSENESCENCE_HUM',nosenescence_hum) 

      !Config Key   = MAX_TURNOVER_TIME
      !Config Desc  = maximum turnover time for grasse
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef,  80.,  80., 80., 80. 
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('MAX_TURNOVER_TIME',max_turnover_time)

      !Config Key   = MIN_TURNOVER_TIME
      !Config Desc  = minimum turnover time for grasse 
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 10., 10., 10., 10. 
      !Config Help  =
      !Config Units = [days]
      CALL getin_p('MIN_TURNOVER_TIME',min_turnover_time)

      !Config Key   = MIN_LEAF_AGE_FOR_SENESCENCE
      !Config Desc  = minimum leaf age to allow senescence g
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, 90., undef, undef, 90., undef, 60., 60., 30., 30., 30., 30.
      !Config Help  =
      !Config Units = [days] 
      CALL getin_p('MIN_LEAF_AGE_FOR_SENESCENCE',min_leaf_age_for_senescence)

      !Config Key   = SENESCENCE_TEMP_C
      !Config Desc  = critical temperature for senescence (C), constant c of aT^2+bT+c, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, 12., undef, 7., 2., -1.375, 5., 5., 10.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('SENESCENCE_TEMP_C',senescence_temp_c)

      !Config Key   = SENESCENCE_TEMP_B
      !Config Desc  = critical temperature for senescence (C), constant b of aT^2+bT+c ,tabulated
      !Config if    = OK_STOMATE 
      !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0., .1, 0., 0., 0.
      !Config Help  =
      !Config Units = [-]
      CALL getin_p('SENESCENCE_TEMP_B',senescence_temp_b)

      !Config Key   = SENESCENCE_TEMP_A
      !Config Desc  = critical temperature for senescence (C), constant a of aT^2+bT+c , tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0.,.00375, 0., 0., 0. 
      !Config Help  =
      !Config Units = [-] 
      CALL getin_p('SENESCENCE_TEMP_A',senescence_temp_a)
      
      !
      ! DGVM
      !

      !Config Key   = RESIDENCE_TIME
      !Config Desc  = residence time of trees
      !Config if    = OK_DGVM and NOT(LPJ_GAP_CONST_MORT)
      !Config Def   = undef, 30.0, 30.0, 40.0, 40.0, 40.0, 80.0, 80.0, 80.0, 0.0, 0.0, 0.0, 0.0 
      !Config Help  =
      !Config Units = [years]
      CALL getin_p('RESIDENCE_TIME',residence_time)

      !Config Key   = TMIN_CRIT
      !Config Desc  = critical tmin, tabulated
      !Config if    = OK_STOMATE
      !Config Def   = undef,  0.0, 0.0, -30.0, -14.0, -30.0, -45.0, -45.0, undef, undef, undef, undef, undef
      !Config Help  = 
      !Config Units = [C]
      CALL getin_p('TMIN_CRIT',tmin_crit)

      !Config Key   = TCM_CRIT
      !Config Desc  = critical tcm, tabulated 
      !Config if    = OK_STOMATE
      !Config Def   = undef, undef, undef, 5.0, 15.5, 15.5, -8.0, -8.0, -8.0, undef, undef, undef, undef
      !Config Help  =
      !Config Units = [C]
      CALL getin_p('TCM_CRIT',tcm_crit)
      
      first_call = .FALSE.
       
   ENDIF !(first_call)
  
 END SUBROUTINE config_stomate_pft_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_clear
!!
!>\BRIEF         This subroutine deallocates memory at the end of the simulation. 
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

 SUBROUTINE pft_parameters_clear
   
   l_first_pft_parameters = .TRUE.
   
   IF (ALLOCATED(pft_to_mtc)) DEALLOCATE(pft_to_mtc)
   IF (ALLOCATED(PFT_name)) DEALLOCATE(PFT_name)
   IF (ALLOCATED(veget_ori_fixed_test_1)) DEALLOCATE(veget_ori_fixed_test_1)   
   IF (ALLOCATED(llaimax)) DEALLOCATE(llaimax)
   IF (ALLOCATED(llaimin)) DEALLOCATE(llaimin)
   IF (ALLOCATED(height_presc)) DEALLOCATE(height_presc)   
   IF (ALLOCATED(type_of_lai)) DEALLOCATE(type_of_lai)
   IF (ALLOCATED(is_tree)) DEALLOCATE(is_tree)
   IF (ALLOCATED(natural)) DEALLOCATE(natural)
   IF (ALLOCATED(is_deciduous)) DEALLOCATE(is_deciduous)
   IF (ALLOCATED(is_evergreen)) DEALLOCATE(is_evergreen)
   IF (ALLOCATED(is_c3)) DEALLOCATE(is_c3)  
   IF (ALLOCATED(is_summergreen)) DEALLOCATE(is_summergreen)
   IF (ALLOCATED(is_needleleaf)) DEALLOCATE(is_needleleaf)
   IF (ALLOCATED(is_broadleaf)) DEALLOCATE(is_broadleaf)
   IF (ALLOCATED(is_tropical)) DEALLOCATE(is_tropical)
   IF (ALLOCATED(is_temperate)) DEALLOCATE(is_temperate)
   IF (ALLOCATED(is_boreal)) DEALLOCATE(is_boreal)
   IF (ALLOCATED(humcste)) DEALLOCATE(humcste)
   IF (ALLOCATED(pref_soil_veg)) DEALLOCATE(pref_soil_veg)
   IF (ALLOCATED(is_c4)) DEALLOCATE(is_c4)  
   IF (ALLOCATED(gsslope)) DEALLOCATE(gsslope)
   IF (ALLOCATED(gsoffset)) DEALLOCATE(gsoffset)
   IF (ALLOCATED(vcmax_fix)) DEALLOCATE(vcmax_fix)
   IF (ALLOCATED(vjmax_fix)) DEALLOCATE(vjmax_fix)
   IF (ALLOCATED(co2_tmin_fix)) DEALLOCATE(co2_tmin_fix)
   IF (ALLOCATED(co2_topt_fix)) DEALLOCATE(co2_topt_fix)
   IF (ALLOCATED(co2_tmax_fix)) DEALLOCATE(co2_tmax_fix) 
   IF (ALLOCATED(ext_coeff)) DEALLOCATE(ext_coeff)
   IF (ALLOCATED(rveg_pft)) DEALLOCATE(rveg_pft)
   IF (ALLOCATED(rstruct_const)) DEALLOCATE(rstruct_const)
   IF (ALLOCATED(kzero)) DEALLOCATE(kzero)
   IF (ALLOCATED(wmax_veg)) DEALLOCATE(wmax_veg)
   IF (ALLOCATED(throughfall_by_pft)) DEALLOCATE(throughfall_by_pft)
   IF (ALLOCATED(snowa_ini)) DEALLOCATE(snowa_ini)
   IF (ALLOCATED(snowa_dec)) DEALLOCATE(snowa_dec)
   IF (ALLOCATED(alb_leaf_vis)) DEALLOCATE(alb_leaf_vis)
   IF (ALLOCATED(alb_leaf_nir)) DEALLOCATE(alb_leaf_nir)   
   IF (ALLOCATED(em_factor_isoprene)) DEALLOCATE(em_factor_isoprene)
   IF (ALLOCATED(em_factor_monoterpene)) DEALLOCATE(em_factor_monoterpene)
   IF (ALLOCATED(em_factor_ORVOC)) DEALLOCATE(em_factor_ORVOC)
   IF (ALLOCATED(em_factor_OVOC)) DEALLOCATE(em_factor_OVOC)
   IF (ALLOCATED(em_factor_MBO)) DEALLOCATE(em_factor_MBO)
   IF (ALLOCATED(em_factor_methanol)) DEALLOCATE(em_factor_methanol)
   IF (ALLOCATED(em_factor_acetone)) DEALLOCATE(em_factor_acetone)
   IF (ALLOCATED(em_factor_acetal)) DEALLOCATE(em_factor_acetal)
   IF (ALLOCATED(em_factor_formal)) DEALLOCATE(em_factor_formal)
   IF (ALLOCATED(em_factor_acetic)) DEALLOCATE(em_factor_acetic)
   IF (ALLOCATED(em_factor_formic)) DEALLOCATE(em_factor_formic)
   IF (ALLOCATED(em_factor_no_wet)) DEALLOCATE(em_factor_no_wet)
   IF (ALLOCATED(em_factor_no_dry)) DEALLOCATE(em_factor_no_dry)
   IF (ALLOCATED(Larch)) DEALLOCATE(Larch)
   IF (ALLOCATED(leaf_tab)) DEALLOCATE(leaf_tab)
   IF (ALLOCATED(sla)) DEALLOCATE(sla)
   IF (ALLOCATED(vcmax_opt)) DEALLOCATE(vcmax_opt)
   IF (ALLOCATED(vjmax_opt)) DEALLOCATE(vjmax_opt)
   IF (ALLOCATED(tphoto_min_a)) DEALLOCATE(tphoto_min_a)
   IF (ALLOCATED(tphoto_min_b)) DEALLOCATE(tphoto_min_b)
   IF (ALLOCATED(tphoto_min_c)) DEALLOCATE(tphoto_min_c)
   IF (ALLOCATED(tphoto_opt_a)) DEALLOCATE(tphoto_opt_a)
   IF (ALLOCATED(tphoto_opt_b)) DEALLOCATE(tphoto_opt_b)
   IF (ALLOCATED(tphoto_opt_c)) DEALLOCATE(tphoto_opt_c)
   IF (ALLOCATED(tphoto_max_a)) DEALLOCATE(tphoto_max_a)
   IF (ALLOCATED(tphoto_max_b)) DEALLOCATE(tphoto_max_b)
   IF (ALLOCATED(tphoto_max_c)) DEALLOCATE(tphoto_max_c)
   IF (ALLOCATED(maint_resp_slope)) DEALLOCATE(maint_resp_slope)
   IF (ALLOCATED(maint_resp_slope_c)) DEALLOCATE(maint_resp_slope_c)
   IF (ALLOCATED(maint_resp_slope_b)) DEALLOCATE(maint_resp_slope_b)
   IF (ALLOCATED(maint_resp_slope_a)) DEALLOCATE(maint_resp_slope_a)
   IF (ALLOCATED(coeff_maint_zero)) DEALLOCATE(coeff_maint_zero)
   IF (ALLOCATED(cm_zero_leaf)) DEALLOCATE(cm_zero_leaf)
   IF (ALLOCATED(cm_zero_sapabove)) DEALLOCATE(cm_zero_sapabove)
   IF (ALLOCATED(cm_zero_sapbelow)) DEALLOCATE(cm_zero_sapbelow)
   IF (ALLOCATED(cm_zero_heartabove)) DEALLOCATE(cm_zero_heartabove)
   IF (ALLOCATED(cm_zero_heartbelow)) DEALLOCATE(cm_zero_heartbelow)
   IF (ALLOCATED(cm_zero_root)) DEALLOCATE(cm_zero_root)
   IF (ALLOCATED(cm_zero_fruit)) DEALLOCATE(cm_zero_fruit)
   IF (ALLOCATED(cm_zero_carbres)) DEALLOCATE(cm_zero_carbres)
   IF (ALLOCATED(flam)) DEALLOCATE(flam)
   IF (ALLOCATED(resist)) DEALLOCATE(resist)
   IF (ALLOCATED(coeff_lcchange_1)) DEALLOCATE(coeff_lcchange_1)
   IF (ALLOCATED(coeff_lcchange_10)) DEALLOCATE(coeff_lcchange_10)
   IF (ALLOCATED(coeff_lcchange_100)) DEALLOCATE(coeff_lcchange_100)
   IF (ALLOCATED(lai_max)) DEALLOCATE(lai_max)
   IF (ALLOCATED(pheno_model)) DEALLOCATE(pheno_model)
   IF (ALLOCATED(pheno_type)) DEALLOCATE(pheno_type)
   IF (ALLOCATED(pheno_gdd_crit_c)) DEALLOCATE(pheno_gdd_crit_c)
   IF (ALLOCATED(pheno_gdd_crit_b)) DEALLOCATE(pheno_gdd_crit_b)
   IF (ALLOCATED(pheno_gdd_crit_a)) DEALLOCATE(pheno_gdd_crit_a)
   IF (ALLOCATED(pheno_gdd_crit)) DEALLOCATE(pheno_gdd_crit)
   IF (ALLOCATED(ngd_crit)) DEALLOCATE(ngd_crit)
   IF (ALLOCATED(ncdgdd_temp)) DEALLOCATE(ncdgdd_temp)
   IF (ALLOCATED(hum_frac)) DEALLOCATE(hum_frac)
   IF (ALLOCATED(lowgpp_time)) DEALLOCATE(lowgpp_time)   
   IF (ALLOCATED(hum_min_time)) DEALLOCATE(hum_min_time)
   IF (ALLOCATED(tau_sap)) DEALLOCATE(tau_sap)
   IF (ALLOCATED(tau_fruit)) DEALLOCATE(tau_fruit)
   IF (ALLOCATED(ecureuil)) DEALLOCATE(ecureuil)
   IF (ALLOCATED(alloc_min)) DEALLOCATE(alloc_min)
   IF (ALLOCATED(alloc_max)) DEALLOCATE(alloc_max)
   IF (ALLOCATED(demi_alloc)) DEALLOCATE(demi_alloc)
   IF (ALLOCATED(leaflife_tab)) DEALLOCATE(leaflife_tab)
   IF (ALLOCATED(leaffall)) DEALLOCATE(leaffall)
   IF (ALLOCATED(leafagecrit)) DEALLOCATE(leafagecrit)
   IF (ALLOCATED(senescence_type)) DEALLOCATE(senescence_type)
   IF (ALLOCATED(senescence_hum)) DEALLOCATE(senescence_hum)
   IF (ALLOCATED(nosenescence_hum)) DEALLOCATE(nosenescence_hum)
   IF (ALLOCATED(max_turnover_time)) DEALLOCATE(max_turnover_time)
   IF (ALLOCATED(min_turnover_time)) DEALLOCATE(min_turnover_time)
   IF (ALLOCATED(min_leaf_age_for_senescence)) DEALLOCATE(min_leaf_age_for_senescence)
   IF (ALLOCATED(senescence_temp_c)) DEALLOCATE(senescence_temp_c)
   IF (ALLOCATED(senescence_temp_b)) DEALLOCATE(senescence_temp_b)
   IF (ALLOCATED(senescence_temp_a)) DEALLOCATE(senescence_temp_a)
   IF (ALLOCATED(senescence_temp)) DEALLOCATE(senescence_temp)
   IF (ALLOCATED(residence_time)) DEALLOCATE(residence_time)
   IF (ALLOCATED(tmin_crit)) DEALLOCATE(tmin_crit)
   IF (ALLOCATED(tcm_crit)) DEALLOCATE(tcm_crit)
   IF (ALLOCATED(lai_initmin)) DEALLOCATE(lai_initmin)
   IF (ALLOCATED(tree)) DEALLOCATE(tree)
   IF (ALLOCATED(bm_sapl)) DEALLOCATE(bm_sapl)
   IF (ALLOCATED(migrate)) DEALLOCATE(migrate)
   IF (ALLOCATED(maxdia)) DEALLOCATE(maxdia)
   IF (ALLOCATED(cn_sapl)) DEALLOCATE(cn_sapl)
   IF (ALLOCATED(leaf_timecst)) DEALLOCATE(leaf_timecst)
   
 END SUBROUTINE pft_parameters_clear

END MODULE pft_parameters
