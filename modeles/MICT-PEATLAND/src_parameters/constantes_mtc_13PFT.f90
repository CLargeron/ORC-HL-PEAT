! =================================================================================================================================
! MODULE       : constantes_mtc
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         This module contains the standard values of the paramters for the 13 metaclasses of vegetation used by ORCHIDEE.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): Didier Solyga : replace default values for humscte at 2 meters soil depth by default values for humcste
!!                   at 4 meters (used for the CMIP simulations). The standard values for 2 meters soil depth are :
!!                   REAL(r_std), PARAMETER, DIMENSION(nvmc) :: humcste_mtc  =  &
!!                   & (/ 5.0,   0.8,   0.8,   1.0,   0.8,   0.8,   1.0,  &
!!                   &    1.0,   0.8,   4.0,   4.0,   4.0,   4.0  /)
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE constantes_mtc

  USE defprec
  USE constantes

  IMPLICIT NONE

  !
  ! METACLASSES CHARACTERISTICS
  !

  INTEGER(i_std), PARAMETER :: nvmc = 13                         !! Number of MTCS fixed in the code (unitless)

  CHARACTER(len=34), PARAMETER, DIMENSION(nvmc) :: MTC_name = &  !! description of the MTC (unitless)
  & (/ 'bare ground                       ', &          !  1
  &    'tropical  broad-leaved evergreen  ', &          !  2
  &    'tropical  broad-leaved raingreen  ', &          !  3
  &    'temperate needleleaf   evergreen  ', &          !  4
  &    'temperate broad-leaved evergreen  ', &          !  5
  &    'temperate broad-leaved summergreen', &          !  6
  &    'boreal    needleleaf   evergreen  ', &          !  7
  &    'boreal    broad-leaved summergreen', &          !  8
  &    'boreal    needleleaf   summergreen', &          !  9
  &    '          C3           grass      ', &          ! 10
  &    '          C4           grass      ', &          ! 11
  &    '          C3           agriculture', &          ! 12
  &    '          C4           agriculture'  /)         ! 13


  !
  ! VEGETATION STRUCTURE
  !
  !-
  !  1. Sechiba
  !
  !   1.1 Labels - Characteristics
  !-
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_tree_mtc  =  &                          !! Is the vegetation type a tree ? 
  & (/  .FALSE.,   .TRUE.,   .TRUE.,    .TRUE.,    .TRUE.,    .TRUE.,   .TRUE., &   !! (true/false)
  &     .TRUE.,    .TRUE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.  /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_deciduous_mtc  =  &                     !! Is PFT deciduous ? (true/false)
  & (/ .FALSE.,   .FALSE.,   .TRUE. ,   .FALSE.,   .FALSE.,   .TRUE.,   .FALSE.,  &
  &    .TRUE. ,   .TRUE. ,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.  /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_evergreen_mtc  =  &                     !! Is PFT evergreen ? (true/false)
  & (/ .FALSE.,   .TRUE.,    .FALSE.,   .TRUE.,    .TRUE.,    .FALSE.,   .TRUE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)        

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_summergreen_mtc  =  &                    !! Is PFT summergreen ? (true/false)
  & (/ .FALSE.,   .FALSE.,    .FALSE.,   .FALSE.,   .FALSE.,   .TRUE.,   .FALSE.,  &
  &    .TRUE.,    .TRUE.,     .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)        

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_needleleaf_mtc  =  &                     !! Is PFT needleleaf ? (true/false)
  & (/ .FALSE.,   .FALSE.,    .FALSE.,   .TRUE.,    .FALSE.,    .FALSE.,   .TRUE.,  &
  &    .FALSE.,   .TRUE.,     .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_broadleaf_mtc  =  &                      !! Is PFT broadleaf ? (true/false)
  & (/ .FALSE.,   .TRUE.,    .TRUE.,    .FALSE.,   .TRUE.,    .TRUE.,   .FALSE.,  &
  &    .TRUE.,    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_tropical_mtc  =  &                       !! Is PFT tropical ? (true/false)
  & (/ .FALSE.,   .TRUE.,    .TRUE.,    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_temperate_mtc  =  &                      !! Is PFT temperate ? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .TRUE.,    .TRUE.,    .TRUE.,   .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_boreal_mtc  =  &                         !! Is PFT boreal ? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .TRUE.,  &
  &    .TRUE.,    .TRUE.,    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_c3_mtc  =  &                             !! Is PFT C3 ? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  .FALSE.,  & 
  &    .FALSE.,   .FALSE.,   .TRUE.,    .FALSE.,   .TRUE.,    .FALSE.  /)

  CHARACTER(LEN=5), PARAMETER, DIMENSION(nvmc) :: type_of_lai_mtc  =  &  !! Type of behaviour of the LAI evolution algorithm
  & (/ 'inter', 'inter', 'inter', 'inter', 'inter',  &                   !! for each vegetation type. (unitless)
  &    'inter', 'inter', 'inter', 'inter', 'inter',  &                   !! Value of type_of_lai : mean or interp
  &    'inter', 'inter', 'inter' /)

  !-
  !  1.2 Prescribed Values
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: veget_ori_fixed_mtc  =  &  !! Value for veget_ori for tests in
  & (/ 0.2,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                !! 0-dim simulations (0-1, unitless)
  &    0.0,   0.0,   0.8,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: llaimax_mtc  =  &          !! laimax for maximum
  & (/ 0.0,   8.0,   8.0,   4.0,   4.5,   4.5,   4.0,  &                !! See also type of lai interpolation (m^2.m^{-2})
  &    4.5,   4.0,   2.0,   2.0,   2.0,   2.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: llaimin_mtc  = &           !! laimin for minimum lai
  & (/ 0.0,   8.0,   0.0,   4.0,   4.5,   0.0,   4.0,  &                !! See also type of lai interpolation (m^2.m^{-2})
  &    0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: height_presc_mtc  =  &     !! prescribed height of vegetation (m)
  & (/  0.0,   30.0,   30.0,   20.0,   20.0,   20.0,   15.0,  &         !! Value for height_presc : one for each vegetation type
  &    15.0,   15.0,    0.5,    0.6,    1.0,    1.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rveg_mtc  =  &   
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &    1.0,   1.0,   1.0,   1.0,   1.0,   1.0   /)

  !-
  !  2. Stomate
  !
  !   2.1 Labels - Characteristics
  !
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: natural_mtc =  &                         !! natural?  (true/false)
  & (/ .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,    .TRUE.,   .TRUE.,  &
  &    .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,   .FALSE.,   .FALSE.  /)

  INTEGER(i_std),PARAMETER, DIMENSION(nvmc) :: leaf_tab_mtc  =  &                 !! leaf type (1-4, unitless)
  & (/  4,   1,   1,   2,   1,   1,   2,   &                                      !! 1=broad leaved tree, 2=needle leaved tree
  &     1,   2,   3,   3,   3,   3   /)                                           !! 3=grass 4=bare ground
  !-
  !   2.2 Prescribed Values
  !-
   REAL(r_std), PARAMETER, DIMENSION(nvmc) :: sla_mtc  =  &                       !! specif leaf area (m^2.gC^{-1})
  & (/ 1.5E-2,   1.53E-2,   2.6E-2,   9.26E-3,     2E-2,   2.6E-2,   9.26E-3,  &
  &    2.6E-2,    1.9E-2,   2.6E-2,    2.6E-2,   2.6E-2,   2.6E-2  /) 


  !
  ! EVAPOTRANSPIRATION (sechiba)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rstruct_const_mtc  =  &  !! Structural resistance. (s.m^{-1}) 
  & (/ 0.0,   25.0,   25.0,   25.0,   25.0,   25.0,   25.0,  &        !! Value for rstruct_const : one for each vegetation type
  &   25.0,   25.0,    2.5,    2.0,    2.0,    2.0   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: kzero_mtc  =  &                  !! A vegetation dependent constant used in the 
  & (/    0.0,   12.E-5,   12.E-5,   12.E-5,   12.E-5,   25.E-5,   12.E-5,  & !! calculation  of the surface resistance. (kg.m^2.s^{-1})
  &    25.E-5,   25.E-5,   30.E-5,   30.E-5,   30.E-5,   30.E-5  /)           !! Value for kzero one for each vegetation type


  !
  ! WATER (sechiba)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: wmax_veg_mtc  =  &        !! Maximum field capacity for each of the 
  & (/ 150.0,   150.0,   150.0,   150.0,   150.0,   150.0,   150.0,  & !! vegetations types (Temporary). (kg.m^{-3})
  &    150.0,   150.0,   150.0,   150.0,   150.0,   150.0  /)          !! Value of wmax_veg : max quantity of water :
                                                                       !! one for each vegetation type

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: humcste_mtc  =  &         !! Root profile description for the different 
  & (/ 5.0,   0.4,   0.4,   1.0,   0.8,   0.8,   1.0,  &               !! vegetations types. (m^{-1})
  &    1.0,   0.8,   4.0,   1.0,   4.0,   1.0  /)                      !! These are the factor in the exponential which gets       
                                                                       !! the root density as a function of depth
                                                                       !! Values for dpu_max = 4.0  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: humcste_cwrr  =  &        !! Root profile description for the different
  & (/ 5.0,   0.8,   0.8,   1.0,   0.8,   0.8,   1.0,  &               !! vegetations types. (m^{-1})
  &    1.0,   0.8,   4.0,   4.0,   4.0,   4.0  /)                      !! These are the factor in the exponential which gets       
                                                                       !! the root density as a function of depth
                                                                       !! Values for dpu_max = 2.0
                                                                       !! (used by using 11 layers hydrology) 


  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: throughfall_by_mtc  =  &  !! Fraction of rain intercepted by the canopy
  & (/ 30.0,   30.0,   30.0,   30.0,   30.0,   30.0,   30.0,  &        !! (0-100, unitless)
  &    30.0,   30.0,   30.0,   30.0,   30.0,   30.0  /)


  !
  ! ALBEDO (sechiba)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: snowa_ini_mtc  =  &     !! Initial snow albedo value for each vegetation type
  & (/ 0.35,    0.0,    0.0,   0.14,   0.14,   0.14,   0.14,  &      !! as it will be used in condveg_snow (unitless)
  &    0.14,   0.14,   0.18,   0.18,   0.18,   0.18  /)              !! Source : Values are from the Thesis of S. Chalita (1992)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: snowa_dec_mtc  =  &     !! Decay rate of snow albedo value for each vegetation type
  & (/ 0.45,    0.0,    0.0,   0.06,   0.06,   0.11,   0.06,  &      !! as it will be used in condveg_snow (unitless)
  &    0.11,   0.11,   0.52,   0.52,   0.52,   0.52  /)              !! Source : Values are from the Thesis of S. Chalita (1992)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alb_leaf_vis_mtc  =  &  !! leaf albedo of vegetation type, visible albedo 
  & (/ 0.00,   0.04,   0.06,   0.06,   0.06,   0.06,   0.06,  &      !! (unitless)
  &    0.06,   0.06,   0.10,   0.10,   0.10,   0.10  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alb_leaf_nir_mtc  =  &  !! leaf albedo of vegetation type, near infrared albedo
  & (/ 0.00,   0.20,   0.22,   0.22,   0.22,   0.22,   0.22,  &      !! (unitless)
  &    0.22,   0.22,   0.30,   0.30,   0.30,   0.30  /)


  !
  ! SOIL - VEGETATION
  !
  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: pref_soil_veg_mtc  =  &       !! The soil tile number for each vegetation
  & (/ 1,   2,   2,   2,   2,   2,   2,  &                                    
  &    2,   2,   3,   3,   3,   3  /)                                         


  !
  ! PHOTOSYNTHESIS
  !
  !-
  ! 1 .CO2
  !-
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_c4_mtc  =  &                            !! flag for C4 vegetation types (true/false)
  & (/ .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  &
  &    .FALSE.,  .FALSE.,   .FALSE.,   .TRUE.,    .FALSE.,   .TRUE.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gsslope_mtc  =  &       !! Slope of the gs/A relation (Ball & al.) (unitless)
  & (/ 0.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,  & 
  &    9.0,   9.0,   9.0,   3.0,   9.0,   3.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gsoffset_mtc  =  &      !! intercept of the gs/A relation (Ball & al.)(unitless)
  & (/  0.0,   0.01,   0.01,   0.01,   0.01,   0.01,   0.01,  &
  &    0.01,   0.01,   0.01,   0.03,   0.01,   0.03  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vcmax_fix_mtc  =  &     !! values used for vcmax when STOMATE is not
  & (/  0.0,   40.0,   50.0,   30.0,   35.0,   40.0,   30.0,  &      !! activated (�mol.m^{-2}.s^{-1})
  &    40.0,   35.0,   60.0,   60.0,   70.0,   70.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vjmax_fix_mtc  =  &     !! values used for vjmax when STOMATE is no
  & (/  0.0,   80.0,   100.0,    60.0,    70.0,    80.0,   60.0,  &  !! activated (�mol.m^{-2}.s^{-1})
  &    80.0,   70.0,   120.0,   120.0,   140.0,   140.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: co2_tmin_fix_mtc  =  &  !! values used for photosynthesis tmin 
  & (/ 0.0,    2.0,    2.0,   -4.0,   -3.0,   -2.0,   -4.0,  &       !! when STOMATE is not activated (C)
  &   -4.0,   -4.0,   -5.0,    6.0,   -5.0,    6.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: co2_topt_fix_mtc  =  &  !! values used for photosynthesis tpopt
  & (/  0.0,   27.5,   27.5,   17.5,   25.0,   20.0,   17.5,  &      !! when STOMATE is not activated (C)
  &    17.5,   17.5,   20.0,   32.5,   20.0,   32.5  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: co2_tmax_fix_mtc  =  &  !! values used for photosynthesis tmax 
  & (/  0.0,   55.0,   55.0,   38.0,   48.0,   38.0,   38.0,  &      !! when STOMATE is not activated (C)
  &    38.0,   38.0,   45.0,   55.0,   45.0,   55.0  /)
  !-
  ! 2 .Stomate
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ext_coeff_mtc  =  &     !! extinction coefficient of the Monsi&Saeki
  & (/ 0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,  &             !! relationship (1953) (unitless)
  &    0.5,   0.5,   0.5,   0.5,   0.5,   0.5  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vcmax_opt_mtc  =  &     !! Maximum rate of carboxylation (�mol.m^{-2}.s^{-1})
  & (/ undef,   65.0,    65.0,    35.0,   45.0,   55.0,   35.0,  &
  &     45.0,   35.0,    70.0,    70.0,   70.0,   70.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vjmax_opt_mtc  =  &     !! Maximum rate of RUbp regeneration (�mol.m^{-2}.s^{-1})
  & (/  undef,   130.0,   130.0,    70.0,    80.0,   110.0,   70.0,  &
  &      90.0,    70.0,   160.0,   160.0,   200.0,   200.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_min_a_mtc  =  &  !! minimum photosynthesis temperature, 
  & (/  undef,   0.0,      0.0,   0.0,   0.0,   0.0,   0.0,  &       !! constant a of ax^2+bx+c (deg C), tabulated
  &       0.0,   0.0,   0.0025,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_min_b_mtc  =  &  !! minimum photosynthesis temperature, 
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &          !! constant b of ax^2+bx+c (deg C), tabulated
  &       0.0,   0.0,   0.1,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_min_c_mtc  =  &  !! minimum photosynthesis temperature, 
  & (/  undef,    2.0,     2.0,   -4.0,   -3.0,   -2.0,   -4.0,  &   !! constant b of ax^2+bx+c (deg C), tabulated
  &      -4.0,   -4.0,   -3.25,   13.0,   -5.0,   13.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_opt_a_mtc  =  &  !! optimum photosynthesis temperature,
  & (/  undef,   0.0,      0.0,   0.0,   0.0,   0.0,   0.0,  &       !! constant a of ax^2+bx+c (deg C), tabulated
  &       0.0,   0.0,   0.0025,   0.0,   0.0,   0.0  /)

  REAL(r_std),  PARAMETER, DIMENSION(nvmc) :: tphoto_opt_b_mtc  =  & !! optimum photosynthesis temperature, 
  & (/  undef,   0.0,    0.0,   0.0,   0.0,   0.0,   0.0,  &         !! constant b of ax^2+bx+c (deg C), tabulated
  &       0.0,   0.0,   0.25,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_opt_c_mtc  =  &  !! optimum photosynthesis temperature, 
  & (/  undef,   37.0,    37.0,   25.0,   32.0,   26.0,   25.0,  &   !! constant c of ax^2+bx+c (deg C), tabulated
  &      25.0,   25.0,   27.25,   36.0,   30.0,   36.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_max_a_mtc  =  &  !! maximum photosynthesis temperature,
  & (/  undef,   0.0,       0.0,   0.0,   0.0,   0.0,   0.0,  &      !! constant a of ax^2+bx+c (deg C), tabulated
  &       0.0,   0.0,   0.00375,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_max_b_mtc  =  &  !! maximum photosynthesis temperature, 
  & (/  undef,   0.0,    0.0,   0.0,   0.0,   0.0,   0.0,  &         !! constant b of ax^2+bx+c (deg C), tabulated
  &       0.0,   0.0,   0.35,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_max_c_mtc  =  &  !! maximum photosynthesis temperature,
  & (/  undef,   55.0,     55.0,   38.0,   48.0,   38.0,   38.0,  &  !! constant c of ax^2+bx+c (deg C), tabulated
  &      38.0,   38.0,   41.125,   55.0,   45.0,   55.0  /)


  !
  ! RESPIRATION (stomate)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_c_mtc  =  &  !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.12,   0.12,   0.16,   0.16,   0.16,   0.16,  &          !! constant c of aT^2+bT+c, tabulated
  &      0.16,   0.16,   0.16,   0.12,   0.16,   0.12  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_b_mtc  =  &  !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.0,        0.0,   0.0,        0.0,   0.0,   0.0,  &      !! constant b of aT^2+bT+c, tabulated
  &       0.0,   0.0,   -0.00133,   0.0,   -0.00133,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_a_mtc  =  &  !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                !! constant a of aT^2+bT+c, tabulated
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_leaf_mtc  =   &                  !! maintenance respiration coefficient
  & (/   undef,   2.35E-3,   2.62E-3,   1.01E-3,   2.35E-3,   2.62E-3,   1.01E-3,  &  !! at 0 deg C,for leaves, tabulated, 
  &    2.62E-3,   2.05E-3,   2.62E-3,   2.62E-3,   2.62E-3,   2.62E-3  /)             !! (gC.gC^{-1}.day^{-1})

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_sapabove_mtc =  &                !! maintenance respiration coefficient 
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for sapwood above,
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! tabulated, (gC.gC^{-1}.day^{-1}) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_sapbelow_mtc  =  &               !! maintenance respiration coefficient
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for sapwood below, 
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! tabulated, (gC.gC^{-1}.day^{-1})  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_heartabove_mtc  =  &             !! maintenance respiration coefficient
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                           !! at 0 deg C, for heartwood above,
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)                                  !! tabulated, (gC.gC^{-1}.day^{-1}) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_heartbelow_mtc  =  &             !! maintenance respiration coefficient
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                           !! at 0 deg C, for heartwood below, 
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)                                  !! tabulated, (gC.gC^{-1}.day^{-1})  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_root_mtc  =  &                   !! maintenance respiration coefficient
  & (/   undef,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,  &  !! at 0 deg C, for roots, tabulated,
  &    1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3,   1.67E-3  /)             !! (gC.gC^{-1}.day^{-1})  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_fruit_mtc  =  &                  !! maintenance respiration coefficient
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for fruits, tabulated,
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! (gC.gC^{-1}.day^{-1})

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cm_zero_carbres_mtc  =  &                !! maintenance respiration coefficient
  & (/   undef,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,  &  !! at 0 deg C, for carbohydrate reserve,
  &    1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4,   1.19E-4  /)             !! tabulated, (gC.gC^{-1}.day^{-1})   


  !
  ! FIRE (stomate)
  !
  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: flam_mtc  =  &         !! flamability: critical fraction of water 
  & (/  undef,   0.15,   0.25,   0.25,   0.25,   0.25,   0.25,  &  !! holding capacity (0-1, unitless)
  &      0.25,   0.25,   0.25,   0.25,   0.35,   0.35  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: resist_mtc  =  &       !! fire resistance (0-1, unitless)
  & (/ undef,   0.95,   0.90,   0.12,   0.50,   0.12,   0.12,  &
  &    0.12,    0.12,    0.0,    0.0,    0.0,    0.0 /) 


  !
  ! FLUX - LUC
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_1_mtc  =  &   !! Coeff of biomass export for the year
  & (/  undef,   0.597,   0.597,   0.597,   0.597,   0.597,   0.597,  &   !! (0-1, unitless)
  &     0.597,   0.597,   0.597,   0.597,   0.597,   0.597  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_10_mtc  =  &  !! Coeff of biomass export for the decade 
  & (/  undef,   0.403,   0.403,   0.299,   0.299,   0.299,   0.299,  &   !! (0-1, unitless)
  &     0.299,   0.299,   0.299,   0.403,   0.299,   0.403  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_100_mtc  =  & !! Coeff of biomass export for the century
  & (/  undef,     0.0,     0.0,   0.104,   0.104,   0.104,   0.104,  &   !! (0-1, unitless)
  &     0.104,   0.104,   0.104,     0.0,   0.104,     0.0  /)


  !
  ! PHENOLOGY
  !
  !-
  ! 1. Stomate
  !-
  REAL(r_std), PARAMETER, DIMENSION (nvmc) :: lai_max_mtc  =  &          !! maximum LAI, PFT-specific (m^2.m^{-2})
  & (/ undef,   7.0,   7.0,   5.0,   5.0,   5.0,   4.5,  &
  &      4.5,   3.0,   2.5,   2.5,   5.0,   5.0  /)

  CHARACTER(len=6), PARAMETER, DIMENSION(nvmc) :: pheno_model_mtc  =  &  !! which phenology model is used? (tabulated) 
  & (/  'none  ',   'none  ',   'moi   ',   'none  ',   'none  ',  &
  &     'ncdgdd',   'none  ',   'ncdgdd',   'ngd   ',   'moigdd',  &
  &     'moigdd',   'moigdd',   'moigdd'  /) 

  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: pheno_type_mtc  =  &     !! type of phenology (0-4, unitless)
  & (/  0,   1,   3,   1,   1,   2,   1,  &                              !! 0=bare ground 1=evergreen,  2=summergreen, 
  &     2,   2,   4,   4,   2,   3  /)                                   !! 3=raingreen,  4=perennial
  !-
  ! 2. Leaf Onset
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_c_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! constant c of aT^2+bT+c
  &     undef,   undef,   270.0,   400.0,   125.0,   400.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_b_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! constant b of aT^2+bT+c
  &     undef,   undef,    6.25,     0.0,     0.0,     0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_a_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,     undef,   undef,   undef,   undef,   undef,  &  !! constant a of aT^2+bT+c
  &     undef,   undef,   0.03125,     0.0,     0.0,     0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ngd_crit_mtc  =  &            !! critical ngd, tabulated. 
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! Threshold -5 degrees (days)
  &     undef,    17.0,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ncdgdd_temp_mtc  =  &         !! critical temperature for the ncd vs. gdd 
  & (/  undef,   undef,   undef,   undef,   undef,     5.0,   undef,  &    !! function in phenology (C)
  &       0.0,   undef,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: hum_frac_mtc  =  &            !! critical humidity (relative to min/max) 
  & (/  undef,   undef,   0.5,   undef,   undef,   undef,   undef, &       !! for phenology (unitless)
  &     undef,   undef,   0.5,     0.5,     0.5,     0.5  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lowgpp_time_mtc  =  &         !! minimum duration of dormance
  & (/  undef,   undef,   30.0,   undef,   undef,   30.0,   undef,  &      !! for phenology (days)
  &      30.0,    30.0,   30.0,    30.0,    30.0,   30.0  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: hum_min_time_mtc  =  &        !! minimum time elapsed since
  & (/  undef,   undef,   50.0,   undef,   undef,   undef,   undef,  &     !! moisture minimum (days)
  &     undef,   undef,   35.0,    35.0,    75.0,    75.0  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tau_sap_mtc  =  &             !! time (days)  
  & (/  undef,   730.0,   730.0,   730.0,   730.0,   730.0,   730.0,  &
  &     730.0,   730.0,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tau_fruit_mtc  =  &           !! fruit lifetime (days)
  & (/  undef,  90.0,    90.0,    90.0,    90.0,   90.0,   90.0,  &
  &      90.0,  90.0,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ecureuil_mtc  =  &            !! fraction of primary leaf and root allocation
  & (/  undef,   0.0,   1.0,   0.0,   0.0,   1.0,   0.0,  &                !! put into reserve (0-1, unitless)
  &       1.0,   1.0,   1.0,   1.0,   1.0,   1.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alloc_min_mtc  =  &           !! NEW - allocation above/below = f(age) 
  & (/  undef,   0.2,     0.2,     0.2,     0.2,    0.2,   0.2,  &         !! - 30/01/04 NV/JO/PF
  &       0.2,   0.2,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alloc_max_mtc  =  &           !! NEW - allocation above/below = f(age) 
  & (/  undef,   0.8,     0.8,     0.8,     0.8,    0.8,   0.8,  &         !! - 30/01/04 NV/JO/PF
  &       0.8,   0.8,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: demi_alloc_mtc  =  &          !! NEW - allocation above/below = f(age) 
  & (/  undef,   5.0,     5.0,     5.0,     5.0,    5.0,   5.0,  &         !! - 30/01/04 NV/JO/PF
  &       5.0,   5.0,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaflife_mtc  =  &            !! leaf longevity, tabulated (??units??)
  & (/  undef,   0.5,   2.0,   0.33,   1.0,   2.0,   0.33,  &
  &       2.0,   2.0,   2.0,   2.0,    2.0,   2.0  /)
  !-
  ! 3. Senescence
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaffall_mtc  =  &             !! length of death of leaves, tabulated (days)
  & (/  undef,   undef,   10.0,   undef,   undef,   10.0,   undef,  &
  &      10.0,    10.0,   10.0,    10.0,    10.0,   10.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leafagecrit_mtc  =  &          !! critical leaf age, tabulated (days)
  & (/  undef,   730.0,   180.0,   910.0,   730.0,   180.0,   910.0,  &
  &     180.0,   180.0,   120.0,   120.0,    90.0,    90.0  /)

  CHARACTER(LEN=6), PARAMETER, DIMENSION(nvmc) :: senescence_type_mtc  =  & !! type of senescence, tabulated (unitless)
  & (/  'none  ',  'none  ',   'dry   ',  'none  ',  'none  ',  &
  &     'cold  ',  'none  ',   'cold  ',  'cold  ',  'mixed ',  &
  &     'mixed ',  'mixed ',   'mixed '            /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_hum_mtc  =  &       !! critical relative moisture availability
  & (/  undef,   undef,   0.3,   undef,   undef,   undef,   undef,  &       !! for senescence (0-1, unitless)
  &     undef,   undef,   0.2,     0.2,     0.3,     0.2  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: nosenescence_hum_mtc  =  &     !! relative moisture availability above which 
  & (/  undef,   undef,   0.8,   undef,   undef,   undef,   undef,  &       !! there is no humidity-related senescence
  &     undef,   undef,   0.3,     0.3,     0.3,     0.3  /)                !! (0-1, unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: max_turnover_time_mtc  =  &    !! maximum turnover time for grasses (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &
  &     undef,   undef,    80.0,    80.0,    80.0,    80.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: min_turnover_time_mtc  =  &    !! minimum turnover time for grasses (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &
  &     undef,   undef,    10.0,    10.0,    10.0,    10.0  /)
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: min_leaf_age_for_senescence_mtc  =  &  !! minimum leaf age to allow 
  & (/  undef,   undef,   90.0,   undef,   undef,   90.0,   undef,  &               !! senescence g (days)
  &      60.0,    60.0,   30.0,    30.0,    30.0,   30.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_temp_c_mtc  =  &    !! critical temperature for senescence (C)
  & (/  undef,   undef,    undef,   undef,   undef,   12.0,   undef,  &     !! constant c of aT^2+bT+c, tabulated
  &       7.0,     2.0,   -1.375,     5.0,    5.0,    10.0  /)              !! (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_temp_b_mtc  =  &    !! critical temperature for senescence (C), 
  & (/  undef,   undef,   undef,   undef,   undef,   0.0,   undef,  &       !! constant b of aT^2+bT+c, tabulated
  &       0.0,     0.0,     0.1,     0.0,     0.0,   0.0  /)                !! (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_temp_a_mtc  =  &    !! critical temperature for senescence (C), 
  & (/  undef,   undef,     undef,   undef,   undef,   0.0,   undef,  &     !! constant a of aT^2+bT+c, tabulated
  &       0.0,     0.0,   0.00375,     0.0,     0.0,   0.0  /)              !! (unitless)


  !
  ! DGVM
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: residence_time_mtc  =  &    !! residence time of trees (years)
  & (/  undef,   30.0,   30.0,   40.0,   40.0,   40.0,   80.0,  &
  &      80.0,   80.0,    0.0,    0.0,    0.0,    0.0  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tmin_crit_mtc  =  &
  & (/  undef,     0.0,     0.0,   -30.0,   -14.0,   -30.0,   -45.0,  &  !! critical tmin, tabulated (C)
  &     -45.0,   undef,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tcm_crit_mtc  =  &
  & (/  undef,   undef,   undef,     5.0,    15.5,    15.5,   -8.0,  &   !! critical tcm, tabulated (C)
  &      -8.0,    -8.0,   undef,   undef,   undef,   undef  /)



  !
  ! Biogenic Volatile Organic Compounds
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_isoprene_mtc = &     !! Isoprene emission factor 
  & (/  0.,    24.,   24.,    8.,   16.,   45.,   8.,  &                    !! (\mu gC.g^{-1}.h^{-1}) 
  &     8.,     8.,   16.,   24.,    5.,    5.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_monoterpene_mtc = &  !! Monoterpene emission factor
  & (/   0.,   0.8,    0.8,   2.4,    1.2,    0.8,    2.4,  &               !! (\mu gC.g^{-1}.h^{-1}) 
  &    2.4,    2.4,    0.8,   1.2,    0.2,     0.2  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_ORVOC_mtc = &        !! ORVOC emissions factor 
  &  (/  0.,    1.5,    1.5,    1.5,    1.5,   1.5,    1.5,  &              !! (\mu gC.g^{-1}.h^{-1}) 
  &     1.5,    1.5,    1.5,    1.5,    1.5,   1.5  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_OVOC_mtc = &         !! OVOC emissions factor 
  &  (/  0.,    1.5,    1.5,    1.5,    1.5,   1.5,    1.5,  &              !! (\mu gC.g^{-1}.h^{-1})
  &     1.5,    1.5,    1.5,    1.5,    1.5,   1.5  /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_MBO_mtc = &          !! MBO emissions factor
  & (/  0.,    0.,    0.,    20.0,    0.,    0.,    0.,  &                  !! (\mu gC.g^{-1}.h^{-1}) 
  &     0.,    0.,    0.,       0.,   0.,    0.  /)  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_methanol_mtc = &     !! Methanol emissions factor 
  & (/  0.,    0.6,   0.6,   1.8,   0.9,   0.6,   1.8,  &                   !! (\mu gC.g^{-1}.h^{-1}) 
  &    1.8,    1.8,   0.6,   0.9,    2.,     2.  /)  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_acetone_mtc = &      !! Acetone emissions factor
  & (/  0.,   0.29,   0.29,   0.87,   0.43,   0.29,   0.87,  &              !! (\mu gC.g^{-1}.h^{-1}) 
  &   0.87,   0.87,   0.29,   0.43,   0.07,   0.07  /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_acetal_mtc = &       !! Acetaldehyde emissions factor 
  & (/  0.,   0.1,    0.1,     0.3,    0.15,   0.1,   0.3,   0.3,   &       !! (\mu gC.g^{-1}.h^{-1})
  &    0.3,   0.1,   0.15,   0.025,   0.025  /)  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_formal_mtc = &       !! Formaldehyde emissions factor
  & (/  0.,   0.07,   0.07,   0.2,     0.1,    0.07,   0.2,  &              !! (\mu gC.g^{-1}.h^{-1})
  &    0.2,    0.2,   0.07,   0.1,   0.017,   0.017  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_acetic_mtc = &       !! Acetic Acid emissions factor
  & (/   0.,   0.002,   0.002,   0.006,   0.003,     0.002,   0.006,   &    !! (\mu gC.g^{-1}.h^{-1})
  &   0.006,   0.006,   0.002,   0.003,   0.0005,   0.0005  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_formic_mtc = &       !! Formic Acid emissions factor
  & (/  0.,   0.01,   0.01,   0.03,     0.015,     0.01,   0.03,  &         !! (\mu gC.g^{-1}.h^{-1})
  &   0.03,   0.03,   0.01,   0.015,   0.0025,   0.0025  /)  

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: em_factor_no_wet_mtc = &        !! NOx emissions factor soil emissions and exponential
  & (/  0.,   2.6,   0.06,   0.03,   0.03,   0.03,   0.03,  &               !! dependancy factor for wet soils (ngN.m^{-2}.s^{-1}) 
  &  0.03,   0.03,   0.36,   0.36,   0.36,   0.36  /)  

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: em_factor_no_dry_mtc = &        !! NOx emissions factor soil emissions and exponential
  & (/  0.,   8.60,   0.40,   0.22,   0.22,   0.22,   0.22,  &              !! dependancy factor for dry soils (ngN.m^{-2}.s^{-1}) 
  &   0.22,   0.22,   2.65,   2.65,   2.65,   2.65  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Larch_mtc = &                  !! Larcher 1991 SAI/LAI ratio (unitless)
  & (/   0.,   0.015,   0.015,   0.003,   0.005,   0.005,   0.003,  &
  &   0.005,   0.003,   0.005,   0.005,   0.008,   0.008  /)  



END MODULE constantes_mtc
