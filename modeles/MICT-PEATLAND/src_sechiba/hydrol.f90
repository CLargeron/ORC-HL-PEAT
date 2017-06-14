

! MODULE        : hydrol
!
! CONTACT       : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module computes the soil moisture processes on continental points. 
!!
!!\n DESCRIPTION : contains hydrol_init, hydrol_var_init, hydrol_waterbal, hydrol_alma,
!!                 hydrol_snow, hydrol_vegupd, hydrol_canop, hydrol_flood, hydrol_soil.
!!                 The assumption in this module is that very high vertical resolution is
!!                 needed in order to properly resolve the vertical diffusion of water in
!!                 the soils. Furthermore we have taken into account the sub-grid variability
!!                 of soil properties and vegetation cover by allowing the co-existence of
!!                 different soil moisture columns in the same grid box.
!!                 This routine was originaly developed by Patricia deRosnay.
!! 
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/hydrol.f90 $
!! $Date: 2013-03-19 18:15:15 +0100 (Tue, 19 Mar 2013) $
!! $Revision: 1223 $
!! \n
!_ ================================================================================================================================

MODULE hydrol

  USE ioipsl
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE sechiba_io
  USE slowproc
  USE grid

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hydrol_main,hydrol_clear 

  !
  ! variables used inside hydrol module : declaration and initialisation
  !
  LOGICAL, SAVE                                   :: l_first_hydrol=.TRUE.   !! Initialisation has to be done one time
  LOGICAL, SAVE                                   :: l_second_hydrol=.TRUE.  !! Initialisation has to be done one time
  !
  LOGICAL, SAVE                                   :: check_cwrr=.FALSE.      !! The check the water balance
  LOGICAL, SAVE                                   :: doponds=.FALSE.         !! Reinfiltration param.
  REAL, SAVE                                      :: max_stagnant            !! Max eau ds stagnant (mm) ! Chloe050813
  !
  CHARACTER(LEN=80) , SAVE                        :: var_name                !! To store variables names for I/O
  !
  REAL(r_std), PARAMETER                          :: drain_rest_cste = 15.0  !! time constant in days to return to free drainage 
                                                                             !! after return flow 
  REAL(r_std), PARAMETER                          :: allowed_err =  2.0E-8_r_std
  REAL(r_std), PARAMETER                          :: EPS1 = EPSILON(un)      !! A small number
  ! one dimension array allocated, computed, saved and got in hydrol module
  ! Values per soil type
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: nvan                !! Van Genuchten coeficients n
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: avan                !! Van Genuchten coeficients a
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcr                 !! Residual humidity [m3/m3]
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs                 !! Saturation humidity [m3/m3]
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: ks                  !! Hydraulic conductivity Saturation for each soil type 
                                                                         !! [mm/day] 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: ds                  !! Hydraulic diffusivity for each soil type and soil level
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pcent               !! Fraction of saturated volumetric soil moisture above 
                                                                         !! which transpir is max 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: free_drain_max      !! Max value of the permeability coeff at the bottom of 
                                                                         !! the soil [mm/day] 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcf                 !! Field capacity [m3/m3]
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcw                 !! Wilting point [m3/m3]
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_awet             !! Vol. wat. cont. above which albedo is cst [m3/m3]
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_adry             !! Vol. wat. cont. below which albedo is cst [m3/m3]
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: psis                !! Matrix potential at saturation [mm]

  ! Values per grid point
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_water_beg    !! Total amount of water at start of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_water_end    !! Total amount of water at end of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_flux         !! Total water flux
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_beg   !! Total amount of water on vegetation at start of time 
                                                                         !! step 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_end   !! Total amount of water on vegetation at end of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_beg  !! Total amount of water in the soil at start of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_end  !! Total amount of water in the soil at end of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_beg         !! Total amount of snow at start of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_end         !! Total amount of snow at end of time step
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delsoilmoist     !! Change in soil moisture
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delintercept     !! Change in interception storage
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delswe           !! Change in SWE
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: swi              !! Soil Wetness Index
  

  ! array allocated, computed, saved and got in hydrol module
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_veget       !! zero/one when veget fraction is zero/higher
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_soiltile    !! zero/one where soil tile is zero/higher 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: humrelv          !! humrel for each soil type
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vegstressv       !! vegstress for each soil type  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:):: us               !! relative humidity 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol         !! Eau tombee sur le sol
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol_ns      !! Eau tombee sur le sol par type de sol
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ae_ns            !! Evaporation du sol nu par type de sol
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: stagnant         !! Reservoir d'eau stagnante Chloe050813
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: evap_bare_lim_ns !! limitation of bare soil evaporation on each soil column 
                                                                         !! (used to deconvoluate vevapnu) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: free_drain_coef  !! Coefficient for free drainage at bottom
 ! Chloe   
!d�plac� dans sechiba.f90
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: wt_soil        !! Water table impose 
 ! 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_bare_ns     !! evaporating bare soil fraction per tile
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: rootsink         !! stress racinaire par niveau et type de sol
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsnowveg       !! Sublimation of snow on vegetation
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: subsnownobio     !! Sublimation of snow on other surface types (ice, lakes, 
                                                                         !! ...) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snowmelt         !! Quantite de neige fondue
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: icemelt          !! Quantite de glace fondue
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsinksoil      !! Excess of sublimation as a sink for the soil
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: vegtot           !! Total  vegetation (veget_max)
  ! The last vegetation map which was used to distribute the reservoirs
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: resdist          !! Distribution of reservoirs
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: vmr              !! variation of veget
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mx_eau_var       !!

  ! arrays used by cwrr scheme
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: nroot            !! nvm * nstm * nslm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kfact_root       !! kjpindex * nslm * nstm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kfact            !! nslm * nscm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: zz               !! nslm+1 * nstm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: dz               !! nslm+1 * nstm

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: mc_lin           !! imin:imax * nscm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: k_lin            !! imin:imax * nslm * nscm
!Isa
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: k_lin_write              !!kjpindex*imin:imax*nscm
!end Isa            !! imin:imax * nslm * nscm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: d_lin            !! imin:imax * nslm * nscm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: a_lin            !! imin:imax * nslm * nscm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: b_lin            !! imin:imax * nslm * nscm

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: humtot           !! (:) Total Soil Moisture, Kg/m2
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: flux             !! (:)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION (:)          :: resolv           !! (:)

!Chloe carte porosit� mcs : 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mcs_njsc         

!! linarization coefficients of hydraulic conductivity K
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: k     !! (:,nslm)
!Isa : conductivit� hydraulique moyenne sur la maille, tenant compte des soiltiles ie des pr�f de la v�g�tation
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kk_moy           !! (:,nslm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)    :: kk           !! (:,nslm, nstm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: a                !! (:,nslm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: b                !!
!! linarization coefficients of hydraulic diffusivity D
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: d                !!

!! matrix coefficients
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: e                !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: f                !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: g1               !!


  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ep               !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: fp               !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: gp               !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: rhs              !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: srhs             !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: gam              !!

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: water2infilt     !! Water to be infiltrated
  !Chloe
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_obj_peat     !! tmc of peat we want to get              

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc              !! (:,nstm) Total moisture content (mm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcr             !! (nstm) Total moisture constent at residual (mm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcs             !! (nstm) Total moisture constent at saturation (mm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter       !! (:,nstm) Total moisture in the litter by soil type
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_mea     !! Total moisture in the litter over the grid

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_wilt  !! (:,nstm) Moisture of litter at wilt pt
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_field !! (:,nstm) Moisture of litter at field cap.
!!! A CHANGER DANS TOUT HYDROL: tmc_litter_res et sat ne devraient pas dependre de ji - tdo
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_res   !! (:,nstm) Moisture of litter at residual moisture.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_sat   !! (:,nstm) Moisture of litter at saturatiion
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_awet  !! (:,nstm) Moisture of litter at mc_awet
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_adry  !! (:,nstm) Moisture of litter at mc_dry
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_wet_mea !! Total moisture in the litter over the grid below which 
                                                                         !! albedo is fixed 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_dry_mea !! Total moisture in the litter over the grid above which 
                                                                         !! albedo is fixed 
  LOGICAL, SAVE                                      :: tmc_init_updated = .FALSE. !! Flag allowing to determine if tmc is initialized.

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: v1               !! (:)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: qflux00          !! flux at the top of the soil column

  !! par type de sol :
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ru_ns            !! ruissellement
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: dr_ns            !! drainage
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tr_ns            !! transpiration 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: cvs_over_veg     !! (:,nvm,nstm) old value of corr_veg_soil/veget_max kept 
                                                                         !! from diag to next split 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: corr_veg_soil    !! (:,nvm,nstm) percentage of each veg. type on each soil 
                                                                         !! of each grid point 
!Chloe
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: mc_obj_peat      !! ideal  mc of peat we want            
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: mc               !! (:,nslm,nstm) \f($m^3 \times m^3$)\f
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soilmoist        !! (:,nslm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: soil_wet         !! soil wetness
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soil_wet_litter  !! soil wetness of the litter
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: qflux            !! fluxes between the soil layers
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: tmat             !!
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: stmat            !!
!Isa pour gel dans le sol 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_hydro_diag
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: profil_froz_hydro
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)    :: profil_froz_hydro_ns
  !par soiltile
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: temp_hydro ! temp profile on hydrological levels
  LOGICAL, SAVE                                     :: ok_freeze_cwrr !! CWRR freezing scheme by I. Gouttevin
  LOGICAL, SAVE                                     :: ok_gel_thd, ok_gel_thermosoil
  LOGICAL, SAVE                                     :: ok_shumdiag_perma !corresponding options
  LOGICAL, SAVE                                     :: ok_shumdiag_interpol 
        ! cle ajoutee par CR: modif d'Isa pour interpoler shumdiag sur l'axe
        ! diag
 REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mcl !! liquid moisture content
 INTEGER(i_std), SAVE                         :: kjitt             !! Time step number

 !Chloe++ Peatland PFT14
 !LOGICAL, SAVE                                     :: ok_routage_peat !! d�fini dans sechiba 
 LOGICAL, SAVE                                     :: ok_sat30cm      !! soil is saturated below 30cm (30-45cm) depth.
 LOGICAL, SAVE                                     :: ok_stagnant     !! Stagnant water reservoir (10cm upper the surface)- Limited Runoff 
 !d�fini dans s�chiba : 
 !LOGICAL, SAVE                                     :: ok_no_drainage  !! No drainage in peatland case
 LOGICAL, SAVE                                     :: ok_reinfilt_peat !! Reinfiltration of soiltile 1 to 3 to soiltile 4
 !Chloe--

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_main
!!
!>\BRIEF         
!!
!! DESCRIPTION :
!! - called only one time for initialisation
!! - called every time step
!! - called one more time at last time step for writing _restart_ file
!!
!! - Choose between initialisation/Restart time/Time Step
!! - 1 Do initialisation ==> hydrol_init
!! - X if check_waterbal ==> hydrol_waterbal
!! - 2 prepares restart file for the next simulation
!! - 3 shared time step
!! - 3.1 computes snow  ==> hydrol_snow
!! - 3.2 computes vegetations reservoirs  ==> hydrol_vegupd
!! - 3.3 computes canopy  ==> hydrol_canop
!! - 3.4 computes surface reservoir  ==> hydrol_flood
!! - 3.5 computes soil hydrologie ==> hydrol_soil
!! - X if check_waterbal ==> hydrol_waterbal
!! - 4 write out file  ==> hydrol_alma/histwrite(*)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, &
       & index, indexveg, indexsoil, indexlayer, control_in, &
       & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc, &
       & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,  &
       & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, &
       & humrel, vegstress, drysoil_frac, evapot, evapot_penm, evap_bare_lim, flood_frac, flood_res, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, reinf_slope,  &
       & rest_id, hist_id, hist2_id,&
       & stempdiag, peatland, ok_routage_peat,ok_no_drainage,water2add_peat, mc_lin_axis_index, wt_soil, wt_soil2 ) !Chloe
!Isa: shumdiag_perma, stempdiag

    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    ! input scalar 
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                         :: rest_id,hist_id  !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                         :: hist2_id         !! _history_ file 2 identifier
    REAL(r_std), INTENT (in)                           :: dtradia          !! Time step in seconds
    LOGICAL, INTENT(in)                                :: ldrestart_read   !! Logical for _restart_ file to read
    LOGICAL, INTENT(in)                                :: ldrestart_write  !! Logical for _restart_ file to write
    !Chloe
    LOGICAL, INTENT(in), DIMENSION(kjpindex)           :: peatland         !! logical if peatland is dominant in a grid cell
    !Chloe
    LOGICAL, INTENT(in)                                :: ok_routage_peat     
    LOGICAL, INTENT(in)                                :: ok_no_drainage 
    ! input fields

    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index            !! Indeces of the points on the map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg        !! Indeces of the points on the 3D map for veg
    INTEGER(i_std),DIMENSION (kjpindex*nstm), INTENT (in):: indexsoil      !! Indeces of the points on the 3D map for soil
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexlayer     !! Indeces of the points on the 3D map for soil layers
    TYPE(control_type), INTENT (in)                    :: control_in       !! Flags that (de)activate parts of the model
    !
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain      !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_snow      !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: returnflow       !! Routed water which comes back into the soil (from the 
                                                                           !! bottom) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinfiltration   !! Routed water which comes back into the soil (at the 
                                                                           !! top) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: irrigation       !! Water from irrigation returning to soil moisture  
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol_new     !! New soil temperature

    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc             !! indexing of PFT to soiltile
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio    !! Total fraction of ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: soilcap          !! Soil capacity
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! fraction of PFT on each soil-hydrology tile
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet         !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type           
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax         !! Maximum water on vegetation for interception
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir         !! Transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: reinf_slope      !! Slope coef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot           !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_penm      !! Soil Potential Evaporation Correction
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac       !! flood fraction

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: humrel           !! Relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vegstress        !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out):: shumdiag         !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: k_litt           !! litter approximate conductivity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: litterhumdiag    !! litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: tot_melt         !! Total melt    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: qsintveg         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)      :: evap_bare_lim    !!
    !Chloe
   REAL(r_std),DIMENSION (kjpindex), INTENT(out)      :: water2add_peat    !!  
   REAL(r_std),DIMENSION (kjpindex,nstm), INTENT(out)      :: wt_soil      !!Water Table position (mm)       
 REAL(r_std),DIMENSION (kjpindex,nstm), INTENT(out)      :: wt_soil2      !!Water Table position (mm)  
!   REAL(r_std),DIMENSION (kjpindex), INTENT(out)      :: wtold        !! Water table pas de temps precedent
    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapnu          !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapsno         !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapflo         !! Floodplain evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: flood_res        !! flood reservoir estimate
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio  !! Water balance on ice, lakes, .. [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio_age !! Snow age on ice, lakes, ...
    !! We consider that any water on the ice is snow and we only peforme a water balance to have consistency.
    !! The water balance is limite to + or - 10^6 so that accumulation is not endless
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: floodout       !! flux out of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: runoff         !! Complete runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: drainage         !! Water on vegetation due to interception
!Isa
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in):: stempdiag        !! diagnostic temp profile from thermosoil
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out):: shumdiag_perma  !! % of porosity filled with water (mc/mcs) used for the thermal computations 

    !! Drainage

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jst, jsl
    REAL(r_std),DIMENSION (kjpindex)                   :: soilwet          !! A temporary diagnostic of soil wetness
    REAL(r_std),DIMENSION (kjpindex)                   :: snowdepth        !! Depth of snow layer
    REAL(r_std),DIMENSION (kjpindex)                   :: njsc_tmp         !! Temporary REAL value for njsc to write it
!Isa
 INTEGER(i_std), DIMENSION(kjpindex*imax)	   :: mc_lin_axis_index
 INTEGER(i_std)                                :: ji, jsc, i

!Chloe reinfiltration runoff des soiltile 1 � 3
   REAL(r_std)     :: reinfilt_watpeat !! Reinfilt runoff(jst=1:3) vers soiltile 4
   !REAL(r_std)     :: fracreinfilt     !! Fraction du ruissellement a reinfiltrer au peat         
    !!_hydrol_main

    !
    !! Choose between initialisation/Restart time/Time Step
    !! 1 Do initialisation ==>hydrol_init
    !
!Isa

kjitt = kjit

    IF (l_first_hydrol) THEN

       IF (long_print) WRITE (numout,*) ' l_first_hydrol : call hydrol_init '

       write(*,*) 'hydrol 331: shumdiag=',shumdiag(1,1)
       CALL hydrol_init (kjit, ldrestart_read, kjpindex, index, rest_id, veget_max, soiltile, humrel,&
            & njsc, vegstress, snow, snow_age, snow_nobio, snow_nobio_age, qsintveg, peatland, ok_routage_peat, wt_soil,wt_soil2) 
!Isa
       CALL hydrol_var_init (kjpindex, veget, veget_max, &
            & soiltile, njsc, mx_eau_var, shumdiag,shumdiag_perma, k_litt, &
            & litterhumdiag, drysoil_frac, evap_bare_lim, ok_routage_peat) 
       !write(*,*) 'hydrol 335: shumdiag=',shumdiag(1,1)

       ! If we check the water balance we first save the total amount of water
    !! X if check_waterbal ==> hydrol_waterbal
    !CHLOE TEST commente if check_waterbal   
    IF (check_waterbal) THEN
          CALL hydrol_waterbal(kjpindex, index, .TRUE., dtradia, veget_max, &
               & totfrac_nobio, qsintveg, snow, snow_nobio,&
               & precip_rain, precip_snow, returnflow, reinfiltration, irrigation, tot_melt, &
               & vevapwet, transpir, vevapnu, vevapsno, vevapflo, floodout, runoff, drainage)
       ENDIF
       !
       IF (almaoutput) THEN
          CALL hydrol_alma(kjpindex, index, .TRUE., qsintveg, snow, snow_nobio, soilwet)
       ENDIF

       RETURN

    ENDIF

    !
    !! 2 prepares restart file for the next simulation
    !
    IF (ldrestart_write) THEN

       IF (long_print) WRITE (numout,*) ' we have to complete restart file with HYDROLOGIC variables '

       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
          WRITE (var_name,"('moistc_',i1)") jst
          CALL restput_p(rest_id, var_name, nbp_glo,  nslm, 1, kjit, mc(:,:,jst), 'scatter',  nbp_glo, index_g)
       END DO
       !
       DO jst=1,nstm
          DO jsl=1,nslm
             ! var_name= "us_1_01" ... "us_3_11"
             WRITE (var_name,"('us_',i1,'_',i2.2)") jst,jsl
             CALL restput_p(rest_id, var_name, nbp_glo,nvm, 1,kjit,us(:,:,jst,jsl),'scatter',nbp_glo,index_g)
          END DO
       END DO
       !
       var_name= 'free_drain_coef'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  free_drain_coef, 'scatter',  nbp_glo, index_g)
       !
!Chloe++
       var_name= 'wt_soil'
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  wt_soil,'scatter',  nbp_glo, index_g)
       !
       var_name= 'wt_soil2'
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  wt_soil2,'scatter',  nbp_glo, index_g)
!Chloe--
       var_name= 'water2infilt'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  water2infilt, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'ae_ns'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  ae_ns, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'stagnant'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  stagnant, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'vegstress'
       CALL restput_p(rest_id, var_name, nbp_glo,   nvm, 1, kjit,  vegstress, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'snow'    
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  snow, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'snow_age'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  snow_age, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'snow_nobio'    
       CALL restput_p(rest_id, var_name, nbp_glo,   nnobio, 1, kjit,  snow_nobio, 'scatter', nbp_glo, index_g)
       !
       var_name= 'snow_nobio_age'
       CALL restput_p(rest_id, var_name, nbp_glo,   nnobio, 1, kjit,  snow_nobio_age, 'scatter', nbp_glo, index_g)
       !
       var_name= 'qsintveg'
       CALL restput_p(rest_id, var_name, nbp_glo, nvm, 1, kjit,  qsintveg, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'resdist'
       CALL restput_p(rest_id, var_name, nbp_glo, nvm, 1, kjit,  resdist, 'scatter',  nbp_glo, index_g)       
       !
       DO jst=1,nstm
          ! var_name= "cvs_over_veg_1" ... "cvs_over_veg_3"
          WRITE (var_name,"('cvs_over_veg_',i1)") jst
          CALL restput_p(rest_id, var_name, nbp_glo,  nvm, 1, kjit, cvs_over_veg(:,:,jst), 'scatter',  nbp_glo, index_g)
       END DO
       !
       IF ( check_waterbal ) THEN
          var_name= 'tot_water_beg'
          CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  tot_water_end, 'scatter', nbp_glo, index_g)
       ENDIF
       !
       RETURN
       !
    END IF

    !
    !!__3 shared time step
    !
    IF (long_print) WRITE (numout,*) 'hydrol pas de temps = ',dtradia

    ! 
    !! 3.1 computes snow  ==>hydrol_snow
    CALL hydrol_snow(kjpindex, dtradia, precip_rain, precip_snow, temp_sol_new, soilcap, &
         & frac_nobio, totfrac_nobio, vevapnu, vevapsno, snow, snow_age, snow_nobio, snow_nobio_age, &
         & tot_melt, snowdepth)

    !
    !! 3.2 computes vegetations reservoirs  ==>hydrol_vegupd
    CALL hydrol_vegupd(kjpindex, veget, veget_max, soiltile, qsintveg,resdist)

    !
    !! 3.3 computes canopy  ==>hydrol_canop
    CALL hydrol_canop(kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, qsintveg,precisol,tot_melt)

    !
    !!___3.4 computes surface reservoir  ==>hydrol_flood
    CALL hydrol_flood(kjpindex, dtradia, vevapflo, flood_frac, flood_res, floodout)

    !
    !! 3.5 computes soil hydrologie ==>hydrol_soil
!Isa : shumdiag_perma, stempdiag
!    write(*,*) 'hydrol 479: shumdiag=',shumdiag(1,1)
    CALL hydrol_soil(kjpindex, dtradia, veget_max, soiltile, njsc, reinf_slope,  &
         & transpir, vevapnu, evapot, evapot_penm, runoff,   &
         & drainage, returnflow, reinfiltration, irrigation, &
         !Isa
         & tot_melt,evap_bare_lim, shumdiag,shumdiag_perma, &
         & k_litt, litterhumdiag, humrel, vegstress, drysoil_frac,&
         & stempdiag,snow, peatland,water2add_peat,&
         &ok_routage_peat,ok_no_drainage, wt_soil,wt_soil2)
         
!    write(*,*) 'hydrol 487: shumdiag=',shumdiag(1,1)
!    write(*,*) 'shumdiag_perma=',shumdiag_perma(1,1)

    ! If we check the water balance we end with the comparison of total water change and fluxes
    !! X if check_waterbal ==> hydrol_waterbal
    IF (check_waterbal) THEN
       CALL hydrol_waterbal(kjpindex, index, .FALSE., dtradia, veget_max, totfrac_nobio, &
            & qsintveg, snow,snow_nobio, precip_rain, precip_snow, returnflow, reinfiltration, &
            & irrigation, tot_melt, vevapwet, transpir, vevapnu, vevapsno, vevapflo, floodout, runoff, drainage)
    ENDIF

    !! 4 write out file  ==> hydrol_alma/histwrite(*)
    !
    ! If we use the ALMA standards
    IF (almaoutput) THEN
       CALL hydrol_alma(kjpindex, index, .FALSE., qsintveg, snow, snow_nobio, soilwet)
    ENDIF


!    write(*,*) 'hydrol 510: almaoutput'
    IF ( .NOT. almaoutput ) THEN
       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
          WRITE (var_name,"('moistc_',i1)") jst
          CALL histwrite(hist_id, trim(var_name), kjit,mc(:,:,jst), kjpindex*nslm, indexlayer)

          ! var_name= "kfactroot_1" ... "kfactroot_3"
          WRITE (var_name,"('kfactroot_',i1)") jst
          CALL histwrite(hist_id, trim(var_name), kjit, kfact_root(:,:,jst), kjpindex*nslm, indexlayer)

          ! var_name= "vegetsoil_1" ... "vegetsoil_3"
          WRITE (var_name,"('vegetsoil_',i1)") jst
          CALL histwrite(hist_id, trim(var_name), kjit,corr_veg_soil(:,:,jst), kjpindex*nvm, indexveg)
       ENDDO


       CALL histwrite(hist_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
       CALL histwrite(hist_id, 'stagnant', kjit, stagnant, kjpindex*nstm, indexsoil)
       CALL histwrite(hist_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite(hist_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite(hist_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
       CALL histwrite(hist_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
       CALL histwrite(hist_id, 'humtot', kjit, humtot, kjpindex, index)
         !Chloe:
       CALL histwrite(hist_id, 'mcs_njsc', kjit, mcs_njsc, kjpindex, index)
        CALL histwrite(hist_id, 'water2add_peat', kjit, water2add_peat, kjpindex, index)  
       CALL histwrite(hist_id, 'wt_soil', kjit, wt_soil, kjpindex*nstm,indexsoil)
       CALL histwrite(hist_id, 'wt_soil2', kjit, wt_soil2, kjpindex*nstm,indexsoil)
        njsc_tmp(:)=njsc(:)
       CALL histwrite(hist_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
       CALL histwrite(hist_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
       CALL histwrite(hist_id, 'drainage', kjit, drainage, kjpindex, index)
       CALL histwrite(hist_id, 'runoff', kjit, runoff, kjpindex, index)
       CALL histwrite(hist_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
       CALL histwrite(hist_id, 'rain', kjit, precip_rain, kjpindex, index)
       CALL histwrite(hist_id, 'snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite(hist_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
       CALL histwrite(hist_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
       CALL histwrite(hist_id, 'SWI', kjit, swi, kjpindex, index)
!Isa
       CALL histwrite(hist_id, 'tot_melt', kjit, tot_melt, kjpindex, index)
       !
       IF ( control_in%do_floodplains ) THEN
          CALL histwrite(hist_id, 'floodout', kjit, floodout, kjpindex, index)
       ENDIF
       !
       IF ( hist2_id > 0 ) THEN
          DO jst=1,nstm
             ! var_name= "mc_1" ... "mc_3"
             WRITE (var_name,"('moistc_',i1)") jst
             CALL histwrite(hist2_id, trim(var_name), kjit,mc(:,:,jst), kjpindex*nslm, indexlayer)

             ! var_name= "kfactroot_1" ... "kfactroot_3"
             WRITE (var_name,"('kfactroot_',i1)") jst
             CALL histwrite(hist2_id, trim(var_name), kjit, kfact_root(:,:,jst), kjpindex*nslm, indexlayer)

             ! var_name= "vegetsoil_1" ... "vegetsoil_3"
             WRITE (var_name,"('vegetsoil_',i1)") jst
             CALL histwrite(hist2_id, trim(var_name), kjit,corr_veg_soil(:,:,jst), kjpindex*nvm, indexveg)
          ENDDO
          CALL histwrite(hist2_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
          CALL histwrite(hist2_id, 'stagnant', kjit, stagnant, kjpindex*nstm, indexsoil) 
          CALL histwrite(hist2_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite(hist2_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite(hist2_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
          CALL histwrite(hist2_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
          CALL histwrite(hist2_id, 'humtot', kjit, humtot, kjpindex, index)
          !Chloe++
          CALL histwrite(hist2_id, 'wt_soil', kjit, wt_soil, kjpindex*nstm,indexsoil)
          CALL histwrite(hist2_id, 'wt_soil2', kjit, wt_soil2, kjpindex*nstm,indexsoil)  
        
          CALL histwrite(hist2_id, 'water2add_peat', kjit, water2add_peat, kjpindex, index)
          !Chloe-- 
            njsc_tmp(:)=njsc(:)
          CALL histwrite(hist2_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
          CALL histwrite(hist2_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
          CALL histwrite(hist2_id, 'drainage', kjit, drainage, kjpindex, index)
          CALL histwrite(hist2_id, 'runoff', kjit, runoff, kjpindex, index)
          IF ( control_in%do_floodplains ) THEN
             CALL histwrite(hist2_id, 'floodout', kjit, floodout, kjpindex, index)
          ENDIF
          CALL histwrite(hist2_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
          CALL histwrite(hist2_id, 'rain', kjit, precip_rain, kjpindex, index)
          CALL histwrite(hist2_id, 'snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite(hist2_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
          CALL histwrite(hist2_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
          CALL histwrite(hist2_id, 'SWI', kjit, swi, kjpindex, index) 
          !
          IF (check_waterbal) THEN
             CALL histwrite(hist2_id, 'TotWater', kjit, tot_water_end, kjpindex, index)
             CALL histwrite(hist2_id, 'TotWaterFlux', kjit, tot_flux, kjpindex, index)
          ENDIF
       ENDIF
    ELSE
       CALL histwrite(hist_id, 'Snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite(hist_id, 'Rainf', kjit, precip_rain, kjpindex, index)
       CALL histwrite(hist_id, 'Qs', kjit, runoff, kjpindex, index)
       CALL histwrite(hist_id, 'Qsb', kjit, drainage, kjpindex, index)
       CALL histwrite(hist_id, 'Qsm', kjit, tot_melt, kjpindex, index)
       CALL histwrite(hist_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
       CALL histwrite(hist_id, 'DelSWE', kjit, delswe, kjpindex, index)
       CALL histwrite(hist_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
       !
       CALL histwrite(hist_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
       CALL histwrite(hist_id, 'SoilWet', kjit, soilwet, kjpindex, index)
       !
       CALL histwrite(hist_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
       CALL histwrite(hist_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
       !
       CALL histwrite(hist_id, 'SnowDepth', kjit, snowdepth, kjpindex, index)
       !
       IF ( hist2_id > 0 ) THEN
          CALL histwrite(hist2_id, 'Snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite(hist2_id, 'Rainf', kjit, precip_rain, kjpindex, index)
          CALL histwrite(hist2_id, 'Qs', kjit, runoff, kjpindex, index)
          CALL histwrite(hist2_id, 'Qsb', kjit, drainage, kjpindex, index)
          CALL histwrite(hist2_id, 'Qsm', kjit, tot_melt, kjpindex, index)
          CALL histwrite(hist2_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
          CALL histwrite(hist2_id, 'DelSWE', kjit, delswe, kjpindex, index)
          CALL histwrite(hist2_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
          !
          CALL histwrite(hist2_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
          CALL histwrite(hist2_id, 'SoilWet', kjit, soilwet, kjpindex, index)
          !
          CALL histwrite(hist2_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
          CALL histwrite(hist2_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
          !
          CALL histwrite(hist2_id, 'SnowDepth', kjit, snowdepth, kjpindex, index)
       ENDIF
    ENDIF
!    write(*,*) 'hydrol 626'
!Isa
     if (ok_freeze_cwrr) then
       CALL histwrite(hist_id, 'profil_froz_hydro', kjit,profil_froz_hydro , kjpindex*nslm, indexlayer)
       DO jst=1,nstm
            WRITE (var_name,"('profil_froz_hydro_',i1)") jst
            CALL histwrite(hist_id, trim(var_name), kjit, profil_froz_hydro_ns(:,:,jst), kjpindex*nslm, indexlayer)
       ENDDO
       CALL histwrite(hist_id, 'temp_hydro', kjit,temp_hydro , kjpindex*nslm, indexlayer)
! Isa pour �criture des k_lin, � faire une fois !
!caveats : ce n'est pas la fa�on la plus rigoureuse de faire...(cf le dernier kjpindex qui devrait �tre kjpij)
!Chloe use k_lin_write d Isa, garder kk_moy non comment�, Comment Debut
!DO ji=1,kjpindex
!	do jst = 1, nscm
!		do i=1, imax
!				mc_lin_axis_index((i-1)*kjpindex+ji) = INDEX(ji) + (i-1)*kjpindex
!				k_lin_write(:,i, jsc) = k_lin(i,:,jsc)
!		enddo
!	enddo
!ENDDO

       !CALL histwrite(hist_id, 'k_lin1', kjit, k_lin_write(:,:,1),kjpindex*imax, mc_lin_axis_index )
       !CALL histwrite(hist_id, 'k_lin2', kjit, k_lin_write(:,:,2),kjpindex*imax, mc_lin_axis_index)
       !CALL histwrite(hist_id, 'k_lin3', kjit, k_lin_write(:,:,3),kjpindex*imax, mc_lin_axis_index)
       CALL histwrite(hist_id, 'kk_moy', kjit, kk_moy,kjpindex*nslm, indexlayer) ! averaged over soiltiles
       CALL histwrite(hist_id, 'kk1', kjit, kk(:,:,1),kjpindex*nslm, indexlayer)
       CALL histwrite(hist_id, 'kk2', kjit, kk(:,:,2),kjpindex*nslm, indexlayer)
       CALL histwrite(hist_id, 'kk3', kjit, kk(:,:,3),kjpindex*nslm, indexlayer)
!end Isa, Chloe, commenter cette zone Comment Fin
     endif ! ok_freeze_cwrr


!TEST A EFFACER
! DO ji=1,kjpindex
!        DO jst=1,nstm
!        tmc(ji,nstm)=600
!        ENDDO
! ENDDO

!DO ji=1,kjpindex
!IF (stagnant(ji,4) .NE. -1.) THEN
!STOP 'stagnant ne -1'
!ENDIF
!ENDDO


    IF (l_second_hydrol) THEN 
       l_second_hydrol=.FALSE.
    ENDIF
    IF (long_print) WRITE (numout,*) ' hydrol_main Done '


  END SUBROUTINE hydrol_main


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_init
!!
!>\BRIEF        Initializations and memory allocation   
!!
!! DESCRIPTION  :
!! - 1 Some initializations
!! - 2 make dynamic allocation with good dimension
!! - 2.1 array allocation for soil texture
!! - 2.2 Soil texture choice
!! - 3 Other array allocation
!! - 4 Open restart input file and read data for HYDROLOGIC process
!! - 5 get restart values if none were found in the restart file
!! - 6 Vegetation array      
!! - 7 set humrelv from us
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!!_ hydrol_init

  SUBROUTINE hydrol_init(kjit, ldrestart_read, kjpindex, index, rest_id, veget_max, soiltile, humrel,&
       &  njsc, vegstress, snow, snow_age, snow_nobio, snow_nobio_age, qsintveg, peatland, ok_routage_peat, wt_soil,wt_soil2)

    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number 
    LOGICAL,INTENT (in)                                 :: ldrestart_read     !! Logical for _restart_ file to read
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier 
    !Chloe+ 131113
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc               !! indexing of PFT to soiltile
    LOGICAL, INTENT(in), DIMENSION(kjpindex)            :: peatland           !! logical peatland if dominant in a grid cell
    LOGICAL, INTENT(in)                                 :: ok_routage_peat     
    
    ! input fields
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max          !! Carte de vegetation max
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)  :: soiltile           !! fraction of PFT on each soil-hydrology tile

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: humrel             !! Stress hydrique, relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: vegstress          !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: snow               !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: snow_age           !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: snow_nobio       !! Snow on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: qsintveg           !! Water on vegetation due to interception
!Chloe
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT(out)  :: wt_soil            !!Water Table position (mm) 
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT(out)  :: wt_soil2            !!Water Table position (mm) 
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ier
    INTEGER(i_std)                                      :: ji, jv, jst, jsl, jsc

    ! initialisation
    IF (l_first_hydrol) THEN 
       l_first_hydrol=.FALSE.
    ELSE 
       WRITE (numout,*) ' l_first_hydrol false . we stop '
       STOP 'hydrol_init'
    ENDIF

    !
    !! 1 Some initializations
    !
    !
    !Config Key   = CHECK_CWRR
    !Config Desc  = Should we check detailed CWRR water balance ?
    !Config Def   = FALSE
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to check
    !Config         the detailed water balance in each time step
    !Config         of CWRR.
    !Config Units = [FLAG]
    !
    check_cwrr = .FALSE.
    CALL getin_p('CHECK_CWRR', check_cwrr)
    !
    !Config Key   = DO_PONDS
    !Config Desc  = Should we include ponds 
    !Config Def   = FALSE
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the ponds and return 
    !Config         the water into the soil moisture. If this is 
    !Config         activated, then there is no reinfiltration 
    !Config         computed inside the hydrol module.
    !Config Units = [FLAG]
    !
    doponds = .FALSE.
    CALL getin_p('DO_PONDS', doponds)

    ! Chloe050813
    max_stagnant = 100.
    CALL getin_p('MAX_STAGNANT', max_stagnant)

     
!Isa
    if ( control%ok_freeze ) then
        write(*,*) 'hydrol 768: ok_freeze active, on active toutes les options d Isabelle'
        ! on active tout par d�faut:
        ok_freeze_cwrr = .true.
        ok_gel_thd = .FALSE.
        ok_gel_thermosoil = .true.
        ok_shumdiag_interpol =.true.
        ok_shumdiag_perma =.true.

    else !if (ok_freeze) then
        write(*,*) 'hydrol 776: ok_freeze desactive, on lit les options 1 par 1'
        ok_freeze_cwrr = .FALSE.
        CALL getin_p('OK_FREEZE_CWRR',ok_freeze_cwrr)
        !calcul de la fraction gelee
        ok_gel_thd = .FALSE.
        CALL getin_p('OK_GEL_THD',ok_gel_thd)
        if (ok_gel_thd) write(*,*) 'ISA : ok_gel_thd'
        ok_gel_thermosoil = .FALSE.
        CALL getin_p('OK_GEL_THERMOSOIL',ok_gel_thermosoil)
        if (ok_gel_thermosoil) write(*,*) 'ISA : ok_gel_thermosoil'
        ok_shumdiag_perma = .FALSE.
        CALL getin_p('OK_SHUMDIAG_PERMA',  ok_shumdiag_perma)
        if ( ok_shumdiag_perma) write(*,*) 'ISA : ok_shumdiag_perma'
        ok_shumdiag_interpol = .FALSE.
        CALL getin_p('OK_SHUMDIAG_INTERPOL',  ok_shumdiag_interpol)
        write(*,*) 'ok_shumdiag_interpol=',ok_shumdiag_interpol
  
    endif !if (ok_freeze) then
     
        !Chloe++

        ! d�fini dans sechiba :
        !ok_routage_peat = .FALSE.
        !CALL getin_p('OK_ROUTAGE_PEAT', ok_routage_peat)
        !if (ok_routage_peat) write(*,*) 'Chlo� : ok_routage_peat'

        ok_sat30cm = .FALSE.
        CALL getin_p('OK_SAT30CM', ok_sat30cm)
        if (ok_sat30cm) write(*,*) 'Chlo� : ok_sat30cm'
        
        ok_stagnant= .FALSE.
        CALL getin_p('OK_STAGNANT', ok_stagnant)
        if (ok_stagnant) write(*,*) 'Chlo� : ok_stagnant'

         !d�fini dans sechiba : 
        ! ok_no_drainage= .FALSE.
        !CALL getin_p('OK_NO_DRAINAGE', ok_no_drainage)
        !if (ok_no_drainage) write(*,*) 'Chlo� : ok_no_drainage'

       ok_reinfilt_peat= .FALSE.
       CALL getin_p('OK_REINFILT_PEAT',  ok_reinfilt_peat)
        if ( ok_reinfilt_peat) write(*,*) 'Chlo� : ok_reinfilt_peat'


        !Chloe--



    ! CR: verif compatibilit� des options:
    if (ok_freeze_cwrr) then
        if (ok_gel_thd) then
            if (ok_gel_thermosoil) then
              write(*,*) 'hydrol 771: options gel sol incompatibles'
              write(*,*) 'ok_gel_thd,ok_gel_thermosoil=',ok_gel_thd,ok_gel_thermosoil
              write(*,*) '1 seule doit etre activee: on stoppe'
              stop
            endif !if (ok_gel_thermosoil) then
        elseif (ok_gel_thermosoil) then
            if (ok_gel_thd) then
              write(*,*) 'hydrol 781: options gel sol incompatibles'
              write(*,*) 'ok_gel_thd,ok_gel_thermosoil=',ok_gel_thd,ok_gel_thermosoil
              write(*,*) '1 seule doit etre activee: on stoppe'
              stop
            endif !if (ok_gel_thermosoil) then
        endif !if (ok_gel_thd) then
    endif !if (ok_freeze_cwrr) then  
    if ( ok_shumdiag_perma) then
        if (ok_shumdiag_interpol) then
        else
            write(*,*) 'hydrol 801: options shumdiag incompatibles'
            write(*,*) 'ok_shumdiag_perma ne marche que si ok_shumdiag_interpol'
            stop
        endif
    endif  

    !! 2 make dynamic allocation with good dimension

    !! 2.1 array allocation for soil texture

    ALLOCATE (nvan(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in nvan allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (avan(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in avan allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcr(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcr allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcs(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcs allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ks(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ks allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ds(nslm,nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ds allocation. We stop. We need nslm*nscm words = ',nslm*nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (pcent(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in pcent allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (free_drain_max(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in free_drain_max allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcf(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcf allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcw(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcw allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF
    
    ALLOCATE (mc_awet(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_awet allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mc_adry(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_adry allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF
       
    ALLOCATE (psis(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in psis allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF
    !!__2.2 Soil texture choose
    SELECTCASE (nscm)
    CASE (3)
              
       nvan = nvan_fao 
       !nvan = nvan_peat  
       avan = avan_fao
       !avan=avan_peat
       mcr = mcr_fao
       !mcr=mcr_peat
       mcs = mcs_fao
       !mcs=mcs_peat
       ks = ks_fao
       !ks=ks_peat
       ds(:,:)=zero
       pcent = pcent_fao
       free_drain_max = free_drain_max_fao
       mcf = mcf_fao
       mcw = mcw_fao
       mc_awet = mc_awet_fao
       mc_adry = mc_adry_fao
       psis = psis_fao
       
        DO ji=1, kjpindex
            IF (peatland(ji) ) THEN 
             !  nvan = nvan_peat       
             !  avan = avan_peat
               mcs = mcs_peat
             !  ks = ks_peat
             ! mcr= mcr_peat
            ENDIF
            
        ENDDO

       
    CASE (5)
              
       DO jsc =1, nscm 
          nvan(jsc) = nvan_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          avan(jsc) = avan_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          mcr(jsc) = mcr_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          mcs(jsc) = mcs_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          ks(jsc) = ks_fao(CEILING(jsc*un/deux))
       ENDDO
       ds(:,:)=zero
       DO jsc =1, nscm 
          pcent(jsc) = pcent_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          free_drain_max(jsc) = free_drain_max_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          mcf(jsc) = mcf_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          mcw(jsc) = mcw_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          mc_awet(jsc) = mc_awet_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          mc_adry(jsc) = mc_adry_fao(CEILING(jsc*un/deux))
       ENDDO
       DO jsc =1, nscm 
          psis(jsc) = psis_fao(CEILING(jsc*un/deux))
       ENDDO
       
    CASE (12)
       
       nvan = nvan_usda
       avan = avan_usda
       mcr = mcr_usda
       mcs = mcs_usda
       ks = ks_usda
       ds(:,:)=zero
       pcent = pcent_usda
       free_drain_max = free_drain_max_usda
       mcf = mcf_usda
       mcw = mcw_usda
       mc_awet = mc_awet_usda
       mc_adry = mc_adry_usda
       psis = psis_usda
       
    CASE DEFAULT
       WRITE (numout,*) 'Unsupported soil type classification. Choose between zobler, fao and usda according to the map'
       STOP 'hydrol_init'
    ENDSELECT

    !! 3 Other array allocation


    ALLOCATE (mask_veget(kjpindex,nvm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mask_veget allocation. We stop. We need kjpindex*nvm words = ',kjpindex*nvm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mask_soiltile(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mask_soiltile allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (humrelv(kjpindex,nvm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in humrelv allocation. We stop. We need kjpindex words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (vegstressv(kjpindex,nvm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vegstressv allocation. We stop. We need kjpindex words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (us(kjpindex,nvm,nstm,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in us allocation. We stop. We need kjpindex words = ',kjpindex*nvm*nstm*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (precisol(kjpindex,nvm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in precisol allocation. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'hydrol_init'
    END IF

!Chloe : d�plac� dans sechiba
!    ALLOCATE (wt_soil(kjpindex,nstm),stat=ier)
!    IF (ier.NE.0) THEN
!       WRITE (numout,*) ' error in wt_soil allocation. We stop. We need kjpindex words = ',kjpindex*nstm
!       STOP 'hydrol_init'
!    END IF


    ALLOCATE (precisol_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in precisol_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (free_drain_coef(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in free_drain_coef allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (frac_bare_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_bare_ns allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (water2infilt(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in water2infilt allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ae_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ae_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ! Reservoir d'eau stagnante
    ALLOCATE (stagnant(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in stagnant allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    ENDIF  
 
    ALLOCATE (evap_bare_lim_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in evap_bare_lim_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (rootsink(kjpindex,nslm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rootsink allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (subsnowveg(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in subsnowveg allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (subsnownobio(kjpindex,nnobio),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in subsnownobio allocation. We stop. We need kjpindex words = ',kjpindex*nnobio
       STOP 'hydrol_init'
    END IF

    ALLOCATE (snowmelt(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowmelt allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (icemelt(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in icemelt allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (subsinksoil(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in subsinksoil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mx_eau_var(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mx_eau_var allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (vegtot(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vegtot allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (resdist(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in resdist allocation. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (vmr(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vmr allocation. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (humtot(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in humtot allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

!Chloe
    ALLOCATE (mcs_njsc(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcs_njsc allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (flux(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in flux allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (resolv(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in resolv allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (k(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in k allocation. We stop. We need kjpindex*nslm words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF
    
!Isa
 if (ok_freeze_cwrr) then
    ALLOCATE (kk_moy(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in kk allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF
    ALLOCATE (kk(kjpindex,nslm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in kk allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF
 endif !if (ok_freeze_cwrr) then

    ALLOCATE (a(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in a allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (b(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in b allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (d(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in d allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (e(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in e allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (f(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in f allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (g1(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in g1 allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ep(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ep allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (fp(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in fp allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (gp(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gp allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (rhs(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rhs allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (srhs(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in srhs allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (gam(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gam allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    !Chloe
    ALLOCATE (tmc_obj_peat(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN 
       WRITE (numout,*) ' error in tmc_obj_peat allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF


    ALLOCATE (tmcs(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmcs allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmcr(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmcr allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litt_mea(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litt_mea allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_res(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_res allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_wilt(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_wilt allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_field(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_field allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_sat(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_sat allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_awet(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_awet allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_adry(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_adry allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litt_wet_mea(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litt_wet_mea allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litt_dry_mea(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litt_dry_mea allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (v1(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in v1 allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (qflux00(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qflux00 allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ru_ns(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ru_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF
    ru_ns(:,:) = zero

    ALLOCATE (dr_ns(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in dr_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF
    dr_ns(:,:) = zero

    ALLOCATE (tr_ns(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tr_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (cvs_over_veg(kjpindex,nvm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in cvs_over_veg allocation. We stop. We need kjpindex*nvm*nstm words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (corr_veg_soil(kjpindex,nvm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in corr_veg_soil allocation. We stop. We need kjpindex*nvm*nstm words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF
!Chloe
    ALLOCATE (mc_obj_peat(kjpindex,nslm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_obj_peat allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mc(kjpindex,nslm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (soilmoist(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soilmoist allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (soil_wet(kjpindex,nslm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soil_wet allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (soil_wet_litter(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soil_wet allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (qflux(kjpindex,nslm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qflux allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmat(kjpindex,nslm,3),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmat allocation. We stop. We need kjpindex words = ',kjpindex*nslm*trois
       STOP 'hydrol_init'
    END IF

    ALLOCATE (stmat(kjpindex,nslm,3),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in stmat allocation. We stop. We need kjpindex words = ',kjpindex*nslm*trois
       STOP 'hydrol_init'
    END IF

    ALLOCATE (nroot(nvm, nstm, nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in nroot allocation. We stop. We need nvm*nstm*nslm words = ',nvm * nstm * nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (kfact_root(kjpindex, nslm, nstm), stat=ier)
    IF (ier .NE. 0) THEN
       WRITE (numout,*) 'error in kfact_root allocation, We stop. We need kjpindex*nslm*nstm words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (kfact(nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in kfact allocation. We stop. We need nslm*nscm words = ',nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (zz(nslm+1, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in zz allocation. We stop. We need (nslm+1)*nstm words = ',(nslm+1) * nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (dz(nslm+1, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in dz allocation. We stop. We need (nslm+1)*nstm words = ',(nslm+1) * nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mc_lin(imin:imax, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_lin allocation. We stop. We need (imax-imin)*nscm words = ',(imax-imin) * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (k_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in k_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF
!Isa

    ALLOCATE (k_lin_write(kjpindex,imin:imax, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in k_lin_write allocation. We stop. We need kjpindex words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF
!end Isa

    ALLOCATE (d_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in d_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (a_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in a_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (b_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in b_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

    !Isa
  if (ok_freeze_cwrr) then
    ALLOCATE (profil_froz_hydro(kjpindex, nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in profil_froz_hydro allocation. We stop. ISA'
       STOP 'hydrol_init'
    END IF
    ALLOCATE (profil_froz_hydro_ns(kjpindex, nslm, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in profil_froz_hydro_ns allocation. We stop. ISA'
       STOP 'hydrol_init'
    END IF
  endif !if (ok_freeze_cwrr) then

    ALLOCATE (temp_hydro(kjpindex, nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in temp_hydro allocation. We stop. ISA'
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcl(kjpindex, nslm, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcl allocation. We stop. ISA'
       STOP 'hydrol_init'
    END IF
    ALLOCATE (frac_hydro_diag(nslm, nbdl),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_hydro_diag allocation. We stop'
       STOP 'calcule_frac_hydro_diag'
    END IF
!end Isa

    !  If we check the water balance we need two more variables
    IF ( check_waterbal ) THEN

       ALLOCATE (tot_water_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_water_beg allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_water_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_water_end allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_flux(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_flux allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

    ENDIF

    ! Soil Wetness Index
    ALLOCATE (swi(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in swi. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    ENDIF  

    !
    !  If we use the almaoutputs we need a few more variables
    !  tdo - they could be allocated only if alma_output, but then they should also be computed only if alma_output
    !
    IF ( almaoutput ) THEN

       ALLOCATE (tot_watveg_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watveg_beg allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_watveg_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watveg_end allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_watsoil_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watsoil_beg allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_watsoil_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watsoil_end allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (delsoilmoist(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in delsoilmoist allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (delintercept(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in delintercept. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (delswe(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in delswe. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       ENDIF

       ALLOCATE (snow_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in snow_beg allocation. We stop. We need kjpindex words =',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (snow_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in snow_end allocation. We stop. We need kjpindex words =',kjpindex
          STOP 'hydrol_init'
       END IF


    ENDIF !IF ( almaoutput ) THEN

    !! 4 Open restart input file and read data for HYDROLOGIC process

    IF (ldrestart_read) THEN

       IF (long_print) WRITE (numout,*) ' we have to read a restart file for HYDROLOGIC variables'

       IF (is_root_prc) CALL ioconf_setatt('UNITS', '-')
       !
       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
           WRITE (var_name,"('moistc_',I1)") jst
           IF (is_root_prc) CALL ioconf_setatt('LONG_NAME',var_name)
           CALL restget_p (rest_id, var_name, nbp_glo, nslm , 1, kjit, .TRUE., mc(:,:,jst), "gather", nbp_glo, index_g)
       END DO
       !
       IF (is_root_prc) CALL ioconf_setatt('UNITS', '-')
       write(*,*) 'hydrol 1541: restget_p us_'
       DO jst=1,nstm
          DO jsl=1,nslm
             ! var_name= "us_1_01" ... "us_3_11"
             WRITE (var_name,"('us_',i1,'_',i2.2)") jst,jsl
             IF (is_root_prc) CALL ioconf_setatt('LONG_NAME',var_name)
             CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., us(:,:,jst,jsl), "gather", nbp_glo, index_g)
          END DO
       END DO
!       write(*,*) 'hydrol 1551'
       !
       var_name= 'free_drain_coef'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','Coefficient for free drainage at bottom of soil')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., free_drain_coef, "gather", nbp_glo, index_g)
       !!Chloe++
       var_name= 'wt_soil'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'mm')
          CALL ioconf_setatt('LONG_NAME','Water table')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., wt_soil, "gather", nbp_glo, index_g)
       !
       var_name= 'wt_soil2'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'mm')
          CALL ioconf_setatt('LONG_NAME','Water table')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., wt_soil2, "gather", nbp_glo, index_g)
       !!Chloe--
       !
       var_name= 'water2infilt'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','Remaining water to be infiltrated on top of the soil')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., water2infilt, "gather", nbp_glo, index_g)
       !
       var_name= 'ae_ns'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'kg/m^2')
          CALL ioconf_setatt('LONG_NAME','Bare soil evap on each soil type')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., ae_ns, "gather", nbp_glo, index_g)
       !
       var_name= 'stagnant'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'kg/m^2')
          CALL ioconf_setatt('LONG_NAME','Stagnant water')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., stagnant, "gather", nbp_glo, index_g)
       !
       var_name= 'snow'        
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'kg/m^2')
          CALL ioconf_setatt('LONG_NAME','Snow mass')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_age'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'd')
          CALL ioconf_setatt('LONG_NAME','Snow age')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow_age, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_nobio'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'kg/m^2')
          CALL ioconf_setatt('LONG_NAME','Snow on other surface types')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_nobio_age'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'd')
          CALL ioconf_setatt('LONG_NAME','Snow age on other surface types')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio_age, "gather", nbp_glo, index_g)
       !
       var_name= 'vegstress'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','Vegetation growth moisture stress')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., vegstress, "gather", nbp_glo, index_g)
       !
       var_name= 'qsintveg'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'kg/m^2')
          CALL ioconf_setatt('LONG_NAME','Intercepted moisture')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., qsintveg, "gather", nbp_glo, index_g)
       !
       var_name= 'resdist'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','Distribution of reservoirs')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., resdist, "gather", nbp_glo, index_g)
       !
       IF (is_root_prc) CALL ioconf_setatt('UNITS', '-')
       DO jst=1,nstm
          ! var_name= "cvs_over_veg_1" ... "cvs_over_veg_3"
          WRITE (var_name,"('cvs_over_veg_',i1)") jst
          IF (is_root_prc) CALL ioconf_setatt('LONG_NAME',var_name)
          CALL restget_p (rest_id, var_name, nbp_glo,  nvm, 1, kjit, .TRUE., cvs_over_veg(:,:,jst), "gather",  nbp_glo, index_g)
       END DO
       !
       IF ( check_waterbal ) THEN
          var_name= 'tot_water_beg'
          IF (is_root_prc) THEN
             CALL ioconf_setatt('UNITS', 'kg/m^2')
             CALL ioconf_setatt('LONG_NAME','Previous Total water')
          ENDIF
          CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., tot_water_beg, "gather", nbp_glo, index_g)
       ENDIF


    !! 5 get restart values if none were found in the restart file
       !
       !Config Key   = HYDROL_MOISTURE_CONTENT
       !Config Desc  = Soil moisture on each soil tile and levels
       !Config If    = HYDROL_CWRR       
       !Config Def   = 0.3
       !Config Help  = The initial value of mc if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       write(*,*) 'hydrol 1681: ok_converge_isaorig=',control%ok_converge_isaorig
       if (control%ok_converge_isaorig) then
          CALL setvar_p (mc, val_exp, 'HYDROL_MOISTURE_CONTENT', 0.33_r_std)
       else !if (ok_converge_isaorig) then 
          CALL setvar_p (mc, val_exp, 'HYDROL_MOISTURE_CONTENT', 0.3_r_std)
       endif!if (ok_converge_isaorig) then
       mcl(:,:,:)=mc(:,:,:)
       write(*,*) 'hydrol 1682'
!Chloe
CALL setvar_p (mc_obj_peat, val_exp, 'HYDROL_MOISTURE_CONTENT_OBJ_PEAT', 0.3_r_std)

       !
       !Config Key   = US_INIT
       !Config Desc  = US_NVM_NSTM_NSLM
       !Config If    = HYDROL_CWRR       
       !Config Def   = 0.0
       !Config Help  = The initial value of us (relative moisture) if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =

!      !Chloe++
        IF(ok_sat30cm) THEN  
            DO ji=1, kjpindex 
                DO jsl=9,nslm
                mc(:,jsl,4)=mcs(njsc(ji)) ! mcs(4)
                 IF (peatland(ji) ) mc(:,jsl,:)=mcs_peat !Chloe
                ENDDO
            ENDDO
        ENDIF

        IF(ok_routage_peat) THEN
            DO ji=1, kjpindex
                DO jsl=1,7
                mc_obj_peat(:,jsl,4)=mc(:,jsl,4)
                ENDDO

                DO jsl=8,nslm
                    mc_obj_peat(:,jsl,4)=mcs(njsc(ji)) ! mcs(4)
                    IF (peatland(ji) ) mc(:,jsl,:)=mcs_peat !Chloe
                ENDDO
            ENDDO
        ENDIF

                !Test Desert mousson :
    !    DO ji=1, kjpindex
    !        DO jsl=1,4
    !        mc(ji,jsl,:)=mcs(njsc(ji))
    !        ENDDO
    !    ENDDO


!     !Chloe--
 

       DO jsl=1,nslm
          CALL setvar_p (us(:,:,:,jsl), val_exp, 'US_INIT', zero)
       ENDDO

       !Config Key  = WT_SOIL
       !Config Desc = Coefficient for free drainage at bottom
       !Config Def  = -1.0, -1.0, -1.0
       !Config Help = The initial value of wt impose
       !Config        in the restart file. This should only be used if the model
       !is 
       !Config        started without a restart file.
       !
       CALL setvar_p (wt_soil, val_exp, 'WT_SOIL', (/ -1.0_r_std, -1.0_r_std,-1.0_r_std, -1.0_r_std /) )
       CALL setvar_p (wt_soil2, val_exp, 'WT_SOIL2', (/ 0.0_r_std, 0.0_r_std,0.0_r_std, 0.0_r_std /) )
       !
       !Config Key   = FREE_DRAIN_COEF
       !Config Desc  = Coefficient for free drainage at bottom
       !Config If    = HYDROL_CWRR       
       !Config Def   = 1.0, 1.0, 1.0
       !Config Help  = The initial value of free drainage if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       !Chloe+ CL081113 ce call ne marche pas, je remet l'ancien (ie NO_KEYWORD)
       ! free_drain_coef(ji,jst) free_drain_max(nscm) or nstm=4 et nscm=3 ne marche plus
       !CALL setvar_p (free_drain_coef, val_exp, 'FREE_DRAIN_COEF', free_drain_max)

       !write(*,*) 'hydrol CLtest 3'

       CALL setvar_p (free_drain_coef, val_exp, 'NO_KEYWORD', 1.0_r_std)

       !Config Key   = WATER_TO_INFILT
       !Config Desc  = Water to be infiltrated on top of the soil
       !Config If    = HYDROL_CWRR    
       !Config Def   = 0.0
       !Config Help  = The initial value of free drainage if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (water2infilt, val_exp, 'WATER_TO_INFILT', zero)
       !
       !Config Key   = EVAPNU_SOIL
       !Config Desc  = Bare soil evap on each soil if not found in restart
       !Config If    = HYDROL_CWRR  
       !Config Def   = 0.0
       !Config Help  = The initial value of bare soils evap if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (ae_ns, val_exp, 'EVAPNU_SOIL', zero)
       !
       !Config Key   = STAGNANT
       !Config Desc  = Stagnant water reservoir
       !Config If    = HYDROL_CWRR
       !Config Def   = 0.0
       !Config Help  = The initial value of stagnant water reservoir if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (stagnant, val_exp, 'STAGNANT', zero)
       !
       !Config Key  = HYDROL_SNOW
       !Config Desc  = Initial snow mass if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow mass if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow, val_exp, 'HYDROL_SNOW', zero)
       !
       !Config Key   = HYDROL_SNOWAGE
       !Config Desc  = Initial snow age if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow age if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow_age, val_exp, 'HYDROL_SNOWAGE', zero)
       !
       !Config Key   = HYDROL_SNOW_NOBIO
       !Config Desc  = Initial snow amount on ice, lakes, etc. if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow_nobio, val_exp, 'HYDROL_SNOW_NOBIO', zero)
       !
       !Config Key   = HYDROL_SNOW_NOBIO_AGE
       !Config Desc  = Initial snow age on ice, lakes, etc. if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow age if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow_nobio_age, val_exp, 'HYDROL_SNOW_NOBIO_AGE', zero)
       !
       !Config Key   = HYDROL_QSV
       !Config Desc  = Initial water on canopy if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of moisture on canopy if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units =
       !
       write(*,*) 'hydrol 1756: setvar_p HYDROL_QSV'
       CALL setvar_p (qsintveg, val_exp, 'HYDROL_QSV', zero)

!Isa
     if (ok_freeze_cwrr) then  
       write(*,*) 'hydrol 1760: setvar_p profil_froz_hydro'
       CALL setvar_p (profil_froz_hydro, val_exp, 'NO_KEYWORD', zero)
       write(*,*) 'hydrol 1762: setvar_p profil_froz_hydro_ns'
       CALL setvar_p (profil_froz_hydro_ns, val_exp, 'NO_KEYWORD', zero)
       write(*,*) 'hydrol 1765: setvar_p kk'
       CALL setvar_p (kk, val_exp, 'NO_KEYWORD', 276.48)
       CALL setvar_p (kk_moy, val_exp, 'NO_KEYWORD', 276.48)
       write(*,*) 'hydrol 1794: setvar_p temp_hydro'
       CALL setvar_p (temp_hydro, val_exp, 'NO_KEYWORD', 280.)
     endif !if (ok_freeze_cwrr) then
!     write(*,*) 'hydrol 1798: setvar_p mcl'
!     CALL setvar_p (mcl, val_exp, 'HYDROL_MOISTURE_CONTENT', 0.3_r_std)
       

    !! 6 Vegetation array      
       !
       ! There is no need to configure the initialisation of resdist. If not available it is the vegetation map
       !
       write(*,*) 'hydrol 1772'
       IF ( MINVAL(resdist) .EQ.  MAXVAL(resdist) .AND. MINVAL(resdist) .EQ. val_exp) THEN
          resdist(:,:) = veget_max(:,:)
          vmr(:,:) = zero
       ELSE
          DO jv = 1, nvm
             DO ji = 1, kjpindex
                IF ( ABS(veget_max(ji,jv)-resdist(ji,jv)) .GT. EPS1 ) THEN
                   vmr(ji,jv) = veget_max(ji,jv) - resdist(ji,jv)
                ELSE
                   vmr(ji,jv) = zero
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !
       !  Remember that it is only frac_nobio + SUM(veget_max(,:)) that is equal to 1. Thus we need vegtot
       !
       DO ji = 1, kjpindex
          vegtot(ji) = SUM(veget_max(ji,:))
       ENDDO
       !
       !
       ! compute the masks for veget


       mask_veget(:,:) = 0
       mask_soiltile(:,:) = 0

       DO jst=1,nstm
          DO ji = 1, kjpindex
             IF(soiltile(ji,jst) .GT. min_sechiba) THEN
                mask_soiltile(ji,jst) = 1
             ENDIF
          END DO
       ENDDO
          
       DO jv = 1, nvm
          DO ji = 1, kjpindex
             IF(veget_max(ji,jv) .GT. min_sechiba) THEN
                mask_veget(ji,jv) = 1
             ENDIF
          END DO
       END DO
          
    !! 7 set humrelv from us

       humrelv(:,:,:) = SUM(us,dim=4)
       vegstressv(:,:,:) = humrelv(:,:,:)
       ! set humrel from humrelv, assuming equi-repartition for the first time step

       humrel(:,:) = zero
       CALL setvar_p (cvs_over_veg, val_exp, 'NO_KEYWORD', un)

       DO jst=1,nstm
          DO jv=1,nvm
             DO ji=1,kjpindex

                vegstress(ji,jv)=vegstress(ji,jv) + vegstressv(ji,jv,jst) * &
                     & soiltile(ji,jst) * cvs_over_veg(ji,jv,jst) * vegtot(ji)

                humrel(ji,jv)=humrel(ji,jv) + humrelv(ji,jv,jst) * & 
                     & soiltile(ji,jst) &
                     & * cvs_over_veg(ji,jv,jst)*vegtot(ji)
                humrel(ji,jv)=MAX(humrel(ji,jv), zero)* mask_veget(ji,jv)           
             END DO
          END DO
       END DO
    ENDIF
    write(*,*) 'hydrol 1841'
    !
    !
    IF (long_print) WRITE (numout,*) ' hydrol_init done '
    !

  
  END SUBROUTINE hydrol_init


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_clear
!!
!>\BRIEF        Deallocate arrays 
!!
!_ ================================================================================================================================
!_ hydrol_clear

  SUBROUTINE hydrol_clear()

    l_first_hydrol=.TRUE.

    ! Allocation for soiltile related parameters
    IF ( ALLOCATED (nvan)) DEALLOCATE (nvan)
    IF ( ALLOCATED (avan)) DEALLOCATE (avan)
    IF ( ALLOCATED (mcr)) DEALLOCATE (mcr)
    IF ( ALLOCATED (mcs)) DEALLOCATE (mcs)
    IF ( ALLOCATED (ks)) DEALLOCATE (ks)
    IF ( ALLOCATED (ds)) DEALLOCATE (ds)
    IF ( ALLOCATED (pcent)) DEALLOCATE (pcent)
    IF ( ALLOCATED (free_drain_max)) DEALLOCATE (free_drain_max)
    IF ( ALLOCATED (mcf)) DEALLOCATE (mcf)
    IF ( ALLOCATED (mcw)) DEALLOCATE (mcw)
    IF ( ALLOCATED (mc_awet)) DEALLOCATE (mc_awet)
    IF ( ALLOCATED (mc_adry)) DEALLOCATE (mc_adry)
    IF ( ALLOCATED (psis)) DEALLOCATE (psis)
    ! Other arrays
    IF (ALLOCATED (mask_veget)) DEALLOCATE (mask_veget)
    IF (ALLOCATED (mask_soiltile)) DEALLOCATE (mask_soiltile)
    IF (ALLOCATED (humrelv)) DEALLOCATE (humrelv)
    IF (ALLOCATED (vegstressv)) DEALLOCATE (vegstressv)
    IF (ALLOCATED (us)) DEALLOCATE (us)
    IF (ALLOCATED  (precisol)) DEALLOCATE (precisol)
    IF (ALLOCATED  (precisol_ns)) DEALLOCATE (precisol_ns)
    !IF (ALLOCATED  (wt_soil)) DEALLOCATE (wt_soil)
    IF (ALLOCATED  (free_drain_coef)) DEALLOCATE (free_drain_coef)
    IF (ALLOCATED  (frac_bare_ns)) DEALLOCATE (frac_bare_ns)
    IF (ALLOCATED  (water2infilt)) DEALLOCATE (water2infilt)
    IF (ALLOCATED  (ae_ns)) DEALLOCATE (ae_ns)
    IF (ALLOCATED  (evap_bare_lim_ns)) DEALLOCATE (evap_bare_lim_ns)
    IF (ALLOCATED  (rootsink)) DEALLOCATE (rootsink)
    IF (ALLOCATED  (subsnowveg)) DEALLOCATE (subsnowveg)
    IF (ALLOCATED  (subsnownobio)) DEALLOCATE (subsnownobio)
    IF (ALLOCATED  (snowmelt)) DEALLOCATE (snowmelt)
    IF (ALLOCATED  (icemelt)) DEALLOCATE (icemelt)
    IF (ALLOCATED  (subsinksoil)) DEALLOCATE (subsinksoil)
    IF (ALLOCATED  (mx_eau_var)) DEALLOCATE (mx_eau_var)
    IF (ALLOCATED  (vegtot)) DEALLOCATE (vegtot)
    IF (ALLOCATED  (vmr)) DEALLOCATE (vmr)
    IF (ALLOCATED  (resdist)) DEALLOCATE (resdist)
    IF (ALLOCATED  (tot_water_beg)) DEALLOCATE (tot_water_beg)
    IF (ALLOCATED  (tot_water_end)) DEALLOCATE (tot_water_end)
    IF (ALLOCATED  (tot_flux)) DEALLOCATE (tot_flux)
    IF (ALLOCATED  (tot_watveg_beg)) DEALLOCATE (tot_watveg_beg)
    IF (ALLOCATED  (tot_watveg_end)) DEALLOCATE (tot_watveg_end)
    IF (ALLOCATED  (tot_watsoil_beg)) DEALLOCATE (tot_watsoil_beg)
    IF (ALLOCATED  (tot_watsoil_end)) DEALLOCATE (tot_watsoil_end)
    IF (ALLOCATED  (delsoilmoist)) DEALLOCATE (delsoilmoist)
    IF (ALLOCATED  (delintercept)) DEALLOCATE (delintercept)
    IF (ALLOCATED  (snow_beg)) DEALLOCATE (snow_beg)
    IF (ALLOCATED  (snow_end)) DEALLOCATE (snow_end)
    IF (ALLOCATED  (delswe)) DEALLOCATE (delswe)
    IF (ALLOCATED  (swi)) DEALLOCATE (swi)
    IF (ALLOCATED  (stagnant)) DEALLOCATE (stagnant)
    ! more allocation for cwrr scheme
    IF (ALLOCATED  (v1)) DEALLOCATE (v1)
    IF (ALLOCATED  (humtot)) DEALLOCATE (humtot)
    IF (ALLOCATED  (flux)) DEALLOCATE (flux)
    IF (ALLOCATED  (resolv)) DEALLOCATE (resolv)
    IF (ALLOCATED  (k)) DEALLOCATE (k)
!Chloe
    IF (ALLOCATED  (mcs_njsc)) DEALLOCATE (mcs_njsc)
!Isa
  if (ok_freeze_cwrr) then
    IF (ALLOCATED  (kk)) DEALLOCATE (kk)
    IF (ALLOCATED  (kk_moy)) DEALLOCATE (kk_moy)
  endif !if (ok_freeze_cwrr) then
    IF (ALLOCATED  (a)) DEALLOCATE (a)
    IF (ALLOCATED  (b)) DEALLOCATE (b)
    IF (ALLOCATED  (d)) DEALLOCATE (d)
    IF (ALLOCATED  (e)) DEALLOCATE (e)
    IF (ALLOCATED  (f)) DEALLOCATE (f)
    IF (ALLOCATED  (g1)) DEALLOCATE (g1)
    IF (ALLOCATED  (ep)) DEALLOCATE (ep)
    IF (ALLOCATED  (fp)) DEALLOCATE (fp)
    IF (ALLOCATED  (gp)) DEALLOCATE (gp)
    IF (ALLOCATED  (rhs)) DEALLOCATE (rhs)
    IF (ALLOCATED  (srhs)) DEALLOCATE (srhs)
    IF (ALLOCATED  (gam)) DEALLOCATE (gam)
!Chloe
    IF (ALLOCATED  (tmc_obj_peat)) DEALLOCATE (tmc_obj_peat)
    IF (ALLOCATED  (tmc)) DEALLOCATE (tmc)
    IF (ALLOCATED  (tmcs)) DEALLOCATE (tmcs)
    IF (ALLOCATED  (tmcr)) DEALLOCATE (tmcr)
    IF (ALLOCATED  (tmc_litter)) DEALLOCATE (tmc_litter)
    IF (ALLOCATED  (tmc_litt_mea)) DEALLOCATE (tmc_litt_mea)
    IF (ALLOCATED  (tmc_litter_res)) DEALLOCATE (tmc_litter_res)
    IF (ALLOCATED  (tmc_litter_wilt)) DEALLOCATE (tmc_litter_wilt)
    IF (ALLOCATED  (tmc_litter_field)) DEALLOCATE (tmc_litter_field)
    IF (ALLOCATED  (tmc_litter_sat)) DEALLOCATE (tmc_litter_sat)
    IF (ALLOCATED  (tmc_litter_awet)) DEALLOCATE (tmc_litter_awet)
    IF (ALLOCATED  (tmc_litter_adry)) DEALLOCATE (tmc_litter_adry)
    IF (ALLOCATED  (tmc_litt_wet_mea)) DEALLOCATE (tmc_litt_wet_mea)
    IF (ALLOCATED  (tmc_litt_dry_mea)) DEALLOCATE (tmc_litt_dry_mea)
    IF (ALLOCATED  (qflux00)) DEALLOCATE (qflux00)
    IF (ALLOCATED  (ru_ns)) DEALLOCATE (ru_ns)
    IF (ALLOCATED  (dr_ns)) DEALLOCATE (dr_ns)
    IF (ALLOCATED  (tr_ns)) DEALLOCATE (tr_ns)
    IF (ALLOCATED  (cvs_over_veg)) DEALLOCATE (cvs_over_veg)
    IF (ALLOCATED  (corr_veg_soil)) DEALLOCATE (corr_veg_soil)
    IF (ALLOCATED  (mc)) DEALLOCATE (mc)
    IF (ALLOCATED  (soilmoist)) DEALLOCATE (soilmoist)
    IF (ALLOCATED  (soil_wet)) DEALLOCATE (soil_wet)
    IF (ALLOCATED  (soil_wet_litter)) DEALLOCATE (soil_wet_litter)
    IF (ALLOCATED  (qflux)) DEALLOCATE (qflux)
    IF (ALLOCATED  (tmat)) DEALLOCATE (tmat)
    IF (ALLOCATED  (stmat)) DEALLOCATE (stmat)
    IF (ALLOCATED  (nroot)) DEALLOCATE (nroot)
    IF (ALLOCATED  (kfact_root)) DEALLOCATE (kfact_root)
    IF (ALLOCATED  (kfact)) DEALLOCATE (kfact)
    IF (ALLOCATED  (zz)) DEALLOCATE (zz)
    IF (ALLOCATED  (dz)) DEALLOCATE (dz)
    IF (ALLOCATED  (mc_lin)) DEALLOCATE (mc_lin)
    IF (ALLOCATED  (k_lin)) DEALLOCATE (k_lin)
    IF (ALLOCATED  (d_lin)) DEALLOCATE (d_lin)
    IF (ALLOCATED  (a_lin)) DEALLOCATE (a_lin)
    IF (ALLOCATED  (b_lin)) DEALLOCATE (b_lin)
!Isa
     IF (ALLOCATED  (frac_hydro_diag)) DEALLOCATE (frac_hydro_diag)
    IF (ALLOCATED  (k_lin_write)) DEALLOCATE (k_lin_write)
    RETURN 

  END SUBROUTINE hydrol_clear

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_tmc_update
!!
!>\BRIEF        This routine updates the soil moisture profiles when the vegetation fraction have changed. 
!!
!! DESCRIPTION  :
!! 
!!    This routine update tmc and mc with variation of veget_max (LAND_USE or DGVM activated)
!! 
!!
!!
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_tmc_update
  SUBROUTINE hydrol_tmc_update ( kjpindex, veget_max, soiltile )
    ! interface description
    ! input scalar
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! domain size
    ! input fields
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! fraction of PFT on each soil-hydrology tile

    ! local declaration
    INTEGER(i_std)                                      :: ji, jv, jst
    REAL(r_std), DIMENSION(nstm)                        :: vmr_soil, soil_old
    REAL(r_std)                                         :: vmr_sum, tmc_dilu, tmc_old

    ! If veget has been updated before restart (with LAND USE)
    DO ji=1,kjpindex

    
        IF ( ANY(ABS(vmr(ji,:)) .GT. EPS1) ) THEN
          vmr_soil(:)=zero
          soil_old(:)=zero

          ! Variation of frac_nobio ( = (1 - sum(veget_max)) - (1 - sum(resdist)) = -sum(vmr) )
          vmr_soil(1)=vmr_soil(1) - SUM(vmr(ji,:))
          soil_old(1)=1-SUM(resdist(ji,:))
          DO jv = 1,nvm
             jst=pref_soil_veg(jv)
             vmr_soil(jst)=vmr_soil(jst)+vmr(ji,jv)
             soil_old(jst)=soil_old(jst)+resdist(ji,jv)
          ENDDO
          ! By definition, we have
          ! sum( vmr_soil(jst) > 0 ) vmr_soil(jst) = - sum(vmr_soil(kst) <= 0) vmr_soil(kst)
          ! veget loss
          vmr_sum=SUM(vmr_soil,MASK=vmr_soil.LT.zero)

          tmc_dilu=zero
          DO jst=1,nstm
             IF ( vmr_soil(jst) < -min_sechiba ) THEN
                ! tmc_dilu = sum(vmr_soil(jst) <= 0) [ vmr_soil(jst) * tmc(jst) ] / {sum(vmr_soil(kst) < 0) vmr_soil(kst)}
                tmc_dilu = tmc_dilu + vmr_soil(jst) / vmr_sum * tmc(ji,jst)
             ENDIF
          ENDDO

          DO jst=1,nstm
             IF ( vmr_soil(jst) > min_sechiba ) THEN
                ! vmr_soil(jst) == 0 => soiltile(jst) = soil_old(jst)
                ! humtot_new = sum( vmr_soil(jst) <= 0 ) tmc(jst) * soiltile(jst)
                !              + sum( vmr_soil(jst) > 0 ) [ tmc(jst) * soil_old(jst) + tmc_dilu * vmr_soil(jst) ]
                !            = sum( jst ) [ tmc(jst) * soil_old(jst) ]
                !              + sum( vmr_soil(jst) <= 0 ) [ tmc(jst) * vmr_soil(jst) ]
                !              + sum( vmr_soil(jst) > 0 ) vmr_soil(jst)
                !                * sum(vmr_soil(lst) <= 0) ( vmr_soil(lst) * tmc(lst) )
                !                / {sum(vmr_soil(kst) <= 0) vmr_soil(kst)}
                !            = sum( jst ) [ tmc(jst) * soil_old(jst) ]
                !              + sum( vmr_soil(jst) <= 0 ) [ tmc(jst) * vmr_soil(jst) ]
                !              - sum(vmr_soil(lst) <= 0) ( vmr_soil(lst) * tmc(lst) )
                !            = humtot_old
                tmc_old = tmc(ji,jst)
                tmc(ji,jst) = ( tmc(ji,jst) * soil_old(jst) + tmc_dilu * vmr_soil(jst) ) / soiltile(ji,jst)
                IF (tmc_old .GT. min_sechiba) THEN
                   mc(ji,:,jst) = mc(ji,:,jst) * tmc(ji,jst) / tmc_old
                ELSE
                   ! create new soil : mean moisture content over all layer to conserve global humidity
                   mc(ji,:,jst) = tmc(ji,jst) / (dpu_max* mille)
!
!??                   mc(ji,:,jst) = mc(ji,:,jst) * soiltile(ji,jst)
!
!??       CALL setvar_p (mc(ji,:,jst), val_exp, 'HYDROL_MOISTURE_CONTENT', 0.3_r_std)
!
!??                   tmc(ji,jst) / dpu_max = mc_
!??                   mc = tmc(ji,jst)
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO

!Chloe
IF (jst .EQ.nstm .AND. tmc(1,4) .GT. tmcs(1,4)) THEN
write(18,*) 'tmc init', tmc(1,4)
ENDIF

    tmc_init_updated = .TRUE.

    IF (long_print) WRITE (numout,*) ' hydrol_tmc_update done '

  END SUBROUTINE hydrol_tmc_update

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_var_init
!!
!>\BRIEF        This routine initializes HYDROLOGIC variables.  
!!
!! DESCRIPTION  :
!! - 1 compute the depths
!! - 2 compute the profile for roots
!! - 3 compute the profile for ksat, a and n Van Genuchten parameter
!! - 4 compute the linearized values of k, a, b and d for the resolution of Fokker Planck equation
!! - 5 water reservoirs initialisation
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_var_init

  SUBROUTINE hydrol_var_init (kjpindex, veget, veget_max, soiltile, njsc,  &
       & mx_eau_var, shumdiag,shumdiag_perma, k_litt, &
       & litterhumdiag, drysoil_frac, evap_bare_lim, ok_routage_peat) 

    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! domain size
    ! input fields
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! max fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget         !! fraction of vegetation type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc          !! indexing of PFT to soiltile 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! fraction of PFT on each soil-hydrology tile
     LOGICAL, INTENT(in)                   :: ok_routage_peat     
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: mx_eau_var    !!
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out) :: shumdiag      !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out) :: shumdiag_perma      !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: k_litt        !! litter cond.
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: litterhumdiag !! litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: drysoil_frac  !! function of litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)       :: evap_bare_lim !! 
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv, jd,jst, jsc, jsl, i
    REAL(r_std)                                         :: m, frac
    REAL(r_std)                                         :: avan_mod, nvan_mod !! modified VG parameters with exponantial profile
    REAL(r_std), DIMENSION(nslm,nscm)                   :: afact, nfact  !! multiplicative factor for decay of a and n 
    ! parameters for "soil densification" with depth
    REAL(r_std)                                         :: dp_comp       !! depth at which the 'compacted value' (Van Genuchten) of 
                                                                         !! ksat is reached 
    REAL(r_std)                                         :: f_ks          !! exponential factor for decay of ksat with depth
    ! Fixed parameters from fitted relationships
    REAL(r_std), PARAMETER                              :: n0 = 0.95     !! fitted value for relation log((n-n0)/(n_ref-n0)) = 
                                                                         !! nk_rel * log(k/k_ref) 
    REAL(r_std), PARAMETER                              :: nk_rel = 0.34 !! fitted value for relation log((n-n0)/(n_ref-n0)) = 
                                                                         !! nk_rel * log(k/k_ref) 
    REAL(r_std), PARAMETER                              :: a0 = 0.00012  !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref) 
    REAL(r_std), PARAMETER                              :: ak_rel = 0.53 !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref) 
    REAL(r_std)                                         :: kfact_max     !! Maximum factor for K due to depth
    REAL(r_std)                                         :: k_tmp, tmc_litter_ratio




!!??Aurelien: Les 3 parametres qui suient pourait peut-�tre mis dans hydrol_init?
    !
    !
    !
    !Config Key   = MAXK_PARAM
    !Config Desc  = Maximum Factor for Ks increase due to vegetation
    !Config Def   = 10.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    !
    kfact_max = 10.0
    CALL getin_p ('MAXK_PARAM', kfact_max)
    !
    !Config Key   = KPROF_PARAM
    !Config Desc  = Factor for Ks decay with depth
    !Config Def   = 2.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =  
    !Config Units = [-]
    !
    f_ks = deux
    CALL getin_p ('KPROF_PARAM', f_ks)
    !
    !Config Key   = KDEPTH_PARAM
    !Config Desc  = Depth for compacted value of Ks 
    !Config Def   = 0.3
    !Config If    = HYDROL_CWRR 
    !Config Help  =  
    !Config Units = [-]
    !
    dp_comp = 0.3
    CALL getin_p ('KDEPTH_PARAM', dp_comp)
    !
    !
    DO jst=1,nstm
       !-
    !! 1 compute the depths
       !-
       zz(1,jst) = zero
       dz(1,jst) = zero
       DO jsl=2,nslm
          zz(jsl,jst) = dpu_max* mille*((2**(jsl-1))-1)/ ((2**(nslm-1))-1)
          dz(jsl,jst) = zz(jsl,jst)-zz(jsl-1,jst)
       ENDDO
       zz(nslm+1,jst) = zz(nslm,jst)
       dz(nslm+1,jst) = zero

       !-
    !! 2 compute the profile for roots
       !-
         !! ??Aurelien: Je ne sais pas dou vient cette formule mais il semble manquer des - devant humcste et il faudrait verifier 
       !! zz ou dz??Aurelien: Reference? 
       DO jv = 1,nvm
          DO jsl = 2, nslm-1
             nroot(jv,jst,jsl) = (EXP(-humcste(jv)*zz(jsl,jst)/mille)) * &
                     & (EXP(humcste(jv)*dz(jsl,jst)/mille/deux) - &
                     & EXP(-humcste(jv)*dz(jsl+1,jst)/mille/deux))/ &
                     & (EXP(-humcste(jv)*dz(2,jst)/mille/deux) &
                     & -EXP(-humcste(jv)*zz(nslm,jst)/mille))
          ENDDO
       ENDDO
       DO jv=1,nvm
          nroot(jv,jst,1) = zero
          nroot(jv,jst,nslm) = (EXP(humcste(jv)*dz(nslm,jst)/mille/deux) -un) * &
                  & EXP(-humcste(jv)*zz(nslm,jst)/mille) / &
                  & (EXP(-humcste(jv)*dz(2,jst)/mille/deux) &
                  & -EXP(-humcste(jv)*zz(nslm,jst)/mille))
       ENDDO
    ENDDO

       ! An additional exponential factor for ks depending on the amount of roots in the soil 
       ! through a geometric average over the vegets
       !!??Aurelien: Pkoi utiliser ks_usda?
    kfact_root(:,:,:) = un
       DO jsl = 1, nslm
          DO jv = 2, nvm
             jst = pref_soil_veg(jv)
             DO ji = 1, kjpindex
                IF(soiltile(ji,jst) .GT. min_sechiba) THEN
                   kfact_root(ji,jsl,jst) = kfact_root(ji,jsl,jst) * &
                        & MAX((MAXVAL(ks_usda)/ks(njsc(ji)))**(- veget(ji,jv) * (humcste(jv)*zz(jsl,jst)/mille - un)/deux), &
                        un)

                ENDIF
             ENDDO
          ENDDO
       ENDDO
    !-
    !! 3 Compute the profile for ksat, a and n
    !-

    ! This is used for FAO2 map (with soil type 2 and 4 being heterogeneous with depth)
    IF (nscm .EQ. 5) THEN
       ! For every soil texture
       DO jsc = 1, nscm
          DO jsl=1,nslm
             frac = MAX(zero, MIN(un, f_ks * (zz(jsl,jsc)/mille - dp_comp)))
             kfact(jsl,jsc) = ( ks(MIN(nscm,jsc+1))/ks(jsc) )**frac
             nfact(jsl,jsc) = ( (nvan(MIN(nscm,jsc+1))-n0)/(nvan(jsc)-n0) )**frac
             afact(jsl,jsc) = ( (avan(MIN(nscm,jsc+1))-a0)/(avan(jsc)-a0) )**frac
          ENDDO
       ENDDO
    ! And this is to put a global soiltype profile 
    ELSE
       ! For every soil texture
       DO jsc = 1, nscm
          DO jsl=1,nslm
             kfact(jsl,jsc) = MIN(MAX(EXP(- f_ks * (zz(jsl,jsc)/mille - dp_comp)), un/kfact_max),un)
             !!??Aurelien: comment nk_rel et ak_rel ont-ils ete choisi? Et ce sont les memes pour chaque type de sol !?!
             nfact(jsl,jsc) = ( kfact(jsl,jsc) )**nk_rel
             afact(jsl,jsc) = ( kfact(jsl,jsc) )**ak_rel
          ENDDO
       ENDDO
    ENDIF
    
    ! For every soil texture
    DO jsc = 1, nscm
       !-
    !! 4 compute the linearized values of k, a, b and d
       !-
    ! Calcul the matrix coef for dublin model:
    ! pice-wise linearised hydraulic conductivity k_lin=alin * mc_lin + b_lin
    ! and diffusivity d_lin in each interval of mc, called mc_lin,
    ! between imin, for residual mcr, 
    ! and imax for saturation mcs.


       mc_lin(imin,jsc)=mcr(jsc)
       mc_lin(imax,jsc)=mcs(jsc)
       
    

       
       DO ji= imin+1, imax-1 
          mc_lin(ji,jsc) = mcr(jsc) + (ji-imin)*(mcs(jsc)-mcr(jsc))/(imax-imin)    
          
       ENDDO




       DO jsl = 1, nslm

          nvan_mod = n0 + (nvan(jsc)-n0) * nfact(jsl,jsc)
          avan_mod = a0 + (avan(jsc)-a0) * afact(jsl,jsc)


          !!??Aurelien: comment n0 et a0 ont-ils ete choisi? Et ce sont les memes pour chaque type de sol !?!
          m = un - un / nvan_mod
          ! Perhaps we may have problems with precision here for some machines (k being very small for ji=imin+1
          ! How can we handle it?
          DO ji = imax,imin+1,-1 
             frac=MIN(un,(mc_lin(ji,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
             k_lin(ji,jsl,jsc) = ks(jsc) * kfact(jsl,jsc) * (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2
          ENDDO



          ! We have to avoid k=0
          k_lin(imin,jsl,jsc) = k_lin(imin+1,jsl,jsc)/mille

          DO ji = imin,imax-1 
             a_lin(ji,jsl,jsc) = (k_lin(ji+1,jsl,jsc)-k_lin(ji,jsl,jsc)) / (mc_lin(ji+1,jsc)-mc_lin(ji,jsc))
             b_lin(ji,jsl,jsc)  = k_lin(ji,jsl,jsc) - a_lin(ji,jsl,jsc)*mc_lin(ji,jsc)


             IF(ji.NE.imin.AND.ji.NE.imax-1) THEN

                frac=MIN(un,(mc_lin(ji,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
                d_lin(ji,jsl,jsc) =(k_lin(ji,jsl,jsc) / (avan_mod*m*nvan_mod)) *  &
                     &  ( (frac**(-un/m))/(mc_lin(ji,jsc)-mcr(jsc)) ) * &
                  &  (  frac**(-un/m) -un ) ** (-m)
                frac=MIN(un,(mc_lin(ji+1,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
                d_lin(ji+1,jsl,jsc) =(k_lin(ji+1,jsl,jsc) / (avan_mod*m*nvan_mod))*&
                     &  ( (frac**(-un/m))/(mc_lin(ji+1,jsc)-mcr(jsc)) ) * &
                  &  (  frac**(-un/m) -un ) ** (-m)
                d_lin(ji,jsl,jsc) = undemi * (d_lin(ji,jsl,jsc)+d_lin(ji+1,jsl,jsc))
                               ELSEIF(ji.EQ.imin) THEN
                d_lin(ji,jsl,jsc) = zero
             ELSEIF(ji.EQ.imax-1) THEN

                frac=MIN(un,(mc_lin(ji,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
                d_lin(ji,jsl,jsc) =(k_lin(ji,jsl,jsc) / (avan_mod*m*nvan_mod)) * &
                     & ( (frac**(-un/m))/(mc_lin(ji,jsc)-mcr(jsc)) ) *  &
                  & (  frac**(-un/m) -un ) ** (-m)
                ds(jsl,jsc) = d_lin(ji,jsl,jsc)
                   
             ENDIF
          ENDDO
       ENDDO
    ENDDO


    !! 5 Water reservoir initialisation
    !
!!$    DO jst = 1,nstm
!!$       DO ji = 1, kjpindex
!!$          mx_eau_var(ji) = mx_eau_var(ji) + soiltile(ji,jst)*&
!!$               &   dpu_max*mille*mcs(njsc(ji))
!!$       END DO
!!$    END DO
!!$    IF (check_CWRR) THEN
!!$       IF ( ANY ( ABS( mx_eau_var(:) - dpu_max*mille*mcs(njsc(:)) ) > min_sechiba ) ) THEN
!!$          ji=MAXLOC ( ABS( mx_eau_var(:) - dpu_max*mille*mcs(njsc(:)) ) , 1)
!!$          WRITE(numout, *) "Erreur formule simplifiée mx_eau_var ! ", mx_eau_var(ji), dpu_max*mille*mcs(njsc(ji))
!!$          WRITE(numout, *) "err = ",ABS(mx_eau_var(ji) - dpu_max*mille*mcs(njsc(ji)))
!!$          STOP 1
!!$       ENDIF
!!$    ENDIF

    mx_eau_var(:) = zero
    mx_eau_var(:) = dpu_max*mille*mcs(njsc(:)) 

    DO ji = 1,kjpindex

       IF (vegtot(ji) .LE. zero) THEN
          mx_eau_var(ji) = mx_eau_eau*deux
          !!??Aurelien: c koi ce deux? a remplacer par dpu_max?
          !!Et �a veut dire quoi vegtot=0? cad frac_nobio=1? Et si 0<frac_nobio<1 ???
       ENDIF

    END DO

    ! Compute the litter humidity, shumdiag and fry
    litterhumdiag(:) = zero
    k_litt(:) = zero
    tmc_litt_mea(:) = zero
    tmc_litt_wet_mea(:) = zero
    tmc_litt_dry_mea(:) = zero
    shumdiag(:,:) = zero
    shumdiag_perma(:,:) = zero ! Isa
    soilmoist(:,:) = zero
    humtot(:) = zero
    tmc(:,:) = zero
    swi(:) = zero

    ! Loop on soil types to compute the variables (ji,jst)
    DO jst=1,nstm 
       DO ji = 1, kjpindex
                   
            tmcs(ji,jst)=dpu_max* mille*mcs(njsc(ji))
            tmcr(ji,jst)=dpu_max* mille*mcr(njsc(ji))

       ENDDO
    ENDDO
      
!Chloe :
DO ji=1,kjpindex
mcs_njsc(ji)=mcs(njsc(ji))
ENDDO

 
    ! The total soil moisture for each soil type:

    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc(ji,jst)= dz(2,jst) * ( trois*mc(ji,1,jst)+ mc(ji,2,jst))/huit
          !Chloe mc_obj_peat
          IF (ok_routage_peat .AND. jst .EQ. 4) THEN 
            tmc_obj_peat(ji,jst)=dz(2,jst) * ( trois*mc_obj_peat(ji,1,jst)+ mc_obj_peat(ji,2,jst))/huit
          ENDIF
       END DO
    ENDDO

    DO jst=1,nstm 
       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) * ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
              !Chloe mc_obj_peat
              IF (ok_routage_peat .AND. jst .EQ. 4) THEN           
                tmc_obj_peat(ji,jst)=tmc_obj_peat(ji,jst) + dz(jsl,jst) * ( trois*mc_obj_peat(ji,jsl,jst) &
& + mc_obj_peat(ji,jsl-1,jst))/huit + dz(jsl+1,jst)*(trois*mc_obj_peat(ji,jsl,jst) + mc_obj_peat(ji,jsl+1,jst))/huit
              ENDIF
          END DO
       END DO
    ENDDO

IF (jst .EQ. nstm .AND. tmc(1,4) .GT. tmcs(1,4)) THEN
write(15,*) 'tmc AV wat2inf,tmcs', tmc(1,4), tmcs(1,4)
ENDIF

    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) +  dz(nslm,jst) * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
          !Chloe mc_obj_peat       
          IF (ok_routage_peat .AND. jst .EQ. 4) THEN     
          tmc_obj_peat(ji,jst) = tmc_obj_peat(ji,jst) + dz(nslm,jst) * (trois * mc_obj_peat(ji,nslm,jst) + &
	 & mc_obj_peat(ji,nslm-1,jst))/huit
          tmc_obj_peat(ji,jst) = tmc_obj_peat(ji,jst) + water2infilt(ji,jst)
          ENDIF
       ENDDO
    END DO


IF (jst .EQ. nstm .AND. tmc(1,4) .GT. tmcs(1,4)) THEN
write(17,*) 'tmc AP wat2inf, wat2inf,tmcs', tmc(1,4), water2infilt(1,4),tmcs(1,4)
ENDIF
    ! If veget has been updated before restart (with LAND USE or DGVM),
    ! tmc and mc must be modified with respect to humtot conservation.
    IF ( control%do_land_use .OR. control%ok_dgvm ) THEN
       CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile )
    ENDIF

    ! The litter variables:
    ! level 1
    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2,jst) * (trois*mc(ji,1,jst)+mc(ji,2,jst))/huit
          tmc_litter_wilt(ji,jst) = dz(2,jst) * mcw(njsc(ji)) / deux
          tmc_litter_res(ji,jst) = dz(2,jst) * mcr(njsc(ji)) / deux
          tmc_litter_field(ji,jst) = dz(2,jst) * mcf(njsc(ji)) / deux
          tmc_litter_sat(ji,jst) = dz(2,jst) * mcs(njsc(ji)) / deux
          tmc_litter_awet(ji,jst) = dz(2,jst) * mc_awet(njsc(ji)) / deux
          tmc_litter_adry(ji,jst) = dz(2,jst) * mc_adry(njsc(ji)) / deux
       ENDDO
    END DO
    ! sum from level 2 to 4
    DO jst=1,nstm 
       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl,jst) * & 
                  & ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
             tmc_litter_wilt(ji,jst) = tmc_litter_wilt(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))*& 
                  & mcw(njsc(ji))/deux
             tmc_litter_res(ji,jst) = tmc_litter_res(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))*& 
                  & mcr(njsc(ji))/deux
             tmc_litter_sat(ji,jst) = tmc_litter_sat(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mcs(njsc(ji))/deux
             tmc_litter_field(ji,jst) = tmc_litter_field(ji,jst) + &
                  & (dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mcf(njsc(ji))/deux
             tmc_litter_awet(ji,jst) = tmc_litter_awet(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mc_awet(njsc(ji))/deux
             tmc_litter_adry(ji,jst) = tmc_litter_adry(ji,jst) + &
                  & (dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mc_adry(njsc(ji))/deux

          END DO
       END DO
    END DO

    ! subsequent calcul of soil_wet_litter (tmc-tmcw)/(tmcf-tmcw)
    DO jst=1,nstm 
       DO ji=1,kjpindex
          soil_wet_litter(ji,jst)=MIN(un, MAX(zero,&
               &(tmc_litter(ji,jst)-tmc_litter_wilt(ji,jst))/&
               & (tmc_litter_field(ji,jst)-tmc_litter_wilt(ji,jst)) ))
       END DO
    ENDDO

    ! Soil wetness profiles (mc-mcw)/(mcs-mcw)
    DO jst=1,nstm 
       DO ji=1,kjpindex
          soil_wet(ji,1,jst) = MIN(un, MAX(zero,&
               &(trois*mc(ji,1,jst) + mc(ji,2,jst) - quatre*mcw(njsc(ji)))&
               & /(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
    !!??Aurelien: a quoi sert cette ligne?

            
          humrelv(ji,1,jst) = zero
       ENDDO
    END DO
    !write(*,*) 'hydrol 2516: soil_wet=',soil_wet(1,1,1)
    !write(*,*) 'mcw(njsc(1))=',mcw(njsc(1))
    !write(*,*) 'mcs(njsc(1))=',mcs(njsc(1))

    DO jst=1,nstm 
       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             soil_wet(ji,jsl,jst) = MIN(un, MAX(zero,&
                  & (trois*mc(ji,jsl,jst) + & 
                  & mc(ji,jsl-1,jst) *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & + mc(ji,jsl+1,jst)*(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & - quatre*mcw(njsc(ji))) / (quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))

          END DO
       END DO
    END DO
    !write(*,*) 'hydrol 2540: soil_wet=',soil_wet(1,2,1)

    DO jst=1,nstm 
       DO ji=1,kjpindex
          soil_wet(ji,nslm,jst) = MIN(un, MAX(zero,&
               & (trois*mc(ji,nslm,jst) &
               & + mc(ji,nslm-1,jst)-quatre*mcw(njsc(ji)))/(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
       ENDDO
    END DO
    !write(*,*) 'hydrol 2551: soil_wet=',soil_wet(1,nslm,1)

    !Now we compute the grid averaged values
    DO jst=1,nstm        
       DO ji=1,kjpindex
          !
          IF ( tmc_litter(ji,jst) < tmc_litter_res(ji,jst)) THEN
             i = imin
          ELSE
             tmc_litter_ratio = (tmc_litter(ji,jst)-tmc_litter_res(ji,jst)) / &
                  & (tmc_litter_sat(ji,jst)-tmc_litter_res(ji,jst))
             i= MAX(MIN(INT((imax-imin)*tmc_litter_ratio)+imin , imax-1), imin)
          ENDIF
          ! k_litt is an averaged conductivity for saturated infiltration in the 'litter' layer
          ! This is used for reinfiltration from surface water
          k_tmp = MAX(k_lin(i,1,njsc(ji))*ks(njsc(ji)), zero)
          k_litt(ji) = k_litt(ji) + soiltile(ji,jst) * SQRT(k_tmp)
       ENDDO
    ENDDO

    DO jst=1,nstm
       DO ji = 1, kjpindex
          humtot(ji) = humtot(ji) + soiltile(ji,jst) * tmc(ji,jst) 

          litterhumdiag(ji) = litterhumdiag(ji) + &
               & soil_wet_litter(ji,jst) * soiltile(ji,jst)

          tmc_litt_wet_mea(ji) =  tmc_litt_wet_mea(ji) + & 
               & tmc_litter_awet(ji,jst)* soiltile(ji,jst)

          tmc_litt_dry_mea(ji) = tmc_litt_dry_mea(ji) + &
               & tmc_litter_adry(ji,jst) * soiltile(ji,jst) 

          tmc_litt_mea(ji) = tmc_litt_mea(ji) + &
               & tmc_litter(ji,jst) * soiltile(ji,jst) 
       ENDDO
    ENDDO
!Isa : correction du calcul de shumdiag pour le mettre sur l'�chelle diag. soilmoist n'est pas concern�.

CALL calcule_frac_hydro_diag

    DO jst=1,nstm 
       DO jd=1,nbdl
!Isa shumdiag_perma
        DO ji=1,kjpindex
                if (ok_shumdiag_interpol) then
		  DO jsl = 1, nslm
             	    shumdiag(ji,jd)= shumdiag(ji,jd) + soil_wet(ji,jsl,jst)  &
                  & *frac_hydro_diag(jsl,jd)* &
                  & ((mcs(njsc(ji))-mcw(njsc(ji))) &
                  & /(mcf(njsc(ji))-mcw(njsc(ji)))) * &
                  & soiltile(ji,jst)


		  ENDDO !DO jsl = 1, nslm
                else !if (ok_shumdiag_interpol) then
                  jsl=jd
                  shumdiag(ji,jsl)= shumdiag(ji,jsl) + soil_wet(ji,jsl,jst) * &
                  & ((mcs(njsc(ji))-mcw(njsc(ji))) &
                  & /(mcf(njsc(ji))-mcw(njsc(ji)))) * &
                  & soiltile(ji,jst)
!                  if (ji.eq.1) then
!                      write(*,*) 'hydrol 2605: shumdiag(ji,jsl)=',shumdiag(ji,jsl)
!                      write(*,*) 'soil_wet(ji,jsl,jst)=',soil_wet(ji,jsl,jst)
!                      write(*,*) 'soiltile(ji,jst)=',soiltile(ji,jst)
!                      write(*,*) 'mcs(njsc(ji))-mcw(njsc(ji))=',mcs(njsc(ji))-mcw(njsc(ji))
!                      write(*,*) 'mcf(njsc(ji))-mcw(njsc(ji))=',mcf(njsc(ji))-mcw(njsc(ji))
!                  endif !if (ji.eq.1) then
                endif !if (ok_shumdiag_interpol) then


                shumdiag(ji,jd) = MAX(MIN(shumdiag(ji,jd), un), zero) 

            if (ok_shumdiag_perma) then
                 DO jsl = 1, nslm    
             	   shumdiag_perma(ji,jd)= shumdiag_perma(ji,jd)  &
                   & + mc(ji,jsl,jst) *frac_hydro_diag(jsl,jd) &
                   & /mcs(njsc(ji))*soiltile(ji,jst)
		           !Chloe+ mettre mcs(njsc(ji)) au lieu de mcs(jst)


                 ENDDO !DO jsl = 1, nslm
            endif !if(ok_shumdiag_perma) then

        END DO !DO ji=1,kjpindex
       END DO !DO jd=1,nbdl
        DO jsl = 1, nslm
             DO ji=1,kjpindex
	             soilmoist(ji,jsl)=soilmoist(ji,jsl) &
                  & +mc(ji,jsl,jst)*soiltile(ji,jst)



	         ENDDO
	    ENDDO !DO jsl = 1, nslm
    END DO !DO jst=1,nstm 
    write(*,*) 'hydrol 2393: shumdiag=',shumdiag(1,1)


if(.NOT.ok_shumdiag_perma) shumdiag_perma(:,:)=shumdiag(:,:)
    !
    !
    DO ji=1,kjpindex
       IF ( tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji) > zero ) THEN
          drysoil_frac(ji) = un + MAX( MIN( (tmc_litt_dry_mea(ji) - tmc_litt_mea(ji)) / &
               & (tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji)), zero), - un)
       ELSE
          drysoil_frac(ji) = zero
       ENDIF
    END DO

    evap_bare_lim = zero

    IF (long_print) WRITE (numout,*) ' hydrol_var_init done '
  END SUBROUTINE hydrol_var_init

 

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_snow
!!
!>\BRIEF        This routine computes snow processes. 
!!
!! DESCRIPTION  :
!! - 0 initialisation
!! - 1 On vegetation
!! - 1.1 Compute snow masse
!! - 1.2 Sublimation 
!! - 1.2.1 Check that sublimation on the vegetated fraction is possible.
!! - 1.3. snow melt only if temperature positive
!! - 1.3.1 enough snow for melting or not
!! - 1.3.2 not enough snow
!! - 1.3.3 negative snow - now snow melt
!! - 1.4 Snow melts only on weight glaciers
!! - 2 On Land ice
!! - 2.1 Compute snow
!! - 2.2 Sublimation 
!! - 2.3 Snow melt only for continental ice fraction
!! - 2.3.1 If there is snow on the ice-fraction it can melt
!! - 2.4 Snow melts only on weight glaciers 
!! - 3 On other surface types - not done yet
!! - 4 computes total melt (snow and ice)
!! - 5 computes snow age on veg and ice (for albedo)
!! - 5.1 Snow age on vegetation
!! - 5.2 Snow age on ice
!! - 6 Diagnose the depth of the snow layer
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_snow

  SUBROUTINE hydrol_snow (kjpindex, dtradia, precip_rain, precip_snow , temp_sol_new, soilcap,&
       & frac_nobio, totfrac_nobio, vevapnu, vevapsno, snow, snow_age, snow_nobio, snow_nobio_age, &
       & tot_melt, snowdepth)

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size
    REAL(r_std), INTENT (in)                                 :: dtradia       !! Time step in seconds
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain   !! Rainfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_snow   !! Snow precipitation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_sol_new  !! New soil temperature
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: soilcap       !! Soil capacity
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: totfrac_nobio !! Total fraction of continental ice+lakes+ ...

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: tot_melt      !! Total melt  
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snowdepth     !! Snow depth

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu       !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapsno      !! Snow evaporation
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow          !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow_age      !! Snow age
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio    !! Ice water balance
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio_age!! Snow age on ice, lakes, ...

    !! 0.4 Local variables

    INTEGER(i_std)                               :: ji, jv
    REAL(r_std), DIMENSION (kjpindex)             :: d_age  !! Snow age change
    REAL(r_std), DIMENSION (kjpindex)             :: xx     !! temporary
    REAL(r_std)                                   :: snowmelt_tmp !! The name says it all !

    !
    ! for continental points
    !

    !
    !!_0 initialisation
    !
    DO jv = 1, nnobio
       DO ji=1,kjpindex
          subsnownobio(ji,jv) = zero
       ENDDO
    ENDDO
    DO ji=1,kjpindex
       subsnowveg(ji) = zero
       snowmelt(ji) = zero
       icemelt(ji) = zero
       subsinksoil(ji) = zero
       tot_melt(ji) = zero
    ENDDO



    !
    !! 1 On vegetation
    !
    DO ji=1,kjpindex
       !
    !! 1.1 Compute snow masse
       !
       snow(ji) = snow(ji) + (un - totfrac_nobio(ji))*precip_snow(ji)
       !
       !
    !! 1.2 Sublimation 
       !      Separate between vegetated and no-veget fractions 
       !      Care has to be taken as we might have sublimation from the
       !      the frac_nobio while there is no snow on the rest of the grid.
       !
       IF ( snow(ji) > snowcri ) THEN
          subsnownobio(ji,iice) = frac_nobio(ji,iice)*vevapsno(ji)
          subsnowveg(ji) = vevapsno(ji) - subsnownobio(ji,iice)
       ELSE
          ! Correction Nathalie - Juillet 2006.
          ! On doit d'abord tester s'il existe un frac_nobio!
          ! Pour le moment je ne regarde que le iice
          IF ( frac_nobio(ji,iice) .GT. min_sechiba) THEN
             subsnownobio(ji,iice) = vevapsno(ji)
             subsnowveg(ji) = zero
          ELSE 
             subsnownobio(ji,iice) = zero
             subsnowveg(ji) = vevapsno(ji)
          ENDIF
       ENDIF
       !
       !
    !! 1.2.1 Check that sublimation on the vegetated fraction is possible.
       !
       IF (subsnowveg(ji) .GT. snow(ji)) THEN
          ! What could not be sublimated goes into soil evaporation
          !         vevapnu(ji) = vevapnu(ji) + (subsnowveg(ji) - snow(ji))
          IF( (un - totfrac_nobio(ji)).GT.min_sechiba) THEN
             subsinksoil (ji) = (subsnowveg(ji) - snow(ji))/ (un - totfrac_nobio(ji))
          END IF
          ! Sublimation is thus limited to what is available
          subsnowveg(ji) = snow(ji)
          snow(ji) = zero
          vevapsno(ji) = subsnowveg(ji) + subsnownobio(ji,iice)
       ELSE
          snow(ji) = snow(ji) - subsnowveg(ji)
       ENDIF
       !
    !! 1.3. snow melt only if temperature positive
       !
       IF (temp_sol_new(ji).GT.tp_00) THEN
          !
          IF (snow(ji).GT.sneige) THEN
             !
             snowmelt(ji) = (un - frac_nobio(ji,iice))*(temp_sol_new(ji) - tp_00) * soilcap(ji) / chalfu0
             !
    !! 1.3.1 enough snow for melting or not
             !
             IF (snowmelt(ji).LT.snow(ji)) THEN
                snow(ji) = snow(ji) - snowmelt(ji)
             ELSE
                snowmelt(ji) = snow(ji)
                snow(ji) = zero
             END IF
             !
          ELSEIF (snow(ji).GE.zero) THEN
             !
    !! 1.3.2 not enough snow
             !
             snowmelt(ji) = snow(ji)
             snow(ji) = zero
          ELSE
             !
    !! 1.3.3 negative snow - now snow melt
             !
             snow(ji) = zero
             snowmelt(ji) = zero
             WRITE(numout,*) 'hydrol_snow: WARNING! snow was negative and was reset to zero. '
             !
          END IF

       ENDIF
    !! 1.4 Snow melts only on weight glaciers
       ! Ice melt only if there is more than a given mass : maxmass_glacier,
       ! Ajouts Edouard Davin / Nathalie de Noblet add extra to melting
       !
       IF ( snow(ji) .GT. maxmass_glacier ) THEN
          snowmelt(ji) = snowmelt(ji) + (snow(ji) - maxmass_glacier)
          snow(ji) = maxmass_glacier
       ENDIF
       !
    END DO
    !
    !! 2 On Land ice
    !
    DO ji=1,kjpindex
       !
    !! 2.1 Compute snow
       !
       !!??Aurelien: pkoi mettre precip_rain en dessous? We considere liquid precipitations becomes instantly snow?  
       snow_nobio(ji,iice) = snow_nobio(ji,iice) + frac_nobio(ji,iice)*precip_snow(ji) + &
            & frac_nobio(ji,iice)*precip_rain(ji)
       !
    !! 2.2 Sublimation 
       !      Was calculated before it can give us negative snow_nobio but that is OK
       !      Once it goes below a certain values (-maxmass_glacier for instance) we should kill
       !      the frac_nobio(ji,iice) !
       !
       snow_nobio(ji,iice) = snow_nobio(ji,iice) - subsnownobio(ji,iice)
       !
    !! 2.3 Snow melt only for continental ice fraction
       !
       snowmelt_tmp = zero
       IF (temp_sol_new(ji) .GT. tp_00) THEN
          !
    !! 2.3.1 If there is snow on the ice-fraction it can melt
          !
          snowmelt_tmp = frac_nobio(ji,iice)*(temp_sol_new(ji) - tp_00) * soilcap(ji) / chalfu0
          !
          IF ( snowmelt_tmp .GT. snow_nobio(ji,iice) ) THEN
             snowmelt_tmp = MAX( zero, snow_nobio(ji,iice))
          ENDIF
          snowmelt(ji) = snowmelt(ji) + snowmelt_tmp
          snow_nobio(ji,iice) = snow_nobio(ji,iice) - snowmelt_tmp
          !
       ENDIF
       !
    !! 2.4 Snow melts only on weight glaciers 
       !      Ice melt only if there is more than a given mass : maxmass_glacier, 
       !
       IF ( snow_nobio(ji,iice) .GT. maxmass_glacier ) THEN
          icemelt(ji) = snow_nobio(ji,iice) - maxmass_glacier
          snow_nobio(ji,iice) = maxmass_glacier
       ENDIF
       !
    END DO

    !
    !! 3 On other surface types - not done yet
    !
    IF ( nnobio .GT. 1 ) THEN
       WRITE(numout,*) 'WE HAVE',nnobio-1,' SURFACE TYPES I DO NOT KNOW'
       WRITE(numout,*) 'CANNOT TREAT SNOW ON THESE SURFACE TYPES'
       STOP 'in hydrol_snow' 
    ENDIF

    !
    !! 4 computes total melt (snow and ice)
    !
    DO ji = 1, kjpindex
       tot_melt(ji) = icemelt(ji) + snowmelt(ji)
    ENDDO

    !
    !! 5 computes snow age on veg and ice (for albedo)
    !
    DO ji = 1, kjpindex
       !
    !! 5.1 Snow age on vegetation
       !
       IF (snow(ji) .LE. zero) THEN
          snow_age(ji) = zero
       ELSE
          snow_age(ji) =(snow_age(ji) + (un - snow_age(ji)/max_snow_age) * dtradia/one_day) &
               & * EXP(-precip_snow(ji) / snow_trans)
       ENDIF
       !
    !! 5.2 Snow age on ice
       !
       ! age of snow on ice: a little bit different because in cold regions, we really
       ! cannot negect the effect of cold temperatures on snow metamorphism any more.
       !
       IF (snow_nobio(ji,iice) .LE. zero) THEN
          snow_nobio_age(ji,iice) = zero
       ELSE
          !
          d_age(ji) = ( snow_nobio_age(ji,iice) + &
               &  (un - snow_nobio_age(ji,iice)/max_snow_age) * dtradia/one_day ) * &
               &  EXP(-precip_snow(ji) / snow_trans) - snow_nobio_age(ji,iice)
          IF (d_age(ji) .GT. min_sechiba ) THEN
             xx(ji) = MAX( tp_00 - temp_sol_new(ji), zero )
             xx(ji) = ( xx(ji) / 7._r_std ) ** 4._r_std
             d_age(ji) = d_age(ji) / (un+xx(ji))
          ENDIF
          snow_nobio_age(ji,iice) = MAX( snow_nobio_age(ji,iice) + d_age(ji), zero )
          !
       ENDIF

    ENDDO



    !
    !! 6 Diagnose the depth of the snow layer
    !

    DO ji = 1, kjpindex
       snowdepth(ji) = snow(ji) /sn_dens
    ENDDO

    IF (long_print) WRITE (numout,*) ' hydrol_snow done '
  END SUBROUTINE hydrol_snow

   
!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_canop
!!
!>\BRIEF        This routine computes canopy processes.
!!
!! DESCRIPTION  :
!! - 1 evaporation off the continents
!! - 1.1 The interception loss is take off the canopy. 
!! - 1.2 precip_rain is shared for each vegetation type
!! - 1.3 Limits the effect and sum what receives soil
!! - 1.4 swap qsintveg to the new value
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_canop

  SUBROUTINE hydrol_canop (kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, &
       & qsintveg,precisol,tot_melt)

    ! 
    ! interface description
    !

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex    !! Domain size
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain !! Rain precipitation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: vevapwet    !! Interception loss
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget_max   !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget       !! Fraction of vegetation type 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: qsintmax    !! Maximum water on vegetation for interception
    REAL(r_std), DIMENSION  (kjpindex), INTENT (in)          :: tot_melt    !! Total melt

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)       :: precisol    !! Water fallen onto the ground (throughfall)

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(inout)     :: qsintveg    !! Water on vegetation due to interception

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv
    REAL(r_std), DIMENSION (kjpindex,nvm)                    :: zqsintvegnew

    ! boucle sur les points continentaux
    ! calcul de qsintveg au pas de temps suivant
    ! par ajout du flux interception loss
    ! calcule par enerbil en fonction
    ! des calculs faits dans diffuco
    ! calcul de ce qui tombe sur le sol
    ! avec accumulation dans precisol
    ! essayer d'harmoniser le traitement du sol nu
    ! avec celui des differents types de vegetation
    ! fait si on impose qsintmax ( ,1) = 0.0
    !
    ! loop for continental subdomain
    !
    !
    !! 1 evaporation off the continents
    !
    !! 1.1 The interception loss is take off the canopy. 
    DO jv=2,nvm
       qsintveg(:,jv) = qsintveg(:,jv) - vevapwet(:,jv)
    END DO

    !     It is raining :
    !! 1.2 precip_rain is shared for each vegetation type
    !     sum (veget (1,nvm)) must be egal to 1-totfrac_nobio.
    !     iniveget computes veget each day
    !
    qsintveg(:,1) = zero
    DO jv=2,nvm
       IF ( ok_throughfall_by_pft ) THEN
          ! Correction Nathalie - Juin 2006 - une partie de la pluie arrivera toujours sur le sol
          ! sorte de throughfall supplementaire
          qsintveg(:,jv) = qsintveg(:,jv) + veget(:,jv) * ((1-throughfall_by_pft(jv))*precip_rain(:))
       ELSE
          qsintveg(:,jv) = qsintveg(:,jv) + veget(:,jv) * precip_rain(:)
       ENDIF
    END DO

    !
    !! 1.3 Limits the effect and sum what receives soil
    !
    precisol(:,1)=veget_max(:,1)*precip_rain(:)
    DO jv=2,nvm
       DO ji = 1, kjpindex
          zqsintvegnew(ji,jv) = MIN (qsintveg(ji,jv),qsintmax(ji,jv)) 
          IF ( ok_throughfall_by_pft ) THEN
             ! correction throughfall Nathalie - Juin 2006
             precisol(ji,jv) = (veget(ji,jv)*throughfall_by_pft(jv)*precip_rain(ji)) + &
                  qsintveg(ji,jv) - zqsintvegnew (ji,jv) + &
                  (veget_max(ji,jv) - veget(ji,jv))*precip_rain(ji)
          ELSE
             precisol(ji,jv) = qsintveg(ji,jv) - zqsintvegnew (ji,jv) + &
                  (veget_max(ji,jv) - veget(ji,jv))*precip_rain(ji)
          ENDIF
       ENDDO
    END DO
    !    
    DO jv=1,nvm
       DO ji = 1, kjpindex
          IF (vegtot(ji).GT.min_sechiba) THEN
             precisol(ji,jv) = precisol(ji,jv)+tot_melt(ji)*veget_max(ji,jv)/vegtot(ji)
          ENDIF
       ENDDO
    END DO



    !   
    !
    !! 1.4 swap qsintveg to the new value
    !
    DO jv=2,nvm
       qsintveg(:,jv) = zqsintvegnew (:,jv)
    END DO

    IF (long_print) WRITE (numout,*) ' hydrol_canop done '
  END SUBROUTINE hydrol_canop


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_vegupd
!!
!>\BRIEF        Vegetation update   
!!
!! DESCRIPTION  :
!!   The vegetation cover has changed and we need to adapt the reservoir distribution 
!!   and the distribution of plants on different soil types.
!!   You may note that this occurs after evaporation and so on have been computed. It is
!!   not a problem as a new vegetation fraction will start with humrel=0 and thus will have no
!!   evaporation. If this is not the case it should have been caught above.
!!
!! - 1 Update of vegetation is it needed?
!! - 2 calculate water mass that we have to redistribute
!! - 3 put it into reservoir of plant whose surface area has grown
!! - 4 Soil tile gestion
!! - 5 update the corresponding masks
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_vegupd

  SUBROUTINE hydrol_vegupd(kjpindex, veget, veget_max, soiltile,qsintveg,resdist)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                            :: kjpindex 
    ! input fields
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)    :: veget            !! New vegetation map
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)     :: veget_max        !! Max. fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)   :: soiltile         !! fraction of PFT on each soil-hydrology tile

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg         !! Water on old vegetation 
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(inout) :: resdist          !! Old vegetation map

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji,jv,jst
    REAL(r_std), DIMENSION (kjpindex,nvm)          :: veget_exist
    REAL(r_std), DIMENSION (kjpindex,nvm)          :: qsintveg2               !! Water on new vegetation
    REAL(r_std), DIMENSION (kjpindex,nvm)          :: vmr                     !! variation of veget
    REAL(r_std), DIMENSION (kjpindex,nvm)          :: qsdq
    REAL(r_std), DIMENSION(kjpindex)               :: vegchtot,vtr, qstr, fra
    REAL(r_std), DIMENSION(kjpindex)               :: tot_corr_veg_soil
    !


    !! 1 If veget has been updated at last time step (with LAND USE or DGVM),
    !! tmc and mc must be modified with respect to humtot conservation.
    IF ( control%do_land_use .OR. control%ok_dgvm .AND. .NOT. tmc_init_updated ) THEN
       CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile )
    ELSE
       IF ( .NOT. tmc_init_updated ) THEN
          ! veget variation vmr has already been computed in hydrol_tmc_update in hydrol_var_init.
          DO jv = 1, nvm
             DO ji = 1, kjpindex
                IF ( ABS(veget_max(ji,jv)-resdist(ji,jv)) .GT. EPS1 ) THEN
                   vmr(ji,jv) = veget_max(ji,jv)-resdist(ji,jv)
                ELSE
                   vmr(ji,jv) = zero
                ENDIF
                !
             ENDDO
          ENDDO
       ENDIF
       tmc_init_updated=.FALSE.
    ENDIF
    !
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF (resdist(ji,jv) .GT. min_sechiba) THEN
             qsintveg2(ji,jv) = qsintveg(ji,jv)/resdist(ji,jv)
          ELSE
             qsintveg2(ji,jv) = zero
          ENDIF
       ENDDO
    ENDDO
    !
    vegchtot(:) = zero
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          vegchtot(ji) = vegchtot(ji) + ABS( vmr(ji,jv) )
       ENDDO
    ENDDO
    !
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF ( vegchtot(ji) .GT. min_sechiba ) THEN
             qsdq(ji,jv) = ABS(vmr(ji,jv)) * qsintveg2(ji,jv)
          ENDIF
       ENDDO
    ENDDO
    !
    !! 2 calculate water mass that we have to redistribute
    !
    qstr(:) = zero
    vtr(:) = zero
    !
    !
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF ( ( vegchtot(ji) .GT. min_sechiba ) .AND. ( vmr(ji,jv) .LT. -min_sechiba ) ) THEN
             qstr(ji) = qstr(ji) + qsdq(ji,jv)
             vtr(ji) = vtr(ji) - vmr(ji,jv)
          ENDIF
       ENDDO
    ENDDO
    !
    !! 3 put it into reservoir of plant whose surface area has grown
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF ( vegchtot(ji) .GT. min_sechiba .AND. ABS(vtr(ji)) .GT. EPSILON(un)) THEN
             fra(ji) = vmr(ji,jv) / vtr(ji)
             IF ( vmr(ji,jv) .GT. min_sechiba)  THEN
                qsintveg(ji,jv) = qsintveg(ji,jv) + fra(ji)* qstr(ji)
             ELSE
                qsintveg(ji,jv) = qsintveg(ji,jv) - qsdq(ji,jv)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
!MM underflow :
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF ( ABS(qsintveg(ji,jv)) < EPS1 ) THEN
             qsintveg(ji,jv) = EPS1
          ENDIF
       ENDDO
    ENDDO



    !! 4 Soil tile gestion
    ! Now that the work is done resdist needs an update !
    resdist(:,:) = veget_max(:,:)

    ! Remember that it is only frac_nobio + SUM(veget_max(,:)) that is equal to 1. Thus we need vegtot
    DO ji = 1, kjpindex
       vegtot(ji) = SUM(veget_max(ji,:))
    ENDDO

    ! Compute the masks for veget
    
    mask_veget(:,:) = 0
    mask_soiltile(:,:) = 0
    
    DO jst=1,nstm
       DO ji = 1, kjpindex
          IF(soiltile(ji,jst) .GT. min_sechiba) THEN
             mask_soiltile(ji,jst) = 1
          ENDIF
       END DO
    ENDDO
          
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF(veget_max(ji,jv) .GT. min_sechiba) THEN
             mask_veget(ji,jv) = 1
          ENDIF
       END DO
    END DO

    ! Distribution of the vegetation depending on the soil type
    veget_exist(:,:) = veget_max(:,:)

    ! Compute corr_veg_soil 
    corr_veg_soil(:,:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji=1,kjpindex
          ! for veget distribution used in sechiba via humrel
          IF (mask_soiltile(ji,jst).GT.0 .AND. vegtot(ji) > min_sechiba) THEN
             corr_veg_soil(ji,jv,jst)=veget_max(ji,jv)/soiltile(ji,jst)
          ENDIF
       ENDDO
    ENDDO

    IF (check_cwrr .AND. l_second_hydrol) THEN
       ! somme(soiltile * corr_veg_soil ) = 1
       tot_corr_veg_soil(:)=zero
       DO jst = 1, nstm
          DO jv = 1,nvm
             DO ji=1,kjpindex
                tot_corr_veg_soil(ji)=tot_corr_veg_soil(ji)+soiltile(ji,jst)*corr_veg_soil(ji,jv,jst)
             ENDDO
          ENDDO
       ENDDO

       DO ji=1,kjpindex
          IF ( ABS( tot_corr_veg_soil(ji) - vegtot(ji) ) > 10*EPS1 ) THEN
             WRITE(numout,*) 'corr_veg_soil SPLIT FALSE:ji=',ji,&
                  tot_corr_veg_soil(ji)
             WRITE(numout,*) 'err',ABS( tot_corr_veg_soil(ji) - vegtot(ji) )
             WRITE(numout,*) 'vegtot',vegtot(ji)
             DO jv=1,nvm
                WRITE(numout,*) 'jv,veget_max,corr_veg_soil',jv,veget_max(ji,jv),corr_veg_soil(ji,jv,:)
             END DO
             STOP 1
          ENDIF
       ENDDO
    ENDIF

    ! Tout dans cette routine est maintenant certainement obsolete (veget_max etant constant) en dehors des lignes suivantes:
    !frac_bare_ns(:,:) = zero
    !DO jst = 1, nstm
    !   DO jv = 1, nvm
    !      DO ji =1, kjpindex
    !         IF(vegtot(ji) .GT. min_sechiba) THEN
    !            frac_bare_ns(ji,jst) = frac_bare_ns(ji,jst) + corr_veg_soil(ji,jv,jst) * frac_bare(ji,jv) / vegtot(ji)
    !         ENDIF
    !      END DO
    !   ENDDO
    !END DO


    ! Tout dans cette routine est maintenant certainement obsolete (veget_max etant constant) en dehors des lignes suivantes:
    frac_bare_ns(:,:) = zero 
    DO jst = 1, nstm 
       DO jv = 1, nvm
          DO ji =1, kjpindex
             IF(vegtot(ji) .GT. min_sechiba) THEN 
                !IF (ok_reinfilt_peat) THEN
                !   frac_bare_ns(ji,jst) = frac_bare_ns(ji,jst) + corr_veg_soil(ji,jv,jst) * frac_bare(ji,jv) / vegtot(ji) * (un/trois)
                !ELSE
                    frac_bare_ns(ji,jst) = frac_bare_ns(ji,jst) + corr_veg_soil(ji,jv,jst) * frac_bare(ji,jv) / vegtot(ji)
                !ENDIF
     
             ENDIF
          END DO
       ENDDO
    END DO

    ! To supress
    vmr(:,:) = zero

    IF (long_print) WRITE (numout,*) ' hydrol_vegupd done '

    RETURN
    !
  END SUBROUTINE hydrol_vegupd


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_flood
!!
!>\BRIEF        This routine computes the evolution of the surface reservoir (floodplain).  
!!
!! DESCRIPTION  :
!! - 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil
!! - 2 Compute the total flux from floodplain floodout (transfered to routing)
!! - 3 Discriminate between precip over land and over floodplain
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_flood

  SUBROUTINE hydrol_flood (kjpindex, dtradia, vevapflo, flood_frac, flood_res, floodout)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !!
    ! input fields
    REAL(r_std), INTENT (in)                                 :: dtradia          !! Time step in seconds
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: flood_frac       !! Fraction of floodplains in grid box

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: floodout         !! Flux to take out from floodplains

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: flood_res        !! Floodplains reservoir estimate
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapflo         !! Evaporation over floodplains

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv           !! Indices
    REAL(r_std), DIMENSION (kjpindex)                        :: temp             !! 

    !- 
    !! 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil 
    !-
    DO ji = 1,kjpindex
       temp(ji) = MIN(flood_res(ji), vevapflo(ji))
    ENDDO
    DO ji = 1,kjpindex
       flood_res(ji) = flood_res(ji) - temp(ji)
       subsinksoil(ji) = subsinksoil(ji) + vevapflo(ji) - temp(ji)
       vevapflo(ji) = temp(ji)
    ENDDO

    !- 
    !! 2 Compute the total flux from floodplain floodout (transfered to routing) 
    !-
    DO ji = 1,kjpindex
       floodout(ji) = vevapflo(ji) - flood_frac(ji) * SUM(precisol(ji,:))
    ENDDO

    !-
    !! 3 Discriminate between precip over land and over floodplain
    !-
    DO jv=1, nvm
       DO ji = 1,kjpindex
          precisol(ji,jv) = precisol(ji,jv) * (1 - flood_frac(ji))
       ENDDO
    ENDDO 



    IF (long_print) WRITE (numout,*) ' hydrol_flood done'

  END SUBROUTINE hydrol_flood


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_soil
!!
!>\BRIEF        This routine computes soil processes with CWRR scheme.
!!
!! DESCRIPTION  :
!! - 0 Arrays initialisation
!! -   for each soil type
!! - 1 We compare water2infilt and water2extract to keep only difference
!! - 1.1 add to the first layer
!! - 1.2 filling layers
!! - 2 Before diffusion scheme
!! - 2.1 Some initialisation necessary for the diffusion scheme to work
!! - 2.2 coefficients are computed for the profile of mc before infiltration:
!! - 3 The infiltration is computed 
!! - 4 Coefficient are recomputed for the profile of mc after infiltration:
!! - 5 Prepar the diffusion scheme
!! - 5.1 Set the values for diffusion scheme
!! - 5.2 verifications for a good soil humidity
!! - 5.3 compute matrix coefficients
!! - 6 Resolve diffusion scheme
!! - 6.1 solve equations assuming atmosphere limiting
!! - 6.2 check if really atmosphere limiting
!! - 6.3 Reset the coefficient for diffusion (only used if resolv(ji) = .TRUE.)
!! - 6.4 resolve the equations with new boundary conditions if necessary
!! - 7 close the water balance
!! - 7.1 compute dr_ns with the bottom boundary condition 
!! - 7.2 compute total soil moisture content
!! - 7.3 deduction of upper flux from soil moisture variation and bottom flux
!! - 7.4 deduction of ae_ns and ru_ns:
!! - 8 Special treatment for the unstable cases
!! - 9 Then compute the temporary surface water and correct the outgoing runoff
!! - 10 smooth again
!! - 11 Optional computation of the fluxes 
!! - 12 We make some useful output
!! - 13 before closing the soil water, we check the water balance of soil
!! - 14 sum 3d variables into 2d variables
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil

  SUBROUTINE hydrol_soil (kjpindex, dtradia, veget_max, soiltile, njsc, reinf_slope, &
       & transpir, vevapnu, evapot, evapot_penm, runoff, drainage,   &
       & returnflow, reinfiltration, irrigation, &
       & tot_melt, evap_bare_lim, shumdiag, shumdiag_perma,   &
       & k_litt, litterhumdiag,  &
       & humrel,vegstress, drysoil_frac, &
!Isa, Chloe
      & stempdiag,snow, peatland, water2add_peat, ok_routage_peat,ok_no_drainage ,wt_soil,wt_soil2)
       
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex 
    ! input fields
    REAL(r_std), INTENT (in)                                 :: dtradia          !! Time step [s]
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: veget_max        !! Map of max vegetation types [-]
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc             !! indexing of PFT to soiltile
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile         !! fraction of PFT on each soil-hydrology tile [-]
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: transpir         !! Transpiration [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: reinf_slope      !! Slope coef
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow       !! Water returning to the deep reservoir [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration   !! Water returning to the top of the soil [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation       !! Irrigation [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot           !! Potential evaporation [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot_penm      !! Potential evaporation Penman [mm] 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt         !! snow melt [mm]
!Isa
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)       :: stempdiag 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: snow          !! Snow mass [Kg/m^2]
    !Chloe
    LOGICAL, INTENT(in), DIMENSION(kjpindex)                 :: peatland      !! logical peatland if dominant in a grid cell
    LOGICAL, INTENT(in)                                      :: ok_routage_peat     
    LOGICAL, INTENT(in)                                      :: ok_no_drainage

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff           !! complete runoff [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage         !! complete drainage [mm]
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: evap_bare_lim    !! limitation of bare soil evaporation on each 
                                                                                 !! soil column [mm] 
    REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (out)     :: shumdiag         !! Relative soil moisture
!Isa shumdiag_perma
    REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (out)     :: shumdiag_perma   
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: k_litt           !! Litter cond.
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: litterhumdiag    !! Litter humidity
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress        !! Veg. moisture stress (only for vegetation 
                                                                                 !! growth) 
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac     !! Function of the litter humidity, that will be 

!Chloe reinfiltration runoff des soiltile 1 � 3
    REAL(r_std)                                              :: reinfilt_watpeat !! Reinfilt runoff(jst=1:3) vers soiltile 4
    REAL(r_std)                                              :: fracreinfilt
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: water2add_peat 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)      :: wt_soil           !!Water Table position (mm)  
 REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)      :: wt_soil2           !!Water Table position (mm)  
!    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: wtold
    ! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu          !! 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)    :: humrel           !! Relative humidity [0-1, dimensionless]

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji, jv, jsl, jst, li       !! indices
    REAL(r_std), PARAMETER                         :: frac_mcs = 0.66            !! temporary depth
    !!??Aurelien: frac_ms inutile?
    REAL(r_std)                                    :: us_tmp,tmp_wt                     !! temporary stress
    REAL(r_std), DIMENSION(kjpindex)               :: temp                       !! temporary value for fluxes
    REAL(r_std), DIMENSION(kjpindex)               :: tmcold, tmcint             !!
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm)     :: moderwilt                  !!
    REAL(r_std), DIMENSION(kjpindex,nslm)          :: mcint                      !! To save mc values for future use
    LOGICAL, DIMENSION(kjpindex)                   :: is_under_mcr               !! Allows under residual soil moisture due to evap 
    LOGICAL, DIMENSION(kjpindex)                   :: is_over_mcs                !! Allows over saturated soil moisture due to 
                                                                                 !! returnflow 
    REAL(r_std), DIMENSION(kjpindex)               :: sum_rootsink               !! Sum of the root sink
    REAL(r_std), DIMENSION(kjpindex)               :: deltahum,diff              !!
    LOGICAL(r_std), DIMENSION(kjpindex)            :: test                       !!
    REAL(r_std), DIMENSION(kjpindex)               :: tsink                      !!
    REAL(r_std), DIMENSION(kjpindex)               :: water2extract              !! Temporary variable [mm]
    REAL(r_std), DIMENSION(kjpindex)               :: returnflow_soil            !! Water from the routing back to the bottom of 
                                                                                 !! the soil [mm] 
    REAL(r_std), DIMENSION(kjpindex)               :: reinfiltration_soil        !! Water from the routing back to the top of the 
                                                                                 !! soil [mm] 
    REAL(r_std), DIMENSION(kjpindex)               :: irrigation_soil            !! Water from irrigation returning to soil 
                                                                                 !! moisture for each soil type [mm] 
    REAL(r_std), DIMENSION(kjpindex)               :: flux_infilt                !!

!chloe
  REAL(r_std), DIMENSION (kjpindex,nstm) :: water2infilt_test
    !! 0 Arrays initialisation

    returnflow_soil(:) = zero
    reinfiltration_soil(:) = zero
    irrigation_soil(:) = zero
    qflux(:,:,:) = zero
    is_under_mcr(:) = .FALSE.
    is_over_mcs(:) = .FALSE.
    flux_infilt(:) = zero
    k(:,:)=zero
!Isa
 if (ok_freeze_cwrr) then
    kk(:,:,:)=zero
    kk_moy(:,:)=zero
  endif !if (ok_freeze_cwrr) then

if (ok_freeze_cwrr) then
!ISA 1
    !
    ! 1.1. calcul de la temp�rature sur l'�chelle hydro (nslm) sur la base du stempdiag (�chelle diagnostique nbdl) 
    !
	CALL calcule_temp_hydro(kjpindex, stempdiag, snow)
!end ISA 1
endif
    !
    ! split 2d variables to 3d variables, per soil type
    !
    CALL hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, transpir, humrel, evap_bare_lim)
    !
    ! Common variables
    !
    DO ji=1,kjpindex
       IF(vegtot(ji).GT.min_sechiba) THEN
          returnflow_soil(ji) = zero
          reinfiltration_soil(ji) = (returnflow(ji) + reinfiltration(ji))/vegtot(ji)
          irrigation_soil(ji) = irrigation(ji)/vegtot(ji)
       ELSE
          returnflow_soil(ji) = zero
          reinfiltration_soil(ji) = zero
          irrigation_soil(ji) = zero
       ENDIF
    ENDDO
       



    !
    !!_  for each soil type
    !
    DO jst = 1,nstm
       !
       !- We compute the sum of the sinks for future check-up
       sum_rootsink(:)=SUM(rootsink(:,:,jst),dim=2)
       DO ji=1,kjpindex
          tsink(ji) = sum_rootsink(ji)+MAX(ae_ns(ji,jst),zero)+subsinksoil(ji)
       ENDDO
       !
       ! The total moisture content (including water2infilt) is saved for balance checks at the end
       tmcold(:) = tmc(:,jst)

       !- The value of mc is kept in mcint, used in the flux computation after diffusion:
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
          ENDDO
       ENDDO

       
       DO ji = 1, kjpindex
          tmcint(ji) = dz(2,jst) * ( trois*mcint(ji,1) + mcint(ji,2) )/huit 
       ENDDO

       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmcint(ji) = tmcint(ji) + dz(jsl,jst) &
                  & * (trois*mcint(ji,jsl)+mcint(ji,jsl-1))/huit &
                  & + dz(jsl+1,jst) * (trois*mcint(ji,jsl)+mcint(ji,jsl+1))/huit
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          tmcint(ji) = tmcint(ji) + dz(nslm,jst) &
               & * (trois * mcint(ji,nslm) + mcint(ji,nslm-1))/huit
       ENDDO
    !! 1 We compare water2infilt and water2extract to keep only difference

       ! The bare soil evaporation is substracted to the soil moisture profile and first to the water available:  
       DO ji = 1, kjpindex
          water2extract(ji) = MIN(water2infilt(ji,jst) + irrigation_soil(ji) + reinfiltration_soil(ji), &
               & MAX(ae_ns(ji,jst),zero) + subsinksoil(ji))
       ENDDO
       ! First we substract from the surface
       DO ji = 1, kjpindex
            water2infilt(ji,jst) = water2infilt(ji,jst) + irrigation_soil(ji) + reinfiltration_soil(ji) - water2extract(ji)
       ENDDO       
       

       ! Then we update the water to extract from the soil
       DO ji = 1, kjpindex
          water2extract(ji) =  MAX(ae_ns(ji,jst),zero) + subsinksoil(ji) - water2extract(ji)
       ENDDO

       !! We add and substract components to the soil
       !! by filling the layers from top to bottom (same for reinfiltration) - this is done by smooth below
       !! 1.1 add to the first layer
       DO ji = 1, kjpindex
          mc(ji,1,jst) = mc(ji,1,jst)  &
      &         - water2extract(ji) * deux / dz(2,jst)
       ENDDO


       

       !!??Aurelien: here, the first layer can be oversaturated, thats why we need hydrol_soil_smooth
    !! 1.2 filling layers
       CALL hydrol_soil_smooth(kjpindex,jst, njsc, is_under_mcr, is_over_mcs)

   
    !! 2 Before diffusion scheme
   
    !! 2.1 Some initialisation necessary for the diffusion scheme to work

       DO ji = 1, kjpindex
          !- We correct rootsink for first two layers so that it is not too low in the first layer
          v1(ji,jst) = dz(2,jst)/huit * (trois * mc(ji,1,jst)+ mc(ji,2,jst))
!Isa :contenu en eau de la couche 1. on transf�re si n�cessaire une partie du rootsink de la couche 1, exc�dentaire par rapport 
!au contenu en eau pr�sent de cette couche,sur le rootsink de la couche 2.
          rootsink(ji,2,jst) = rootsink(ji,2,jst) + MAX(rootsink(ji,1,jst)-v1(ji,jst), zero) 
          rootsink(ji,1,jst) = MIN(rootsink(ji,1,jst),v1(ji,jst))
          !- estimate maximum evaporation flux in mm/step, assuming the water is available
          flux(ji) = zero
          IF(vegtot(ji).GT.min_sechiba .AND. .NOT. is_under_mcr(ji) ) THEN
             !- Flux = evapot_penm   if frac_bare_ns > min_sechiba
             !-      = zero          else
             flux(ji) = evapot_penm(ji) * &
                  & AINT(frac_bare_ns(ji,jst)+un-min_sechiba)    ! fonction PARTIE ENTIERE         
          ENDIF
       ENDDO

       ! Then we prepare the infiltration (first the irrigation and below mcr fill up)
       ! Initialise the flux to be infiltrated 
      
       
       !Chloe 21/03/14 et 30/07/14
       !d�ja sur une boucle de jst
       reinfilt_watpeat = 0. 
       IF (ok_reinfilt_peat .AND. (jst .EQ. nstm) ) THEN
           DO ji = 1, kjpindex
              IF ( soiltile(ji,4) .GT. zero) THEN
                  !Reinfiltrate runoff to peat soiltile 4 in function of peat fraction (gerhard's 2003 methods)
                  !fracreinfilt = soiltile(ji,4) ! Comme Krinner 2003
                  ! fracreinfilt = 1. ! la totalite
                  !reinfilt_watpeat= 0
                  ! ru_ns(ji,1:3) = ru_ns(ji,1:3) * (1.-fracreinfilt)                
                  reinfilt_watpeat = ru_ns(ji,1)*soiltile(ji,1) + ru_ns(ji,2)*soiltile(ji,2) + ru_ns(ji,3)*soiltile(ji,3)
                  !IF (soiltile(ji,4) .LT. 0.01) THEN
                  !    fracreinfilt = soiltile(ji,4)
                  !    ru_ns(ji,1:3) = ru_ns(ji,1:3) * (1.-fracreinfilt) 
                  !    water2infilt(ji,4)= water2infilt(ji,4) + (reinfilt_watpeat*fracreinfilt)/soiltile(ji,4)
                  !ELSE
                      ru_ns(ji,1:3)=0
                      water2infilt(ji,4)= water2infilt(ji,4) + reinfilt_watpeat/soiltile(ji,4)
                  !ENDIF
              ENDIF

           ENDDO
       ENDIF
       !
        !
       !Chloe++ R�infiltration de l'eau pr�sente dans le r�servoir stagnant
       DO ji = 1, kjpindex
      
         IF ( ok_stagnant ) THEN
           water2infilt(ji,4) = water2infilt(ji,4) + stagnant(ji,4)
           stagnant(ji,4)=0.
         ENDIF 
         flux_infilt(ji) = water2infilt(ji,jst) + precisol_ns(ji,jst) 
         !- The incoming flux is also first dedicated to fill the soil up to mcr (in case needed)



        ENDDO
       !Chloe--


       DO jsl = 1, nslm
          WHERE (is_under_mcr(:))
             mc(:,jsl,jst) = mc(:,jsl,jst) + flux_infilt(:) / (dpu_max*mille)
          ENDWHERE
       END DO
       WHERE (is_under_mcr(:))
          flux_infilt(:) = zero
       ENDWHERE



    !! 2.2 coefficients are computed for the profile of mc before infiltration:
       CALL hydrol_soil_coef(kjpindex,jst,njsc)
!do ji=1, kjpindex
!enddo
    !! 3 The infiltration is computed 
       CALL hydrol_soil_infilt(kjpindex, jst, dtradia, njsc, flux_infilt)
!do ji=1, kjpindex
!enddo
    !! 4 Coefficient are recomputed for the profile of mc after infiltration:
       CALL hydrol_soil_coef(kjpindex,jst,njsc)


    !! 5 Prepar the diffusion scheme
 
    !! 5.1 Set the values for diffusion scheme
       CALL hydrol_soil_setup(kjpindex,jst,dtradia)

    !! 5.2 verifications for a good soil humidity
    ! We only run the scheme in case we are not under mcr after precip and various reinfiltrations and not over mcs
       resolv(:) = (mask_soiltile(:,jst) .GT. 0) .AND. & 
               & (.NOT. is_under_mcr(:)) .AND. (.NOT. is_over_mcs(:))    

       ! In oversaturated case, we first take the evaporation from the outgoing water (the rest will be taken from the soil)
       sum_rootsink(:)=SUM(rootsink(:,:,jst),dim=2)
       WHERE (is_over_mcs(:))
          mc(:,1,jst) = mc(:,1,jst) - flux(:) * deux / dz(2,jst)
       ENDWHERE
       DO jsl=1, nslm
          WHERE (is_over_mcs(:))
             mc(:,jsl,jst) = mc(:,jsl,jst) - sum_rootsink(:) / (dpu_max*mille) 
          ENDWHERE
       ENDDO
       ! In under residual case, we equally spread the transpiration over the layers
       DO jsl = 1, nslm
          WHERE (is_under_mcr(:))
             mc(:,jsl,jst) = mc(:,jsl,jst) - sum_rootsink(:) / (dpu_max*mille) 
          ENDWHERE
       ENDDO
!Isa
if (ok_freeze_cwrr) then
	do ji =1, kjpindex
	    do jsl = 1, nslm
		mcl(ji,jsl,:)= MIN(mc(ji,jsl,:),mcr(njsc(ji))+(1-profil_froz_hydro_ns(ji,jsl, :))*(mc(ji,jsl,:)-mcr(njsc(ji))))
                !borne inf�rieure mc pour les cas under_mcr
        enddo
	enddo
else !if (ok_freeze_cwrr) then
   mcl(:,:,:)=mc(:,:,:)
endif !if (ok_freeze_cwrr) then
!Isa : in the following, mc -> mcl


       !Chloe++ mcs(njsc(ji))?
        if (ok_sat30cm ) then
            do ji=1, kjpindex
                do jsl=9,nslm
                  mc(ji,jsl,4)=mcs(njsc(ji))
                IF (peatland(ji) ) mc(:,jsl,:)=mcs_peat !Chloe

                enddo
            enddo
        endif

        if (ok_routage_peat ) then
            do ji=1, kjpindex
                DO jsl=1,7
                    mc_obj_peat(ji,jsl,4)=mc(ji,jsl,4)
                ENDDO 
             
                do jsl=8,nslm
                  mc_obj_peat(ji,jsl,4)=mcs(njsc(ji))
                IF (peatland(ji) ) mc(:,jsl,:)=mcs_peat !Chloe

                enddo
            enddo
        endif

        !Test Desert mousson :
        !DO ji=1, kjpindex
        !    DO jsl=1,4
        !    mc(ji,jsl,:)=mcs(njsc(ji))
        !    ENDDO
        !ENDDO




       !Chloe--



    !! 5.3 compute matrix coefficients
       !- First layer
       DO ji = 1, kjpindex
          ! We directly use the value of k (a and b are only used un hydrol_soil_coef to compute k)
          !-        First layer
          !-        Isa : ici flux=Epot, en mm/tstep
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               &  -flux(ji) - (b(ji,1)+b(ji,2))/deux*(dtradia/one_day) - rootsink(ji,1,jst)
       ENDDO

       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  & +  gp(ji,jsl) * mcl(ji,jsl+1,jst) & 
                  & + (b(ji,jsl-1) - b(ji,jsl+1))  & 
                  & * (dtradia/one_day) / deux & 
                  & - rootsink(ji,jsl,jst) 
          ENDDO
       ENDDO
       
       !- Last layer
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
!Isa : introduction du free_drain_coef, mis � 1, en attente de correction
 		& + (b(ji,jsl-1)+b(ji,jsl) & 
                & *(1.-2.*free_drain_coef(ji,jst))) & 
                & /deux * (dtradia/one_day) &
                & - rootsink(ji,jsl,jst)

       ENDDO

       !- store the equations in case needed again
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3) 
          ENDDO
       ENDDO

    !! 6 Resolve diffusion scheme

    !! 6.1 solve equations assuming atmosphere limiting
    ! (sufficient soil moisture for Evapot)
       CALL hydrol_soil_tridiag(kjpindex,jst)


    !! 6.2 check if really atmosphere limiting
       DO ji = 1, kjpindex
          !
          !- Prepare to rerun in case of under residual with evaporation 
          !- Isa : pr�cision : avec les mc non modifi�s par le hydrol_soil_tridiag pr�c�dent !
          !ISA 3: mc -> mcl
          !-
          resolv(ji) = (mcl(ji,1,jst).LT.(mcr(njsc(ji))).AND.flux(ji).GT.min_sechiba)
          !end ISA 3
       ENDDO

    !! 6.3 Reset the coefficient for diffusion (only used if resolv(ji) = .TRUE.)
       DO jsl=1,nslm
          !- The new condition is to put the upper layer at residual soil moisture
          DO ji = 1, kjpindex
                rhs(ji,jsl) = srhs(ji,jsl)
                tmat(ji,jsl,1) = stmat(ji,jsl,1)
                tmat(ji,jsl,2) = stmat(ji,jsl,2)
                tmat(ji,jsl,3) = stmat(ji,jsl,3)
             END DO
          END DO

       DO ji = 1, kjpindex
          tmat(ji,1,2) = un
          tmat(ji,1,3) = zero
          rhs(ji,1) = mcr(njsc(ji))

       ENDDO

    !! 6.4 resolve the equations with new boundary conditions if necessary
       CALL hydrol_soil_tridiag(kjpindex,jst)
!       write(*,*) 'mc(1,37,1)=',mcl(1,37,1)
!Isa


if (ok_freeze_cwrr) then
	do ji =1, kjpindex
	    do jsl = 1, nslm
		mc(ji,jsl,:)=MAX(mcl(ji,jsl,:), mcl(ji,jsl,:)+profil_froz_hydro_ns(ji,jsl,:)*(mc(ji,jsl,:)-mcr(njsc(ji))))
        enddo
	enddo
else
   mc(:,:,:)=mcl(:,:,:)
endif


!       !Chloe++
        if (ok_sat30cm ) then
            do ji=1, kjpindex
                do jsl=9,nslm
                  mc(ji,jsl,4)=mcs(njsc(ji))
                 IF (peatland(ji) ) mc(:,jsl,:)=mcs_peat !Chloe
                enddo
            enddo
        endif

        if (ok_routage_peat ) then
            do ji=1, kjpindex
                DO jsl=1,7
                    mc_obj_peat(ji,jsl,4)=mc(ji,jsl,4)
                ENDDO                

                DO jsl=8,nslm
                  mc_obj_peat(ji,jsl,4)=mcs(njsc(ji))
                  IF (peatland(ji) ) mc(:,jsl,:)=mcs_peat !Chloe
                ENDDO
            enddo
        endif

!       !Chloe--




    !! 7 close the water balance
    !! 7.1 compute dr_ns with the bottom boundary condition 
    !initialize qflux at bottom of diffusion and avoid over saturated or under residual soil moisture 

       DO ji = 1, kjpindex
          dr_ns(ji,jst)=zero
          jsl=nslm
          IF (.NOT. is_under_mcr(ji)) THEN
           ! IF ( .NOT. ok_no_drainage .AND. .NOT. jst .EQ. 4 ) THEN ! Chloe+ pas de drainage pour les peatland 04/07/13, only soiltype 4
!Isa : correction pr ajout de free_drain_coef
             dr_ns(ji,jst) = mask_soiltile(ji,jst)*k(ji,jsl) * (dtradia/one_day)*free_drain_coef(ji,jst)
           ! ENDIF
          ENDIF
          IF (ok_no_drainage) THEN
              dr_ns(ji,4)=zero
          ENDIF
       ENDDO

       !Isa precision : si on est en under_mcr avant la r�sol pr�c�dente on ne
       !'tridiag' pas donc pas de calcul de drainage.          

    !! 7.2 compute total soil moisture content
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit 
       ENDDO
!       write(*,*) 'hydrol 3865: tmc(37,1)=',tmc(37,1)
!       write(*,*) 'mc(37,1:2,1)=',mc(37,1:2,1)
!       write(*,*) 'dz(2,1)=',dz(2,1)
    

       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) &
                  & * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO


       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm,jst) &
               & * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO



       DO ji = 1, kjpindex
    !! 7.3 deduction of upper flux from soil moisture variation and bottom flux
          qflux00(ji,jst) = mask_soiltile(ji,jst) * &
               & (MIN(tmcs(ji,jst),tmc(ji,jst))-tmcint(ji)+SUM(rootsink(ji,:,jst))+dr_ns(ji,jst)-returnflow_soil(ji))

    !! 7.4 deduction of ae_ns and ru_ns:
    ! ae_ns+ru_ns=precisol_ns+irrigation-q0  
    ! (returnflow introduit par le haut) !!??Aurelien: are you sure
          ! deduction of ae_ns and ru_ns:
          ! ae_ns+ru_ns=precisol_ns+irrigation-q0 
          ! Isa: qflux00 = tt ce qui est rentr� dans le sol entre la fin et le
          ! d�but de la routine. Donc:
          ! qflux00+ae_ns+runoff=precisol_ns+irrigation_soil+reinfiltration_soil

          ae_ns(ji,jst) = MAX(MIN((precisol_ns(ji,jst) &
               & +water2infilt(ji,jst)-water2extract(ji)-qflux00(ji,jst)) ,evapot_penm(ji)),zero) * mask_soiltile(ji,jst)
          ru_ns(ji,jst) = (precisol_ns(ji,jst)  &
              & +water2infilt(ji,jst)-water2extract(ji)-qflux00(ji,jst)-ae_ns(ji,jst)) * mask_soiltile(ji,jst)

       ENDDO !DO ji = 1, kjpindex

!Chloe Gel du sol = pas de runoff :
!IF (ok_freeze_cwrr) THEN
!    WHERE (temp_hydro(:,1) .LT. (tzero - deux))
!        ru_ns(:,jst)=zero
!    ENDWHERE
!ENDIF

!!! Chloe 10:10/14 
!!! On impose le passage du reservoir avant de ruisseler
!! Mais apr�s avoir calcul� le ruissellement quand meme, �vident non ?
    ! Chloe++ Stockage du r�servoir d'eau stagnante                                                                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DO ji = 1, kjpindex
!       !We already are in loop soiltile
!       IF ( ok_stagnant .AND. jst .EQ. nstm  ) THEN
!            stagnant(ji,jst) = stagnant(ji,jst) + ru_ns(ji,jst)
!            ru_ns(ji,jst) = 0.
!    
!         IF ( stagnant(ji,jst) .GT. max_stagnant ) THEN
!            ru_ns(ji,jst) = stagnant(ji,jst) - max_stagnant
!            stagnant(ji,jst) = max_stagnant
!         ENDIF
!
!       ELSE
!          stagnant(:,jst) = zero ! Valeur non physique
!       ENDIF
!   ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!
             
    !! 8 Special treatment for the unstable cases
    ! when boundary condition mc1=mcr leads to negative runoff

!    write(*,*) 'hydrol 3599: ru_ns=',ru_ns(1,1)
!    write(*,*) 'ru_ns(37,1)=',ru_ns(37,1)
!    write(*,*) 'precisol_ns(37,1)=',precisol_ns(37,1)
!    write(*,*) 'water2infilt(37,1)=',water2infilt(37,1)
!    write(*,*) 'water2extract(37)=',water2extract(37)
!    write(*,*) 'qflux00(37,1)=',qflux00(37,1)
!    write(*,*) 'ae_ns(37,1)=',ae_ns(37,1)
!    write(*,*) 'mask_soiltile(37,1)=',mask_soiltile(37,1)
!    write(*,*) 'tmcs(37,1)=',tmcs(37,1)
!    write(*,*) 'tmc(37,1)=',tmc(37,1)
!    write(*,*) 'tmcint(37)=',tmcint(37)
!    write(*,*) 'rootsink(37,:,1)=',rootsink(37,:,1)
!    write(*,*) 'dr_ns(37,1)=',dr_ns(37,1)
!    write(*,*) 'returnflow_soil(37)=',returnflow_soil(37)
       IF (long_print) THEN
          DO ji = 1, kjpindex
             IF (ru_ns(ji,jst).LT.-min_sechiba) THEN
                WRITE (numout,*) 'Negative runoff corrected', ji,jst,ru_ns(ji,jst), mc(ji,1,jst), tmc(ji,jst)
             ENDIF
          ENDDO
       ENDIF
!       write(*,*) 'hydrol 3910'

       DO ji = 1, kjpindex
          temp(ji) = MIN(ru_ns(ji,jst),zero)
       ENDDO
       DO ji = 1, kjpindex
          ru_ns(ji,jst) = ru_ns(ji,jst) - temp(ji)
          ! We correct this by taking water from the whole soil
          qflux00(ji,jst) = qflux00(ji,jst) + temp(ji)
       ENDDO
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mc(ji,jsl,jst) = mc(ji,jsl,jst) + temp(ji) / (dpu_max*mille)
          ENDDO
       ENDDO



       ! Avoid under-precision value for the 3 outward flux
       DO ji = 1, kjpindex  
          IF (ABS(ae_ns(ji,jst)).LT.min_sechiba) THEN
             ae_ns(ji,jst) = zero
          ENDIF

          IF(ABS(ru_ns(ji,jst)).LT.min_sechiba) THEN
             ru_ns(ji,jst) = zero
          ENDIF

          IF(ABS(dr_ns(ji,jst)).LT.min_sechiba) THEN
             dr_ns(ji,jst) = zero
          ENDIF
       ENDDO




    !! 9 Then compute the temporary surface water and correct the outgoing runoff
       !IF ( .NOT. doponds ) THEN
        !Chloe : 
        ! if stagnant : wat2inf=0 ici, water2infilt prend les quantit�s de runoff
        ! total auparavant sans la prise en compte de la pente.
        ! Cette m�thode permet de ne pas avoir un tmc > tmcs = 860kg/m2 si
        ! mcs=0.43
        IF ( .NOT. doponds .AND. .NOT. ok_stagnant) THEN
          DO ji = 1, kjpindex
             water2infilt(ji,jst) = reinf_slope(ji) * ru_ns(ji,jst)
          ENDDO
       ELSE
          DO ji = 1, kjpindex           
             water2infilt(ji,jst) = zero
          ENDDO
       ENDIF
       !
       DO ji = 1, kjpindex
       !Chloe add :          
       !IF (ru_ns(ji,jst) .GT. zero) THEN 
          ru_ns(ji,jst) = ru_ns(ji,jst) - water2infilt(ji,jst)
       !ENDIF   
          water2infilt(ji,jst) = water2infilt(ji,jst) + ae_ns(ji,jst)
      END DO

!Chloe : mettre le ok_stagnant apr�s les calculs et corrections appliqu�s au
!runoff : 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !DO ji = 1, kjpindex
!  !     !We already are in loop soiltile
   !   IF ( ok_stagnant .AND. jst .EQ. nstm  ) THEN
   !         stagnant(ji,jst) = stagnant(ji,jst) + ru_ns(ji,jst)
   !         ru_ns(ji,jst) = zero
   !  
   !      IF ( stagnant(ji,jst) .GT. max_stagnant .AND. jst .EQ. nstm ) THEN
   !            ru_ns(ji,jst) = stagnant(ji,jst) - max_stagnant
   !            stagnant(ji,jst) = max_stagnant
   !      ENDIF
   !            !stagnant(ji,4) = MIN(stagnant(ji,4) + ru_ns(ji,4),max_stagnant)
   !            !ru_ns(ji,4)= MAX(stagnant(ji,4) - max_stagnant, zero)
   !   ENDIF
   !ENDDO

DO ji = 1, kjpindex
!       !We already are in loop soiltile
      IF ( ok_stagnant .AND. jst .EQ. nstm  ) THEN
         stagnant(ji,jst) = stagnant(ji,jst) + ru_ns(ji,jst)
         ru_ns(ji,jst)= MAX(stagnant(ji,jst) - max_stagnant, zero)
         stagnant(ji,jst) = MIN(stagnant(ji,jst),max_stagnant)
    ENDIF
ENDDO

!DO ji=1,kjpindex
!IF ( ok_stagnant .AND. jst .EQ. nstm  ) THEN
!    IF (stagnant(1,4) .LT. max_stagnant .AND. ru_ns(1,4) .GT. zero) THEN
!       write(*,*) 'Chloe 4 BUG BUG', stagnant(1,4), ru_ns(1,4)
!    ENDIF 
!    IF (stagnant(ji,jst) .LT. max_stagnant .AND. ru_ns(ji,jst) .GT. zero) THEN
!       write(*,*) 'Chloe ji,jst BUG BUG', stagnant(ji,jst), ru_ns(ji,jst)
!    ENDIF
!    
!ENDIF 
!ENDDO

!!!!!!!CHLOE 30 JUILLET : Mauvaise place : doit etre avant call
!hydrol_soil_infilt
       !Chloe 210314
       !
     !  reinfilt_watpeat = zero
     !  IF (ok_reinfilt_peat .AND. (jst .EQ. nstm) ) THEN
     !      DO ji = 1, kjpindex
     !            IF ( soiltile(ji,4) .GT. zero) THEN
     !               !Reinfiltrate runoff to peat soiltile 4
     !               ! in function of peat fraction (gerhard's 2003 methods)
     !             reinfilt_watpeat = ru_ns(ji,1)*soiltile(ji,1) + ru_ns(ji,2)*soiltile(ji,2) + ru_ns(ji,3)*soiltile(ji,3) 
     !              ru_ns(ji,1:3)=0
     !             water2infilt(ji,4)= water2infilt(ji,4) + reinfilt_watpeat/soiltile(ji,4)
     ! 
     !           ENDIF
     !      ENDDO
     !  ENDIF
     !  !
     !  !Chloe end
  

    !! 10 smooth again
    !!??Aurelien: Why?
!       write(*,*) 'hydrol 3960'
       CALL hydrol_soil_smooth(kjpindex, jst, njsc, is_under_mcr, is_over_mcs)


    !! 11 Optional computation of the fluxes 
       IF ( check_cwrr ) THEN
          CALL hydrol_soil_flux(kjpindex,jst,mcint,returnflow_soil)
       ENDIF

    !! 12 We make some useful output
    !- Total soil moisture, soil moisture at litter levels, soil wetness...

       !-total soil moisture:
       DO ji=1,kjpindex
          tmc(ji,jst)= dz(2,jst) * (trois*mc(ji,1,jst) + mc(ji,2,jst))/huit
       END DO

       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) * ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO

       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) +  dz(nslm,jst) * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       END DO
      
       ! the litter is the 4 top levels of the soil
       ! we compute various field of soil moisture for the litter (used for stomate and for albedo)

       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst)+ mc(ji,2,jst))/huit
          tmc_litter(ji,jst) = tmc_litter(ji,jst) 
       END DO

       ! sum from level 1 to 4

       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl,jst) * & 
                  & ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO

       
       ! subsequent calcul of soil_wet_litter (tmc-tmcw)/(tmcf-tmcw)

       DO ji=1,kjpindex
          soil_wet_litter(ji,jst) = MIN(un, MAX(zero,&
               & (tmc_litter(ji,jst)-tmc_litter_wilt(ji,jst)) / &
               & (tmc_litter_field(ji,jst)-tmc_litter_wilt(ji,jst)) ))
       END DO

       ! Soil wetness profiles (mc-mcw)/(mcs-mcw)
       ! soil_wet is the ratio of soil moisture to available soil moisture for plant
       ! (ie soil moisture at saturation minus soil moisture at wilting point).

       DO ji=1,kjpindex
          soil_wet(ji,1,jst) = MIN(un, MAX(zero,&
               & (trois*mc(ji,1,jst) + mc(ji,2,jst) - quatre*mcw(njsc(ji)))&
               & /(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
          humrelv(ji,1,jst) = zero

       END DO

       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             soil_wet(ji,jsl,jst) = MIN(un, MAX(zero,&
                  & (trois*mc(ji,jsl,jst) + & 
                  & mc(ji,jsl-1,jst) *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & + mc(ji,jsl+1,jst)*(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & - quatre*mcw(njsc(ji))) / (quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
              

          END DO
       END DO

       DO ji=1,kjpindex
          soil_wet(ji,nslm,jst) = MIN(un, MAX(zero,&
               & (trois*mc(ji,nslm,jst) &
               & + mc(ji,nslm-1,jst)-quatre*mcw(njsc(ji)))/(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
            

       END DO

       ! we compute the moderation of transpiration due to wilting point
       ! moderwilt is a factor which is zero if soil moisture is below the wilting point
       ! and is un if soil moisture is above the wilting point.


       DO jsl=1,nslm
          DO ji=1,kjpindex
             moderwilt(ji,jsl,jst) = INT( MAX(soil_wet(ji,jsl,jst), zero) + un - min_sechiba )
          END DO
       END DO

       ! we compute the new humrelv to use in sechiba:
       ! loop on each vegetation type
!       write(*,*) 'hydrol 4054'
       humrelv(:,1,jst) = zero   

       ! calcul of us for each layer and vegetation type.
       DO jv = 2,nvm
          DO ji=1,kjpindex
             !- Here we make the assumption that roots do not take water from the 1st layer. 
             !- Comment the us=0 if you want to change this.
!             us(ji,jv,jst,1) = moderwilt(ji,1,jst)*MIN(un,((trois*mc(ji,1,jst) + mc(ji,2,jst)) &
!                  & /(quatre*mcs(jst)*pcent(jst))) )* (un-EXP(-humcste(jv)*dz(2,jst)/mille/deux)) &
!                  & /(un-EXP(-humcste(jv)*zz(nslm,jst)/mille))
             us(ji,jv,jst,1) = zero
             humrelv(ji,jv,jst) = MAX(us(ji,jv,jst,1),zero)
          END DO
       ENDDO

       DO jsl = 2,nslm-1
          DO jv = 2, nvm
             DO ji=1,kjpindex
                !
!                us_tmp= (trois*mc(ji,jsl,jst) + &
!                     & mc(ji,jsl-1,jst)*(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))+ &
!                     & mc(ji,jsl+1,jst)*(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))) &
!                     & /(quatre*mcs(njsc(ji))*pcent(njsc(ji)))
                ! us_tmp should not be negative
                !!??Aurelien: Useless, here we are sure mc between mcr and mcs
!                us_tmp=MAX(us_tmp, zero)
                ! us is computed with a SQRT in order for it to grow more rapidly with soil moisture.
                ! it is not essential
!Isa : influence du gel dans le sol

                     if (ok_freeze_cwrr) then

                       us_tmp= (trois*(1-profil_froz_hydro_ns(ji,jsl, jst))*mc(ji,jsl,jst)+ &
                     & (1-profil_froz_hydro_ns(ji,jsl-1, jst))*mc(ji,jsl-1,jst) &
                     &  *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))+ &
                     & (1-profil_froz_hydro_ns(ji,jsl+1, jst))*mc(ji,jsl+1,jst) &
                     &  *(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))) &
                     & /(quatre*mcs(njsc(ji))*pcent(njsc(ji)))


                     else !if (ok_freeze_cwrr) then

                        us_tmp= (trois*mc(ji,jsl,jst) + &
                     & mc(ji,jsl-1,jst) &
                     &  *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))+ &
                     & mc(ji,jsl+1,jst) &
                     &  *(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))) &
                     & /(quatre*mcs(njsc(ji))*pcent(njsc(ji)))


                     endif !if (ok_freeze_cwrr) then

                ! us_tmp should not be negative
                !!??Aurelien: Useless, here we are sure mc between mcr and mcs
                us_tmp=MAX(us_tmp, zero)
                ! us is computed with a SQRT in order for it to grow more rapidly with soil moisture.
                ! it is not essential
                us(ji,jv,jst,jsl) = moderwilt(ji,jsl,jst) * MIN( un, SQRT(us_tmp) ) * nroot(jv,jst,jsl)
                !

                !Chloe++
                if(ok_sat30cm) then
                    if (jsl.gt.8) then
                      us(ji,jv,4,jsl)=1.
                    endif
                endif


                !Chloe--

                us(ji,jv,jst,jsl) = MAX(us (ji,jv,jst,jsl), zero)

                humrelv(ji,jv,jst) = MAX((humrelv(ji,jv,jst) + us(ji,jv,jst,jsl)),zero)
             END DO
          END DO
       ENDDO



    ! for jsl = nslm :
       DO jv = 2, nvm
          DO ji=1,kjpindex

             if (ok_freeze_cwrr) then
                us_tmp = (trois*(1-profil_froz_hydro_ns(ji,nslm, jst))*mc(ji,nslm,jst)  &
                  & + (1-profil_froz_hydro_ns(ji,nslm-1, jst))*mc(ji,nslm-1,jst))  &
                  & / (quatre*mcs(njsc(ji))*pcent(njsc(ji)))


             else !if (ok_freeze_cwrr) then  
                us_tmp = (trois*mc(ji,nslm,jst)  &
                  & + mc(ji,nslm-1,jst))  &
                  & / (quatre*mcs(njsc(ji))*pcent(njsc(ji)))



             endif !if (ok_freeze_cwrr) then
             !
             ! us_tmp should not be negative
             us_tmp=MAX(us_tmp, zero)
             !
             us(ji,jv,jst,nslm) =moderwilt(ji,nslm,jst)*  &
                  & MIN( un, SQRT(us_tmp) ) * nroot(jv,jst,nslm)

             us(ji,jv,jst,nslm) = MAX(us(ji,jv,jst,nslm), zero)



            !Chloe++
             if(ok_sat30cm) then
                us(ji,jv,4,nslm)=1.
             endif

            !Chloe--


             humrelv(ji,jv,jst) = MAX(zero,MIN(un, humrelv(ji,jv,jst) &
                  & + us(ji,jv,jst,nslm)))
             vegstressv(ji,jv,jst) = humrelv(ji,jv,jst)
          END DO
       END DO

       DO jv = 2, nvm
          DO ji = 1, kjpindex
             IF (corr_veg_soil(ji,jv,jst) .LT. min_sechiba) THEN
                humrelv(ji,jv,jst) = zero  
             ENDIF
          END DO
       END DO
   
    !! 13 before closing the soil water, we check the water balance of soil

       IF(check_cwrr) THEN
          DO ji = 1,kjpindex
             deltahum(ji) = (tmc(ji,jst) - tmcold(ji))
             diff(ji)     = precisol_ns(ji,jst)-ru_ns(ji,jst)-dr_ns(ji,jst)-tsink(ji) &
                  & + irrigation_soil(ji) + returnflow_soil(ji) + reinfiltration_soil(ji)

             test(ji) = (abs(deltahum(ji)-diff(ji))*mask_soiltile(ji,jst) .GT. allowed_err)
          ENDDO
             
          DO ji = 1,kjpindex
             IF(test(ji)) THEN
                
                WRITE (numout,*)'CWRR pat: bilan non nul',ji,jst,njsc(ji),deltahum(ji)-diff(ji)
                WRITE (numout,*)'tmc,tmcold,diff',tmc(ji,jst),tmcold(ji),deltahum(ji)
                WRITE(numout,*) 'evapot,evapot_penm,ae_ns',evapot(ji),evapot_penm(ji),ae_ns(ji,jst)
                WRITE (numout,*)'flux,ru_ns,qdrain,tsink,q0,precisol',flux(ji),ru_ns(ji,jst), &
                     &      dr_ns(ji,jst),tsink(ji),qflux00(ji,jst),precisol_ns(ji,jst)
                WRITE (numout,*)'water2infilt',water2infilt(ji,jst)
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                WRITE (numout,*)'irrigation, returnflow, reinfiltration', &
                     & irrigation_soil(ji),returnflow_soil(ji),reinfiltration_soil(ji)
                WRITE (numout,*)'mc',mc(ji,:,jst)
                WRITE (numout,*)'qflux',qflux(ji,:,jst)
                WRITE (numout,*)'veget_max', veget_max(ji,:)
                WRITE (numout,*)'k', k(ji,:)
                waterbal_error=.TRUE.
                CALL ipslerr(2, 'hydrol_soil', 'We will STOP after hydrol_waterbal.',&
                     & 'CWRR water balance check','')
                
             ENDIF
          ENDDO
!          write(*,*) 'hydrol 4206'

          DO ji = 1,kjpindex
             
             IF(MINVAL(mc(ji,:,jst)).LT.-min_sechiba) THEN
                WRITE (numout,*)'CWRR MC NEGATIVE', &
                     & ji,lalo(ji,:),MINLOC(mc(ji,:,jst)),jst,mc(ji,:,jst)
                WRITE (numout,*)'evapot,evapot_penm,ae_ns',evapot(ji),evapot_penm(ji),ae_ns(ji,jst)
                WRITE (numout,*)'flux,ru_ns,qdrain,tsink,q0,precisol',flux(ji),ru_ns(ji,jst), &
                     &      dr_ns(ji,jst),tsink(ji),qflux00(ji,jst),precisol_ns(ji,jst)
                WRITE (numout,*)'water2infilt',water2infilt(ji,jst)
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                WRITE (numout,*)'irrigation, returnflow, reinfiltration', &
                     & irrigation_soil(ji),returnflow_soil(ji),reinfiltration_soil(ji)
                WRITE (numout,*)'mc',mc(ji,:,jst)
                WRITE (numout,*)'qflux',qflux(ji,:,jst)
                WRITE (numout,*)'veget_max', veget_max(ji,:)
                WRITE (numout,*)'k', k(ji,:)
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                waterbal_error=.TRUE.
                CALL ipslerr(2, 'hydrol_soil', 'We will STOP after hydrol_waterbal.',&
                     & 'CWRR MC NEGATIVE','')
             ENDIF
          END DO

          DO ji=1,kjpindex
             IF (ru_ns(ji,jst)*soiltile(ji,jst).LT.-min_sechiba) THEN
                WRITE (numout,*) 'Negative runoff', ji,jst, mask_soiltile(ji,jst) 
                WRITE (numout,*) 'mc1, mc2', mc(ji,1,jst), mc(ji,2,jst)
                WRITE (numout,*) 'mcint1, mcint2', mcint(ji,1), mcint(ji,2)
                WRITE (numout,*) 'qflux1, flux', qflux(ji,nslm,jst), flux(ji)
                WRITE (numout,*) 'is_over_mcs, is_under_mcr, test', &
                     & is_over_mcs(ji), is_under_mcr(ji), tmc(ji,jst)-tmcint(ji)+qflux(ji,nslm,jst)+SUM(rootsink(ji,:,jst))
                WRITE (numout,*)'mc', mc(ji,:,jst)
                WRITE (numout,*)'mcint', mcint(ji,:)
                WRITE (numout,*)'qflux', qflux(ji,:,jst)
                WRITE (numout,*)'rootsink1,evapot_penm,vegtot', rootsink(ji,1,jst), evapot_penm(ji), vegtot(ji)
                WRITE (numout,*) 'ae_ns, tsink, returnflow, reinfiltration, precisol_ns, irrigation, qflux0, ru_ns', &
                     & ae_ns(ji,jst), tsink(ji), returnflow_soil(ji), reinfiltration_soil(ji), &
                     & precisol_ns(ji,jst), irrigation_soil(ji), qflux00(ji,jst), ru_ns(ji,jst)
                waterbal_error=.TRUE.
                CALL ipslerr(2, 'hydrol_soil', 'We will STOP after hydrol_waterbal.',&
                     & 'Negative runoff, non-saturated soil','')
             ENDIF
          ENDDO
       ENDIF

       IF (long_print) WRITE (numout,*) ' hydrol_soil done for jst =', jst     
!Isa
if (ok_freeze_cwrr) then
do ji = 1, kjpindex
kk_moy(ji,:) =kk_moy(ji,:)+soiltile(ji,jst)*k(ji,:) 
kk(ji,:,jst)=k(ji,:)
enddo
endif !if (ok_freeze_cwrr) then


    END DO  ! end of loop on soiltile begin Line 3859
!    write(*,*) 'hydrol 4263'

    !
    !! 14 sum 3d variables into 2d variables
    !
    CALL hydrol_diag_soil (kjpindex, veget_max, soiltile, njsc, runoff, drainage, &
         & evap_bare_lim, evapot, vevapnu, returnflow, reinfiltration, irrigation, &
         & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, &
         & drysoil_frac,tot_melt, water2add_peat, ok_routage_peat)
!    write(*,*) 'hydrol 4271'


!!!Chloe++
            DO ji=1,kjpindex
                DO jst=1,nstm
                   wt_soil2(ji,jst)=0.0_r_std
                      DO jsl=1,nslm
                        !tmp_wt=MIN(mc(ji,jsl,jst)/(0.98*mcs(njsc(ji))),un)
                         tmp_wt=MIN(MAX(zero,(mc(ji,jsl,jst)-mcr(njsc(ji)))/(mcs(njsc(ji))-mcr(njsc(ji)))),un)
                         wt_soil2(ji,jst)=wt_soil2(ji,jst) + tmp_wt*dz(jsl,jst)
                      ENDDO
                      IF (stagnant(ji,jst) .GT. zero) THEN
                            wt_soil2(ji,jst)=wt_soil2(ji,jst)+stagnant(ji,jst)
                      ENDIF
                          wt_soil2(ji,jst)= deux*mille - wt_soil2(ji,jst)
!                         wt_soil2(ji,jst)= -(wt_soil2(ji,jst)-deux*mille)
                 ENDDO
                    !!! TEST CHLOE POUR FLUX METHANE ZERO :
                    !wt_soil2(ji,4)=100
            ENDDO



   
!!!!!!!!!!! Chlo� !!!!!!!!!!
            !wt_soil(:,:)=-1.0_r_std
            !wtold(:)=wt_soil(:,4) 
            DO ji=1,kjpindex
                DO jst=1,nstm
                   wt_soil(ji,jst)=-1.0_r_std
                    DO jsl=1,nslm
                        IF (mc(ji,jsl,jst) .GE. 0.39 .AND. ABS(wt_soil(ji,jst)+1.0_r_std) .LT. EPSILON(1.0_r_std) ) THEN
                            !wt_soil = depth du haut de la couche qui est satur� en eau
                            wt_soil(ji,jst) =  zz(jsl,jst) - dz(jsl,jst)/2 
                        ENDIF
                
                    ENDDO
                       IF (wt_soil(ji,jst) .EQ. -1.0_r_std .AND. mc(ji,11,jst) .LT. 0.39 .AND. mc(ji,1,jst) .LT. 0.39) THEN
                            wt_soil(ji,jst)=2100.
                       ENDIF
                       IF (jst .EQ. nstm .AND. wt_soil(ji,jst) .EQ. zero .AND. stagnant(ji,jst) .GT. zero) THEN
                            wt_soil(ji,jst)= -un*stagnant(ji,jst)
                       ENDIF
                ENDDO
            ENDDO


!Chloe v2 17/10/14
!DO ji=1,kjpindex
!    DO jst=1,nstm
!        DO jsl=1,nslm
!             IF (mcs(njsc(ji)) .EQ. 0.41) THEN
!                IF(mc(ji,jsl,jst) .GE. 0.3895 .AND. wt_soil(ji,jst) .GT. -1.0_r_std) THEN  
!                    wt_soil(ji,jst) =  zz(jsl,jst) - dz(jsl,jst)/2
!                ENDIF
!             ENDIF
!             IF  (mcs(njsc(ji)) .EQ. 0.43) THEN
!                IF(mc(ji,jsl,jst) .GE. 0.4085 .AND. wt_soil(ji,jst) .GT. -1.0_r_std) THEN
!                    wt_soil(ji,jst) =  zz(jsl,jst) - dz(jsl,jst)/2
!                ENDIF
!             ENDIF
!            IF (wt_soil(ji,jst) .EQ. -1.0_r_std) THEN
!                wt_soil(ji,jst)=2100
!            ENDIF
!        ENDDO
!    ENDDO
!ENDDO
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 RETURN



  END SUBROUTINE hydrol_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_infilt
!!
!>\BRIEF        Infiltration
!!
!! DESCRIPTION  :
!! - 1 First layer
!! - 2 Infiltration layer by layer 
!! - 2.1 Initialisation
!! - 2.2 Infiltrability of each layer if under a saturated one
!! - 2.3 We compute the mean rate at which water actually infiltrate:
!! - 2.4 From which we deduce the time it takes to fill up the layer or to end the time step...
!! - 2.5 The water enters in the layer
!! - 3 Verification
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_infilt

  SUBROUTINE hydrol_soil_infilt(kjpindex, ins, dtradia, njsc, flux_infilt)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! GLOBAL (in or inout)
    INTEGER(i_std), INTENT(in)                        :: kjpindex        !! Domain size
    REAL(r_std), INTENT (in)                          :: dtradia         !! Time step in seconds
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc            !! indexing of PFT to soiltile
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)    :: flux_infilt     !! Water to infiltrate


    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jsl, ins        !! Indices
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf_pot         !! infiltrable water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf             !! infiltrated water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: dt_tmp              !! time remaining before the end of the time step
    REAL(r_std), DIMENSION (kjpindex)             :: dt_inf              !! the time it takes to complete the infiltration in the 
                                                                         !! layer 
    REAL(r_std)                                   :: k_m                 !! the mean conductivity used for the saturated front
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tmp          !! infiltration rate for the considered layer 
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tot          !! total infiltration 
    REAL(r_std), DIMENSION (kjpindex)             :: flux_tmp            !! rate at which precip hits the ground

    ! If data (or coupling with GCM) was available, a parameterization for subgrid rainfall could be performed


    DO ji = 1, kjpindex
    !-
    !_ 1 First layer
       !-
       ! First we fill up the first layer (about 1mm) without any resistance and quasi-immediately
       wat_inf_pot(ji) = MAX((mcs(njsc(ji))-mc(ji,1,ins)) * dz(2,ins) / deux, zero)
       wat_inf(ji) = MIN(wat_inf_pot(ji), flux_infilt(ji))
       mc(ji,1,ins) = mc(ji,1,ins) + wat_inf(ji) * deux / dz(2,ins)
       !
    ENDDO
    !-
    !! 2 Infiltration layer by layer 
    !! 2.1 Initialisation
    ! Initialize a countdown for infiltration during the time-step and the value of potential runoff
    dt_tmp(:) = dtradia / one_day
    infilt_tot(:) = wat_inf(:)
    ! Compute the rate at which water will try to infiltrate each layer
    flux_tmp(:) = (flux_infilt(:)-wat_inf(:)) / dt_tmp(:)
    !
    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
    !! 2.2 Infiltrability of each layer if under a saturated one
          ! This is computed by an simple arithmetic average because 
          ! the time step (30min) is not appropriate for a geometric average (advised by Haverkamp and Vauclin)
          k_m = (k(ji,jsl) + ks(njsc(ji))*kfact(jsl-1,njsc(ji))*kfact_root(ji,jsl,ins)) / deux

!Isa gel Valdai
!        if (ok_freeze_cwrr.and.temp_hydro(ji, jsl) .lt. tzero)  k_m = k(ji,jsl)
! modif CR pour eviter temp_hydro non init
        if (ok_freeze_cwrr) then
            if (temp_hydro(ji, jsl) .lt. tzero) then
                k_m = k(ji,jsl)
            endif
        endif


    !! 2.3 We compute the mean rate at which water actually infiltrate:
          !- Subgrid: Exponential distribution of  k around i_m, but average p directly used 
          !!??Aurelien: A big black cloud for me
          infilt_tmp(ji) = k_m * (un - EXP(- flux_tmp(ji) / k_m)) 
    !! 2.4 From which we deduce the time it takes to fill up the layer or to end the time step...
          wat_inf_pot(ji) =  MAX((mcs(njsc(ji))-mc(ji,jsl,ins)) * (dz(jsl,ins) + dz(jsl+1,ins)) / deux, zero)

          IF ( infilt_tmp(ji) > min_sechiba) THEN
             dt_inf(ji) =  MIN(wat_inf_pot(ji)/infilt_tmp(ji), dt_tmp(ji))
             ! The water infiltration TIME has to limited by what is still available for infiltration.
             IF ( dt_inf(ji) * infilt_tmp(ji) > flux_infilt(ji)-infilt_tot(ji) ) THEN
                dt_inf(ji) = MAX(flux_infilt(ji)-infilt_tot(ji),zero)/infilt_tmp(ji)
             ENDIF
          ELSE
             dt_inf(ji) = dt_tmp(ji)
          ENDIF



    !! 2.5 The water enters in the layer
          wat_inf(ji) = dt_inf(ji) * infilt_tmp(ji)
          ! bviously the moisture content
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + &
               & wat_inf(ji) * deux / (dz(jsl,ins) + dz(jsl+1,ins))
          ! the time remaining before the next time step
          dt_tmp(ji) = dt_tmp(ji) - dt_inf(ji)
          ! and finally the infilt_tot (which is just used to check if there is a problem, below) 
          infilt_tot(ji) = infilt_tot(ji) + infilt_tmp(ji) * dt_inf(ji)
   ! IF (infilt_tot(ji) .LT. -min_sechiba) THEN
          !write(*,*) 'infilt_tot', infilt_tot(ji)
   !       write(*,*) 'k,ks(njsc(ji)) ji, , jsl, jst, mc', k(ji,1:2),ks(njsc(ji)), ji, jsl,ins, mc(ji,1,ins)
   !       write(*,*) 'Chloe infilt_tmp, dt_inf', infilt_tmp(ji), dt_inf(ji)
   !       write(*,*) 'Chloe k_m', k_m
   !       write(*,*) 'Chloe flux_tmp(ji)', flux_tmp(ji)
   !       write(*,*) 'flux_infilt(:), wat_inf(:)', flux_infilt(ji), wat_inf(ji)
   !     waterbal_error=.TRUE.
   !       CALL ipslerr(3, 'hydrol_soil_infilt', 'We will STOP after hydrol_soil_infilt.','','')
   ! ENDIF  

       ENDDO
    ENDDO
   


    !! 3 Verification
    DO ji = 1, kjpindex
        IF (infilt_tot(ji) .LT. -min_sechiba .OR. infilt_tot(ji) .GT. flux_infilt(ji) + min_sechiba) THEN
   WRITE (numout,*) 'Error in the calculation of infilt tot', infilt_tot(ji)
          WRITE (numout,*) 'k, ji, jst, mc', k(ji,1:2), ji, ins, mc(ji,1,ins)
    !      write(*,*) 'flux_infilt(:), wat_inf(:)', flux_infilt(ji), wat_inf(ji)
          waterbal_error=.TRUE.
          CALL ipslerr(3, 'hydrol_soil_infilt', 'We will STOP after hydrol_soil_infilt.','','')
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE hydrol_soil_infilt


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_smooth
!!
!>\BRIEF        Smooth soil moisture values: avoid over-saturation or under-residual values.
!!
!! DESCRIPTION  :
!! - 1 Avoid over-saturation values
!! - 1.1 top to bottom
!! - 1.2 bottom to top
!! - 2 Avoid below residual values
!! - 2.1 top to bottom
!! - 2.2 bottom to top
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth

  SUBROUTINE hydrol_soil_smooth(kjpindex, ins, njsc, is_under_mcr, is_over_mcs)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables



    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! number of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc            !! indexing of PFT to soiltile


    !! 0.2 Output variables

    LOGICAL, DIMENSION(kjpindex), INTENT(out)          :: is_under_mcr    !! Allows under residual soil moisture due to evap 
    LOGICAL, DIMENSION(kjpindex), INTENT(out)          :: is_over_mcs     !! Allows over saturated soil moisture due to returnflow 

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                       :: ji,jsl
    REAL(r_std)                          :: excess
    REAL(r_std), DIMENSION(kjpindex)     :: excessji

       
    !-
    !! 1 Avoid over-saturation values
    !-
    
    ! in case of over-saturation we put the water where it is possible

    !! 1.1 top to bottom
    DO jsl = 1, nslm-2
       DO ji=1, kjpindex
          excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) + excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl+1,ins)+dz(jsl+2,ins))
      
       ENDDO

    ENDDO


    jsl = nslm-1
    DO ji=1, kjpindex

       excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) + excess * &
            &  (dz(jsl,ins)+dz(jsl+1,ins))/dz(jsl+1,ins)
    ENDDO
   

    jsl = nslm
    DO ji=1, kjpindex

       excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) + excess * &
            &  dz(jsl,ins)/(dz(jsl-1,ins)+dz(jsl,ins))
    ENDDO


    !! 1.2 bottom to top
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) + excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl-1,ins)+dz(jsl,ins))

       ENDDO
    ENDDO

    DO ji=1, kjpindex
       excessji(ji) = mask_soiltile(ji,ins) * MAX(mc(ji,1,ins)-mcs(njsc(ji)),zero)
    ENDDO

    DO ji=1, kjpindex
       mc(ji,1,ins) = mc(ji,1,ins) - excessji(ji)

       is_over_mcs(ji) = (excessji(ji) .GT. min_sechiba)
    ENDDO
    

    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excessji(ji) * dz(2,ins) / (deux * dpu_max*mille)
       ENDDO
    ENDDO
    

    !-
    !! 2 Avoid below residual values
    !-
       
    ! Smooth the profile to avoid negative values of punctual soil moisture

   
    ! 2.1 top to bottom
    DO jsl = 1,nslm-2

       DO ji=1, kjpindex
          excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl+1,ins)+dz(jsl+2,ins))
 
       ENDDO

    ENDDO

    jsl = nslm-1
    DO ji=1, kjpindex
       excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
            &  (dz(jsl,ins)+dz(jsl+1,ins))/dz(jsl+1,ins)

    ENDDO


    jsl = nslm
    DO ji=1, kjpindex
       excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
            &  dz(jsl,ins)/(dz(jsl-1,ins)+dz(jsl,ins))
    ENDDO


    !! 2.2 bottom to top
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl-1,ins)+dz(jsl,ins))

        ENDDO
    ENDDO

    DO ji=1, kjpindex
         excessji(ji) = mask_soiltile(ji,ins) * MAX(mcr(njsc(ji))-mc(ji,1,ins),zero)

    ENDDO
   
    DO ji=1, kjpindex
       mc(ji,1,ins) = mc(ji,1,ins) + excessji(ji)
       is_under_mcr(ji) = (excessji(ji) .GT. min_sechiba)

    ENDDO

        DO jsl = 1, nslm
            DO ji=1, kjpindex
               mc(ji,jsl,ins) = mc(ji,jsl,ins) - excessji(ji) * dz(2,ins) / (deux * dpu_max*mille)
            ENDDO

        ENDDO
   
    ! We just get sure that mc remains at 0 where soiltile=0
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mask_soiltile(ji,ins) * mc(ji,jsl,ins)
        ENDDO

    ENDDO
 



    RETURN
  END SUBROUTINE hydrol_soil_smooth


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_flux
!!
!>\BRIEF        This subroutine computes the hydrological fluxes between the different soil layers.
!!
!! DESCRIPTION  :
!! - 1 Initialize qflux from the bottom, with dr_ns
!! - 2 between layer nslm and nslm-1   
!! - 3 we go to top and deduct qflux(1:nslm-2)   
!! - 4 Water balance verification
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_flux

  SUBROUTINE hydrol_soil_flux(kjpindex,ins,mcint,returnflow_soil)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! index of soil type
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: mcint           !! mc values at the beginning of the time step
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)      :: returnflow_soil !! returnflow

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jsl,ji
    REAL(r_std), DIMENSION(kjpindex)                   :: temp
    
    !- Compute the flux at every level from bottom to top (using mc and sink values)
    DO ji = 1, kjpindex
    !! 1 Initialize qflux from the bottom, with dr_ns
       jsl = nslm
       qflux(ji,jsl,ins) = dr_ns(ji,ins) - returnflow_soil(ji)
    !!_ between layer nslm and nslm-1   
       jsl = nslm-1
       qflux(ji,jsl,ins) = qflux(ji,jsl+1,ins) & 
            &  + (mc(ji,jsl,ins)-mcint(ji,jsl) &
            &  + trois*mc(ji,jsl+1,ins) - trois*mcint(ji,jsl+1)) &
            &  * (dz(jsl+1,ins)/huit) &
            &  + rootsink(ji,jsl+1,ins) 
    ENDDO
    !! 3 we go to top and deduct qflux(1:nslm-2)   
    DO jsl = nslm-2,1,-1
       DO ji = 1, kjpindex
          qflux(ji,jsl,ins) = qflux(ji,jsl+1,ins) & 
               &  + (mc(ji,jsl,ins)-mcint(ji,jsl) &
               &  + trois*mc(ji,jsl+1,ins) - trois*mcint(ji,jsl+1)) &
               &  * (dz(jsl+1,ins)/huit) &
               &  + rootsink(ji,jsl+1,ins) &
               &  + (dz(jsl+2,ins)/huit) &
               &  * (trois*mc(ji,jsl+1,ins) - trois*mcint(ji,jsl+1) &
               &  + mc(ji,jsl+2,ins)-mcint(ji,jsl+2)) 
       END DO
    ENDDO
    
    !! 4 Water balance verification  
    DO ji = 1, kjpindex
       temp(ji) =  qflux(ji,1,ins) + (dz(2,ins)/huit) &
            &  * (trois* (mc(ji,1,ins)-mcint(ji,1)) + (mc(ji,2,ins)-mcint(ji,2))) &
            &  + rootsink(ji,1,ins)
    ENDDO

    DO ji = 1, kjpindex
       IF (ABS(qflux00(ji,ins)-temp(ji)).GT. deux*min_sechiba) THEN
          WRITE(numout,*) 'Problem in the water balance, qflux computation', qflux00(ji,ins),temp(ji)
          WRITE (numout,*) 'returnflow_soil', returnflow_soil(ji)
          WRITE(numout,*) 'ji', ji, 'jsl',jsl,'ins',ins
          WRITE(numout,*) 'mcint', mcint(ji,:)
          WRITE(numout,*) 'mc', mc(ji,:,ins)
          WRITE (numout,*) 'rootsink', rootsink(ji,1,ins)
          waterbal_error=.TRUE.
          CALL ipslerr(3, 'hydrol_soil_flux', 'We will STOP after hydrol_soil_flux.',&
               & 'Problem in the water balance, qflux computation','')
       ENDIF
    ENDDO



    RETURN
  END SUBROUTINE hydrol_soil_flux


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_tridiag
!!
!>\BRIEF        This subroutine solves a set of linear equations which has a tridiagonal coefficient matrix. 
!!
!! DESCRIPTION  : None
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_tridiag 

  SUBROUTINE hydrol_soil_tridiag(kjpindex,ins)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! number of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     ::  ji,jsl
    REAL(r_std), DIMENSION(kjpindex)                   :: bet

!Isa
    REAL(r_std)                                         :: mc_beg

    DO ji = 1, kjpindex

       IF (resolv(ji)) THEN
          bet(ji) = tmat(ji,1,2)
          mcl(ji,1,ins) = rhs(ji,1)/bet(ji)
       ENDIF
    ENDDO

    DO jsl = 2,nslm
       DO ji = 1, kjpindex
          
          IF (resolv(ji)) THEN

             gam(ji,jsl) = tmat(ji,jsl-1,3)/bet(ji)
             bet(ji) = tmat(ji,jsl,2) - tmat(ji,jsl,1)*gam(ji,jsl)
             mcl(ji,jsl,ins) = (rhs(ji,jsl)-tmat(ji,jsl,1)*mcl(ji,jsl-1,ins))/bet(ji)
          ENDIF

       ENDDO
    ENDDO

    DO ji = 1, kjpindex
       IF (resolv(ji)) THEN
          DO jsl = nslm-1,1,-1
             mcl(ji,jsl,ins) = mcl(ji,jsl,ins) - gam(ji,jsl+1)*mcl(ji,jsl+1,ins)
          ENDDO
       ENDIF
    ENDDO



    RETURN

  END SUBROUTINE hydrol_soil_tridiag
 
!Isa
 SUBROUTINE hydrol_soil_coef1(kjpindex,ins,njsc) !pr infiltration comme en sol non gel�
    !
    IMPLICIT NONE
    !
    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              ! index of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc            !! Indeces of the soiltype

    ! local 
    INTEGER(i_std) :: jsl,ji,i


    DO jsl=1,nslm
       DO ji=1,kjpindex 
          i= MAX(MIN(INT((imax-imin)*(mc(ji,jsl,ins)-mcr(njsc(ji)))&
               &         / (mcs(njsc(ji))-mcr(njsc(ji))))+imin , imax-1), imin)
          
          a(ji,jsl) = a_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          b(ji,jsl) = b_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          d(ji,jsl) = d_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          k(ji,jsl) = MAX(k_lin(imin+1,jsl,njsc(ji)), &
               & a_lin(i,jsl,njsc(ji)) * mc(ji,jsl,ins) + b_lin(i,jsl,njsc(ji)))
                 ENDDO ! loop on grid
    ENDDO



    RETURN
  END SUBROUTINE hydrol_soil_coef1


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_coef
!!
!>\BRIEF        Computes coef for the linearised hydraulic conductivity 
!! k_lin=a_lin mc_lin+b_lin and the linearised diffusivity d_lin. 
!!
!! DESCRIPTION  :
!! First, we identify the interval i in which the current value of mc is located.
!! Then, we give the values of the linearized parameters to compute 
!! conductivity and diffusivity as K=a*mc+b and d.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_coef
 
  SUBROUTINE hydrol_soil_coef(kjpindex,ins,njsc)

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              !! Index of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc             !! indexing of PFT to soiltile


    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: jsl,ji,i
    INTEGER :: mc_entier_max, mc_entier
    REAL(r_std)                                       :: mc_ratio
!Isa
    REAL(r_std)              :: mc_used,x,m, mc_incremente ! contenu en eau liquide r�el 

    !-first, we identify the interval i in which the current value of mc is located
    !-then, we give the values of the linearized parameters to compute 
    ! conductivity and diffusivity as K=a*mc+b and d
    !
    !Isa : Van Genuchten parameter for thermodynamical calculation
    !
!Chloe+ nvan(ins) chang� en nvan(njsc)   
    !
    DO jsl=1,nslm
       DO ji=1,kjpindex 


       m = 1.-1./nvan(njsc(ji))

if (.NOT.ok_freeze_cwrr) then
	  i= MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc(ji,jsl,ins)-mcr(njsc(ji)))/(mcs(njsc(ji))-mcr(njsc(ji))))))
          a(ji,jsl) = a_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          b(ji,jsl) = b_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          d(ji,jsl) = d_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          k(ji,jsl) = MAX(k_lin(imin+1,jsl,njsc(ji)), &
               & a_lin(i,jsl,njsc(ji)) *mc(ji,jsl,ins) + b_lin(i,jsl,njsc(ji)))


else !if (.NOT.ok_freeze_cwrr) then
        if (ok_gel_thermosoil.OR.(mc(ji,jsl, ins).lt.(mcr(njsc(ji))+min_sechiba))) then
                if (temp_hydro(ji, jsl).ge.273._r_std) then
                        x=1._r_std
                else if (273._r_std.gt.temp_hydro(ji, jsl).AND.temp_hydro(ji, jsl).ge.271._r_std) then 
                        x=(temp_hydro(ji, jsl)-271._r_std)/2._r_std
                else 
                        x=0._r_std
                endif
        else if (ok_gel_thd) then
                if (temp_hydro(ji, jsl).ge.273._r_std) then
                        x=1._r_std
                else if (273._r_std.gt.temp_hydro(ji, jsl).AND.temp_hydro(ji, jsl).ge.271._r_std) then 
                        x=MIN(((mcs(njsc(ji))-mcr(njsc(ji))) &
                  & *((1000.*avan(njsc(ji))*(tzero-temp_hydro(ji, jsl)) &
                  & *lhf/tzero/10.)**nvan(njsc(ji))+1.)**(-m))/(mc(ji,jsl, ins) &
                  & -mcr(njsc(ji))),1._r_std)                
          else 
                        x=0._r_std
                endif

        endif
	profil_froz_hydro_ns(ji, jsl,ins)=1._r_std-x



!x=liquid saturation degree/residual=(mcl-mcr)/(mcs-mcr)
!1-x=frozen saturation degree/residua=(mcf-mcr)/(mcs-mcr)
!(1-x)==profil_froz_hydro

	mc_used = mcr(njsc(ji))+x*(mc(ji,jsl,ins)-mcr(njsc(ji))) 

 
    ! utilis� pour le calcul des ppt�s hydro uniquement => pas de pb si
        ! sous-saturation ni sursaturation (x limit� par 1).
!
! calcul de k based on mc_liq
!
!        if (temp_hydro(ji, jsl).ge.271._r_std) then	
	  i= MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc_used-mcr(njsc(ji)))/(mcs(njsc(ji))-mcr(njsc(ji))))))

	  a(ji,jsl) = a_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          b(ji,jsl) = b_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          d(ji,jsl) = d_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
          k(ji,jsl) = MAX(k_lin(imin+1,jsl,njsc(ji)), &
               & a_lin(i,jsl,njsc(ji)) * mc_used + b_lin(i,jsl,njsc(ji)))

  
! Chloe 030314 :
       !    mc_entier_max=FLOOR((mcs(njsc(ji))-mcr(njsc(ji)))/0.05)
       ! 
       !   DO mc_entier= 1,mc_entier_max
       !    mc_incremente=mcr(njsc(ji))+(mc_entier-1.)*0.05
       !    mc_used = mcr(njsc(ji))+x*(mc_incremente-mcr(njsc(ji))) 
       !    i=  MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc_used-mcr(njsc(ji)))/(mcs(njsc(ji))-mcr(njsc(ji))))))
       !    k(ji,jsl) = MAX(k_lin(imin+1,jsl,njsc(ji)), &
       !        & a_lin(i,jsl,njsc(ji)) * mc_used + b_lin(i,jsl,njsc(ji)))
       !    IF (jsl .EQ. nslm) then 
       !        write(16,*) 'k', k(1,:)
       !        write(17,*) 'mc_used k', mc_used
       !        !IF (mc_entier .EQ. mc_entier_max) STOP
       !    ENDIF
       !    ENDDO


!	else
!	  d(ji,jsl)=zero
!	  k(ji,jsl)=zero
!          a(ji,jsl)=zero
!          b(ji,jsl)=zero
!	endif

endif !if (.NOT.ok_freeze_cwrr) then
         ! GK040314
         !write(15,*) '--',i,jsl,njsc(ji)
         !write(15,*) '  ', a_lin(i,jsl,njsc(ji)),b_lin(i,jsl,njsc(ji))
         !write(15,*) '  ',k_lin(imin+1,jsl,njsc(ji)) 
         !write(15,*) '  ',k(ji,jsl),mc(ji,jsl,ins),mcr(njsc(ji)),mcs(njsc(ji)), mc_used

       ENDDO ! loop on grid
         ! GK040314
         !write(15,*) 
    ENDDO ! loop on hydro layers
    


    RETURN
  END SUBROUTINE hydrol_soil_coef


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_setup
!!
!>\BRIEF        This subroutine computes the matrix coef.  
!!
!! DESCRIPTION  : None 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : matrix coef
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_soil_setup(kjpindex,ins,dtradia)


    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL(r_std), INTENT (in)                           :: dtradia          !! Time step in seconds
    ! parameters
    INTEGER(i_std), INTENT(in)                        :: kjpindex          !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins               !! index of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: jsl,ji
    REAL(r_std)                        :: temp3, temp4

    !-we compute tridiag matrix coefficients (LEFT and RIGHT) 
    ! of the system to solve [LEFT]*mc_{t+1}=[RIGHT]*mc{t}+[add terms]: 
    ! e(nslm),f(nslm),g1(nslm) for the [left] vector
    ! and ep(nslm),fp(nslm),gp(nslm) for the [right] vector

    ! w_time=1 (in constantes_soil) indicates implicit computation for diffusion 
    temp3 = w_time*(dtradia/one_day)/deux
    temp4 = (un-w_time)*(dtradia/one_day)/deux

    ! Passage to arithmetic means for layer averages also in this subroutine : Aurelien 11/05/10

    !- coefficient for first layer
    DO ji = 1, kjpindex
       e(ji,1) = zero
       f(ji,1) = trois * dz(2,ins)/huit  + temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))+a(ji,1))
       g1(ji,1) = dz(2,ins)/(huit)       - temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))-a(ji,2))
       ep(ji,1) = zero
       fp(ji,1) = trois * dz(2,ins)/huit - temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))+a(ji,1))
       gp(ji,1) = dz(2,ins)/(huit)       + temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))-a(ji,2))
    ENDDO

    !- coefficient for medium layers

    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
          e(ji,jsl) = dz(jsl,ins)/(huit)                        - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins))+a(ji,jsl-1))

          f(ji,jsl) = trois * (dz(jsl,ins)+dz(jsl+1,ins))/huit  + temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins)) )

          g1(ji,jsl) = dz(jsl+1,ins)/(huit)                     - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins))-a(ji,jsl+1))

          ep(ji,jsl) = dz(jsl,ins)/(huit)                       + temp4 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins))+a(ji,jsl-1))

          fp(ji,jsl) = trois * (dz(jsl,ins)+dz(jsl+1,ins))/huit - temp4 &
               & * ( (d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins)) )

          gp(ji,jsl) = dz(jsl+1,ins)/(huit)                     + temp4 &
               & *((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins))-a(ji,jsl+1))
       ENDDO
    ENDDO

    !- coefficient for last layer
    DO ji = 1, kjpindex
       e(ji,nslm) = dz(nslm,ins)/(huit)        - temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm,ins))+a(ji,nslm-1))
       f(ji,nslm) = trois * dz(nslm,ins)/huit  + temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) / (dz(nslm,ins)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,1)))
       g1(ji,nslm) = zero
       ep(ji,nslm) = dz(nslm,ins)/(huit)       + temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm,ins))+a(ji,nslm-1))
       fp(ji,nslm) = trois * dz(nslm,ins)/huit - temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm,ins)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,1)))
       gp(ji,nslm) = zero
    ENDDO


    RETURN
  END SUBROUTINE hydrol_soil_setup


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_split_soil
!!
!>\BRIEF        Splits 2d variables into 3d variables, per soil type. 
!!
!! DESCRIPTION  :
!! - 1 split 2d variables into 3d variables, per soil type
!! - 1.1 precipitation
!! - 1.2 evaporation
!! - 1.2.1 vevapnu_old
!! - 1.2.2 ae_ns new
!! - 1.3 transpiration
!! - 1.4 root sink
!! - 2 Verification
!! - 2.1 Check of mc
!! - 2.2 Check if the deconvolution is correct and conserves the fluxes
!! - 2.2.1 the precisol and evapnu
!! - 2.2.2 the transpiration and root sink
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_split_soil

  SUBROUTINE hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, transpir, humrel,evap_bare_lim)
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)       :: veget_max        !! max Vegetation map 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile         !! fraction of PFT on each soil-hydrology tile
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: vevapnu          !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: transpir         !! Transpiration
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: humrel           !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evap_bare_lim    !!   

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jv, jsl, jst
    REAL(r_std), Dimension (kjpindex)             :: vevapnu_old
    REAL(r_std), Dimension (kjpindex)             :: tmp_check1
    REAL(r_std), Dimension (kjpindex)             :: tmp_check2
    REAL(r_std), DIMENSION (kjpindex,nstm)        :: tmp_check3
    REAL(r_std)                                   :: test
    !
    !
    !! 1 split 2d variables into 3d variables, per soil type
    !
    !
    !! 1.1 precipitation
    precisol_ns(:,:)=zero
    DO jv=1,nvm
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF(veget_max(ji,jv).GT.min_sechiba) THEN
                precisol_ns(ji,jst)=precisol_ns(ji,jst)+precisol(ji,jv)* &
                     & corr_veg_soil(ji,jv,jst) /vegtot(ji) / veget_max(ji,jv)
             ENDIF
          END DO
       END DO
    END DO
    !
    !
    !! 1.2 evaporation
    !! 1.2.1 vevapnu_old
    vevapnu_old(:)=zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF ( vegtot(ji) .GT. min_sechiba) THEN
             vevapnu_old(ji)=vevapnu_old(ji)+ &
                  & ae_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          ENDIF
       END DO
    END DO
    !
    !
    !
    !! 1.2.2 ae_ns new
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF (vevapnu_old(ji).GT.min_sechiba) THEN   
             IF(evap_bare_lim(ji).GT.min_sechiba) THEN       
                ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
             ELSE
                IF(vevapnu_old(ji).GT.min_sechiba) THEN  
                   ae_ns(ji,jst)=ae_ns(ji,jst) * vevapnu(ji)/vevapnu_old(ji)
                ELSE
                   ae_ns(ji,jst)=zero
                ENDIF
             ENDIF
          ELSEIF(frac_bare_ns(ji,jst).GT.min_sechiba) THEN
             IF(evap_bare_lim(ji).GT.min_sechiba) THEN  
                ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
             ELSE
                IF(tot_bare_soil(ji).GT.min_sechiba) THEN  
                   ae_ns(ji,jst) = vevapnu(ji) * frac_bare_ns(ji,jst)/tot_bare_soil(ji)
                ELSE
                   ae_ns(ji,jst) = zero
                ENDIF
             ENDIF
          ENDIF
          precisol_ns(ji,jst)=precisol_ns(ji,jst)+MAX(-ae_ns(ji,jst),zero)
       END DO
    END DO
    !
    !
    !! 1.3 transpiration
    tr_ns(:,:)=zero
    DO jv=1,nvm
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF (humrel(ji,jv).GT.min_sechiba) THEN 
                tr_ns(ji,jst)=tr_ns(ji,jst)+ cvs_over_veg(ji,jv,jst)*humrelv(ji,jv,jst)* & 
                     & transpir(ji,jv)/humrel(ji,jv)
             ENDIF
          END DO
       END DO
    END DO

    !
    !
    !! 1.4 root sink
    rootsink(:,:,:)=zero
    DO jv=1,nvm
       DO jsl=1,nslm
          DO jst=1,nstm
             DO ji=1,kjpindex
                IF (humrel(ji,jv).GT.min_sechiba) THEN 
                   rootsink(ji,jsl,jst) = rootsink(ji,jsl,jst) &
                        & + cvs_over_veg(ji,jv,jst)* (transpir(ji,jv)*us(ji,jv,jst,jsl))/ &
                        & humrel(ji,jv)
                END IF
             END DO
          END DO
       END DO
    END DO

    !! 2 Verification
    !! 2.1 Check of mc
    IF(check_cwrr) THEN
       DO jsl=1,nslm
          DO jst=1,nstm
             DO ji=1,kjpindex
                IF(mc(ji,jsl,jst).LT.-0.05) THEN
                   WRITE(numout,*) 'CWRR split-----------------------------------------------'
                   WRITE(numout,*) 'ji,jst,jsl',ji,jst,jsl
                   WRITE(numout,*) 'mc',mc(ji,jsl,jst)
                   WRITE(numout,*) 'rootsink,us',rootsink(ji,:,jst),us(ji,:,jst,jsl)
                   WRITE(numout,*) 'corr_veg_soil',corr_veg_soil(ji,:,jst)
                   WRITE(numout,*) 'transpir',transpir(ji,:)
                   WRITE(numout,*) 'veget_max',veget_max(ji,:)
                   WRITE(numout,*) 'cvs_over_veg',cvs_over_veg(ji,:,jst)
                   WRITE(numout,*) 'humrel',humrel(ji,:)
                   WRITE(numout,*) 'humrelv (pour ce jst)',humrelv(ji,:,jst)
                   WRITE(numout,*) 'ae_ns',ae_ns(ji,jst)
                   WRITE(numout,*) 'tr_ns',tr_ns(ji,jst)
                   WRITE(numout,*) 'vevapnuold',vevapnu_old(ji)
                ENDIF
             END DO
          END DO
       END DO
    ENDIF


    !! 2.2 Check if the deconvolution is correct and conserves the fluxes

    IF (check_cwrr) THEN


       tmp_check1(:)=zero
       tmp_check2(:)=zero  

    !! 2.2.1 the precisol and evapnu

       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + &
                  & (precisol_ns(ji,jst)-MAX(-ae_ns(ji,jst),zero))* &
                  & soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + precisol(ji,jv)
          END DO
       END DO


       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji)- tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'PRECISOL SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)

             DO jv=1,nvm
                WRITE(numout,'(a,i2.2,"|",F13.4,"|",F13.4,"|",3(F9.6))') 'jv,veget_max, precisol, corr_veg_soil ',&
                     jv,veget_max(ji,jv),precisol(ji,jv),corr_veg_soil(ji,jv,:)
             END DO

             DO jst=1,nstm
                WRITE(numout,*) 'jst,precisol_ns',jst,precisol_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             waterbal_error=.TRUE.
             CALL ipslerr(2, 'hydrol_split_soil', 'We will STOP after hydrol_split_soil.',&
                  & 'check_CWRR','PRECISOL SPLIT FALSE')
          ENDIF

       END DO


       tmp_check1(:)=zero

       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + ae_ns(ji,jst)* &
                  & soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji)- vevapnu(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'VEVAPNU SPLIT FALSE:ji, Sum(ae_ns), vevapnu =',ji,tmp_check1(ji),vevapnu(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- vevapnu(ji))
             WRITE(numout,*) 'ae_ns',ae_ns(ji,:)
             WRITE(numout,*) 'vegtot',vegtot(ji)
             WRITE(numout,*) 'evap_bare_lim, evap_bare_lim_ns',evap_bare_lim(ji), evap_bare_lim_ns(ji,:)
             WRITE(numout,*) 'tot_bare_soil,frac_bare_ns',tot_bare_soil(ji),frac_bare_ns(ji,:)
             WRITE(numout,*) 'vevapnu_old',vevapnu_old(ji)
             DO jst=1,nstm
                WRITE(numout,*) 'jst,ae_ns',jst,ae_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
                WRITE(numout,*) 'veget_exist/soiltile', veget_max(ji,:)/vegtot(ji)/soiltile(ji,jst)
                WRITE(numout,*) "corr_veg_soil",corr_veg_soil(ji,:,jst)
             END DO
             waterbal_error=.TRUE.
             CALL ipslerr(2, 'hydrol_split_soil', 'We will STOP after hydrol_split_soil.',&
                  & 'check_CWRR','VEVAPNU SPLIT FALSE')
          ENDIF
       ENDDO

    !! 2.2.2 the transpiration and root sink

       tmp_check1(:)=zero
       tmp_check2(:)=zero  


       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + tr_ns(ji,jst)* &
                  & soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + transpir(ji,jv)
          END DO
       END DO

       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji)- tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'TRANSPIR SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)

             DO jv=1,nvm
                WRITE(numout,*) 'jv,veget_max, transpir',jv,veget_max(ji,jv),transpir(ji,jv)
                DO jst=1,nstm
                   WRITE(numout,*) 'corr_veg_soil:ji,jv,jst',ji,jv,jst,corr_veg_soil(ji,jv,jst)
                END DO
             END DO

             DO jst=1,nstm
                WRITE(numout,*) 'jst,tr_ns',jst,tr_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             waterbal_error=.TRUE.
             CALL ipslerr(2, 'hydrol_split_soil', 'We will STOP after hydrol_split_soil.',&
                  & 'check_CWRR','TRANSPIR SPLIT FALSE')
          ENDIF

       END DO


       tmp_check3(:,:)=zero

       DO jst=1,nstm
          DO jsl=1,nslm
             DO ji=1,kjpindex
                tmp_check3(ji,jst)=tmp_check3(ji,jst) + rootsink(ji,jsl,jst)
             END DO
          END DO
       ENDDO

       DO jst=1,nstm
          DO ji=1,kjpindex
             IF(ABS(tmp_check3(ji,jst)- tr_ns(ji,jst)).GT.allowed_err) THEN
                WRITE(numout,*) 'ROOTSINK SPLIT FALSE:ji,jst=', ji,jst,&
                     & tmp_check3(ji,jst),tr_ns(ji,jst)
                WRITE(numout,*) 'err',ABS(tmp_check3(ji,jst)- tr_ns(ji,jst))
                WRITE(numout,*) 'HUMREL(jv=1:13)',humrel(ji,:)
                WRITE(numout,*) 'TRANSPIR',transpir(ji,:)
                DO jv=1,nvm 
                   WRITE(numout,*) 'jv=',jv,'us=',us(ji,jv,jst,:)
                ENDDO
                waterbal_error=.TRUE.
                CALL ipslerr(2, 'hydrol_split_soil', 'We will STOP after hydrol_split_soil.',&
                  & 'check_CWRR','ROOTSINK SPLIT FALSE')
             ENDIF
          END DO
       END DO

    ENDIF



    RETURN 

  END SUBROUTINE hydrol_split_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_diag_soil
!!
!>\BRIEF        ??
!!
!! DESCRIPTION  :
!! - 1 Apply mask_soiltile
!! - 2 sum 3d variables in 2d variables with fraction of vegetation per soil type
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_diag_soil

  SUBROUTINE hydrol_diag_soil (kjpindex, veget_max, soiltile, njsc, runoff, drainage, &
       & evap_bare_lim, evapot, vevapnu, returnflow, reinfiltration, irrigation, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, drysoil_frac, tot_melt, water2add_peat, ok_routage_peat)
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max       !! Max. vegetation type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc            !! indexing of PFT to soiltile
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile        !! fraction of PFT on each soil-hydrology tile
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot          !! 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow      !! Water returning to the deep reservoir
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration  !! Water returning to the top of the soil
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation      !! Water from irrigation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt        !!
    !Chloe
    LOGICAL, INTENT(in)                                      :: ok_routage_peat     

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac    !! Function of litter wetness
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff          !! complete runoff
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage        !! Drainage
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: evap_bare_lim   !! 
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out)      :: shumdiag        !! relative soil moisture       !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out)      :: shumdiag_perma ! Isa
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: k_litt          !! litter cond.
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: litterhumdiag   !! litter humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: humrel          !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress       !! Veg. moisture stress (only for vegetation growth)
!Chloe
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: water2add_peat !! water we need to add for peat saturation  

    !! 0.3 Modified variables 

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu         !!

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv, jsl, jst, i
    REAL(r_std), DIMENSION (kjpindex)                                    :: mask_vegtot
    REAL(r_std)                                              :: k_tmp, tmc_litter_ratio
    !
    ! Put the prognostics variables of soil to zero if soiltype is zero

!Isa
    INTEGER(i_std)				             :: jd

!    write(*,*) 'hydrol 5422'
    !! 1 Apply mask_soiltile
    DO jst=1,nstm 

       DO ji=1,kjpindex

       !!??Aurelien: who and why coment this line : if et endif
 !         IF(soiltile(ji,jst).EQ.zero) THEN

             ae_ns(ji,jst) = ae_ns(ji,jst) * mask_soiltile(ji,jst)
             dr_ns(ji,jst) = dr_ns(ji,jst) * mask_soiltile(ji,jst)
             ru_ns(ji,jst) = ru_ns(ji,jst) * mask_soiltile(ji,jst)
             tmc(ji,jst) =  tmc(ji,jst) * mask_soiltile(ji,jst)
             !Chloe
             IF (ok_routage_peat .AND. jst .EQ. 4) tmc_obj_peat(ji,jst)=tmc_obj_peat(ji,jst) * mask_soiltile(ji,jst)
               !Chloe
               !write(*,*) 'Chloe ji, tmc_obj_peat, tmc', ji, tmc_obj_peat(ji,4), tmc(ji,4)          
             if (ok_freeze_cwrr) then
                profil_froz_hydro_ns(ji,:,jst)=profil_froz_hydro_ns(ji,:,jst)*mask_soiltile(ji,jst)
             endif !if (ok_freeze_cwrr) then

             DO jv=1,nvm
                humrelv(ji,jv,jst) = humrelv(ji,jv,jst) * mask_soiltile(ji,jst)
                DO jsl=1,nslm
                   us(ji,jv,jst,jsl) = us(ji,jv,jst,jsl)  * mask_soiltile(ji,jst)
                END DO
             END DO

             DO jsl=1,nslm          
                mc(ji,jsl,jst) = mc(ji,jsl,jst)  * mask_soiltile(ji,jst)
             END DO

  !        ENDIF

       END DO
    END DO

    runoff(:) = zero
    drainage(:) = zero
    humtot(:) = zero
    evap_bare_lim(:) = zero
    evap_bare_lim_ns(:,:) = zero
    shumdiag(:,:)= zero
    shumdiag_perma(:,:)=zero
    k_litt(:) = zero
    litterhumdiag(:) = zero
    tmc_litt_mea(:) = zero
    soilmoist(:,:) = zero
    humrel(:,:) = zero
    vegstress(:,:) = zero
    swi(:) = zero
    if (ok_freeze_cwrr) then
       profil_froz_hydro(:,:)=zero
    endif !if (ok_freeze_cwrr) then 
!    write(*,*) 'hydrol 5470'
    !
    !! 2 sum 3d variables in 2d variables with fraction of vegetation per soil type

    DO ji = 1, kjpindex
       mask_vegtot(ji) = 0
       IF(vegtot(ji) .GT. min_sechiba) THEN
          mask_vegtot(ji) = 1
       ENDIF
    END DO
    
    DO ji = 1, kjpindex 
       ! Here we weight ae_ns by the fraction of bare evaporating soil. 
       ! This is given by frac_bare_ns, taking into account bare soil under vegetation
       ae_ns(ji,:) = mask_vegtot(ji) * ae_ns(ji,:) * frac_bare_ns(ji,:)
    END DO


!!!!!Schizophr�nie : 10/10/14 pourquoi tu met ca la chlo� ?
!!! 1410 : parce que la c'est juste avant le calcul du runoff, pardi !
    ! Chloe++ Stockage du r�servoir d'eau stagnante
 !  DO ji = 1, kjpindex 
 !   DO jst = 1, nstm
 !      IF ( ok_stagnant .AND. jst .EQ. nstm  ) THEN
 !           stagnant(ji,jst) = stagnant(ji,jst) + ru_ns(ji,jst) 
 !           ru_ns(ji,jst) = 0.
 !        
 !        IF ( stagnant(ji,jst) .GT. max_stagnant ) THEN
 !           ru_ns(ji,jst) = stagnant(ji,jst) - max_stagnant
 !           stagnant(ji,jst) = max_stagnant
 !        ENDIF

!          IF (ok_freeze_cwrr) THEN
!             IF (temp_hydro(ji,1) .LT. ZeroCelsius) THEN
!                water2infilt(ji,jst) = water2infilt(ji,jst) + ru_ns(ji,jst)
!                ru_ns(:,jst) =  zero
!             ENDIF
!          ENDIF

  !     ELSE
  !        stagnant(:,jst) = zero ! Valeur non physique
  !     ENDIF
  !  ENDDO
  ! ENDDO
    ! Chloe--
      

!##### GARDE FOUUUU #####
!DO ji=1,kjpindex
!IF (stagnant(ji,4) .NE. -1.) THEN
!STOP 'stagnant ne -1'
!ENDIF
!ENDDO

    ! Chloe--


    DO jst = 1, nstm
    DO ji = 1, kjpindex 
          drainage(ji) = mask_vegtot(ji) * (drainage(ji) + vegtot(ji)*soiltile(ji,jst) * dr_ns(ji,jst))
          runoff(ji) = mask_vegtot(ji) *  (runoff(ji) +   vegtot(ji)*soiltile(ji,jst) * ru_ns(ji,jst)) &
               & + (1 - mask_vegtot(ji)) * (tot_melt(ji) + irrigation(ji) + returnflow(ji) + reinfiltration(ji))
          humtot(ji) = mask_vegtot(ji) * (humtot(ji) + soiltile(ji,jst) * tmc(ji,jst))
          if (ok_freeze_cwrr) then 
	    profil_froz_hydro(ji,:)=mask_vegtot(ji) * (profil_froz_hydro(ji,:) + soiltile(ji,jst) * profil_froz_hydro_ns(ji,:, jst))
          endif !if (ok_freeze_cwrr) then
         !Isa
       END DO
    END DO
   
   
    ! Chloe
    ! Quantity of water needed to get saturation in peatland (used in routing module)
    DO jst=1,nstm
        IF (ok_routage_peat .AND. jst .EQ. 4) THEN   
            DO ji= 1,kjpindex 
                water2add_peat(ji)= tmc_obj_peat(ji,4) - tmc(ji,4)
             ! write(*,*) 'Chloe 2 ji tmc_obj, tmc', ji, tmc_obj_peat(ji,4), tmc(ji,4), water2add_peat(ji)
            ENDDO
        ENDIF
    ENDDO


    DO jst=1,nstm
       DO ji=1,kjpindex
          IF ((evapot(ji).GT.min_sechiba) .AND. &
               & (tmc_litter(ji,jst).GT.(tmc_litter_wilt(ji,jst)))) THEN
             evap_bare_lim_ns(ji,jst) = ae_ns(ji,jst) / evapot(ji)
          ELSEIF((evapot(ji).GT.min_sechiba).AND. &
               & (tmc_litter(ji,jst).GT.(tmc_litter_res(ji,jst)))) THEN
             evap_bare_lim_ns(ji,jst) =  (un/deux) * ae_ns(ji,jst) / evapot(ji)
          END IF
       END DO
    END DO
!    write(*,*) 'hydrol 5513'

     DO ji = 1, kjpindex
         evap_bare_lim(ji) =   SUM(evap_bare_lim_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
        IF(evap_bare_lim(ji).GT.un + min_sechiba) THEN
           WRITE(numout,*) 'CWRR DIAG EVAP_BARE_LIM TOO LARGE', ji, &
                & evap_bare_lim(ji),evap_bare_lim_ns(ji,:)
        ENDIF 
     ENDDO
    ! we add the excess of snow sublimation to vevapnu

    DO ji = 1,kjpindex
       vevapnu(ji) = vevapnu (ji) + subsinksoil(ji)*vegtot(ji)
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             IF(veget_max(ji,jv).GT.min_sechiba) THEN
                vegstress(ji,jv)=vegstress(ji,jv)+vegstressv(ji,jv,jst)*soiltile(ji,jst) &
                     & * corr_veg_soil(ji,jv,jst) *vegtot(ji)/veget_max(ji,jv)
                vegstress(ji,jv)= MAX(vegstress(ji,jv),zero)
             ENDIF
          END DO
       END DO
    END DO

    cvs_over_veg(:,:,:) = zero
    DO jv=1,nvm
       DO ji=1,kjpindex
          IF(veget_max(ji,jv).GT.min_sechiba) THEN
             DO jst=1,nstm
                cvs_over_veg(ji,jv,jst) = corr_veg_soil(ji,jv,jst)/vegtot(ji) / veget_max(ji,jv)
             ENDDO
          ENDIF
       END DO
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             humrel(ji,jv)=humrel(ji,jv)+humrelv(ji,jv,jst)*soiltile(ji,jst) &
                  & * cvs_over_veg(ji,jv,jst)*vegtot(ji)
             humrel(ji,jv)=MAX(humrel(ji,jv),zero)
          END DO
       END DO
    END DO

    DO jst=1,nstm
       DO ji=1,kjpindex
          ! We compute here a mean k for the 'litter' used for reinfiltration from floodplains of ponds
          !
          IF ( tmc_litter(ji,jst) < tmc_litter_res(ji,jst)) THEN
             i = imin
          ELSE
             tmc_litter_ratio = (tmc_litter(ji,jst)-tmc_litter_res(ji,jst)) / &
                  & (tmc_litter_sat(ji,jst)-tmc_litter_res(ji,jst))
             i= MAX(MIN(INT((imax-imin)*tmc_litter_ratio)+imin, imax-1), imin)
          ENDIF
          !
          !
          k_tmp = MAX(k_lin(i,1,njsc(ji))*ks(njsc(ji)), zero)
          k_litt(ji) = k_litt(ji) + soiltile(ji,jst) * SQRT(k_tmp)

       ENDDO
    ENDDO
!   write(*,*) 'hydrol 5578'

    DO jst=1,nstm        

       DO ji=1,kjpindex
          litterhumdiag(ji) = litterhumdiag(ji) + &
               & soil_wet_litter(ji,jst) * soiltile(ji,jst)

          tmc_litt_mea(ji) = tmc_litt_mea(ji) + &
               & tmc_litter(ji,jst) * soiltile(ji,jst) 

       END DO

       DO jd=1,nbdl
          DO ji=1,kjpindex
             if (ok_shumdiag_interpol) then
               DO jsl=1,nslm  
                  shumdiag(ji,jd)= shumdiag(ji,jd) + soil_wet(ji,jsl,jst)  &
                  & *frac_hydro_diag(jsl,jd)* &
                  & ((mcs(njsc(ji))-mcw(njsc(ji)))/(mcf(njsc(ji))-mcw(njsc(ji)))) * &
                  & soiltile(ji,jst)

               enddo !DO jd=1,nbdl 
             else !if (ok_shumdiag_interpol) then
                 jsl=jd
                 shumdiag(ji,jsl)= shumdiag(ji,jsl) + soil_wet(ji,jsl,jst) * &
                  & ((mcs(njsc(ji))-mcw(njsc(ji)))/(mcf(njsc(ji))-mcw(njsc(ji)))) * &
                  & soiltile(ji,jst)
             endif !if (ok_shumdiag_interpol) then
             jsl=jd
             soilmoist(ji,jsl)=soilmoist(ji,jsl)+mc(ji,jsl,jst)*soiltile(ji,jst)

             if (ok_shumdiag_perma) then
               DO jsl=1,nslm
                shumdiag_perma(ji,jd)= shumdiag_perma(ji,jd)  &
                  & + mc(ji,jsl,jst) *frac_hydro_diag(jsl,jd) &
                  & /mcs(njsc(ji))*soiltile(ji,jst)
               enddo !DO jsl=1,nslm
             endif !if (ok_shumdiag_perma) then
!             shumdiag(ji,jsl) = MAX(MIN(shumdiag(ji,jsl), un), zero) 
          END DO !DO ji=1,kjpindex
       END DO !DO jd=1,nbdl

    END DO !DO jst=1,nstm  

    ! First we compute swi (ALMIP requirement) - we assume here that dz is independant of jst (jst=1)
!Isa : correction nslm -> nbdl
    jst=1
    DO ji=1,kjpindex
       swi(ji) = swi(ji) + shumdiag(ji,1) * (dz(2,jst))/(deux*dpu_max*mille)
       DO jsl=2,nbdl-1 
          swi(ji) = swi(ji) + shumdiag(ji,jsl) * (dz(jsl,jst)+dz(jsl+1,jst))/(deux*dpu_max*mille)
       ENDDO
       jsl = nbdl
       swi(ji) = swi(ji) + shumdiag(ji,jsl) * (dz(jsl,jst))/(deux*dpu_max*mille)
    END DO

    ! For stomate we need shumdiag to be bounded by 0 and 1
    shumdiag(:,:) = MIN(shumdiag(:,:), un)
    shumdiag(:,:) = MAX(shumdiag(:,:), zero) 
if (.NOT.ok_shumdiag_perma) shumdiag_perma(:,:)=shumdiag(:,:)

    DO ji=1,kjpindex
       IF ( tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji) > zero ) THEN
          drysoil_frac(ji) = un + MAX( MIN( (tmc_litt_dry_mea(ji) - tmc_litt_mea(ji)) / &
               & (tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji)), zero), - un)
       ELSE
          drysoil_frac(ji) = zero
       ENDIF
    END DO
  
!    write(*,*) 'hydrol 5532'
  END SUBROUTINE hydrol_diag_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_waterbal 
!!
!>\BRIEF        Checks the water balance.
!!
!! DESCRIPTION  :
!! This routine checks the water balance. First it gets the total
!! amount of water and then it compares the increments with the fluxes.
!! The computation is only done over the soil area as over glaciers (and lakes?)
!! we do not have water conservation.
!! This verification does not make much sense in REAL*4 as the precision is the same as some
!! of the fluxes
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_waterbal

  SUBROUTINE hydrol_waterbal (kjpindex, index, first_call, dtradia, veget_max, totfrac_nobio, &
       & qsintveg, snow,snow_nobio, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, tot_melt, &
       & vevapwet, transpir, vevapnu, vevapsno, vevapflo, floodout, runoff, drainage)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                        :: kjpindex     !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index        !! Indeces of the points on the map
    LOGICAL, INTENT (in)                               :: first_call   !! At which time is this routine called ?
    REAL(r_std), INTENT (in)                           :: dtradia      !! Time step in seconds
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max    !! Max Fraction of vegetation type 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio!! Total fraction of continental ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg     !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow         !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio !!Ice water balance
    !
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain  !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_snow  !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: returnflow   !! Water to the bottom
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinfiltration !! Water to the top
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: irrigation   !! Water from irrigation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: tot_melt     !! Total melt
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet     !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir     !! Transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapnu      !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapsno     !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapflo     !! Floodplains evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: floodout     !! flow out of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: runoff       !! complete runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: drainage     !! Drainage
   ! REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: stagnant     !! stagnant



    !! 0.2 Output variables


    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: ji
    REAL(r_std) :: watveg, delta_water
    !
    !
    !
    IF ( ALL( tot_water_beg(:) == val_exp ) ) THEN

       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))
          tot_water_beg(ji) = humtot(ji)*vegtot(ji) + watveg + snow(ji)&
               & + SUM(snow_nobio(ji,:))
       ENDDO

       tot_water_end(:) = tot_water_beg(:)
       tot_flux(:) = zero
       RETURN
    ELSE IF ( first_call ) THEN
       tot_water_end(:) = tot_water_beg(:)
       tot_flux(:) = zero


       RETURN
    ENDIF

    tot_water_end(:) = zero
    tot_flux(:) = zero

    !
    DO ji = 1, kjpindex
       !
       ! If the fraction of ice, lakes, etc. does not complement the vegetation fraction then we do not
       ! need to go any further
       !
       IF ( ABS(un - (totfrac_nobio(ji) + vegtot(ji))) .GT. allowed_err ) THEN
          WRITE(numout,*) 'HYDROL problem in vegetation or frac_nobio on point ', ji
          WRITE(numout,*) 'totfrac_nobio : ', totfrac_nobio(ji)
          WRITE(numout,*) 'vegetation fraction : ', vegtot(ji)
          waterbal_error=.TRUE.
          CALL ipslerr(2, 'hydrol_waterbal', 'We will STOP after hydrol_waterbal.','','')
       ENDIF
    ENDDO

    IF ( .NOT. waterbal_error ) THEN
       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))

          !Chloe++ Water balance verification with total mass water with runoff peatland processes
          IF (ok_stagnant) THEN
                tot_water_end(ji) = humtot(ji)*vegtot(ji) + watveg + &
                & snow(ji) + SUM(snow_nobio(ji,:)) + SUM(stagnant(ji,:)*mask_soiltile(ji,:))*vegtot(ji)
          ELSE
                tot_water_end(ji) = humtot(ji)*vegtot(ji) + watveg + &
                & snow(ji) + SUM(snow_nobio(ji,:))
          ENDIF

          tot_flux(ji) =  precip_rain(ji) + precip_snow(ji) + irrigation(ji) - &
               & SUM(vevapwet(ji,:)) - SUM(transpir(ji,:)) - vevapnu(ji) - vevapsno(ji) - vevapflo(ji) + &
               & floodout(ji) - runoff(ji) - drainage(ji) + returnflow(ji) + reinfiltration(ji)

           !write(*,*) 'Check_Waterbal activ�'
           !write(40,*) '',tot_flux(ji)
       ENDDO

      
       
       DO ji = 1, kjpindex
          !
          delta_water = tot_water_end(ji) - tot_water_beg(ji)
          !
          !
          !  Set some precision ! This is a wild guess and corresponds to what works on an IEEE machine
          !  under double precision (REAL*8).
          !

          !
          IF ( ABS(delta_water-tot_flux(ji)) .GT. deux*allowed_err ) THEN
             WRITE(numout,*) '------------------------------------------------------------------------- '
             WRITE(numout,*) 'HYDROL does not conserve water. The erroneous point is : ', ji
             WRITE(numout,*) 'Coord erroneous point', lalo(ji,:)
             WRITE(numout,*) 'The error in mm/s is :', (delta_water-tot_flux(ji))/dtradia, ' and in mm/dt : ', &
                  & delta_water-tot_flux(ji)
             WRITE(numout,*) 'delta_water : ', delta_water, ' tot_flux : ', tot_flux(ji)
             WRITE(numout,*) 'Actual and allowed error : ', ABS(delta_water-tot_flux(ji)), allowed_err
             WRITE(numout,*) 'vegtot : ', vegtot(ji)
             WRITE(numout,*) 'precip_rain : ', precip_rain(ji)
             WRITE(numout,*) 'precip_snow : ',  precip_snow(ji)
             WRITE(numout,*) 'Water from routing. Reinfiltration/returnflow/irrigation : ', reinfiltration(ji), &
                  & returnflow(ji),irrigation(ji)
             WRITE(numout,*) 'Total water in soil humtot:',  humtot(ji)
             WRITE(numout,*) 'mc:' , mc(ji,:,:)
             WRITE(numout,*) 'Water on vegetation watveg:', watveg
             WRITE(numout,*) 'Snow mass snow:', snow(ji)
             WRITE(numout,*) 'Snow mass on ice snow_nobio:', SUM(snow_nobio(ji,:))
             WRITE(numout,*) 'Melt water tot_melt:', tot_melt(ji)
             WRITE(numout,*) 'evapwet : ', vevapwet(ji,:)
             WRITE(numout,*) 'transpir : ', transpir(ji,:)
             WRITE(numout,*) 'evapnu, evapsno, evapflo: ', vevapnu(ji), vevapsno(ji), vevapflo(ji)
             WRITE(numout,*) 'drainage,runoff,floodout : ', drainage(ji),runoff(ji),floodout(ji)
             waterbal_error=.TRUE.
             CALL ipslerr(2, 'hydrol_waterbal', 'We will STOP after hydrol_waterbal.','','')
          ENDIF
          !
       ENDDO
       !
       ! Transfer the total water amount at the end of the current timestep top the begining of the next one.
       !
       tot_water_beg = tot_water_end

     
    ENDIF

  END SUBROUTINE hydrol_waterbal


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_alma 
!!
!>\BRIEF        This routine computes the changes in soil moisture and interception storage for the ALMA outputs.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_alma

  SUBROUTINE hydrol_alma (kjpindex, index, first_call, qsintveg, snow, snow_nobio, soilwet)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                        :: kjpindex     !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index        !! Indeces of the points on the map
    LOGICAL, INTENT (in)                              :: first_call   !! At which time is this routine called ?
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg     !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow         !! Snow water equivalent
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Water balance on ice, lakes, .. [Kg/m^2]

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: soilwet     !! Soil wetness

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: ji
    REAL(r_std) :: watveg
    !
    !
    !


    IF ( first_call ) THEN

       tot_watveg_beg(:) = zero
       tot_watsoil_beg(:) = zero
       snow_beg(:)        = zero
       !
       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))
          tot_watveg_beg(ji) = watveg
          tot_watsoil_beg(ji) = humtot(ji)
          snow_beg(ji)        = snow(ji)+ SUM(snow_nobio(ji,:))
       ENDDO
       !
       tot_watveg_end(:) = tot_watveg_beg(:)
       tot_watsoil_end(:) = tot_watsoil_beg(:)
       snow_end(:)        = snow_beg(:)

       RETURN

    ENDIF
    !
    ! Calculate the values for the end of the time step
    !
    tot_watveg_end(:) = zero
    tot_watsoil_end(:) = zero
    snow_end(:) = zero
    delintercept(:) = zero
    delsoilmoist(:) = zero
    delswe(:) = zero
    !
    DO ji = 1, kjpindex
       watveg = SUM(qsintveg(ji,:))
       tot_watveg_end(ji) = watveg
       tot_watsoil_end(ji) = humtot(ji)
       snow_end(ji) = snow(ji)+ SUM(snow_nobio(ji,:))
       !
       delintercept(ji) = tot_watveg_end(ji) - tot_watveg_beg(ji)
       delsoilmoist(ji) = tot_watsoil_end(ji) - tot_watsoil_beg(ji)
       delswe(ji)       = snow_end(ji) - snow_beg(ji)
       !
       !
    ENDDO
    !
    !
    ! Transfer the total water amount at the end of the current timestep top the begining of the next one.
    !
    tot_watveg_beg = tot_watveg_end
    tot_watsoil_beg = tot_watsoil_end
    snow_beg(:) = snow_end(:)
    !
    DO ji = 1,kjpindex
       IF ( mx_eau_var(ji) > 0 ) THEN
          soilwet(ji) = tot_watsoil_end(ji) / mx_eau_var(ji)
       ELSE
          soilwet(ji) = zero
       ENDIF
    ENDDO
   

  END SUBROUTINE hydrol_alma
  !
  SUBROUTINE calcule_temp_hydro(kjpindex, stempdiag, snow)

        INTEGER(i_std), INTENT(in)                               :: kjpindex 
	REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in) :: stempdiag
	REAL(r_std),DIMENSION (kjpindex), INTENT (in) :: snow
!locals
	INTEGER jh, jd, ji
 	REAL(r_std) :: snow_h, profil_froz, m, mc_used
    	REAL(r_std)  :: lev_diag, prev_diag, lev_prog, prev_prog
    	REAL(r_std), DIMENSION(nslm,nbdl) :: intfactt
    !
    !

do ji=1,kjpindex

 	snow_h = snow(ji)/sn_dens
	intfactt(:,:)=0.
        !
        prev_diag = snow_h
        DO jh = 1, nslm
	    if (jh.EQ.1) then
          	lev_diag = zz(2,nstm)/1000./2.+snow_h
	    elseif (jh.EQ.nslm) then
		lev_diag = zz(nslm,nstm)/1000.+snow_h

	    else
		lev_diag = zz(jh, nstm)/1000. &
                  & +(zz(jh+1,nstm)-zz(jh,nstm))/1000./2.+snow_h

	    endif
          prev_prog = 0.0
          DO jd = 1, nbdl
                   lev_prog = diaglev(jd)
                   if ((lev_diag.GT.diaglev(nbdl).and. &
                  & prev_diag.lt.diaglev(nbdl)-min_sechiba)) then
                        lev_diag=diaglev(nbdl)          
                   endif               
	           intfactt(jh,jd) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog),&
                             & 0.0)/(lev_diag-prev_diag)
            	   prev_prog = lev_prog
          ENDDO
          if (lev_diag.GT.diaglev(nbdl).and. &
                  & prev_diag.ge.diaglev(nbdl)-min_sechiba) intfactt(jh,nbdl)=1.
           prev_diag = lev_diag
        ENDDO
enddo
! 
!            WRITE(numout,*) 'temp_hydro -- temp_hydro --  temp_hydro' 
!            DO jh = 1, nslm
!               WRITE(numout,*) jh, '-', intfactt(jh,1:nbdl)
!            ENDDO
!            WRITE(numout,*) "SUM -- SUM -- SUM SUM -- SUM -- SUM"
!            DO jh = 1, nslm
!               WRITE(numout,*) jh, '-', SUM(intfactt(jh,1:nbdl))
!            ENDDO
!            WRITE(numout,*) 'temp_hydro -- temp_hydro --  temp_hydro' 

        !

   temp_hydro(:,:)=0.

    DO jd= 1, nbdl
      DO jh= 1, nslm
        DO ji = 1, kjpindex
          temp_hydro(ji,jh) = temp_hydro(ji,jh) + stempdiag(ji,jd)*intfactt(jh,jd)
        ENDDO
      ENDDO  
    ENDDO


END SUBROUTINE calcule_temp_hydro


SUBROUTINE calcule_frac_hydro_diag
!output : frac_hydro_diag
!locals
	INTEGER(i_std) :: jd, jh
        REAL(r_std)    :: prev_hydro, next_hydro, prev_diag, next_diag
!
		frac_hydro_diag(:,:)=0.
		prev_diag = 0.0

		DO jd = 1, nbdl 

			next_diag = diaglev(jd)
			prev_hydro = 0.0
			DO jh = 1, nslm
				if (jh.EQ.1) then
					next_hydro = zz(2,nstm)/1000./2.
				elseif (jh.EQ.nslm) then
					next_hydro = zz(nslm,nstm)/1000.
				else
					next_hydro = zz(jh, nstm)/1000.+(zz(jh+1,nstm)-zz(jh,nstm))/1000./2.
				endif
			    frac_hydro_diag(jh,jd) = MAX(MIN(next_hydro, next_diag)-MAX(prev_hydro, prev_diag), 0.)/(next_diag - prev_diag)
			    prev_hydro=next_hydro
			ENDDO

			prev_diag = next_diag
	    ENDDO

! do jd =1,nbdl
! write(*,*) "SUM(frac_hydro_diag(,jd))",jd, SUM(frac_hydro_diag(:,jd))
! enddo
! write(*,*) "last test", frac_hydro_diag(1,:)


  END SUBROUTINE calcule_frac_hydro_diag

  !
END MODULE hydrol
