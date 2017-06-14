! =================================================================================================================================
! MODULE 	: stomate_data
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "stomate_data" module defines the values about the PFT parameters. It will print
!! the values of the parameters for STOMATE in the standard outputs. 
!!
!!\n DESCRIPTION: None 
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months))
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_data.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_data

  ! modules used:

  USE constantes
  USE pft_parameters
  USE defprec
  

  IMPLICIT NONE

  INTEGER(i_std),PARAMETER :: ibare_sechiba = 1   !! Index for bare soil in Sechiba (unitless)

  INTEGER(i_std),SAVE :: bavard=1                 !! Level of online diagnostics in STOMATE (0-4, unitless)

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: hori_index     !! Move to Horizontal indices

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: horipft_index  !! Horizontal + PFT indices

  ! Land cover change

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip10_index   !! Horizontal + P10 indices
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip100_index  !! Horizontal + P100 indice
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip11_index   !! Horizontal + P11 indices
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip101_index  !! Horizontal + P101 indices
!Chloe:
 INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip171_index   !! horizontal + p171indices


  INTEGER(i_std),SAVE :: itime                 !! time step
  INTEGER(i_std),SAVE :: hist_id_stomate       !! STOMATE history file ID
  INTEGER(i_std),SAVE :: hist_id_stomate_IPCC  !! STOMATE history file ID for IPCC output
  INTEGER(i_std),SAVE :: rest_id_stomate       !! STOMATE restart file ID

  REAL(r_std),PARAMETER :: adapted_crit = 1. - ( 1. / euler ) !! critical value for being adapted (1-1/e) (unitless)
  REAL(r_std),PARAMETER :: regenerate_crit = 1. / euler       !! critical value for being regenerative (1/e) (unitless)


  ! private & public routines

  PUBLIC data

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: data
!!
!>\BRIEF         This routine defines the values of the PFT parameters. It will print the values of the parameters for STOMATE
!!               in the standard outputs of ORCHIDEE. 
!!
!! DESCRIPTION : This routine defines PFT parameters. It initializes the pheno_crit structure by tabulated parameters.\n
!!               Some initializations are done for parameters. The SLA is calculated according *to* Reich et al (1992).\n
!!               Another formulation by Reich et al(1997) could be used for the computation of the SLA.
!!               The geographical coordinates might be used for defining some additional parameters
!!               (e.g. frequency of anthropogenic fires, irrigation of agricultural surfaces, etc.). \n
!!               For the moment, this possibility is not used. \n
!!               The specifc leaf area (SLA) is calculated according Reich et al, 1992 by :
!!               \latexonly
!!               \input{stomate_data_SLA.tex}
!!               \endlatexonly
!!               The sapling (young) biomass for trees and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_tree.tex}
!!               \endlatexonly
!!               The sapling biomass for grasses and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_grass.tex}
!!               \endlatexonly
!!               The critical stem diameter is given by the following formula :
!!               \latexonly
!!               \input{stomate_data_stem_diameter.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months)) 
!!                   \latexonly
!!                   \input{stomate_data_SLA_Reich_97.tex}
!!                   \endlatexonly
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!! - Reich PB, Walters MB, Ellsworth DS, (1992), Leaf life-span in relation to leaf, plant and 
!! stand characteristics among diverse ecosystems. Ecological Monographs, Vol 62, pp 365-392.
!! - Reich PB, Walters MB, Ellsworth DS (1997) From tropics to tundra: global convergence in plant 
!!  functioning. Proc Natl Acad Sci USA, 94:13730 13734
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE data (npts, lalo)


    !! 0. Variables and parameter declaration


    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                   :: npts    !! [DISPENSABLE] Domain size (unitless)
    REAL(r_std),DIMENSION (npts,2), INTENT (in)  :: lalo    !! [DISPENSABLE] Geographical coordinates (latitude,longitude)

    !! 0.4 Local variables

    INTEGER(i_std)                               :: j       !! Index (unitless)
    REAL(r_std)                                  :: alpha   !! alpha's : (unitless)
    REAL(r_std)                                  :: dia     !! stem diameter (m)
    REAL(r_std)                                  :: csa_sap !! Crown specific area sapling @tex $(m^2.ind^{-1})$ @endtex

!_ ================================================================================================================================

    IF ( bavard .GE. 1 ) WRITE(numout,*) 'data: PFT characteristics'

    !- pheno_gdd_crit
    pheno_gdd_crit(:,1) = pheno_gdd_crit_c(:)
    pheno_gdd_crit(:,2) = pheno_gdd_crit_b(:)         
    pheno_gdd_crit(:,3) = pheno_gdd_crit_a(:) 
    !
    !- senescence_temp
    senescence_temp(:,1) = senescence_temp_c(:)
    senescence_temp(:,2) = senescence_temp_b(:)
    senescence_temp(:,3) = senescence_temp_a(:)
    !
    !- maint_resp_slope
    maint_resp_slope(:,1) = maint_resp_slope_c(:)              
    maint_resp_slope(:,2) = maint_resp_slope_b(:)
    maint_resp_slope(:,3) = maint_resp_slope_a(:)
    !
    !-coeff_maint_zero
    coeff_maint_zero(:,ileaf) = cm_zero_leaf(:)
    coeff_maint_zero(:,isapabove) = cm_zero_sapabove(:)
    coeff_maint_zero(:,isapbelow) = cm_zero_sapbelow(:)
    coeff_maint_zero(:,iheartabove) = cm_zero_heartabove(:)
    coeff_maint_zero(:,iheartbelow) = cm_zero_heartbelow(:)
    coeff_maint_zero(:,iroot) = cm_zero_root(:)
    coeff_maint_zero(:,ifruit) = cm_zero_fruit(:)
    coeff_maint_zero(:,icarbres) = cm_zero_carbres(:)


    IF ( bavard .GE. 1 ) WRITE(numout,*) 'data: PFT characteristics'

    DO j = 2,nvm ! Loop over # PFTS 

       IF ( bavard .GE. 1 ) WRITE(numout,'(a,i3,a,a)') '    > PFT#',j,': ', PFT_name(j)

       !
       ! 1 tree? (true/false)
       !

       IF ( leaf_tab(j) .LE. 2 ) THEN
          tree(j) = .TRUE.
       ELSE
          tree(j) = .FALSE.
       ENDIF !( leaf_tab(j) .LE. 2 )

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       tree: (::tree) ', tree(j)

       !
       ! 2 flamability (0-1, unitless)
       !

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       litter flamability (::flam) :', flam(j)

       !
       ! 3 fire resistance (unitless)
       !

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       fire resistance (::resist) :', resist(j)

       !
       ! 4 specific leaf area per mass carbon = 2 * sla / dry mass (m^2.gC^{-1})
       !

       ! SZ: Reich et al, 1992 find no statistically significant differences between broadleaved and coniferous
       ! forests, specifically, the assumption that grasses grow needles is not justified. Replacing the function
       ! with the one based on Reich et al. 1997. Given that sla=100cm2/gDW at 9 months, sla is:
       ! sla=exp(5.615-0.46*ln(leaflon in months))

       ! Oct 2010 : sla values are prescribed by values given by N.Viovy 

       ! includes conversion from 
       !!       sla(j) = 2. * 1e-4 * EXP(5.615 - 0.46 * log(12./leaflife_tab(j)))
       !!\latexonly
       !!\input{stomate_data_SLA.tex}
       !!\endlatexonly
!       IF ( leaf_tab(j) .EQ. 2 ) THEN
!
!          ! needle leaved tree
!          sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!
!       ELSE
!
!          ! broad leaved tree or grass (Reich et al 1992)
!          sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!
!       ENDIF

!$      IF ( leaf_tab(j) .EQ. 1 ) THEN
!$
!$        ! broad leaved tree
!$
!$        sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!$
!$      ELSE
!$
!$        ! needle leaved or grass (Reich et al 1992)
!$
!$        sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!$
!$      ENDIF
!$
!$      IF ( ( leaf_tab(j) .EQ. 2 ) .AND. ( pheno_type_tab(j) .EQ. 2 ) ) THEN
!$
!$        ! summergreen needle leaf
!$
!$        sla(j) = 1.25 * sla(j)
!$
!$      ENDIF

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       specific leaf area (m**2/gC) (::sla):', sla(j), 12./leaflife_tab(j)

       !
       ! 5 sapling characteristics
       !

       IF ( tree(j) ) THEN

          !> 5.1 trees

          !!\latexonly
          !!\input{stomate_data_sapl_tree.tex}
          !!\endlatexonly

          alpha = alpha_tree

          bm_sapl(j,ileaf) = &
               &     ((bm_sapl_leaf(1)*pipe_tune1*(mass_ratio_heart_sap *bm_sapl_leaf(2)*sla(j)/(pi*pipe_k1)) & 
               &     **bm_sapl_leaf(3))/sla(j))**bm_sapl_leaf(4)

          IF ( pheno_type(j) .NE. 1 ) THEN
             ! not evergreen
             bm_sapl(j,icarbres) = bm_sapl_carbres * bm_sapl(j,ileaf)
          ELSE
             bm_sapl(j,icarbres) = zero
          ENDIF ! (pheno_type_tab(j) .NE. 1 )

          csa_sap = bm_sapl(j,ileaf) / ( pipe_k1 / sla(j) )

          dia = (mass_ratio_heart_sap * csa_sap * dia_coeff(1) / pi ) ** dia_coeff(2)

          bm_sapl(j,isapabove) = &
               bm_sapl_sapabove * pipe_density * csa_sap * pipe_tune2 * dia ** pipe_tune3
          bm_sapl(j,isapbelow) = bm_sapl(j,isapabove)

          bm_sapl(j,iheartabove) = bm_sapl_heartabove * bm_sapl(j,isapabove)
          bm_sapl(j,iheartbelow) = bm_sapl_heartbelow * bm_sapl(j,isapbelow)

       ELSE

          !> 5.2 grasses

          !!\latexonly
          !!\input{stomate_data_sapl_grass.tex}
          !!\endlatexonly

          alpha = alpha_grass

          IF ( natural(j) ) THEN
             bm_sapl(j,ileaf) = init_sapl_mass_leaf_nat / sla(j)
          ELSE
             bm_sapl(j,ileaf) = init_sapl_mass_leaf_agri / sla(j)
          ENDIF

          bm_sapl(j,icarbres) = init_sapl_mass_carbres *bm_sapl(j,ileaf)

          bm_sapl(j,isapabove) = zero
          bm_sapl(j,isapbelow) = zero

          bm_sapl(j,iheartabove) = zero
          bm_sapl(j,iheartbelow) = zero

       ENDIF !( tree(j) )

       bm_sapl(j,iroot) = init_sapl_mass_root * (1./alpha) * bm_sapl(j,ileaf)

       bm_sapl(j,ifruit) = init_sapl_mass_fruit  * bm_sapl(j,ileaf)

       IF ( bavard .GE. 1 ) THEN
          WRITE(numout,*) '       sapling biomass (gC):'
          WRITE(numout,*) '         leaves: (::bm_sapl(j,ileaf))',bm_sapl(j,ileaf)
          WRITE(numout,*) '         sap above ground: (::bm_sapl(j,ispabove)):',bm_sapl(j,isapabove)
          WRITE(numout,*) '         sap below ground: (::bm_sapl(j,isapbelow))',bm_sapl(j,isapbelow)
          WRITE(numout,*) '         heartwood above ground: (::bm_sapl(j,iheartabove))',bm_sapl(j,iheartabove)
          WRITE(numout,*) '         heartwood below ground: (::bm_sapl(j,iheartbelow))',bm_sapl(j,iheartbelow)
          WRITE(numout,*) '         roots: (::bm_sapl(j,iroot))',bm_sapl(j,iroot)
          WRITE(numout,*) '         fruits: (::bm_sapl(j,ifruit))',bm_sapl(j,ifruit)
          WRITE(numout,*) '         carbohydrate reserve: (::bm_sapl(j,icarbres))',bm_sapl(j,icarbres)
       ENDIF !( bavard .GE. 1 ) 

       !
       ! 6 migration speed (m/year)
       !

       IF ( tree(j) ) THEN

          migrate(j) = migrate_tree

       ELSE

          ! can be any value as grasses are, per *definition*, everywhere (big leaf).
          migrate(j) = migrate_grass

       ENDIF !( tree(j) )

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       migration speed (m/year): (::migrate(j))', migrate(j)

       !
       ! 7 critical stem diameter: beyond this diameter, the crown area no longer
       !     increases (m)
       !

       IF ( tree(j) ) THEN

          !!\latexonly
          !!\input{stomate_data_stem_diameter.tex}
          !!\endlatexonly

          maxdia(j) = ( ( pipe_tune4 / ((pipe_tune2*pipe_tune3)/(maxdia_coeff(1)**pipe_tune3)) ) &
               ** ( un / ( pipe_tune3 - un ) ) ) * maxdia_coeff(2)
          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year

       ELSE

          maxdia(j) = undef
          cn_sapl(j)=1

       ENDIF !( tree(j) )

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       critical stem diameter (m): (::maxdia(j))', maxdia(j)

       !
       ! 8 Coldest tolerable temperature (K)
       !

       IF ( ABS( tmin_crit(j) - undef ) .GT. min_stomate ) THEN
          tmin_crit(j) = tmin_crit(j) + ZeroCelsius
       ELSE
          tmin_crit(j) = undef
       ENDIF 

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       coldest tolerable temperature (K): (::tmin_crit(j))', tmin_crit(j)

       !
       ! 9 Maximum temperature of the coldest month: need to be below this temperature
       !      for a certain time to regrow leaves next spring *(vernalization)* (K)
       !

       IF ( ABS ( tcm_crit(j) - undef ) .GT. min_stomate ) THEN
          tcm_crit(j) = tcm_crit(j) + ZeroCelsius
       ELSE
          tcm_crit(j) = undef
       ENDIF

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       vernalization temperature (K): (::tcm_crit(j))', tcm_crit(j)

       !
       ! 10 critical values for phenology
       !

       ! 10.1 model used

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       phenology model used: (::pheno_model(j)) ',pheno_model(j)

       ! 10.2 growing degree days. What kind of gdd is meant (i.e. threshold 0 or -5 deg C
       !        or whatever), depends on how this is used in stomate_phenology.


       IF ( ( bavard .GE. 1 ) .AND. ( ALL(pheno_gdd_crit(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical GDD is a function of long term T (C): (::gdd)'
          WRITE(numout,*) '          ',pheno_gdd_crit(j,1), &
               ' + T *',pheno_gdd_crit(j,2), &
               ' + T^2 *',pheno_gdd_crit(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( pheno_model(j) .EQ. 'moigdd' ) .OR. &
            ( pheno_model(j) .EQ. 'humgdd' )       ) .AND. &
            ( ANY(pheno_gdd_crit(j,:) .EQ. undef) )                      ) THEN
          STOP 'problem with phenology parameters, critical GDD. (::pheno_model)'
       ENDIF

       ! 10.3 number of growing days

       IF ( ( bavard .GE. 1 ) .AND. ( ngd_crit(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical NGD: (::ngd_crit(j))', ngd_crit(j)

       ! 10.4 critical temperature for ncd vs. gdd function in phenology (C)

       IF ( ( bavard .GE. 1 ) .AND. ( ncdgdd_temp(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical temperature for NCD vs. GDD (C): (::ncdgdd_temp(j))', &
            ncdgdd_temp(j)

       ! 10.5 humidity fractions (0-1, unitless)

       IF ( ( bavard .GE. 1 ) .AND. ( hum_frac(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical humidity fraction: (::hum_frac(j))', &
            &  hum_frac(j)

       ! 10.6 minimum time during which there was no photosynthesis (days)

       IF ( ( bavard .GE. 1 ) .AND. ( lowgpp_time(j) .NE. undef ) ) &
            WRITE(numout,*) '         minimum dormance duration (d): (::lowgpp_time(j))', &
       & lowgpp_time(j)

       ! 10.7 minimum time elapsed since moisture minimum (days)

       IF ( ( bavard .GE. 1 ) .AND. ( hum_min_time(j) .NE. undef ) ) &
            WRITE(numout,*) '         time to wait after moisture min (d): (::hum_min_time(j))', &
        &    hum_min_time(j)

       !
       ! 11 critical values for senescence
       !

       ! 11.1 type of senescence

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       type of senescence: (::senescence_type(j))',senescence_type(j)

       ! 11.2 critical temperature for senescence (C)

       IF ( ( bavard .GE. 1 ) .AND. ( ALL(senescence_temp(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical temperature for senescence (C) is'
          WRITE(numout,*) '          a function of long term T (C): (::senescence_temp)'
          WRITE(numout,*) '          ',senescence_temp(j,1), &
               ' + T *',senescence_temp(j,2), &
               ' + T^2 *',senescence_temp(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'cold' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )      ) .AND. &
            ( ANY(senescence_temp(j,:) .EQ. undef ) )           ) THEN
          STOP 'problem with senescence parameters, temperature. (::senescence_type)'
       ENDIF

       ! 11.3 critical relative moisture availability for senescence

       IF ( ( bavard .GE. 1 ) .AND. ( senescence_hum(j) .NE. undef ) ) &
            WRITE(numout,*)  ' max. critical relative moisture availability for' 
            WRITE(numout,*)  ' senescence: (::senescence_hum(j))',  &
            & senescence_hum(j)

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( senescence_hum(j) .EQ. undef )                   ) THEN
          STOP 'problem with senescence parameters, humidity.(::senescence_type)'
       ENDIF

       ! 14.3 relative moisture availability above which there is no moisture-related
       !      senescence (0-1, unitless)

       IF ( ( bavard .GE. 1 ) .AND. ( nosenescence_hum(j) .NE. undef ) ) &
            WRITE(numout,*) '         relative moisture availability above which there is' 
            WRITE(numout,*) '             no moisture-related senescence: (::nosenescence_hum(j))', &
            &  nosenescence_hum(j)

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( nosenescence_hum(j) .EQ. undef )                   ) THEN
          STOP 'problem with senescence parameters, humidity. (::senescence_type)'
       ENDIF

       !
       ! 12 sapwood -> heartwood conversion time (days)
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       sapwood -> heartwood conversion time (d): (::tau_sap(j))', tau_sap(j)

       !
       ! 13 fruit lifetime (days)
       !

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       fruit lifetime (d): (::tau_fruit(j))', tau_fruit(j)

       !
       ! 14 length of leaf death (days)
       !      For evergreen trees, this variable determines the lifetime of the leaves.
       !      Note that it is different from the value given in leaflife_tab.
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       length of leaf death (d): (::leaffall(j))', leaffall(j)

       !
       ! 15 maximum lifetime of leaves (days)
       !

       IF ( ( bavard .GE. 1 ) .AND. ( leafagecrit(j) .NE. undef ) ) &
            WRITE(numout,*) '       critical leaf age (d): (::leafagecrit(j))', leafagecrit(j)

       !
       ! 16 time constant for leaf age discretisation (days)
       !

       leaf_timecst(j) = leafagecrit(j) / REAL( nleafages,r_std )

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       time constant for leaf age discretisation (d): (::leaf_timecst(j))', &
            leaf_timecst(j)

       !
       ! 17 minimum lai, initial (m^2.m^{-2})
       !

       IF ( tree(j) ) THEN
          lai_initmin(j) = lai_initmin_tree
       ELSE
          lai_initmin(j) = lai_initmin_grass
       ENDIF !( tree(j) )

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       initial LAI: (::lai_initmin(j))', lai_initmin(j)

       !
       ! 19 maximum LAI (m^2.m^{-2})
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       critical LAI above which no leaf allocation: (::lai_max(j))', lai_max(j)

       !
       ! 20 fraction of primary leaf and root allocation put into reserve (0-1, unitless)
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       reserve allocation factor: (::ecureuil(j))', ecureuil(j)

       !
       ! 21 maintenance respiration coefficient (g/g/day) at 0 deg C
       !

       IF ( bavard .GE. 1 ) THEN

          WRITE(numout,*) '       maintenance respiration coefficient (g/g/day) at 0 deg C:'
          WRITE(numout,*) '         . leaves: (::coeff_maint_zero(j,ileaf))',coeff_maint_zero(j,ileaf)
          WRITE(numout,*) '         . sapwood above ground: (::coeff_maint_zero(j,isapabove)) ',&
                        & coeff_maint_zero(j,isapabove)
          WRITE(numout,*) '         . sapwood below ground: (::coeff_maint_zero(j,isapbelow))  ',&
                       & coeff_maint_zero(j,isapbelow)
          WRITE(numout,*) '         . heartwood above ground: (::coeff_maint_zero(j,iheartabove)) ',&
                       & coeff_maint_zero(j,iheartabove)
          WRITE(numout,*) '         . heartwood below ground: (::coeff_maint_zero(j,iheartbelow)) ',&
                       & coeff_maint_zero(j,iheartbelow)
          WRITE(numout,*) '         . roots: (::coeff_maint_zero(j,iroot))',coeff_maint_zero(j,iroot)
          WRITE(numout,*) '         . fruits: (::coeff_maint_zero(j,ifruit)) ',coeff_maint_zero(j,ifruit)
          WRITE(numout,*) '         . carbohydrate reserve: (::coeff_maint_zero(j,icarbres)) ',&
                       & coeff_maint_zero(j,icarbres)

       ENDIF !( bavard .GE. 1 )

       !
       ! 22 parameter for temperature sensitivity of maintenance respiration
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       temperature sensitivity of maintenance respiration (1/K) is'
       WRITE(numout,*) '          a function of long term T (C): (::maint_resp_slope)'
       WRITE(numout,*) '          ',maint_resp_slope(j,1),' + T *',maint_resp_slope(j,2), &
            ' + T^2 *',maint_resp_slope(j,3)

       !
       ! 23 natural ?
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       Natural: (::natural(j))', natural(j)

       !
       ! 24 Vcmax et Vjmax (umol.m^{-2}.s^{-1}) 
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       Maximum rate of carboxylation: (::vcmax_opt(j))', vcmax_opt(j)

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       Maximum rate of RUbp regeneration: (::vjmax_opt(j))', vjmax_opt(j)

       !
       ! 25 constants for photosynthesis temperatures
       !

       IF ( bavard .GE. 1 ) THEN
          WRITE(numout,*) '       min. temperature for photosynthesis as a function of long term T (C):'
          WRITE(numout,*) '       (::tphoto_min_c(j))'
          WRITE(numout,*) '          ',tphoto_min_c(j), &
               ' + T*',tphoto_min_b(j), &
               ' + T^2*',tphoto_min_a(j)
          WRITE(numout,*) '       opt. temperature for photosynthesis as a function of long term T (C):' 
          WRITE(numout,*) '       (::tphoto_opt_c(j))'
          WRITE(numout,*) '          ',tphoto_opt_c(j), &
               ' + T*',tphoto_opt_b(j), &
               ' + T^2*',tphoto_opt_a(j)
          WRITE(numout,*) '       max. temperature for photosynthesis as a function of long term T (C):' 
          WRITE(numout,*)'        (::tphoto_max_c(j))'
          WRITE(numout,*) '          ',tphoto_max_c(j), &
               ' + T*',tphoto_max_b(j), &
               ' + T^2*',tphoto_max_a(j)


          !
          ! 26 Properties
          !

          WRITE(numout,*) '       Slope of the gs/A relation: (::gsslope(j))', gsslope(j)
          WRITE(numout,*) '       Intercept of the gs/A relation: (::gsoffset(j))', gsoffset(j)
          WRITE(numout,*) '       C4 photosynthesis: (::is_c4(j))', is_c4(j)
          WRITE(numout,*) '       Depth constant for root profile (m): (::1./humcste(j))', 1./humcste(j)

       ENDIF !( bavard .GE. 1 ) 

       !
       ! 27 extinction coefficient of the Monsi and Saeki (1953) relationship 
       !
       IF ( bavard .GE. 1 ) THEN
          WRITE(numout,*) '       extinction coefficient: (::ext_coeff(j))', ext_coeff(j)
       ENDIF !( bavard .GE. 1 )

       !
       ! 28 check coherence between tree definitions
       !      this is not absolutely necessary (just security)
       !

       IF ( tree(j) .NEQV. is_tree(j) ) THEN
          STOP 'Definition of tree/not tree not coherent (::tree)'
       ENDIF

    ENDDO ! Loop over # PFTS 

    !
    ! 29 time scales for phenology and other processes (in days)
    !

    tau_longterm = coeff_tau_longterm * one_year

    IF ( bavard .GE. 1 ) THEN

       WRITE(numout,*) '   > time scale for ''monthly'' moisture availability (d): (::tau_hum_month)', &
            tau_hum_month
       WRITE(numout,*) '   > time scale for ''weekly'' moisture availability (d): (::tau_hum_week)', &
           tau_hum_week
       WRITE(numout,*) '   > time scale for ''monthly'' 2 meter temperature (d): (::tau_t2m_month)', &
            tau_t2m_month
       WRITE(numout,*) '   > time scale for ''weekly'' 2 meter temperature (d): (::tau_t2m_week)', &
            tau_t2m_week
       WRITE(numout,*) '   > time scale for ''weekly'' GPP (d): (::tau_gpp_week)', &
            tau_gpp_week
       WRITE(numout,*) '   > time scale for ''monthly'' soil temperature (d): (::tau_tsoil_month)', &
            tau_tsoil_month
       WRITE(numout,*) '   > time scale for ''monthly'' soil humidity (d): (::tau_soilhum_month)', &
            tau_soilhum_month
       WRITE(numout,*) '   > time scale for vigour calculations (y): (::tau_longterm / one_year)', &
            tau_longterm / one_year

    ENDIF !( bavard .GE. 1 ) 

    !
    ! 30 fraction of allocatable biomass which is lost as growth respiration (0-1, unitless)
    !

    IF ( bavard .GE. 1 ) &
         WRITE(numout,*) '   > growth respiration fraction: (::frac_growthresp)', frac_growthresp

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving data'

  END SUBROUTINE data

END MODULE stomate_data
