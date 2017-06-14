! =================================================================================================================================
! MODULE       : stomate_turnover.f90
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module manages the end of the growing season and calculates herbivory and turnover of leaves, fruits, fine roots.
!! 
!!\n DESCRIPTION: This subroutine calculates leaf senescence due to climatic conditions or as a 
!! function of leaf age and new LAI, and subsequent turnover of the different plant biomass compartments (sections 1 to 6), 
!! herbivory (section 7), fruit turnover for trees (section 8) and sapwood conversion (section 9). 
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_turnover.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_turnover

  ! modules used:

  USE ioipsl
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC turn, turn_clear

  LOGICAL, SAVE                          :: firstcall = .TRUE.           !! first call (true/false)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : turn_clear
!!
!>\BRIEF        Set flag ::firstcall to .TRUE., and therefore activate section 1 
!!              of subroutine turn which writes a message to the output.
!!                
!_ ================================================================================================================================

  SUBROUTINE turn_clear
    firstcall=.TRUE.
  END SUBROUTINE turn_clear


!! ================================================================================================================================
!! SUBROUTINE    : turn
!!
!>\BRIEF         Calculate turnover of leaves, roots, fruits and sapwood due to aging or climatic 
!!               induced senescence. Calculate herbivory.
!!
!! DESCRIPTION : This subroutine determines the turnover of leaves and fine roots (and stems for grasses)
!! and simulates following processes:
!! 1. Mean leaf age is calculated from leaf ages of separate leaf age classes. Should actually 
!!    be recalculated at the end of the routine, but it does not change too fast. The mean leaf 
!!    age is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_lma_update_eqn1.tex}
!!    \endlatexonly
!!    \n 
!! 2. Meteorological senescence: the detection of the end of the growing season and shedding 
!!    of leaves, fruits and fine roots due to unfavourable meteorological conditions.
!!    The model distinguishes three different types of "climatic" leaf senescence, that do not 
!!    change the age structure: sensitivity to cold temperatures, to lack of water, or both. 
!!    If meteorological conditions are fulfilled, a flag ::senescence is set to TRUE. Note
!!    that evergreen species do not experience climatic senescence.
!!    Climatic senescence is triggered by sensitivity to cold temperatures where the critical 
!!    temperature for senescence is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_temp_crit_eqn2.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to lack of water availability where the 
!!    moisture availability critical level is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_moist_crit_eqn3.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to temperature or to lack of water where
!!    critical temperature and moisture availability are calculated as above.\n
!!    Trees in climatic senescence lose their fine roots at the same rate as they lose their leaves. 
!!    The rate of biomass loss of both fine roots and leaves is presribed through the equation:
!!    \latexonly
!!    \input{turnover_clim_senes_biomass_eqn4.tex}
!!    \endlatexonly
!!    \n
!!    with ::leaffall(j) a PFT-dependent time constant which is given in 
!!    ::stomate_constants. In grasses, leaf senescence is extended to the whole plant 
!!    (all carbon pools) except to its carbohydrate reserve.    
!! 3. Senescence due to aging: the loss of leaves, fruits and  biomass due to aging
!!    At a certain age, leaves fall off, even if the climate would allow a green plant
!!    all year round. Even if the meteorological conditions are favorable for leaf maintenance,
!!    plants, and in particular, evergreen trees, have to renew their leaves simply because the 
!!    old leaves become inefficient. Roots, fruits (and stems for grasses) follow leaves. 
!!    The ??senescence?? rate varies with leaf age. Note that plant is not declared senescent 
!!    in this case (wchich is important for allocation: if the plant loses leaves because of 
!!    their age, it can renew them). The leaf turnover rate due to aging of leaves is calculated
!!    using the following equation:
!!    \latexonly
!!    \input{turnover_age_senes_biomass_eqn5.tex}
!!    \endlatexonly
!!    \n
!!    Drop all leaves if there is a very low leaf mass during senescence. After this, the biomass 
!!    of different carbon pools both for trees and grasses is set to zero and the mean leaf age 
!!    is reset to zero. Finally, the leaf fraction and leaf age of the different leaf age classes 
!!    is set to zero. For deciduous trees: next to leaves, also fruits and fine roots are dropped.
!!    For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected:
!! 4. Update the leaf biomass, leaf age class fraction and the LAI
!!    Older leaves will fall more frequently than younger leaves and therefore the leaf age 
!!    distribution needs to be recalculated after turnover. The fraction of biomass in each 
!!    leaf class is updated using the following equation:
!!    \latexonly
!!    \input{turnover_update_LeafAgeDistribution_eqn6.tex}
!!    \endlatexonly
!!    \n
!! 5. Simulate herbivory activity and update leaf and fruits biomass. Herbivore activity 
!!    affects the biomass of leaves and fruits as well as stalks (only for grasses).
!!    However, herbivores do not modify leaf age structure.
!! 6. Calculates fruit turnover for trees. Trees simply lose their fruits with a time 
!!    constant ::tau_fruit(j), that is set to 90 days for all PFTs in ::stomate_constants 
!! 7. Convert sapwood to heartwood for trees and update heart and softwood above and 
!!    belowground biomass. Sapwood biomass is converted into heartwood biomass 
!!    with a time constant tau ::tau_sap(j) of 1 year. Note that this biomass conversion 
!!    is not added to "turnover" as the biomass is not lost. For the updated heartwood, 
!!    the sum of new heartwood above and new heartwood below after converting sapwood to 
!!    heartwood, is saved as ::hw_new(:). Creation of new heartwood decreases the age of 
!!    the plant ??carbon?? with a factor that is determined by: old heartwood ::hw_old(:) 
!!    divided by the new heartwood ::hw_new(:)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLES: ::Biomass of leaves, fruits, fine roots and sapwood above (latter for grasses only), 
!!                        ::Update LAI, ::Update leaf age distribution with new leaf age class fraction 
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199. 
!! - McNaughton, S. J., M. Oesterheld, D. A. Frank and K. J. Williams (1989), 
!! Ecosystem-level patterns of primary productivity and herbivory in terrestrial habitats, 
!! Nature, 341, 142-144, 1989. 
!! - Sitch, S., C. Huntingford, N. Gedney, P. E. Levy, M. Lomas, S. L. Piao, , Betts, R., Ciais, P., Cox, P., 
!! Friedlingstein, P., Jones, C. D., Prentice, I. C. and F. I. Woodward : Evaluation of the terrestrial carbon  
!! cycle, future plant geography and climate-carbon cycle feedbacks using 5 dynamic global vegetation 
!! models (dgvms), Global Change Biology, 14(9), 2015–2039, 2008. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale=0.5]{turnover_flowchart_1.png}
!! \includegraphics[scale=0.5]{turnover_flowchart_2.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE turn (npts, dt, PFTpresent, &
       herbivores, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       moiavail_week, moiavail_month, tlong_ref, t2m_month, t2m_week, veget_max, &
       leaf_age, leaf_frac, age, lai, biomass, &
       turnover, senescence,turnover_time)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 

    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size - number of grid cells 
                                                                                       !! (unitless) 
    REAL(r_std), INTENT(in)                                    :: dt                   !! time step (dt_days)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                   :: PFTpresent           !! PFT exists (true/false)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! last year's maximum moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! last year's minimum moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "weekly" moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "monthly" moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tlong_ref            !! "long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "weekly" 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: veget_max            !! "maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground (unitless) 

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out)       :: turnover             !! Turnover @tex ($gC m^{-2}$) @endtex
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                  :: senescence           !! is the plant senescent? (true/false) 
                                                                                       !! (interesting only for deciduous trees: 
                                                                                       !! carbohydrate reserve) 
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age             !! age of the leaves (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac            !! fraction of leaves in leaf age class 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! age (years)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: lai                  !! leaf area index @tex ($m^2 m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)     :: biomass              !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! turnover_time of grasses (days)

    !! 0.4 Local  variables

    REAL(r_std), DIMENSION(npts,nvm)                           :: leaf_meanage         !! mean age of the leaves (days)
    REAL(r_std), DIMENSION(npts)                               :: dturnover            !! Intermediate variable for turnover ??
                                                                                       !! @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts)                               :: moiavail_crit        !! critical moisture availability, function 
                                                                                       !! of last year's moisture availability 
                                                                                       !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                               :: tl                   !! long term annual mean temperature, (C)
    REAL(r_std), DIMENSION(npts)                               :: t_crit               !! critical senescence temperature, function 
                                                                                       !! of long term annual temperature (K) 
    LOGICAL, DIMENSION(npts)                                   :: shed_rest            !! shed the remaining leaves? (true/false)
    REAL(r_std), DIMENSION(npts)                               :: sapconv              !! Sapwood conversion @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts)                               :: hw_old               !! old heartwood mass @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts)                               :: hw_new               !! new heartwood mass @tex ($gC m^{-2}$) 
                                                                                       !! @endtex 
    REAL(r_std), DIMENSION(npts)                               :: lm_old               !! old leaf mass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nleafages)                     :: delta_lm             !! leaf mass change for each age class @tex 
                                                                                       !! ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts)                               :: turnover_rate        !! turnover rate (unitless) 
    REAL(r_std), DIMENSION(npts,nvm)                           :: leaf_age_crit        !! critical leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm)                           :: new_turnover_time    !! instantaneous turnover time (days)
    INTEGER(i_std)                                             :: j,m                  !! Index (unitless)

!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering turnover'

    !! 1. first call - output messages

    IF ( firstcall ) THEN

       WRITE(numout,*) 'turnover:'

       WRITE(numout,*) ' > minimum mean leaf age for senescence (days) (::min_leaf_age_for_senescence) : '&
            ,min_leaf_age_for_senescence

       firstcall = .FALSE.


    ENDIF

    !! 2. Initializations 

    !! 2.1 set output to zero
    turnover(:,:,:) = zero
    new_turnover_time(:,:) = zero
    senescence(:,:) = .FALSE.

    !! 2.2 Recalculate mean leaf age
    !      Mean leaf age is recalculated from leaf ages of separate leaf age classes. Should actually be recalculated at the 
    !      end of this routine, but it does not change too fast.
    !      The mean leaf age is calculated using the following equation:
    !      \latexonly
    !      \input{turnover_lma_update_eqn1.tex}
    !      \endlatexonly
    !      \n
    leaf_meanage(:,:) = zero

    DO m = 1, nleafages
       leaf_meanage(:,:) = leaf_meanage(:,:) + leaf_age(:,:,m) * leaf_frac(:,:,m)
    ENDDO

    !! 3. Climatic senescence

    !     Three different types of "climatic" leaf senescence,
    !     that do not change the age structure. 
    DO j = 2,nvm ! Loop over # PFTs

       !! 3.1 Determine if there is climatic senescence. 
       !      The climatic senescence can be of three types:
       !      sensitivity to cold temperatures, to lack of water, or both. If meteorological conditions are 
       !      fulfilled, a flag senescence is set to TRUE.
       !      Evergreen species do not experience climatic senescence.

       SELECT CASE ( senescence_type(j) )

       CASE ( 'cold' )

          !! 3.1.1 Summergreen species: Climatic senescence is triggered by sensitivity to cold temperatures
          !        Climatic senescence is triggered by sensitivity to cold temperatures as follows: 
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the monthly temperature is low enough (i.e. when monthly temperature ::t2m_month(:) is below a critical 
          !        temperature ::t_crit(:), which is calculated in this module),
          !        AND the temperature tendency is negative (i.e. when weekly temperatures ::t2m_week(:) are lower than monthly 
          !        temperatures ::t2m_month(:))
          !        If these conditions are met, senescence is set to TRUE.
          !
          !        The critical temperature for senescence is calculated using the following equation:
          !        \latexonly
          !        \input{turnover_temp_crit_eqn2.tex}
          !        \endlatexonly
          !        \n
          !
          ! Critical temperature for senescence may depend on long term annual mean temperature
          tl(:) = tlong_ref(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(j,1) + &
               tl(:) * senescence_temp(j,2) + &
               tl(:)*tl(:) * senescence_temp(j,3)

          WHERE ( ( biomass(:,j,ileaf) .GT. zero ) .AND. &
               ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( t2m_month(:) .LT. t_crit(:) ) .AND. ( t2m_week(:) .LT. t2m_month(:) ) )


             senescence(:,j) = .TRUE.

          ENDWHERE

       CASE ( 'dry' )

          !! 3.1.2 Raingreen species: Climatic senescence is triggered by sensitivity to lack of water availability 
          !        Climatic senescence is triggered by sensitivity to lack of water availability as follows:  
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the moisture availability drops below a critical level (i.e. when weekly moisture availability 
          !        ::moiavail_week(:,j) is below a critical moisture availability ::moiavail_crit(:),
          !        which is calculated in this module), 
          !        If these conditions are met, senescence is set to TRUE.
          !
          !        The moisture availability critical level is calculated using the following equation:
          !        \latexonly
          !        \input{turnover_moist_crit_eqn3.tex}
          !        \endlatexonly
          !        \n
          moiavail_crit(:) = &
               MIN( MAX( minmoiavail_lastyear(:,j) + hum_frac(j) * &
               ( maxmoiavail_lastyear(:,j) - minmoiavail_lastyear(:,j) ), &
               senescence_hum(j) ), &
               nosenescence_hum(j) )

          WHERE ( ( biomass(:,j,ileaf) .GT. zero ) .AND. &
               ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
               ( moiavail_week(:,j) .LT. moiavail_crit(:) ) )

             senescence(:,j) = .TRUE.

          ENDWHERE

       CASE ( 'mixed' )

          !! 3.1.3 Mixed criterion: Climatic senescence is triggered by sensitivity to temperature or to lack of water  
          !        Climatic senescence is triggered by sensitivity to temperature or to lack of water availability as follows:
          !        If biomass is large enough (i.e. when it is greater than zero), 
          !        AND (i.e. when leaf mean age is above a certain PFT-dependent treshold ::min_leaf_age_for_senescence,
          !        which is given in ::stomate_constants),      
          !        AND the moisture availability drops below a critical level (i.e. when weekly moisture availability 
          !        ::moiavail_week(:,j) is below a critical moisture availability ::moiavail_crit(:), calculated in this module), 
          !        OR 
          !        the monthly temperature is low enough (i.e. when monthly temperature ::t2m_month(:) is below a critical 
          !        temperature ::t_crit(:), calculated in this module),
          !        AND the temperature tendency is negative (i.e. when weekly temperatures ::t2m_week(:) are lower than 
          !        monthly temperatures ::t2m_month(:)).
          !        If these conditions are met, senescence is set to TRUE.
          moiavail_crit(:) = &
               MIN( MAX( minmoiavail_lastyear(:,j) + hum_frac(j) * &
               (maxmoiavail_lastyear(:,j) - minmoiavail_lastyear(:,j) ), &
               senescence_hum(j) ), &
               nosenescence_hum(j) )

          tl(:) = tlong_ref(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(j,1) + &
               tl(:) * senescence_temp(j,2) + &
               tl(:)*tl(:) * senescence_temp(j,3)

          IF ( tree(j) ) THEN
             ! critical temperature for senescence may depend on long term annual mean temperature
             WHERE ( ( biomass(:,j,ileaf) .GT. zero ) .AND. &
                  ( leaf_meanage(:,j) .GT. min_leaf_age_for_senescence(j) ) .AND. &
                  ( ( moiavail_week(:,j) .LT. moiavail_crit(:) ) .OR. &
                  ( ( t2m_month(:) .LT. t_crit(:) ) .AND. ( t2m_week(:) .LT. t2m_month(:) ) ) ) )
                senescence(:,j) = .TRUE.
             ENDWHERE
          ELSE
             new_turnover_time(:,j)=max_turnover_time(j)+ new_turnover_time_ref
             WHERE ((moiavail_week(:,j) .LT. moiavail_month(:,j))&
                  .AND. (moiavail_week(:,j) .LT. nosenescence_hum(j)))
                new_turnover_time(:,j)=max_turnover_time(j) * &
                     (un - nosenescence_hum(j)+moiavail_week(:,j)) * &
                     (un - nosenescence_hum(j)+moiavail_week(:,j)) + &
                     min_turnover_time(j)
             ENDWHERE

             WHERE (new_turnover_time(:,j) .GT. turnover_time(:,j)*1.1)
                new_turnover_time(:,j)=max_turnover_time(j)+ new_turnover_time_ref
             ENDWHERE

             turnover_time(:,j)=(turnover_time(:,j)*dt_turnover_time/dt+new_turnover_time(:,j))/(dt_turnover_time/dt + un)

          ENDIF


       !! Evergreen species do not experience climatic senescence
       CASE ( 'none' )

          
       !! In case no climatic senescence type is recognized.
       CASE default

          WRITE(numout,*) '  turnover: don''t know how to treat this PFT.'
          WRITE(numout,*) '  number (::j) : ',j
          WRITE(numout,*) '  senescence type (::senescence_type(j)) : ',senescence_type(j)

          STOP

       END SELECT

       !! 3.2 Drop leaves and roots, plus stems and fruits for grasses

       IF ( tree(j) ) THEN

          !! 3.2.1 Trees in climatic senescence lose their fine roots at the same rate as they lose their leaves. 
          !        The rate of biomass loss of both fine roots and leaves is presribed through the equation:
          !        \latexonly
          !        \input{turnover_clim_senes_biomass_eqn4.tex}
          !        \endlatexonly
          !        \n
          !         with ::leaffall(j) a PFT-dependent time constant which is given in ::stomate_constants),
          WHERE ( senescence(:,j) )

             turnover(:,j,ileaf) = biomass(:,j,ileaf) * dt / leaffall(j)
             turnover(:,j,iroot) = biomass(:,j,iroot) * dt / leaffall(j)

             biomass(:,j,ileaf) = biomass(:,j,ileaf) - turnover(:,j,ileaf)
             biomass(:,j,iroot) = biomass(:,j,iroot) - turnover(:,j,iroot)

          ENDWHERE

       ELSE

          !! 3.2.2 In grasses, leaf senescence is extended to the whole plant 
          !        In grasses, leaf senescence is extended to the whole plant (all carbon pools) except to its
          !        carbohydrate reserve.      
          WHERE (turnover_time(:,j) .LT. max_turnover_time(j)) 
             turnover(:,j,ileaf) = biomass(:,j,ileaf) * dt / turnover_time(:,j)
             turnover(:,j,isapabove) = biomass(:,j,isapabove) * dt / turnover_time(:,j)
             turnover(:,j,iroot) = biomass(:,j,iroot) * dt / turnover_time(:,j) 
             turnover(:,j,ifruit) = biomass(:,j,ifruit) * dt / turnover_time(:,j)
          ELSEWHERE
             turnover(:,j,ileaf)= zero
             turnover(:,j,isapabove) = zero
             turnover(:,j,iroot) = zero
             turnover(:,j,ifruit) = zero
          ENDWHERE
          biomass(:,j,ileaf) = biomass(:,j,ileaf) - turnover(:,j,ileaf)
          biomass(:,j,isapabove) = biomass(:,j,isapabove) - turnover(:,j,isapabove)
          biomass(:,j,iroot) = biomass(:,j,iroot) - turnover(:,j,iroot)
          biomass(:,j,ifruit) = biomass(:,j,ifruit) - turnover(:,j,ifruit)

       ENDIF      ! tree/grass

    ENDDO        ! loop over PFTs

    !! 4. Leaf fall
    !     At a certain age, leaves fall off, even if the climate would allow a green plant
    !     all year round. Even if the meteorological conditions are favorable for leaf maintenance,
    !     plants, and in particular, evergreen trees, have to renew their leaves simply because the 
    !     old leaves become inefficient.   
    !     Roots, fruits (and stems) follow leaves. The decay rate varies with leaf age.
    !     Note that plant is not declared senescent in this case (wchich is important for allocation:
    !     if the plant loses leaves because of their age, it can renew them).
    !
    !     The leaf turnover rate due to aging of leaves is calculated using the following equation:
    !     \latexonly
    !     \input{turnover_age_senes_biomass_eqn5.tex}
    !     \endlatexonly
    !     \n
    DO j = 2,nvm ! Loop over # PFTs

       !! save old leaf mass
       lm_old(:) = biomass(:,j,ileaf)

       !! initialize leaf mass change in age class
       delta_lm(:,:) = zero

       IF ( tree(j) ) THEN

          !! 4.1 Trees: leaves, roots, fruits roots and fruits follow leaves.

          !! 4.1.1 Critical age: prescribed for trees
          leaf_age_crit(:,j) = leafagecrit(j)

          !! 4.1.2 Loop over leaf age classes
          DO m = 1, nleafages
             turnover_rate(:) = zero
             WHERE ( leaf_age(:,j,m) .GT. leaf_age_crit(:,j)/2. )

                turnover_rate(:) =  &
                     MIN( 0.99_r_std, dt / ( leaf_age_crit(:,j) * &
                     ( leaf_age_crit(:,j) / leaf_age(:,j,m) )**quatre ) )

                dturnover(:) = biomass(:,j,ileaf) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,ileaf) = turnover(:,j,ileaf) + dturnover(:)
                biomass(:,j,ileaf) = biomass(:,j,ileaf) - dturnover(:)

                ! save leaf mass change
                delta_lm(:,m) = - dturnover(:)

                dturnover(:) = biomass(:,j,iroot) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,iroot) = turnover(:,j,iroot) + dturnover(:)
                biomass(:,j,iroot) = biomass(:,j,iroot) - dturnover(:)

                dturnover(:) = biomass(:,j,ifruit) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,ifruit) = turnover(:,j,ifruit) + dturnover(:)
                biomass(:,j,ifruit) = biomass(:,j,ifruit) - dturnover(:)

             ENDWHERE

          ENDDO

       ELSE

          !! 4.2 Grasses: leaves, roots, fruits, sap follow leaves.

          !! 4.2.1 Critical leaf age depends on long-term temperature
          !        Critical leaf age depends on long-term temperature
          !        generally, lower turnover in cooler climates.
          leaf_age_crit(:,j) = &
               MIN( leafagecrit(j) * leaf_age_crit_coeff(1) , &
               MAX( leafagecrit(j) * leaf_age_crit_coeff(2) , &
               leafagecrit(j) - leaf_age_crit_coeff(3) * &
               ( tlong_ref(:)-ZeroCelsius - leaf_age_crit_tref ) ) )

          ! 4.2.2 Loop over leaf age classes
          DO m = 1, nleafages

             WHERE ( leaf_age(:,j,m) .GT. leaf_age_crit(:,j)/2. )

                turnover_rate(:) =  &
                     MIN( 0.99_r_std, dt / ( leaf_age_crit(:,j) * &
                     ( leaf_age_crit(:,j) / leaf_age(:,j,m) )**quatre ) )

                dturnover(:) = biomass(:,j,ileaf) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,ileaf) = turnover(:,j,ileaf) + dturnover(:)
                biomass(:,j,ileaf) = biomass(:,j,ileaf) - dturnover(:)

                ! save leaf mass change
                delta_lm(:,m) = - dturnover(:)

                dturnover(:) = biomass(:,j,isapabove) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,isapabove) = turnover(:,j,isapabove) + dturnover(:)
                biomass(:,j,isapabove) = biomass(:,j,isapabove) - dturnover(:)

                dturnover(:) = biomass(:,j,iroot) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,iroot) = turnover(:,j,iroot) + dturnover(:)
                biomass(:,j,iroot) = biomass(:,j,iroot) - dturnover(:)

                dturnover(:) = biomass(:,j,ifruit) * leaf_frac(:,j,m) * turnover_rate(:)
                turnover(:,j,ifruit) = turnover(:,j,ifruit) + dturnover(:)
                biomass(:,j,ifruit) = biomass(:,j,ifruit) - dturnover(:)

             ENDWHERE
          ENDDO
       ENDIF       ! tree/grass ?

       !! 4.3 Recalculate the fraction of leaf biomass in each leaf age class.
       !      Older leaves will fall more fast than younger leaves and therefore 
       !      the leaf age distribution needs to be recalculated after turnover. 
       !      The fraction of biomass in each leaf class is updated using the following equation:
       !      \latexonly
       !      \input{turnover_update_LeafAgeDistribution_eqn6.tex}
       !      \endlatexonly
       !      \n
       !
       !      new fraction = new leaf mass of that fraction / new total leaf mass
       !                   = (old fraction*old total leaf mass ::lm_old(:) + biomass change of that fraction ::delta_lm(:,m)  ) /
       !                     new total leaf mass ::biomass(:,j,ileaf
       DO m = 1, nleafages

          WHERE ( biomass(:,j,ileaf) .GT. zero )
             leaf_frac(:,j,m) = ( leaf_frac(:,j,m)*lm_old(:) + delta_lm(:,m) ) / biomass(:,j,ileaf)
          ELSEWHERE
             leaf_frac(:,j,m) = zero
          ENDWHERE

       ENDDO

    ENDDO         ! loop over PFTs

    !! 5. New (provisional) LAI 
    !     ::lai(:,j) is determined from the leaf biomass ::biomass(:,j,ileaf) and the 
    !     specific leaf surface :: sla(j) (m^2 gC^{-1})
    !     The leaf area index is updated using the following equation:
    !     \latexonly
    !     \input{turnover_update_LAI_eqn7.tex}
    !     \endlatexonly
    !     \n

    !    lai(:,ibare_sechiba) = zero
    !    DO j = 2, nvm ! Loop over # PFTs
    !        lai(:,j) = biomass(:,j,ileaf) * sla(j)
    !    ENDDO

    !! 6. Definitely drop all leaves if there is a very low leaf mass during senescence.

    !     Both for deciduous trees and grasses same conditions are checked:
    !     If biomass is large enough (i.e. when it is greater than zero), 
    !     AND when senescence is set to true
    !     AND the leaf biomass drops below a critical minimum biomass level (i.e. when it is lower than half
    !     the minimum initial LAI ::lai_initmin(j) divided by the specific leaf area ::sla(j),
    !     ::lai_initmin(j) is set to 0.3 in stomate_data.f90 and sla is a constant that is set to 0.015366 m2/gC), 
    !     If these conditions are met, the flag ::shed_rest(:) is set to TRUE.
    !
    !     After this, the biomass of different carbon pools both for trees and grasses is set to zero
    !     and the mean leaf age is reset to zero.
    !     Finally, the leaf fraction and leaf age of the different leaf age classes is set to zero.
    DO j = 2,nvm ! Loop over # PFTs

       shed_rest(:) = .FALSE.

       !! 6.1 For deciduous trees: next to leaves, also fruits and fine roots are dropped 
       !      For deciduous trees: next to leaves, also fruits and fine roots are dropped: fruit ::biomass(:,j,ifruit) 
       !      and fine root ::biomass(:,j,iroot) carbon pools are set to zero.
       IF ( tree(j) .AND. ( senescence_type(j) .NE. 'none' ) ) THEN

          ! check whether we shed the remaining leaves
          WHERE ( ( biomass(:,j,ileaf) .GT. zero ) .AND. senescence(:,j) .AND. &
               ( biomass(:,j,ileaf) .LT. (lai_initmin(j) / 2.)/sla(j) )             )

             shed_rest(:) = .TRUE.

             turnover(:,j,ileaf)  = turnover(:,j,ileaf) + biomass(:,j,ileaf)
             turnover(:,j,iroot)  = turnover(:,j,iroot) + biomass(:,j,iroot)
             turnover(:,j,ifruit) = turnover(:,j,ifruit) + biomass(:,j,ifruit)

             biomass(:,j,ileaf)  = zero
             biomass(:,j,iroot)  = zero
             biomass(:,j,ifruit) = zero

             ! reset leaf age
             leaf_meanage(:,j) = zero

          ENDWHERE

       ENDIF

       !! 6.2 For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected: 
       !      For grasses: all aboveground carbon pools, except the carbohydrate reserves are affected: 
       !      fruit ::biomass(:,j,ifruit), fine root ::biomass(:,j,iroot) and sapwood above 
       !      ::biomass(:,j,isapabove) carbon pools are set to zero. 
       IF ( .NOT. tree(j) ) THEN

          ! Shed the remaining leaves if LAI very low.
          WHERE ( ( biomass(:,j,ileaf) .GT. zero ) .AND. senescence(:,j) .AND. &
               (  biomass(:,j,ileaf) .LT. (lai_initmin(j) / 2.)/sla(j) ))

             shed_rest(:) = .TRUE.

             turnover(:,j,ileaf) = turnover(:,j,ileaf) + biomass(:,j,ileaf)
             turnover(:,j,isapabove) = turnover(:,j,isapabove) + biomass(:,j,isapabove)
             turnover(:,j,iroot) = turnover(:,j,iroot) + biomass(:,j,iroot)
             turnover(:,j,ifruit) = turnover(:,j,ifruit) + biomass(:,j,ifruit)

             biomass(:,j,ileaf) = zero
             biomass(:,j,isapabove) = zero
             biomass(:,j,iroot) = zero
             biomass(:,j,ifruit) = zero

             ! reset leaf age
             leaf_meanage(:,j) = zero

          ENDWHERE

       ENDIF

       !! 6.3 Reset the leaf age structure: the leaf fraction and leaf age of the different leaf age classes is set to zero.
      
       DO m = 1, nleafages

          WHERE ( shed_rest(:) )

             leaf_age(:,j,m) = zero
             leaf_frac(:,j,m) = zero

          ENDWHERE

       ENDDO

    ENDDO          ! loop over PFTs
    
    !! 7. Herbivore activity: elephants, cows, gazelles but no lions.
 
    !     Herbivore activity affects the biomass of leaves and fruits as well 
    !     as stalks (only for grasses). Herbivore activity does not modify leaf 
    !     age structure. Herbivores ::herbivores(:,j) is the time constant of 
    !     probability of a leaf to be eaten by a herbivore, and is calculated in 
    !     ::stomate_season. following Mc Naughton et al. [1989].

    IF ( ok_herbivores ) THEN

       ! If the herbivore activity is allowed (if ::ok_herbivores is true, which is set in run.def), 
       ! remove the amount of biomass consumed by herbivory from the leaf biomass ::biomass(:,j,ileaf) and 
       ! the fruit biomass ::biomass(:,j,ifruit).
       ! The daily amount consumed equals the biomass multiplied by 1 day divided by the time constant ::herbivores(:,j).
       DO j = 2,nvm ! Loop over # PFTs

          IF ( tree(j) ) THEN

             !! For trees: only the leaves and fruit carbon pools are affected

             WHERE (biomass(:,j,ileaf) .GT. zero)
                ! added by shilong
                WHERE (herbivores(:,j).GT. zero)
                   dturnover(:) = biomass(:,j,ileaf) * dt / herbivores(:,j)
                   turnover(:,j,ileaf) = turnover(:,j,ileaf) + dturnover(:)
                   biomass(:,j,ileaf) = biomass(:,j,ileaf) - dturnover(:)

                   dturnover(:) = biomass(:,j,ifruit) * dt / herbivores(:,j)
                   turnover(:,j,ifruit) = turnover(:,j,ifruit) + dturnover(:)
                   biomass(:,j,ifruit) = biomass(:,j,ifruit) - dturnover(:)
                ENDWHERE
             ENDWHERE

          ELSE

             ! For grasses: all aboveground carbon pools are affected: leaves, fruits and sapwood above
             WHERE ( biomass(:,j,ileaf) .GT. zero )
                ! added by shilong
                WHERE (herbivores(:,j) .GT. zero)
                   dturnover(:) = biomass(:,j,ileaf) * dt / herbivores(:,j)
                   turnover(:,j,ileaf) = turnover(:,j,ileaf) + dturnover(:)
                   biomass(:,j,ileaf) = biomass(:,j,ileaf) - dturnover(:)

                   dturnover(:) = biomass(:,j,isapabove) * dt / herbivores(:,j)
                   turnover(:,j,isapabove) = turnover(:,j,isapabove) + dturnover(:)
                   biomass(:,j,isapabove) = biomass(:,j,isapabove) - dturnover(:)

                   dturnover(:) = biomass(:,j,ifruit) * dt / herbivores(:,j)
                   turnover(:,j,ifruit) = turnover(:,j,ifruit) + dturnover(:)
                   biomass(:,j,ifruit) = biomass(:,j,ifruit) - dturnover(:)
                ENDWHERE

             ENDWHERE

          ENDIF  ! tree/grass?

       ENDDO    ! loop over PFTs

    ENDIF ! end herbivores

    !! 8. Fruit turnover for trees

    !     Fruit turnover for trees: trees simply lose their fruits with a time constant ::tau_fruit(j), 
    !     that is set to 90 days for all PFTs in ::stomate_constants
    DO j = 2,nvm ! Loop over # PFTs

       IF ( tree(j) ) THEN

          dturnover(:) = biomass(:,j,ifruit) * dt / tau_fruit(j)
          turnover(:,j,ifruit) = turnover(:,j,ifruit) + dturnover(:)
          biomass(:,j,ifruit) = biomass(:,j,ifruit) - dturnover(:)


       ENDIF

    ENDDO       ! loop over PFTs

    !! 9 Conversion of sapwood to heartwood both for aboveground and belowground sapwood and heartwood.

    !   Following LPJ (Sitch et al., 2003), sapwood biomass is converted into heartwood biomass 
    !   with a time constant tau ::tau_sap(j) of 1 year.
    !   Note that this biomass conversion is not added to "turnover" as the biomass is not lost!
    DO j = 2,nvm ! Loop over # PFTs

       IF ( tree(j) ) THEN

          !! For the recalculation of age in 9.2 (in case the vegetation is not dynamic ie. ::control%ok_dgvm is FALSE), 
          !! the heartwood above and below is stored in ::hw_old(:).
          IF ( .NOT. control%ok_dgvm ) THEN
             hw_old(:) = biomass(:,j,iheartabove) + biomass(:,j,iheartbelow)
          ENDIF

          !! 9.1 Calculate the rate of sapwood to heartwood conversion 
          !      Calculate the rate of sapwood to heartwood conversion with the time constant ::tau_sap(j) 
          !      and update aboveground and belowground sapwood ::biomass(:,j,isapabove) and ::biomass(:,j,isapbelow)
          !      and heartwood ::biomass(:,j,iheartabove) and ::biomass(:,j,iheartbelow).

          ! Above the ground
          sapconv(:) = biomass(:,j,isapabove) * dt / tau_sap(j)
          biomass(:,j,isapabove) = biomass(:,j,isapabove) - sapconv(:)
          biomass(:,j,iheartabove) =  biomass(:,j,iheartabove) + sapconv(:)

          ! Below the ground
          sapconv(:) = biomass(:,j,isapbelow) * dt / tau_sap(j)
          biomass(:,j,isapbelow) = biomass(:,j,isapbelow) - sapconv(:)
          biomass(:,j,iheartbelow) =  biomass(:,j,iheartbelow) + sapconv(:)

          !! 9.2 If the vegetation is not dynamic, the age of the plant is decreased. 
          !      The updated heartwood, the sum of new heartwood above and new heartwood below after 
          !      converting sapwood to heartwood, is saved as ::hw_new(:) .
          !      Creation of new heartwood decreases the age of the plant with a factor that is determined by: 
          !      old heartwood ::hw_old(:) divided by the new heartwood ::hw_new(:)
          IF ( .NOT. control%ok_dgvm ) THEN

             hw_new(:) = biomass(:,j,iheartabove) + biomass(:,j,iheartbelow)

             WHERE ( hw_new(:) .GT. zero )

                age(:,j) = age(:,j) * hw_old(:)/hw_new(:)

             ENDWHERE

          ENDIF

       ENDIF

    ENDDO       ! loop over PFTs


    ! Write mean leaf age and time constant of probability of a leaf to be eaten by a herbivore 
    ! to the stomate output file.
    CALL histwrite (hist_id_stomate, 'LEAF_AGE', itime, &
         leaf_meanage, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'HERBIVORES', itime, &
         herbivores, npts*nvm, horipft_index)

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving turnover'

  END SUBROUTINE turn

END MODULE stomate_turnover
