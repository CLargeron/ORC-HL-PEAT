! =================================================================================================================================
! MODULE       : lpj_establish
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Establish pft's
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: 
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!! - Haxeltine, A. and I. C. Prentice (1996), BIOME3: An equilibrium
!!        terrestrial biosphere model based on ecophysiological constraints, 
!!        resource availability, and competition among plant functional types,
!!        Global Biogeochemical Cycles, 10(4), 693-709.\n
!! - Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/lpj_establish.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE lpj_establish

  ! modules used:
  USE ioipsl
  USE stomate_data
  USE constantes

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC establish,establish_clear

  LOGICAL, SAVE                          :: firstcall = .TRUE.           !! first call
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : fire_clear 
!!
!>\BRIEF       Set the firstcall flag to .TRUE. and activate initialization
!_ ================================================================================================================================

  SUBROUTINE establish_clear
    firstcall = .TRUE.
  END SUBROUTINE establish_clear


! =================================================================================================================================
! SUBROUTINE   : establish 
!
!>\BRIEF       Calculate sstablishment of new woody PFT and herbaceous PFTs
!!
!! DESCRIPTION : Establishments of new woody and herbaceous PFT are simulated. 
!! Maximum establishment rate (0.12) declines due to competition for light (space).
!! There are two establishment estimates: one for the for DGVM and one for the 
!! static cases.\n
!! In the case of DGVM, competitive process of establishment for the area of 
!! available space is represented using more detailed description compared with static 
!! one. Biomass and distribution of plant age are updated on the basis of changes 
!! in number of individuals. Finally excess sapwood of is converted to heartwood.
!!
!! \latexonly
!! \input{equation_lpj_establish.tex}
!! \endlatexonly
!! \n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)    :
!! Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!    dynamics in the modelling of terrestrial ecosystems: comparing two
!!    contrasting approaches within European climate space,
!!    Global Ecology and Biogeography, 10, 621-637.
!!
!! FLOWCHART       : 
!! \latexonly
!! \includegraphics[scale = 0.7]{establish.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE establish (npts, dt, PFTpresent, regenerate, &
       neighbours, resolution, need_adjacent, herbivores, &
       precip_annual, gdd0, lm_lastyearmax, &
       cn_ind, lai, avail_tree, avail_grass,  npp_longterm, &
       leaf_age, leaf_frac, &
       ind, biomass, age, everywhere, co2_to_bm,veget_max, woodmass_ind)
 
    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                                :: npts            !! Domain size - number of pixels (dimensionless)    
    REAL(r_std), INTENT(in)                                   :: dt              !! Time step of vegetation dynamics for stomate 
                                                                                 !! (days)            
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: PFTpresent      !! Is pft there (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: regenerate      !! Winter sufficiently cold (unitless)   
    INTEGER(i_std), DIMENSION(npts,8), INTENT(in)             :: neighbours      !! indices of the 8 neighbours of each grid point 
                                                                                 !! (unitless);  
                                                                                 !! [1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW]  
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                :: resolution      !! resolution at each grid point (m); 1=E-W, 2=N-S     
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: need_adjacent   !! in order for this PFT to be introduced, does it
                                                                                 !! have to be present in an adjacent grid box?  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: herbivores      !! time constant of probability of a leaf to 
                                                                                 !! be eaten by a herbivore (days)     
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: precip_annual   !! annual precipitation (mm year^{-1}) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: gdd0            !! growing degree days (degree C)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lm_lastyearmax  !! last year's maximum leaf mass for each PFT 
                                                                                 !! (gC m^{-2 })
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: cn_ind          !! crown area of individuals (m^2)        
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lai             !! leaf area index OF an individual plant 
                                                                                 !! (m^2 m^{-2})           
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: avail_tree      !! space availability for trees (unitless)     
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: avail_grass     !! space availability for grasses (unitless)  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: npp_longterm    !! longterm NPP, for each PFT (gC m^{-2})   
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: veget_max       !! "maximal" coverage fraction of a PFT 
                                                                                 !! (LAI -> infinity) on ground (unitless)
    !! 0.2 Output variables
    
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_age        !! leaf age (days)     
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac       !! fraction of leaves in leaf age class (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind             !! Number of individuals (individuals m^{-2})          
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)    :: biomass         !! biomass (gC m^{-2 })     
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age             !! mean age (years)       
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere      !! is the PFT everywhere in the grid box or very 
                                                                                 !! localized (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm       !! biomass up take for establishment i.e. 
                                                                                 !! pseudo-photosynthesis (gC m^{-2} day^{-1}) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: woodmass_ind    !! woodmass of the individual, needed to calculate
                                                                                 !! crownarea in lpj_crownarea (gC m^{-2 })  

    !! 0.4 Local variables

    REAL(r_std)                                               :: tau_eatup       !! time during which a sapling can be entirely 
                                                                                 !! eaten by herbivores (days) 
    REAL(r_std), DIMENSION(npts,nvm)                          :: fpc_nat         !! new fpc, foliage projective cover: fractional
                                                                                 !! coverage (unitless)       
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_climate_tree  !! maximum tree establishment rate, 
                                                                                              !! based on climate only (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_climate_grass !! maximum grass establishment rate,
                                                                                              !! based on climate only (unitless)
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_tree          !! maximum tree establishment rate, 
                                                                                              !! based on climate and fpc 
                                                                                              !! (unitless) 
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_grass         !! maximum grass establishment rate,
                                                                                              !! based on climate and fpc 
                                                                                              !! (unitless) 
    REAL(r_std), DIMENSION(npts)                              :: sumfpc          !! total natural fpc (unitless)
    REAL(r_std), DIMENSION(npts)                              :: fracnat         !! total fraction occupied by natural 
                                                                                 !! vegetation (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: sumfpc_wood     !! total woody fpc (unitless)    
    REAL(r_std), DIMENSION(npts)                              :: spacefight_tree !! for trees, measures the total concurrence for
                                                                                 !! available space (unitless)      
    REAL(r_std), DIMENSION(npts)                              :: spacefight_grass!! for grasses, measures the total concurrence 
                                                                                 !! for available space (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                          :: d_ind           !! change in number of individuals per time step 
                                                                                 !! (individuals m^{-2} day{-1})         
    REAL(r_std), DIMENSION(npts)                              :: bm_new          !! biomass increase (gC m^{-2 })        
    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)    
    REAL(r_std), DIMENSION(npts)                              :: b1              !! temporary variable           
    REAL(r_std), DIMENSION(npts)                              :: sm2             !! new sap mass (gC m^{-2 })          
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass of an individual (gC m^{-2})  
    REAL(r_std), DIMENSION(npts)                              :: leaf_mass_young !! carbon mass in youngest leaf age class 
                                                                                 !! (gC m^{-2})
    REAL(r_std), DIMENSION(npts)                              :: sm_at           !! ratio of hw(above) to total hw, sm(above) to 
                                                                                 !! total sm (unitless)         
    REAL(r_std), DIMENSION(npts)                              :: factor          !! reduction factor for establishment if many 
                                                                                 !! trees or grasses are present (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: total_bm_c      !! Total carbon mass for all pools (gC m^{-2})    
    REAL(r_std), DIMENSION(npts)                              :: total_bm_sapl   !! Total sappling biomass for all pools 
                                                                                 !! (gC m^{-2})  
    INTEGER(i_std)                                            :: nfrontx         !! from how many sides is the grid box invaded
                                                                                 !! (unitless?)   
    INTEGER(i_std)                                            :: nfronty         !! from how many sides is the grid box invaded
                                                                                 !! (unitless?)   
   !LOGICAL, DIMENSION(npts)                                  :: many_new        !! daily establishment rate is large compared to 
                                                                                 !! present number of individuals
    REAL(r_std), DIMENSION(npts)                              :: vn              !! flow due to new individuals veget_max after 
                                                                                 !! establishment, to get a proper estimate of 
                                                                                 !! carbon and nitrogen 
    REAL(r_std), DIMENSION(npts)                              :: lai_ind         !! lai on each PFT surface (m^2 m^{-2})   
    INTEGER(i_std)                                            :: i,j,k,m         !! indices (unitless)       
!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering establish'

  !! 1. messages and initialization

    ! Assumption: time durasion that sapling is completely eaten by hervioures is a half year?   
    ! No reference
    tau_eatup = one_year/2.

    !! 1.1 First call only
    IF ( firstcall ) THEN

       WRITE(numout,*) 'establish:'

       WRITE(numout,*) '   > time during which a sapling can be entirely eaten by herbivores (d): ', &
            tau_eatup

       firstcall = .FALSE.

    ENDIF

  !! 2. recalculate fpc

    IF (control%ok_dgvm) THEN
       fracnat(:) = un

       !! 2.1 Only natural part of the grid cell
       do j = 2,nvm ! Loop over # PFTs
          
          IF ( .NOT. natural(j) ) THEN
             fracnat(:) = fracnat(:) - veget_max(:,j)
          ENDIF
       ENDDO ! Loop over # PFTs
       
       sumfpc(:) = zero

       !! 2.2 Total natural fpc on grid
       !      The overall fractional coverage of a PFT in a grid is calculated here.
       !      FPC is related to mean individual leaf area index by the Lambert-Beer law.
       !      See Eq. (1) in tex file.\n
       DO j = 2,nvm ! Loop over # PFTs
          IF ( natural(j) ) THEN
             WHERE(fracnat(:).GT.min_stomate)
                WHERE (lai(:,j) == val_exp) 
                   fpc_nat(:,j) = cn_ind(:,j) * ind(:,j) / fracnat(:)
                ELSEWHERE
                   fpc_nat(:,j) = cn_ind(:,j) * ind(:,j) / fracnat(:) & 
                        * ( un - exp( - lm_lastyearmax(:,j) * sla(j) * ext_coeff(j) ) )
                ENDWHERE
             ENDWHERE

             WHERE ( PFTpresent(:,j) )
                sumfpc(:) = sumfpc(:) + fpc_nat(:,j)
             ENDWHERE
          ELSE

             fpc_nat(:,j) = zero

          ENDIF

       ENDDO ! Loop over # PFTs
       
       !! 2.3 Total woody fpc on grid and number of regenerative tree pfts
       !      Total woody FPC increases by adding new FPC.
       !      Under the condition that temperature in last winter is higher than a threshold, 
       !      woody plants is exposed in higher competitive environment.
       sumfpc_wood(:) = zero
       spacefight_tree(:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( tree(j) .AND. natural(j) ) THEN

             ! total woody fpc
             WHERE ( PFTpresent(:,j) )
                sumfpc_wood(:) = sumfpc_wood(:) + fpc_nat(:,j)
             ENDWHERE

             ! how many trees are competing? Count a PFT fully only if it is present
             ! on the whole grid box.
             WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit ) )
                spacefight_tree(:) = spacefight_tree(:) + everywhere(:,j)
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

       !! 2.4 Total number of natural grasses on grid\n
       !     Grass increment equals 'everywhere'\n
       spacefight_grass(:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( .NOT. tree(j) .AND. natural(j) ) THEN

             ! Count a PFT fully only if it is present on a grid.
             WHERE ( PFTpresent(:,j) )
                spacefight_grass(:) = spacefight_grass(:) + everywhere(:,j)
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

       !! 2.5 Maximum establishment rate, based on climate only\n
       WHERE ( ( precip_annual(:) .GE. precip_crit ) .AND. ( gdd0(:) .GE. gdd_crit_estab ) )

          estab_rate_max_climate_tree(:) = estab_max_tree ! 'estab_max_*'; see 'stomate_constants.f90'
          estab_rate_max_climate_grass(:) = estab_max_grass

       ELSEWHERE

          estab_rate_max_climate_tree(:) = zero
          estab_rate_max_climate_grass(:) = zero

       ENDWHERE

       !! 2.6 Reduce maximum tree establishment rate if many trees are present.
       !      In the original DGVM, this is done using a step function which yields a
       !      reduction by factor 4 if sumfpc_wood(i) .GT.  fpc_crit - 0.05.
       !      This can lead to small oscillations (without consequences however).
       !      Here, a steady linear transition is used between fpc_crit-0.075 and
       !      fpc_crit-0.025.
       !      factor(:) = 1. - 15. * ( sumfpc_wood(:) - (fpc_crit-.075))
       !      factor(:) = MAX( 0.25_r_std, MIN( 1._r_std, factor(:)))
       !      SZ modified according to Smith et al. 2001
       !      See Eq. (2) in header
       factor(:)=(un - exp(- establish_scal_fact * (un - sumfpc_wood(:))))*(un - sumfpc_wood(:))
       estab_rate_max_tree(:) = estab_rate_max_climate_tree(:) * factor(:)

       !! 2.7 Modulate grass establishment rate.
       !      If canopy is not closed (fpc < fpc_crit-0.05), normal establishment.
       !      If canopy is closed, establishment is reduced by a factor 4.
       !      Factor is linear between these two bounds.
       !      This is different from the original DGVM where a step function is
       !      used at fpc_crit-0.05 (This can lead to small oscillations,
       !      without consequences however).
       !      factor(:) = 1. - 15. * ( sumfpc(:) - (fpc_crit-.05))
       !      factor(:) = MAX( 0.25_r_std, MIN( 1._r_std, factor(:)))
       !      estab_rate_max_grass(:) = estab_rate_max_climate_grass(:) * factor(:)
       !      SZ modified to true LPJ formulation, grasses are only allowed in the
       !      fpc fraction not occupied by trees..., 080806
       !      estab_rate_max_grass(:)=MAX(0.98-sumfpc(:),zero)
       !      See Eq. (3) in header
       estab_rate_max_grass(:) = MAX(MIN(estab_rate_max_climate_grass(:), max_tree_coverage - sumfpc(:)),zero)

       !! 2.8 Longterm grass NPP for competition between C4 and C3 grasses
       !      to avoid equal veget_max, the idea is that more reestablishment
       !      is possible for the more productive PFT
       factor(:) = min_stomate
       DO j = 2,nvm ! Loop over # PFTs
          IF ( natural(j) .AND. .NOT.tree(j)) & 
               factor(:) = factor(:) + npp_longterm(:,j) * &
               lm_lastyearmax(:,j) * sla(j)
       ENDDO ! Loop over # PFTs

       !! 2.9 Establish natural PFTs
       d_ind(:,:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( natural(j) ) THEN ! only for natural PFTs

             !! 2.9.1 PFT expansion across the grid box. Not to be confused with areal coverage.
             IF ( treat_expansion ) THEN

                ! only treat plants that are regenerative and present and still can expand
                DO i = 1, npts ! Loop over # pixels - domain size

                   IF ( PFTpresent(i,j) .AND. &
                        ( everywhere(i,j) .LT. un ) .AND. &
                        ( regenerate(i,j) .GT. regenerate_crit ) ) THEN

                      ! from how many sides is the grid box invaded (separate x and y directions
                      ! because resolution may be strongly anisotropic)
                      ! For the moment we only look into 4 direction but that can be expanded (JP) 
                      nfrontx = 0
                      IF ( neighbours(i,3) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,3),j) .GT. 1.-min_stomate ) nfrontx = nfrontx+1
                      ENDIF
                      IF ( neighbours(i,7) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,7),j) .GT. 1.-min_stomate ) nfrontx = nfrontx+1
                      ENDIF

                      nfronty = 0
                      IF ( neighbours(i,1) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,1),j) .GT. 1.-min_stomate ) nfronty = nfronty+1
                      ENDIF
                      IF ( neighbours(i,5) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,5),j) .GT. 1.-min_stomate ) nfronty = nfronty+1
                      ENDIF
                      
                      everywhere(i,j) = &
                           everywhere(i,j) + migrate(j) * dt/one_year * &
                           ( nfrontx / resolution(i,1) + nfronty / resolution(i,2) )
                      
                      IF ( .NOT. need_adjacent(i,j) ) THEN
                         
                         ! in that case, we also assume that the PFT expands from places within
                         ! the grid box (e.g., oasis).
                         ! What is this equation? No reference.
                         everywhere(i,j) = &
                              everywhere(i,j) + migrate(j) * dt/one_year * &
                              2. * SQRT( pi*everywhere(i,j)/(resolution(i,1)*resolution(i,2)) )

                      ENDIF

                      everywhere(i,j) = MIN( everywhere(i,j), un )

                   ENDIF

                ENDDO ! Loop over # pixels - domain size

             ENDIF ! treat expansion?

             !! 2.9.2 Establishment rate
             !      - Is lower if the PFT is only present in a small part of the grid box
             !        (after its introduction), therefore multiplied by "everywhere".
             !      - Is divided by the number of PFTs that compete ("spacefight").
             !      - Is modulated by space availability (avail_tree, avail_grass).

             !! 2.9.2.1 present and regenerative trees
             IF ( tree(j) ) THEN

                WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit ) )
                   
                   
                   d_ind(:,j) = estab_rate_max_tree(:)*everywhere(:,j)/spacefight_tree(:) * &
                        avail_tree(:) * dt/one_year

                ENDWHERE

             !! 2.9.2.2 present and regenerative grasses
             ELSE

                WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit )  & 
                     .AND.factor(:).GT.min_stomate .AND. spacefight_grass(:).GT. min_stomate) 
                   
                   d_ind(:,j) = estab_rate_max_grass(:)*everywhere(:,j)/spacefight_grass(:) * &
                        MAX(min_stomate,npp_longterm(:,j)*lm_lastyearmax(:,j)*sla(j)/factor(:)) * fracnat(:) * dt/one_year
                   
                ENDWHERE
                
             ENDIF  ! tree/grass
             
          ENDIF ! if natural
       ENDDO  ! Loop over # PFTs
       
  !! 3. Lpj establishment in static case 

    !     Lpj establishment in static case, SZ 080806, account for real LPJ dynamics in
    !     prescribed vegetation, i.e. population dynamics within a given area of the grid cell.
    ELSE 

       d_ind(:,:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          WHERE(ind(:,j)*cn_ind(:,j).GT.min_stomate)
             lai_ind(:) = sla(j) * lm_lastyearmax(:,j)/(ind(:,j)*cn_ind(:,j))
          ELSEWHERE
             lai_ind(:) = zero
          ENDWHERE

          !! 3.1 For natural woody PFTs
          IF ( natural(j) .AND. tree(j)) THEN 

             ! See Eq. (4) in tex file.            
             fpc_nat(:,j) =  MIN(un, cn_ind(:,j) * ind(:,j) * & 
                  MAX( ( un - exp( - ext_coeff(j) * lai_ind(:) ) ), min_cover ) )


             WHERE (veget_max(:,j).GT.min_stomate.AND.ind(:,j).LE.2.)

                !! 3.1.1 Only establish into growing stands 
                !        Only establish into growing stands, ind can become very
                !        large in the static mode because LAI is very low in poor 
                !        growing conditions, favouring continuous establishment. 
                !        To avoid this a maximum IND is set. BLARPP: This should be
                !        replaced by a better stand density criteria.
                factor(:)=(un - exp(-establish_scal_fact * (un - fpc_nat(:,j))))*(un - fpc_nat(:,j))

                estab_rate_max_tree(:) = estab_max_tree * factor(:) 

                !! 3.1.2 do establishment for natural PFTs\n
                d_ind(:,j) = MAX( zero, estab_rate_max_tree(:) * dt/one_year)

             ENDWHERE

             !SZ: quickfix: to simulate even aged stand, uncomment the following lines...
             !where (ind(:,j) .LE. min_stomate)
             !d_ind(:,j) = 0.1 !MAX( 0.0, estab_rate_max_tree(:) * dt/one_year)
             WHERE (veget_max(:,j).GT.min_stomate .AND. ind(:,j).EQ.zero)
                d_ind(:,j) = ind_0_estab
             ENDWHERE

          !! 3.2 For natural grass PFTs
          ELSEIF ( natural(j) .AND. .NOT.tree(j)) THEN 

             WHERE (veget_max(:,j).GT.min_stomate)

                fpc_nat(:,j) =  cn_ind(:,j) * ind(:,j) * & 
                     MAX( ( un - exp( - ext_coeff(j) * lai_ind(:) ) ), min_cover )

                d_ind(:,j) = MAX(zero , (un - fpc_nat(:,j)) * dt/one_year )

             ENDWHERE

             WHERE (veget_max(:,j).GT.min_stomate .AND. ind(:,j).EQ. zero)
                d_ind(:,j) = ind_0_estab 
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

    ENDIF ! DGVM OR NOT

  !! 4. Biomass calculation

    DO j = 2,nvm ! Loop over # PFTs

       IF ( natural(j) ) THEN ! only for natural PFTs

          !! 4.1 Herbivores reduce establishment rate
          !      We suppose that saplings are vulnerable during a given time after establishment.
          !      This is taken into account by preventively reducing the establishment rate.
          IF ( ok_herbivores ) THEN

             d_ind(:,j) = d_ind(:,j) * EXP( - tau_eatup/herbivores(:,j) )

          ENDIF

          !! 4.2 Total biomass.
          !      Add biomass only if d_ind, over one year, is of the order of ind.
          !      save old leaf mass to calculate leaf age
          leaf_mass_young(:) = leaf_frac(:,j,1) * biomass(:,j,ileaf)

          ! total biomass of existing PFT to limit biomass added from establishment
          total_bm_c(:) = zero

          DO k = 1, nparts
             total_bm_c(:)=total_bm_c(:)+biomass(:,j,k)
          ENDDO
          IF(control%ok_dgvm) THEN
             vn(:) = veget_max(:,j)
          ELSE
             vn(:) = un
          ENDIF
          total_bm_sapl(:)=zero
          DO k = 1, nparts
             WHERE(d_ind(:,j).GT.min_stomate.AND.vn(:).GT.min_stomate)

                total_bm_sapl(:) = total_bm_sapl(:) + & 
                     bm_sapl(j,k) * d_ind(:,j) / vn(:)
             ENDWHERE
          ENDDO

          !! 4.3 Woodmass calculation

          !! 4.3.1 with DGVM
          IF(control%ok_dgvm) THEN

             ! SZ calculate new woodmass_ind and veget_max after establishment (needed for correct scaling!)
             ! essential correction for MERGE!
             IF(tree(j))THEN
                DO i=1,npts ! Loop over # pixels - domain size
                   IF((d_ind(i,j)+ind(i,j)).GT.min_stomate) THEN

                      IF((total_bm_c(i).LE.min_stomate) .OR. (veget_max(i,j) .LE. min_stomate)) THEN

                         ! new wood mass of PFT
                         woodmass_ind(i,j) = &
                              & (((biomass(i,j,isapabove) + biomass(i,j,isapbelow) &
                              & + biomass(i,j,iheartabove) + biomass(i,j,iheartbelow))*veget_max(i,j)) &
                              & + (bm_sapl(j,isapabove) + bm_sapl(j,isapbelow) &
                              & + bm_sapl(j,iheartabove) + bm_sapl(j,iheartbelow))*d_ind(i,j))/(ind(i,j) + d_ind(i,j))

                      ELSE
 
                         ! new biomass is added to the labile pool, hence there is no change 
                         ! in CA associated with establishment
                         woodmass_ind(i,j) = &
                              & (biomass(i,j,isapabove) + biomass(i,j,isapbelow) &
                              & +biomass(i,j,iheartabove) + biomass(i,j,iheartbelow))*veget_max(i,j) &
                              & /(ind(i,j) + d_ind(i,j))

                      ENDIF

                      ! new diameter of PFT
                      dia(i) = (woodmass_ind(i,j)/(pipe_density*pi/4.*pipe_tune2)) &
                           &                **(1./(2.+pipe_tune3))
                      vn(i) = (ind(i,j) + d_ind(i,j))*pipe_tune1*MIN(dia(i),maxdia(j))**pipe_tune_exp_coeff

                   ENDIF
                ENDDO ! Loop over # pixels - domain size
             ELSE ! for grasses, cnd=1, so the above calculation cancels
                vn(:) = ind(:,j) + d_ind(:,j)
             ENDIF

          !! 4.3.2 without DGVM (static)\n
          ELSE 
             DO i=1,npts ! Loop over # pixels - domain size
                IF(tree(j).AND.(d_ind(i,j)+ind(i,j)).GT.min_stomate) THEN
                   IF(total_bm_c(i).LE.min_stomate) THEN

                      ! new wood mass of PFT
                      woodmass_ind(i,j) = &
                           & (((biomass(i,j,isapabove) + biomass(i,j,isapbelow) &
                           & + biomass(i,j,iheartabove) + biomass(i,j,iheartbelow))) &
                           & + (bm_sapl(j,isapabove) + bm_sapl(j,isapbelow) &
                           & + bm_sapl(j,iheartabove) + bm_sapl(j,iheartbelow))*d_ind(i,j))/(ind(i,j)+d_ind(i,j))

                   ELSE
 
                      ! new biomass is added to the labile pool, hence there is no change 
                      ! in CA associated with establishment
                      woodmass_ind(i,j) = &
                           & (biomass(i,j,isapabove) + biomass(i,j,isapbelow) &
                           & + biomass(i,j,iheartabove) + biomass(i,j,iheartbelow)) &
                           & /(ind(i,j) + d_ind(i,j))

                   ENDIF
                ENDIF
             ENDDO ! Loop over # pixels - domain size

             vn(:) = un ! cannot change in static!, and veget_max implicit in d_ind

          ENDIF

          !! 4.4 total biomass of PFT added by establishment defined over veget_max ...
          total_bm_sapl(:) = zero
          DO k = 1, nparts ! Loop over # litter tissues (nparts=8); see 'stomate_constants.f90'
             WHERE(d_ind(:,j).GT.min_stomate.AND.total_bm_c(:).GT.min_stomate.AND.vn(:).GT.min_stomate)

                total_bm_sapl(:) = total_bm_sapl(:) + & 
                     bm_sapl(j,k) * d_ind(:,j) / vn(:)
             ENDWHERE
          ENDDO ! Loop over # litter tissues

          !! 4.5 Update biomass at each component
          DO k = 1, nparts ! Loop over # litter tissues

             bm_new(:) = zero

             ! first ever establishment, C flows
             WHERE( d_ind(:,j).GT.min_stomate .AND. &
                  total_bm_c(:).LE.min_stomate .AND. &
                  vn(:).GT.min_stomate)
                ! WHERE ( many_new(:) )

                !bm_new(:) = d_ind(:,j) * bm_sapl(j,k) / veget_max (:,j)
                bm_new(:) = d_ind(:,j) * bm_sapl(j,k) / vn(:)

                biomass(:,j,k) = biomass(:,j,k) + bm_new(:)

                co2_to_bm(:,j) = co2_to_bm(:,j) + bm_new(:) / dt

             ENDWHERE

             ! establishment into existing population, C flows
             WHERE(d_ind(:,j).GT.min_stomate.AND.total_bm_c(:).GT.min_stomate)

                bm_new(:) = total_bm_sapl(:) * biomass(:,j,k) / total_bm_c(:)

                biomass(:,j,k) = biomass(:,j,k) + bm_new(:)
                co2_to_bm(:,j) = co2_to_bm(:,j) + bm_new(:) / dt

             ENDWHERE
          ENDDO ! Loop over # litter tissues

          

          !! 4.6 Decrease leaf age in youngest class if new leaf biomass is higher than old one.
          WHERE ( d_ind(:,j) * bm_sapl(j,ileaf) .GT. min_stomate )
 
             ! reset leaf ages. Should do a real calculation like in the npp routine, 
             ! but this case is rare and not worth messing around.
             ! SZ 080806, added real calculation now, because otherwise leaf_age/leaf_frac
             ! are not initialised for the calculation of vmax, and hence no growth at all.
             ! logic follows that of stomate_npp.f90, just that it's been adjusted for the code here
             leaf_age(:,j,1) = leaf_age(:,j,1) * leaf_mass_young(:) / &
                  ( leaf_mass_young(:) + d_ind(:,j) * bm_sapl(j,ileaf) )

          ENDWHERE

          leaf_mass_young(:) = leaf_mass_young(:) + d_ind(:,j) * bm_sapl(j,ileaf)   

          !! 4.7 Youngest class: new mass in youngest class divided by total new mass
          WHERE ( biomass(:,j,ileaf) .GT. min_stomate )
             ! new age class fractions (fraction in youngest class increases)
             leaf_frac(:,j,1) = leaf_mass_young(:) / biomass(:,j,ileaf)

          ENDWHERE

          !! 4.8 Other classes: old mass in leaf age class divided by new mass
          DO m = 2, nleafages

             WHERE ( biomass(:,j,ileaf) .GT. min_stomate )

                leaf_frac(:,j,m) = leaf_frac(:,j,m) * & 
                     ( biomass(:,j,ileaf) + d_ind(:,j) * bm_sapl(j,ileaf) ) /  biomass(:,j,ileaf)

             ENDWHERE

          ENDDO

          !! 4.9 Update age and number of individuals
          WHERE ( d_ind(:,j) .GT. min_stomate )

             age(:,j) = age(:,j) * ind(:,j) / ( ind(:,j) + d_ind(:,j) )

             ind(:,j) = ind(:,j) + d_ind(:,j)

          ENDWHERE

          !! 4.10 Convert excess sapwood to heartwood
          IF ( tree(j) ) THEN

             sm2(:) = biomass(:,j,isapabove) + biomass(:,j,isapbelow)

             WHERE ( ( d_ind(:,j) .GT. min_stomate ) .AND. &
                  ( biomass(:,j,isapabove) + biomass(:,j,isapbelow) ) .GT. sm2(:) ) !!?? l.h.s. equals r.h.s ?

                sm_at(:) = biomass(:,j,isapabove) / &
                   ( biomass(:,j,isapabove) + biomass(:,j,isapbelow) )

                !SZ to clarify with Gerhard Krinner: This is theoretically inconsistent because 
                ! the allocation to sapwood and leaves do not follow the LPJ logic in stomate_alloc.f90
                ! hence imposing this here not only solves for the uneveness of age (mixing new and average individual)
                ! but also corrects for the discrepancy between SLAVE and LPJ logic of allocation, thus leads to excess 
                ! heartwood and thus carbon accumulation!
                ! should be removed.
                biomass(:,j,iheartabove) = biomass(:,j,iheartabove) + &
                     ( biomass(:,j,isapabove) - sm2(:) * sm_at(:) )
                biomass(:,j,isapabove) = sm2(:) * sm_at(:)

                biomass(:,j,iheartbelow) = biomass(:,j,iheartbelow) + &
                     ( biomass(:,j,isapbelow) - sm2(:) * (un - sm_at) )
                biomass(:,j,isapbelow) = sm2(:) * (un - sm_at(:))

             ENDWHERE

          ENDIF ! tree

       ENDIF ! natural

    ENDDO ! Loop over # PFTs

  !! 5. history

    d_ind = d_ind / dt

    CALL histwrite (hist_id_stomate, 'IND_ESTAB', itime, d_ind, npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'ESTABTREE', itime, estab_rate_max_tree, npts, hori_index)
    CALL histwrite (hist_id_stomate, 'ESTABGRASS', itime, estab_rate_max_grass, npts, hori_index)

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving establish'

  END SUBROUTINE establish

END MODULE lpj_establish
