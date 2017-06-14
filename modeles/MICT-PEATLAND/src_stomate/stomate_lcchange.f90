! =================================================================================================================================
! MODULE       : stomate_lcchange
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Impact of land cover change on carbon stocks
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_lcchange.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================


MODULE stomate_lcchange

  ! modules used:
  
  USE ioipsl
  USE stomate_data
  USE pft_parameters
  USE constantes
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC lcchange_main
  
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : lcchange_main
!!
!>\BRIEF        Impact of land cover change on carbon stocks
!!
!! DESCRIPTION  : This subroutine is activated by 'LAND_COVER_CHANGE = y' in configuration file.
!! In the case of 'LAND_COVER_CHANGE = y', the impact of land cover change on 
!! carbon stocks is computed in this subroutine. The land cover change is written
!! by the difference of new and previous "maximal" coverage fraction of a PFT. 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',
!! and increase/decrease of each component, are estimated.\n
!!
!! Main structure of lpj_establish.f90 is:
!! 1. Initialization
!! 2. Calculation of changes in carbon stocks and biomass by land cover change
!! 3. Update 10 year- and 100 year-turnover pool contents
!! 4. History
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::prod10, ::prod100, ::flux10, ::flux100,
!!   :: cflux_prod10 and :: cflux_prod100 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \latexonly 
!!     \includegraphics[scale=0.5]{lcchange.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  
  SUBROUTINE lcchange_main ( npts, dt_days, veget_max, veget_max_new,&
       biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &        
       co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, tree, cn_ind,flux10,flux100, &
       prod10,prod100,&
       convflux,&
       cflux_prod10,cflux_prod100, leaf_frac,&
       npp_longterm, lm_lastyearmax, litter, carbon)
    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt_days          !! Time step of vegetation dynamics for stomate
                                                                                  !! (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(INOUT)           :: veget_max_new    !! new "maximal" coverage fraction of a PFT (LAI
                                                                                  !! -> infinity) on ground
    REAL(r_std) , DIMENSION (nvm, nparts), INTENT(in)         :: bm_sapl          !! biomass of sapling 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    LOGICAL, DIMENSION(nvm), INTENT(in)                       :: tree             !! Does the pft contains trees?

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts), INTENT(out)                 :: convflux         !! release during first year following land cover
                                                                                  !! change
    REAL(r_std), DIMENSION(npts), INTENT(out)                 :: cflux_prod10     !! total annual release from the 10 year-turnover
                                                                                  !! pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts), INTENT(out)                 :: cflux_prod100    !! total annual release from the 100 year-
                                                                                  !! turnover pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out)      :: turnover_daily   !! Turnover rates 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex

    !! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: veget_max        !! "maximal" coverage fraction of a PFT (LAI ->
                                                                                  !! infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)    :: biomass          !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind              !! Number of individuals @tex ($m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
                                                                                  !! growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm        !! biomass uptaken 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)    :: bm_to_litter     !! conversion of biomass to litter 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind           !! crown area of individuals 
                                                                                  !! @tex ($m^{2}$) @endtex
    REAL(r_std), DIMENSION(npts,0:10), INTENT(inout)          :: prod10           !! products remaining in the 10 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (10 + 1 : input from year of land
                                                                                  !! cover change)
    REAL(r_std), DIMENSION(npts,0:100), INTENT(inout)         :: prod100          !! products remaining in the 100 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (100 + 1 : input from year of land
                                                                                  !! cover change)
    REAL(r_std), DIMENSION(npts,10), INTENT(inout)            :: flux10           !! annual release from the 10/100 year-turnover 
                                                                                  !! pool compartments
    REAL(r_std), DIMENSION(npts,100), INTENT(inout)           :: flux100          !! annual release from the 10/100 year-turnover
                                                                                  !! pool compartments
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs), INTENT(inout):: litter           !! metabolic and structural litter, above and 
                                                                                  !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon           !! carbon pool: active, slow, or passive 
                                                                                  !! @tex ($gC m^{-2}$) @endtex

    !! 0.4 Local variables

    INTEGER(i_std)                                            :: i, j, k, l, m    !! indices (unitless)
    REAL(r_std)                                               :: bm_new           !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nparts)                        :: biomass_loss     !! biomass loss @tex ($gC m^{-2}$) @endtex
    REAL(r_std)                                               :: above            !! aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nlevs)                   :: dilu_lit         !! Litter dilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ncarb)                         :: dilu_soil_carbon !! Soil Carbondilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(nvm)                                :: delta_veg        !! changes in "maximal" coverage fraction of PFT 
    REAL(r_std)                                               :: delta_veg_sum    !! sum of delta_veg
    REAL(r_std),DIMENSION(npts,nvm)                           :: delta_ind        !! change in number of individuals  
!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering lcchange_main'
    
  !! 1. initialization
    
    prod10(:,0)         = zero
    prod100(:,0)        = zero   
    above               = zero
    convflux(:)         = zero
    cflux_prod10(:)     = zero
    cflux_prod100(:)    = zero
    delta_ind(:,:)      = zero
    delta_veg(:)        = zero
    
  !! 2. calculation of changes in carbon stocks and biomass by land cover change\n
    
    DO i = 1, npts ! Loop over # pixels - domain size
       
       !! 2.1 initialization of carbon stocks\n
       delta_veg(:) = veget_max_new(i,:)-veget_max(i,:)
       delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.0.)
       
       dilu_lit(i,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:) = zero
       
       !! 2.2 if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       DO j=2, nvm
          IF ( delta_veg(j) < -min_stomate ) THEN 
             dilu_lit(i,:,:) = dilu_lit(i,:,:) + delta_veg(j)*litter(i,:,j,:) / delta_veg_sum
             dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
             biomass_loss(i,:) = biomass_loss(i,:) + biomass(i,j,:)*delta_veg(j) / delta_veg_sum
          ENDIF
       ENDDO
       
       !! 2.3 
       DO j=2, nvm ! Loop over # PFTs

          !! 2.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN

             !! 2.3.1.1 Initial setting of new establishment
             IF (veget_max(i,j) .LT. min_stomate) THEN 
                IF (tree(j)) THEN

                   ! cn_sapl(j)=0.5; stomate_data.f90
                   cn_ind(i,j) = cn_sapl(j) 
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

                ! large_value = 1.E33_r_std
                when_growthinit(i,j) = large_value 
                leaf_frac(i,j,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf) * ind(i,j)
             ENDIF
             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j) 
             ENDIF
             
             !! 2.3.1.2 Update of biomass in each each carbon stock component 
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n
             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
                bm_new = delta_ind(i,j) * bm_sapl(j,k) 
                IF (veget_max(i,j) .GT. min_stomate) THEN

                   ! in the case that bm_new is overestimated compared with biomass?
                   IF ((bm_new/delta_veg(j)) > biomass(i,j,k)) THEN
                      bm_new = biomass(i,j,k)*delta_veg(j)
                   ENDIF
                ENDIF
                biomass(i,j,k) = ( biomass(i,j,k) * veget_max(i,j) + bm_new )  / veget_max_new(i,j)
                co2_to_bm(i,j) = co2_to_bm(i,j)+  (bm_new* dt_days) / (one_year * veget_max_new(i,j))
             ENDDO ! loop over # carbon stock components
             
             !! 2.3.1.3 Calculation of dilution in litter, soil carbon, and  input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?

             ! Litter
             litter(i,:,j,:)=(litter(i,:,j,:) * veget_max(i,j) + &
                  dilu_lit(i,:,:) * delta_veg(j)) / veget_max_new(i,j)
            
             ! Soil carbon
             carbon(i,:,j)=(carbon(i,:,j) * veget_max(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_max_new(i,j)

             ! Litter input
             bm_to_litter(i,j,isapbelow) = bm_to_litter(i,j,isapbelow) + &
                  &                         biomass_loss(i,isapbelow)*delta_veg(j) / veget_max_new(i,j)
             bm_to_litter(i,j,iheartbelow) = bm_to_litter(i,j,iheartbelow) + biomass_loss(i,iheartbelow) *delta_veg(j) &
                  &  / veget_max_new(i,j)
             bm_to_litter(i,j,iroot) = bm_to_litter(i,j,iroot) + biomass_loss(i,iroot)*delta_veg(j) / veget_max_new(i,j)
             bm_to_litter(i,j,ifruit) = bm_to_litter(i,j,ifruit) + biomass_loss(i,ifruit)*delta_veg(j) / veget_max_new(i,j)
             bm_to_litter(i,j,icarbres) = bm_to_litter(i,j,icarbres) + &
                  &                         biomass_loss(i,icarbres)   *delta_veg(j) / veget_max_new(i,j)
             bm_to_litter(i,j,ileaf) = bm_to_litter(i,j,ileaf) + biomass_loss(i,ileaf)*delta_veg(j) / veget_max_new(i,j)
             age(i,j)=age(i,j)*veget_max(i,j)/veget_max_new(i,j)
             
          !! 2.3.2 The case that vegetation coverage of PFTj is no change or decreases
          ELSE 
 
             !! 2.3.2.1 Biomass export
             ! coeff_lcchange_*:  Coeff of biomass export for the year, decade, and century
             above = biomass(i,j,isapabove) + biomass(i,j,iheartabove)
             convflux(i)  = convflux(i)  - ( coeff_lcchange_1(j) * above * delta_veg(j) ) 
             prod10(i,0)  = prod10(i,0)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0) = prod100(i,0) - ( coeff_lcchange_100(j) * above * delta_veg(j) )
             
             !! 2.3.2.2 Total reduction
             IF ( veget_max_new(i,j) .LT. min_stomate ) THEN 
                
                veget_max_new(i,j)= zero
                ind(i,j) = zero
                biomass(i,j,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                carbon(i,:,j) = zero
                litter(i,:,j,:) = zero
                bm_to_litter(i,j,:) = zero
                turnover_daily(i,j,:) = zero
                
             ENDIF
             
          ENDIF ! End if PFT's coverage reduction
          
       ENDDO ! Loop over # PFTs
       
       !! 2.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i) =  cflux_prod10(i) + flux10(i,m)
          prod10(i,m)     =  prod10(i,m-1)   - flux10(i,m-1)
          flux10(i,m)     =  flux10(i,m-1)
          
          IF (prod10(i,m) .LT. 1.0) prod10(i,m) = zero
       ENDDO
       
       cflux_prod10(i) = cflux_prod10(i) + flux10(i,1) 
       flux10(i,1)     = 0.1 * prod10(i,0)
       prod10(i,1)     = prod10(i,0)
       
       !! 2.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i)  =  cflux_prod100(i) + flux100(i,m)
          prod100(i,m)      =  prod100(i,m-1)   - flux100(i,m-1)
          flux100(i,m)      =  flux100(i,m-1)
          
          IF (prod100(i,m).LT.1.0) prod100(i,m) = zero
       ENDDO
       
       cflux_prod100(i)  = cflux_prod100(i) + flux100(i,1) 
       flux100(i,1)      = 0.01 * prod100(i,0)
       prod100(i,1)      = prod100(i,0)
       prod10(i,0)        = zero
       prod100(i,0)       = zero 
       
    ENDDO ! Loop over # pixels - domain size
    
  !! 3. history
    
    veget_max(:,:) = veget_max_new(:,:)
    convflux        = convflux/one_year*dt_days
    cflux_prod10    = cflux_prod10/one_year*dt_days
    cflux_prod100   = cflux_prod100/one_year*dt_days
    
    IF (bavard.GE.4) WRITE(numout,*) 'Leaving lcchange_main'
    
  END SUBROUTINE lcchange_main
  
END MODULE stomate_lcchange
