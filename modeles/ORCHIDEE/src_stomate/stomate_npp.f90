! =================================================================================================================================
! MODULE          : stomate_npp
!
! CONTACT         : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE         : IPSL (2006)
!                 This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF          This modules calculates NPP: Maintenance and growth respiration
!!
!!\n DESCRIPTION: We calculate first the maintenance respiration. This is substracted from the
!!                allocatable biomass (and from the present biomass if the GPP is too low).\n
!!                Of the rest, a part is lost as growth respiration, while the other part is
!!                effectively allocated.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_npp.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_npp

  ! modules used:

  USE ioipsl
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC npp_calc,npp_calc_clear

  LOGICAL, SAVE                                              :: firstcall = .TRUE.         !! first call

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: npp_calc_clear
!!
!>\BRIEF        : Set the flag ::firstcall to .TRUE. and as such activate section 
!! 1.1 of the subroutine npp_calc (see below).\n
!_ ================================================================================================================================

  SUBROUTINE npp_calc_clear
    firstcall=.TRUE.
  END SUBROUTINE npp_calc_clear





!! ================================================================================================================================
!! SUBROUTINE	: npp_calc
!!
!>\BRIEF        Calculate NPP as the difference between GPP and respiration (= growth + maintenance respiration).
!!              Update biomass of all compartments after calculating respiration and allocation.
!!
!!
!! DESCRIPTION  : NPP is calculated from three components: Gross Primary Productivity (GPP), maintenance respiration 
!! and growth respiration (all in @tex $ gC.m^{-2}dt^{-1} $ @endtex), following the convention that positive fluxes denote 
!! fluxes plants to the atmosphere. GPP is the input variable from which, in the end, NPP or total allocatable biomass 
!! @tex $(gC.m^{-2}dt^{-1}))$ @endtex is calculated. Net primary production is then calculated as:\n	
!! NPP = GPP - growth_resp - maint-resp   [eq. 1]\n   
!!	
!! The calculation of maintenance respiration is done in routine stomate_resp.f90. Maintenance respiration is calculated for 
!! the whole plant and is therefore removed from the total allocatable biomass. In order to prevent all allocatable biomass 
!! from being used for maintenance respiration, a limit fraction of total allocatable biomass, tax_max, is defined (in 
!! variables declaration). If maintenance respiration exceeds tax_max (::bm_tax_max), the maximum allowed allocatable biomass
!! will be respired and the remaining respiration, required in excess of tax_max, is taken out from tissues already present in
!! the plant (biomass).\n  
!! 
!! After total allocatable biomass has been updated by removing maintenance respiration, total allocatable biomass is distributed 
!! to all plant compartments according to the f_alloc fractions calculated in stomate_alloc.f90.\n
!!
!! Growth respiration is calculated as a fraction of allocatable biomass for each part of the plant. The fraction coefficient 
!! ::frac_growth_resp is defined in stomate_constants.f90 and is currently set to be the same for all plant compartments. 
!! Allocatable biomass of all plant compartments are updated by removing what is lost through growth respiration. Net allocatable
!! biomass (total allocatable biomass after maintenance and growth respiration) is added to the current biomass for  each plant 
!! compartment.
!!
!! Finally, leaf age and plant age are updated. Leaf age is described with the concept of "leaf age classes". A number of leaf 
!! age classes (nleafages) is defined in stomate_constants.f90. Each leaf age class contains a fraction (::leaf_frac) of the 
!! total leaf biomass. When new biomass is added to leaves, the age of the biomass in the youngest leaf age class is decreased. 
!! The fractions of leaves in the other leaf ages classes are also updated as the total biomass has increased. Plant age is 
!! updated first by increasing the age of the previous biomass by one time step, and then by adjusting this age as the average 
!! of the ages of the previous and the new biomass.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::npp
!!
!! REFERENCE(S)	: 
!! - F.W.T.Penning De Vries, A.H.M. Brunsting, H.H. Van Laar. 1974. Products, requirements and efficiency of biosynthesis a 
!! quantitative approach. Journal of Theoretical Biology, Volume 45, Issue 2, June 1974, Pages 339-377.
!! 
!! FLOWCHART : 
!! \latexonly
!! \includegraphics[scale=0.14]{stomate_npp_flow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE npp_calc (npts, dt, &
       PFTpresent, &
       tlong_ref, t2m, tsoil, lai, rprof, &
       gpp, f_alloc, bm_alloc, resp_maint_part,&
       biomass, leaf_age, leaf_frac, age, &
       resp_maint, resp_growth, npp)
   
!! 0 Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt               !! Time step (days)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: PFTpresent       !! PFT exists (true/false)
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: tlong_ref        !! Long term annual mean 2 meter reference 
                                                                                  !! Temperature (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: t2m              !! Temperature at 2 meter (K)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)             :: tsoil            !! Soil temperature of each soil layer (K)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lai              !! PFT leaf area index (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: rprof            !! PFT root depth as calculated in stomate.f90
                                                                                  !! from root profile parameter humcste (m) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: gpp              !! PFT gross primary productivity 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)       :: f_alloc          !! Fraction of total allocatable biomass that 
                                                                                  !! goes into each plant part (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)       :: resp_maint_part  !! Maintenance respiration of different plant 
                                                                                  !! parts @tex $(gC.m^{-2}dt^{-1})$ @endtex
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: resp_maint       !! PFT maintenance respiration 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: resp_growth      !! PFT growth respiration 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex				
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: npp              !! PFT net primary productivity 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(out)      :: bm_alloc         !! PFT biomass increase, i.e. NPP per plant part 
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)    :: biomass          !! PFT total biomass of each plant part 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_age         !! PFT age of different leaf age classes (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! PFT fraction of total leaves in leaf age 
                                                                                  !! class (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! PFT age (years)

    !! 0.4 Local variables

    REAL(r_std), SAVE, DIMENSION(0:nbdl)                      :: z_soil           !! Soil levels  representing soil depth (m)				
    REAL(r_std), DIMENSION(npts,nvm)                          :: t_root           !! Root temperature (convolution of root and 
                                                                                  !! soil temperature profiles)(K)
    REAL(r_std), DIMENSION(npts,nvm,nparts)                   :: coeff_maint      !! PFT maintenance respiration coefficients of 
                                                                                  !! different plant compartments at 0 deg C 
                                                                                  !! @tex $(g.g^{-1}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nparts)                       :: t_maint          !! Temperature which is pertinent for maintenance
                                                                                  !! respiration, which is air/root temperature for
                                                                                  !! above/below-ground compartments (K)
    REAL(r_std), DIMENSION(npts)                              :: rpc              !! Scaling factor for integrating vertical soil 
                                                                                  !! profiles (unitless)
    REAL(r_std), DIMENSION(npts)                              :: tl               !! Long term annual mean temperature (C)
    REAL(r_std), DIMENSION(npts)                              :: slope            !! Slope of maintenance respiration coefficient
                                                                                  !! (1/K)
    REAL(r_std), DIMENSION(npts,nparts)                       :: resp_growth_part !! Growth respiration of different plant parts
                                                                                  !! @tex $(gC.m^{-2}dt^{-1})$ @endtex		
    REAL(r_std), DIMENSION(npts,nvm)                          :: bm_alloc_tot     !! Allocatable biomass for the whole plant
                                                                                  !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                              :: bm_add           !! Biomass increase @tex $(gC.m^{-2})$ @endtex		
    REAL(r_std), DIMENSION(npts)                              :: bm_new           !! New biomass @tex $(gC.m^{-2})$ @endtex	
    REAL(r_std), DIMENSION(npts,nvm)                          :: leaf_mass_young  !! Leaf mass in youngest age class 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex		
    REAL(r_std), DIMENSION(npts,nvm)                          :: lm_old           !! Leaf mass after maintenance respiration 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex			
    REAL(r_std), DIMENSION(npts,nvm)                          :: bm_create        !! Biomass created when biomass<0 because of dark
                                                                                  !! respiration @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                              :: bm_tax_max       !! Maximum part of allocatable biomass used for 
                                                                                  !! respiration @tex $(gC.m^{-2})$ @endtex	
    REAL(r_std), DIMENSION(npts)                              :: bm_pump          !! Biomass that remains to be taken away 
                                                                                  !! @tex $(gC.m^{-2})$ @endtex
    INTEGER(i_std)                                            :: i,j,k,l,m        !! Indeces(unitless)

!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering npp'
    
 !! 1. Initializations
    
    !! 1.1 First call
    IF ( firstcall ) THEN

       !! 1.1.1 Soil levels
       !  Get the depth of the different soil layers (number of layers=nbdl) 
       !  previously calculated as variable diaglev in routines sechiba.f90 and slowproc.f90  
       z_soil(0) = zero
       z_soil(1:nbdl) = diaglev(1:nbdl)

       !! 1.1.2 Output message
       !  Write message including value used for tax_max 	
       WRITE(numout,*) 'npp:'

       WRITE(numout,*) '   > max. fraction of allocatable biomass used for'// &
            ' maint. resp.:', tax_max

       firstcall = .FALSE.

    ENDIF ! End if first call

    !! 1.2 Set output variables to zero
    bm_alloc(:,:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    npp(:,:) = zero

!![BEGIN DISPENSABLE] between 'begin dispensable' and 'end dispensable', section can be removed. coefficients and temperatures are
!! calculated here but not used anymore since maintenance respiration temperature is now done in stomate_resp.f90
    
    !! 1.3 Root temperature
    !  Calculate root temperature as the convolution of root and temperature profiles
    !  assuming an exponential root profile.
    DO j = 2,nvm	! Loop over # PFTs

       !! 1.3.1 Calculate the integration constant ::rpc 
       !  rpc is an integration constant such that the integral of the root profile is 'one' (scaling factor)
       !  z_soil(nbdl) are the depths of the different layers
       !  rprof is the root depth. For now it is prescribed for each pft in routine stomate.f90 as the inverse of the 
       !  root profile parameter humcste (defined in constantes_veg.f90)
       rpc(:) = un / ( un - EXP( -z_soil(nbdl) / rprof(:,j) ) )

       ! 1.3.2 integrate over the # of soil layers (nbdl levels)

       t_root(:,j) = zero

       DO l = 1, nbdl	! Loop over # of soil layers 

          t_root(:,j) = &
               t_root(:,j) + tsoil(:,l) * rpc(:) * &
               ( EXP( -z_soil(l-1)/rprof(:,j) ) - EXP( -z_soil(l)/rprof(:,j) ) )

       ENDDO	! End of loop over # of soil layers

    ENDDO	! End of loop over # PFTs

    !! 1.4 Total allocatable biomass
    ! total allocatable biomass during this time step determined from GPP.
    ! GPP was calculated as CO2 assimilation in enerbil.f90
    bm_alloc_tot(:,:) = gpp(:,:) * dt

    
 !! 2. Calculate maintenance respiration coefficients  	

    DO j = 2,nvm ! Loop over # of PFTs

       
       !! 2.1 Temperatures used for calculating the maintenance respiration for the different plant parts.
       !      - for above-ground parts, we use 2-meters air temperature, (::t2m)
       !      - for below-ground parts, we use root temperature calculated in section 1.3 of this routine 
       !        as a convolution of root and temperature profiles 

       !! 2.1.1 Aboveground biomass
       t_maint(:,ileaf) = t2m(:)
       t_maint(:,isapabove) = t2m(:)
       t_maint(:,ifruit) = t2m(:)

       !! 2.1.2 Belowground biomass
       t_maint(:,isapbelow) = t_root(:,j)
       t_maint(:,iroot) = t_root(:,j)

       !! 2.1.3 Heartwood biomass
       !        Heartwood does does not respire. Any temperature could have been set [code cleaning: set t(heartwood) to undef to 'undef']
       t_maint(:,iheartbelow) = t_root(:,j)
       t_maint(:,iheartabove) = t2m(:)

       !! 2.1.4 Reserve biomass
       !        Use aboveground temperature for trees and belowground temeperature for grasses
       IF ( tree(j) ) THEN
          t_maint(:,icarbres) = t2m(:)
       ELSE
          t_maint(:,icarbres) = t_root(:,j)
       ENDIF

       
       !! 2.2 Calculate slope and coeff_maint
       !      Maintenance respiration is a fraction of biomass defined by the coefficient 
       !      coeff_maint. Coeff_maint is defined through a linear relationship with temperature
       !      which slope is the coefficient 'slope' and which intercept is 'coeff_maint_zero'.
       !      Coeff_maint_zero is defined in stomate_data.f90. 
       !      Slope is calculated here through a second-degree polynomial equation that makes it
       !     dependent on the long term temperature (to represent adaptation of the ecosystem).
       !      in stomate_resp.f90 both coefficients (slope and coeff_maint_resp) depend on:
       !      - the long term 2-m temperature (tlong_ref) converted to celsius 
       !      - the relevant temperatures used for the different plant parts as defined in section 2.1
       !      - parameters maint_resp_slope and coef_maint_zero defined in stomate_constants.f90 
       !        (cm_zero_plantpartname)
       tl(:) = tlong_ref(:) - ZeroCelsius
       slope(:) = maint_resp_slope(j,1) + tl(:) * maint_resp_slope(j,2) + &
            tl(:)*tl(:) * maint_resp_slope(j,3)

       DO k = 1, nparts

          coeff_maint(:,j,k) = &
               MAX( coeff_maint_zero(j,k) * &
               ( un + slope(:) * (t_maint(:,k)-ZeroCelsius) ), zero )

       ENDDO	! Loop over # of plant parts

    ENDDO	! Loop over # of PFTs

![END DISPENSABLE]  
 
 !! 3. Calculate maintenance and growth respiration

    ! First, total maintenance respiration for the whole plant is calculated by summing maintenance 
    ! respiration of the different plant compartments. Then, maintenance respiration is subtracted 
    ! from whole-plant allocatable biomass (up to a maximum fraction of the total allocatable biomass). 
    ! Growth respiration is then calculated for each plant compartment as a fraction of remaining 
    ! allocatable biomass for this compartment. NPP is calculated by substracting total autotrophic 
    ! respiration from GPP i.e. NPP = GPP - maintenance resp - growth resp.
    DO j = 2,nvm	! Loop over # of PFTs

![BEGIN DISPENSABLE]     
       !! 3.1 Maintenance respiration of the different plant compartments
![END DISPENSABLE]       

       !! 3.1 Maintenance respiration of the different plant parts
       !      Maintenance respiration of the different plant parts is calculated in 
       !      stomate_resp.f90 as a function of the plant's temperature, 
       !      the long term temperature and plant coefficients 
       !      VPP killer:
       ! [BEGIN DISPENSABLE]      resp_maint(:,j) = SUM( resp_maint_part(:,:), DIM=2 ) [END DISPENSABLE]
       resp_maint(:,j) = zero

       !  Following the calculation of hourly maintenance respiration, verify that 
       !  the PFT has not been killed after calcul of resp_maint_part in stomate.
       DO k= 1, nparts
          WHERE (PFTpresent(:,j))
             resp_maint(:,j) = resp_maint(:,j) + resp_maint_part(:,j,k)
          ENDWHERE
       ENDDO
       
       !! 3.2 Substract maintenance respiration from allocatable biomass
       !      The total maintenance respiration calculated in 3.2 is substracted  from the newly 
       !      produced allocatable biomass (bm_alloc_tot). However, ensure that not all allocatable 
       !      biomass is removed by setting a maximum to the fraction of allocatable biomass used 
       !      for maintenance respiration: tax_max. If the maintenance respiration is larger than 
       !      tax_max,the amount tax_max is taken from allocatable biomass, and the remaining of 
       !      maintenance respiration is taken from the tissues themselves (biomass). We suppose 
       !      that respiration is not dependent on leaf age -> therefore the leaf age structure is 
       !      not changed. 
       !      The maximum fraction of allocatable biomass used for respiration is defined as tax_max. 
       !      The value of tax_max is set in the declarations section (0.4 Local variables) of this
       !      routine
       bm_tax_max(:) = tax_max * bm_alloc_tot(:,j)

       DO i = 1, npts	! Loop over # of pixels

          ! If there is enough allocatable biomass to cover maintenance respiration, 
	  ! then biomass associated with maintenance respiration is removed from allocatable biomass
          IF ( ( bm_alloc_tot(i,j) .GT. zero ) .AND. &
               ( ( resp_maint(i,j) * dt ) .LT. bm_tax_max(i) ) ) THEN	
	
             bm_alloc_tot(i,j) = bm_alloc_tot(i,j) - resp_maint(i,j) * dt

             ! [BEGIN DISPENSABLE] Same as current line except that zero became min_stomate 
             ! Shilong        ELSEIF ( resp_maint(i,j) .GT. zero ) THEN [END DISPENSABLE] 

          ! If there is not enough allocatable biomass to cover maintenance respiration, the  
          ! - maximum allowed allocatable biomass (bm_tax_max) is removed from allocatable biomass.
          ELSEIF ( resp_maint(i,j) .GT. min_stomate ) THEN
	     
             bm_alloc_tot(i,j) = bm_alloc_tot(i,j) - bm_tax_max(i)

             ! ::bm_pump is the amount of maintenance respiration that exceeds the maximum allocatable biomass
	     ! This amount of biomass still needs to be respired and will be removed from tissues biomass of each
	     ! plant compartment
             bm_pump(i) = resp_maint(i,j) * dt - bm_tax_max(i)

             ! The biomass is removed from each plant compartment tissues as the ratio of the maintenance		
             ! respiration of the plant compartment to the total maintenance respiration (resp_maint_part/resp_maint)
             biomass(i,j,ileaf) = biomass(i,j,ileaf) - &
                  bm_pump(i) * resp_maint_part(i,j,ileaf) / resp_maint(i,j)
             biomass(i,j,isapabove) = biomass(i,j,isapabove) - &
                  bm_pump(i) * resp_maint_part(i,j,isapabove) / resp_maint(i,j)
             biomass(i,j,isapbelow) = biomass(i,j,isapbelow) - &
                  bm_pump(i) * resp_maint_part(i,j,isapbelow) / resp_maint(i,j)
             biomass(i,j,iroot) = biomass(i,j,iroot) - &
                  bm_pump(i) * resp_maint_part(i,j,iroot) / resp_maint(i,j)
             biomass(i,j,ifruit) = biomass(i,j,ifruit) - &
                  bm_pump(i) * resp_maint_part(i,j,ifruit) / resp_maint(i,j)
             biomass(i,j,icarbres) = biomass(i,j,icarbres) - &
                  bm_pump(i) * resp_maint_part(i,j,icarbres) / resp_maint(i,j)

          ENDIF	! End if there is enough allocatable biomass to cover maintenance respiration

       ENDDO   ! Fortran95: WHERE - ELSEWHERE construct

       
       !! 3.3 Allocate allocatable biomass to different plant compartments.
       !      The amount of allocatable biomass of each compartment is a fraction according f_alloc of total 
       !      allocatable biomass (the f_alloc of the different plant parts are calculated in stomate_alloc.f90)
       DO k = 1, nparts
          bm_alloc(:,j,k) = f_alloc(:,j,k) * bm_alloc_tot(:,j)
       ENDDO

       
       !! 3.4 Calculate growth respiration of each plant compartment.
       !      Growth respiration of a plant compartment is a fraction of the allocatable biomass remaining after
       !      maintenance respiration losses have been taken into account. The fraction of allocatable biomass 
       !      removed for growth respiration is the same for all plant compartments and is defined by the parameter
       !      frac_growth_resp in stomate_constants.f90. Allocatable biomass ::bm_alloc is updated as a result of 
       !      the removal of growth resp.
       resp_growth_part(:,:) = frac_growthresp * bm_alloc(:,j,:) / dt
       bm_alloc(:,j,:) = ( un - frac_growthresp ) * bm_alloc(:,j,:)

       
       !! 3.5 Total growth respiration 
       !      Calculate total growth respiration of the plant as the sum of growth respiration of all plant parts	
       resp_growth(:,j) = zero

       DO k = 1, nparts
          resp_growth(:,j) = resp_growth(:,j) + resp_growth_part(:,k)
       ENDDO

    ENDDO ! # End Loop over # of PFTs

    
 !! 4. Update the biomass with newly allocated biomass after respiration
 
    !  Save the old leaf biomass for later. "old" leaf mass is leaf mass after maintenance respiration in the case 
    !  where maintenance respiration has required taking biomass from tissues in section 3.3 
    lm_old(:,:) = biomass(:,:,ileaf)
    biomass(:,:,:) = biomass(:,:,:) + bm_alloc(:,:,:)

    
 !! 5. Deal with negative biomasses
    
    !  Biomass can become negative in some rare cases, as the GPP can be negative. This corresponds to very 
    !  situations that can be seen as the 'creation' of a seed ('virtual photosynthesis'). In this case, we set
    !  biomass to some small value min_stomate. For carbon budget to remain balanced, this creation of matter (carbon) 
    !  is taken into account by decreasing the autotrophic respiration by the same amount that has been added to biomass 
    !  for it to become positive. In this case, maintenance respiration can become negative in extreme cases (deserts�)!!

    DO k = 1, nparts	! Loop over # of plant parts

       DO j = 2,nvm	! Loop over # of PFTs

          WHERE ( biomass(:,j,k) .LT. zero )

             bm_create(:,j) = min_stomate - biomass(:,j,k)
             
             biomass(:,j,k) = biomass(:,j,k) + bm_create(:,j)
             
             resp_maint(:,j) = resp_maint(:,j) - bm_create(:,j) / dt

          ENDWHERE

       ENDDO	! Loop over # of PFTs

    ENDDO	! Loop over # plant parts

    
 !! 6. Calculate NPP (See Eq 1 in header)
    
    !  Calculate the NPP @tex $(gC.m^{-2}dt^{-1})$ @endtex as the difference between GPP
    !  and autotrophic respiration (maintenance and growth respirations)
    DO j = 2,nvm	! Loop over # PFTs
       npp(:,j) = gpp(:,j) - resp_growth(:,j) - resp_maint(:,j)
    ENDDO	! Loop over # PFTs

    
 !! 7. Update leaf age

    !  Leaf age is needed for calculation of turnover and vmax in stomate_turnover.f90 and stomate_vmax.f90 routines. 
    !  Leaf biomass is distributed according to its age into several "age classes" with age class=1 representing the
    !  youngest class, and consisting of the most newly allocated leaf biomass 
    
    !! 7.1 Update quantity and age of the leaf biomass in the youngest class
    !      The new amount of leaf biomass in the youngest age class (leaf_mass_young) is the sum of :
    !      - the leaf biomass that was already in the youngest age class (leaf_frac(:,j,1) * lm_old(:,j)) with the 
    !        leaf age given in leaf_age(:,j,1) 
    !      - and the new biomass allocated to leaves (bm_alloc(:,j,ileaf)) with a leaf age of zero.
    DO j = 2,nvm
       leaf_mass_young(:,j) = leaf_frac(:,j,1) * lm_old(:,j) + bm_alloc(:,j,ileaf)
    ENDDO

    ! The age of the updated youngest age class is the average of the ages of its 2 components: bm_alloc(leaf) of age
    ! '0', and leaf_frac*lm_old(=leaf_mass_young-bm_alloc) of age 'leaf_age(:,j,1)' 
    DO j = 2,nvm
       WHERE ( ( bm_alloc(:,j,ileaf) .GT. zero ) .AND. &
         ( leaf_mass_young(:,j) .GT. zero ) )

          leaf_age(:,j,1) = MAX ( zero, &
               & leaf_age(:,j,1) * &
               & ( leaf_mass_young(:,j) - bm_alloc(:,j,ileaf) ) / &
               & leaf_mass_young(:,j) )
          
       ENDWHERE
    ENDDO

    !! 7.2 Update leaf age
    !      Update fractions of leaf biomass in each age class (fraction in youngest class increases)

    !! 7.2.1 Update age of youngest leaves
    !        For age class 1 (youngest class), because we have added biomass to the youngest class, we need to update
    !        the fraction of total leaf biomass that belongs to the youngest age class : updated mass in class divided
    !        by new total leaf mass
    DO j = 2,nvm
       WHERE ( biomass(:,j,ileaf) .GT. min_stomate )

          leaf_frac(:,j,1) = leaf_mass_young(:,j) / biomass(:,j,ileaf)

       ENDWHERE
    ENDDO

    !! 7.2.2 Update age of other age classes
    !        Because the total leaf biomass has changed, we need to update the fraction of leaves in each age class:
    !        mass in leaf age class (from previous fraction of leaves in this class and previous total leaf biomass) 
    !        divided by new total mass
    DO m = 2, nleafages	! Loop over # leaf age classes

       DO j = 2,nvm	! Loop over # PFTs
          WHERE ( biomass(:,j,ileaf) .GT. min_stomate )

             leaf_frac(:,j,m) = leaf_frac(:,j,m) * lm_old(:,j) / biomass(:,j,ileaf)

          ENDWHERE
       ENDDO

    ENDDO	! Loop over # leaf age classes

 !! 8. Update whole-plant age 
    
    !! 8.1 PFT age
    !      At every time step, increase age of the biomass that was already present at previous time step. 
    !      Age is expressed in years, and the time step 'dt' in days so age increase is: dt divided by number 
    !      of days in a year.
    WHERE ( PFTpresent(:,:) )

       age(:,:) = age(:,:) + dt/one_year

    ELSEWHERE

       age(:,:) = zero

    ENDWHERE

    !! 8.2 Age of grasses and crops
    !  For grasses and crops, biomass with age 0 has been added to the whole plant with age 'age'. New biomass is the sum of 
    !  the current total biomass in all plant parts (bm_new), bm_new(:) = SUM( biomass(:,j,:), DIM=2 ). The biomass that has 
    !  just been added is the sum of the allocatable biomass of all plant parts (bm_add), its age is zero. bm_add(:) = 
    !  SUM( bm_alloc(:,j,:), DIM=2 ). Before allocation, the plant biomass is bm_new-bm_add, its age is "age(:,j)". The age of
    !  the new biomass is the average of the ages of previous and added biomass.
    !  For trees, age is treated in "establish" if vegetation is dynamic, and in turnover routines if it is static (in this 
    !  case, only the age of the heartwood is accounted for).
    DO j = 2,nvm

       IF ( .NOT. tree(j) ) THEN

          bm_new(:) = biomass(:,j,ileaf) + biomass(:,j,isapabove) + &
               biomass(:,j,iroot) + biomass(:,j,ifruit)
          bm_add(:) = bm_alloc(:,j,ileaf) + bm_alloc(:,j,isapabove) + &
               bm_alloc(:,j,iroot) + bm_alloc(:,j,ifruit)

          WHERE ( ( bm_new(:) .GT. zero ) .AND. ( bm_add(:) .GT. zero ) )
             age(:,j) = age(:,j) * ( bm_new(:) - bm_add(:) ) / bm_new(:)
          ENDWHERE

       ENDIF

    ENDDO

 !! 9. Write history files

    ! Save in history file the variables describing the biomass allocated to the plant parts
    CALL histwrite (hist_id_stomate, 'BM_ALLOC_LEAF', itime, &
         bm_alloc(:,:,ileaf), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'BM_ALLOC_SAP_AB', itime, &
         bm_alloc(:,:,isapabove), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'BM_ALLOC_SAP_BE', itime, &
         bm_alloc(:,:,isapbelow), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'BM_ALLOC_ROOT', itime, &
         bm_alloc(:,:,iroot), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'BM_ALLOC_FRUIT', itime, &
         bm_alloc(:,:,ifruit), npts*nvm, horipft_index)
    CALL histwrite (hist_id_stomate, 'BM_ALLOC_RES', itime, &
         bm_alloc(:,:,icarbres), npts*nvm, horipft_index)

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving npp'

  END SUBROUTINE npp_calc

END MODULE stomate_npp
