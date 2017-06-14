! =================================================================================================================================
! MODULE       : stomate_soilcarbon
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Calculate soil dynamics largely following the Century model
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_soilcarbon.f90 $ 
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_soilcarbon

  ! modules used:

  USE ioipsl
  USE stomate_data
  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC soilcarbon,soilcarbon_clear

  ! Variables shared by all subroutines in this module
  
  LOGICAL, SAVE    :: firstcall = .TRUE.   !! Is this the first call? (true/false)

CONTAINS


!! ================================================================================================================================
!!  SUBROUTINE   : soilcarbon_clear
!!
!>\BRIEF        Set the flag ::firstcall to .TRUE. and as such activate sections 1.1.2 and 1.2 of the subroutine soilcarbon 
!! (see below).
!! 
!_ ================================================================================================================================
  
  SUBROUTINE soilcarbon_clear
    firstcall=.TRUE.
  ENDSUBROUTINE soilcarbon_clear


!! ================================================================================================================================
!!  SUBROUTINE   : soilcarbon
!!
!>\BRIEF        Computes the soil respiration and carbon stocks, essentially 
!! following Parton et al. (1987).
!!
!! DESCRIPTION	: The soil is divided into 3 carbon pools, with different 
!! characteristic turnover times : active (1-5 years), slow (20-40 years) 
!! and passive (200-1500 years).\n
!! There are three types of carbon transferred in the soil:\n
!! - carbon input in active and slow pools from litter decomposition,\n
!! - carbon fluxes between the three pools,\n
!! - carbon losses from the pools to the atmosphere, i.e., soil respiration.\n
!!
!! The subroutine performs the following tasks:\n
!!
!! Section 1.\n
!! The flux fractions (f) between carbon pools are defined based on Parton et 
!! al. (1987). The fractions are constants, except for the flux fraction from
!! the active pool to the slow pool, which depends on the clay content,\n
!! \latexonly
!! \input{soilcarbon_eq1.tex}
!! \endlatexonly\n
!! In addition, to each pool is assigned a constant turnover time.\n
!!
!! Section 2.\n
!! The carbon input, calculated in the stomate_litter module, is added to the 
!! carbon stock of the different pools.\n
!!
!! Section 3.\n
!! First, the outgoing carbon flux of each pool is calculated. It is 
!! proportional to the product of the carbon stock and the ratio between the 
!! iteration time step and the residence time:\n
!! \latexonly
!! \input{soilcarbon_eq2.tex}
!! \endlatexonly
!! ,\n
!! Note that in the case of crops, the additional multiplicative factor 
!! integrates the faster decomposition due to tillage (following Gervois et 
!! al. (2008)).
!! In addition, the flux from the active pool depends on the clay content:\n
!! \latexonly
!! \input{soilcarbon_eq3.tex}
!! \endlatexonly
!! ,\n
!! Each pool is then cut from the carbon amount corresponding to each outgoing
!! flux:\n
!! \latexonly
!! \input{soilcarbon_eq4.tex}
!! \endlatexonly\n
!! Second, the flux fractions lost to the atmosphere is calculated in each pool
!! by subtracting from 1 the pool-to-pool flux fractions. The soil respiration 
!! is then the summed contribution of all the pools,\n
!! \latexonly
!! \input{soilcarbon_eq5.tex}
!! \endlatexonly\n
!! Finally, each carbon pool accumulates the contribution of the other pools:
!! \latexonly
!! \input{soilcarbon_eq6.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUTS VARIABLE(S): carbon, resp_hetero_soil
!!
!! REFERENCE(S)   :
!! - Parton, W.J., D.S. Schimel, C.V. Cole, and D.S. Ojima. 1987. Analysis of 
!! factors controlling soil organic matter levels in Great Plains grasslands. 
!! Soil Sci. Soc. Am. J., 51, 1173-1179.
!! - Gervois, S., P. Ciais, N. de Noblet-Ducoudre, N. Brisson, N. Vuichard, 
!! and N. Viovy (2008), Carbon and water balance of European croplands 
!! throughout the 20th century, Global Biogeochem. Cycles, 22, GB2022, 
!! doi:10.1029/2007GB003018.
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics[scale=0.5]{soilcarbon_flowchart.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE soilcarbon (npts, dt, clay, &
       soilcarbon_input, control_temp, control_moist, &
       carbon, resp_hetero_soil)


!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                            :: npts             !! Domain size (unitless)
    REAL(r_std), INTENT(in)                               :: dt               !! Time step \f$(dtradia one_day^{-1})$\f !Chloe : Non c'est faux c'est one_day/dtradia !
    REAL(r_std), DIMENSION(npts), INTENT(in)              :: clay             !! Clay fraction (unitless, 0-1) 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(in)    :: soilcarbon_input !! Amount of carbon going into the carbon pools from litter decomposition \f$(gC m^{-2} day^{-1})$\f
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(in)        :: control_temp     !! Temperature control of heterotrophic respiration (unitless: 0->1)
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(in)        :: control_moist    !! Moisture control of heterotrophic respiration (unitless: 0.25->1)
! Chloe :
!    REAL(r_std), DIMENSION(npts), INTENT(in)              :: wpro_conso_carb  !! Amount of Carbon consumed by methanogenese (Peatland emissions)
    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)         :: resp_hetero_soil !! Soil heterotrophic respiration \f$(gC m^{-2} (dtradia one_day^{-1})^{-1})$\f

    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout) :: carbon           !! Soil carbon pools: active, slow, or passive, \f$(gC m^{2})$\f

    !! 0.4 Local variables
    REAL(r_std), SAVE, DIMENSION(ncarb)                   :: carbon_tau       !! Residence time in carbon pools (days)
    REAL(r_std), DIMENSION(npts,ncarb,ncarb)              :: frac_carb        !! Flux fractions between carbon pools 
                                                                              !! (second index=origin, third index=destination) 
                                                                              !! (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,ncarb)                    :: frac_resp        !! Flux fractions from carbon pools to the atmosphere (respiration) (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,ncarb)                    :: fluxtot          !! Total flux out of carbon pools \f$(gC m^{2})$\f
    REAL(r_std), DIMENSION(npts,ncarb,ncarb)              :: flux             !! Fluxes between carbon pools \f$(gC m^{2})$\f
    CHARACTER(LEN=7), DIMENSION(ncarb)                    :: carbon_str       !! Name of the carbon pools for informative outputs (unitless)
    INTEGER(i_std)                                        :: k,kk,m,i           !! Indices (unitless)

!_ ================================================================================================================================

    !! bavard is the level of diagnostic information, 0 (none) to 4 (full)
    IF (bavard.GE.3) WRITE(numout,*) 'Entering soilcarbon' 

!! 1. Initializations
    
    !! 1.1 Get soil "constants"
    !! 1.1.1 Flux fractions between carbon pools: depend on clay content, recalculated each time
    ! From active pool: depends on clay content
    frac_carb(:,iactive,iactive) = zero
    frac_carb(:,iactive,ipassive) = frac_carb_ap
    frac_carb(:,iactive,islow) = un - (metabolic_ref_frac - active_to_pass_clay_frac*clay(:)) - frac_carb(:,iactive,ipassive)

    ! 1.1.1.2 from slow pool

    frac_carb(:,islow,islow) = zero
    frac_carb(:,islow,iactive) = frac_carb_sa
    frac_carb(:,islow,ipassive) = frac_carb_sp

    ! From passive pool
    frac_carb(:,ipassive,ipassive) = zero
    frac_carb(:,ipassive,iactive) = frac_carb_pa
    frac_carb(:,ipassive,islow) = frac_carb_ps

    IF ( firstcall ) THEN

        !! 1.1.2 Residence times in carbon pools (days)
        carbon_tau(iactive) = carbon_tau_iactive * one_year       ! 1.5 years. This is same as CENTURY. But, in Parton et al. (1987), it's weighted by moisture and temperature dependences.
        carbon_tau(islow) = carbon_tau_islow * one_year          ! 25 years. This is same as CENTURY. But, in Parton et al. (1987), it's weighted by moisture and temperature dependences.
        carbon_tau(ipassive) = carbon_tau_ipassive * one_year       ! 1000 years. This is same as CENTURY. But, in Parton et al. (1987), it's weighted by moisture and temperature dependences.
        
        !! 1.2 Messages : display the residence times  
        carbon_str(iactive) = 'active'
        carbon_str(islow) = 'slow'
        carbon_str(ipassive) = 'passive'
        
        WRITE(numout,*) 'soilcarbon:'
        
        WRITE(numout,*) '   > minimal carbon residence time in carbon pools (d):'
        DO k = 1, ncarb ! Loop over carbon pools
          WRITE(numout,*) '(1, ::carbon_str(k)):',carbon_str(k),' : (1, ::carbon_tau(k)):',carbon_tau(k)
        ENDDO
        
        WRITE(numout,*) '   > flux fractions between carbon pools: depend on clay content'
        
        firstcall = .FALSE.
        
    ENDIF

    !! 1.3 Set soil respiration to zero
    resp_hetero_soil(:,:) = zero

!! 2. Update the carbon stocks with the different soil carbon input
    carbon(:,:,:) = carbon(:,:,:) + soilcarbon_input(:,:,:) * dt

!! 3. Fluxes between carbon reservoirs, and to the atmosphere (respiration) \n

    !! 3.1. Determine the respiration fraction : what's left after
    ! subtracting all the 'pool-to-pool' flux fractions
    ! Diagonal elements of frac_carb are zero
    !    VPP killer:
    !     frac_resp(:,:) = 1. - SUM( frac_carb(:,:,:), DIM=3 )
    frac_resp(:,:) = un - frac_carb(:,:,iactive) - frac_carb(:,:,islow) - &
         frac_carb(:,:,ipassive) 

    !! 3.2. Calculate fluxes

    DO m = 2,nvm ! Loop over # PFTs

      !! 3.2.1. Flux out of pools

      DO k = 1, ncarb ! Loop over carbon pools from which the flux comes
        
        ! Determine total flux out of pool
        ! S.L. Piao 2006/05/05 - for crop multiply tillage factor of decomposition
        ! Not crop
!! Chloe++
!! On n'applique pas si peat :
        ! IF ( natural(m) ) THEN
          IF (natural(m) .AND. .NOT. m .EQ. 14) THEN
!Chloe--
            fluxtot(:,k) = dt/carbon_tau(k) * carbon(:,k,m) * &
                  control_moist(:,ibelow) * control_temp(:,ibelow)
         ! C3 crop
          ELSEIF ( (.NOT. natural(m)) .AND. (.NOT. is_c4(m)) ) THEN
             fluxtot(:,k) = dt/carbon_tau(k) * carbon(:,k,m) * &
                  control_moist(:,ibelow) * control_temp(:,ibelow) * flux_tot_coeff(1)
          ! C4 Crop   
          ELSEIF ( (.NOT. natural(m)) .AND. is_c4(m) ) THEN
             fluxtot(:,k) = dt/carbon_tau(k) * carbon(:,k,m) * &
                  control_moist(:,ibelow) * control_temp(:,ibelow) * flux_tot_coeff(2)
          ENDIF
        ! END - S.L. Piao 2006/05/05 - for crop multiply tillage factor of decomposition
        
!! Chloe++ Decomp Limitee
        !! On applique un facteur de decomposition control_moist=0.35 (Cf Wania et al 2009)
        !! pour le cas du peat :
        IF (m .EQ. 14) THEN
             fluxtot(:,k) = dt/carbon_tau(k) * carbon(:,k,m) * 0.35 * control_temp(:,ibelow)            
   !          fluxtot(:,k) = dt/carbon_tau(k) * carbon(:,k,m) * control_moist(:,ibelow) * control_temp(:,ibelow)
        ENDIF
!!Chloe--


        ! Carbon flux from active pools depends on clay content
        IF ( k .EQ. iactive ) THEN
            fluxtot(:,k) = fluxtot(:,k) * ( un - flux_tot_coeff(3) * clay(:) )
        ENDIF

        ! Update the loss in each carbon pool
        carbon(:,k,m) = carbon(:,k,m) - fluxtot(:,k)


        !Chloe :
        !IF ( k .EQ. iactive .AND. m .EQ. 14) THEN 
        !write(*,*) 'Chloe AV conso',carbon(:,k,m)
        !carbon(:,k,m) = carbon(:,k,m) - wpro_conso_carb(:)*dt
        !write(*,*) 'Chloe AP conso',carbon(:,k,m), wpro_conso_carb(:)*dt
        !ENDIF
        !Chloe --

        ! Fluxes towards the other pools (k -> kk)
        DO kk = 1, ncarb ! Loop over the carbon pools where the flux goes
          flux(:,k,kk) = frac_carb(:,k,kk) * fluxtot(:,k)
        ENDDO
        
      ENDDO ! End of loop over carbon pools
      
      !! 3.2.2 respiration
      !BE CAREFUL: Here resp_hetero_soil is divided by dt to have a value which corresponds to
      ! the sechiba time step but then in stomate.f90 resp_hetero_soil is multiplied by dt.
      ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day
      ! time step and avoid some constructions that could create bug during future developments.
      !
      !       VPP killer:
      !       resp_hetero_soil(:,m) = SUM( frac_resp(:,:) * fluxtot(:,:), DIM=2 ) / dt
      
      resp_hetero_soil(:,m) = &
         ( frac_resp(:,iactive) * fluxtot(:,iactive) + &
         frac_resp(:,islow) * fluxtot(:,islow) + &
         frac_resp(:,ipassive) * fluxtot(:,ipassive)  ) / dt
      
      !! 3.2.3 add fluxes to active, slow, and passive pools
      !       VPP killer:
      !       carbon(:,:,m) = carbon(:,:,m) + SUM( flux(:,:,:), DIM=2 )
    
      DO k = 1, ncarb ! Loop over carbon pools
        carbon(:,k,m) = carbon(:,k,m) + &
           flux(:,iactive,k) + flux(:,ipassive,k) + flux(:,islow,k)
      ENDDO ! Loop over carbon pools

!! Chloe TEST VALIDATION GUETTE :
!! 217.8 = moyenne de clabile du peat des HL de la moyenne annuelle de 2013
!! 82.15 = moyenne clabile des HL de la moy annuelle 2013 si fpeat > 10%
!        carbon(:,1,14)=217.8
!         carbon(:,1,14)=82.15
!         carbon(:,1,14)=130 
     
    ENDDO ! End loop over PFTs
   
    !!!!! Chloe ++ :
     ! write(*,*) 'Chloe dt,one_day',dt,one_day
       !!!!! Conservation of carbon mass when we use  peatland emission model :
       !DO k=1,ncarb
       !IF (m .EQ. 14 .AND. k .EQ. 1) THEN      
        !carbon(:,1,14) = carbon(:,1,14) - wpro_conso_carb(:)*dt
        !    DO i= 1,npts
        !    IF (carbon(i,1,14) .LT. 0) carbon(:,1,14)=zero
        !    ENDDO
        !IF (wpro_conso_carb(1) .GT. 0) THEN
        !    write(*,*) 'carbon,wpro,dt',carbon(1,1,14),wpro_conso_carb(1),dt
        !ENDIF
       !ENDIF
       !ENDDO
    !!! Chloe --

    
    IF (bavard.GE.4) WRITE(numout,*) 'Leaving soilcarbon'
    
  END SUBROUTINE soilcarbon

END MODULE stomate_soilcarbon
