! =================================================================================================================================
! MODULE       : stomate_prescribe
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Initialize and update density, crown area.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_prescribe.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_prescribe

  ! modules used:

  USE ioipsl
  USE stomate_data
  USE pft_parameters
  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC prescribe,prescribe_clear

    ! first call
    LOGICAL, SAVE                                              :: firstcall = .TRUE.

CONTAINS

! =================================================================================================================================
!! SUBROUTINE   : prescribe_clear
!!
!>\BRIEF        : Set the firstcall flag back to .TRUE. to prepare for the next simulation.
!_=================================================================================================================================

  SUBROUTINE prescribe_clear
    firstcall=.TRUE.
  END SUBROUTINE prescribe_clear

!! ================================================================================================================================
!! SUBROUTINE   : prescribe
!!
!>\BRIEF         Works only with static vegetation and agricultural PFT. Initialize biomass,
!!               density, crown area in the first call and update them in the following.
!!
!! DESCRIPTION (functional, design, flags): \n
!! This module works only with static vegetation and agricultural PFT.
!! In the first call, initialize density of individuals, biomass, crown area,
!! and leaf age distribution to some reasonable value. In the following calls,
!! these variables are updated.
!!
!! To fulfill these purposes, pipe model are used:
!! \latexonly 
!!     \input{prescribe1.tex}
!!     \input{prescribe2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLES(S): ::ind, ::cn_ind, ::leaf_frac
!!
!! REFERENCES   :
!! - Krinner, G., N. Viovy, et al. (2005). "A dynamic global vegetation model 
!!   for studies of the coupled atmosphere-biosphere system." Global 
!!   Biogeochemical Cycles 19: GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!   plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!   global vegetation model, Global Change Biology, 9, 161-185.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE prescribe (npts, &
                        veget_max, PFTpresent, everywhere, when_growthinit, &
                        biomass, leaf_frac, ind, cn_ind)

!! INTERFACE DESCRIPTION

!! INPUT SCALAR
    INTEGER(i_std), INTENT(in)                                :: npts            !! Domain size (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_max       !! "maximal" coverage fraction of a PFT (LAI -> infinity) on ground (unitless;0-1)

!! MODIFIED FIELD
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent      !! PFT present (0 or 1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere      !! is the PFT everywhere in the grid box or very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit !! how many days ago was the beginning of the growing season (days)
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)    :: biomass         !! biomass (gC/(m^2 of ground))
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac       !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind             !! density of individuals (1/(m^2 of ground))
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind          !! crown area per individual (m^2)

!! OUTPUT FIELDS 

!! LOCAL VARIABLES
    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass (gC/(m^2 of ground))
    REAL(r_std), DIMENSION(npts)                              :: woodmass_ind    !! woodmass of an individual (gC)
    INTEGER(i_std)                                            :: i,j             !! index (unitless)

!_ ================================================================================================================================

    DO j = 2,nvm

      ! only when the DGVM is not activated or agricultural PFT.

      IF ( ( .NOT. control%ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) ) THEN

        !
        !! 1.Update crown area
        !

        cn_ind(:,j) = zero

        IF ( tree(j) ) THEN

          !
          !! 1.1 treat for trees
          !

          dia(:) = zero

          DO i = 1, npts ! loop over grid points

            IF ( veget_max(i,j) .GT. zero ) THEN

              !! 1.1.1 calculate wood mass on an area basis, which include sapwood and heartwood aboveground and belowground.

              woodmass(i) = (biomass(i,j,isapabove) + biomass(i,j,isapbelow) + &
                   biomass(i,j,iheartabove) + biomass(i,j,iheartbelow)) * veget_max(i,j)  

              IF ( woodmass(i) .GT. min_stomate ) THEN

                !! 1.1.2 calculate critical individual density
!?? the logic for 1.1.3 and 1.1.2 is strange, it should be the case that first to calculate critical woodmass per individual,
!?? then calculate critical density.


                ! how to derive the following equation:
                ! first, TreeHeight=pipe_tune2 * Diameter^{pipe_tune3}
                ! we assume the tree is an ideal cylinder, so it volume is: Volume = pi*(Dia/2)^2*Height = pi/4 * Dia * pipe_tune2*Dia^{pipe_tune3} 
                !                                                                  = pi/4 * pipe_tune2 * Dia^{2+pipe_tune3}
                ! last, the woodmass_per_individual = pipe_density * Volume = pipe_density*pi/4.*pipe_tune2 * Dia^{2+pipe_tune3}             
                ind(i,j) = woodmass(i) / &
                           ( pipe_density*pi/4.*pipe_tune2 * maxdia(j)**(2.+pipe_tune3) )

                !! 1.1.3 individual biomass corresponding to this critical density of individuals

                woodmass_ind(i) = woodmass(i) / ind(i,j)

                !! 1.1.4 calculate stem diameter per individual tree

                dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                         ( un / ( 2. + pipe_tune3 ) )

                !! 1.1.5 calculate provisional tree crown area for per individual tree

                ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

                !! 1.1.6 If total tree crown area for this tree PFT exceeds its veget_max, tree density is recalculated.

                IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_max(i,j) ) THEN

                  ind(i,j) = veget_max(i,j) / cn_ind(i,j)

                ELSE

                   ind(i,j) = ( veget_max(i,j) / &
                        &     ( pipe_tune1 * (woodmass(i)/(pipe_density*pi/4.*pipe_tune2)) &
                        &     **(pipe_tune_exp_coeff/(2.+pipe_tune3)) ) ) &
                        &     ** (1./(1.-(pipe_tune_exp_coeff/(2.+pipe_tune3))))
                   

                  woodmass_ind(i) = woodmass(i) / ind(i,j)

                  dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                           ( un / ( 2. + pipe_tune3 ) )

                  ! final crown area
                  cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

                ENDIF

              ELSE !woodmas=0  => impose some value

                dia(:) = maxdia(j)

                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

              ENDIF ! IF ( woodmass(i) .GT. min_stomate )

            ENDIF    ! veget_max .GT. 0.

          ENDDO      ! loop over grid points

        ELSE !grass

          !
          !! 1.2 grasses: crown area always set to 1m**2
          !

          WHERE ( veget_max(:,j) .GT. zero )
            cn_ind(:,j) = un
          ENDWHERE

        ENDIF   !IF ( tree(j) )

        !
        !! 2 density of individuals
        !
        
        WHERE ( veget_max(:,j) .GT. zero )

          ind(:,j) = veget_max(:,j) / cn_ind(:,j)  

        ELSEWHERE

          ind(:,j) = zero

        ENDWHERE

      ENDIF     ! IF ( ( .NOT. control%ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) )

    ENDDO       ! loop over PFTs

    !
    !!? it's better to move the code for first call at the beginning of the module.
    !! 2 If it's the first call for this module, 
    !

    IF ( firstcall ) THEN

      WRITE(numout,*) 'prescribe:'

      ! impose some biomass if zero and PFT prescribed

      WRITE(numout,*) '   > Imposing initial biomass for prescribed trees, '// &
                      'initial reserve mass for prescribed grasses.'
      WRITE(numout,*) '   > Declaring prescribed PFTs present.'

      DO j = 2,nvm ! loop over PFTs
        DO i = 1, npts ! loop over grid points

          ! is vegetation static or PFT agricultural?
          ! Static vegetation or agricultural PFT
          IF ( ( .NOT. control%ok_dgvm ) .OR. &
               ( ( .NOT. natural(j) ) .AND. ( veget_max(i,j) .GT. min_stomate ) ) ) THEN

            !
            !! 2.1 if tree biomass is extremely small, prescribe the biomass by assuming they have sapling biomass, which is a constant in the model.
            !!     then set all the leaf age as 1.
            !
            ! if tree PFT and biomass too small, prescribe the biomass to a value.
            IF ( tree(j) .AND. &
                 ( veget_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:) ) .LE. min_stomate ) ) THEN
               !!? here the code is redundant, as "veget_max(i,j) .GT. min_stomate" is already met in the above if condition.
               IF (veget_max(i,j) .GT. min_stomate) THEN
                  biomass(i,j,:) = (bm_sapl_rescale * bm_sapl(j,:) * ind(i,j)) / veget_max(i,j)
               ELSE
                  biomass(i,j,:) = zero
               ENDIF

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

              ! seasonal trees: no leaves at beginning

              IF ( pheno_model(j) .NE. 'none' ) THEN

                biomass(i,j,ileaf) = zero
                leaf_frac(i,j,1) = zero

              ENDIF

            ENDIF

            !
            !! 2.2 for grasses, set only the carbon reserve pool to "sapling" carbon reserve pool.
            !!     and set all leaf age to 1.

            IF ( ( .NOT. tree(j) ) .AND. &
                 ( veget_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:) ) .LE. min_stomate ) ) THEN

              biomass(i,j,icarbres) = bm_sapl(j,icarbres) * ind(i,j) / veget_max(i,j)

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

            ENDIF

            !
            !! 2.3 declare all PFTs with positive veget_max as present everywhere in that grid box
            !

            IF ( veget_max(i,j) .GT. min_stomate ) THEN
              PFTpresent(i,j) = .TRUE.
              everywhere(i,j) = un
            ENDIF

          ENDIF   ! not control%ok_dgvm  or agricultural

        ENDDO ! loop over grid points
      ENDDO ! loop over PFTs

      firstcall = .FALSE.

    ENDIF

  END SUBROUTINE prescribe

END MODULE stomate_prescribe
