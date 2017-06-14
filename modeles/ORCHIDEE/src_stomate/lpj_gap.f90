! =================================================================================================================================
! MODULE       : lpj_gap
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Simulate mortality of individuals and update biomass, litter and 
!! stand density of the PFT
!!
!!\n DESCRIPTION : Simulate mortality of individuals and update biomass, litter and 
!! stand density of the PFT. This module differs from lpj_kill.f90 in that this 
!! module kills individuals within a PFT and lpj_kill.f90 removes a PFT from a 
!! gridbox
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : 
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!         plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!         global vegetation model, Global Change Biology, 9, 161-185.\n
!! - Waring, R. H. (1983). "Estimating forest growth and efficiency in relation 
!!         to canopy leaf area." Advances in Ecological Research 13: 327-354.\n
!! 
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/lpj_gap.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE lpj_gap

  ! modules used:

  USE ioipsl
  USE stomate_data
  USE pft_parameters
  USE parallel
  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC gap,gap_clear
  
  ! Variable declaration 

  LOGICAL, SAVE             :: firstcall = .TRUE.                !! first call flag

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : gap_clear
!!
!>\BRIEF        Set the firstcall flag back to .TRUE. to prepare for the next simulation.
!_ ================================================================================================================================
  
  SUBROUTINE gap_clear
    firstcall = .TRUE.
  END SUBROUTINE gap_clear


!! ================================================================================================================================
!! SUBROUTINE   : gap
!!
!>\BRIEF        Simulate tree and grass mortality, transfer dead biomass to litter and update stand density
!!
!! DESCRIPTION  : Calculate mortality of trees and grasses, transfer the dead biomass to litter pool,
!! and update biomass pool and number of individuals. To get tree mortality, it's possible to choose either a 
!! constant mortality rate; or to calculate the tree mortality rate based on it's growth efficiency, which is 
!! defined as this year's net biomass increment per leaf area.\n
!!
!! When using growth efficiency mortality, first calculate the net biomass increment for the last year, then 
!! calculate the growth efficiency and finally calculate the growth efficiency mortality.\n
!!
!! Eqation to calculate growth efficiency:
!! \latexonly 
!!     \input{gap1.tex}
!! \endlatexonly
!! Where $greff$ is growth efficiency, $\Delta$ is net biomass increment, 
!! $C_{leaf} is last year's leaf biomass, and $SLA$ the specific leaf area. 
!! 
!!
!! Eqation to calculate growth efficiency mortality:
!! \latexonly 
!!     \input{gap2.tex}
!! \endlatexonly
!! Where $mort_{greff}$ is the growth efficiency mortality, $greff$ is growth  
!! efficiency, $k_{mort1}$ is asymptotic maximum mortality rate.
!! 
!! The name for variable ::availability is not well chosen. Actually the meaning of the variable is mortailty 
!! rate derived from growth efficiency related mortality. ?? Suggestion: change the name "availability" to 
!! "mortgref", which means "mortality caused by ".\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::biomass; biomass, ::ind density of individuals, ::bm_to_litter biomass transfer 
!! to litter and ::mortality mortality (fraction of trees that is dying per time step)
!!
!! REFERENCE(S)   :
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!         plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!         global vegetation model, Global Change Biology, 9, 161-185.
!! - Waring, R. H. (1983). "Estimating forest growth and efficiency in relation 
!!         to canopy leaf area." Advances in Ecological Research 13: 327-354.
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
  
  SUBROUTINE gap (npts, dt, &
       npp_longterm, turnover_longterm, lm_lastyearmax, &
       PFTpresent, biomass, ind, bm_to_litter, mortality)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                              :: npts                    !! Domain size (-)
    REAL(r_std), INTENT(in)                                 :: dt                      !! Time step (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: npp_longterm            !! "Long term" (default 3-year) net primary
                                                                                       !! productivity  
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)     :: turnover_longterm       !! "Long term" (default 3-year) turnover  
                                                                                       !! rate @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: lm_lastyearmax          !! Last year's maximum leaf mass, for each 
                                                                                       !! PFT @tex $(gC m^{-2})$ @endtex
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                :: PFTpresent              !! Is the pft present in the pixel

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm),INTENT(out)            :: mortality               !! Mortality (fraction of trees that is 
                                                                                       !! dying per time step)

    !! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)  :: biomass                 !! Biomass @tex $(gC m^{-2}) $@endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)         :: ind                     !! Number of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(inout)  :: bm_to_litter            !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                            :: delta_biomass           !! Net biomass increase for the previous 
                                                                                       !! year @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts)                            :: dmortality              !! Dead biomass caused by mortality 
                                                                                       !! @tex $(gC m^{-2}) $@endtex
    REAL(r_std), DIMENSION(npts)                            :: vigour                  !! Growth efficiency, an indicator of tree 
                                                                                       !! vitality, used to calculate mortality
    REAL(r_std), DIMENSION(npts)                            :: availability            !! Mortality rate derived by growth 
                                                                                       !! efficiency @tex $(year^{-1})$ @endtex 
    INTEGER(i_std)                                          :: j,k,m                   !! Indices

!_ ================================================================================================================================

   IF ( firstcall ) THEN

       firstcall = .FALSE.

    ENDIF

    IF (bavard.GE.3) WRITE(numout,*) 'Entering gap',lpj_gap_const_mort

    mortality(:,:) = zero

    ! loop over #PFT
    DO j = 2,nvm

 !! 1. Tree mortality

       IF ( tree(j) ) THEN 

          !! 1.1 Use growth efficiency or constant mortality?
          IF ( .NOT.  lpj_gap_const_mort  ) THEN

             !! 1.1.1 Estimate net biomass increment
             !        To calculate growth efficiency mortality, first estimate net biomass increment by 
             !        subtracting turnover from last year's NPP.
             WHERE ( PFTpresent(:,j) .AND. ( lm_lastyearmax(:,j) .GT. min_stomate ) )

            !??! the following should be removed
            ! note that npp_longterm is now actually longterm growth efficiency (NPP/LAI)
            ! to be fair to deciduous trees
            ! calculate net biomass increment
             delta_biomass(:) = MAX( npp_longterm(:,j) - ( turnover_longterm(:,j,ileaf) + &
                  turnover_longterm(:,j,iroot) + turnover_longterm(:,j,ifruit) + & 
                  turnover_longterm(:,j,isapabove) + turnover_longterm(:,j,isapbelow) ) ,zero)

            !! 1.1.2 Calculate growth efficiency 
            !        Calculate growth efficiency by dividing net biomass increment by last year's 
            !        maximum LAI. (corresponding actually to the maximum LAI of the previous year)
             vigour(:) = delta_biomass(:) / (lm_lastyearmax(:,j)*sla(j))

             ELSEWHERE

                vigour(:) = zero

             ENDWHERE

             !! 1.1.3 Calculate growth efficiency mortality rate
             WHERE ( PFTpresent(:,j) )
                
                availability(:) = availability_fact / ( un + ref_greff * vigour(:) )
                ! Scale mortality by timesteps per year
                mortality(:,j) = MAX(min_avail,availability(:))  * dt/one_year  

             ENDWHERE

          ELSE  ! .NOT. lpj_gap_const_mort

             !! 1.2 Use constant mortality accounting for the residence time of each tree PFT
             WHERE ( PFTpresent(:,j) )

                mortality(:,j) = dt/(residence_time(j)*one_year)

             ENDWHERE

          ENDIF ! .NOT.  lpj_gap_const_mort

          !! 1.3 Mortality in DGVM
          !      If the long term NPP is zero, all trees are killed
          !??! This is this only applied in the DGVM maybe in order to make the DGVM respond faster and thus make the vegetation dynamics more dynamic?
          !??! the link here with lpj_kill.f90 is still not clear and so would leave to who especially working on this. 
          IF ( control%ok_dgvm ) THEN

             WHERE ( PFTpresent(:,j) .AND. ( npp_longterm(:,j) .LE. min_stomate ) )

                mortality(:,j) = un

             ENDWHERE

          ENDIF

          !! 1.4 Update biomass and litter pools 
          !    Update biomass and litter pool after dying and transfer recently died biomass to litter
          DO k = 1, nparts

             WHERE ( PFTpresent(:,j) )

                dmortality(:) =  mortality(:,j) * biomass(:,j,k)
                bm_to_litter(:,j,k) = bm_to_litter(:,j,k) + dmortality(:)
                
                biomass(:,j,k) = biomass(:,j,k) - dmortality(:)

             ENDWHERE

          ENDDO

          !! 1.5 In case of dynamic vegetation, update tree individual density
          IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN

             WHERE ( PFTpresent(:,j) )

                ind(:,j) = ind(:,j) * ( un - mortality(:,j) )

             ENDWHERE

          ENDIF
       ELSE 

 !! 2. Grasses mortality

          ! For grasses, if last year's NPP is very small (less than 10 gCm^{-2}year{-1})
          ! the grasses completely die
          IF ( .NOT.control%ok_dgvm .AND. .NOT.lpj_gap_const_mort) THEN

             WHERE ( PFTpresent(:,j) .AND. ( npp_longterm(:,j) .LE. npp_longterm_init ) )

                mortality(:,j) = un

             ENDWHERE

             ! Update biomass and litter pools
             DO k = 1, nparts

                WHERE ( PFTpresent(:,j) )

                   dmortality(:) =  mortality(:,j) * biomass(:,j,k)
                   
                   bm_to_litter(:,j,k) = bm_to_litter(:,j,k) + dmortality(:)
                   
                   biomass(:,j,k) = biomass(:,j,k) - dmortality(:)

                ENDWHERE
             ENDDO
             
          ENDIF

       ENDIF   !IF ( tree(j) )

    ENDDO      !loop over pfts

 !! 3. Write to history files

    ! output in fraction of trees that dies/day.
    mortality = mortality / dt

    CALL histwrite (hist_id_stomate, 'MORTALITY', itime, &
         mortality, npts*nvm, horipft_index)

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving gap'

  END SUBROUTINE gap

END MODULE lpj_gap
