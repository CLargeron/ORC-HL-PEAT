! =================================================================================================================================
! MODULE 	: stomate_assimtemp
!
! CONTACT	: orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculate the photosynthesis temperatures.
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_assimtemp.f90 $ 
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ =================================================================================================================================

MODULE stomate_assimtemp

  ! modules used:

  USE pft_parameters
  USE constantes 

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC assim_temp

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: assim_temp
!!
!>\BRIEF        This subroutine calculates the minimal, maximal and optimal 
!! temperatures for photosynthesis. 
!!
!! \n DESCRIPTION (definitions, functional, design, flags): The temperatures are 
!! calculated as second order polynomials of the long term reference temperature tlong_ref. 
!!
!! A climatological annual temperature is given as an input file of ORCHIDEE, at a 1.125°(lon)x1.121°(lat) spatial resolution. 
!! This temperature was derived from a IIASA database (Leemans and Cramer, 1991).
!! This initial climatological annual temperature is first interpolated at the spatial grid of the
!! meteorological forcing fields in stomate_io.f90::get_reftemp (\f$tlong\_ref\f$).
!! This long term reference temperature is then updated along the ORCHIDEE run at each stomate time step in stomate_season.f90.
!! 
!! The polynomial coefficients depend on the PFT and their values are defined in 
!! stomate_constants.f90. Coefficients for degrees 1 and 2 of the polynomials are 
!! zero, meaning the temperatures are just predefined constants (Wullschleger, 1993), except for C3
!! grasses, which has been designed specifically for ORCHIDEE (no reference publication other than Krinner et al. (2005)).
!!
!! This routine is called once at the beginning by stomate_var_init and then at each stomate time step by stomateLpj.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): minimal (t_photo_min), maximal (t_photo_max) and 
!! optimal temperatures (t_photo_opt) for photosynthesis
!!
!! REFERENCE(S)	: 
!! - Krinner, G., Viovy, N., Noblet-Ducoudre, N. de, Ogee, J., Polcher, J., 
!! Friedlingstein, P., Sitch, S., and Prentice, I. C (2005). A dynamic global 
!! vegetation model for studies of the coupled atmosphere-biosphere system. 
!! Global Biogeochem. Cycles, 19, GB1015
!! - Leemans, R. and Cramer,W. (1991). The IIASA Database for Mean Monthly Values
!! of Temperature, Precipitation, and Cloudiness on a Global Terrestrial Grid.
!! Research Report, INTERNATIONAL INSTITUTE FOR APPLIED SYSTEMS ANALYSIS Laxenburg, Austria.
!! International Standard Book Number 3-7045-0113-1. 
!! - Wullschleger, S.D. (1993). Biochemical limitations to carbon assimilation in C\f$_3\f$ plants: 
!! a retrospective analysis of the A/C\f$_i\f$ curves from 109 species. Journal of Experimental Botany, 44 (262), 907-920
!!
!! FLOWCHART    : None
!!
!! REVISION(S)	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE assim_temp (npts, tlong_ref, t2m_month, t_photo_min, t_photo_opt, t_photo_max)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                                  :: npts        !! Domain size    
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: tlong_ref   !! "long term" 2 meter reference temperatures (K)    
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: t2m_month   !! "monthly" 2-meter temperatures (K)

    !
    !! 0.2 Output variables 
    !
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)               :: t_photo_min !! Minimum temperature for photosynthesis (K)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)               :: t_photo_opt !! Optimum temperature for photosynthesis (K)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)               :: t_photo_max !! Maximum temperature for photosynthesis (K)

    !
    !! 0.3 Modified variables
    !

    !
    !! 0.4 Local variables
    !
    REAL(r_std), DIMENSION(npts)                                :: tl          !! "Long term" 2 meter reference temperatures 
                                                                               !! (degree C)    
    INTEGER(i_std)                                              :: j           !! Index (unitless)

!_ ================================================================================================================================

    tl(:) = tlong_ref(:) - ZeroCelsius

    DO j = 2,nvm	! Loop over # PFTs 

       !
       !! 1. Normal case
       !

       t_photo_min(:,j) = tphoto_min_c(j) + tphoto_min_b(j)*tl(:) + tphoto_min_a(j)*tl(:)*tl(:) + ZeroCelsius
       t_photo_opt(:,j) = tphoto_opt_c(j) + tphoto_opt_b(j)*tl(:) + tphoto_opt_a(j)*tl(:)*tl(:) + ZeroCelsius
       t_photo_max(:,j) = tphoto_max_c(j) + tphoto_max_b(j)*tl(:) + tphoto_max_a(j)*tl(:)*tl(:) + ZeroCelsius

       !
       !! 2. If the monthly temperature is too low, we set tmax < tmin.
       !!    Therefore, photosynthesis will not be possible (we need tmin < t < tmax).
       !     t2m_month is calculated at every stomate time step using a linear relaxation method, in stomate_season.f90 
       !     (see section [27] in Krinner et al. (2005)).
       !

!! min_stomate is a constant defined in stomate_constants.f90 (= 1.E-8_r_std).
       WHERE ( t2m_month(:) .LT. t_photo_min(:,j) )
          t_photo_max(:,j) = t_photo_min(:,j) - min_stomate
       ENDWHERE

    ENDDO	! Loop over # PFTs

  END SUBROUTINE assim_temp

END MODULE stomate_assimtemp
