! =================================================================================================================================
! MODULE 	: qsat_moisture
!
! CONTACT       : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE       : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "qsat_moisture" module contains public tools functions like qsat, dev_qsat.   
!!
!!\n DESCRIPTION: This module is the result of the splitting of constantes_veg.\n
!!                As the subroutines qsatcalc, dev_qsatcalc are used only by enerbil and diffuco, they are part of SECHIBA 
!!                component. 
!!
!! REFERENCE(S)	: 
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2012-03-08 09:35:49 +0100 (Thu, 08 Mar 2012) $
!! $Revision: 737 $
!! \n
!_ ================================================================================================================================

MODULE qsat_moisture

  USE defprec
  USE constantes
  USE IOIPSL

!-
  IMPLICIT NONE
!-

  LOGICAL,SAVE :: l_qsat_first=.TRUE.                !! First call to qsat subroutines and functions (true/false)


  INTEGER(i_std),PARAMETER :: max_temp=370           !! Maximum temperature for saturated humidity (K) and also used as 
                                                     !! the size of local array to keep saturated humidity (unitless)

  INTEGER(i_std),PARAMETER :: min_temp=100           !! Minimum temperature for saturated humidity (K)

  REAL(r_std),DIMENSION(max_temp),SAVE :: qsfrict    !! Array to keep saturated humidity at each temperature level 
                                                     !! (kg of water/kg of air)

  CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: qsatcalc
!! 
!>\BRIEF          This routine calculates the saturated humidity using the pressure
!!                and the temperature for all pixels. 
!!
!! DESCRIPTION : This routine interpolates qsat between temperatures by the following formula :
!!               \latexonly
!!               \input{qsatcalc.tex}
!!               \endlatexonly
!!               \n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S) : qsat_out
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

    SUBROUTINE qsatcalc (kjpindex,temp_in,pres_in,qsat_out)

      IMPLICIT NONE

      !! 0. Variables and parameters declaration
      
      !! 0.1 Input variables
      
      INTEGER(i_std),INTENT(in) :: kjpindex                   !! Domain size (unitless)
      REAL(r_std),DIMENSION(kjpindex),INTENT(in)  :: temp_in  !! Temperature in degre Kelvin (K)
      REAL(r_std),DIMENSION(kjpindex),INTENT(in)  :: pres_in  !! Pressure (Pa)
      
      !! 0.2 Output variables 
      
      REAL(r_std),DIMENSION(kjpindex),INTENT(out) :: qsat_out !! Saturated humidity at the surface (kg of water/kg of air)
      
      !! 0.4 Local variables
      
      INTEGER(i_std), DIMENSION(kjpindex) :: jt               !! Temporary array stocking the truncated temperatures in Kelvin 
                                                              !!(converted into integers)
      INTEGER(i_std)                      :: ji               !! indices (unitless)
      REAL(r_std),DIMENSION(kjpindex)     :: zz_a, zz_b, zz_f !! Temporary variables
      INTEGER(i_std)                      :: nbad             !! Number of points where the temperature is too high or too low
      INTEGER(i_std),DIMENSION(1)         :: lo               !! Temporary vector to mark the position of the highest temperature
                                                              !! or the lowest temperature over all the pixels in jt (unitless)
!_ ================================================================================================================================

      !-
      !! 1.Initialize qsfrict array if needed
      !-
      IF (l_qsat_first) THEN
         !-
         CALL qsfrict_init
         l_qsat_first = .FALSE.
         !-
      ENDIF !(l_qsat_first)
      
      !-
      !! 2. Computes qsat interpolation into two successive temperature
      !-
      jt = INT(temp_in(:))
       
      !! 2.1 Diagnostic pixels where the temperature is too high 
      nbad = COUNT(jt(:) >= max_temp-1)
      
      IF (nbad > 0) THEN
         WRITE(numout,*) ' qsatcalc: temperature too high at ', &
              &    nbad, ' points.'
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(2,'qsatcalc','diffuco', '', &
                 &                     'temperature incorect.') ! Warning message
         ELSE
            lo = MAXLOC(temp_in(:))
            WRITE(numout,*) &
                 &     'Maximum temperature ( ',MAXVAL(temp_in),') found at ',lo(1)
            WHERE (jt(:) >= max_temp-1)   jt(:) = max_temp-1
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF ! (nbad > 0)
 
      
      !! 2.2 Diagnostic pixels where the temperature is too low
      nbad = COUNT(jt(:) <= min_temp)
      
      IF (nbad > 0) THEN
         WRITE(numout,*) ' qsatcalc: temperature too low at ', &
              &    nbad, ' points.'
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(2,'qsatcalc','diffuco', '', &
                 &                   'temperature incorect.') ! Warning message
         ELSE
            lo = MINLOC(temp_in(:))
            WRITE(numout,*) &
                 &     'Minimum temperature ( ',MINVAL(temp_in),') found at ',lo(1)
            WHERE (jt(:) <= min_temp)   jt(:) = min_temp
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF! (nbad > 0)
 
      !! 2.3 Temporary variables needed for interpolation
      DO ji = 1, kjpindex ! Loop over # pixels
         
         zz_f(ji) = temp_in(ji)-FLOAT(jt(ji))
         zz_a(ji) = qsfrict(jt(ji))
         zz_b(ji) = qsfrict(jt(ji)+1)
         
      ENDDO ! Loop over # pixels
      
      !-
      !! 3. Interpolation between these two values
      !-
      DO ji = 1, kjpindex ! Loop over # pixels
         
         qsat_out(ji) = ((zz_b(ji)-zz_a(ji))*zz_f(ji)+zz_a(ji))/pres_in(ji)
      
         !write(*,*) 'Chloe qsat_out, pres_in', qsat_out, pres_in          
   
      ENDDO  ! Loop over # pixels


    END SUBROUTINE qsatcalc
 
!! ================================================================================================================================
!! FUNCTION 	: [DISPENSABLE] qsat
!!
!>\BRIEF         This function computes deviation the saturated humidity with the pressure
!! and the temperature for a scalar. 
!!
!! DESCRIPTION : This routine is obsolete : replaced by the subroutine qsatcalc. \n
!!               qsat is interpolated by : \n
!!               \latexonly
!!               \input{qsat.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): None\n
!!
!! RETURN VALUE : qsat_result
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================    

    FUNCTION qsat (temp_in,pres_in) RESULT (qsat_result)

      IMPLICIT NONE

      !! 0. Variables and parameters declaration

      !! 0.1 Input variables

      REAL(r_std),INTENT(in) :: temp_in   !! Temperature (K)
      REAL(r_std),INTENT(in) :: pres_in   !! Pressure    (Pa)

      !! 0.2 Result
      
      REAL(r_std) :: qsat_result          !! Saturated humidity calculated at the surface (kg/kg)
      
      !! 0.4 Local variables
      
      INTEGER(i_std)   :: jt              !! Temporary scalar stocking the truncated temperature in Kelvin 
                                          !! (converted into integer)
      REAL(r_std)      :: zz_a,zz_b,zz_f  !! Temporary scalar variables      
 
!_ ================================================================================================================================

      !-
      !! 1.Initialize qsfrict array if needed
      !-
      IF (l_qsat_first) THEN
         !-
         CALL qsfrict_init
         l_qsat_first = .FALSE.
         !-
      ENDIF

      !-
      !! 2. Computes qsat interpolation into two successive temperatures
      !-
      jt = INT(temp_in)

      !! 2.1 Is the temperature too high ?
      IF (jt >= max_temp-1) THEN
         WRITE(numout,*) &
              &   ' We stop. temperature too BIG : ',temp_in, &
              &   ' approximation for : ',jt
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(2,'qsat','', '',&
                 &                   'temperature incorect.') ! Warning message
         ELSE
            qsat_result = 999999.
            RETURN
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF !(jt >= max_temp-1)

      !! 2.2 Is the temperature too low ?
      IF (jt <= min_temp ) THEN
         WRITE(numout,*) &
              &   ' We stop. temperature too SMALL : ',temp_in, &
              &   ' approximation for : ',jt
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(2,'qsat','', '',&
                 &                   'temperature incorect.')
         ELSE
            qsat_result = -999999.
            RETURN
         ENDIF!(.NOT.diag_qsat)
         !-
      ENDIF !(jt <= min_temp )

      !! 2.3 Temporary variables needed for interpolation
      zz_f = temp_in-FLOAT(jt)
      zz_a = qsfrict(jt)
      zz_b = qsfrict(jt+1)

      !! 3. Interpolates between these two values

      qsat_result = ((zz_b-zz_a)*zz_f+zz_a)/pres_in
!      write(*,*) 'Chloe qsat toot high qsat, pres_in', qsat_result, pres_in

    END FUNCTION qsat


!! ================================================================================================================================
!! SUBROUTINE 	: dev_qsatcalc
!!
!>\BRIEF         This routine calculates the deviation of the saturated humidity qsat.
!!
!! DESCRIPTION : The deviation of qsat is calculated by :
!!               \latexonly
!!               \input{dev_qsatcalc.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S) : dev_qsat_out
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART	: None
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics[scale = 1]{pheno_moigdd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================ 

    SUBROUTINE dev_qsatcalc (kjpindex,temp_in,pres_in,dev_qsat_out)

      IMPLICIT NONE

      !! 0. Variables and parameters declaration
      
      !! 0.1 Input variables

      INTEGER(i_std),INTENT(in)                   :: kjpindex       !! Domain size (unitless)
      REAL(r_std),DIMENSION(kjpindex),INTENT(in)  :: temp_in        !! Temperature (K)
      REAL(r_std),DIMENSION(kjpindex),INTENT(in)  :: pres_in        !! Pressure


      !! 0.2 Output variables 
      
      REAL(r_std),DIMENSION(kjpindex),INTENT(out) :: dev_qsat_out   !! Result (??units??) 

      
      !! 0.4 Local variables
      
      INTEGER(i_std),DIMENSION(kjpindex) :: jt                      !! Temporary array stocking the truncated temperatures
                                                                    !! in Kelvin (converted into integers) 
      INTEGER(i_std)                     :: ji                      !! Indice (unitless)
      REAL(r_std),DIMENSION(kjpindex)    :: zz_a, zz_b, zz_c, zz_f  !! Temporary vector variables
      INTEGER(i_std)                     :: nbad                    !! Number of points where the temperature is too high or too low

!_ ================================================================================================================================

      !-
      !! 1.Initialize qsfrict array if needed
      !-
      IF (l_qsat_first) THEN
         !-
         CALL qsfrict_init
         l_qsat_first = .FALSE.
         !-
      ENDIF

      !-
      !! 2. Compute qsat interpolation into two successive temperature
      !-
      jt = INT(temp_in(:)+undemi)
      
      !! 2.1 Pixels where the temperature is too high 
      nbad = COUNT( jt(:) >= max_temp-1 )
      
      IF (nbad > 0) THEN
         WRITE(numout,*) &
              &   ' dev_qsatcalc: temperature too high at ',nbad,' points.'
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(3,'dev_qsatcalc','', '', &
                 &                   'temperature incorect.') ! Fatal error
         ELSE
            WHERE (jt(:) >= max_temp-1)   jt(:) = max_temp-1
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF !(nbad > 0)
      
      !! 2.2 Pixels where the temperature is too low
      nbad = COUNT( jt(:) <= min_temp )
      
      IF (nbad > 0) THEN
         WRITE(numout,*) &
              &   ' dev_qsatcalc: temperature too low at ',nbad,' points.'
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(3,'dev_qsatcalc', '', '',&
                 &                   'temperature incorect.') ! Fatal error
         ELSE
            WHERE (jt(:) <= min_temp)   jt(:) = min_temp
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF !(nbad > 0)

      !! 2.3 Temporary variables needed for interpolation
      DO ji=1,kjpindex  ! Loop over # pixels
         
         zz_f(ji) = temp_in(ji)+undemi-FLOAT(jt(ji))
         zz_a(ji) = qsfrict(jt(ji)-1)
         zz_b(ji) = qsfrict(jt(ji))
         zz_c(ji) = qsfrict(jt(ji)+1)
         
      ENDDO ! Loop over # pixels
      
      !-
      !! 3. Interpolates between these two values
      !-
      DO ji = 1, kjpindex ! Loop over # pixels
         
         dev_qsat_out(ji) = &
              &   ((zz_c(ji)-deux*zz_b(ji)+zz_a(ji))*(zz_f(ji)-un) + &
              &                         zz_c(ji)-zz_b(ji))/pres_in(ji)
         
      ENDDO ! Loop over # pixels


    END SUBROUTINE dev_qsatcalc

!! ================================================================================================================================
!! FUNCTION 	: [DISPENSABLE] dev_qsat
!!
!>\BRIEF         This function computes deviation of qsat. 
!!
!! DESCRIPTION : The deviation of qsat is calculated by :
!!               \latexonly
!!               \input{dev_qsat.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : dev_qsat_result
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

    FUNCTION dev_qsat (temp_in,pres_in) RESULT (dev_qsat_result)

      IMPLICIT NONE

      !! 0. Variables and parameters declaration 
   
      !! 0.1 Input variables

      REAL(r_std),INTENT(in)  :: pres_in          !! Pressure (Pa)
      REAL(r_std),INTENT(in)  :: temp_in          !! Temperture (K)
      
      !! 0.2 Result

      REAL(r_std) :: dev_qsat_result              !! (??units??) !!
      
      !! 0.4 Local variables

      INTEGER(i_std)  :: jt                       !! Index (unitless)
      REAL(r_std)     :: zz_a, zz_b, zz_c, zz_f   !! Temporary scalars   

!_ ================================================================================================================================
 
      !-
      !! 1.Initialize qsfrict array if needed
      !-
      IF (l_qsat_first) THEN
         !-
         CALL qsfrict_init
         l_qsat_first = .FALSE.
         !-
      ENDIF

      !-
      !! 2. computes qsat deviation interpolation
      !!    into two successive temperature
      !-
      jt = INT(temp_in+undemi)

      !! 2.1 Is the temperature too high ?
      IF (jt >= max_temp-1) THEN
         !-
         WRITE(numout,*) &
              &   ' We stop. temperature too HIGH : ',temp_in, &
              &   ' approximation for : ',jt
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(3,'dev_qsat','', '',&
                 &                   'temperature incorect.') ! Fatal error
         ELSE
            dev_qsat_result = 999999.
            RETURN
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF !(jt >= max_temp-1)
      !-
      !! 2.2 Is the temperature too low ?
      IF (jt <= min_temp ) THEN
         WRITE(numout,*) &
              &   ' We stop. temperature too LOW : ',temp_in, &
              &   ' approximation for : ',jt
         !-
         IF (.NOT.diag_qsat) THEN
            CALL ipslerr(3,'dev_qsat','', '',&
                 &                    'temperature incorect.')
         ELSE
            dev_qsat_result = -999999.
            RETURN
         ENDIF !(.NOT.diag_qsat)
         !-
      ENDIF !(jt <= min_temp )

      !! 2.3  Temporary variables for interpolation
      zz_f = temp_in+undemi-FLOAT(jt)
      zz_a = qsfrict(jt-1)
      zz_b = qsfrict(jt)
      zz_c = qsfrict(jt+1)

      !-
      !! 3. Interpolate
      !-
      dev_qsat_result=((zz_c-deux*zz_b+zz_a)*(zz_f-un)+zz_c-zz_b)/pres_in

    END FUNCTION dev_qsat
 

!! ================================================================================================================================
!! SUBROUTINE 	: qsfrict_init
!!
!>\BRIEF         The qsfrict_init routine initialises qsfrict array to store 
!!               precalculated values for qsat by using Goff-Gratch equations.  
!!
!! DESCRIPTION : This routine calculates the specific humidity qsat as a function of temperature in 
!!               Kelvin by using the modified Goff-Gratch equations(1946): \n
!!               \latexonly
!!               \input{goff_gratch.tex}
!!               \endlatexonly
!!               qsfrict is initialized by the following formulas : \n
!!               \latexonly
!!               \input{qsfrict_init.tex}
!!               \endlatexonly
!!               These values are used by the subroutines qsatcalc, dev_qsat. \n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::qsfrict
!!
!! REFERENCE(S)	:
!! - Algorithme d'un ensemble de paramétrisation physique (1998), 
!! Note de Laurent Li décrivant les paramétrisations physiques incluses dans le modèle (pdf),
!! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques 
!! - Goff, J. A., and S. Gratch (1946) Low-pressure properties of water from −160 to 212 °F, in Transactions of the
!! American Society of Heating and Ventilating Engineers, pp 95–122, presented at the 52nd annual meeting of the
!! American Society of Heating and Ventilating Engineers, New York, 1946.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

    SUBROUTINE qsfrict_init
      
      IMPLICIT NONE

      !! 0. Variables and parameters declaration
      
      !! 0.4 Local variables

      INTEGER(i_std) :: ji                             !! Indice(unitless)
      REAL(r_std)    :: zrapp,zcorr,ztemperature,zqsat !! Temporary vector variables      

!_ ================================================================================================================================

      !! 1. Initialisation
      zrapp = msmlr_h2o/msmlr_air
      zcorr = 0.00320991_r_std
      
      !! 2. Computes saturated humidity one time and store in qsfrict local array
      DO ji=100,max_temp ! Loop over size(qsfrict) : each position of qsfrict matches a temperature 
         
         ztemperature = FLOAT(ji)
         !-
         IF (ztemperature < 273._r_std) THEN 
            zqsat = zrapp*10.0_r_std**(2.07023_r_std-zcorr*ztemperature &
                 &             -2484.896/ztemperature+3.56654*LOG10(ztemperature)) ! Equilibrium water vapor - solid
         ELSE 
            zqsat = zrapp*10.0**(23.8319-2948.964/ztemperature & 
                 &              -5.028*LOG10(ztemperature) &
                 &              -29810.16*EXP(-0.0699382*ztemperature) &
                 &              +25.21935*EXP(-2999.924/ztemperature)) ! Equilibrium water vapor - liquid
         ENDIF !(ztemperature < 273._r_std)
         !-
         qsfrict (ji) = zqsat
         
      ENDDO ! Loop over size(qsfrict) 

      !! 3. Set to zero the non-computed values
      qsfrict(1:100) = zero 
      !-
      IF (long_print) WRITE (numout,*) ' qsfrict_init done'


    END SUBROUTINE qsfrict_init


END MODULE qsat_moisture
