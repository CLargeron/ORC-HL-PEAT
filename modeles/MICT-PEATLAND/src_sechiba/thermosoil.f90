! =================================================================================================================================
! MODULE       : thermosoil
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculates the soil temperatures by solving the heat
!! diffusion equation within the soil
!!
!!\n DESCRIPTION : General important informations about the numerical scheme and
!!                 the soil vertical discretization:\n
!!               - the soil is divided into "ngrnd" (=7 by default) layers, reaching to as
!!                 deep as 5.5m down within the soil, with thiscknesses
!!                 following a geometric series of ration 2.\n
!!               - "jg" is usually used as the index going from 1 to ngrnd to describe the
!!                  layers, from top (jg=1) to bottom (jg=ngrnd)\n
!!               - the thermal numerical scheme is implicit finite differences.\n
!!                 -- When it is resolved in thermosoil_profile at the present timestep t, the
!!                 dependancy from the previous timestep (t-1) is hidden in the
!!                 integration coefficients cgrnd and dgrnd, which are therefore
!!                 calculated at the very end of thermosoil_main (call to
!!                 thermosoil_coef) for use in the next timestep.\n
!!                 -- At timestep t, the system becomes :\n 
!!
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!
!!                 (the bottom boundary condition has been used to obtained this equation).\n
!!                 To solve it, the uppermost soil temperature T(1) is required.
!!                 It is obtained from the surface temperature Ts, which is
!!                 considered a linear extrapolation of T(1) and T(2)\n
!!
!!                           Ts=(1-lambda)*T(1) -lambda*T(2) \n 
!!                                      -- EQ2--\n
!!
!!                 -- caveat 1 : Ts is called 'temp_soil_new' in this routine,
!!                 don' t act.\n
!!                 -- caveat 2 : actually, the surface temperature at time t Ts
!!                 depends on the soil temperature at time t through the
!!                 ground heat flux. This is again implicitly solved, with Ts(t)
!!                 expressed as :\n
!!
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflux+otherfluxes(Ts(t))\n 
!!                                      -- EQ3 --\n
!!
!!                 and the dependency from the previous timestep is hidden in
!!                 soilcap and soilflux (apparent surface heat capacity and heat
!!                 flux respectively). Soilcap and soilflux are therefore
!!                 calculated at the previsou timestep, at the very end of thermosoil
!!                 (final call to thermosoil_coef) and stored to be used at the next time step.
!!                 At timestep t, EQ3 is solved for Ts in enerbil, and Ts
!!                 is used in thermosoil to get T(1) and solve EQ1.\n
!!
!! - lambda is the @tex $\mu$ @endtex of F. Hourdin' s PhD thesis, equation (A28); ie the
!! coefficient of the linear extrapolation of Ts (surface temperature) from T1 and T2 (ptn(jg=1) and ptn(jg=2)), so that:\n
!! Ts= (1+lambda)*T(1)-lambda*T(2) --EQ2-- \n
!! lambda = (zz_coeff(1))/((zz_coef(2)-zz_coef(1))) \n
!!
!! - cstgrnd is the attenuation depth of the diurnal temperature signal
!! (period : one_day) as a result of the heat conduction equation
!! with no coefficients :
!!\latexonly
!!\input{thermosoil_var_init0.tex}
!!\endlatexonly
!!  -- EQ4 --\n
!! This equation results from the change of variables :
!! z' =z*sqrt(Cp/K) where z' is the new depth (homogeneous
!! to sqrt(time) ), z the real depth (in m), Cp and K the soil heat
!! capacity and conductivity respectively.\n
!!
!! the attenuation depth of a diurnal thermal signal for EQ4 is therefore homogeneous to sqrt(time) and
!! equals : \n
!! cstgrnd = sqrt(oneday/Pi)
!!
!! - lskin is the attenuation depth of the diurnal temperature signal
!! (period : one_day) within the soil for the complete heat conduction equation
!! (ie : with coefficients)
!!\latexonly
!!\input{thermosoil_var_init00.tex}
!!\endlatexonly
!! -- EQ5 --  \n
!! it can be retrieved from cstgrnd using the change of variable  z' =z*sqrt(Cp/K):\n
!! lskin = sqrt(K/Cp)*cstgrnd =  sqrt(K/Cp)*sqrt(oneday//Pi)\n
!! 
!! In thermosoil, the ratio lskin/cstgrnd is frequently used as the
!! multiplicative factor to go from
!!'adimensional' depths (like z' ) to real depths (z). z' is not really
!! adimensional but is reffered to like this in the code.
!!
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/thermosoil.f90 $
!! $Date: 2013-06-20 15:26:25 +0200 (Thu, 20 Jun 2013) $
!! $Revision: 1322 $
!! \n
!_ ================================================================================================================================

MODULE thermosoil

  USE ioipsl
 
  ! modules used :
  USE constantes
  USE constantes_soil
  USE sechiba_io
  USE grid
  USE parallel

  IMPLICIT NONE

  !private and public routines :
  PRIVATE
  PUBLIC :: thermosoil_main,thermosoil_clear 

  LOGICAL, SAVE                                   :: l_first_thermosoil=.TRUE.!! does the initialisation of the routine 
                                                                              !! (true/false)
  CHARACTER(LEN=80) , SAVE                        :: var_name                 !! To store variables names for the 
                                                                              !! input-outputs dealt with by IOIPSL
  REAL(r_std), SAVE                               :: lambda, cstgrnd, lskin   !! See Module description
  REAL(r_std), SAVE                               :: fz1, zalph               !! usefull constants for diverse use
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: ptn                      !! vertically discretized 
                                                                              !! soil temperatures @tex ($K$) @endtex. 
  REAL(r_std), SAVE, DIMENSION (ngrnd)            :: zz                       !! depths of the soil thermal numerical nodes. 
                                                                              !! Caveats: they are not exactly the centers of the
                                                                              !! thermal layers, see the calculation in 
                                                                              !! ::thermosoil_var_init  @tex ($m$) @endtex.
  REAL(r_std), SAVE, DIMENSION (ngrnd)            :: zz_coef                  !! depths of the boundaries of the thermal layers,
                                                                              !! see the calculation in 
                                                                              !! thermosoil_var_init  @tex ($m$) @endtex.
  REAL(r_std), SAVE, DIMENSION (ngrnd)            :: dz1                      !! numerical constant used in the thermal numerical
                                                                              !! scheme  @tex ($m^{-1}$) @endtex. ; it corresponds
                                                                              !! to the coefficient  @tex $d_k$ @endtex of equation
                                                                              !! (A.12) in F. Hourdin PhD thesis.
  REAL(r_std), SAVE, DIMENSION (ngrnd)            :: dz2                      !! thicknesses of the thermal layers  @tex ($m$)
                                                                              !! @endtex; typically: 
                                                                              !! dz2(jg)=zz_coef(jg+1)-zz_coef(jg); calculated once 
                                                                              !! and for all in thermosoil_var_init
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: z1                       !! constant of the numerical scheme; it is an 
                                                                              !! intermediate buffer for the calculation of the 
                                                                              !! integration coefficients cgrnd and dgrnd.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: cgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: dgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa                    !! volumetric vertically discretized soil heat 
                                                                              !! capacity  @tex ($J K^{-1} m^{-3}$) @endtex. 
                                                                              !! It depends on the soil
                                                                              !! moisture content (wetdiag) and is calculated at 
                                                                              !! each time step in thermosoil_coef.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa                   !! vertically discretized soil thermal conductivity 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex. Same as pcapa.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: zdz1                     !! numerical constant of the numerical scheme; it is
                                                                              !! an intermediate buffer for the calculation of the 
                                                                              !! integration coefficients cgrnd and dgrnd 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: zdz2                     !! numerical constant of the numerical scheme; it is 
                                                                              !! an intermediate buffer for the calculation of the 
                                                                              !! integration coefficients cgrnd and dgrnd
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_en                 !! heat capacity used for surfheat_incr and 
                                                                              !! coldcont_incr 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: ptn_beg                  !! vertically discretized temperature at the 
                                                                              !! beginning of the time step  @tex ($K$) @endtex; 
                                                                              !! is used in 
                                                                              !! thermosoil_energy for energy-related diagnostic of
                                                                              !! the routine.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: temp_sol_beg             !! Surface temperature at the beginning of the 
                                                                              !! timestep  @tex ($K$) @endtex
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: surfheat_incr            !! Change in soil heat content during the timestep 
                                                                              !!  @tex ($J$) @endtex.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: coldcont_incr            !! Change in snow heat content  @tex ($J$) @endtex.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: wetdiag                  !! Soil wetness on the thermodynamical levels 
                                                                              !! (1, ngrnd) (0-1, dimensionless). corresponds to the
                                                                              !! relative soil humidity to the wilting point when 
                                                                              !! the 11-layers hydrology (hydrol) is used, see more
                                                                              !! precisions in thermosoil_humlev.

 !Isa------------------------------------------------------------------------
  !  Permafrost-related constants: 
  !
  ! temperature range over which soil freezing takes place (K)
  REAL(r_std), PARAMETER        	 :: fr_dT =2

  ! Isa : the new porosity is a correction, old value was 0.15
  REAL(r_std), PARAMETER                :: poros = .41 ! (moyenne des mcs sur classif USDA)

! heat capacity of ice-filled saturated soil (J m**-3 K**-1, at porosity of 0.41)
!Isa : the new heat capacity for saturated frozen soil is a correction ; old
!value was 1.81e6 (corresponding to poros = 0.15), the new value is computed
!based on : so_capa_ice = so_capa_dry + poros*capa_ice*rho_ice, with
!capa_ice=2.06 J/Kg/K at 0°C.
! this is a divergence from Gouttevin et al., 2012, where so_capa_ice = 2.3e6
! was considered (mistake...)
 
  INTEGER, PARAMETER	   	 :: bavard = 4
  INTEGER, SAVE   	         :: flioid	  !! flio output file ID


!Isa 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: profil_froz
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: wetdiaglong          !! Long-term soil humidity (for permafrost) if ok_freeze_thermix ; wetdiag sinon.
  LOGICAL, SAVE    				   :: ok_freeze_thermix
!Isa E
! CR: option supl pour converg avec trunc
        LOGICAL, SAVE    :: ok_thermix_trunc
        LOGICAL, SAVE    :: ok_wetdiaglong
!        LOGICAL, SAVE    :: ok_converge_isaorig ! on le met dans constantes

    REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:)    :: pcappa_supp
     REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:)    :: E_sol_lat_couche
    LOGICAL, SAVE     				   :: ok_Ecorr
!Isa permafrost_map
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:) :: overburden
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:) :: excess_ice
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:) :: permafrost
    LOGICAL, SAVE     				   :: read_permafrost_map
!Isa reftemp
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:) :: reftemp
    LOGICAL, SAVE     				   :: read_reftemp
    ! end isa

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_main
!!
!>\BRIEF        Thermosoil_main computes the soil thermal properties and dynamics, ie solves
!! the heat diffusion equation within the soil. The soil temperature profile is
!! then interpolated onto the diagnostic axis.
!!
!! DESCRIPTION : The resolution of the soil heat diffusion equation 
!! relies on a numerical finite-difference implicit scheme
!! fully described in the reference and in the header of the thermosoil module.
!! - The dependency of the previous timestep hidden in the 
!! integration coefficients cgrnd and dgrnd (EQ1), calculated in thermosoil_coef, and 
!! called at the end of the routine to prepare for the next timestep.
!! - The effective computation of the new soil temperatures is performed in thermosoil_profile. 
!!
!! - The calling sequence of thermosoil_main is summarized in the flowchart below.
!! - Thermosoil_init and thermosoil_var_init initialize the variables from
!! restart files or with default values; they also set up
!! the vertical discretization for the numerical scheme.
!! - thermosoil_coef calculates the coefficients for the numerical scheme for the very first iteration of thermosoil;
!! after that, thermosoil_coef is called only at the end of the module to calculate the coefficients for the next timestep.
!! - thermosoil_profile solves the numerical scheme.\n
!!
!! - Flags : one unique flag : THERMOSOIL_TPRO (to be set to the desired initial soil in-depth temperature in K; by default 280K)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): vertically discretized soil temperatures ptn, soil
!! thermal properties (pcapa, pkappa), apparent surface heat capacity (soilcap)
!! and heat flux (soilflux) to be used in enerbil at the next timestep to solve
!! the surface energy balance.
!!
!! REFERENCE(S) : 
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!!  Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin' s PhD thesis relative to the thermal
!!  integration scheme has been scanned and is provided along with the documentation, with name : 
!!  Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{thermosoil_flowchart.png}
!! \endlatexonly
!! 
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, indexgrnd, &
!Isa
       & indexnbdl, control_in, &
       ! GK
        & temp_sol_new, snow, soilcap, soilflx,  &
        & shumdiag_perma, stempdiag, ptnlev1, rest_id, hist_id, hist2_id)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id,hist_id  !! Restart_ file and history file identifier 
                                                                              !! (unitless)
    INTEGER(i_std),INTENT (in)                            :: hist2_id         !! history file 2 identifier (unitless)
    REAL(r_std), INTENT (in)                              :: dtradia          !! model iteration time step in seconds (s)
    LOGICAL, INTENT(in)                                   :: ldrestart_read   !! Logical for restart files to be read 
                                                                              !! (true/false)
    LOGICAL, INTENT(in)                                   :: ldrestart_write  !! Logical for restart files to be writen 
                                                                              !! (true/false)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: index            !! Indeces of the points on the map (unitless)
    INTEGER(i_std),DIMENSION (kjpindex*ngrnd), INTENT (in):: indexgrnd        !! Indeces of the points on the 3D map (vertical 
                                                                              !! dimension towards the ground) (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new     !! Surface temperature at the present time-step,
                                                                              !! Ts @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass @tex ($kg$) @endtex.
                                                                              !! Caveat: when there is snow on the
                                                                              !! ground, the snow is integrated into the soil for
                                                                              !! the calculation of the thermal dynamics. It means
                                                                              !! that the uppermost soil layers can completely or 
                                                                              !! partially consist in snow. In the second case, zx1
                                                                              !! and zx2 are the fraction of the soil layer 
                                                                              !! consisting in snow and 'normal' soil, respectively
                                                                              !! This is calculated in thermosoil_coef.
!    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)    :: shumdiag         !! Relative soil humidity on the diagnostic axis 
                                                                              !! (0-1, unitless). Caveats: when "hydrol" 
                                                                              !! (the 11-layers hydrology) 
                                                                              !! is used, this humidity is 
                                                                              !! calculated with respect to the wilting point: 
                                                                              !! shumdiag= (mc-mcw)/(mcs-mcw), with mc : moisture 
                                                                              !! content; mcs : saturated soil moisture content; 
                                                                              !! mcw: soil moisture content at the wilting point. 
                                                                              !! When the 2-layers hydrology "hydrolc" is used, 
                                                                              !! shumdiag is just a soil wetness index, from 0 to 1
                                                                              !! but cannot direcly be linked to a soil moisture 
                                                                              !! content.
    ! Isa:                                                                          
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in) :: shumdiag_perma         !! Diagnostic of relative humidity
        ! ça remplace shumdiag

    ! specific Isa:
    INTEGER(i_std),DIMENSION (kjpindex*nbdl), INTENT (in) :: indexnbdl       !! Indeces of the points on the 3D map
    ! CR: peut-être à enlever, je ne suis pas sure de bien comprendre la gestion
    ! des flags
    TYPE(control_type), INTENT (in)                    :: control_in       !! Flags that (de)activate parts of the model
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilcap          !! apparent surface heat capacity
                                                                              !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilflx          !! apparent soil heat flux @tex ($W m^{-2}$) @endtex
                                                                              !! , positive 
                                                                              !! towards the soil, writen as Qg (ground heat flux) 
                                                                              !! in the history files, and computed at the end of 
                                                                              !! thermosoil for the calculation of Ts in enerbil, 
                                                                              !! see EQ3.
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (inout) :: stempdiag        !! diagnostic temperature profile @tex ($K$) @endtex
                                                                              !! , eg on the 
                                                                              !! diagnostic axis (levels:1:nbdl). The soil 
                                                                              !! temperature is put on this diagnostic axis to be
                                                                              !! used by other modules (slowproc.f90; routing.f90;
                                                                              !! hydrol or hydrolc when a frozen soil 
                                                                              !! parametrization is used..)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: ptnlev1           !! 1st level soil temperature   

    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex,ngrnd)                :: temp             !! buffer
    REAL(r_std),DIMENSION (kjpindex,ngrnd-1)              :: temp1            !! buffer
    REAL(r_std),DIMENSION (kjpindex)                      :: temp2            !! buffer
    CHARACTER(LEN=80)                                     :: var_name         !! To store variables names for I/O


!_ ================================================================================================================================

  !! 1. do initialisation

  IF (l_first_thermosoil) THEN

        IF (long_print) WRITE (numout,*) ' l_first_thermosoil : call thermosoil_init '

        
        !! 1.1. Allocate and initialize soil temperatures variables
        !! by reading restart files or using default values.
        CALL thermosoil_init (kjit, ldrestart_read, kjpindex, index, rest_id)

       
        !! 1.2.Computes physical constants and arrays; initializes soil thermal properties; produces the first stempdiag
        !!  Computes some physical constants and arrays depending on the soil vertical discretization 
        !! (lskin, cstgrnd, zz, zz_coef, dz1, dz2); get the vertical humidity onto the thermal levels, and 
        !! initializes soil thermal properties (pkappa, pcapa); produces the first temperature diagnostic stempdiag.
        
        CALL thermosoil_var_init (kjpindex, zz, zz_coef, dz1, dz2, pkappa, pcapa, pcapa_en, &
        &        shumdiag_perma, stempdiag, snow, &
!Isa
        &        profil_froz)


        !
        !! 1.3. Computes cgrd, dgrd, soilflx and soilcap coefficients from restart values or initialisation values.
        ! computes cgrd and dgrd coefficient from previous time step (restart)
        !
        CALL thermosoil_coef (kjpindex, dtradia, temp_sol_new, snow, ptn, soilcap, soilflx, zz, dz1, dz2, z1, zdz1,&
           & zdz2, cgrnd, dgrnd, pcapa, pcapa_en, pkappa,&
!Isa
           & profil_froz, pcappa_supp, stempdiag)
        
        
        !! 1.4. call to thermosoil_energy, if you wish to perform some checks (?)
        !!?? the usefulness of this routine seems questionable.
        CALL thermosoil_energy (kjpindex, temp_sol_new, soilcap, .TRUE.) 

    
        ! ajout CR from trunk
        !! 1.5. read restart files for other variables than ptn.
        !!?? mind the use of ok_var here.
        !!?? ok_var is a function of sechiba_io_p.f90, documented as follows :
        !!!! pour déclancher les restarts rajoutés avec un paramètre externe
           !!FUNCTION ok_var ( varname )
           !!CHARACTER(LEN=*), INTENT(IN) :: varname
           !!LOGICAL ok_var
           !!ok_var=.FALSE.
           !!CALL getin_p(varname, ok_var)
           !!END FUNCTION ok_var
         !!
         !! from what we understand, it looks for the chain varname in
         !!run.def; if absent, returns .FALSE., and the variable named
         !!'varname' is not searched for in the restart. This looks like a
         !!trick to read variables in restart files when they are not read
         !!there by default. For all variables in the following sequence, ok_var
         !!is by default false, so don' t bother about this.
         !! this is also logical as those variables have been initialized
         !!above.
         !!?? so maybe this part of the code could be deleted to add clarity.
        
        
!    write(*,*) 'thermosoil 223: ok_converge_isaorig=',control%ok_converge_isaorig 
!    if (control%ok_converge_isaorig) then
!             ! il n'y avait pas d'init???
!         else ! if (ok_converge_isaorig) then
         IF (ldrestart_read) THEN
           IF (long_print) WRITE (numout,*) ' we have to READ a restart file for THERMOSOIL variables'

           var_name= 'cgrnd'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','Cgrnd coefficient.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, .TRUE., temp1, "gather", nbp_glo, index_g)
              IF (MINVAL(temp1) < MAXVAL(temp1) .OR. MAXVAL(temp1) .NE. val_exp) THEN
                 cgrnd(:,:)=temp1(:,:)
              ENDIF
           ENDIF

           var_name= 'dgrnd'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','Dgrnd coefficient.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, .TRUE., temp1, "gather", nbp_glo, index_g)
              IF (MINVAL(temp1) < MAXVAL(temp1) .OR. MAXVAL(temp1) .NE. val_exp) THEN
                 dgrnd(:,:)=temp1(:,:)
              ENDIF
           ENDIF

           var_name= 'z1'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp2, "gather", nbp_glo, index_g)
              IF (MINVAL(temp2) < MAXVAL(temp2) .OR. MAXVAL(temp2) .NE. val_exp) THEN
                 z1(:)=temp2(:)
              ENDIF
           ENDIF

           var_name= 'pcapa'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 pcapa(:,:)=temp(:,:)
              ENDIF
           ENDIF

           var_name= 'pcapa_en'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 pcapa_en(:,:)=temp(:,:)
              ENDIF
           ENDIF

           var_name= 'pkappa'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 pkappa(:,:)=temp(:,:)
              ENDIF
           ENDIF

           var_name= 'zdz1'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, .TRUE., temp1, "gather", nbp_glo, index_g)
              IF (MINVAL(temp1) < MAXVAL(temp1) .OR. MAXVAL(temp1) .NE. val_exp) THEN
                 zdz1(:,:)=temp1(:,:)
              ENDIF
           ENDIF

           var_name= 'zdz2'
           CALL ioconf_setatt('UNITS', '-')
           CALL ioconf_setatt('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 zdz2(:,:)=temp(:,:)
              ENDIF
           ENDIF

           var_name='temp_sol_beg'
           CALL ioconf_setatt('UNITS', 'K')
           CALL ioconf_setatt('LONG_NAME','Old Surface temperature')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp2, "gather", nbp_glo, index_g)
              IF (MINVAL(temp2) < MAXVAL(temp2) .OR. MAXVAL(temp2) .NE. val_exp) THEN
                 temp_sol_beg(:) = temp2(:)
              ENDIF
           ENDIF

        ENDIF !ldrestart_read
!        endif !if (OK_restart_pkappa) then
!        ! end CR from trunk

        RETURN 

    ENDIF !l_first_thermosoil

    
  !! 2. Prepares the restart files for the next simulation

    !!?? do all the coefficients (cgrnd, dgrnd...) be put in the restart file
    !! as they are by default not read there, but calculated in
    !!thermosoil_var_init from the restart or initial temperature ?
    !! exceptions are soilcap and soilflx, used in enerbil, and of course ptn.
    IF (ldrestart_write) THEN

        IF (long_print) WRITE (numout,*) ' we have to complete restart file with THERMOSOIL variables'

        var_name= 'ptn'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, ptn, 'scatter', nbp_glo, index_g)

        if (ok_wetdiaglong) then
	var_name= 'wetdiaglong'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, wetdiaglong, 'scatter', nbp_glo, index_g)
        endif !if (ok_wetdiaglong) then
        var_name= 'E_sol_lat_couche'
        CALL restput_p(rest_id, var_name, nbp_glo,ngrnd , 1, kjit, E_sol_lat_couche, 'scatter', nbp_glo, index_g)

         ! CR: ajout from trunc
        var_name= 'cgrnd'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, cgrnd, 'scatter', nbp_glo, index_g)
        var_name= 'dgrnd'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, dgrnd, 'scatter', nbp_glo, index_g)

        var_name= 'z1'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, z1, 'scatter', nbp_glo, index_g)

        var_name= 'pcapa'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, pcapa, 'scatter', nbp_glo, index_g)

        var_name= 'pcapa_en'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, pcapa_en, 'scatter', nbp_glo, index_g)

        var_name= 'pkappa'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, pkappa, 'scatter', nbp_glo, index_g)

        var_name= 'zdz1'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, zdz1, 'scatter', nbp_glo, index_g)

        var_name= 'zdz2'
        CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, zdz2, 'scatter', nbp_glo, index_g)

        var_name= 'temp_sol_beg'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, temp_sol_beg, 'scatter', nbp_glo, index_g)

        var_name= 'soilcap'  
        CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  soilcap, 'scatter',  nbp_glo, index_g)
        
        var_name= 'soilflx'  
        CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  soilflx, 'scatter',  nbp_glo, index_g)

        ! read in enerbil
        var_name= 'temp_sol_new'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, temp_sol_new, 'scatter', nbp_glo, index_g)

        RETURN 

    END IF !ldrestart_write


  !! 3. Put the soil wetness diagnostic on the levels of the soil temperature

    !!?? this could logically be put just before the last call to
    !!thermosoil_coef, as the results are used there...    !
    CALL thermosoil_humlev(kjpindex, shumdiag_perma, snow)
    
    ! Compute long-term soil humidity (for permafrost)
    !    
    if (ok_wetdiaglong) then
        CALL thermosoil_wlupdate( kjpindex, dtradia, ptn, wetdiag, wetdiaglong )
    else
        wetdiaglong(:,:)=wetdiag(:,:)
    endif

    !  !! 4. Effective computation of the soil temperatures profile, using the cgrd and dgrd coefficients from previsou tstep
    ! computes profile with previous cgrd and dgrd coefficient
    !
!    if (control%ok_converge_isaorig) then
!    if (1.eq.0) then
!       CALL thermosoil_profile_isaorig (kjpindex, temp_sol_new, ptn)
!    else
       CALL thermosoil_profile (kjpindex, temp_sol_new, ptn, stempdiag)
!    endif
! on appelle thermosoil_profile comme dans le trunk, car ça ne change rien

  !! 5. Call to thermosoil_energy, still to be clarified..

    CALL thermosoil_energy (kjpindex, temp_sol_new, soilcap, .FALSE.)
    
  !! 6. Writing the history files according to the ALMA standards (or not..)

    !in only one file (hist2_id <=0) or in 2 different files (hist2_id >0).
    IF ( .NOT. almaoutput ) THEN
      CALL histwrite(hist_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      
      IF ( control_in%hydrol_cwrr ) THEN
      CALL histwrite(hist_id, 'ptn_beg', kjit, ptn_beg, kjpindex*ngrnd, indexgrnd)
!Isa
      CALL histwrite(hist_id, 'profil_froz', kjit, profil_froz, kjpindex*ngrnd, indexgrnd)
      CALL histwrite(hist_id, 'pkappa', kjit, pkappa, kjpindex*ngrnd, indexgrnd)
      CALL histwrite(hist_id, 'pcapa', kjit, pcapa, kjpindex*ngrnd, indexgrnd)
       CALL histwrite(hist_id, 'pcappa_supp', kjit, pcappa_supp, kjpindex*ngrnd, indexgrnd)
      CALL histwrite(hist_id, 'wetdiag', kjit, wetdiag(:,:), kjpindex*ngrnd, indexgrnd)
      CALL histwrite(hist_id, 'wetdiaglong', kjit, wetdiaglong(:,:), kjpindex*ngrnd, indexgrnd)
      CALL histwrite(hist_id, 'shumdiag_perma', kjit, shumdiag_perma, kjpindex*nbdl, indexnbdl)
      CALL histwrite(hist_id, 'stempdiag', kjit, stempdiag, kjpindex*nbdl, indexnbdl)
!Isa
      endif !IF ( control_flags%hydrol_cwrr ) THEN
      !write(*,*) 'thermosoil 225: write Qg, not almaoutput'
      CALL histwrite(hist_id, 'Qg', kjit, soilflx, kjpindex, index)

    ELSE !IF ( .NOT. almaoutput ) THEN
      CALL histwrite(hist_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      CALL histwrite(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite(hist_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
      CALL histwrite(hist_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
    ENDIF  !IF ( .NOT. almaoutput ) THEN
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite(hist2_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
       ELSE
          CALL histwrite(hist2_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
          CALL histwrite(hist2_id, 'Qg', kjit, soilflx, kjpindex, index)
          CALL histwrite(hist2_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
          CALL histwrite(hist2_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
       ENDIF
    ENDIF
    
  !! 7. A last final call to thermosoil_coef
 
    !! A last final call to thermosoil_coef, which calculates the different
    !!coefficients (cgrnd, dgrnd, dz1, z1, zdz2, soilcap, soilflx) from this time step to be
    !!used at the next time step, either in the surface temperature calculation
    !!(soilcap, soilflx) or in the soil thermal numerical scheme.
    !
    !Isa
    CALL thermosoil_coef (kjpindex, dtradia, temp_sol_new, snow, ptn, soilcap, soilflx, zz,&
           & dz1, dz2, z1, zdz1,zdz2, cgrnd, dgrnd, pcapa, pcapa_en, pkappa, profil_froz, pcappa_supp, stempdiag)

    ! CR: ajout du trunk
    ptnlev1(:) = ptn(:,1)

 !   ptn_beg(:,:)      = ptn(:,:) ! seulement dans isa: 
 ! pas necessaire dans trunk car deja fait dans thermosoil_energy
    IF (long_print) WRITE (numout,*) ' thermosoil_main done '

  END SUBROUTINE thermosoil_main

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_init
!!
!>\BRIEF        Allocates local and global arrays; initializes soil temperatures using either restart files
!! or a fixed value set by the flag THERMOSOIL_TPRO.
!!		  
!! DESCRIPTION  : flag : THERMOSOIL_TPRO (to be set to the desired initial temperature in K; by default 280K).
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_init(kjit, ldrestart_read, kjpindex, index, rest_id)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number (unitless) 
    LOGICAL,INTENT (in)                                 :: ldrestart_read     !! Logical for restart file to read (true/false)
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map (unitless)
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! Restart file identifier (unitless)
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ier, i
    CHARACTER(LEN=80)                                  :: var_name            !! To store variables names for I/O

!_ ================================================================================================================================

  !! 1. Initialisation

    !! Initialisation has to be done only one time, so the logical
    !! logical l_first_thermosoil has to be set to .FALSE. now..
    IF (l_first_thermosoil) THEN 
        l_first_thermosoil=.FALSE.
    ELSE 
        WRITE (numout,*) ' l_first_thermosoil false . we stop '
        STOP 'thermosoil_init'
    ENDIF
! Isa ++
    if (control%ok_freeze) then
        ! on active toutes les options d'Isabelle:
        write(*,*) 'thermosoil 487: ok_freeze active, on active toutes les options d Isabelle'
        !CL mettre ok_ecorr a yes pour diminuer les erreurs du à la fonte de la neige
        ok_Ecorr = .true.
        read_permafrost_map = .false.
        read_reftemp = .false.
        ok_freeze_thermix=.true.
        ok_thermix_trunc=.false.
        ok_wetdiaglong = .true.        
    else !if (ok_freeze) then
        write(*,*) 'thermosoil 495: ok_freeze desactive, on lit les options 1 par 1'
        ! on active les options comme on veut
        ok_Ecorr = .FALSE.
        CALL getin_p ('OK_ECORR',ok_Ecorr)
        if (ok_Ecorr) write(*,*)'ISA : ok_Ecorr!'
        read_permafrost_map = .false.
        CALL getin_p ('READ_PERMAFROST_MAP',read_permafrost_map)
        if (read_permafrost_map) write(*,*)'ISA : read_permafrost_map !'
        read_reftemp = .false.
        CALL getin_p ('READ_REFTEMP',read_reftemp)
        if (read_reftemp) write(*,*)'ISA : read_reftemp!'
        ok_freeze_thermix = .FALSE.
        CALL getin_p ('OK_FREEZE_THERMIX',ok_freeze_thermix)
        if (ok_freeze_thermix) write(*,*)'ISA : ok_freeze_thermix !'
        ! Isa --
        ! CR: option en plus pour converger avec le trunk:
        if (ok_freeze_thermix) then
        else
                ok_thermix_trunc = .TRUE.
                CALL getin_p ('OK_THERMIX_trunc',ok_thermix_trunc)
                write(*,*) 'ok_thermix_trunc=',ok_thermix_trunc
        endif
        ok_wetdiaglong = .FALSE.
        CALL getin_p ('OK_WETDIAGLONG',ok_wetdiaglong)
    endif !!if (ok_freeze) then
!    ok_converge_isaorig=.TRUE.
!    CALL getin_p ('ok_converge_isaorig',ok_converge_isaorig)

  !! 2. Arrays allocations

    ALLOCATE (ptn(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ptn allocation. We stop. We need ',kjpindex,' fois ',ngrnd,' words = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (z1(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in z1 allocation. We STOP. We need ',kjpindex,' words '
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (cgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in cgrnd allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*(ngrnd-1)
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (dgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in dgrnd allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*(ngrnd-1)
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (pcapa(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (pkappa(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pkappa allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (zdz1(kjpindex,ngrnd-1),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zdz1 allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*(ngrnd-1)
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (zdz2(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zdz2 allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (surfheat_incr(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in surfheat_incr allocation. We STOP. We need ',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (coldcont_incr(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in coldcont_incr allocation. We STOP. We need ',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (pcapa_en(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa_en allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (ptn_beg(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ptn_beg allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (temp_sol_beg(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in temp_sol_beg allocation. We STOP. We need ',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (wetdiag(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in wetdiag allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    wetdiag(:,:)=val_exp
    ALLOCATE (wetdiaglong(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in wetdiaglong allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

!Isa ++
    ALLOCATE (profil_froz(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa_en allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    
     ALLOCATE (pcappa_supp(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa_supp allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF
    
     ALLOCATE (E_sol_lat_couche(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in E_sol_sat allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (permafrost(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in permafrost allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (excess_ice(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in excess_ice allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (overburden(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in overburden allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (reftemp(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in reftemp allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    !Isa --
    !
  !! 3. Reads restart files for soil temperatures only 
    
    !! Reads restart files for soil temperatures only. If no restart file is
    !! found,  the initial soil temperature is by default set to 280K at all depths. The user
    !! can decide to initialize soil temperatures at an other value, in which case he should set the flag THERMOSOIL_TPRO
    !! to this specific value in the run.def.

    IF (ldrestart_read) THEN
        IF (long_print) WRITE (numout,*) ' we have to READ a restart file for THERMOSOIL variables'

        var_name= 'ptn'
        CALL ioconf_setatt('UNITS', 'K')
        CALL ioconf_setatt('LONG_NAME','Soil Temperature profile')
        CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., ptn, "gather", nbp_glo, index_g)
        !
        !Config Key   = THERMOSOIL_TPRO
        !Config Desc  = Initial soil temperature profile if not found in restart
        !Config Def   = 280.
        !Config If    = OK_SECHIBA
        !Config Help  = The initial value of the temperature profile in the soil if 
        !Config         its value is not found in the restart file. This should only 
        !Config         be used if the model is started without a restart file. Here
        !Config         we only require one value as we will assume a constant 
        !Config         throughout the column.
        !Config Units = Kelvin [K]
!
	if (read_permafrost_map) then 
		CALL read_permafrostmap(kjpindex,lalo,overburden,excess_ice,permafrost)
		DO i=1, ngrnd
		  where (ptn(:,i).eq.val_exp.AND.permafrost(:).eq.1) 
			ptn(:,i)=272._r_std
		  elsewhere 
			ptn(:,i)=280._r_std
		  endwhere
		ENDDO
	else if (read_reftemp) then
		CALL read_reftempfile(kjpindex,lalo,reftemp)
       		CALL setvar_p (ptn, val_exp, 'NO_KEYWORD' ,reftemp)
	else
       		CALL setvar_p (ptn, val_exp,'THERMOSOIL_TPRO',280._r_std)
	endif
 
        ptn_beg(:,:) = ptn(:,:)
     
	

 !      if (ok_wetdiaglong) then
 ! on initialise toujours, ça ne peut pas faire de mal et ca permet convergence
 ! avec isaorig
       var_name= 'wetdiaglong'
       CALL ioconf_setatt('UNITS', '-')
       CALL ioconf_setatt('LONG_NAME','Long-term soil humidity')
       CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., wetdiaglong, "gather", nbp_glo, index_g)

        CALL restget_p (rest_id, 'E_sol_lat_couche', nbp_glo, ngrnd, 1, kjit, .TRUE., E_sol_lat_couche, "gather", nbp_glo, index_g)
        CALL setvar_p (E_sol_lat_couche, val_exp,'NO_KEYWORD',zero)
!        endif !if (ok_wetdiaglong) then


   ENDIF


    IF (long_print) WRITE (numout,*) ' thermosoil_init done '


  END SUBROUTINE thermosoil_init

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_clear
!!
!>\BRIEF        Sets the flag l_first_thermosoil to true and desallocates the allocated arrays.
!! ??!! the call of thermosoil_clear originates from sechiba_clear but the calling sequence and 
!! its purpose require further investigation.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
 SUBROUTINE thermosoil_clear()

        l_first_thermosoil=.TRUE.
 
        IF ( ALLOCATED (ptn)) DEALLOCATE (ptn)
        IF ( ALLOCATED (z1)) DEALLOCATE (z1) 
        IF ( ALLOCATED (cgrnd)) DEALLOCATE (cgrnd) 
        IF ( ALLOCATED (dgrnd)) DEALLOCATE (dgrnd) 
        IF ( ALLOCATED (pcapa)) DEALLOCATE (pcapa)
        IF ( ALLOCATED (pkappa))  DEALLOCATE (pkappa)
        IF ( ALLOCATED (zdz1)) DEALLOCATE (zdz1)
        IF ( ALLOCATED (zdz2)) DEALLOCATE (zdz2)
        IF ( ALLOCATED (pcapa_en)) DEALLOCATE (pcapa_en)
        IF ( ALLOCATED (ptn_beg)) DEALLOCATE (ptn_beg)
        IF ( ALLOCATED (temp_sol_beg)) DEALLOCATE (temp_sol_beg)
        IF ( ALLOCATED (surfheat_incr)) DEALLOCATE (surfheat_incr)
        IF ( ALLOCATED (coldcont_incr)) DEALLOCATE (coldcont_incr)
        IF ( ALLOCATED (wetdiag)) DEALLOCATE (wetdiag)
!Isa
        IF ( ALLOCATED (profil_froz)) DEALLOCATE (profil_froz)
        IF ( ALLOCATED (wetdiaglong)) DEALLOCATE (wetdiaglong)

  END SUBROUTINE thermosoil_clear
  !!
  !!
  FUNCTION fz(rk) RESULT (fz_result)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    REAL(r_std), INTENT(in)                        :: rk
    
    !! 0.2 Output variables

    REAL(r_std)                                    :: fz_result
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    fz_result = fz1 * (zalph ** rk - un) / (zalph - un)

  END FUNCTION fz


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_var_init
!!
!>\BRIEF        Define and initializes the soil thermal parameters
!!		  
!! DESCRIPTION	: This routine\n
!! 1. Defines the parameters ruling the vertical grid of the thermal scheme (fz1, zalpha).\n
!! 2. Defines the scaling coefficients for adimensional depths (lskin, cstgrnd, see explanation in the 
!!    variables description of thermosoil_main). \n
!! 3. Calculates the vertical discretization of the soil (zz, zz_coef, dz2) and the constants used
!!    in the numerical scheme and which depend only on the discretization (dz1, lambda).\n
!! 4. Initializes the soil thermal parameters (capacity, conductivity) based on initial soil moisture content.\n
!! 5. Produces a first temperature diagnostic based on temperature initialization.\n
!!
!! The scheme comprizes ngrnd=7 layers by default.
!! The layer' s boundaries depths (zz_coef) follow a geometric series of ratio zalph=2 and first term fz1.\n
!! zz_coef(jg)=fz1.(1-zalph^jg)/(1-zalph) \n
!! The layers' boudaries depths are first calculated 'adimensionally', ie with a
!! discretization adapted to EQ5. This discretization is chosen for its ability at
!! reproducing a thermal signal with periods ranging from days to centuries. (see
!! Hourdin, 1992). Typically, fz1 is chosen as : fz1=0.3*cstgrnd (with cstgrnd the
!! adimensional attenuation depth). \n
!! The factor lskin/cstgrnd is then used to go from adimensional depths to
!! depths in m.\n
!! zz(real)=lskin/cstgrnd*zz(adimensional)\n
!! Similarly, the depths of the numerical nodes are first calculated
!! adimensionally, then the conversion factor is applied.\n
!! the numerical nodes (zz) are not exactly the layers' centers : their depths are calculated as follows:\n
!! zz(jg)=fz1.(1-zalph^(jg-1/2))/(1-zalph)\n
!! The values of zz and zz_coef used in the default thermal discretization are in the following table.
!! \latexonly
!! \includegraphics{thermosoil_var_init1.jpg}
!! \endlatexonly\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCE(S)	:
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of 
!! planetary atmospheres, Ph.D. thesis, Paris VII University.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_var_init(kjpindex, zz, zz_coef, dz1, dz2, pkappa, pcapa, pcapa_en, &
  &     shumdiag_perma, stempdiag, snow, &
!Isa
  &     profil_froz)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size (unitless)
    REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (in)      :: shumdiag_perma    !! Relative soil humidity on the diagnostic axis 
                                                                                  !! (unitless), [0,1]. (see description of the 
                                                                                  !! variables of thermosoil_main for more 
                                                                                  !! explanations) 
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)     	     :: snow              !! Snow quantity
    
    !! 0.2 Output variables

    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: zz                !! depths of the layers'numerical nodes 
                                                                                  !! @tex ($m$)@endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: zz_coef		  !! depths of the layers'boundaries 
                                                                                  !! @tex ($m$)@endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: dz1               !! numerical constant depending on the vertical
                                                                                  !! thermal grid only @tex  ($m^{-1}$) @endtex. 
                                                                                  !! (see description
                                                                                  !! of the variables of thermosoil_main for more
                                                                                  !! explanations)
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: dz2               !! thicknesses of the soil thermal layers 
                                                                                  !! @tex ($m$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(out)     :: pcapa             !! volumetric vertically discretized soil heat 
                                                                                  !! capacity @tex ($J K^{-1} m^{-3}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(out)     :: pcapa_en	  !! volumetric vertically discretized heat 
                                                                                  !! capacity used in thermosoil_energy
                                                                                  !! @tex ($J K^{-1} m^{-3}$) @endtex ;
                                                                                  !! usefulness still to be clarified.
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(out)     :: pkappa            !! vertically discretized soil thermal 
                                                                                  !! conductivity @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (out)     :: stempdiag         !! Diagnostic temperature profile @tex ($K$)
                                                                                  !! @endtex                                                                                  
!Isa++
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(out)     :: profil_froz   

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ier, ji, jg
    REAL(r_std)                                              :: sum
!_ ================================================================================================================================

  !! 1. Initialization of the parameters of the vertical discretization and of the attenuation depths
    
    cstgrnd=SQRT(one_day / pi)
    lskin = SQRT(so_cond / so_capa * one_day / pi)
    fz1 = 0.3_r_std * cstgrnd
    zalph = deux
    
  !! 2.  Computing the depth of the thermal levels (numerical nodes) and the layers boundaries
   
    !! Computing the depth of the thermal levels (numerical nodes) and 
    !! the layers boundariesusing the so-called
    !! adimentional variable z' = z/lskin*cstgrnd (with z in m)
    
    !! 2.1 adimensional thicknesses of the layers
    DO jg=1,ngrnd

    !!?? code simplification hopefully possible here with up-to-date compilers !
    !!! This needs to be solved soon. Either we allow CPP options in SECHIBA or the VPP
    !!! fixes its compiler 
    !!!#ifdef VPP5000
      dz2(jg) = fz(REAL(jg,r_std)-undemi+undemi) - fz(REAL(jg-1,r_std)-undemi+undemi)
    !!!#else
    !!!      dz2(jg) = fz(REAL(jg,r_std)) - fz(REAL(jg-1,r_std))
    !!!#endif
    ENDDO
    
    !! 2.2 adimentional depth of the numerical nodes and layers' boudaries
    DO jg=1,ngrnd
      zz(jg)      = fz(REAL(jg,r_std) - undemi)
      zz_coef(jg) = fz(REAL(jg,r_std)-undemi+undemi)
    ENDDO

    !! 2.3 Converting to meters
    DO jg=1,ngrnd
      zz(jg)      = zz(jg) /  cstgrnd * lskin
      zz_coef(jg) = zz_coef(jg) / cstgrnd * lskin 
      dz2(jg)     = dz2(jg) /  cstgrnd * lskin
!       write(*,*) 'Chloe ngrnd,zz',jg,zz(jg)
    ENDDO

    !! 2.4 Computing some usefull constants for the numerical scheme
    DO jg=1,ngrnd-1
      dz1(jg)  = un / (zz(jg+1) - zz(jg))
    ENDDO
    lambda = zz(1) * dz1(1)
    !! 2.5 Get the wetness profile on the thermal vertical grid from the diagnostic axis
    CALL thermosoil_humlev(kjpindex, shumdiag_perma, snow)
    !
    ! Compute long-term soil humidity (for permafrost)
    CALL setvar_p (wetdiaglong, val_exp,'NO_KEYWORD',wetdiag(:,:))
    ! cette routine veut dire que wetdiaglong=wetdiag si wetdiaglong=val_exp

    !! 2.6 Thermal conductivity at all levels
!Isa

    write(*,*) 'ok_freeze_thermix=',ok_freeze_thermix
if (ok_freeze_thermix) then
	CALL thermosoil_getdiff( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en,profil_froz, pcappa_supp)
else !if (ok_freeze_thermix) then
    if (ok_thermix_trunc) then
        ! pour convergence avec le trunc
        CALL thermosoil_getdiff_old_thermix_trunc2( kjpindex, pkappa, pcapa, pcapa_en )
    else !if (ok_thermix_trunc) then
	CALL thermosoil_getdiff_old_thermix_with_snow( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en,profil_froz )
    endif !if (ok_thermix_trunc) then
endif !if (ok_freeze_thermix) then

  !! 3. Diagnostics : consistency checks on the vertical grid.
    WRITE (numout,*) 'diagnostic des niveaux dans le sol' !!?? to be changed,
    WRITE (numout,*) 'niveaux intermediaires et pleins'
    sum = zero
    DO jg=1,ngrnd
      sum = sum + dz2(jg)
      WRITE (numout,*) zz(jg),sum
    ENDDO

  !! 4. Compute a first diagnostic temperature profile

    CALL thermosoil_diaglev(kjpindex, stempdiag)

    IF (long_print) WRITE (numout,*) ' thermosoil_var_init done '

  END SUBROUTINE thermosoil_var_init


  

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_coef
!!
!>\BRIEF        Calculate soil thermal properties, integration coefficients, apparent heat flux,
!! surface heat capacity,  
!!
!! DESCRIPTION	: This routine computes : \n
!!		1. the soil thermal properties. \n 
!!		2. the integration coefficients of the thermal numerical scheme, cgrnd and dgrnd,
!!              which depend on the vertical grid and on soil properties, and are used at the next 
!!              timestep.\n
!!              3. the soil apparent heat flux and surface heat capacity soilflux
!!              and soilcap, used by enerbil to compute the surface temperature at the next
!!              timestep.\n
!!             -  The soil thermal properties depend on water content (wetdiag) and on the presence 
!!              of snow : snow is integrated into the soil for the thermal calculations, ie if there 
!!              is snow on the ground, the first thermal layer(s) consist in snow, depending on the 
!!              snow-depth. If a layer consists out of snow and soil, wheighed soil properties are 
!!              calculated\n
!!             - The coefficients cgrnd and dgrnd are the integration
!!              coefficients for the thermal scheme \n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!              They correspond respectively to $\beta$ and $\alpha$ from F. Hourdin\'s thesis and 
!!              their expression can be found in this document (eq A19 and A20)
!!             - soilcap and soilflux are the apparent surface heat capacity and flux
!!               used in enerbil at the next timestep to solve the surface
!!               balance for Ts (EQ3); they correspond to $C_s$ and $F_s$ in F.
!!               Hourdin\'s PhD thesis and are expressed in eq. A30 and A31. \n
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflux+otherfluxes(Ts(t)) \n
!!                                      -- EQ3 --\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cgrnd, dgrnd, pcapa, pkappa, soilcap, soilflx
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_coef (kjpindex, dtradia, temp_sol_new, snow, ptn, soilcap, soilflx, zz, dz1, dz2, z1, zdz1,&
           & zdz2, cgrnd, dgrnd, pcapa, pcapa_en, pkappa, &
!Isa
	& profil_froz, pcappa_supp, stempdiag)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: kjpindex     !! Domain size (unitless)
    REAL(r_std), INTENT (in)                               :: dtradia      !! Time step in seconds @tex ($s$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new !! soil surface temperature @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: snow         !! snow mass @tex ($Kg$) @endtex
    REAL(r_std), DIMENSION (ngrnd), INTENT(in)             :: zz           !! depths of the soil thermal numerical nodes 
                                                                           !! @tex ($m$) @endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(in)             :: dz1          !! numerical constant depending on the vertical 
                                                                           !! thermal grid only @tex ($m^{-1}$) @endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(in)             :: dz2          !! thicknesses of the soil thermal layers
                                                                           !! @tex ($m$) @endtex 
    ! Isa: ptn devient inout alors qu'il était in dans le trunk                                                                       
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT (inout)   :: ptn          !! vertically discretized soil temperatures  
                                                                           !! @tex ($K$) @endtex
    
    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilcap      !! surface heat capacity
                                                                           !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilflx      !! surface heat flux @tex ($W m^{-2}$) @endtex,
                                                                           !! positive towards the 
                                                                           !! soil, writen as Qg (ground heat flux) in the history 
                                                                           !! files.
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: z1           !! numerical constant @tex ($W m^{-1} K^{-1}$) @endtex

    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: cgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (beta in F. Hourdin thesis)
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: dgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (alpha in F. Hourdin thesis)
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: zdz1         !! numerical (buffer) constant 
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex

    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(out)   :: zdz2         !! numerical (buffer) constant  
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex
!Isa
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(out)     :: profil_froz
!Isa E
     REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcappa_supp
     REAL(r_std),DIMENSION(kjpindex,nbdl),INTENT(out)    :: stempdiag


    !! 0.3 Modified variable

    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(inout) :: pcapa        !! volumetric vertically discretized soil heat capacity
                                                                           !! @tex ($J K^{-1} m^{-3}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(inout) :: pcapa_en     !! volumetric vertically discretized heat capacity used 
                                                                           !! to calculate surfheat_incr
                                                                           !! @tex ($J K^{-1} m^{-3}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(inout) :: pkappa       !! vertically discretized soil thermal conductivity 
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex

    !! 0.4 Local variables

    INTEGER(i_std)                                         :: ji, jg
    REAL(r_std), DIMENSION(kjpindex)                       :: snow_h       !! snow_h is the snow height @tex ($m$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)                       :: zx1, zx2     !! zx1 and zx2 are the layer fraction consisting in snow
                                                                           !! and soil respectively.
!_ ================================================================================================================================

  !! 1. Computation of the soil thermal properties
   
    ! Computation of the soil thermal properties; snow properties are also accounted for


    ! Isa
if (ok_freeze_thermix) then
    CALL thermosoil_getdiff( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en, profil_froz, pcappa_supp)
else !if (ok_freeze_thermix) then
    if (ok_thermix_trunc) then
      CALL thermosoil_getdiff_old_thermix_trunc( kjpindex, snow, pkappa, pcapa, pcapa_en )
    else !if (ok_thermix_trunc) then
      CALL thermosoil_getdiff_old_thermix_with_snow( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en, profil_froz) 
    endif !if (ok_thermix_trunc) then

endif !if (ok_freeze_thermix) then
!    write(*,*) 'thermosoil 980: pkappa=',pkappa(1,1)
!    write(*,*) 'pcapa,pcapa_en=',pcapa(1,1),pcapa_en(1,1)

!Isa
if (ok_Ecorr) then
    CALL thermosoil_readjust(kjpindex, ptn)
endif

!if (control%ok_converge_isaorig) then
!if (1.eq.0) then
!    write(*,*) 'thermosoil 1042: call thermosoil_diaglev'
!    CALL thermosoil_diaglev(kjpindex, stempdiag)
!endif !if (control%ok_converge_isaorig) then
! sinon, on le fait plutôt à la fin de thermosoil_profile
! en fait, ça ne change rien, donc même quand on teste convergence avec isa, on
! appelle thermosoil_diaglev a la fin de thermosoil_profile

    !! 2. computation of the coefficients of the numerical integration scheme

    ! cgrnd, dgrnd

    !! 2.1.  some "buffer" values

    DO jg=1,ngrnd
      DO ji=1,kjpindex
        zdz2(ji,jg)=pcapa(ji,jg) * dz2(jg)/dtradia
      ENDDO
    ENDDO
    
    DO jg=1,ngrnd-1
      DO ji=1,kjpindex
        zdz1(ji,jg) = dz1(jg) * pkappa(ji,jg)
      ENDDO
    ENDDO !DO jg=1,ngrnd-1

    
    !! 2.2.  the coefficients ! 
    DO ji = 1,kjpindex
      z1(ji) = zdz2(ji,ngrnd) + zdz1(ji,ngrnd-1)
      cgrnd(ji,ngrnd-1) = zdz2(ji,ngrnd) * ptn(ji,ngrnd) / z1(ji)
      dgrnd(ji,ngrnd-1) = zdz1(ji,ngrnd-1) / z1(ji)
    ENDDO

    DO jg = ngrnd-1,2,-1
      DO ji = 1,kjpindex
        z1(ji) = un / (zdz2(ji,jg) + zdz1(ji,jg-1) + zdz1(ji,jg) * (un - dgrnd(ji,jg)))
        cgrnd(ji,jg-1) = (ptn(ji,jg) * zdz2(ji,jg) + zdz1(ji,jg) * cgrnd(ji,jg)) * z1(ji)
        dgrnd(ji,jg-1) = zdz1(ji,jg-1) * z1(ji)
      ENDDO
    ENDDO

  !! 3. Computation of the apparent ground heat flux 
    
    !! Computation of the apparent ground heat flux (> towards the soil) and
    !! apparent surface heat capacity, used at the next timestep by enerbil to
    !! compute the surface temperature.
    DO ji = 1,kjpindex
      soilflx(ji) = zdz1(ji,1) * (cgrnd(ji,1) + (dgrnd(ji,1)-1.) * ptn(ji,1))
      soilcap(ji) = (zdz2(ji,1) * dtradia + dtradia * (un - dgrnd(ji,1)) * zdz1(ji,1))
      z1(ji) = lambda * (un - dgrnd(ji,1)) + un
      soilcap(ji) = soilcap(ji) / z1(ji)
      soilflx(ji) = soilflx(ji) + &
         & soilcap(ji) * (ptn(ji,1) * z1(ji) - lambda * cgrnd(ji,1) - temp_sol_new(ji)) / dtradia 
    ENDDO

    IF (long_print) WRITE (numout,*) ' thermosoil_coef done '

  END SUBROUTINE thermosoil_coef

  !! Computation of : the ground temperature evolution
  !!  
  SUBROUTINE thermosoil_profile_isaorig (kjpindex, temp_sol_new, ptn)

    ! interface description
    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex       !! Domain size
    ! input fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new   !! New soil temperature
    ! modified fields
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (inout)   :: ptn            !! Different levels soil temperature
    ! output fields
    ! local declaration
    INTEGER(i_std)                                          :: ji, jg

    !   ------
    !
    !   Method: implicit time integration
    !   -------
    !   Consecutives ground temperatures are related by:
    !           T(k+1) = C(k) + D(k)*T(k)  (1)
    !   the coefficients C and D are computed at the t-dt time-step.
    !   Routine structure:
    !   -new temperatures are computed  using (1)
    !           
    !
    !    surface temperature
    DO ji = 1,kjpindex
      ptn(ji,1) = (lambda * cgrnd(ji,1) + temp_sol_new(ji)) / (lambda * (un - dgrnd(ji,1)) + un)
    ENDDO

    !   other temperatures
    DO jg = 1,ngrnd-1
      DO ji = 1,kjpindex
        ptn(ji,jg+1) = cgrnd(ji,jg) + dgrnd(ji,jg) * ptn(ji,jg)
      ENDDO
    ENDDO


    IF (long_print) WRITE (numout,*) ' thermosoil_profile done '

  END SUBROUTINE thermosoil_profile_isaorig

  !! ================================================================================================================================
!! SUBROUTINE   : thermosoil_profile
!!
!>\BRIEF        In this routine solves the numerical soil thermal scheme, ie calculates the new soil temperature profile; 
!! This profile is then exported onto the diagnostic axis (call thermosoil_diaglev)
!!
!! DESCRIPTION	: The calculation of the new soil temperature profile is based on
!! the cgrnd and dgrnd values from the previous timestep and the surface temperature Ts aka temp_sol_new. (see detailed
!! explanation in the header of the thermosoil module or in the reference).\n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k)\n
!!                                      -- EQ1 --\n
!!                           Ts=(1-lambda)*T(1) -lambda*T(2)\n 
!!                                      -- EQ2--\n
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

 SUBROUTINE thermosoil_profile (kjpindex, temp_sol_new, ptn, stempdiag)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex       !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new   !! Surface temperature at the present time-step 
                                                                               !! @tex ($K$) @endtex
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out)      :: stempdiag      !! diagnostic temperature profile 
                                                                               !! @tex ($K$) @endtex

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (inout)   :: ptn            !! vertically discretized soil temperatures 
                                                                               !! @tex ($K$) @endtex


    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jg
!_ ================================================================================================================================
    
  !! 1. Computes the soil temperatures ptn.

    !! 1.1. ptn(jg=1) using EQ1 and EQ2
    DO ji = 1,kjpindex
      ptn(ji,1) = (lambda * cgrnd(ji,1) + temp_sol_new(ji)) / (lambda * (un - dgrnd(ji,1)) + un)
    ENDDO

    !! 1.2. ptn(jg=2:ngrnd) using EQ1.
    DO jg = 1,ngrnd-1
      DO ji = 1,kjpindex
        ptn(ji,jg+1) = cgrnd(ji,jg) + dgrnd(ji,jg) * ptn(ji,jg)
      ENDDO
    ENDDO

  !! 2. Put the soil temperatures onto the diagnostic axis 
  
    !! Put the soil temperatures onto the diagnostic axis for convenient
    !! use in other routines (stomate..)
    CALL thermosoil_diaglev(kjpindex, stempdiag)

    IF (long_print) WRITE (numout,*) ' thermosoil_profile done '

  END SUBROUTINE thermosoil_profile
!!
!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_diaglev
!!
!>\BRIEF        Interpolation of the soil in-depth temperatures onto the diagnostic profile.
!!
!! DESCRIPTION  : This is a very easy linear interpolation, with intfact(jd, jg) the fraction
!! the thermal layer jg comprised within the diagnostic layer jd. The depths of
!! the diagnostic levels are diaglev(1:nbdl), computed in slowproc.f90.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_diaglev(kjpindex, stempdiag)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                          :: kjpindex       !! Domain size (unitless)
   
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out) :: stempdiag      !! Diagnostic soil temperature profile @tex ($K$) @endtex
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jd, jg
    REAL(r_std)                                         :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)      :: intfact
    LOGICAL, PARAMETER                                  :: check=.FALSE.
!_ ================================================================================================================================
    
  !! 1. Computes intfact(jd, jg)

    !! Computes intfact(jd, jg), the fraction
    !! the thermal layer jg comprised within the diagnostic layer jd.

    IF ( .NOT. ALLOCATED(intfact)) THEN
        
        ALLOCATE(intfact(nbdl, ngrnd))
        
        prev_diag = zero
        DO jd = 1, nbdl
          lev_diag = diaglev(jd)
          prev_prog = zero
          DO jg = 1, ngrnd
             IF ( jg == ngrnd .AND. (prev_prog + dz2(jg)) < lev_diag ) THEN
                lev_prog = lev_diag
             ELSE
                lev_prog = prev_prog + dz2(jg)
             ENDIF
            intfact(jd,jg) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog),&
                        & zero)/(lev_diag-prev_diag)
            prev_prog = lev_prog
          ENDDO
           prev_diag = lev_diag
        ENDDO

        IF ( check ) THEN
           WRITE(numout,*) 'thermosoil_diaglev -- thermosoil_diaglev -- thermosoil_diaglev --' 
           DO jd = 1, nbdl
              WRITE(numout,*) jd, '-', intfact(jd,1:ngrnd)
           ENDDO
           WRITE(numout,*) "SUM -- SUM -- SUM SUM -- SUM -- SUM"
           DO jd = 1, nbdl
              WRITE(numout,*) jd, '-', SUM(intfact(jd,1:ngrnd))
           ENDDO
           WRITE(numout,*) 'thermosoil_diaglev -- thermosoil_diaglev -- thermosoil_diaglev --' 
        ENDIF
        
    ENDIF

 !! 2. does the interpolation

    stempdiag(:,:) = zero
    DO jg = 1, ngrnd
      DO jd = 1, nbdl
        DO ji = 1, kjpindex
          stempdiag(ji,jd) = stempdiag(ji,jd) + ptn(ji,jg)*intfact(jd,jg)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE thermosoil_diaglev

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_humlev
!!
!>\BRIEF        Interpolates the diagnostic soil humidity profile shumdiag(nbdl, diagnostic axis) onto 
!!              the thermal axis, which gives wetdiag(ngrnd, thermal axis).
!!
!! DESCRIPTION  : Same as in thermosoil_diaglev : This is a very easy linear interpolation, with intfactw(jd, jg) the fraction
!! the thermal layer jd comprised within the diagnostic layer jg. 
!!?? I would think wise to change the indeces here, to keep jD for Diagnostic
!!?? and jG for Ground thermal levels...
!! 
!! The depths of the diagnostic levels are diaglev(1:nbdl), computed in slowproc.f90.
!! Recall that when the 11-layer hydrology is used,
!! wetdiag and shumdiag are with reference to the moisture content (mc)
!! at the wilting point mcw : wetdiag=(mc-mcw)/(mcs-mcw).
!! with mcs the saturated soil moisture content.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): wetdiag (soil soil humidity profile on the thermal axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoil_humlev(kjpindex, shumdiag, snow)
  
  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                            :: kjpindex    !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)    :: shumdiag    !! Relative soil humidity on the diagnostic axis. 
                                                                         !! (0-1, unitless). Caveats : when "hydrol" (the 11-layers
                                                                         !! hydrology) is used, this humidity is calculated with 
                                                                         !! respect to the wilting point : 
                                                                         !! shumdiag= (mc-mcw)/(mcs-mcw), with mc : moisture 
                                                                         !! content; mcs : saturated soil moisture content; mcw: 
                                                                         !! soil moisture content at the wilting point. when the 2-layers
                                                                         !! hydrology "hydrolc" is used, shumdiag is just
                                                                         !! a diagnostic humidity index, with no real physical 
                                                                         !! meaning.
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow 
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                       :: ji, jd, jg
    REAL(r_std)                                          :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL(r_std), DIMENSION(ngrnd,nbdl)                   :: intfactw     !! fraction of each diagnostic layer (jd) comprized within
                                                                         !! a given thermal layer (jg)(0-1, unitless) 
    REAL(r_std), DIMENSION(kjpindex)               :: snow_h
    LOGICAL, PARAMETER :: check=.FALSE.

!_ ================================================================================================================================
    
  !! 1. computes intfactw(jd,jg), the fraction of each diagnostic layer (jg) comprized within a given thermal layer (jd)
    IF ( check ) &
         WRITE(numout,*) 'thermosoil_humlev --' 

    ! Snow height
    snow_h(:)=snow(:)/sn_dens
    !
    wetdiag(:,:) = zero
    DO ji=1,kjpindex
       prev_diag = zero
       DO jd = 1, ngrnd
          lev_diag = prev_diag + dz2(jd)
          prev_prog = snow_h(ji)
          DO jg = 1, nbdl
             IF ( jg == nbdl .AND. diaglev(jg)+snow_h(ji) < lev_diag ) THEN
                lev_prog = lev_diag+snow_h(ji)
             ELSE
                lev_prog = diaglev(jg)+snow_h(ji)
             ENDIF
             intfactw(jd,jg) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog), zero)/(lev_diag-prev_diag)
             prev_prog = lev_prog
          ENDDO
          prev_diag = lev_diag
       ENDDO

      ! write(*,*) 'THERMOSOIL HUMLEV 1 ', wetdiag(1,2) !CLTEST
       
       DO jg = 1, nbdl
          DO jd = 1, ngrnd
            wetdiag(ji,jd) = wetdiag(ji,jd) &
     &          + shumdiag(ji,jg)*intfactw(jd,jg)
          ENDDO !DO jd = 1, ngrnd
       ENDDO !DO jg = 1, nbdl
       !
       IF ( check ) THEN
          DO jd = 1, ngrnd
             WRITE(numout,*) ji,jd, '-', intfactw(jd,1:nbdl),'-sum-', SUM(intfactw(jd,1:nbdl))
          ENDDO
       ENDIF
    ENDDO !do ji=1,kjpindex
    
     !write(*,*) 'THERMOSOIL HUMLEV 2 ', wetdiag(1,2) !CLTEST

     IF ( check ) &
        WRITE(numout,*) 'thermosoil_humlev --' 
     ! write(*,*) 'thermosoil_humlev 1205: wetdiag=',wetdiag(1,1)

  END SUBROUTINE thermosoil_humlev


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_energy
!!
!>\BRIEF         Energy check-up.
!!
!! DESCRIPTION  : I didn\'t comment this routine since at do not understand its use, please
!! ask initial designers (Jan ? Nathalie ?).
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ??
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_energy(kjpindex, temp_sol_new, soilcap, first_call)

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! Domain size (unitless)
    LOGICAL, INTENT (in)                           :: first_call   !! First call (true/false)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: temp_sol_new !! Surface temperature at the present time-step, Ts 
                                                                   !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: soilcap      !! Apparent surface heat capacity 
                                                                   !! @tex ($J m^{-2} K^{-1}$) @endtex, 
                                                                   !! see eq. A29 of F. Hourdin\'s PhD thesis.
    
    !! 0.2 Output variables

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    
    INTEGER(i_std)                                 :: ji, jg
!_ ================================================================================================================================

    IF (first_call) THEN

     DO ji = 1, kjpindex
      surfheat_incr(ji) = zero
      coldcont_incr(ji) = zero
      temp_sol_beg(ji)  = temp_sol_new(ji)
      
      DO jg = 1, ngrnd
       ptn_beg(ji,jg)   = ptn(ji,jg)
      ENDDO
      
     ENDDO
    
     RETURN

    ENDIF

     DO ji = 1, kjpindex
      surfheat_incr(ji) = zero
      coldcont_incr(ji) = zero
     ENDDO
    
    !  Sum up the energy content of all layers in the soil.
    DO ji = 1, kjpindex
   
       IF (pcapa_en(ji,1) .LE. sn_capa) THEN
          
          ! Verify the energy conservation in the surface layer
          coldcont_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          surfheat_incr(ji) = zero
       ELSE
          
          ! Verify the energy conservation in the surface layer
          surfheat_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          coldcont_incr(ji) = zero
       ENDIF
    ENDDO

    ptn_beg(:,:)      = ptn(:,:)
    temp_sol_beg(:)   = temp_sol_new(:)

  END SUBROUTINE thermosoil_energy


!Isa++
  SUBROUTINE thermosoil_readjust(kjpindex, ptn)
    ! interface description
    ! input scalar
    INTEGER(i_std), INTENT(in)                             :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(inout)    :: ptn
    INTEGER(i_std)  :: ji, jg
    REAL(r_std) :: ptn_tmp
    !WARNING : pcapa est en J/K/m; pcappa_supp est en J/K
    !local
    REAL(r_std), DIMENSION(kjpindex, ngrnd, ngrnd) :: intfact_soil_th

    do jg=1, ngrnd
        do ji=1, kjpindex
        !Isa : here, all soil latent energy is put into E_sol_lat_couche(ji, 1)
        !because the variable soil layers make it difficult to keep track of all
        !layers in this version
        E_sol_lat_couche(ji, 1)=E_sol_lat_couche(ji, 1)+pcappa_supp(ji,jg)*(ptn(ji,jg)-ptn_beg(ji,jg))

        enddo
   enddo

!Isa test
    do ji=1, kjpindex

        if (E_sol_lat_couche(ji,1).GT.min_sechiba.AND.MINVAL(ptn(ji,:)).GT.tzero+fr_dT/2.) then
                !1. normalement, à cette température, il n'y a plus de neige
                !donc pas la peine d'utiliser x1 et x2
                !2. répartition de l'excess d'énergie sur 2.7m = 6 premiers
                !niveaux..
         !plus d'énergie au dégel qu'au gel; cette énergie aurait dû chauffer au lieu de dégeler 
		!=> on monte les températures en se limitant à 0.5°C à chaque pas de temps..
                do jg=1,6
		ptn_tmp=ptn(ji,jg)

		ptn(ji,jg)=ptn(ji,jg)+min(E_sol_lat_couche(ji,1)/pcapa(ji,jg)/zz_coef(6), 0.5)
		E_sol_lat_couche(ji,1)=E_sol_lat_couche(ji,1)-(ptn(ji,jg)-ptn_tmp)*pcapa(ji,jg)*dz2(jg)
                enddo
		else if (E_sol_lat_couche(ji,1).LT.-min_sechiba.AND.MINVAL(ptn(ji,:)).GT.tzero+fr_dT/2.) then
		!pas assez d'énergie lors du dégel. Cette énergie aurait dû dégeler au lieu de chauffer;
		!=> on rabaisse les températures en conséquence
                do jg=1,6
                ptn_tmp=ptn(ji,jg)

		ptn(ji,jg)=max(tzero+fr_dT/2., ptn_tmp+E_sol_lat_couche(ji,1)/pcapa(ji,jg)/zz_coef(6))
		E_sol_lat_couche(ji,1)=E_sol_lat_couche(ji,1)+(ptn_tmp-ptn(ji,jg))*pcapa(ji,jg)*dz2(jg)
               enddo
	endif 
    enddo

  END SUBROUTINE thermosoil_readjust
! Isa --
  !Isa ajout Bruno
  !-------------------------------------------------------------------

  SUBROUTINE thermosoil_wlupdate( kjpindex, dt, ptn, hsd, hsdlong )
  ! Updates the long-term soil humidity
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),INTENT(in)			 	:: dt
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)   	:: hsd
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(inout) 	:: hsdlong
    ! CR: debugage pour voir où c'est actif
    integer(i_std) ::  ji,jg,nactif

    !
     WHERE ( ptn(:,:) .GT. tzero + fr_dT/2. )
      hsdlong(:,:) = ( hsd(:,:) * dt + hsdlong(:,:) * (tau_freezesoil-dt) ) / tau_freezesoil
!Isa try.. could be used one day...
!hsdlong(:,:) = hsd(:,:)
     ENDWHERE 

    ! verif si c'est actif:
    nactif=0
        do ji=1,kjpindex
          do jg=1,ngrnd
            if (ptn(ji,jg) .GT. tzero + fr_dT/2. ) then
                nactif=nactif+1
            endif
          enddo
        enddo
!    write(*,*) 'thermosoil 1382: ok_wetdiaglong,nactif=',ok_wetdiaglong,nactif
!     if (kjpindex.ge.9) then
!       write(*,*) 'thermosoil 1390: hsd(9,1)=',hsd(9,1)
!       write(*,*) 'hsdlong(9,1)=',hsdlong(9,1)
!     endif

   END SUBROUTINE thermosoil_wlupdate
   
!-------------------------------------------------------------------

!Isa ajout from thermosoil_bruno
  SUBROUTINE thermosoil_getdiff( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en,profil_froz, pcappa_supp)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)				:: kjpindex
!Isa E_corr : in -> inout
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(inout)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)   	:: wetdiaglong
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pkappa
!Isa
     REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcappa_supp
!Isa modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: profil_froz
    !    
    REAL						:: x, E_supp, Epaiss !Isa
    REAL(r_std), DIMENSION(kjpindex)             	:: snow_h
    REAL(r_std), DIMENSION(kjpindex,ngrnd) 		:: zx1, zx2    
    ! Heat capacity of ice/water mixture
    REAL			   			:: cap_iw
    ! Thermal conductivity for saturated soil
    REAL						:: csat
    INTEGER						:: ji,jg


    DO ji = 1,kjpindex
      
      ! 1. Determine the fractions of snow and soil
      
      snow_h(ji) = snow(ji) / sn_dens
      
      !
      !  1.1. The first level
      !
      IF ( snow_h(ji) .GT. zz_coef(1) ) THEN

          ! the 1st level is in the snow => the 1st layer is entirely snow
          zx1(ji,1) = 1.
	  zx2(ji,1) = 0.
	  	  
      ELSE IF ( snow_h(ji) .GT. zero ) THEN      
      
          ! the 1st level is beyond the snow and the snow is present
          zx1(ji,1) = snow_h(ji) / zz_coef(1)
          zx2(ji,1) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)	
	 
      ELSE
      
          ! there is no snow at all, quoi ;-)
          zx1(ji,1) = 0.
	  zx2(ji,1) = 1.       
	  
      ENDIF
      
      !
      !  1.2. The other levels except the two last (too deep to be accounted for by the current hydrology??)
      !

      DO jg = 2, ngrnd !- 2
        IF ( snow_h(ji) .GT. zz_coef(jg) ) THEN

            ! the current level is in the snow => the current layer is entirely snow
            zx1(ji,jg) = 1.
	    zx2(ji,jg) = 0.
	  	  
        ELSE IF ( snow_h(ji) .GT. zz_coef(jg-1) ) THEN
	
   	    ! the current layer is partially snow
            zx1(ji,jg) = (snow_h(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
            zx2(ji,jg) = ( zz_coef(jg) - snow_h(ji)) / (zz_coef(jg) - zz_coef(jg-1))

        ELSE
	
	    ! both levels are out of snow => the current layer is entirely soil	  
            zx1(ji,jg) = 0.
	    zx2(ji,jg) = 1.       
	  	
        ENDIF
      ENDDO
      DO jg = 1, ngrnd !-2 
         !
         ! 2. Calculate heat capacity with allowance for permafrost
         !    For the moment we don't take into account porosity changes related mx_eau_eau changes
         !    which would be reflected in so_capa_wet. 
         !    so_capa_ice implies porosity of 0.15 (mx_eau_eau=150) -> not anymore :
         !    Isa : changed to account for poros = 0.4
         !    Isa : changed wetdiag to have the real fraction of porosity filled with water

         ! 2.1. soil heat capacity depending on temperature and humidity
	 
         IF (ptn(ji,jg) .LT. tzero-fr_dT/2.) THEN

	    ! frozen soil
             pcapa(ji,jg) = so_capa_dry + wetdiaglong(ji,jg)*(so_capa_ice - so_capa_dry)!Isa : old version, proved to be correct
	    ! pcapa(ji,jg) = 2.E6 !DKtest - compar. w/ theor. sol.

	    profil_froz(ji,jg) = 1.
  !Isa E
 	    pcappa_supp(ji,jg)= 0.

     	 ELSEIF (ptn(ji,jg) .GT. tzero+fr_dT/2.) THEN
	    ! unfrozen soil	 
     	    pcapa(ji,jg) =  so_capa_dry + wetdiaglong(ji,jg)*(so_capa_wet - so_capa_dry)!Isa : old version, proved to be correct
    	     ! pcapa(ji,jg) = 4.E6 !DKtest - compar. w/ theor. sol.

	    profil_froz(ji,jg) = 0.
 !Isa E
 	    pcappa_supp(ji,jg)= 0.

     	 ELSE

     	   ! x is the unfrozen fraction of soil water	   	   
     	   x = (ptn(ji,jg)-(tzero-fr_dT/2.)) / fr_dT
	   
           profil_froz(ji,jg) = (1. - x)

	   ! net heat capacity of the ice/water mixture
     	   cap_iw = x * so_capa_wet + (1.-x) * so_capa_ice
     	    ! cap_iw = x * 4.E6 + (1.-x) * 2.E6 !DKtest - compar. w/ theor. sol. 

	   pcapa(ji,jg) = so_capa_dry + wetdiaglong(ji,jg)*(cap_iw-so_capa_dry) + wetdiaglong(ji,jg)*poros*lhf*rho_water/fr_dT

!Isa E
 	   pcappa_supp(ji,jg)= wetdiaglong(ji,jg)*poros*lhf*rho_water/fr_dT*zx2(ji,jg)*dz2(jg)
	   
 	   IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, hcap_trans: ', jg, zz(jg), pcapa(ji,jg)	 
	   
         ENDIF

         !
	 ! 2.2. Take into account the snow and soil fractions in the layer
	 !
	 
         pcapa(ji,jg) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg)

	 !
	 ! 2.3. Calculate the heat capacity for energy conservation check 
	 !        (??, does not influence other results, just written to history file)
	 
	 IF ( zx1(ji,jg).GT.0. ) THEN
            pcapa_en(ji,jg) = sn_capa
	 ELSE
            pcapa_en(ji,jg) = pcapa(ji,jg)
	 ENDIF
 
         !
         ! 3. Calculate the heat conductivity with allowance for permafrost (Farouki, 1981, Cold Reg. Sci. Technol.)
         !
	 
	 ! 3.1. unfrozen fraction
	 
         x = (ptn(ji,jg)-(tzero-fr_dT/2.)) / fr_dT * poros
         x = MIN( poros, MAX( 0., x ) )

         ! 3.2. saturated conductivity
	 
      	 csat = cond_solid**(1.-poros) * cond_ice**(poros-x) * cond_water**x

     	
       	 ! 3.3. unsaturated conductivity

     	 pkappa(ji,jg) =(csat - so_cond_dry)*wetdiaglong(ji,jg) + so_cond_dry

         IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, cond: ', jg, zz(jg), pkappa(ji,jg)	 

         !
	 ! 3.4. Take into account the snow and soil fractions in the layer
	 	 
         pkappa(ji,jg) = un / ( zx1(ji,jg) / sn_cond + zx2(ji,jg) / pkappa(ji,jg) )

         IF (bavard.GE.9) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, cond t/a snow:', jg, zz(jg), pkappa(ji,jg)

	 
      END DO            
    ENDDO   


   
   END SUBROUTINE thermosoil_getdiff


!Isa
     SUBROUTINE thermosoil_getdiff_old_thermix_with_snow( kjpindex, ptn, wetdiag, snow, pkappa, pcapa, pcapa_en,profil_froz)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)   	:: wetdiag
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pkappa
!modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: profil_froz
    !    
    REAL						:: x
    REAL(r_std), DIMENSION(kjpindex)             	:: snow_h
    REAL(r_std), DIMENSION(kjpindex,ngrnd) 		:: zx1, zx2    
    INTEGER						:: ji,jg


    DO ji = 1,kjpindex
      
      ! 1. Determine the fractions of snow and soil
      
      snow_h(ji) = snow(ji) / sn_dens
      
      !
      !  1.1. The first level
      !
      IF ( snow_h(ji) .GT. zz_coef(1) ) THEN

          ! the 1st level is in the snow => the 1st layer is entirely snow
          zx1(ji,1) = 1.
	  zx2(ji,1) = 0.
	  	  
      ELSE IF ( snow_h(ji) .GT. zero ) THEN      
      
          ! the 1st level is beyond the snow and the snow is present
          zx1(ji,1) = snow_h(ji) / zz_coef(1)
          zx2(ji,1) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)	
	 
      ELSE
      
          ! there is no snow at all, quoi ;-)
          zx1(ji,1) = 0.
	  zx2(ji,1) = 1.       
	  
      ENDIF
      
      !
      !  1.2. The other levels except the two last (too deep to be accounted for by the current hydrology??)
      !

      DO jg = 2, ngrnd !- 2
        IF ( snow_h(ji) .GT. zz_coef(jg) ) THEN

            ! the current level is in the snow => the current layer is entirely snow
            zx1(ji,jg) = 1.
	    zx2(ji,jg) = 0.
	  	  
        ELSE IF ( snow_h(ji) .GT. zz_coef(jg-1) ) THEN
	
   	    ! the current layer is partially snow
            zx1(ji,jg) = (snow_h(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
            zx2(ji,jg) = ( zz_coef(jg) - snow_h(ji)) / (zz_coef(jg) - zz_coef(jg-1))

        ELSE
	
	    ! both levels are out of snow => the current layer is entirely soil	  
            zx1(ji,jg) = 0.
	    zx2(ji,jg) = 1.       
	  	
        ENDIF
      ENDDO
            
      DO jg = 1, ngrnd !-2 
         !
         ! 2. Calculate frozen profile for hydrolc.f90
	 !
	 IF (bavard.GE.8) write(*,'(A24,I3,F5.2,G14.7)')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg)
	 
         IF (ptn(ji,jg) .LT. tzero-fr_dT/2.) THEN
	    profil_froz(ji,jg) = 1.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg)	 
	    
     	 ELSEIF (ptn(ji,jg) .GT. tzero+fr_dT/2.) THEN
	    profil_froz(ji,jg) = 0.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, : Soil_temp', jg, zz(jg), ptn(ji,jg)	 
	   
     	 ELSE
	 
     	   ! x is the unfrozen fraction of soil water	   	   
     	   x = (ptn(ji,jg)-(tzero-fr_dT/2.)) / fr_dT	   
           profil_froz(ji,jg) = (1. - x)
	   
         ENDIF

         ! 3. heat capacity calculation
	 !
         ! 3.0 old heat capacity calculation
        pcapa(ji,jg) = so_capa_dry + wetdiag(ji,jg)*(so_capa_wet - so_capa_dry)

	 ! 3.1. Still some improvement from the old_version : Take into account the snow and soil fractions in the layer

         pcapa(ji,jg) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg)

	 ! 3.2. Calculate the heat capacity for energy conservation check 
	 !        (??, does not influence other results, just written to history file)
	 
	 IF ( zx1(ji,jg).GT.0. ) THEN
            pcapa_en(ji,jg) = sn_capa
	 ELSE
            pcapa_en(ji,jg) = pcapa(ji,jg)
	 ENDIF
         !
         !4. heat conductivity calculation
	 !
         !4.0 old heat conductivity calculation
         pkappa(ji,jg) = so_cond_dry + wetdiag(ji,jg)*(so_cond_wet - so_cond_dry)

         !4.0 Still some improvement from the old_version : Take into account the snow and soil fractions in the layer

         pkappa(ji,jg) = un / ( zx1(ji,jg) / sn_cond + zx2(ji,jg) / pkappa(ji,jg) )


	 
      END DO            
      !
    ENDDO   
    
   
   END SUBROUTINE thermosoil_getdiff_old_thermix_with_snow

!Isa
     SUBROUTINE thermosoil_getdiff_old_thermix( kjpindex, ptn, wetdiag, snow, pkappa, pcapa, pcapa_en,profil_froz)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)   	:: wetdiag
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pkappa
!modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: profil_froz
    !    
    REAL						:: x 
    INTEGER						:: ji,jg

    
    DO ji = 1,kjpindex            
      DO jg = 1, ngrnd 
         !
         ! 2. Calculate frozen profile for hydrolc.f90
	 !
	 IF (bavard.GE.8) write(*,'(A24,I3,F5.2,G14.7)')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg)
	 
         IF (ptn(ji,jg) .LT. tzero-fr_dT/2.) THEN
	    profil_froz(ji,jg) = 1.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg)	 
	    
     	 ELSEIF (ptn(ji,jg) .GT. tzero+fr_dT/2.) THEN
	    profil_froz(ji,jg) = 0.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, : Soil_temp', jg, zz(jg), ptn(ji,jg)	 
	   
     	 ELSE
	 
     	   ! x is the unfrozen fraction of soil water	   	   
     	   x = (ptn(ji,jg)-(tzero-fr_dT/2.)) / fr_dT	   
           profil_froz(ji,jg) = (1. - x)
	   
         ENDIF

         ! 3. heat capacity calculation
	 !

         ! 3.1 old heat capacity calculation
        pcapa(ji,jg) = so_capa_dry + wetdiag(ji,jg)*(so_capa_wet - so_capa_dry)

	 ! 3.2. Calculate the heat capacity for energy conservation check 
	 !        (??, does not influence other results, just written to history file)
	pcapa_en(ji,1) = so_capa_dry + wetdiag(ji,1)*(so_capa_wet - so_capa_dry)

         !4. heat conductivity calculation
	 !
         !4.0 old heat conductivity calculation
         pkappa(ji,jg) = so_cond_dry + wetdiag(ji,jg)*(so_cond_wet - so_cond_dry)

	 
      END DO            
      !
    ENDDO   
   
   END SUBROUTINE thermosoil_getdiff_old_thermix

   SUBROUTINE thermosoil_getdiff_old_thermix_trunc( kjpindex, snow, pkappa, pcapa, pcapa_en )

    INTEGER(i_std), INTENT(in) :: kjpindex
    
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pkappa
    INTEGER						:: ji,jg
    REAL(r_std), DIMENSION(kjpindex)                       :: snow_h       !! snow_h is the snow height @tex ($m$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)                       :: zx1, zx2     !! zx1 and zx2 
                             !! are the layer fraction consisting in snowand soil respectively.

     
    ! Computation of the soil thermal properties; snow properties are also accounted for

    DO ji = 1,kjpindex
      snow_h(ji) = snow(ji) / sn_dens

      IF ( snow_h(ji) .GT. zz_coef(1) ) THEN
          pcapa(ji,1) = sn_capa
          pcapa_en(ji,1) = sn_capa
          pkappa(ji,1) = sn_cond
      ELSE IF ( snow_h(ji) .GT. zero ) THEN
          pcapa_en(ji,1) = sn_capa
          zx1(ji) = snow_h(ji) / zz_coef(1)
          zx2(ji) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)
          pcapa(ji,1) = zx1(ji) * sn_capa + zx2(ji) * so_capa_wet
          pkappa(ji,1) = un / ( zx1(ji) / sn_cond + zx2(ji) / so_cond_wet )
      ELSE
          pcapa(ji,1) = so_capa_dry + wetdiag(ji,1)*(so_capa_wet - so_capa_dry)
          pkappa(ji,1) = so_cond_dry + wetdiag(ji,1)*(so_cond_wet - so_cond_dry)
          pcapa_en(ji,1) = so_capa_dry + wetdiag(ji,1)*(so_capa_wet - so_capa_dry)
      ENDIF
      !
      DO jg = 2, ngrnd - 2
        IF ( snow_h(ji) .GT. zz_coef(jg) ) THEN
            pcapa(ji,jg) = sn_capa
            pkappa(ji,jg) = sn_cond
            pcapa_en(ji,jg) = sn_capa
        ELSE IF ( snow_h(ji) .GT. zz_coef(jg-1) ) THEN
            zx1(ji) = (snow_h(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
            zx2(ji) = ( zz_coef(jg) - snow_h(ji)) / (zz_coef(jg) - zz_coef(jg-1))
            pcapa(ji,jg) = zx1(ji) * sn_capa + zx2(ji) * so_capa_wet
            pkappa(ji,jg) = un / ( zx1(ji) / sn_cond + zx2(ji) / so_cond_wet )
            pcapa_en(ji,jg) = sn_capa
        ELSE
            pcapa(ji,jg) = so_capa_dry + wetdiag(ji,jg)*(so_capa_wet - so_capa_dry)
            pkappa(ji,jg) = so_cond_dry + wetdiag(ji,jg)*(so_cond_wet - so_cond_dry)
            pcapa_en(ji,jg) = so_capa_dry + wetdiag(ji,jg)*(so_capa_wet - so_capa_dry)
        ENDIF
      ENDDO
    
    ENDDO ! DO ji = 1,kjpindex
 
    END SUBROUTINE thermosoil_getdiff_old_thermix_trunc

    SUBROUTINE thermosoil_getdiff_old_thermix_trunc2( kjpindex,  pkappa, pcapa, pcapa_en )

    INTEGER(i_std), INTENT(in) :: kjpindex
    
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(out)    :: pkappa
    INTEGER						:: ji,jg

     
    DO jg = 1,ngrnd
      DO ji = 1,kjpindex
        pkappa(ji,jg) = so_cond_dry + wetdiag(ji,jg)*(so_cond_wet - so_cond_dry)
        pcapa(ji,jg) = so_capa_dry + wetdiag(ji,jg)*(so_capa_wet - so_capa_dry)
        pcapa_en(ji,jg) = so_capa_dry + wetdiag(ji,jg)*(so_capa_wet - so_capa_dry)
      ENDDO
    ENDDO

    END SUBROUTINE thermosoil_getdiff_old_thermix_trunc2

  SUBROUTINE read_permafrostmap(kjpindex,lalo,overburden,excess_ice,permafrost)
    
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout) :: overburden
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout) :: excess_ice
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout) :: permafrost
    
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: xx,yy, permafrost_file, continuous_file, discontinuous_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: sporadic_file, isolated_file, overburden_file, excess_ice_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:) :: x,y
    REAL(r_std) :: lev(1), date, dt
    INTEGER(i_std) :: itau(1)
    REAL(r_std),DIMENSION(kjpindex) :: tref
    
    ! plus bas, on prend la temperature lue dans un fichier climato si celui-ci existe
    filename = 'NONE'
    CALL getin('PERMAFROST_MAP_FILE',filename)
    IF ( filename .EQ. "NONE" ) THEN
    ELSE
       CALL flininfo(filename,iml, jml, lml, tml, fid)
       ALLOCATE (yy(iml,jml), stat=ier)
       ALLOCATE (xx(iml,jml), stat=ier)
       ALLOCATE (x(iml),y(jml), stat=ier)
       ALLOCATE (continuous_file(iml,jml), stat=ier)
       ALLOCATE (discontinuous_file(iml,jml), stat=ier)
       ALLOCATE (sporadic_file(iml,jml), stat=ier)
       ALLOCATE (isolated_file(iml,jml), stat=ier)
       ALLOCATE (overburden_file(iml,jml), stat=ier)
       ALLOCATE (excess_ice_file(iml,jml), stat=ier)
       ALLOCATE (permafrost_file(iml,jml), stat=ier)
       CALL flinopen (filename, .FALSE., iml, jml, lml, &
            xx, yy, lev, tml, itau, date, dt, fid)
       CALL flinget (fid, 'continuous_permafrost', iml, jml, lml, tml, &
            1, 1, continuous_file)
       CALL flinget (fid, 'discontinuous_permafrost', iml, jml, lml, tml, &
            1, 1, discontinuous_file)
       CALL flinget (fid, 'sporadic_permafrost', iml, jml, lml, tml, &
            1, 1, sporadic_file)
       CALL flinget (fid, 'isolated_permafrost', iml, jml, lml, tml, &
            1, 1, isolated_file)
       CALL flinget (fid, 'thick_overburden', iml, jml, lml, tml, &
            1, 1, overburden_file)
       CALL flinget (fid, 'high_ground_ice_content', iml, jml, lml, tml, &
            1, 1, excess_ice_file)
       CALL flinclo (fid)
       ! On suppose que le fichier est regulier.
       ! Si ce n'est pas le cas, tant pis. Les temperatures seront mal
       ! initialisees et puis voila. De toute maniere, il faut avoir
       ! l'esprit mal tourne pour avoir l'idee de faire un fichier de
       ! climatologie avec une grille non reguliere.
       permafrost_file(:,:) = continuous_file + discontinuous_file + sporadic_file + isolated_file
       x(:) = xx(:,1)
       y(:) = yy(1,:)
       ! prendre la valeur la plus proche
       DO ip = 1, kjpindex
          dlonmin = HUGE(1.)
          DO ix = 1,iml
             dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
             IF ( dlon .LT. dlonmin ) THEN
                imin = ix
                dlonmin = dlon
             ENDIF
          ENDDO
          dlatmin = HUGE(1.)
          DO iy = 1,jml
             dlat = ABS(lalo(ip,1)-y(iy))
             IF ( dlat .LT. dlatmin ) THEN
                jmin = iy
                dlatmin = dlat
             ENDIF
          ENDDO
          permafrost(ip) = permafrost_file(imin,jmin)
          overburden(ip) = overburden_file(imin,jmin)
          excess_ice(ip) = excess_ice_file(imin,jmin)
       ENDDO
       DEALLOCATE (yy)
       DEALLOCATE (xx)
       DEALLOCATE (x)
       DEALLOCATE (continuous_file)
       DEALLOCATE (discontinuous_file)
       DEALLOCATE (sporadic_file)
       DEALLOCATE (isolated_file)
       DEALLOCATE (overburden_file)
       DEALLOCATE (excess_ice_file)
       DEALLOCATE (permafrost_file)
    ENDIF
    WRITE(*,*) 'cdk: #points permafrost', sum(permafrost)
    WRITE(*,*) 'cdk: #points overburden', sum(overburden)
    WRITE(*,*) 'cdk: #points excess_ice', sum(excess_ice)
    
  END SUBROUTINE read_permafrostmap

 SUBROUTINE read_reftempfile(kjpindex,lalo,reftemp)
    
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo
    REAL(r_std), DIMENSION(kjpindex, ngrnd), INTENT(inout) :: reftemp

    
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: xx,yy
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: reftemp_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:) :: x,y
    REAL(r_std) :: lev(1), date, dt
    INTEGER(i_std) :: itau(1)
    REAL(r_std),DIMENSION(kjpindex) :: tref
    
    ! plus bas, on prend la temperature lue dans un fichier climato si celui-ci existe
    filename = 'reftemp.nc'
    CALL getin('REFTEMP_FILE',filename)

       CALL flininfo(filename,iml, jml, lml, tml, fid)
       ALLOCATE (yy(iml,jml), stat=ier)
       ALLOCATE (xx(iml,jml), stat=ier)
       ALLOCATE (x(iml),y(jml), stat=ier)
       ALLOCATE (reftemp_file(iml,jml), stat=ier)

       CALL flinopen (filename, .FALSE., iml, jml, lml, &
            xx, yy, lev, tml, itau, date, dt, fid)
       CALL flinget (fid, 'temperature', iml, jml, lml, tml, &
            1, 1, reftemp_file)
       CALL flinclo (fid)
       ! On suppose que le fichier est regulier.
       ! Si ce n'est pas le cas, tant pis. Les temperatures seront mal
       ! initialisees et puis voila. De toute maniere, il faut avoir
       ! l'esprit mal tourne pour avoir l'idee de faire un fichier de
       ! climatologie avec une grille non reguliere.
       x(:) = xx(:,1)
       y(:) = yy(1,:)
       ! prendre la valeur la plus proche
       DO ip = 1, kjpindex
          dlonmin = HUGE(1.)
          DO ix = 1,iml
             dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
             IF ( dlon .LT. dlonmin ) THEN
                imin = ix
                dlonmin = dlon
             ENDIF
          ENDDO
          dlatmin = HUGE(1.)
          DO iy = 1,jml
             dlat = ABS(lalo(ip,1)-y(iy))
             IF ( dlat .LT. dlatmin ) THEN
                jmin = iy
                dlatmin = dlat
             ENDIF
          ENDDO
          reftemp(ip, :) = reftemp_file(imin,jmin)+273.15
       ENDDO
       DEALLOCATE (yy)
       DEALLOCATE (xx)
       DEALLOCATE (x)
       DEALLOCATE (reftemp_file)


    
  END SUBROUTINE read_reftempfile



END MODULE thermosoil
