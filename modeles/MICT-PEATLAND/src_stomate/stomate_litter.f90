! =================================================================================================================================
! MODULE       : stomate_litter
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Update litter and lignine content after litter fall and 
!! calculating litter decomposition.      
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_stomate/stomate_litter.f90 $
!! $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!! $Revision: 947 $
!! \n
!_ ================================================================================================================================

MODULE stomate_litter

  ! modules used:

  USE ioipsl
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC littercalc,littercalc_clear, deadleaf

  LOGICAL, SAVE                        :: firstcall = .TRUE.       !! first call

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : littercalc_calc
!!
!!\BRIEF        Set the flag ::firstcall to .TRUE. and as such activate section
!! 1.1 of the subroutine littercalc (see below).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE littercalc_clear
    firstcall =.TRUE.
  END SUBROUTINE littercalc_clear


!! ================================================================================================================================
!! SUBROUTINE   : littercalc
!!
!!\BRIEF        Calculation of the litter decomposition and therefore of the 
!! heterotrophic respiration from litter following Parton et al. (1987).
!!
!! DESCRIPTION  : The littercal routine splits the litter in 4 pools: 
!! aboveground metaboblic, aboveground structural, belowground metabolic and 
!! belowground structural. the fraction (F) of plant material going to metabolic 
!! and structural is defined following Parton et al. (1987)
!! \latexonly
!! \input{littercalc1.tex}
!! \endlatexonly
!! \n
!! where L is the lignin content of the plant carbon pools considered and CN 
!! its CN ratio. L and CN are fixed parameters for each plant carbon pools,
!! therefore it is the ratio between each plant carbon pool within a PFT, which
!! controlled the part of the total litter, that will be considered as
!! recalcitrant (i.e. structural litter) or labile (i.e. metabolic litter).\n  
!! 
!! The routine calculates the fraction of aboveground litter which is metabolic
!! or structural (the litterpart variable) which is then used in lpj_fire.f90.\n 
!! 
!! In the section 2, the routine calculate the new plant material entering the
!! litter pools by phenological death of plants organs (corresponding to the
!! variable turnover) and by fire, herbivory and others non phenological causes
!! (variable bm_to_litter). This calculation is first done for each PFT and then
!! the values calculated for each PFT are added up. Following the same approach
!! the lignin content of the total structural litter is calculated and will be
!! then used as a factor control of the decomposition of the structural litter
!! (lignin_struc) in the section 5.1.2. A test is performed to avoid that we add
!! more lignin than structural litter. Finally, the variable litterpart is
!! updated.\n
!! 
!! In the section 3 and 4 the temperature and the moisture controlling the
!! decomposition are calculated for above and belowground. For aboveground
!! litter, air temperature and litter moisture are calculated in sechiba and used 
!! directly. For belowground, soil temperature and moisture are also calculated 
!! in sechiba but are modulated as a function of the soil depth. The modulation 
!! is a multiplying factor exponentially distributed between 0 (in depth) and 1
!! in surface.\n  
!! 
!! Then, in the section 5, the routine calculates the structural litter decomposition 
!! (C) following first order kinetics following Parton et al. (1987).
!! \latexonly
!! \input{littercalc2.tex}
!! \endlatexonly
!! \n
!! with k the decomposition rate of the structural litter. 
!! k corresponds to
!! \latexonly
!! \input{littercalc3.tex}
!! \endlatexonly
!! \n
!! with littertau the turnover rate, T a function of the temperature and M a function of
!! the moisture described below.\n
!!  
!! Then, the fraction of dead leaves (DL) composed by aboveground structural litter is
!! calculated as following
!! \latexonly
!! \input{littercalc4.tex}
!! \endlatexonly
!! \n
!! with k the decomposition rate of the structural litter previously
!! described.\n
!!
!! In the section 5.1, the fraction of decomposed structural litter
!! incorporated to the soil (Input) and its associated heterotrophic respiration are
!! calculated. For structural litter, the C decomposed could go in the active
!! soil carbon pool or in the slow carbon, as described in 
!! stomate_soilcarbon.f90.\n
!! \latexonly
!! \input{littercalc5.tex}
!! \endlatexonly
!! \n
!! with f a parameter describing the fraction of structural litter incorporated
!! into the considered soil carbon pool, C the amount of litter decomposed and L 
!! the amount of lignin in the litter. The litter decomposed which is not
!! incorporated into the soil is respired.\n
!!
!! In the section 5.2, the fraction of decomposed metabolic litter
!! incorporated to the soil and its associated heterotrophic respiration are
!! calculated with the same approaches presented for 5.1 but no control factor
!! depending on the lignin content are used.\n
!! 
!! In the section 6 the dead leaf cover is calculated through a call to the 
!! deadleaf subroutine presented below.\n
!!
!! MAIN OUTPUT VARIABLES: ::deadleaf_cover, ::resp_hetero_litter, ::soilcarbon_input, 
!! ::control_temp, ::control_moist
!!
!! REFERENCES:
!! - Parton, WJ, Schimel, DS, Cole, CV, and Ojima, DS. 1987. Analysis
!! of factors controlling soil organic matter levels in Great Plains
!! grasslands. Soil Science Society of America journal (USA)
!! (51):1173-1179.
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics(scale=0.5){littercalcflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE littercalc (npts, dt, &
       turnover, bm_to_litter, &
       veget_max, tsurf, tsoil, soilhum, litterhum, &
       litterpart, litter, dead_leaves, lignin_struc, &
       deadleaf_cover, resp_hetero_litter, &
       soilcarbon_input, control_temp, control_moist,wpro_conso_carb)
  
    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

!Chloe ++ :
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: wpro_conso_carb    !! quantity of methane consumed by methano (walter'smodel)
!Chloe --        

    INTEGER(i_std), INTENT(in)                                  :: npts               !! Domain size - number of grid pixels
    REAL(r_std), INTENT(in)                                     :: dt                 !! Time step of the simulations for stomate
                                                                                      !! @tex $(dtradia one\_day^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)         :: turnover           !! Turnover rates of plant biomass 
                                                                                      !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)         :: bm_to_litter       !! Conversion of biomass to litter 
                                                                                      !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)                  :: veget_max          !! PFT "Maximal" coverage fraction of a PFT 
                                                                                      !! defined in the input vegetation map 
                                                                                      !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: tsurf              !! Temperature (K) at the surface
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)               :: tsoil              !! Soil temperature (K)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)               :: soilhum            !! Daily soil humidity of each soil layer 
                                                                                      !! (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: litterhum          !! Daily litter humidity (unitless)

    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts), INTENT(out)                   :: deadleaf_cover     !! Fraction of soil covered by dead leaves 
                                                                                      !! over all PFTs (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)               :: resp_hetero_litter !! Litter heterotrophic respiration. The unit
                                                                                      !! is given by m^2 of ground.  
                                                                                      !! @tex $(gC dtradia one\_day^{-1}) m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(out)         :: soilcarbon_input   !! Quantity of carbon going into carbon pools
                                                                                      !! from litter decomposition. The unit is  
                                                                                      !! given by m^2 of ground 
                                                                                      !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(out)             :: control_temp       !! Temperature control of heterotrophic 
                                                                                      !! respiration, above and below (0-1, 
                                                                                      !! unitless)
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(out)             :: control_moist      !! Moisture control of heterotrophic 
                                                                                      !! respiration (0.25-1, unitless)
   
    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)       :: litterpart         !! Fraction of litter above the ground 
                                                                                      !! belonging to the different PFTs (0-1, 
                                                                                      !! unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs), INTENT(inout) :: litter             !! Metabolic and structural litter,above and
                                                                                      !! below ground. The unit is given by m^2 of 
                                                                                      !! ground @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)       :: dead_leaves        !! Dead leaves per ground unit area, per PFT,
                                                                                      !! metabolic and structural in 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)       :: lignin_struc       !! Ratio Lignin content in structural litter,
                                                                                      !! above and below ground, (0-1, unitless)
    
    !! 0.4 Local variables
 
    REAL(r_std), SAVE, DIMENSION(nparts,nlitt)                  :: litterfrac         !! The fraction of leaves, wood, etc. that 
                                                                                      !! goes into metabolic and structural 
                                                                                      !! litterpools (0-1, unitless)
    REAL(r_std), SAVE, DIMENSION(0:nbdl)                        :: z_soil             !! Soil levels (m)
    REAL(r_std), DIMENSION(npts)                                :: rpc                !! Integration constant for vertical root 
                                                                                      !! profiles (unitless)
    REAL(r_std), SAVE, DIMENSION(nlitt)                         :: litter_tau         !! Turnover time in litter pools (days)
    REAL(r_std), SAVE, DIMENSION(nlitt,ncarb,nlevs)             :: frac_soil          !! Fraction of litter that goes into soil 
                                                                                      !! (litter -> carbon, above and below). The
                                                                                      !! remaining part goes to the atmosphere
    REAL(r_std), DIMENSION(npts)                                :: tsoil_decomp       !! Temperature used for decompostition in 
                                                                                      !! soil (K)
    REAL(r_std), DIMENSION(npts)                                :: soilhum_decomp     !! Humidity used for decompostition in soil
                                                                                      !! (unitless)
    REAL(r_std), DIMENSION(npts)                                :: fd                 !! Fraction of structural or metabolic litter
                                                                                      !! decomposed (unitless)
    REAL(r_std), DIMENSION(npts)                                :: qd                 !! Quantity of structural or metabolic litter
                                                                                      !! decomposed @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlevs)                      :: old_struc          !! Old structural litter, above and below 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlitt,nlevs)                :: litter_inc_PFT     !! Increase of litter, per PFT, metabolic and
                                                                                      !! structural, above and below ground. The 
                                                                                      !! unit is given by m^2 of ground.  
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs)                :: litter_inc         !! Increase of metabolic and structural 
                                                                                      !! litter, above and below ground. The unit 
                                                                                      !! is given by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlevs)                      :: lignin_struc_inc   !! Lignin increase in structural litter, 
                                                                                      !! above and below ground. The unit is given 
                                                                                      !! by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlitt)                      :: litter_pft         !! Metabolic and structural litter above the 
                                                                                      !! ground per PFT
    REAL(r_std), DIMENSION(npts)                                :: zdiff_min          !! Intermediate field for looking for minimum
                                                                                      !! of what? this is not used in the code. 
                                                                                      !! [??CHECK] could we delete it?
    CHARACTER(LEN=10), DIMENSION(nlitt)                         :: litter_str         !! Messages to write output information about
                                                                                      !! the litter
    CHARACTER(LEN=22), DIMENSION(nparts)                        :: part_str           !! Messages to write output information about
                                                                                      !! the plant
    CHARACTER(LEN=7), DIMENSION(ncarb)                          :: carbon_str         !! Messages to write output information about
                                                                                      !! the soil carbon
    CHARACTER(LEN=5), DIMENSION(nlevs)                          :: level_str          !! Messages to write output information about
                                                                                      !! the level (aboveground or belowground litter)
    INTEGER(i_std)                                              :: i,j,k,l,m          !! Indices (unitless)
!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering littercalc'

  !! 1. Initialisations of the different fields during the first call of the routine 

    IF ( firstcall ) THEN

       !! 1.1.3 litter fractions:
       !!   what fraction of leaves, wood, etc. goes into metabolic and structural litterpools
       DO k = 1, nparts

          litterfrac(k,imetabolic) = metabolic_ref_frac - metabolic_LN_ratio * LC(k) * CN(k)
          litterfrac(k,istructural) = un - litterfrac(k,imetabolic)

       ENDDO

       !! 1.1.4 residence times in litter pools (days)
       litter_tau(imetabolic) = tau_metabolic * one_year      ! .5 years
       litter_tau(istructural) = tau_struct * one_year     ! 3 years

       !! 1.1.5 decomposition flux fraction that goes into soil
       !       (litter -> carbon, above and below)
       !       1-frac_soil goes into atmosphere
       frac_soil(:,:,:) = zero

       ! structural litter: lignin fraction goes into slow pool + respiration,
       !                    rest into active pool + respiration
       frac_soil(istructural,iactive,iabove) = frac_soil_struct_aa
       frac_soil(istructural,iactive,ibelow) = frac_soil_struct_ab
       frac_soil(istructural,islow,iabove) = frac_soil_struct_sa
       frac_soil(istructural,islow,ibelow) = frac_soil_struct_sb

       ! metabolic litter: all goes into active pool + respiration.
       !   Nothing into slow or passive pool.
       frac_soil(imetabolic,iactive,iabove) = frac_soil_metab_aa
       frac_soil(imetabolic,iactive,ibelow) = frac_soil_metab_ab

       
       !! 1.2 soil levels
       z_soil(0) = zero
       z_soil(1:nbdl) = diaglev(1:nbdl)

       
       !! 1.3 messages
       litter_str(imetabolic) = 'metabolic'
       litter_str(istructural) = 'structural'

       carbon_str(iactive) = 'active'
       carbon_str(islow) = 'slow'
       carbon_str(ipassive) = 'passive'

       level_str(iabove) = 'above'
       level_str(ibelow) = 'below'

       part_str(ileaf) = 'leaves'
       part_str(isapabove) = 'sap above ground'
       part_str(isapbelow) = 'sap below ground'
       part_str(iheartabove) = 'heartwood above ground'
       part_str(iheartbelow) = 'heartwood below ground'
       part_str(iroot) = 'roots'
       part_str(ifruit) = 'fruits'
       part_str(icarbres) = 'carbohydrate reserve'

       WRITE(numout,*) 'litter:'

       WRITE(numout,*) '   > C/N ratios: '
       DO k = 1, nparts
          WRITE(numout,*) '       ', part_str(k), ': ',CN(k)
       ENDDO

       WRITE(numout,*) '   > Lignine/C ratios: '
       DO k = 1, nparts
          WRITE(numout,*) '       ', part_str(k), ': ',LC(k)
       ENDDO

       WRITE(numout,*) '   > fraction of compartment that goes into litter: '
       DO k = 1, nparts
          DO m = 1, nlitt
             WRITE(numout,*) '       ', part_str(k), '-> ',litter_str(m), ':',litterfrac(k,m)
          ENDDO
       ENDDO

       WRITE(numout,*) '   > scaling depth for decomposition (m): ',z_decomp

       WRITE(numout,*) '   > minimal carbon residence time in litter pools (d):'
       DO m = 1, nlitt
          WRITE(numout,*) '       ',litter_str(m),':',litter_tau(m)
       ENDDO

       WRITE(numout,*) '   > litter decomposition flux fraction that really goes '
       WRITE(numout,*) '     into carbon pools (rest into the atmosphere):'
       DO m = 1, nlitt
          DO l = 1, nlevs
             DO k = 1, ncarb
                WRITE(numout,*) '       ',litter_str(m),' ',level_str(l),' -> ',&
                     carbon_str(k),':', frac_soil(m,k,l)
             ENDDO
          ENDDO
       ENDDO

       firstcall = .FALSE.

    ENDIF

    
    !! 1.3 litter above the ground per PFT.
    DO j = 2, nvm ! Loop over # PFTs

       DO k = 1, nlitt !Loop over litter pool
          litter_pft(:,j,k) = litterpart(:,j,k) * litter(:,k,j,iabove)
       ENDDO

    ENDDO

    
    !! 1.4 set output to zero
    deadleaf_cover(:) = zero
    resp_hetero_litter(:,:) = zero
    soilcarbon_input(:,:,:) = zero

    
  !! 2. Add biomass to different litterpools (per m^2 of ground)
    
    !! 2.1 first, save old structural litter (needed for lignin fractions).
    !     above/below
    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 2,nvm !Loop over PFTs

          old_struc(:,m,l) = litter(:,istructural,m,l)

       ENDDO
    ENDDO

    
    !! 2.2 update litter, dead leaves, and lignin content in structural litter
    litter_inc(:,:,:,:) = zero
    lignin_struc_inc(:,:,:) = zero

    DO j = 2,nvm !Loop over PFTs

       !! 2.2.1 litter
       DO k = 1, nlitt    !Loop over litter pools (metabolic and structural)

          !! 2.2.2 calculate litter increase (per m^2 of ground).
          !       Only a given fracion of fruit turnover is directly coverted into litter.
          !       Litter increase for each PFT, structural and metabolic, above/below
               litter_inc_PFT(:,j,k,iabove) = &
               litterfrac(ileaf,k) * bm_to_litter(:,j,ileaf) + &
               litterfrac(isapabove,k) * bm_to_litter(:,j,isapabove) + &
               litterfrac(iheartabove,k) * bm_to_litter(:,j,iheartabove) + &
               litterfrac(ifruit,k) * bm_to_litter(:,j,ifruit) + &
               litterfrac(icarbres,k) * bm_to_litter(:,j,icarbres) + &
               litterfrac(ileaf,k) * turnover(:,j,ileaf) + &
               litterfrac(isapabove,k) * turnover(:,j,isapabove) + &
               litterfrac(iheartabove,k) * turnover(:,j,iheartabove) + &
               litterfrac(ifruit,k) * turnover(:,j,ifruit) + &
               litterfrac(icarbres,k) * turnover(:,j,icarbres)

               litter_inc_PFT(:,j,k,ibelow) = &
               litterfrac(isapbelow,k) * bm_to_litter(:,j,isapbelow) + &
               litterfrac(iheartbelow,k) * bm_to_litter(:,j,iheartbelow) + &
               litterfrac(iroot,k) * bm_to_litter(:,j,iroot) + &
               litterfrac(isapbelow,k) * turnover(:,j,isapbelow) + &
               litterfrac(iheartbelow,k) * turnover(:,j,iheartbelow) + &
               litterfrac(iroot,k) * turnover(:,j,iroot)

          ! litter increase, met/struct, above/below
              litter_inc(:,k,j,iabove) = litter_inc(:,k,j,iabove) + litter_inc_PFT(:,j,k,iabove)
              litter_inc(:,k,j,ibelow) = litter_inc(:,k,j,ibelow) + litter_inc_PFT(:,j,k,ibelow)

          !! 2.2.3 dead leaves, for soil cover.
               dead_leaves(:,j,k) = &
               dead_leaves(:,j,k) + &
               litterfrac(ileaf,k) * ( bm_to_litter(:,j,ileaf) + turnover(:,j,ileaf) )

          !! 2.2.4 lignin increase in structural litter
          IF ( k .EQ. istructural ) THEN

             lignin_struc_inc(:,j,iabove) = &
                  lignin_struc_inc(:,j,iabove) + &
                  LC(ileaf) * bm_to_litter(:,j,ileaf) + &
                  LC(isapabove) * bm_to_litter(:,j,isapabove) + &
                  LC(iheartabove) * bm_to_litter(:,j,iheartabove) + &
                  LC(ifruit) * bm_to_litter(:,j,ifruit) + &
                  LC(icarbres) * bm_to_litter(:,j,icarbres) + &
                  LC(ileaf) * turnover(:,j,ileaf) + &
                  LC(isapabove) * turnover(:,j,isapabove) + &
                  LC(iheartabove) * turnover(:,j,iheartabove) + &
                  LC(ifruit) * turnover(:,j,ifruit) + &
                  LC(icarbres) * turnover(:,j,icarbres)

                  lignin_struc_inc(:,j,ibelow) = &
                  lignin_struc_inc(:,j,ibelow) + &
                  LC(isapbelow) * bm_to_litter(:,j,isapbelow) + &
                  LC(iheartbelow) * bm_to_litter(:,j,iheartbelow) + &
                  LC(iroot) * bm_to_litter(:,j,iroot) + &
                  LC(isapbelow)*turnover(:,j,isapbelow) + &
                  LC(iheartbelow)*turnover(:,j,iheartbelow) + &
                  LC(iroot)*turnover(:,j,iroot)

          ENDIF

       ENDDO
    ENDDO

    !! 2.2.5 add new litter (struct/met, above/below)
    litter(:,:,:,:) = litter(:,:,:,:) + litter_inc(:,:,:,:)

    !! 2.2.6 for security: can't add more lignin than structural litter (above/below)
    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 2,nvm !Lopp over PFTs

          lignin_struc_inc(:,m,l) = &
               MIN( lignin_struc_inc(:,m,l), litter_inc(:,istructural,m,l) )

       ENDDO
    ENDDO

    !! 2.2.7 new lignin content: add old lignin and lignin increase, divide by 
    !!       total structural litter (above/below)
    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 2,nvm !Loop over PFTs
          WHERE( litter(:,istructural,m,l) .GT. min_stomate )

       !MM : Soenke modif
       ! Best vectorization ?
!!$       lignin_struc(:,:,:) = &
!!$            ( lignin_struc(:,:,:)*old_struc(:,:,:) + lignin_struc_inc(:,:,:) ) / &
!!$            litter(:,istructural,:,:,icarbon)

             lignin_struc(:,m,l) = lignin_struc(:,m,l) * old_struc(:,m,l)
             lignin_struc(:,m,l) = lignin_struc(:,m,l) + lignin_struc_inc(:,m,l)
             lignin_struc(:,m,l) = lignin_struc(:,m,l) / litter(:,istructural,m,l)
          ELSEWHERE
             lignin_struc(:,m,l) = zero
          ENDWHERE
       ENDDO
    ENDDO

    
    !! 2.3 new litter fraction per PFT (for structural and metabolic litter, above
    !!       the ground).
    DO j = 2,nvm !Loop over PFTs

       WHERE ( litter(:,:,j,iabove) .GT. min_stomate )

          litterpart(:,j,:) = &
               ( litter_pft(:,j,:) + litter_inc_PFT(:,j,:,iabove) ) / litter(:,:,j,iabove)

       ELSEWHERE

          litterpart(:,j,:) = zero

       ENDWHERE

    ENDDO

    
  !! 3. Temperature control on decay: Factor between 0 and 1

    !! 3.1 above: surface temperature
    control_temp(:,iabove) = control_temp_func (npts, tsurf)

    
    !! 3.2 below: convolution of temperature and decomposer profiles
    !!            (exponential decomposer profile supposed)
   
    !! 3.2.1 rpc is an integration constant such that the integral of the root profile is 1.
    rpc(:) = un / ( un - EXP( -z_soil(nbdl) / z_decomp ) )

    !! 3.2.2 integrate over the nbdl levels
    tsoil_decomp(:) = zero

    DO l = 1, nbdl

       tsoil_decomp(:) = &
            tsoil_decomp(:) + tsoil(:,l) * rpc(:) * &
            ( EXP( -z_soil(l-1)/z_decomp ) - EXP( -z_soil(l)/z_decomp ) )

    ENDDO

    control_temp(:,ibelow) = control_temp_func (npts, tsoil_decomp)

  !! 4. Moisture control. Factor between 0 and 1
    
    !! 4.1 above the ground: litter humidity
    control_moist(:,iabove) = control_moist_func (npts, litterhum)

    !
    !! 4.2 below: convolution of humidity and decomposer profiles
    !            (exponential decomposer profile supposed)

    !! 4.2.1 rpc is an integration constant such that the integral of the root profile is 1.
    rpc(:) = un / ( un - EXP( -z_soil(nbdl) / z_decomp ) )

    !! 4.2.2 integrate over the nbdl levels
    soilhum_decomp(:) = zero

    DO l = 1, nbdl !Loop over soil levels

       soilhum_decomp(:) = &
            soilhum_decomp(:) + soilhum(:,l) * rpc(:) * &
            ( EXP( -z_soil(l-1)/z_decomp ) - EXP( -z_soil(l)/z_decomp ) )

    ENDDO

    control_moist(:,ibelow) = control_moist_func (npts, soilhum_decomp)

  !! 5. fluxes from litter to carbon pools and respiration

    DO l = 1, nlevs !Loop over litter levels (above and below ground)
       DO m = 2,nvm !Loop over PFTs

          !! 5.1 structural litter: goes into active and slow carbon pools + respiration

          !! 5.1.1 total quantity of structural litter which is decomposed
          fd(:) = dt/litter_tau(istructural) * &
               control_temp(:,l) * control_moist(:,l) * exp( -litter_struct_coef * lignin_struc(:,m,l) )

          qd(:) = litter(:,istructural,m,l) * fd(:)

          litter(:,istructural,m,l) = litter(:,istructural,m,l) - qd(:)

          !! 5.1.2 decompose same fraction of structural part of dead leaves. Not exact
          !!       as lignine content is not the same as that of the total structural litter.
          ! to avoid a multiple (for ibelow and iabove) modification of dead_leaves,
          ! we do this test to do this calcul only ones in 1,nlev loop
          if (l == iabove)  dead_leaves(:,m,istructural) = dead_leaves(:,m,istructural) * ( un - fd(:) )

          !! 5.1.3 non-lignin fraction of structural litter goes into
          !!       active carbon pool + respiration
          soilcarbon_input(:,iactive,m) = soilcarbon_input(:,iactive,m) + &
               frac_soil(istructural,iactive,l) * qd(:) * ( 1. - lignin_struc(:,m,l) ) / dt

      !BE CAREFUL: Here resp_hetero_litter is divided by dt to have a value which corresponds to
      ! the sechiba time step but then in stomate.f90 resp_hetero_litter is multiplied by dt.
      ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day
      ! time step and avoid some constructions that could create bug during future developments.
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(istructural,iactive,l) ) * qd(:) * &
               ( 1. - lignin_struc(:,m,l) ) / dt

          !! 5.1.4 lignin fraction of structural litter goes into
          !!       slow carbon pool + respiration
          soilcarbon_input(:,islow,m) = soilcarbon_input(:,islow,m) + &
               frac_soil(istructural,islow,l) * qd(:) * lignin_struc(:,m,l) / dt

      !BE CAREFUL: Here resp_hetero_litter is divided by dt to have a value which corresponds to
      ! the sechiba time step but then in stomate.f90 resp_hetero_litter is multiplied by dt.
      ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day
      ! time step and avoid some constructions that could create bug during future developments.
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(istructural,islow,l) ) * qd(:) * lignin_struc(:,m,l) / dt

          
          !! 5.2 metabolic litter goes into active carbon pool + respiration
         
          !! 5.2.1 total quantity of metabolic litter that is decomposed
          fd(:) = dt/litter_tau(imetabolic) * control_temp(:,l) * control_moist(:,l)

          qd(:) = litter(:,imetabolic,m,l) * fd(:)

          litter(:,imetabolic,m,l) = litter(:,imetabolic,m,l) - qd(:)

          !! 5.2.2 decompose same fraction of metabolic part of dead leaves.
          !  to avoid a multiple (for ibelow and iabove) modification of dead_leaves,
          !  we do this test to do this calcul only ones in 1,nlev loop
          if (l == iabove)  dead_leaves(:,m,imetabolic) = dead_leaves(:,m,imetabolic) * ( 1. - fd(:) )


          !! 5.2.3 put decomposed litter into carbon pool + respiration
          soilcarbon_input(:,iactive,m) = soilcarbon_input(:,iactive,m) + &
               frac_soil(imetabolic,iactive,l) * qd(:) / dt

           !!!!Chloe ++
           !IF (m .EQ. 14) THEN 
           ! soilcarbon_input(:,iactive,m) = soilcarbon_input(:,iactive,m) - wpro_conso_carb(:)*dt
           !ENDIF
           !!!!Chloe --

      !BE CAREFUL: Here resp_hetero_litter is divided by dt to have a value which corresponds to
      ! the sechiba time step but then in stomate.f90 resp_hetero_litter is multiplied by dt.
      ! Perhaps it could be simplified. Moreover, we must totally adapt the routines to the dtradia/one_day
      ! time step and avoid some constructions that could create bug during future developments.
          resp_hetero_litter(:,m) = resp_hetero_litter(:,m) + &
               ( 1. - frac_soil(imetabolic,iactive,l) ) * qd(:) / dt

       ENDDO
    ENDDO

    
  !! 6. calculate fraction of total soil covered by dead leaves

    CALL deadleaf (npts, veget_max, dead_leaves, deadleaf_cover)

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving littercalc'

  END SUBROUTINE littercalc


!! ==============================================================================================================================\n
!! SUBROUTINE   : deadleaf
!!
!>\BRIEF        This routine calculates the deadleafcover. 
!!
!! DESCRIPTION  : It first calculates the lai corresponding to the dead leaves (LAI) using 
!! the dead leaves carbon content (DL) the specific leaf area (sla) and the 
!! maximal coverage fraction of a PFT (vegetmax) using the following equations:
!! \latexonly
!! \input{deadleaf1.tex}
!! \endlatexonly
!! \n
!! Then, the dead leaf cover (DLC) is calculated as following:\n
!! \latexonly
!! \input{deadleaf2.tex}
!! \endlatexonly
!! \n
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE: ::deadleaf_cover
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE deadleaf (npts, veget_max, dead_leaves, deadleaf_cover)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: npts           !! Domain size - number of grid pixels (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)  :: dead_leaves    !! Dead leaves per ground unit area, per PFT, 
                                                                          !! metabolic and structural  
                                                                          !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)          :: veget_max      !! PFT "Maximal" coverage fraction of a PFT defined in 
                                                                          !! the input vegetation map 
                                                                          !! @tex $(m^2 m^{-2})$ @endtex 
    
    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts), INTENT(out)           :: deadleaf_cover !! Fraction of soil covered by dead leaves over all PFTs
                                                                          !! (0-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                        :: dead_lai       !! LAI of dead leaves @tex $(m^2 m^{-2})$ @endtex
    INTEGER(i_std)                                      :: j              !! Index (unitless)
!_ ================================================================================================================================
    
  !! 1. LAI of dead leaves
  
    dead_lai(:) = zero

    DO j = 2,nvm !Loop over PFTs
       dead_lai(:) = dead_lai(:) + ( dead_leaves(:,j,imetabolic) + dead_leaves(:,j,istructural) ) * sla(j) &
            * veget_max(:,j)
    ENDDO

  !! 2. fraction of soil covered by dead leaves

    deadleaf_cover(:) = un - exp( - 0.5 * dead_lai(:) )

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving deadleaf'

  END SUBROUTINE deadleaf


!! ================================================================================================================================
!! FUNCTION     : control_moist_func
!!
!>\BRIEF        Calculate moisture control for litter and soild C decomposition
!!
!! DESCRIPTION  : Calculate moisture control factor applied
!! to litter decomposition and to soil carbon decomposition in
!! stomate_soilcarbon.f90 using the following equation: \n
!! \latexonly
!! \input{control_moist_func1.tex}
!! \endlatexonly
!! \n
!! with M the moisture control factor and soilmoisutre, the soil moisture 
!! calculated in sechiba.
!! Then, the function is ranged between 0.25 and 1:\n
!! \latexonly
!! \input{control_moist_func2.tex}
!! \endlatexonly
!! \n
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE : ::moistfunc_result
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  FUNCTION control_moist_func (npts, moist_in) RESULT (moistfunc_result)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
          
    INTEGER(i_std), INTENT(in)               :: npts                !! Domain size - number of grid pixel (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in) :: moist_in            !! relative humidity (unitless)

    !! 0.2 Output variables
   
    REAL(r_std), DIMENSION(npts)             :: moistfunc_result    !! Moisture control factor (0.25-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    moistfunc_result(:) = -moist_coeff(1) * moist_in(:) * moist_in(:) + moist_coeff(2)* moist_in(:) - moist_coeff(3)
    moistfunc_result(:) = MAX( 0.25_r_std, MIN( un, moistfunc_result(:) ) )

  END FUNCTION control_moist_func


!! ================================================================================================================================
!! FUNCTION     : control_temp_func
!!
!>\BRIEF        Calculate temperature control for litter and soild C decomposition
!!
!! DESCRIPTION  : Calculate temperature control factor applied
!! to litter decomposition and to soil carbon decomposition in
!! stomate_soilcarbon.f90 using the following equation: \n
!! \latexonly
!! \input{control_temp_func1.tex}
!! \endlatexonly
!! \n
!! with T the temperature control factor, temp the temperature in Kelvin of 
!! the air (for aboveground litter) or of the soil (for belowground litter 
!! and soil)
!! Then, the function is limited in its maximal range to 1:\n
!! \latexonly
!! \input{control_temp_func2.tex}
!! \endlatexonly
!! \n
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE: ::tempfunc_result
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  FUNCTION control_temp_func (npts, temp_in) RESULT (tempfunc_result)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                 :: npts            !! Domain size - number of land pixels (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)   :: temp_in         !! Temperature (K)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts)               :: tempfunc_result !! Temperature control factor (0-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    tempfunc_result(:) = exp( soil_Q10 * ( temp_in(:) - (ZeroCelsius+tsoil_ref)) / Q10 )
    tempfunc_result(:) = MIN( un, tempfunc_result(:) )

  END FUNCTION control_temp_func

END MODULE stomate_litter
