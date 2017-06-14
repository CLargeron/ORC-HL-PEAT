!! This module define solar 
!!
!! @call sechiba_main
!! @Version : $Revision: 861 $, $Date: 2012-05-04 16:28:38 +0200 (Fri, 04 May 2012) $
!! 
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_global/solar.f90 $
!< $Date: 2012-05-04 16:28:38 +0200 (Fri, 04 May 2012) $
!< $Author: nicolas.vuichard $
!< $Revision: 861 $
!!
!! @author Marie-Alice Foujols, Jan Polcher and Martial Mancip
!! 
!!
MODULE solar

  USE defprec
  USE constantes
  USE parallel

  IMPLICIT NONE



CONTAINS
  !
  !f90doc CONTAINS
  ! 
  !
  !
  !-
  !=====================================================================
  SUBROUTINE solarang (julian, julian0, iim, jjm, lon, lat, csang)
    !---------------------------------------------------------------------
    !- This subroutine computes the solar angle according to the method
    !- used by GSWP and developed by J.C. Morill.
    !- See for further details :
    !- http://www.atmo.arizona.edu/~morrill/swrad.html
    !---------------------------------------------------------------------
    !
    USE calendar
    !
    IMPLICIT NONE
    !-
    INTEGER,INTENT(in) :: iim, jjm
    REAL,INTENT(in)  :: julian
    REAL,INTENT(in)  :: julian0
    REAL,DIMENSION(iim,jjm), INTENT(in)  :: lon, lat
    REAL,DIMENSION(iim,jjm), INTENT(out) :: csang
    !-
    REAL :: gamma, dec, decd
    REAL :: et, gmt, le, ls, lcorr, latime, omega, omegad
    REAL :: llatd, llat
    INTEGER :: igmt, ilon, ilat
    INTEGER,SAVE,ALLOCATABLE :: zone(:)
    REAL,SAVE,ALLOCATABLE :: lhour(:)
    !
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    !
    IF (check) WRITE(numout,*) 'We get the right calendar information'
    !-
    !- 1) Day angle gamma
    !-
    !   gamma = 2.*pi*MOD(julian,one_year)/one_year
    gamma = 2.*pi*(julian-julian0)/one_year
    !-
    !- 2) Solar declination (assumed constant for a 24 hour period)  page 7
    !-    in radians
    !-
    IF (check) WRITE(numout,*) 'Solar declination'
    !
    dec = ( 0.006918-0.399912*COS(gamma)+0.070257*SIN(gamma) &
         &       -0.006758*COS(2*gamma)+0.000907*SIN(2*gamma)      &
         &       -0.002697*COS(3*gamma)+0.00148*SIN(3*gamma))
    decd = dec*(180/pi)
    !-
    !- maximum error 0.0006 rad (<3'), leads to error
    !- of less than 1/2 degree in zenith angle
    !-
    IF (check) WRITE(numout,*) 'Equation of time'
    !- 3)  Equation of time  page 11
    !-
    et = ( 0.000075+0.001868*COS(gamma)-0.032077*SIN(gamma)&
         &      -0.014615*COS(2*gamma)-0.04089*SIN(2*gamma))*229.18
    !-
    !- Get the time zones for the current time
    !-
    gmt = 24.*(julian-INT(julian))
    IF (.NOT.ALLOCATED(zone))  ALLOCATE(zone(iim))
    IF (.NOT.ALLOCATED(lhour)) ALLOCATE(lhour(iim))
    !-
    !igmt = NINT(gmt)
    IF (check) WRITE(numout,*) 'Get time zone'
    CALL time_zone(gmt, iim, jjm, lon, zone, lhour)
    !-
    !- Start a loop through the grid
    !-
    IF (check) WRITE(numout,*) 'Start a loop through the grid'
    DO ilon=1,iim
       !---
       !--- 4) Local apparent time  page 13
       !---
       !--- ls     standard longitude (nearest 15 degree meridian)
       !--- le     local longtitude
       !--- lhour  local standard time
       !--- latime local apparent time
       !--- lcorr  longitudunal correction (minutes)
       !---
       le = lon(ilon,1)
       ls = ((zone(ilon)-1)*15)-180.
       lcorr = 4.*(ls-le)*(-1)
       latime = lhour(ilon)+lcorr/60.+et/60.
       IF (latime <  0.) latime = latime+24
       IF (latime > 24.) latime = latime-24
       !---
       !--- 5) Hour angle omega  page 15
       !---
       !--- hour angle is zero at noon, positive in the morning
       !--- It ranges from 180 to -180
       !---
       !--- omegad is omega in degrees, omega is in radians
       !---
       omegad = (latime-12.)*(-15.)
       omega  = omegad*pi/180.
       !---
       DO ilat=1,jjm
          !-----
          !----- 6)  Zenith angle  page 15
          !-----
          !----- csang cosine of zenith angle (radians)
          !----- llatd =  local latitude in degrees
          !----- llat  =  local latitude in radians
          !-----
          llatd = lat(1,ilat)
          llat  = llatd*pi/180.
          csang(ilon,ilat) = &
               &  MAX(zero,SIN(dec)*SIN(llat)+COS(dec)*COS(llat)*COS(omega))
       ENDDO
    ENDDO
    !----------------------
  END SUBROUTINE solarang
  !-
  !=====================================================================
  SUBROUTINE time_zone (gmt, iim, jjm, lon, zone, lhour)
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: iim, jjm, zone(iim)
    REAL :: lon(iim,jjm), gmt, lhour(iim)
    !-
    INTEGER :: deg
    !!??   REAL :: deg
    INTEGER :: i, ilon, change
    !---------------------------------------------------------------------
    DO ilon=1,iim
       !---
       !--- Convert longitude index to longtitude in degrees
       !---
       deg = lon(ilon,1)
       !---
       !--- Determine into which time zone (15 degree interval) the
       !--- longitude falls.
       !---
       DO i=1,25
          IF (deg < (-187.5+(15*i))) THEN
             zone(ilon) = i
             IF (zone(ilon) == 25)   zone(ilon) = 1
             EXIT
          ENDIF
       ENDDO
       !---
       !--- Calculate change (in number of hours) from GMT time to
       !--- local hour.  Change will be negative for zones < 13 and
       !--- positive for zones > 13.
       !---
       !--- There is also a correction for lhour < 0 and lhour > 23
       !--- to lhour between 0 and 23.
       !---
       IF (zone(ilon) < 13) THEN
          change = zone(ilon)-13
          lhour(ilon) = gmt+change
       ENDIF
       !---
       IF (zone(ilon) == 13) THEN
          lhour(ilon) = gmt
       ENDIF
       !---
       IF (zone(ilon) > 13) THEN
          change = zone(ilon)-13
          lhour(ilon) = gmt+change
       ENDIF
       IF (lhour(ilon) <  0) lhour(ilon) = lhour(ilon)+24
       IF (lhour(ilon) >= 24) lhour(ilon) = lhour(ilon)-24
       !---
    ENDDO
    !-----------------------
  END SUBROUTINE time_zone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END MODULE solar
