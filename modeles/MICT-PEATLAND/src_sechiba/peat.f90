  !! ================================================================================================================================
!! SUBROUTINE   : slowproc_peat
!!
!>\BRIEF         Interpolate the peat map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): Chloe Largeron
!!
!! MAIN OUTPUT VARIABLE(S): :: peatland 

!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_peat(nbpt, lalo, neighbours, resolution, contfrac, peatland )
   
    !   This subroutine should read the Peat map and interpolate to the model grid. The method
    !   is to get fraction of peat soil for each grid box.

    !!  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)    :: nbpt                   !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)       :: lalo(nbpt,2)           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)    :: neighbours(nbpt,8)     !! Vector of neighbours for each grid point 
                                                            !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)       :: resolution(nbpt,2)     !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)       :: contfrac(nbpt)         !! Fraction of land in each grid box.
    !

    !  0.2 OUTPUT
    !
    LOGICAL, INTENT(out)      :: peatland(nbpt)     !! The fraction of clay as used by STOMATE


    !  0.3 LOCAL
    !
    CHARACTER(LEN=80) :: filepeat
   
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, fopt, ilf, nbexp
    REAL(r_std) :: lev(1), date, dt
    INTEGER(i_std) :: itau(1)
    INTEGER(i_std)    :: idi, idi_last, nbvmax
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: soiltext, soiltext2
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_rel, lon_rel
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: peatfrac
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: peatsoil

    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)  :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: resol_lu
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: solt, solt2
    REAL(r_std) ::  sgn, coslat
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std) :: nix, njx
   
    INTEGER                  :: ALLOC_ERR
    LOGICAL ::           ok_interpol = .FALSE. ! optionnal return of aggregate_2d

! ==============================================================================
    !
    !
    !  Needs to be a configurable variable
    !
    !
    !Config Key   = PEAT_FILE
    !Config Desc  = Name of file from peat soil is read
    !Config Def   = Global_Peat.nc
    !Config Units = [FILE]
    !
    
    filepeat = 'Global_Peat.nc'
    CALL getin_p('PEAT_FILE', filepeat)

    IF (is_root_prc) THEN
       CALL flininfo(filepeat,iml, jml, lml, tml, fid)
       CALL flinclo(fid)

    ENDIF
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    
    ! Global_Peat file is .05° peatland fraction file.
    !
    ALLOC_ERR=-1
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(peatfrac(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of peatfrac : ",ALLOC_ERR
      STOP 
    ENDIF
     ALLOC_ERR=-1
    ALLOCATE(peatsoil(ib), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of peatsoil : ",ALLOC_ERR
      STOP 
    ENDIF


    ALLOC_ERR=-1
    ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of resol_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    IF (is_root_prc) CALL flinopen(filepeat, .FALSE., iml, jml, lml, lon_rel, lat_rel, lev, tml, itau, date, dt, fid)
    CALL bcast(lon_rel)
    CALL bcast(lat_rel)
    CALL bcast(itau)
    CALL bcast(date)
    CALL bcast(dt)
    
    !
    IF (is_root_prc) CALL flinget(fid, 'Peatlands', iml, jml, lml, tml, 1, 1, peatfrac)
    CALL bcast(soiltext)
    !
    !
    IF (is_root_prc) CALL flinclo(fid)
    !
    nbexp = 0
    !
    !
    ! Mask of permitted variables.
    !
    mask(:,:) = zero
    DO ip=1,iml
       DO jp=1,jml
          !
          IF (peatfrac(ip,jp) .GT. -1.) THEN
             mask(ip,jp) = un
          ENDIF
          !
          ! Resolution in longitude
          !
          coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )     
          IF ( ip .EQ. 1 ) THEN
             resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip,jp) ) * pi/180. * R_Earth * coslat
          ELSEIF ( ip .EQ. iml ) THEN
             resol_lu(ip,jp,1) = ABS( lon_rel(ip,jp) - lon_rel(ip-1,jp) ) * pi/180. * R_Earth * coslat
          ELSE
             resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
          ENDIF
          !
          ! Resolution in latitude
          !
          IF ( jp .EQ. 1 ) THEN
             resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp) - lat_rel(ip,jp+1) ) * pi/180. * R_Earth
          ELSEIF ( jp .EQ. jml ) THEN
             resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp) ) * pi/180. * R_Earth
          ELSE
             resol_lu(ip,jp,2) =  ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp+1) )/2. * pi/180. * R_Earth
          ENDIF
          !
       ENDDO
    ENDDO
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
       njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
       nbvmax = nix*njx
    ENDIF
    CALL bcast(nbvmax)
    !
    callsign = "Peatlands"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
    WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx
       !
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_rel, lat_rel, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          !
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO

  nbexp = 0
    DO ib = 1, nbpt
      !-
      !- Peatland fraction
      !- 
      fopt = COUNT(sub_area(ib,:) > zero)
      !
      !    Check that we found some points
      !
      IF ( fopt .EQ. 0) THEN
         nbexp = nbexp + 1
         peatsoil(ib) = 0. 
      ELSE
         peatsoil(ib) = zero
         ! Initialize last index to the highest possible 
         idi_last=nbvmax
         DO idi=1, nbvmax
            ! Leave the do loop if all sub areas are treated, sub_area <= 0
            IF ( sub_area(ib,idi) <= zero ) THEN
               ! Set last index to the last one used
               idi_last=idi-1
               ! Exit do loop
               EXIT
            END IF

            ip = sub_index(ib,idi,1)
            jp = sub_index(ib,idi,2)

            peatsoil(ib) = peatsoil(ib) + peatfrac(ip,jp) * sub_area(ib,idi)
         ENDDO

         IF ( idi_last >= 1 ) THEN ! devrait tjs etre le cas
            peatsoil(ib) = peatsoil(ib) / SUM(sub_area(ib,1:idi_last)) 
         ELSE
            peatsoil(ib) = 0.
         ENDIF
       ENDIF
       !
    ENDDO
    !
    WRITE(numout,*) 'Interpolation Done'
    !
    IF ( nbexp .GT. 0 ) THEN
       WRITE(numout,*) 'slowproc_peat : The interpolation of the peatland fraction had ', nbexp
       WRITE(numout,*) 'slowproc_peat : points without data. This are either coastal points or'
       WRITE(numout,*) 'slowproc_peat : ice covered land.'
       WRITE(numout,*) 'slowproc_peat : The problem was solved by setting the peatland fraction to 0.'
    ENDIF
    !
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)
    !
    DEALLOCATE (resol_lu)
    DEALLOCATE (peatfrac)
    DEALLOCATE (peatsoil)
    !
    RETURN
 
 
  END SUBROUTINE slowproc_peat


