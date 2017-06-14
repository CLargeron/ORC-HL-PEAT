! Definition and allocation of parallel datas.
! Initialization of parallel or sequentiel IOs.
! Definition of Load Balancing functions.

!-
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parallel/data_para.f90 $ 
!< $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!< $Author: didier.solyga $
!< $Revision: 947 $
!-

MODULE data_para

!-
  USE defprec
  USE ioipsl
!-
#include "src_parallel.h"
!-
!-
! Unit for output messages
  INTEGER(i_std), SAVE :: numout = 6

  INTEGER, SAVE :: mpi_size                                           !! Number of parallel processes
  INTEGER, SAVE :: mpi_rank                                           !! my rank num
  INTEGER, SAVE :: root_prc                                           !! rank of root proc
  LOGICAL, SAVE :: is_root_prc                                        !! Only root proc is true
  INTEGER, SAVE :: nbp_loc                                            !! number of local continental points
  INTEGER, SAVE :: nbp_glo                                            !! number of global continental points
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: nbp_para_nb
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: nbp_para_begin
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: nbp_para_end
  
  INTEGER,SAVE :: iim_g, jjm_g                                  !! Dimension of global fields
  ! i x j 2D points (not land points) index
  INTEGER,SAVE,allocatable,dimension(:) :: ij_para_nb           ! Number of 2D points for each mpi_rank block
  INTEGER,SAVE,allocatable,dimension(:) :: ij_para_begin        ! First 2D point for each mpi_rank block
  INTEGER,SAVE,allocatable,dimension(:) :: ij_para_end          ! Last 2D point for each mpi_rank block
  ! i 2D index 
  INTEGER,SAVE,allocatable,dimension(:) :: ii_para_begin        ! First i index of 2D point for each mpi_rank block
  INTEGER,SAVE,allocatable,dimension(:) :: ii_para_end          ! Last i index of 2D point for each mpi_rank block
  ! j 2D index
  INTEGER,SAVE,allocatable,dimension(:) :: jj_para_nb           ! Number of complete j lines for each mpi_rank block
  INTEGER,SAVE,allocatable,dimension(:) :: jj_para_begin        ! First j index of 2D point for each mpi_rank block
  INTEGER,SAVE,allocatable,dimension(:) :: jj_para_end          ! Last j index of 2D point for each mpi_rank block
  
  INTEGER,SAVE :: ii_begin
  INTEGER,SAVE :: ii_end
  INTEGER,SAVE :: jj_begin
  INTEGER,SAVE :: jj_end
  INTEGER,SAVE :: jj_nb
  INTEGER,SAVE :: ij_begin
  INTEGER,SAVE :: ij_end
  INTEGER,SAVE :: ij_nb
 
  
  !! Global array used by stomate and sechiba
  !-
  !! index of land points on the 2D map
  INTEGER(i_std),ALLOCATABLE,DIMENSION(:),SAVE   :: index_g
  !-
  !! indices of the 4 neighbours of each grid point (1=N, 2=E, 3=S, 4=W)
  INTEGER(i_std),ALLOCATABLE,DIMENSION(:,:),SAVE :: neighbours_g
  !-
  !! resolution at each grid point in m (1=E-W, 2=N-S)
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:),SAVE    :: resolution_g 
  REAL(r_std),ALLOCATABLE,DIMENSION(:),SAVE    :: area_g 
  !-
  !! Geographical coordinates
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:),SAVE    :: lalo_g
  ! Global grid, for all process
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE     :: lon_g, lat_g, zlev_g
  !-
  !! Fraction of continents
  REAL(r_std),ALLOCATABLE,DIMENSION(:),SAVE      :: contfrac_g
 
  INTEGER, SAVE :: MPI_COMM_ORCH
  INTEGER, SAVE :: MPI_REAL_ORCH
  INTEGER, SAVE :: MPI_INT_ORCH
  LOGICAL, SAVE :: cpl_lmdz
  
  INTEGER, SAVE :: orch_domain_id
  
CONTAINS
  
  SUBROUTINE init_para(cpl_lmdz_x, communicator)

    IMPLICIT NONE
    
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    LOGICAL :: cpl_lmdz_x
    INTEGER,OPTIONAL :: communicator
    INTEGER :: ierr
    CHARACTER(LEN=20) :: filename
    INTEGER :: i, div
    CHARACTER(LEN=4) :: num
    LOGICAL, PARAMETER :: check = .FALSE.

    cpl_lmdz=cpl_lmdz_x
      
#ifdef CPP_PARA
    ! Orchidee communicator :
    IF (.NOT. cpl_lmdz) THEN
       CALL MPI_INIT(ierr)
       IF (ierr /= 0) THEN
          WRITE(*,*) 'INIT_PARA : MPI_INIT failed with ',ierr
          STOP "INIT_PARA"          
       ENDIF
       MPI_COMM_ORCH=MPI_COMM_WORLD
    ELSE
       MPI_COMM_ORCH=communicator
    ENDIF
	  
    
    ! Int and real precision
    IF (MPI_VERSION == 1) THEN
       ! Version MPI 1.x
       IF (i_std==i_4) THEN
          MPI_INT_ORCH=MPI_INTEGER
       ELSEIF (i_std==i_8) THEN
          MPI_INT_ORCH=MPI_INTEGER
       ELSE
          MPI_INT_ORCH=MPI_INTEGER
       ENDIF
       
       IF (r_std==r_4) THEN
          MPI_REAL_ORCH=MPI_REAL
       ELSEIF (r_std==r_8) THEN
          MPI_REAL_ORCH=MPI_DOUBLE_PRECISION
       ELSE
          MPI_REAL_ORCH=MPI_REAL
       ENDIF
    ELSE
       ! Others MPI
       IF (i_std==i_4) THEN
          MPI_INT_ORCH=MPI_INTEGER4
       ELSEIF (i_std==i_8) THEN
          MPI_INT_ORCH=MPI_INTEGER8
       ELSE
          MPI_INT_ORCH=MPI_INTEGER
       ENDIF
         
       IF (r_std==r_4) THEN
          MPI_REAL_ORCH=MPI_REAL4
       ELSEIF (r_std==r_8) THEN
          MPI_REAL_ORCH=MPI_REAL8
       ELSE
          MPI_REAL_ORCH=MPI_REAL
       ENDIF
    ENDIF
      
    CALL MPI_COMM_SIZE(MPI_COMM_ORCH,mpi_size,ierr)
    IF (ierr /= 0) THEN
       WRITE(*,*) 'INIT_PARA : MPI_COMM_SIZE failed with ',ierr
       STOP "INIT_PARA"
    ENDIF
    CALL MPI_COMM_RANK(MPI_COMM_ORCH,mpi_rank,ierr)
    IF (ierr /= 0) THEN
       WRITE(*,*) 'INIT_PARA : MPI_COMM_RANK failed with ',ierr
       STOP "INIT_PARA"
    ENDIF
      
#else
    mpi_rank=0
    mpi_size=1    
#endif
    root_prc=0
      
    IF (mpi_rank==0) THEN
       is_root_prc=.TRUE.
    ELSE
       is_root_prc=.FALSE.
    ENDIF

    ! Open mpi_rank outputs or select stdout 
    IF (mpi_size > 1) THEN
       write(num,'(I4.4)') mpi_rank
       
       numout = 100
       
       filename =  'out_orchidee_'//num
       
       OPEN(UNIT=numout,FILE=TRIM(filename),ACTION='write',STATUS='unknown',FORM='formatted',IOSTAT=ierr) 
       IF (ierr /= 0) THEN
#ifdef CPP_PARA
          CALL MPI_FINALIZE(ierr)
#endif
          WRITE(*,*) "Erreur can't open file ", filename
          STOP "INIT_PARA"
       ENDIF
    ELSE
       numout = 6
    ENDIF

    IF (check) THEN
#ifdef CPP_PARA
       WRITE(numout,*) 'version MPI ', MPI_VERSION, MPI_SUBVERSION
       WRITE(numout,*) 'INTEGERS ',MPI_INTEGER4,MPI_INTEGER8,MPI_INTEGER,MPI_INT_ORCH
       WRITE(numout,*) 'REALS ',MPI_REAL4,MPI_REAL8,MPI_REAL,MPI_REAL_ORCH
#endif
       WRITE(numout,*) 'RANK ',mpi_rank,' SIZE ',mpi_size
       WRITE(numout,*) "Am I root process ?",is_root_prc
       WRITE(numout,*) "Init_para : For process number ",mpi_rank, "output file = ",filename
    ENDIF
      
  END SUBROUTINE init_para
 
  SUBROUTINE init_data_para(iim,jjm,nbpoints,index_x)

    IMPLICIT NONE
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif

    INTEGER,INTENT(IN) :: iim
    INTEGER,INTENT(IN) :: jjm
    INTEGER,INTENT(IN) :: nbpoints
    INTEGER,DIMENSION(:),INTENT(IN) :: index_x

    INTEGER,ALLOCATABLE,DIMENSION(:) :: index_l
    INTEGER,ALLOCATABLE,DIMENSION(:) :: displs
#ifdef CPP_PARA
    INTEGER :: ierr
#endif  
    INTEGER :: i
    LOGICAL, PARAMETER :: check=.FALSE.


    IF (check) WRITE(numout,*) 'INIT_DATA_PARA',iim,jjm,nbpoints,index_x
    ALLOCATE(nbp_para_nb(0:mpi_size-1))
    ALLOCATE(nbp_para_begin(0:mpi_size-1))
    ALLOCATE(nbp_para_end(0:mpi_size-1))
    ALLOCATE(jj_para_nb(0:mpi_size-1))
    ALLOCATE(jj_para_begin(0:mpi_size-1))
    ALLOCATE(jj_para_end(0:mpi_size-1))
    ALLOCATE(ii_para_begin(0:mpi_size-1))
    ALLOCATE(ii_para_end(0:mpi_size-1))
    ALLOCATE(ij_para_nb(0:mpi_size-1))
    ALLOCATE(ij_para_begin(0:mpi_size-1))
    ALLOCATE(ij_para_end(0:mpi_size-1))

    IF (cpl_lmdz) THEN
#ifdef CPP_PARA
       CALL MPI_GATHER(nbpoints,1,MPI_INT_ORCH,nbp_para_nb,1,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
#else
       nbp_para_nb(0)=nbpoints
#endif

       IF (is_root_prc) THEN
          nbp_glo=sum(nbp_para_nb)
          ALLOCATE(displs(0:mpi_size-1))
	  ALLOCATE(index_l(nbp_glo))
          displs(0)=0
          DO i=1,mpi_size-1
             displs(i)=displs(i-1)+nbp_para_nb(i-1)
          ENDDO
       ENDIF
#ifdef CPP_PARA
       CALL MPI_BCAST(nbp_glo,1,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
       CALL MPI_GATHERV(index_x,nbpoints,MPI_INT_ORCH,index_l,nbp_para_nb,displs,&
            MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)  
#else
       ! if not parallized, nb_glo=nbpoints=nbp_para_nb
       IF (is_root_prc) index_l(:)=index_x(1:nbp_glo)
#endif
    ELSE
       IF (is_root_prc) THEN
          nbp_glo=nbpoints
          ALLOCATE(index_l(nbp_glo))
          index_l(:)=index_x(1:nbp_glo)
          CALL Read_Load_Balance(nbpoints,nbp_para_nb)
       ENDIF
    ENDIF

    IF (is_root_prc) THEN

       IF (check) WRITE(numout,*) '==== DISTRIB ====',nbpoints,nbp_para_nb

       nbp_para_begin(0)=1
       nbp_para_end(0)=nbp_para_nb(0)

       DO i=1,mpi_size-1
          nbp_para_begin(i)=nbp_para_end(i-1)+1
          nbp_para_end(i)=nbp_para_begin(i)+nbp_para_nb(i)-1
       ENDDO

       nbp_loc=nbp_para_nb(mpi_rank)

       iim_g=iim
       jjm_g=jjm


       ij_para_begin(0)=1
       ij_para_end(0)=index_l(nbp_para_end(0))
       ij_para_nb(0)=ij_para_end(0)-ij_para_begin(0)+1

       DO i=1,mpi_size-1
          ij_para_begin(i)=ij_para_end(i-1)+1
          ij_para_end(i)=index_l(nbp_para_end(i))
	  ij_para_nb(i)=ij_para_end(i)-ij_para_begin(i)+1
       ENDDO

       ij_para_end(mpi_size-1)=iim*jjm
       ij_para_nb(mpi_size-1)=ij_para_end(mpi_size-1)-ij_para_begin(mpi_size-1)+1


       DO i=0,mpi_size-1
          jj_para_begin(i)=(ij_para_begin(i)-1)/iim + 1
	  jj_para_end(i)=(ij_para_end(i)-1)/iim + 1
          jj_para_nb(i)=jj_para_end(i)-jj_para_begin(i)+1

	  ii_para_begin(i)=MOD(ij_para_begin(i)-1,iim)+1
	  ii_para_end(i)=MOD(ij_para_end(i)-1,iim)+1
       ENDDO

       IF (check) THEN
          WRITE(numout,*) '==== DECOUP ===='
          WRITE(numout,*) 'nbp_para_begin=',nbp_para_begin
          WRITE(numout,*) 'nbp_para_end  =',nbp_para_end
          WRITE(numout,*) 'nbp_loc=',nbp_loc

          WRITE(numout,*) 'ij_para_begin=',ij_para_begin
          WRITE(numout,*) 'ij_para_end=',ij_para_end
          WRITE(numout,*) 'ij_para_nb=',ij_para_nb
          WRITE(numout,*) 'jj_para_begin=',jj_para_begin
          WRITE(numout,*) 'jj_para_end=',jj_para_end
          WRITE(numout,*) 'jj_para_nb=',jj_para_nb	
          WRITE(numout,*) 'ii_para_begin=',ii_para_begin
          WRITE(numout,*) 'ii_para_end=',ii_para_end
       ENDIF

       !
       !- Root need the global variables
       !-
       ALLOCATE(resolution_g(nbp_glo,2),area_g(nbp_glo),lalo_g(nbp_glo,2))
       ALLOCATE(neighbours_g(nbp_glo,8),contfrac_g(nbp_glo),index_g(nbp_glo))
       ALLOCATE(lon_g(iim_g, jjm_g), lat_g(iim_g, jjm_g), zlev_g(iim_g, jjm_g))

       index_g(:)=index_l(1:nbp_glo)
    ELSE
       ALLOCATE(resolution_g(0,0),area_g(0),lalo_g(0,0))
       ALLOCATE(neighbours_g(0,0),contfrac_g(0),index_g(0))
    ENDIF

#ifdef CPP_PARA
    IF (is_root_prc) WRITE(numout,*) 'nbp_para_nb =',nbp_para_nb
    CALL MPI_BCAST(nbp_para_nb,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(nbp_para_begin,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(nbp_para_end,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(jj_para_nb,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(jj_para_begin,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(jj_para_end,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(ii_para_begin,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(ii_para_end,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(ij_para_nb,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(ij_para_begin,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(ij_para_end,mpi_size,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(iim_g,1,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(jjm_g,1,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    CALL MPI_BCAST(nbp_glo,1,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
#endif

    nbp_loc=nbp_para_nb(mpi_rank)

    ij_nb=ij_para_nb(mpi_rank)
    ij_begin=ij_para_begin(mpi_rank)
    ij_end=ij_para_end(mpi_rank)

    jj_nb=jj_para_nb(mpi_rank)
    jj_begin=jj_para_begin(mpi_rank)
    jj_end=jj_para_end(mpi_rank)

    ii_begin=ii_para_begin(mpi_rank)
    ii_end=ii_para_end(mpi_rank)

    IF (check) &
         WRITE(numout,*) 'Init_io_para'
    CALL Init_io_para

    IF (is_root_prc ) THEN 
       DEALLOCATE(index_l)
       IF ( cpl_lmdz) THEN
          DEALLOCATE(displs)
       ENDIF
    ENDIF

    IF (check) &
         WRITE(numout,*)  'DATA PARA',nbp_loc,nbp_glo,jj_begin,jj_end,ii_begin,ii_end       
  END SUBROUTINE init_data_para
  

  SUBROUTINE Init_io_para
    USE ioipsl
    IMPLICIT NONE
    
    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 

    ddid=(/ 1,2 /)
    dsg=(/ iim_g, jjm_g /)
    dsl=(/ iim_g, jj_nb /)
    dpf=(/ 1,jj_begin /)
    dpl=(/ iim_g, jj_end /)
    dhs=(/ ii_begin-1,0 /)
    if (mpi_rank==mpi_size-1) then
      dhe=(/0,0/)
    else
      dhe=(/ iim_g-ii_end,0 /)  
    endif
    
    call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'APPLE',orch_domain_id)
  END SUBROUTINE Init_io_para
  
  SUBROUTINE Read_Load_balance(NbPoints,Nbpoints_loc)

    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: NbPoints
    INTEGER,INTENT(OUT) :: Nbpoints_loc(0:mpi_size-1)
#ifdef CPP_PARA  
    INTEGER :: unit_number=10
    CHARACTER(len=255)  :: filename='Load_balance_orchidee.dat'
    INTEGER :: j
#endif
    INTEGER :: i,s,ierr
    
    Nbpoints_loc(:) = 0

#ifdef CPP_PARA
    OPEN(UNIT=unit_number,FILE=trim(filename),STATUS='old',FORM='formatted',IOSTAT=ierr) 
#else
    ierr=1
#endif   
    s=0

#ifdef CPP_PARA  
    IF (ierr==0) THEN
       i=0
       !- Reading for any balancing file (even with a bad structure)
       DO WHILE (i < mpi_size .AND. ierr == 0) 
          READ (unit_number,*,IOSTAT=ierr) j,Nbpoints_loc(i)
          s=s+Nbpoints_loc(i)
          i=i+1
       ENDDO
       CLOSE(unit_number)
    ENDIF
#endif   
    
    !- Correction of bad balancing file (or an empty file) => same nb of points for each procs
    IF (ierr/=0 .OR. s/=Nbpoints) THEN
       DO i=0,mpi_size-1
          Nbpoints_loc(i)=Nbpoints/mpi_size
          IF (MOD(Nbpoints,mpi_size) > i) Nbpoints_loc(i)=Nbpoints_loc(i)+1
       ENDDO
    ENDIF
    
  END SUBROUTINE Read_Load_balance
  
  SUBROUTINE Write_Load_balance(times)
    IMPLICIT NONE
    REAL(r_std),INTENT(IN) :: times
  
#ifdef CPP_PARA  
    CHARACTER(len=255)  :: filename='Load_balance_orchidee.dat'
    INTEGER :: unit_number=10
    INTEGER :: i,ierr
    REAL :: All_Times(0:mpi_size-1)
    REAL :: average
    REAL :: efficiency
    INTEGER :: dp,S
    INTEGER :: New_nbpoints(0:mpi_size-1)
#endif
    
    WRITE(numout,*) 'time',times
#ifndef CPP_PARA
    RETURN
#else

    CALL MPI_GATHER(times,1,MPI_REAL_ORCH,All_times,1,MPI_REAL_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    IF (is_root_prc) WRITE(numout,*) 'ALL_times',All_times

    IF (is_root_prc) THEN
     
       OPEN(UNIT=unit_number,FILE=trim(filename),STATUS='replace',FORM='formatted',IOSTAT=ierr)
       
       average=sum(All_times(:))/mpi_size
       DO i=0,mpi_size-1
          efficiency=All_times(i)/Nbp_para_nb(i)
          New_nbpoints(i)=Nbp_para_nb(i)-(All_times(i)-average)/efficiency
       ENDDO
       
       S=sum(new_nbpoints(:))
       dp=nbp_glo-S
       
       IF ( dp > 0 ) THEN
          DO WHILE ( dp > 0 )
             New_nbpoints(MOD(dp,mpi_size))=New_nbpoints(MOD(dp,mpi_size))+1
             dp=dp-1
          ENDDO
       ELSE
          dp=-dp
          DO WHILE ( dp > 0 )
             New_nbpoints(MOD(dp,mpi_size))=New_nbpoints(MOD(dp,mpi_size))-1
             dp=dp-1
          ENDDO
       ENDIF
       
       ! If this algorithm diverge, we use previous repartition.
       IF ( ANY(New_nbpoints(:) .LE. 0) ) THEN
          New_nbpoints(:)=Nbp_para_nb(:)
       ENDIF
       
       DO i=0,mpi_size-1
          WRITE(Unit_number,*) i,New_nbpoints(i)
       ENDDO
       CLOSE(Unit_number)
    ENDIF
#endif
 
  END SUBROUTINE Write_Load_Balance
  
  
  
END MODULE data_para

#include "mpi_dummy.h"
