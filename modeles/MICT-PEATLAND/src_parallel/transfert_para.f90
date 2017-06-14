! Low level parallel communication encapsulations for ORCHIDEE.

!-
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parallel/transfert_para.f90 $ 
!< $Date: 2012-07-19 15:12:52 +0200 (Thu, 19 Jul 2012) $
!< $Author: didier.solyga $
!< $Revision: 947 $
!-

MODULE transfert_para

  USE data_para
  USE timer
!-
  IMPLICIT NONE
!-
#include "src_parallel.h"
!-

  INTERFACE bcast
    MODULE PROCEDURE bcast_c, bcast_c1,                           &
                     bcast_i,bcast_i1,bcast_i2,bcast_i3,bcast_i4, &
                     bcast_r,bcast_r1,bcast_r2,bcast_r3,bcast_r4, &
		     bcast_l,bcast_l1,bcast_l2,bcast_l3,bcast_l4
  END INTERFACE

  INTERFACE scatter
    MODULE PROCEDURE scatter_i,scatter_i1,scatter_i2,scatter_i3, &
                     scatter_r,scatter_r1,scatter_r2,scatter_r3, &
		     scatter_l,scatter_l1,scatter_l2,scatter_l3
  END INTERFACE

!!$  INTERFACE gather_s
!!$    MODULE PROCEDURE gather_is, &
!!$                     gather_rs, &
!!$		     gather_ls
!!$  END INTERFACE
  
  INTERFACE gather
    MODULE PROCEDURE gather_i,gather_i1,gather_i2,gather_i3, &
                     gather_r,gather_r1,gather_r2,gather_r3, &
		     gather_l,gather_l1,gather_l2,gather_l3  
  END INTERFACE
  
  INTERFACE scatter2D
    MODULE PROCEDURE scatter2D_i,scatter2D_i1,scatter2D_i2,scatter2D_i3, &
                     scatter2D_r0,scatter2D_r,scatter2D_r1,scatter2D_r2,scatter2D_r3, &
		     scatter2D_l,scatter2D_l1,scatter2D_l2,scatter2D_l3
  END INTERFACE

  INTERFACE gather2D
    MODULE PROCEDURE gather2D_i,gather2D_i1,gather2D_i2,gather2D_i3, &
                     gather2D_r0,gather2D_r,gather2D_r1,gather2D_r2,gather2D_r3, &
		     gather2D_l,gather2D_l1,gather2D_l2,gather2D_l3
  END INTERFACE 
  
  INTERFACE reduce_sum
    MODULE PROCEDURE reduce_sum_i,reduce_sum_i1,reduce_sum_i2,reduce_sum_i3,reduce_sum_i4, &
                     reduce_sum_r,reduce_sum_r1,reduce_sum_r2,reduce_sum_r3,reduce_sum_r4
  END INTERFACE 
     
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Broadcast --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Les chaine de charactère -- !!

  SUBROUTINE bcast_c(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var
    CHARACTER(LEN=len(Var)),DIMENSION(1) :: Var1    
#ifndef CPP_PARA
    RETURN
#else
    IF (is_root_prc) &
         Var1(1)=Var
    CALL bcast_cgen(Var1,1)
    Var=Var1(1)
#endif
  END SUBROUTINE bcast_c

  SUBROUTINE bcast_c1(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var(:)

#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_cgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_c1


!! -- Les entiers -- !!

  SUBROUTINE bcast_i(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var
    INTEGER,DIMENSION(1) :: Var1

#ifndef CPP_PARA
    RETURN
#else
   IF (is_root_prc) &
         Var1(1)=Var 
    CALL bcast_igen(Var1,1)
    Var=Var1(1)
#endif
  END SUBROUTINE bcast_i

  SUBROUTINE bcast_i1(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_igen(Var,size(Var))
#endif
  END SUBROUTINE bcast_i1

  SUBROUTINE bcast_i2(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_igen(Var,size(Var))
#endif
  END SUBROUTINE bcast_i2

  SUBROUTINE bcast_i3(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_igen(Var,size(Var))
#endif
  END SUBROUTINE bcast_i3

  SUBROUTINE bcast_i4(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_igen(Var,size(Var))
#endif
  END SUBROUTINE bcast_i4


!! -- Les reels -- !!

  SUBROUTINE bcast_r(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var
    REAL,DIMENSION(1) :: Var1   
#ifndef CPP_PARA
    RETURN
#else
    IF (is_root_prc) &
         Var1(1)=Var
    CALL bcast_rgen(Var1,1)
    Var=Var1(1)
#endif
  END SUBROUTINE bcast_r

  SUBROUTINE bcast_r1(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_rgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_r1

  SUBROUTINE bcast_r2(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_rgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_r2

  SUBROUTINE bcast_r3(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_rgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_r3

  SUBROUTINE bcast_r4(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_rgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_r4
  
!! -- Les booleans -- !!

  SUBROUTINE bcast_l(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var
    LOGICAL,DIMENSION(1) :: Var1
   
#ifndef CPP_PARA
    RETURN
#else
    IF (is_root_prc) &
         Var1(1)=Var
    CALL bcast_lgen(Var1,1)
    Var=Var1(1)
#endif
  END SUBROUTINE bcast_l

  SUBROUTINE bcast_l1(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_lgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_l1

  SUBROUTINE bcast_l2(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_lgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_l2

  SUBROUTINE bcast_l3(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_lgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_l3

  SUBROUTINE bcast_l4(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:,:)
   
#ifndef CPP_PARA
    RETURN
#else
    CALL bcast_lgen(Var,size(Var))
#endif
  END SUBROUTINE bcast_l4
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter_i(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(nbp_loc) :: VarOut

    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

     IF (is_root_prc) THEN
      CALL scatter_igen(VarIn,Varout,1)
     ELSE
      CALL scatter_igen(dummy,Varout,1)
    ENDIF
    
#endif
  END SUBROUTINE scatter_i

  SUBROUTINE scatter_i1(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
        
#ifdef CPP_PARA
    INTEGER :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_igen(VarIn,Varout,Size(VarOut,2))
    ELSE
      CALL scatter_igen(dummy,Varout,Size(VarOut,2))
    ENDIF
    
#endif
  END SUBROUTINE scatter_i1
  
  SUBROUTINE scatter_i2(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3))
    ELSE
      CALL scatter_igen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3))
    ENDIF
#endif
  END SUBROUTINE scatter_i2

  SUBROUTINE scatter_i3(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ELSE
      CALL scatter_igen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ENDIF
  
#endif
  END SUBROUTINE scatter_i3


  SUBROUTINE scatter_r(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_rgen(VarIn,Varout,1)
    ELSE
      CALL scatter_rgen(dummy,Varout,1)
    ENDIF
  
#endif
  END SUBROUTINE scatter_r

  SUBROUTINE scatter_r1(VarIn, VarOut)

  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
        
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_rgen(VarIn,Varout,Size(VarOut,2))
    ELSE
      CALL scatter_rgen(dummy,Varout,Size(VarOut,2))      
    ENDIF
  
#endif
  END SUBROUTINE scatter_r1
  
  SUBROUTINE scatter_r2(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3))
    ELSE
      CALL scatter_rgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3))
    ENDIF
  
#endif
  END SUBROUTINE scatter_r2

  SUBROUTINE scatter_r3(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ELSE
      CALL scatter_rgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ENDIF
  
#endif
  END SUBROUTINE scatter_r3


  SUBROUTINE scatter_l(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA    
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_lgen(VarIn,Varout,1)
    ELSE
      CALL scatter_lgen(dummy,Varout,1)
    ENDIF
    
#endif
  END SUBROUTINE scatter_l

  SUBROUTINE scatter_l1(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_lgen(VarIn,Varout,Size(VarOut,2))
    ELSE
      CALL scatter_lgen(dummy,Varout,Size(VarOut,2))      
    ENDIF
  
#endif
  END SUBROUTINE scatter_l1
  
  SUBROUTINE scatter_l2(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3))
    ELSE
      CALL scatter_lgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3))
    ENDIF
  
#endif
  END SUBROUTINE scatter_l2

  SUBROUTINE scatter_l3(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else
    IF (is_root_prc) THEN
      CALL scatter_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ELSE
      CALL scatter_lgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ENDIF
  
#endif
  END SUBROUTINE scatter_l3  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Gather   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$  SUBROUTINE gather_is(VarIn, VarOut)
!!$    USE data_para
!!$    USE timer
!!$
!!$    IMPLICIT NONE
!!$  
!!$#ifdef CPP_PARA
!!$    INCLUDE 'mpif.h'
!!$#endif
!!$    
!!$    INTEGER,INTENT(IN) :: VarIn
!!$    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
!!$  
!!$#ifdef CPP_PARA
!!$    INTEGER :: nb,i,index_para,rank
!!$    INTEGER :: ierr
!!$    LOGICAL :: flag=.FALSE.
!!$    LOGICAL, PARAMETER :: check=.FALSE.
!!$#endif
!!$
!!$#ifndef CPP_PARA
!!$    VarOut(:)=VarIn
!!$    RETURN
!!$#else
!!$
!!$    IF (timer_state(timer_mpi)==running) THEN
!!$      flag=.TRUE.
!!$    ELSE
!!$      flag=.FALSE.
!!$    ENDIF
!!$    
!!$    IF (flag) CALL suspend_timer(timer_mpi)
!!$
!!$    IF (check) &
!!$         WRITE(numout,*) "gather_rgen VarIn=",VarIn    
!!$
!!$#ifdef CPP_PARA
!!$    CALL MPI_GATHER(VarIn,1,MPI_INT_ORCH,VarOut,1,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
!!$#endif
!!$
!!$    IF (check) &
!!$         WRITE(numout,*) "gather_rgen VarOut=",VarOut
!!$    IF (flag) CALL resume_timer(timer_mpi)
!!$#endif
!!$  END SUBROUTINE gather_is
!!$
!!$  SUBROUTINE gather_rs(VarIn, VarOut)
!!$    USE data_para
!!$    USE timer
!!$
!!$    IMPLICIT NONE
!!$  
!!$#ifdef CPP_PARA
!!$    INCLUDE 'mpif.h'
!!$#endif
!!$
!!$    REAL,INTENT(IN) :: VarIn
!!$    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
!!$  
!!$#ifdef CPP_PARA
!!$    INTEGER :: nb,i,index_para,rank
!!$    INTEGER :: ierr
!!$    LOGICAL :: flag=.FALSE.
!!$    LOGICAL, PARAMETER :: check=.FALSE.
!!$#endif
!!$
!!$#ifndef CPP_PARA
!!$    VarOut(:)=VarIn
!!$    RETURN
!!$#else
!!$
!!$    IF (timer_state(timer_mpi)==running) THEN
!!$      flag=.TRUE.
!!$    ELSE
!!$      flag=.FALSE.
!!$    ENDIF
!!$    
!!$    IF (flag) CALL suspend_timer(timer_mpi)
!!$
!!$    IF (check) &
!!$         WRITE(numout,*) "gather_rgen VarIn=",VarIn    
!!$
!!$#ifdef CPP_PARA
!!$    CALL MPI_GATHER(VarIn,1,MPI_REAL_ORCH,VarOut,1,MPI_REAL_ORCH,root_prc,MPI_COMM_ORCH,ierr)
!!$#endif
!!$
!!$    IF (check) &
!!$         WRITE(numout,*) "gather_rgen VarOut=",VarOut
!!$
!!$    IF (flag) CALL resume_timer(timer_mpi)
!!$#endif
!!$  END SUBROUTINE gather_rs
!!$
!!$  SUBROUTINE gather_ls(VarIn, VarOut)
!!$    USE data_para
!!$    USE timer
!!$
!!$    IMPLICIT NONE
!!$  
!!$#ifdef CPP_PARA
!!$    INCLUDE 'mpif.h'
!!$#endif
!!$    
!!$    LOGICAL,INTENT(IN) :: VarIn
!!$    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
!!$  
!!$#ifdef CPP_PARA
!!$    INTEGER :: nb,i,index_para,rank
!!$    INTEGER :: ierr
!!$    LOGICAL :: flag=.FALSE.
!!$    LOGICAL, PARAMETER :: check=.FALSE.
!!$#endif
!!$
!!$#ifndef CPP_PARA
!!$    VarOut(:)=VarIn
!!$    RETURN
!!$#else
!!$
!!$    IF (timer_state(timer_mpi)==running) THEN
!!$      flag=.TRUE.
!!$    ELSE
!!$      flag=.FALSE.
!!$    ENDIF
!!$    
!!$    IF (flag) CALL suspend_timer(timer_mpi)
!!$
!!$    IF (check) &
!!$         WRITE(numout,*) "gather_rgen VarIn=",VarIn    
!!$
!!$#ifdef CPP_PARA
!!$    CALL MPI_GATHER(VarIn,1,MPI_LOGICAL,VarOut,1,MPI_LOGICAL,root_prc,MPI_COMM_ORCH,ierr)
!!$#endif
!!$
!!$    IF (check) &
!!$         WRITE(numout,*) "gather_rgen VarOut=",VarOut
!!$    IF (flag) CALL resume_timer(timer_mpi)
!!$#endif
!!$  END SUBROUTINE gather_ls

!!!!! --> Les entiers

  SUBROUTINE gather_i(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

!    if (SIZE(VarIn) == 1) call stopit
    IF (is_root_prc) THEN
      CALL gather_igen(VarIn,VarOut,1)
    ELSE
      CALL gather_igen(VarIn,dummy,1)
    ENDIF
  
#endif
  END SUBROUTINE gather_i

!!!!!

  SUBROUTINE gather_i1(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

!    if (SIZE(VarIn) == 1) stop
    IF (is_root_prc) THEN
      CALL gather_igen(VarIn,VarOut,Size(VarIn,2))
    ELSE
      CALL gather_igen(VarIn,dummy,Size(VarIn,2))
    ENDIF
  
#endif
  END SUBROUTINE gather_i1

!!!!!
  
  SUBROUTINE gather_i2(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_igen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3))
    ELSE
      CALL gather_igen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3))
    ENDIF
  
#endif
  END SUBROUTINE gather_i2

!!!!!

  SUBROUTINE gather_i3(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_igen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ELSE
      CALL gather_igen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ENDIF
  
#endif
  END SUBROUTINE gather_i3

!!!!! --> Les reels

  SUBROUTINE gather_r(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

!    if (SIZE(VarIn) == 1) call stopit
    IF (is_root_prc) THEN
      CALL gather_rgen(VarIn,VarOut,1)
    ELSE
      CALL gather_rgen(VarIn,dummy,1)
    ENDIF
  
#endif
  END SUBROUTINE gather_r

!!!!!

  SUBROUTINE gather_r1(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_rgen(VarIn,VarOut,Size(VarIn,2))
    ELSE
      CALL gather_rgen(VarIn,dummy,Size(VarIn,2))
    ENDIF
  
#endif
  END SUBROUTINE gather_r1

!!!!!
  
  SUBROUTINE gather_r2(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_rgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3))
    ELSE
      CALL gather_rgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3))      
    ENDIF
  
#endif
  END SUBROUTINE gather_r2

!!!!!

  SUBROUTINE gather_r3(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_rgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ELSE
      CALL gather_rgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ENDIF
  
#endif
  END SUBROUTINE gather_r3

!!!!! --> Les booleen

  SUBROUTINE gather_l(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

!    if (SIZE(VarIn) == 1) call stopit
    IF (is_root_prc) THEN
      CALL gather_lgen(VarIn,VarOut,1)
    ELSE
      CALL gather_lgen(VarIn,dummy,1)      
    ENDIF
  
#endif
  END SUBROUTINE gather_l

!!!!!

  SUBROUTINE gather_l1(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_lgen(VarIn,VarOut,Size(VarIn,2))
    ELSE
      CALL gather_lgen(VarIn,dummy,Size(VarIn,2))
    ENDIF
  
#endif
  END SUBROUTINE gather_l1

!!!!!
  
  SUBROUTINE gather_l2(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_lgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3))
    ELSE
      CALL gather_lgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3))
    ENDIF
  
#endif
  END SUBROUTINE gather_l2

!!!!!

  SUBROUTINE gather_l3(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather_lgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ELSE
      CALL gather_lgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))      
    ENDIF
  
#endif
  END SUBROUTINE gather_l3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter2D   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter2D_i(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_igen(VarIn,VarOut,1)
    ELSE
      CALL scatter2D_igen(dummy,VarOut,1)
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_i

  SUBROUTINE scatter2D_i1(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_igen(VarIn,VarOut,SIZE(VarOut,3))
    ELSE
      CALL scatter2D_igen(dummy,VarOut,SIZE(VarOut,3))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_i1

  SUBROUTINE scatter2D_i2(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_igen(VarIn,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4))
    ELSE
      CALL scatter2D_igen(dummy,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_i2
  
  SUBROUTINE scatter2D_i3(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:,:)=VarIn(:,:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_igen(VarIn,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4)*SIZE(VarOut,5))
    ELSE
      CALL scatter2D_igen(dummy,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4)*SIZE(VarOut,5))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_i3


  SUBROUTINE scatter2D_r0(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_rgen(VarIn,VarOut,1)
    ELSE
      CALL scatter2D_rgen(dummy,VarOut,1)      
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_r0

  SUBROUTINE scatter2D_r(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_rgen(VarIn,VarOut,1)
    ELSE
      CALL scatter2D_rgen(dummy,VarOut,1)      
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_r

  SUBROUTINE scatter2D_r1(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_rgen(VarIn,VarOut,SIZE(VarOut,3))
    ELSE
      CALL scatter2D_rgen(dummy,VarOut,SIZE(VarOut,3))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_r1

  SUBROUTINE scatter2D_r2(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_rgen(VarIn,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4))
    ELSE
      CALL scatter2D_rgen(dummy,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_r2
  
  SUBROUTINE scatter2D_r3(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:,:)=VarIn(:,:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_rgen(VarIn,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4)*SIZE(VarOut,5))
    ELSE
      CALL scatter2D_rgen(dummy,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4)*SIZE(VarOut,5))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_r3  
  
  
  SUBROUTINE scatter2D_l(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_lgen(VarIn,VarOut,1)
    ELSE
      CALL scatter2D_lgen(dummy,VarOut,1)
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_l

  SUBROUTINE scatter2D_l1(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA    
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_lgen(VarIn,VarOut,SIZE(VarOut,3))
    ELSE
      CALL scatter2D_lgen(dummy,VarOut,SIZE(VarOut,3))
    ENDIF
  

#endif
  END SUBROUTINE scatter2D_l1

  SUBROUTINE scatter2D_l2(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_lgen(VarIn,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4))
    ELSE
      CALL scatter2D_lgen(dummy,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4))
    ENDIF
  
#endif
  END SUBROUTINE scatter2D_l2
  
  SUBROUTINE scatter2D_l3(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    LOGICAL :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:,:,:,:)=VarIn(:,:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL scatter2D_lgen(VarIn,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4)*SIZE(VarOut,5))
    ELSE
      CALL scatter2D_lgen(dummy,VarOut,SIZE(VarOut,3)*SIZE(VarOut,4)*SIZE(VarOut,5))
    ENDIF

#endif
  END SUBROUTINE scatter2D_l3 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Gather2D   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE gather2D_i(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_igen(VarIn,VarOut,1)
    ELSE
      CALL gather2D_igen(VarIn,dummy,1)
    ENDIF

#endif
  END SUBROUTINE gather2D_i

  SUBROUTINE gather2D_i1(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_igen(VarIn,VarOut,SIZE(VarIn,3))
    ELSE
      CALL gather2D_igen(VarIn,dummy,SIZE(VarIn,3))
    ENDIF

#endif
  END SUBROUTINE gather2D_i1

  SUBROUTINE gather2D_i2(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_igen(VarIn,VarOut,SIZE(VarIn,3)*SIZE(VarIn,4))
    ELSE
      CALL gather2D_igen(VarIn,dummy,SIZE(VarIn,3)*SIZE(VarIn,4))
    ENDIF

#endif
  END SUBROUTINE gather2D_i2
  
  SUBROUTINE gather2D_i3(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:,:)=VarIn(:,:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_igen(VarIn,VarOut,SIZE(VarIn,3)*SIZE(VarIn,4)*SIZE(VarIn,5))
    ELSE
      CALL gather2D_igen(VarIn,dummy,SIZE(VarIn,3)*SIZE(VarIn,4)*SIZE(VarIn,5))
    ENDIF

#endif
  END SUBROUTINE gather2D_i3


  SUBROUTINE gather2D_r0(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_rgen(VarIn,VarOut,1)
    ELSE
      CALL gather2D_rgen(VarIn,dummy,1)
    ENDIF

#endif
  END SUBROUTINE gather2D_r0

  SUBROUTINE gather2D_r(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_rgen(VarIn,VarOut,1)
    ELSE
      CALL gather2D_rgen(VarIn,dummy,1)
    ENDIF

#endif
  END SUBROUTINE gather2D_r

  SUBROUTINE gather2D_r1(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_rgen(VarIn,VarOut,SIZE(VarIn,3))
    ELSE
      CALL gather2D_rgen(VarIn,dummy,SIZE(VarIn,3))
    ENDIF

#endif
  END SUBROUTINE gather2D_r1

  SUBROUTINE gather2D_r2(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_rgen(VarIn,VarOut,SIZE(VarIn,3)*SIZE(VarIn,4))
    ELSE
      CALL gather2D_rgen(VarIn,dummy,SIZE(VarIn,3)*SIZE(VarIn,4))
    ENDIF

#endif
  END SUBROUTINE gather2D_r2
  
  SUBROUTINE gather2D_r3(VarIn, VarOut)

    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:,:)=VarIn(:,:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_rgen(VarIn,VarOut,SIZE(VarIn,3)*SIZE(VarIn,4)*SIZE(VarIn,5))
    ELSE
      CALL gather2D_rgen(VarIn,dummy,SIZE(VarIn,3)*SIZE(VarIn,4)*SIZE(VarIn,5))
    ENDIF
  

#endif
  END SUBROUTINE gather2D_r3  
  
  
  SUBROUTINE gather2D_l(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

#ifdef CPP_PARA    
    LOGICAL :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_lgen(VarIn,VarOut,1)
    ELSE
      CALL gather2D_lgen(VarIn,dummy,1)
    ENDIF
  

#endif
  END SUBROUTINE gather2D_l

  SUBROUTINE gather2D_l1(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA    
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_lgen(VarIn,VarOut,SIZE(VarIn,3))
    ELSE
      CALL gather2D_lgen(VarIn,dummy,SIZE(VarIn,3))
    ENDIF
  

#endif
  END SUBROUTINE gather2D_l1

  SUBROUTINE gather2D_l2(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

#ifdef CPP_PARA    
    LOGICAL :: dummy
#endif

#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_lgen(VarIn,VarOut,SIZE(VarIn,3)*SIZE(VarIn,4))
    ELSE
      CALL gather2D_lgen(VarIn,dummy,SIZE(VarIn,3)*SIZE(VarIn,4))
    ENDIF
  

#endif
  END SUBROUTINE gather2D_l2
  
  SUBROUTINE gather2D_l3(VarIn, VarOut)

    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
#ifdef CPP_PARA    
    LOGICAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:,:)=VarIn(:,:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL gather2D_lgen(VarIn,VarOut,SIZE(VarIn,3)*SIZE(VarIn,4)*SIZE(VarIn,5))
    ELSE
      CALL gather2D_lgen(VarIn,dummy,SIZE(VarIn,3)*SIZE(VarIn,4)*SIZE(VarIn,5))
    ENDIF
  

#endif
  END SUBROUTINE gather2D_l3  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des reduce_sum   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE reduce_sum_i(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  :: VarIn
    INTEGER,INTENT(OUT) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut=VarIn
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_igen(VarIn,Varout,1)
    ELSE
      CALL reduce_sum_igen(VarIn,dummy,1)
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_i

  SUBROUTINE reduce_sum_i1(VarIn, VarOut)

    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_i1

  SUBROUTINE reduce_sum_i2(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_i2

  SUBROUTINE reduce_sum_i3(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_i3

  SUBROUTINE reduce_sum_i4(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    INTEGER :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_i4                  
  
  
  SUBROUTINE reduce_sum_r(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut=VarIn
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_rgen(VarIn,Varout,1)
    ELSE
      CALL reduce_sum_rgen(VarIn,dummy,1)
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_r

  SUBROUTINE reduce_sum_r1(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_r1

  SUBROUTINE reduce_sum_r2(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_r2

  SUBROUTINE reduce_sum_r3(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_r3

  SUBROUTINE reduce_sum_r4(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
#ifdef CPP_PARA
    REAL :: dummy
#endif
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#else

    IF (is_root_prc) THEN
      CALL reduce_sum_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
#endif
  END SUBROUTINE reduce_sum_r4 
  
                            
END MODULE transfert_para    

#ifdef CPP_PARA

  SUBROUTINE bcast_cgen(var,nb)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    CHARACTER(LEN=*),DIMENSION(nb),INTENT(INOUT) :: Var
    INTEGER,INTENT(IN) :: nb
    
    INCLUDE 'mpif.h'

    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (check) &
         WRITE(numout,*) "bcast_cgen before bcast Var",Var
    IF (flag) CALL suspend_timer(timer_mpi)
    CALL MPI_BCAST(Var,nb*LEN(Var(1)),MPI_CHARACTER,root_prc,MPI_COMM_ORCH,ierr)
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
         WRITE(numout,*) "bcast_cgen after bcast Var",Var
        
  END SUBROUTINE bcast_cgen
      
  SUBROUTINE bcast_igen(var,nb)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    INTEGER,DIMENSION(nb),INTENT(INOUT) :: var
    INTEGER,INTENT(IN) :: nb
    
    INCLUDE 'mpif.h'

    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
   
    IF (check) &
         WRITE(numout,*) "bcast_igen before bcast Var",Var
    CALL MPI_BCAST(Var,nb,MPI_INT_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
         WRITE(numout,*) "bcast_igen after bcast Var",Var    
        
  END SUBROUTINE bcast_igen
  
  SUBROUTINE bcast_rgen(var,nb)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    REAL,DIMENSION(nb),INTENT(INOUT) :: var
    INTEGER,INTENT(IN) :: nb
    
    INCLUDE 'mpif.h'

    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (check) &
         WRITE(numout,*) "bcast_rgen before bcast Var",Var
    IF (flag) CALL suspend_timer(timer_mpi)    
    CALL MPI_BCAST(Var,nb,MPI_REAL_ORCH,root_prc,MPI_COMM_ORCH,ierr)
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
         WRITE(numout,*) "bcast_rgen after bcast Var",Var
    
  END SUBROUTINE bcast_rgen
  
  SUBROUTINE bcast_lgen(var,nb)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    LOGICAL,DIMENSION(nb),INTENT(INOUT) :: Var
    INTEGER,INTENT(IN) :: nb
    
    INCLUDE 'mpif.h'

    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.


    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (check) &
         WRITE(numout,*) "bcast_lgen before bcast Var",Var
    IF (flag) CALL suspend_timer(timer_mpi)    
    CALL MPI_BCAST(Var,nb,MPI_LOGICAL,root_prc,MPI_COMM_ORCH,ierr)
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
         WRITE(numout,*) "bcast_lgen after bcast Var",Var

  END SUBROUTINE bcast_lgen

  
  SUBROUTINE scatter_igen(VarIn, VarOut, dimsize)
    USE data_para
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(nbp_glo,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(nbp_loc,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    INTEGER,DIMENSION(dimsize*nbp_glo) :: VarTmp
    
    INTEGER :: nb,i,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
    
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp(Index_Para:Index_Para+nb-1)=VarIn(nbp_para_begin(rank):nbp_para_end(rank),i)
          Index_Para=Index_Para+nb
        ENDDO
      ENDDO
      IF (check) THEN
         WRITE(numout,*) "scatter_igen VarIn",VarIn
         WRITE(numout,*) "scatter_igen VarTmp",VarTmp
      ENDIF
    ENDIF
      
    CALL MPI_SCATTERV(VarTmp,counts,displs,MPI_INT_ORCH,VarOut,nbp_loc*dimsize,   &
                      MPI_INT_ORCH,root_prc, MPI_COMM_ORCH,ierr)
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "scatter_igen VarOut",VarOut

  END SUBROUTINE scatter_igen

  SUBROUTINE scatter_rgen(VarIn, VarOut, dimsize)
    USE data_para
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(nbp_glo,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(nbp_loc,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    REAL,DIMENSION(dimsize*nbp_glo) :: VarTmp
    
    INTEGER :: nb,i,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
    
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp(Index_Para:Index_Para+nb-1)=VarIn(nbp_para_begin(rank):nbp_para_end(rank),i)
          Index_Para=Index_Para+nb
        ENDDO
      ENDDO
      IF (check) THEN
         WRITE(numout,*) "scatter_rgen VarIn",VarIn
         WRITE(numout,*) "scatter_rgen VarTmp",VarTmp
      ENDIF
    ENDIF
      
    CALL MPI_SCATTERV(VarTmp,counts,displs,MPI_REAL_ORCH,VarOut,nbp_loc*dimsize,   &
                      MPI_REAL_ORCH,root_prc, MPI_COMM_ORCH,ierr)

    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "scatter_rgen VarOut",VarOut

  END SUBROUTINE scatter_rgen
  
  SUBROUTINE scatter_lgen(VarIn, VarOut, dimsize)
    USE data_para
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(nbp_glo,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(nbp_loc,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    LOGICAL,DIMENSION(dimsize*nbp_glo) :: VarTmp
    
    INTEGER :: nb,i,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
    
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp(Index_Para:Index_Para+nb-1)=VarIn(nbp_para_begin(rank):nbp_para_end(rank),i)
          Index_Para=Index_Para+nb
        ENDDO
      ENDDO
      IF (check) THEN
         WRITE(numout,*) "scatter_lgen VarIn",VarIn
         WRITE(numout,*) "scatter_lgen VarTmp",VarTmp
      ENDIF
    ENDIF
      
    CALL MPI_SCATTERV(VarTmp,counts,displs,MPI_LOGICAL,VarOut,nbp_loc*dimsize,   &
                      MPI_LOGICAL,root_prc, MPI_COMM_ORCH,ierr)
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "scatter_lgen VarOut",VarOut

  END SUBROUTINE scatter_lgen  

  SUBROUTINE gather_igen(VarIn, VarOut, dimsize)
    USE data_para
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(nbp_loc,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(nbp_glo,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'
    
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    INTEGER,DIMENSION(dimsize*nbp_glo) :: VarTmp
    
    INTEGER :: nb,i,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)

    IF (is_root_prc) THEN
      Index_Para=1
      IF (check) &
           WRITE(numout,*) "gather_igen mpi_size, dimsize, nbp_glo",mpi_size, dimsize, nbp_glo
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
	Index_Para=Index_Para+nb*dimsize
      ENDDO
       IF (check) &
            WRITE(numout,*) "gather_igen nbp_para_nb, displs, counts,Index_Para-1",nbp_para_nb, displs, counts,Index_Para-1
     
    ENDIF
    
    IF (check) &
         WRITE(numout,*) "gather_igen VarIn=",VarIn    
    CALL MPI_GATHERV(VarIn,nbp_loc*dimsize,MPI_INT_ORCH,VarTmp,counts,displs,   &
         MPI_INT_ORCH,root_prc, MPI_COMM_ORCH,ierr)

    IF (check) &
         WRITE(numout,*) "gather_igen dimsize,VarTmp=",dimsize,VarTmp
		          
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        DO i=1,dimsize
          VarOut(nbp_para_begin(rank):nbp_para_end(rank),i)=VarTmp(Index_Para:Index_Para+nb-1)
	  Index_Para=Index_Para+nb
        ENDDO
      ENDDO
    ENDIF
    IF (check) &
         WRITE(numout,*) "gather_igen VarOut=",VarOut
    IF (flag) CALL resume_timer(timer_mpi)

  END SUBROUTINE gather_igen  

  SUBROUTINE gather_rgen(VarIn, VarOut, dimsize)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(nbp_loc,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(nbp_glo,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'
  
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    REAL,DIMENSION(dimsize*nbp_glo) :: VarTmp
    
    INTEGER :: nb,i,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)

    IF (is_root_prc) THEN
      Index_Para=1
      IF (check) &
           WRITE(numout,*) "gather_rgen mpi_size, dimsize, nbp_glo",mpi_size, dimsize, nbp_glo
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
	Index_Para=Index_Para+nb*dimsize
      ENDDO
      IF (check) &
           WRITE(numout,*) "gather_rgen nbp_para_nb, displs, counts,Index_Para-1",nbp_para_nb, displs, counts,Index_Para-1
     
    ENDIF
    
    IF (check) &
         WRITE(numout,*) "gather_rgen VarIn=",VarIn    
    CALL MPI_GATHERV(VarIn,nbp_loc*dimsize,MPI_REAL_ORCH,VarTmp,counts,displs,   &
                      MPI_REAL_ORCH,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
         WRITE(numout,*) "gather_rgen dimsize,VarTmp=",dimsize,VarTmp
		          
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        DO i=1,dimsize
          VarOut(nbp_para_begin(rank):nbp_para_end(rank),i)=VarTmp(Index_Para:Index_Para+nb-1)
	  Index_Para=Index_Para+nb
        ENDDO
      ENDDO
    ENDIF
    IF (check) &
         WRITE(numout,*) "gather_rgen VarOut=",VarOut
    IF (flag) CALL resume_timer(timer_mpi)

  END SUBROUTINE gather_rgen  

  SUBROUTINE gather_lgen(VarIn, VarOut, dimsize)
    USE data_para
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(nbp_loc,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(nbp_glo,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    LOGICAL,DIMENSION(dimsize*nbp_glo) :: VarTmp
    
    INTEGER :: nb,i,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.


    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)

    IF (is_root_prc) THEN
      Index_Para=1
      IF (check) &
           WRITE(numout,*) "gather_lgen mpi_size, dimsize, nbp_glo",mpi_size, dimsize, nbp_glo
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
	Index_Para=Index_Para+nb*dimsize
      ENDDO
      IF (check) &
           WRITE(numout,*) "gather_lgen nbp_para_nb, displs, counts,Index_Para-1",nbp_para_nb, displs, counts,Index_Para-1
    ENDIF
    
    IF (check) &
         WRITE(numout,*) "gather_lgen VarIn=",VarIn    
    CALL MPI_GATHERV(VarIn,nbp_loc*dimsize,MPI_LOGICAL,VarTmp,counts,displs,   &
                      MPI_LOGICAL,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
         WRITE(numout,*) "gather_lgen dimsize,VarTmp=",dimsize,VarTmp
		          
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=nbp_para_nb(rank)
        DO i=1,dimsize
          VarOut(nbp_para_begin(rank):nbp_para_end(rank),i)=VarTmp(Index_Para:Index_Para+nb-1)
	  Index_Para=Index_Para+nb
        ENDDO
      ENDDO
    ENDIF
    IF (check) &
         WRITE(numout,*) "gather_lgen VarOut=",VarOut
    IF (flag) CALL resume_timer(timer_mpi)

  END SUBROUTINE gather_lgen
  

  SUBROUTINE scatter2D_igen(VarIn, VarOut, dimsize)
    USE data_para, iim=>iim_g,jjm=>jjm_g
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(iim*jjm,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(iim*jj_nb,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1)   :: displs
    INTEGER,DIMENSION(0:mpi_size-1)   :: counts
    INTEGER,DIMENSION(dimsize*iim*jjm)   :: VarTmp1
    INTEGER,DIMENSION(ij_nb,dimsize)     :: VarTmp2
    
    INTEGER :: nb,i,ij,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
    
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp1(Index_Para:Index_Para+nb-1)=VarIn(ij_para_begin(rank):ij_para_end(rank),i)
          Index_Para=Index_Para+nb
        ENDDO
      ENDDO
      IF (check) THEN
         WRITE(numout,*) "scatter2D_igen VarIn",VarIn
         WRITE(numout,*) "scatter2D_igen VarTmp1",VarTmp1
      ENDIF
    ENDIF
      
    CALL MPI_SCATTERV(VarTmp1,counts,displs,MPI_INT_ORCH,VarTmp2,ij_nb*dimsize,   &
                      MPI_INT_ORCH,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
         WRITE(numout,*) "scatter2D_igen VarTmp2",VarTmp2
   
    DO i=1,dimsize
      DO ij=1,ij_nb
        VarOut(ij+ii_begin-1,i)=VarTmp2(ij,i)
      ENDDO
    ENDDO
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "scatter2D_igen VarOut",VarOut

  END SUBROUTINE scatter2D_igen
  
  
  SUBROUTINE scatter2D_rgen(VarIn, VarOut, dimsize)
    USE data_para, iim=>iim_g,jjm=>jjm_g
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(iim*jjm,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(iim*jj_nb,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1)   :: displs
    INTEGER,DIMENSION(0:mpi_size-1)   :: counts
    REAL,DIMENSION(dimsize*iim*jjm)   :: VarTmp1
    REAL,DIMENSION(ij_nb,dimsize)     :: VarTmp2
    REAL,DIMENSION(iim*jj_nb,dimsize) :: VarOut_bis
    
    INTEGER :: nb,i,ij,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
    
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp1(Index_Para:Index_Para+nb-1)=VarIn(ij_para_begin(rank):ij_para_end(rank),i)
          Index_Para=Index_Para+nb
        ENDDO
      ENDDO
      IF (check) THEN
         WRITE(numout,*) "scatter2D_rgen VarIn",VarIn
         WRITE(numout,*) "scatter2D_rgen VarTmp1",VarTmp1
      ENDIF
    ENDIF
    nb=ij_nb*dimsize
    IF (check) &
         WRITE(numout,*) "ij_nb*dimsize",ij_nb*dimsize
      
    CALL MPI_SCATTERV(VarTmp1,counts,displs,MPI_REAL_ORCH,VarTmp2,nb,   &
                      MPI_REAL_ORCH,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
         WRITE(numout,*) "scatter2D_rgen VarTmp2",VarTmp2

    DO i=1,dimsize
      DO ij=1,ij_nb
        VarOut(ij+ii_begin-1,i)=VarTmp2(ij,i)
      ENDDO
    ENDDO

    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "scatter2D_rgen VarOut",VarOut

  END SUBROUTINE scatter2D_rgen

  SUBROUTINE scatter2D_lgen(VarIn, VarOut, dimsize)
    USE data_para, iim=>iim_g,jjm=>jjm_g
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(iim*jjm,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(iim*jj_nb,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1)   :: displs
    INTEGER,DIMENSION(0:mpi_size-1)   :: counts
    LOGICAL,DIMENSION(dimsize*iim*jjm)   :: VarTmp1
    LOGICAL,DIMENSION(ij_nb,dimsize)     :: VarTmp2
    
    INTEGER :: nb,i,ij,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)
    
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp1(Index_Para:Index_Para+nb-1)=VarIn(ij_para_begin(rank):ij_para_end(rank),i)
          Index_Para=Index_Para+nb
        ENDDO
      ENDDO
      IF (check) THEN
         WRITE(numout,*) "scatter2D_lgen VarIn",VarIn
         WRITE(numout,*) "scatter2D_lgen VarTmp1",VarTmp1
      ENDIF
    ENDIF
      
    CALL MPI_SCATTERV(VarTmp1,counts,displs,MPI_LOGICAL,VarTmp2,ij_nb*dimsize,   &
                      MPI_LOGICAL,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
         WRITE(numout,*) "scatter2D_lgen VarTmp2",VarTmp2
   
    DO i=1,dimsize
      DO ij=1,ij_nb
        VarOut(ij+ii_begin-1,i)=VarTmp2(ij,i)
      ENDDO
    ENDDO
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "scatter2D_lgen VarOut",VarOut

  END SUBROUTINE scatter2D_lgen


  SUBROUTINE gather2D_igen(VarIn, VarOut, dimsize)
    USE data_para, iim=>iim_g,jjm=>jjm_g
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(iim*jj_nb,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(iim*jjm,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    INTEGER,DIMENSION(ij_nb,dimsize)   :: VarTmp1
    INTEGER,DIMENSION(dimsize*iim*jjm) :: VarTmp2
    
    INTEGER :: nb,i,ij,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL,PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)

    IF (is_root_prc) THEN
      Index_Para=1
      IF (check) &
           WRITE(numout,*) "gather2D_igen mpi_size, dimsize, nbp_glo",mpi_size, dimsize, nbp_glo
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
	Index_Para=Index_Para+nb*dimsize
      ENDDO
      IF (check) &
           WRITE(numout,*) "gather2D_igen nbp_para_nb, displs, counts,Index_Para-1",nbp_para_nb, displs, counts,Index_Para-1
    ENDIF
    DO i=1,dimsize
       DO ij=1,ij_nb
          VarTmp1(ij,i)=VarIn(ij+ii_begin-1,i)
       ENDDO
    ENDDO
    
    IF (check) THEN
       WRITE(numout,*) "gather2D_igen VarIn=",VarIn    
       WRITE(numout,*) "gather2D_igen VarTmp1=",VarTmp1
    ENDIF
    CALL MPI_GATHERV(VarTmp1,ij_nb*dimsize,MPI_INT_ORCH,VarTmp2,counts,displs,   &
                     MPI_INT_ORCH,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
       WRITE(numout,*) "gather2D_igen VarTmp2=",VarTmp2
		          
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        DO i=1,dimsize
          VarOut(ij_para_begin(rank):ij_para_end(rank),i)=VarTmp2(Index_Para:Index_Para+nb-1)
	  Index_Para=Index_Para+nb
        ENDDO
      ENDDO
    ENDIF
    
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "gather2D_igen VarOut=",VarOut

  END SUBROUTINE gather2D_igen   


  
  SUBROUTINE gather2D_rgen(VarIn, VarOut, dimsize)
    USE data_para, iim=>iim_g,jjm=>jjm_g
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(iim*jj_nb,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(iim*jjm,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'
  
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    REAL,DIMENSION(ij_nb,dimsize)   :: VarTmp1
    REAL,DIMENSION(dimsize*iim*jjm) :: VarTmp2
    
    INTEGER :: nb,i,ij,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL,PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)

    IF (is_root_prc) THEN
      Index_Para=1
      IF (check) &
           WRITE(numout,*) "gather2D_rgen mpi_size, dimsize, nbp_glo",mpi_size, dimsize, nbp_glo
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
	Index_Para=Index_Para+nb*dimsize
      ENDDO
      IF (check) &
           WRITE(numout,*) "gather2D_rgen nbp_para_nb, displs, counts,Index_Para-1",nbp_para_nb, displs, counts,Index_Para-1
    ENDIF
    
    DO i=1,dimsize
      DO ij=1,ij_nb
        VarTmp1(ij,i)=VarIn(ij+ii_begin-1,i)
      ENDDO
    ENDDO

    IF (check) THEN
       WRITE(numout,*) "gather2D_rgen VarIn=",VarIn    
       WRITE(numout,*) "gather2D_rgen VarTmp1=",VarTmp1
    ENDIF
    CALL MPI_GATHERV(VarTmp1,ij_nb*dimsize,MPI_REAL_ORCH,VarTmp2,counts,displs,   &
                     MPI_REAL_ORCH,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
       WRITE(numout,*) "gather2D_rgen VarTmp2=",VarTmp2

    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        DO i=1,dimsize
          VarOut(ij_para_begin(rank):ij_para_end(rank),i)=VarTmp2(Index_Para:Index_Para+nb-1)
	  Index_Para=Index_Para+nb
        ENDDO
      ENDDO
    ENDIF
    
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "gather2D_rgen VarOut=",VarOut

  END SUBROUTINE gather2D_rgen   

  SUBROUTINE gather2D_lgen(VarIn, VarOut, dimsize)
    USE data_para, iim=>iim_g,jjm=>jjm_g
    USE timer

    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(iim*jj_nb,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(iim*jjm,dimsize) :: VarOut
  
    INCLUDE 'mpif.h'
  
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    LOGICAL,DIMENSION(ij_nb,dimsize)   :: VarTmp1
    LOGICAL,DIMENSION(dimsize*iim*jjm) :: VarTmp2
    
    INTEGER :: nb,i,ij,index_para,rank
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL,PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (flag) CALL suspend_timer(timer_mpi)

    IF (is_root_prc) THEN
      Index_Para=1
      IF (check) &
           WRITE(numout,*) "gather2D_lgen mpi_size, dimsize, nbp_glo",mpi_size, dimsize, nbp_glo
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        displs(rank)=Index_Para-1
        counts(rank)=nb*dimsize
	Index_Para=Index_Para+nb*dimsize
      ENDDO
      IF (check) &
           WRITE(numout,*) "gather2D_lgen nbp_para_nb, displs, counts,Index_Para-1",nbp_para_nb, displs, counts,Index_Para-1
    ENDIF
    
    DO i=1,dimsize
      DO ij=1,ij_nb
        VarTmp1(ij,i)=VarIn(ij+ii_begin-1,i)
      ENDDO
    ENDDO
    
    IF (check) THEN
       WRITE(numout,*) "gather2D_lgen VarIn=",VarIn    
       WRITE(numout,*) "gather2D_lgen VarTmp1=",VarTmp1
    ENDIF
    CALL MPI_GATHERV(VarTmp1,ij_nb*dimsize,MPI_LOGICAL,VarTmp2,counts,displs,   &
                     MPI_LOGICAL,root_prc, MPI_COMM_ORCH,ierr)
    IF (check) &
       WRITE(numout,*) "gather2D_lgen VarTmp2=",VarTmp2
		          
    IF (is_root_prc) THEN
      Index_Para=1
      DO rank=0,mpi_size-1
        nb=ij_para_nb(rank)
        DO i=1,dimsize
          VarOut(ij_para_begin(rank):ij_para_end(rank),i)=VarTmp2(Index_Para:Index_Para+nb-1)
	  Index_Para=Index_Para+nb
        ENDDO
      ENDDO
    ENDIF
    
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "gather2D_lgen VarOut=",VarOut

  END SUBROUTINE gather2D_lgen   

  SUBROUTINE reduce_sum_igen(VarIn,VarOut,nb)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    INTEGER,DIMENSION(nb),INTENT(IN) :: VarIn
    INTEGER,DIMENSION(nb),INTENT(OUT) :: VarOut    
    INTEGER,INTENT(IN) :: nb
    
    INCLUDE 'mpif.h'

    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (check) &
       WRITE(numout,*) "reduce_sum_igen VarIn",VarIn
    IF (flag) CALL suspend_timer(timer_mpi)
   
    CALL MPI_REDUCE(VarIn,VarOut,nb,MPI_INT_ORCH,MPI_SUM,root_prc,MPI_COMM_ORCH,ierr)
            
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "reduce_sum_igen VarOut",VarOut

  END SUBROUTINE reduce_sum_igen
  
  SUBROUTINE reduce_sum_rgen(VarIn,VarOut,nb)
    USE data_para
    USE timer

    IMPLICIT NONE
    
    REAL,DIMENSION(nb),INTENT(IN) :: VarIn
    REAL,DIMENSION(nb),INTENT(OUT) :: VarOut    
    INTEGER,INTENT(IN) :: nb

    INCLUDE 'mpif.h'
    
    INTEGER :: ierr
    LOGICAL :: flag=.FALSE.
    LOGICAL, PARAMETER :: check=.FALSE.

    IF (timer_state(timer_mpi)==running) THEN
      flag=.TRUE.
    ELSE
      flag=.FALSE.
    ENDIF
    
    IF (check) &
       WRITE(numout,*) "reduce_sum_rgen VarIn",VarIn
    IF (flag) CALL suspend_timer(timer_mpi)
   
    CALL MPI_REDUCE(VarIn,VarOut,nb,MPI_REAL_ORCH,MPI_SUM,root_prc,MPI_COMM_ORCH,ierr)
        
    IF (flag) CALL resume_timer(timer_mpi)
    IF (check) &
       WRITE(numout,*) "reduce_sum_rgen VarOut",VarOut

  END SUBROUTINE reduce_sum_rgen

  subroutine stopit
    USE ioipsl
    call MPI_FINALIZE

    CALL ipslerr (3,'transfert_para : gather', &
         &          'A gather function was called with a VarIn',&
         &          ' which size is only one.', &
         &          '(must be strickly greater than one )')
  end subroutine stopit
#endif
