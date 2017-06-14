! Overlap of IOIPSL functions for specific parallel use in ORCHIDEE.

!-
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parallel/ioipsl_para.f90 $ 
!< $Date: 2012-02-20 14:58:25 +0100 (Mon, 20 Feb 2012) $
!< $Author: didier.solyga $
!< $Revision: 720 $
!-

MODULE ioipsl_para
  USE ioipsl
  USE data_para
  USE transfert_para
!-
  IMPLICIT NONE
!-
#include "src_parallel.h"
!-
  INTERFACE getin_p
    MODULE PROCEDURE getin_p_c,getin_p_c1,   &
         getin_p_i,getin_p_i1,getin_p_i2,&
         getin_p_r,getin_p_r1,getin_p_r2,&
         getin_p_l,getin_p_l1,getin_p_l2
  END INTERFACE
!-
  INTERFACE restput_p
     MODULE PROCEDURE &
          restput_p_r3d, restput_p_r2d, restput_p_r1d, &
          restput_p_opp_r2d, restput_p_opp_r1d
  END INTERFACE
!-
  INTERFACE restget_p
     MODULE PROCEDURE &
          restget_p_r3d, restget_p_r2d, restget_p_r1d, &
          restget_p_opp_r2d, restget_p_opp_r1d
  END INTERFACE

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Definition des getin -> bcast      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Les chaines de caracteres -- !!
  
  SUBROUTINE getin_p_c(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    CHARACTER(LEN=*),INTENT(INOUT) :: VarOut    

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_c  


  SUBROUTINE getin_p_c1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    CHARACTER(LEN=*),INTENT(INOUT) :: VarOut(:)    

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_c1 

!! -- Les entiers -- !!
  
  SUBROUTINE getin_p_i(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut    

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_i

  SUBROUTINE getin_p_i1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut(:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_i1

  SUBROUTINE getin_p_i2(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    INTEGER,INTENT(INOUT) :: VarOut(:,:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_i2

!! -- Les flottants -- !!
  
  SUBROUTINE getin_p_r(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_r

  SUBROUTINE getin_p_r1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut(:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_r1

  SUBROUTINE getin_p_r2(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    REAL,INTENT(INOUT) :: VarOut(:,:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_r2

!! -- Les Booleens -- !!
  
  SUBROUTINE getin_p_l(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_l

  SUBROUTINE getin_p_l1(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut(:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_l1

  SUBROUTINE getin_p_l2(VarIn,VarOut)
    IMPLICIT NONE    
    CHARACTER(LEN=*),INTENT(IN) :: VarIn
    LOGICAL,INTENT(INOUT) :: VarOut(:,:)

    IF (is_root_prc) CALL getin(VarIn,VarOut)
    CALL bcast(VarOut)
  END SUBROUTINE getin_p_l2
!-
!-----------------------------
!-----------------------------
!-----------------------------
!-
  SUBROUTINE restget_p_opp_r1d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
! DO NOT USE THIS FUNCTION WITH NON GRID VARIABLE !
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim*jjm*llm) )
       CALL restget &
            (fid, vname_q, iim, jjm, llm, itau, def_beha, &
            temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    CALL scatter(temp_g,var)
    IF (is_root_prc) DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_r1d
!-
!===
!-
  SUBROUTINE restget_p_opp_r2d &
  (fid, vname_q, iim, jjm, llm, itau, def_beha, &
   var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
    !-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm) )
       CALL restget &
            (fid, vname_q, iim, jjm, llm, itau, def_beha, &
            temp_g, MY_OPERATOR, nbindex, ijndex)
    ENDIF
    CALL scatter(temp_g,var)
    IF (is_root_prc) DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_opp_r2d
!-
!===
!-
  SUBROUTINE restget_p_r1d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
! DO NOT USE THIS FUNCTION WITH NON GRID VARIABLE !
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL :: def_beha
    REAL :: var(:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim*jjm*llm) )
       CALL restget &
            (fid,vname_q,iim,jjm,llm,itau,def_beha,temp_g)
    ENDIF
    CALL scatter(temp_g,var)
    IF (is_root_prc) DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_r1d
!-
!===
!-
  SUBROUTINE restget_p_r2d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL :: def_beha
    REAL :: var(:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm) )
       CALL restget &
            (fid,vname_q,iim,jjm,llm,itau,def_beha,temp_g)
    ENDIF
    CALL scatter(temp_g,var)
    IF (is_root_prc) DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_r2d
!-
!===
!-
  SUBROUTINE restget_p_r3d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    LOGICAL def_beha
    REAL :: var(:,:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g

    IF (is_root_prc) THEN 
       ALLOCATE( temp_g(iim,jjm,llm) )
       CALL restget &
            (fid,vname_q,iim,jjm,llm,itau,def_beha,temp_g)
    ENDIF
    CALL scatter(temp_g,var)
    IF (is_root_prc) DEALLOCATE(temp_g)
  END SUBROUTINE restget_p_r3d
!-
!-----------------------------
!-----------------------------
!-
  SUBROUTINE restput_p_opp_r1d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) ALLOCATE( temp_g(iim*jjm*llm) )
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)

       DEALLOCATE( temp_g )
    ENDIF
          
  END SUBROUTINE restput_p_opp_r1d
!-
!===
!-
  SUBROUTINE restput_p_opp_r2d &
  (fid, vname_q, iim, jjm, llm, itau, var, MY_OPERATOR, nbindex, ijndex)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:)
    CHARACTER(LEN=*) :: MY_OPERATOR
    INTEGER :: nbindex, ijndex(nbindex)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) ALLOCATE( temp_g(iim,jjm) )
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput &
            (fid, vname_q, iim, jjm, llm, itau, temp_g, MY_OPERATOR, nbindex, ijndex)
       DEALLOCATE( temp_g )
    ENDIF
          
  END SUBROUTINE restput_p_opp_r2d
!-
!===
!-
  SUBROUTINE restput_p_r1d (fid,vname_q,iim,jjm,llm,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:)
    !-----------------------------
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_g

    IF (is_root_prc) ALLOCATE( temp_g(iim*jjm*llm) )
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput (fid,vname_q,iim,jjm,llm,itau,temp_g)
       DEALLOCATE( temp_g )
    ENDIF
          
  END SUBROUTINE restput_p_r1d
!-
!===
!-
  SUBROUTINE restput_p_r2d (fid,vname_q,iim,jjm,llm,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_g

    IF (is_root_prc) ALLOCATE( temp_g(iim,jjm) )
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput (fid,vname_q,iim,jjm,llm,itau,temp_g)
       DEALLOCATE( temp_g )
    ENDIF
          
  END SUBROUTINE restput_p_r2d
!-
!===
!-
  SUBROUTINE restput_p_r3d (fid,vname_q,iim,jjm,llm,itau,var)
    IMPLICIT NONE
!-
    INTEGER :: fid
    CHARACTER(LEN=*) :: vname_q
    INTEGER :: iim, jjm, llm, itau
    REAL :: var(:,:,:)
    !-------------------------
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_g

    IF (is_root_prc) ALLOCATE( temp_g(iim,jjm,llm) )
    CALL gather(var,temp_g)
    IF (is_root_prc) THEN
       CALL restput (fid,vname_q,iim,jjm,llm,itau,temp_g)
       DEALLOCATE( temp_g )
    ENDIF
          
  END SUBROUTINE restput_p_r3d

END MODULE ioipsl_para
