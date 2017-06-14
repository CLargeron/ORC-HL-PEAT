! Yann Meurdesoif functions for parallel tests.

!-
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parallel/orch_write_field_p.f90 $ 
!< $Date: 2012-02-20 14:58:25 +0100 (Mon, 20 Feb 2012) $
!< $Author: didier.solyga $
!< $Revision: 720 $
!-

MODULE Orch_Write_field_p
  
  interface WriteField_p
    module procedure WriteField_4d_p,WriteField_3d_p,WriteField_2d_p
  end interface WriteField_p
 
  interface WriteFieldI_p
    module procedure WriteFieldI_3d_p,WriteFieldI_2d_p,WriteFieldI_1d_p
  end interface WriteFieldI_p  
  
  
CONTAINS

  SUBROUTINE init_WriteField_p(index_Write_Field_p)
  USE parallel
  USE Write_Field, only : Init_WriteField
  IMPLICIT NONE
    INTEGER,INTENT(in) :: index_Write_Field_p(nbp_loc)

    IF (is_root_prc) CALL Init_WriteField(iim_g,jjm_g,nbp_glo,index_g)
    
  END SUBROUTINE init_WriteField_p

  SUBROUTINE WriteField_4d_p(name,Field)
    USE parallel
    USE Write_field, only : WriteField
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:,:) :: Field 
      INTEGER, DIMENSION(4) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: Field_g
      
      
      Dim=shape(Field)
      
      ALLOCATE(Field_g(iim_g,jjm_g,Dim(3),Dim(4)))
      CALL Gather2D(Field,Field_g)

      IF (is_root_prc) CALL WriteField(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteField_4d_p
    
  SUBROUTINE WriteField_3d_p(name,Field)
    USE parallel
    USE Write_field, only : WriteField
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:) :: Field 
      INTEGER, DIMENSION(3) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Field_g
      
      
      Dim=shape(Field)
      
      ALLOCATE(Field_g(iim_g,jjm_g,Dim(3)))
      CALL Gather2D(Field,Field_g)

      IF (is_root_prc) CALL WriteField(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteField_3d_p

  SUBROUTINE WriteField_2d_p(name,Field)
    USE parallel
    USE Write_field, only : WriteField
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:) :: Field 
      INTEGER, DIMENSION(2) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:) :: Field_g
      
      
      Dim=shape(Field)
      
      ALLOCATE(Field_g(iim_g,jjm_g))
      CALL Gather2D(Field,Field_g)

      IF (is_root_prc) CALL WriteField_gen(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteField_2d_p

  SUBROUTINE WriteFieldI_3d_p(name,Field)
    USE parallel
    USE Write_field, only : WriteFieldI
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:,:) :: Field 
      INTEGER, DIMENSION(3) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Field_g
      
      
      Dim=shape(Field)
      
      ALLOCATE(Field_g(nbp_glo,Dim(2),Dim(3)))
      CALL gather(Field,Field_g)
      
      IF (is_root_prc) CALL WriteFieldI(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteFieldI_3d_p

  SUBROUTINE WriteFieldI_2d_p(name,Field)
    USE parallel
    USE Write_field, only : WriteFieldI
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:,:) :: Field 
      INTEGER, DIMENSION(2) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:,:) :: Field_g
      
      
      Dim=shape(Field)
      
      ALLOCATE(Field_g(nbp_glo,Dim(2)))
      CALL gather(Field,Field_g)
      
      IF (is_root_prc) CALL WriteFieldI(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteFieldI_2d_p    

  SUBROUTINE WriteFieldI_1d_p(name,Field)
    USE parallel
    USE Write_field, only : WriteFieldI
    IMPLICIT NONE
      CHARACTER(len=*) :: name
      REAL, DIMENSION(:) :: Field 
      INTEGER, DIMENSION(1) :: Dim
      
      REAL, ALLOCATABLE, DIMENSION(:) :: Field_g
      
      
      Dim=shape(Field)
      
      ALLOCATE(Field_g(nbp_glo))
      CALL gather(Field,Field_g)
      
      IF (is_root_prc) CALL WriteFieldI(name,Field_g)  
      
      DEALLOCATE(Field_g)
  END SUBROUTINE WriteFieldI_1d_p    
    
END MODULE Orch_Write_field_p
