! Parallel tools : Barrier and Finalize.

!-
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_parallel/tools_para.f90 $ 
!< $Date: 2012-02-20 14:58:25 +0100 (Mon, 20 Feb 2012) $
!< $Author: didier.solyga $
!< $Revision: 720 $
!-

MODULE tools_para
!-
  USE ioipsl
  USE defprec
  USE timer
  USE data_para
!-
#include "src_parallel.h"
!-
CONTAINS
  SUBROUTINE stop_mpi()
#ifdef CPP_PARA
    CALL MPI_COMM_FREE(MPI_COMM_ORCH,ierr)
    CALL MPI_FINALIZE(ierr)
#endif
    CALL ipslerr(3,'STOP_MPI','MPI has been stopped in ORCHIDEE.',&
 &                 "Don't know the reason","Please verify output.")
  ENd subroutine stop_mpi

  SUBROUTINE barrier_para()
#ifdef CPP_PARA
    CALL MPI_BARRIER(MPI_COMM_ORCH,ierr)
#endif
  END SUBROUTINE barrier_para

  SUBROUTINE finalize_para(timer_global,timer_mpi)
    INTEGER, INTENT(IN) :: timer_global, timer_mpi
    DOUBLE PRECISION :: cpu_time_mpi

    cpu_time_mpi = Get_cpu_time(timer_mpi)
    WRITE(numout,*) '*********************************************************'
    WRITE(numout,*) ' TEMPS GLOBAL   ---> REAL TIME :',Get_real_time(timer_global)
    WRITE(numout,*) ' TEMPS GLOBAL   ---> CPU TIME  :',Get_cpu_time(timer_global)
    WRITE(numout,*) ' TEMPS HORS MPI ---> REAL TIME :',Get_real_time(timer_mpi)
    WRITE(numout,*) ' TEMPS HORS MPI ---> CPU TIME  :',cpu_time_mpi
    WRITE(numout,*) '*********************************************************' 
    WRITE(numout,*) 'END OF RUN.'

    CALL Write_Load_Balance(REAL(cpu_time_mpi, r_std))
#ifdef CPP_PARA
!    CALL MPI_COMM_FREE(MPI_COMM_ORCH,ierr)
    CALL MPI_FINALIZE(ierr)
#endif
  END SUBROUTINE finalize_para

END MODULE tools_para
