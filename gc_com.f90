!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! GC - General Communication primitives package. For use on multiprocessor
! shared memory and message passing systems.
!
!
! LICENSING TERMS
!
!  GC is provided free of charge. Unless otherwise agreed with SINTEF, use and
!  redistribution in source and binary forms are permitted provided that
!
!      (1) source distributions retain all comments appearing within this file
!          header, and
!
!      (2) distributions including binaries display the following
!          acknowledgement:
!
!              "This product includes software developed by SINTEF.",
!
!          in the documentation or other materials provided with the
!          distribution and in all advertising materials mentioning features or
!          use of this software.
!
!  The name of SINTEF may not be used to endorse or promote products derived
!  from this software without specific prior written permission.  SINTEF
!  disclaims any warranty that this software will be fit for any specific
!  purposes. In no event shall SINTEF be liable for any loss of performance or
!  for indirect or consequential damage or direct or indirect injury of any
!  kind. In no case shall SINTEF be liable for any representation or warranty
!  make to any third party by the users of this software.
!
! USAGE
!
!  Define zero or ONE of the INTERFACE flags below with -D<name> and compile
!  using the C (or Fortran) preprocessor.
!
! INTERFACE FLAGS
!
!  SHM_SRC  -  CRAY MPP systems CRI using shared memory (SHMEM)
!  NX2_SRC  -  Intel Paragon using Intel/NX2 message passing
!  PVM_SRC  -  Public ORNL PVM/Cray PVP  message passing
!  1  -  The Message Passing interface
!  MPL_SRC  -  The IBM Message Passing Language
!
!
! INTERFACE CONFIGURATION FLAGS
!
!  Define zero or more of the flags below with -D<name> and compile using the
!  C (or Fortran) preprocessor.
!
!  FLP_32B  -  use 32 bit floating point precision. Enabled by default for
!              all systems except CRI PVM & MPP.
!  1  -  use 64 bit floating point precision. Enabled by default on CRI
!              systems.
!  PVM_T3D  -  enable CRAY MPP PVM specifics (for backwards compatibility only,
!              this functionality is always enabled on CRI MPP systems using
!              PVM_SRC).
!  PVM_V33  -  has PVM v3.3 or higher, with psend/precv.
!  PVM_NBS  -  PVM No Byte Swap - disable byte sex swapping. Byte swapping
!              of INTEGER and REAl data is on by default.
!  PVM_DIP  -  PVM Data In Place, enabled with ORNL PVM v3.3 and higher.
!              Do not copy message data from user space into PVM buffers
!              on send. Data must NOT be modified before leaving the sending
!              node (before arriving receiving node with CRI MPP PVM).
!
!  NX2_32B, NX2_64B, PVM_32B, PVM_64B, MPI_32B, MPI_64B are supported for
!  backward compatibility with prerelease versions only. Use the FLP_nnB flags.
!
!
! REVISION HISTORY
!
!  1.0 - Initial version based on the GCOM interface by R.Skaalin and
!        the PAR interface by J.Amundsen. The MPL part is due to Z. Christidis,
!        IBM T. J. Watson Research Center.
!
!-------------------------------------------------------------------------------
! $Id: gc_com.f90,v 1.1 2006-10-27 14:54:41 mifapw Exp $
! (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!===  GC dimensioning. The user may change this by preprocessor directives
!     On the command line when compiling GC. The default setting corresponds
!     to -DMAX_PROC=1024 -DMAX_COLL=1024 -DMAX_PT2PT=16384

!     Maximum number of processors



!     Maximum size of the arrays used in collective operations



!     Maximum size of the arrays sent/received in point-to-point communication




!=== End of GC dimensioning. Do NOT change anything below this point !


!     Default architecture flags setup


!     Consistency checks. Beware: Check C and util source before changing !


!     Size of sync work arrays - used internally only







!     The node id used for IO


!     Message tags reserved by GC
!     BEWARE: limits are duplicated in the header files !




















!     INFO tags reserved by GC, used in SHMEM only



!     Derived quantities

!     PVM byte swap support (not supported with CRI systems)



      SUBROUTINE GC_INIT (PATH, ME, NPROC)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Initialize all (machine dependent) variables used in the
!     *  communication. Currently this only applies to a work array used
!     *  in the reduction function of shmem.
!     *
!     * Input:
!     *  NPROC   - number of nodes (PVM_SRC only)
!     *  PATH    - path name of executable to be spawned (PVM_SRC only)
!     *            used only if first character is non-blank.
!     *
!     * Output:
!     *  NPROC   - number of nodes 
!     *          - -1 if error occured (PVM_SRC only)
!     *  ME      - my node ID, 0..NPROC-1
!     *
!     * NOTES:       
!     *
!     * Initialize BCAST_SYNC_WRK, REDUCESYNC (SHMEM only).
!     ******************************************************************

      IMPLICIT NONE
      CHARACTER*(*) PATH
      INTEGER ME, NPROC
	integer GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED

      INCLUDE 'mpif.h'
      INTEGER INFO


      INTEGER I

      IF (INITED) RETURN


      CALL MPI_INIT(INFO)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, INFO)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, INFO)


      GC__NPROC = NPROC
      INITED = .TRUE.
      END


      SUBROUTINE GC_EXIT()
!     ******************************************************************
!     * Purpose:
!     *
!     *  Controlled cleanup of the parallel system.
!     *
!     * Input:
!     *
!     * Output:
!     *
!     * NOTES:
!     *  This is a dummy routine for SHMEM, NX and MPL.
!     *  
!     ******************************************************************
      IMPLICIT NONE


      INTEGER INFO


      CALL MPI_FINALIZE(INFO)


      END


      SUBROUTINE GC_ABORT (ME, NPROC, MESG)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Aborts program execution on all processors and clean up
!     *  the parallel system.
!     *
!     * Input:
!     *  ME      - my node ID, 0..NPROC-1
!     *  NPROC   - number of nodes (all except PVM_SRC)
!     *  MESG    - abort message to be printed on stdout
!     *
!     * Output:
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER ME, NPROC
      INTEGER I, INFO
      CHARACTER*(*) MESG

      INCLUDE 'mpif.h'


! Flush has different name with xlf compiler.
! Do we need it at all?
!      CALL FLUSH_(6)

      WRITE(*,*) 'gc_abort: ', MESG


      CALL MPI_ABORT(MPI_COMM_WORLD,9,INFO)


!     We should never reach this point ...
      CALL EXIT(1)
      END


      SUBROUTINE GC_RSEND (MSG, LEN, RECI, INFO, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send a real array from this processor to processor RECI.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  RECI    - receiver of the message
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  RARR    - name of the array on recieving processor
!     *            (SHM_SRC only)
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI,INFO
      real*8 RARR(LEN), SARR(LEN)
      INCLUDE 'mpif.h'


      CALL MPI_SEND(SARR, 8*LEN, MPI_BYTE, RECI, MSG,&
          MPI_COMM_WORLD, INFO)

      END

! NON-BLOCKING and WAIT routines added by Heiko Klein and Peter Wind

      SUBROUTINE GC_RSEND_NONBLOCK (REQ,MSG,LEN,RECI,INFO,RARR,SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send a real array from this processor to processor RECI.
!     *  with non-blocking IO, only implemented for MPI!!!
!     *
!     * Input:
!     *  REQ     - request
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  RECI    - receiver of the message
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  RARR    - name of the array on recieving processor
!     *            (SHM_SRC only)
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI,INFO
      INTEGER REQ(1)
      real*8 RARR(LEN), SARR(LEN)

      INCLUDE 'mpif.h'

      CALL MPI_ISEND(SARR, 8*LEN, MPI_BYTE, RECI, MSG,&
          MPI_COMM_WORLD, REQ, INFO)

      END


      SUBROUTINE GC_RRECV_NONBLOCK (REQ,MSG,LEN,SEND,INFO,RARR,SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Receive a real array from processor SEND. non-blocking
!     *  non-blocking only implemented for mpi!
!     *
!     * Input:
!     *  REQ     - request handler (integer(1))
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message (SEND = -1 means any processor)
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  SARR    - name of the array on the sending processor
!     *            (SHM_SRC only)
!     *
!     * Output:
!     *  RARR    - array to be received
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO
      INTEGER REQ(1)
      real*8 RARR(LEN), SARR(LEN)

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_IRECV(RARR, 8*LEN, MPI_BYTE, SEND, MSG,&
          MPI_COMM_WORLD, STATUS, REQ, INFO)


      END


      SUBROUTINE GC_WAIT (REQ,INFO)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Waiting for corresponding non-blocking send/receive to 
!     *  finish.
!     *  Only implemented for MPI!!!
!     *
!     * Output:
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

      INTEGER INFO
      INTEGER REQ(1)


      CALL MPI_WAIT(REQ, STATUS, INFO)


      END


      SUBROUTINE GC_RRECV (MSG, LEN, SEND, INFO, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Receive a real array from processor SEND.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message (SEND = -1 means any processor)
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  SARR    - name of the array on the sending processor
!     *            (SHM_SRC only)
!     *
!     * Output:
!     *  RARR    - array to be received
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO
      real*8 RARR(LEN), SARR(LEN)


      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)


      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, 8*LEN, MPI_BYTE, SEND, MSG,&
          MPI_COMM_WORLD, STATUS, INFO)


      END


      SUBROUTINE GC_ISEND (MSG, LEN, RECI, INFO, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send a integer array from this processor to processor RECI.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  RECI    - receiver of the message
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  RARR    - name of the array on recieving processor
!     *            (SHM_SRC only)
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  INFO    - status of send. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI, INFO, RARR(LEN), SARR(LEN)

      INCLUDE 'mpif.h'


      CALL MPI_SEND(SARR, 4*LEN, MPI_BYTE, RECI, MSG,&
          MPI_COMM_WORLD, INFO)

      END


      SUBROUTINE GC_IRECV (MSG, LEN, SEND, INFO, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Receive a integer array from processor SEND.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message (SEND = -1 means any processor)
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  SARR    - name of the array on the sending processor
!     *            (SHM_SRC only)
!     *
!     * Output:
!     *  RARR    - array to be received
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO, RARR(LEN), SARR(LEN)

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, 4*LEN, MPI_BYTE, SEND, MSG,&
          MPI_COMM_WORLD, STATUS, INFO)


      END


      SUBROUTINE GC_BSEND (MSG, LEN, RECI, INFO, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send a byte array from this processor to processor RECI.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of BYTES in message
!     *  RECI    - receiver of the message
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  RARR    - name of the array on recieving processor
!     *            (SHM_SRC only)
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI, INFO&
             ,RARR(LEN/4+1), SARR(LEN/4+1)

      INCLUDE 'mpif.h'


      CALL MPI_SEND(SARR, LEN, MPI_BYTE, RECI, MSG, MPI_COMM_WORLD,&
          INFO)


      END


      SUBROUTINE GC_BRECV (MSG, LEN, SEND, INFO, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Receive a byte array from processor SEND.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of BYTES in message
!     *  SEND    - sender of the message (SEND = -1 means any processor)
!     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
!     *            should be used (SHM_SRC only)
!     *  SARR    - name of the array on the sending processor
!     *            (SHM_SRC only)
!     *
!     * Output:
!     *  RARR    - array to be received
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO, RARR(LEN), SARR(LEN/4+1)

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)


      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, LEN, MPI_BYTE, SEND, MSG, MPI_COMM_WORLD,&
          STATUS, INFO)


      END


      SUBROUTINE GC_GSYNC (NPROC, INFO)                  
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Synchronize the processors. Mainly used in front of
!     *  (asynchronous) communication and in connection with timing.
!     *
!     * Input:
!     *  NPROC   - number of nodes 
!     *  INFO    - flag deciding if this is a syncronization in front
!     *            of a SHMEM_PUT (SHM_SRC only)
!     *
!     * Output:
!     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     * 
!     *  No node can continue execution before everybody have reached
!     *  this point.
!     *  On CRAY MPP the cache is flushed if INFO == -9998.
!     *  This feature should be used directly after a SHM_PUT.
!     *  
!     ******************************************************************

      IMPLICIT NONE
      INTEGER NPROC,INFO
      INCLUDE 'mpif.h'


      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)


      END


      SUBROUTINE GC_SSYNC (NPROC, INFO)                  
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Synchronize the processors. Only on the CRAY T3D if SHM_SRC.
!     *
!     * Input:
!     *  NPROC   - number of nodes 
!     *  INFO    - flag deciding if this is a syncronization in front
!     *            of a SHMEM_PUT (SHM_SRC only)
!     *
!     * Output:
!     *
!     * NOTES:       
!     * 
!     *  The cache is flushed if INFO == -9998.
!     *  This feature should be used directly after a SHM_PUT.
!     *  
!     ******************************************************************

      IMPLICIT NONE
      INTEGER NPROC, INFO


      END


      SUBROUTINE GC_RBCAST (MSG, LEN, SEND, NPROC, INFO, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Broadcast a real array to every processor.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message
!     *  NPROC   - Number of processors
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  SARR    - array to be received (on nodes != SEND)
!     *  INFO    - status of bcast. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO
      real*8 SARR(LEN)

      INCLUDE 'mpif.h'


      CALL MPI_BCAST(SARR, 8*LEN, MPI_BYTE, SEND,&
          MPI_COMM_WORLD, INFO)


      END


      SUBROUTINE GC_IBCAST (MSG, LEN, SEND, NPROC, INFO, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Broadcast an integer array to every processor.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message
!     *  NPROC   - Number of processors
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  SARR    - array to be received (on nodes != SEND)
!     *  INFO    - status of bcast. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO, SARR(LEN)

      INCLUDE 'mpif.h'


      CALL MPI_BCAST(SARR, 4*LEN, MPI_BYTE, SEND,&
          MPI_COMM_WORLD, INFO)



      END


      SUBROUTINE GC_BBCAST (MSG, LEN, SEND, NPROC, INFO, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Broadcast a byte array to every processor.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of BYTES in message
!     *  SEND    - sender of the message
!     *  NPROC   - Number of processors
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  SARR    - array to be received (on nodes != SEND)
!     *  INFO    - status of bcast. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO, SARR(LEN/4+1)
      INCLUDE 'mpif.h'


      CALL MPI_BCAST(SARR, LEN, MPI_BYTE, SEND, MPI_COMM_WORLD, INFO)

      END


      SUBROUTINE GC_RSUM (LEN, NPROC, INFO, SSUM)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate the real sum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  SSUM    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  SSUM    - array containing the sums across the nodes
!     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      real*8 SSUM(LEN)
      real*8 REDUCE_DATA_WRK(4096)

      INCLUDE 'mpif.h'
      INTEGER I

      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SSUM(I)
      ENDDO

      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SSUM, LEN, &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)


      END


      SUBROUTINE GC_RMIN (LEN, NPROC, INFO, SMIN)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Finds the real minimum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  SMIN    - array with elements of which the elementwise minimum
!     *            across the nodes is to be found
!     *
!     * Output:
!     *  SMIN    - array containing the minimums across the nodes
!     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      real*8 SMIN(LEN)
      real*8 REDUCE_DATA_WRK(4096)

      INCLUDE 'mpif.h'
      INTEGER I

	DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMIN(I)
      ENDDO

      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMIN, LEN, &
          MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO)

      END


      SUBROUTINE GC_RMAX (LEN, NPROC, INFO, SMAX)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Finds the real maximum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  SMAX    - array with elements of which the elementwise maximum
!     *            across the nodes is to be found
!     *
!     * Output:
!     *  SMAX    - array containing the maximums across the nodes
!     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      real*8 SMAX(LEN)
      real*8 REDUCE_DATA_WRK(4096)

      INCLUDE 'mpif.h'
      INTEGER I


      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMAX(I)
      ENDDO

      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMAX, LEN, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, INFO)

      END


      SUBROUTINE GC_ISUM (LEN, NPROC, INFO, ISUM)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate the integer sum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  ISUM    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  ISUM    - array containing the sums across the nodes
!     *  INFO    - status of isum. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, ISUM(LEN)
      INTEGER REDUCE_DATA_IWRK(4096)

      INCLUDE 'mpif.h'
      INTEGER I

      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = ISUM(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, ISUM, LEN, MPI_INTEGER,&
          MPI_SUM, MPI_COMM_WORLD, INFO)


      END


      SUBROUTINE GC_IMIN (LEN, NPROC, INFO, IMIN)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Finds the real minimum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  IMIN    - array with elements of which the elementwise minimum
!     *            across the nodes is to be found
!     *
!     * Output:
!     *  IMIN    - array containing the minimums across the nodes
!     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, IMIN(LEN)
      INTEGER REDUCE_DATA_IWRK(4096)
      INCLUDE 'mpif.h'
      INTEGER I


      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMIN(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, IMIN, LEN, MPI_INTEGER,&
          MPI_MIN, MPI_COMM_WORLD, INFO)

      END


      SUBROUTINE GC_IMAX (LEN, NPROC, INFO, IMAX)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Finds the integer maximum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  IMAX    - array with elements of which the elementwise maximum
!     *            across the nodes is to be found
!     *
!     * Output:
!     *  IMAX    - array containing the maximums across the nodes
!     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, IMAX(LEN)
      INTEGER REDUCE_DATA_IWRK(4096)

      INCLUDE 'mpif.h'
      INTEGER I

      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMAX(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, IMAX, LEN, MPI_INTEGER,&
          MPI_MAX, MPI_COMM_WORLD, INFO)

      END


      INTEGER FUNCTION GC_RSIZE ()
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Return the number of bytes in a real variable
!     *
!     * Input:
!     *
!     * Output:
!     *
!     * NOTES:       
!     * 
!     ******************************************************************

      IMPLICIT NONE


      GC_RSIZE = 8
    
      END


      INTEGER FUNCTION GC_ISIZE ()
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Return the number of bytes in an integer variable
!     *
!     * Input:
!     *
!     * Output:
!     *
!     * NOTES:       
!     * 
!     ******************************************************************

      IMPLICIT NONE


      GC_ISIZE = 4

    
      END


      INTEGER FUNCTION GC_COMLEN (LTYPE, LLEN, LAST, FIRST)
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Return the number of bytes in a common block
!     *
!     * Input:
!     *  LTYPE     - data type of the last variable, 'I' or 'R'
!     *  LLEN      - number of elements in the last variable
!     *  LAST      -  variable of the COMMON BLOCK
!     *  FIRST     - first variable of the COMMON BLOCK
!     *
!     * Output:
!     *  GC_COMLEN - length of COMMON block in bytes, or 0 if error
!     *              in LTYPE
!     *
!     * NOTES:
!     *  - This routine does NOT work with CRI MPP systems where the
!     *    LAST or FIRST arguments are CHARACTER variables. Further on,
!     *    it is subject to corrupt the CRI MPP call stack if used this
!     *    way.
!     ******************************************************************

      IMPLICIT NONE
      CHARACTER*(*) LTYPE
      INTEGER LLEN, LAST(1), FIRST(1)

      INTEGER LOC, NBYTES



!     Systems with UNIX Fortran character argument passing
      NBYTES = LOC(LAST) - LOC(FIRST)



      IF (LTYPE(1:1) .EQ. 'r' .OR. LTYPE(1:1) .EQ. 'R') THEN
         NBYTES = NBYTES + 8*LLEN
      ELSE IF (LTYPE(1:1) .EQ. 'i' .OR. LTYPE(1:1) .EQ. 'I') THEN
         NBYTES = NBYTES + 4*LLEN
      ELSE IF (LTYPE(1:1) .EQ. 'c' .OR. LTYPE(1:1) .EQ. 'C') THEN
         NBYTES = NBYTES + LLEN
      ELSE
!        Error in LTYPE
         NBYTES = 0
      ENDIF
      
      GC_COMLEN = NBYTES
      END


      SUBROUTINE GC_CONFIG (MXPROC, MXCOLL, MXPT2PT, INTF)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Return information about the GC configuration.
!     *
!     * Output:
!     *  MXPROC    - maximum numbers of processors compiled into the
!     *              interface
!     *  MXCOLL    - maximum number of elements for collective
!     *              operations
!     *  MXPT2PT   - maximum number of elements for point to point
!     *              operations
!     *  INTF      - name of interface selected at compile time
!     *
!     * NOTES:
!     *    
!     ******************************************************************

      IMPLICIT NONE
      INTEGER MXPROC, MXCOLL, MXPT2PT
      CHARACTER*(*) INTF


      MXPROC = 64
      MXCOLL = 4096
      MXPT2PT = 16384
      INTF = 'GC_MPI'
      END


      BLOCK DATA

      INTEGER GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED
      DATA GC__NPROC/-1/, INITED/.FALSE./

      END


!===============================================================================
!  This part contains INTERNAL routines to be used within the GC
!  interface ONLY.
!===============================================================================

      CHARACTER*(*) FUNCTION GC__SUFFIX(PATH)
!
!     Support function returning file name part of a slash separated path
!
      INTEGER I, PEND
      CHARACTER*(*) PATH

      PEND = LEN(PATH)
      I = PEND
 100  CONTINUE
        IF (PATH(I:I) .EQ. '/') THEN
           GC__SUFFIX = PATH(I+1:PEND)
           RETURN
        ENDIF

        I = I - 1
      IF ( I .GT. 1) GOTO 100

      GC__SUFFIX = PATH(1:PEND)
      END


      INTEGER FUNCTION GC__FCONFIG(NPROC, INTFID)
!
!     Support function for mixed C and Fortran GC library calls
!
      IMPLICIT NONE
      INTEGER NPROC, INTFID, GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED


      IF (.NOT. INITED) THEN
         GC__FCONFIG = -1
         RETURN
      ENDIF
      NPROC = GC__NPROC
      INTFID = 4
      GC__FCONFIG = 0
      END


      INTEGER FUNCTION GC__FPVMCONFIG(NPROC, FTIDS)
!
!     Support function for mixed C and Fortran GC library calls
!
      IMPLICIT NONE
      INTEGER NPROC, FTIDS(0:2*NPROC-1), GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED

      IF (.NOT. INITED) THEN
         GC__FPVMCONFIG = -1
         RETURN
      ENDIF

      GC__FPVMCONFIG = 0
      END


      INTEGER FUNCTION GC__C2F(NPROC, INTFID, CTIDS)
!
!     Support function for mixed C and Fortran GC library calls
!
      IMPLICIT NONE
      INTEGER NPROC, INTFID, CTIDS(0:NPROC-1)
      INTEGER GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED

      IF (INITED) THEN
         GC__C2F = 0
         RETURN
      ELSE IF (NPROC .GT. 64&
             .OR. INTFID .NE. 4) THEN
         GC__C2F = -1
         RETURN
      ELSE
         GC__C2F = 0
      ENDIF
      GC__NPROC = NPROC

      INITED = .TRUE.
      END

