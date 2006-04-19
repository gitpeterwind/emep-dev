# 1 "gc_com.F"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GC - General Communication primitives package. For use on multiprocessor
C shared memory and message passing systems.
C
C
C LICENSING TERMS
C
C  GC is provided free of charge. Unless otherwise agreed with SINTEF, use and
C  redistribution in source and binary forms are permitted provided that
C
C      (1) source distributions retain all comments appearing within this file
C          header, and
C
C      (2) distributions including binaries display the following
C          acknowledgement:
C
C              "This product includes software developed by SINTEF.",
C
C          in the documentation or other materials provided with the
C          distribution and in all advertising materials mentioning features or
C          use of this software.
C
C  The name of SINTEF may not be used to endorse or promote products derived
C  from this software without specific prior written permission.  SINTEF
C  disclaims any warranty that this software will be fit for any specific
C  purposes. In no event shall SINTEF be liable for any loss of performance or
C  for indirect or consequential damage or direct or indirect injury of any
C  kind. In no case shall SINTEF be liable for any representation or warranty
C  make to any third party by the users of this software.
C
C USAGE
C
C  Define zero or ONE of the INTERFACE flags below with -D<name> and compile
C  using the C (or Fortran) preprocessor.
C
C INTERFACE FLAGS
C
C  SHM_SRC  -  CRAY MPP systems CRI using shared memory (SHMEM)
C  NX2_SRC  -  Intel Paragon using Intel/NX2 message passing
C  PVM_SRC  -  Public ORNL PVM/Cray PVP  message passing
C  1  -  The Message Passing interface
C  MPL_SRC  -  The IBM Message Passing Language
C
C
C INTERFACE CONFIGURATION FLAGS
C
C  Define zero or more of the flags below with -D<name> and compile using the
C  C (or Fortran) preprocessor.
C
C  FLP_32B  -  use 32 bit floating point precision. Enabled by default for
C              all systems except CRI PVM & MPP.
C  1  -  use 64 bit floating point precision. Enabled by default on CRI
C              systems.
C  PVM_T3D  -  enable CRAY MPP PVM specifics (for backwards compatibility only,
C              this functionality is always enabled on CRI MPP systems using
C              PVM_SRC).
C  PVM_V33  -  has PVM v3.3 or higher, with psend/precv.
C  PVM_NBS  -  PVM No Byte Swap - disable byte sex swapping. Byte swapping
C              of INTEGER and REAl data is on by default.
C  PVM_DIP  -  PVM Data In Place, enabled with ORNL PVM v3.3 and higher.
C              Do not copy message data from user space into PVM buffers
C              on send. Data must NOT be modified before leaving the sending
C              node (before arriving receiving node with CRI MPP PVM).
C
C  NX2_32B, NX2_64B, PVM_32B, PVM_64B, MPI_32B, MPI_64B are supported for
C  backward compatibility with prerelease versions only. Use the FLP_nnB flags.
C
C
C REVISION HISTORY
C
C  1.0 - Initial version based on the GCOM interface by R.Skaalin and
C        the PAR interface by J.Amundsen. The MPL part is due to Z. Christidis,
C        IBM T. J. Watson Research Center.
C
C-------------------------------------------------------------------------------
C $Id: gc_com.f,v 1.1 2006-04-19 08:10:09 mifapw Exp $
C (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C===  GC dimensioning. The user may change this by preprocessor directives
C     On the command line when compiling GC. The default setting corresponds
C     to -DMAX_PROC=1024 -DMAX_COLL=1024 -DMAX_PT2PT=16384

C     Maximum number of processors



C     Maximum size of the arrays used in collective operations



C     Maximum size of the arrays sent/received in point-to-point communication




C=== End of GC dimensioning. Do NOT change anything below this point !


C     Default architecture flags setup
# 106

# 112

# 115

# 118


# 122



# 127


C     Backwards compability
# 132

# 135

# 138

# 141

# 144


C     Consistency checks. Beware: Check C and util source before changing !



# 157

# 165

# 173


# 177





# 189

# 194

# 197



C     Size of sync work arrays - used internally only







C     The node id used for IO


C     Message tags reserved by GC
C     BEWARE: limits are duplicated in the header files !




















C     INFO tags reserved by GC, used in SHMEM only



C     Derived quantities
# 240







# 250

      
# 256






# 266






C     PVM byte swap support (not supported with CRI systems)




# 279



      SUBROUTINE GC_INIT (PATH, ME, NPROC)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Initialize all (machine dependent) variables used in the
C     *  communication. Currently this only applies to a work array used
C     *  in the reduction function of shmem.
C     *
C     * Input:
C     *  NPROC   - number of nodes (PVM_SRC only)
C     *  PATH    - path name of executable to be spawned (PVM_SRC only)
C     *            used only if first character is non-blank.
C     *
C     * Output:
C     *  NPROC   - number of nodes 
C     *          - -1 if error occured (PVM_SRC only)
C     *  ME      - my node ID, 0..NPROC-1
C     *
C     * NOTES:       
C     *
C     * Initialize BCAST_SYNC_WRK, REDUCESYNC (SHMEM only).
C     ******************************************************************

      IMPLICIT NONE
      CHARACTER*(*) PATH
      INTEGER ME, NPROC
	integer GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED
# 321

# 325

# 339


      INCLUDE 'mpif.h'
      INTEGER INFO

# 347


      INTEGER I

      IF (INITED) RETURN

# 366


# 371


# 485



      CALL MPI_INIT(INFO)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, INFO)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, INFO)


# 496


# 501


      GC__NPROC = NPROC
      INITED = .TRUE.
      END


      SUBROUTINE GC_EXIT()
C     ******************************************************************
C     * Purpose:
C     *
C     *  Controlled cleanup of the parallel system.
C     *
C     * Input:
C     *
C     * Output:
C     *
C     * NOTES:
C     *  This is a dummy routine for SHMEM, NX and MPL.
C     *  
C     ******************************************************************
      IMPLICIT NONE


      INTEGER INFO



# 531



      CALL MPI_FINALIZE(INFO)


      END


      SUBROUTINE GC_ABORT (ME, NPROC, MESG)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Aborts program execution on all processors and clean up
C     *  the parallel system.
C     *
C     * Input:
C     *  ME      - my node ID, 0..NPROC-1
C     *  NPROC   - number of nodes (all except PVM_SRC)
C     *  MESG    - abort message to be printed on stdout
C     *
C     * Output:
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER ME, NPROC
      INTEGER I, INFO
      CHARACTER*(*) MESG
# 564

# 569


      INCLUDE 'mpif.h'



      CALL FLUSH(6)

      WRITE(*,*) 'gc_abort: ', MESG

# 582


# 586


# 594



      CALL MPI_ABORT(MPI_COMM_WORLD,9,INFO)


# 602


C     We should never reach this point ...
      CALL EXIT(1)
      END


      SUBROUTINE GC_RSEND (MSG, LEN, RECI, INFO, RARR, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Send a real array from this processor to processor RECI.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  RECI    - receiver of the message
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  RARR    - name of the array on recieving processor
C     *            (SHM_SRC only)
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI,INFO
      real*8 RARR(LEN), SARR(LEN)
# 638

# 643


      INCLUDE 'mpif.h'

# 650



# 655


# 659


# 669



      CALL MPI_SEND(SARR, 8*LEN, MPI_BYTE, RECI, MSG,
     *     MPI_COMM_WORLD, INFO)


# 678


      END

C NON-BLOCKING and WAIT routines added by Heiko Klein and Peter Wind

      SUBROUTINE GC_RSEND_NONBLOCK (REQ,MSG,LEN,RECI,INFO,RARR,SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Send a real array from this processor to processor RECI.
C     *  with non-blocking IO, only implemented for MPI!!!
C     *
C     * Input:
C     *  REQ     - request
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  RECI    - receiver of the message
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  RARR    - name of the array on recieving processor
C     *            (SHM_SRC only)
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI,INFO
      INTEGER REQ(1)
      real*8 RARR(LEN), SARR(LEN)
# 716

# 721


      INCLUDE 'mpif.h'

# 728



# 733


# 737


# 747



      CALL MPI_ISEND(SARR, 8*LEN, MPI_BYTE, RECI, MSG,
     *     MPI_COMM_WORLD, REQ, INFO)


# 756


      END


      SUBROUTINE GC_RRECV_NONBLOCK (REQ,MSG,LEN,SEND,INFO,RARR,SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Receive a real array from processor SEND. non-blocking
C     *  non-blocking only implemented for mpi!
C     *
C     * Input:
C     *  REQ     - request handler (integer(1))
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  SEND    - sender of the message (SEND = -1 means any processor)
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  SARR    - name of the array on the sending processor
C     *            (SHM_SRC only)
C     *
C     * Output:
C     *  RARR    - array to be received
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO
      INTEGER REQ(1)
      real*8 RARR(LEN), SARR(LEN)
# 793

# 801


      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

# 809



# 814


# 818


# 848



      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_IRECV(RARR, 8*LEN, MPI_BYTE, SEND, MSG,
     *     MPI_COMM_WORLD, STATUS, REQ, INFO)


# 859


      END


      SUBROUTINE GC_WAIT (REQ,INFO)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Waiting for corresponding non-blocking send/receive to 
C     *  finish.
C     *  Only implemented for MPI!!!
C     *
C     * Output:
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

      INTEGER INFO
      INTEGER REQ(1)


      CALL MPI_WAIT(REQ, STATUS, INFO)


      END


      SUBROUTINE GC_RRECV (MSG, LEN, SEND, INFO, RARR, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Receive a real array from processor SEND.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  SEND    - sender of the message (SEND = -1 means any processor)
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  SARR    - name of the array on the sending processor
C     *            (SHM_SRC only)
C     *
C     * Output:
C     *  RARR    - array to be received
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO
      real*8 RARR(LEN), SARR(LEN)
# 924

# 932


      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

# 940



# 945


# 949


# 979



      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, 8*LEN, MPI_BYTE, SEND, MSG,
     *     MPI_COMM_WORLD, STATUS, INFO)


# 990


      END


      SUBROUTINE GC_ISEND (MSG, LEN, RECI, INFO, RARR, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Send a integer array from this processor to processor RECI.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  RECI    - receiver of the message
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  RARR    - name of the array on recieving processor
C     *            (SHM_SRC only)
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  INFO    - status of send. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI, INFO, RARR(LEN), SARR(LEN)
# 1023

# 1028


      INCLUDE 'mpif.h'

# 1035



# 1040


# 1044


# 1054



      CALL MPI_SEND(SARR, 4*LEN, MPI_BYTE, RECI, MSG,
     *     MPI_COMM_WORLD, INFO)


# 1063


      END


      SUBROUTINE GC_IRECV (MSG, LEN, SEND, INFO, RARR, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Receive a integer array from processor SEND.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  SEND    - sender of the message (SEND = -1 means any processor)
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  SARR    - name of the array on the sending processor
C     *            (SHM_SRC only)
C     *
C     * Output:
C     *  RARR    - array to be received
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO, RARR(LEN), SARR(LEN)
# 1096

# 1104


      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

# 1112



# 1117


# 1121


# 1151



      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, 4*LEN, MPI_BYTE, SEND, MSG,
     *     MPI_COMM_WORLD, STATUS, INFO)


# 1162


      END


      SUBROUTINE GC_BSEND (MSG, LEN, RECI, INFO, RARR, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Send a byte array from this processor to processor RECI.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of BYTES in message
C     *  RECI    - receiver of the message
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  RARR    - name of the array on recieving processor
C     *            (SHM_SRC only)
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI, INFO
     *        ,RARR(LEN/4+1), SARR(LEN/4+1)
# 1196

# 1201


      INCLUDE 'mpif.h'

# 1208



# 1214


# 1218


# 1228



      CALL MPI_SEND(SARR, LEN, MPI_BYTE, RECI, MSG, MPI_COMM_WORLD,
     *     INFO)


# 1237


      END


      SUBROUTINE GC_BRECV (MSG, LEN, SEND, INFO, RARR, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Receive a byte array from processor SEND.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of BYTES in message
C     *  SEND    - sender of the message (SEND = -1 means any processor)
C     *  INFO    - flag deciding if SHMEM_GET (default) or SHMEM_PUT 
C     *            should be used (SHM_SRC only)
C     *  SARR    - name of the array on the sending processor
C     *            (SHM_SRC only)
C     *
C     * Output:
C     *  RARR    - array to be received
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO, RARR(LEN), SARR(LEN/4+1)
# 1270

# 1278


      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)

# 1286



# 1292


# 1296


# 1306



      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, LEN, MPI_BYTE, SEND, MSG, MPI_COMM_WORLD,
     *     STATUS, INFO)


# 1317


      END


      SUBROUTINE GC_GSYNC (NPROC, INFO)                  
C     ******************************************************************
C     * Purpose:
C     *  
C     *  Synchronize the processors. Mainly used in front of
C     *  (asynchronous) communication and in connection with timing.
C     *
C     * Input:
C     *  NPROC   - number of nodes 
C     *  INFO    - flag deciding if this is a syncronization in front
C     *            of a SHMEM_PUT (SHM_SRC only)
C     *
C     * Output:
C     *  INFO    - status of send 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     * 
C     *  No node can continue execution before everybody have reached
C     *  this point.
C     *  On CRAY MPP the cache is flushed if INFO == -9998.
C     *  This feature should be used directly after a SHM_PUT.
C     *  
C     ******************************************************************

      IMPLICIT NONE
      INTEGER NPROC,INFO
# 1354


      INCLUDE 'mpif.h'

# 1361


# 1366


# 1370


# 1394



      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)


# 1402


      END


      SUBROUTINE GC_SSYNC (NPROC, INFO)                  
C     ******************************************************************
C     * Purpose:
C     *  
C     *  Synchronize the processors. Only on the CRAY T3D if SHM_SRC.
C     *
C     * Input:
C     *  NPROC   - number of nodes 
C     *  INFO    - flag deciding if this is a syncronization in front
C     *            of a SHMEM_PUT (SHM_SRC only)
C     *
C     * Output:
C     *
C     * NOTES:       
C     * 
C     *  The cache is flushed if INFO == -9998.
C     *  This feature should be used directly after a SHM_PUT.
C     *  
C     ******************************************************************

      IMPLICIT NONE
      INTEGER NPROC, INFO

# 1433


      END


      SUBROUTINE GC_RBCAST (MSG, LEN, SEND, NPROC, INFO, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Broadcast a real array to every processor.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  SEND    - sender of the message
C     *  NPROC   - Number of processors
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  SARR    - array to be received (on nodes != SEND)
C     *  INFO    - status of bcast. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO
      real*8 SARR(LEN)
# 1470

# 1474

# 1483


      INCLUDE 'mpif.h'

# 1490



# 1497


# 1505


# 1529



      CALL MPI_BCAST(SARR, 8*LEN, MPI_BYTE, SEND,
     *     MPI_COMM_WORLD, INFO)


# 1538


      END


      SUBROUTINE GC_IBCAST (MSG, LEN, SEND, NPROC, INFO, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Broadcast an integer array to every processor.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of elements in message
C     *  SEND    - sender of the message
C     *  NPROC   - Number of processors
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  SARR    - array to be received (on nodes != SEND)
C     *  INFO    - status of bcast. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO, SARR(LEN)
# 1574

# 1578

# 1587


      INCLUDE 'mpif.h'

# 1594



# 1601


# 1609


# 1633



      CALL MPI_BCAST(SARR, 4*LEN, MPI_BYTE, SEND,
     *     MPI_COMM_WORLD, INFO)


# 1642


      END


      SUBROUTINE GC_BBCAST (MSG, LEN, SEND, NPROC, INFO, SARR)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Broadcast a byte array to every processor.
C     *
C     * Input:
C     *  MSG     - message tag
C     *  LEN     - number of BYTES in message
C     *  SEND    - sender of the message
C     *  NPROC   - Number of processors
C     *  SARR    - array to be sent
C     *
C     * Output:
C     *  SARR    - array to be received (on nodes != SEND)
C     *  INFO    - status of bcast. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO, SARR(LEN/4+1)
# 1678

# 1682

# 1691


      INCLUDE 'mpif.h'

# 1698



# 1705


# 1713


# 1734



      CALL MPI_BCAST(SARR, LEN, MPI_BYTE, SEND, MPI_COMM_WORLD, INFO)


# 1742


      END


      SUBROUTINE GC_RSUM (LEN, NPROC, INFO, SSUM)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Calculate the real sum across all processors and distribute
C     *  the result to all the processors.
C     *
C     * Input:
C     *  LEN     - number of elements in message
C     *  NPROC   - number of processors
C     *  SSUM    - array with elements to be added up across the nodes
C     *
C     * Output:
C     *  SSUM    - array containing the sums across the nodes
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      real*8 SSUM(LEN)
      real*8 REDUCE_DATA_WRK(4096)
# 1778

# 1781

# 1787

# 1790


      INCLUDE 'mpif.h'
      INTEGER I

# 1799



# 1806


# 1810


# 1872



      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SSUM(I)
      ENDDO

      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SSUM, LEN, 
     *     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
# 1884



# 1893

      END


      SUBROUTINE GC_RMIN (LEN, NPROC, INFO, SMIN)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Finds the real minimum across all processors and distribute
C     *  the result to all the processors.
C     *
C     * Input:
C     *  LEN     - number of elements in message
C     *  NPROC   - number of processors
C     *  SMIN    - array with elements of which the elementwise minimum
C     *            across the nodes is to be found
C     *
C     * Output:
C     *  SMIN    - array containing the minimums across the nodes
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      real*8 SMIN(LEN)
      real*8 REDUCE_DATA_WRK(4096)
# 1929

# 1932

# 1938

# 1941


      INCLUDE 'mpif.h'
      INTEGER I

# 1950



# 1957


# 1961


# 2023



	DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMIN(I)
      ENDDO

      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMIN, LEN, 
     *     MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO)
# 2035



# 2044


      END


      SUBROUTINE GC_RMAX (LEN, NPROC, INFO, SMAX)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Finds the real maximum across all processors and distribute
C     *  the result to all the processors.
C     *
C     * Input:
C     *  LEN     - number of elements in message
C     *  NPROC   - number of processors
C     *  SMAX    - array with elements of which the elementwise maximum
C     *            across the nodes is to be found
C     *
C     * Output:
C     *  SMAX    - array containing the maximums across the nodes
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      real*8 SMAX(LEN)
      real*8 REDUCE_DATA_WRK(4096)
# 2081

# 2084

# 2090

# 2093


      INCLUDE 'mpif.h'
      INTEGER I

# 2102



# 2109


# 2113


# 2175



      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMAX(I)
      ENDDO

      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMAX, LEN, 
     *     MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, INFO)
# 2187



# 2196


      END


      SUBROUTINE GC_ISUM (LEN, NPROC, INFO, ISUM)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Calculate the integer sum across all processors and distribute
C     *  the result to all the processors.
C     *
C     * Input:
C     *  LEN     - number of elements in message
C     *  NPROC   - number of processors
C     *  ISUM    - array with elements to be added up across the nodes
C     *
C     * Output:
C     *  ISUM    - array containing the sums across the nodes
C     *  INFO    - status of isum. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, ISUM(LEN)
      INTEGER REDUCE_DATA_IWRK(4096)
# 2231

# 2234

# 2240

# 2243


      INCLUDE 'mpif.h'
      INTEGER I

# 2252



# 2259


# 2263


# 2325



      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = ISUM(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, ISUM, LEN, MPI_INTEGER,
     *     MPI_SUM, MPI_COMM_WORLD, INFO)


# 2341

      END


      SUBROUTINE GC_IMIN (LEN, NPROC, INFO, IMIN)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Finds the real minimum across all processors and distribute
C     *  the result to all the processors.
C     *
C     * Input:
C     *  LEN     - number of elements in message
C     *  NPROC   - number of processors
C     *  IMIN    - array with elements of which the elementwise minimum
C     *            across the nodes is to be found
C     *
C     * Output:
C     *  IMIN    - array containing the minimums across the nodes
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, IMIN(LEN)
      INTEGER REDUCE_DATA_IWRK(4096)
# 2376

# 2379

# 2385

# 2388


      INCLUDE 'mpif.h'
      INTEGER I

# 2397



# 2404


# 2408


# 2470



      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMIN(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, IMIN, LEN, MPI_INTEGER,
     *     MPI_MIN, MPI_COMM_WORLD, INFO)


# 2486


      END


      SUBROUTINE GC_IMAX (LEN, NPROC, INFO, IMAX)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Finds the integer maximum across all processors and distribute
C     *  the result to all the processors.
C     *
C     * Input:
C     *  LEN     - number of elements in message
C     *  NPROC   - number of processors
C     *  IMAX    - array with elements of which the elementwise maximum
C     *            across the nodes is to be found
C     *
C     * Output:
C     *  IMAX    - array containing the maximums across the nodes
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and 1 only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, IMAX(LEN)
      INTEGER REDUCE_DATA_IWRK(4096)
# 2522

# 2525

# 2531

# 2534


      INCLUDE 'mpif.h'
      INTEGER I

# 2543



# 2550


# 2554


# 2616



      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMAX(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, IMAX, LEN, MPI_INTEGER,
     *     MPI_MAX, MPI_COMM_WORLD, INFO)


# 2632


      END


      INTEGER FUNCTION GC_RSIZE ()
C     ******************************************************************
C     * Purpose:
C     *  
C     *  Return the number of bytes in a real variable
C     *
C     * Input:
C     *
C     * Output:
C     *
C     * NOTES:       
C     * 
C     ******************************************************************

      IMPLICIT NONE


      GC_RSIZE = 8
# 2657

    
      END


      INTEGER FUNCTION GC_ISIZE ()
C     ******************************************************************
C     * Purpose:
C     *  
C     *  Return the number of bytes in an integer variable
C     *
C     * Input:
C     *
C     * Output:
C     *
C     * NOTES:       
C     * 
C     ******************************************************************

      IMPLICIT NONE

# 2680

      GC_ISIZE = 4

    
      END


      INTEGER FUNCTION GC_COMLEN (LTYPE, LLEN, LAST, FIRST)
C     ******************************************************************
C     * Purpose:
C     *  
C     *  Return the number of bytes in a common block
C     *
C     * Input:
C     *  LTYPE     - data type of the last variable, 'I' or 'R'
C     *  LLEN      - number of elements in the last variable
C     *  LAST      -  variable of the COMMON BLOCK
C     *  FIRST     - first variable of the COMMON BLOCK
C     *
C     * Output:
C     *  GC_COMLEN - length of COMMON block in bytes, or 0 if error
C     *              in LTYPE
C     *
C     * NOTES:
C     *  - This routine does NOT work with CRI MPP systems where the
C     *    LAST or FIRST arguments are CHARACTER variables. Further on,
C     *    it is subject to corrupt the CRI MPP call stack if used this
C     *    way.
C     ******************************************************************

      IMPLICIT NONE
      CHARACTER*(*) LTYPE
      INTEGER LLEN, LAST(1), FIRST(1)
# 2715

      INTEGER LOC, NBYTES


# 2729

# 2742

C     Systems with UNIX Fortran character argument passing
      NBYTES = LOC(LAST) - LOC(FIRST)



      IF (LTYPE(1:1) .EQ. 'r' .OR. LTYPE(1:1) .EQ. 'R') THEN
         NBYTES = NBYTES + 8*LLEN
      ELSE IF (LTYPE(1:1) .EQ. 'i' .OR. LTYPE(1:1) .EQ. 'I') THEN
         NBYTES = NBYTES + 4*LLEN
      ELSE IF (LTYPE(1:1) .EQ. 'c' .OR. LTYPE(1:1) .EQ. 'C') THEN
         NBYTES = NBYTES + LLEN
      ELSE
C        Error in LTYPE
         NBYTES = 0
      ENDIF
      
      GC_COMLEN = NBYTES
      END


      SUBROUTINE GC_CONFIG (MXPROC, MXCOLL, MXPT2PT, INTF)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Return information about the GC configuration.
C     *
C     * Output:
C     *  MXPROC    - maximum numbers of processors compiled into the
C     *              interface
C     *  MXCOLL    - maximum number of elements for collective
C     *              operations
C     *  MXPT2PT   - maximum number of elements for point to point
C     *              operations
C     *  INTF      - name of interface selected at compile time
C     *
C     * NOTES:
C     *    
C     ******************************************************************

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


C===============================================================================
C  This part contains INTERNAL routines to be used within the GC
C  interface ONLY.
C===============================================================================

      CHARACTER*(*) FUNCTION GC__SUFFIX(PATH)
C
C     Support function returning file name part of a slash separated path
C
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
C
C     Support function for mixed C and Fortran GC library calls
C
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
C
C     Support function for mixed C and Fortran GC library calls
C
      IMPLICIT NONE
      INTEGER NPROC, FTIDS(0:2*NPROC-1), GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED
# 2862



      IF (.NOT. INITED) THEN
         GC__FPVMCONFIG = -1
         RETURN
      ENDIF
# 2874


      GC__FPVMCONFIG = 0
      END


      INTEGER FUNCTION GC__C2F(NPROC, INTFID, CTIDS)
C
C     Support function for mixed C and Fortran GC library calls
C
      IMPLICIT NONE
      INTEGER NPROC, INTFID, CTIDS(0:NPROC-1)
      INTEGER GC__NPROC
      LOGICAL INITED
      COMMON /GC__INIT/ GC__NPROC, INITED
# 2897

# 2901



      IF (INITED) THEN
         GC__C2F = 0
         RETURN
      ELSE IF (NPROC .GT. 64
     *        .OR. INTFID .NE. 4) THEN
         GC__C2F = -1
         RETURN
      ELSE
         GC__C2F = 0
      ENDIF
      GC__NPROC = NPROC

# 2927


# 2934


      INITED = .TRUE.
      END


# 3004

