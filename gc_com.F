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
C  MPI_SRC  -  The Message Passing interface
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
C  FLP_64B  -  use 64 bit floating point precision. Enabled by default on CRI
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
C $Id: gc_com.F,v 1.4 2002-09-13 08:55:39 mifads Exp $
C (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C===  GC dimensioning. The user may change this by preprocessor directives
C     On the command line when compiling GC. The default setting corresponds
C     to -DMAX_PROC=1024 -DMAX_COLL=1024 -DMAX_PT2PT=16384

C     Maximum number of processors
#if !defined(MAX_PROC)
#define MAX_PROC 64
#endif
C     Maximum size of the arrays used in collective operations
#if !defined(MAX_COLL)
#define MAX_COLL 4096      
#endif
C     Maximum size of the arrays sent/received in point-to-point communication
#if !defined(MAX_PT2PT)
#define MAX_PT2PT 16384
#endif

C=== End of GC dimensioning. Do NOT change anything below this point !


C     Default architecture flags setup
#if defined(_CRAYMPP) && defined(PVM_SRC) && !defined(PVM_T3D)
C     PVM_T3D is defined for backwards compatibility only
#define PVM_T3D
#endif
#if defined(_CRAY)
#define FLP_64B
#if defined(FLP_32B)
#undef FLP_32B
#endif
#else
#if defined(FLP_64B) && defined(FLP_32B)
#undef FLP_32B
#endif
#if !defined(FLP_64B) && !defined(FLP_32B)
#define FLP_32B
#endif
#endif
#if defined (PVM_DIP)
#define GC__PVMENCODING PVMINPLACE
#else
#define GC__PVMENCODING PVMDEFAULT
#endif
#if defined(_AIX)
#define FLUSH FLUSH_
#endif

C     Backwards compability
#if defined(NX2_64B) && defined(NX2_32B)
#undef NX2_32B
#endif
#if defined(PVM_64B) && defined(PVM_32B)
#undef PVM_32B
#endif
#if defined(MPI_64B) && defined(MPI_32B)
#undef MPI_32B
#endif
#if defined(NX2_32B) || defined(PVM_32B) || defined(MPI_32B)
#define FLP_32B
#endif
#if defined(NX2_64B) || defined(PVM_64B) || defined(MPI_64B)
#define FLP_64B
#endif

C     Consistency checks. Beware: Check C and util source before changing !
#undef GC__INTERFACE_ID
#undef GC__INTERFACE
#undef GC__ERROR
#if defined(SHM_SRC)
#if defined(GC__INTERFACE)
#define GC__ERROR
#else
#define GC__INTERFACE_ID 1
#define GC__INTERFACE 'GC_SHM'
#endif
#endif
#if defined(NX2_SRC)
#if defined(GC__INTERFACE)
#define GC__ERROR
#else
#define GC__INTERFACE_ID 2
#define GC__INTERFACE 'GC_NX2'
#endif
#endif
#if defined(PVM_SRC)
#if defined(GC__INTERFACE)
#define GC__ERROR
#else
#define GC__INTERFACE_ID 3
#define GC__INTERFACE 'GC_PVM'
#endif
#endif
#if defined(MPI_SRC)
#if defined(GC__INTERFACE)
#define GC__ERROR
#else
#define GC__INTERFACE_ID 4
#define GC__INTERFACE 'GC_MPI'
#endif
#endif
#if defined(MPL_SRC)
#if defined(GC__INTERFACE)
#define GC__ERROR
#else
#define GC__INTERFACE_ID 5
#define GC__INTERFACE 'GC_MPL'
#endif
#endif
#if !defined(GC__INTERFACE)
#define GC__INTERFACE_ID 0
#define GC__INTERFACE 'GC_NONE'
#define GC__SERIAL
#endif
#if defined(GC__ERROR)
ERROR GC: multiple INTERFACES defined.
#endif


C     Size of sync work arrays - used internally only
#if !defined(GC__BCAST_SYNC)
#define GC__BCAST_SYNC  128
#endif
#if !defined(GC__COLL_SYNC)
#define GC__COLL_SYNC 36
#endif

C     The node id used for IO
#define GC__IONODE 0

C     Message tags reserved by GC
C     BEWARE: limits are duplicated in the header files !
#define GCID__TASK  999999999
#define GCID__SYNC0 999999998
#define GCID__SYNC1 999999997
#define GCID__BCAST 999999996
#define GCID__RSUM0 999999995
#define GCID__RSUM1 999999994
#define GCID__RMIN0 999999993
#define GCID__RMIN1 999999992
#define GCID__RMAX0 999999991
#define GCID__RMAX1 999999990
#define GCID__ISUM0 999999989
#define GCID__ISUM1 999999988
#define GCID__IMIN0 999999987
#define GCID__IMIN1 999999986
#define GCID__IMAX0 999999985
#define GCID__IMAX1 999999984
#define GCID__ARCH  999999983
#define GCID__FIRST GCID__ARCH
#define GCID__LAST  GCID__TASK

C     INFO tags reserved by GC, used in SHMEM only
#define GC__SHM_GET  -9999
#define GC__SHM_PUT  -9998

C     Derived quantities
#if defined(_CRAY) || defined(SHM_SRC)
#define GC__ISIZE   8
#else
#define GC__ISIZE   4
#endif

#if defined(FLP_64B)
#define GC__RSIZE   8
#define REAL real*8
#else
#define GC__RSIZE   4
#define REAL real*4
#endif
      
#if defined(NX2_SRC) && defined(FLP_64B)
#define GC__NX_GSUM GDSUM
#define GC__NX_GHIGH GDHIGH
#define GC__NX_GLOW GDLOW
#else
#define GC__NX_GSUM GSSUM
#define GC__NX_GHIGH GSHIGH
#define GC__NX_GLOW GSLOW
#endif

#if defined(MPL_SRC) && defined(FLP_64B)
#define GC__MPL_SUM D_VADD
#define GC__MPL_MIN D_VMIN
#define GC__MPL_MAX D_VMAX
#else
#define GC__MPL_SUM S_VADD
#define GC__MPL_MIN S_VMIN
#define GC__MPL_MAX S_VMAX
#endif

C     PVM byte swap support (not supported with CRI systems)

#if (!defined(PVM_SRC) || defined(_CRAY)) && !defined(PVM_NBS)
#define PVM_NBS
#endif
#if !defined(PVM_NBS)
#define GC__MYARCH 1
#endif


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
#if defined(SHM_SRC)
      INCLUDE 'mpp/shmem.fh'
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
      INTEGER MY_PE
      INTRINSIC MY_PE
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
      INTEGER MYNODE, NUMNODES
#endif
#if defined(PVM_SRC)
#if !defined(PVM_NBS) && (defined(__alpha) || defined(ultrix))
      INTEGER dec_rsize, dec_isize
      common /gc__dec/dec_rsize, dec_isize
#endif
#if defined(_CRAYMPP)
      INTEGER IME, ITID
#endif
      INCLUDE 'fpvm3.h'
      INTEGER MYTID, MYPTID, MYARC, INFO
      CHARACTER*256 ARGBUF, FILE, GC__SUFFIX
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER INFO
#endif
#if defined(MPL_SRC)
      INTEGER INFO, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif

      INTEGER I

      IF (INITED) RETURN

#if defined(SHM_SRC)
C
C     Initialize the synchronization workspace 
C     with SHMEM_SYNC_VALUE from  'mpp/shmem.h'.
C
      DO I = 1,GC__BCAST_SYNC
         BCAST_SYNC_WRK(I) = SHMEM_SYNC_VALUE
      ENDDO
      DO I = 1,GC__COLL_SYNC
         REDUCE_SYNC_WRK(I) = SHMEM_SYNC_VALUE
      ENDDO
      ME = MY_PE()
      NPROC = N$PES
#endif

#if defined(NX2_SRC)
      ME = MYNODE()
      NPROC = NUMNODES()
#endif

#if defined(PVM_SRC) 
#if defined(_CRAYMPP)
      CALL PVMFMYTID(MYTID)
      CALL PVMFGETPE(MYTID, ME)
      CALL PVMFGSIZE(PVMALL, NPROC)
      IF (ME .NE. 0) THEN
         CALL PVMFINITSEND(GC__PVMENCODING, INFO)
         CALL PVMFPACK(BYTE1, ME, GC__ISIZE, 1, INFO)
         CALL PVMFPACK(BYTE1, MYTID, GC__ISIZE, 1, INFO)
         CALL PVMFSEND(0, GCID__TASK, INFO)
         CALL PVMFRECV(0, GCID__BCAST, INFO)
         CALL PVMFUNPACK(BYTE1, TIDS, GC__ISIZE*(NPROC+1), 1, INFO)
      ELSE
         TIDS(-1) = -1
         TIDS(0) = MYTID
         DO I = 1,NPROC-1
            CALL PVMFRECV(-1, GCID__TASK, INFO)
            CALL PVMFUNPACK(BYTE1, IME, GC__ISIZE, 1, INFO)
            CALL PVMFUNPACK(BYTE1, ITID, GC__ISIZE, 1, INFO)
            TIDS(IME) = ITID
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING, INFO)
         CALL PVMFPACK(BYTE1, TIDS, GC__ISIZE*(NPROC+1), 1, INFO)
         CALL PVMFBCAST(PVMALL, GCID__BCAST, INFO)
      ENDIF
#else
C
C     Spawn NPROC-1 slaves, broadcast NPROC and PVM task IDs
C     If the user has supplied a non-blank path argument, use this
C     unconditionally. Otherwise use file name part of argv[0].
C
#if !defined(PVM_NBS) && (defined(__alpha) || defined(ultrix))
      dec_rsize = GC__RSIZE
      dec_isize = GC__ISIZE
#endif
      CALL PVMFMYTID(MYTID)
      CALL PVMFPARENT(MYPTID)
#if !defined(PVM_NBS) && !defined(_CRAY)
      MYARC = GC__MYARCH
#endif
      IF (MYPTID .EQ. PVMNOPARENT) THEN
         TIDS(-1) = -1
         TIDS(0) = MYTID
         IF (NPROC .GT. 1) THEN
            IF (PATH(1:1) .NE. ' ') THEN
               FILE = PATH
            ELSE
               CALL GETARG(0, ARGBUF)
               FILE = GC__SUFFIX(ARGBUF)
            ENDIF
            CALL PVMFSPAWN(FILE, GC__PVMENCODING, '*', NPROC-1, TIDS(1), I)
            IF (I .NE. NPROC-1) THEN
               NPROC = -1
               RETURN
            ENDIF
            CALL PVMFINITSEND(GC__PVMENCODING, I)
#if !defined(PVM_NBS) && !defined(_CRAY)
            CALL PVMFPACK(BYTE1, MYARC, GC__ISIZE, 1, I)
#endif
            CALL PVMFPACK(BYTE1, NPROC, GC__ISIZE, 1, I)
            CALL PVMFPACK(BYTE1, TIDS, GC__ISIZE*(NPROC+1), 1, I)
            CALL PVMFMCAST(NPROC-1, TIDS(1), GCID__TASK, I)
         ENDIF
      ELSE
         CALL PVMFRECV(MYPTID, GCID__TASK, I)
#if !defined(PVM_NBS) && !defined(_CRAY)
         CALL PVMFUNPACK(BYTE1, ARCS(0), GC__ISIZE, 1, I)
#endif
         CALL PVMFUNPACK(BYTE1, NPROC, GC__ISIZE, 1, I)
#if !defined(PVM_NBS) && !defined(_CRAY)
         IF (ARCS(0) .NE. MYARC) CALL GC__ISWAP(1, NPROC)
#endif
         CALL PVMFUNPACK(BYTE1, TIDS, GC__ISIZE*(NPROC+1), 1, I)
#if !defined(PVM_NBS) && !defined(_CRAY)
         IF (ARCS(0) .NE. MYARC) CALL GC__ISWAP(NPROC+1, TIDS)
#endif
      ENDIF
      DO I = 0, NPROC-1
         IF (MYTID .EQ. TIDS(I)) ME = I
      ENDDO
#if !defined(PVM_NBS) && !defined(_CRAY)
C     Exchange architectures
      ARCS(ME) = MYARC
      IF (ME .EQ. 0) THEN
         DO I = 1, NPROC-1
            CALL PVMFRECV(TIDS(I), GCID__ARCH, INFO)
            CALL PVMFUNPACK(BYTE1, ARCS(I), GC__ISIZE, 1, INFO)
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING, INFO)
         CALL PVMFPACK(BYTE1, ARCS, GC__ISIZE*NPROC, 1, INFO)
         CALL PVMFMCAST(NPROC-1, TIDS(1), GCID__ARCH, INFO)
      ELSE
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(0), GCID__ARCH, ARCS(ME), GC__ISIZE,
     $                  BYTE1, INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING, INFO)
         CALL PVMFPACK(BYTE1, ARCS(ME), GC__ISIZE, 1, INFO)
         CALL PVMFSEND(TIDS(0), GCID__ARCH, INFO)
#endif
         CALL PVMFRECV(TIDS(0), GCID__ARCH, INFO)
         CALL PVMFUNPACK(BYTE1, ARCS, GC__ISIZE*NPROC, 1, INFO)
      ENDIF
      DO I = 0, NPROC-1
         IF (ARCS(I) .EQ. MYARC) THEN
            ARCS(I) = 0
         ELSE
            ARCS(I) = -1
         ENDIF
      ENDDO
#endif
#endif
#endif

#if defined(MPI_SRC)
      CALL MPI_INIT(INFO)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, INFO)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, INFO)
#endif

#if defined(MPL_SRC)
      CALL MP_ENVIRON(NPROC,ME)
      CALL MP_TASK_QUERY (MPLBUF, 4, 3)
#endif

#if defined(GC__SERIAL)
      ME = 0
      NPROC = 1
#endif

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

#if defined(PVM_SRC) || defined(MPI_SRC)
      INTEGER INFO
#endif


#if defined(PVM_SRC)
      CALL PVMFEXIT(INFO)
#endif

#if defined(MPI_SRC)
      CALL MPI_FINALIZE(INFO)
#endif

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
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif

#if !defined(_CRAY) 
      CALL FLUSH(6)
#endif
      WRITE(*,*) 'gc_abort: ', MESG

#if defined(_CRAYMPP)
C     On CRAY MPP we always use EXIT().
      CALL EXIT(1)
#endif

#if defined(NX2_SRC)
      CALL KILL(0,9)
#endif

#if defined(PVM_SRC) 
      DO I = 0, NPROC-1
         IF (I .NE. ME) CALL PVMFKILL(TIDS(I),INFO)
      ENDDO
      CALL PVMFEXIT(INFO)
      CALL EXIT(1)
#endif

#if defined(MPI_SRC)
      CALL MPI_ABORT(MPI_COMM_WORLD,9,INFO)
#endif

#if defined(MPL_SRC)
      CALL MP_STOPALL(INFO)
#endif

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
C     *  INFO    - status of send 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI,INFO
      REAL RARR(LEN), SARR(LEN)
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      IF (INFO .EQ. GC__SHM_PUT) CALL SHMEM_PUT(RARR, SARR, LEN, RECI)
#endif

#if defined(NX2_SRC)
      CALL CSEND(MSG, SARR, GC__RSIZE*LEN, RECI, 0)
#endif

#if defined(PVM_SRC)
#if defined(PVM_V33)
      CALL PVMFPSEND(TIDS(RECI), MSG, SARR, GC__RSIZE*LEN, BYTE1, INFO)       
#else
      CALL PVMFINITSEND(GC__PVMENCODING,INFO)
      CALL PVMFPACK(BYTE1, SARR, GC__RSIZE*LEN, 1, INFO)
      CALL PVMFSEND(TIDS(RECI), MSG, INFO)
#endif
#endif

#if defined(MPI_SRC)
      CALL MPI_SEND(SARR, GC__RSIZE*LEN, MPI_BYTE, RECI, MSG,
     $     MPI_COMM_WORLD, INFO)
#endif

#if defined(MPL_SRC)
      CALL MP_BSEND(SARR,GC__RSIZE*LEN,RECI,MSG)
#endif

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
C     *  INFO    - status of send 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO
      REAL RARR(LEN), SARR(LEN)
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1), IBUF
      COMMON /GC__TIDS/ TIDS, ARCS
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAY)
      INTEGER ASEND, AMSG, ALEN, I
#endif
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      IF (INFO .NE. GC__SHM_PUT) CALL SHMEM_GET(RARR, SARR, LEN, SEND)
#endif

#if defined(NX2_SRC)
      CALL CRECV(MSG, RARR, GC__RSIZE*LEN)
#endif

#if defined(PVM_SRC)
#if defined(PVM_V33)
      CALL PVMFPRECV(TIDS(SEND), MSG, RARR, GC__RSIZE*LEN, BYTE1, 
     $     ASEND, AMSG, ALEN, INFO)
#else
      CALL PVMFRECV(TIDS(SEND), MSG, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAY)
      IF (SEND .EQ. -1) CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
      CALL PVMFUNPACK(BYTE1, RARR, GC__RSIZE*LEN, 1, INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAY)
C     We have to figure out the sender if using the -1 wildcard
      IF (SEND .EQ. -1) THEN
         I = -1
100      CONTINUE
            I = I + 1
         IF (I .LT. MAX_PROC .AND. ASEND .NE. TIDS(I)) GOTO 100
         IF (I .EQ. MAX_PROC) THEN
            INFO = -1
            RETURN
         ENDIF
         ASEND = I
      ELSE
         ASEND = SEND
      ENDIF
      IF (ARCS(ASEND) .NE. 0) CALL GC__RSWAP(LEN, RARR)
#endif
#endif

#if defined(MPI_SRC)
      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, GC__RSIZE*LEN, MPI_BYTE, SEND, MSG,
     $     MPI_COMM_WORLD, STATUS, INFO)
#endif

#if defined(MPL_SRC)
      IF (SEND .EQ. -1) SEND = MPLBUF(1)
      CALL MP_BRECV(RARR, GC__RSIZE*LEN, SEND, MSG, INFO)
#endif

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
C     *  INFO    - status of send. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI, INFO, RARR(LEN), SARR(LEN)
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      IF (INFO .EQ. GC__SHM_PUT) CALL SHMEM_PUT(RARR, SARR, LEN, RECI)
#endif

#if defined(NX2_SRC)
      CALL CSEND(MSG, SARR, GC__ISIZE*LEN, RECI, 0)
#endif

#if defined(PVM_SRC)
#if defined(PVM_V33)
      CALL PVMFPSEND(TIDS(RECI), MSG, SARR, GC__ISIZE*LEN, BYTE1, INFO)
#else
      CALL PVMFINITSEND(GC__PVMENCODING,INFO)
      CALL PVMFPACK(BYTE1, SARR, GC__ISIZE*LEN, 1, INFO)
      CALL PVMFSEND(TIDS(RECI), MSG, INFO)
#endif
#endif

#if defined(MPI_SRC)
      CALL MPI_SEND(SARR, GC__ISIZE*LEN, MPI_BYTE, RECI, MSG,
     $     MPI_COMM_WORLD, INFO)
#endif

#if defined(MPL_SRC)
      CALL MP_BSEND(SARR,GC__ISIZE*LEN,RECI,MSG)
#endif

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
C     *  INFO    - status of send 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO, RARR(LEN), SARR(LEN)
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1), IBUF
      COMMON /GC__TIDS/ TIDS, ARCS
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAY)
      INTEGER ASEND, AMSG, ALEN, I
#endif
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      IF (INFO .NE. GC__SHM_PUT) CALL SHMEM_GET(RARR, SARR, LEN, SEND)
#endif

#if defined(NX2_SRC)
      CALL CRECV(MSG, RARR, GC__ISIZE*LEN)
#endif

#if defined(PVM_SRC)
#if defined(PVM_V33)
      CALL PVMFPRECV(TIDS(SEND), MSG, RARR, GC__ISIZE*LEN, BYTE1, 
     $     ASEND, AMSG, ALEN, INFO)
#else
      CALL PVMFRECV(TIDS(SEND), MSG, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAY)
      IF (SEND .EQ. -1) CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
      CALL PVMFUNPACK(BYTE1, RARR, GC__ISIZE*LEN, 1, INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAY)
C     We have to figure out the sender if using the -1 wildcard
      IF (SEND .EQ. -1) THEN
         I = -1
100      CONTINUE
            I = I + 1
         IF (I .LT. MAX_PROC .AND. ASEND .NE. TIDS(I)) GOTO 100
         IF (I .EQ. MAX_PROC) THEN
            INFO = -1
            RETURN
         ENDIF
         ASEND = I
      ELSE
         ASEND = SEND
      ENDIF
      IF (ARCS(ASEND) .NE. 0) CALL GC__ISWAP(LEN, RARR)
#endif
#endif

#if defined(MPI_SRC)
      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, GC__ISIZE*LEN, MPI_BYTE, SEND, MSG,
     $     MPI_COMM_WORLD, STATUS, INFO)
#endif

#if defined(MPL_SRC)
      IF (SEND .EQ. -1) SEND = MPLBUF(1)
      CALL MP_BRECV(RARR, GC__ISIZE*LEN, SEND, MSG, INFO)
#endif

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
C     *  INFO    - status of send 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, RECI, INFO
     $        ,RARR(LEN/GC__ISIZE+1), SARR(LEN/GC__ISIZE+1)
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      IF (INFO .EQ. GC__SHM_PUT) 
     $    CALL SHMEM_PUT(RARR, SARR, LEN/GC__ISIZE+1, RECI)
#endif

#if defined(NX2_SRC)
      CALL CSEND(MSG, SARR, LEN, RECI, 0)
#endif

#if defined(PVM_SRC)
#if defined(PVM_V33)
      CALL PVMFPSEND(TIDS(RECI), MSG, SARR, LEN, BYTE1, INFO)
#else
      CALL PVMFINITSEND(GC__PVMENCODING,INFO)
      CALL PVMFPACK(BYTE1, SARR, LEN, 1, INFO)
      CALL PVMFSEND(TIDS(RECI), MSG, INFO)
#endif
#endif

#if defined(MPI_SRC)
      CALL MPI_SEND(SARR, LEN, MPI_BYTE, RECI, MSG, MPI_COMM_WORLD,
     $     INFO)
#endif

#if defined(MPL_SRC)
      CALL MP_BSEND(SARR,LEN,RECI,MSG)
#endif

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
C     *  INFO    - status of send 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, INFO, RARR(LEN), SARR(LEN/GC__ISIZE+1)
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#if defined(PVM_V33)
      INTEGER ASEND, AMSG, ALEN
#endif
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      IF (INFO .NE. GC__SHM_PUT) 
     $    CALL SHMEM_GET(RARR, SARR, LEN/GC__ISIZE+1, SEND)
#endif

#if defined(NX2_SRC)
      CALL CRECV(MSG, RARR, LEN)
#endif

#if defined(PVM_SRC)
#if defined(PVM_V33)
      CALL PVMFPRECV(TIDS(SEND), MSG, RARR, LEN, BYTE1, ASEND, 
     $     AMSG, ALEN, INFO)
#else
      CALL PVMFRECV(TIDS(SEND), MSG, INFO)
      CALL PVMFUNPACK(BYTE1, RARR, LEN, 1, INFO)
#endif
#endif

#if defined(MPI_SRC)
      IF (SEND .EQ. -1) SEND = MPI_ANY_SOURCE
      CALL MPI_RECV(RARR, LEN, MPI_BYTE, SEND, MSG, MPI_COMM_WORLD,
     $     STATUS, INFO)
#endif

#if defined(MPL_SRC)
      IF (SEND .EQ. -1) SEND = MPLBUF(1)
      CALL MP_BRECV(RARR, LEN, SEND, MSG, INFO)
#endif

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
C     *  INFO    - status of send 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     * 
C     *  No node can continue execution before everybody have reached
C     *  this point.
C     *  On CRAY MPP the cache is flushed if INFO == GC__SHM_PUT.
C     *  This feature should be used directly after a SHM_PUT.
C     *  
C     ******************************************************************

      IMPLICIT NONE
      INTEGER NPROC,INFO
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, ISYNC, MYTID
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif

#if defined(SHM_SRC)
      IF (INFO .EQ. GC__SHM_PUT) CALL SHMEM_UDCFLUSH()
      CALL BARRIER()
#endif

#if defined(NX2_SRC)
      CALL GSYNC()
#endif

#if defined(PVM_SRC)
#if defined(_CRAYMPP)
      CALL PVMFBARRIER(PVMALL, -1, INFO)
#else
      ISYNC = 1
      CALL PVMFMYTID(MYTID)
      IF (MYTID .EQ. TIDS(GC__IONODE)) THEN
         DO I = 1, NPROC-1
            CALL PVMFRECV(-1, GCID__SYNC0, INFO)
            CALL PVMFUNPACK(BYTE1, ISYNC, GC__ISIZE, 1, INFO)
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING, INFO)
         CALL PVMFPACK(BYTE1, ISYNC, GC__ISIZE, 1, INFO)
         CALL PVMFMCAST(NPROC-1, TIDS(1), GCID__SYNC1, INFO)
      ELSE
         CALL PVMFINITSEND(GC__PVMENCODING, INFO)
         CALL PVMFPACK(BYTE1, ISYNC, GC__ISIZE, 1, INFO)
         CALL PVMFSEND(TIDS(GC__IONODE), GCID__SYNC0, INFO)
         CALL PVMFRECV(TIDS(GC__IONODE), GCID__SYNC1, INFO)
         CALL PVMFUNPACK(BYTE1, ISYNC, GC__ISIZE, 1, INFO)
      ENDIF
#endif
#endif

#if defined(MPI_SRC)
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
#endif

#if defined(MPL_SRC)
      CALL MP_SYNC(MPLBUF(4))
#endif

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
C     *  The cache is flushed if INFO == GC__SHM_PUT.
C     *  This feature should be used directly after a SHM_PUT.
C     *  
C     ******************************************************************

      IMPLICIT NONE
      INTEGER NPROC, INFO

#if defined(SHM_SRC)
      IF (INFO .EQ. GC__SHM_PUT) CALL SHMEM_UDCFLUSH()
      CALL BARRIER()
#endif

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
C     *  INFO    - status of bcast. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO
      REAL SARR(LEN)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
      INTEGER MYNODE
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER MYTID
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#if defined(PVM_V33)
      INTEGER ASEND, AMSG, ALEN
#endif
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_BROADCAST(SARR,SARR,LEN,SEND,0,0,NPROC,
     $     BCAST_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      IF (SEND .EQ. MYNODE()) THEN
         CALL CSEND(MSG, SARR, GC__RSIZE*LEN, -1, 0)
      ELSE
         CALL CRECV(MSG, SARR, GC__RSIZE*LEN)
      ENDIF
#endif

#if defined(PVM_SRC)
      CALL PVMFMYTID(MYTID)
      IF (TIDS(SEND) .EQ. MYTID) THEN
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1, SARR, GC__RSIZE*LEN, 1, INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL, MSG, INFO)
#else
         CALL PVMFMCAST(NPROC, TIDS(0), MSG, INFO)
#endif
      ELSE
#if defined(PVM_V33)
         CALL PVMFPRECV(TIDS(SEND), MSG, SARR, GC__RSIZE*LEN, BYTE1, 
     $        ASEND, AMSG, ALEN, INFO)
#else
         CALL PVMFRECV(TIDS(SEND), MSG, INFO)
         CALL PVMFUNPACK(BYTE1, SARR, GC__RSIZE*LEN, 1, INFO)
#if !defined(PVM_NBS) && !defined(_CRAY)
         IF (ARCS(SEND) .NE. 0) CALL GC__RSWAP(LEN, SARR)
#endif
#endif
      ENDIF
#endif

#if defined (MPI_SRC)
      CALL MPI_BCAST(SARR, GC__RSIZE*LEN, MPI_BYTE, SEND,
     $     MPI_COMM_WORLD, INFO)
#endif

#if defined (MPL_SRC)
      CALL MP_BCAST(SARR, GC__RSIZE*LEN, SEND, MPLBUF(4))
#endif

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
C     *  INFO    - status of bcast. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO, SARR(LEN)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
      INTEGER MYNODE
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER MYTID
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#if defined(PVM_V33)
      INTEGER ASEND, AMSG, ALEN
#endif
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_BROADCAST(SARR,SARR,LEN,SEND,0,0,NPROC,
     $     BCAST_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      IF (SEND .EQ. MYNODE()) THEN
         CALL CSEND(MSG, SARR, GC__ISIZE*LEN, -1, 0)
      ELSE
         CALL CRECV(MSG, SARR, GC__ISIZE*LEN)
      ENDIF
#endif

#if defined(PVM_SRC)
      CALL PVMFMYTID(MYTID)
      IF (TIDS(SEND) .EQ. MYTID) THEN
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1, SARR, GC__ISIZE*LEN, 1, INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL, MSG, INFO)
#else
         CALL PVMFMCAST(NPROC, TIDS(0), MSG, INFO)
#endif
      ELSE
#if defined(PVM_V33)
         CALL PVMFPRECV(TIDS(SEND), MSG, SARR, GC__ISIZE*LEN, BYTE1, 
     $        ASEND, AMSG, ALEN, INFO)
#else
         CALL PVMFRECV(TIDS(SEND), MSG, INFO)
         CALL PVMFUNPACK(BYTE1, SARR, GC__ISIZE*LEN, 1, INFO)
#if !defined(PVM_NBS) && !defined(_CRAY)
         IF (ARCS(SEND) .NE. 0) CALL GC__ISWAP(LEN, SARR)
#endif
#endif
      ENDIF
#endif

#if defined (MPI_SRC)
      CALL MPI_BCAST(SARR, GC__ISIZE*LEN, MPI_BYTE, SEND,
     $     MPI_COMM_WORLD, INFO)
#endif

#if defined (MPL_SRC)
      CALL MP_BCAST(SARR, GC__ISIZE*LEN, SEND, MPLBUF(4))
#endif

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
C     *  INFO    - status of bcast. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER MSG, LEN, SEND, NPROC, INFO, SARR(LEN/GC__ISIZE+1)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
      INTEGER MYNODE
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER MYTID
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#if defined(PVM_V33)
      INTEGER ASEND, AMSG, ALEN
#endif
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
#endif
#if defined(MPL_SRC)
      INTEGER MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_BROADCAST(SARR,SARR,LEN/GC__ISIZE+1,SEND,0,0,NPROC,
     $     BCAST_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      IF (SEND .EQ. MYNODE()) THEN
         CALL CSEND(MSG, SARR, LEN, -1, 0)
      ELSE
         CALL CRECV(MSG, SARR, LEN)
      ENDIF
#endif

#if defined(PVM_SRC)
      CALL PVMFMYTID(MYTID)
      IF (TIDS(SEND) .EQ. MYTID) THEN
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1, SARR, LEN, 1, INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL, MSG, INFO)
#else
         CALL PVMFMCAST(NPROC, TIDS(0), MSG, INFO)
#endif
      ELSE
#if defined(PVM_V33)
         CALL PVMFPRECV(TIDS(SEND), MSG, SARR, LEN, BYTE1, ASEND, 
     $        AMSG, ALEN, INFO)
#else
         CALL PVMFRECV(TIDS(SEND), MSG, INFO)
         CALL PVMFUNPACK(BYTE1, SARR, LEN, 1, INFO)
#endif
      ENDIF
#endif

#if defined (MPI_SRC)
      CALL MPI_BCAST(SARR, LEN, MPI_BYTE, SEND, MPI_COMM_WORLD, INFO)
#endif

#if defined (MPL_SRC)
      CALL MP_BCAST(SARR, LEN, SEND, MPLBUF(4))
#endif

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
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      REAL SSUM(LEN)
      REAL REDUCE_DATA_WRK(MAX_COLL)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, L, MYTID, IBUF
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAY)
      INTEGER ASEND, AMSG, ALEN, J
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER I
#endif
#if defined(MPL_SRC)
      INTEGER I, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
      EXTERNAL GC__MPL_SUM
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_REAL8_SUM_TO_ALL(SSUM,SSUM,LEN,0,0,NPROC,
     &     REDUCE_DATA_WRK,REDUCE_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      CALL GC__NX_GSUM(SSUM,LEN,REDUCE_DATA_WRK)
#endif

#if defined(PVM_SRC)
C
C      The PVM source should use a binary tree, but that is
C      complicated to code in case of non-power of two processors.
C
      CALL PVMFMYTID(MYTID)
      IF (MYTID .NE. TIDS(GC__IONODE)) THEN
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(GC__IONODE),GCID__RSUM0,SSUM,GC__RSIZE*LEN,
     $        BYTE1,INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,SSUM,GC__RSIZE*LEN,1,INFO)
         CALL PVMFSEND(TIDS(GC__IONODE),GCID__RSUM0,INFO)
#endif
         CALL PVMFRECV(TIDS(GC__IONODE),GCID__RSUM1,INFO)
         CALL PVMFUNPACK(BYTE1,SSUM,GC__RSIZE*LEN,1,INFO)
#if !defined(PVM_NBS) && !defined(_CRAY)
         IF (ARCS(GC__IONODE) .NE. 0) CALL GC__RSWAP(LEN, SSUM)
#endif
      ELSE
         DO I = 1,NPROC-1
#if defined(PVM_V33)
            CALL PVMFPRECV(-1,GCID__RSUM0,REDUCE_DATA_WRK,GC__RSIZE*LEN,
     $           BYTE1,ASEND,AMSG,ALEN,INFO)
#else
            CALL PVMFRECV(-1,GCID__RSUM0, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
            CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
            CALL PVMFUNPACK(BYTE1,REDUCE_DATA_WRK,GC__RSIZE*LEN,1,INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
C           We could have specified a receive tid of TIDS(I) and avoided
C           this, but taking into account the high comp to comm ratio of
C           PVM this approach is probably faster.
            J = 0
100         CONTINUE
               J = J + 1
            IF (J .LT. NPROC .AND. ASEND .NE. TIDS(J)) GOTO 100
            IF (J .EQ. NPROC) THEN
               INFO = -1
               RETURN
            ENDIF
            ASEND = J
            IF (ARCS(ASEND) .NE. 0) CALL GC__RSWAP(LEN,
     $         REDUCE_DATA_WRK)
#endif
            DO L=1,LEN
               SSUM(L) = SSUM(L) + REDUCE_DATA_WRK(L)
            ENDDO
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,SSUM,GC__RSIZE*LEN,1,INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL,GCID__RSUM1,INFO)
#else
         CALL PVMFMCAST(NPROC-1,TIDS(1),GCID__RSUM1,INFO)
#endif
      ENDIF
#endif

#if defined(MPI_SRC)
      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SSUM(I)
      ENDDO
#if defined(FLP_64B)
      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SSUM, LEN, 
     $     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
#else
      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SSUM, LEN, MPI_REAL, MPI_SUM,
     $     MPI_COMM_WORLD, INFO)
#endif
#endif

#if defined(MPL_SRC)
      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SSUM(I)
      END DO
      CALL MP_COMBINE(REDUCE_DATA_WRK,SSUM,GC__RSIZE*LEN,GC__MPL_SUM,
     $                MPLBUF(4))
#endif
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
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      REAL SMIN(LEN)
      REAL REDUCE_DATA_WRK(MAX_COLL)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, L, MYTID, IBUF
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAYMPP)
      INTEGER ASEND, AMSG, ALEN, J
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER I
#endif
#if defined(MPL_SRC)
      INTEGER I, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
      EXTERNAL GC__MPL_MIN
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_REAL8_MIN_TO_ALL(SMIN,SMIN,LEN,0,0,NPROC,
     &     REDUCE_DATA_WRK,REDUCE_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      CALL GC__NX_GLOW(SMIN,LEN,REDUCE_DATA_WRK)
#endif

#if defined(PVM_SRC)
C
C      The PVM source should use a binary tree, but that is
C      complicated to code in case of non-power of two processors.
C
      CALL PVMFMYTID(MYTID)
      IF (MYTID .NE. TIDS(GC__IONODE)) THEN
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(GC__IONODE),GCID__RMIN0,SMIN,GC__RSIZE*LEN,
     $        BYTE1,INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,SMIN,GC__RSIZE*LEN,1,INFO)
         CALL PVMFSEND(TIDS(GC__IONODE),GCID__RMIN0,INFO)
#endif
         CALL PVMFRECV(TIDS(GC__IONODE),GCID__RMIN1,INFO)
         CALL PVMFUNPACK(BYTE1,SMIN,GC__RSIZE*LEN,1,INFO)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
         IF (ARCS(GC__IONODE) .NE. 0) CALL GC__RSWAP(LEN, SMIN)
#endif
      ELSE
         DO I = 1,NPROC-1
#if defined(PVM_V33)
            CALL PVMFPRECV(-1,GCID__RMIN0,REDUCE_DATA_WRK,GC__RSIZE*LEN,
     $           BYTE1,ASEND,AMSG,ALEN,INFO)
#else
            CALL PVMFRECV(-1, GCID__RMIN0, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
            CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
            CALL PVMFUNPACK(BYTE1,REDUCE_DATA_WRK,GC__RSIZE*LEN,1,INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
C           We could have specified a receive tid of TIDS(I) and avoided
C           this, but taking into account the high comp to comm ratio of
C           PVM this approach is probably faster.
            J = 0
100         CONTINUE
               J = J + 1
            IF (J .LT. NPROC .AND. ASEND .NE. TIDS(J)) GOTO 100
            IF (J .EQ. NPROC) THEN
               INFO = -1
               RETURN
            ENDIF
            ASEND = J
            IF (ARCS(ASEND) .NE. 0) CALL GC__RSWAP(LEN,
     $         REDUCE_DATA_WRK)
#endif
            DO L=1,LEN
               SMIN(L) = MIN(SMIN(L),REDUCE_DATA_WRK(L))
            ENDDO
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,SMIN,GC__RSIZE*LEN,1,INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL,GCID__RMIN1,INFO)
#else
         CALL PVMFMCAST(NPROC-1,TIDS(1),GCID__RMIN1,INFO)
#endif
      ENDIF
#endif

#if defined(MPI_SRC)
	DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMIN(I)
      ENDDO
#if defined(FLP_64B)
      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMIN, LEN, 
     $     MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO)
#else
      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMIN, LEN, MPI_REAL, MPI_MIN,
     $     MPI_COMM_WORLD, INFO)
#endif
#endif

#if defined(MPL_SRC)
      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMIN(I)
      END DO
      CALL MP_COMBINE(REDUCE_DATA_WRK,SMIN,GC__RSIZE*LEN,GC__MPL_MIN,
     $                MPLBUF(4))
#endif

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
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO
      REAL SMAX(LEN)
      REAL REDUCE_DATA_WRK(MAX_COLL)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, L, MYTID, IBUF
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAYMPP)
      INTEGER ASEND, AMSG, ALEN, J
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER I
#endif
#if defined(MPL_SRC)
      INTEGER I, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
      EXTERNAL GC__MPL_MAX
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_REAL8_MAX_TO_ALL(SMAX,SMAX,LEN,0,0,NPROC,
     &     REDUCE_DATA_WRK,REDUCE_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      CALL GC__NX_GHIGH(SMAX,LEN,REDUCE_DATA_WRK)
#endif

#if defined(PVM_SRC)
C
C      The PVM source should use a binary tree, but that is
C      complicated to code in case of non-power of two processors.
C
      CALL PVMFMYTID(MYTID)
      IF (MYTID .NE. TIDS(GC__IONODE)) THEN
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(GC__IONODE),GCID__RMAX0,SMAX,GC__RSIZE*LEN,
     $        BYTE1,INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,SMAX,GC__RSIZE*LEN,1,INFO)
         CALL PVMFSEND(TIDS(GC__IONODE),GCID__RMAX0,INFO)
#endif
         CALL PVMFRECV(TIDS(GC__IONODE),GCID__RMAX1,INFO)
         CALL PVMFUNPACK(BYTE1,SMAX,GC__RSIZE*LEN,1,INFO)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
         IF (ARCS(GC__IONODE) .NE. 0) CALL GC__RSWAP(LEN, SMAX)
#endif
      ELSE
         DO I = 1,NPROC-1
#if defined(PVM_V33)
            CALL PVMFPRECV(-1,GCID__RMAX0,REDUCE_DATA_WRK,GC__RSIZE*LEN,
     $           BYTE1,ASEND,AMSG,ALEN,INFO)
#else
            CALL PVMFRECV(-1, GCID__RMAX0, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
            CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
            CALL PVMFUNPACK(BYTE1,REDUCE_DATA_WRK,GC__RSIZE*LEN,1,INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
C           We could have specified a receive tid of TIDS(I) and avoided
C           this, but taking into account the high comp to comm ratio of
C           PVM this approach is probably faster.
            J = 0
100         CONTINUE
               J = J + 1
            IF (J .LT. NPROC .AND. ASEND .NE. TIDS(J)) GOTO 100
            IF (J .EQ. NPROC) THEN
               INFO = -1
               RETURN
            ENDIF
            ASEND = J
            IF (ARCS(ASEND) .NE. 0) CALL GC__RSWAP(LEN,
     $         REDUCE_DATA_WRK)
#endif
            DO L=1,LEN
               SMAX(L) = MAX(SMAX(L),REDUCE_DATA_WRK(L))
            ENDDO
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,SMAX,GC__RSIZE*LEN,1,INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL,GCID__RMAX1,INFO)
#else
         CALL PVMFMCAST(NPROC-1,TIDS(1),GCID__RMAX1,INFO)
#endif
      ENDIF
#endif

#if defined(MPI_SRC)
      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMAX(I)
      ENDDO
#if defined(FLP_64B)
      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMAX, LEN, 
     $     MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, INFO)
#else
      CALL MPI_ALLREDUCE(REDUCE_DATA_WRK, SMAX, LEN, MPI_REAL, MPI_MAX,
     $     MPI_COMM_WORLD, INFO)
#endif
#endif

#if defined(MPL_SRC)
      DO I = 1,LEN
         REDUCE_DATA_WRK(I) = SMAX(I)
      ENDDO
      CALL MP_COMBINE(REDUCE_DATA_WRK,SMAX,GC__RSIZE*LEN,GC__MPL_MAX,
     $                MPLBUF(4))
#endif

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
C     *  INFO    - status of isum. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, ISUM(LEN)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, L, MYTID, IBUF
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAYMPP)
      INTEGER ASEND, AMSG, ALEN, J
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER I
#endif
#if defined(MPL_SRC)
      INTEGER I, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
      EXTERNAL I_VADD
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_INT8_SUM_TO_ALL(ISUM,ISUM,LEN,0,0,NPROC,
     &     REDUCE_DATA_IWRK,REDUCE_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      CALL GISUM(ISUM,LEN,REDUCE_DATA_IWRK)
#endif

#if defined(PVM_SRC)
C
C      The PVM source should use a binary tree, but that is
C      complicated to code in case of non-power of two processors.
C
      CALL PVMFMYTID(MYTID)
      IF (MYTID .NE. TIDS(GC__IONODE)) THEN
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(GC__IONODE),GCID__ISUM0,ISUM,GC__ISIZE*LEN,
     $        BYTE1,INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,ISUM,GC__ISIZE*LEN,1,INFO)
         CALL PVMFSEND(TIDS(GC__IONODE),GCID__ISUM0,INFO)
#endif
         CALL PVMFRECV(TIDS(GC__IONODE),GCID__ISUM1,INFO)
         CALL PVMFUNPACK(BYTE1,ISUM,GC__ISIZE*LEN,1,INFO)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
         IF (ARCS(GC__IONODE) .NE. 0) CALL GC__ISWAP(LEN, ISUM)
#endif
      ELSE
         DO I = 1,NPROC-1
#if defined(PVM_V33)
            CALL PVMFPRECV(-1,GCID__ISUM0,REDUCE_DATA_IWRK,GC__ISIZE*LEN
     $          ,BYTE1,ASEND,AMSG,ALEN,INFO)
#else
            CALL PVMFRECV(-1, GCID__ISUM0, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
            CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
            CALL PVMFUNPACK(BYTE1,REDUCE_DATA_IWRK,GC__ISIZE*LEN,1,INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
C           We could have specified a receive tid of TIDS(I) and avoided
C           this, but taking into account the high comp to comm ratio of
C           PVM this approach is probably faster.
            J = 0
100         CONTINUE
               J = J + 1
            IF (J .LT. NPROC .AND. ASEND .NE. TIDS(J)) GOTO 100
            IF (J .EQ. NPROC) THEN
               INFO = -1
               RETURN
            ENDIF
            ASEND = J
            IF (ARCS(ASEND) .NE. 0) CALL GC__ISWAP(LEN,
     $         REDUCE_DATA_IWRK)
#endif
            DO L=1,LEN
               ISUM(L) = ISUM(L) + REDUCE_DATA_IWRK(L)
            ENDDO
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,ISUM,GC__ISIZE*LEN,1,INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL,GCID__ISUM1,INFO)
#else
         CALL PVMFMCAST(NPROC-1,TIDS(1),GCID__ISUM1,INFO)
#endif
      ENDIF
#endif

#if defined(MPI_SRC)
      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = ISUM(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, ISUM, LEN, MPI_INTEGER,
     $     MPI_SUM, MPI_COMM_WORLD, INFO)
#endif

#if defined(MPL_SRC)
      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = ISUM(I)
      ENDDO
      CALL MP_COMBINE(REDUCE_DATA_IWRK,ISUM,GC__ISIZE*LEN,I_VADD,
     $                MPLBUF(4))
#endif
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
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, IMIN(LEN)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, L, MYTID, IBUF
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAYMPP)
      INTEGER ASEND, AMSG, ALEN, J
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER I
#endif
#if defined(MPL_SRC)
      INTEGER I, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
      EXTERNAL I_VMIN
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_INT8_MIN_TO_ALL(IMIN,IMIN,LEN,0,0,NPROC,
     &     REDUCE_DATA_IWRK,REDUCE_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      CALL GILOW(IMIN,LEN,REDUCE_DATA_IWRK)
#endif

#if defined(PVM_SRC)
C
C      The PVM source should use a binary tree, but that is
C      complicated to code in case of non-power of two processors.
C
      CALL PVMFMYTID(MYTID)
      IF (MYTID .NE. TIDS(GC__IONODE)) THEN
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(GC__IONODE),GCID__IMIN0,IMIN,GC__ISIZE*LEN,
     $        BYTE1,INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,IMIN,GC__ISIZE*LEN,1,INFO)
         CALL PVMFSEND(TIDS(GC__IONODE),GCID__IMIN0,INFO)
#endif
         CALL PVMFRECV(TIDS(GC__IONODE),GCID__IMIN1,INFO)
         CALL PVMFUNPACK(BYTE1,IMIN,GC__ISIZE*LEN,1,INFO)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
         IF (ARCS(GC__IONODE) .NE. 0) CALL GC__ISWAP(LEN, IMIN)
#endif
      ELSE
         DO I = 1,NPROC-1
#if defined(PVM_V33)
            CALL PVMFPRECV(-1,GCID__IMIN0,REDUCE_DATA_IWRK,GC__ISIZE*LEN
     $          ,BYTE1,ASEND,AMSG,ALEN,INFO)
#else
            CALL PVMFRECV(-1,GCID__IMIN0,IBUF)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
            CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, IBUF)
#endif
            CALL PVMFUNPACK(BYTE1,REDUCE_DATA_IWRK,GC__ISIZE*LEN,1,INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
C           We could have specified a receive tid of TIDS(I) and avoided
C           this, but taking into account the high comp to comm ratio of
C           PVM this approach is probably faster.
            J = 0
100         CONTINUE
               J = J + 1
            IF (J .LT. NPROC .AND. ASEND .NE. TIDS(J)) GOTO 100
            IF (J .EQ. NPROC) THEN
               INFO = -1
               RETURN
            ENDIF
            ASEND = J
            IF (ARCS(ASEND) .NE. 0) CALL GC__ISWAP(LEN,
     $         REDUCE_DATA_IWRK)
#endif
            DO L=1,LEN
               IMIN(L) = MIN(IMIN(L),REDUCE_DATA_IWRK(L))
            ENDDO
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,IMIN,GC__ISIZE*LEN,1,INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL,GCID__IMIN1,INFO)
#else
         CALL PVMFMCAST(NPROC-1,TIDS(1),GCID__IMIN1,INFO)
#endif
      ENDIF
#endif

#if defined(MPI_SRC)
      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMIN(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, IMIN, LEN, MPI_INTEGER,
     $     MPI_MIN, MPI_COMM_WORLD, INFO)
#endif

#if defined(MPL_SRC)
      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMIN(I)
      END DO
      CALL MP_COMBINE(REDUCE_DATA_IWRK,IMIN,GC__ISIZE*LEN,I_VMIN,
     $                MPLBUF(4))
#endif

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
C     *  INFO    - status of rsum. 0 is OK (PVM_SRC and MPI_SRC only),
C     *            refer to the header files for nonzero status codes
C     *
C     * NOTES:       
C     *
C     ******************************************************************

      IMPLICIT NONE
      INTEGER LEN, NPROC, INFO, IMAX(LEN)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
#if defined(SHM_SRC)
      INTEGER BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(NX2_SRC)
      INCLUDE 'fnx.h'
#endif
#if defined(PVM_SRC)
      INCLUDE 'fpvm3.h'
      INTEGER I, L, MYTID, IBUF
      INTEGER TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif
#if (defined(PVM_V33) || !defined(PVM_NBS)) && !defined(_CRAYMPP)
      INTEGER ASEND, AMSG, ALEN, J
#endif
#if defined(MPI_SRC)
      INCLUDE 'mpif.h'
      INTEGER I
#endif
#if defined(MPL_SRC)
      INTEGER I, MPLBUF(4)
      COMMON /MPLCOM/MPLBUF
      EXTERNAL I_VMAX
#endif


#if defined(SHM_SRC)
      CALL BARRIER()
      CALL SHMEM_INT8_MAX_TO_ALL(IMAX,IMAX,LEN,0,0,NPROC,
     &     REDUCE_DATA_IWRK,REDUCE_SYNC_WRK)
#endif

#if defined(NX2_SRC)
      CALL GIHIGH(IMAX,LEN,REDUCE_DATA_IWRK)
#endif

#if defined(PVM_SRC)
C
C      The PVM source should use a binary tree, but that is
C      complicated to code in case of non-power of two processors.
C
      CALL PVMFMYTID(MYTID)
      IF (MYTID .NE. TIDS(GC__IONODE)) THEN
#if defined(PVM_V33)
         CALL PVMFPSEND(TIDS(GC__IONODE),GCID__IMAX0,IMAX,GC__ISIZE*LEN,
     $        BYTE1,INFO)
#else
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,IMAX,GC__ISIZE*LEN,1,INFO)
         CALL PVMFSEND(TIDS(GC__IONODE),GCID__IMAX0,INFO)
#endif
         CALL PVMFRECV(TIDS(GC__IONODE),GCID__IMAX1,INFO)
         CALL PVMFUNPACK(BYTE1,IMAX,GC__ISIZE*LEN,1,INFO)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
         IF (ARCS(GC__IONODE) .NE. 0) CALL GC__ISWAP(LEN, IMAX)
#endif
      ELSE
         DO I = 1,NPROC-1
#if defined(PVM_V33)
            CALL PVMFPRECV(-1,GCID__IMAX0,REDUCE_DATA_IWRK,GC__ISIZE*LEN
     $          ,BYTE1,ASEND,AMSG,ALEN,INFO)
#else
            CALL PVMFRECV(-1, GCID__IMAX0, IBUF)
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
            CALL PVMFBUFINFO(IBUF, ALEN, AMSG, ASEND, INFO)
#endif
            CALL PVMFUNPACK(BYTE1,REDUCE_DATA_IWRK,GC__ISIZE*LEN,1,INFO)
#endif
#if !defined(PVM_NBS) && !defined(_CRAYMPP)
C           We could have specified a receive tid of TIDS(I) and avoided
C           this, but taking into account the high comp to comm ratio of
C           PVM this approach is probably faster.
            J = 0
100         CONTINUE
               J = J + 1
            IF (J .LT. NPROC .AND. ASEND .NE. TIDS(J)) GOTO 100
            IF (J .EQ. NPROC) THEN
               INFO = -1
               RETURN
            ENDIF
            ASEND = J
            IF (ARCS(ASEND) .NE. 0) CALL GC__ISWAP(LEN,
     $         REDUCE_DATA_IWRK)
#endif
            DO L=1,LEN
               IMAX(L) = MAX(IMAX(L),REDUCE_DATA_IWRK(L))
            ENDDO
         ENDDO
         CALL PVMFINITSEND(GC__PVMENCODING,INFO)
         CALL PVMFPACK(BYTE1,IMAX,GC__ISIZE*LEN,1,INFO)
#if defined(_CRAYMPP)
         CALL PVMFBCAST(PVMALL,GCID__IMAX1,INFO)
#else
         CALL PVMFMCAST(NPROC-1,TIDS(1),GCID__IMAX1,INFO)
#endif
      ENDIF
#endif

#if defined(MPI_SRC)
      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMAX(I)
      ENDDO
      CALL MPI_ALLREDUCE(REDUCE_DATA_IWRK, IMAX, LEN, MPI_INTEGER,
     $     MPI_MAX, MPI_COMM_WORLD, INFO)
#endif

#if defined(MPL_SRC)
      DO I = 1,LEN
         REDUCE_DATA_IWRK(I) = IMAX(I)
      END DO
      CALL MP_COMBINE(REDUCE_DATA_IWRK,IMAX,GC__ISIZE*LEN,I_VMAX,
     $                MPLBUF(4))
#endif

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

#if defined(FLP_64B)
      GC_RSIZE = 8
#else
      GC_RSIZE = 4
#endif
    
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

#if defined(_CRAY) || defined (SHM_SRC)
      GC_ISIZE = 8
#else
      GC_ISIZE = 4
#endif
    
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
#if defined(_CRAYMPP)
      INTEGER LASTADDR, FIRSTADDR
#endif
      INTEGER LOC, NBYTES


#if defined(_CRAY) && !defined(_CRAYMPP)
#if defined(_ADDR32) || defined(_CRAY2)
C     32-bit adressable CRI PVP systems.
      NBYTES = 8 * (AND(LOC(LAST),X'FFFFFFFF') -
     $              AND(LOC(FIRST),X'FFFFFFFF'))
#else
C     24-bit adressable CRI PVP systems.
      NBYTES = 8 * (AND(LOC(LAST),X'FFFFFF') -
     $              AND(LOC(FIRST),X'FFFFFF'))
#endif
#else
#if defined(_CRAYMPP)
C     CRI MPP systems
#define GC__MINADDR X'FFFFFF'
      LASTADDR = LOC(LAST)
      FIRSTADDR = LOC(FIRST)
  
C     Check if FIRST holds a CHARACTER length and not an address
      IF (FIRSTADDR .LE. GC__MINADDR) THEN
         GC_COMLEN = 0
         RETURN
      ENDIF
      NBYTES = LASTADDR - FIRSTADDR
#else
C     Systems with UNIX Fortran character argument passing
      NBYTES = LOC(LAST) - LOC(FIRST)
#endif
#endif

      IF (LTYPE(1:1) .EQ. 'r' .OR. LTYPE(1:1) .EQ. 'R') THEN
         NBYTES = NBYTES + GC__RSIZE*LLEN
      ELSE IF (LTYPE(1:1) .EQ. 'i' .OR. LTYPE(1:1) .EQ. 'I') THEN
         NBYTES = NBYTES + GC__ISIZE*LLEN
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


      MXPROC = MAX_PROC
      MXCOLL = MAX_COLL
      MXPT2PT = MAX_PT2PT
      INTF = GC__INTERFACE
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
      INTFID = GC__INTERFACE_ID
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
#if defined(PVM_SRC)
      INTEGER I, TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif


      IF (.NOT. INITED) THEN
         GC__FPVMCONFIG = -1
         RETURN
      ENDIF
#if defined(PVM_SRC)
      DO I = 0, NPROC-1
         FTIDS(I) = TIDS(I)
         FTIDS(I+NPROC) = ARCS(I)
      ENDDO
#endif

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
#if defined(SHM_SRC)
      INCLUDE 'mpp/shmem.fh'
      INTEGER I, BCAST_SYNC_WRK(GC__BCAST_SYNC)
      INTEGER REDUCE_SYNC_WRK(GC__COLL_SYNC)    
      REAL REDUCE_DATA_WRK(MAX_COLL)
      INTEGER REDUCE_DATA_IWRK(MAX_COLL)
      EQUIVALENCE (REDUCE_DATA_WRK,REDUCE_DATA_IWRK)
      COMMON /GC__SHM/ BCAST_SYNC_WRK, REDUCE_SYNC_WRK, REDUCE_DATA_WRK   
#endif
#if defined(PVM_SRC)
      INTEGER I, TIDS(-1:MAX_PROC-1), ARCS(0:MAX_PROC-1)
      COMMON /GC__TIDS/ TIDS, ARCS
#endif


      IF (INITED) THEN
         GC__C2F = 0
         RETURN
      ELSE IF (NPROC .GT. MAX_PROC
     $        .OR. INTFID .NE. GC__INTERFACE_ID) THEN
         GC__C2F = -1
         RETURN
      ELSE
         GC__C2F = 0
      ENDIF
      GC__NPROC = NPROC

#if defined(SHM_SRC)
C
C     Initialize the synchronization workspace 
C     with SHMEM_SYNC_VALUE from  'mpp/shmem.h'.
C
      DO I = 1,GC__BCAST_SYNC
         BCAST_SYNC_WRK(I) = SHMEM_SYNC_VALUE
      ENDDO
      DO I = 1,GC__COLL_SYNC
         REDUCE_SYNC_WRK(I) = SHMEM_SYNC_VALUE
      ENDDO
#endif

#if defined(PVM_SRC)
      TIDS(-1) = -1
      DO I = 0, NPROC-1
         TIDS(I) = CTIDS(I)
      ENDDO
#endif

      INITED = .TRUE.
      END


#if !defined(PVM_NBS)

#if defined(__alpha) || defined(ultrix)
C     Hack to avoid faulty compiler optiomization on DEC systems
#define gcrsize dec_rsize
#define gcisize dec_isize
#else
#define gcrsize GC__RSIZE
#define gcisize GC__ISIZE
#endif

      SUBROUTINE GC__RSWAP(LEN, ARR)
C
C     Support routine for byte swapping GC real arrays
C
      IMPLICIT NONE
      INTEGER LEN
      real*GC__RSIZE ARR(LEN)

      INTEGER I, J
      real*GC__RSIZE RSRC, RDEST
      character*GC__RSIZE CSRC, CDEST
      EQUIVALENCE (RSRC, CSRC), (RDEST, CDEST)
#if defined(__alpha) || defined(ultrix)
      INTEGER dec_rsize, dec_isize
      common /gc__dec/dec_rsize, dec_isize
#endif


      DO I = 1, LEN
         RSRC = ARR(I)
         DO J = 1, gcrsize
            CDEST(gcrsize+1-J:gcrsize+1-J) = CSRC(J:J)
         ENDDO
         ARR(I) = RDEST
      ENDDO
      END


      SUBROUTINE GC__ISWAP(LEN, ARR)
C
C     Support routine for byte swapping GC integer arrays
C
      IMPLICIT NONE
      INTEGER LEN, ARR(LEN)

      INTEGER I, J
      integer*GC__ISIZE ISRC, IDEST
      character*GC__ISIZE CSRC, CDEST
      EQUIVALENCE (ISRC, CSRC), (IDEST, CDEST)
#if defined(__alpha) || defined(ultrix)
      INTEGER dec_rsize, dec_isize
      common /gc__dec/dec_rsize, dec_isize
#endif


      DO I = 1, LEN
         ISRC = ARR(I)
         DO J = 1, gcisize
            CDEST(gcisize+1-J:gcisize+1-J) = CSRC(J:J)
         ENDDO
         ARR(I) = IDEST
      ENDDO
      END
#endif
