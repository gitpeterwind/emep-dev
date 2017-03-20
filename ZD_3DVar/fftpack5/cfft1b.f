CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: cfft1b.f,v 1.1 2007/11/27 13:53:26 mkahnert Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CFFT1B (N, INC, C, LENC, WSAVE, LENSAV,
     1                  WORK, LENWRK, IER)
      INTEGER  N, INC, LENC, LENSAV, LENWRK, IER
      COMPLEX       C(LENC)
      REAL     WSAVE(LENSAV)     ,WORK(LENWRK)
C
      IER = 0
C
      IF (LENC .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('CFFT1B ', 6)
      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N))) + 4) THEN
        IER = 2
        CALL XERFFT ('CFFT1B ', 8)
      ELSEIF (LENWRK .LT. 2*N) THEN
        IER = 3
        CALL XERFFT ('CFFT1B ', 10)
      ENDIF
C
      IF (N .EQ. 1) RETURN
C
      IW1 = N+N+1
      CALL C1FM1B (N,INC,C,WORK,WSAVE,WSAVE(IW1),
     1                           WSAVE(IW1+1))
      RETURN
      END
