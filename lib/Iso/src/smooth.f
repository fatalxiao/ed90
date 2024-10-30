      SUBROUTINE SMOOTH(NROW, NCOL, NDIM, X, W, A, B, NCYCLE, ICYCLE, G,
     *   EPS1, EPS2, IFAULT, FX, PW, W1, WT, NW)
C
C     ALGORITHM AS 206  APPL. STATIST. (1984) VOL.33, NO.3
C
C     Subroutine to order a two-dimensional array using an algorithm of
C     Dykstra & Robertson (1982).   The ordering is done so that the
C     regression function is increasing in each independent variable.
C
C     Incorporates corrections from Applied Statistics vol.35(3), 
C     vol.36(1) and vol.40(1).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NROW,NCOL), W(NROW,NCOL), A(NROW,NCOL,4)
      DIMENSION B(NDIM,5), G(NROW,NCOL)
      DIMENSION FX(NDIM), PW(NDIM), W1(NDIM), WT(NDIM)
      INTEGER NW(NDIM)
      DATA ZERO/0.0d0/, DELTA/0.00001d0/, FRACT/0.5d0/
C
      IFAULT = 0
C Rprintf("epsilon =%f12.6",EPS)
C Rprintf("zero =%f12.6",ZERO)
C Rprintf("delta =%f12.6",DELTA)
C Rprintf("fract =%f12.6",FRACT)
C
C     Check that there are at least 2 rows and columns
C
      IF (NROW .LT. 2 .OR. NCOL .LT. 2) GO TO 120
C
C     Check that weights are positive or zero
C
      WSUM = ZERO
      WXSUM = ZERO
      WMIN = 1.0d+08
      DO 3 I = 1, NROW
        DO 4 J = 1, NCOL
          WW = W(I,J)
          IF (WW .LT. ZERO) GO TO 110
          IF (WW .LT. DELTA) GO TO 3
          WSUM = WSUM + WW
          WXSUM = WXSUM + WW * X(I,J)
          IF (WW .LT. WMIN) WMIN = WW
    4 CONTINUE
    3 CONTINUE
      IF (WSUM .LT. DELTA) GO TO 130
      WMEAN = WXSUM / WSUM
C
      DO 5 I = 1, NROW
      DO 6 J = 1, NCOL
        WW = W(I,J)
        A(I,J,3) = WW
        A(I,J,4) = X(I,J)
        IF (WW .GE. DELTA) GO TO 6
        A(I,J,3) = FRACT * WMIN
        A(I,J,4) = WMEAN
        ICT = ICT + 1
        IFAULT = 4
    6 CONTINUE
    5 CONTINUE
C
C     Initialize R and C to zero, and set up workspace
C
      IFLAG = 0
      DELR = ZERO
      DELC = ZERO
      ITIC = 0
    8 ITIC = ITIC + 1
      DO 10 I = 1, NROW
      DO 11 J = 1, NCOL
        G(I,J) = A(I,J,4)
        A(I,J,2) = ZERO
        A(I,J,1) = ZERO
   11 CONTINUE
   10 CONTINUE
C
C     Initialize counter for number of cycles
C
      ICOUNT = 0
      IF (IFLAG .EQ. 1) GO TO 55
      IF (ITIC .EQ. 3 .AND. DELC .GT. DELR) GO TO 55
C
C     Smooth over rows
C
   25 JCOUNT = 0
      DO 50 I = 1, NROW
      DO 30 J = 1, NCOL
        B(J,1) = G(I,J) - A(I,J,1)
        B(J,2) = A(I,J,3)
   30   CONTINUE
C
      CALL PAV(NCOL, NDIM, EPS2, B, 1, B(1,2), B(1,3), FX,
     * PW, W1, WT, NW)
C
      KCOUNT = 0
      DO 40 J = 1, NCOL
        ORD = B(J,3)
        A(I,J,1) = ORD - B(J,1)
        IF (ABS(ORD - G(I,J)) .LT. EPS1) KCOUNT = KCOUNT + 1
        G(I,J) = ORD
   40   CONTINUE
C
C     Determine if there is no change in the Ith row from the previous
C     iteration
C
      IF (KCOUNT .EQ. NCOL) JCOUNT = JCOUNT + 1
   50 CONTINUE
C
      ICOUNT = ICOUNT + 1
      IF (ICOUNT .EQ. 2 .AND. IFLAG .EQ. 1)
     *    CALL DIST(A(1,1,1), NROW, NCOL, DELR, IFLAG)
      IF (ICOUNT .EQ. 2 .AND. IFLAG .EQ. 2 .AND. ITIC .EQ. 2) GO TO 8
      IF (ICOUNT .EQ. 1) GO TO 55
      IF (NCYCLE .EQ. ICOUNT) GO TO 100
C
C     Determine if there has been no change in all rows from the
C     previous iteration
C
      IF (JCOUNT .EQ. NROW) GO TO 90
C
C     Smooth over columns
C
   55 LCOUNT = 0
      DO 80 J = 1, NCOL
      DO 60 I = 1, NROW
        A(I,J,2) = G(I,J) - A(I,J,2)
   60  CONTINUE
C
      CALL PAV(NROW, NDIM, EPS2, A(1,J,2), 1, A(1,J,3), B(1,3), FX,
     * PW, W1, WT, NW)
C
      MCOUNT = 0
      DO 70 I = 1, NROW
        ORD = B(I,3)
        A(I,J,2) = ORD - A(I,J,2)
        IF (ABS(ORD - G(I,J)) .LT. EPS1) MCOUNT = MCOUNT + 1
        G(I,J) = ORD
   70   CONTINUE
C
C     Determine if there is no change in the Jth column from the
C     previous iteration
C
      IF (MCOUNT .EQ. NROW) LCOUNT = LCOUNT + 1
   80 CONTINUE
C
      ICOUNT = ICOUNT + 1
      IF (ICOUNT .EQ. 2 .AND. IFLAG .EQ. 0)
     *    CALL DIST(A(1,1,2), NROW, NCOL, DELC, IFLAG)
      IF (ICOUNT .EQ. 2 .AND. IFLAG .EQ. 1) GO TO 8
      IF (ICOUNT .EQ. 1) GO TO 25
C
C     Determine if there is has been no change in any column from the
C     previous iteration
C
      IF (LCOUNT .EQ. NCOL) GO TO 90
C
C     Check if number of cycles has been reached
C
      IF (NCYCLE .EQ. ICOUNT) GO TO 100
      GO TO 25
   90 ICYCLE = ICOUNT
      RETURN
C
  100 ICYCLE = ICOUNT
      IFAULT = IFAULT + 1
      RETURN
  110 IFAULT = 2
      RETURN
  120 IFAULT = 3
      RETURN
  130 IFAULT = 6
      RETURN
      END


      SUBROUTINE PAV(K, NDIM, EPS2, X, IORDER, W, FINALX,
     *FX, PW, W1, WT, NW)
C
C     ALGORITHM AS 206.1  APPL. STATIST. (1984) VOL.33, NO.3
C
C     Apply pool adjacent violators theorem
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NW(NDIM)
      dimension X(NDIM), W(NDIM), FINALX(NDIM), FX(NDIM), PW(NDIM)
      dimension W1(NDIM), WT(NDIM)
C
C     DATA EPS/1.0d-06/
C
C     Set up workspace
C
      NWC = K
      DO 10 I = 1, K
      NW(I) = 1
      FX(I) = X(I)
      IF (IORDER .EQ. 0) FX(I) = -FX(I)
      WT(I) = W(I)
      PW(I) = WT(I) * FX(I)
      W1(I) = W(I)
   10 CONTINUE
      IBEL = K - 1
   20 I = 0
   30 I = I + 1
   35 IF (I .GT. IBEL) GO TO 50
      I1 = I + 1
C
C     Determine if pooling is required
C
      IF (FX(I) - FX(I1) .LE. EPS2) GO TO 30
C
C     Pool the adjacent values
C
      PW(I) = PW(I) + PW(I1)
      W1(I) = W1(I) + W1(I1)
      FX(I) = PW(I) / W1(I)
      NW(I) = NW(I) + NW(I1)
      NWC = NWC - 1
      IF (I1 .GT. IBEL) GO TO 45
      DO 40 J = I1, IBEL
      J1 = J + 1
      PW(J) = PW(J1)
      W1(J) = W1(J1)
      FX(J) = FX(J1)
      NW(J) = NW(J1)
   40 CONTINUE
   45 IBEL = IBEL - 1
      GO TO 35
   50 ICOUNT = 0
      IF (IBEL .LE. 0) GO TO 70
C
C     Determine if all values are ordered
C
      DO 60 L = 1, IBEL
      IF (FX(L) - FX(L+1) .LE. EPS2) ICOUNT = ICOUNT + 1
   60 CONTINUE
      IF (ICOUNT .NE. IBEL) GO TO 20
C
C     Recover final ordered values
C
   70 J = 1
      JL = 1
      JU = NW(1)
   80 DO 90 L = JL,JU
      FINALX(L) = FX(J)
   90 CONTINUE
      J = J + 1
      IF (J .GT. NWC) GO TO 100
      JL = JU + 1
      JU = JU + NW(J)
      GO TO 80
  100 IF (IORDER .EQ. 1) RETURN
      DO 110 I = 1, K
      FINALX(I) = -FINALX(I)
  110 CONTINUE
      RETURN
      END


      SUBROUTINE DIST(A1, NROW, NCOL, DEL, IFLAG)
C
C     ALGORITHM AS 206.2  APPL. STATIST. (1984) VOL.33, NO.3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension A1(NROW,NCOL)
C
      DATA ZERO/0.0d0/
C
      DEL = ZERO
      DO 20 I = 1, NROW
      DO 30 J = 1, NCOL
        DEL = DEL + A1(I,J) * A1(I,J)
   30 CONTINUE
   20 CONTINUE
      IFLAG = IFLAG + 1
      RETURN
      END
