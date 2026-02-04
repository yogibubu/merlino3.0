C=======================================================================
C  MAIN DRIVER - LOCAL SYMMETRY / COLLECTIVE INTERNAL COORDINATES
C=======================================================================

      PROGRAM MAIN_LOCAL_COORDS
      IMPLICIT NONE

C=========================== LIMITS ====================================
      INTEGER MAXAT, MAXBND, MAXDEG, MAXQ
      INTEGER MAXGRP, MAXORD, MAXRING, MAXRSZ, MAXCOL

      PARAMETER (MAXAT   = 200)
      PARAMETER (MAXBND  = 400)
      PARAMETER (MAXDEG  = 12)
      PARAMETER (MAXQ    = 3000)
      PARAMETER (MAXGRP  = 1500)
      PARAMETER (MAXORD  = 20)
      PARAMETER (MAXRING = 200)
      PARAMETER (MAXRSZ  = 12)
      PARAMETER (MAXCOL  = 4000)

C=========================== MOLECULE ==================================
      INTEGER NAT, NBOND, IERR
      INTEGER ZAT(MAXAT)
      INTEGER BND(2,MAXBND)
      DOUBLE PRECISION X(MAXAT), Y(MAXAT), ZC(MAXAT)

C=========================== ATOM CLASSES ==============================
      INTEGER ACLASS(MAXAT), NCLASS

C=========================== PRIMITIVES ================================
      INTEGER NQ
      INTEGER QTYPE(MAXQ), QSUB(MAXQ), QCENTER(MAXQ)
      INTEGER QATOM(4,MAXQ)

C=========================== GROUPING ==================================
      INTEGER QSIG(6,MAXQ)
      INTEGER NGRP, GRP_SIZE(MAXGRP)
      INTEGER GRP_Q(MAXORD,MAXGRP)

C=========================== RINGS =====================================
      INTEGER NRING
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)

      INTEGER NCLUST
      INTEGER CLUST_SIZE(MAXRING)
      INTEGER CLUST_RING(MAXRING,MAXRING)

C=========================== COLLECTIVE ================================
      INTEGER NCOL
      INTEGER COL_NTERM(MAXCOL)
      INTEGER COL_Q(MAXORD,MAXCOL)
      DOUBLE PRECISION COL_COEFF(MAXORD,MAXCOL)
      DOUBLE PRECISION COL_NORM(MAXCOL)

C=========================== INPUT =====================================
      CALL READ_XYZ('geom.xyz', NAT, ZAT, X, Y, ZC, IERR)
      IF (IERR .NE. 0) STOP 'ERROR reading geom.xyz'

      CALL READ_CONNECTIVITY('conn.dat', NBOND, BND, IERR)
      IF (IERR .NE. 0) STOP 'ERROR reading conn.dat'

C=========================== STEP 1 ====================================
      CALL BUILD_ATOM_CLASSES(
     &   MAXAT, MAXBND, MAXDEG,
     &   NAT, ZAT, NBOND, BND,
     &   ACLASS, NCLASS, IERR)

C=========================== STEP 2 ====================================
      CALL BUILD_PRIMITIVE_COORDS(
     &   MAXAT, MAXBND, MAXQ, MAXDEG,
     &   NAT, ZAT, NBOND, BND,
     &   X, Y, ZC,
     &   NQ, QTYPE, QSUB, QATOM, QCENTER,
     &   IERR)

C=========================== STEP 3 ====================================
      CALL BUILD_SIGNATURES(
     &   MAXQ,
     &   NQ, QTYPE, QSUB, QATOM, QCENTER,
     &   ACLASS,
     &   QSIG)

      CALL GROUP_PRIMITIVES(
     &   MAXQ, MAXGRP, MAXORD,
     &   NQ, QSIG,
     &   NGRP, GRP_SIZE, GRP_Q)

C=========================== STEP 3.5 ==================================
      CALL FIND_RINGS(
     &   MAXAT, MAXBND, MAXRING, MAXRSZ,
     &   NAT, NBOND, BND,
     &   NRING, RING_SIZE, RING_ATOM)

      CALL NORMALIZE_RINGS(
     &   MAXRING, MAXRSZ,
     &   NRING, RING_SIZE, RING_ATOM)

      CALL REMOVE_DUPLICATE_RINGS(
     &   MAXRING, MAXRSZ,
     &   NRING, RING_SIZE, RING_ATOM)

      CALL BUILD_RING_CLUSTERS(
     &   MAXRING,
     &   NRING, RING_SIZE, RING_ATOM,
     &   NCLUST, CLUST_SIZE, CLUST_RING)

C=========================== STEP 4 ====================================
      CALL BUILD_ALL_COLLECTIVE_COORDS(
     &   MAXQ, MAXORD, MAXCOL, MAXRSZ, MAXRING,
     &   NQ, QTYPE, QSUB, QATOM, QCENTER,
     &   NGRP, GRP_SIZE, GRP_Q,
     &   NRING, RING_SIZE, RING_ATOM,
     &   NCLUST, CLUST_SIZE, CLUST_RING,
     &   NCOL, COL_NTERM, COL_Q, COL_COEFF, COL_NORM)

C=========================== DIAGNOSTIC PRINT ==========================
      WRITE(*,*)
      WRITE(*,*) '===== STRUCTURAL SUMMARY ====='
      WRITE(*,*) 'NATOMS        = ', NAT
      WRITE(*,*) 'NBONDS        = ', NBOND
      WRITE(*,*) 'PRIMITIVES    = ', NQ
      WRITE(*,*) 'GROUPS        = ', NGRP
      WRITE(*,*) 'RINGS         = ', NRING
      WRITE(*,*) 'RING CLUSTERS = ', NCLUST
      WRITE(*,*) 'COLLECTIVE    = ', NCOL
      WRITE(*,*)

      STOP
      END

