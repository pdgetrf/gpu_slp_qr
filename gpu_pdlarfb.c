/* pdlarfb.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

http://www.netlib.org/f2c/libf2c.zip
 */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "f2c.h"
#include "util_gpu.h"


/* Table of constant values */

static doublereal dzero = 0.;
static doublereal done = 1.;
static doublereal mone = -1.;
static int c_n1 = -1;
int ione = 1;

/* Subroutine */ int gpu_pdlarfb_(char *side, char *trans, char *direct, char *
		storev, int *m, int *n, int *k, doublereal *v, int *
		iv, int *jv, int *descv, doublereal *t, doublereal *c__, 
		int *ic, int *jc, int *descc, doublereal *work, ftnlen 
		side_len, ftnlen trans_len, ftnlen direct_len, ftnlen storev_len,
		double *A2, int *descA2, double *W, double *V, double *pinnbuf)
{
	/* System generated locals */
	int i__1, i__2, i__3, i__4;

	/* Local variables */
	static int ilastcol, ilastrow;
	extern /* Subroutine */ int pb_topget__(int *, char *, char *, char *,
			ftnlen, ftnlen, ftnlen);
	static int ii, jj, kp, kq, lv, lw, ldc, iic, jjc, mpc, nqc, mbv, ldv, 
			   nbv, iiv, jjv, ipt, ipv, ipw, mqv, npv, mpc0, nqc0, ipw1, mqv0, 
			   npv0, ioff, wide, itop;
	static char uplo[1];
	static int iibeg, ibase, jjbeg;
	extern int iceil_(int *, int *);
	static int ioffc, iiend, iccol, jjend;
	extern /* Subroutine */ int dgemm_(char *, char *, int *, int *, 
			int *, doublereal *, doublereal *, int *, doublereal *, 
			int *, doublereal *, doublereal *, int *, ftnlen, ftnlen);
	extern logical lsame_(char *, char *, ftnlen, ftnlen);
	static int ileft, ioffv, npcol, ivcol, icrow, mycol;
	extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
			int *, int *, doublereal *, doublereal *, int *, 
			doublereal *, int *, ftnlen, ftnlen, ftnlen, ftnlen);
	static int ictxt, iinxt, jjnxt, nprow, ivrow, myrow, icoffc, height, 
			   iroffc, icoffv;
	extern /* Subroutine */ int dlacpy_(char *, int *, int *, 
			doublereal *, int *, doublereal *, int *, ftnlen), 
		   dlaset_(char *, int *, int *, doublereal *, doublereal *, 
				   doublereal *, int *, ftnlen);
	static int iright, iroffv;
	extern /* Subroutine */ int blacs_gridinfo__(int *, int *, 
			int *, int *, int *);
	extern int numroc_(int *, int *, int *, int *, 
			int *);
	static int mydist;
	extern /* Subroutine */ int dgebr2d_(int *, char *, char *, int *,
			int *, doublereal *, int *, int *, int *, ftnlen,
			ftnlen);
	static char transt[1];
	extern /* Subroutine */ int dgebs2d_(int *, char *, char *, int *,
			int *, doublereal *, int *, ftnlen, ftnlen), infog1l_(
				int *, int *, int *, int *, int *, int *, 
				int *), infog2l_(int *, int *, int *, int *, 
					int *, int *, int *, int *, int *, int *, 
					int *), dtrbr2d_(int *, char *, char *, char *, char *, 
						int *, int *, doublereal *, int *, int *, int 
						*, ftnlen, ftnlen, ftnlen, ftnlen), dgsum2d_(int *, char *, 
							char *, int *, int *, doublereal *, int *, int *, 
							int *, ftnlen, ftnlen), dtrbs2d_(int *, char *, char *, 
								char *, char *, int *, int *, doublereal *, int *, 
								ftnlen, ftnlen, ftnlen, ftnlen), pbdtran_(int *, char *, char 
									*, int *, int *, int *, doublereal *, int *, 
									doublereal *, doublereal *, int *, int *, int *, 
									int *, int *, doublereal *, ftnlen, ftnlen);
	static char colbtop[1];
	static logical forward;
	static char rowbtop[1];

	static cudaStream_t fstream;
	cudaStreamCreate(&fstream);


	/*  -- ScaLAPACK auxiliary routine (version 1.7) -- */
	/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
	/*     and University of California, Berkeley. */
	/*     May 1, 1997 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  PDLARFB applies a real block reflector Q or its transpose Q**T to a */
	/*  real distributed M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) */
	/*  from the left or the right. */

	/*  Notes */
	/*  ===== */

	/*  Each global data object is described by an associated description */
	/*  vector.  This vector stores the information required to establish */
	/*  the mapping between an object element and its corresponding process */
	/*  and memory location. */

	/*  Let A be a generic term for any 2D block cyclicly distributed array. */
	/*  Such a global array has an associated description vector DESCA. */
	/*  In the following comments, the character _ should be read as */
	/*  "of the global array". */

	/*  NOTATION        STORED IN      EXPLANATION */
	/*  --------------- -------------- -------------------------------------- */
	/*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case, */
	/*                                 DTYPE_A = 1. */
	/*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating */
	/*                                 the BLACS process grid A is distribu- */
	/*                                 ted over. The context itself is glo- */
	/*                                 bal, but the handle (the int */
	/*                                 value) may vary. */
	/*  M_A    (global) DESCA( M_ )    The number of rows in the global */
	/*                                 array A. */
	/*  N_A    (global) DESCA( N_ )    The number of columns in the global */
	/*                                 array A. */
	/*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute */
	/*                                 the rows of the array. */
	/*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute */
	/*                                 the columns of the array. */
	/*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first */
	/*                                 row of the array A is distributed. */
	/*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the */
	/*                                 first column of the array A is */
	/*                                 distributed. */
	/*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local */
	/*                                 array.  LLD_A >= MAX(1,LOCr(M_A)). */

	/*  Let K be the number of rows or columns of a distributed matrix, */
	/*  and assume that its process grid has dimension p x q. */
	/*  LOCr( K ) denotes the number of elements of K that a process */
	/*  would receive if K were distributed over the p processes of its */
	/*  process column. */
	/*  Similarly, LOCc( K ) denotes the number of elements of K that a */
	/*  process would receive if K were distributed over the q processes of */
	/*  its process row. */
	/*  The values of LOCr() and LOCc() may be determined via a call to the */
	/*  ScaLAPACK tool function, NUMROC: */
	/*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ), */
	/*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ). */
	/*  An upper bound for these quantities may be computed by: */
	/*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A */
	/*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A */

	/*  Arguments */
	/*  ========= */

	/*  SIDE    (global input) CHARACTER */
	/*          = 'L': apply Q or Q**T from the Left; */
	/*          = 'R': apply Q or Q**T from the Right. */

	/*  TRANS   (global input) CHARACTER */
	/*          = 'N':  No transpose, apply Q; */
	/*          = 'T':  Transpose, apply Q**T. */

	/*  DIRECT  (global input) CHARACTER */
	/*          Indicates how Q is formed from a product of elementary */
	/*          reflectors */
	/*          = 'F': Q = H(1) H(2) . . . H(k) (Forward) */
	/*          = 'B': Q = H(k) . . . H(2) H(1) (Backward) */

	/*  STOREV  (global input) CHARACTER */
	/*          Indicates how the vectors which define the elementary */
	/*          reflectors are stored: */
	/*          = 'C': Columnwise */
	/*          = 'R': Rowwise */

	/*  M       (global input) INTEGER */
	/*          The number of rows to be operated on i.e the number of rows */
	/*          of the distributed submatrix sub( C ). M >= 0. */

	/*  N       (global input) INTEGER */
	/*          The number of columns to be operated on i.e the number of */
	/*          columns of the distributed submatrix sub( C ). N >= 0. */

	/*  K       (global input) INTEGER */
	/*          The order of the matrix T (= the number of elementary */
	/*          reflectors whose product defines the block reflector). */

	/*  V       (local input) DOUBLE PRECISION pointer into the local memory */
	/*          to an array of dimension ( LLD_V, LOCc(JV+K-1) ) if */
	/*          STOREV = 'C', ( LLD_V, LOCc(JV+M-1)) if STOREV = 'R' and */
	/*          SIDE = 'L', ( LLD_V, LOCc(JV+N-1) ) if STOREV = 'R' and */
	/*          SIDE = 'R'. It contains the local pieces of the distributed */
	/*          vectors V representing the Householder transformation. */
	/*          See further details. */
	/*          If STOREV = 'C' and SIDE = 'L', LLD_V >= MAX(1,LOCr(IV+M-1)); */
	/*          if STOREV = 'C' and SIDE = 'R', LLD_V >= MAX(1,LOCr(IV+N-1)); */
	/*          if STOREV = 'R', LLD_V >= LOCr(IV+K-1). */

	/*  IV      (global input) INTEGER */
	/*          The row index in the global array V indicating the first */
	/*          row of sub( V ). */

	/*  JV      (global input) INTEGER */
	/*          The column index in the global array V indicating the */
	/*          first column of sub( V ). */

	/*  DESCV   (global and local input) INTEGER array of dimension DLEN_. */
	/*          The array descriptor for the distributed matrix V. */

	/*  T       (local input) DOUBLE PRECISION array, dimension MB_V by MB_V */
	/*          if STOREV = 'R' and NB_V by NB_V if STOREV = 'C'. The trian- */
	/*          gular matrix T in the representation of the block reflector. */

	/*  C       (local input/local output) DOUBLE PRECISION pointer into the */
	/*          local memory to an array of dimension (LLD_C,LOCc(JC+N-1)). */
	/*          On entry, the M-by-N distributed matrix sub( C ). On exit, */
	/*          sub( C ) is overwritten by Q*sub( C ) or Q'*sub( C ) or */
	/*          sub( C )*Q or sub( C )*Q'. */

	/*  IC      (global input) INTEGER */
	/*          The row index in the global array C indicating the first */
	/*          row of sub( C ). */

	/*  JC      (global input) INTEGER */
	/*          The column index in the global array C indicating the */
	/*          first column of sub( C ). */

	/*  DESCC   (global and local input) INTEGER array of dimension DLEN_. */
	/*          The array descriptor for the distributed matrix C. */

	/*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK) */
	/*          If STOREV = 'C', */
	/*            if SIDE = 'L', */
	/*              LWORK >= ( NqC0 + MpC0 ) * K */
	/*            else if SIDE = 'R', */
	/*              LWORK >= ( NqC0 + MAX( NpV0 + NUMROC( NUMROC( N+ICOFFC, */
	/*                         NB_V, 0, 0, NPCOL ), NB_V, 0, 0, LCMQ ), */
	/*                         MpC0 ) ) * K */
	/*            end if */
	/*          else if STOREV = 'R', */
	/*            if SIDE = 'L', */
	/*              LWORK >= ( MpC0 + MAX( MqV0 + NUMROC( NUMROC( M+IROFFC, */
	/*                         MB_V, 0, 0, NPROW ), MB_V, 0, 0, LCMP ), */
	/*                         NqC0 ) ) * K */
	/*            else if SIDE = 'R', */
	/*              LWORK >= ( MpC0 + NqC0 ) * K */
	/*            end if */
	/*          end if */

	/*          where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ), */

	/*          IROFFV = MOD( IV-1, MB_V ), ICOFFV = MOD( JV-1, NB_V ), */
	/*          IVROW = INDXG2P( IV, MB_V, MYROW, RSRC_V, NPROW ), */
	/*          IVCOL = INDXG2P( JV, NB_V, MYCOL, CSRC_V, NPCOL ), */
	/*          MqV0 = NUMROC( M+ICOFFV, NB_V, MYCOL, IVCOL, NPCOL ), */
	/*          NpV0 = NUMROC( N+IROFFV, MB_V, MYROW, IVROW, NPROW ), */

	/*          IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ), */
	/*          ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ), */
	/*          ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ), */
	/*          MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ), */
	/*          NpC0 = NUMROC( N+ICOFFC, MB_C, MYROW, ICROW, NPROW ), */
	/*          NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ), */

	/*          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions; */
	/*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling */
	/*          the subroutine BLACS_GRIDINFO. */

	/*  Alignment requirements */
	/*  ====================== */

	/*  The distributed submatrices V(IV:*, JV:*) and C(IC:IC+M-1,JC:JC+N-1) */
	/*  must verify some alignment properties, namely the following */
	/*  expressions should be true: */

	/*  If STOREV = 'Columnwise' */
	/*    If SIDE = 'Left', */
	/*      ( MB_V.EQ.MB_C .AND. IROFFV.EQ.IROFFC .AND. IVROW.EQ.ICROW ) */
	/*    If SIDE = 'Right', */
	/*      ( MB_V.EQ.NB_C .AND. IROFFV.EQ.ICOFFC ) */
	/*  else if STOREV = 'Rowwise' */
	/*    If SIDE = 'Left', */
	/*      ( NB_V.EQ.MB_C .AND. ICOFFV.EQ.IROFFC ) */
	/*    If SIDE = 'Right', */
	/*      ( NB_V.EQ.NB_C .AND. ICOFFV.EQ.ICOFFC .AND. IVCOL.EQ.ICCOL ) */
	/*  end if */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Quick return if possible */

	/* Parameter adjustments */
	--work;
	--descc;
	--c__;
	--t;
	--descv;
	--v;

	/* Function Body */
	if (*m <= 0 || *n <= 0 || *k <= 0) {
		return 0;
	}

	/*     Get grid parameters */

	ictxt = descc[2];
	blacs_gridinfo__(&ictxt, &nprow, &npcol, &myrow, &mycol);

	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
		*(unsigned char *)transt = 'T';
	} else {
		*(unsigned char *)transt = 'N';
	}
	forward = lsame_(direct, "F", (ftnlen)1, (ftnlen)1);
	if (forward) {
		*(unsigned char *)uplo = 'U';
	} else {
		*(unsigned char *)uplo = 'L';
	}

	infog2l_(iv, jv, &descv[1], &nprow, &npcol, &myrow, &mycol, &iiv, &jjv, &
			ivrow, &ivcol);
	infog2l_(ic, jc, &descc[1], &nprow, &npcol, &myrow, &mycol, &iic, &jjc, &
			icrow, &iccol);
	ldc = descc[9];
	ldv = descv[9];
	iic = min(iic,ldc);
	iiv = min(iiv,ldv);
	iroffc = (*ic - 1) % descc[5];
	icoffc = (*jc - 1) % descc[6];
	mbv = descv[5];
	nbv = descv[6];
	iroffv = (*iv - 1) % mbv;
	icoffv = (*jv - 1) % nbv;
	i__1 = *m + iroffc;
	mpc = numroc_(&i__1, &descc[5], &myrow, &icrow, &nprow); // mpc: rows I own
	i__1 = *n + icoffc;
	nqc = numroc_(&i__1, &descc[6], &mycol, &iccol, &npcol); // nqc: columns I own
	if (mycol == iccol) {
		nqc -= icoffc;
	}
	if (myrow == icrow) {
		mpc -= iroffc;
	}
	/* Computing MIN */
	/* Computing MAX */
	i__3 = 1, i__4 = jjc + nqc - 1;
	i__1 = jjc, i__2 = max(i__3,i__4);
	jjc = min(i__1,i__2);
	/* Computing MIN */
	/* Computing MAX */
	i__3 = 1, i__4 = numroc_(&descv[4], &nbv, &mycol, &descv[8], &npcol);
	i__1 = jjv, i__2 = max(i__3,i__4);
	jjv = min(i__1,i__2);
	ioffc = iic + (jjc - 1) * ldc;
	ioffv = iiv + (jjv - 1) * ldv;

	if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) 
	{
		/*        V is stored columnwise */

		if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) 
		{

			/*           Form  Q*sub( C )  or  Q'*sub( C  ) */
			/*           Locally V( IOFFV ) is MPV x K, C( IOFFC ) is MPC x NQC */
			/*           WORK( IPV ) is MPC x K = V( IOFFV ), MPC = MPV */
			/*           WORK( IPW ) is NQC x K = C( IOFFC )' * V( IOFFV ) */

			ipv = 1;
			ipw = ipv + mpc * *k;
			lv = max(1,mpc);
			lw = max(1,nqc);

			/*           Broadcast V to the other process columns. */

			pb_topget__(&ictxt, "Broadcast", "Rowwise", rowbtop, (ftnlen)9, (
						ftnlen)7, (ftnlen)1);
			
			if (mycol == ivcol) 
			{
				dgebs2d_(&ictxt, "Rowwise", rowbtop, &mpc, k, &v[ioffv], &ldv,
						(ftnlen)7, (ftnlen)1);
				if (myrow == ivrow) 
				{
					dtrbs2d_(&ictxt, "Rowwise", rowbtop, uplo, "Non unit", k, 
							k, &t[1], &nbv, (ftnlen)7, (ftnlen)1, (ftnlen)1, (
								ftnlen)8);
				}
				dlacpy_("All", &mpc, k, &v[ioffv], &ldv, &work[ipv], &lv, (
							ftnlen)3);
			} 
			else 
			{
				dgebr2d_(&ictxt, "Rowwise", rowbtop, &mpc, k, &work[ipv], &lv,
						&myrow, &ivcol, (ftnlen)7, (ftnlen)1);
				if (myrow == ivrow) 
				{
					dtrbr2d_(&ictxt, "Rowwise", rowbtop, uplo, "Non unit", k, 
							k, &t[1], &nbv, &myrow, &ivcol, (ftnlen)7, (
								ftnlen)1, (ftnlen)1, (ftnlen)8);
				}
			}

			if (forward) 
			{
				/*              WORK(IPV) = ( V1 ) where V1 is unit lower triangular, */
				/*                          ( V2 ) zeroes upper triangular part of V1 */
				mydist = (myrow - ivrow + nprow) % nprow;
				/* Computing MAX */
				i__1 = 0, i__2 = mydist * mbv - iroffv;
				itop = max(i__1,i__2);
				iibeg = iiv;
				iiend = iibeg + mpc - 1;
				/* Computing MIN */
				i__1 = iceil_(&iibeg, &mbv) * mbv;
				iinxt = min(i__1,iiend);

L10:
				if (*k - itop > 0) {
					i__1 = iinxt - iibeg + 1;
					i__2 = *k - itop;
					dlaset_("Upper", &i__1, &i__2, &dzero, &done, &work[ipv 
							+ iibeg - iiv + itop * lv], &lv, (ftnlen)5);
					mydist += nprow;
					itop = mydist * mbv - iroffv;
					iibeg = iinxt + 1;
					/* Computing MIN */
					i__1 = iinxt + mbv;
					iinxt = min(i__1,iiend);
					goto L10;
				}

			} 
			else 
			{

				/*              WORK(IPV) = ( V1 ) where V2 is unit upper triangular, */
				/*                          ( V2 ) zeroes lower triangular part of V2 */

				jj = jjv;
				ioff = (*iv + *m - *k - 1) % mbv;
				i__1 = *iv + *m - *k;
				infog1l_(&i__1, &mbv, &nprow, &myrow, &descv[7], &ii, &
						ilastrow);
				i__1 = *k + ioff;
				kp = numroc_(&i__1, &mbv, &myrow, &ilastrow, &nprow);
				if (myrow == ilastrow) {
					kp -= ioff;
				}
				mydist = (myrow - ilastrow + nprow) % nprow;
				itop = mydist * mbv - ioff;
				/* Computing MIN */
				i__1 = itop + mbv;
				ibase = min(i__1,*k);
				/* Computing MIN */
				i__1 = max(0,itop);
				itop = min(i__1,*k);

L20:
				if (jj <= jjv + *k - 1) {
					height = ibase - itop;
					i__1 = itop - jj + jjv;
					dlaset_("All", &kp, &i__1, &dzero, &dzero, &work[ipv + ii 
							- iiv + (jj - jjv) * lv], &lv, (ftnlen)3);
					dlaset_("Lower", &kp, &height, &dzero, &done, &work[ipv 
							+ ii - iiv + itop * lv], &lv, (ftnlen)5);
					/* Computing MAX */
					i__1 = 0, i__2 = kp - height;
					kp = max(i__1,i__2);
					ii += height;
					jj = jjv + ibase;
					mydist += nprow;
					itop = mydist * mbv - ioff;
					/* Computing MIN */
					i__1 = itop + mbv;
					ibase = min(i__1,*k);
					itop = min(itop,*k);
					goto L20;
				}
			}

			/*
			 * xxx
			 */

			/*           WORK( IPW ) = C( IOFFC )' * V  (NQC x MPC x K) -> NQC x K */

			if (mpc > 0) 
			{
#define GPU
#ifdef GPU
				if (nqc > 0)
				{
					// copy V to GPU
					cublasSetMatrix(mpc, *k, sizeof(double), &work[ipv], lv, V, mpc);

					int iiA2, jjA2;
					infog2l_(ic, ic, descA2, &nprow, &npcol, &myrow, &mycol, 
							&iiA2, &jjA2, &icrow, &iccol);
					iiA2--; jjA2--;

					/*
					if (myrow==2 && mycol==0)
						printf ("\n(%d,%d): k=%d, mpc=%d, nqc=%d, ic=%d, jc=%d, iiA2=%d, jjA2=%d, icrow=%d, iccol=%d, ldc=%d\n", 
								myrow, mycol, *k, mpc, nqc, *ic, *jc, 
								iiA2, jjA2, icrow, iccol, ldc);
								*/

					// perform GEMM
					cublasDgemm(MagmaTrans, MagmaNoTrans, nqc, *k, mpc, done, 
							&A2[jjA2*ldc+iiA2], ldc,
							V, mpc, dzero, W, nqc);
					/*
					double *A3;
					TESTING_DEVALLOC (A3, double, mpc*nqc);
					cublasSetMatrix(mpc, nqc, sizeof(double), &c__[ioffc], ldc, A3, mpc);
					cublasDgemm(MagmaTrans, MagmaNoTrans, nqc, *k, mpc, done, 
							A3, mpc,
							V, mpc, dzero, W, nqc);
					*/

					// copy W back
					cublasGetMatrix(nqc, *k, sizeof(double), W, nqc, &work[ipw], lw);
				}

#else
				dgemm_("Transpose", "No transpose", &nqc, k, &mpc, &done, &
						c__[ioffc], &ldc, &work[ipv], &lv, &dzero, &work[ipw],
						&lw, (ftnlen)9, (ftnlen)12);
#endif
			} 
			else 
			{
				dlaset_("All", &nqc, k, &dzero, &dzero, &work[ipw], &lw, (ftnlen)3);
			}

			dgsum2d_(&ictxt, "Columnwise", " ", &nqc, k, &work[ipw], &lw, &
					ivrow, &mycol, (ftnlen)10, (ftnlen)1);


			if (myrow == ivrow) 
			{

				/*              WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T */

				dtrmm_("Right", uplo, transt, "Non unit", &nqc, k, &done, &t[
						1], &nbv, &work[ipw], &lw, (ftnlen)5, (ftnlen)1, (
							ftnlen)1, (ftnlen)8);
				dgebs2d_(&ictxt, "Columnwise", " ", &nqc, k, &work[ipw], &lw, 
						(ftnlen)10, (ftnlen)1);
			} 
			else 
			{
				dgebr2d_(&ictxt, "Columnwise", " ", &nqc, k, &work[ipw], &lw, 
						&ivrow, &mycol, (ftnlen)10, (ftnlen)1);
			}

			/*               C            C      -     V       *     W' */
			/*           C( IOFFC ) = C( IOFFC ) - WORK( IPV ) * WORK( IPW )' */
			/*                        MPC x NQC    MPC x K         K x NQC */
#ifdef GPU
			infog2l_(&ione, jc, &descc[1], &nprow, &npcol, &myrow, &mycol, 
					&iic, &jjc, &icrow, &iccol);

			if (mpc>0 && nqc>0)
			{
				cublasSetMatrix(nqc, *k, sizeof(double), &work[ipw], lw, W, nqc);

				int iiA2, jjA2;
				int jnc = *jc-nbv;
				infog2l_(&ione, &jnc, descA2, &nprow, &npcol, &myrow, &mycol, 
						&iiA2, &jjA2, &icrow, &iccol);
				iiA2--; jjA2--;

				if (mycol==iccol)// i have the next panel
				{
					// send the next panel to host
#define ASYN
#ifdef ASYN
					cudaMemcpyAsync	(pinnbuf, &A2[jjA2*ldc+iiA2],
							ldc*nbv*sizeof(double), cudaMemcpyDeviceToHost, fstream);
#endif
					
					//printf ("(%d,%d) has the panel when ic=%d, jc=%d, copying A2(%d,%d) of size %dx%d, mpc=%d, npc=%d\n"
					//, myrow, mycol, *ic, *jc, iiA2, jjA2, ldc, nbv, mpc, nqc);

					// GPU performs the update 
					if (nqc>*k)
					{
						int iiA3, jjA3;
						infog2l_(ic, jc, descA2, &nprow, &npcol, &myrow, &mycol, 
								&iiA3, &jjA3, &icrow, &iccol);
						iiA3--;		jjA3--;
						cublasDgemm(MagmaNoTrans, MagmaTrans, mpc, nqc-*k, *k, 
								mone, V, mpc,
								W+*k, nqc,	done, A2+jjA3*ldc+iiA3, ldc);
					}


					// CPU performs the update 
#ifdef ASYN
					cudaStreamSynchronize	(fstream); 	
					memcpy (&c__[(jjc-1)*ldc+(iic-1)+1], pinnbuf, ldc*nbv*sizeof(double));
#else
					cublasGetMatrix(ldc, nbv, sizeof(double), &A2[jjA2*ldc+iiA2], ldc, &c__[(jjc-1)*ldc+(iic-1)+1], ldc);
#endif
					dgemm_("No transpose", "Transpose", &mpc, k, k, &mone, 
							&work[ipv], &lv, &work[ipw], &lw, &done, &c__[ioffc], &ldc, 
							(ftnlen)12, (ftnlen)9);
				}
				else
				{
					// perform GEMM
					infog2l_(ic, &jnc, descA2, &nprow, &npcol, &myrow, &mycol, 
							&iiA2, &jjA2, &icrow, &iccol);
					iiA2--;		jjA2--;
					//printf ("\n(%d,%d): k=%d, mpc=%d, nqc=%d, ic=%d, jc=%d, iiA2=%d, jjA2=%d, icrow=%d, iccol=%d, ldc=%d\n", 
					//		myrow, mycol, *k, mpc, nqc, *ic, *jc, 
					//		iiA2, jjA2, icrow, iccol, ldc);
					cublasDgemm(MagmaNoTrans, MagmaTrans, mpc, nqc, *k, mone, V, mpc,
							W, nqc,	done, &A2[jjA2*ldc+iiA2], ldc);
				}
					
				cudaThreadSynchronize();

				// copy W back
				//cublasGetMatrix(mpc, nqc, sizeof(double), A2, mpc, &c__[ioffc], ldc);
			}
			else if (mycol==iccol) // i'm in the next panel but I don't have things for the trailing update
			{
				int iiA2, jjA2;
				int jnc = *jc-nbv;
				infog2l_(&ione, &jnc, descA2, &nprow, &npcol, &myrow, &mycol, 
						&iiA2, &jjA2, &icrow, &iccol);
				iiA2--; jjA2--;

				// send the next panel to host
				cudaMemcpyAsync	(pinnbuf, &A2[jjA2*ldc+iiA2],
						ldc*nbv*sizeof(double), cudaMemcpyDeviceToHost, fstream);

				cudaStreamSynchronize	(fstream); 	
				memcpy (&c__[(jjc-1)*ldc+(iic-1)+1], pinnbuf, ldc*nbv*sizeof(double));
			}
#else
			dgemm_("No transpose", "Transpose", &mpc, &nqc, k, &mone, &work[
					ipv], &lv, &work[ipw], &lw, &done, &c__[ioffc], &ldc, 
					(ftnlen)12, (ftnlen)9);
#endif

		} 
		else 
		{

			printf ("not supported\n");

			/*           Form sub( C )*Q or sub( C )*Q' */

			/*           ICOFFC = IROFFV is required by the current transposition */
			/*           routine PBDTRAN */

			i__1 = *n + iroffv;
			npv0 = numroc_(&i__1, &mbv, &myrow, &ivrow, &nprow);
			if (myrow == ivrow) {
				npv = npv0 - iroffv;
			} else {
				npv = npv0;
			}
			if (mycol == iccol) {
				nqc0 = nqc + icoffc;
			} else {
				nqc0 = nqc;
			}

			/*           Locally V( IOFFV ) is NPV x K C( IOFFC ) is MPC x NQC */
			/*           WORK( IPV ) is K x NQC0 = [ . V( IOFFV ) ]' */
			/*           WORK( IPW ) is NPV0 x K = [ . V( IOFFV )' ]' */
			/*           WORK( IPT ) is the workspace for PBDTRAN */

			ipv = 1;
			ipw = ipv + *k * nqc0;
			ipt = ipw + npv0 * *k;
			lv = max(1,*k);
			lw = max(1,npv0);

			if (mycol == ivcol) {
				if (myrow == ivrow) {
					dlaset_("All", &iroffv, k, &dzero, &dzero, &work[ipw], &
							lw, (ftnlen)3);
					ipw1 = ipw + iroffv;
					dlacpy_("All", &npv, k, &v[ioffv], &ldv, &work[ipw1], &lw,
							(ftnlen)3);
				} else {
					ipw1 = ipw;
					dlacpy_("All", &npv, k, &v[ioffv], &ldv, &work[ipw1], &lw,
							(ftnlen)3);
				}

				if (forward) {

					/*                 WORK(IPW) = ( . V1' V2' )' where V1 is unit lower */
					/*                 triangular, zeroes upper triangular part of V1 */

					mydist = (myrow - ivrow + nprow) % nprow;
					/* Computing MAX */
					i__1 = 0, i__2 = mydist * mbv - iroffv;
					itop = max(i__1,i__2);
					iibeg = iiv;
					iiend = iibeg + npv - 1;
					/* Computing MIN */
					i__1 = iceil_(&iibeg, &mbv) * mbv;
					iinxt = min(i__1,iiend);

L30:
					if (*k - itop > 0) {
						i__1 = iinxt - iibeg + 1;
						i__2 = *k - itop;
						dlaset_("Upper", &i__1, &i__2, &dzero, &done, &work[
								ipw1 + iibeg - iiv + itop * lw], &lw, (ftnlen)
								5);
						mydist += nprow;
						itop = mydist * mbv - iroffv;
						iibeg = iinxt + 1;
						/* Computing MIN */
						i__1 = iinxt + mbv;
						iinxt = min(i__1,iiend);
						goto L30;
					}

				} else {

					/*                 WORK( IPW ) = ( . V1' V2' )' where V2 is unit upper */
					/*                 triangular, zeroes lower triangular part of V2. */

					jj = jjv;
					i__1 = *iv + *n - *k;
					infog1l_(&i__1, &mbv, &nprow, &myrow, &descv[7], &ii, &
							ilastrow);
					ioff = (*iv + *n - *k - 1) % mbv;
					i__1 = *k + ioff;
					kp = numroc_(&i__1, &mbv, &myrow, &ilastrow, &nprow);
					if (myrow == ilastrow) {
						kp -= ioff;
					}
					mydist = (myrow - ilastrow + nprow) % nprow;
					itop = mydist * mbv - ioff;
					/* Computing MIN */
					i__1 = itop + mbv;
					ibase = min(i__1,*k);
					/* Computing MIN */
					i__1 = max(0,itop);
					itop = min(i__1,*k);

L40:
					if (jj <= jjv + *k - 1) {
						height = ibase - itop;
						i__1 = itop - jj + jjv;
						dlaset_("All", &kp, &i__1, &dzero, &dzero, &work[ipw1 
								+ ii - iiv + (jj - jjv) * lw], &lw, (ftnlen)3)
							;
						dlaset_("Lower", &kp, &height, &dzero, &done, &work[
								ipw1 + ii - iiv + itop * lw], &lw, (ftnlen)5);
						/* Computing MAX */
						i__1 = 0, i__2 = kp - height;
						kp = max(i__1,i__2);
						ii += height;
						jj = jjv + ibase;
						mydist += nprow;
						itop = mydist * mbv - ioff;
						/* Computing MIN */
						i__1 = itop + mbv;
						ibase = min(i__1,*k);
						itop = min(itop,*k);
						goto L40;
					}
				}
			}

			i__1 = *n + iroffv;
			pbdtran_(&ictxt, "Columnwise", "Transpose", &i__1, k, &mbv, &work[
					ipw], &lw, &dzero, &work[ipv], &lv, &ivrow, &ivcol, &c_n1,
					&iccol, &work[ipt], (ftnlen)10, (ftnlen)9);

			/*           WORK( IPV ) = ( . V' ) -> WORK( IPV ) = V' is K x NQC */

			if (mycol == iccol) {
				ipv += icoffc * lv;
			}

			/*           WORK( IPW ) becomes MPC x K = C( IOFFC ) * V */
			/*           WORK( IPW ) = C( IOFFC ) * V  (MPC x NQC x K) -> MPC x K */

			lw = max(1,mpc);

			if (nqc > 0) 
			{
				dgemm_("No transpose", "Transpose", &mpc, k, &nqc, &done, &
						c__[ioffc], &ldc, &work[ipv], &lv, &dzero, &work[ipw],
						&lw, (ftnlen)12, (ftnlen)9);
			} else 
			{
				dlaset_("All", &mpc, k, &dzero, &dzero, &work[ipw], &lw, 
						(ftnlen)3);
			}

			dgsum2d_(&ictxt, "Rowwise", " ", &mpc, k, &work[ipw], &lw, &myrow,
					&ivcol, (ftnlen)7, (ftnlen)1);

			/*           WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T */

			if (mycol == ivcol) {
				if (myrow == ivrow) {

					/*                 Broadcast the block reflector to the other rows. */

					dtrbs2d_(&ictxt, "Columnwise", " ", uplo, "Non unit", k, 
							k, &t[1], &nbv, (ftnlen)10, (ftnlen)1, (ftnlen)1, 
							(ftnlen)8);
				} else {
					dtrbr2d_(&ictxt, "Columnwise", " ", uplo, "Non unit", k, 
							k, &t[1], &nbv, &ivrow, &mycol, (ftnlen)10, (
								ftnlen)1, (ftnlen)1, (ftnlen)8);
				}
				dtrmm_("Right", uplo, trans, "Non unit", &mpc, k, &done, &t[
						1], &nbv, &work[ipw], &lw, (ftnlen)5, (ftnlen)1, (
							ftnlen)1, (ftnlen)8);

				dgebs2d_(&ictxt, "Rowwise", " ", &mpc, k, &work[ipw], &lw, (
							ftnlen)7, (ftnlen)1);
			} else {
				dgebr2d_(&ictxt, "Rowwise", " ", &mpc, k, &work[ipw], &lw, &
						myrow, &ivcol, (ftnlen)7, (ftnlen)1);
			}

			/*               C            C      -     W       *     V' */
			/*           C( IOFFC ) = C( IOFFC ) - WORK( IPW ) * WORK( IPV ) */
			/*                        MPC x NQC    MPC x K         K x NQC */

			dgemm_("No transpose", "No transpose", &mpc, &nqc, k, &mone, &
					work[ipw], &lw, &work[ipv], &lv, &done, &c__[ioffc], &
					ldc, (ftnlen)12, (ftnlen)12);
		}
	} 
	else 
	{

		printf ("not supported\n");
		/*        V is stored rowwise */

		if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) 
		{

			/*           Form Q*sub( C ) or Q'*sub( C ) */

			/*           IROFFC = ICOFFV is required by the current transposition */
			/*           routine PBDTRAN */

			i__1 = *m + icoffv;
			mqv0 = numroc_(&i__1, &nbv, &mycol, &ivcol, &npcol);
			if (mycol == ivcol) {
				mqv = mqv0 - icoffv;
			} else {
				mqv = mqv0;
			}
			if (myrow == icrow) {
				mpc0 = mpc + iroffc;
			} else {
				mpc0 = mpc;
			}

			/*           Locally V( IOFFV ) is K x MQV, C( IOFFC ) is MPC x NQC */
			/*           WORK( IPV ) is MPC0 x K = [ . V( IOFFV ) ]' */
			/*           WORK( IPW ) is K x MQV0 = [ . V( IOFFV ) ] */
			/*           WORK( IPT ) is the workspace for PBDTRAN */

			ipv = 1;
			ipw = ipv + mpc0 * *k;
			ipt = ipw + *k * mqv0;
			lv = max(1,mpc0);
			lw = max(1,*k);

			if (myrow == ivrow) {
				if (mycol == ivcol) {
					dlaset_("All", k, &icoffv, &dzero, &dzero, &work[ipw], &
							lw, (ftnlen)3);
					ipw1 = ipw + icoffv * lw;
					dlacpy_("All", k, &mqv, &v[ioffv], &ldv, &work[ipw1], &lw,
							(ftnlen)3);
				} else {
					ipw1 = ipw;
					dlacpy_("All", k, &mqv, &v[ioffv], &ldv, &work[ipw1], &lw,
							(ftnlen)3);
				}

				if (forward) {

					/*                 WORK( IPW ) = ( . V1 V2 ) where V1 is unit upper */
					/*                 triangular, zeroes lower triangular part of V1 */

					mydist = (mycol - ivcol + npcol) % npcol;
					/* Computing MAX */
					i__1 = 0, i__2 = mydist * nbv - icoffv;
					ileft = max(i__1,i__2);
					jjbeg = jjv;
					jjend = jjv + mqv - 1;
					/* Computing MIN */
					i__1 = iceil_(&jjbeg, &nbv) * nbv;
					jjnxt = min(i__1,jjend);

L50:
					if (*k - ileft > 0) {
						i__1 = *k - ileft;
						i__2 = jjnxt - jjbeg + 1;
						dlaset_("Lower", &i__1, &i__2, &dzero, &done, &work[
								ipw1 + ileft + (jjbeg - jjv) * lw], &lw, (
									ftnlen)5);
						mydist += npcol;
						ileft = mydist * nbv - icoffv;
						jjbeg = jjnxt + 1;
						/* Computing MIN */
						i__1 = jjnxt + nbv;
						jjnxt = min(i__1,jjend);
						goto L50;
					}

				} else {

					/*                 WORK( IPW ) = ( . V1 V2 ) where V2 is unit lower */
					/*                 triangular, zeroes upper triangular part of V2. */

					ii = iiv;
					i__1 = *jv + *m - *k;
					infog1l_(&i__1, &nbv, &npcol, &mycol, &descv[8], &jj, &
							ilastcol);
					ioff = (*jv + *m - *k - 1) % nbv;
					i__1 = *k + ioff;
					kq = numroc_(&i__1, &nbv, &mycol, &ilastcol, &npcol);
					if (mycol == ilastcol) {
						kq -= ioff;
					}
					mydist = (mycol - ilastcol + npcol) % npcol;
					ileft = mydist * nbv - ioff;
					/* Computing MIN */
					i__1 = ileft + nbv;
					iright = min(i__1,*k);
					/* Computing MIN */
					i__1 = max(0,ileft);
					ileft = min(i__1,*k);

L60:
					if (ii <= iiv + *k - 1) {
						wide = iright - ileft;
						i__1 = ileft - ii + iiv;
						dlaset_("All", &i__1, &kq, &dzero, &dzero, &work[ipw1 
								+ ii - iiv + (jj - jjv) * lw], &lw, (ftnlen)3)
							;
						dlaset_("Upper", &wide, &kq, &dzero, &done, &work[
								ipw1 + ileft + (jj - jjv) * lw], &lw, (ftnlen)
								5);
						/* Computing MAX */
						i__1 = 0, i__2 = kq - wide;
						kq = max(i__1,i__2);
						ii = iiv + iright;
						jj += wide;
						mydist += npcol;
						ileft = mydist * nbv - ioff;
						/* Computing MIN */
						i__1 = ileft + nbv;
						iright = min(i__1,*k);
						ileft = min(ileft,*k);
						goto L60;
					}
				}
			}

			/*           WORK( IPV ) = WORK( IPW )' (replicated) is MPC0 x K */

			i__1 = *m + icoffv;
			pbdtran_(&ictxt, "Rowwise", "Transpose", k, &i__1, &nbv, &work[
					ipw], &lw, &dzero, &work[ipv], &lv, &ivrow, &ivcol, &
					icrow, &c_n1, &work[ipt], (ftnlen)7, (ftnlen)9);

			/*           WORK( IPV ) = ( . V )' -> WORK( IPV ) = V' is MPC x K */

			if (myrow == icrow) {
				ipv += iroffc;
			}

			/*           WORK( IPW ) becomes NQC x K = C( IOFFC )' * V' */
			/*           WORK( IPW ) = C( IOFFC )' * V'  (NQC x MPC x K) -> NQC x K */

			lw = max(1,nqc);

			if (mpc > 0) {
				dgemm_("Transpose", "No transpose", &nqc, k, &mpc, &done, &
						c__[ioffc], &ldc, &work[ipv], &lv, &dzero, &work[ipw],
						&lw, (ftnlen)9, (ftnlen)12);
			} else {
				dlaset_("All", &nqc, k, &dzero, &dzero, &work[ipw], &lw, (
							ftnlen)3);
			}

			dgsum2d_(&ictxt, "Columnwise", " ", &nqc, k, &work[ipw], &lw, &
					ivrow, &mycol, (ftnlen)10, (ftnlen)1);

			/*           WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T */

			if (myrow == ivrow) {
				if (mycol == ivcol) {

					/*                 Broadcast the block reflector to the other columns. */

					dtrbs2d_(&ictxt, "Rowwise", " ", uplo, "Non unit", k, k, &
							t[1], &mbv, (ftnlen)7, (ftnlen)1, (ftnlen)1, (
								ftnlen)8);
				} else {
					dtrbr2d_(&ictxt, "Rowwise", " ", uplo, "Non unit", k, k, &
							t[1], &mbv, &myrow, &ivcol, (ftnlen)7, (ftnlen)1, 
							(ftnlen)1, (ftnlen)8);
				}
				dtrmm_("Right", uplo, transt, "Non unit", &nqc, k, &done, &t[
						1], &mbv, &work[ipw], &lw, (ftnlen)5, (ftnlen)1, (
							ftnlen)1, (ftnlen)8);

				dgebs2d_(&ictxt, "Columnwise", " ", &nqc, k, &work[ipw], &lw, 
						(ftnlen)10, (ftnlen)1);
			} else {
				dgebr2d_(&ictxt, "Columnwise", " ", &nqc, k, &work[ipw], &lw, 
						&ivrow, &mycol, (ftnlen)10, (ftnlen)1);
			}

			/*               C            C      -     V'      *     W' */
			/*           C( IOFFC ) = C( IOFFC ) - WORK( IPV ) * WORK( IPW )' */
			/*                        MPC x NQC    MPC x K         K x NQC */

			dgemm_("No transpose", "Transpose", &mpc, &nqc, k, &mone, &work[
					ipv], &lv, &work[ipw], &lw, &done, &c__[ioffc], &ldc, (
						ftnlen)12, (ftnlen)9);

		} 
		else 
		{

			/*           Form Q*sub( C ) or Q'*sub( C ) */

			/*           Locally V( IOFFV ) is K x NQV, C( IOFFC ) is MPC x NQC */
			/*           WORK( IPV ) is K x NQV = V( IOFFV ), NQV = NQC */
			/*           WORK( IPW ) is MPC x K = C( IOFFC ) * V( IOFFV )' */

			ipv = 1;
			ipw = ipv + *k * nqc;
			lv = max(1,*k);
			lw = max(1,mpc);

			/*           Broadcast V to the other process rows. */

			pb_topget__(&ictxt, "Broadcast", "Columnwise", colbtop, (ftnlen)9,
					(ftnlen)10, (ftnlen)1);
			if (myrow == ivrow) {
				dgebs2d_(&ictxt, "Columnwise", colbtop, k, &nqc, &v[ioffv], &
						ldv, (ftnlen)10, (ftnlen)1);
				if (mycol == ivcol) {
					dtrbs2d_(&ictxt, "Columnwise", colbtop, uplo, "Non unit", 
							k, k, &t[1], &mbv, (ftnlen)10, (ftnlen)1, (ftnlen)
							1, (ftnlen)8);
				}
				dlacpy_("All", k, &nqc, &v[ioffv], &ldv, &work[ipv], &lv, (
							ftnlen)3);
			} else {
				dgebr2d_(&ictxt, "Columnwise", colbtop, k, &nqc, &work[ipv], &
						lv, &ivrow, &mycol, (ftnlen)10, (ftnlen)1);
				if (mycol == ivcol) {
					dtrbr2d_(&ictxt, "Columnwise", colbtop, uplo, "Non unit", 
							k, k, &t[1], &mbv, &ivrow, &mycol, (ftnlen)10, (
								ftnlen)1, (ftnlen)1, (ftnlen)8);
				}
			}

			if (forward) {

				/*              WORK(IPW) = ( V1 V2 ) where V1 is unit upper */
				/*              triangular, zeroes lower triangular part of V1 */

				mydist = (mycol - ivcol + npcol) % npcol;
				/* Computing MAX */
				i__1 = 0, i__2 = mydist * nbv - icoffv;
				ileft = max(i__1,i__2);
				jjbeg = jjv;
				jjend = jjv + nqc - 1;
				/* Computing MIN */
				i__1 = iceil_(&jjbeg, &nbv) * nbv;
				jjnxt = min(i__1,jjend);

L70:
				if (*k - ileft > 0) {
					i__1 = *k - ileft;
					i__2 = jjnxt - jjbeg + 1;
					dlaset_("Lower", &i__1, &i__2, &dzero, &done, &work[ipv 
							+ ileft + (jjbeg - jjv) * lv], &lv, (ftnlen)5);
					mydist += npcol;
					ileft = mydist * nbv - icoffv;
					jjbeg = jjnxt + 1;
					/* Computing MIN */
					i__1 = jjnxt + nbv;
					jjnxt = min(i__1,jjend);
					goto L70;
				}

			} else {

				/*              WORK( IPW ) = ( . V1 V2 ) where V2 is unit lower */
				/*              triangular, zeroes upper triangular part of V2. */

				ii = iiv;
				i__1 = *jv + *n - *k;
				infog1l_(&i__1, &nbv, &npcol, &mycol, &descv[8], &jj, &
						ilastcol);
				ioff = (*jv + *n - *k - 1) % nbv;
				i__1 = *k + ioff;
				kq = numroc_(&i__1, &nbv, &mycol, &ilastcol, &npcol);
				if (mycol == ilastcol) {
					kq -= ioff;
				}
				mydist = (mycol - ilastcol + npcol) % npcol;
				ileft = mydist * nbv - ioff;
				/* Computing MIN */
				i__1 = ileft + nbv;
				iright = min(i__1,*k);
				/* Computing MIN */
				i__1 = max(0,ileft);
				ileft = min(i__1,*k);

L80:
				if (ii <= iiv + *k - 1) {
					wide = iright - ileft;
					i__1 = ileft - ii + iiv;
					dlaset_("All", &i__1, &kq, &dzero, &dzero, &work[ipv + ii 
							- iiv + (jj - jjv) * lv], &lv, (ftnlen)3);
					dlaset_("Upper", &wide, &kq, &dzero, &done, &work[ipv + 
							ileft + (jj - jjv) * lv], &lv, (ftnlen)5);
					/* Computing MAX */
					i__1 = 0, i__2 = kq - wide;
					kq = max(i__1,i__2);
					ii = iiv + iright;
					jj += wide;
					mydist += npcol;
					ileft = mydist * nbv - ioff;
					/* Computing MIN */
					i__1 = ileft + nbv;
					iright = min(i__1,*k);
					ileft = min(ileft,*k);
					goto L80;
				}

			}

			/*           WORK( IPV ) is K x NQC = V = V( IOFFV ) */
			/*           WORK( IPW ) = C( IOFFC ) * V'  (MPC x NQC x K) -> MPC x K */

			if (nqc > 0) {
				dgemm_("No Transpose", "Transpose", &mpc, k, &nqc, &done, &
						c__[ioffc], &ldc, &work[ipv], &lv, &dzero, &work[ipw],
						&lw, (ftnlen)12, (ftnlen)9);
			} else {
				dlaset_("All", &mpc, k, &dzero, &dzero, &work[ipw], &lw, (
							ftnlen)3);
			}

			dgsum2d_(&ictxt, "Rowwise", " ", &mpc, k, &work[ipw], &lw, &myrow,
					&ivcol, (ftnlen)7, (ftnlen)1);

			/*           WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T */

			if (mycol == ivcol) {
				dtrmm_("Right", uplo, trans, "Non unit", &mpc, k, &done, &t[
						1], &mbv, &work[ipw], &lw, (ftnlen)5, (ftnlen)1, (
							ftnlen)1, (ftnlen)8);
				dgebs2d_(&ictxt, "Rowwise", " ", &mpc, k, &work[ipw], &lw, (
							ftnlen)7, (ftnlen)1);
			} else {
				dgebr2d_(&ictxt, "Rowwise", " ", &mpc, k, &work[ipw], &lw, &
						myrow, &ivcol, (ftnlen)7, (ftnlen)1);
			}

			/*               C            C      -     W       *     V */
			/*           C( IOFFC ) = C( IOFFC ) - WORK( IPW ) * WORK( IPV ) */
			/*                        MPC x NQC    MPC x K         K x NQC */

			dgemm_("No transpose", "No transpose", &mpc, &nqc, k, &mone, &
					work[ipw], &lw, &work[ipv], &lv, &done, &c__[ioffc], &
					ldc, (ftnlen)12, (ftnlen)12);

		}
	}


	cudaStreamDestroy(fstream);

	return 0;

	/*     End of PDLARFB */

} /* pdlarfb_ */

