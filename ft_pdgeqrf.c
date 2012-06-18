/* pdgeqrf.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "util_ft.h"
#include "f2c.h"

extern int err_step;
extern int err_block;
extern int real_err_step;
extern int err_XX;	// myrow
extern int err_YY;	// mycol, both are 0-based

extern int numroc_(int * N, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );
extern void descset_ (int * desc, int *M, int *N, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld);

extern int cc;

extern int ft_pdlarfb(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, doublereal *v, int *
	iv, int *jv, int *descv, doublereal *t, doublereal *c__, 
	int *ic, int *jc, int *descc, doublereal *work, ftnlen 
	side_len, ftnlen trans_len, ftnlen direct_len, ftnlen storev_len, t_checksuite *cs, t_Grid *grid, int *ReEntry);

/* Table of constant values */

static int c__1 = 1;
static int c__2 = 2;
static int c__6 = 6;

/* Subroutine */ int ft_pdgeqrf(int *m, int *n, doublereal *a, int *
	ia, int *ja, int *descA, doublereal *tau, doublereal *work, 
	int *lwork, int *info, t_Grid *grid, t_checksuite *cs, int ReEntry)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static int i__, j, k;
    extern /* Subroutine */ int pb_topget__(int *, char *, char *, char *,
	     ftnlen, ftnlen, ftnlen), pb_topset__(int *, char *, char *, 
	    char *, ftnlen, ftnlen, ftnlen);
    static int jb, jn, mp0, nq0, ipw, idum1[1], idum2[1];
    extern int iceil_(int *, int *);
    static int icoff, iacol, iinfo, npcol, iarow, mycol, lwmin, ictxt, 
	    nprow, myrow;
    extern /* Subroutine */ int blacs_gridinfo__(int *, int *, 
	    int *, int *, int *);
    extern int numroc_(int *, int *, int *, int *, 
	    int *);
    static logical lquery;
    extern /* Subroutine */ int chk1mat_(int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *),
	     pdgeqr2_(int *, int *, doublereal *, int *, int *
	    , int *, doublereal *, doublereal *, int *, int *);
    extern int indxg2p_(int *, int *, int *, int *, 
	    int *);
    extern /* Subroutine */ int pdlarfb_(char *, char *, char *, char *, 
	    int *, int *, int *, doublereal *, int *, int 
	    *, int *, doublereal *, doublereal *, int *, int *, 
	    int *, doublereal *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    pdlarft_(char *, char *, int *, int *, doublereal *, 
	    int *, int *, int *, doublereal *, doublereal *, 
	    doublereal *, ftnlen, ftnlen), pxerbla_(int *, char *, 
	    int *, ftnlen);
    static char colbtop[1], rowbtop[1];
    extern /* Subroutine */ int pchk1mat_(int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *);


/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 25, 2001 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  PDGEQRF computes a QR factorization of a real distributed M-by-N */
/*  matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = Q * R. */

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

/*  M       (global input) INTEGER */
/*          The number of rows to be operated on, i.e. the number of rows */
/*          of the distributed submatrix sub( A ). M >= 0. */

/*  N       (global input) INTEGER */
/*          The number of columns to be operated on, i.e. the number of */
/*          columns of the distributed submatrix sub( A ). N >= 0. */

/*  A       (local input/local output) DOUBLE PRECISION pointer into the */
/*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)). */
/*          On entry, the local pieces of the M-by-N distributed matrix */
/*          sub( A ) which is to be factored.  On exit, the elements on */
/*          and above the diagonal of sub( A ) contain the min(M,N) by N */
/*          upper trapezoidal matrix R (R is upper triangular if M >= N); */
/*          the elements below the diagonal, with the array TAU, */
/*          represent the orthogonal matrix Q as a product of elementary */
/*          reflectors (see Further Details). */

/*  IA      (global input) INTEGER */
/*          The row index in the global array A indicating the first */
/*          row of sub( A ). */

/*  JA      (global input) INTEGER */
/*          The column index in the global array A indicating the */
/*          first column of sub( A ). */

/*  DESCA   (global and local input) INTEGER array of dimension DLEN_. */
/*          The array descriptor for the distributed matrix A. */

/*  TAU     (local output) DOUBLE PRECISION array, dimension */
/*          LOCc(JA+MIN(M,N)-1). This array contains the scalar factors */
/*          TAU of the elementary reflectors. TAU is tied to the */
/*          distributed matrix A. */

/*  WORK    (local workspace/local output) DOUBLE PRECISION array, */
/*                                                     dimension (LWORK) */
/*          On exit, WORK(1) returns the minimal and optimal LWORK. */

/*  LWORK   (local or global input) INTEGER */
/*          The dimension of the array WORK. */
/*          LWORK is local input and must be at least */
/*          LWORK >= NB_A * ( Mp0 + Nq0 + NB_A ), where */

/*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ), */
/*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ), */
/*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ), */
/*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ), */
/*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ), */

/*          and NUMROC, INDXG2P are ScaLAPACK tool functions; */
/*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling */
/*          the subroutine BLACS_GRIDINFO. */

/*          If LWORK = -1, then LWORK is global input and a workspace */
/*          query is assumed; the routine only calculates the minimum */
/*          and optimal size for all work arrays. Each of these */
/*          values is returned in the first entry of the corresponding */
/*          work array, and no error message is issued by PXERBLA. */

/*  INFO    (global output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  If the i-th argument is an array and the j-entry had */
/*                an illegal value, then INFO = -(i*100+j), if the i-th */
/*                argument is a scalar and had an illegal value, then */
/*                INFO = -i. */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(ja) H(ja+1) . . . H(ja+k-1), where k = min(m,n). */

/*  Each H(i) has the form */

/*     H(j) = I - tau * v * v' */

/*  where tau is a real scalar, and v is a real vector with v(1:i-1) = 0 */
/*  and v(i) = 1; v(i+1:m) is stored on exit in A(ia+i:ia+m-1,ja+i-1), */
/*  and tau in TAU(ja+i-1). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Get grid parameters */

    /* Parameter adjustments */
    --work;
    --tau;
    --descA;
    --a;

	int ShouldRecover = ReEntry;

    /* Function Body */
    ictxt = descA[2];
    blacs_gridinfo__(&ictxt, &nprow, &npcol, &myrow, &mycol);

/*     Test the input parameters */

    *info = 0;
    if (nprow == -1) {
	*info = -602;
    } else {
	chk1mat_(m, &c__1, n, &c__2, ia, ja, &descA[1], &c__6, info);
	if (*info == 0) {
	    icoff = (*ja - 1) % descA[6];
	    iarow = indxg2p_(ia, &descA[5], &myrow, &descA[7], &nprow);
	    iacol = indxg2p_(ja, &descA[6], &mycol, &descA[8], &npcol);
	    i__1 = *m + (*ia - 1) % descA[5];
	    mp0 = numroc_(&i__1, &descA[5], &myrow, &iarow, &nprow);
	    i__1 = *n + icoff;
	    nq0 = numroc_(&i__1, &descA[6], &mycol, &iacol, &npcol);
	    lwmin = descA[6] * (mp0 + nq0 + descA[6]);

		if (!ReEntry)
			work[1] = (doublereal) lwmin;

	    lquery = *lwork == -1;
	    if (*lwork < lwmin && ! lquery) {
		*info = -9;
	    }
	}
	if (*lwork == -1) {
	    idum1[0] = -1;
	} else {
	    idum1[0] = 1;
	}
	idum2[0] = 9;
	pchk1mat_(m, &c__1, n, &c__2, ia, ja, &descA[1], &c__6, &c__1, idum1, 
		idum2, info);
    }

    if (*info != 0) {
	i__1 = -(*info);
	pxerbla_(&ictxt, "PDGEQRF", &i__1, (ftnlen)7);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    k = min(*m,*n);
    ipw = descA[6] * descA[6] + 1;
    pb_topget__(&ictxt, "Broadcast", "Rowwise", rowbtop, (ftnlen)9, (ftnlen)7,
	     (ftnlen)1);
    pb_topget__(&ictxt, "Broadcast", "Columnwise", colbtop, (ftnlen)9, (
	    ftnlen)10, (ftnlen)1);
    pb_topset__(&ictxt, "Broadcast", "Rowwise", "I-ring", (ftnlen)9, (ftnlen)
	    7, (ftnlen)6);
    pb_topset__(&ictxt, "Broadcast", "Columnwise", " ", (ftnlen)9, (ftnlen)10,
	     (ftnlen)1);
	
	/*	   
	 *	   Initial Checksum of [A] -> [A, C, C]
	 *				      [R      ]	<--- both uninitialized 
	 *				      [R      ] <--+
	 */
	int realn = *n-(cs->Nc);
	int N_1 = realn+1;
	cc=0;
	int mn = min(*m,*n);
	/*
	{
		int one = 1;
		printmatrix_ (&a[1], &descA[1], "A", &one);
	}
	*/
#ifdef CHECKPOINT
	if (!ReEntry)
		checkpointing ('r', &a[1], *m, realn, 1, 1, &descA[1], &a[1], 1, N_1, &descA[1], cs, grid);
#endif

/*     Handle the first block of columns separately */

/* Computing MIN */
    i__1 = iceil_(ja, &descA[6]) * descA[6], i__2 = *ja + k - 1;
    jn = min(i__1,i__2);
    jb = jn - *ja + 1;
	
	/* make a local copy & checkpoint localcopy */
	if (!ReEntry)
	{
		local_copy (&a[1], &descA[1], 1, 1, cs, grid);
	}

/*     Compute the QR factorization of the first block A(ia:ia+m-1,ja:jn) */

	if (!ReEntry)
		pdgeqr2_(m, &jb, &a[1], ia, ja, &descA[1], &tau[1], &work[1], lwork, &iinfo);

    if (*ja + jb <= *ja + *n - 1) 
	{

		/*        Form the triangular factor of the block reflector */
		/*        H = H(ja) H(ja+1) . . . H(jn) */

		if (!ReEntry)
			pdlarft_("Forward", "Columnwise", m, &jb, &a[1], ia, ja, &descA[1], &
					tau[1], &work[1], &work[ipw], (ftnlen)7, (ftnlen)10);

		/*        Apply H' to A(ia:ia+m-1,ja+jb:ja+n-1) from the left */

		i__1 = *n - jb;
		i__2 = *ja + jb;

		ft_pdlarfb("Left", "Transpose", "Forward", "Columnwise", m, &i__1, &jb, 
				&a[1], ia, ja, &descA[1], &work[1], &a[1], ia, &i__2, &descA[1], 
				&work[ipw], (ftnlen)4, (ftnlen)9, (ftnlen)7, (ftnlen)10, cs, grid, &ReEntry);
    }

/*     Loop over the remaining blocks of columns */

    i__1 = *ja + k - 1;
    i__2 = descA[6];
	
	int jj = 1; 
	int dec = 0;
	int jc=1;
    for (j = jn + 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2, jj++) 
	{
		if (jj%(grid->npcol)==0)
		{
			if (!ReEntry)
			{
				local_copy (&a[1], &descA[1], j, j, cs, grid);
				/*
				if (j==31)
				{
					int one = 1;
					printmatrix_ (cs->localcopy, cs->descL, "B", &one);
				}
				*/
			}
			jc=j;
		}

		/* Computing MIN */
		i__3 = k - j + *ja;
		jb = min(i__3,descA[6]);
		i__ = *ia + j - *ja;

		/*        Compute the QR factorization of the current block */
		/*        A(i:ia+m-1,j:j+jb-1) */

		i__3 = *m - j + *ja;
		if (!ReEntry)
			pdgeqr2_(&i__3, &jb, &a[1], &i__, &j, &descA[1], &tau[1], &work[1], 
					lwork, &iinfo);

		/* 
		 *  Checkpointing L		  
		 */
#ifdef CHECKPOINT
		if (((jj+1)%(grid->npcol))==0 || j+i__2>realn)
		{
			/*
			if (grid->myrank_mpi==0)
				printf ("j=%d, n=%d, jx=%d\n", j, mn-j+*ja+2*jb, j-(grid->npcol-1)*jb);
				*/
			if (!ReEntry)
				checkpointing ('d', &a[1], mn-j+*ja+2*jb, realn, jc, jc, &descA[1], &a[1], *m+1, j, &descA[1], cs, grid);
		}
#endif

		dec += (((jj+1)%(grid->npcol))==0 || j+i__2>realn)?(jb*2):0;
		if (j + jb <= *ja + *n - 1 - dec) 
		{

			/*           Form the triangular factor of the block reflector */
			/*           H = H(j) H(j+1) . . . H(j+jb-1) */

			i__3 = *m - j + *ja;

			if (!ReEntry)
				pdlarft_("Forward", "Columnwise", &i__3, &jb, &a[1], &i__, &j, &
						descA[1], &tau[1], &work[1], &work[ipw], (ftnlen)7, (
							ftnlen)10);

			/*           Apply H' to A(i:ia+m-1,j+jb:ja+n-1) from the left */

			i__3 = *m - j + *ja;
			i__4 = *n - j - jb + *ja - dec;
			i__5 = j + jb;

			int ret=0;
			ret = ft_pdlarfb("Left", "Transpose", "Forward", "Columnwise", &i__3, 
					&i__4, &jb, &a[1], &i__, &j, &descA[1], &work[1], &a[1], 
					&i__, &i__5, &descA[1], &work[ipw], (ftnlen)4, (ftnlen)9, (ftnlen)7, 
					(ftnlen)10, cs, grid, &ReEntry);

			if (ret==-1)
			{
				return -1;
			}
		}
		
		/*
		 * inject fault & recover checksum and zone 1 ----------------------------------------------------------
		 */
		//#define INJECT
		#ifdef INJECT
		int nb=grid->nb;
		real_err_step=err_block*nb*npcol+nb*err_step+1;
		
		if (j==real_err_step)
		//if ((j-1)/i__2==real_err_step)
		{
			/*
			{
				int one = 1;
				printmatrix_ (&a[1], &descA[1], "A", &one);
			}
			*/

			MPI_Barrier (MPI_COMM_WORLD);

			if (!ShouldRecover)
			{
				/*
				{
					int one = 1;
					printmatrix_ (&a[1], &descA[1], "A", &one);
				}
				*/

				// kill proce
				/*
				if (grid->myrow==err_XX && grid->mycol==err_YY)
				{
					kill_proc (&a[1], &descA[1], grid);

					memset (cs->localcopy, 0, (cs->np_L)*(cs->nq_L)*sizeof(double));	
				}
				else
				{
					dump_to_disk (&a[1], &descA[1], tau, work, *lwork, grid, cs);
				}

				MPI_Barrier (MPI_COMM_WORLD);
				return -1;
				*/
			}
			else
			{
				//load_from_disk (&a[1], &descA[1], tau, work, *lwork, grid, cs);
				//MPI_Barrier (MPI_COMM_WORLD);

				// recover
				recover (CHECKSUM, ZONE1, ZONE2, ZONE3, &a[1], *m, realn, &descA[1], err_XX, err_YY, j, grid, cs);
				/*
				{
					int one = 1;
					printmatrix_ (&a[1], &descA[1], "B", &one);
				}
				*/

				/*
				{
					int one = 1;
					//printmatrix_ (cs->localcopy, cs->descL, "B", &one);
					printmatrix_ (&a[1], &descA[1], "B", &one);
				}
				*/
			}


			MPI_Barrier (MPI_COMM_WORLD);
		}
		#endif //----------------------------------------------------------

		/* L10: */
	}

#ifdef VERIFY_CHK
	verify_checkpointing (*m, *n, &a[1], &descA[1], grid, cs);
#endif

    pb_topset__(&ictxt, "Broadcast", "Rowwise", rowbtop, (ftnlen)9, (ftnlen)7,
	     (ftnlen)1);
    pb_topset__(&ictxt, "Broadcast", "Columnwise", colbtop, (ftnlen)9, (
	    ftnlen)10, (ftnlen)1);

    work[1] = (doublereal) lwmin;

    return 0;

/*     End of PDGEQRF */

} /* pdgeqrf_ */

