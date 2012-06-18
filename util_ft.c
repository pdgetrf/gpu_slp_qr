#include "util_ft.h"
#include "slp.h"
#define NB 80
#define DTYPE_	0
#define CTXT_	1
#define M_		2
#define N_		3
#define MB_		4
#define NB_		5
#define RSRC_	6
#define CSRC_	7
#define LLD_	8
#define TAG_L	1021
#define TAG_LL	904	
#define TAG 99

extern int real_err_step;
extern int err_step;
extern int err_XX;	// myrow
extern int err_YY;	// mycol, both are 0-based
extern int cc;

///////////////////////////////////////////////
//    print out matrix for check			 //
///////////////////////////////////////////////
void printmatrix_ (double * A, int * descA, char * NAME, int * len)
{
	double * work2;
	int NOUT = 6;
	int ione = 1, izero = 0;
	int m = descA[2], n = descA[3];

	printf ("\n");

	work2 = (double *)malloc (m*m*sizeof (double));
	pdlaprnt_ (&m, &n, A, &ione, &ione, descA, &izero, &izero, NAME, &NOUT, work2, *len);
	free (work2);

	printf ("\n");

}

void fprintmatrix_ (double * A, int m, int n, int * descA, char * NAME, int len)
{
	double * work2;
	int NOUT = 6;
	int ione = 1, izero = 0;

	printf ("\n");

	work2 = (double *)malloc (m*m*sizeof (double));
	pdlaprnt_ (&m, &n, A, &ione, &ione, descA, &izero, &izero, NAME, &NOUT, work2, len);
	free (work2);

	printf ("\n");

}

double verifyQR(double *Aorg, double * A, int M, int N, int * descA, double * tau, t_Grid *grid)
{

	double resid;
	double *work=NULL;
	int ione = 1, izero = 0;;
	double done=1.0, mone=-1.0;

	//----- non-ft -----//
	int Mp0, Nq0, iarow, iacol, lwork;

	// grid parameters
	int	nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	//----------- allocate workspace ----------// -_-||

	//LWORK = NB_A * ( 2*Mp0 + Nq0 + NB_A ), where
	//Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ) * NB_A,
	//Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ) * MB_A,
	//IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
	//IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),NPROW ),
	//IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),NPCOL )

	iarow = indxg2p_ (&ione, &nb, &myrow, &izero, &nprow);
	iacol = indxg2p_ (&ione, &nb, &mycol, &izero, &npcol);

	Mp0 = numroc_(&M , &nb, &myrow, &iarow, &nprow )*nb;
	Nq0 = numroc_(&N , &nb, &mycol, &iacol, &npcol )*nb;

	lwork = nb*(2*Mp0 + Nq0 + nb);
	work = (double *)malloc(lwork*sizeof(double));

	pdgeqrrv_ (&M, &N, A, &ione, &ione, descA, tau, work);

	pdmatadd_ ( &M, &N, &done, Aorg, &ione, &ione, descA, &mone, A, &ione, &ione, descA);
	resid = pdlange_("F", &M, &N, A, &ione, &ione, descA, work)/pdlange_("F", &M, &N, Aorg, &ione, &ione, descA, work)/M;

	free (work);

	return resid;
}




/*
 * produce distributed matrix,  
 */
void distr_matrix (bool fill,	double **A, int *descA, 
		int M, int N, t_Grid *grid,
		int *np_A, int *nq_A) 
{
	int info;
	int izero = 0;

	// grid parameters
	int	nb = grid->nb;
	int	ictxt = grid->ictxt;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	// allocate the generator matrix and check matrix
	int np_iA = numroc_( &M, &nb, &myrow, &izero, &nprow );
	int nq_iA = numroc_( &N, &nb, &mycol, &izero, &npcol );

	if (np_iA*nq_iA!=0)
	{
		*A = (double *)malloc(np_iA*nq_iA*sizeof(double)) ;
		if (*A == NULL) ABORT;
		memset (*A, 0, np_iA*nq_iA*sizeof(double));
	}

	if (descA != NULL)
	{
		int itemp = MAX( 1, np_iA );
		descinit_( descA, &M, &N, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
		if (info != 0) ABORT;
	}

	if (fill)
	{
		// fill in random numbers
		int seed = 800;
		pdmatgen_ (&ictxt, "N", "N", &M, &N, &nb, &nb, *A, 
				descA+8, descA+6, descA+7, 
				&seed, &izero, &np_iA, &izero, &nq_iA, 
				&myrow, &mycol, &nprow, &npcol);
	}

	/* set np and nq */
	if (np_A != NULL)
		*np_A = np_iA;
	if (nq_A != NULL)
		*nq_A = nq_iA;
}


void set_grid (t_Grid *grid, int ictxt, int nprocs_mpi, int myrank_mpi, int myrow, int mycol, int nprow, int npcol, int nb)
{
	grid->ictxt = ictxt;
	grid->myrank_mpi = myrank_mpi;
	grid->nprocs_mpi = nprocs_mpi;
	grid->myrow = myrow;
	grid->mycol = mycol;
	grid->nprow = nprow;
	grid->npcol = npcol;
	grid->nb = nb;
}

