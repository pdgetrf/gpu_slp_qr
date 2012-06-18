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
extern int ifcleanfile;

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

void set_FTKit (double *A, t_checksuite *cs, double *tau, double *work, int lwork)
{
	cs->kit.A = A;
	cs->kit.localcopy = cs->localcopy;
	cs->kit.tau = tau;
	cs->kit.work = work;
	cs->kit.lwork = lwork;
}

void load_from_disk (double *A, int *descA, double *tau, double *work, int lwork, t_Grid *grid, t_checksuite *cs)
{
	int i;
	// grid parameters
	int	nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	int izero=0;

	int np_A = numroc_( &descA[M_], &nb, &myrow, &izero, &nprow );
	int nq_A = numroc_( &descA[N_], &nb, &mycol, &izero, &npcol );

	FILE * pFile;
	char filename[256];
	sprintf (filename, "/ssd/dump_%d", grid->myrank_mpi);
	pFile = fopen ( filename , "rb" );
	
	if (pFile)
	{
		//printf ("(%d,%d): loading from %s\n", myrow, mycol, filename);
		// load A
		fread (A, sizeof(double), np_A*nq_A, pFile );

		// read local copy
		fread (cs->localcopy, sizeof(double), (cs->np_L)*(cs->nq_L), pFile );

		// load tau 
		fread (tau, sizeof(double), cs->MaxLocalCol, pFile );
	
		// load work 
		fread (work, sizeof(double), lwork, pFile );

		fclose (pFile);

		ifcleanfile = 1;
	}
}

void dump_to_disk (double *A, int *descA, double *tau, double *work, int lwork, t_Grid *grid, t_checksuite *cs)
{
	int i;
	// grid parameters
	int	nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	int izero=0;

	int np_A = numroc_( &descA[M_], &nb, &myrow, &izero, &nprow );
	int nq_A = numroc_( &descA[N_], &nb, &mycol, &izero, &npcol );
	
	FILE * pFile;
	char filename[256];
	sprintf (filename, "/ssd/dump_%d", grid->myrank_mpi);
	//printf ("(%d,%d): dumping to %s\n", myrow, mycol, filename);
	pFile = fopen ( filename , "wb" );
	if (pFile == NULL)
	{
		printf ("can't open dumpfile %s\n", filename);
		exit(-1);
	}
	
	// dump A
	fwrite (A, sizeof(double), np_A*nq_A, pFile );

	// dump local copy
	fwrite (cs->localcopy, sizeof(double), (cs->np_L)*(cs->nq_L), pFile );

	// dump tau 
	fwrite (tau, sizeof(double), cs->MaxLocalCol, pFile );

	// dump work 
	fwrite (work, sizeof(double), lwork, pFile );

	fclose (pFile);
}

void inject_errors (double *C, int *descC, int *err_XX, int *err_YY, t_Grid *grid)
{
	int i;
	// grid parameters
	int	nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	int izero=0;
	int np_C = numroc_( &descC[M_], &nb, &myrow, &izero, &nprow );
	int nq_C = numroc_( &descC[N_], &nb, &mycol, &izero, &npcol );

	if (myrow == err_XX[0] && mycol == err_YY[0])
	{
		for (i=0; i<np_C*nq_C; i++)
			C[i]++;
		printf ("proc (%d,%d) is dead\n", myrow, mycol);
	}
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



#ifdef old_verify
double verifyLU(double * A, int M, int * descA, int * ipiv, double *Acpy, t_Grid *grid)
{
	int ione=1;
	double done = 1.0, mone=-1.0e0;
	double * work = NULL;

	pdmatadd_ ( &M, &M, &done, Acpy, &ione, &ione, descA, &mone, A, &ione, &ione, descA);
	return pdlange_("F", &M, &M, A, &ione, &ione, descA, work)/pdlange_("F", &M, &M, Acpy, &ione, &ione, descA, work)/M;
}
#endif

#ifdef REASSEMBLE	//<----- not being used
double verifyLU(double *A, double * A, int * descA, int * ipiv, double *Acpy, t_Grid *grid)
{
	int m, n, m_nb, n_nb;
	int ione=1, izero = 0;
	double done = 1.0, mone=(-1.0e0), pone=(1.0e0), dzero = 0;
	double * work = NULL;
	int i;
	int Aseed = 100;
	double resid;
	int Mp0, Nq0, iarow, iacol, lwork;

	// grid parameters
	m = 
	int	nb = grid->nb;
	int	ictxt = grid->ictxt;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	//----------- allocate workspace ----------// -_-||

	//*          LWORK >= Mp0 * NB_A + Nq0 * MB_A, where
	//*
	//      *          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
	//      *          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
	//      *          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
	//      *          Mp0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
	//      *          Nq0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),

	iarow = indxg2p_ (&ione, &nb, &myrow, &izero, &nprow);
	iacol = indxg2p_ (&ione, &nb, &mycol, &izero, &npcol);

	Mp0 = numroc_(&m , &nb, &myrow, &iarow, &nprow );
	Nq0 = numroc_(&n , &nb, &mycol, &iacol, &npcol );

	lwork = Mp0 * nb + Nq0 * nb;
	work = (double *)malloc(lwork*sizeof(double));

	//----------- A - A' -----------------//
	pdgetrrv_ (&m, &n, A, &ione, &ione, descA, ipiv, work );

	pdmatadd_ ( &m, &n, &done, Acpy, &ione, &ione, descA, &mone, A, &ione, &ione, descA);
	resid = pdlange_("F", &m, &n, A, &ione, &ione, descA, work);

	return resid;
}
#endif

/*
 * recovery the broken panel factorizations
 */
void recover_zone2 (double *A, int MA, int NA, int *descA, int errx, int erry, int errstep, int vic_id, t_Grid *grid, t_checksuite *cs)
{
	int izero=0;
	int i, j, k;
	int lrindx, lcindx, rsrc, csrc;
	MPI_Comm part_comm;	// participating communicator 
	int lindx, rocsrc, nchm[2];
	int gi, gj;
	double *Aoff, *CC, *T;
	Aoff = CC = T = NULL;
	double *Achk_old, *Achk_new;
	int descAchk[9], npp, nqq;

	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	int iamvic = (grid->myrank_mpi==vic_id);

	MPI_Group world_group, part_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	// determine local data size (non-checksum)
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );

	int NAA = NA+cs->Nc;
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
	int nqq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );

	// roll back to the start status of this Q-panel factorization
	int jj = errstep/(npcol*nb)*npcol*nb+1;
	/*
	if (grid->myrank_mpi==0)
		printf ("resetting panel starting (%d,%d)\n", jj, jj);
		*/
	panel_reset (A, descA, jj, jj, cs, grid);

	// re-factorize  
	int MM = MA-jj+1;
	int NN = npcol*nb; 
	
	/*
	if (grid->myrank_mpi==0)
		printf ("re-factorizing panels starting (%d,%d) to (%d,%d) of size %d x %d\n", jj, jj, errstep, errstep, MM, NN);
		*/

	int info;
	int lwork = -1;
	double lazywork;
	ft_pdgeqrf0 (&MM, &NN, A, &jj, &jj, descA, NULL, &lazywork, &lwork, &info, grid, cs, errstep);
	lwork = (int)lazywork;
	double *work = (double*)malloc(lwork*sizeof(double));

	int nchkc = numroc_( &NN, &nb, &mycol, &izero, &npcol ); //LOCc(N_A) 
	double *tau = (double*)malloc(NA*sizeof(double));

	ft_pdgeqrf0 (&MM, &NN, A, &jj, &jj, descA, tau, work, &lwork, &info, grid, cs, errstep);
	free (tau);
	free (work);
}

void recover_local_fillin (double *A, int MA, int NA, int *descA, int errx, int erry, int errstep, int vic_id, t_Grid *grid, t_checksuite *cs)
{
	int izero=0;
	int i, j, k;
	int lrindx, lcindx, rsrc, csrc;
	MPI_Comm part_comm;	// participating communicator 
	int lindx, rocsrc, nchm[2];
	int gi, gj;
	double *Aoff, *CC, *T;
	Aoff = CC = T = NULL;
	double *Achk_old, *Achk_new;
	int descAchk[9], npp, nqq;

	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	int iamvic = (grid->myrank_mpi==vic_id);

	MPI_Group world_group, part_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	// determine local data size (non-checksum)
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );

	int NAA = NA+nb*2;
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
	int nqq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );

	distr_matrix (false, &Achk_old, descAchk, MA, nb*2, grid, &npp, &nqq);
	distr_matrix (false, &Achk_new, descAchk, MA, nb*2, grid, &npp, &nqq);

	// copy the old checksum into buffer
	Cpdgemr2d (MA, nb*2, A, 1, NA+1, descA, Achk_old, 1, 1, descAchk, grid->ictxt);

	// re-calculate the checksum
	hor_checkpointing_local (A, descA, 1, 1, 1, cs, grid);

	// copy the new checksum into buffer
	Cpdgemr2d (MA, nb*2, A, 1, NA+1, descA, Achk_new, 1, 1, descAchk, grid->ictxt);

	int N_M = nb*2, ione=1;
	double done = 1, mone = -1;
	pdmatadd_ ( &MA, &N_M, &done, Achk_old, &ione, &ione, descAchk, &mone, Achk_new, &ione, &ione, descAchk);
	Cpdgemr2d (MA, nb*2, Achk_new, 1, 1, descAchk, A, 1, NA+1, descA, grid->ictxt);

	// copy resut back in place
	int i_c, j_c;
	if (iamvic)
		j_c = nq_A;
	MPI_Bcast( &j_c, 1, MPI_INT, vic_id, MPI_COMM_WORLD );

	for (j=0; j<j_c; j+=nb)
	{
		int gi, gj, g;
		if (iamvic)
		{
			// source matrix column (global)
			int j_1 = j + 1;
			g = indxl2g_( &j_1, &nb, &mycol, &izero, &npcol );

			gj = NA+1;

			// source matrix row (global)
			int ione = 1;
			gi = indxl2g_( &ione, &nb, &myrow, &izero, &nprow );
		}
		MPI_Bcast( &gi, 1, MPI_INT, vic_id, MPI_COMM_WORLD );
		MPI_Bcast( &gj, 1, MPI_INT, vic_id, MPI_COMM_WORLD );
			
		infog2l_ (&gi, &gj, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
		//printf ("(%d,%d) <--- (%d,%d)\n", rsrc, vic_id%npcol, rsrc, csrc);

		// copy the new checksum into buffer
		if (csrc==(vic_id%npcol))
		{
			// local copy
			if (iamvic)
			{
				//printf ("(%d,%d) local copy when j=%d\n", myrow, mycol, j);
				memcpy (A+j*np_A, A+(lcindx-1)*np_A, np_A*nb*sizeof(double));
			}
		}
		else
		{
			// MPI copy
			if (myrow == rsrc)
			{
				MPI_Status status;
				if (iamvic)
				{
					//printf ("(%d,%d) reciving from (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, csrc, j);
					MPI_Recv(A+j*np_A, np_A*nb, MPI_DOUBLE, csrc, 1, cs->rowcomm, &status);
				}
				else if (mycol==csrc)
				{
					//printf ("(%d,%d) sending to (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, vic_id%npcol, j);
					MPI_Send(A+(lcindx-1)*np_A, np_A*nb, MPI_DOUBLE, vic_id%npcol, 1, cs->rowcomm);
				}
			}
		}
	}

	// copy the correct checksum back //
	Cpdgemr2d (MA, nb*2, Achk_old, 1, 1, descAchk, A, 1, NA+1, descA, grid->ictxt);

	if (npp*nqq>0)
	{
		free (Achk_new);
		free (Achk_old);
	}
}

void recover_zone1 (double *A, int MA, int NA, int *descA, int errx, int erry, int errstep, int vic_id, t_Grid *grid, t_checksuite *cs)
{
	int izero=0;
	int i, j, k;
	int lrindx, lcindx, rsrc, csrc;
	MPI_Comm part_comm;	// participating communicator 
	int lindx, rocsrc, nchm[2];
	int gi, gj;
	double *Aoff, *CC, *T;
	Aoff = CC = T = NULL;
	double *Achk_old, *Achk_new;
	int descAchk[9], npp, nqq;

	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	int iamvic = (grid->myrank_mpi==vic_id);

	MPI_Group world_group, part_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	// determine local data size (non-checksum)
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );

	int NAA = NA+cs->Nc;
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
	int nqq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );

	distr_matrix (false, &Achk_old, descAchk, MA, NAA-MA, grid, &npp, &nqq);
	distr_matrix (false, &Achk_new, descAchk, MA, NAA-MA, grid, &npp, &nqq);

	// copy the old checksum into buffer
	Cpdgemr2d (MA, NAA-MA, A, 1, MA+1, descA, Achk_old, 1, 1, descAchk, grid->ictxt);

	// re-calculate the checksum
	checkpointing('r', A, MA, MA, 1, 1, descA, A, 1, MA, descA, cs, grid);
	//checkpointing_short ('r', A, MA, MA, 1, 1, descA, A, 1, MA, descA, cs, grid, errstep);

	// copy the new checksum into buffer
	Cpdgemr2d (MA, NAA-MA, A, 1, MA+1, descA, Achk_new, 1, 1, descAchk, grid->ictxt);

	int N_M = NAA-MA, ione=1;
	double done = 1, mone = -1;
	pdmatadd_ ( &MA, &N_M, &done, Achk_old, &ione, &ione, descAchk, &mone, Achk_new, &ione, &ione, descAchk);
	Cpdgemr2d (MA, NAA-MA, Achk_new, 1, 1, descAchk, A, 1, NA+1, descA, grid->ictxt);

	// copy resut back in place
	int i_c, j_c;
	if (iamvic)
		j_c = nq_A;
	MPI_Bcast( &j_c, 1, MPI_INT, vic_id, MPI_COMM_WORLD );

	for (j=0; j<j_c; j+=nb)
	{
		int gi, gj, g;
		if (iamvic)
		{
			// source matrix column (global)
			int j_1 = j + 1;
			g = indxl2g_( &j_1, &nb, &mycol, &izero, &npcol );

			gj = NAA-((int)(g/(npcol*nb))+1)*nb*2+1;

			// source matrix row (global)
			int ione = 1;
			gi = indxl2g_( &ione, &nb, &myrow, &izero, &nprow );
		}
		MPI_Bcast( &gi, 1, MPI_INT, vic_id, MPI_COMM_WORLD );
		MPI_Bcast( &gj, 1, MPI_INT, vic_id, MPI_COMM_WORLD );
			
		infog2l_ (&gi, &gj, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
		//printf ("(%d,%d) <--- (%d,%d)\n", rsrc, vic_id%npcol, rsrc, csrc);

		// copy the new checksum into buffer
		if (csrc==(vic_id%npcol))
		{
			// local copy
			if (iamvic)
			{
				//printf ("(%d,%d) local copy when j=%d\n", myrow, mycol, j);
				memcpy (A+j*np_A, A+(lcindx-1)*np_A, np_A*nb*sizeof(double));
			}
		}
		else
		{
			// MPI copy
			if (myrow == rsrc)
			{
				MPI_Status status;
				if (iamvic)
				{
					//printf ("(%d,%d) reciving from (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, csrc, j);
					MPI_Recv(A+j*np_A, np_A*nb, MPI_DOUBLE, csrc, 1, cs->rowcomm, &status);
				}
				else if (mycol==csrc)
				{
					//printf ("(%d,%d) sending to (%d,%d) for MPI_COPY when j=%d\n", myrow, mycol, myrow, vic_id%npcol, j);
					MPI_Send(A+(lcindx-1)*np_A, np_A*nb, MPI_DOUBLE, vic_id%npcol, 1, cs->rowcomm);
				}
			}
		}
	}

	// copy the correct checksum back //
	Cpdgemr2d (MA, NAA-MA, Achk_old, 1, 1, descAchk, A, 1, MA+1, descA, grid->ictxt);

	if (npp*nqq>0)
	{
		free (Achk_new);
		free (Achk_old);
	}

#ifdef OLD_RECOVER_ZONE1
	// range of recover
	if (iamvic)
	{
		infog1l_ (&errstep, &nb, &npcol, &mycol, &izero, &lindx, &rocsrc);
		nchm[0] = npp_A;
		nchm[1] = lindx-1;
		if (mycol != rocsrc)
			nchm[1]--;
	}
	MPI_Bcast( nchm, 2, MPI_INT, vic_id, MPI_COMM_WORLD );

	for (j=0; j<=nchm[1]; j+=nb)	//col
	{
		if (iamvic)
		{
			int j_1 = j+1;
			gj = indxl2g_( &j_1, &nb, &mycol, &izero, &npcol );
		}
		MPI_Bcast( &gj, 1, MPI_INT, vic_id, MPI_COMM_WORLD );

		// offset A
		infog1l_ (&gj, &nb, &npcol, &mycol, &izero, &lcindx, &csrc); 
		if (mycol == csrc)
			Aoff = A+(lcindx-1)*np_A;

		for (i=0; i<nchm[0]; i+=nb)	//row
		{
			// determine who owns the targeting checksum
			gi = i+MA+1;

			infog2l_ (&gi, &gj, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
			int root = (rsrc == myrow && csrc == mycol)?(grid->myrank_mpi):(-1);

			MPI_Allreduce ( MPI_IN_PLACE, &root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			// determine who is in the computation
			int in =	((mycol==csrc) &&
					((i<npp_A)||(root==grid->myrank_mpi)))
				*(root+1);												// root+1 because root can be 0

			MPI_Comm_split(MPI_COMM_WORLD, in, myrow, &part_comm); 

			// translate root to part_comm
			int part_root;
			MPI_Comm_group(part_comm, &part_group);
			MPI_Group_translate_ranks (world_group, 1, &root, part_group, &part_root);

			// calculate the recovery 
			double *S=NULL, *R=NULL;
			if (in)
			{
				//printf ("(%d,%d) is in when i=%d, j=%d, root=%d\n", myrow, mycol, i, j, grid->myrank_mpi==root);

				// copy participants into seperate buffers
				S = (double*)malloc(nb*nb*sizeof(double));
				if (iamvic || i>=npp_A)
					memset (S, 0, nb*nb*sizeof(double));
				else
					for (k=0; k<nb; k++)
						memcpy (S+k*nb, Aoff+i+k*np_A, nb*sizeof(double));

				if (grid->myrank_mpi==root)
				{
					R = (double*)malloc(nb*nb*sizeof(double));
					memset (R, 0, nb*nb*sizeof(double));
				}

				MPI_Reduce (S, R, nb*nb, MPI_DOUBLE, MPI_SUM, part_root, part_comm);

				if (grid->myrank_mpi==root)
				{
					CC = A+(lcindx-1)*np_A+(lrindx-1);
					int ii,jj;
					for (ii=0; ii<nb; ii++)
						for (jj=0; jj<nb; jj++)
							R[jj*nb+ii] = *(CC+jj*np_A+ii) - R[jj*nb+ii]; 
				}

				free (S);
			}

			//MPI_Barrier (MPI_COMM_WORLD);

			// fix the holes 
			MPI_Request req1,req2;
			MPI_Status status;
			if (grid->myrank_mpi == root)	// root sends the patch
				MPI_Isend(R, nb*nb, MPI_DOUBLE, vic_id, TAG, MPI_COMM_WORLD, &req2);

			if (iamvic)	//victim receives the patch
			{
				S = (double*)malloc(nb*nb*sizeof(double));
				MPI_Irecv(S, nb*nb, MPI_DOUBLE, root, TAG, MPI_COMM_WORLD, &req1);
			}

			if (root==grid->myrank_mpi)
			{
				MPI_Wait(&req2, &status);
				free (R);
			}

			if (iamvic)
			{
				MPI_Wait(&req1, &status);

				// victim puts patch in place
				CC = A+j*np_A+i;
				for (k=0; k<nb; k++)
					memcpy (CC+k*np_A, S+k*nb, nb*sizeof(double));

				free (S);
			}

			MPI_Comm_free (&part_comm);
			R=S=NULL;
		}
	}
#endif
}

void recover_local (char side, double *A, int MA, int NA, int *descA, int errx, int erry, int errstep, int vic_id, t_Grid *grid, t_checksuite *cs)
{
	int izero = 0, i, j, k;
	int tt[2], g[2];
	int lindx, rocsrc, nchm[2];

	// grid parameters
	int	nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	int iamvic = (grid->myrank_mpi==vic_id);

	// get nchkc
	int np_A = cs->np_L; 

	NA = npcol*nb;
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );

	int nqq_A = cs->nq_L; 

	// horizontal, and only horizontal
	if (side == 'b' || side == 'r')
	{
		if (nqq_A<nq_A)
		{
			printf ("(%d,%d): quick return from recovering checksum\n",myrow, mycol);
			return;		// quick return is victim doesn't carry checksum
		}

		int i_c, j_c;
		if (iamvic)
		{
			i_c = np_A;
			j_c = (nqq_A-nq_A);
		}
		MPI_Bcast( &i_c, 1, MPI_INT, vic_id, MPI_COMM_WORLD );
		MPI_Bcast( &j_c, 1, MPI_INT, vic_id, MPI_COMM_WORLD );

		for (i=0; i<i_c; i+=nb)
			for (j=0, k=0; j<j_c; j+=nb, k++)
			{
				if (iamvic)
				{
					int i_1 = i+1;
					g[0] = indxl2g_( &i_1, &nb, &myrow, &izero, &nprow );
					int j_1 = nq_A + j + 1;
					g[1] = indxl2g_( &j_1, &nb, &mycol, &izero, &npcol );

					tt[1]=(((g[1]-NA-1)/nb)%2==0)?(g[1]+nb):(g[1]-nb);
					tt[0]=g[0];
					//printf ("(%d,%d): block (%d,%d) <------ (%d,%d)\n", myrow, mycol, g[0], g[1], tt[0],tt[1]);
				}
				MPI_Bcast( tt, 2, MPI_INT, vic_id, MPI_COMM_WORLD );
				MPI_Bcast( g, 2, MPI_INT, vic_id, MPI_COMM_WORLD );

				Cpdgemr2d (nb, nb, A, tt[0], tt[1], descA, A, g[0], g[1], descA, grid->ictxt);
			}

		recover_local_fillin (A, MA, NA, descA, errx, erry, errstep, vic_id, grid, cs);
	}
}

void recover_checksum (char side, double *A, int MA, int NA, int *descA, int errx, int erry, int errstep, int vic_id, t_Grid *grid, t_checksuite *cs)
{
	int izero = 0, i, j, k;
	int tt[2], g[2];
	int lindx, rocsrc, nchm[2];

	// grid parameters
	int	nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	int iamvic = (grid->myrank_mpi==vic_id);

	// get nchkc
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );

	int NAA = NA+cs->Nc;
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
	int nqq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );

	// horizontal, and only horizontal
	if (side == 'b' || side == 'r')
	{
		/*
		if (iamvic)
			printf ("(%d,%d): recovering checksum\n",myrow, mycol);
			*/

		if (nqq_A<nq_A)
		{
			printf ("(%d,%d): quick return from recovering checksum\n",myrow, mycol);
			return;		// quick return is victim doesn't carry checksum
		}

		int i_c, j_c;
		if (iamvic)
		{
			i_c = np_A;
			j_c = (nqq_A-nq_A);
		}
		MPI_Bcast( &i_c, 1, MPI_INT, vic_id, MPI_COMM_WORLD );
		MPI_Bcast( &j_c, 1, MPI_INT, vic_id, MPI_COMM_WORLD );

		for (i=0; i<i_c; i+=nb)
			for (j=0, k=0; j<j_c; j+=nb, k++)
			{
				if (iamvic)
				{
					int li = i+1;
					g[0] = indxl2g_( &li, &nb, &myrow, &izero, &nprow );
					int j_1 = nq_A + j + 1;
					g[1] = indxl2g_( &j_1, &nb, &mycol, &izero, &npcol );

					tt[1]=(((g[1]-NA-1)/nb)%2==0)?(g[1]+nb):(g[1]-nb);
					tt[0]=g[0];
					//printf ("(%d,%d): block (%d,%d) <------ (%d,%d)\n", myrow, mycol, g[0], g[1], tt[0],tt[1]);
				}
				MPI_Bcast( tt, 2, MPI_INT, vic_id, MPI_COMM_WORLD );
				MPI_Bcast( g, 2, MPI_INT, vic_id, MPI_COMM_WORLD );

				Cpdgemr2d (nb, nb, A, tt[0], tt[1], descA, A, g[0], g[1], descA, grid->ictxt);
			}
	}
}


void recover (int CHK, int Z1, int Z2, int Z3, double *A, int m, int n, int *descA, int errx, int erry, int errstep, t_Grid *grid, t_checksuite *cs)
{
	// broadcast victim's id
	int vic_id = (grid->myrow==errx && grid->mycol==erry)?(grid->myrank_mpi):-1;
	MPI_Allreduce ( MPI_IN_PLACE, &vic_id, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	// recover checksum
	if (CHK)
	{
		recover_checksum ('b', A, m, n, descA, errx, erry, errstep, vic_id, grid, cs);
		//ROOTSAY ("checksum recovered");
	}

	MPI_Barrier (MPI_COMM_WORLD);

	// recover non-Q panels 
	if (Z1)
	{
		recover_zone1 (A, m, n, descA, errx, erry, errstep, vic_id, grid, cs);
		//ROOTSAY ("ABFT region recovered");
	}

	MPI_Barrier (MPI_COMM_WORLD);
	
	if (Z3)
	{
		recover_local ('b', cs->localcopy, m, n, cs->descL, errx, erry, errstep, vic_id, grid, cs);
		//ROOTSAY ("localcopy recovered");
	}

	if (Z2)
	{
		recover_zone2 (A, m, n, descA, errx, erry, errstep, vic_id, grid, cs);
		MPI_Barrier (MPI_COMM_WORLD);
		//ROOTSAY ("checkpointing region recovered");
	}

	MPI_Barrier (MPI_COMM_WORLD);

}


void kill_proc(double *A, int *descA, t_Grid *grid)
{
	int izero=0;
	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	int np_A = numroc_( &descA[M_], &nb, &myrow, &izero, &nprow );
	int nq_A = numroc_( &descA[N_], &nb, &mycol, &izero, &npcol );

	//printf ("\nzeroing out proc (%d,%d) of size %d x %d at step %d\n", myrow, mycol, np_A, nq_A, real_err_step);

	memset (A, 0, np_A*nq_A*sizeof(double)); 
}

int gLDA;
void SumBlock( double *in, double *inout, int *len, MPI_Datatype *dptr ) 
{ 
	int i, j; 

	// sum
	for (i=0; i<NB; i++)
		for (j=0; j<NB; j++)
			inout[i*gLDA+j] += in[i*gLDA+j];
} 

void panel_reset(double *A, int *descA, int IA, int JA, t_checksuite *cs, t_Grid *grid)
{
	int i, j, k;
	int izero = 0;
	double *Aoff;
	int MaxLocalRow = cs->MaxLocalRow;
	int MaxLocalCol = cs->MaxLocalCol;

	int nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	
	int MA = descA[2]; 
	int NA = descA[3]-cs->Nc; 
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
	
	// determine number of rows to checkpoint locally 
	int lrindx, lcindx, rsrc, csrc;
	infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);

	int len = (np_A-lrindx+1)*sizeof(double);
	if (len>0 && lcindx<nq_A && lrindx<np_A)
	{
		// make a local_copy if demanded
		//printf ("(%d,%d) copying %d row %d col at (%d,%d)\n", myrow, mycol, (np_A-lrindx+1), nb, IA, JA);
		Aoff = A+(lcindx-1)*np_A+(lrindx-1);
		for (i=0; i<nb; i++)
			memcpy (Aoff+i*np_A, cs->localcopy+i*np_A, len);
	}
}
void local_copy(double *A, int *descA, int IA, int JA, t_checksuite *cs, t_Grid *grid)
{
	int i, j, k;
	int izero = 0;
	double *Aoff;
	int MaxLocalRow = cs->MaxLocalRow;
	int MaxLocalCol = cs->MaxLocalCol;

	int nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	
	int MA = descA[2]; 
	int NA = descA[3]-cs->Nc; 
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
						
	// determine number of rows to checkpoint locally 
	int lrindx, lcindx, rsrc, csrc;
	infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);

	int len = (np_A-lrindx+1)*sizeof(double);
	if (len>0 && lcindx<nq_A && lrindx<np_A)
	{
		// make a local_copy if demanded
//		printf ("(%d,%d) copying %d row %d col at (%d,%d)\n", myrow, mycol, (np_A-lrindx+1), nb, IA, JA);
		Aoff = A+(lcindx-1)*np_A+(lrindx-1);
		for (i=0; i<nb; i++)
			memcpy (cs->localcopy+i*np_A, Aoff+i*np_A, len);
	}

	// checkpoint local copy
	hor_checkpointing_local (cs->localcopy, cs->descL, 1, 1, 1, cs, grid);
}

/*
 * perform checkpointing for the local copy
 * assuming JA is on process column 0
 */
void hor_checkpointing_local (double *A, int *descA, int IA, int JA, int iscopy, t_checksuite *cs, t_Grid *grid)
{
	int i, j, k;
	int izero = 0;
	double *Aoff;

	int nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	
	int MA = descA[2]; 
	int np_A = cs->np_L; 
	
	int NA = descA[3]-2*nb; 
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
	
	// temp buffers
	double *S = cs->S;
	double *R = cs->R;

	// determine number of rows to checkpoint locally 
	int lrindx, lcindx, rsrc, csrc;
	infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);
	int rows = np_A-lrindx+1;

	// determine target locatoin
	int tc = nb*npcol+1;
	//printf ("(%d,%d) JA=%d, tc=%d\n", myrow, mycol, JA, tc);	

	// find root
	int root, rootcr, rootcc;
	infog2l_ (&IA, &tc, descA, &nprow, &npcol, &myrow, &mycol, &rootcr, &rootcc, &rsrc, &csrc);
	root = (mycol == csrc);
//	if (root)
//		printf ("(%d,%d) is root (csrc=%d) when JA=%d\n", myrow, mycol, csrc, JA);

	// checkpointing
	for (i=0; i<rows; i+=nb)
	{
		// determine participants 
		int in = ((lrindx+i)<np_A);
		if (in)
		{
			Aoff = A+(lcindx-1)*np_A+(lrindx-1)+i;
			for (k=0; k<nb; k++)
				memcpy (S+k*nb, Aoff+k*np_A, nb*sizeof(double));
		}
		else
			memset (S, 0, nb*nb*sizeof(double));

		if ((lrindx+i)<np_A)
		{
			MPI_Reduce (S, R, nb*nb, MPI_DOUBLE, MPI_SUM, csrc, cs->rowcomm);

			if (root)
			{
				Aoff = A+(rootcc-1)*np_A+(rootcr-1)+i;  
				for (k=0; k<nb; k++)
					memcpy (Aoff+k*np_A, R+k*nb, nb*sizeof(double));
			}
		}
	}
		
	if (iscopy)
	{
		//make a copy
		Cpdgemr2d (MA-IA+1, nb, A, IA, tc, descA, A, IA, tc+nb, descA, grid->ictxt);
	}
}

/*
 * perform one step of checkpointing
 * assuming JA is on process column 0
 */
void hor_checkpointing (double *A, int *descA, int IA, int JA, int iscopy, t_checksuite *cs, t_Grid *grid)
{
	int i, j, k;
	int izero = 0;
	double *Aoff;
	int MaxLocalRow = cs->MaxLocalRow;
	int MaxLocalCol = cs->MaxLocalCol;

	int nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;
	
	int MA = descA[2]; 
	int NA = descA[3]-cs->Nc; 
	int np_A = numroc_( &MA, &nb, &myrow, &izero, &nprow );
	int nq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );
						
	
	// temp buffers
	//double *S = (double*)malloc(nb*nb*sizeof(double));
	//double *R = (double*)malloc(nb*nb*sizeof(double));
	double *S = cs->S;
	double *R = cs->R;

	// determine number of rows to checkpoint locally 
	int lrindx, lcindx, rsrc, csrc;
	infog2l_ (&IA, &JA, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc);
	int rows = np_A-lrindx+1;

	/*
	// make a local_copy if demanded
	if (iscopy && (lcindx<nq_A) && (lrindx+i)<np_A)
	{
		printf ("(%d,%d) copying %d row %d col at (%d,%d)\n", myrow, mycol, (np_A-lrindx), nb, IA, JA);
		int len = (np_A-lrindx)*sizeof(double);
		Aoff = A+(lcindx-1)*np_A+(lrindx-1);
		for (i=0; i<nb; i++)
			memcpy (cs->localcopy+i*np_A, Aoff+i*np_A, len);
	}
	*/

	// determine target locatoin
	int tc = ((NA-JA+1)/(nb*npcol)+((NA-JA+1)%(nb*npcol)!=0)-1)*2*nb+NA+1;
	//printf ("(%d,%d) JA=%d, tc=%d\n", myrow, mycol, JA, tc);	

	// find root
	int root, rootcr, rootcc;
	infog2l_ (&IA, &tc, descA, &nprow, &npcol, &myrow, &mycol, &rootcr, &rootcc, &rsrc, &csrc);
	root = (mycol == csrc);
//	if (root)
//		printf ("(%d,%d) is root (csrc=%d) when JA=%d\n", myrow, mycol, csrc, JA);

	// checkpointing
	for (i=0; i<rows; i+=nb)
	{
		// determine participants 
		int in = ((lcindx<nq_A) && (lrindx+i)<np_A);
		if (in)
		{
			//printf ("(%d,%d) IA=%d, JA=%d, lrindx=%d, lcindx=%d, rows=%d, i=%d\n", myrow, mycol, IA, JA, lrindx, lcindx, rows, i);	
			Aoff = A+(lcindx-1)*np_A+(lrindx-1)+i;
			for (k=0; k<nb; k++)
				memcpy (S+k*nb, Aoff+k*np_A, nb*sizeof(double));
		}
		else
			memset (S, 0, nb*nb*sizeof(double));

		if ((lrindx+i)<np_A)
		{
			MPI_Reduce (S, R, nb*nb, MPI_DOUBLE, MPI_SUM, csrc, cs->rowcomm);

			if (root)
			{
				Aoff = A+(rootcc-1)*np_A+(rootcr-1)+i;  
				for (k=0; k<nb; k++)
					memcpy (Aoff+k*np_A, R+k*nb, nb*sizeof(double));
			}
		}
	}
		
	//free (S);
	//free (R);

	if (iscopy)
	{
		//make a copy
		Cpdgemr2d (MA-IA+1, nb, A, IA, tc, descA, A, IA, tc+nb, descA, grid->ictxt);
	}
}

/*
 * checkpoint A using G into Ac on left or bottom or both
 */
void checkpointing_short (char side,
		double *A, int MA, int NA, int IA, int JA, int *descA, 
		double *Acd, int IAcd, int JAcd, int *descAcd,
		t_checksuite *cs, t_Grid *grid, int errstep)
{
	// quick return if possible
	if (side != 'r' && side != 'b' && side != 'd')	ABORT;

	static int clear=1; 
	int izero=0;
	int i, j, k, r;
	int lrindx, lcindx, rsrc, csrc;
	double *R=NULL, *S=NULL;

	int MaxLocalRow = cs->MaxLocalRow;
	int MaxLocalCol = cs->MaxLocalCol;

	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	MPI_Group world_group, part_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	if (side=='d' || side=='b')
	{
		hor_checkpointing (A, descA, IA, JA, COPY, cs, grid); 
	}
	else if (side=='r' || side=='b')
	{
		// determine local data size (non-checksum)
		int MAA = MA+cs->Mc;
		int np_A = numroc_( &MAA, &nb, &myrow, &izero, &nprow );
		int npp_A = numroc_( &MA, &nb, &myrow, &izero, &nprow ); // no checksum np_A 

		int NAA = NA+cs->Nc;
		int nq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );
		int nqq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );

		gLDA = np_A;

		// zero out the initial checksum
		if (nq_A > nqq_A)
			memset (A+np_A*nqq_A, 0, (nq_A-nqq_A)*np_A*sizeof(double));

		// checkpointing
		/*
		   printf ("%d: MaxLocalRow=%d, MaxLocalCol=%d\n", grid->myrank_mpi, MaxLocalRow, MaxLocalCol);
		   printf ("%d: npp_A=%d, nqq_A=%d\n", grid->myrank_mpi, npp_A, nqq_A);
		 */

		int Qnb=npcol*nb;
		for (i=0; i<NA; i+=Qnb)
		{
			if ((i+1)<errstep && (i+1+Qnb)>errstep)
			{
				if (grid->myrank_mpi==0)
					printf ("skipping step %d in re-checkpointing\n", i);
			}
			else
			{
				// creating checksum
				hor_checkpointing (A, descA, IA, JA+i, COPY, cs, grid); 
			}
		}
	}
}

/*
 * checkpoint A using G into Ac on left or bottom or both
 */
void checkpointing (char side,
		double *A, int MA, int NA, int IA, int JA, int *descA, 
		double *Acd, int IAcd, int JAcd, int *descAcd,
		t_checksuite *cs, t_Grid *grid)
{
	// quick return if possible
	if (side != 'r' && side != 'b' && side != 'd')	ABORT;

	static int clear=1; 
	int izero=0;
	int i, j, k, r;
	int lrindx, lcindx, rsrc, csrc;
	double *R=NULL, *S=NULL;

	int MaxLocalRow = cs->MaxLocalRow;
	int MaxLocalCol = cs->MaxLocalCol;

	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	MPI_Group world_group, part_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	if (side=='d' || side=='b')
	{
		hor_checkpointing (A, descA, IA, JA, COPY, cs, grid); 
	}
	else if (side=='r' || side=='b')
	{
		// determine local data size (non-checksum)
		int MAA = MA+cs->Mc;
		int np_A = numroc_( &MAA, &nb, &myrow, &izero, &nprow );
		int npp_A = numroc_( &MA, &nb, &myrow, &izero, &nprow ); // no checksum np_A 

		int NAA = NA+cs->Nc;
		int nq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );
		int nqq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );

		gLDA = np_A;

		// zero out the initial checksum
		if (nq_A > nqq_A)
			memset (A+np_A*nqq_A, 0, (nq_A-nqq_A)*np_A*sizeof(double));

		// checkpointing
		/*
		   printf ("%d: MaxLocalRow=%d, MaxLocalCol=%d\n", grid->myrank_mpi, MaxLocalRow, MaxLocalCol);
		   printf ("%d: npp_A=%d, nqq_A=%d\n", grid->myrank_mpi, npp_A, nqq_A);
		 */

		int Qnb=npcol*nb;
		for (i=0; i<NA; i+=Qnb)
		{
			// creating checksum
			hor_checkpointing (A, descA, IA, JA+i, COPY, cs, grid); 
		}
	}
}

/*
 * verify checkpointing after computing
 */
void verify_checkpointing (int M, int N, double *A, int *descA, t_Grid *grid, t_checksuite *cs)
{

	double *Achk_old, *Achk_new;
	int descAchk[9], npp, nqq;

	distr_matrix (false, &Achk_old, descAchk, M, N-M, grid, &npp, &nqq);
	distr_matrix (false, &Achk_new, descAchk, M, N-M, grid, &npp, &nqq);

	// copy the old checksum into buffer
	Cpdgemr2d (M, N-M, A, 1, M+1, descA, Achk_old, 1, 1, descAchk, grid->ictxt);

	// re-calculate the checksum
	checkpointing ('r', A, M, M, 1, 1, descA, A, 1, M, descA, cs, grid);

	// copy the new checksum into buffer
	Cpdgemr2d (M, N-M, A, 1, M+1, descA, Achk_new, 1, 1, descAchk, grid->ictxt);

	// cross check
	int N_M = N-M, ione=1;
	double done = 1, mone = -1;
	pdmatadd_ ( &M, &N_M, &done, Achk_new, &ione, &ione, descAchk, &mone, Achk_old, &ione, &ione, descAchk);
	double resid = pdlange_("F", &M, &N_M, Achk_old, &ione, &ione, descAchk, NULL)/pdlange_("F", &M, &N_M, Achk_new, &ione, &ione, descAchk, NULL)/M;
	if (grid->myrank_mpi==0)
		printf ("resid=%e\n", resid);

	if (npp*nqq>0)
	{
		free (Achk_new);
		free (Achk_old);
	}
}

#ifdef old
/*
 * checkpoint A using G into Ac on left or bottom or both
 */
void checkpointing (char side,
		double *A, int MA, int NA, int IA, int JA, int *descA, 
		double *Acd, int IAcd, int JAcd, int *descAcd,
		t_checksuite *cs, t_Grid *grid)
{
	// quick return if possible
	if (side != 'r' && side != 'b' && side != 'd')	ABORT;

	static int clear=1; 
	int izero=0;
	int i, j, k, r;
	int lrindx, lcindx, rsrc, csrc;
	double *R=NULL, *S=NULL;

	int MaxLocalRow = cs->MaxLocalRow;
	int MaxLocalCol = cs->MaxLocalCol;

	int nb=grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	MPI_Group world_group, part_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	if (side=='r' || side=='b')
	{
		// determine local data size (non-checksum)
		int MAA = MA+cs->Mc;
		int np_A = numroc_( &MAA, &nb, &myrow, &izero, &nprow );
		int npp_A = numroc_( &MA, &nb, &myrow, &izero, &nprow ); // no checksum np_A 

		int NAA = NA+cs->Nc;
		int nq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );
		int nqq_A = numroc_( &NA, &nb, &mycol, &izero, &npcol );

		gLDA = np_A;

		// zero out init checksum
		if (nq_A > nqq_A)
			memset (A+np_A*nqq_A, 0, (nq_A-nqq_A)*np_A*sizeof(double));

		// checkpointing
		/*
		   printf ("%d: MaxLocalRow=%d, MaxLocalCol=%d\n", grid->myrank_mpi, MaxLocalRow, MaxLocalCol);
		   printf ("%d: npp_A=%d, nqq_A=%d\n", grid->myrank_mpi, npp_A, nqq_A);
		 */

		int offset=npcol*nb;
		for (i=0; i<MaxLocalRow; i+=nb)
		{
			for (r=0; r<offset; r+=nb)
			{
				for (j=r; j<MaxLocalCol; j+=offset)
				{
					// determine who owns the targeting checksum
					int i_1 = i+1;
					int gi = indxl2g_( &i_1, &nb, &myrow, &izero, &nprow );
					int gj = JAcd + j;
					infog2l_ (&gi, &gj, descAcd, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );

					int root = rsrc*npcol+csrc;

					// determine who is in the computation
					int in = (((i<npp_A) && (j<nqq_A))||(root==grid->myrank_mpi && i<npp_A))*(root+1);	// root+1 because root can be 0
					// & i<npp_A because root might not be within data range
					if (myrow==rsrc && cs->hcomm[j/nb].created == 0)
					{
						MPI_Comm_split(cs->rowcomm, in, mycol, &(cs->hcomm[j/nb].part_comm)); 

						// translate root to part_comm
						MPI_Comm_group(cs->hcomm[j/nb].part_comm, &part_group);
						MPI_Group_translate_ranks (cs->rowgroup, 1, &csrc, part_group, &(cs->hcomm[j/nb].part_root));
						cs->hcomm[j/nb].created = 1;
					}

					// do the sum
					if (in)
					{

						//printf ("%d(%d,%d): (%d,%d) -> %d is root when i=%d, j=%d, lrindx=%d, lcindx=%d, in=%d, part_root=%d\n", 
						//			grid->myrank_mpi, myrow, mycol, gi, gj, root, i, j, lrindx, lcindx, in, part_root);

						S = (double*)malloc(nb*nb*sizeof(double));
						memset (S, 0, nb*nb*sizeof(double));

						if (grid->myrank_mpi==root)
						{
							R = (double*)malloc(nb*nb*sizeof(double));
							memset (R, 0, nb*nb*sizeof(double));
						}

						if (j<nqq_A && i<npp_A)
							for (k=0; k<nb; k++)
								memcpy (S+k*nb, A+j*np_A+i+k*np_A, nb*sizeof(double));

						MPI_Reduce (S, R, nb*nb, MPI_DOUBLE, MPI_SUM, cs->hcomm[j/nb].part_root, cs->hcomm[j/nb].part_comm);

						if (grid->myrank_mpi==root)
						{
							for (k=0; k<nb; k++)
								memcpy (Acd+(lcindx-1)*np_A+(lrindx-1)+k*np_A, R+k*nb, nb*sizeof(double));
							free (R);
						}

						free(S); 
					}

					R=S=NULL;
				}
			}
		}

		// make a copy
		int nchkc = cs->MaxLocalCol;
		Cpdgemr2d (MA, nchkc, Acd, IAcd, JAcd, descA, Acd, IAcd, JAcd+nchkc+(cs->IsNskew)*nb, descA, grid->ictxt);
	}

	if (side=='d' || side=='b')
	{
		// determine local data size (non-checksum)
		int MAA = MA+cs->Mc;
		int np_A = numroc_( &MAA, &nb, &myrow, &izero, &nprow );
		int npp_A = numroc_( &MA, &nb, &myrow, &izero, &nprow ); // no checksum np_A 

		int NAA = NA+cs->Nc;
		int nq_A = numroc_( &NAA, &nb, &mycol, &izero, &npcol );

		gLDA = np_A;

		// zero out init checksum (only do this during the 1st call)
		if (clear)
		{
			if (np_A > npp_A)
				for (j=0; j<nq_A; j++)
					memset (A+npp_A+np_A*j, 0, (np_A-npp_A)*sizeof(double)); 
			clear = 0;
		}

		// offset A to (IA, JA)
		infog1l_ (&JA, &nb, &npcol, &mycol, &izero, &lcindx, &csrc); 
		if (mycol == csrc)
			A += (lcindx-1)*np_A;

		// checkpointing
		if (mycol==csrc)
		{
			/*
			int cut;
			infog1l_ (&IA, &nb, &nprow, &myrow, &izero, &cut, &rsrc); 
			cut--;
			if (myrow<rsrc)
				cut-=nb;
				*/

			int offset=nprow*nb;
			//printf ("(%d,%d) JA=%d, cut=%d, offset=%d\n", myrow, mycol, JA, cut, offset);
			for (r=0; r<offset; r+=nb)
			{
				for (i=r; i<MaxLocalRow; i+=offset)
				{
					// determine who owns the targeting checksum
					int gi = i+IAcd;
					int gj = JAcd;
					infog2l_ (&gi, &gj, descA, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );

					int root = rsrc*npcol+csrc;

					// determine who is in the computation
					int i_1 = i+1;
					int myg = indxl2g_(&i_1, &nb, &myrow, &izero, &nprow);

					// do the sum
					if (i<npp_A && myg>=IA)
					{
						for (k=0; k<nb; k++)
							memcpy (cs->S+k*nb, A+i+k*np_A, nb*sizeof(double));
			//			printf ("(%d,%d) offering (%d,%d) --> (%d,%d)\n", myrow, mycol, myg, gj, gi, gj);
					}
					else
						memset (cs->S, 0, nb*nb*sizeof(double));

					MPI_Reduce (cs->S, cs->R, nb*nb, MPI_DOUBLE, MPI_SUM, rsrc, cs->colcomm);
					cc++;

					if (grid->myrank_mpi==root)
					{
						for (k=0; k<nb; k++)
							memcpy (Acd+(lcindx-1)*np_A+(lrindx-1)+k*np_A, cs->R+k*nb, nb*sizeof(double));
					}
				}
			}
		}

		// make a copy
		int nchkr = cs->MaxLocalRow;
		Cpdgemr2d (nchkr, nb, Acd, IAcd, JAcd, descA, Acd, IAcd+nchkr+(cs->IsMskew)*nb, JAcd, descA, grid->ictxt);
	}
}
#endif

#ifdef JUNK
/*
 * checkpoint A using G into Ac on left or bottom or both
 */
void checkpointing (char side,
		double *A, int MA, int NA, int IA, int JA, int *descA, 
		double *Acd, int IAcd, int JAcd, int *descAcd,
		t_checksuite *cs, t_Grid *grid)
{

	// quick return if possible
	if (side != 'r' && side != 'b' && side != 'd')	ABORT;

	int M, N, K;
	int ione = 1;
	double pone = 1.0;
	double dzero = 0.0;
	double *T;
	int descT[9];
	int np_T, nq_T;

	if (side=='r' || side=='b')
	{
		M = MA;
		K = NA;
		N = NGr;

		// checkpointing on the right
		pdgemm_ ("N", "N", &M, &N, &K, &pone,	A, &IA, &JA, descA,
				Gr, &IGr, &JGr, descGr,
				&dzero,	Acr, &IAcr, &JAcr, descAcr); 

		side = (side=='r')?'x':side;
	}


	if (side=='d' || side=='b')
	{
		distr_matrix (false, &T, descT, MGd, NA, grid, &np_T, &nq_T);
		pdlacpy_ ("A", &MGd, &NA, Gd, &ione, &ione, descGd, T, &ione, &ione, descT);
		pdtrmm_ ("R", "L", "N", "U", &MGd, &NA, &pone, A, &IA, &JA, descA, 
				T, &ione, &ione, descT);

		int INA = IA+NA;
		int MNA = NGd-MA;
		int t_1 = JGd+MA;
		pdgemm_ ("N", "N", &MGd, &NA, &MNA, &pone,	Gd, &IGd, &t_1, descGd,
				A, &INA, &JA, descA,
				&pone,	T, &ione, &ione, descT); 

		//pdlacpy_ ("A", &MGd, &NA, T, &ione, &ione, descT, Acd, &IAcd, &JAcd, descAcd);	//	<---- pdlacpy does NOT work between matrices with 
		//		  different distribution
		Cpdgemr2d (MGd, NA, T, ione, ione, descT, Acd, IAcd, JAcd, descAcd, grid->ictxt);

		if (np_T*nq_T>0)
			free (T);
	}
}
#endif

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

/*
 * produce generator matrix, and then check matrix from the generator matrix
 * M, N are the sizes of the generator matrix
 */
void generator_check (char side, char type, double **G, int *descG, 
		double **H, int *descH, 
		int M, int N, t_Grid *grid,
		int *np_G, int *nq_G,
		int *np_H, int *nq_H) 
{
	int info;
	int i, j;
	int izero = 0;
	int ione = 1;

	// grid parameters
	int	nb = grid->nb;
	int	ictxt = grid->ictxt;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;
	int	nprow = grid->nprow;
	int	npcol = grid->npcol;

	// allocate the generator matrix and check matrix
	int MM;
	int NN;
	if (side=='r')
	{	MM = M+1; NN = N;	}
	else
	{	MM = M; NN = N+1;	}

	int np_iG = numroc_( &M, &nb, &myrow, &izero, &nprow );
	int nq_iG = numroc_( &N, &nb, &mycol, &izero, &npcol );

	int np_iH = numroc_( &MM, &nb, &myrow, &izero, &nprow );
	int nq_iH = numroc_( &NN, &nb, &mycol, &izero, &npcol );

	if (np_iH*nq_iH!=0)
	{
		*H = (double *)malloc(np_iH*nq_iH*sizeof(double)) ;
		if (*H == NULL) ABORT;
		memset (*H, 0, np_iH*nq_iH*sizeof(double));
	}
	if (np_iG*nq_iG!=0)
	{
		*G = (double *)malloc(np_iG*nq_iG*sizeof(double)) ;
		if (*G == NULL) ABORT;
		memset (*G, 0, np_iG*nq_iG*sizeof(double));
	}

	int itemp = MAX( 1, np_iH );
	descinit_( descH, &MM, &NN, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
	if (info != 0) ABORT;

	itemp = MAX( 1, np_iG );
	descinit_( descG, &M, &N, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
	if (info != 0) ABORT;

	if (type == 'w')
	{
		// fill in random numbers
		int seed = (type=='r')?900:500;
		pdmatgen_ (&ictxt, "N", "N", &M, &N, &nb, &nb, *G, 
				descG+8, descG+6, descG+7, 
				&seed, &izero, &np_iG, &izero, &nq_iG, 
				&myrow, &mycol, &nprow, &npcol);
	}
	else
	{
		for (i=0; i<np_iG; i++)
			for (j=0; j<nq_iG; j++)
				(*G)[j*np_iG+i] = 1;
	}

	pdlacpy_ ("A", &M, &N, *G, &ione, &ione, descG, *H, &ione, &ione, descH);

	int lrindx, lcindx, rsrc, csrc; 
	if (side == 'd')
	{
		// put -I to the end of Hd
		for (i=1; i<=M; i++)
		{
			j = N + i;

			infog2l_ (&i, &j, descH, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
			if (myrow == rsrc && mycol == csrc)
				(*H)[(lcindx-1)*np_iH+lrindx-1] = -1;
		}
	}
	else
	{
		// put -I to the end of Hr
		for (i=1; i<=N; i++)
		{
			j = M + i;

			infog2l_ (&j, &i, descH, &nprow, &npcol, &myrow, &mycol, &lrindx, &lcindx, &rsrc, &csrc );
			if (myrow == rsrc && mycol == csrc)
				(*H)[(lcindx-1)*np_iH+lrindx-1] = -1;
		}
	}

	/* set np and nq */
	*np_H = np_iH;
	*nq_H = nq_iH;
	*np_G = np_iG;
	*nq_G = nq_iG;
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


/*
 * get a healthy node to do recovery on
 */
int get_recovery_node (int *err_cords, int noe)
{
	return 0;	//for now we'
}


void get_checksuite (int Nc, int Mc, int IsNskew, int IsMskew, int MaxLocalCol, int MaxLocalRow, t_checksuite *cs, t_Grid *grid)
{
	int nb = grid->nb;
	int	myrow = grid->myrow;
	int	mycol = grid->mycol;

	cs->Nc=Nc;
	cs->Mc=Mc;
	cs->IsNskew=IsNskew;
	cs->IsMskew=IsMskew;
	cs->MaxLocalCol=MaxLocalCol;
	cs->MaxLocalRow=MaxLocalRow;

	cs->R = (double*)calloc(nb*nb, sizeof(double));
	cs->S = (double*)calloc(nb*nb, sizeof(double));

	int i;
	/*
	// generate vcomm
	cs->vcomm = (t_Comm *)malloc(MaxLocalRow/nb*sizeof(t_Comm));
	if (cs->vcomm == NULL)
	{
	printf ("error allocating vertical communicator array\n");
	exit(1);
	}
	for (i=0; i<MaxLocalRow/nb; i++)
	{
	(cs->vcomm[i]).part_root=-1;
	(cs->vcomm[i]).created=0;
	}
	 */

	// generate hcomm
	cs->hcomm = (t_Comm *)malloc(MaxLocalCol/nb*sizeof(t_Comm));
	if (cs->hcomm == NULL)
	{
		printf ("error allocating vertical communicator array\n");
		exit(1);
	}
	for (i=0; i<MaxLocalCol/nb; i++)
	{
		(cs->hcomm[i]).part_root=-1;
		(cs->hcomm[i]).created=0;
	}

	MPI_Comm_split(MPI_COMM_WORLD, myrow, mycol, &(cs->rowcomm)); 
	MPI_Comm_split(MPI_COMM_WORLD, mycol, myrow, &(cs->colcomm)); 
	MPI_Comm_group(cs->colcomm, &(cs->colgroup));
	MPI_Comm_group(cs->rowcomm, &(cs->rowgroup));
}

void del_checksuite (t_checksuite *cs, int nb)
{
	cs->Nc=-1;
	cs->Mc=-1;
	cs->IsNskew=-1;
	cs->IsMskew=-1;

	int i;
	/*
	   for (i=0; i<cs->MaxLocalRow/nb; i++)
	   {
	   if (cs->vcomm[i].created == 1)
	   MPI_Comm_free (&(cs->vcomm[i].part_comm));
	   }
	   free (cs->vcomm);
	 */

	for (i=0; i<cs->MaxLocalCol/nb; i++)
	{
		if (cs->hcomm[i].created == 1)
			MPI_Comm_free (&(cs->hcomm[i].part_comm));
	}
	free (cs->hcomm);

	free (cs->S);
	free (cs->R);
}


