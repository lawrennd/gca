//extern "C" void dgetri(
//	const int *n,			/* (input) */
//	double *a,			/* a[n][lda] (input/output) */
//	const int *lda,			/* (input) */
//	const int *ipiv,		/* ipiv[n] (input) */
//	double *work,			/* work[lwork] (workspace/output) */
//	const int *lwork,		/* (input) */
//	int *info			/* (output) */
//	);
extern "C" void dgetri(
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	const int *ipiv,		// ipiv[n] (input)
	double *work,			// work[lwork] (workspace/output)
	const int &lwork,		// (input)
	int &info			// (output)
	);
