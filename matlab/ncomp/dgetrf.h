//extern "C" void dgetrf(
//	const int *m,			/* (input) */
//	const int *n,			/* (input) */
//	double *a,			/* a[n][lda] (input/output) */
//	const int *lda,			/* (input) */
//	int *ipiv,			/* ipiv[min(m,n)] (output) */
//	int *info			/* (output) */
//	);

extern "C" void dgetrf(
	const int &m,			// (input)
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	int *ipiv,			// ipiv[min(m,n)] (output)
	int &info			// (output)
	);

/*extern "C" void dgetrf(
	const int &m,			// (input)
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	int *ipiv,			// ipiv[min(m,n)] (output)
	int &info			// (output)
  );*/

