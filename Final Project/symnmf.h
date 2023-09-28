double** sym_compute(double **vectorX, int n, int d);
double** ddg_compute(double **A_similarity_matrix, int n);
double** norm_compute(double **A_similarity_matrix, int n, double **D_diagonal_degree_matrix);
double** compute_symnmf(double **H, double **W, int n, int k);
int* analysis_compute(double **H, int n, int k);