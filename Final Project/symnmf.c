#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"

enum GoalType
{
    SYM,
    DDG,
    NORM,
    INVALID
};

enum GoalType get_goal_type(const char *goal)
{
    if (strcmp(goal, "sym") == 0)
    {
        return SYM;
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        return DDG;
    }
    else if (strcmp(goal, "norm") == 0)
    {
        return NORM;
    }
    else
    {
        return INVALID;
    }
}

double **build_mat(int rows, int cols)
{
    double **m = (double **)calloc(rows ,sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        m[i] = (double *)calloc(cols , sizeof(double));
        for (int j = 0; j < cols; j++)
        {
            m[i][j] = 0.0;
        }
    }
    return m;
}

double **matrixMultiply(double **A, double **B, int rowsA, int colsA, int colsB)
{
    double **result = (double **)calloc(rowsA , sizeof(double *));
    for (int i = 0; i < rowsA; i++)
    {
        result[i] = (double *)calloc(colsB , sizeof(double));
        for (int j = 0; j < colsB; j++)
        {
            result[i][j] = 0.0;
            for (int k = 0; k < colsA; k++)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

double **transposeMatrix(double **A, int rows, int cols)
{
    double **transposed = (double **)calloc(cols , sizeof(double *));
    for (int i = 0; i < cols; i++)
    {
        transposed[i] = (double *)calloc(rows , sizeof(double));
        for (int j = 0; j < rows; j++)
        {
            transposed[i][j] = A[j][i];
        }
    }
    return transposed;
}

double **update_H(double **H, double **W, int n, int k)
{
    double b = 0.5;

    double **Hnext = build_mat(n, k);
    double **HT = build_mat(n, k);
    double **Mone = build_mat(n, k);
    double **mechne = build_mat(n, k);
    double **temp_mechne = build_mat(n, n);


    HT = transposeMatrix(H, n, k);
    Mone = matrixMultiply(W, H, n, n, k);
    temp_mechne = matrixMultiply(H, HT, n, k, n);
    mechne = matrixMultiply(temp_mechne, H, n, n, k);
    


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < k; j++)
        {
            Hnext[i][j] = H[i][j] * (1 - b + b * (Mone[i][j] / mechne[i][j]));
        }
    }
    for (int i = 0; i < n; i++)
    {
        free(Mone[i]);
        free(mechne[i]);
        free(temp_mechne[i]);
    }

    for (int i = 0; i < k; i++)
    {
        free(HT[i]);
    }
    
    free(HT);
    free(Mone);
    free(mechne);
    free(temp_mechne);

    return Hnext;
}

double squaredFrobeniusNorm(double **A, double **B, int rows, int cols)
{
    double sum = 0.0;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            double diff = A[i][j] - B[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}




double **compute_symnmf(double **H, double **W, int n, int k)
{

    double eps = 0.0001;
    double **NewH = build_mat(n, k);
    int round = 0;
    while (round < 300)
    {
        round = round + 1;
        NewH = update_H(H, W, n, k);
        double diff = squaredFrobeniusNorm(H, NewH, n, k);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                H[i][j] = NewH[i][j];
            }
        }
        if (diff < eps)
        {
            break;
        }

    }

    for (int i = 0; i < n; i++)
    {
        free(NewH[i]);
    }

    free(NewH);

    return H;
}


int* analysis_compute(double **H, int n, int k){
    int *res = malloc(n * sizeof(int));
    int i, j;

    for (i = 0; i < n; i++) {
        double max = H[i][0];  
        int maxindex = 0;

        for (j = 1; j < k; j++) {
            if (H[i][j] > max) {
                max = H[i][j];
                maxindex = j;
            }
        }

        res[i] = maxindex;
    }

    return res;
}


double **sym_compute(double **vectorX, int n, int d)
{

    double **A_similarity_matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        A_similarity_matrix[i] = (double *)malloc(n * sizeof(double));
    }


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                double diff_sum = 0.0;


                for (int k = 0; k < d; k++)
                {
                    double diff = vectorX[i][k] - vectorX[j][k];
                    diff_sum += diff * diff;
                }


                A_similarity_matrix[i][j] = exp(-diff_sum / 2.0);
            }
            else
            {
                A_similarity_matrix[i][j] = 0.0;
            }
        }
    }
    return A_similarity_matrix;
}

double **ddg_compute(double **A_similarity_matrix, int n)
{

    double **D_diagonal_Degree_Matrix =  (double **)calloc(n, sizeof(double));
    for (int i = 0; i < n; i++)
    {
        D_diagonal_Degree_Matrix[i] =  (double *)calloc(n, sizeof(double));
    }

    for (int i = 0; i < n; i++)
    {
        double degree = 0.0;


        for (int j = 0; j < n; j++)
        {
            degree += A_similarity_matrix[i][j];
        }


        D_diagonal_Degree_Matrix[i][i] = degree;
    }




    return D_diagonal_Degree_Matrix;
}

double **norm_compute(double **A_similarity_matrix, int n, double **D_diagonal_degree_matrix)
{

    double **W_graph_laplacian_matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        W_graph_laplacian_matrix[i] = (double *)malloc(n * sizeof(double));
    }


    double **D_sqrt_inverse = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        D_sqrt_inverse[i] = (double *)malloc(n * sizeof(double));
        double sqrt_degree = sqrt(D_diagonal_degree_matrix[i][i]);
        for (int j = 0; j < n; j++)
        {
            D_sqrt_inverse[i][j] = (i == j) ? (1.0 / sqrt_degree) : 0.0;
        }
    }


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            W_graph_laplacian_matrix[i][j] = D_sqrt_inverse[i][i] * A_similarity_matrix[i][j] * D_sqrt_inverse[j][j];
        }
    }


    for (int i = 0; i < n; i++)
    {
        free(D_sqrt_inverse[i]);
    }
    free(D_sqrt_inverse);

    return W_graph_laplacian_matrix;
}



void print_matrix(double **matrix, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%.4lf", matrix[i][j]);
            if (j < cols - 1) {
                    printf(",");
                }
            
        }
        printf("\n");
    }
}


int main(int argc, char *argv[])
{

    int n = 0;
    int d = 0; 
    double **vectorX;

    if (argc != 3)
    {
        printf("An Error Has Occurred\n");
        return 1;
    }

    const char *goal = argv[1];
    const char *file_name = argv[2];

    enum GoalType goalType = get_goal_type(goal);

    if (goalType == INVALID)
    {
        printf("An Error Has Occurred\n");
        return 1;
    }
    else
    {
        int i, j, temp; 

        FILE *file = fopen(file_name, "r");
    if (file == NULL) {
        return 1;
    }

    while ((temp = fgetc(file)) != EOF) {
        if (temp == '\n') {
            n++;
        } else if ((temp == ',') && (n == 0)) {
            d++;
        }
    }
    d++;

    rewind(file);

    vectorX = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        vectorX[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            if (fscanf(file, "%lf,", &vectorX[i][j]) != 1) {
                return 1;
            }
        }
    }

    fclose(file);




        if (goalType == SYM)
        {
            double **A_similarity_matrix = sym_compute(vectorX, n, d);


            print_matrix(A_similarity_matrix, n, n);


            for (int i = 0; i < n; i++)
            {
                free(A_similarity_matrix[i]);
            }
            free(A_similarity_matrix);
        }

        else {if (goalType == DDG)
        {
            double **A_similarity_matrix = sym_compute(vectorX, n, d);
            double **D_diagonal_degree_matrix = ddg_compute(A_similarity_matrix, n);


            print_matrix(D_diagonal_degree_matrix, n, n);


            for (int i = 0; i < n; i++)
            {
                free(A_similarity_matrix[i]);
            }
            free(A_similarity_matrix);


            for (int i = 0; i < n; i++)
            {
                free(D_diagonal_degree_matrix[i]);
            }
            free(D_diagonal_degree_matrix);
        }

        else {if (goalType == NORM)
        {
            double **A_similarity_matrix = sym_compute(vectorX, n, d);
            double **D_diagonal_degree_matrix = ddg_compute(A_similarity_matrix, n);
            double **W_graph_laplacian_matrix = norm_compute(A_similarity_matrix, n, D_diagonal_degree_matrix);


            print_matrix(W_graph_laplacian_matrix, n, n);


            for (int i = 0; i < n; i++)
            {
                free(A_similarity_matrix[i]);
            }
            free(A_similarity_matrix);


            for (int i = 0; i < n; i++)
            {
                free(D_diagonal_degree_matrix[i]);
            }
            free(D_diagonal_degree_matrix);


            for (int i = 0; i < n; i++)
            {
                free(W_graph_laplacian_matrix[i]);
            }
            free(W_graph_laplacian_matrix);
        }}}


        for (int i = 0; i < n; i++)
        {
            free(vectorX[i]);
        }
        free(vectorX);
    }

    return 0;
}