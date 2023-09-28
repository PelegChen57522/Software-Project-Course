#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"


static PyObject *symnmf(PyObject *self, PyObject *args)
{
    PyObject* vectorX_pyObj;  
    PyObject* H_init_pyObj; 
    int n, k, d;         
    
    if (!PyArg_ParseTuple(args, "OOiii", &vectorX_pyObj,&H_init_pyObj,&k,&n,&d))
    {
        return NULL;
    }

    
    double **vectorX = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        vectorX[i] = (double *)malloc(d * sizeof(double));

        PyObject *datapoint_pyObj = PyList_GetItem(vectorX_pyObj, i);
        for (int j = 0; j < d; j++)
        {
            PyObject *coordinate_pyObj = PyList_GetItem(datapoint_pyObj, j);
            vectorX[i][j] = PyFloat_AsDouble(coordinate_pyObj);
        }
    }

    double **H_init = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        H_init[i] = (double *)malloc(k * sizeof(double));

        PyObject *H_init_line_pyObj = PyList_GetItem(H_init_pyObj, i);
        for (int j = 0; j < k; j++)
        {
            PyObject *coordinate_obj = PyList_GetItem(H_init_line_pyObj, j);
            H_init[i][j] = PyFloat_AsDouble(coordinate_obj);
        }
    }



    double **A_similarity_matrix = sym_compute(vectorX, n, d);


    double **D_diagonal_degree_matrix = ddg_compute(A_similarity_matrix, n);


    double **W_graph_laplacian_matrix = norm_compute(A_similarity_matrix, n, D_diagonal_degree_matrix);

    
    double** Final_H = compute_symnmf(H_init, W_graph_laplacian_matrix,n,k);


    if (Final_H == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;
    }



    PyObject *result = PyList_New(n);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        PyObject *row = PyList_New(k);
        if (!row) return NULL;

        for (int j = 0; j < k; j++) {
            PyObject *value = PyFloat_FromDouble(Final_H[i][j]);
            if (!value) return NULL;
            
            PyList_SetItem(row, j, value);
        }
        
        PyList_SetItem(result, i, row);
    }




    for (int i = 0; i < n; i++)
    {
        free(vectorX[i]);
        free(A_similarity_matrix[i]);
        free(D_diagonal_degree_matrix[i]);
        free(W_graph_laplacian_matrix[i]);
        free(Final_H[i]);
    }
    free(vectorX);
    free(A_similarity_matrix);
    free(D_diagonal_degree_matrix);
    free(W_graph_laplacian_matrix);
    free(Final_H);

    return result;
}





static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *vectorX_pyObj;

    if (!PyArg_ParseTuple(args, "O", &vectorX_pyObj))
    {
        return NULL;
    }

    int n = PyList_Size(vectorX_pyObj);
    int d = PyList_Size(PyList_GetItem(vectorX_pyObj, 0));

    double **vectorX = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        vectorX[i] = (double *)malloc(d * sizeof(double));

        PyObject *datapoint_pyObj = PyList_GetItem(vectorX_pyObj, i);
        for (int j = 0; j < d; j++)
        {
            PyObject *coordinate_pyObj = PyList_GetItem(datapoint_pyObj, j);
            vectorX[i][j] = PyFloat_AsDouble(coordinate_pyObj);
        }
    }

    double **A_similarity_matrix = sym_compute(vectorX, n, d);

    if (A_similarity_matrix == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;
    }

    PyObject *result = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *row = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            PyList_SetItem(row, j, PyFloat_FromDouble(A_similarity_matrix[i][j]));
        }
        PyList_SetItem(result, i, row);
    }


    for (int i = 0; i < n; i++)
    {
        free(vectorX[i]);
        free(A_similarity_matrix[i]);
    }
    free(vectorX);
    free(A_similarity_matrix);

    return result;
}


static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *vectorX_pyObj;

    if (!PyArg_ParseTuple(args, "O", &vectorX_pyObj))
    {
        return NULL;
    }

    int n = PyList_Size(vectorX_pyObj);
    int d = PyList_Size(PyList_GetItem(vectorX_pyObj, 0));

    double **vectorX = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        vectorX[i] = (double *)malloc(d * sizeof(double));

        PyObject *datapoint_pyObj = PyList_GetItem(vectorX_pyObj, i);
        for (int j = 0; j < d; j++)
        {
            PyObject *coordinate_pyObj = PyList_GetItem(datapoint_pyObj, j);
            vectorX[i][j] = PyFloat_AsDouble(coordinate_pyObj);
        }
    }

    double **A_similarity_matrix = sym_compute(vectorX, n, d);

    double **D_diagonal_degree_matrix = ddg_compute(A_similarity_matrix, n);


    if (D_diagonal_degree_matrix == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;
    }

    PyObject *result = PyList_New(n);
    if (!result) return NULL;

    for (int i = 0; i < n; i++) {
        PyObject *row = PyList_New(n);
        if (!row) return NULL;

        for (int j = 0; j < n; j++) {
            PyObject *value = PyFloat_FromDouble(D_diagonal_degree_matrix[i][j]);
            if (!value) return NULL;
            
            PyList_SetItem(row, j, value);
        }
        
        PyList_SetItem(result, i, row);
    }




    for (int i = 0; i < n; i++)
    {
        free(vectorX[i]);
        free(A_similarity_matrix[i]);
        free(D_diagonal_degree_matrix[i]);
    }
    free(vectorX);
    free(A_similarity_matrix);
    free(D_diagonal_degree_matrix);

    return result;
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *vectorX_pyObj;

    if (!PyArg_ParseTuple(args, "O", &vectorX_pyObj))
    {
        return NULL;
    }

    int n = PyList_Size(vectorX_pyObj);
    int d = PyList_Size(PyList_GetItem(vectorX_pyObj, 0));

    double **vectorX = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        vectorX[i] = (double *)malloc(d * sizeof(double));

        PyObject *datapoint_pyObj = PyList_GetItem(vectorX_pyObj, i);
        for (int j = 0; j < d; j++)
        {
            PyObject *coordinate_pyObj = PyList_GetItem(datapoint_pyObj, j);
            vectorX[i][j] = PyFloat_AsDouble(coordinate_pyObj);
        }
    }

    double **A_similarity_matrix = sym_compute(vectorX, n, d);

    double **D_diagonal_degree_matrix = ddg_compute(A_similarity_matrix, n);

    double **W_graph_laplacian_matrix = norm_compute(A_similarity_matrix, n, D_diagonal_degree_matrix);

    if (W_graph_laplacian_matrix == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;
    }

    PyObject *result = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *row = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            PyList_SetItem(row, j, PyFloat_FromDouble(W_graph_laplacian_matrix[i][j]));
        }
        PyList_SetItem(result, i, row);
    }

    for (int i = 0; i < n; i++)
    {
        free(vectorX[i]);
        free(A_similarity_matrix[i]);
        free(D_diagonal_degree_matrix[i]);
        free(W_graph_laplacian_matrix[i]);
    }
    free(vectorX);
    free(A_similarity_matrix);
    free(D_diagonal_degree_matrix);
    free(W_graph_laplacian_matrix);

    return result;
}



static PyObject *analysis(PyObject *self, PyObject *args)
{
    PyObject *H_pyObj;
    int n,k;

    if (!PyArg_ParseTuple(args, "Oii", &H_pyObj, &n, &k))
    {
        return NULL;
    }


    double **H = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        H[i] = (double *)malloc(k * sizeof(double));

        PyObject *row_pyObj = PyList_GetItem(H_pyObj, i);
        for (int j = 0; j < k; j++)
        {
            PyObject *coordinate_pyObj = PyList_GetItem(row_pyObj, j);
            H[i][j] = PyFloat_AsDouble(coordinate_pyObj);
        }
    }

    int *res_analysis = analysis_compute(H,n,k);

    if (res_analysis == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;
    }

    PyObject *result = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyList_SetItem(result, i, PyLong_FromLong(res_analysis[i]));
    }

    for (int i = 0; i < n; i++)
    {
        free(H[i]);
    }
    free(H);
    free(res_analysis);

    return result;
}



static PyMethodDef methods[] = {
    {"symnmf", symnmf, METH_VARARGS, "symnmf"},
    {"sym", sym, METH_VARARGS, "Compute the similarity matrix"},
    {"ddg", ddg, METH_VARARGS, "Compute the diagonal degree matrix"},
    {"norm", norm, METH_VARARGS, "Compute the normalized similarity matrix"},
    {"analysis", analysis, METH_VARARGS, "Compute analysis"},
    {NULL, NULL, 0, NULL} 
};

static struct PyModuleDef symnmf_module = {
    PyModuleDef_HEAD_INIT,
    "symnmfModule",
    NULL,
    -1,
    methods};

PyMODINIT_FUNC PyInit_symnmfModule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmf_module);
    if (!m)
    {
        return NULL;
    }
    return m;
}