// import the Python header
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// HW1 imports
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void update_centroids(double **centroids, double **data, int *cluster_assignment, int n, int k, int d)
{
    double *cluster_sizes = (double *)malloc(k * sizeof(double));
    int i, j;

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            centroids[i][j] = 0;
        }
        cluster_sizes[i] = 0;
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            centroids[cluster_assignment[i]][j] += data[i][j];
        }
        cluster_sizes[cluster_assignment[i]]++;
    }

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            centroids[i][j] /= cluster_sizes[i];
        }
    }

    free(cluster_sizes);
}

double euclidean_distance(double *a, double *b, int d)
{
    double distance = 0.0;
    int i;
    for (i = 0; i < d; i++)
    {
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(distance);
}

// assigns data points to the nearest cluster centroid in a k-means clustering algorithm.
void assign_clusters(double **centroids, double **data, int *cluster_assignment, int n, int k, int d)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        double min_distance = 1e100;
        int min_cluster_id = -1;
        for (j = 0; j < k; j++)
        {
            double distance = euclidean_distance(data[i], centroids[j], d);
            if (distance < min_distance)
            {
                min_distance = distance;
                min_cluster_id = j;
            }
        }
        cluster_assignment[i] = min_cluster_id;
    }
}

PyObject * kmeans(double **centroids, double **data, int n, int k, int d, int max_iter, double epsilon)
{
    int *cluster_assignment = (int *)malloc(n * sizeof(int));
    int i, j, t;
    double **old_centroids;
    double max_change;

    for (t = 0; t < max_iter; t++)
    {

        assign_clusters(centroids, data, cluster_assignment, n, k, d);

        old_centroids = (double **)malloc(k * sizeof(double *));
        for (i = 0; i < k; i++)
        {
            old_centroids[i] = (double *)malloc(d * sizeof(double));
            for (j = 0; j < d; j++)
            {
                old_centroids[i][j] = centroids[i][j];
            }
        }

        update_centroids(centroids, data, cluster_assignment, n, k, d);

        max_change = 0.0;
        for (i = 0; i < k; i++)
        {
            double dist = euclidean_distance(centroids[i], old_centroids[i], d);
            if (dist > max_change)
            {
                max_change = dist;
            }
        }

        for (i = 0; i < k; i++)
        {
            free(old_centroids[i]);
        }
        free(old_centroids);

        if (max_change <= epsilon)
        {
            break;
        }
    }


    free(cluster_assignment);

    PyObject *final_centroids = PyList_New(k);
    for (i = 0; i < k; i++) {
        PyObject *centroid = PyList_New(d);
        for (j = 0; j < d; j++) {
            PyObject *coordinate = PyFloat_FromDouble(centroids[i][j]);
            PyList_SetItem(centroid, j, coordinate);
        }
        PyList_SetItem(final_centroids, i, centroid);
    }

    return final_centroids;

}




static PyObject *fit(PyObject *self, PyObject *args)
{

    int k ;
    int max_iter;
    int i;
    double epsilon;


    PyObject *datapoints_pyObj;
    PyObject *centroids_pyObj;

    if (!PyArg_ParseTuple(args, "OOiid", &datapoints_pyObj, &centroids_pyObj, &k, &max_iter, &epsilon))
    {
        return NULL;
    }


    int n = PyList_Size(datapoints_pyObj);
    int d = PyList_Size(PyList_GetItem(datapoints_pyObj, 0));

    double **datapoints = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        datapoints[i] = (double *)malloc(d * sizeof(double));

        PyObject *datapoint_pyObj = PyList_GetItem(datapoints_pyObj, i);
        for (int j = 0; j < d; j++)
        {
            PyObject *coordinate_pyObj = PyList_GetItem(datapoint_pyObj, j);
            datapoints[i][j] = PyFloat_AsDouble(coordinate_pyObj);
        }
    }


    double **centroids = (double **)malloc(k * sizeof(double *));
    for (int i = 0; i < k; i++)
    {
        centroids[i] = (double *)malloc(d * sizeof(double));

        PyObject *centroid_pyObj = PyList_GetItem(centroids_pyObj, i);
        for (int j = 0; j < d; j++)
        {
            PyObject *coordinate_obj = PyList_GetItem(centroid_pyObj, j);
            centroids[i][j] = PyFloat_AsDouble(coordinate_obj);
        }
    }



    PyObject *final_centroids = kmeans(centroids, datapoints, n, k, d, max_iter, epsilon);
    if (final_centroids == NULL) {
        printf("An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        free(datapoints[i]);
    }
    free(datapoints);

    for (i = 0; i < k; i++)
    {
        free(centroids[i]);
    }
    free(centroids);

    return Py_BuildValue("O", final_centroids);

}

// PyMethodDef
// In order to call the methods defined in your module, we will need to tell the Python interpreter about them first so we wil use PyMethodDef
static PyMethodDef kmeans_methods[] = {
    {"fit", // The name the user would write to invoke this particular function in Python.
     fit,
     METH_VARARGS,                              // flags indicating parameters accepted for this function
     PyDoc_STR("Perform k-means clustering.")}, // The docstring for the function
    {NULL, NULL, 0, NULL}                       // The last entry must be all NULL as shown to act as a sentinel. Python looks for this entry to know that all of the functions for the module have been defined.
};

// PyModuleDef
// This initiates the module using the above definitions.
static struct PyModuleDef kmeans_module = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    "K-means clustering module",
    -1,
    kmeans_methods};

// PyMODINIT_FUNC
// When a Python program imports your module for the first time, it will call PyInit_mykmeanssp(void):
PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeans_module);
    if (!m)
    {
        return NULL;
    }
    return m;
}
