#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void update_centroids(double **centroids, double **data, int *cluster_assignment, int n, int k, int d) {
    double *cluster_sizes = (double *)malloc(k * sizeof(double));
    int i, j;

    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            centroids[i][j] = 0;
        }
        cluster_sizes[i] = 0;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            centroids[cluster_assignment[i]][j] += data[i][j];
        }
        cluster_sizes[cluster_assignment[i]]++;
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            centroids[i][j] /= cluster_sizes[i];
        }
    }

    free(cluster_sizes);
}

double euclidean_distance(double *a, double *b, int d) {
    double distance = 0.0;
    int i;
    for (i = 0; i < d; i++) {
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(distance);
}

void assign_clusters(double **centroids, double **data, int *cluster_assignment, int n, int k, int d) {
    int i, j;
    for (i = 0; i < n; i++) {
        double min_distance = 1e100;
        int min_cluster_id = -1;
        for (j = 0; j < k; j++) {
            double distance = euclidean_distance(data[i], centroids[j], d);
            if (distance < min_distance) {
                min_distance = distance;
                min_cluster_id = j;
            }
        }
        cluster_assignment[i] = min_cluster_id;
    }
}

void kmeans(double **centroids, double **data, int n, int k, int d, int max_iter, double epsilon) {
    int *cluster_assignment = (int *)malloc(n * sizeof(int));
    int i, j, t;
    double **old_centroids;
    double max_change;

    for (t = 0; t < max_iter; t++) {


        assign_clusters(centroids, data, cluster_assignment, n, k, d);

        old_centroids = (double **)malloc(k * sizeof(double *));
        for (i = 0; i < k; i++) {
            old_centroids[i] = (double *)malloc(d * sizeof(double));
            for (j = 0; j < d; j++) {
                old_centroids[i][j] = centroids[i][j];
            }
        }

        update_centroids(centroids, data, cluster_assignment, n, k, d);

        max_change = 0.0;
        for (i = 0; i < k; i++) {
            double dist = euclidean_distance(centroids[i], old_centroids[i], d);
            if (dist > max_change) {
                max_change = dist;
            }
        }

        for (i = 0; i < k; i++) {
            free(old_centroids[i]);
        }
        free(old_centroids);

        if (max_change <= epsilon) {
            break;
        }
    }

    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            printf("%.4f", centroids[i][j]);
            if (j < d - 1) {
                printf(",");
            }
        }
        printf("\n");
    }

    free(cluster_assignment);
}


int main(int argc, char *argv[]) {
    
    int k, max_iter, n, d, i, j;
    double epsilon;
    double **data, **centroids;
    char ch;
    char *endptr;
    char *endptr2;

    if (argc < 2 || argc > 3) {
        fprintf(stderr, "An Error Has Occurred\n");
        return 1;
    }

    k = strtol(argv[1], &endptr, 10);
    if (*endptr != '\0' || k <= 1) {
        fprintf(stderr, "Invalid number of clusters!\n");
        return 1;
    }

    if (argc == 3) {
        max_iter = strtol(argv[2], &endptr2, 10);
        if (*endptr2 != '\0' || max_iter <= 1 || max_iter >= 1000) {
            fprintf(stderr, "Invalid maximum iteration!\n");
            return 1;
        }
    } else {
        max_iter = 200;
    }



    epsilon = 0.001;

    /* Calculate dimensions and number of data points */
    n = 0;
    d = 0;
    while ((ch = getchar()) != EOF) {
        if (ch == '\n') {
            n++;
        } else if (n == 0 && ch == ',') {
            d++;
        }
    }
    d++; /* Add 1 because the last dimension is not followed by a comma */


    rewind(stdin); /* Reset file pointer to the beginning of the input */

    data = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        data[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            if (scanf("%lf,", &data[i][j]) != 1) {
                fprintf(stderr, "An Error Has Occurred\n");
                return 1;
            }
        }
    }

    centroids = (double **)malloc(k * sizeof(double *));
    for (i = 0; i < k; i++) {
        centroids[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            centroids[i][j] = data[i][j];
        }
    }

    kmeans(centroids, data, n, k, d, max_iter, epsilon);

    for (i = 0; i < n; i++) {
        free(data[i]);
    }
    free(data);

    for (i = 0; i < k; i++) {
        free(centroids[i]);
    }
    free(centroids);

    return 0;
}
