// main.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/fbm.h"

int main(void) {
    
    srand((unsigned)time(NULL));

    double H = 0.7;
    int    n = 128;
    double T = 128.0;
    double *path = NULL;
    double *times = malloc((n + 1) * sizeof *times);

    path = simulate_fBm(H, n, T);

    for (int i = 0; i <= n; ++i) {
        double t = T * (double)i / (double)n;
        times[i] = (double)(t);
    }

    // print times array elements
    for (int i = 0; i <= n; ++i) {
        printf("times[%d] = %.17g\n", i, times[i]);
    }

    // Save full path to CSV (columns: t, B_H(t))
    FILE *fp = fopen("fbm.csv", "w");
    if (!fp) { perror("fopen"); free(path); return 1; }
    fprintf(fp, "t,BH\n");
    for (int i = 0; i <= n; ++i) {
        double t = T * (double)i / (double)n;
        fprintf(fp, "%.17g,%.17g\n", t, path[i]);
    }
    fclose(fp);

    free(path);
    free(times);
    return 0;
}
