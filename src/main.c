// main.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/fbm.h"

int main(void) {
    
    srand((unsigned)time(NULL));

    char filename[50];

    int target_path_count = 5;

    double H = 0.7;
    int    n = 128;
    double T = 128.0;
    int counter = 1;
    double *path = NULL;
    double *times = malloc((n + 1) * sizeof *times);

    for (int i = 0; i <= n; ++i) {
        double t = T * (double)i / (double)n;
        times[i] = (double)(t);
    }

    while (counter <= target_path_count) {
        
        path = simulate_fBm(H, n, T);
        // Save full path to CSV (columns: t, B_H(t))
        // create a variable holding the file name
        sprintf(filename, "data/fbm_%d.csv", counter);
        FILE *fp = fopen(filename, "w");
        if (!fp) { perror("fopen"); free(path); free(times); return 1; }
        fprintf(fp, "t,BH\n");
        for (int i = 0; i <= n; ++i) {
            double t = T * (double)i / (double)n;
            fprintf(fp, "%.17g,%.17g\n", t, path[i]);
        }
        fclose(fp);

        counter++;

    }

    free(path);
    free(times);
    return 0;
}
