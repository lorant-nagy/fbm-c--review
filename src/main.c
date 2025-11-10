// main.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../include/fbm_patched.h"

int main(void) {
    
    srand((unsigned)time(NULL));

    char filename[100];

    int target_path_count = 1000;

    double H = 0.8;
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
        char dir[100];
        sprintf(dir, "data/fbm_H%0.1f_T%0.1f_n%d", H, T, n);
        
        // Create directory if it doesn't exist
        struct stat st = {0};
        if (stat(dir, &st) == -1) {
            if (mkdir(dir, 0755) != 0) {
                perror("mkdir");
                free(path);
                free(times);
                return 1;
            }
        }
        
        sprintf(filename, "data/fbm_H%0.1f_T%0.1f_n%d/fbm_%d.csv", H, T, n, counter);
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
