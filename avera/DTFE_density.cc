#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include<DTFE.h>
#include "global_variables.h"

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

//density field estimation with DTFE method
void DTFE_density(REAL** x)
{
    printf("DTFE starting...\n");
    fflush(stdout); // Kikényszerítjük a kiírást!

    //időmérés
    REAL start_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
    REAL omp_start_time = omp_get_wtime();

    char NO_grid[64]; sprintf(NO_grid, "%i", DENSITY_CELLS);
    char L_str[64]; sprintf(L_str, "%.13f", L);
    char box_min_str[64]; sprintf(box_min_str, "0.0");
    char grid_str[64]; sprintf(grid_str, "--grid");
    char periodic_str[64]; sprintf(periodic_str, "--periodic");
    char box_str[64]; sprintf(box_str, "--box");
    char verbose_str[64]; sprintf(verbose_str, "--verbose");
    char DTFE_str[64]; sprintf(DTFE_str, "DTFE");
    char one_str[64]; sprintf(one_str, "1");

    int argc;
    char *argv_DTFE[] = {DTFE_str, grid_str, NO_grid, periodic_str, box_str, box_min_str, L_str, box_min_str, L_str, box_min_str, L_str, verbose_str, one_str};
    char sph_str[64]; sprintf(sph_str, "--SPH");
    char neighbors_str[64]; sprintf(neighbors_str, "40");
    char *argv_SPH[] = {DTFE_str, grid_str, NO_grid, periodic_str, box_str, box_min_str, L_str, box_min_str, L_str, box_min_str, L_str, verbose_str, one_str, sph_str, neighbors_str};

    User_options userOptions;
    
    printf("[DEBUG 1] Opciok beolvasasa...\n"); 
    fflush(stdout);

    if(NONISOTROPIC_EXPANSION == 3) {
        argc = 13;
        userOptions.readOptions(argc, argv_DTFE, false, true);
    } else {
        printf("SPH density estimation...\n"); fflush(stdout);
        sprintf(one_str, "1");
        argc = 15;
        userOptions.readOptions(argc, argv_SPH, false, true);
    }

    printf("[DEBUG 2] Opciok rendben. Vektorok elokeszitese...\n"); 
    fflush(stdout);

    std::vector<Particle_data> particles;
    std::vector<Sample_point> samplingCoordinates;
    Quantities uquantities;
    Quantities aquantities;

    unsigned N_un = (unsigned) N;
    printf("[DEBUG 3] Reszecskek szama (N): %u\n", N_un); 
    fflush(stdout);

    particles.reserve(N);
    double M_new = (double) M/rho_crit*pow(a_max/a_start, 3.0);

    printf("[DEBUG 4] Reszecskek masolasa indul az x tombbul...\n"); 
    fflush(stdout);

    for (size_t i=0; i<N_un; ++i) {
        Particle_data temp;
        for (int j=0; j<3; ++j) {
            temp.position(j) = (double) x[i][j];
        }
        temp.weight() = M_new;
        particles.push_back(temp);
        
        // Csekkoljuk az elsőt és az utolsót:
        if (i == 0) { printf("[DEBUG 4.1] Elso reszecske atmasolva.\n"); fflush(stdout); }
        if (i == N_un - 1) { printf("[DEBUG 4.2] Utolso reszecske atmasolva.\n"); fflush(stdout); }
    }

    printf("[DEBUG 5] Vektor epites kesz. DTFE FOFUGGVENY HIVASA...\n"); 
    fflush(stdout);

    // ITT A NAGY UGRÁS
    DTFE(&particles, samplingCoordinates, userOptions, &uquantities, &aquantities);

    printf("[DEBUG 6] DTFE fofuggveny tuléltük! Suruseg szamolasa...\n"); 
    fflush(stdout);

    int DENSITY_CELLS3 = pow(DENSITY_CELLS,3);
    int rho_min = 0;
    int rho_max = 0;
    for(int i=0; i<DENSITY_CELLS3; i++) {
        RHO[i] = ((double) aquantities.density[i]);
        if(RHO[rho_min] > RHO[i]) rho_min = i;
        if(RHO[rho_max] < RHO[i]) rho_max = i;
    }

    REAL end_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
    REAL omp_end_time = omp_get_wtime();

    printf("The minimal and the maximal density (in Omega_m):\n Rho_min = %g\t Rho_max = %g\t\n", RHO[rho_min],RHO[rho_max]);
    printf("...DTFE done.\n"); fflush(stdout);

    if(NONISOTROPIC_EXPANSION == 3) {
        printf("DTFE CPU time = %lfs\n", end_time-start_time);
        printf("DTFE RUN time = %lfs\n", omp_end_time-omp_start_time);
    } else {
        printf("SPH CPU time = %lfs\n", end_time-start_time);
        printf("SPH RUN time = %lfs\n", omp_end_time-omp_start_time);
    }
    return;
}