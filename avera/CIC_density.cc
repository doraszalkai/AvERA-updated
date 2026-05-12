#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "global_variables.h"

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

extern int empty_cells;
extern int non_empty_cells;

// A fv neve marad DTFE_density, hogy a step.cc megkapja, de valójában CIC-et számol
void DTFE_density(REAL** x)
{
    printf("CIC density estimation starting...\n");
    REAL start_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
    REAL omp_start_time = omp_get_wtime();

    unsigned N_un = (unsigned) N;
    double M_new = (double) M/rho_crit*pow(a_max/a_start, 3.0);

    // --- SZAPUDI CIC KÓDJA KEZDŐDIK ---
    int DENSITY_CELLS3 = pow(DENSITY_CELLS, 3);
    double cell_size = L / DENSITY_CELLS;

    // Initialize density grid
    for(int i = 0; i < DENSITY_CELLS3; i++)
        RHO[i] = 0.0;

    // CIC: assign each particle's mass to 8 nearest grid points
    for(size_t i = 0; i < N_un; i++)
    {
        double x_norm = x[i][0] / cell_size;
        double y_norm = x[i][1] / cell_size;
        double z_norm = x[i][2] / cell_size;
        
        int ix = (int)floor(x_norm);
        int iy = (int)floor(y_norm);
        int iz = (int)floor(z_norm);
        
        double dx = x_norm - ix;
        double dy = y_norm - iy;
        double dz = z_norm - iz;
        
        // Distribute mass to 8 corners of the cell
        for(int jx = 0; jx < 2; jx++)
        {
            for(int jy = 0; jy < 2; jy++)
            {
                for(int jz = 0; jz < 2; jz++)
                {
                    int nx = (ix + jx + DENSITY_CELLS) % DENSITY_CELLS;
                    int ny = (iy + jy + DENSITY_CELLS) % DENSITY_CELLS;
                    int nz = (iz + jz + DENSITY_CELLS) % DENSITY_CELLS;
                    
                    double wx = (jx == 0) ? (1.0 - dx) : dx;
                    double wy = (jy == 0) ? (1.0 - dy) : dy;
                    double wz = (jz == 0) ? (1.0 - dz) : dz;
                    
                    int cell_idx = nx * DENSITY_CELLS * DENSITY_CELLS + ny * DENSITY_CELLS + nz;
                    RHO[cell_idx] += M_new * wx * wy * wz;
                }
            }
        }
    }

    // Convert to density (mass per cell volume)
    double cell_volume = pow(cell_size, 3);
    int rho_min = 0;
    int rho_max = 0;

    for(int i = 0; i < DENSITY_CELLS3; i++) {
        RHO[i] /= cell_volume;
        if(RHO[rho_min] > RHO[i]) rho_min = i;
        if(RHO[rho_max] < RHO[i]) rho_max = i;
    }
    // --- SZAPUDI CIC KÓDJA VÉGE ---

    REAL end_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
    REAL omp_end_time = omp_get_wtime();

    printf("The minimal and the maximal density (in Omega_m):\n Rho_min = %g\t Rho_max = %g\t\n", RHO[rho_min],RHO[rho_max]);
    printf("...CIC done.\n");
    printf("CIC CPU time = %lfs\n", end_time-start_time);
    printf("CIC RUN time = %lfs\n", omp_end_time-omp_start_time);
    
    // Számolás Csabai kérésére
    empty_cells = 0;
    non_empty_cells = 0;

    for(int i = 0; i < DENSITY_CELLS3; i++) {
        if(RHO[i] > 1e-10) { // Egy nagyon kicsi küszöb a numerikus zaj ellen
            non_empty_cells++;
        } else {
            empty_cells++;
        }
    }
    
    double empty_percent = (double)empty_cells / DENSITY_CELLS3 * 100.0;
    
    // Kiírás a standard outputra (ami a szimulációs logba kerül)
    printf("--------------------------------------------------\n");
    printf("CIC CELL STATISTICS:\n");
    printf("Total cells: %d\n", DENSITY_CELLS3);
    printf("Empty cells: %d (%.2f%%)\n", empty_cells, empty_percent);
    printf("Non-empty cells: %d\n", non_empty_cells);
    printf("--------------------------------------------------\n");
    fflush(stdout); // Kikényszerítjük, hogy azonnal írja ki
    
    return;
}