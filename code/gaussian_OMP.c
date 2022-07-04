
// The size of linear system is defined by N_cons, N_var in gaussianLib.h
// module load compilers/intel && module load mpi/intel-mpi/2019u5/bin && cd COMP328/Assignment
// icc -qopenmp -O0 gaussian_OMP.c gaussianLib.c -o gaussian_OMP.exe && ./gaussian_OMP.exe 2


#pragma warning(disable : 1368)

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "gaussianLib.h"

// variables
// char** x[N_var][2];
// The coefficient matrix
double A[N_cons][N_var];
// The constant values
double b[N_var][1];
// The augmented matrix
double M[N_cons][N_var + 1];

int rank_M=0, rank_A=0;

double start, end, duration;


int main(int argc, char* argv[]){
    
    if(argc==2){
        debug = atoi(argv[1]);
    }
    
    load_linear_system(option, A, b);
    form_M(M,A,b);
    Gaussian_OMP(M, &rank_A, &rank_M);    
    show_solution(M, N_var, rank_A, rank_M);
    
    duration = end-start;
    if(debug>=-1){
        printf("Run time for Gaussian elimination openMP: %f milliseconds\n", duration * 1000.0);
    }
    
}

void Gaussian_OMP(double M[N_row][N_col], int* rank_A, int* rank_M){
/*  
    find the left-most non-zero elements(pivot),
    swap pivot row to the top of unprocessed rows (the line under the last processed Row-Echelon part), 
    reduce pivot to one (pivot row devided by pivot), 
    set other elements (except pivot) in pivot row to 0, all other rows minus some multiple of pivot row,
    repeat above until the last column has been checked.
        
*/
    
    int row_flag = 0, col = 0, row = 0, j = 0, i = 0;
    
    // clock_gettime(CLOCK_REALTIME, &start2); 
    start=omp_get_wtime();

    // find the left-most non-zero elements(pivot), swap rows, reduce pivot to one, same col other ele=0
    for(col = 0; col < N_col; ++col){
        for(row = row_flag; row < N_row; ++row){ // column by column, each column top to down
            if(M[row][col] != 0.0) { // find the left-most non-zero elements under the row_flag row
                
                // swap the target pivot row to the row_flag row
                swap_two_rows(M, row, row_flag); 
                
                // change the target pivot to 1
                double pivot = M[row_flag][col];
                // #pragma omp parallel for num_threads(20) default(none) private(j) shared(M, N_col, col, row_flag, pivot) // delete num_threads when using bash //越并越多
                for(j=N_col-1; j>=col; --j){ // pivot itself cannot be changed to one before other elements being divided by it
                    M[row_flag][j] = M[row_flag][j] / pivot;
                    // printf("%d thread %d handles j=%d\n", row_flag, omp_get_thread_num(), j);
                }
                
                // printf("\n");
                // change the elements in the target pivot column to 0
                #pragma omp parallel for default(none) private(i, j) shared(M, N_col, col, row_flag, N_row)
                // #pragma omp parallel for num_threads(2) default(none) private(i, j) shared(M, N_col, col, row_flag, N_row)
                for(i=0; i<N_row; ++i){
                    if(i!=row_flag &&  M[i][col]!=0){ // no need to change target pivot row, or rows with target pivot column element being 0 
                        for(j=N_col-1; j>=col; --j){ // M[i][col] cannot be changed before being reduced by other elements of the same row
                          // printf("%d thread %d: M[%d][%d] = M[%d][%d] - M[%d][%d]\n", row_flag, omp_get_thread_num(), i,j, i,col, row_flag,col); 
                          M[i][j] = M[i][j] -  M[i][col] * M[row_flag][j];
                        }
                    }
                }
                
                ++ row_flag;
                *rank_M=row_flag;
                *rank_A= col!=N_col-1 ? row_flag:row_flag-1;
                // return;
                break;
            }
        }
        
    }
   
    end=omp_get_wtime();
    
    if(debug>=2){
        printf("\n");
        printf("Gaussian Elimination to the reduced row-echolen form:\n");
        print_M(M);
    }
}