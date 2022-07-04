

// The size of linear system is defined by N_cons, N_var in gaussianLib.h
// module load compilers/intel && module load mpi/intel-mpi/2019u5/bin && cd COMP328/Assignment/compare
// mpiicc -O0 gaussian_MPI.c gaussianLib.c -o gaussian.exe && mpirun -np 2 ./gaussian.exe 2


#pragma warning(disable : 1368)
#pragma warning(disable : 167)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "gaussianLib.h"

double A[N_cons][N_var];
// The constant values
double b[N_var][1];
// The augmented matrix
double M[N_cons][N_var + 1]; // [N_row][N_col] may cause error

int my_rank=0, comm_size=0;

int root = 0;

double start, finish, elapsed;


void Gaussian_MPI(double M[N_row][N_col], int* rank_A, int* rank_M);

int main(int argc, char* argv[]){

    if(argc==2){
        debug = atoi(argv[1]);
    }

    int rank_M=0, rank_A=0;
        
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
    if(N_cons % comm_size != 0){
        if(my_rank == root){
            printf("!! N_cons %d is indivisible by comm_size %d , please change the number of processors !!\n", N_cons, comm_size);
        }
        MPI_Finalize();
        return -1;
    }
    

    if(my_rank==root){
        load_linear_system(option, A, b);
        form_M(M,A,b);
    }
   
    Gaussian_MPI(M, &rank_A, &rank_M);
    
    if(my_rank==root){
        show_solution(M, N_var, rank_A, rank_M);
        if(debug>=-1){
            printf("Run time for Gaussian elimination of by MPI: %f milliseconds\n", (finish-start)*1000.0);
        }
    }
    
    MPI_Finalize();
    
    return 0;
}


void Gaussian_MPI(double M[N_row][N_col], int* rank_A, int* rank_M){
/*  
    find the left-most non-zero elements(pivot),
    swap pivot row to the top of unprocessed rows (the line under the last processed Row-Echelon part), 
    reduce pivot to one (pivot row devided by pivot), 
    set other elements (except pivot) in pivot row to 0, all other rows minus some multiple of pivot row,
    repeat above until the last column has been checked.
        
*/
    
    int part_size = N_cons * (N_var+1) /comm_size;
    int row_range = N_cons / comm_size;
    int start_row = my_rank * row_range;
    int end_row = (my_rank+1) * row_range - 1; // my rows in [start_row, end_row]
    
    double M_seg[N_cons][N_var+1];

    double pivot_line[N_var+1];
    
    double temp_line[N_var+1];
    
    int row_flag = 0;
    
    int local_flag = 0, local_row = 0;
    
    int find_pivot = 0;
    
    int target_processor = 0, row_processor = 0;
    
    int col, row, i,j;
    
    
    start = MPI_Wtime();
    
    // root scatter M to all the procesors's M_seg
    MPI_Scatter(M, part_size, MPI_DOUBLE, M_seg, part_size, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    // column by column (left to right), each column top to down, to find the pivot of current stage
    for(col = 0; col < N_col; ++col){
        for(row = row_flag; row < N_row; ++row){ 
            
            // the processor having the "row"th row
            row_processor = (int) floor(((double)row / (double)row_range));
            // printf("col %d, row_flag: %d row %d -- row_processor: %d\n", col, row_flag, row, row_processor);
        
            // whether the target pivot is in this row (of the col)
            find_pivot = 0;
            
            // if I have the row (implemented by only one processor), check whether the element is 0, not 0 -> find pivot
            if(my_rank == row_processor){
                // map the index in M to the index in M_seg
                local_row = row % row_range;
                // check wether it is the pivot role (<=> pivot nonzero)
                if(M_seg[local_row][col] != 0){
                    // change the find_pivot flag 
                    find_pivot = 1; 
                }
            }
            
            // everybody updates the find_pivot flag;
            MPI_Bcast(&find_pivot, 1, MPI_INT, row_processor, MPI_COMM_WORLD);
            
            // If we find the pivot role, do the computation, and update row_flag, other wise continue to search the nex row
            if(find_pivot == 1){
                // the processor that the pivot row should be in (after swapping)
                target_processor = (int) floor(((double)row_flag / (double)row_range));
                // map the index in M to the index in M_seg
                local_flag = row_flag % row_range;
                
                // swap the "row"th row and the "row_flag"th row
                // exchange data using sendrecv
                if(row != row_flag){
                    if(target_processor == row_processor){
                        if(my_rank == target_processor){
                            swap_two_rows(M_seg, local_flag, local_row);
                            // print_M(M_seg);
                        }
                    }else{
                        if(my_rank == target_processor){
                            // MPI_Send(M_seg[local_flag], N_col, MPI_DOUBLE, row_processor, 0,  MPI_COMM_WORLD); // tag 0: "row_flag"th row
                            // MPI_Recv(M_seg[local_flag], N_col, MPI_DOUBLE, row_processor, 1,  MPI_COMM_WORLD, MPI_STATUS_IGNORE); // tag 1: "row"th row
                            
                             MPI_Sendrecv(M_seg[local_flag], N_col, MPI_DOUBLE, row_processor, 0,
                                          M_seg[local_flag], N_col, MPI_DOUBLE, row_processor, 1,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            
                        }
                        if(my_rank == row_processor){
                            // MPI_Send(M_seg[local_row], N_col, MPI_DOUBLE, target_processor, 1,  MPI_COMM_WORLD);
                            // MPI_Recv(M_seg[local_row], N_col, MPI_DOUBLE, target_processor, 0,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            
                            MPI_Sendrecv(M_seg[local_row], N_col, MPI_DOUBLE, target_processor, 1,
                                         M_seg[local_row], N_col, MPI_DOUBLE, target_processor, 0,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                    }
                    
                }
                    
                // change the pivot to one
                if(my_rank == target_processor){
                    for(j=N_col-1; j>=col; --j){ // pivot itself cannot be changed to one before other elements being divided by it
                        M_seg[local_flag][j] = M_seg[local_flag][j] / M_seg[local_flag][col];
                    }
                    memcpy(pivot_line, M_seg[local_flag], N_col * sizeof(double));
                }
                // broad cast pivot line to everybody's pivot_line
                MPI_Bcast(pivot_line, N_col, MPI_DOUBLE, target_processor, MPI_COMM_WORLD); 
                
                // change the other elements in pivot column to 0
                for(i=0; i<row_range; ++i){
                    if((my_rank * row_range + i != row_flag) && (M_seg[i][col] != 0)){ // no need to change target pivot row, or rows with target pivot column element being 0 
                        for(j=N_col-1; j>=col; --j){ // M[i][col] cannot be changed before being reduced by other elements of the same row
                            // printf("M[%d][%d] = M[%d][%d] - M[%d][%d]", i,j, i,col, row_flag,col); 
                            M_seg[i][j] = M_seg[i][j] -  M_seg[i][col] * pivot_line[j];
                        }
                    }
                }
                
                // update row_flag, rank of the Augmented Matrix M, rank of the Coefficient Matrix A
                ++ row_flag;
                *rank_M=row_flag;
                *rank_A= col!=N_col-1 ? row_flag:row_flag-1;
                 
                break;
                
            }     
        }
    }
    
    // root collects M_seg to form M in reduced row-echolen form
    MPI_Gather(M_seg, part_size, MPI_DOUBLE, M, part_size, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    finish = MPI_Wtime();
    
    if(my_rank == root && debug>=2){
        printf("\n");
        printf("Gaussian Elimination to the reduced row-echolen form:\n");
        print_M(M);
    }

}






















