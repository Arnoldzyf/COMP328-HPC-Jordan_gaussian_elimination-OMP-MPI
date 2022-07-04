#pragma warning(disable : 1368)
#pragma warning(disable : 167)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gaussianLib.h"

const int N_row = N_cons, N_col = N_var + 1;

int option=0;

int debug=-1;

void print_M_seg(double M_seg[][N_col], int N_row){
    int i, j;
	  for(i=0; i<N_row; ++i){ //row
		    for(j=0; j<N_col; ++j){ //column
			      printf("\t %.2f", M_seg[i][j]);
		    }
        printf("\n");
	  }
}

void print_M(double M[N_row][N_col]){ // !!!!  double **X  may cause core dump //can only handle one size??
	  int i, j;
	  for(i=0; i<N_row; ++i){ //row
		    for(j=0; j<N_col; ++j){ //column
			      printf("\t %.2f", M[i][j]);
		    }
        printf("\n");
	  }
}


void fill_A(double A[N_cons][N_var]){
    srand(4);
	  int i,j;
	  for(i=0; i<N_cons; ++i){
        for(j=0; j<N_var; ++j){
		        A[i][j] = round(1000.0 * rand() / RAND_MAX);
        } 
	  }
}

void fill_b(double b[N_var][1]){
    srand(2);
    int i,j;
	  for(i=0; i<N_cons; ++i){
        b[i][0] = round(1000.0 * rand() / RAND_MAX);
	  }
}

void form_M(double M[N_row][N_col], double A[N_cons][N_var], double b[N_var][1]){
    int i,j;
    for(i=0; i<N_row; ++i){
        for(j=0; j<N_col; ++j){
            if(j!=N_col-1){
                M[i][j]=A[i][j];
            }else{
                M[i][j]=b[i][0];
            }
        }
    }
    if(debug>=2){
        printf("\n");
        printf("The resulting Augmented Matrix M:\n");
        print_M(M);
    }
    
}

void load_linear_system(int option, double A[N_cons][N_var], double b[N_var][1]){
    if(debug>=0){
        // printf("Load linear system Ax=b of case %d...\n", option);
        printf("size: N_cons, N_var = %d, %d\n", N_cons, N_var);
    }
   
    
    switch(option) {
        case 0:{ // random 4 4
            fill_A(A);
            fill_b(b);
            break;
        }
        case 1:{ // https://www.digitmath.com/gaussian-elimination.html // 4 4 unique: -1 2 1 3
            // printf("Unique solutions\n");
            double A_[N_cons][N_var] = {  // name as A will lead to problems
                {0, 1, 1, -2},
                {1, -4, -7, -1},
                {1, 2, -1, 0},
                {2, 4, 1, -3}
            };
            memcpy(A, A_, sizeof(A_)); 
            double b_[N_cons][1] = {-3, -19, 2, -2};
            memcpy(b, b_, sizeof(b_));
            break;
        }
        case 2:{ // gaussian.pdf // 4 4 infinite
            double A_[N_cons][N_var] = {
                {1, 1, -1, 4},
                {2, 2, -1, 7},
                {0, 0, 0, 0}, 
                {3, 3, -2, 11}
            };
            memcpy(A, A_, sizeof(A_));
            double b_[N_cons][1] = {2, 1, 0, 3};
            memcpy(b, b_, sizeof(b_));
            break;
        }
        case 3:{ // https://www2.math.upenn.edu/~rimmer/math240/gausselim.pdf  // 3 3 infeasible
            double A_[N_cons][N_var] = {
                {3, 6, 6},
                {3, -6, -3},
                {3, -2, 0}, 
            };
            memcpy(A, A_, sizeof(A_)); // !size of A will go wrong
            double b_[N_cons][1] = {5, 2, 1};
            memcpy(b, b_, sizeof(b_));
            break;
        }
        
    }
    
    if(debug>=2){
        printf("\n");
        int i,j;
        // print A
        printf("Coefficient Matrix A:\n");
        for(i=0; i<N_cons; ++i){
            for(j=0; j<N_var; ++j){
    		        printf("\t %.2f", A[i][j]);
            }
            printf("\n");
    	  }
        // print x
        printf("Variables x.T:\n");
        for(i=0; i<N_cons; ++i){
           printf("\t x%d", i+1);
    	  }
        printf("\n");
        // print b
        printf("Constants b.T:\n");
        for(i=0; i<N_cons; ++i){
           printf("\t %.2f", b[i][0]);
    	  }
        printf("\n");
    }
    
}

void swap_two_rows(double M[][N_col], int row_i, int row_j){
    if(row_i == row_j) {return;} // no need swapping
    int k = 0;
    double temp[N_col];
    for(k = 0; k < N_col; ++k){ 
        temp[k] = M[row_i][k]; // store the (row_i+1) row into temp
        M[row_i][k] = M[row_j][k]; // move the (row_j+1) row to the (row_i+1) row
        M[row_j][k] = temp[k]; // move the original (row_i+1) row to the (row_j+1) row
    }
}


// M should be reduced row-echolen form
void show_solution(double M[N_row][N_col], int n_var, int rank_A, int rank_M){ // N_var connot be a parameter
    if(debug>=1){
        printf("\n");
        if(rank_M!=rank_A){
            printf("The linear system is infeasible, ");
            printf("since rank(A)=%d != rank(M)=%d.\n", rank_A, rank_M);
        }else if(rank_M==rank_A && rank_M < n_var){
            printf("The linear system has infinitely many solutions, ");
            printf("because rank(A)=%d == rank(M)=%d < N_var=%d.\n", rank_A, rank_M, n_var);
        }else if(rank_M==rank_A && rank_M == n_var){
            printf("The linear system has a unique solution: ");
            int i;
            for(i=0; i<n_var; ++i){
                printf("x%d = %.3f, ", i+1, M[i][N_col-1]);
            }
            printf("\n");
        }else{
            printf("Oops, sth. went wrong!\n");
        }
    }
}


