
#define N_cons 60
#define N_var 60

// size of the augmented matrix
extern const int N_row, N_col;

// choose test cases, 0: randomly assign values
extern int option;

// whether to print out results, -2: Nothing, -1: time, 0: + size, 1: + value, 2: + input + gaussian matrix
extern int debug;

void print_M_seg(double M_seg[][N_col], int N_row);

void print_M(double M[N_row][N_col]);

void form_M(double M[N_row][N_col], double A[N_cons][N_var], double b[N_var][1]);

void load_linear_system(int option, double A[N_cons][N_var], double b[N_var][1]);

void swap_two_rows(double M[][N_col], int row_i, int row_j);

void show_solution(double M[N_row][N_col], int n_var, int rank_A, int rank_M);



void Gaussian_OMP(double M[N_row][N_col], int* rank_A, int* rank_M);