#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gurobi_c.h"
#include "omp.h"

// #define NUM_THREADS 3
// #define N_ITER 50
// #define ALPHA 0.1

int main (int argc, char *argv[]) 
{   
    FILE *input, *output; 
    output = fopen("./output/decomp_lsq.csv", "wb+");
    int nthreads;
    double start, end;
    
    int N_ITER = 100;
    double ALPHA = 0.2;
    int NUM_THREADS = 1;
    char *eptr;
    if (argc > 1) {
        // set alpha runtime
        ALPHA = strtod(argv[1], &eptr);
    }

    /* Load matrix A and initialize */
    input = fopen("./input/A.matrix", "r");
    size_t len = 0;
    char *line = NULL;
    ssize_t read = getline(&line, &len, input);

    /** Get size of A matrix **/
    int rows, cols;
    int i = 0;
    char *tok;
    for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n")) {
      if (i == 0) {
        rows = atoi(tok);
        printf("ROWS: %d\n", rows);
      } else if (i == 1) {
        cols = atoi(tok);
        printf("COLS: %d\n", cols);
      }
      i++;
    }
    
    double A[rows][cols];
    /* Init number of threads as size of x */
    NUM_THREADS = cols;

    /** Load values from A.matrix to A **/
   
    int row = 0;
    while ( (read = getline(&line, &len, input)) != -1 ) {
      i = 0;
      for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n")) {
        char *eptr;
        A[row][i] = strtod(tok, &eptr);
        i++;
      }
      row++;
    }
	
    fclose(input);

    /* Load matrix B and initialize */
    input = fopen("./input/B.matrix", "r");

    double b[cols];
    
    /** Load values from b.matrix to b **/
    row = 0;
    while ( (read = getline(&line, &len, input)) != -1 ) {
      for (tok = strtok(line, ";"); tok && *tok; tok = strtok(NULL, ";\n")) {
        char *eptr;
        b[row] = strtod(tok, &eptr);
      }
      row++;
    }
	
    fclose(input);

    omp_set_num_threads(NUM_THREADS);

    // intialize x (solution) and the dual variables (as initial guesses)
    double x[cols];
    double y[rows];
    for (int i = 0; i < cols; i++) {
      x[i] = 0;
    }
    for (int i = 0; i < rows; i++) {
      y[i] = 0;
    }


    start = omp_get_wtime();
    for (int i = 0; i < N_ITER; i++) {

        #pragma omp parallel
        {
            int id, nthrds;

            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();

            /* Based on thread ID */
            nthreads = nthrds;
            GRBenv   *env   = NULL;
            GRBmodel *model = NULL;
            int       error = 0;
            double    sol[1];
            int       qrow[1];
            int       qcol[1];
            double    qval[1];
            int       optimstatus;
            double    objval;

            /* Create environment */
            char env_log_name[15];
            sprintf(env_log_name, "decomp%d.log", id);
            error = GRBloadenv(&env, env_log_name);
            if (error) goto QUIT;

            /* Create an empty model */

            error = GRBnewmodel(env, &model, "decomp", 0, NULL, NULL, NULL, NULL, NULL);
            if (error) goto QUIT;

            /* Add variable */
            error = GRBaddvars(model, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                NULL);
            if (error) goto QUIT;

            /* Quadratic objective terms */

            qrow[0] = 0; 
            qcol[0] = 0; 
            qval[0] = 1;

            error = GRBaddqpterms(model, 1, qrow, qcol, qval);
            if (error) goto QUIT;

            /* Linear objective term */
            /* Calculating coeff */
            double coeff = 0;
            for (int j = 0; j < rows; j++) {
                coeff = coeff + y[j]*A[j][id];
            }
            error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 0, coeff);
            if (error) goto QUIT;

            /* Optimize model */

            error = GRBoptimize(model);
            if (error) goto QUIT;

            /* Write model to 'qp.lp' */
            char env_model_name[15];
            sprintf(env_model_name, "decomp%d.lp", id);
            error = GRBwrite(model, env_model_name);
            if (error) goto QUIT;

            /* Capture solution information */

            error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
            if (error) goto QUIT;

            error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
            if (error) goto QUIT;

            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 1, sol);
            if (error) goto QUIT;

            printf("\nOptimization complete\n");
            if (optimstatus == GRB_OPTIMAL) {
                printf("Optimal objective: %.4e\n", objval);

                printf("  x%d = %.4f \n", id, sol[0]);
                x[id] = sol[0];
            } else if (optimstatus == GRB_INF_OR_UNBD) {
                printf("Model is infeasible or unbounded\n");
            } else {
                printf("Optimization was stopped early\n");
            }
            
            QUIT:

                /* Error reporting */

                if (error) {
                    printf("ERROR: %s\n", GRBgeterrormsg(env));
                    exit(1);
                }

                /* Free model */

                GRBfreemodel(model);

                /* Free environment */

                GRBfreeenv(env);
            
        }
        printf("  x1=%.4f, x2=%.4f, x3=%.4f, x4=%.4f, x5=%.4f, x6=%.4f, x7=%.4f, x8=%.4f, x9=%.4f, x10=%.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9] );

        
        fprintf(output, "%d, %f, %f, %f \n", i, x[0], x[1], x[2]);
        
        /* Update y */
        double ax;
        for (int j = 0; j < rows; j++) {
            ax = 0;
            for (int k = 0; k < cols; k++) {
                ax = ax + x[k]*A[j][k];
            }
            y[j] = y[j] + ALPHA * (ax-b[j]);
        }

    }

    end = omp_get_wtime();
    printf("Time to complete: %.4f\n\n", end - start);
    printf("Number of threads: %d\n", nthreads);
    fclose(output);
}