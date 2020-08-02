#include <stdlib.h>
#include <stdio.h>
#include "gurobi_c.h"
#include "omp.h"

#define NUM_THREADS 3
// #define N_ITER 50
// #define ALPHA 0.1

int main (int argc, char *argv[]) 
{   
    FILE *fp = fopen("./output/decomp_lsq.csv", "wb+");
    int nthreads;
    double start, end;
    
    int N_ITER = 50;
    double ALPHA = 0.2;
    char *eptr;
    if (argc > 1) {
        // set alpha runtime
        ALPHA = strtod(argv[1], &eptr);
    }

    omp_set_num_threads(NUM_THREADS);

    double A[2][3] = {{1, 2, 3},{1, 1, 0}};
    double b[] = {4, 2};

    // initiall guess for dual variables
    double y[] = {1, 1};

    double x[] = {0,0,0};
    start = omp_get_wtime();
    for (int i = 0; i < N_ITER; i++) {

        #pragma omp parallel
        {
            int id, nthrds;

            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();

            
            if (id == 0) {
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

                error = GRBloadenv(&env, "decomp1.log");
                if (error) goto QUIT;

                /* Create an empty model */

                error = GRBnewmodel(env, &model, "decomp1", 0, NULL, NULL, NULL, NULL, NULL);
                if (error) goto QUIT;

                /* Add variables */
                char *names[] = {"x"};
                error = GRBaddvars(model, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                    names);
                if (error) goto QUIT;

                /* Quadratic objective terms */

                qrow[0] = 0; 
                qcol[0] = 0; 
                qval[0] = 1;

                error = GRBaddqpterms(model, 1, qrow, qcol, qval);
                if (error) goto QUIT;

                /* Linear objective term */

                error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 0, y[0]*A[0][0] + y[1]*A[1][0]);
                if (error) goto QUIT;

                /* Optimize model */

                error = GRBoptimize(model);
                if (error) goto QUIT;

                /* Write model to 'qp.lp' */

                error = GRBwrite(model, "decomp1.lp");
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

                    printf("  x=%.4f \n", sol[0]);
                    x[0] = sol[0];
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
            // Solve problem 2
            if (id == 1) {
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

                error = GRBloadenv(&env, "decomp2.log");
                if (error) goto QUIT;

                /* Create an empty model */

                error = GRBnewmodel(env, &model, "decomp2", 0, NULL, NULL, NULL, NULL, NULL);
                if (error) goto QUIT;

                /* Add variables */
                char *names[] = {"y"};
                error = GRBaddvars(model, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                    names);
                if (error) goto QUIT;

                /* Quadratic objective terms */

                qrow[0] = 0; 
                qcol[0] = 0; 
                qval[0] = 1;

                error = GRBaddqpterms(model, 1, qrow, qcol, qval);
                if (error) goto QUIT;

                /* Linear objective term */

                error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 0, y[0]*A[0][1] + y[1]*A[1][1]);
                if (error) goto QUIT;

                /* Optimize model */

                error = GRBoptimize(model);
                if (error) goto QUIT;

                /* Write model to 'qp.lp' */

                error = GRBwrite(model, "decomp2.lp");
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

                    printf("  y=%.4f \n", sol[0]);
                    x[1] = sol[0];
                } else if (optimstatus == GRB_INF_OR_UNBD) {
                    printf("Model is infeasible or unbounded\n");
                } else {
                    printf("Optimization was stopped early\n");
                }
                
            }
            // Solve problem 3
            if (id == 2) {
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

                error = GRBloadenv(&env, "decomp3.log");
                if (error) goto QUIT;

                /* Create an empty model */

                error = GRBnewmodel(env, &model, "decomp3", 0, NULL, NULL, NULL, NULL, NULL);
                if (error) goto QUIT;

                /* Add variables */
                char *names[] = {"z"};
                error = GRBaddvars(model, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                    names);
                if (error) goto QUIT;

                /* Quadratic objective terms */

                qrow[0] = 0; 
                qcol[0] = 0; 
                qval[0] = 1;

                error = GRBaddqpterms(model, 1, qrow, qcol, qval);
                if (error) goto QUIT;

                /* Linear objective term */

                error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 0, y[0]*A[0][2] + y[1]*A[1][2]);
                if (error) goto QUIT;

                /* Optimize model */

                error = GRBoptimize(model);
                if (error) goto QUIT;

                /* Write model to 'qp.lp' */

                error = GRBwrite(model, "decomp3.lp");
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

                    printf("  z=%.4f \n", sol[0]);
                    x[2] = sol[0];
                } else if (optimstatus == GRB_INF_OR_UNBD) {
                    printf("Model is infeasible or unbounded\n");
                } else {
                    printf("Optimization was stopped early\n");
                }
                
            }
            
        }
        printf("X RESULT: %f, %f, %f \n", x[0], x[1], x[2]);
        
        fprintf(fp, "%d, %f, %f, %f \n", i, x[0], x[1], x[2]);
        

        y[0] = y[0] + ALPHA * (A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2] - b[0]);
        y[1] = y[1] + ALPHA * (A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2] - b[1]);
    }

    end = omp_get_wtime();
    printf("Time to complete: %.4f\n\n", end - start);
    printf("Number of threads: %d\n", nthreads);
    fclose(fp);
}