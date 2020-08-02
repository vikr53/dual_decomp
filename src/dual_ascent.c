/* This example formulates and solves the following simple QP model using Dual Acent:

     minimize    x^2 +  y^2 + z^2
     subject to  x + 2 y + 3 z = 4
                 x +   y       = 2
                 x, y, z non-negative
*/

#include <stdlib.h>
#include <stdio.h>
#include "gurobi_c.h"
#include "omp.h"

int main(int argc, char *argv[]) {
    FILE *fp;
    fp=fopen("./output/ascent_lsq.csv", "wb+");
    double start, end;
    double N_ITER = 50;
    double ALPHA = 0.2;
    double A[2][3] = {{1, 2, 3},{1, 1, 0}};
    double b[] = {4, 2};

    // initiall guess for dual variables
    double y[] = {1, 1};

    double x[] = {0,0,0};
    start = omp_get_wtime();
    for (int i = 0; i < N_ITER; i++) {
        GRBenv   *env   = NULL;
        GRBmodel *model = NULL;
        int       error = 0;
        double    sol[3];
        int       ind[3];
        double    val[3];
        int       qrow[5];
        int       qcol[5];
        double    qval[5];
        char      vtype[3];
        int       optimstatus;
        double    objval;
        double    start, end;

        /* Create environment */

        error = GRBloadenv(&env, "ascent.log");
        if (error) goto QUIT;

        /* Create an empty model */

        error = GRBnewmodel(env, &model, "ascent", 0, NULL, NULL, NULL, NULL, NULL);
        if (error) goto QUIT;

        /* Add variables */
        char *names[] = {"x", "y", "z"};
        error = GRBaddvars(model, 3, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                            names);
        if (error) goto QUIT;

        /* Quadratic objective terms */

        qrow[0] = 0; qrow[1] = 1; qrow[2] = 1; qrow[3] = 1; qrow[4] = 2;
        qcol[0] = 0; qcol[1] = 1; qcol[2] = 1; qcol[3] = 2; qcol[4] = 2;
        qval[0] = 1; qval[1] = 1; qval[2] = 0; qval[3] = 0; qval[4] = 1;

        error = GRBaddqpterms(model, 5, qrow, qcol, qval);
        if (error) goto QUIT;

        /* Linear objective term */
        /* x */
        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 0, y[0]*A[0][0] + y[1]*A[1][0]);
        if (error) goto QUIT;
        /* y */
        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 1, y[0]*A[0][1] + y[1]*A[1][1]);
        if (error) goto QUIT;
        /* z */
        error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, 2, y[0]*A[0][2] + y[1]*A[1][2]);
        if (error) goto QUIT;

        /* Optimize model */

        error = GRBoptimize(model);
        if (error) goto QUIT;

        /* Write model to 'qp.lp' */

        error = GRBwrite(model, "dual_ascent.lp");
        if (error) goto QUIT;

        /* Capture solution information */

        error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
        if (error) goto QUIT;

        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
        if (error) goto QUIT;

        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 3, sol);
        if (error) goto QUIT;

        if (optimstatus == GRB_OPTIMAL) {
            printf("Optimal objective: %.4e\n", objval);
            printf(" x=%.4f, y=%.4f, z=%.4f\n", sol[0], sol[1], sol[2]);
            x[0] = sol[0];
            x[1] = sol[1];
            x[2] = sol[2];
        } else if (optimstatus == GRB_INF_OR_UNBD) {
            printf("Model is infeasible or unbounded\n");
        } else {
            printf("Optimization was stopped early\n");
        }
        
        fprintf(fp, "%d, %f, %f, %f \n", i, x[0], x[1], x[2]);

        y[0] = y[0] + ALPHA * (A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2] - b[0]);
        y[1] = y[1] + ALPHA * (A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2] - b[1]);

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
    
    end = omp_get_wtime();
    printf("Time to complete: %.4f\n\n", end - start);
}