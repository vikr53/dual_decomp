/* This example formulates and solves the following simple QP model:

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

    start = omp_get_wtime();
    /* Create environment */

    error = GRBloadenv(&env, "qp.log");
    if (error) goto QUIT;

    /* Create an empty model */

    error = GRBnewmodel(env, &model, "qp", 0, NULL, NULL, NULL, NULL, NULL);
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

    /* First constraint: x + 2 y + 3 z <= 4 */

    ind[0] = 0; ind[1] = 1; ind[2] = 2;
    val[0] = 1; val[1] = 2; val[2] = 3;

    error = GRBaddconstr(model, 3, ind, val, GRB_EQUAL, 4.0, "c0");
    if (error) goto QUIT;

    /* Second constraint: x + y >= 1 */

    ind[0] = 0; ind[1] = 1;
    val[0] = 1; val[1] = 1;

    error = GRBaddconstr(model, 2, ind, val, GRB_EQUAL, 2.0, "c1");
    if (error) goto QUIT;

    /* Optimize model */

    error = GRBoptimize(model);
    if (error) goto QUIT;

    /* Write model to 'qp.lp' */

    error = GRBwrite(model, "lsq.lp");
    if (error) goto QUIT;

    /* Capture solution information */

    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;

    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;

    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 3, sol);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");
    if (optimstatus == GRB_OPTIMAL) {
        printf("Optimal objective: %.4e\n", objval);
        end = omp_get_wtime();
        printf("Time to complete: %.10f\n\n", end - start);
        printf("  x=%.4f, y=%.4f, z=%.4f\n", sol[0], sol[1], sol[2]);
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

    return 0;
    

}