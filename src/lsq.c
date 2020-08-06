/* This example formulates and solves the following simple QP model:

     minimize    x^2 +  y^2 + z^2
     subject to  x + 2 y + 3 z = 4
                 x +   y       = 2
                 x, y, z non-negative
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gurobi_c.h"
#include "omp.h"

int main(int argc, char *argv[]) {
    FILE *input;

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

    printf("A %f, %f, %f\n", A[0][0], A[0][1], A[0][2]);

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
    
    GRBenv   *env   = NULL;
    GRBmodel *model = NULL;
    int       error = 0;
    double    sol[cols];
    int       ind[cols];
    double    val[cols];
    int       qrow[cols];
    int       qcol[cols];
    double    qval[cols];
    char      vtype[cols];
    int       optimstatus;
    double    objval;
    double    start, end;
    double    x[cols];

    start = omp_get_wtime();
    /* Create environment */

    error = GRBloadenv(&env, "qp.log");
    if (error) goto QUIT;

    /* Create an empty model */

    error = GRBnewmodel(env, &model, "qp", 0, NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;

    /* Add variables */
    error = GRBaddvars(model, cols, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL);
    if (error) goto QUIT;

    /* Quadratic objective terms */
    for (int i = 0; i < cols; i++) {
	  qrow[i] = i; 
	  qcol[i] = i; 
	  qval[i] = 1; 
	}

    error = GRBaddqpterms(model, cols, qrow, qcol, qval);
    if (error) goto QUIT;

    /* First constraint: x + 2 y + 3 z <= 4 */
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            ind[j]=j;
            val[j]=A[i][j];
        }
        error = GRBaddconstr(model, cols, ind, val, GRB_EQUAL, b[i], NULL);
        if (error) goto QUIT;
    }
    
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

    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cols, sol);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");
    if (optimstatus == GRB_OPTIMAL) {
        printf("Optimal objective: %.4e\n", objval);
        end = omp_get_wtime();
        printf("Time to complete: %.10f\n\n", end - start);
        printf("  x1=%.4f, x2=%.4f, x3=%.4f, x4=%.4f, x5=%.4f, x6=%.4f, x7=%.4f, x8=%.4f, x9=%.4f, x10=%.4f, x11=%.4f\n", sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6], sol[7], sol[8], sol[9], sol[10] );
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