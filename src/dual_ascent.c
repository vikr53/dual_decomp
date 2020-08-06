/* This example formulates and solves the following simple QP model using Dual Acent:

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
    FILE *input, *output;
    output = fopen("./output/ascent_lsq.csv", "wb+");
    double start, end;
    double N_ITER = 100;
    double ALPHA = 0.2;

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
    
    //printf("A: %f %f %f %f", A[1][0], A[2][1], A[3][2], A[4][3]);
    //printf("b: %f %f %f %f", b[0], b[1], b[2], b[9]);
     

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
        GRBenv   *env   = NULL;
        GRBmodel *model = NULL;
        int       error = 0;
        double    sol[cols];
        int       ind[cols];
        double    val[cols];
        int       qrow[cols];
        int       qcol[cols];
        double    qval[cols];
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
        error = GRBaddvars(model, cols, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                            NULL);
        if (error) goto QUIT;

        /* Quadratic objective terms */

        for (int j = 0; j < (sizeof(x)/sizeof(x[0])); j++) {
            qrow[j] = j; 
            qcol[j] = j; 
            qval[j] = 1; 
        }

        error = GRBaddqpterms(model, cols, qrow, qcol, qval);
        if (error) goto QUIT;

        /* Linear objective term */
        /* Calculating coefficients - y^TA */
        double coeff;
        for (int j = 0; j < cols; j++) {
            coeff = 0;
            for (int k = 0; k < rows; k++) {
                coeff = coeff + y[k]*A[k][j];
            }
            error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, j, coeff);
            if (error) goto QUIT;
        }

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

        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cols, sol);
        if (error) goto QUIT;

        if (optimstatus == GRB_OPTIMAL) {
            printf("Optimal objective: %.4e\n", objval);
            for (int j = 0; j < cols; j++) {
                x[j] = sol[j];
            }
        } else if (optimstatus == GRB_INF_OR_UNBD) {
            printf("Model is infeasible or unbounded\n");
        } else {
            printf("Optimization was stopped early\n");
        }
        
        printf("  x1=%.4f, x2=%.4f, x3=%.4f, x4=%.4f, x5=%.4f, x6=%.4f, x7=%.4f, x8=%.4f, x9=%.4f, x10=%.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9] );
        fprintf(output, "%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", i, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        /* Update y */
        double ax;
        for (int j = 0; j < rows; j++) {
            ax = 0;
            for (int k = 0; k < cols; k++) {
                ax = ax + x[k]*A[j][k];
            }
            y[j] = y[j] + ALPHA * (ax-b[j]);
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

    fclose(output);
    end = omp_get_wtime();
    printf("Time to complete: %.4f\n\n", end - start);
}
