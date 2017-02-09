#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#define BUFFSIZE 1000
#define ZERO 0.0
#define COEFS 1000

typedef struct matrix {
    int constraints[COEFS + 1];
    int variables[COEFS + 1];
    double coefs[COEFS + 1];
    int size;
} matrix;

matrix* newMatrix() {
    matrix* m = (matrix*) malloc(sizeof(matrix));
    m->size = 1;
    return m;
}

void addCoef(matrix* m, int constraint, int variable, double coef) {
    m->constraints[m->size] = constraint;
    m->variables[m->size] = variable;
    m->coefs[m->size] = coef;
    m->size = m->size + 1;
}

void loadMatrix(glp_prob* lp, matrix* m) {
    glp_load_matrix(lp, m->size - 1, m->constraints, m->variables, m->coefs);
}

void addAuxVariables(glp_prob* lp, double* coefs, int start, int end, int bound) {
    glp_add_rows(lp, (end - start));
    for (int i = start; i <= end; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "aux_%d", i);
        glp_set_row_name(lp, i, name);
        glp_set_row_bnds(lp, i, bound, coefs[i], ZERO);
    }
}

void addPositivityConstraints(glp_prob* lp, int start, int end) {
    glp_add_cols(lp, (end - start));
    for (int i = start; i <= end; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "positivity_%d", i);
        glp_set_col_name(lp, i, name);
        glp_set_col_bnds(lp, i, GLP_LO, ZERO, ZERO);
        glp_set_obj_coef(lp, i, 1.0);
    }
}

void generateColumn(int i, int j) {
    /* TODO */
}

int main(int argc, char const *argv[]) {

    return 0;
}
