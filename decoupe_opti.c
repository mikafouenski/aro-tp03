#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#define BUFFSIZE 1000
#define ZERO 0.0
#define COEFS 1000
#define NBCONSTRAINTS 4

typedef struct matrix {
    int constraints[COEFS + 1];
    int variables[COEFS + 1];
    double coefs[COEFS + 1];
    int next;
} matrix;

matrix* newMatrix() {
    matrix* m = (matrix*) malloc(sizeof(matrix));
    m->next = 1;
    return m;
}

void addCoef(matrix* m, int constraint, int variable, double coef) {
    m->constraints[m->next] = constraint;
    m->variables[m->next] = variable;
    m->coefs[m->next] = coef;
    m->next = m->next + 1;
}

void initMatrix(matrix* m) {
    double init[5][4] = {{2, 0, 0, 0},
                         {1, 1, 0, 1},
                         {1, 1, 0, 0},
                         {1, 0, 1, 1},
                         {1, 0, 1, 0}};

    for (int j = 1; j <= 5; ++j)
        for (int i = 1; i <= NBCONSTRAINTS; ++i)
            addCoef(m, i, j, init[j][i]);
}

void loadMatrix(glp_prob* lp, matrix* m) {
    glp_load_matrix(lp, m->next - 1, m->constraints, m->variables, m->coefs);
}

void addConstraints(glp_prob* lp, double* coefs, int start, int end, int bound) {
    glp_add_rows(lp, end);
    for (int i = start; i <= end; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "aux_%d", i);
        glp_set_row_name(lp, i, name);
        glp_set_row_bnds(lp, i, bound, coefs[i - 1], ZERO);
    }
}

void addVarrables(glp_prob* lp, int start, int end) {
    glp_add_cols(lp, end);
    for (int i = start; i <= end; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "x_%d", i);
        glp_set_col_name(lp, i, name);
        glp_set_col_bnds(lp, i, GLP_LO, ZERO, ZERO);
        glp_set_obj_coef(lp, i, 1.0);
    }
}

int main(int argc, char const *argv[]) {
    double coefs[NBCONSTRAINTS] = {97, 610, 395, 211};

    glp_prob *lp = glp_create_prob();
    matrix* m = (matrix*) malloc (sizeof(matrix));
    glp_set_prob_name(lp, "decoupe");
    glp_set_obj_dir(lp, GLP_MIN);

    addConstraints(lp, coefs, 1, NBCONSTRAINTS, GLP_LO);
    addVarrables(lp, 1, 2);
    printf("Ajout OK.\n");

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.meth = GLP_DUAL;

    initMatrix(m);

    loadMatrix(lp, m);

    glp_simplex(lp, &parm);

    glp_write_lp(lp, NULL ,"Prog.txt");

    // double result[4];
    // for (int i = 0; i < 4; ++i){
    //     result[i] = glp_get_row_dual(lp,i);
    // } 


    return 0;
}
