#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>

#define BUFFSIZE 1000
#define COEFS 1000
#define NBCONSTRAINTS 4
#define NBINIT 5

typedef struct matrix {
    int ia[COEFS + 1];
    int ja[COEFS + 1];
    double ar[COEFS + 1];
    int next;
} matrix;

matrix* newMatrix() {
    matrix* m = (matrix*) malloc(sizeof(matrix));
    m->next = 1;
    return m;
}

void addCoef(matrix* m, int i, int j, double coef) {
    m->ia[m->next] = i;
    m->ja[m->next] = j;
    m->ar[m->next] = coef;
    m->next = m->next + 1;
}

void load(glp_prob* lp, matrix* m) {
    glp_load_matrix(lp, m->next - 1, m->ia, m->ja, m->ar);
}

void addConstraints(glp_prob* lp, double* orders, int start, int end) {
    glp_add_rows(lp, end);
    for (int i = start; i <= end; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "aux_%d", i);
        glp_set_row_name(lp, i, name);
        glp_set_row_bnds(lp, i, GLP_LO, orders[i - 1], 0);
    }
}

void addVarrables(glp_prob* lp, int start, int end) {
    glp_add_cols(lp, end);
    for (int i = start; i <= end; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "x_%d", i);
        glp_set_col_name(lp, i, name);
        glp_set_col_bnds(lp, i, GLP_LO, 0, 0);
        glp_set_obj_coef(lp, i, 1);
    }
}

void init(glp_prob* lp, matrix* m, double* orders) {
    addConstraints(lp, orders, 1, NBCONSTRAINTS);
    addVarrables(lp, 1, 5);
    // 5 première colonnes
    double init[5][NBCONSTRAINTS] = {
        {2, 0, 0, 0},
        {1, 1, 0, 1},
        {1, 1, 0, 0},
        {1, 0, 1, 1},
        {1, 0, 1, 0}
    };

    for (int j = 1; j <= 5; ++j)
        for (int i = 1; i <= NBCONSTRAINTS; ++i)
            addCoef(m, i, j, init[j - 1][i - 1]);
}

void simplexDual(glp_prob* lp, double* sol) {
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.meth = GLP_DUAL;
    parm.msg_lev = GLP_MSG_OFF;
    glp_simplex(lp, &parm);
    for (int i = 1; i <= NBCONSTRAINTS; ++i)
        sol[i - 1] = glp_get_row_dual(lp, i);
}

void knapsack(double* dual) {
    double sizes[NBCONSTRAINTS] = {45, 36, 31, 14};
    glp_prob *knapsack = glp_create_prob();
    matrix* m = newMatrix();
    glp_set_prob_name(knapsack, "Knapsack");
    glp_set_obj_dir(knapsack, GLP_MAX);

    glp_add_rows(knapsack, 1);
    glp_set_row_name(knapsack, 1, "size");
    glp_set_row_bnds(knapsack, 1, GLP_UP, 0, 100);

    glp_add_cols(knapsack, NBCONSTRAINTS);
    for (int i = 1; i <= NBCONSTRAINTS; ++i) {
        char name[BUFFSIZE];
        sprintf(name, "x_%d", i);
        glp_set_col_name(knapsack, i, name);
        glp_set_col_bnds(knapsack, i, GLP_LO, 0, 0);
        glp_set_obj_coef(knapsack, i, dual[i - 1]);
        glp_set_col_kind(knapsack, i, GLP_IV);
        addCoef(m, 1, i, sizes[i]);
    }

    load(knapsack, m);

    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    glp_intopt(knapsack, &parm);

    glp_print_sol(knapsack, "a");

    free(m);
    glp_delete_prob(knapsack);
}

int main(int argc, char const *argv[]) {
    double orders[NBCONSTRAINTS] = {97, 610, 395, 211};

    glp_prob *lp = glp_create_prob();
    matrix* m = newMatrix();
    glp_set_prob_name(lp, "Cut");
    glp_set_obj_dir(lp, GLP_MIN);
    printf("Creation OK.\n");

    init(lp, m, orders);
    printf("Init OK.\n");

    load(lp, m);
    printf("Load OK.\n");

    double sol[NBCONSTRAINTS];
    simplexDual(lp, sol);
    printf("Simplex OK.\n");

    knapsack(sol);

    // glp_print_sol(lp, "a");
    // glp_write_lp(lp, NULL, "a");





    // free
    free(m);
    glp_delete_prob(lp);
    glp_free_env();

    return 0;
}
