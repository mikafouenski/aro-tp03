#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */

int main(){
	int i,j,k = 1;
	double matrix_coef[20] = {2.0,1.0,1.0,1.0,1.0, 0.0,1.0,1.0,0.0,0.0, 0.0,0.0,0.0,1.0,1.0, 0.0,1.0,0.0,1.0,0.0};

    glp_prob *lp;

	int ia[1+1000], ja[1+1000];

	double ar[1+1000], z, x1, x2, x3, x4, x5, y1, y2, y3, y4;

	lp = glp_create_prob();
	glp_set_prob_name(lp, "découpe");
	glp_set_obj_dir(lp, GLP_MIN);
	glp_add_rows(lp, 4);
	glp_set_row_name(lp, 1, "R1");
	glp_set_row_bnds(lp, 1, GLP_LO, 97.0, 0.0);
	glp_set_row_name(lp, 2, "R2");
	glp_set_row_bnds(lp, 2, GLP_LO, 610.0, 0.0);
	glp_set_row_name(lp, 3, "R3");
	glp_set_row_bnds(lp, 3, GLP_LO, 395.0, 0.0);
	glp_set_row_name(lp, 4, "R4");
	glp_set_row_bnds(lp, 4, GLP_LO, 211.0, 0.0);

	glp_add_cols(lp, 5);
	glp_set_col_name(lp, 1, "P1");
	glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 1, 1.0);
	glp_set_col_name(lp, 2, "P2");
	glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 2, 1.0);
	glp_set_col_name(lp, 3, "P3");
	glp_set_col_bnds(lp, 3, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 3, 1.0);
	glp_set_col_name(lp, 4, "P4");
	glp_set_col_bnds(lp, 4, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 4, 1.0);
	glp_set_col_name(lp, 5, "P5");
	glp_set_col_bnds(lp, 5, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 5, 1.0);

	for(i=1;i<5;i++){
		for(j=1;j<6;j++){
			ia[k] = i; ja[k] = j; ar[k] = matrix_coef[k-1];
			k++;
		}
	}

	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.meth = GLP_DUAL;

	glp_load_matrix(lp, 20, ia, ja, ar);
	glp_simplex(lp, NULL);
	glp_simplex(lp, &parm);
	z = glp_get_obj_val(lp);
	x1 = glp_get_col_prim(lp, 1);
	x2 = glp_get_col_prim(lp, 2);
	x3 = glp_get_col_prim(lp, 3);
	x4 = glp_get_col_prim(lp, 4);
	x5 = glp_get_col_prim(lp, 5);

	y1 = glp_get_row_dual(lp,1);
	y2 = glp_get_row_dual(lp,2);
	y3 = glp_get_row_dual(lp,3);
	y4 = glp_get_row_dual(lp,4);

	/*SAC A DOS*/

	double matrix_coef2[4] = {45,36,31,14};

	glp_prob *lp2;
	int ia2[1+1000], ja2[1+1000];
	double ar2[1+1000];

	lp2 = glp_create_prob();
	glp_set_prob_name(lp2, "Sac a dos");
	glp_set_obj_dir(lp2, GLP_MAX);
	glp_add_rows(lp2, 1);
	glp_set_row_name(lp2, 1, "SOMME");
	glp_set_row_bnds(lp2, 1, GLP_UP, 0, 100);

	glp_add_cols(lp2, 4);
	glp_set_col_name(lp2, 1, "P1");
	glp_set_col_bnds(lp2, 1, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp2, 1, y1);
	glp_set_col_kind(lp2,1,GLP_CV);

	glp_set_col_name(lp2, 2, "P2");
	glp_set_col_bnds(lp2, 2, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp2, 2, y2);
	glp_set_col_kind(lp2,2,GLP_CV);

	glp_set_col_name(lp, 3, "P3");
	glp_set_col_bnds(lp2, 3, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp2, 3, y3);
	glp_set_col_kind(lp2,3,GLP_CV);

	glp_set_col_name(lp2, 4, "P4");
	glp_set_col_bnds(lp2, 4, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp2, 4, y4);
	glp_set_col_kind(lp2,4,GLP_CV);

	k=1;
	for(j=1;j<5;j++){
			ia2[k] = 1; ja2[k] = j; ar2[k] = matrix_coef2[k-1];
			k++;
	}
	glp_load_matrix(lp2, 4, ia2, ja2, ar2);
	glp_simplex(lp2, NULL);

	/*FIN SAC A DOS */

	glp_write_lp(lp,NULL,"Prog.txt");
	glp_write_lp(lp2,NULL,"SacADos.txt");
	printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g; x4= %g; x5= %g; y1 = %g, y2 = %g, y3 = %g, y4 = %g \n",z, x1, x2, x3, x4, x5,y1,y2,y3,y4);
	glp_delete_prob(lp);
	glp_delete_prob(lp2);
}
