#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */

void sacADos(int y1, int y2, int y3, int y4,double* resultat){
    int i,j,k;
    double matrix_coef[4] = {45,36,31,14};
    double matrix_resultat[4] = {y1,y2,y3,y4};
    char name[10] ;
    glp_prob *lp;
    int ia[1+1000], ja[1+1000];
    double ar[1+1000];

    lp = glp_create_prob();
    glp_set_prob_name(lp, "Sac a dos");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, 1);
    glp_set_row_name(lp, 1, "SOMME");
    glp_set_row_bnds(lp, 1, GLP_UP, 0, 100);
    
    glp_add_cols(lp, 4);
    for(i=0;i<4;i++){
        snprintf(name,10, "P%d", i+1);
        glp_set_col_name(lp, i+1, name);
        glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, i+1, matrix_resultat[i]);
        glp_set_col_kind(lp,i+1,GLP_CV);
    }
    
    k=1;
    for(j=1;j<=4;j++){
            ia[k] = 1; ja[k] = j; ar[k] = matrix_coef[k-1];
            k++;
    }
    glp_load_matrix(lp, 4, ia, ja, ar);
    glp_simplex(lp, NULL);
    resultat[0] =  glp_get_obj_val(lp);  
    for(i =i ;i<5;i++){
        resultat[i]=glp_get_col_prim(lp, i-1);
    }                 
    glp_delete_prob(lp);
}

void newMatCoef( double* matrix_coef, double* matrice_sac, int taille){
    int i;
    double newMat[taille+4];
    for(i=0;i<taille;i++){
        newMat[i]=matrix_coef[i];
    }
    for(i=0;i<4;i++){
        newMat[taille+i]=matrice_sac[i];
    }
    matrix_coef=newMat;
}

double* resolutionGLPK(double* matrix_coef,double* matrix_commande, int nbContrainte, int nbVariable, char* nom, int MIN_MAX){
    int i,j,k = 1;
    int tailleMat = 20;
    double y1, y2, y3, y4;
    double ResultatSacADos[5];
    glp_prob *lp;
    char name[10];
    double ar[1+1000];
    int ia[1+1000], ja[1+1000];
   
    for(;;){
        lp = glp_create_prob();
        glp_set_obj_dir(lp, MIN_MAX);
        glp_set_prob_name(lp, nom);
        glp_add_rows(lp, nbContrainte);

        for(i=0;i<nbContrainte;i++){
            snprintf(name,10, "R%d", i+1);
            glp_set_row_name(lp, i+1,name);
            glp_set_row_bnds(lp, i+1, GLP_LO, matrix_commande[i], 0.0);
        }
        glp_add_cols(lp, nbVariable);
        
        for(i=0;i<nbVariable;i++){
            snprintf(name,10, "P%d", i+1);
            glp_set_col_name(lp, i+1, name);
            glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(lp, i+1, 1.0);
        }

        for(j=1;j<=nbVariable;j++){
            for(i=1;i<=nbContrainte;i++){
                ia[k] = i; ja[k] = j; ar[k] = matrix_coef[k-1];
                k++;
            }
        }

        glp_load_matrix(lp, nbContrainte*nbVariable, ia, ja, ar);
        glp_simplex(lp, NULL);
        glp_write_lp(lp,NULL,"Prog.txt");
        

        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.meth = GLP_DUAL;
        glp_simplex(lp, &parm);

        y1 = glp_get_row_dual(lp,1);
        y2 = glp_get_row_dual(lp,2);
        y3 = glp_get_row_dual(lp,3);
        y4 = glp_get_row_dual(lp,4);
        sacADos(y1,y2,y3,y4,ResultatSacADos);
        glp_delete_prob(lp);
        //if (ResultatSacADos[0]==0){break;}

        //newMatCoef(matrix_coef,ResultatSacADos,tailleMat);
        tailleMat += 4;
    }

}

void resolution(double* matrix_coef,double* matrix_commande, int nbContrainte, int nbVariable){
    
   //double z, x1, x2, x3, x4, x5, y1, y2, y3, y4, resultatSacADos;
    resolutionGLPK(matrix_coef,matrix_commande,nbContrainte,nbVariable,"Decoupe",GLP_MIN);


    //z = glp_get_obj_val(lp);
    /*x1 = glp_get_col_prim(lp, 1);
    x2 = glp_get_col_prim(lp, 2);
    x3 = glp_get_col_prim(lp, 3);
    x4 = glp_get_col_prim(lp, 4);
    x5 = glp_get_col_prim(lp, 5);
    */

    //printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g; x4= %g; x5= %g; y1 = %g, y2 = %g, y3 = %g, y4 = %g, resultatSacADos = %g \n",z, x1, x2, x3, x4, x5,y1,y2,y3,y4,resultatSacADos);

}
 

int main(){
    double matrix_coef[20] = {2,0,0,0,1,1,0,1,1,1,0,0,1,0,1,1,1,0,1,0}; //les cinqs premieres colonnes
    double matrix_commande[4] = {97.0,610.0,395.0,211.0};
    resolution(matrix_coef,matrix_commande,4,5);
}
