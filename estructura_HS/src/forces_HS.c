/* Forces from structure factor */

#include "libraries.h"

int main(int argc, char **argv)
{
    /* I will calculate just the first 10 time lapses */
    int n_t = 10;
    double t = 0.0; 
    double delta_t = 0.000001;
    double k = 7.1;
    int i = 0;
    double eta = 0.25;
    double ck = 0.0;
    double sk = 0.0;
    double f = 0.0;
    /*double fself = 0.0;*/
    FILE *fp1;
    /*FILE *fp2;*/

    ck = calculate_ck(eta, k);
    sk = calculate_sk(eta, ck);

    fp1 = fopen("data/force.dat", "w");
    /*fp2 = fopen("data/force_self.dat", "w");*/

    for(i = 0; i < n_t; i++)
    {
        t += delta_t;

    	f = exp((-1.0) * k * k * t / sk);
    	/*fself = exp(- k * k * t);*/

    	fprintf(fp1, "%f\t%f\n", t, f);
    	/* fprintf(fp2, "%f\t%f\n", t, fself); */
    }

    fclose(fp1);
    /* fclose(fp2); */

return(1);
}