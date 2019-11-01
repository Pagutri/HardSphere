/* Structure factor for hard sphere */

#include "libraries.h"

int main(int argc, char **argv)
{
    double eta = 0.25; /* Volume fraction */
    double k = 0.001; 
    double delta_k = 0.01;
    double ck = 0.0; /* Structure */
    double sk = 0.0; /* Structure factor */
    int n_ks = 4100;
    int i = 0;

    FILE *fp;
    fp = fopen("data/my_data.dat", "w");

    for(i = 0; i < n_ks; i++)
    {
        ck = calculate_ck(eta, k);
        sk = calculate_sk(eta, ck);

        fprintf(fp, "%f\t%f\n", k, sk);

        k += delta_k;
    }

    fclose(fp);
    
return(1);
}
