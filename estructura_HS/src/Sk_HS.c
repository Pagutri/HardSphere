/*Factor de estructura de esfera dura*/

#include <stdio.h>
#include <math.h>

double calcular_ck(double eta, double k);
double calcular_sk(double eta, double ck);

int main(int argc, char **argv)
{
	double eta = 0.25; /* Fraccion de volumen */
	double k = 0.001; 
	double delta_k = 0.01;
	double ck = 0.0; /* Estructura */
	double sk = 0.0; /* Factor de estructura */
	int n_ks = 4100;
	int i = 0;

	FILE *fp;
	fp = fopen("data/my_data.dat", "w");

	for(i = 0; i < n_ks; i++)
	{
		ck = calcular_ck(eta, k);
		/*printf("ck: %lf  ", ck);*/
		sk = calcular_sk(eta, ck);
		/*printf("sk: %lf\n", sk);*/

		fprintf(fp, "%f\t%f\n", k, sk);

		k += delta_k;
	}

	fclose(fp);
	
return(1);
}

double calcular_ck(double eta, double k)
{
	double alpha = 0.0;
	double delta = 0.0;
	double gamma = 0.0;
	double denom = 0.0;
	double ck = 0.0;
	double pi = 4.0 * atan(1.0);

	denom = pow(1.0 - eta, 4.0);

	alpha = -pow(1.0 + 2.0 * eta, 2.0) / denom;
	delta = 6.0 * eta * pow(1.0 + 0.5 * eta, 2.0) / denom;
	gamma = -0.5 * eta * alpha;
	/*printf("alpha: %lf, delta: %lf, gamma: %lf\n", alpha, delta, gamma);

	ck = 4.0 * pi * (alpha * pow(k, 3.0) * (sin(k) - k * cos(k)) + 24.0 * gamma + \
	                 delta * pow(k, 2.0) * (2.0 * k * sin(k) - (pow(k, 2.0) - 2.0) * cos(k) - 2.0) + \
		             gamma * ((pow(k, 3.0) - 24.0 * k) * sin(k) - \
			                  (pow(k, 4.0) - 12.0 * pow(k, 2.0) + 24.0) * cos(k))) / pow(k, 6.0);*/

	ck = 4.0 * pi * (alpha * pow(k, 3.0) * (sin(k) - k * cos(k)) + 24.0 * gamma + \
	                 delta * pow(k, 2.0) * (2.0 * k * sin(k) - (pow(k, 2.0) - 2.0) * cos(k) - 2.0) + \
		             gamma * ((pow(k, 3.0) - 24.0 * k) * sin(k) - \
			                  (pow(k, 4.0) - 12.0 * pow(k, 2.0) + 24.0) * cos(k))) / pow(k, 6.0);

	return(ck);
}

double calcular_sk(double eta, double ck)
{
	double rho = 0.0;
	double sk = 0.0;
	double pi = 4.0 * atan(1.0);

	rho = 6.0 * eta / pi;
	sk = 1.0 / (1.0 - rho * ck);

	return(sk);
}