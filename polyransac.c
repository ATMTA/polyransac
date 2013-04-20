#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <time.h>

const double PI = 3.14159265359;

typedef struct Params{
	int polyDegree;	//max degree of out polynom
	int maxTrial;	//number of RANSAC iterations
	int maxMPts; 	//Max. #pts to fit polynomial
	int size;		//length of total sample
	double tolerance; // Pt. to model thr. for inliers
} Params;

typedef struct Result{
	double inliners;	//relative number of inliners
	double* coef;	//coefficients of model
} Result;

double d_pow(double, int);
double d_max(double, double);
double d_min(double, double);
int count_inliners(double*, double*, double*, Params*);
void solveRLS(double*, double*, Params*, Result*);
void solveLS(double*, double*, int, int, double*, double*);
void printInFile(double*, int, const char*);

void solveRLS(double* x, double* y, Params* par, Result* result){	
	int max_pts = par->maxMPts;
	int size = par->size;	
	int degree = par->polyDegree;
	int max_trial = par->maxTrial;
	
	result->inliners = -1.0;	
	
	double* V = (double*)malloc((degree+1)*max_pts*sizeof(double));	
	double* b = (double*)malloc(max_pts*sizeof(double));

	double* hs_v = (double*)malloc(max_pts*sizeof(double));
	double* hs_w = (double*)malloc((degree+1)*sizeof(double));
	
	int current_index, i, k;
	double current_x, temp, pNoOutliners;	
	const double p = 0.99;
	
	//begin ransac
	for(i=0; i<max_pts; i++)
		V[i] = 1.0;
	
	int trialCount = 0; int N = 1;
	while(N>trialCount){
	
		for(i=0; i<max_pts;){
			current_index = rand()%size;
			current_x = x[current_index];
			k=0;
			while(k<i && V[max_pts+k] != current_x)
				k++;
			if(k!=i)
				continue;
			temp = current_x;
			for(k=1; k<=degree; k++, temp *= current_x)
				V[i+k*max_pts] = temp;				
			b[i] = y[current_index];
			i++;
		}
		
		solveLS(V, b, max_pts, degree+1, hs_v, hs_w);
		
		k = count_inliners(x,y,b, par);	
		if(1.0*k/size > result->inliners){
			result->inliners = 1.0*k/size;
			for(i=0; i<=degree; i++)
				result->coef[i] = b[i];			
			pNoOutliners = 1.0 - d_pow(1.0*k/size, max_pts);
			pNoOutliners = d_max(1e-9, pNoOutliners);
			pNoOutliners = d_min(1-1e-9, pNoOutliners);
			N = log(1-p)/log(pNoOutliners);			
		}

		trialCount++;
		if(trialCount>max_trial)
			break;
	}
}

int count_inliners(double* x, double* y, double* coef, Params* par){
	int res = 0;
	int size = par->size;
	int degree = par->polyDegree;	
	double tolerance = par->tolerance;
	int k;
	double poly_val, current_x;
	for(int i=0; i<size; i++){
		current_x = x[i];			
		poly_val = 0.0;
		//вычисление значения полинома
		for(k=degree; k>=0; k--)
			poly_val = poly_val*current_x+coef[k];
		//вычисление отклонения
		if((poly_val-y[i])*(poly_val-y[i]) < tolerance)
			res++;		
	}
	return res;
	//return 15;
}

double d_max(double x, double y){
	return x>y ? x:y;
}

double d_min(double x, double y){
	return x<y ? x:y;
}

double d_pow(double x, int y){
	double res = 1.0;
	int i = 0;
	while(i<y){
		res *= x;
		i++;
	}
	return res;
}

// решение полноранговой (m>n) задачи наименьших квадратов
// через QR-разложение (алгоритм Хаусхолдера)
// матрица хранится в column-major order
// на выходе в матрице в правом верхнем углу записана R
// в нижней части хранятся использованные векторы 
// Хаусхолдера с нормированной 1-й компонентой (v=1)
// решение хранится в pb(1:n)
void solveLS(double* pA, double* pb, int m, int n, double* pv, double* pw){
	double mu, sign;	
	for(int j=0; j<n; j++)
	{
		//matrix by column-major a(i,j) = A->pmatrix(i+j*rows)
		//вычисление очередного вектора Хаусхолдера
		mu = 0.0;
		for(int l=j; l<m; l++)
			mu += pA[j*m+l]*pA[j*m+l];

		if(mu > 1e-10){
			mu = sqrt(mu);		
			if(abs(pA[j+j*m])<1e-10)
				sign = 0;
			else
				sign = pA[j+j*m]>0 ? 1 : -1;
			mu = pA[j+j*m] + sign*mu;
		} else {
			mu = 1.0;
		}
		pv[j] = 1.0;
		for(int l=j+1; l<m; l++)
			pv[l] = pA[j*m+l]/mu;
		
		// модификация матрицы
		mu = 0.0;
		for(int l=j; l<m; l++)
			mu += pv[l]*pv[l];
		mu = -2/mu;			
		for(int l=j; l<n; l++)
		{
			pw[l] = 0.0;
			for(int k=j; k<m; k++)
				pw[l] += pA[k+l*m]*pv[k];
			pw[l] *= mu;
		}				

		for(int k=j; k<n; k++)
			for(int i=j; i<m; i++)
				pA[i+k*m] += pv[i]*pw[k];

		// запись нетривиальной части вектора Хаусхолдера
		for(int k= j+1; k<m; k++)
			pA[k+j*m] = pv[k];		

		// модификация правой части
		for(int l=j; l<m; l++)
		{
			pv[l] = pb[j];
			for(int k=j+1; k<m; k++)
				pv[l] += pb[k]*pA[k+j*m];
			pv[l] *= mu;
		}

		pb[j] += pv[j];
		for(int l=j+1; l<m; l++)
			pb[l] += pA[l+j*m]*pv[l];
	}
	
	// решение верхнетреугольной системы Rx=Q'b
	for(int j=n-1; j>0; j--){
		pb[j] /= pA[j+j*m];
		for(int k=j-1; k>-1; k--)
			pb[k] -= pA[k+j*m]*pb[j];
	}
	pb[0] /= pA[0];
}

int main(int argc, char** argv){
	
	// параметры функционирования алгоритма
	// (размер исходной выборки, размер подвыборки для построения модели,
	// максимальное число итераций, максимальная степень полинома,
	// порог для инлайнеров)

	// maxMPts >= polyDegree+1
	Params par;
	par.size = 5000;
	par.maxMPts = 17;
	par.maxTrial = 100;
	par.polyDegree = 16;
	par.tolerance = 1.0;
	
	// максимальное число "неинлайнеров"
	double outliners = 0.8;	
	
	Result res;
	res.coef = (double*)malloc((par.polyDegree+1)*sizeof(double));
	
	double* x = (double*)malloc(sizeof(double)*par.size);
	double* y = (double*)malloc(sizeof(double)*par.size);

	double u_noise1 =0.0;
	double u_noise2 =0.0;
	double n_noise =0.0;	
	
	//генерация выборки с гауссовским шумом
	for(int i=0; i<par.size; i++)
	{
		x[i] = i;
		//y[i] = 2*i*i+0.5;
		if(i%2)
			y[i] = 2*i*i+0.5+n_noise;
		else {
			u_noise1 = 1.0*rand()/RAND_MAX;
			u_noise2 = 1.0*rand()/RAND_MAX;
			n_noise = cos(2*PI*u_noise1)*sqrt(-2*log(u_noise2));			
			y[i] = 2*i*i+0.5+n_noise;
			n_noise = sin(2*PI*u_noise1)*sqrt(-2*log(u_noise2));
		}
	}

	//добавление выбросов
	for(int i=0; i<outliners*par.size; i++)
		y[rand()%par.size] += 10-(rand()%20);	
	
	//запись в файлы
	printInFile(x, par.size, "C:\\1\\x.txt");
	printInFile(y, par.size, "C:\\1\\y.txt");

	//выполнение RANSAC
	clock_t ev_time = clock();
	solveRLS(x,y, &par, &res);	
	ev_time = clock()-ev_time;
	
	printf("Total time:\t%.6f\n\n", 1.0*ev_time/CLOCKS_PER_SEC);	
	printf("Inliners = %f\n", par.size*res.inliners);	
	printf("\n\nCoefficients:\n\n");
	for(int i=par.polyDegree; i>-1; i--)
		printf("x^%d:\t%.8f\n", i, res.coef[i]);
	printf("\n");

	system("pause");
	return 0;
}

void printInFile(double *mt, int size, const char* fname){	
	FILE* pf = fopen(fname, "wb");	
	for(int i=0; i<size; i++)		
		fprintf(pf, "%.8f\n", mt[i]);
	fclose(pf);
}