# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# define f(i,j) (double)i+j+1
# define EPS 0,0000000001

int gauss(int, double* A, double* b, double* x);
double* prodMatrixVec(int, double* A, double* b);
int fill_matrix(const char* s, const int n, double** A, double** b);
double norm(int n, double* x);
double* getColumn(int n, int m, double* A, int j, double* a);
int Print_matrix (int n, const double* A);
int changeColumn(int n, int m, double* A, int j, double* a);
int constructU(int n, int m, double* x, double** U);
double* prodUVec(int n, int m, double* x, double* b);
int prodUMatrix(int n, double* b,int m, double* x, double* A);
int create_matrix (const char* s, int n);


int gauss(int n, double* A, double* b, double* x)
{
    int i;
    int k;
    double s=0;
    
    if (A[n*n-1] == 0) { printf("THE MATRIX IS SINGULAR!\n"); return 0; } 
      
    for(i=n-1; i>=0; i--)
    {
    	for(k=i+1; k<n; k++) s+=x[k]*A[i*n+k];
	x[i]=(b[i]-s)/A[i*n+i];
	s=0;
    }
    
return 1;
}

double* prodMatrixVec(int n, double* A, double* b)
{
    int i;
    int j;
    double s=0;
    double* y;
    y=(double*) malloc(n* sizeof(double));
    for(i=0; i<n; i++)
    {
    	for(j=0; j<n; j++) s+=A[i*n+j]*b[j];
	    y[i]=s;
     	s=0;
    }
    for(i=0; i<n; i++) b[i]=y[i];
    if(!y) free(y);
    return b;
}


int fill_matrix(const char* s, const int n, double** A, double** b)
{
    double* a;
    double* B;
    int i, j;
    int k;
    FILE* f;
    f=fopen(s, "r");
    if (f==NULL) return 2;
    a=(double* )malloc(n*n* sizeof(double));
    if (a==NULL){fclose(f); return 3;}
    B=(double* )malloc(n* sizeof(double));
    if (B==NULL){free(a); fclose(f); return 4;}
    for (i=0; i<n; i++)
	    {
        for(j=0; j<n; j++)
           		{
                      if ((k=fscanf(f, "%lf", a+i*n+j)) != 1) { fclose(f); free(a); free(B); return 5;}
                }
            if ((k=fscanf(f, "%lf", B+i))!=1){fclose(f); free(a); free(B); return 6;}
        }
    *A=a;
    *b=B;
    fclose(f);
//printf("FILL\n");
    return 1;
}


double norm(int n, double* x)
{
       double s=0;
       int i;
       for(i=0; i<n; i++) s+=x[i]*x[i];
return sqrt(s);
}

double* getColumn(int n, int m, double* A, int j, double* a)
{
    int i;
    for (i=m; i<n; i++) a[i]=A[i*n+j];
return a;
}


int Print_matrix (int n, const double* A){
int i, j;
for(i=0; i<n; i++)
	{
   for(j=0; j<n; j++)
            {
            printf("%lf ", A[i*n+j]);
            }
    printf("\n");
    }
return 0;
}

int changeColumn(int n, int m, double* A, int j, double* a)
{
    int i;
    for (i=m; i<n; i++) A[i*n+j]=a[i];
return 0;
}

/*
int constructU(int n, int m, double* x, double** U)
{
if(!(*U))free(*U);
double* u;
u=(double*) malloc(n*n*sizeof(double));
if (u==NULL) return 1;
int i, j;
for(i=0; i<n; i++)
{
	for(j=0; j<n; j++)
   {
   if(i<m || j<m) u[i*n+j]=0;
   else u[i*n+j]=-2*x[i]*x[j];
   if(i==j) u[i*n+j]+=1;
   }
}
*U=u;
return 0;
}
*/

double* prodUVec(int n, int m, double* x, double* b)
{
        double c=0;
        int i;
        for(i=m; i<n; i++) c+=x[i]*b[i];
        for(i=m; i<n; i++) b[i]=b[i]-2*c*x[i];
return b;
}


int prodUMatrix(int n, double* b, int m, double* x, double* A)
{        
        int i;
        for(i=m; i<n; i++)
        {
                 changeColumn(n, m, A, i, prodUVec(n, m, x, getColumn(n, m, A, i, b)));
        }
        return 0;
}

int create_matrix (const char* s, int n)
{
   int i, j;
   
   FILE* f;
   f=fopen(s, "w");
   if (f==NULL) return 2;
   for(i=0; i<n; i++)
   {
            for(j=0; j<n; j++)
            {
                     if(j<=i) fprintf(f, "%f ", f(i,j));
                     else fprintf(f, "%f ", (double)0);
                     //printf("put ");
            }
            fprintf(f, "%f ", (double)1);
            fprintf(f, "\n");
   }
   fclose(f);
   //printf("END");
return 1;
} 

int method_otrajenia(int n, double* A, double* b, double** x1)
{
	int i;
	int m;
	double temp;
	double* b1;
	double* x;
	x=(double*) malloc(n* sizeof(double));
	b1=(double*) malloc(n* sizeof(double));
	for(m=0; m<n-1; m++)
    {
//printf("I'M IN CIKL\n");
//Print_matrix(n, A);
             getColumn(n, m, A, m, x);
//printf("x on step %d\n", m+1);
//for(i=0; i<n; i++) printf("%lf\n", x[i]);
             x[m]=x[m]-norm(n-m, x+m);
             temp=norm(n-m, x+m);
//printf("x[%d] = %lf\n", m, x[m]);
             if (abs(x[m])>EPS)
                {
                for(i=m; i<n; i++) x[i]=x[i]/temp;
                prodUMatrix(n, b1, m, x, A);
                prodUVec(n, m, x, b);
//printf("m = %d\n", m+1);
//Print_matrix(n, A);
//printf("\n");
                }
    }
//printf("start print\n ");
//Print_matrix(n, A);
//printf("end printf\n\n");
    if(gauss(n, A, b, x) != 1) { printf("  THE MATRIX IS SINGULAR!\n"); return 0;};
    *x1=x;
    if(!b1) free(b1);
return 1;
}

int PrintNevyaska(int n, double* A, double* b, double* x, char* s)
{
    int i;
    if(fill_matrix(s, n, &A, &b) != 1) { printf("PROGRAM STOPPED BECAUSE PROBLEMS WITH READING MATRIX IN COUNTING NEVYASKA!!\n"); return 2;}
    prodMatrixVec(n, A, x);
    for(i=0; i<n; i++) x[i]=x[i]-b[i];
    printf("Nevyaska = %lf", norm(n, x));
    printf("\n"); 
    return 1;   
}






int main(void){
    time_t timeStart, timeEnd;
    int n;
    int i;
    int m;
    double* A;
    double* b;
    double* x;
    char* s=(char*)"inputFile.txt";
    char* s1=(char*)"inputFileRandom.txt";
    printf("Input dimension ");
    if(scanf("%d", &n) !=1) printf("WRONG INPUT DATA!\n");
    printf("\n");
    printf("chose");
    if(scanf("%d", &m) != 1) printf("WRONG INPUT DATA!\n");
    timeStart=clock();
    if (m==1)
    {
    if ((i=fill_matrix(s, n, &A, &b))!=1) { printf("wrong input data!\n"); return 1;}
    }
    else {
          if(create_matrix(s1, n) != 1) { printf("FILE NAME ERROR!!\n"); return 2;}
          fill_matrix(s1, n, &A, &b);
          //Print_matrix(n, A);
         }
    //Print_matrix(n, A); 
    if (method_otrajenia(n, A, b, &x) != 1) { printf("system was not solve!\n"); return 4; }
    //for (i=0; i<n*n; i++) printf("do%d", i);
      //  A1[i]=A[i]; b1[i]=b[i];}
    
    for(i=0; i<n; i++) printf("%lf\n", x[i]);
    //for(i=0; i<n; i++) a[i]=x[i];
    //Print_matrix(n, A1);
    
    if (m==1) 
    { 
    	if(PrintNevyaska(n ,A, b, x, s) != 1) {printf("program stopped while nevyaska!\n"); free(A); free(b); free(x); return 3;} 
    }
    else
    {
    	if(PrintNevyaska(n ,A, b, x, s1) != 1) {printf("program stopped while nevyaska!\n"); free(A); free(b); free(x); return 3;}
    }
    
    //Print_matrix(n, A);
    //for(i=0; i<n; i++) printf("\n %.5lf", x[i]);
    if(!s)free(s);
    if(!s1)free(s1);
    if(!A)free(A);
    if(!b)free(b);
    if(!x)free(x);
    timeEnd=clock();
    printf("\nTIME = %lf\n", difftime(timeEnd, timeStart)/CLOCKS_PER_SEC);
    // scanf("%d", &n);
return 0;
}
