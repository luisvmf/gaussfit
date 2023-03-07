#include <stdio.h>
#include <stdlib.h>
#include "lbfgs.c"
#include "lbfgs.h"
#include "gaussfitconfig.h"
//Copyright (c) 2023 Luís Victor Muller Fabris. Apache License.


//Data structures for the gaussian linked list.
typedef struct{
	double scale;
	double A;
	double sigma;
	double mu;
	int type; //0 for gaussian 1 for lorentzian
}gaussfit_listdata_private;

typedef gaussfit_listdata_private* listelement;

struct gaussfit_listnode_private{
	listelement theelement;
	struct gaussfit_listnode_private *next;
	struct gaussfit_listnode_private *previous;
};

typedef struct gaussfit_listnode_private gaussfit_listnode_privatetype;

typedef struct gaussfit_listnode_private *gaussfit_gaussfitdata_private;

//Data structure for the gaussian fit process.
typedef struct{
	double *x0scale;
	double *type; ////0 for gaussian 1 for lorentzian
	double *initialparam;
	double *x;
	double *y;
	int Ngauss;
	gaussfit_gaussfitdata_private gaussianInitialGuess;
	int Npoints;
}gaussfitparam;


void initlist(gaussfit_gaussfitdata_private *list){
	*list=(gaussfit_listnode_privatetype*)malloc(sizeof(gaussfit_listnode_privatetype));
	//The type list is a pointer, when passed by reference we have a pointer-pointer.
	(*list)->next=*list;
	(*list)->previous=*list;
}
//Returns 1 if the list is empty.
int isempty(gaussfit_gaussfitdata_private *list){
	return ((*list)->next==*list);
}
//Inserts listelement theelement into the list.
void insertlist(gaussfit_gaussfitdata_private *list,listelement theelement){
	(*list)->previous->next=(gaussfit_listnode_privatetype*)malloc(sizeof(gaussfit_listnode_privatetype));
	(*list)->previous->next->previous=(*list)->previous;
	(*list)->previous=(*list)->previous->next;
	(*list)->previous->next=*list;
	(*list)->previous->theelement=theelement;
}
//Removes the element gaussfit_listnode_privateaux from the list.
void removelist(gaussfit_gaussfitdata_private *list, gaussfit_listnode_privatetype *gaussfit_listnode_privateaux){
	if(gaussfit_listnode_privateaux==NULL){
		return;
	}
	if(gaussfit_listnode_privateaux==*list){
		puts("Error on removelist(gaussfit_gaussfitdata_private *list, gaussfit_listnode_privatetype *gaussfit_listnode_privateaux), list is empty!");
		return;
	}
	gaussfit_listnode_privateaux->previous->next=gaussfit_listnode_privateaux->next;
	gaussfit_listnode_privateaux->next->previous=gaussfit_listnode_privateaux->previous;
	free(gaussfit_listnode_privateaux->theelement);
	free(gaussfit_listnode_privateaux);
}
//Returns the list element at position pos.
gaussfit_listnode_privatetype *getlistelement(gaussfit_gaussfitdata_private *list, int pos){
	if(isempty(list)==1){
		return NULL;
	}
	gaussfit_listnode_privatetype *aux=(*list)->next;
	int i=0;
	while(aux!=*list){
		if(i==pos){
			return aux;
		}
		aux=aux->next;
		i++;
	}
	return NULL;
}

//Get number of added initial guesses or fitted gaussians.
int getNumberFittedGaussians(gaussfitparam *myparam){
	if(myparam->Ngauss<=0){
		int i=0;
		gaussfit_listnode_privatetype *aux=getlistelement(&(myparam->gaussianInitialGuess),i);
		while(aux!=NULL){
			i=i+1;
			aux=getlistelement(&(myparam->gaussianInitialGuess),i);
		}
		return i;
	}else{
		return myparam->Ngauss;
	}
}

//Get result of gaussian fit.
gaussfit_listdata_private* getInitialGuess(gaussfitparam *myparam,int i){
	gaussfit_gaussfitdata_private *list=&(myparam->gaussianInitialGuess);
	if(i>=getNumberFittedGaussians(myparam)){
		printf("Error! In getFitResult(&myparam,i) 'i' must be a integer between 0 and getNumberFittedGaussians(&myparam)-1.\n");
		exit(EXIT_FAILURE);
	}
	if(i<0){
		printf("Error! In getFitResult(&myparam,i) 'i' must be a integer between 0 and getNumberFittedGaussians(&myparam)-1.\n");
		exit(EXIT_FAILURE);
	}
	if(isempty(list)){
		puts("getInitialGuess(gaussfit_gaussfitdata_private *list): list is empty!");
		exit(EXIT_FAILURE);
	}
	gaussfit_listdata_private* dataret=(gaussfit_listdata_private *)malloc(sizeof(gaussfit_listdata_private));
	gaussfit_listnode_privatetype *aux=getlistelement(list,i);
	dataret->A=aux->theelement->A;
	dataret->sigma=aux->theelement->sigma;
	dataret->mu=aux->theelement->mu;
	return dataret;
}

//Sort sort_x and sort_y in ascending order using sort_x as index. Both sort_x and sort_y start at index zero and must have the same size. N_shell is the size of this arrays.
void shellsort(double *sort_x, double *sort_y, int N_shell){
    int i, j, k;
	float aux_shell;
    for (i = N_shell / 2; i > 0; i = i / 2){
        for (j = i; j < N_shell; j++){
            for(k = j - i; k >= 0; k = k - i){
				if(k<0)
					return;
				if(k>=N_shell)
					return;
				if((k+i)<0)
					return;
				if((k+i)>=N_shell)
					return;
                if (sort_x[k+i] >= sort_x[k]){
                    break;
                }else{
                    aux_shell = sort_x[k];
                    sort_x[k] = sort_x[k+i];
                    sort_x[k+i] = aux_shell;
                    aux_shell = sort_y[k];
                    sort_y[k] = sort_y[k+i];
                    sort_y[k+i] = aux_shell;
                }
            }
        }
    }
}

//Differentiate the spectrum.
void differentiate(double *data,double *datax,int datalen){
	int i=0;
	double *vecdiff=NULL;
	vecdiff=(double *)malloc((datalen+100)*sizeof(double));
	datalen=datalen-1;
	if(vecdiff==NULL){
		perror("Malloc failed");
		return;
	}
	while(i<datalen){
		vecdiff[i]=((data[i+1])-(data[i]))/((datax[i+1])-(datax[i]));
		i=i+1;
	}
	vecdiff[i]=vecdiff[i-1];
	i=0;
	while(i<=datalen){
		data[i]=vecdiff[i];
		i=i+1;
	}
	free(vecdiff);
}

void private_boxcarint(double *dados,int dadoslen,int loopselect){
	int contador1=1;
	double *dadosb;
	dadosb=(double *)malloc((dadoslen+50100)*sizeof(double));
	if(dadosb==NULL){
		perror("Malloc failed");
		return;
	}
	dadosb[0]=(dados[0]+dados[1])/2;
	if(loopselect==0){
		while(contador1<(dadoslen-1)){
			dadosb[contador1]=((dados[contador1-1])+(dados[contador1])+(dados[contador1+1]))/3;
			contador1=contador1+1;
		}
	}else{
		while(contador1<(dadoslen-1)){
			dadosb[contador1]=((dados[contador1])+(dados[contador1+1]))/2;
			contador1=contador1+1;
		}
	}
	dadosb[dadoslen-1]=(dados[dadoslen-1]+dados[dadoslen-2])/2;
	contador1=0;
	while(contador1<=dadoslen){
		dados[contador1]=dadosb[contador1];
		contador1=contador1+1;
	}
	free(dadosb);
}

//Applies the boxcar to the spectrum. 
void boxcar(gaussfitparam *myparam, int boxcarsize){
	double *data=myparam->y;
	int datalen=myparam->Npoints;
	shellsort(myparam->x,myparam->y,myparam->Npoints);
	int cnt_boxcar=0;
	int loopselect=0;
	while(cnt_boxcar<boxcarsize){
		private_boxcarint(data,datalen,loopselect);
		cnt_boxcar=cnt_boxcar+1;
	}
}

float private_log(float a){
	if(a==0){
		return -pow(10.0,200.0);//-Infinity
	}
	return log(a);
}


//Calculate gaussian y value at point x.
double gaussian(double x,double A,double mu,double sigma){
	double y=A*(pow(2.71828,-(pow((x-mu),2.0) /(2.0*sigma*sigma))));
	return y;
}

//Calculate lorentzian y value at point x.
double lorentzian(double x,double A,double mu,double sigma){
	double y=A*4*0.25*sigma*sigma/(4*0.25*sigma*sigma+(x-mu)*(x-mu));
	return y;
}

//Calculate (0.25*sigma*sigma+(x-mu)*(x-mu)) value at point x.
double lorentzianhelper(double x,double A,double mu,double sigma){
	double y=1.0/(4*0.25*sigma*sigma+(x-mu)*(x-mu));
	return y;
}

//Calculate (pow(2.71828,(-pow((x-mu),2.0) /(2.0*sigma*sigma)))) exponential y value at point x.
double exphelper(double x,double A,double mu,double sigma){
	double y=(pow(2.71828,(-pow((x-mu),2.0) /(2.0*sigma*sigma))));
	return y;
}

int private_fit(double* x,double* y,int count,int order, long double* coef){
	long double **A=NULL;
	long double **xpower=NULL;
	A = (long double **)malloc((count+order+300)*sizeof(long double*));
	for(int imat=0;imat<(count+order+300);imat++) A[imat]=(long double *)malloc(2*(count+order+300)*sizeof(long double));
	xpower = (long double **)malloc((2*(count+order+300))*sizeof(long double*));
	for(int imat=0;imat<(2*(count+order+300));imat++) xpower[imat]=(long double *)malloc((2*(count+order+300))*sizeof(long double));
	long double ratio=0;
	long double aux;
	int i,j,k,l;
	if(count<=order){
		for(int imat=0;imat<(2*(count+order+300));imat++) free(xpower[imat]);
		free(xpower);
		for(int imat=0;imat<(count+order+300);imat++) free(A[imat]);
		free(A);
		return -1;
	}
	order=order+1;
	for(i=0;i<count;i++){
		xpower[i][0]=1;
		for(j=1;j<=2*order;j++){
			xpower[i][j]=xpower[i][j-1]*x[i];
		}
	}
	for(i=0;i<order;i++){
		A[i][order]=0;
		for(k=0;k<count;k++){
			A[i][order]=A[i][order]+xpower[k][i]*y[k];
		}
		for(j=0;j<order;j++){
			A[i][j]=0;
			for(k=0;k<count;k++){
				A[i][j]=A[i][j]+xpower[k][j+i];
			}
		}	
	}
	int oldcount=count;
	count=order;
	i=j=k=l=0;
	for (i=0;i<count;i++){
		if(A[i][i]==0){
			l=1;
			while(A[i+l][i]==0&&(i+l)<count){
				l++;
			}
			if ((i+l)==count){
				break;
			}
			for (j=i,k=0;k<=count;k++){
				aux=A[j][k];
				A[j][k]=A[j+l][k];
				A[j+l][k]=aux;
			}
		}
		for (j=0;j<count;j++){
			if (i!=j){
				if(A[i][i]==0){
					for(int imat=0;imat<(2*(count+order+300));imat++) free(xpower[imat]);
					free(xpower);
					for(int imat=0;imat<(count+order+300);imat++) free(A[imat]);
					free(A);
					return -1;
				}
				ratio=A[j][i]/A[i][i];
				for (k=0;k<=count;k++)
					A[j][k]=A[j][k]-(A[i][k])*ratio;
			}
		}
	}
	for(int i=0;i<=order;i++){
		coef[i]=A[i][count]/A[i][i];
	}
	for(int imat=0;imat<(2*(count+order+300));imat++) free(xpower[imat]);
	free(xpower);
	for(int imat=0;imat<(count+order+300);imat++) free(A[imat]);
	free(A);
	return 0;
}

//Removes the spectrum baseline.
void removebaseline(gaussfitparam *myparam, int order){
	double *yb=myparam->y;
	int N=myparam->Npoints;
	shellsort(myparam->x,myparam->y,myparam->Npoints);
	if(order>=0){
		long double *coefbaseline=NULL;
		int i=0;
		int j=0;
		N=N-1;
		coefbaseline=(long double *)malloc(50000*sizeof(long double));
		double *y=(double *)malloc((2*N+500)*sizeof(double));
		for(i=0;i<=N;i++){
			y[i]=(double) i;
		}
		i=0;
		int res=private_fit(y,yb,N,order,coefbaseline);
		double aux=0;
		double auxb=0;
		for (i=0;i<=N;i++){
			aux=0;
			auxb=1;
			for (j=0;j<=order;j++){
				aux=aux+auxb*coefbaseline[j];
				auxb=auxb*y[i];
			}
			yb[i]=yb[i]-aux;
		}
		aux=yb[0];
		for(i=0;i<=N;i++){
			if(yb[i]<aux){
				aux=yb[i];
			}
		}
		if(aux<0){
			for(i=0;i<=N;i++){
				yb[i]=yb[i]-aux;
			}
		}
		free(coefbaseline);
		free(y);
	}
}

//Iterate fit.
static lbfgsfloatval_t evaluate(void *instanceb,const lbfgsfloatval_t *xhelper,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step){
    int i,j;
	gaussfitparam *instance=(gaussfitparam *) instanceb;
    lbfgsfloatval_t fx = 0.0;
	double gausssum;
	int k=0;
	for (j=0;j<instance->Ngauss;j++){
		g[k+0]=0;
		g[k+1]=0;
		g[k+2]=0;
		k=k+3;
	}
	double exphelpersum[3];
	double x[3];
    for (i = 0;i < instance->Npoints;i++) {
		gausssum=0;
		k=0;
		for (j=0;j<instance->Ngauss;j++){
			x[0]=(xhelper[k+0]);
			x[1]=(xhelper[k+1]);
			x[2]=(xhelper[k+2]);
			if(instance->type[j]==0){//Gaussian
				gausssum=gausssum+gaussian(instance->x[i],x[0]*(instance->x0scale[j]),x[1],x[2]);
			}else{
				//Lorentzian
				gausssum=gausssum+lorentzian(instance->x[i],x[0]*(instance->x0scale[j]),x[1],x[2]);
			}
			k=k+3;
		}
		k=0;
		for (j=0;j<instance->Ngauss;j++){
			x[0]=(xhelper[k+0]);
			x[1]=(xhelper[k+1]);
			x[2]=(xhelper[k+2]);
			if(instance->type[j]==0){//Gaussian
				exphelpersum[0]=exphelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2]);
				exphelpersum[1]=x[0]*instance->x0scale[j]*(instance->x[i]-x[1])*exphelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2])/(x[2]*x[2]);
				exphelpersum[2]=x[0]*instance->x0scale[j]*(instance->x[i]-x[1])*(instance->x[i]-x[1])*exphelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2])/(x[2]*x[2]*x[2]);
				g[k+0]=g[k+0]-2*(instance->y[i]-gausssum)*exphelpersum[0];
				g[k+1]=g[k+1]-2*(instance->y[i]-gausssum)*exphelpersum[1];
				g[k+2]=g[k+2]-2*(instance->y[i]-gausssum)*exphelpersum[2];
			}else{
				//Lorentzian
				exphelpersum[0]=0.25*4*x[2]*x[2]*lorentzianhelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2]);
				exphelpersum[1]=2*x[0]*instance->x0scale[j]*x[2]*x[2]*(instance->x[i]-x[1])*lorentzianhelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2])*lorentzianhelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2]);
				exphelpersum[2]=2*x[0]*instance->x0scale[j]*x[2]*lorentzianhelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2])-2*x[0]*instance->x0scale[j]*x[2]*x[2]*x[2]*lorentzianhelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2])*lorentzianhelper(instance->x[i],x[0]*instance->x0scale[j],x[1],x[2]);
				g[k+0]=g[k+0]-2*(instance->y[i]-gausssum)*exphelpersum[0];
				g[k+1]=g[k+1]-2*(instance->y[i]-gausssum)*exphelpersum[1];
				g[k+2]=g[k+2]-2*(instance->y[i]-gausssum)*exphelpersum[2];
			}
			k=k+3;
		}
        fx=fx+(instance->y[i]-gausssum)*(instance->y[i]-gausssum);
    }
    return fx;
}

//Print progress message if Debugiter==1.
static int progress(void *instanceb,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,int n,int k,int ls){
	gaussfitparam *instance=(gaussfitparam *) instanceb;
	if(Debugiter==1){
		printf("Iteration %d:\n", k);
		int i;
		for (i = 0;i<instance->Ngauss;i++) {
			printf("     Values %i: A=%f  σ=%f  μ=%f\n",i,(x[3*i]),(x[3*i+1]),(x[3*i+2]));
		}
		printf("  Gradient=(%f,%f,%f)\n", g[0], g[1],g[2]);
		printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
		printf("\n");
	}
    return 0;
}

//Start gaussians fit.
double fitgaussians(gaussfitparam *myparam){
	//Initialize arrays and variables according to gaussian list.
	int i=0;
	gaussfit_listnode_privatetype *aux=getlistelement(&(myparam->gaussianInitialGuess),i);
	while(aux!=NULL){
		i=i+1;
		aux=getlistelement(&(myparam->gaussianInitialGuess),i);
	}
	myparam->Ngauss=i;
	myparam->initialparam=(double *)malloc(3*myparam->Ngauss*sizeof(double));
	myparam->x0scale=(double *)malloc(myparam->Ngauss*sizeof(double));
	myparam->type=(double *)malloc(myparam->Ngauss*sizeof(double));
	shellsort(myparam->x,myparam->y,myparam->Npoints);
	//Copy data from gaussian list to array.
	int jj=0;
	int kk=0;
	while(jj<i){
		aux=getlistelement(&(myparam->gaussianInitialGuess),jj);
		myparam->x0scale[jj]=aux->theelement->scale;
		myparam->type[jj]=aux->theelement->type;
        myparam->initialparam[kk] = aux->theelement->A;
        myparam->initialparam[kk+1] = aux->theelement->mu;
        myparam->initialparam[kk+2] = aux->theelement->sigma;
		jj=jj+1;
		kk=kk+3;
	}
	//Initialize fit with liblbfgs.
	double *initialparam=myparam->initialparam;
    i=0;
	int ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(3*myparam->Ngauss);
    lbfgs_parameter_t param;
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }
	for (i = 0;i < 3*myparam->Ngauss;i++) {
		x[i]=initialparam[i];
	}
    lbfgs_parameter_init(&param);
    param.linesearch = gaussfit_linesearch;
	param.gtol=gaussfit_gtol;
	param.ftol=gaussfit_ftol;
	param.m=gaussfit_m;
	param.past=gaussfit_past;
	param.delta=gaussfit_delta;
	param.epsilon=gaussfit_epsilon;
	param.max_linesearch=gaussfit_max_linesearch;
	param.max_iterations=gaussfit_max_iterations;
	param.min_step=gaussfit_min_step;
	param.max_step=gaussfit_max_step;
	param.wolfe=gaussfit_wolfe;
	param.xtol=gaussfit_xtol;
	//Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
    ret = lbfgs(3*myparam->Ngauss, x, &fx, evaluate, progress, myparam, &param);
    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d, %s Final ", ret,lbfgs_strerror(ret));
    printf("fx = %f.\n", fx);
	for (i = 0;i < 3*myparam->Ngauss;i++) {
		initialparam[i]=(x[i]);
	}
	lbfgs_free(x);
}


//Get result of gaussian fit.
gaussfit_listdata_private* getFitResult(gaussfitparam *myparam,int i){
	double *initialparam=myparam->initialparam;
	if(i>=myparam->Ngauss){
		printf("Error! In getFitResult(&myparam,i) 'i' must be a integer between 0 and getNumberFittedGaussians(&myparam)-1.\n");
		exit(EXIT_FAILURE);
	}
	if(i<0){
		printf("Error! In getFitResult(&myparam,i) 'i' must be a integer between 0 and getNumberFittedGaussians(&myparam)-1.\n");
		exit(EXIT_FAILURE);
	}
	gaussfit_listdata_private* dataret=(gaussfit_listdata_private *)malloc(sizeof(gaussfit_listdata_private));
	dataret->A=(initialparam[3*i]);
	dataret->mu=(initialparam[3*i+1]);
	dataret->sigma=(initialparam[3*i+2]);
	dataret->type=myparam->type[i];
	return dataret;
}

//Show fited data and curve using gnuplot.
void showfitresult_gnuplot(gaussfitparam *myparam){
	FILE *fd = fopen("/dev/shm/comluisvmftestearqfitgauss.dat", "w"); 
	FILE *fdb = fopen("/dev/shm/comluisvmftestearqfitgaussl.dat", "w"); 
	for(int i=0;i<myparam->Npoints;i++){
		float gausssum=0;
		int k=0;
		float x[3];
		fprintf(fd,"%f ",myparam->x[i]);
		fprintf(fdb,"%f ",myparam->x[i]);
		for (int j=0;j<myparam->Ngauss;j++){
			x[0]=(myparam->initialparam[k+0]);
			x[1]=(myparam->initialparam[k+1]);
			x[2]=(myparam->initialparam[k+2]);
			if(myparam->type[j]==0){//Gaussian
				gausssum=gausssum+gaussian(myparam->x[i],x[0]*(myparam->x0scale[j]),x[1],x[2]);
			}else{
				//Lorentzian
				gausssum=gausssum+lorentzian(myparam->x[i],x[0]*(myparam->x0scale[j]),x[1],x[2]);
			}
			k=k+3;
		}
		fprintf(fd,"%f %f",myparam->y[i],gausssum);
		fprintf(fdb,"%f %f",myparam->y[i],gausssum);
		k=0;
		for (int j=0;j<myparam->Ngauss;j++){
			x[0]=(myparam->initialparam[k+0]);
			x[1]=(myparam->initialparam[k+1]);
			x[2]=(myparam->initialparam[k+2]);
			if(myparam->type[j]==0){//Gaussian
				fprintf(fd," %f",gaussian(myparam->x[i],x[0]*myparam->x0scale[j],x[1],x[2]));
			}else{
				//Lorentzian
				fprintf(fdb," %f",lorentzian(myparam->x[i],x[0]*myparam->x0scale[j],x[1],x[2]));
			}
			k=k+3;
		}
		fprintf(fd,"\n");
		fprintf(fdb,"\n");
	}
	fclose(fd);
	fclose(fdb);
	system("gnuplot -p -e \"plot '/dev/shm/comluisvmftestearqfitgauss.dat' using 1:2 w lines lw 3.5 lc 'gray' title 'data','/dev/shm/comluisvmftestearqfitgauss.dat' using 1:3 w lines lw 3 lc 'black' title 'Fit', for [i=4:*] '/dev/shm/comluisvmftestearqfitgauss.dat' using 1:i with lines lw 1.1 title 'Gaussian '.(i-3),  for [i=4:*] '/dev/shm/comluisvmftestearqfitgaussl.dat' using 1:i with lines lw 1.1 title 'Lorentzian '.(i-3)\"");
}

//Data structures initializations.
void initgaussfit(gaussfitparam *myparam){
	initlist(&(myparam->gaussianInitialGuess));
	myparam->Ngauss=-1;
}

//Add gaussian initial guess to fit.
void addgaussian(gaussfitparam *myparam, double A, double sigma, double mu){
	myparam->Ngauss=-1;
	gaussfit_listdata_private *mydata=(gaussfit_listdata_private *)malloc(sizeof(gaussfit_listdata_private));
	mydata->scale=1;//Scale of variable A.
	mydata->A=A;
	mydata->sigma=sigma;
	mydata->mu=mu;
	mydata->type=0;
	insertlist(&(myparam->gaussianInitialGuess),mydata);
}

//Add lorentzian initial guess to fit.
void addlorentzian(gaussfitparam *myparam, double A, double sigma, double mu){
	myparam->Ngauss=-1;
	gaussfit_listdata_private *mydata=(gaussfit_listdata_private *)malloc(sizeof(gaussfit_listdata_private));
	mydata->scale=1;//Scale of variable A.
	mydata->A=A;
	mydata->sigma=sigma;
	mydata->mu=mu;
	mydata->type=1;
	insertlist(&(myparam->gaussianInitialGuess),mydata);
}

int private_findpeaks(double *xb,double *yb,double* Apeaks, double *sigmapeaks, double* mupeaks, float risingthreshold,int N,float fallafterpeak,float riseaftermin){
	double *xbb=(double *)malloc((2*N+100)*sizeof(double));
	double *ybb=(double *)malloc((2*N+100)*sizeof(double));
	double *adjustsection=(double *)malloc((2*N+100)*sizeof(double));
	double *adjustsectionx=(double *)malloc((2*N+100)*sizeof(double));
	int *zeroderivativepoints=(int *)malloc((2*N+100)*sizeof(double));
	int j=0;
	double valmaxderivada=yb[0];
	int zeroderivativepointscount=0;
	int realpeakscount=0;
	//-------------------------------------------
	//-------------------------------------------
	//Here we update the processed spectrum to the auxiliar arrays "xbb" and "ybb".
	for (j=0;j<N;j++){
		ybb[j]=yb[j];
		xbb[j]=xb[j];
	}
	//-------------------------------------------
	//-------------------------------------------
	differentiate(ybb,xbb,N);
	for (j=0;j<N;j++){
		if(yb[j]>valmaxderivada){
			valmaxderivada=yb[j];
		}
	}
	double toplast;
	double minok;
	int init;
	if(xbb[0]>xbb[1]){//This if is used to handle inverted spectra correctly.
		init=-1;
		toplast=yb[0];
		minok=yb[0];
		for (j=4;j<N-6;j++){
			if(yb[j]<minok){
				minok=yb[j];
			}
			if(ybb[j]<0){
				if(ybb[j+1]>0){
					//Test concavity.
					if(yb[j]<=yb[j+1]){
						if(yb[j+2]<=yb[j+1]){
							if(yb[j+1]>risingthreshold*valmaxderivada){
								if(init==-1){
									toplast=yb[j];
									init=1;
								}
								if(minok<=fallafterpeak*toplast){
									if(yb[j]>riseaftermin*minok){
										minok=yb[j];
										toplast=yb[j];
										zeroderivativepoints[zeroderivativepointscount]=j;
										zeroderivativepointscount=zeroderivativepointscount+1;
									}
								}
							}
						}
					}
				}
			}
		}
	}else{//This else is used to handle inverted spectra correctly.
		init=-1;
		toplast=yb[0];
		minok=yb[0];
		for (j=4;j<N-6;j++){
			if(yb[j]<minok){
				minok=yb[j];
			}
			if(ybb[j]>0){
				if(ybb[j+1]<0){
					//Test concavity.
					if(yb[j]<=yb[j+1]){
						if(yb[j+2]<=yb[j+1]){
							if(yb[j+1]>risingthreshold*valmaxderivada){
								if(init==-1){
									toplast=yb[j];
									init=1;
								}
								if(minok<=fallafterpeak*toplast){
									if(yb[j]>riseaftermin*minok){
										minok=yb[j];
										toplast=yb[j];
										zeroderivativepoints[zeroderivativepointscount]=j;
										zeroderivativepointscount=zeroderivativepointscount+1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for (j=0;j<zeroderivativepointscount;j++){
		double highpos=yb[zeroderivativepoints[j]-3];
		int jkk=-3;
		int thebigpeakindex=-3;
		while(jkk<4){
			if(yb[zeroderivativepoints[j]+jkk]>highpos){
				highpos=yb[zeroderivativepoints[j]+jkk];
				thebigpeakindex=jkk;
			}
			jkk=jkk+1;
		}
		adjustsection[0]=yb[zeroderivativepoints[j]+thebigpeakindex-1];
		adjustsection[1]=yb[zeroderivativepoints[j]+thebigpeakindex];
		adjustsection[2]=yb[zeroderivativepoints[j]+thebigpeakindex+1];
		//The spectrum should have only positive values for the gaussian fit.
		while(1){
			int changedfdhdhfdfh=1;
			if(adjustsection[0]<=0){
				adjustsection[0]=adjustsection[0]+50000;
				adjustsection[1]=adjustsection[1]+50000;
				adjustsection[2]=adjustsection[2]+50000;
				changedfdhdhfdfh=0;
			}
			if(adjustsection[1]<=0){
				adjustsection[0]=adjustsection[0]+50000;
				adjustsection[1]=adjustsection[1]+50000;
				adjustsection[2]=adjustsection[2]+50000;
				changedfdhdhfdfh=0;
			}
			if(adjustsection[2]<=0){
				adjustsection[0]=adjustsection[0]+50000;
				adjustsection[1]=adjustsection[1]+50000;
				adjustsection[2]=adjustsection[2]+50000;
				changedfdhdhfdfh=0;
			}
			if(changedfdhdhfdfh==1){
				break;
			}
		}
		adjustsectionx[0]=xb[zeroderivativepoints[j]+thebigpeakindex-1];
		adjustsectionx[1]=xb[zeroderivativepoints[j]+thebigpeakindex];
		adjustsectionx[2]=xb[zeroderivativepoints[j]+thebigpeakindex+1];
		double x1=adjustsectionx[0];
		double x2=adjustsectionx[1];
		double x3=adjustsectionx[2];
		double y1=adjustsection[0];
		double y2=adjustsection[1];
		double y3=adjustsection[2];
		int ok=1;
		double mu_guess=1.0;
		double k=1;
		if(y2==0){
			k=pow(10.0,200.0);//Infinity
			ok=0;
		}
		if(y3==0){
			k=pow(10.0,200.0);//Infinity
			ok=0;
		}
		if(ok==1){
			if(log(y1/y2)==0){
				k=pow(10.0,200.0);//Infinity
				ok=0;
			}
		}
		if(ok==1){
			k=(private_log(y1/y3)/private_log(y1/y2));
		}
		if((-2*x2*k+2*x1*k+2*x3-2*x1)!=0){
			mu_guess=(x3*x3-x1*x1+(k*(x1*x1-x2*x2)))/(-2*x2*k+2*x1*k+2*x3-2*x1);
		}
		ok=1;
		if(y2==0){
			k=pow(10.0,200.0);//Infinity
			ok=0;
		}
		if(ok==1){
			k=private_log(y1/y2);
		}
		double sigma_guess=0;
		if(k!=0){
			sigma_guess=sqrt(((x2-mu_guess)*(x2-mu_guess)-(x1-mu_guess)*(x1-mu_guess))/(2*k));
		}else{
			sigma_guess=pow(10.0,200.0);//Infinity
		}
		double A_guess=0;
		if(sigma_guess!=0){
			A_guess=y3*pow(2.71828,(x3-mu_guess)*(x3-mu_guess)/(2*sigma_guess*sigma_guess));
		}else{
			A_guess=0;
		}
		Apeaks[realpeakscount]=A_guess;
		mupeaks[realpeakscount]=mu_guess;
		sigmapeaks[realpeakscount]=sigma_guess;
		realpeakscount=realpeakscount+1;
	}
	free(xbb);
	free(ybb);
	free(adjustsection);
	free(adjustsectionx);
	free(zeroderivativepoints);
	return realpeakscount;
}


//Done, documment improve noise tolerance with some other threshold, ex porcentage of fall after peak, noisetol does a boxcar without affecting the spectra outside of the function-> added float fallafterpeak,float riseaftermin.
void automaticguess(gaussfitparam *myparam, float risingthreshold, int noisetol,float fallafterpeak,float riseaftermin,int type){
	int N=myparam->Npoints;
	double *Apeaks=(double *)malloc((2*N+100)*sizeof(double));
	double *sigmapeaks=(double *)malloc((2*N+100)*sizeof(double));
	double *mupeaks=(double *)malloc((2*N+100)*sizeof(double));
	double *x_aux=(double *)malloc((2*N+100)*sizeof(double));
	double *y_aux=(double *)malloc((2*N+100)*sizeof(double));
	for (int j=0;j<N;j++){
		x_aux[j]=myparam->x[j];
		y_aux[j]=myparam->y[j];
	}
	double *data=y_aux;
	int datalen=myparam->Npoints;
	shellsort(x_aux,y_aux,N);
	int cnt_boxcar=0;
	int loopselect=0;
	while(cnt_boxcar<noisetol){
		private_boxcarint(data,datalen,loopselect);
		cnt_boxcar=cnt_boxcar+1;
	}
	int peakscount=private_findpeaks(x_aux, y_aux, Apeaks, sigmapeaks, mupeaks, risingthreshold, N, fallafterpeak, riseaftermin);
	for(int i=0;i<peakscount;i++){
		if(type==0){
			addgaussian(myparam,Apeaks[i],sigmapeaks[i],mupeaks[i]);
		}else{
			addlorentzian(myparam,Apeaks[i],sigmapeaks[i],mupeaks[i]);
		}
	}
	free(Apeaks);
	free(sigmapeaks);
	free(mupeaks);
}
