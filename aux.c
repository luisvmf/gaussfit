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
		while(contador1<dadoslen){
			dadosb[contador1]=((dados[contador1-1])+(dados[contador1])+(dados[contador1+1]))/3;
			contador1=contador1+1;
		}
	}else{
		while(contador1<dadoslen){
			dadosb[contador1]=((dados[contador1])+(dados[contador1+1]))/2;
			contador1=contador1+1;
		}
	}
	dadosb[dadoslen]=(dados[dadoslen]+dados[dadoslen-1])/2;
	contador1=0;
	while(contador1<=dadoslen){
		dados[contador1]=dadosb[contador1];
		contador1=contador1+1;
	}
	free(dadosb);
}
void boxcar(double *dados, int boxcarsize,int dadoslen){
	int cnt_boxcar=0;
	int loopselect=0;
	while(cnt_boxcar<boxcarsize){
		private_boxcarint(dados,dadoslen,loopselect);
		cnt_boxcar=cnt_boxcar+1;
	}
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
void fitbaseline(double *y,double *yb, int order,int N){
	if(order>0){
		long double *coefbaseline=NULL;
		int i=0;
		int j=0;
		coefbaseline=(long double *)malloc(50000*sizeof(long double));		
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
		free(coefbaseline);
	}
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
