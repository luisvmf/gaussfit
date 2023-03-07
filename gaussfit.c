#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "libgaussfit.c"
//Copyright (c) 2023 Luís Victor Muller Fabris. Apache License.

typedef struct{
	double *xdatapoints;
	double *ydatapoints;
	int Nlines;
}xydatastructmain;

void replaceDoubleCharacter(char *str,char rep){
	char *dest = str;
	while (*str != '\0'){
		while (*str == rep && *(str + 1) == rep){
			str++;
		}
		*dest++ = *str++;
    }
    *dest = '\0';
}

int hasDoubleCharacter(char *str,char rep){
	char *dest = str;
	while (*str != '\0'){
		if(*str == rep && *(str + 1) == rep){
			return 1;
		}
		*dest++ = *str++;
    }
	return 0;
}

void replaceAll(char *str, char orig, char rep) {
	int i=1;
	char *strb;
	while(i>0){
		i=0;
		while((strb = strchr(str, orig)) != NULL) {
			*strb = rep;
			i++;
		}
	}
}

int countlines(char *filename){
	FILE* fd = fopen(filename, "r");
	int i=0;
	size_t len=0;
	char *line;
	while ((getline(&line, &len, fd)) != -1){
		i=i+1;
	}
	fclose(fd);
	return i;
}

void printusage(){
	printf("Usage:\n gaussfit [-p] [-m] [-a] [-g A=XX,sigma=XX,mu=XX] [-g A=XX,sigma=XX,mu=XX] [-f file]\n\n");
	printf("   -w X  File save mode.\n             X=0-> Don't save results (default).\n             X=1-> Saves final fit parameters\n                   (filename.gaussfit).\n             X=2-> Saves final fit parameters\n                   and spectrum with boxcar and\n                   baselined removed (filename.gaussfit\n                   and filename.proc).\n");
	printf("   -b X  Baseline polynomial order (X is an integer).\n");
	printf("   -s X  Boxcar smoothing (X is an integer).\n");
	printf("   -p    Plot graph and gaussians with gnuplot.\n         Gnuplot must be installed and on path\n         for this option to be used.\n");
	printf("   -m    Manualy specify gaussians/Lorentzians initial\n         guess with option -g or -l. Must not be used\n         together with option -a.\n");
	printf("   -a    Automaticaly set gaussians initial guess.\n");
	printf("   -g A=XX,sigma=XX,mu=XX\n         Set gaussians initial guess, where XX are\n         float point values. Options -g or -l must\n         be used if option -m is used. \n");
	printf("   -l A=XX,sigma=XX,mu=XX\n         Set Lorentzian initial guess, where XX are\n         float point values. Options -g or -l must\n         be used if option -m is used. \n");
	printf("   -f file\n         Set file containing data to\n         be fitted. Must be used.\n");
	printf("\n\nEquations:\n                       2                             2\n");
	printf("               (x - μᵢ)                           A σᵢ\n");
	printf("              ___________           lᵢ(x)= ___________________\n");
	printf("                    2                          2           2\n");
	printf("                 2 σᵢ                        σᵢ  + (x - μᵢ)\n");
	printf("gᵢ(x) =    Aᵢ e                     \n");
	printf("\n\n\nf(x) =  Σ (gᵢ(x)+lᵢ(x))\n");
	printf("        ᵢ\n");
}

xydatastructmain *readfile(char *filename){
	FILE* fd = fopen(filename, "r");
	if(fd==NULL){
		perror("Error, can't open file specified in argument -f");
		printusage();
		exit(EXIT_FAILURE);
	}
	char *line;
	size_t len=0;
	int i;
	int N=countlines(filename);
	int j=0;
	int k=0;
	xydatastructmain *returndatareadfile=(xydatastructmain *)malloc(sizeof(xydatastructmain));
	returndatareadfile->xdatapoints=(double *)malloc(N*sizeof(double));
	returndatareadfile->ydatapoints=(double *)malloc(N*sizeof(double));
	returndatareadfile->Nlines=N;
	while ((getline(&line, &len, fd)) != -1){
		replaceAll(line,',','.');
		replaceAll(line,'\r','\n');
		replaceDoubleCharacter(line,'\n');
		replaceDoubleCharacter(line,' ');
		replaceAll(line,' ','\t');
		replaceDoubleCharacter(line,'\t');
		//Check if line contains invalid characters (not number, ., \t +, -, \r, \0 or \n). XXX XXX XXX This is not 100% to prevent invalid lines, but will get most of normal comment line cases.
		char allowedcharacters[17]={'0','1','2','3','4','5','6','7','8','9','.','\t','-','+','\n','\r','\0'};
		int okchar;
		int iscommentline=0;
		if(hasDoubleCharacter(line,'+')==1){
			iscommentline=1;
		}
		if(hasDoubleCharacter(line,'-')==1){
			iscommentline=1;
		}
		if(hasDoubleCharacter(line,'.')==1){
			iscommentline=1;
		}
		for(int ii=0;ii<strlen(line);ii++){
			okchar=0;
			for(int jj=0;jj<strlen(allowedcharacters);jj++){
				if(line[ii]==allowedcharacters[jj]){
					okchar=1;
				}
			}
			if(okchar==0){
				iscommentline=1;
			}
		}
		//Check if line is empty.
		if(line[0]=='\n'){
			iscommentline=1;
		}
		if(line[0]=='\r'){
			iscommentline=1;
		}
		//End check, continue only if it doesn't contain invalid characters, if it contains discard line (it is a comment line).
		if(iscommentline!=1){
			char* tmp = strdup(line);
			char* str= strtok(tmp, "\t");
			i=0;
			while(str!=NULL){
				if(i==0){
					if(j>=N){
						perror("File changed while reading.\n");
						exit(EXIT_FAILURE);
					}
					returndatareadfile->xdatapoints[j]=atof(str);
					j=j+1;
				}else{
					if(k>=N){
						perror("File changed while reading.\n");
						exit(EXIT_FAILURE);
					}
					returndatareadfile->ydatapoints[k]=atof(str);
					k=k+1;
				}
				str= strtok(NULL, "\t");
				i++;
				if(i>1){
					break;
				}
			}
			free(tmp);
			if(k>j){
				k=k-1;
				returndatareadfile->Nlines=returndatareadfile->Nlines-1;
			}
			if(j>k){
				j=j-1;
				returndatareadfile->Nlines=returndatareadfile->Nlines-1;
			}
		}else{
			returndatareadfile->Nlines=returndatareadfile->Nlines-1;
		}
    }
	return returndatareadfile;
}

int FileExists(char *file){
	FILE *fd;
	if (fd = fopen(file, "r")){
		fclose(fd);
		return 1;
	}
	return 0;
}



//XXX TODO automatic initial guess;
int main(int argc, char **argv){
	gaussfitparam myparam;
	initgaussfit(&myparam);
	int c;
	enum{
		A=0,
		sigma=1,
		mu=2,
		THE_END
	};
	const char *mount_opts[]={
		[A] = "A",
		[sigma] = "sigma",
		[mu] = "mu",
		[THE_END] = NULL
	};
	enum{
		rth=0,
		nt=1,
		faf=2,
		risaf=3,
		tyf=4,
		THE_ENDb
	};
	const char *mount_optsauto[]={
		[rth] = "rth",
		[nt] = "nt",
		[faf] = "faf",
		[risaf]="risaf",
		[tyf]="tyf",
		[THE_ENDb] = NULL
	};
	float rthv=0.2;
	int ntv=95;
	float fafv=0.95;
	int tyfv=1;
	float risafv=1.1;
	char *value=NULL;
	char *valueauto=NULL;
	int manualgaussian=-1;
	float curgaussvalues[4];
	int curaddedgaussvalues[4];
	char *file=NULL;
	int addedgauss=0;
	int plot=0;
	int baselineorder=-1;
	int boxcarorder=0;
	int writemode=0;
	while ((c = getopt (argc, argv, "b:s:a:mf:pg:w:l:")) != -1){ 
			switch (c){
				case 'b':
					baselineorder=atoi(optarg);
					break;
				case 'w':
					writemode=atoi(optarg);
					if(writemode!=0){
						if(writemode!=1){
							if(writemode!=2){
								printf("Error: Invalid -w value.\n");
								printusage();
								return -1;				
							}
						}
					}
					break;
				case 's':
					boxcarorder=atoi(optarg);
					if(boxcarorder<0){
						printf("Error: Boxcar smoothing must be positive.\n");
						printusage();
						return -1;
					}
					break;
				case 'm':
					if(manualgaussian==0){
						printf("Error: Option -m and -a can't be used at the same time.\n");
						printusage();
						return -1;
					}
					manualgaussian=1;
					break;
				case 'p':
					plot=1;
					break;
				case 'a':
					if(manualgaussian==1){
						printf("Error: Option -m and -a can't be used at the same time.\n");
						printusage();
						return -1;
					}
					manualgaussian=0;
						while (*optarg != '\0'){
							valueauto=NULL;
							switch(getsubopt (&optarg, (char *const *) mount_optsauto, &valueauto)){
								case rth:
									if(valueauto==NULL){
										break;
									}
									rthv=atof(valueauto);
									break;
								case nt:
									if(valueauto==NULL){
										break;
									}
									ntv=atoi(valueauto);
									break;
								case faf:
									if(valueauto==NULL){
										break;
									}
									fafv=atof(valueauto);
									break;
								case risaf:
									if(valueauto==NULL){
										break;
									}
									risafv=atof(valueauto);
									break;
								case tyf:
									if(valueauto==NULL){
										break;
									}
									tyfv=atoi(valueauto);
									break;
								default:
									break;
							}
						}
					break;
				case 'f':
					if(file!=NULL){
						printf("Error: Option -f used more than once.\n");
						printusage();
						return -1;
					}
					file=optarg;
					break;
				case 'g':
					curaddedgaussvalues[0]=-1;
					curaddedgaussvalues[1]=-1;
					curaddedgaussvalues[2]=-1;
					if(manualgaussian==1){
						while (*optarg != '\0'){
							value=NULL;
							switch(getsubopt (&optarg, (char *const *) mount_opts, &value)){
								case A:
									if(value==NULL){
										break;
									}
									curgaussvalues[1]=atof(value);
									curaddedgaussvalues[0]=0;
									break;
								case sigma:
									if(value==NULL){
										break;
									}
									curgaussvalues[2]=atof(value);
									curaddedgaussvalues[1]=0;
									break;
								case mu:
									if(value==NULL){
										break;
									}
									curgaussvalues[3]=atof(value);
									curaddedgaussvalues[2]=0;
									break;
								default:
									break;
							}
						}
						for(int ii=0;ii<=2;ii++){
							if(curaddedgaussvalues[ii]==-1){
								printf("Error: Invalid parameters for gaussian specified.\n");
								printusage();
								return -1;
							}
						}
						addgaussian(&myparam,curgaussvalues[1],curgaussvalues[2],curgaussvalues[3]);
						addedgauss=1;
					}else{
						printf("Error: Option -m is required to be used before option -g.\n");
						printusage();
						return -1;
					}
					break;
				case 'l':
					curaddedgaussvalues[0]=-1;
					curaddedgaussvalues[1]=-1;
					curaddedgaussvalues[2]=-1;
					if(manualgaussian==1){
						while (*optarg != '\0'){
							value=NULL;
							switch(getsubopt (&optarg, (char *const *) mount_opts, &value)){
								case A:
									if(value==NULL){
										break;
									}
									curgaussvalues[1]=atof(value);
									curaddedgaussvalues[0]=0;
									break;
								case sigma:
									if(value==NULL){
										break;
									}
									curgaussvalues[2]=atof(value);
									curaddedgaussvalues[1]=0;
									break;
								case mu:
									if(value==NULL){
										break;
									}
									curgaussvalues[3]=atof(value);
									curaddedgaussvalues[2]=0;
									break;
								default:
									break;
							}
						}
						for(int ii=0;ii<=2;ii++){
							if(curaddedgaussvalues[ii]==-1){
								printf("Error: Invalid parameters for lorentzian specified.\n");
								printusage();
								return -1;
							}
						}
						addlorentzian(&myparam,curgaussvalues[1],curgaussvalues[2],curgaussvalues[3]);
						addedgauss=1;
					}else{
						printf("Error: Option -m is required to be used before option -l.\n");
						printusage();
						return -1;
					}
					break;
				case '?':
					printusage();
					return -1;
				default:
					abort();
			}
		}
		if(file==NULL){
			printf("Error: File not specified.\n");
			printusage();
			return -1;
		}
		if(manualgaussian==-1){
			printf("Error: Option -m or -a must be used.\n");
			printusage();
			return -1;
		}
		if(manualgaussian==1){
			if(addedgauss==0){
				printf("Error: Option -m used without option -g or -l.\n");
				printusage();
				return -1;
			}
		}
		xydatastructmain *returndatareadfile=readfile(file);
		myparam.x=returndatareadfile->xdatapoints;
		myparam.y=returndatareadfile->ydatapoints;
		myparam.Npoints=returndatareadfile->Nlines;
		boxcar(&myparam,boxcarorder);
		removebaseline(&myparam,baselineorder);
		if(manualgaussian==0){
			printf("\x1b[0m\x1b[36mObtaining automatic initial guess with parameters:\x1b[0m\n\x1b[32m (rth,nt,faf,riseaf)=(%f,%i,%f,%f)\n",rthv,ntv,fafv,risafv);
			automaticguess(&myparam,rthv,ntv,fafv,risafv,tyfv);
		}
		printf("\x1b[36mGaussians initial guess:\x1b[0m\n\x1b[32m");
		int ijj=0;
		while(ijj<getNumberFittedGaussians(&myparam)){
			printf("(A,sigma,mu)=(%f,%f,%f)\n",getInitialGuess(&myparam,ijj)->A,getInitialGuess(&myparam,ijj)->sigma,getInitialGuess(&myparam,ijj)->mu);
			ijj=ijj+1;
		}
		printf("\x1b[0m\x1b[36mInitiating iterations:\x1b[0m\n\x1b[32m");
		fitgaussians(&myparam);
		printf("Iterations complete.\n");
		printf("\x1b[0m\x1b[36mResult:\x1b[0m\n");
		printf("                       2                             2\n");
		printf("              -(x - μᵢ)                           A σᵢ\n");
		printf("              ___________           lᵢ(x)= ___________________\n");
		printf("                    2                          2           2\n");
		printf("                 2 σᵢ                        σᵢ  + (x - μᵢ)\n");
		printf("gᵢ(x) =    Aᵢ e                     \n");
		printf("\n\n\nf(x) =  Σ (gᵢ(x)+lᵢ(x))\n");
		printf("        ᵢ");
		printf("\n     ┌───────────────────────────────────────────────────────────────────┐\n     │                     ...:::Final Results:::...                     │\n     ├──────────────┬─────────────────┬────────────────┬─────────────────┤\n");
		int i=0;
		for (i = 0;i<getNumberFittedGaussians(&myparam)-1;i++) {
			if(getFitResult(&myparam,i)->type==0){
				printf("     │Gaussian   %03i│ A=%+013.5f │σ=%+013.5f │μ=%+013.5f  │\n     ├──────────────┼─────────────────┼────────────────┼─────────────────┤\n",i,getFitResult(&myparam,i)->A,getFitResult(&myparam,i)->sigma,getFitResult(&myparam,i)->mu);
			}else{
				printf("     │Lorentzian %03i│ A=%+013.5f │σ=%+013.5f │μ=%+013.5f  │\n     ├──────────────┼─────────────────┼────────────────┼─────────────────┤\n",i,getFitResult(&myparam,i)->A,getFitResult(&myparam,i)->sigma,getFitResult(&myparam,i)->mu);
			}
		}
		if(getFitResult(&myparam,i)->type==0){
			printf("     │Gaussian   %03i│ A=%+013.5f │σ=%+013.5f │μ=%+013.5f  │\n",i,getFitResult(&myparam,i)->A,getFitResult(&myparam,i)->sigma,getFitResult(&myparam,i)->mu);
		}else{
			printf("     │Lorentzian %03i│ A=%+013.5f │σ=%+013.5f │μ=%+013.5f  │\n",i,getFitResult(&myparam,i)->A,getFitResult(&myparam,i)->sigma,getFitResult(&myparam,i)->mu);
		}
		printf("     └──────────────┴─────────────────┴────────────────┴─────────────────┘\n\n");
		if(writemode==1){
			char *filefinal=(char *)malloc((strlen(file)+11)*sizeof(char));
			sprintf(filefinal,"%s.gaussfit",file);
			if(FileExists(filefinal)==0){
				FILE *fd = fopen(filefinal, "w");
				fprintf(fd,"#n\tA\tσ\tμ\n");
				for (int j=0;j<getNumberFittedGaussians(&myparam);j++){
					fprintf(fd,"%i %f %f %f\n",j,getFitResult(&myparam,j)->A,getFitResult(&myparam,j)->sigma,getFitResult(&myparam,j)->mu);
				}
				fclose(fd);
			}else{
				printf("Error! File %s exists, skipping write.\n",filefinal);
			}
			free(filefinal);
		}
		if(writemode==2){
			char *filefinal=(char *)malloc((strlen(file)+11)*sizeof(char));
			sprintf(filefinal,"%s.gaussfit",file);
			if(FileExists(filefinal)==0){
				FILE *fd = fopen(filefinal, "w");
				fprintf(fd,"#n\tA\tσ\tμ\ttype(0=Gaussian,1=Lorentzian)\n");
				for (int j=0;j<getNumberFittedGaussians(&myparam);j++){
					fprintf(fd,"%i %f %f %f %i\n",j,getFitResult(&myparam,j)->A,getFitResult(&myparam,j)->sigma,getFitResult(&myparam,j)->mu,getFitResult(&myparam,i)->type);
				}
				fclose(fd);
			}else{
				printf("Error! File %s exists, skipping write.\n",filefinal);
			}
			sprintf(filefinal,"%s.proc",file);
			if(FileExists(filefinal)==0){
				FILE *fd = fopen(filefinal, "w");
				fprintf(fd,"#X\tY\n");
				for (int j=0;j<myparam.Npoints;j++){
					fprintf(fd,"%f %f\n",myparam.x[j],myparam.y[j]);
				}
				fclose(fd);
			}else{
				printf("Error! File %s exists, skipping write.\n",filefinal);
			}
			free(filefinal);
		}
		if(plot==1){
			showfitresult_gnuplot(&myparam);
		}
		return 0;
}
