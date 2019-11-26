#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "struct.h"


int chkcmdline( int argc, char **argv,CMD *cmd){
        char **p;
	int num=0;
	void errmsg();
        if(argc < 2){
		errmsg();
                return(FALSE);
        }
	//default values
        p=argv;
	cmd->map_t=0.00;
	cmd->Nthr=2;
	cmd->dreso=2.00;
	cmd->MergeDist=0.50;
	cmd->Filter=0.10;
	cmd->Mode=0;
	cmd->LocalR=10.0;
	cmd->Dkeep=0.5;
	cmd->Nround=5000;
	cmd->Nnb=30;
	cmd->Ntabu=100;
	cmd->Nsim=10;
	cmd->Allow=1.01;
	cmd->Nbeam=20;
	cmd->Sym=false;
	cmd->CopyMode=false;
        cmd->Mode=2; //set MST mode as default
        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->filename,*(++p));
                	--argc; break;
		 case 'Y':
			strcpy(cmd->symfilename,*(++p));
			cmd->Sym=true;
                	--argc; break;
		 case 'C':
			cmd->CopyMode=true;
                	break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->map_t=atof(*(++p)); 
			--argc; break;
		 case 'g':
			cmd->dreso=atof(*(++p)); 
			--argc; break;
		 case 'f':
                        cmd->Filter=atof(*(++p));
                        --argc; break;
		 case 'm':
                        cmd->MergeDist=atof(*(++p));
                        --argc; break;
		 case 'R':
                        cmd->LocalR=atof(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'r':
                        cmd->Nround=atoi(*(++p));
                        --argc; break;
		 case 'b':
                        cmd->Nnb=atoi(*(++p));
                        --argc; break;
		 case 'l':
                        cmd->Ntabu=atoi(*(++p));
                        --argc; break;
		  case 's':
                        cmd->Nsim=atoi(*(++p));
                        --argc; break;
		 case 'a':
                        cmd->Allow=atoi(*(++p));
                        --argc; break;
		 case 'T':
                        cmd->Mode=0;
                        break;
		 case 'G':
                        cmd->Mode=1;
                        break;
		 case 'M':
                        cmd->Mode=2;
                        break;
		 case 'L':
                        cmd->Mode=3;
                        break;
		 case 'V':
                        cmd->Mode=4;
                        break;
		 case 'W':
                        cmd->Mode=5;
                        break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	  }
	 }
        }
	if(cmd->Sym==false){
	 puts("##Missing Rotation file!");
			errmsg(); 
			return(FALSE);
	}
	//option check
	printf("#Number of threads= %d\n",cmd->Nthr);
	printf("#Map Threshold= %f\n",cmd->map_t);
	printf("#Band Width= %f\n",cmd->dreso);
	printf("#Merge Dist= %f\n",cmd->MergeDist);
	printf("#Filtering Dens Cut= %f\n",cmd->Filter);
	printf("#Mode= %d 0:Tabu, 1:Graph, 2: MST\n",cmd->Mode);
	printf("#Num of Round= %d\n",cmd->Nround);
	printf("#Num of Tabu List= %d\n",cmd->Ntabu);
	printf("#Num of Neighbors= %d \n",cmd->Nnb);
	printf("#Num of Rounds= %d \n",cmd->Nround);
        return(TRUE);
}

void errmsg(){
	puts("Usage: MainmastSeg -i [MAP.mrc] -Y [Rotation Matrix] [(option)]");
	puts("Fast MainmastSeg C-lang & multi thread version");
	puts("---Mode---");
	puts("-L : Fast LDP search mode");
	puts("-M : Minimum Spanning Tree Mode (default)");
	puts("-G : Graph Mode");
	puts("-W : Generating MRC file Mode. File name segment\%d.mrc");
	puts("-V : Movie Mode");
	//puts("-T : Tabu Search Mode (Default)");
	puts("---Options---");
	printf("-c [int  ] :Number of cores for threads def=%d\n",2);
	printf("-t [float] :Threshold of density map def=%.3f\n",0.00);
	printf("-g [float] : bandwidth of the gaussian filter\n");
        printf("             def=2.0, sigma = 0.5*[float]\n");
	printf("-f [float] :Filtering for representative points def=%.3f\n",0.10);
	printf("-m [float] :After MeanShifting merge LDPs where d<[float] def=%.3f\n",0.50);
	printf("-R [float] :Radius of Local MST def=%.3f\n",10.0);
	printf("-k [float] :keep edges where d<[float] def=%.3f\n",0.5);
/*
	puts("---Options for Tabu Search---");
	printf("-r [int  ] :Number of rounds def=%d\n",5000);
	printf("-b [int  ] :Number of neighbors def=%d\n",30);
	printf("-l [int  ] :Size of Tabu-list def=%d\n",100);
	printf("-s [int  ] :Number of Models def=%d\n",10);
	printf("-a [float] :Total(Tree) < Total(MST)*[float] def=%f\n",1.01);
*/
	//puts("---Options for Symmetry Search---");
	//puts("-Y [Symmetry File]: For Modelling Challenge 2019!!");
	//puts("-C : Coords copy Mode");
	printf("Thi is Ver %.3f\n",VER);
}
