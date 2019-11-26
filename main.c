/*
caldep + fragment + sphere
*/
//#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "struct.h"
#include "func.h"
#include "mrc.h"
#include "sym.h"
//#include "scoring.h"

#define PDB_STRLEN 55

void malloc_error(char *a){
 fprintf(stderr,"malloc error in %s\n",a);
 exit(0);
}
double gettimeofday_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int readlist(char *fname,char **list){
 int num=0;
 FILE *fp;
 int len;
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(list[num],LIN,fp)!=NULL){
  len=strlen(list[num]);
  list[num][len-1]='\0';//ignore terminal \n
  num++;
 }
 fclose(fp);
 return TRUE;
}

int line_num(char *fname){
 int num=0;
 FILE *fp;
 char line[LIN];
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(line,LIN,fp)!=NULL){
  num++;
 }
 fclose(fp);
 return num;
}


CMD cmd;
bool *KenFlag;

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 POINTS pt;
 MRC mrc;
 GRAPH g;
 TREE mst;
 //Get Options
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }
 
 if(readmrc(&mrc,cmd.filename))
  return(0);


 SYM s;
 if(cmd.Sym==true){
  if(readsym(&s,cmd.symfilename))
	  return (0);
 
 }else{
  //put data
  s.mtx[0][0][0]=1.00; s.mtx[0][0][1]=0.00; s.mtx[0][0][2]=0.00;
  s.mtx[0][1][0]=0.00; s.mtx[0][1][1]=1.00; s.mtx[0][1][2]=0.00;
  s.mtx[0][2][0]=0.00; s.mtx[0][2][1]=0.00; s.mtx[0][2][2]=1.00;
  s.Nsym=1;
 }

 //upsampling
 if(upsampling(&mrc,cmd.map_t))
  return(0);

 printf("#Nact= %d\n",mrc.Nact);



 if(cmd.Mode==3){//fast LDPs only
  if(fastLDP(&mrc,&pt))
   return(0);
  //if(MergePoints(&mrc,&pt))
  // return(0);
  ShowLDP(&mrc,&pt);
  t4=gettimeofday_sec();
  printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
  return 0;
 }

 //Mean Shifting
 if(meanshift(&mrc,&pt))
  return(0);

 if(MergePoints(&mrc,&pt))
  return(0);

 //New for Symmetry
 int *Ctabu,Ntb;

 if(cmd.Sym==true){
	 printf("#Size...%d\n",pt.Ncd*s.Nsym);
  if((Ctabu=(int *)malloc(sizeof(int)*pt.Ncd*4*s.Nsym*s.Nsym))==NULL){
   printf("#Memory overflow...Ncd= %d * Nsym= %d\n",pt.Ncd,s.Nsym);
   return true;
  }
  if((s.tbl=(int **)malloc(sizeof(int*)*pt.Ncd*s.Nsym*2))==NULL)
	  return true;
  for(int i=0;i<pt.Ncd*s.Nsym*2;i++)
   if((s.tbl[i]=(int *)malloc(sizeof(int)*(s.Nsym)))==NULL)
	  return true;

  if(FindCorrPointsSym(&pt,&mrc,&s,Ctabu,&Ntb,cmd.CopyMode))
   return(0);
 }


 //Set Up Graph and MST
 if(SetUpGraphSym(&pt,&g,&mrc,&mst,Ctabu,Ntb,&s))
  return(0);

 //return 0;

 if(cmd.Mode==1){//Graph
  //ShowModel(&mrc,&pt);
  ShowModelChainCIF(&mrc,&pt,&g);
  ShowGraph(&g);
  return 0;
 }
 if(cmd.Mode==2){//MST
  //ShowModel(&mrc,&pt);
  ShowModelChainCIF(&mrc,&pt,&g);
  ShowTree(&g,&mst);
  return 0;
 }
 if(cmd.Mode==4){//Movie mode
  ShowModelChainCIF_Movie(&mrc,&pt,&g);
  return 0;
 }
 if(cmd.Mode==5){//Movie mode
  GenerateSegmentMrc(&mrc,&pt,&g);
  ShowModelChainCIF(&mrc,&pt,&g);//Show MSTs
  ShowTree(&g,&mst);
  return 0;
 }

  t4=gettimeofday_sec();
  printf("#FINISHED TOTAL TIME= %f\n",t4-t1);

  return(0);//STOP!!

 //Optimize Graph by Tabu seatch
 //TREE results[100];//Max 100 simu
 TREE *results;//Max 100 simu

 if(s.Nsym!=g.Nchain){
	 printf("#CHECK!! Nsym %d != Nchain %d\n",s.Nsym,g.Nchain);
	 //Update
	 g.Nchain=s.Nsym;
	 //return 0;
 }
 results=(TREE *)malloc(sizeof(TREE)*s.Nsym*cmd.Nsim*2);
 //if((KenFlag=(bool *)malloc(sizeof(bool)*g.Nnode))==NULL)
 // return 0;

 if(TabuSym(&g,&mst,results))
  return 0;

 //Volume and density data
 ShowPathSym(&mrc,&pt,&g,results,cmd.Nsim,s.Nsym);

 t4=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;

}

