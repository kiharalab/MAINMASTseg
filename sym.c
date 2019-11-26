//Symmetry Operation
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "struct.h"
#include "mrc.h"
#include "sym.h"

extern CMD cmd;

bool readsym(SYM *s,char *filename){
	int i,j; 
        FILE *fpin; 
        char line[LIN], buf[LIN]; 
        float x,y,z;
	float t;
	int a,mod;
	int cnt,pos;
        
        if((fpin=fopen(filename,"r")) == NULL){ 
                fprintf(stderr,"Can't open %s\n",filename); 
                return(-1); 
        }

	cnt=0;
        while(fgets(line,LIN,fpin)){
	 if(!strncmp(line,"REMARK",3)){
   	  sscanf(line,"REMARK %d %s %d %f %f %f %f",&a,buf,&mod,&x,&y,&z,&t);
	  pos=cnt%3;
	  printf("#MOD %d pos=%d %f %f %f %f\n",mod,pos,x,y,z,t);
	  s->mtx[mod-1][pos][0]=x;
	  s->mtx[mod-1][pos][1]=y;
	  s->mtx[mod-1][pos][2]=z;
	  s->t[mod-1][pos]=t;
	  cnt++;
 	 }
	}
	s->Nsym=mod;
	fclose(fpin);

	printf("#Nsym= %d\n",s->Nsym);
	return false;
}

bool readCIF(CIF *out,char *filename){
	int i,j; 
        FILE *fpin; 
        char line[LIN], buf[LIN]; 
        float x,y,z;
	float t;
	int a,mod;
	int cnt,pos;
        
        if((fpin=fopen(filename,"r")) == NULL){ 
                fprintf(stderr,"Can't open %s\n",filename); 
                return(-1); 
        }

	cnt=0;
        while(fgets(line,LIN,fpin)){
	 if(!strncmp(line,"REMARK",3)){
   	  sscanf(line,"REMARK %d %s %d %f %f %f %f",&a,buf,&mod,&x,&y,&z,&t);
 	 }
	}
	fclose(fpin);

	return false;
}


bool MergePointsSym(MRC *m,POINTS *p, int Nsym){
 double dcut=cmd.MergeDist/m->widthx;
 double d2cut=dcut*dcut;
 double rdcut=cmd.Filter;
 double dmax,dmin,drange;
 bool *stock;
 double sym_cd[10][3];
 int *tmp_member,*tmp_member2;
 double ang_unit=2*PI/(double)Nsym;

 if((stock=(bool *)malloc(sizeof(bool)*p->Ncd))==NULL)
  return true;
 if((tmp_member=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 if((tmp_member2=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 dmax=0;dmin=999999.99;

 for(int i=0;i<p->Ncd;i++){
  if(p->dens[i]<dmin)
   dmin=p->dens[i];
  if(p->dens[i]>dmax)
   dmax=p->dens[i];
 }
 drange=dmax-dmin;
 double rv_range=1.00/drange;
 printf("#C%d Symmetry\n",Nsym);
 printf("#dmax= %f dmin= %f\n",dmax,dmin);
 //init member
 if((p->member=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 for(int i=0;i<p->Ncd;i++)
  p->member[i]=i;
 for(int i=0;i<p->Ncd;i++)
  stock[i]=true;

 //#pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<p->Ncd-1;i++){
  double tmp[3],d2;
  if((p->dens[i]-dmin)*rv_range < rdcut)
   stock[i]=false;
  if(stock[i]==false)
   continue;

  //Generate Symmetry COORDs
  double ang=0;
  for(int sym=0;sym<Nsym;sym++){
   //rotate around Z-axis
   ang=ang_unit*(double)sym;
   sym_cd[sym][0]=cos(ang)*(p->cd[i][0]-m->xdim*0.5)
		 -sin(ang)*(p->cd[i][1]-m->ydim*0.5)
		 +m->xdim*0.5;// X
   sym_cd[sym][1]=sin(ang)*(p->cd[i][0]-m->xdim*0.5)
		 +cos(ang)*(p->cd[i][1]-m->ydim*0.5)
		 +m->ydim*0.5;// Y
   sym_cd[sym][2]=p->cd[i][2];//Keep Z
   //printf("%f %f %f\n",sym_cd[sym][0],sym_cd[sym][1],sym_cd[sym][2]);
  }
  //return false;

  for(int j=i+1;j<p->Ncd;j++){
   if(stock[j]==false)
    continue;


	for(int sym=0;sym<Nsym;sym++){
	 tmp[0]=sym_cd[sym][0]-p->cd[j][0];
  	 tmp[1]=sym_cd[sym][1]-p->cd[j][1];
  	 tmp[2]=sym_cd[sym][2]-p->cd[j][2];
  	 d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];

	 if(d2<d2cut){
    	  //Keep high dens
    	  if(p->dens[i]>p->dens[j]){
    	   stock[j]=false;
    	   p->member[j]=i;
    	  }else{
    	   stock[i]=false;
    	   p->member[i]=j;
    	   break;
    	  }
   	 }
	 if(stock[j]==false||stock[i]==false)
	  break;
	}
	if(stock[i]==false)
	  break;
/*
   tmp[0]=p->cd[i][0]-p->cd[j][0];
   tmp[1]=p->cd[i][1]-p->cd[j][1];
   tmp[2]=p->cd[i][2]-p->cd[j][2];
   d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];


   if(d2<d2cut){
    //Keep high dens
    if(p->dens[i]>p->dens[j]){
     stock[j]=false;
     p->member[j]=i;
    }else{
     stock[i]=false;
     p->member[i]=j;
     break;
    }
   }
*/
  }
 }
 //printf("62809; %d\n",p->member[62809]);
 //Update member data
 for(int i=0;i<p->Ncd;i++){
  int now=p->member[i];
  for(int j=0;j<p->Ncd;j++){
   if(now==p->member[now])
    break;
   now=p->member[now];
  }
  p->member[i]=now;
 }
 //printf("62809; %d\n",p->member[62809]);
 //Copy
 int Nmerge=0;
 for(int i=0;i<p->Ncd;i++){
  if(stock[i]){
   p->cd[Nmerge][0]=p->cd[i][0];
   p->cd[Nmerge][1]=p->cd[i][1];
   p->cd[Nmerge][2]=p->cd[i][2];
   p->dens[Nmerge]=p->dens[i];
   tmp_member[i]=Nmerge;
   Nmerge++;
  }else{
   tmp_member[i]=-1;
  }
 }
 for(int i=0;i<p->Ncd;i++)
  tmp_member2[i]=tmp_member[p->member[i]];
 for(int i=0;i<p->Ncd;i++)
  p->member[i]=tmp_member2[i];
 
 //printf("62809 tmp; %d\n",tmp_member[62809]);
 //printf("62809; %d\n",p->member[62809]);
 printf("#After Merge: %d\n",Nmerge);



 //Copy and Rotate cd
 #pragma omp parallel for schedule(dynamic,5)
 for(int sym=1;sym<Nsym;sym++){
  int ii;
  double ang=0;
 	for(int i=0;i<Nmerge;i++){
	 ii=Nmerge+i;
	 ang=ang_unit*(double)sym;
  	 p->cd[ii][0]=cos(ang)*(p->cd[i][0]-m->xdim*0.5)
		 -sin(ang)*(p->cd[i][1]-m->ydim*0.5)
		 +m->xdim*0.5;// X
  	 p->cd[ii][1]=sin(ang)*(p->cd[i][0]-m->xdim*0.5)
		 +cos(ang)*(p->cd[i][1]-m->ydim*0.5)
		 +m->ydim*0.5;// Y
  	 p->cd[ii][2]=p->cd[i][2];//Keep Z

	 p->dens[ii]=p->dens[i];

 	}
 }
 p->Ncd=Nmerge*Nsym;
 printf("#After Symmtetry: %d\n",p->Ncd);
 return false;
}

bool FindCorrPointsCn(POINTS *p,MRC *mrc,int Nsym, int *tabu,int *Ntb){
 int i,j,dif;
 double dcut=cmd.MergeDist/mrc->widthx;
 double d2cut=dcut*dcut;
 int Ntabu=0;
 puts("#Chek Circle Symmetry");

 //Find Correspondig LDPs by Cn symmetry
 double ang_unit=2*PI/(double)Nsym;
 int Nc=0;

 //if((tabu=(int *)malloc(sizeof(int)*p->Ncd*p->Ncd*2))==NULL)
 // return true;


 for(int i=0;i<p->Ncd;i++){
  double sym_cd[10][3],d2;
  double ang=0;
  double tmp[3];

  for(int sym=0;sym<Nsym;sym++){
   //rotate around Z-axis
   ang=ang_unit*(double)sym;
   sym_cd[sym][0]=cos(ang)*(p->cd[i][0]-mrc->xdim*0.5)
		 -sin(ang)*(p->cd[i][1]-mrc->ydim*0.5)
		 +mrc->xdim*0.5;// X
   sym_cd[sym][1]=sin(ang)*(p->cd[i][0]-mrc->xdim*0.5)
		 +cos(ang)*(p->cd[i][1]-mrc->ydim*0.5)
		 +mrc->ydim*0.5;// Y
   sym_cd[sym][2]=p->cd[i][2];//Keep Z
  }
 	for(int j=i+1;j<p->Ncd;j++){
	 for(int sym=0;sym<Nsym;sym++){
	  tmp[0]=sym_cd[sym][0]-p->cd[j][0];
  	  tmp[1]=sym_cd[sym][1]-p->cd[j][1];
  	  tmp[2]=sym_cd[sym][2]-p->cd[j][2];
  	  d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
	  if(d2<d2cut){
	   tabu[Ntabu*2]=i;
	   tabu[Ntabu*2+1]=j;
	   Ntabu++;
	   break;
    	  }
   	 }
  	}
 }
 *Ntb=Ntabu;
 printf("#Nchain-tabu= %d\n",Ntabu/2);
 return false;
}

bool FindCorrPointsDn(POINTS *p,MRC *mrc,int Nsym, int *tabu,int *Ntb){
 int i,j,dif;
 double dcut=cmd.MergeDist/mrc->widthx;
 double d2cut=dcut*dcut;
 int Ntabu=0;
 puts("#Chek Dihedral Symmetry");

 //Find Correspondig LDPs by Cn symmetry
 double ang_unit=2*PI/(double)Nsym;
 int Nc=0;

 //if((tabu=(int *)malloc(sizeof(int)*p->Ncd*p->Ncd*2))==NULL)
 // return true;


 for(int i=0;i<p->Ncd;i++){
  double sym_cd[20][3],d2;
  double ang=0;
  double tmp[3];

  for(int sym=0;sym<Nsym;sym++){
   //rotate around Z-axis
   ang=ang_unit*(double)sym;
   sym_cd[sym][0]=cos(ang)*(p->cd[i][0]-mrc->xdim*0.5)
		 -sin(ang)*(p->cd[i][1]-mrc->ydim*0.5)
		 +mrc->xdim*0.5;// X
   sym_cd[sym][1]=sin(ang)*(p->cd[i][0]-mrc->xdim*0.5)
		 +cos(ang)*(p->cd[i][1]-mrc->ydim*0.5)
		 +mrc->ydim*0.5;// Y
   sym_cd[sym][2]=p->cd[i][2];//Keep Z

   //Flip around X-axis
   sym_cd[Nsym+sym][0]= sym_cd[sym][0];
   sym_cd[Nsym+sym][1]=-sym_cd[sym][1]+mrc->ydim;
   sym_cd[Nsym+sym][2]=-sym_cd[sym][2]+mrc->zdim;
  }
 	for(int j=i+1;j<p->Ncd;j++){
	 for(int sym=0;sym<Nsym*2;sym++){
	  tmp[0]=sym_cd[sym][0]-p->cd[j][0];
  	  tmp[1]=sym_cd[sym][1]-p->cd[j][1];
  	  tmp[2]=sym_cd[sym][2]-p->cd[j][2];
  	  d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
	  if(d2<d2cut){
	   tabu[Ntabu*2]=i;
	   tabu[Ntabu*2+1]=j;
	   Ntabu++;
	   break;
    	  }
   	 }
  	}
 }
 *Ntb=Ntabu;
 printf("#Nchain-tabu= %d\n",Ntabu/2);
 return false;
}

bool FindCorrPointsSym(POINTS *p,MRC *mrc,SYM *s, int *tabu,int *Ntb,bool cmode){
 int i,j,dif;
 int Nsym=s->Nsym;
 double dcut=cmd.MergeDist/mrc->widthx;
 double d2cut=dcut*dcut;
 int Ntabu=0;
 int *used;
 if((used=(int *)calloc(sizeof(int),Nsym*p->Ncd))==NULL)
	 return true;
 if((p->ignore=(bool *)calloc(sizeof(bool),Nsym*p->Ncd))==NULL)
	 return true;

 puts("#Chek Symmetry");

 //Find Correspondig LDPs by Symmetry Matrix

 for(int i=0;i<p->Ncd;i++){
  p->ignore[i]=false;
	 if(used[i]==1)
	  continue;

  double sym_cd[100][3],d2;//up to 100-mer
  double ang=0;
  double tmp[3];

  //make coordinates
  for(int sym=0;sym<Nsym;sym++){
   //rotate based on mtx
	  
   tmp[0]=p->cd[i][0]*mrc->widthx+mrc->orgxyz[0];
   tmp[1]=p->cd[i][1]*mrc->widthx+mrc->orgxyz[1];
   tmp[2]=p->cd[i][2]*mrc->widthx+mrc->orgxyz[2];
	  
   sym_cd[sym][0]=tmp[0]*s->mtx[sym][0][0]
	   +tmp[1]*s->mtx[sym][0][1]
	   +tmp[2]*s->mtx[sym][0][2]
	   +s->t[sym][0];
   sym_cd[sym][1]=tmp[0]*s->mtx[sym][1][0]
	   +tmp[1]*s->mtx[sym][1][1]
	   +tmp[2]*s->mtx[sym][1][2]
	   +s->t[sym][1];
   sym_cd[sym][2]=tmp[0]*s->mtx[sym][2][0]
	   +tmp[1]*s->mtx[sym][2][1]
	   +tmp[2]*s->mtx[sym][2][2]
	   +s->t[sym][2];
   
   sym_cd[sym][0]=(sym_cd[sym][0]-mrc->orgxyz[0])/mrc->widthx;
   sym_cd[sym][1]=(sym_cd[sym][1]-mrc->orgxyz[1])/mrc->widthx;
   sym_cd[sym][2]=(sym_cd[sym][2]-mrc->orgxyz[2])/mrc->widthx;

  }
  //Ignore center 
  bool cent_flag=false;
  for(int sym1=0;sym1<Nsym;sym1++){
  	for(int sym2=sym1+1;sym2<Nsym;sym2++){
	 tmp[0]=sym_cd[sym1][0]-sym_cd[sym2][0];
	 tmp[1]=sym_cd[sym1][1]-sym_cd[sym2][1];
	 tmp[2]=sym_cd[sym1][2]-sym_cd[sym2][2];
	 d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
	 if(d2<d2cut){
	 	 cent_flag=true;
		 break;
	 }
	}      	
  }
  if(cent_flag==true){
	  p->ignore[i]=true;
	  continue;
  }
  
  //puts("#CHECKING.. SYM_CHECK");
  for(int sym=1;sym<Nsym;sym++){
	bool flag=false;
	double d2min=d2cut;
	int minid=-1;
   //printf("#Sym%d i=%d / %d\n",sym,i,p->Ncd*s->Nsym*2);
	s->tbl[i][sym-1]=-1; //<<<<<<<HERE
   //printf("#*Sym%d i=%d / %d\n",sym,i,p->Ncd*s->Nsym*2);
	  //find closest LDPs
	
   	for(int j=0;j<p->Ncd;j++){
	   //if(used[j]==1)//used
	   //	   continue;

  	 tmp[0]=sym_cd[sym][0]-p->cd[j][0];
	 if(tmp[0]>dcut||tmp[0]<-dcut)
	  continue;
  	 tmp[1]=sym_cd[sym][1]-p->cd[j][1];
  	 tmp[2]=sym_cd[sym][2]-p->cd[j][2];

  	 d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];

	 if(d2<d2min){
		 d2min=d2;
		 minid=j;
	 }
    	}
	 //printf("#minid= %d\n",minid);
	if(minid!=-1){//Find close point
	 tabu[Ntabu*2]=i;
	 tabu[Ntabu*2+1]=minid;
	 Ntabu++;
	 s->tbl[i][sym-1]=minid;

		 if(cmode==true){
		  //copy coords
		  p->cd[minid][0]=sym_cd[sym][0];
		  p->cd[minid][1]=sym_cd[sym][1];
		  p->cd[minid][2]=sym_cd[sym][2];
		 }
	 //printf("#Ntabu= %d\n",Ntabu);
	}else if(cmode==true){//No same points, copy coords
    	 //printf("%d sym %d NO MEMBER!!\n",i,sym);
    	 //Add
    	 j=p->Ncd;
    	 p->Ncd++;
	 
	 p->cd[j]=(double*)malloc(sizeof(double)*3);
    	 p->cd[j][0]=sym_cd[sym][0];
         p->cd[j][1]=sym_cd[sym][1];
         p->cd[j][2]=sym_cd[sym][2];
	
	 tabu[Ntabu*2]=i;
         tabu[Ntabu*2+1]=j;
         //used[j]=1;
	 s->tbl[i][sym-1]=j;

         Ntabu++;
	 //printf("#Ntabu= %d\n",Ntabu);
   	}
  }
 }
 puts("#FIN SYM_CHECK");
 //Complete Table
 for(int i=0;i<p->Ncd;i++){
  //if(used[i]==0)
  // continue;
  if(p->ignore[i]==true)
   continue;

  double sym_cd[100][3],d2;//up to 100mer
  double ang=0;
  double tmp[3];

  for(int sym=0;sym<Nsym;sym++){
   tmp[0]=p->cd[i][0]*mrc->widthx+mrc->orgxyz[0];
   tmp[1]=p->cd[i][1]*mrc->widthx+mrc->orgxyz[1];
   tmp[2]=p->cd[i][2]*mrc->widthx+mrc->orgxyz[2];
	  
   sym_cd[sym][0]=tmp[0]*s->mtx[sym][0][0]
	   +tmp[1]*s->mtx[sym][0][1]
	   +tmp[2]*s->mtx[sym][0][2]
	   +s->t[sym][0];
   sym_cd[sym][1]=tmp[0]*s->mtx[sym][1][0]
	   +tmp[1]*s->mtx[sym][1][1]
	   +tmp[2]*s->mtx[sym][1][2]
	   +s->t[sym][1];
   sym_cd[sym][2]=tmp[0]*s->mtx[sym][2][0]
	   +tmp[1]*s->mtx[sym][2][1]
	   +tmp[2]*s->mtx[sym][2][2]
	   +s->t[sym][2];
   
   sym_cd[sym][0]=(sym_cd[sym][0]-mrc->orgxyz[0])/mrc->widthx;
   sym_cd[sym][1]=(sym_cd[sym][1]-mrc->orgxyz[1])/mrc->widthx;
   sym_cd[sym][2]=(sym_cd[sym][2]-mrc->orgxyz[2])/mrc->widthx;
  }
  for(int sym=1;sym<Nsym;sym++){
  	if(s->tbl[i][sym-1]!=-1)//have data?
  	 continue;
	  bool flag=false;
	  double dmin=d2cut;
	  int minid=-1;
   	for(int j=0;j<p->Ncd;j++){
  	 tmp[0]=sym_cd[sym][0]-p->cd[j][0];
  	 tmp[1]=sym_cd[sym][1]-p->cd[j][1];
  	 tmp[2]=sym_cd[sym][2]-p->cd[j][2];
  	 d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
	 if(dmin>d2){
		minid=j;
		dmin=d2;
	 }
   	}
	if(minid!=-1){
	 tabu[Ntabu*2]=i;
	 tabu[Ntabu*2+1]=minid;
	 Ntabu++;
	 flag=true;
	 s->tbl[i][sym-1]=minid;
	}
   	//if(minid==-1)//No cd???
	   //printf("##Cannot find SYM COORDs %d sym%d dmin= %f\n",i,sym,dmin);
  	}
 }
 *Ntb=Ntabu;
 printf("#Nchain-tabu= %d\n",Ntabu/2);
 return false;
}


int cmp_simple(const void *c1, const void *c2){

 SIMPLE_TBL a=*(SIMPLE_TBL *)c1;
 SIMPLE_TBL b=*(SIMPLE_TBL *)c2;

 if(a.Nmem<b.Nmem) return 1;
 if(a.Nmem>b.Nmem) return -1;
 return 0;
}


//Make Graph and MST
bool SetUpGraphSym(POINTS *p, GRAPH *g,MRC *mrc,TREE *mst,int *tabu,int Ntabu,SYM *s){

 int i,j,dif;
 double dcut=cmd.LocalR/mrc->widthx;
 double d2cut=dcut*dcut;
 int *cid,*tmp_tb;
 int Nsym=s->Nsym;
 //printf("Ntabu= %d\n",Ntabu);

 //MALLOC GRAPH------
 
 if((g->adj=(bool **)malloc(sizeof(bool *)*p->Ncd))==NULL)
  return true;
 for(i=0;i<p->Ncd;i++){
  if((g->adj[i]=(bool *)malloc(sizeof(bool)*p->Ncd))==NULL)
   return true;
  for(j=0;j<p->Ncd;j++)
   g->adj[i][j]=false;
 }
 
 if((g->cid=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 if((cid=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;

 if((g->node=(NODE *)malloc(sizeof(NODE)*p->Ncd))==NULL)
  return true;

 if((tmp_tb=(int *)malloc(sizeof(int)*Ntabu*2))==NULL)
  return true;

 puts("#Fin Malloc");
 //END---------------


 g->Nnode=p->Ncd;//*******
 int Ne=0;
 
 #pragma omp parallel for reduction(+:Ne) schedule(dynamic,5)
 for(int ii=0;ii<p->Nori;ii++){
  int m1=p->member[ii];
  int m2;

  if(m1==-1)
   continue;
  if(p->ignore[m1]==true)
	  continue;

  for(int jj=ii+1;jj<p->Nori;jj++){
   m2=p->member[jj];
   if(m2==-1)
    continue;
   if(m1==m2)
    continue;
   if(g->adj[m1][m2])
    continue;
  if(p->ignore[m2]==true)
	  continue;
   //check
   if((p->origrid[ii][0]-p->origrid[jj][0])*(p->origrid[ii][0]-p->origrid[jj][0])>1)
    continue;
   if((p->origrid[ii][1]-p->origrid[jj][1])*(p->origrid[ii][1]-p->origrid[jj][1])>1)
    continue;
   if((p->origrid[ii][2]-p->origrid[jj][2])*(p->origrid[ii][2]-p->origrid[jj][2])>1)
    continue;
   
   g->adj[m1][m2]=true;
   g->adj[m2][m1]=true;
   Ne++;
   //printf("%d m1:%d %d m2:%d\n",ii,m1,jj,m2);
   //Update with symmetry data
   	for(int sym=0;sym<s->Nsym-1;sym++){
	 int id1=s->tbl[m1][sym];
	 int id2=s->tbl[m2][sym];

	 	if(id1==-1 || id2==-1)
		 continue;

	 g->adj[id1][id2]=true;
  	 g->adj[id2][id1]=true;
  	 Ne++;
	
	}
  }
 }
 
 //return true;
 //Ne=p->Ncd*p->Ncd*0.5;
 printf("#Fin checking connect Ne= %d\n",Ne);
 if((g->edge=(EDGE *)malloc(sizeof(EDGE)*Ne))==NULL)
  return true;

 int Nc=0;
 for(i=0;i<p->Ncd;i++)
  g->cid[i]=i;

 //copy
 for(i=0;i<Ntabu*2;i++)
  tmp_tb[i]=tabu[i];

 Ne=0;
 for(i=0;i<p->Ncd;i++){
  double d=0;
  double tmp[3];

  if(p->ignore[i]==true)
	  continue;
  //g->cid[i]=i;

  //for(j=i+1;j<g->Nnode;j++){
  for(j=i+1;j<p->Ncd;j++){
   if(p->ignore[j]==true)
	  continue;
   if(g->adj[i][j]==false)
    continue;
   tmp[0]=p->cd[i][0]-p->cd[j][0];
   tmp[1]=p->cd[i][1]-p->cd[j][1];
   tmp[2]=p->cd[i][2]-p->cd[j][2];

   d=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
   //if(d>d2cut)
//	   continue;
   g->edge[Ne].d=d;//map based
   g->edge[Ne].id1=i;
   g->edge[Ne].id2=j;
   //printf("*%d %f %d %d\n",Ne,g->edge[Ne].d,g->edge[Ne].id1,g->edge[Ne].id2);
   Ne++;
  }
 }
 printf("#Nnode= %d Ntabu= %d Ne= %d\n",g->Nnode,Ntabu,Ne);

 //very slow but....
 //edge table for sym
 if((g->etbl=(int **)malloc(sizeof(int *)*Ne))==NULL)
  return true;
 for(i=0;i<Ne;i++)
  if((g->etbl[i]=(int *)malloc(sizeof(int)*(Nsym-1)))==NULL)
	  return true;

 //Search...slow
 for(int i=0;i<Ne;i++){
	 int id11,id12,id21,id22;
	 int id1=g->edge[i].id1;
	 int id2=g->edge[i].id2;
	for(int sym=0;sym<Nsym-1;sym++){
	 int find=-1;
	 id11=s->tbl[id1][sym];
         id12=s->tbl[id2][sym];
	 if(id11==-1||id12==-1){
	  g->etbl[i][sym]=-1;
	  continue;
	 }
	 if(id11==id12){
	  g->etbl[i][sym]=-1;
	  continue;
	 }

	 	for(int j=0;j<Ne;j++){
			id21=g->edge[j].id1;
			id22=g->edge[j].id2;
			if((id21==id11 && id22==id12)||(id21==id12 && id22==id11)){
			 find=j;
				break;
			}
		
		}
		if(find==-1){
		 printf("#No SYM Edge:etbl %d ori: %d %d sym %d %d:%d\n",i,id1,id2,sym,id11,id12);
		 //return true;
		}
		g->etbl[i][sym]=find;
	}
 
 }
 
 

 //New
 g->Ne=Ne;
 //sort
 qsort(g->edge,Ne,sizeof(EDGE),cmp_edge_d);

 //MST for sym
 int v1,v2,tmp_cid;
 int Nt=0;
 EDGE **tree;
 if((tree=(EDGE **)malloc(sizeof(EDGE *)*Ne))==NULL)
  return true;

 int MaxCid=0;
 for(i=0;i<Ne;i++)
  g->edge[i].mst=false;
 for(i=0;i<Ne;i++){
	 if(g->edge[i].mst==true)
		 continue;
 	
  v1=g->edge[i].id1;
  v2=g->edge[i].id2;

  g->edge[i].mst=false;

  if(g->cid[v1]==g->cid[v2])
   continue;

  //check tabu list
  bool flag=false;
  for(int t=0;t<Ntabu;t++){
   if(tmp_tb[t*2]==g->cid[v1] && tmp_tb[t*2+1]==g->cid[v2]){
    flag=true;
    break;
   }
   if(tmp_tb[t*2]==g->cid[v2] && tmp_tb[t*2+1]==g->cid[v1]){
    flag=true;
    break;
   }
  }
  if(flag==true)
   continue;

  /*
  //check tabu for sym
  for(int sym=0;sym<Nsym-1;sym++){
   int I=g->etbl[i][sym];
   int V1=g->edge[I].id1;
   int V2=g->edge[I].id2;
   if(g->edge[I].mst==true){
	   printf("#Already Used %d\n",I);
    flag=true;
   }
   	//????????
	if(g->cid[V1]==g->cid[V2]){
   	 flag=true;
	 //printf("##cid %d %d = %d %d\n",V1,g->cid[V1],V2,g->cid[V2]);
	 break;
     	}
	//--------------
	for(int t=0;t<Ntabu;t++){
   	 if(tmp_tb[t*2]==g->cid[V1] && tmp_tb[t*2+1]==g->cid[V2]){
   	  flag=true;
   	  break;
   	 }
   	 if(tmp_tb[t*2]==g->cid[V2] && tmp_tb[t*2+1]==g->cid[V1]){
   	  flag=true;
   	  break;
   	 }
  	}
	if(flag==true)
		break;
  }
  if(flag==true)
   continue;
*/

  //NEW
  g->edge[i].mst=true;

  tree[Nt]=&(g->edge[i]);
  tmp_cid=g->cid[v2];
  Nt++;
  if(MaxCid<tmp_cid)
   MaxCid=tmp_cid;

  //update cid in the tree
  for(j=0;j<Nt;j++){
   if(g->cid[tree[j]->id1]==tmp_cid)
    g->cid[tree[j]->id1]=g->cid[v1];
   if(g->cid[tree[j]->id2]==tmp_cid)
    g->cid[tree[j]->id2]=g->cid[v1];
  }

  //update tmp_tb
  for(j=0;j<Ntabu;j++){
   if(tmp_tb[j*2]==tmp_cid)
    tmp_tb[j*2]=g->cid[v1];
   if(tmp_tb[j*2+1]==tmp_cid)
    tmp_tb[j*2+1]=g->cid[v1];
  }

  //Update Sym
  /*
  for(int sym=0;sym<Nsym-1;sym++){
   int I=g->etbl[i][sym];
   g->edge[I].mst=true;
   tree[Nt]=&(g->edge[I]);
   v1=g->edge[I].id1;
   v2=g->edge[I].id2;
   tmp_cid=g->cid[v2];
   Nt++;
   if(Nt>Ne)
	   printf("%d sym:%d %d/%d\n",i,sym,Nt,Ne);
   
   if(MaxCid<tmp_cid)
    MaxCid=tmp_cid;

   //update cid in the tree
   for(j=0;j<Nt;j++){
    if(g->cid[tree[j]->id1]==tmp_cid)
     g->cid[tree[j]->id1]=g->cid[v1];
    if(g->cid[tree[j]->id2]==tmp_cid)
     g->cid[tree[j]->id2]=g->cid[v1];
   }

   //update tmp_tb
   for(j=0;j<Ntabu;j++){
    if(tmp_tb[j*2]==tmp_cid)
     tmp_tb[j*2]=g->cid[v1];
    if(tmp_tb[j*2+1]==tmp_cid)
     tmp_tb[j*2+1]=g->cid[v1];
   }
  }
  */
  //printf("Nt=%d\n",Nt);
 }
 g->tree=tree;
 g->Nt=Nt;
 printf("#Nt= %d / Ne= %d\n",Nt,Ne);

//**************


 //Remove isolated tree
 printf("#MaxCid= %d\n",MaxCid);
 int *Ncid;
 //if((Ncid=(int *)calloc(sizeof(int),(g->Nnode)))==NULL)
 if((Ncid=(int *)calloc(sizeof(int),(p->Ncd)))==NULL)
  return true;
 //for(i=0;i<g->Nnode;i++)
 for(i=0;i<p->Ncd;i++)
  Ncid[g->cid[i]]++;
 int UseCid=-1;
 int Nuse=0;
 for(i=0;i<=MaxCid;i++){
  if(Nuse<Ncid[i]){
   UseCid=i;
   Nuse=Ncid[i];
  }
 }
 printf("#UseCid= %d N= %d/%d\n",UseCid,Nuse,g->Nnode);
 if(UseCid==-1)
  return true;
 
 g->Nchain=0;
 SIMPLE_TBL simtbl[100];
 for(i=0;i<=MaxCid;i++){
  if(Ncid[i] > 20){//Ignore too small fragments
   printf("##Chain%d Ncd= %d / %d %f\n",i,Ncid[i],g->Nnode,(double)Ncid[i]/(double)g->Nnode);

   simtbl[g->Nchain].cid=i;
   simtbl[g->Nchain].Nmem=Ncid[i];;

   //g->chain[g->Nchain]=i;
   //g->Nchain_mem[g->Nchain]=Ncid[i];
   g->Nchain++;
  }
 }
 qsort(simtbl,g->Nchain,sizeof(SIMPLE_TBL),cmp_simple);
 //Get only top Nsym
 for(i=0;i<g->Nchain;i++){
  g->chain[i]=simtbl[i].cid;
  g->Nchain_mem[i]=simtbl[i].Nmem;
  printf("##SORT Chain%d Ncd= %d / %d %f\n",i,g->Nchain_mem[i],g->Nnode,(double)g->Nchain_mem[i]/(double)g->Nnode);
 }

//*******************

 //return 0;
/*
 //clean edge shift
 int Ntmp=0;
 for(i=0;i<Ne;i++){
  if(UseCid==g->cid[g->edge[i].id1]){
   g->edge[Ntmp]=g->edge[i];
   Ntmp++;
  }
 }
 g->Ne=Ntmp;
 //sort
 qsort(g->edge,g->Ne,sizeof(EDGE),cmp_edge_d);
 
 //
 Nt=0;
 //init cid
 for(i=0;i<g->Nnode;i++)
  g->cid[i]=i;
 for(i=0;i<g->Ne;i++){
  v1=g->edge[i].id1;
  v2=g->edge[i].id2;
  g->edge[i].mst=false;
  g->edge[i].local=false;
  //assigin edge-id
  g->edge[i].eid=i;//!!!!!!
  if(g->cid[v1]==g->cid[v2])
   continue;

  tree[Nt]=&(g->edge[i]);
  tmp_cid=g->cid[v2];
  g->edge[i].mst=true;//Used in MST
  Nt++;
  if(MaxCid<tmp_cid)
   MaxCid=tmp_cid;
  //update cid in the tree
  for(j=0;j<Nt;j++){
   if(g->cid[tree[j]->id1]==tmp_cid)
    g->cid[tree[j]->id1]=g->cid[v1];
   if(g->cid[tree[j]->id2]==tmp_cid)
    g->cid[tree[j]->id2]=g->cid[v1];
  }
 }
 //g->tree=tree;
 g->Nt=Nt;
 
 printf("#After cleaning.. Nt= %d Ne= %d\n",Nt,Ntmp);
*/

 //Set Edge dens
 #pragma omp parallel for schedule(dynamic,5)
 for(int ii=0;ii<g->Ne;ii++){
  double cd1[3],vec[3],dens,MinDens;
  int v1=g->edge[ii].id1;
  int v2=g->edge[ii].id2;
 
  //v1->v2
  vec[0]=p->cd[v2][0]-p->cd[v1][0];
  vec[1]=p->cd[v2][1]-p->cd[v1][1];
  vec[2]=p->cd[v2][2]-p->cd[v1][2];
  
  MinDens=99999999;
  for(int jj=1;jj<11;jj++){
   cd1[0]=p->cd[v1][0]+vec[0]*0.1*(double)(jj);
   cd1[1]=p->cd[v1][1]+vec[1]*0.1*(double)(jj);
   cd1[2]=p->cd[v1][2]+vec[2]*0.1*(double)(jj);
   dens=meanshift_pos(mrc,cd1);
   if(dens<MinDens)
    MinDens=dens;
  }
  //printf("#dens %d %f / %f v1:%d v2:%d\n",ii,dens,0.5*(p->dens[v1]+p->dens[v2]),v1,v2);
  g->edge[ii].dens=MinDens*g->edge[ii].d;
  //printf("#MinDens %d = %f\n",ii,g->edge[ii].dens);
 }
 puts("#FIN Set Edge dens data");


 //Local MST
 #pragma omp parallel for schedule(dynamic,5)
 for(int ii=0;ii<g->Nnode;ii++){
  double vec[3],d2;
  int tmpid;
  int cid[32000];//max 30k coordinates 31733 < emd-6555
  //init cid
  for(int jj=0;jj<g->Nnode;jj++)
   cid[jj]=jj;

  	for(int jj=0;jj<g->Ne;jj++){
  	 int v1=g->edge[jj].id1;
  	 int v2=g->edge[jj].id2;
  	 if(cid[v1]==cid[v2])
  	  continue;

	 //dist
	 vec[0]=p->cd[ii][0]-p->cd[v1][0];
 	 vec[1]=p->cd[ii][1]-p->cd[v1][1];
 	 vec[2]=p->cd[ii][2]-p->cd[v1][2];
	 d2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
	 if(d2>d2cut)
	  continue;

	 vec[0]=p->cd[ii][0]-p->cd[v2][0];
 	 vec[1]=p->cd[ii][1]-p->cd[v2][1];
 	 vec[2]=p->cd[ii][2]-p->cd[v2][2];
	 d2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
	 if(d2>d2cut)
	  continue;

	 g->edge[jj].local=true;

	 tmpid=cid[v2];
	 for(int kk=0;kk<g->Nnode;kk++){
	  if(cid[kk]==tmpid)//update cid
	   cid[kk]=cid[v1];
	 }
  	}
 }
 //End Local MST
 puts("#End Local MST");
 //Set Keep Flag
 double dkeep=cmd.Dkeep/mrc->widthx;
 for(int jj=0;jj<g->Ne;jj++){
  g->edge[jj].keep=false;

  if(g->edge[jj].mst && g->edge[jj].d<dkeep)
   g->edge[jj].keep=true;
 }


 //init
 for(int i=0;i<g->Nnode;i++)
  g->node[i].N=0;
 //input GRAPH
 for(int i=0;i<g->Ne;i++){
  if(g->edge[i].local==false && g->edge[i].mst==false)
   continue;
  int id=g->edge[i].id1;

  //chain!!
  if(g->cid[g->edge[i].id1]!=g->cid[g->edge[i].id2] && g->edge[i].mst==true){
  	printf("Edge %d !! %d != %d\n",i,g->cid[g->edge[i].id1],g->cid[g->edge[i].id2]);
	  continue;
  }

  g->node[id].e[g->node[id].N]=&(g->edge[i]);
  g->node[id].N++;

  //if(g->node[id].N>10)
  //printf("id %d : N= %d\n",id,g->node[id].N);
  id=g->edge[i].id2;
  g->node[id].e[g->node[id].N]=&(g->edge[i]);
  g->node[id].N++;
  if(g->node[id].N>20)
   printf("#**id %d : N= %d\n",id,g->node[id].N);
 }


 //input MST
 puts("#Input MST**");
 mst->Nnode=g->Nnode;
 mst->Ne=g->Nt; //#Edges in the MSTree
 mst->Etotal=g->Ne;//All Edges

 if(InitTree(mst,true)){
  printf("##WARRNING!! MALLOC ERROR!!\n");
  return true;
 }

 //return true;

 puts("#Input MST*1*");
 for(int i=0;i<g->Ne;i++){
  g->edge[i].eid=i;
  if(g->edge[i].mst==false)
   continue;
  int id=g->edge[i].id1;
  //int eid=g->edge[i].eid;
  int eid=g->edge[i].eid;

  mst->len+=g->edge[i].d;
  mst->St=id;//Starting Point

  mst->ActE[eid]=true;//Active Edge = MST
 }




 puts("#End SetUp");
 return false;
}

//Simple Tabu search for Symmetry
bool TabuSym(GRAPH *g,TREE *init_t,TREE *Tr){
 int Nround=cmd.Nround;
 int Nnb=cmd.Nnb;
 int Ntabu=cmd.Ntabu;
 int Nsim=cmd.Nsim;
 double AllowLen=init_t->len*cmd.Allow;
 MOVE *move;
 int tabu[1000];
 int Ntb=0;
 int tbp=0;
 int Nsym=g->Nchain;
 TREE Tnow,Tbest,Tnb[500];

 move=(MOVE *)malloc(sizeof(MOVE)*Nnb);

 printf("#MST Len= %f\n",init_t->len);


 double InitSco;
 TREE t[100];
 //Neighbors
 for(int i=0;i<Nnb;i++){
  CopyTree(init_t,&Tnb[i],true);//copy& malloc
 }
 puts("#Fin Malloc for Tnb");
 //Results
 for(int i=0;i<Nsim*Nsym;i++)
  CopyTree(init_t,&Tr[i],true);
 //CopyTree(init_t,&Tr[0],true);
 puts("#Fin Malloc for Tr");

 CopyTree(init_t,&Tnow,true);
 puts("#Fin Malloc for Tnow");
 CopyTree(init_t,&Tbest,true);
 puts("#Fin Malloc for Tbest");


 for(int sym=0;sym<Nsym;sym++){
  CopyTree(init_t,&Tnow,false);
  CopyTree(init_t,&Tbest,false);
  int cid=g->chain[sym];//chain ID
  double score=QualityTreeChain(g,&Tnow,cid);//score only for cid
  Tnow.score=score;
  Tbest.score=score;
  printf("#chain:%d Initscore= %f\n",sym,score);
 
 	for(int i=0;i<Nsim;i++){
	 int mod_id=i+Nsim*sym;
 	 CopyTree(init_t,&Tnow,false);

  	 tbp=0;
  	 Ntb=0;
  	 double pre_best=0;

	 printf(" score: %f\n",move[0].score);
  		for(int n=0;n<Nround;n++){
   		 //****SetCutTbl(g,&Tnow,tabu,Ntb,AllowLen);
   		 SetCutTblChain(g,&Tnow,tabu,Ntb,AllowLen,cid);
	 	 if(Tnow.Nmv==0){
	 	  printf("#Reset tabu Nmv= %d\n",Tnow.Nmv);
	 	  //reset
	 	  Ntb=0;
	 	  tbp=0;
	 	  //n--;
	 	  continue;
	 	 }

		puts("START Nb search");
	 		#pragma omp parallel for schedule(dynamic,5)
	 		for(int nb=0;nb<Nnb;nb++){
 	 		 if(nb>=Tnow.Nmv){
	 		  move[nb].score=0;
	 		  Tnb[nb].score=0;
	 		  continue;
	 		 }

	  		 CopyTree(&Tnow,&Tnb[nb],false);

	  		 MoveTree(g,&Tnb[nb],Tnow.mv[nb].cut_id,Tnow.mv[nb].add_id);//Cut and Add

	  		 move[nb].score=QualityTreeChain(g,&Tnb[nb],cid);

	  		 //New!!
	  		 move[nb].cut_id=Tnow.mv[nb].cut_id;
	  		 move[nb].add_id=Tnow.mv[nb].add_id;

	  		 Tnb[nb].score=move[nb].score;

				//printf("nb= %d insco= %f\n",nb,move[nb].score);
	  		 double ken;
	  		 if(i>0){
	  		  for(int j=0;j<i;j++){
				  //int mod_id=i+Nsim*sym;
			   int mod_id2=j+Nsim*sym;
			   //printf("mod_id2= %d / %d nb= %d/%d\n",mod_id2,Nsim*Nsym,nb,Nnb);
			   //printf("Tr[mod_id2].Lpath= %d\n",Tr[mod_id2].Lpath);
	  		   //ken=Kendall(Tnb[nb].Path,Tnb[nb].Lpath,Tr[j].Path,Tr[j].Lpath);
	  		   ken=Kendall(Tnb[nb].Path,Tnb[nb].Lpath,Tr[mod_id2].Path,Tr[mod_id2].Lpath);
	  		   //ken=Kendall(Tr[mod_id2].Path,Tr[mod_id2].Lpath,Tr[mod_id2].Path,Tr[mod_id2].Lpath);
	  		   //ken=Kendall(Tnb[nb].Path,Tnb[nb].Lpath,Tnb[nb].Path,Tnb[nb].Lpath);
	    		   //printf("Kendall= %f\n",ken);
	    		 	if(fabs(ken)>0.99){
	    		 	 move[nb].score=0;
	    		 	 Tnb[nb].score=0;
	    		 	 break;
	    	 		}
	   	 	  }
	  	 	 }
			}
		puts("FIN Nb search");
	 	//check best movement
	 	double best_sco=0;
	 	int best_mv=0;
	 		for(int nb=0;nb<Nnb;nb++){
	 		 if(best_sco<move[nb].score && move[nb].score !=pre_best ){
	 		  best_sco=move[nb].score;
	 		  best_mv=nb;
	 		 }
	 		}
		puts("FIN best search");
	 	 //Update Sim Best
	 	 //if(n==0||Tr[i].score < Tnb[best_mv].score){
	 	 if(n==0||Tr[mod_id].score < Tnb[best_mv].score){
			 //printf("bestmv= %d %d/%d\n",best_mv,mod_id,Nsym*Nsim);
	 	  CopyTree(&Tnb[best_mv],&Tr[mod_id],false);
	 	 }
		puts("FIN copy results");
	 	 //Update Best
	 	 if(best_sco > Tbest.score){
	 	  CopyTree(&Tnb[best_mv],&Tbest,false);
	 	  Tbest.score=best_sco;
	 	 }
		//puts("FIN copy best");
	 	 printf("#ch%d Sim: %3d/%3d Round= %4d/%4d Nmv= %6d RoundBestSco= %.1f TotalBest= %.1f\n",sym,i+1,Nsim,n+1,Nround,Tnow.Nmv,best_sco,Tbest.score);
	 	 //Update Now
	 	 CopyTree(&Tnb[best_mv],&Tnow,false);

	 
	 	 pre_best=best_sco;

	 	 //Update Tabu List
	 	 tabu[tbp]=move[best_mv].cut_id;
	 	 if(Ntb<Ntabu) Ntb++;

	 	 tbp++;
	 	 if(tbp>=Ntb) tbp=0;

	 	 tabu[tbp]=move[best_mv].add_id;
	 	 if(Ntb<Ntabu) Ntb++;
	 	 tbp++;
	 	 if(tbp>=Ntb) tbp=0;
	 	 //printf("tbp=%d Ntb=%d\n",tbp,Ntb);
  		}//END OF ROUND
	}//END OF SIM
 }//END OF SYM


 return false;
}

double QualityTreeChain(GRAPH *g,TREE *t,int cid){
 int i,j;

 //DFS
 int st=t->St;
 int Nst=1;
 int v,w,maxi,maxi2;
 double maxd=0;
 double sco=0;
 int pre_st=-11;
 //find cid
 for(int i=0;i<g->Nnode;i++){
	 if(g->cid[i]==cid){
		 st=i;
		 break;
	 }
 }
 //printf("####St= %d cid= %d\n",st,cid);
 
 while(1){
  t->stock[0]=st;
  Nst=1;
  t->nextv[st]=-1;
  //init ActN
  //printf("Total Node= %d\n",t->Nnode);
  for(i=0;i<t->Nnode;i++){
   t->cost[i]=0.00000;
   t->ActN[i]=false;
  }
  maxd=0;
  maxi=-1;
        while(1){
         if(Nst==0)
          break;
         v=t->stock[Nst-1];//pop
         Nst--;
         if(t->ActN[v])
          continue;
         //printf("#pop %d %d/%d Act? %d\n",Nst,v,t->Nnode,t->ActN[v]);
         t->ActN[v]=true;
         //branch
         //printf("Nbra= %d\n",g->node[v].N);
         for(j=0;j<g->node[v].N;j++){
          int eid=g->node[v].e[j]->eid;
          //printf("check eid %d %d\n",eid,t->ActE[eid]);
          if(t->ActE[eid]==false) continue;

          w=g->node[v].e[j]->id1;
          //printf("w1= %d or %d\n",w,g->node[v].e[j]->id2);
          if(t->ActN[w]) continue;
          if(g->cid[w]!=cid) continue;
          //push
          t->stock[Nst]=w;
          Nst++;

          t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
          t->nextv[w]=v;//path

          //printf("next1 %d -> %d or %d\n",v,w,g->node[v].e[j]->id2);
          if(maxd<t->cost[w]){
           maxd=t->cost[w];
           maxi=w;
          }
         }
         for(j=0;j<g->node[v].N;j++){
          int eid=g->node[v].e[j]->eid;
          //printf("eid: %d %d\n",eid,t->ActE[eid]);
          if(t->ActE[eid]==false) continue;

          w=g->node[v].e[j]->id2;
          //printf("w2= %d or %d\n",w,g->node[v].e[j]->id1);
          if(t->ActN[w]) continue;

	  //check chain!!
	  if(g->cid[w]!=cid) continue;
          //push
          t->stock[Nst]=w;
          Nst++;

          t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
          t->nextv[w]=v;//path
          //printf("next2 %d -> %d\n",v,w);

          if(maxd<t->cost[w]){
           maxd=t->cost[w];
           maxi=w;
          }
 }

        }
  //printf("#Maxi=  %d -> %d pre= %d dist=%f\n",st,maxi,pre_st,maxd);
  if(maxi==pre_st){
   break;
  }
  pre_st=st;
  st=maxi;
 }
 if(maxi==-1)
  return 0;

 //Set Start&End
 t->St=maxi;
 t->Ed=st;//!!! Aug-24

 sco=maxd*maxd;

 //Set 1st path
 for(i=0;i<t->Nnode;i++)
  t->ActN[i]=false;

 t->ActN[maxi]=true;
 int now=maxi;
 t->Path[0]=now;
 t->Lpath=1;
 while(1){
  now=t->nextv[now];
  t->ActN[now]=true;

  t->Path[t->Lpath]=now;
  t->Lpath++;

  //printf("#set true %d\n",now);
  if(now==st)
   break;
 }
int besti;
 double bestc,diffc;

 //for(int b=2;b<t->Nnode;b++){
 for(int b=2;b<100;b++){
  besti=-1;
  bestc=0;

        for(int i=0;i<t->Nnode;i++){
         if(t->ActN[i])
          continue;

        ////////
         if(g->node[i].N==0)
          continue;
        ///////
	//chain id
	 if(g->cid[i]!=cid)
		 continue;

         //trace back
         now=i;
         while(1){
          //printf("now= %d -> %d %d\n",now,t->nextv[now],t->node[now].N);
          now=t->nextv[now];
          if(t->ActN[now])
           break;
         }
         diffc=t->cost[i]-t->cost[now];
         if(bestc<diffc){
          bestc=diffc;
          besti=i;
         }
        }
  if(besti==-1)
   break;

  sco+=bestc*bestc;
  //printf("#%d= %f\n",b,sco);
  //trace back again
  now=besti;
         while(1){
          //printf("now=%d\n",now);
          now=t->nextv[now];
          if(t->ActN[now])
           break;
          t->ActN[now]=true; //!!!!!! New Aug-24
         }

 }
 t->score=sco;
 return sco;
}


bool SetCutTblChain(GRAPH *g,TREE *t,int *tabu,int Ntb,double LenCutoff,int cid){
 int i,j,n=0;
 
 for(i=0;i<g->Ne;i++){

  int id1=g->edge[i].id1;
  int id2=g->edge[i].id2;

  if(g->node[id1].N==1) continue;
  if(g->node[id2].N==1) continue;

  //chain!!
  if(g->cid[id1]!=cid) continue;
  if(g->cid[id2]!=cid) continue;

  //Active & not keep
  if(t->ActE[i] && g->edge[i].keep==false){
   bool flag=false;

   for(j=0;j<Ntb;j++){
    if(tabu[j]==i){
     flag=true;
     break;
    }
   }
   if(flag==true)
    continue;

   double len=t->len-g->edge[i].d;
	//Set Addtbl
	SplitChain(g,t,i);
	
 	for(j=0;j<g->Ne;j++){
	 int p1=g->edge[j].id1;
  	 int p2=g->edge[j].id2;

 	 if(t->ActE[j])
 	  continue;
 	 if(t->cid[g->edge[j].id1]*t->cid[g->edge[j].id2]!=-1)
 	  continue;
	 if(g->edge[j].local==false && g->edge[j].mst==false)
	  continue;

	 //chain!!
	 //Only considering the same initial-tree, MST
	 if(g->cid[p1]!=cid) continue;
  	 if(g->cid[p2]!=cid) continue;

  	//Total length restraints
  	 if(len+g->edge[j].d>LenCutoff)
  	  continue;

	 //tabu list
	 bool flag=false;
  	 for(int k=0;k<Ntb;k++){
  	  if(tabu[k]==j){
  	   flag=true;
  	   break;
  	  }
  	 }
  	 if(flag==true)
   	  continue;
	 
	 t->mv[n].cut_id=i;
	 t->mv[n].add_id=j;
	 n++;
 	}
  }
 }
 t->Nmv=n;
 //printf("Nmv= %d\n",t->Nmv);
 if(n>2)
  ShuffleMv(t->mv,n);

 return false;
}

void ShowPathSym(MRC *m,POINTS *p,GRAPH *g, TREE *t, int n,int Nsym){
 //Show path
 //For density*volume, using K-clustring data
 int i,j,k,now,ind;
 int Natm=1;
 double tmp[3];
 int st,ed;
 double fmax=0;
 int xydim=m->xdim*m->ydim;
 //path data
 char code[61]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789";
 
 for(int sym=0;sym<Nsym;sym++){

 //for(int mid=0;mid<n;mid++){
 for(int mod=0;mod<n;mod++){
  int mid=mod+n*sym;
  st=t[mid].St;
  ed=t[mid].Ed;
  //printf("st= %d ed= %d\n",st,ed);
  fmax=0;

	//Assign Volume and Density Data
	for(i=0;i<t[mid].Lpath;i++)
	 t[mid].cost[i]=0;
//puts("START");
	for(int j=0;j<p->Nori;j++){
  	 //printf("1st= %d ed= %d j= %d len= %d\n",st,ed,j,t[mid].Lpath);

	 int x=p->origrid[j][0];
	 int y=p->origrid[j][1];
	 int z=p->origrid[j][2];


	 ind=xydim*z+m->xdim*y+x;
	 //if(p->mask[j]==0.00||m->dens[ind]==0.00)
	 if(m->dens[ind]==0.00)
	  continue;
  	 //map origin xyz and xwidth based coordinates
  	 tmp[0]=(double)x; tmp[1]=(double)y; tmp[2]=(double)z;

  	 //printf("1st= %d ed= %d j= %d\n",st,ed,j);
	 int minid=-1;
	 double mindis=(10.0*m->widthx)*(10.0*m->widthx);//set 10A
	 double d2;
		//Search closest path point
		for(int rnum=0;rnum<t[mid].Lpath;rnum++){
		 i=t[mid].Path[rnum];
  	 //printf("1st= %d ed= %d i= %d\n",st,ed,i);
		 d2=(p->cd[i][0]-tmp[0])*(p->cd[i][0]-tmp[0])
		   +(p->cd[i][1]-tmp[1])*(p->cd[i][1]-tmp[1])
		   +(p->cd[i][2]-tmp[2])*(p->cd[i][2]-tmp[2]);
  		 if(d2<mindis){
		  minid=rnum;
		  mindis=d2;
		 }
		}

	 if(minid!=-1)
	  t[mid].cost[minid]+=m->dens[ind];
	 //printf("minid= %d d= %f\n",minid,sqrt(d2));
	}
	//}}}

	for(i=0;i<t[mid].Lpath;i++)
	 if(fmax<t[mid].cost[i])
	  fmax=t[mid].cost[i];



  printf("#SCORE: %f LEN=%d fmax=%f sym= %d\n",t[mid].score,t[mid].Lpath,fmax,sym);
  printf("MODEL %d\n",mid+1);



	Natm=1;
	for(int rnum=0;rnum<t[mid].Lpath;rnum++){
	 i=t[mid].Path[rnum];
	 tmp[0]=p->cd[i][0]*m->widthx+m->orgxyz[0];
 	 tmp[1]=p->cd[i][1]*m->widthx+m->orgxyz[1];
 	 tmp[2]=p->cd[i][2]*m->widthx+m->orgxyz[2];
 	 printf("ATOM  %5d  CA  ALA %c%4d    ",Natm,code[sym],Natm);
 	 //printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->dens[i]);
 	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,t[mid].cost[rnum]/fmax);
 	 Natm++;
	}
  printf("TER\nENDMDL\n");
 }
 }
}
