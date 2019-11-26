
typedef struct{
	double mtx[100][3][3];
	double t[100][3];
	int Nsym;
	int **tbl;
} SYM;


typedef struct{
	float **cd;
	int Ncd,Nch;
	int *cid;
	char *chain;
} CIF;

typedef struct{
	int Nmem;
	int cid;
} SIMPLE_TBL;



bool readsym(SYM *,char *);
bool MergePointsSym(MRC *,POINTS *,int);
bool SetUpGraphSym(POINTS *, GRAPH *,MRC *,TREE *,int *,int,SYM *);
bool FindCorrPointsCn(POINTS *,MRC *,int, int *,int *);
bool FindCorrPointsDn(POINTS *,MRC *,int, int *,int *);
bool FindCorrPointsSym(POINTS *,MRC *,SYM *, int *,int *,bool);

bool TabuSym(GRAPH *,TREE *,TREE *);

double QualityTreeChain(GRAPH *,TREE *,int);
bool SetCutTblChain(GRAPH *,TREE *,int *,int,double,int);
void ShowPathSym(MRC *,POINTS *,GRAPH *,TREE *,int,int);

bool readCIF_simple(CIF *,char *);
