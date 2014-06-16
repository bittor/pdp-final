#define FINOFFSET 1
#define COUNTOFFSET 33
#define NOFIGHTOFFSET 47
#define BRD_TOTAL_SIZE 48

#define MAX_EATLIST 60
#define MAX_MOVELIST 60
#define MAX_FLIPLIST 32

#define EatRate 1
#define MoveRate 2
#define FlipRate 7

#define FIN_K 0
#define	FIN_G 1
#define FIN_M 2
#define FIN_R 3
#define FIN_N 4
#define FIN_C 5
#define FIN_P 6
#define FIN_k 7
#define FIN_g 8
#define FIN_m 9 
#define FIN_r 10
#define FIN_n 11
#define FIN_c 12
#define FIN_p 13
#define FIN_X 14
#define FIN_E 15 

typedef int CLR;

__constant int ADJ[32][4]={
	{ 1,-1,-1, 4},{ 2,-1, 0, 5},{ 3,-1, 1, 6},{-1,-1, 2, 7},
	{ 5, 0,-1, 8},{ 6, 1, 4, 9},{ 7, 2, 5,10},{-1, 3, 6,11},
	{ 9, 4,-1,12},{10, 5, 8,13},{11, 6, 9,14},{-1, 7,10,15},
	{13, 8,-1,16},{14, 9,12,17},{15,10,13,18},{-1,11,14,19},
	{17,12,-1,20},{18,13,16,21},{19,14,17,22},{-1,15,18,23},
	{21,16,-1,24},{22,17,20,25},{23,18,21,26},{-1,19,22,27},
	{25,20,-1,28},{26,21,24,29},{27,22,25,30},{-1,23,26,31},
	{29,24,-1,-1},{30,25,28,-1},{31,26,29,-1},{-1,27,30,-1}
};


typedef enum lvl {
	LVL_K=0,
	LVL_G=1,
	LVL_M=2,
	LVL_R=3,
	LVL_N=4,
	LVL_C=5,
	LVL_P=6
} LVL;




CLR GetColor(int f) {
	if(f<FIN_X) return f/7;
	return -1;
}

LVL GetLevel(int f) {
	switch (f%7){
		case 0: return LVL_K;
		case 1: return LVL_G;
		case 2: return LVL_M;
		case 3: return LVL_R;
		case 4: return LVL_N;
		case 5: return LVL_C;
		case 6: return LVL_P;
	}
	return LVL_P;
}

int myrand(long *seed){
	*seed = (*seed ^ 0x5DEECE66DL) & ((1L << 48) - 1);
	return *seed >> 16;
}

int ChkEats(int fa, int fb) {
	if(fa>=FIN_X)return false;
	if(fb==FIN_X)return false;
	if(fb==FIN_E)return true ;
	if(GetColor(fb)==GetColor(fa))return false;

	const LVL la=GetLevel(fa);
	if(la==LVL_C)return true ;

	const LVL lb=GetLevel(fb);
	if(la==LVL_K)return lb!=LVL_P;
	if(la==LVL_P)return lb==LVL_P||lb==LVL_K;

	return la<=lb;
}

int bittuhEatGen(int *brd, int *movelst){
	if(brd[0] == -1) return 0;
    movelst[0] = 0;

    for(int st = 0; st < 32; st++){
        int stFin = brd[FINOFFSET+st];
		if(GetColor(stFin)!=brd[0])continue;  //not my pieces
        LVL stLevel=GetLevel(stFin);

        if(stLevel != LVL_C){
            for(int z=0;z<4;z++) {          //neighbor
                int ed=ADJ[st][z];
                if(ed == -1) continue;
                int edFin = brd[FINOFFSET+ed];
                if(ChkEats(stFin, edFin) && edFin != FIN_E){
                    movelst[2*movelst[0]+1] = st;
					movelst[2*movelst[0]+2] = ed;
					movelst[0]++;
				}
			}
        }
		else{
            for(int z=0;z<4;z++) {
                int jump=0;
                for(int ed = ADJ[st][z]; ed != -1 && jump < 2; ed = ADJ[ed][z]) {
                    int edFin = brd[FINOFFSET+ed];
                    if(edFin != FIN_E) jump++;
                    //destination is not Empty or Sealed, and is opponent
                    if(edFin < FIN_X && GetColor(edFin) != brd[0]){
						movelst[2*movelst[0]+1] = st;
						movelst[2*movelst[0]+2] = ed;
						movelst[0]++;
						break;
					}
                }
            }
        }
    }
	return movelst[0];
}

int bittuhMoveGen( int *brd, int *movelst){
    if(brd[0]==-1)return false;
	movelst[0]=0;
	for(int st=0;st<32;st++){
        const int stFin=brd[FINOFFSET+st];
		if(GetColor(stFin)!=brd[0])continue;  //not my pieces
		for(int z=0;z<4;z++) {          //neighbor
            const int ed=ADJ[st][z];
            if(ed == -1)continue;
            if(brd[FINOFFSET+ed] == FIN_E){
                movelst[2*movelst[0]+1] = st;
				movelst[2*movelst[0]+2] = ed;
				movelst[0]++;
			}
        }
	}
	return movelst[0];
}

int bittuhFlipGen(int *brd, int *movelst){
    if(brd[0]==-1)return false;
	movelst[0]=0;
	for(int p=0;p<32;p++){
        int pFin = brd[FINOFFSET+p];
        if(pFin == FIN_X){
			movelst[2*movelst[0]+1] = p;
			movelst[2*movelst[0]+2] = p;
			movelst[0]++;
		}
	}
	return movelst[0];
}

int myChkLose(int *brd){	//only detect whether any live piece(no matter it's moveable) exist and any unrevealed pieces
	if(brd[0]==-1) return -1;
	
	//if any my piece still unrevealed?
	for(int i=7*brd[0];i<7*brd[0]+7;i++) {
		if(brd[COUNTOFFSET+i]!=0) return 0;
	}
	// if any live piece is mine?
	for(int p=0;p<32;p++)
		if(GetColor(brd[FINOFFSET+p])==brd[0])
			return 0;
	
	return 1;
}

void Flip(int *brd, int pos, long *seed){
	int f;
	if(brd[FINOFFSET+pos] == FIN_X){
		int i, sum = 0;
		for(i=0;i<14;i++) sum += brd[COUNTOFFSET+i];
		sum = myrand(seed)%sum;
		for(i=0;i<14;i++) if((sum-=brd[COUNTOFFSET+i])<0) break;
		f = i;
	}
	brd[FINOFFSET+pos] = f;
	brd[COUNTOFFSET+f]--;
	if(brd[0] == -1) brd[0]=GetColor(f);
	brd[0]^=1;
}

void Move(int *brd, int st, int ed, long *seed){
	if(st != ed){
		if(brd[FINOFFSET+ed] != FIN_E) brd[NOFIGHTOFFSET] = 0;
		else brd[NOFIGHTOFFSET]++;
		brd[FINOFFSET+ed] = brd[FINOFFSET+st];
		brd[FINOFFSET+st] = FIN_E;
		brd[0] ^= 1;
	}
	else{
		Flip(brd, st, seed);
		brd[NOFIGHTOFFSET] = 0;
	}
}

__kernel
void simulate(const __constant int *brd, __global int *result){
	int depth = 0;
	const CLR player = brd[0];
	int tempBrd[BRD_TOTAL_SIZE];
	int eatlist[MAX_EATLIST];
	int movelist[MAX_MOVELIST];
	int fliplist[MAX_FLIPLIST];
	int idx = get_global_id(0);
	long seed = (long)idx * 0x5a7d9e654;


	for(int i = 0; i < BRD_TOTAL_SIZE; i++){
		tempBrd[i] = brd[i];
		//result[i] = tempBrd[i];
	}


	while(bittuhFlipGen(tempBrd, fliplist) + bittuhMoveGen(tempBrd, movelist) != 0 &&
			myChkLose(tempBrd) != 1 && 
			tempBrd[NOFIGHTOFFSET] < 40 && depth < 1000)
	{
		bittuhEatGen(tempBrd, eatlist); 
		int randn = myrand(&seed)%(EatRate*eatlist[0]+MoveRate*movelist[0]+FlipRate*fliplist[0]);
		if(randn < EatRate*eatlist[0]){
			randn = randn%eatlist[0];
			Move(tempBrd, eatlist[1+2*randn], eatlist[2+2*randn], &seed);
		}
		else if((randn -= EatRate*eatlist[0]) < MoveRate*movelist[0]){
			randn = randn%movelist[0];
			Move(tempBrd, movelist[1+2*randn], movelist[2+2*randn], &seed);
		}
		else{
			randn -= MoveRate*movelist[0];
			randn = randn%fliplist[0];
			Move(tempBrd, fliplist[1+2*randn], fliplist[1+2*randn], &seed);
		}
		depth++;
	}
	
	
	int lose = myChkLose(tempBrd);
	if(lose != 1 && bittuhFlipGen(tempBrd, fliplist) + bittuhMoveGen(tempBrd, movelist) > 0)
		result[idx*3+2]++; //draw
	else if(player == tempBrd[0])
		result[idx*3]++; //opponent win
	else
		result[idx*3+1]++;	//opponent lose
	result[idx*3+3] = depth;
	result[idx*3+4] = tempBrd[NOFIGHTOFFSET];
	for(int i = 0; i < BRD_TOTAL_SIZE; i++){
		result[idx*3+5+i] = tempBrd[i];
	}
	//result[idx*3+53] = 
		
}

