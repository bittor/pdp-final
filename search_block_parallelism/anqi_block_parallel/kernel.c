#define FINOFFSET 1
#define COUNTOFFSET 33
#define NOFIGHTOFFSET 47
#define BRD_TOTAL_SIZE 49

#define MAX_EATLIST 120
#define MAX_MOVELIST 120
#define MAX_FLIPLIST 65

#define RESULT_PER_WORKITEM 5

#define EatRate 6
#define MoveRate 2
#define FlipRate 2

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

#pragma OPENCL EXTENSION cl_amd_printf : enable
typedef char CLR;

__constant int adj[32][4]={
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
	if(f <= FIN_P) return 0;
	if(f <= FIN_p) return 1;
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

int myrand(int *seed){
	unsigned int hi, lo;
	hi = 16807 * (*seed >> 16);
	lo = 16807 * (*seed & 0xFFFF);
	lo += (hi & 0x7FFF) << 16;
	lo += hi >> 15;
	if(lo > 2147483647)
		lo -= 2147483647;
	*seed = lo;
	return *seed;
}

char ChkEats(char fa, char fb) {
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

char bittuhEatGen(char *brd, char *movelst, char ADJ[32][4]){
	if(brd[0] == -1) return 0;
	movelst[0] = 0;

	for(char st = 0; st < 32; st++){
		char stFin = brd[FINOFFSET+st];
		if(GetColor(stFin)!=brd[0])continue;  //not my pieces
		LVL stLevel=GetLevel(stFin);

		if(stLevel != LVL_C){
			for(char z=0;z<4;z++) {          //neighbor
				char ed=ADJ[st][z];
				if(ed == -1) continue;
				char edFin = brd[FINOFFSET+ed];
				if(ChkEats(stFin, edFin) && edFin != FIN_E){
					movelst[2*movelst[0]+1] = st;
					movelst[2*movelst[0]+2] = ed;
					movelst[0]++;
				}
			}
		}
		else{
			for(char z=0;z<4;z++) {
				char jump=0;
				for(char ed = ADJ[st][z]; ed != -1 && jump < 2; ed = ADJ[ed][z]) {
					char edFin = brd[FINOFFSET+ed];
					if(edFin == FIN_E || ++jump != 2){
						continue;
					}
					//destination is not Empty or Sealed, and is opponent
					if(edFin < FIN_X && GetColor(edFin) != brd[0] && jump > 0){
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

char bittuhMoveGen( char *brd, char *movelst, char ADJ[32][4]){
	char c = 1;
	if(brd[0]==-1)return 255;
	movelst[0]=0;
	for(char st=0;st<32;st++){
		char stFin=brd[FINOFFSET+st];
		c++;
		if(GetColor(stFin)!=brd[0])continue;  //not my pieces
		c++;
		for(char z=0;z<4;z++) {          //neighbor
			const char ed=ADJ[st][z];
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

char bittuhFlipGen(char *brd, char *movelst){
	if(brd[0]==-1)return false;
	movelst[0]=0;
	for(char p=0;p<32;p++){
		char pFin = brd[FINOFFSET+p];
		if(pFin == FIN_X){
			movelst[2*movelst[0]+1] = p;
			movelst[2*movelst[0]+2] = p;
			movelst[0]++;
		}
	}
	return movelst[0];
}
/*
int myChkLose(int *brd){	//only detect whether any live piece(no matter it's moveable) exist and any unrevealed pieces
	int c = 0;
	if(brd[0]==-1) return -1;

	//if any my piece still unrevealed?
	for(int i=7*brd[0];i<7*brd[0]+7;i++) {
		if(brd[COUNTOFFSET+i]!=0) return 2;
		c++;
	}
	
	// if any live piece is mine?
	for(int p=0;p<32;p++){
	if(GetColor(brd[FINOFFSET+p])==brd[0])
	return 1;
	}

	return 0;
}*/

void Flip(char *brd, char pos, int *seed){
	char f;
	if(brd[FINOFFSET+pos] == FIN_X){
		char i, sum = 0;
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

void Move(char *brd, char st, char ed, int *seed){
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
	const char ADJ[32][4]={
	{ 1,-1,-1, 4},{ 2,-1, 0, 5},{ 3,-1, 1, 6},{-1,-1, 2, 7},
	{ 5, 0,-1, 8},{ 6, 1, 4, 9},{ 7, 2, 5,10},{-1, 3, 6,11},
	{ 9, 4,-1,12},{10, 5, 8,13},{11, 6, 9,14},{-1, 7,10,15},
	{13, 8,-1,16},{14, 9,12,17},{15,10,13,18},{-1,11,14,19},
	{17,12,-1,20},{18,13,16,21},{19,14,17,22},{-1,15,18,23},
	{21,16,-1,24},{22,17,20,25},{23,18,21,26},{-1,19,22,27},
	{25,20,-1,28},{26,21,24,29},{27,22,25,30},{-1,23,26,31},
	{29,24,-1,-1},{30,25,28,-1},{31,26,29,-1},{-1,27,30,-1}};

	int depth = 0;
	int highestdepth = 0;
	const CLR player = brd[0];
	char tempBrd[BRD_TOTAL_SIZE];
	char localBrd[BRD_TOTAL_SIZE];
	char eatlist[MAX_EATLIST];
	char movelist[MAX_MOVELIST];
	char fliplist[MAX_FLIPLIST];
	int idx = get_global_id(0);
	int seed = (int)(idx+1) * brd[48];
	int numSimulate = brd[49];

	for(char i = 0; i < BRD_TOTAL_SIZE; i++){
		localBrd[i] = (char)brd[i];
	}

	int win, lose, draw;
	win = 0;
	lose = 0;
	draw = 0;
	for(int i = 0; i < numSimulate; i++){
		depth = 0;
		//copy board buffer to private address space
		for(char j = 0; j < BRD_TOTAL_SIZE; j++){
			tempBrd[j] = localBrd[j];
		}

		while(tempBrd[NOFIGHTOFFSET] < 40)
		{
			char k = bittuhFlipGen(tempBrd, fliplist);
			char l = bittuhMoveGen(tempBrd, movelist, ADJ);
			char m = bittuhEatGen(tempBrd, eatlist, ADJ);
			if(k+l+m == 0) break;

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

		if(tempBrd[NOFIGHTOFFSET] >= 40)
			draw++; //draw
		else if(player == tempBrd[0])
			win++; //opponent win
		else{
			lose++;	//opponent lose
		}

		

		//result[idx*RESULT_PER_WORKITEM+3] += depth;
		if(highestdepth < depth) highestdepth = depth;
	}
	result[idx*RESULT_PER_WORKITEM] = draw;
	result[idx*RESULT_PER_WORKITEM+1] = win;
	result[idx*RESULT_PER_WORKITEM+2] = lose;
	result[idx*RESULT_PER_WORKITEM+4] = highestdepth;


	/*
	//profiling for only one work item
	for(int i = 0; i < BRD_TOTAL_SIZE; i++){
	result[5+i] = tempBrd[i];
	}
	result[53] = depth;
	result[54] = bittuhEatGen(tempBrd, eatlist);
	for(int i = 0; i < eatlist[0]; i++){
	result[55+2*i] = eatlist[1+2*i];
	result[56+2*i] = eatlist[2+2*i];
	}
	 */




}

