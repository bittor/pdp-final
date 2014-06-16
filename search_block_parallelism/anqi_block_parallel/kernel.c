#define FINOFFSET 1

typedef int CLR;



enum LVL {
	LVL_K=0, // �ӱN King
	LVL_G=1, // �K�h Guard
	LVL_M=2, // �۶H Minister
	LVL_R=3, // �Ϩ� Rook     // BIG5 �S���H������
	LVL_N=4, // �X�� kNight
	LVL_C=5, // ���� Cannon
	LVL_P=6  // �L�� Pawn
};

enum FIN {
	FIN_K= 0 /* 'K' �� */ , FIN_k= 7 /* 'k' �N */ , FIN_X=14 /* 'X' ��½ */ ,
	FIN_G= 1 /* 'G' �K */ , FIN_g= 8 /* 'g' �h */ , FIN_E=15 /* '-' �Ů� */ ,
	FIN_M= 2 /* 'M' �� */ , FIN_m= 9 /* 'm' �H */ ,
	FIN_R= 3 /* 'R' �� */ , FIN_r=10 /* 'r' �� */ ,
	FIN_N= 4 /* 'N' �X */ , FIN_n=11 /* 'n' �� */ ,
	FIN_C= 5 /* 'C' �� */ , FIN_c=12 /* 'c' �� */ ,
	FIN_P= 6 /* 'P' �L */ , FIN_p=13 /* 'p' �� */
};

CLR GetColor(FIN f) {
	return f<FIN_X?f/7:-1;
}

LVL GetLevel(FIN f) {
	assert(f<14);
	return LVL(f%7);
}

int ChkEats(FIN fa,FIN fb) {
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

int bittuhEatGen(__local int *brd, __local int *movelst, __constant int *adjacent){
	if(brd[0] == -1) return 0;
    lst[0] = 0;

    for(int st = 0; st < 32; st++){
        FIN stFin=brd[FINOFFSET+st];
		if(GetColor(stFin)!=brd.who)continue;  //not my pieces
        LVL stLevel=GetLevel(stFin);

        if(stLevel != LVL_C){
            for(int z=0;z<4;z++) {          //neighbor
                int ed=adjacent[st*4+z];
                if(ed == -1) continue;
                FIN edFin = brd[FINOFFSET+ed];
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
                for(int ed = adjacent[st*4+z]; ed != -1 && jump < 2; ed = adjacent[ed*4+z]) {
                    FIN edFin = brd[FINOFFSET+ed];
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
    
}

int bittuhMoveGen(__local int *brd, __local int *movelst, __constant int *adjacent){
    if(brd[0]==-1)return false;
	movelst[0]=0;
	for(POS st=0;st<32;st++){
        const FIN stFin=brd[FINOFFSET+st];
		if(GetColor(stFin)!=brd[0])continue;  //not my pieces
		for(int z=0;z<4;z++) {          //neighbor
            const POS ed=adjacent[st*4+z];
            if(ed == -1)continue;
            if(brd[FINOFFSET+ed] == FIN_E){
                movelst[2*movelst[0]+1] = st;
				movelst[2*movelst[0]+2] = ed;
				movelst[0]++;
			}
        }
	}
}

int bittuhFlipGen(__local int *brd, __local int *movelst){
    if(brd[0]==-1)return false;
	movelst[0]=0;
	for(POS p=0;p<32;p++){
        FIN pf = brd[FINOFFSET+p];
        if(pf == FIN_X){
			movelst[2*movelst[0]+1] = st;
			movelst[2*movelst[0]+2] = ed;
			movelst[0]++;
		}

	}

}


