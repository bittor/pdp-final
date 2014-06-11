/*****************************************************************************\
 * Theory of Computer Games: Fall 2012
 * Chinese Dark Chess Search Engine Template by You-cheng Syu
 *
 * This file may not be used out of the class unless asking
 * for permission first.
\*****************************************************************************/
#include<cstdio>
#include<cstdlib>
#include<windows.h>
#include"anqi.hh"

#define SIMULATION_ITERATION 20000
#define SIMULATION_STEP 20

typedef  int SCORE;
static const SCORE INF=1000001;
static const SCORE WIN=1000000;
SCORE SearchMax(const BOARD&,int,int);
SCORE SearchMin(const BOARD&,int,int);

DWORD Tick;     // �}�l�ɨ�
int   TimeOut;  // �ɭ�
MOV   BestMove; // �j�X�Ӫ��̨ε۪k

bool TimesUp() {
	return GetTickCount()-Tick>=TimeOut;
}

// �@�ӭ��q�����誺�f�����
SCORE Eval(const BOARD &B) {
	int cnt[2]={0,0};
	for(POS p=0;p<32;p++){const CLR c=GetColor(B.fin[p]);if(c!=-1)cnt[c]++;}
	for(int i=0;i<14;i++)cnt[GetColor(FIN(i))]+=B.cnt[i];
	return cnt[B.who]-cnt[B.who^1];
}

SCORE RandomSimulation(const BOARD &B){
    BOARD N(B);
    for(int i=0; i<SIMULATION_STEP; i++){
        if(N.ChkLose())
            return i%2==0 ? +WIN: -WIN;

        MOVLST lst;
        if(TimesUp() || N.MoveGen(lst)==0)
            return Eval(N);

        int p = rand()%lst.num;
        N.Move(lst.mov[p]);
    }
    return Eval(N);
}


MOV Play(const BOARD &B) {
	POS p; int c=0;

	// �s�C���H�H��½�l
	if(B.who==-1){p=rand()%32;return MOV(p,p);}

	// �Y�j�X�Ӫ����G�|��{�b�n�N�ηj�X�Ӫ����k
	MOVLST lst;
    SCORE best_score = -INT_MAX;
    B.MoveGen(lst);

    for(int i=0; i<lst.num; i++){
        SCORE current = 0;
        BOARD N(B);
        N.Move(lst.mov[i]);

        for(int j=0; j<SIMULATION_ITERATION; j++){
            current += RandomSimulation(N);
        }

        if(current > best_score)
            BestMove = lst.mov[i];
    }

	// �_�h�H�K½�@�Ӧa�� ���p�ߥi��w�g�S�a��½�F
	for(p=0;p<32;p++)if(B.fin[p]==FIN_X)c++;
	if(c==0)return BestMove;
	if(lst.num!=0 && rand()%2==0) return BestMove;

	c=rand()%c;
	for(p=0;p<32;p++)if(B.fin[p]==FIN_X&&--c<0)break;
	return MOV(p,p);
}

int main() {
	srand(Tick=GetTickCount());

	BOARD B;
	TimeOut=(B.LoadGame("board.txt")-3)*1000;
	if(!B.ChkLose())Output(Play(B));

	return 0;
}
