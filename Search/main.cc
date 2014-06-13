/*****************************************************************************\
 * Theory of Computer Games: Fall 2012
 * Chinese Dark Chess Search Engine Template by You-cheng Syu
 *
 * This file may not be used out of the class unless asking
 * for permission first.
 \*****************************************************************************/
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<windows.h>
#include"anqi.hh"


#define FIRST_LEVEL_SIMULATION 2000
#define OTHER_LEVEL_SIMULATION 100
#define ENABLE_PROFILING

#define PRINT_FIRST_LEVEL 1
#define NOT_PRINT_FIRST_LEVEL 0
#define EXPLORE_PARA 0.1


#ifdef ENABLE_PROFILING
DWORD ENABLE_PROFILING_Tick;
DWORD simulateTick;
DWORD searchTick;
int height;
int drawCount;
int highestSimulateDepth;
int pruned;
#endif

DWORD Tick;     // 開始時刻
int   TimeOut;  // 時限
double piii = 1;
int totalPlay;

bool TimesUp() {
	return GetTickCount()-Tick>=TimeOut;
}

/********* bittuh's simulation ***************************/
/* 0 if player of current board lose, 1 for win, 2 for draw */
int bittuhSimulate(BOARD *brd) {
    ENABLE_PROFILING_Tick = GetTickCount();
	const CLR player=brd->who;
	int randn;
	MOVLST lst;
	BOARD tempBoard;
	memcpy(&tempBoard, brd, sizeof(BOARD));
    int simulateDepth = 0;

	while(!tempBoard.ChkLose()&&tempBoard.noFight<40){  // randomly play until win/lose/draw
#ifdef ENABLE_PROFILING
		simulateDepth++;
		if(highestSimulateDepth < simulateDepth) highestSimulateDepth = simulateDepth;
#endif

		bittuhEatGen(tempBoard, lst);   //Eat as the highest priority
		if(lst.num == 0){
                bittuhNotEatGen(tempBoard, lst);    //if no Eat, use other kinds
		}

		randn = rand()%lst.num;
		tempBoard.Move(lst.mov[randn]);

	}

#ifdef ENABLE_PROFILING
    if(!tempBoard.ChkLose()){
        drawCount++;
	}
	simulateTick += GetTickCount()- ENABLE_PROFILING_Tick;
#endif // DEBUG


	if(tempBoard.ChkLose()) return tempBoard.who == player ? 0 : 1; //win or lose
	else return 2; //draw

}

/************** bittuh's expansion ********************/
void explore(NODE *parent, const int numSimulation, int *parent_WIN, int *parent_LOSE){ //return number of wins of parent node board
	int i, j;
	NODE *curNode, *newNode;
	MOVLST legal_moves;

	bittuhAllMoveGen(*(parent->posi), legal_moves);

    *parent_WIN = 0;
    *parent_LOSE = 0;
	// for all legal moves, add a new node to the game tree
	for(i = 0; i < legal_moves.num; i++){
		newNode = (NODE*)malloc(sizeof(NODE));
		newNode->posi = (BOARD*)malloc(sizeof(BOARD));
		memcpy(newNode->posi, parent->posi, sizeof(BOARD));
		newNode->parent = parent;
		newNode->child = NULL;
		newNode->siblg = NULL;
		newNode->Depth = parent->Depth + 1;
#ifdef ENABLE_PROFILING
		if(newNode->Depth > height){
                height = newNode->Depth;
		}
#endif
		newNode->premove = legal_moves.mov[i];
		newNode->posi->Move(legal_moves.mov[i]);
		newNode->W = 0;
		newNode->L = 0;
		newNode->D = 0;

		/*************simulate****************/
		for(j = 0; j < numSimulation; j++){
            int simulatResult = bittuhSimulate(newNode->posi);
			if(simulatResult == 0){
				newNode->L++;
				*parent_WIN+=1;
			}
			else if(simulatResult == 1){
                newNode->W++;
                *parent_LOSE+=1;
			}
			else
                newNode->D++;
			totalPlay++;
		}
		/**************************************/
		if(i == 0){
			parent->child = newNode;
			curNode = parent->child;
		}
		else{
			curNode->siblg = newNode;
			curNode = curNode->siblg;
		}
	}
}
/*********************************************************/

double UCB(NODE *node, double c){
   //int W = node->Depth%2 ? node->L:node->W;
   //int L = node->Depth%2 ? node->W:node->L;
   int W = node->L;
   int L = node->W;
   int D = node->D;
   //int W = node->W;
   //int L = node->L;
   int N = W+L+D;
   return (double)W/N+(double)c*(double)pow(log(totalPlay)/N, 0.5);
}

double standard_devia(NODE *node){
   int L = node->Depth%2? node->W:node->L;
   int W = node->Depth%2? node->L:node->W;
   int D = node->D;
   int Ni = W+L+D;
   return pow(((double)W - (double)Ni * pow((double)W/Ni, 2))/Ni, 0.5);
}


/******************* bittuh's selection **********************/
NODE *find_best_child(NODE *parent, double c, int profile_flag){
	int i = 0, k = 0, mm, nn, m, n;
	double LEO, REO;
	double ucb;

	NODE *curNode;
	NODE *preNode, *nextNode;
	double best_score = UCB(parent->child, c);
	NODE *best_child = parent->child;
	for(curNode = parent->child; curNode!= NULL ;curNode = curNode->siblg){
        if(profile_flag == PRINT_FIRST_LEVEL)
            fprintf(stderr, "%d->%d: %d/%d/%d, <%f>\n", curNode->premove.st, curNode->premove.ed, curNode->W, curNode->L, curNode->D, UCB(curNode, c));

		ucb = UCB(curNode, c);
		if(best_score <= ucb){
			best_score = ucb;
			best_child = curNode;
		}
	}


	curNode = parent->child;
	preNode = parent;
	nextNode = curNode->siblg;

	// PRUNE!!!!!!

    if(parent->W + parent->L + parent->D > 10000){
        while(curNode!= NULL){
	       if(curNode->W + curNode->L + curNode->D > 1000 && best_child->W + best_child->L + best_child->D > 1000){
                LEO = (double)best_child->W/(best_child->W + best_child->L + best_child->D) - piii * standard_devia(best_child);
                REO = (double)curNode->W/(curNode->W + curNode->L + curNode->D) + piii * standard_devia(curNode);
                if(LEO > REO){  //prune current node
                    pruned++;
                    if(preNode == parent) preNode->child = curNode->siblg;
                    else preNode->siblg = curNode->siblg;
                    free(curNode->posi);
                    free(curNode);
                    if(preNode == parent) curNode == parent->child;
                    else curNode == preNode->siblg;
                    }
                else{   //keep current node and move on
                    if(preNode == parent) preNode = preNode->child;
                    else preNode = preNode->siblg;
                    curNode = preNode->siblg;
                    }
            }
            else{
                if(preNode == parent) preNode = preNode->child;
                else preNode = preNode->siblg;
                curNode = preNode->siblg;
                }
        }
    }


	return best_child;
}
/*********************************************************/


/*********** bittuh's play **************/
NODE* bittuhPlay(BOARD *brd, double c){
	NODE *root = (NODE*)malloc(sizeof(NODE));
	NODE *curNode = root;
	MOVLST lst;
	int tempWin, tempLose, tempTotal;
	int i, old_totalPlay, k = 0;
	root->parent = NULL;
	root->siblg = NULL;
	root->child = NULL;
	root->posi = brd;
	root->Depth = 0;
	explore(root, FIRST_LEVEL_SIMULATION, &tempWin, &tempLose);
	//if(tempWin == 0 || tempLose == 0) fprintf(stderr, "fuck\n");
	root->W = tempWin;
	root->L = tempLose;
	root->D = totalPlay - tempWin - tempLose;

	// tree-growing loop
	while(GetTickCount()-Tick < TimeOut){
		/********* selection ***********/
		curNode = root;
		ENABLE_PROFILING_Tick = GetTickCount();
		while(curNode->child != NULL){
			curNode = find_best_child(curNode, c, NOT_PRINT_FIRST_LEVEL);
		}
		searchTick += GetTickCount() - ENABLE_PROFILING_Tick;

		old_totalPlay = totalPlay;
		/******* expansion & simulation *******/
		explore(curNode, OTHER_LEVEL_SIMULATION, &tempWin, &tempLose);
		//if(tempWin == 0 || tempLose == 0) fprintf(stderr, "fuck\n");
		//fprintf(stderr, "~~curNode: %d, %d\n", curNode->W, curNode->W+curNode->L+curNode->D);
		/********* back propogation **********/
		tempTotal = totalPlay - old_totalPlay;
		for(k = 0; curNode->parent != NULL; k++){
			if(k%2 == 0){
				curNode->W += tempWin;
				curNode->L += tempLose;
			}
			else{
				curNode->L += tempWin;
				curNode->W += tempLose;
			}
			curNode->D += tempTotal - tempWin - tempLose;
            //fprintf(stderr, "~~curNode: %d, %d\n", curNode->W, curNode->W+curNode->L+curNode->D);
			curNode = curNode->parent;
		}
		/**************************************/
	}

	//give out final answer
	curNode = find_best_child(root, c, PRINT_FIRST_LEVEL);
	return curNode;

}
/**************************************************/
int main() {
	srand(Tick=GetTickCount());

	BOARD B;
	NODE* bestNode;
	MOV mymove;
	TimeOut=(B.LoadGame("board.txt")-3)*1000;
	totalPlay = 0;

#ifdef ENABLE_PROFILING
	simulateTick = 0;
	drawCount = 0;
	height = 0;
	highestSimulateDepth = 0;
	searchTick = 0;
	pruned = 0;
#endif

	B.Display();

	if(!B.ChkLose()){
        bestNode = bittuhPlay(&B,EXPLORE_PARA);
        Output(mymove = bestNode->premove);
	}
#ifdef ENABLE_PROFILING
    B.Move(mymove);
    fprintf(stderr, "=================================\n");
    fprintf(stderr, "play done, total simulate %d times, draw %d.\n", totalPlay, drawCount);
    fprintf(stderr, "tree height:%d\n", height);
    fprintf(stderr, "pruned %d nodes\n", pruned);
    fprintf(stderr, "highest simulate moves in one iteration: %d\n", highestSimulateDepth);
    fprintf(stderr, "spend %ld in %ld ticks to simulate\n", simulateTick, GetTickCount()-Tick);
    fprintf(stderr, "spend %ld in %ld ticks to search\n", searchTick, GetTickCount()-Tick);
    fprintf(stderr,"\n");
    B.Display();
#endif // ENABLE_PROFILING

    if(mymove.st != mymove.ed)fprintf(stderr,"best: %d->%d",mymove.st, mymove.ed);
    else fprintf(stderr,"best: flip %d",mymove.st);
	scanf("%d");

	return 0;
}
