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
#include<CL/cl.h>
#include"anqi.hh"

#define MAX_SRCSIZE 8000

#define FIRST_LEVEL_SIMULATION 1000
#define OTHER_LEVEL_SIMULATION 100
#define ENABLE_PROFILING

#define GPU_WORKITEM 1
#define RESULT_PER_WORKITEM 5

#define BRDBUFFER_SIZE 48

#define PRINT_FIRST_LEVEL 1
#define NOT_PRINT_FIRST_LEVEL 0
#define EXPLORE_PARA 0.1

static const int adj[32][4]={
	{ 1,-1,-1, 4},{ 2,-1, 0, 5},{ 3,-1, 1, 6},{-1,-1, 2, 7},
	{ 5, 0,-1, 8},{ 6, 1, 4, 9},{ 7, 2, 5,10},{-1, 3, 6,11},
	{ 9, 4,-1,12},{10, 5, 8,13},{11, 6, 9,14},{-1, 7,10,15},
	{13, 8,-1,16},{14, 9,12,17},{15,10,13,18},{-1,11,14,19},
	{17,12,-1,20},{18,13,16,21},{19,14,17,22},{-1,15,18,23},
	{21,16,-1,24},{22,17,20,25},{23,18,21,26},{-1,19,22,27},
	{25,20,-1,28},{26,21,24,29},{27,22,25,30},{-1,23,26,31},
	{29,24,-1,-1},{30,25,28,-1},{31,26,29,-1},{-1,27,30,-1}
};


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
   return (double)(W+D/2)/N+(double)c*(double)pow(log(totalPlay)/N, 0.5);
}

double standard_devia(NODE *node){
   int L = node->Depth%2? node->W:node->L;
   int W = node->Depth%2? node->L :node->W;
   int D = node->D;
   int Ni = W+L+D;
   return pow(((double)(W+D/2) - (double)Ni * pow((double)(W+D/2)/Ni, 2))/Ni, 0.5);
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
            fprintf(stderr, "%d->%d: %d/%d/%d, %.2f/%.2f/%.2f <%.3f>\n",
                    curNode->premove.st, curNode->premove.ed,
                    curNode->W, curNode->L, curNode->D,
                    (double)curNode->W/(double)(curNode->W+curNode->L+curNode->D),
                    (double)curNode->L/(double)(curNode->W+curNode->L+curNode->D),
                    (double)curNode->D/(double)(curNode->W+curNode->L+curNode->D),
                    UCB(curNode, c));

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
		for(k = 0; curNode != NULL; k++){
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

/********** play with block parallelism ****************/
NODE* bittuhPlay2(BOARD *brd, double c){
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

    NODE *subroot;
    for(subroot = root->child; subroot != NULL; subroot = subroot->siblg){
        old_totalPlay = totalPlay;
        explore(subroot, OTHER_LEVEL_SIMULATION, &tempWin, &tempLose);
        tempTotal = totalPlay - old_totalPlay;
        subroot->W += tempWin;
        subroot->L += tempLose;
        subroot->D += tempTotal - tempWin - tempLose;
        root->W += tempLose;
        root->L += tempWin;
        root->D += tempTotal - tempWin - tempLose;
    }
	// tree-growing loop
	while(GetTickCount()-Tick < TimeOut){

        for(subroot = root->child; subroot != NULL; subroot = subroot->siblg){
            /********* selection ***********/

            curNode = subroot;
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
            for(k = 0; curNode != NULL; k++){
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
        }
		/**************************************/
	}

	//give out final answer
	curNode = find_best_child(root, c, PRINT_FIRST_LEVEL);
	return curNode;

}


int main() {

    FILE *fp;
    fp = fopen("kernel.c", "r");
    size_t srcsize;

    char *kernelsrc = (char*)malloc(MAX_SRCSIZE*sizeof(char));
    srcsize = fread(kernelsrc, sizeof(char), MAX_SRCSIZE, fp);
    fclose(fp);
    fprintf(stderr, "...%d\n", srcsize);



    //OpenCL setup

    cl_int status;

	/* get platform */
	cl_uint numPlatforms = 0;
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get platform number\n");
	cl_platform_id *platforms = NULL;
	platforms = (cl_platform_id*)malloc(numPlatforms*sizeof(cl_platform_id));
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get platform ID\n");

    /* get device */
	cl_uint numDevices = 0;
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get device number\n");
	cl_device_id *devices;
	devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get device ID\n");

    /* context */
    cl_context context;
	context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &status);
	if(status != CL_SUCCESS) fprintf(stderr, "context create err\n");

    /* command queue */
	cl_command_queue cmdQueue;
	cmdQueue = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &status);
	if(status != CL_SUCCESS) fprintf(stderr, "command queue create err\n");


    /* program */
    char msg[4096];
    size_t len;
    cl_program program = clCreateProgramWithSource(context, 1, (const char**)&kernelsrc, (const size_t*)&srcsize, &status);

	if(status != CL_SUCCESS) fprintf(stderr, "create program error\n");
	status = clBuildProgram(program, numDevices, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS){
        fprintf(stderr, "build error %d\n", status);
        status = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 4096*sizeof(char), msg, &len);
        fprintf(stderr, "%d\n%s\n", status, msg);
    }

    cl_mem bufBoard = clCreateBuffer(context, CL_MEM_READ_ONLY, BRDBUFFER_SIZE*sizeof(int), NULL, &status);
    if(status != CL_SUCCESS) fprintf(stderr, "fuck1!\n");
    cl_mem bufResult = clCreateBuffer(context, CL_MEM_READ_WRITE, (GPU_WORKITEM * RESULT_PER_WORKITEM+200)*sizeof(int), NULL, &status);
    if(status != CL_SUCCESS) fprintf(stderr, "fuck2!\n");
    //cl_mem bufADJ = clCreateBuffer(context, CL_MEM_READ_WRITE, 32*4*sizeof(int), NULL, &status);
    //if(status != CL_SUCCESS) fprintf(stderr, "fuck2_2!\n");


    BOARD A;
    A.LoadGame("board4.txt");
    A.Display();
    int brdarray[48];
    brdarray[0] = A.who;
    for(int i = 0; i < 32; i++){
        brdarray[1+i] = A.fin[i];
    }
    for(int i = 0; i < 14; i++){
        brdarray[33+i] = A.cnt[i];
    }
    brdarray[47] = A.noFight;

    status = clEnqueueWriteBuffer(cmdQueue, bufBoard, CL_TRUE, 0, 48*sizeof(int), brdarray, 0, NULL, NULL);
    if(status != CL_SUCCESS) fprintf(stderr, "fuck3!\n");
    //status = clEnqueueWriteBuffer(cmdQueue, bufADJ, CL_TRUE, 0, 32*4*sizeof(int), adj, 0, NULL, NULL);
    //if(status != CL_SUCCESS) fprintf(stderr, "fuck3_1!\n");


    cl_kernel kernel = clCreateKernel(program, "simulate", &status);
    if(status != CL_SUCCESS) fprintf(stderr, "fuck4! %d\n", status);
    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufBoard);
    if(status != CL_SUCCESS) fprintf(stderr, "fuck5! %d\n", status);
    status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufResult);
    if(status != CL_SUCCESS) fprintf(stderr, "fuck6! %d\n", status);
    //status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufADJ);
    //if(status != CL_SUCCESS) fprintf(stderr, "fuck6! %d\n", status);

    size_t globalWorkSize[1];
	globalWorkSize[0] = GPU_WORKITEM;
	fprintf(stderr, "enqueu work\n");
	for(int i = 0; i < 1; i++){
        status = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
        if(status != CL_SUCCESS) fprintf(stderr, "fuck7! %d\n", status);

	}
	status = clFinish(cmdQueue);


	int result[GPU_WORKITEM * RESULT_PER_WORKITEM+200] = {0};
	clEnqueueReadBuffer(cmdQueue, bufResult, CL_TRUE, 0, (RESULT_PER_WORKITEM * GPU_WORKITEM + 200)*sizeof(int), result, 0, NULL, NULL);
	if(status != CL_SUCCESS) fprintf(stderr, "fuck8!\n");

    int draw = 0;
    int win = 0;
    int lose = 0;
    for(int i = 0; i < GPU_WORKITEM; i++){
        draw += result[RESULT_PER_WORKITEM*i];
        win += result[RESULT_PER_WORKITEM*i+1];
        lose += result[RESULT_PER_WORKITEM*i+2];
        fprintf(stderr, "<%d,", result[RESULT_PER_WORKITEM*i]);
        fprintf(stderr, "%d,", result[RESULT_PER_WORKITEM*i+1]);
        fprintf(stderr, "%d>\n", result[RESULT_PER_WORKITEM*i+2]);
    }
    fprintf(stderr, "%d/%d/%d\n", draw, win, lose);

    fprintf(stderr, "who: %d\n", result[5]);
    for(int i=1; i < 48; i++){
        fprintf(stderr, "%d: %d..\n", i, result[5+i]);
    }
    fprintf(stderr, "last random number: %d\n", result[3]);
    fprintf(stderr,"simulate depth: %d\n", result[53]);
    fprintf(stderr,"still has %d XXXmove...\n", result[54]);
    for(int i = 0; i < result[54]; i++){
        fprintf(stderr, "[%d->%d]\n", result[55+2*i], result[56+2*i]);
    }

	BOARD C;
    C.who = result[5];
    for(int i = 0; i < 32; i++){
            C.fin[i] = FIN(result[6+i]);
    }
    for(int i = 0; i < 14; i++){
        C.cnt[i] = result[38+i];
    }

    C.Display();


/*
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
        bestNode = bittuhPlay2(&B,EXPLORE_PARA);
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
*/	scanf("%d");

	return 0;
}
