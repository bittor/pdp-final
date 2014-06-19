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
#include<sys/time.h>
#include<omp.h>

#define MAX_SRCSIZE 8000



#define USING_GPU_SIMULATE
#define USING_CPU_SIMULATE

#define USING_BASE_MCTS
#define USING_ROOTPARA_MCTS
#define USING_OUR_MCTS

#define FIRST_LEVEL_SIMULATION 2000
#define OTHER_LEVEL_SIMULATION 500

#define WORKITEM 16
#define WORKITEM_SIMULATE 64
#define FIRST_LVL_GPU_SIMULATE 64
#define OTHER_LVL_GPU_SIMULATE 16
#define MAX_GROUP_NUM 1
#define RESULT_PER_WORKITEM 5

#define NUM_THREADS 4

#define BRDBUFFER_SIZE 53

#define EXPLORE_PARA 0.6

#define ENABLE_PROFILING
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

cl_mem bufResult[NUM_THREADS];
cl_command_queue *cmdQueue;
cl_program program;
cl_context context;
cl_uint numPlatforms;
cl_platform_id *platforms;
cl_uint numDevices;
cl_device_id *devices;
int *queueCounter;


#ifdef ENABLE_PROFILING
DWORD ENABLE_PROFILING_Tick;
DWORD simulateTick;
DWORD searchTick;
int tree_height;
int drawCount;
int highestSimulateDepth;
long totalSimulateDepth;
int pruned;
NODE *root;
#endif

DWORD Tick;     // 開始時刻
int   TimeOut;  // 時限
double piii = 1;
int totalPlay;

bool TimesUp() {
	return GetTickCount()-Tick>=TimeOut;
}
void clCmdProfilingStub(cl_event e){
    cl_int status;
    cl_ulong tSubmit, tStart, tEnd;
    status = clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &tSubmit, NULL);
    status |= clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &tStart, NULL);
    status |= clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &tEnd, NULL);
    if(status != CL_SUCCESS){
        fprintf(stderr, "$:%d profile command execution time error in GPUSimulate...\n", status);
    }
    else{
        fprintf(stderr, "* submit to start: %lu ms\n", (tStart-tSubmit)*1e-6);
        fprintf(stderr, "* start to end   : %lu ms\n", (tEnd-tStart)*1e-6);
    }
}
/********* bittuh's simulation ***************************/

void bittuhSimulate(BOARD *brd, const int numSimulation, int *win, int *lose, int *draw) {

#ifdef ENABLE_PROFILING
    ENABLE_PROFILING_Tick = GetTickCount();
#endif // ENABLE_PROFILING

	const CLR player=brd->who;
	int randn;
	MOVLST lst;
	BOARD tempBoard;
    int simulateDepth;

    *win = 0;
    *lose = 0;
    *draw = 0;

    for(int i = 0; i < numSimulation; i++){
        simulateDepth = 0;
        memcpy(&tempBoard, brd, sizeof(BOARD));
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
        simulateTick += GetTickCount()- ENABLE_PROFILING_Tick;
#endif // ENABLE_PROFILING

        if(!tempBoard.ChkLose()) *draw += 1;
        else if(tempBoard.who == player) *win += 1;
        else *lose += 1;
    }

}


void bittuhGPUSimulate(BOARD *brd, int *win, int *lose, int *draw, int workitem, int numSim){

    //copy brd to GPU buffer
    int tid = omp_get_thread_num();
    char tempArray[BRDBUFFER_SIZE];
    tempArray[0] = brd->who;
    for(int i = 0; i < 32; i++) tempArray[1+i] = brd->fin[i];
    for(int i = 0; i < 14; i++) tempArray[33+i] = brd->cnt[i];
    tempArray[47] = brd->noFight;
    *(int*)(&tempArray[49]) = rand();
    tempArray[48] = numSim;

    cl_int status;
    cl_mem bufBoard = clCreateBuffer(context, CL_MEM_READ_WRITE, BRDBUFFER_SIZE * sizeof(int), NULL, &status);
    //bufResult[tid] = clCreateBuffer(context, CL_MEM_READ_WRITE, (WORKITEM * RESULT_PER_WORKITEM)*sizeof(int), NULL, &status);
    status = clEnqueueWriteBuffer(cmdQueue[1], bufBoard, CL_TRUE, 0, BRDBUFFER_SIZE*sizeof(char), tempArray, 0, NULL, NULL);
    if(status != CL_SUCCESS){
        fprintf(stderr, "$:%d Write buffer in GPUSimulate fail...\n", status);
        exit(1);
    }

    cl_kernel kernel = clCreateKernel(program, "simulate", &status);
    status |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufBoard);
    status |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufResult[tid]);
    if(status != CL_SUCCESS){
        fprintf(stderr, "$:%d Kernel creation error in GPUSimulate...\n", status);
        exit(2);
    }

#ifdef ENABLE_PROFILING
    ENABLE_PROFILING_Tick = GetTickCount();
#endif // ENABLE_PROFILING


    size_t globalWorkSize[1];
    cl_event event;
	globalWorkSize[0] = workitem;
	for(int i = 0; i < MAX_GROUP_NUM; i++){
        status = clEnqueueNDRangeKernel(cmdQueue[1], kernel, 1, NULL, globalWorkSize, NULL, 0, NULL, &event);
        if(status != CL_SUCCESS){
            fprintf(stderr, "$:%d enqueue kernel error in GPUSimulate...\n", status);
            exit(3);
        }

	}
	status = clWaitForEvents(1, &event);
	if(status != CL_SUCCESS){
        fprintf(stderr, "$:%d fail to wait for result in GPUSimulate...\n", status);
        exit(4);
	}

	int *result = (int*)malloc(sizeof(int) * WORKITEM * RESULT_PER_WORKITEM);
	clEnqueueReadBuffer(cmdQueue[1], bufResult[0], CL_TRUE, 0, (RESULT_PER_WORKITEM * WORKITEM )*sizeof(int), result, 0, NULL, &event);
	if(status != CL_SUCCESS){
        fprintf(stderr, "$:%d read result buffer error in GPUSimulate...\n", status);
        exit(5);
	}

#ifdef ENABLE_PROFILING
	simulateTick += GetTickCount()- ENABLE_PROFILING_Tick;
#endif // DEBUG

	*draw = 0;
    *win = 0;
    *lose = 0;
    int depthsum = 0;
    for(int i = 0; i < WORKITEM; i++){
        *draw += result[RESULT_PER_WORKITEM*i];
        *win += result[RESULT_PER_WORKITEM*i+1];
        *lose += result[RESULT_PER_WORKITEM*i+2];
        #ifdef ENABLE_PROFILING
        if(highestSimulateDepth < result[RESULT_PER_WORKITEM*i+4])highestSimulateDepth = result[RESULT_PER_WORKITEM*i+4];
        #endif // ENABLE_PROFILING
    }
    clReleaseMemObject(bufBoard);
    clReleaseKernel(kernel);
}

NODE *createNode(NODE* parent, MOV mv){
        NODE *newNode = (NODE*)malloc(sizeof(NODE));
		newNode->posi = (BOARD*)malloc(sizeof(BOARD));
		memcpy(newNode->posi, parent->posi, sizeof(BOARD));
		newNode->parent = parent;
		newNode->child = NULL;
		newNode->siblg = NULL;
		newNode->Depth = parent->Depth + 1;
#ifdef ENABLE_PROFILING
		if(newNode->Depth > tree_height){
                tree_height = newNode->Depth;
		}
#endif
		newNode->premove = mv;
		newNode->posi->Move(mv);
		newNode->W = 0;
		newNode->L = 0;
		newNode->D = 0;

		return newNode;
}

/************** bittuh's expansion ********************/
void explore(NODE *parent, const int numSimulation, int *parent_WIN, int *parent_LOSE, int *parent_DRAW){
	int i, j;
	NODE *curNode, *newNode;
	MOVLST legal_moves;
    NODE *temp = parent;

	bittuhAllMoveGen(*(parent->posi), legal_moves);

    *parent_WIN = 0;
    *parent_LOSE = 0;
    *parent_DRAW = 0;

	// for all legal moves, add a new node to the game tree
	for(i = 0; i < legal_moves.num; i++){

        newNode = createNode(parent, legal_moves.mov[i]);

		/*************simulate****************/
		int win, lose, draw;
		//bittuhSimulate(newNode->posi, numSimulation, &win, &lose, &draw);
		bittuhGPUSimulate(newNode->posi, &win, &lose, &draw, WORKITEM, numSimulation);
		newNode->L += win;
		*parent_WIN += win;
		newNode->W += lose;
		*parent_LOSE += lose;
		newNode->D += draw;
		*parent_DRAW += draw;
		totalPlay += win+lose+draw;

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
    Node *child;
    int N = 0;
    double x_mean_value = 0, x_square = 0;

    child = node->child;
    while(child != NULL){
        double ucb = UCB(child, EXPLORE_PARA);

        x_mean_value += ucb;
        x_square += ucb*ucb;

        N++;
        child = child->siblg;
    }

    return sqrt((x_square/N)-(x_mean_value/N)*(x_mean_value/N));
//   int L = node->Depth%2? node->W:node->L;
//   int W = node->Depth%2? node->L :node->W;
//   int D = node->D;
//   int Ni = W+L+D;
//   return pow(((double)(W+D/2) - (double)Ni * pow((double)(W+D/2)/Ni, 2))/Ni, 0.5);
}

void progressivePruning(NODE *parent)
{
    if(parent->W + parent->L + parent->D <= 1000)
        return;

    NODE *curNode, *prevNode;
    double worstRightExpected = 100000000;

    curNode = parent->child;
    while(curNode != NULL){
        if(curNode->child == NULL){
            curNode = curNode->siblg;
            continue;
        }

        Node *child = curNode->child;
        int N = 0;
        double x_mean_value = -1;

        while(child != NULL){
            x_mean_value += UCB(child, EXPLORE_PARA);

            N++;
            child = child->siblg;
        }

        double rightExpected = x_mean_value/N+2*standard_devia(curNode);
        if(rightExpected < worstRightExpected){
            worstRightExpected = rightExpected;
            //curNode->posi->Display();
            //fprintf(stderr, "%d \@ W/D/L: %d/%d/%d\n", curNode->Depth, curNode->W, curNode->D, curNode->L);
        }

        curNode = curNode->siblg;
    }

    prevNode = NULL;
    curNode = parent->child;
    while(curNode != NULL){
        if(curNode->child == NULL){
            prevNode = curNode;
            curNode = curNode->siblg;
            continue;
        }

        Node *child = curNode->child;
        int N = 0;
        double x_mean_value = -1;

        while(child != NULL){
            x_mean_value += UCB(child, EXPLORE_PARA);

            N++;
            child = child->siblg;
        }

        double leftExpected = x_mean_value/N-2*standard_devia(curNode);

        if(leftExpected > worstRightExpected){
            if(prevNode != NULL)
                prevNode->siblg = curNode->siblg;
            else
                parent->child = curNode->siblg;
            if(curNode->premove.st == 9 && curNode->premove.ed == 13){
            fprintf(stderr, "Pruned :\n");
            curNode->posi->Display();
            fprintf(stderr, "%d \@ W/D/L: %d/%d/%d\n", curNode->Depth, curNode->W, curNode->D, curNode->L);

            Node *tmp = curNode->child;
            while(tmp != NULL){
            fprintf(stderr, "curNode: W/D/L %d %d %d", tmp->W, tmp->D, tmp->L);
            tmp = tmp->siblg;
            }
            int a;
            scanf("%d", &a);
            }
#ifdef ENABLE_PROFILING
                pruned++;
#endif // ENABLE_PROFILING
        }

        prevNode = curNode;
        curNode = curNode->siblg;
    }
}

/******************* bittuh's selection **********************/
NODE *find_best_child(NODE *parent, double c){
	int i = 0, k = 0, mm, nn, m, n;
	double ucb;

	NODE *curNode;
	NODE *preNode, *nextNode;
	double best_score = UCB(parent->child, c);
	NODE *best_child = parent->child;
	for(curNode = parent->child; curNode!= NULL ;curNode = curNode->siblg){
		ucb = UCB(curNode, c);
		if(best_score <= ucb){
			best_score = ucb;
			best_child = curNode;
		}
	}

	curNode = parent->child;
	preNode = parent;
	nextNode = curNode->siblg;

    //progressivePruning(parent);
	// PRUNE!!!!!!



	return best_child;
}
/***************************************************************/


/************************* bittuh's play ***********************/
NODE* baseMCTSPlay(BOARD *brd, double c){
	root = (NODE*)malloc(sizeof(NODE));
	NODE *curNode = root;
	MOVLST lst;
	int tempWin, tempLose, tempDraw;
	int i, old_totalPlay, k = 0;
	root->parent = NULL;
	root->siblg = NULL;
	root->child = NULL;
	root->posi = brd;
	root->Depth = 0;
#ifdef USING_GPU_SIMULATE
	explore(root, FIRST_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
    explore(root, FIRST_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif
	root->W = tempWin;
	root->L = tempLose;
	root->D = tempDraw;

	// tree-growing loop
	while(GetTickCount()-Tick < TimeOut){
		/********* selection ***********/
		curNode = root;
#ifdef ENABLE_PROFILING
		ENABLE_PROFILING_Tick = GetTickCount();
#endif // ENABLE_PROFILING
		while(curNode->child != NULL){
			curNode = find_best_child(curNode, c);
		}
#ifdef ENABLE_PROFILING
		searchTick += GetTickCount() - ENABLE_PROFILING_Tick;
#endif // ENABLE_PROFILING

		/******* expansion & simulation *******/

#ifdef USING_GPU_SIMULATE
        explore(curNode, OTHER_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
        explore(root, OTHER_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif
		/********* back propogation **********/
		for(k = 0; curNode != NULL; k++){
			if(k%2 == 0){
				curNode->W += tempWin;
				curNode->L += tempLose;
			}
			else{
				curNode->L += tempWin;
				curNode->W += tempLose;
			}
			curNode->D += tempDraw;
			curNode = curNode->parent;
		}
		/**************************************/
	}

	//give out final answer

	return find_best_child(root, c);


}
/**************************************************************/

/****************** play with multi root ****************/
NODE* rootParallelMCTSPlay(BOARD *brd, double c){
	root = (NODE*)malloc(sizeof(NODE));
	MOVLST lst;
	int tempWin, tempLose, tempDraw;
	int i, old_totalPlay, k = 0;
	root->parent = NULL;
	root->siblg = NULL;
	root->child = NULL;
	root->posi = brd;
	root->Depth = 0;
#ifdef USING_GPU_SIMULATE
	explore(root, FIRST_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
    explore(root, FIRST_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif
	root->W = tempWin;
	root->L = tempLose;
	root->D = tempDraw;
    fprintf(stderr, "first level: %d/%d/%d\n", tempWin, tempLose, tempDraw);

    int rootCount = 0;
    NODE *temproot;
    for(temproot = root->child; temproot != NULL; temproot = temproot->siblg) rootCount++;
    NODE **subroot = (NODE**)malloc(sizeof(NODE*)*rootCount);

    temproot = root->child;
    for(int i = 0; i < rootCount; i++){
        subroot[i] = temproot;
        temproot = temproot->siblg;
    }

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tmpWin, tmpLose, tmpDraw;
        #pragma omp for private(i) schedule(dynamic,1)
        for(int i = 0; i < rootCount; i++){
#ifdef USING_GPU_SIMULATE
            explore(subroot[i], OTHER_LVL_GPU_SIMULATE, &tmpWin, &tmpLose, &tmpDraw);
#elif defined USING_CPU_SIMULATE
            explore(subroot[i], OTHER_LEVEL_SIMULATION, &tmpWin, &tmpLose, &tmpDraw);
#endif

            subroot[i]->W += tmpWin;
            subroot[i]->L += tmpLose;
            subroot[i]->D += tmpDraw;
            #pragma omp critical
            {
                fprintf(stderr, "%d/%d/%d\n", tmpWin, tmpLose, tmpDraw);
                root->W += tmpLose;
                root->L += tmpWin;
                root->D += tmpDraw;
            }
        }
    }
    fprintf(stderr, "~~~~\\\\~~~~~~~\n");
	// tree-growing loop
	#pragma omp parallel num_threads(NUM_THREADS)
	{
        #pragma omp for private(tempWin, tempLose, tempDraw, i) schedule(dynamic,1)
        for(int i = 0; i < rootCount; i++){
            NODE *curNode;
            while(GetTickCount()-Tick < TimeOut){
                curNode = subroot[i];
                #ifdef ENABLE_PROFILING
                ENABLE_PROFILING_Tick = GetTickCount();
                #endif // ENABLE_PROFILING

                while(curNode->child != NULL){
                    curNode = find_best_child(curNode, c);
                }
                #ifdef ENABLE_PROFILING
                searchTick += GetTickCount() - ENABLE_PROFILING_Tick;
                #endif // ENABLE_PROFILING


#ifdef USING_GPU_SIMULATE
                explore(curNode, OTHER_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
                explore(curNode, OTHER_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif
                #pragma omp critical
                {
                    fprintf(stderr, "%d/%d/%d\n", tempWin, tempLose, tempDraw);
                    for(k = 0; curNode != NULL; k++){
                        if(k%2 == 0){
                            curNode->W += tempWin;
                            curNode->L += tempLose;
                        }
                        else{
                            curNode->L += tempWin;
                            curNode->W += tempLose;
                        }
                        curNode->D += tempDraw;
                        curNode = curNode->parent;
                    }
                    //fprintf(stderr, "root: %d/%d/%d\n", root->W, root->L, root->D);
                }

            }
        }
	}

	//give out final answer
	return find_best_child(root, c);

}

/***************** play with multi root and sort ****************/
NODE* ourParallelMCTSPlay(BOARD *brd, double c){
	root = (NODE*)malloc(sizeof(NODE));
	MOVLST lst;
	int tempWin, tempLose, tempDraw;
	int i, old_totalPlay, k = 0;
	root->parent = NULL;
	root->siblg = NULL;
	root->child = NULL;
	root->posi = brd;
	root->Depth = 0;
#ifdef USING_GPU_SIMULATE
	explore(root, FIRST_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
    explore(root, FIRST_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif
	root->W = tempWin;
	root->L = tempLose;
	root->D = tempDraw;


    int rootCount = 0;
    NODE *temproot;
    for(temproot = root->child; temproot != NULL; temproot = temproot->siblg) rootCount++;
    NODE **subroot = (NODE**)malloc(sizeof(NODE*)*rootCount);
    omp_lock_t *subrootLock = (omp_lock_t*)malloc(sizeof(omp_lock_t)*rootCount);

    temproot = root->child;
    for(int i = 0; i < rootCount; i++){
        subroot[i] = temproot;
        temproot = temproot->siblg;
    }

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for private(tempWin, tempLose, tempDraw, i) schedule(dynamic,1)
        //for(temproot = root->child; temproot != NULL; temproot = temproot->siblg){
        for(int i = 0; i < rootCount; i++){
#ifdef USING_GPU_SIMULATE
	explore(subroot[i], OTHER_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
    explore(subroot[i], OTHER_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif


            subroot[i]->W += tempWin;
            subroot[i]->L += tempLose;
            subroot[i]->D += tempDraw;
            #pragma omp critical
            {
                root->W += tempLose;
                root->L += tempWin;
                root->D += tempDraw;
            }
        }
    }

    for(int i = 0; i < rootCount; i++)
        omp_init_lock(&(subrootLock[i]));
	// tree-growing loop
	#pragma omp parallel num_threads(NUM_THREADS) private(tempWin, tempLose, tempDraw)
	{
	    NODE *curNode;
            while(GetTickCount()-Tick < TimeOut){
                int index;
                do{

                    double highestUCB = 0;
                    for(int i=0; i < rootCount; i++){
                        double tmp;
                        if((tmp = UCB(subroot[i], c)) > highestUCB){
                            highestUCB = tmp;
                            index = i;
                        }
                    }
                }
                while(!omp_test_lock(&(subrootLock[index])));
                curNode = subroot[index];

                #ifdef ENABLE_PROFILING
                ENABLE_PROFILING_Tick = GetTickCount();
                #endif // ENABLE_PROFILING
                int level = 0;
                while(curNode->child != NULL){
                    level++;
                    curNode = find_best_child(curNode, c);
                }
                #ifdef ENABLE_PROFILING
                searchTick += GetTickCount() - ENABLE_PROFILING_Tick;
                #endif // ENABLE_PROFILING

#ifdef USING_GPU_SIMULATE
                explore(curNode, OTHER_LVL_GPU_SIMULATE, &tempWin, &tempLose, &tempDraw);
#elif defined USING_CPU_SIMULATE
                explore(curNode, OTHER_LEVEL_SIMULATION, &tempWin, &tempLose, &tempDraw);
#endif

                #pragma omp critical
                {
                    for(k = 0; curNode != NULL; k++){
                        if(k%2 == 0){
                            curNode->W += tempWin;
                            curNode->L += tempLose;
                        }
                        else{
                            curNode->L += tempWin;
                            curNode->W += tempLose;
                        }
                        curNode->D += tempDraw;
                        curNode = curNode->parent;
                    }
                }
                omp_unset_lock(&(subrootLock[index]));
            }
	}

	//give out final answer
	return find_best_child(root, c);

}

#ifdef ENABLE_PROFILING
void initProfilingInfo(){
	simulateTick = 0;
	drawCount = 0;
	tree_height = 0;
	highestSimulateDepth = 0;
	totalSimulateDepth = 0;
	searchTick = 0;
	pruned = 0;
}

void showProfilingInfo(){
    NODE *curNode;
    fprintf(stderr, "================== Profling Info ===============\n");
    if(root != NULL){
    for(curNode = root->child; curNode!= NULL ;curNode = curNode->siblg){
        //fprintf(stderr, "%d/%d/%d\n", curNode->W, curNode->L, curNode->D);


        fprintf(stderr, "%d->%d: %d/%d/%d=%d, %.2f/%.2f/%.2f <%.3f>\n",
            curNode->premove.st, curNode->premove.ed,
            curNode->W, curNode->L, curNode->D, curNode->W+curNode->L+curNode->D,
            (double)curNode->W/(double)(curNode->W+curNode->L+curNode->D),
            (double)curNode->L/(double)(curNode->W+curNode->L+curNode->D),
            (double)curNode->D/(double)(curNode->W+curNode->L+curNode->D),
            UCB(curNode, EXPLORE_PARA));

    }
    }
    fprintf(stderr, "total simulate %d times, <%d>.\n", totalPlay, root->W+root->L+root->D);
    fprintf(stderr, "tree height:%d\n", tree_height);
    fprintf(stderr, "pruned %d nodes\n", pruned);
    fprintf(stderr, "avg simulate moves: %ld\n", totalSimulateDepth/totalPlay);
    fprintf(stderr, "highest simulate moves in one iteration: %d\n", highestSimulateDepth);
    fprintf(stderr, "spend %ld to simulate\n", simulateTick);
    fprintf(stderr, "spend %ld to search\n", searchTick);
    fprintf(stderr,"\n");
    fprintf(stderr, "================================================\n");
}
#endif // ENABLE_PROFILING

void initOpenCL(){
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
    numPlatforms = 0;
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get platform number\n");
	platforms = NULL;
	platforms = (cl_platform_id*)malloc(numPlatforms*sizeof(cl_platform_id));
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get platform ID\n");

    /* get device */

	numDevices = 0;
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get device number\n");

	devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	if(status != CL_SUCCESS) fprintf(stderr, "$:%d error when get device ID...\n", status);

    for(int i = 0; i < numDevices; i++){

        char *dvcInfo;
        cl_uint maxComputeUnits, maxDim;
        cl_ulong globmemsize, localmemsize, maxConstantBuf;
        size_t valueSize, maxWorkGroupSize, maxWorkItemSize[3];
        // print device name
        clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 0, NULL, &valueSize);
        dvcInfo = (char*) malloc(valueSize);
        clGetDeviceInfo(devices[i], CL_DEVICE_NAME, valueSize, dvcInfo, NULL);
        printf("## Device [%s]\n", dvcInfo);
        free(dvcInfo);

        clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, 0, NULL, &valueSize);
        dvcInfo = (char*) malloc(valueSize);
        clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, valueSize, dvcInfo, NULL);
        printf("   Vendored by \"%s\"\n", dvcInfo);
        free(dvcInfo);

        // print parallel compute units
        clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS,
        sizeof(maxComputeUnits), &maxComputeUnits, NULL);
        printf(" + Parallel compute units: %u\n", maxComputeUnits);

        // print hardware device version
        clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(globmemsize), &globmemsize, NULL);
        printf(" + Global Memory Size: %lu bytes\n", globmemsize);

        // print hardware device version
        clGetDeviceInfo(devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localmemsize), &localmemsize, NULL);
        printf(" + Local Memory Size: %lu bytes\n", localmemsize);

        // print hardware device version
        clGetDeviceInfo(devices[i], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(maxConstantBuf), &maxConstantBuf, NULL);
        printf(" + Max Constant Memory Size: %lu bytes\n", maxConstantBuf);

        // print hardware device version
        clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(maxWorkGroupSize), &maxWorkGroupSize, NULL);
        printf(" + Max Work Group Size: %zu \n", maxWorkGroupSize);

        // print hardware device version
        clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(maxDim), &maxDim, NULL);
        clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, maxDim*sizeof(size_t), &maxWorkItemSize, NULL);
        printf(" + Max Work Items Size in all %u dimension: ", maxDim);
        for(int j = 0; j < maxDim; j++) printf("%zu ", maxWorkItemSize[j]);
        printf("\n");



    }
    /* context */

	context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &status);
	if(status != CL_SUCCESS) fprintf(stderr, "context create err\n");

    /* command queue */
    cmdQueue = (cl_command_queue*)malloc(sizeof(cl_command_queue) * numDevices);
    queueCounter = (int*)malloc(sizeof(int) * numDevices);
    for(int i = 0; i < numDevices; i++){
        cmdQueue[i] = clCreateCommandQueue(context, devices[i], CL_QUEUE_PROFILING_ENABLE, &status);
        if(status != CL_SUCCESS) fprintf(stderr, "command queue create err\n");
    }


    /* program */
    char msg[4096];
    size_t len;
    program = clCreateProgramWithSource(context, 1, (const char**)&kernelsrc, (const size_t*)&srcsize, &status);

	if(status != CL_SUCCESS) fprintf(stderr, "create program error\n");
	status = clBuildProgram(program, numDevices, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS){
        fprintf(stderr, "build error %d\n", status);
        status = clGetProgramBuildInfo(program, devices[1], CL_PROGRAM_BUILD_LOG, 4096*sizeof(char), msg, &len);
        fprintf(stderr, "%d\n%s\n", status, msg);
    }

    //create CL memory buffer object

    if(status != CL_SUCCESS) fprintf(stderr, "fuck1!\n");
    for(int i = 0; i < NUM_THREADS; i++){
        bufResult[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, (WORKITEM * RESULT_PER_WORKITEM)*sizeof(int), NULL, &status);
        if(status != CL_SUCCESS){
            fprintf(stderr, "fuck2!\n");
            return;
        }
    }
}

int main() {

    initOpenCL();

#ifdef ENABLE_PROFILING
    initProfilingInfo();
#endif // ENABLE_PROFILING

	srand(Tick=GetTickCount());

	BOARD BBB;

	MOV mymove;
	TimeOut=(BBB.LoadGame("board3.txt")-3)*1000;
	totalPlay = 0;


    struct timeval  tv1, tv2;
    gettimeofday(&tv1, NULL);
	BBB.Display();

    NODE* bestNode;
	if(!BBB.ChkLose()){

        //bestNode = baseMCTSPlay(&BBB,EXPLORE_PARA);

        //bestNode = rootParallelMCTSPlay(&BBB,EXPLORE_PARA);

        bestNode = ourParallelMCTSPlay(&BBB,EXPLORE_PARA);


        Output(mymove = bestNode->premove);
	}

#ifdef ENABLE_PROFILING
    showProfilingInfo();
#endif // ENABLE_PROFILING

    if(mymove.st != mymove.ed)fprintf(stderr,"best: %d->%d\n",mymove.st, mymove.ed);
    else fprintf(stderr,"best: flip %d\n",mymove.st);
    BBB.Display();

	return 0;
}
