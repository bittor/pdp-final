#include <stdio.h>
#include <stdlib.h>

#include <CL/cl.h>

const char *src =
"#pragma OPENCL EXTENSION cl_amd_printf : enable                                \n"
"__kernel                                                                       \n"
"void add                                                                       \n"
"(__global char *array, __global int *result, __global char *result1)                                   \n"
"{                                                                              \n"
"   int idx = get_global_id(0);                                                 \n"
"   result[idx] = (int)array[idx]+(int)array[idx+1];                            \n"
"   result1[idx] = array[idx];                                                 \n"
"   result1[idx] = array[idx+1];                                                \n"
"}                                                                              \n";

const char* programSource1 =
"__kernel                                                                       \n"
"void reduce_array																\n"
"(__global int *array, __global int *rrow, const int n)							\n"
"{																				\n"
"	int matrix_idx = get_global_id(0);											\n"
"	int row_idx = get_global_id(1);												\n"
"	int i;																		\n"
"	rrow[matrix_idx * n + row_idx] = 0;											\n"
"	for(i = 0; i < n; i++)														\n"
"		rrow[matrix_idx * n + row_idx] += array[matrix_idx*n*n + row_idx*n + i];\n"
"}																				\n"
;

const char* programSource2 =
"__kernel																		\n"
"void sum_v1																	\n"
"(__global int *array, __global int *rrow, __global int *sum, const int n)		\n"
"{																				\n"
"	int matrix_idx = get_global_id(0);											\n"
"	int i;																		\n"
"	sum[matrix_idx] = 0;														\n"
"	for(i = 0; i < n*n; i++)													\n"
"		sum[matrix_idx] += array[matrix_idx*n*n + i] * rrow[matrix_idx*n + i%n];\n"
"}																				\n"
;


int main(int argc, char *argv[]){
	int i, m, n;
	int *array;
	int *sum;
	FILE *fp;
    char msg[2048];

	/*
	fp = fopen(argv[1], "rb");
	fread(&m, sizeof(int), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	array = (int*)malloc(sizeof(int) * m * n * n);
	sum = (int*)malloc(sizeof(int) * m);
	fread(array, sizeof(int), m * n * n, fp);
	fclose(fp);
	*/


	cl_int status;


	cl_uint numPlatforms = 0;
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get platform number\n");
	cl_platform_id *platforms = NULL;
	platforms = (cl_platform_id*)malloc(numPlatforms*sizeof(cl_platform_id));
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get platform ID\n");


	cl_uint numDevices = 0;
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get device number\n");
	cl_device_id *devices;
	devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	if(status != CL_SUCCESS)
		fprintf(stderr, "error when get device ID\n");




	cl_context context;
	context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &status);
	if(status != CL_SUCCESS) fprintf(stderr, "context create err\n");


	cl_command_queue cmdQueue;
	cmdQueue = clCreateCommandQueue(context, devices[0], 0, &status);
	if(status != CL_SUCCESS) fprintf(stderr, "command queue create err\n");

	/*
	cl_mem bufarray;
	bufarray = clCreateBuffer(context, CL_MEM_READ_ONLY, m*n*n*sizeof(int), NULL, &status);
	cl_mem bufrrow;
	bufrrow = clCreateBuffer(context, CL_MEM_READ_WRITE, m*n*sizeof(int), NULL, &status);
	cl_mem bufsum;
	bufsum = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
	*/
	int myarray[2] = {0<<24+1<<16+2<<8+3, 4<<24+5<<16+6<<8+7};
	int result[4];

	cl_mem bufarray;
	bufarray = clCreateBuffer(context, CL_MEM_READ_ONLY, 8*sizeof(int8_t), NULL, &status);
	cl_mem bufaaa;
	bufaaa = clCreateBuffer(context, CL_MEM_READ_WRITE, 8*sizeof(int8_t), NULL, &status);
	cl_mem bufResult;
	bufResult = clCreateBuffer(context, CL_MEM_READ_WRITE, 4*sizeof(int), NULL, &status);

	status = clEnqueueWriteBuffer(cmdQueue, bufarray, CL_TRUE, 0, 2*sizeof(int), myarray, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(cmdQueue, bufarray, CL_TRUE, 0, m*n*n*sizeof(int), array, 0, NULL, NULL);

	cl_program program = clCreateProgramWithSource(context, 1, (const char**)&src, NULL, &status);
	if(status != CL_SUCCESS) fprintf(stderr, "create program error\n");
	status = clBuildProgram(program, numDevices, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS){
		fprintf(stderr, "create program error %d\n", status);
		int len;
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 2048*sizeof(char), msg, &len);
        fprintf(stderr, "%s\n", msg);
    }

	/*
	cl_program program1 = clCreateProgramWithSource(context, 1, (const char**)&programSource1, NULL, &status);
	cl_program program2 = clCreateProgramWithSource(context, 1, (const char**)&programSource2, NULL, &status);
	status = clBuildProgram(program1, numDevices, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS) fprintf(stderr, "create program1 error\n");
	status = clBuildProgram(program2, numDevices, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS) fprintf(stderr, "create program1 error\n");
	*/
	cl_kernel kernel;
	kernel = clCreateKernel(program, "add", &status);
	if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~%d\n", status);
    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufarray);
    if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~%d\n", status);
    status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufResult);
    if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~%d\n", status);
    status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufaaa);
    if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~%d\n", status);

	/*
	cl_kernel kernel1;
	kernel1 = clCreateKernel(program1, "reduce_array", &status);
	status |= clSetKernelArg(kernel1, 0, sizeof(cl_mem), &bufarray);
	if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~\n");
	status = clSetKernelArg(kernel1, 1, sizeof(cl_mem), &bufrrow);
	if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~~\n");
	status = clSetKernelArg(kernel1, 2, sizeof(int), &n);
	if(status != CL_SUCCESS) fprintf(stderr, "create kernel1 error~~~\n");


	cl_kernel kernel2;
	kernel2 = clCreateKernel(program2, "sum_v1", &status);
	status |= clSetKernelArg(kernel2, 0, sizeof(cl_mem), &bufarray);
	status |= clSetKernelArg(kernel2, 1, sizeof(cl_mem), &bufrrow);
	status |= clSetKernelArg(kernel2, 2, sizeof(cl_mem), &bufsum);
	status |= clSetKernelArg(kernel2, 3, sizeof(int), &n);
	if(status != CL_SUCCESS) fprintf(stderr, "create kernel2 error\n");
	*/

	size_t globalWorkSize[1];
	globalWorkSize[0] = 4;
	status = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);

	status = clFinish(cmdQueue);
	/*
	size_t globWorkSize1[2];
	globWorkSize1[0] = m;
	globWorkSize1[1] = n;
	status = clEnqueueNDRangeKernel(cmdQueue, kernel1, 2, NULL, globWorkSize1, NULL, 0, NULL, NULL);
	status = clFinish(cmdQueue);


	size_t globWorkSize2[1];
	globWorkSize2[0] = m;
	status = clEnqueueNDRangeKernel(cmdQueue, kernel2, 1, NULL, globWorkSize2, NULL, 0, NULL, NULL);
	status = clFinish(cmdQueue);
	*/
	/*
	clEnqueueReadBuffer(cmdQueue, bufsum, CL_TRUE, 0, m*sizeof(int), sum, 0, NULL, NULL);

	for(i = 0; i < m; i++){
		printf("%d\n", sum[i]);
		if(m <15) fprintf(stderr,"%d\n",sum[i]);
		}
		*/
		int aaarray[2];
    clEnqueueReadBuffer(cmdQueue, bufResult, CL_TRUE, 0, 4*sizeof(int), result, 0, NULL, NULL);
    clEnqueueReadBuffer(cmdQueue, bufaaa, CL_TRUE, 0, 2*sizeof(int), aaarray, 0, NULL, NULL);
    for(i = 0; i < 4; i++)
        printf("%d\n", result[i]);
        printf("[%d]\n", aaarray[0]);
	return 0;

}
