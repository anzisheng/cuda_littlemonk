#include <cstdio>
#include <stdio.h>
#include <math.h>
//                     
//                   b00 b01 b02 b03
//                   b10 b11 b12 b13
//                   b20 b21 b22 b23
//                   b30 b32 b32 b33
//
//a00 a01 a02 a03    c00 c01 c02 c03
//a10 a11 a12 a13    c10 c12 c13 c14
//a20 a21 a22 a23    c20 c21 c22 c23
//a30 a31 a32 a33    c30 c31 c32 c33
//
//c21 = a20 *b01 + a21*b11 + c21 = c20 *b01 + c21*b11 + aa

__golbal__ void gpu_matrix(int* a, int* b, int* c, int M, int N, int K);
 void cpu_matrix(int* a, int *b, int * c, int M, int N, int K);
#define M 1000
#define N 500
#define K 1000

#define BLOCK_SIZE 16

__managed__ int a[M*N];
__managed__ int b[N*K];
__managed__ int c_gpu[M*K];
__managed__ int c_cpu[M*K];
__global__ void gpu_matrix(int* a, int* b, int* c_gpu, int m, int n, int k)
{
    __shared__ int sub_a[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ int sub_b[BLOCK_SIZE][BLOCK_SIZE];
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int tmp = 0;
    int idx = 0;
    for(int step = 0; step < N/BLOCK_SIZE; step++)
    {
        int step_x =  step * BLOCK_SIZE + threadIdx.x;
        int step_y = y;
        idx = step_y * N + step_x;
        if(step_x > N || step_y > M)
        {
            sub_a[threadIdx.y][threadIdx.x] = 0;
        }
        else {
            sub_a[threadIdx.y][threadIdx.x] = a.[idx];

        }
        step_x = x;
        step_y = step * BLOCK_SIZE+ threadIdx.y;
        idx = step*k+step_x;
        
    }

}
void cpu_matrix(int* a, int*b, int* c_cpu, int m, int n, int k )
{
    for(int y = 0; y < m; y++)
    {
        for(int x = 0; x < k; x++)
        {
            int  temp = 0;
            for(int step = 0; step < n; step++)
            {
                temp+= a[y*n + step] * b[step*k+ x];

            } 
            c_cpu[y*k + x] = temp;
        }
    }


}


 int main()
 {
    for(int y=0; y < M ; y++)
    {
        for(int x= 0; x <N; x++)
        {
            a[y*N+x] = rand()%1024;
        }
    }

    for(int y=0; y < N; y++)
    {
        for(int x = 0; x < K; x++)
        b[y*K+x] = rand()%1024;
    }

// define thread dimensions
    unsigned int grid_x = (K + BLOCK_SIZE-1)/BLOCK_SIZE;
    unsigned int grid_y = (M + BLOCK_SIZE-1)/BLOCK_SIZE;

    dim3 dimGrid(grid_x, grid_y);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

    gpu_matrix(a, b, c_gpu, M, N, K);
    cpu_matrix(a, b, c_cpu, M,N,K);
    bool errors = false;
    for(int y = 0; y < M; y++)
    {
        for(int x = 0; x < K; x++)
        {
            if(fabs(c_cpu[y*K +x] - c_gpu[y*K+x]) > 1.0e-10)
            {
                errors = true;

            }
        }
    }
    printf("result: %s \n", errors?"errors":"passed");

    return 0; 
    
    
 }