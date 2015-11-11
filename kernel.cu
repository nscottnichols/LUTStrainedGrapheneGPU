
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <algorithm>
#include "gputimer.h"

#define NUMX 256
#define NUMY 256
#define STRAIN 0.00
#define POISSON 0.165
#define A0 1.42
#define NUM_THREADS 262144
#define BLOCK_WIDTH 1024
#define RESOLUTION 100
#define ZRES 100
#define SMALL_Z 0.00001
#define Z_MAX 20.0



struct atom_t {
	float atompos, vLJ;

	friend std::istream& operator >> (std::istream& ins, atom_t& r);
};

std::istream& operator >> (std::istream& ins, atom_t& r) {
	ins >> r.atompos >> r.vLJ;
	return ins;
};

template<class T>
void print_array(T *array, int size)
{
	printf("{ ");
	for (int i = 0; i < size; i++)  { printf("%e ", array[i]); }
	printf("}\n");
}

__global__ void potentialLJ(
	float epsilon,
	float sigma,
	int x,
	int y,
	int z,
	float *pos_carbonx,
	float *pos_carbony,
	float *pos_atomx,
	float *pos_atomy,
	float *pos_atomz,
	float *potLJ,
	int numCarbons) 
{
	const int2 thread_2D_pos = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y);

	const int thread_1D_pos = thread_2D_pos.y * numCarbons + thread_2D_pos.x;
	float cx = pos_carbonx[thread_1D_pos];
	float cy = pos_carbony[thread_1D_pos];

	float rx = pos_atomx[x] - cx;
	float ry = pos_atomy[y] - cy;
	float rz = pos_atomz[z];
	float invr = rsqrtf(rx * rx + ry * ry + rz * rz);
	float sor6 = sigma * sigma * sigma * sigma * sigma * sigma * invr * invr * invr * invr * invr * invr;
	float vLJ = 4.0 * epsilon * sor6 * (sor6 - 1.0);
	// accumulate effect of all particles
	atomicAdd(potLJ, vLJ);
	
}

int main(int argc, char **argv)
{
	GpuTimer timer;
	// declare variables
	int numCarbons = NUM_THREADS;
	int numx = NUMX;
	int numy = NUMY;
	
	float epsilon = 16.2463;
	float sigma = 2.74;
	float smallz = SMALL_Z;
	float zmax = Z_MAX;
	//////////////////////////////////////////////////////////////////////////////	
	// Create graphene lattice
	//////////////////////////////////////////////////////////////////////////////
	float ax = A0 * (1.00 + (STRAIN * (1.0 - (3.0 * POISSON)) / 4.0)); // New Way
	float ay = A0 * (1.00 + STRAIN);
	float ay2 = (3.0/8.0)*A0*(4.0+(3.0*STRAIN)-(STRAIN*POISSON));
	float d0 = sqrt((float)3.0) * ax;
	// Create a vector object that contains numx elements.
	std::vector<float> transx;
	for (int i = 0; i < numx; ++i) {
		transx.push_back((i - numx / 2)*d0);
	}
	std::vector<float> transx2;
	for (int i = 0; i < numx; ++i) {
		transx2.push_back((d0 / 2) + (i - numx / 2)*d0);
	}
	std::vector<float> x(((transx.size() + transx2.size()) * 2 * numy));
	for (int k = 0; k < numy; ++k) {
		for (int j = 0; j < numy; ++j) {
			for (int i = k * 4 * numy; i < x.size(); i += 4) {
				x[i] = transx[k];
				x[i + 1] = transx[k];
				x[i + 2] = transx2[k];
				x[i + 3] = transx2[k];
			}
		}
	}
	// Create a vector object that contains numy elements.
	std::vector<float> transy;
	for (int i = 0; i < numy; ++i) {
		transy.push_back((((ay + ax) / 2) + (ax + ay + ay)*(i - numy / 2)));
	}
	std::vector<float> transy2;
	for (int i = 0; i < numy; ++i) {
		transy2.push_back(((-(ay + ax) / 2) + (ax + ay + ay)*(i - numy / 2)));
	}
	std::vector<float> transy3;
	for (int i = 0; i < numy; ++i) {
		transy3.push_back(((ay / 2) + (ax + ay + ay)*(i - numy / 2)));
	}
	std::vector<float> transy4;
	for (int i = 0; i < numy; ++i) {
		transy4.push_back(((-ay / 2) + (ax + ay + ay)*(i - numy / 2)));
	}
	std::vector<float> tempy(4 * numx);
	for (int i = 0, k = 0; k < numy; i += 4, ++k) {
		tempy[i] = transy[k];
		tempy[i + 1] = transy2[k];
		tempy[i + 2] = transy3[k];
		tempy[i + 3] = transy4[k];
	}
	std::vector<float> y;
	for (int i = 0; i < numx; ++i) {
		copy(tempy.begin(), tempy.end(), back_inserter(y));
	}




	//////////////////////////////////////////////////////////////////////////////	
	// Test positions
	//////////////////////////////////////////////////////////////////////////////

	std::vector<float> pos_atomx;
	for (int i = 0; i < RESOLUTION + 1; ++i) {
		pos_atomx.push_back(i*d0/(2*(RESOLUTION)));
	}

	std::vector<float> pos_atomy;
	for (int i = 0; i < RESOLUTION + 1; ++i) {
		pos_atomy.push_back(i*ay2/(RESOLUTION));
	}

	std::vector<float> pos_atomz;
	for (int i = 0; i < ZRES + 1; ++i) {
		pos_atomz.push_back(smallz + i*zmax/(RESOLUTION));
	}


	//////////////////////////////////////////////////////////////////////////////
	// declare and allocate host memory
	//////////////////////////////////////////////////////////////////////////////
	float* h_pos_carbonx = &x[0];
	float* h_pos_carbony = &y[0];
	float* h_pos_atomx = &pos_atomx[0];
	float* h_pos_atomy = &pos_atomy[0];
	float* h_pos_atomz = &pos_atomz[0];
	float h_potLJ[RESOLUTION*RESOLUTION*ZRES];
	float h_potLJtemp[0];
	const int CARBON_BYTES = NUM_THREADS * sizeof(float);
	const int RES_BYTES = RESOLUTION * sizeof(float);
	const int ZRES_BYTES = ZRES * sizeof(float);
	const int POTLJ_BYTES = RESOLUTION * RESOLUTION * ZRES * sizeof(float);
	

	memset(h_potLJ, 0.0, POTLJ_BYTES);
	memset(h_potLJtemp, 0.0, sizeof(float));


	// declare, allocate, and zero out GPU memory
	float *d_pos_carbonx;
	float *d_pos_carbony;
	float *d_pos_atomx;
	float *d_pos_atomy;
	float *d_pos_atomz;
	float *d_potLJ;
	cudaMalloc((void **)&d_pos_carbonx, CARBON_BYTES);
	cudaMalloc((void **)&d_pos_carbony, CARBON_BYTES);
	cudaMalloc((void **)&d_pos_atomx, RES_BYTES);
	cudaMalloc((void **)&d_pos_atomy, RES_BYTES);
	cudaMalloc((void **)&d_pos_atomz, ZRES_BYTES);
	cudaMalloc((void **)&d_potLJ, sizeof(float));

	// now copy data from host memory to device memory
	cudaMemcpy((void *)d_pos_carbonx, (void *)h_pos_carbonx, CARBON_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_pos_carbony, (void *)h_pos_carbony, CARBON_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_pos_atomx, (void *)h_pos_atomx, RES_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_pos_atomy, (void *)h_pos_atomy, RES_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_pos_atomz, (void *)h_pos_atomz, ZRES_BYTES, cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_potLJ, (void *)h_potLJtemp, sizeof(float), cudaMemcpyHostToDevice);


	//////////////////////////////////////////////////////////////////////////////
	// launch the kernel
	//////////////////////////////////////////////////////////////////////////////
	printf("%d total threads in %d blocks writing into %d array elements\n",
		NUM_THREADS, NUM_THREADS / BLOCK_WIDTH, NUM_THREADS);	
	timer.Start();
	
	for (int i = 0; i < RESOLUTION; ++i) {
		for (int j = 0; j < RESOLUTION; ++j){
			for (int k = 0; k < ZRES; ++k){

	potentialLJ << <NUM_THREADS / BLOCK_WIDTH, BLOCK_WIDTH >> >(
		epsilon,
		sigma,
		i,
		j,
		k,
		d_pos_carbonx,
		d_pos_carbony,
		d_pos_atomx,
		d_pos_atomy,
		d_pos_atomz,
		d_potLJ,
		numCarbons);
	

	// copy back the array of sums from GPU and print
	cudaMemcpy(h_potLJtemp, d_potLJ, sizeof(float), cudaMemcpyDeviceToHost);
	h_potLJ[i * RESOLUTION * ZRES + j * ZRES + k] = h_potLJtemp[0];

	memset(h_potLJtemp, 0.0, sizeof(float));
	cudaMemcpy((void *)d_potLJ, (void *)h_potLJtemp, sizeof(float), cudaMemcpyHostToDevice);
			}
		}
	}


	// End LOOP

	
	//std::ofstream xout("carbonx.txt");
	//for (int i = 0; i < numCarbons; i++)
	//{
	//	xout << h_pos_carbonx[i] << std::endl; //writing ith character of array in the file
	//}
	//std::ofstream yout("carbony.txt");
	//for (int i = 0; i < numCarbons; i++)
	//{
	//	yout << h_pos_carbony[i] << std::endl; //writing ith character of array in the file
	//}

	std::ofstream fout("out.txt");
	for (int i = 0; i < RESOLUTION*RESOLUTION*ZRES; i++)
	{
		fout << h_potLJ[i] << std::endl; //writing ith character of array in the file
	}
	timer.Stop();
	

	printf("Time elapsed = %g ms\n", timer.Elapsed());
	//print_array(&pos_atomx[0], RESOLUTION+1);
	//print_array(&pos_atomy[0], RESOLUTION+1);
	//print_array(&pos_atomz[0], ZRES+1);

	// free GPU memory allocation and exit
	cudaFree(d_pos_carbonx);
	cudaFree(d_pos_carbony);
	cudaFree(d_pos_atomx);
	cudaFree(d_pos_atomy);
	cudaFree(d_pos_atomz);
	cudaFree(d_potLJ);
	return 0;
}


