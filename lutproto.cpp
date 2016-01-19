

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <algorithm>

#define STRAIN 0.00
#define POISSON 0.165
#define STRETCH 3
#define RESOLUTION 3
#define ZRES 3
#define FILENAME "out.txt"
#define A0 1.42
#define ZBOUND 20.0

struct potential_t {
	float atompos, vLJ;

	friend std::istream& operator >> (std::istream& ins, potential_t& r);
};

std::istream& operator >> (std::istream& ins, potential_t& r) {
	ins >> r.vLJ;
	return ins;
};


template<class T>
void print_arraye(T *array, int size)
{
	printf("{ ");
	for (int i = 0; i < size; i++)  { printf("%e ", array[i]); }
	printf("}\n");
}

template<class T>
void print_arrayf(T *array, int size)
{
	printf("{ ");
	for (int i = 0; i < size; i++)  { printf("%f ", array[i]); }
	printf("}\n");
}

template<class T>
void print_arrayd(T *array, int size)
{
	printf("{ ");
	for (int i = 0; i < size; i++)  { printf("%d ", array[i]); }
	printf("}\n");
}



int main(int argc, char **argv)
{
	float ax = (sqrt(3.0)/8.0)*A0*(4.0+STRAIN-(3.0*STRAIN*POISSON));
	float ay = (3.0/8.0)*A0*(4.0+(3.0*STRAIN)-(STRAIN*POISSON));

	float xb = STRETCH*ax;
	float yb = STRETCH*ay;
	float zb = ZBOUND;


	//Read vLJ from file
	std::vector<float> potential_vLJ;
	potential_t potential_datum;
	std::ifstream f(FILENAME);
	while (f) {
		f >> potential_datum;
		potential_vLJ.push_back(potential_datum.vLJ);
		}

	//xy indices for LUT to LUT
	std::vector<int> xy_ind;
	if(STRETCH%2){
		for (int i = 0; i < (2*RESOLUTION*STRETCH+1); ++i) {
			xy_ind.push_back(abs(i%((RESOLUTION)*2)-(RESOLUTION)));
		}
	}
	else{
		for (int i = 0; i < (2*RESOLUTION*STRETCH+1); ++i) {
			xy_ind.push_back(abs(abs(i%((RESOLUTION)*2)-(RESOLUTION))-(RESOLUTION)));
		}
	}
	
	std::vector<int> z_ind;
	for (int i = 0; i < (2*ZRES+1); ++i) {
			z_ind.push_back(abs(i-ZRES));
		}

	

	//x,y, and z positions
	std::vector<float> pos_atomx;
	for (float i = 0; i < (2*RESOLUTION*STRETCH+1); ++i) {
		pos_atomx.push_back(xb * ((i/(RESOLUTION*STRETCH))-1.0));
	}

	std::vector<float> pos_atomy;
	for (float i = 0; i < (2*RESOLUTION*STRETCH+1); ++i) {
		pos_atomy.push_back(yb * ((i/(RESOLUTION*STRETCH))-1.0));
	}

	std::vector<float> pos_atomz;
	for (float i = 0; i < (2*ZRES+1); ++i) {
		pos_atomz.push_back(zb * ((i/(ZRES))-1.0));
	}
	
	#ifdef DEBUG
	print_arrayd(&xy_ind[0], (2*RESOLUTION*STRETCH+1));
	print_arrayd(&z_ind[0], (2*ZRES+1));
	print_arrayf(&pos_atomx[0], (2*RESOLUTION*STRETCH+1));
	print_arrayf(&pos_atomy[0], (2*RESOLUTION*STRETCH+1));
	print_arrayf(&pos_atomz[0], (2*ZRES+1));
	#endif
	return 0;
}


