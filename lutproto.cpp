

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <algorithm>

#define STRAIN 0.00
#define POISSON 0.165
#define STRETCH 4
#define RESOLUTION 100
#define ZRES 100
#define FILENAME "out.txt"


struct potential_t {
	float atompos, vLJ;

	friend std::istream& operator >> (std::istream& ins, potential_t& r);
};

std::istream& operator >> (std::istream& ins, potential_t& r) {
	ins >> r.vLJ;
	return ins;
};


template<class T>
void print_array(T *array, int size)
{
	printf("{ ");
	for (int i = 0; i < size; i++)  { printf("%e ", array[i]); }
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
	std::vector<float> xy_ind;
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
	print_array(&xy_ind[0], (2*RESOLUTION*STRETCH+1)	);
	return 0;
}


