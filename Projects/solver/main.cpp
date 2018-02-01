#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "global.h"

enum SIMULATION_TYPE
{
	PARTICLES = 0,
	LINES = 1,
	MESH = 2
};

void simulate(int type, int numParticles)
{
	std::string name, ext;
	switch(type)
	{
		case PARTICLES:
			name = "output_files/point_cloud/particle";
			ext = ".bgeo";

			generatePointCloud(name, ext, numParticles);

			break;
		case LINES:
			name = "output_files/segment_mesh/curve";
			ext = ".poly";

			generateSegmentMesh(name, ext, numParticles);

			break;
		case MESH:
			name = "output_files/3d_mesh/mesh";
			ext = ".obj";
			
			generateTriangleMesh(name, ext, numParticles);			
	
			break;
		default:
			break;
	}
}

int main(int argc, char* argv[])
{
	int sim = std::atoi(argv[1]);

	if(sim) {
		int simType = std::atoi(argv[2]);
		int numParticles = std::atoi(argv[3]);
		simulate(simType, numParticles);
	}
	else {
		std::string input = argv[2];
		std::string output = argv[3];
		parse(input, output);
	}

	return 0;
}
