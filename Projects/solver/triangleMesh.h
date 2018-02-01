#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename T, int dim>
class TriangleMesh
{
public:
	TriangleMesh<T, dim>(int n) : numParticles(n), particles(n)
	{
		std::cout << "Number of Particles: " << n << std::endl;

		if(numParticles < 3) {
			std::cout << "ERROR" << std::endl;
		}
		else if(numParticles == 3) {
			std::vector<int> idx;
			idx.push_back(0);
			idx.push_back(1);
			idx.push_back(2);
			indices.push_back(idx);
		}
		else {
			for(int i = 0; i < numParticles; i++) {
				std::vector<int> index;
				int idx1 = i;
				int idx2, idx3;

				if(i + 2 >= numParticles) {
					idx3 = 1;
				}
				else {
					idx3 = i + 2;
				}

				if(i + 1 >= numParticles) {
					idx2 = 0;
				}
				else {
					idx2 = i + 1;
				}

				index.push_back(idx1);
				index.push_back(idx2);
				index.push_back(idx3);
				indices.push_back(index);
			}
		}
	}

	int numParticles;
	Particles<T, dim> particles;
	std::vector<std::vector<int>> indices;

	void update(float t)
	{
		int count = 0;
		for(int i = 0; i < particles.numParticles; i++) {
			if(count == 2) {
				count = 0;
				particles.positions.at(i)[0] = particles.positions.at(i)[0] + particles.velocities.at(i)[0] * count; 
				particles.positions.at(i)[1] = particles.positions.at(i)[1] + particles.velocities.at(i)[1] * count; 
				particles.positions.at(i)[2] = particles.positions.at(i)[2] + particles.velocities.at(i)[2] * count; 
			}
			else {
				particles.positions.at(i)[0] = particles.positions.at(i)[0] - particles.velocities.at(i)[0] * count; 
				particles.positions.at(i)[1] = particles.positions.at(i)[1] - particles.velocities.at(i)[1] * count; 
				particles.positions.at(i)[2] = particles.positions.at(i)[2] - particles.velocities.at(i)[2] * count; 
			}

			count++;
		}
	}
};
