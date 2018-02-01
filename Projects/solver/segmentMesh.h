#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>

enum options
{
	POSITIONS,
	VELOCITIES,
	MASSES
};

template<typename T, int dim>
class SegmentMesh
{
public:
	SegmentMesh<T, dim>() : k(0), gamma(0), restLen(0) {}

	SegmentMesh<T, dim>(int n) : 	numParticles(n), particles(n),
									k(0), gamma(0), restLen(0)
	{
		updateIndices();	
	}

	int numParticles;
	Particles<T, dim> particles;	
	std::vector<std::vector<int>> indices;

	T k;		// Young's Modulus
	T gamma;	// Damping Coefficient
	T restLen;	// Rest Length

	// Connects particles at [i, i + 1] or [(n - 1), 0]
	void updateIndices()
	{
		// TODO
		// Come up with a better way to connect points
		
		for(int i = 0; i < particles.numParticles; i++) {
			std::vector<int> poly;
			
			// Last particle
			if(i + 1 >= particles.numParticles) {
				poly.push_back(i);
				poly.push_back(0);
			}
			else {
				poly.push_back(i);
				poly.push_back(i + 1);
			}

			indices.push_back(poly);
		}
	}

	void printParticles(int option)
	{
		switch(option)
		{
			case POSITIONS:
				std::cout << "Positions" << std::endl;
				particles.printPositions();
				break;
			case VELOCITIES:
				std::cout << "Velocities" << std::endl;
				particles.printVelocities();
				break;	
			case MASSES:
				std::cout << "Masses" << std::endl;
				particles.printMasses();
				break;
			default:
				particles.printProperties();
				break;				
		}
	}

	void printIndices()
	{
		std::cout << "Indices" << std::endl;
		for(size_t i = 0; i < indices.size(); i++) {
			std::cout << "Poly " << i << ":";
			for(size_t j = 0; j < indices.at(i).size(); j++) {
				std::cout << " " << indices.at(i).at(j);
			}
			std::cout << std::endl;
		}
	}

	void printYoungsMod()
	{
		std::cout << "Young's Modulus: ";
		std::cout << k << std::endl;
	}

	void printDampingCoeff()
	{
		std::cout << "Damping Coefficient: ";
		std::cout << gamma << std::endl;
	}

	void printRestLength()
	{
		std::cout << "Rest Length: ";
		std::cout << restLen << std::endl;
	}
};
