#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stdlib.h>
#include <math.h>

template<typename T, int dim>
class Particles
{
public: 
	Particles<T, dim>() : numParticles(0) {}

	Particles<T, dim>(int n) : numParticles(n) 
	{
		// https://stackoverflow.com/questions/41868817/fill-dynamic-vector-using-eigen-library
		Eigen::VectorXd massVector = Eigen::VectorXd(numParticles);
	
		for(int i = 0; i < numParticles; i++) {
			float mass;
			Eigen::Matrix<T, 3, 1> pos, vel, goal;

			mass = (T)1;
			// -n to n
			pos = Eigen::Matrix<T, 3, 1>(rand() % n - 1, rand() % n - 1, rand() % n - 1);
			// -3 to 3
			vel = Eigen::Matrix<T, 3, 1>(sin(rand() % 7 - 3), sin(rand() % 7 - 3), sin(rand() % 7 - 3));

			massVector[i] = mass;

			positions.push_back(pos);
			velocities.push_back(vel);
		}	

		masses = massVector.array().matrix().asDiagonal();
	}

	int numParticles; 
	// https://stackoverflow.com/questions/19821144/eigen-dynamic-matrix-initialization
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> masses;		
	std::vector<Eigen::Matrix<T, dim, 1>> positions;
	std::vector<Eigen::Matrix<T, dim, 1>> velocities;
	std::vector<Eigen::Matrix<T, dim, 1>> forces;
	std::vector<Eigen::Matrix<T, dim, 1>> dampForces;

	void updatePositions(T t)
	{
		// bool pingpong = false;
		for(int i = 0; i < numParticles; i++) {
			positions.at(i)[0] = positions.at(i)[0] + t * velocities.at(i)[0];
			positions.at(i)[1] = positions.at(i)[1] + t * velocities.at(i)[1];
			positions.at(i)[2] = positions.at(i)[2] + t * velocities.at(i)[2];
		}
	}

	void updateMasses()
	{
		Eigen::VectorXd massVector = Eigen::VectorXd(numParticles);
		for(int i = 0; i < numParticles; i++) {
			massVector[i] = (T)1;
		}		

		masses = massVector.array().matrix().asDiagonal();
	}

	void addParticle(T x, T y, T z) 
	{
		Eigen::Matrix<T, 3, 1> pos(x, y, z);
		positions.push_back(pos);
		++numParticles;
	}

	void addParticle(Eigen::Matrix<T, dim, 1> pos) 
	{
		positions.push_back(pos);
		++numParticles;
	}

	void addParticle(Eigen::Matrix<T, dim, 1> pos, Eigen::Matrix<T, dim, 1> vel) 
	{
		positions.push_back(pos);
		velocities.push_back(vel);
		++numParticles;
	}

	void printPositions()
	{
		for(int i = 0; i < numParticles; i++) {
			std::cout << "Particle " << i << ": <" << positions.at(i)[0] << ", " <<  positions.at(i)[1] << ", " << positions.at(i)[2] << ">" << std::endl; 
		}
	}

	void printVelocities()
	{
		for(int i = 0; i < numParticles; i++) {
			std::cout << "Particle " << i << ": <" << velocities.at(i)[0] << ", " <<  velocities.at(i)[1] << ", " << velocities.at(i)[2] << ">" << std::endl; 
		}	
	}

	void printMasses()
	{
		for(int i = 0; i < numParticles; i++) {
			std::cout << "Particle " << i << ":  " << masses.diagonal()[i] << std::endl;
		}
	}

	void printProperties()
	{
		for(int i = 0; i < numParticles; i++) {
			std::cout << "Particle " << i << std::endl;
			std::cout << "Pos:  " << positions.at(i)[0] << ", " <<  positions.at(i)[1] << ", " << positions.at(i)[2] << ">" << std::endl; 
			std::cout << "Vel:  " << velocities.at(i)[0] << ", " <<  velocities.at(i)[1] << ", " << velocities.at(i)[2] << ">" << std::endl; 
			std::cout << "Mass: " << masses.diagonal()[i] << std::endl;
			std::cout << "-------------------------" << std::endl;
		}
	}
};


