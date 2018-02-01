/*  
 * This file handles the file IO of the solver.
*/

#include <Partio.h>

#include "particles.h"
#include "segmentMesh.h"
#include "triangleMesh.h"

#define NUMFRAMES 120
#define OUTPUT_DIRECTORY "output_files/"

using T = double;
constexpr int dim = 3;

// Writes to .bgeo or .pda files
template <class T, int dim>
void writePartio(Particles<T, dim> particles, const std::string& particleFile)
{
	Partio::ParticlesDataMutable* parts = Partio::create();
	Partio::ParticleAttribute posH, vH, mH;
	mH = parts->addAttribute("m", Partio::VECTOR, 1);
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);

	for(int i = 0; i < particles.numParticles; i++) {
		// Get the index of the particle we're working on
		int idx = parts->addParticle();
 		
		// Get the address of the mass, position, velocity
		float* m = parts->dataWrite<float>(mH, idx);
		float* p = parts->dataWrite<float>(posH, idx);
		float* v = parts->dataWrite<float>(vH, idx);
	
		// Set mass
		// m[0] = particles.masses.at(i);
		m[0] = particles.masses.diagonal()[i];
		// Set position and velocity
		for(int j = 0; j < 3; j++) {
			p[j] = (T)particles.positions.at(i)[j];
			v[j] = (T)particles.velocities.at(i)[j];
		}		
	}

	Partio::write(particleFile.c_str(), *parts);
	parts->release();
}

void generatePointCloud(const std::string &name, const std::string &ext, int numParticles)
{
	// Create particles 
	Particles<T, dim> particles(numParticles);

	// Create a bgeo file for each frame
	for(int frame = 0; frame < NUMFRAMES; frame++) {
		std::string file = name + std::to_string(frame + 1) + ext;
        writePartio<T,dim>(particles, file);
	}
}

void generateSegmentMesh(const std::string &fileName, const std::string &ext, int numParticles)
{
	SegmentMesh<T, dim> curve(numParticles);

	for(int frame = 0; frame < NUMFRAMES; frame++) {
		std::string f = fileName + std::to_string(frame + 1) + ext;
		curve.particles.updatePositions(frame);
		
		std::ofstream file(f);
		if(file.is_open()) {
			file << "POINTS\n";
			for(int i = 0; i < curve.numParticles; i++) {
				std::string point = std::to_string(i + 1) + ": ";
		
				for(int j = 0; j < 3; j++) {
					point += std::to_string(curve.particles.positions.at(i)[j]) + " ";
				}
		
				file << point;
				file << "\n";
			}	

			file << "POLYS\n";
			for(int i = 0; (size_t)i < curve.indices.size(); i++) {
				file << std::to_string(i + 1) + ": ";
				
				for(int j = 0; (size_t)j < curve.indices.at(i).size(); j++) {
					file << curve.indices.at(i).at(j) << " ";
				}

				file << "\n";
			}

			file << "END\n";
	
			file.close();
		} 
		else {
			std::cout << "Couldn't open " << fileName << std::endl;
		}

		std::string bgeo = "segment_mesh/frame" + std::to_string(frame + 1) + ".bgeo"; 
		writePartio<T, dim>(curve.particles, bgeo);
	}
}

void generateTriangleMesh(const std::string &fileName, const std::string &ext, int numParticles)
{
	TriangleMesh<T, dim> mesh(numParticles);
	Particles<T, dim> particles(numParticles);
	mesh.particles = particles;

	std::ofstream file("3d_mesh/triangleMesh.obj");
	if(file.is_open()) {
		for(int i = 0; i < mesh.numParticles; i++) {
			// Eigen::Matrix<T, 3, 1> vertex = mesh.particles.positions.at(i);
			Eigen::Matrix<T, 3, 1> vertex = mesh.particles.positions.at(i);

			// Check for duplicates
			for(int j = 0; j < mesh.numParticles; j++) {
				if(i == j) {
					continue;
				}

				if(vertex == mesh.particles.positions.at(j)) {
					std::cout << "DUPLICATE" << std::endl;
				}
			}

			// std::cout << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;	

			file << "v";
			for(int j = 0; j < 3; j++) {
				file << " " << vertex[j];
			}	
			file << "\n";
		}

		file << "\n";
		
		for(int i = 0; (size_t)i < mesh.indices.size(); i++) {
			int count = 0;
			for(int j = 0; (size_t)j < mesh.indices.at(i).size(); j++) {
				++count;

				if(j % 3 == 0 && count == 1) {
					file << "f ";
				}

				file << mesh.indices.at(i).at(j);
			
				if(j % 3 != 2) {
					file << " ";
				}
				else {
					file << " ";
					if(count == 4) {
						file << "\n";
					}
					else {
						file << " ";
					}
				}
			} 
			file << "\n";
		} 
	
		file.close();
	}
	else {
		std::cout << "Couldn't open " << fileName << std::endl;
	}	

	for(int frame = 0; frame < NUMFRAMES; frame++) {
		mesh.particles.updatePositions((T).1);

		std::string bgeo = "3d_mesh/frame" + std::to_string(frame + 1) + ".bgeo";
		writePartio(mesh.particles, bgeo);
	}
}

T computeLength(Eigen::Matrix<T, 3, 1> &p1, Eigen::Matrix<T, 3, 1> &p2)
{
	float xSqr = (p1[0] - p2[0]) * (p1[0] - p2[0]);
	float ySqr = (p1[1] - p2[1]) * (p1[1] - p2[1]);
	float zSqr = (p1[2] - p2[2]) * (p1[2] - p2[2]);

	float length = sqrt(xSqr + ySqr + zSqr);

	return length;
}

void parse(const std::string &filename, const std::string &output)
{
	// Read file
	std::ifstream file(filename);
	// Open file for writing
	std::ofstream outputFile(output);

	if(file.is_open() && outputFile.is_open()) {
		std::string line = "";
		// int testcase = 0;
		// Loop through each line in the file
		while(std::getline(file, line)) {
			// std::cout << "Testcase " << ++testcase << std::endl;

            // Create a new segment for each line
			SegmentMesh<T, dim> segment;
			T k, gamma, L;

			bool youngsMod 	= false;
			bool dampCoeff 	= false;
			bool restLength	= false;

			int count = 0;
			T x, y, z;

			char c;					// Use this to iterate through each character 
			std::string num = "";	// Concatenate c if not a space and convert to number after
			std::vector<Eigen::Matrix<T, 3, 1>> vectors;	// Stores all the positions and velocities

			// Iterate through the line and extract numbers
			for(size_t i = 0; i < line.length(); i++) {
				c = line.at(i);

				// Found a number
				if(c == ' ') {
					T n = (T)std::stof(num);
					num = "";

					// First three numbers are the spring components
					if(!youngsMod) {
						k = n;					// Set Young's Modulus
						youngsMod = true;		// Flag to get the next component
						continue;
					} 
					else if(!dampCoeff) {
						gamma = n;				// Set Damping Coefficient
						dampCoeff = true;		// Flag to get the next component
						continue;
					}
					else if(!restLength) {
						L = n;					// Set Rest Length
						restLength = true;		// Flag to start getting positions
						continue;
					}

					++count;	

					if(count % 3 == 1) {
						x = n;					// Set x component
					} 
					else if(count % 3 == 2) {
						y = n;					// Set y component
					}
					else if(count % 3 == 0) {
						z = n;					// Set z component

						vectors.push_back(Eigen::Matrix<T, 3, 1>(x, y, z));
						count = 0;				// Keep numbers low for debugging				
					}
				} 
				// Concat c to string until we find a space to convert to number
				else {
					num += c;
				}
			}

			// Extract positions and velocities from vectors
			for(size_t i = 0; i < vectors.size(); i++) {
				// Position and Velocity should align at halfway point
				// n = 8
				// Indices ----
				// Pos: 0 1 2 3 
				// Vel: 4 5 6 7 
				if(i < vectors.size() / 2) {
					Eigen::Matrix<T, 3, 1> pos = vectors.at(i);
					Eigen::Matrix<T, 3, 1> vel = vectors.at(i + vectors.size() / 2);
					segment.particles.addParticle(pos, vel);
				}
				// Already processed last particle
				else {
					break;
				}
			}

			// Update mass matrix
			segment.particles.updateMasses();

			// Set spring properties
			segment.k 		= k;
			segment.gamma 	= gamma;
			segment.restLen = L;

			// segment.particles.printProperties();
			// segment.printYoungsMod();
			// segment.printDampingCoeff();
			// segment.printRestLength();

			// TODO
			// segment number of particles is wrong

			Eigen::Matrix<T, 3, 1> p0 = segment.particles.positions.at(0);
			Eigen::Matrix<T, 3, 1> p1 = segment.particles.positions.at(1);
			Eigen::Matrix<T, 3, 1> p2 = segment.particles.positions.at(2);
			Eigen::Matrix<T, 3, 1> p3 = segment.particles.positions.at(3);

			Eigen::Matrix<T, 3, 1> v0 = segment.particles.velocities.at(0);
			Eigen::Matrix<T, 3, 1> v1 = segment.particles.velocities.at(1);
			Eigen::Matrix<T, 3, 1> v2 = segment.particles.velocities.at(2);
			Eigen::Matrix<T, 3, 1> v3 = segment.particles.velocities.at(3);

			// Calculate f0
			float l01 = (p0 - p1).norm();//computeLength(p0, p1);
			float l03 = (p0 - p3).norm();//computeLength(p0, p3);

			Eigen::Matrix<T, 3, 1> n01 = (p0 - p1) / l01;
			Eigen::Matrix<T, 3, 1> n03 = (p0 - p3) / l03;

			Eigen::Matrix<T, 3, 1> f01 = -segment.k * ((l01 / segment.restLen) - 1) * n01;
			Eigen::Matrix<T, 3, 1> f03 = -segment.k * ((l03 / segment.restLen) - 1) * n03;

			Eigen::Matrix<T, 3, 1> f0 = f01 + f03;
			segment.particles.forces.push_back(f0);
			// std::cout << "f0\n" << f0 << std::endl;

			// Calculate fd0
			Eigen::Matrix<T, 3, 1> fd01 = -segment.gamma * (n01 * n01.transpose()) * (v0 - v1);
			Eigen::Matrix<T, 3, 1> fd03 = -segment.gamma * (n03 * n03.transpose()) * (v0 - v3);

			Eigen::Matrix<T, 3, 1> fd0 = fd01 + fd03;
			segment.particles.dampForces.push_back(fd0);
			// std::cout << "fd0\n" << fd0 << std::endl;

			// Calculate f1 
			float l10 = (p1 - p0).norm();//computeLength(p1, p0);
			float l13 = (p1 - p3).norm();//computeLength(p1, p3);
			float l12 = (p1 - p2).norm();//computeLength(p1, p2);

			Eigen::Matrix<T, 3, 1> n10 = (p1 - p0) / l10;
			Eigen::Matrix<T, 3, 1> n13 = (p1 - p3) / l13;
			Eigen::Matrix<T, 3, 1> n12 = (p1 - p2) / l12;
			
			Eigen::Matrix<T, 3, 1> f10 = -segment.k * ((l10 / segment.restLen) - 1) * n10;
			Eigen::Matrix<T, 3, 1> f13 = -segment.k * ((l13 / segment.restLen) - 1) * n13;
			Eigen::Matrix<T, 3, 1> f12 = -segment.k * ((l12 / segment.restLen) - 1) * n12;

			Eigen::Matrix<T, 3, 1> f1 = f10 + f13 + f12;
			segment.particles.forces.push_back(f1);
			// std::cout << "f1\n" << f1 << std::endl;

			// Calculate fd1
			Eigen::Matrix<T, 3, 1> fd10 = -segment.gamma * (n10 * n10.transpose()) * (v1 - v0);
			Eigen::Matrix<T, 3, 1> fd13 = -segment.gamma * (n13 * n13.transpose()) * (v1 - v3);
			Eigen::Matrix<T, 3, 1> fd12 = -segment.gamma * (n12 * n12.transpose()) * (v1 - v2);

			Eigen::Matrix<T, 3, 1> fd1 = fd10 + fd13 + fd12;
			segment.particles.dampForces.push_back(fd1);
			// std::cout << "fd1\n" << fd1 << std::endl;

			// Calculate f2
			float l21 = (p2 - p1).norm();//computeLength(p2, p1);
			float l23 = (p2 - p3).norm();//computeLength(p2, p3);

			Eigen::Matrix<T, 3, 1> n21 = (p2 - p1) / l21;
			Eigen::Matrix<T, 3, 1> n23 = (p2 - p3) / l23;

			Eigen::Matrix<T, 3, 1> f21 = -segment.k * ((l21 / segment.restLen) - 1) * n21;
			Eigen::Matrix<T, 3, 1> f23 = -segment.k * ((l23 / segment.restLen) - 1) * n23;

			Eigen::Matrix<T, 3, 1> f2 = f21 + f23;
			segment.particles.forces.push_back(f2);
			// std::cout << "f2\n" << f2 << std::endl;

			// Calculate fd2
			Eigen::Matrix<T, 3, 1> fd21 = -segment.gamma * (n21 * n21.transpose()) * (v2 - v1);
			Eigen::Matrix<T, 3, 1> fd23 = -segment.gamma * (n23 * n23.transpose()) * (v2 - v3);

			Eigen::Matrix<T, 3, 1> fd2 = fd21 + fd23;
			segment.particles.dampForces.push_back(fd2);
			// std::cout << "fd2\n" << fd2 << std::endl;

			// Calculate f3
			float l30 = (p3 - p0).norm();//computeLength(p3, p0);
			float l32 = (p3 - p2).norm();//computeLength(p3, p2);
			float l31 = (p3 - p1).norm();//computeLength(p3, p1);

			Eigen::Matrix<T, 3, 1> n30 = (p3 - p0) / l30;
			Eigen::Matrix<T, 3, 1> n31 = (p3 - p1) / l31;
			Eigen::Matrix<T, 3, 1> n32 = (p3 - p2) / l32;
			
			Eigen::Matrix<T, 3, 1> f30 = -segment.k * ((l30 / segment.restLen) - 1) * n30;
			Eigen::Matrix<T, 3, 1> f31 = -segment.k * ((l31 / segment.restLen) - 1) * n31;
			Eigen::Matrix<T, 3, 1> f32 = -segment.k * ((l32 / segment.restLen) - 1) * n32;

			Eigen::Matrix<T, 3, 1> f3 = f30 + f31 + f32;
			segment.particles.forces.push_back(f3);
			// std::cout << "f3\n" << f3 << std::endl;

			// Calculate fd1
			Eigen::Matrix<T, 3, 1> fd30 = -segment.gamma * (n30 * n30.transpose()) * (v3 - v0);
			Eigen::Matrix<T, 3, 1> fd31 = -segment.gamma * (n31 * n31.transpose()) * (v3 - v1);
			Eigen::Matrix<T, 3, 1> fd32 = -segment.gamma * (n32 * n32.transpose()) * (v3 - v2);

			Eigen::Matrix<T, 3, 1> fd3 = fd30 + fd31 + fd32;
			segment.particles.dampForces.push_back(fd3);
			// std::cout << "fd3\n" << fd3 << std::endl;

			outputFile << f0[0] << " " << f0[1] << " " << f0[2] << " ";
			outputFile << f1[0] << " " << f1[1] << " " << f1[2] << " ";
			outputFile << f2[0] << " " << f2[1] << " " << f2[2] << " ";
			outputFile << f3[0] << " " << f3[1] << " " << f3[2] << " ";
			outputFile << fd0[0] << " " << fd0[1] << " " << fd0[2] << " ";
			outputFile << fd1[0] << " " << fd1[1] << " " << fd1[2] << " ";
			outputFile << fd2[0] << " " << fd2[1] << " " << fd2[2] << " ";
			outputFile << fd3[0] << " " << fd3[1] << " " << fd3[2] << "\n";

			// std::cout << "--------------" << std::endl;
		}

		// std::cout << "DONE" << std::endl;
	}
	else {
		std::cout << "Couldn't read " << filename << std::endl;
		return;
	}
}

// Solve Ax = b
// void solveLinearSystem()
// {
// 	Eigen::VectorXd x = Eigen::VectorXd(3 * segment.particles.numParticles);
// 	for(int i = 0; i < segment.particles.numParticles; i++) {
// 		if(	i >= 3 * segment.particles.numParticles || 
// 			i + 1 >= 3 * segment.particles.numParticles || 
// 			i + 2 >= 3 * segment.particles.numParticles) 
// 		{
// 			break;
// 		}

// 		Eigen::VectorXd part = segment.particles.positions.at(i);
// 		// std::cout << part[0] << " " << part[1] << " " << part[2] << std::endl;

// 		x[3 * i] 		= part[0];	// x
// 		x[3 * i + 1] 	= part[1];	// y
// 		x[3 * i + 2] 	= part[2];	// z
// 	}

// 	Eigen::VectorXd b = segment.particles.masses * x;
// 	std::cout << b << std::endl;
// }