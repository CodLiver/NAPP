
// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2018.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>

#include <time.h>//test
#include <omp.h>


double t = 0;
double tFinal = 0;
double tPlot = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0001;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
	NumberOfBodies = (argc - 2) / 7;

	x = new double*[NumberOfBodies];
	v = new double*[NumberOfBodies];
	mass = new double[NumberOfBodies];

	int readArgument = 1;

	tPlotDelta = std::stof(argv[readArgument]); readArgument++;
	tFinal = std::stof(argv[readArgument]); readArgument++;

	for (int i = 0; i < NumberOfBodies; i++) {
		x[i] = new double[3];
		v[i] = new double[3];

		x[i][0] = std::stof(argv[readArgument]); readArgument++;
		x[i][1] = std::stof(argv[readArgument]); readArgument++;
		x[i][2] = std::stof(argv[readArgument]); readArgument++;

		v[i][0] = std::stof(argv[readArgument]); readArgument++;
		v[i][1] = std::stof(argv[readArgument]); readArgument++;
		v[i][2] = std::stof(argv[readArgument]); readArgument++;

		mass[i] = std::stof(argv[readArgument]); readArgument++;

		if (mass[i] <= 0.0) {
			std::cerr << "invalid mass for body " << i << std::endl;
			exit(-2);
		}
	}

	std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

	if (tPlotDelta <= 0.0) {
		std::cout << "plotting switched off" << std::endl;
		tPlot = tFinal + 1.0;
	}
	else {
		std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
		tPlot = 0.0;
	}
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
	videoFile.open("result.pvd");
	videoFile << "<?xml version=\"1.0\"?>" << std::endl
		<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
		<< "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
	videoFile << "</Collection>"
		<< "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
	static int counter = -1;
	counter++;
	std::stringstream filename;
	filename << "result-" << counter << ".vtp";
	std::ofstream out(filename.str().c_str());
	out << "<VTKFile type=\"PolyData\" >" << std::endl
		<< "<PolyData>" << std::endl
		<< " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
		<< "  <Points>" << std::endl
		<< "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";

	for (int i = 0; i < NumberOfBodies; i++) {
		out << x[i][0]
			<< " "
			<< x[i][1]
			<< " "
			<< x[i][2]
			<< " ";
	}

	out << "   </DataArray>" << std::endl
		<< "  </Points>" << std::endl
		<< " </Piece>" << std::endl
		<< "</PolyData>" << std::endl
		<< "</VTKFile>" << std::endl;

	videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
 void updateBody() {
 	maxV = 0.0;
 	minDx = std::numeric_limits<double>::max();

 	double acceleration[3*NumberOfBodies] = { 0 };//ax,ay,az,...

  	//lessons I learnt during the imp of my first Deep L. n.

 	#pragma omp parallel for
 	for (int k = 0; k < NumberOfBodies; k++)
 	{
 		int aPosK = k * 3;
 		for (int i = k+1; i < NumberOfBodies; i++) {
 			int aPosI = i * 3;

 			double dists0 = x[i][0] - x[k][0];//x
 			double dists1 = x[i][1] - x[k][1];//y
 			double dists2 = x[i][2] - x[k][2];//z

 			double distance = sqrt(
 				dists0 * dists0 +
 				dists1 * dists1 +
 				dists2 * dists2
 			);

 			double divver = distance * distance * distance;
 			double mIdiv = mass[i] / divver;
 			double mKdiv = mass[k] / divver;

 			acceleration[aPosK] += dists0 * mIdiv;
 			acceleration[aPosK + 1] += dists1 * mIdiv;
 			acceleration[aPosK + 2] += dists2 * mIdiv;

 			acceleration[aPosI] -= dists0 * mKdiv;
 			acceleration[aPosI + 1] -= dists1 * mKdiv;
 			acceleration[aPosI + 2] -= dists2 * mKdiv;

 			minDx = std::min(minDx, distance);
 		}

 	}
  	#pragma omp parallel for
 	for (int k = 0; k < NumberOfBodies; k++)
 	{
	    int aPosK = k * 3;
	    x[k][0] += timeStepSize * v[k][0];
	    x[k][1] += timeStepSize * v[k][1];
	    x[k][2] += timeStepSize * v[k][2];

	    v[k][0] += timeStepSize * acceleration[aPosK];
	    v[k][1] += timeStepSize * acceleration[aPosK + 1];
	    v[k][2] += timeStepSize * acceleration[aPosK + 2];

	  maxV = std::max(maxV,std::sqrt(v[k][0] * v[k][0] + v[k][1] * v[k][1] + v[k][2] * v[k][2]));
	  }

 	t += timeStepSize;
 }

//./assignment1 0.01 100.0  7.0 4.0 10.0 4.0 4.0 4.0 10.0 7.0 3.0 11.0 5.0 2.0 2.0 6.0 9.0 15.0 1.0 4.0 3.0 5.0 9.0 14.0 11.0 9.0 1.0 3.0 3.0 7.0 10.0 8.0 3.0 1.0 4.0 2.0 6.0 8.0 11.0 10.0 0.0 3.0 4.0 1.0 11.0 14.0 5.0 3.0 0.0 1.0 5.0 14.0 9.0 15.0 3.0 4.0 4.0 1.0 9.0 11.0 3.0 2.0 4.0 0.0 10.0 13.0 12.0 3.0 1.0 0.0 5.0 8.0
/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
	if (argc == 1) {
		std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time objects" << std::endl
			<< "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
			<< "  final-time      simulated time (greater 0)" << std::endl
			<< std::endl
			<< "Examples:" << std::endl
			<< "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
			<< "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
			<< "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
			<< std::endl
			<< "In this naive code, only the first body moves" << std::endl;

		return -1;
	}
	else if ((argc - 3) % 7 != 0) {
		std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
		return -2;
	}

	setUp(argc, argv);

	//openParaviewVideoFile(); //test

	int snapshotCounter = 0;
	if (t > tPlot) {
		//printParaviewSnapshot(); //test
		std::cout << "plotted initial setup" << std::endl;
		tPlot = tPlotDelta;
	}

	struct timespec start, end;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	int timeStepCounter = 0;
	while (t <= tFinal) {
		updateBody();
		timeStepCounter++;
		if (t >= tPlot) {
			// printf("cur step %d\n",timeStepCounter);
			//printParaviewSnapshot(); //test
			/*std::cout << "plot next snapshot"
				<< ",\t time step=" << timeStepCounter
				<< ",\t t=" << t
				<< ",\t dt=" << timeStepSize
				<< ",\t v_max=" << maxV
				<< ",\t dx_min=" << minDx
				<< std::endl;*/

			tPlot += tPlotDelta;
		}
	}

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

	printf("size: %f\t nsec: %09zu\n", tFinal/tPlotDelta, (end.tv_sec - start.tv_sec) * 1000000000 + end.tv_nsec - start.tv_nsec);
	//closeParaviewVideoFile(); //test

	return 0;
}
