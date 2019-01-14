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

double lennardJones(double distance, double x) {

	/*if (distance)
	{
		//printf("d %f\n", distance);
		double e = 3.4e-10;
		double s = 1.65e-21;
		double a = pow(s / distance, 6);
		double b = (24 * e*a) / distance;
		//printf("res %f\n", 2 * b*a - b);
		return 2 * b*a - b;//b*(2*a-1);
	}
	else {
		return 0;
	}*/
	double distDx = x / (distance*distance);
	//double e = 1.65e-21;
	double mult = 3.4e-10*distDx;
	double mult6 = mult*mult*mult*mult*mult*mult;
	return (7.92e-20*mult6*distDx)*(mult6 - 0.5);
		
	//return (-4 * 3.4e-10*x)*(-12 * pow((1.65e-21*x/ (distance*distance)), 12) + 6 * pow((1.65e-21 *x / (distance*distance)), 6)) / (distance*distance);


}

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
	maxV = 0.0;
	minDx = std::numeric_limits<double>::max();

	double force[3*NumberOfBodies] = { 0 };//ax,ay,az,...


	for (int k = 0; k < NumberOfBodies; k++)
	{
		int aPosK = k * 3;
		double dists[6];
		
		for (int i = k+1; i < NumberOfBodies; i++) {
			int aPosI = i * 3;


			dists[0] = x[i][0] - x[k][0];//x
			dists[1] = x[i][1] - x[k][1];//y
			dists[2] = x[i][2] - x[k][2];//z

			const double distance = sqrt(
				dists[0] * dists[0] +
				dists[1] * dists[1] +
				dists[2] * dists[2]
			);


			dists[3] = lennardJones(distance, dists[0]);
			dists[4] = lennardJones(distance, dists[1]);
			dists[5] = lennardJones(distance, dists[2]);

			for (l = 0; l < 3; ++l)
			{
				force[aPosK+l] += dists[3+l];
				force[aPosI + l] -= dists[3+l];
			}

			minDx = std::min(minDx, distance);
		}

		x[k][0] += timeStepSize * v[k][0];
		x[k][1] += timeStepSize * v[k][1];
		x[k][2] += timeStepSize * v[k][2];

		v[k][0] += timeStepSize * force[aPosK] / mass[k];
		v[k][1] += timeStepSize * force[aPosK + 1] / mass[k];
		v[k][2] += timeStepSize * force[aPosK + 2] / mass[k];

		maxV = std::max(maxV, std::sqrt(v[k][0] * v[k][0] + v[k][1] * v[k][1] + v[k][2] * v[k][2]));
		//printf("reached\n");
	}
	/*for (int i = 0; i < 3 * NumberOfBodies; i++)
	{
		if (i%3==0 && i>0)
		{
			printf("\n");
			printf("%f,", acceleration[i]);
		}
		else {
			printf("%f,", acceleration[i]);
		}
		
	}
	printf("done\n");*/

	t += timeStepSize;
}

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

	openParaviewVideoFile(); //test

	int snapshotCounter = 0;
	if (t > tPlot) {
		printParaviewSnapshot(); //test
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
			printParaviewSnapshot(); //test
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
	printf("size: %f\t nsec: %09zu\n", tFinal / tPlotDelta, (end.tv_sec - start.tv_sec) * 1000000000 + end.tv_nsec - start.tv_nsec);
	closeParaviewVideoFile(); //test

	return 0;
}
