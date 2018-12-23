// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2018.c -o assignment
//
// Run it with
//
// ./assignment (linux) or assigment.exe .. var
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
//
// installed mingw, g++ in D/. for anything: open mingw compiler go to ./NAPP by c\ run the above command.
// this is a project under repo of vs. the other is under NAPP
// to add arg: go to properties -> debug
//
// before submit anything, test with linux, the given code etc. read assignment. MinGW for G++ which linux already has, ParaView for modeling
//
// (C) 2018 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>

#include <algorithm>//added

//global data
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
 * Set up scenario from the command line. scenario parser: timeInc, maxtime, X-xyz V:xyz M
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
	NumberOfBodies = (argc - 2) / 7;

	x = new double*[NumberOfBodies];
	v = new double*[NumberOfBodies];
	mass = new double[NumberOfBodies];

	int readArgument = 1;

	//timeStepSize = std::stof(argv[readArgument]);
	tPlotDelta = std::stof(argv[readArgument]); readArgument++;

	//printf("setup tplotdelta %f \n",tPlotDelta);
	tFinal = std::stof(argv[readArgument]); readArgument++;
	//printf("setup tfinal %f \n", tFinal);

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
 * opens paraview file
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
	videoFile.open("result.pvd");
	videoFile << "<?xml version=\"1.0\"?>" << std::endl
		<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
		<< "<Collection>";
}





/**
 * closes paraview file
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
	videoFile << "</Collection>"
		<< "</VTKFile>" << std::endl;
}


/**
 * paraview instance maker with xml
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
		<< "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

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

double lennardJones(double distance) {

	if (distance)
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
	}
	
}

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
	//declares vars
	maxV = 0.0;
	minDx = std::numeric_limits<double>::max();

	//declares 3 way forces, x,y,z
	double force[1];
	force[0] = 0.0;
	force[1] = 0.0;
	force[2] = 0.0;

	//for each body, calcs distance to all and updates 3way force with formula.
	/*Working area*/

	for (int k = 0; k < NumberOfBodies; k++) {

		for (int i = 0; i < NumberOfBodies && i != k; i++) {
			const double distance = sqrt(
				(x[k][0] - x[i][0]) * (x[k][0] - x[i][0]) +
				(x[k][1] - x[i][1]) * (x[k][1] - x[i][1]) +
				(x[k][2] - x[i][2]) * (x[k][2] - x[i][2])
			);

			//printf("d0 %f\n", x[i][0] - x[k][0]);
			force[0] += lennardJones(x[i][0] - x[k][0]); //* mass[i] * mass[k] / distance / distance / distance;
			//printf("res0 %f\n", lennardJones(x[i][0] - x[k][0]));
			//printf("d1 %f\n", x[i][1] - x[k][1]);
			force[1] += lennardJones(x[i][1] - x[k][1]); //* mass[i] * mass[k] / distance / distance / distance;
			//printf("res1 %f\n", lennardJones(x[i][1] - x[k][1]));
			//printf("d2 %f\n", x[i][2] - x[k][2]);
			force[2] += lennardJones(x[i][2] - x[k][2]); //* mass[i] * mass[k] / distance / distance / distance;
			//printf("res2 %f\n", lennardJones(x[i][2] - x[k][2]));

			minDx = std::min(minDx, distance);

		}

		//updates dist with time*v aka x=V*t
		x[k][0] = x[k][0] + tPlotDelta * v[k][0];//timeStepSize
		x[k][1] = x[k][1] + tPlotDelta * v[k][1];
		x[k][2] = x[k][2] + tPlotDelta * v[k][2];

		//updates velo with time-force aka v=F*t/m  = ma*t/m = x/t^2*t
		v[k][0] = v[k][0] + tPlotDelta * force[0] / mass[k];
		v[k][1] = v[k][1] + tPlotDelta * force[1] / mass[k];
		v[k][2] = v[k][2] + tPlotDelta * force[2] / mass[k];

		// takes sqrt of all way v's to get directional v.
		
		//maxV = v[k][0];

	}
	maxV = std::sqrt(v[k][0] * v[k][0] + v[k][1] * v[k][1] + v[k][2] * v[k][2]);
	//inc of time 
	//printf("isnide t: %f += tPlotDelta: %f \n", t, tPlotDelta);
	t += tPlotDelta;//timeStepSize;
}


/**
 * Main routine.
 * result.pnv maker
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

	openParaviewVideoFile();

	int snapshotCounter = 0;
	if (t > tPlot) {
		printParaviewSnapshot();
		std::cout << "plotted initial setup" << std::endl;
		tPlot = tPlotDelta;
	}
	
	int timeStepCounter = 0;
	while (t <= tFinal) {
		updateBody();
		timeStepCounter++;
		if (t >= tPlot) {
			printParaviewSnapshot();
			std::cout << "plot next snapshot"
				<< ",\t time step=" << timeStepCounter
				<< ",\t t=" << t
				<< ",\t dt=" << timeStepSize
				<< ",\t v_max=" << maxV
				<< ",\t dx_min=" << minDx
				<< std::endl;

			tPlot += tPlotDelta;
			//printf("tPlot += tPlotDelta -> %f+=%f,,,,,timeStepSize %f \n", tPlot,tPlotDelta,timeStepSize);
		}
	}

	closeParaviewVideoFile();

	return 0;

}
