// SPADExp tools for OpenMX
// preproc_main

// include libraries from include path
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstring>
#define _USE_MATH_DEFINES
#include <cmath>

// include libraries from current directory
#include "inputTools.hpp"
#include "kfold.hpp"

// namespace
using namespace std;

void outProd(double* a, double* b, double* c);
double inProd(double* a, double* b);
double calcKp(double p, double q, double r_X, double r_Y, double r_O, double XY,
						double* origin_frac, double* x_frac, double* y_frac, bool curved, double* kp);
double orth_threshold=0.001;

int main(int argc, const char** argv){
	cout << "Starting preproc tool for OpenMX..." << endl;
	bool use_first_bz_only=true;
	const char* input_path=nullptr;
	const char* output_path=nullptr;
	if(argc==3){
		input_path=argv[1];
		output_path=argv[2];
	}else if(argc==4 && strcmp(argv[1], "--no-fold")==0){
		use_first_bz_only=false;
		input_path=argv[2];
		output_path=argv[3];
	}else{
		printf("Usage: %s [--no-fold] input_file output_file\n", argv[0]);
		return 1;
	}
	printf("K-point mode: %s\n", use_first_bz_only
		? "first-BZ calculation with folding"
		: "full requested-zone calculation without folding");

	ifstream input(input_path);
	if(!input){
		printf("Error: cannot open %s\n", input_path);
		return 1;
	}
	load_file(&input);
	input.close();

	// variables
	int i, j;
	bool curved;
	int dimension;
	double origin_frac[3]={0};
	double x_frac[3]={0};
	double y_frac[3]={0};
	double origin_orth[3];
	double x_orth[3];
	double y_orth[3];
	double xMin, xMax;
	double yMin, yMax;
	int x_count;
	int y_count;
	double cell[3][3];
	spadexp::kfold::ReciprocalCell rec_cell;

	if(load_int("SPADExp.dimension", &dimension)==1){
		printf("SPADExp.dimension is %d\n", dimension);
		if(!(1<=dimension && dimension<=2)){
			printf("Error: invalid dimension\n");
			return 0;
		}	
	}else{
		return 0;
	}
	
	if(load_logical("SPADExp.curved", &curved)==1){
		printf("SPADExp.curved is %s\n", (curved ? "true": "false"));
	}else{
		return 0;
	}

	if(load_doublev("SPADExp.origin", 3, origin_frac)==1){
		// printf("Origin is (%.3f, %.3f, %.3f) in fractional coordinates\n", origin_frac[0], origin_frac[1], origin_frac[2]);
	}else{
		return 0;
	}

	int line_range_start=find_str("<SPADExp.range");
	int line_range_end=find_str("SPADExp.range>");

	if(line_range_start>=0 && line_range_end>=0){
		if(line_range_end-line_range_start-1==dimension){
			char* range=get_line(line_range_start+1);
			if(sscanf(range, "%lf %lf %lf %lf %lf %d", &x_frac[0], &x_frac[1], &x_frac[2], &xMin, &xMax, &x_count)!=6){
				printf("Error in parsing x range");
				return 0;
			}
			if(x_count<1){
				printf("Error: x_count must be positive\n");
				return 0;
			}
			range=get_line(line_range_start+2);
			if(dimension>1){
				if(sscanf(range, "%lf %lf %lf %lf %lf %d", &y_frac[0], &y_frac[1], &y_frac[2], &yMin, &yMax, &y_count)!=6){
					printf("Error in parsing y range");
					return 0;
				}
				if(y_count<1){
					printf("Error: y_count must be positive\n");
					return 0;
				}
			}
		}else{
			printf("Dimension error\n");
			return 0;
		}
	}else{
			return 0;
	}

	// printf("(%.3f, %.3f, %.3f) [%.3f, %.3f] %d\n", x_frac[0], x_frac[1], x_frac[2], xMin, xMax, x_count);
	// printf("(%.3f, %.3f, %.3f) [%.3f, %.3f] %d\n", y_frac[0], y_frac[1], y_frac[2], yMin, yMax, y_count);

	int line_bandCell_start=find_str("<Band.KPath.UnitCell");
	int line_bandCell_end=find_str("Band.KPath.UnitCell>");

	int line_cell_start=find_str("<Atoms.UnitVectors");
	int line_cell_end=find_str("Atoms.UnitVectors>");

	if(line_bandCell_start>=0 && line_bandCell_end>=0){
		line_cell_start=line_bandCell_start;
		line_cell_end=line_bandCell_end;
	}

	if(line_cell_end-line_cell_start-1!=3){
		printf("Cell error");
		return 0;
	}

	for(i=0; i<3; i++){
		char* cell_str=get_line(line_cell_start+1+i);
		if(sscanf(cell_str, "%lf %lf %lf", &cell[i][0], &cell[i][1], &cell[i][2])!=3){
			printf("Error in parsing cell");
			return 0;
		}
	}

	printf("Unit vectors (real space)\n");
	for(i=0; i<3; i++){
	 	printf("a_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), cell[i][0], cell[i][1], cell[i][2]);
	}

	// calculate the reciprocal unit cell
	double vBuf[3];
	outProd(cell[1], cell[2], vBuf);
	double det=inProd(cell[0], vBuf);

	outProd(cell[1], cell[2], vBuf);
	for(i=0; i<3; i++){
		rec_cell[0][i]=2*M_PI*vBuf[i]/det;
	}
	outProd(cell[2], cell[0], vBuf);
	for(i=0; i<3; i++){
		rec_cell[1][i]=2*M_PI*vBuf[i]/det;
	}
	outProd(cell[0], cell[1], vBuf);
	for(i=0; i<3; i++){
		rec_cell[2][i]=2*M_PI*vBuf[i]/det;
	}
	
	printf("Unit vectors (reciprocal space)\n");
	for(i=0; i<3; i++){
	 	printf("b_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), rec_cell[i][0], rec_cell[i][1], rec_cell[i][2]);
	}
	printf("\n");

	for(i=0; i<3; i++){
		origin_orth[i]=0;
		x_orth[i]=0;
		y_orth[i]=0;
		for(j=0; j<3; j++){
			origin_orth[i]+=rec_cell[j][i]*origin_frac[j];
			x_orth[i]+=rec_cell[j][i]*x_frac[j];
			y_orth[i]+=rec_cell[j][i]*y_frac[j];
		}
	}

	double r_O=sqrt(inProd(origin_orth, origin_orth));
	double r_X=sqrt(inProd(x_orth, x_orth));
	double r_Y=sqrt(inProd(y_orth, y_orth));

	printf("O = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", origin_orth[0], origin_orth[1], origin_orth[2], r_O);
	printf("X = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", x_orth[0], x_orth[1], x_orth[2], r_X);
	printf("Y = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", y_orth[0], y_orth[1], y_orth[2], r_Y);
	double OX=inProd(origin_orth, x_orth);
	double XY=inProd(x_orth, y_orth);
	double YO=inProd(y_orth, origin_orth);
	if(abs(OX)>orth_threshold){
		printf("Error: O and X are not orthogonal\n");
		return 0;
	}
	if(abs(YO)>orth_threshold){
		printf("Error: Y and O are not orthogonal\n");
		return 0;
	}

	if(curved){
		printf("Curved mesh accepted (X and Y may be non-orthogonal)\n\n");
	}else{
		printf("Planar mesh accepted (X and Y may follow non-orthogonal reciprocal vectors)\n\n");
	}

	// K=pX+qY+rO
	// r=1 if not curved, else so that |K| = r_O
	double p, q;
	double dp, dq;
	double kp[3];
	double r_Kp;
	if(x_count<2){
		dp=0;
	}else{
		dp=(xMax-xMin)/(x_count-1);
	}
	if(y_count<2){
		dq=0;
	}else{
		dq=(yMax-yMin)/(y_count-1);
	}
	vector<array<double, 3> > requested_kpoints;
	vector<array<double, 3> > reduced_kpoints;
	vector<array<double, 3> > unique_kpoints;
	vector<array<int, 3> > reciprocal_shifts;
	vector<int> unique_index;
	const double unique_tolerance=1.0e-9;
	auto append_kpoint = [&](double p_value, double q_value){
		array<double, 3> requested;
		array<double, 3> reduced;
		array<int, 3> shift;
		r_Kp=calcKp(p_value, q_value, r_X, r_Y, r_O, XY,
						 origin_frac, x_frac, y_frac, curved, requested.data());
		int found=-1;
		if(use_first_bz_only){
			const spadexp::kfold::FoldResult folded_point =
				spadexp::kfold::fold_to_first_bz(requested, rec_cell);
			reduced=folded_point.reduced;
			shift=folded_point.shift;
			for(size_t iu=0; iu<unique_kpoints.size(); iu++){
				double distance2=0.0;
				for(int axis=0; axis<3; axis++){
					double delta=reduced[axis]-unique_kpoints[iu][axis];
					distance2+=delta*delta;
				}
				if(distance2<unique_tolerance*unique_tolerance){
					found=(int)iu;
					break;
				}
			}
			if(found<0){
				found=(int)unique_kpoints.size();
				unique_kpoints.push_back(reduced);
			}
		}else{
			reduced=requested;
			shift={0, 0, 0};
			found=(int)unique_kpoints.size();
			unique_kpoints.push_back(requested);
		}
		requested_kpoints.push_back(requested);
		reduced_kpoints.push_back(reduced);
		reciprocal_shifts.push_back(shift);
		unique_index.push_back(found);
		if(use_first_bz_only){
			printf("%11.6f %11.6f %11.6f -> %11.6f %11.6f %11.6f  G=(%d,%d,%d) # r=%11.5f\n",
				requested[0], requested[1], requested[2], reduced[0], reduced[1], reduced[2],
				shift[0], shift[1], shift[2], r_Kp);
		}else{
			printf("%11.6f %11.6f %11.6f (direct) # r=%11.5f\n",
				requested[0], requested[1], requested[2], r_Kp);
		}
	};

	cout << (use_first_bz_only
		? "requested -> first-BZ kpoints in fractional coordinates"
		: "requested kpoints in fractional coordinates (no folding)") << endl;
	if(dimension==2){
		for(j=0; j<y_count; j++){
			for(i=0; i<x_count; i++){
				p=xMin+dp*i;
				q=yMin+dq*j;
				append_kpoint(p, q);
			}
		}
	}else{
		q=0;
		for(i=0; i<x_count; i++){
			p=xMin+dp*i;
			append_kpoint(p, q);
		}
	}

	// prepare the output file
	// neglect MO.fileout, MO.Nkpoint, and <MO.kpoint> block
	ofstream output(output_path);
	copyInput(&output);
	output << endl;
	output << "#" << endl;
	output << (use_first_bz_only
		? "# MO output folded into the first Brillouin zone"
		: "# MO output at every requested k point; no Brillouin-zone folding") << endl;
	output << "# written by openmx_preproc" << endl;
	output << "#" << endl << endl;
	output << "MO.fileout on" << endl;
	output << "MO.Nkpoint " << unique_kpoints.size() << endl;
	output << "<MO.kpoint" << endl;
	output.precision(12);
	output.setf(ios::fixed, ios::floatfield);
	for(const auto& point: unique_kpoints){
		output << point[0] << " " << point[1] << " " << point[2] << endl;
	}
	output << "MO.kpoint>" << endl;

	if(use_first_bz_only){
		// OpenMX ignores this SPADExp-specific block.  The postprocessor uses it
		// to expand the unique first-BZ eigenvectors back to the requested grid.
		output << "SPADExp.KFold on" << endl;
		output << "SPADExp.KFold.RequestedCount " << requested_kpoints.size() << endl;
		output << "SPADExp.KFold.UniqueCount " << unique_kpoints.size() << endl;
		output << "<SPADExp.KFold.Map" << endl;
		for(size_t ir=0; ir<requested_kpoints.size(); ir++){
			output << requested_kpoints[ir][0] << " " << requested_kpoints[ir][1] << " "
			       << requested_kpoints[ir][2] << " " << reduced_kpoints[ir][0] << " "
			       << reduced_kpoints[ir][1] << " " << reduced_kpoints[ir][2] << " "
			       << reciprocal_shifts[ir][0] << " " << reciprocal_shifts[ir][1] << " "
			       << reciprocal_shifts[ir][2] << " " << unique_index[ir] << endl;
		}
		output << "SPADExp.KFold.Map>" << endl;
		printf("First-BZ folding: %zu requested k points -> %zu unique OpenMX k points\n",
			requested_kpoints.size(), unique_kpoints.size());
	}else{
		printf("No folding: %zu requested k points -> %zu direct OpenMX k points\n",
			requested_kpoints.size(), unique_kpoints.size());
	}

	cout << "Preparation finished." << endl;
	
}


void outProd(double* a, double* b, double* c){
	// c=a\times b
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
}

double inProd(double* a, double* b){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double calcKp(double p, double q, double r_X, double r_Y, double r_O, double XY,
						double* origin_frac, double* x_frac, double* y_frac, bool curved, double* kp){
	double r;
	if(!curved){
		r=1.0;
	}else{
		// |pX+qY|^2+(r|O|)^2=|O|^2.  Keep the X.Y cross term so
		// reciprocal-periodic b1/b2 meshes work for non-orthogonal cells.
		double transverse2=pow(p*r_X, 2)+pow(q*r_Y, 2)+2.0*p*q*XY;
		double radial2=pow(r_O, 2)-transverse2;
		// Points outside the propagating sphere are retained at kz=0 so the
		// rectangular data shape is preserved; the GUI masks them physically.
		r=(radial2>0.0 && r_O>0.0) ? sqrt(radial2)/r_O : 0.0;
	}
	double r_Kp=sqrt(pow(p*r_X, 2)+pow(q*r_Y, 2)+2.0*p*q*XY
					+pow(r*r_O, 2));
	int i;
	for(i=0; i<3; i++){
		kp[i]=p*x_frac[i]+q*y_frac[i]+r*origin_frac[i];
	}
	return r_Kp;
	
}
