// calcPSF main program

// include libraries from include path
#include <iostream>
#include <unistd.h>
#include <cstring>
#include <cstdlib>

// include libraries from current directory
#include "variables.hpp"
#include "setup.hpp"
#include "load_input.hpp"
#include "validation_Thomas_Fermi.hpp"
#include "calculation_Thomas_Fermi.hpp"
#include "validation_atomic_wfn.hpp"
#include "calculation_atomic_wfn.hpp"
#include "validation_scf_atom.hpp"
#include "calculation_scf_atom.hpp"
#include "validation_phase_shift.hpp"
#include "calculation_phase_shift.hpp"
#include "validation_pad.hpp"
#include "calculation_pad.hpp"

#include "log.hpp"

//namespace
using namespace std;

void report_received_electric_field(){
	char message[512];
	const char* labels[3]={"Ex", "Ey", "Ez"};
	cout << "----SPADExp received electric field----" << endl;
	if(Log_file_set){
		fprintf(Log_file_obj, "----SPADExp received electric field----\n");
	}
	for(int component=0; component<3; component++){
		snprintf(message, sizeof(message),
			"%s = %+.16e %+.16e*i (Cartesian)", labels[component],
			PA_electric_field[component][0], PA_electric_field[component][1]);
		cout << message << endl;
		if(Log_file_set){
			fprintf(Log_file_obj, "%s\n", message);
		}
	}
	snprintf(message, sizeof(message),
		"Electric_field raw order [Ex.re Ex.im Ey.re Ey.im Ez.re Ez.im] = "
		"%.16e %.16e %.16e %.16e %.16e %.16e",
		PA_electric_field[0][0], PA_electric_field[0][1],
		PA_electric_field[1][0], PA_electric_field[1][1],
		PA_electric_field[2][0], PA_electric_field[2][1]);
	cout << message << endl;
	if(Log_file_set){
		fprintf(Log_file_obj, "%s\n", message);
		fflush(Log_file_obj);
	}
	cout << endl;
}

// main function
int main(int argc, const char** argv){
	cout << "Starting SPADExp program..." << endl;
	int exit_code=EXIT_FAILURE;
	
	// stdin, stdout, stderr are terminal or not (file)
	int stdin_is_terminal=isatty(fileno(stdin));
	int stdout_is_terminal=isatty(fileno(stdout));
	int stderr_is_terminal=isatty(fileno(stderr));

	cout << "----IO information----" << endl;
	cout << "stdin is " << (stdin_is_terminal ? "" : "NOT ") << "a terminal" << endl;
	cout << "stdout is " << (stdout_is_terminal ? "" : "NOT ") << "a terminal" << endl;
	cout << "stderr is " << (stderr_is_terminal ? "" : "NOT ") << "a terminal" << endl;
	cout << endl;

	// stdin should not be a terminal
	if(stdin_is_terminal==1){
		cerr << "SPADExp error: stdin should be an input file" << endl;
		printf("Usage: %s < input_file\n", argv[0]);
		return EXIT_FAILURE;
	}

	// initialize: load default values and constants
	initialize();
	
	// load input
	cout << "----Load input----" << endl;
	if(load_input()!=1){
		goto FINALIZATION;
	}
	cout << endl;

	// open output file
	if(Output_file_set==false){
		cout << "Error: invalid Output_file" << endl;
		goto FINALIZATION;
	}
	Output_file_obj=fopen(Output_file, "w");
	if(Output_file_obj==NULL){
		cout << "Error: could not open the output file" << endl;
		goto FINALIZATION;
	}
	
	// open log file if necessary
	if(Log_file_set){
		Log_file_obj=fopen(Log_file, "w");
		if(Log_file_obj==NULL){
			cout << "Error: could not open the log file" << endl;
			goto FINALIZATION;
		}
	}

	// Report exactly what the parser handed to the PAD calculation.  This is
	// deliberately emitted before validation and before any matrix-element
	// conversion so input problems can be diagnosed independently.
	if(strcmp(Calculation, "PAD")==0){
		report_received_electric_field();
	}

	// perform calculation
	cout << "----Perform calculation----" << endl;
	int validation_result;
	if(strcmp(Calculation, "Thomas-Fermi")==0){
		// Thomas-Fermi potential calculation
		write_log((char*)"----Thomas-Fermi potential calculation----");
		validation_result=validation_Thomas_Fermi();
		if(validation_result==1){
			if(TF_test){
				test_Thomas_Fermi();
			}else{
				calc_Thomas_Fermi();
			}
			exit_code=EXIT_SUCCESS;
		}
	}else if(strcmp(Calculation, "Atomic-wfn")==0){
		// Atomic wavefunction calculation
		write_log((char*)"----Atomic wavefunction calculation----");
		validation_result=validation_atomic_wfn();
		if(validation_result==1){
			sequence_atomic_wfn();
			exit_code=EXIT_SUCCESS;
		}
	}else if(strcmp(Calculation, "SCF-atom")==0){
		// self-consistent field calculation of an atom
		write_log((char*)"----SCF calculation of an atom----");
		validation_result=validation_scf_atom();
		if(validation_result==1){
			scf_calc_atom();
			exit_code=EXIT_SUCCESS;
		}
	}else if(strcmp(Calculation, "Phase-shift")==0){
		// phase shift calculation within an atomic potential
		write_log((char*)"----Phase shift calculation----");
		validation_result=validation_phase_shift();
		if(validation_result==1){
			calc_phase_shift();
			exit_code=EXIT_SUCCESS;
		}
	}else if(strcmp(Calculation, "PAD")==0){
		// Photoemission angular distribution calculation
		write_log((char*)"----Photoemission angular distribution calculation----");
		validation_result=validation_PAD();
		if(validation_result==1){
			calculate_PAD();
			exit_code=EXIT_SUCCESS;
		}
	}else{
		cout << "Error: invalid calculation" << endl;
		goto FINALIZATION;
	}

	goto FINALIZATION;
	
 FINALIZATION:
	if(exit_code==EXIT_FAILURE){
		cerr << "SPADExp failed. Review the messages above and the configured log file." << endl;
	}else{
		cout << "SPADExp completed successfully." << endl;
	}
	// finalize: free memory allocated by "new" in initialize()
	finalize();
	return exit_code;
}
