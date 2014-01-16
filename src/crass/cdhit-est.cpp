// =============================================================================
// CD-HI-EST
// http://cd-hit.org/
// Cluster Database at High Identity (EST version)
// modified from CD-HI
//
// program written by 
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// Modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include "cdhit-common.h"
//over-write some defs in cd-hi.h for est version
#undef MAX_UAA
#define MAX_UAA 4

//over-write some defs in cd-hi-init.h for est version

void setaa_to_na();
void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> & Comp_AAN_idx);
void make_comp_iseq(int len, char *iseq_comp, char *iseq);

Options options;
SequenceDB seq_db;

////////////////////////////////////  MAIN /////////////////////////////////////
int cdHitMain(std::string& inputFile, std::string& outputFile, float identityThreshold, int threads)
{

	options.cluster_thd = 0.95;
	options.NAA = 8;
	options.NAAN = NAA8;
	seq_db.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();
	mat.set_to_na(); //mat.set_gap(-6,-1);
    
    //BEGIN MODIFICATIONS
    options.input = inputFile;
    options.output = outputFile;
    options.cluster_thd = identityThreshold;
    options.useIdentity = true;
    options.isEST = true;
    //options.backupFile = true;
#ifndef NO_OPENMP
    int cpu = omp_get_num_procs();
    if( threads > cpu ){
        options.threads = cpu;
        printf( "Warning: total number of CPUs in the system is %i\n", cpu );
    } 
    if( threads == 0 ){
        options.threads = cpu;
        printf( "total number of CPUs in the system is %i\n", cpu );
    } else {
        options.threads = threads;
    }
    if( options.threads != threads ) printf( "Actual number of CPUs to be used: %i\n\n", threads );
//#else
    //printf( "threads is ignored: multi-threading with OpenMP is NOT enabled!\n" );
#endif
    // END MODIFICATIONS
    
	// ***********************************    parse command line and open file
	//if (options.SetOptions( argc, argv, false, true ) == 0) return 1;
	options.Validate();

	InitNAA( MAX_UAA );
	seq_db.NAAN = NAAN_array[options.NAA];

	if ( options.option_r ) {
		Comp_AAN_idx.resize( seq_db.NAAN );
		make_comp_short_word_index(options.NAA, NAAN_array, Comp_AAN_idx);
	}

	seq_db.Read( inputFile.c_str(), options );
	//cout << "total seq: " << seq_db.sequences.size() << endl;
	seq_db.SortDivide( options );
	seq_db.DoClustering( options );

	seq_db.WriteClusters( inputFile.c_str(), outputFile.c_str(), options );

	// write a backup clstr file in case next step crashes
	seq_db.WriteExtra1D( options );

    CleanUpTempFiles();
	return 0;
}
