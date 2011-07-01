//
//  SmithWaterman.cpp
//  crass
#include "SmithWaterman.h"


//////////////////////////////////////////////
// Simple Ends-Free Smith-Waterman Algorithm
//
// You will be prompted for input sequences
// Penalties and match scores are hard-coded
//
// Program does not perform multiple tracebacks if 
// it finds several alignments with the same score
//
// By Nikhil Gopal
// Similar implementation here: https://wiki.uni-koeln.de/biologicalphysics/index.php/Implementation_of_the_Smith-Waterman_local_alignment_algorithm
//////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <map>
using namespace std;

int ind; 
double penalty=-4;



double similarityScore(char a, char b)
{
	double result;
	if(a==b)
	{
		result=1;
	}
	else
	{
		result=penalty;
	}
	return result;
}

double findMax(double array[], int length)
{
	double max = array[0];
	ind = 0;
    
	for(int i=1; i<length; i++)
	{
		if(array[i] > max)
		{
			max = array[i];
			ind=i;
		}
	}
	return max;
}

stringPair smithWaterman(std::string seqA, std::string seqB )
{

	// initialize some variables
	int lengthSeqA = (int)seqA.length();
	int lengthSeqB = seqB.length();
	
	// initialize matrix
	double matrix[lengthSeqA+1][lengthSeqB+1];
	for(int i=0;i<=lengthSeqA;i++)
	{
		for(int j=0;j<=lengthSeqB;j++)
		{
			matrix[i][j]=0;
		}
	}
    
	double traceback[4];
	int I_i[lengthSeqA+1][lengthSeqB+1];
	int I_j[lengthSeqA+1][lengthSeqB+1];
    
    
	//start populating matrix
	double matrix_max = 0;
    int i_max = 0, j_max = 0;
	for (int i=1;i<=lengthSeqA;i++)
	{
		for(int j=0;j<=lengthSeqB;j++)
        {
			traceback[0] = matrix[i-1][j-1]+similarityScore(seqA[i-1],seqB[j-1]);
			traceback[1] = matrix[i-1][j]+penalty;
			traceback[2] = matrix[i][j-1]+penalty;
			traceback[3] = 0;
			matrix[i][j] = findMax(traceback,4);
			if(matrix[i][j]> matrix_max)
			{
				matrix_max = matrix[i][j];
				i_max = i;
				j_max = j;
			}
			switch(ind)
			{
				case 0:
					I_i[i][j] = i-1;
					I_j[i][j] = j-1;
					break;
				case 1:
					I_i[i][j] = i-1;
                    I_j[i][j] = j;
                    break;
				case 2:
					I_i[i][j] = i;
                    I_j[i][j] = j-1;
                    break;
				case 3:
					I_i[i][j] = i;
                    I_j[i][j] = j;
                    break;
			}
        }
	}
    
	
	int current_i = i_max,current_j = j_max;
	
    int next_i = I_i[current_i][current_j];
	int next_j = I_j[current_i][current_j];
	
    int tick = 0;
    
	while((next_j!=0) && (next_i!=0) && ((current_i!=next_i) || (current_j!=next_j)))
	{

        
		current_i = next_i;
		current_j = next_j;
		next_i = I_i[current_i][current_j];
		next_j = I_j[current_i][current_j];
		tick++;
	}
    
    std::pair<std::string, std::string> ret_pair(seqA.substr(current_i, i_max - current_i), seqB.substr(current_j, j_max - current_j));
    
    return ret_pair;
}
