#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <string.h>
#include <vector>
#include <queue>
using namespace std;

#include "mpi.h"



int multimerge(int * start[], const int lengths[], const int Number, int newArray[], const int newArrayLength);
int compare_ints(const void *a, const void *b);
bool WRITE_TO_FILE = true;


int main( int argc, char* argv[])
{
   	// processor rank, and total number of processors
    int myid, numprocs;
    
    // for timing used by root processor (#0)
    double startwtime = 0.0, endwtime;

    MPI::Init(argc, argv);

    // Collect information about the number of
    // processors and rank of this processor
    numprocs = MPI::COMM_WORLD.Get_size();
    myid = MPI::COMM_WORLD.Get_rank();

	int randomSeed = 1000;
	int myDataSize = 4000;
	for(int i=0;i<argc;i++)
	{
		// search for special arguments
		if (strcmp(argv[i],"-DS")==0)
		{
			myDataSize = atoi(argv[i+1]); i++;
		}
		else if (strcmp(argv[i],"-SR")==0)
		{
			randomSeed = atoi(argv[i+1]); i++;
		}
	}


	int myData[myDataSize];
	int myDataLengths[numprocs];
	int myDataStarts[numprocs];
	
	// communication buffer used for determination of pivot values
	int pivotbuffer[numprocs*numprocs];
	int pivotbufferSize;

	// Compute the individual lengths of mydata[] to be distributed
	// to the numproc processors.  The last processor gets the "extras".
	for(int i=0;i<numprocs;i++)
	{
		myDataLengths[i] = myDataSize/numprocs;
		myDataStarts[i]= i*myDataSize/numprocs;
	}
	myDataLengths[numprocs-1]+=(myDataSize%numprocs);

    // Root node initializes the testing data, and also starts
    // the timer
    if (myid == 0)
    {
		// set random seed and randomly generate testing data
		srandom(randomSeed);
		for(int index=0; index<myDataSize; index++)
		{
			myData[index] = random()% 900;
			if (myDataSize < 50)
				cout << myData[index] << " ";
		}
		
		cout << endl;
		

		startwtime = MPI::Wtime();
    }


	//  The data is scattered to all processors from the root processor (#0)

	if (myid==0)
	{
		MPI::COMM_WORLD.Scatterv(myData,myDataLengths,myDataStarts,MPI::INT,
				MPI_IN_PLACE,myDataLengths[myid],MPI::INT,0);
	}
	else
	{			
		MPI::COMM_WORLD.Scatterv(myData,myDataLengths,myDataStarts,MPI::INT,
				myData,myDataLengths[myid],MPI::INT,0);
	}
      
	// All processors sort their piece of the data using cstdlib::quicksort
	qsort(myData,myDataLengths[myid], sizeof(int), compare_ints);

	// All processors collect regular samples from sorted list
	
	for(int index=0;index<numprocs;index++)
	{
		pivotbuffer[index]= myData[index*myDataLengths[myid]/numprocs];
	}


	// The root processor gathers all pivot candidates from the processors

	if (myid==0)
	{
		MPI::COMM_WORLD.Gather(MPI_IN_PLACE,numprocs,MPI::INT,
			pivotbuffer,numprocs,MPI::INT,0);
	}
	else
	{
		MPI::COMM_WORLD.Gather(pivotbuffer,numprocs,MPI::INT,
			pivotbuffer,numprocs,MPI::INT,0);
	}

	//  Root processor multimerges the lists together and then selects final pivot values to broadcast

	if (myid == 0)
	{
		// multimerge the numproc sorted lists into one
		int *starts[numprocs];  // array of lists
		int lengths[numprocs];  // array of lengths of lists
		for(int i=0;i<numprocs;i++)
		{
			starts[i]=&pivotbuffer[i*numprocs];
			lengths[i]=numprocs;
		}
		int tempbuffer[numprocs*numprocs];  // merged list
		multimerge(starts,lengths,numprocs,tempbuffer,numprocs*numprocs);

		// regularly select numprocs-1 of pivot candidates to broadcast as partition pivot values for myData

		for(int i=0; i<numprocs-1; i++)
		{
			pivotbuffer[i] = tempbuffer[(i+1)*numprocs];
		}				
	}

	// Root processor broadcasts the partition values
	MPI::COMM_WORLD.Bcast(pivotbuffer,numprocs-1,MPI::INT,0);

	// Partition information for myData[]: 
	// 		index of beginning of ith class is classStart[i],
	//		length of ith class is classLength[i], and
	// 		members of ith class, myData[j], have the property
	//   		pivotbuffer[i-1]<= myData[j] < pivotbuffer[i]
	int classStart[numprocs];
	int classLength[numprocs];
	
	// need for each processor to partition its list using the values of pivotbuffer

	int dataindex=0;
	for(int classindex=0; classindex<numprocs-1; classindex++)
	{
		classStart[classindex] = dataindex;
		classLength[classindex]=0;

		// as long as dataindex refers to data in the current class
		while((dataindex< myDataLengths[myid]) 
			&& (myData[dataindex]<=pivotbuffer[classindex]))
		{
			classLength[classindex]++;
			dataindex++;
		}		
	}
	// set Start and Length for last class
	classStart[numprocs-1] = dataindex;
	classLength[numprocs-1] = myDataLengths[myid] - dataindex;
	
	
	int recvbuffer[myDataSize];    // buffer to hold all members of class i
	int recvLengths[numprocs];     // on myid, lengths of each myid^th class
	int recvStarts[numprocs];      // indices of where to start the store from 0, 1, ...

	// processor iprocessor functions as the root and gathers from the
	// other processors all of its sorted values in the iprocessor^th class.  
	for(int iprocessor=0; iprocessor<numprocs; iprocessor++)
	{	
		// Each processor, iprocessor gathers up the numproc lengths of the sorted
		// values in the iprocessor class
		MPI::COMM_WORLD.Gather(&classLength[iprocessor], 1, MPI::INT, 
			recvLengths,1,MPI::INT,iprocessor);
	

		// From these lengths the myid^th class starts are computed on
		// processor myid
		if (myid == iprocessor)
		{
			recvStarts[0]=0;
			for(int i=1;i<numprocs; i++)
			{
				recvStarts[i] = recvStarts[i-1]+recvLengths[i-1];
			}
		}

		// each iprocessor gathers up all the members of the iprocessor^th 
		// classes from the other nodes
		MPI::COMM_WORLD.Gatherv(&myData[classStart[iprocessor]],
			classLength[iprocessor],MPI::INT,
			recvbuffer,recvLengths,recvStarts,MPI::INT,iprocessor);
	}
		
	
	// multimerge these numproc lists on each processor
	int *mmStarts[numprocs]; // array of list starts
	for(int i=0;i<numprocs;i++)
	{
		mmStarts[i]=recvbuffer+recvStarts[i];
	}
	multimerge(mmStarts,recvLengths,numprocs,myData,myDataSize);
	
	int mysendLength = recvStarts[numprocs-1] + recvLengths[numprocs-1];


	int sendLengths[numprocs]; // lengths of consolidated classes
	int sendStarts[numprocs];  // starting points of classes
	// Root processor gathers up the lengths of all the data to be gathered
	MPI::COMM_WORLD.Gather(&mysendLength,1,MPI::INT,
		sendLengths,1,MPI::INT,0);

	// The root processor compute starts from lengths of classes to gather
	if (myid == 0)
	{
		sendStarts[0]=0;
		for(int i=1; i<numprocs; i++)
		{
			sendStarts[i] = sendStarts[i-1]+sendLengths[i-1];
		}	
	}

	// Now we let processor #0 gather the pieces and glue them together in
	// the right order
	int sortedData[myDataSize];
	MPI::COMM_WORLD.Gatherv(myData,mysendLength,MPI::INT,
		sortedData,sendLengths,sendStarts,MPI::INT,0);

	// the root processor prints the elapsed clock time
    if (myid == 0)
	{
		endwtime = MPI::Wtime();
        //cout << "wall clock time (seconds) = " 
		cout <<  setprecision(4) << endwtime-startwtime << endl;

		// cout << "Data set sorted." << endl;

		if (myDataSize < 50)
		{
			for(int index=0; index<myDataSize; index++)
			{
				cout << myData[index] << " ";	   
			}  
		cout << endl;
		}

		if (WRITE_TO_FILE) {
        	FILE *f = fopen("time.csv", "a");
			if (f == NULL)
			{
			    printf("Error opening file!\n");
			    exit(1);
			}
			fprintf(f, "%d,%d,%.4f \n", numprocs, myDataSize, endwtime-startwtime);
			fclose(f);
        }

	}
        
    MPI::Finalize();
    return 0;
}


// ________________________________________________________

struct mmdata 
{
	int stindex;
	int index;
	int stvalue;

	mmdata(int st=0, int id=0, int stv = 0):stindex(st),index(id),stvalue(stv){}

};


int multimerge(int * starts[], const int lengths[], const int Number, 
			   int newArray[], const int newArrayLength)
{
 	// Create priority queue.  There will be at most one item in the priority queue
 	// for each of the Number lists.
 	priority_queue< mmdata> priorities;

 	// Examine each of the Number start[] lists, place the first location into 
	// the priority 	queue if the list is not empty
 	for(int i=0; i<Number;i++)
 	{
		if (lengths[i]>0)
		{
			priorities.push(mmdata(i,0,starts[i][0]));
		}
	}


	// As long as priorities is not empty, pull off the top member (the smallest 
	//value from list i), push it into the newArray, and place the next element from 
	// list i in the priority queue
	int newArrayindex = 0;  // index into the merged array
	while (!priorities.empty() && (newArrayindex<newArrayLength))
	{
		// grab the smallest element, and remove it from the priority queue
		mmdata xxx = priorities.top();
		priorities.pop();

		// insert this smallest element into the merged array
		newArray[newArrayindex++] = starts[xxx.stindex][xxx.index];

		// if start[xxx.stindex] is not empty, place the next member into priority
		if ( lengths[xxx.stindex]>(xxx.index+1))
		{
			priorities.push(mmdata(xxx.stindex, xxx.index+1, 
								starts[xxx.stindex][xxx.index+1]));
		}
}

// return the logical size of the merged array
return newArrayindex;
}

bool operator<( const mmdata & One, const mmdata & Two)
{
	return One.stvalue > Two.stvalue;
}

int compare_ints(const void *a, const void *b)
{
	int myint1 = *reinterpret_cast<const int *>(a);
	int myint2 = *reinterpret_cast<const int *>(b);
	if (myint1<myint2) return -1;
	if (myint1>myint2) return 1;
	return 0;
}
