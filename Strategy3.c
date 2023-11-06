// Online C compiler to run C program online

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <string.h>

#define MATRIX_D 100
#define MATRIX_P 250
#define MINR_VALUE 0
#define MAXR_VALUE 1
#define MIN_VALUE -30.0
#define MAX_VALUE 30.0
#define N 10 //cannot be less than K
#define M 5000
#define K 10
#define KROWS 25

float global_best_fit;
float global_worst_fit;
float global_best_x[MATRIX_D];
float global_worst_x[MATRIX_D];

int counter=0;

int threadsnum= K;
int mergingsize= 1;

pthread_mutex_t counter_mutex = PTHREAD_MUTEX_INITIALIZER;


///////////////////////// STRUCTURES ///////////////////////////////////

typedef struct ReturnValues
{

    float best_f;

    float worst_f;

    float best_x[MATRIX_D];

    float worst_x[MATRIX_D];

} ReturnValues;

typedef struct JayaParams
{

    float (*threadmatrix)[MATRIX_D];

    float *fit;
     
    int rowsnum;
    
    ReturnValues *rv;
    
    int i;

} JayaParams;

/////////////////////// Supporting functions///////////////////////////////////

// Rastrigin function

float rastrigin(float x[MATRIX_D], int D) {

float sum = 0.0;

	for (int j = 0; j < D; j++) {

	sum += (pow(x[j], 2) - 10 * cos(2 * M_PI * (x[j]*(M_PI/180.0))) + 10);

	}

	return sum;

}



float sphere(float arr[], int D){

        float sum = 0;

    for (int i = 0; i < D; i++) {

        sum += arr[i] * arr[i];

    }

    return sum;



}



float Rosenbrock(float arr[], int D)
{

    float sum = 0.0;

    for (int j = 0; j < D - 1; j++)
    {

        float term1 = 100 * pow((arr[j + 1] - pow(arr[j], 2)), 2);

        float term2 = pow(1 - arr[j], 2);

        sum += term1 + term2;
    }

    return sum;
}

void *fitness(float threadmatrix[][MATRIX_D], float fit[], int rnum)
{
   float row[MATRIX_D];

    for (int i = 0; i < rnum; i++)
    {

        for (int j = 0; j < MATRIX_D; j++)
        {

            row[j] = threadmatrix[i][j];
        }

        float rv = sphere(row, MATRIX_D);

        fit[i] = rv;
    }

}

///////////////////////////// Jayas////////////////////////////////////



void *Jaya1(JayaParams *jp)
{
	int rnum = jp->rowsnum;

    // calculate fitness

    fitness(jp->threadmatrix, jp->fit,rnum);


    // get the best row and worst row, best fit worst fit values of this thread

    float best_fit = jp->fit[0];

    int best_index = 0;

    float worst_fit = jp->fit[0];

    int worst_index = 0;

    for (int i = 0; i < rnum; i++)
    {

        float fi = (jp->fit[i]);

        if (fi <= best_fit)
        {

            best_fit = fi;

            best_index = i;
        }

        if ((fi) >= worst_fit)
        {

            worst_fit = fi;

            worst_index = i;
        }
    }

    (jp->rv)->best_f = best_fit;

    (jp->rv)->worst_f = worst_fit;

    //get best row and worst row and save them in struct

    float thbe = 0;

    float thb = 0;

    for (int i = 0; i < MATRIX_D; i++)
    {

        thbe = jp->threadmatrix[best_index][i];

        (jp->rv)->best_x[i] = thbe;

       
    }

 

    for (int i = 0; i < MATRIX_D; i++)
    {

        thb = jp->threadmatrix[worst_index][i];

        (jp->rv)->worst_x[i] = thb;

        
    }

   

    


}

void *Jaya2(JayaParams *jp)
{


	int rnum = jp->rowsnum;
    float newMatrix[rnum][MATRIX_D]; // x'

    float newFit[rnum]; // f'

    float r1[MATRIX_D];

    float r2[MATRIX_D];

    float xi[MATRIX_D]; // current row



    //X' 

    for (int j = 0; j < rnum; j++)
    {

        

        for (int k = 0; k < MATRIX_D; k++)
        {
        xi[k] = jp->threadmatrix[j][k];
        	
        r1[k] = ((float)rand() / RAND_MAX) * (MAXR_VALUE - MINR_VALUE) + MINR_VALUE;

        r2[k] = ((float)rand() / RAND_MAX) * (MAXR_VALUE - MINR_VALUE) + MINR_VALUE;

            newMatrix[j][k] = xi[k] + (r1[k] * ((jp->rv)->best_x[k] - xi[k])) - (r2[k] * ((jp->rv)->worst_x[k] - xi[k]));

        }

    }



    // calculate new fitness

    fitness(newMatrix, newFit, rnum);
    	
    //update the submatrix of this thread
    for (int i = 0; i < rnum; i++)
    {



        if (jp->fit[i] >= newFit[i])
        {

            for (int j = 0; j < MATRIX_D; j++)
            {

                jp->threadmatrix[i][j] = newMatrix[i][j];

            }
        }
    }

    
    pthread_mutex_lock(&counter_mutex);
    counter++;
    pthread_mutex_unlock(&counter_mutex);

    
}

void *fullJaya(JayaParams *jp){
	for(int n=0; n<M;n++){
	pthread_mutex_lock(&counter_mutex);
    	if((counter+1) < M){
    		pthread_mutex_unlock(&counter_mutex);
		Jaya1(jp);
		Jaya2(jp);
		}
	else{
	   pthread_mutex_unlock(&counter_mutex);
	   break;
	}
	}

	pthread_exit(NULL);

}

///////////////////////////// Global comparisons /////////////////////////////////

void *compareBestAndWorst(float bestfs[], float bestxs[][MATRIX_D], float worstfs[], float worstxs[][MATRIX_D])
{
    float bestf = bestfs[0];

    int bestfi = 0;

    float bestx[MATRIX_D];

    // find the minimum value of bestfs and its corresponding row from bestxs

    for (int i = 0; i < threadsnum; i++)
    {

        if (bestfs[i] <= bestf)
        {

            bestf = bestfs[i];

            for (int j = 0; j < MATRIX_D; j++)
            {

                bestx[j] = bestxs[i][j];
            }
        }
    }

    float worstf = worstfs[0];

    float worstx[MATRIX_D];

    // find the maximum value of worstfs and its corresponding row from worstxs

    for (int i = 0; i < threadsnum; i++)
    {

        if (worstfs[i] >= worstf)
        {

            worstf = worstfs[i];

            for (int j = 0; j < MATRIX_D; j++)
            {

                worstx[j] = worstxs[i][j];
            }
        }
    }

    // changing global values
    global_best_fit = bestf;

    for (int j = 0; j < MATRIX_D; j++)
    {

        global_best_x[j] = bestx[j];

    }



    global_worst_fit = worstf;

    for (int j = 0; j < MATRIX_D; j++)
    {

        global_worst_x[j] = worstx[j];
    }

    
}

////////////////////////////// Partition function/////////////////////////////////////

void* partition (int *startRows, int *endRows){
		
	int* rowsNumForEachThread = malloc(threadsnum * sizeof(int)); //array for the numbers of rows for each threads (since it keeps changing)
	for(int i=0;i<threadsnum;i++){
		if(i==0)
	           rowsNumForEachThread[i] = mergingsize * KROWS; //because we're merging the first 2 always so the size of submatrix of first thread is always changing
	         else
	            rowsNumForEachThread[i]= KROWS;
	}
	int ind=0;
	int i=0;
	//deciding the start index and end index of each submatrix
	while(i<MATRIX_P){
	   
	   startRows[ind] = i;
	   endRows[ind] = i + rowsNumForEachThread[ind];
	   
	   i = i + rowsNumForEachThread[ind];
	   ind++;
	
	}
	free(rowsNumForEachThread);
}

int main()
{

    

    float matrix[MATRIX_P][MATRIX_D]; // original matrix

    srand(time(NULL));

    printf("First Population");

    // Generate random values for the matrix (first population
   FILE *fp;
   fp = fopen("input.txt", "r");
    if (fp == NULL) {
        printf("Failed to open file.");
        return 1;
    }

    for (int i = 0; i < MATRIX_P; i++) {
        for (int j = 0; j < MATRIX_D; j++) {
            fscanf(fp, "%f", &matrix[i][j]);
        }
    }

    fclose(fp);
    
   pthread_t threads[K]; // array of threads
   JayaParams *jp[K]; // array of structure
    

	printf("\n");

   clock_t start_time = clock();
    for( int n=0; n<N;n++){
    	int* startRows = malloc(threadsnum * sizeof(int));
        int * endRows = malloc(threadsnum * sizeof(int));
	partition(startRows, endRows);
	int totalNumOfRows = mergingsize * KROWS;
	
	for (int i = 0; i < threadsnum; i++)
    	{

		jp[i] = malloc(sizeof(JayaParams));
		
		if(i==0){ //for the merged thread
			
			jp[i]->threadmatrix = malloc(totalNumOfRows * sizeof(float[MATRIX_D]));
			jp[i]->fit = malloc(totalNumOfRows * sizeof(float));
			jp[i] -> rowsnum = totalNumOfRows;
		}else{
		         jp[i]->threadmatrix = malloc(KROWS * sizeof(float[MATRIX_D]));
		         jp[i]->fit = malloc(KROWS * sizeof(float));
		         jp[i] -> rowsnum = KROWS;
		}

		jp[i]->rv = malloc(sizeof(ReturnValues)); // the return values from the thread for Jaya 
		
		jp[i] -> i= i;
    	}
	
	//dividing the matrix to parts for each thread to take (a part can have different size from the others)
	    printf("\n");

    for (int rowi = 0; rowi < threadsnum; rowi++)
    {
        int arri = 0;

        for (int j = startRows[rowi]; j < endRows[rowi]; j++)
        {

            for (int h = 0; h < MATRIX_D; h++)
            {

                jp[rowi]->threadmatrix[arri][h] = matrix[j][h];
            }

            arri++;
        }
        

    } // rowi
       for (int rowi = 0; rowi < threadsnum; rowi++)

            pthread_create(&threads[rowi], NULL, (void *)fullJaya, jp[rowi]);

        // waiting for all threads to come back from jaya 1

        for (int i = 0; i < K; i++)
        {
            pthread_join(threads[i], NULL);

        } // i
        
        counter=0;
        //getting best/worst fitness of all best/worst fitnesses and its corresponding x
        float bestfs[threadsnum]; // array for the best fs from all threads

        float bestxs[threadsnum][MATRIX_D]; // array for the x row of corresponding best f

        float worstfs[threadsnum]; // array for the worst fs from all threads

        float worstxs[threadsnum][MATRIX_D]; // array for the x row of corresponding worst f

        for (int i = 0; i < threadsnum; i++)
        {

            bestfs[i] = (jp[i]->rv)->best_f;

            worstfs[i] = (jp[i]->rv)->worst_f;

            for (int j = 0; j < MATRIX_D; j++)
            {

                bestxs[i][j] = (jp[i]->rv)->best_x[j];

                worstxs[i][j] = (jp[i]->rv)->worst_x[j];

            } // j

        } // i

        // to place the best/worst of best/worst fitness and its corresponding x in the global variables

        compareBestAndWorst(bestfs, bestxs, worstfs, worstxs);
        
        //update original matrix 
        for (int k = 0; k < threadsnum; k++)
        {

            int sr = startRows[k];

            int er = endRows[k];

            int arrind1 = 0;

            for (int i = sr; i < er; i++)
            {

               // printf("\n");

                for (int j = 0; j < MATRIX_D; j++)
                {

                    matrix[i][j] = jp[k]->threadmatrix[arrind1][j];

                    //printf("%f ", matrix[i][j]);

                } // j

                arrind1++;

            } // i

        } // k
        
 
        
    	printf("/////////////////////////////////\n");
    	
    	 //free mallocs because sizes will change for the next iteration
    	free(startRows);
    	free(endRows);
    	for (int i = 0; i < threadsnum; i++) {
	    free(jp[i]->rv);
	    free(jp[i]->fit);
	    free(jp[i]->threadmatrix);
	    free(jp[i]);
     	}
     	if((threadsnum-1)==0)
    		threadsnum=1;
    	else{
    	     threadsnum--;
    	     mergingsize++;
    	  }

    	
    
    }//N
    clock_t end_time = clock();

    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nExecution time: %f seconds\n", execution_time);

    printf("\n last best f %f\n", global_best_fit);

    printf("last best x\n");

    for (int i = 0; i < MATRIX_D; i++)
    {

        printf("%f ", global_best_x[i]);
    }

    printf("\n last worst f %f\n", global_worst_fit);

    printf("last worst x\n");

    for (int i = 0; i < MATRIX_D; i++)
    {

        printf("%f ", global_worst_x[i]);
    }
    


    return 0;
}

