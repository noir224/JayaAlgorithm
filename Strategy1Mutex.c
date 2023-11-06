// Online C compiler to run C program online

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <string.h>

#define MATRIX_D 100
#define MATRIX_P 100
#define MINR_VALUE 0
#define MAXR_VALUE 1
#define MIN_VALUE -30.0
#define MAX_VALUE 30.0
#define N 100000
#define K 5
#define KROWS 50

pthread_mutex_t best_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t worst_mutex = PTHREAD_MUTEX_INITIALIZER;

float global_best_fit;
float global_worst_fit;
float global_best_x[MATRIX_D];
float global_worst_x[MATRIX_D];

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

    float best_f;

    float worst_f;

    float best_x[MATRIX_D];

    float worst_x[MATRIX_D];

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

void *fitness(float threadmatrix[KROWS][MATRIX_D], float fit[KROWS])
{

    float row[MATRIX_D];

    for (int i = 0; i < KROWS; i++)
    {

        for (int j = 0; j < MATRIX_D; j++)
        {

            row[j] = threadmatrix[i][j];
        }

        float rv = Rosenbrock(row, MATRIX_D);

        fit[i] = rv;
    }

    // return fit;
}

///////////////////////////// Jayas////////////////////////////////////

void *Jaya1(JayaParams *jp)
{

    // calculate fitness

    fitness(jp->threadmatrix, jp->fit); 
    
    // get the best row and worst row, best fit worst fit values of this thread

    float best_fit = jp->fit[0];

    int best_index = 0;

    float worst_fit = jp->fit[0];

    int worst_index = 0;

    for (int i = 0; i < KROWS; i++)
    {

        float fi = (jp->fit[i]);

        if (fi < best_fit)
        {

            best_fit = fi;

            best_index = i;
        }

        if ((fi) > worst_fit)
        {

            worst_fit = fi;

            worst_index = i;
        }
    }

    (jp->rv)->best_f = best_fit; //Best fit for this thread

    (jp->rv)->worst_f = worst_fit; //Worst fit for this thread


    // get best row and worst row and save them in struct

    float best_row[MATRIX_D]; // best solution

    float worst_row[MATRIX_D]; // worst solution

    float thbe = 0;

    float thb = 0;

    for (int i = 0; i < MATRIX_D; i++)
    {

        thbe = jp->threadmatrix[best_index][i];

        best_row[i] = thbe;

        (jp->rv)->best_x[i] = thbe;

    }


    for (int i = 0; i < MATRIX_D; i++)
    {

        thb = jp->threadmatrix[worst_index][i];

        worst_row[i] = thb;

        (jp->rv)->worst_x[i] = thb;

       
    }

    //change the global fitnesses and their corresponding xs 
    pthread_mutex_lock(&best_mutex);

    if (global_best_fit >= best_fit)
    {

        global_best_fit = best_fit;

        for (int i = 0; i < MATRIX_D; i++)

            global_best_x[i] = best_row[i];
    }

    pthread_mutex_unlock(&best_mutex);

    pthread_mutex_lock(&worst_mutex);

    if (global_worst_fit <= worst_fit)
    {

        global_worst_fit = worst_fit;

        for (int i = 0; i < MATRIX_D; i++)

            global_worst_x[i] = worst_row[i];
    }

    pthread_mutex_unlock(&worst_mutex);

    pthread_exit(NULL);

}

void *Jaya2(JayaParams *jp)
{


    float newMatrix[KROWS][MATRIX_D]; // x'

    float newFit[KROWS]; // f'

    float r1[MATRIX_D];

    float r2[MATRIX_D];

    float xi[MATRIX_D]; // current row

    // generate new population

  

    // X'

    for (int j = 0; j < KROWS; j++)
    {

        for (int k = 0; k < MATRIX_D; k++)
        {
         xi[k] = jp->threadmatrix[j][k];
	r1[k] = ((float)rand() / RAND_MAX) * (MAXR_VALUE - MINR_VALUE) + MINR_VALUE;

        r2[k] = ((float)rand() / RAND_MAX) * (MAXR_VALUE - MINR_VALUE) + MINR_VALUE;
            newMatrix[j][k] = xi[k] + (r1[k] * (global_best_x[k] - xi[k])) - (r2[k] * (global_worst_x[k] - xi[k]));
        }

      
    }


    // calculate new fitness F'

    fitness(newMatrix, newFit);
    //New Population X

    for (int i = 0; i < KROWS; i++)
    {


        if (jp->fit[i] >= newFit[i])
        {

            for (int j = 0; j < MATRIX_D; j++)
            {

                jp->threadmatrix[i][j] = newMatrix[i][j];

            }
        }
    }

    pthread_exit(NULL);
}

///////////////////////////// Global comparisons /////////////////////////////////



int main()
{

    pthread_t threads[K]; // array of threads
    float matrix[MATRIX_P][MATRIX_D]; // original matrix

    srand(time(NULL));
    // Generate random values for the matrix (first population taken from last stage)
   /*FILE *fp;
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

*/
		    for (int i = 0; i < MATRIX_P; i++)
    {

        //printf("\n");

        for (int j = 0; j < MATRIX_D; j++)
        {

            matrix[i][j] = ((float)rand() / RAND_MAX) * (MAX_VALUE - MIN_VALUE) + MIN_VALUE;

           // printf("%d %f  ", i, matrix[i][j]);
        }
    }

    int startRows[K]; //indecies of starting rows for each thread

    int endRows[K]; //indecies of ending rows for each thread

    // to divide the matrix to parts for threads

    int ind = 0;

    for (int i = 0; i < MATRIX_P; i += KROWS)
    {

        startRows[ind] = i;

        endRows[ind] = i + KROWS;

        ind++;
    }

    // passed parameters to jaya for each thread
    JayaParams *jp[K];

    for (int i = 0; i < K; i++)
    {

        jp[i] = malloc(sizeof(JayaParams));

        jp[i]->threadmatrix = malloc(KROWS * sizeof(float[MATRIX_D]));

        jp[i]->fit = malloc(KROWS * sizeof(float));

        jp[i]->rv = malloc(sizeof(ReturnValues)); // the return values from the thread for Jaya 1

        jp[i]->i = i;
    }

    
    // filling the part of matrix that each matrix should take
    for (int rowi = 0; rowi < K; rowi++)
    {

        float threadmatrix[KROWS][MATRIX_D];

        int arri = 0;

        for (int j = startRows[rowi]; j < endRows[rowi]; j++)
        {

            for (int h = 0; h < MATRIX_D; h++)
            {

                threadmatrix[arri][h] = matrix[j][h];

                jp[rowi]->threadmatrix[arri][h] = threadmatrix[arri][h];
            }

            arri++;
        }

    } // rowi


    clock_t start_time = clock();
	//creating threads
    for (int n = 0; n < N; n++)
    {

        global_best_fit = INFINITY;

        global_worst_fit = -INFINITY;

        for (int rowi = 0; rowi < K; rowi++)

            pthread_create(&threads[rowi], NULL, (void *)Jaya1, jp[rowi]);

        // waiting for all threads to come back from jaya 1

        for (int i = 0; i < K; i++)
        {
            pthread_join(threads[i], NULL);

        } // i

	// Jaya 2 threads

        for (int i = 0; i < K; i++)
        {

            pthread_create(&threads[i], NULL, (void *)Jaya2, jp[i]);

        } // k

        for (int i = 0; i < K; i++)
        {

            pthread_join(threads[i], NULL);

        } // i

        // final update on original population

        for (int k = 0; k < K; k++)
        {

            int sr = startRows[k];

            int er = endRows[k];

            int arrind1 = 0;

            for (int i = sr; i < er; i++)
            {
                for (int j = 0; j < MATRIX_D; j++)
                {

                    matrix[i][j] = jp[k]->threadmatrix[arrind1][j];

                } // j

                arrind1++;

            } // i

        } // k

    } // N

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
    
    for (int i = 0; i < K; i++) {
	    free(jp[i]->rv);
	    free(jp[i]->fit);
	    free(jp[i]->threadmatrix);
	    free(jp[i]);
     }


    return 0;
}

