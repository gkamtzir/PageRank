/**
* FILE: gauss_seidel.c
* THMMY, 7th semester, Parallel and Distributed Systems: 4th assignment
* Imperative Page Rank algorithm using Gauss-Seidel approach without
* compression. This algorithm is the most generic one and it's used as a 
* reference for the validity of the results.
* Author:
*   Kamtziridis Georgios, 8542, gkamtzir@auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int get_total_number_of_nodes();
void initialize_matrix(int number_of_nodes);
void populate_matrix(int number_of_nodes);

void replace_dangling_nodes(int number_of_nodes);
void calculate_A_matrix(int number_of_nodes, double alpha);
double calculate_error(int number_of_nodes);
void gauss_seidel(int number_of_nodes, int max_number_of_iterations);

//Global variables.
double **A, *b, *page_rank, *page_rank_previous;
int *number_of_outbound_links;
double desired_error;
char folder_name[128];

int main(int argc, char **argv) {

	if (argc != 2) {
        printf("You have to provide a number that indicates which dataset to choose\n"
            "1: Gun Control (2955 nodes)\n2: Search Engines (11659 nodes)\n");
        exit(1);
    }

    if (atoi(argv[1]) == 1) {

        sprintf(folder_name, "gun_control");

    } else {

        sprintf(folder_name, "search_engines");

    }

	int i, j, number_of_nodes, max_number_of_iterations;
	double p, alpha;
    
	//Initializing variables.
	number_of_nodes = get_total_number_of_nodes();
	max_number_of_iterations = 50;
	desired_error = 0.000001;
	alpha = 0.85;
	p = 1.0 / number_of_nodes;

	number_of_outbound_links = (int *)malloc(sizeof(int) * number_of_nodes);
	page_rank = (double *)malloc(sizeof(double) * number_of_nodes);
	page_rank_previous = (double *)malloc(sizeof(double) * number_of_nodes);
	b = (double *)malloc(sizeof(double) * number_of_nodes);
	A = (double **)malloc(sizeof(double *) * number_of_nodes);

	if (number_of_outbound_links == NULL || page_rank == NULL || b == NULL 
		|| page_rank_previous == NULL || A == NULL) {
		printf("ERROR: could not allocate memory.\n");
		exit(1);
	}

	for (i = 0; i < number_of_nodes; i++) {

		A[i] = (double *)malloc(sizeof(double) * number_of_nodes);

		if (A[i] == NULL) {
			printf("ERROR: could not allocate memory.\n");
			exit(1);
		}
		
		//Initializing the page rank vectors.
		page_rank[i] = 1.0 / number_of_nodes;
		page_rank_previous[i] = 1.0 / number_of_nodes;

		//Initializing the 'b' vector.
    	b[i] = (1-alpha) * p;
	}

	initialize_matrix(number_of_nodes);
	populate_matrix(number_of_nodes);

	//Replacing the dangling nodes.
	replace_dangling_nodes(number_of_nodes);

	//Calcluating the 'A' matrix.
	calculate_A_matrix(number_of_nodes, alpha);

	// //Page rank calculation.
    gauss_seidel(number_of_nodes, max_number_of_iterations);
    
	//Storing the results in a file.
	FILE *file;

    file = fopen("./page_rank_simple", "w");

    for (i = 0; i < number_of_nodes; i++)
    {
		fprintf(file, "%f\n", page_rank[i]);
    }

	fclose(file);

	//Freeing memory.
	free(number_of_outbound_links);
	free(page_rank);
	free(b);

	free(page_rank_previous);

	for (i = 0; i < number_of_nodes; i++) {
		free(A[i]);
	}
	free(A);

	return 0;

}

//Fetches the total number of nodes from the 'node' file.
int get_total_number_of_nodes() {

	FILE *file_nodes;
	int number_of_nodes;
	char filename[128];

    sprintf(filename, "./dataset/%s/nodes", folder_name);

	file_nodes = fopen(filename,"r");

	if (file_nodes == NULL) {	
		printf("ERROR: could not open file 'nodes'\n");
		exit(1);
	}

	if (fscanf(file_nodes,"%d", &number_of_nodes) != 1) {
		printf("Could not read value\n");
		fclose(file_nodes);
		exit(1);
	}
	fclose(file_nodes);

	return number_of_nodes;
}

void initialize_matrix(int number_of_nodes) {

	int i, j;

	for (i = 0; i < number_of_nodes; i ++) {	
		number_of_outbound_links[i] = 0;
		for (j = 0; j < number_of_nodes; j ++) {
			A[i][j] = 0.0;
		}
	}
}

void populate_matrix(int number_of_nodes) {

	FILE *file_list;
	int i, j;
	char filename[128];

    sprintf(filename, "./dataset/%s/inv_adj_list", folder_name);

	file_list = fopen(filename,"r");

	if (file_list == NULL) {
		printf("ERROR: could not open file 'inv_adj_list'\n");
		exit(1);
	}

	for (i = 0; i < number_of_nodes; i ++) {

		if (fscanf(file_list,"%*d: %d",&j) != 1) {
			printf("Could not read value\n");
			fclose(file_list);
			exit(1);
		}
		
		while (j != -1) {

			A[i][j] = 1.0;
			number_of_outbound_links[j]++;

			if (fscanf(file_list,"%d",&j) != 1) {
				printf("Could not read value\n");
				fclose(file_list);
				exit(1);
			}
		
		}
	}

	fclose(file_list);
}

//Calculating the H = S + p*d matrix and setting up the dangling points.
void replace_dangling_nodes(int number_of_nodes) {

	int i, j, dangling_node;

	//Calculating the H = S + p*d matrix.
	for ( j = 0; j < number_of_nodes; j++) {

		dangling_node = 1;
		for (i = 0; i < number_of_nodes; i++) {

			if (A[i][j] > 0.0) {
				dangling_node = 0;
				break;
			}

		}

		if (dangling_node == 0) {
			continue;
		}

		for (i = 0; i < number_of_nodes; i++) {
			A[i][j] = 1.0;
		}

		number_of_outbound_links[j] = number_of_nodes;

	}
}

//Calculating the (I - a * H) matrix.
void calculate_A_matrix(int number_of_nodes, double alpha) {
	
	int i, j;

	for (j = 0; j < number_of_nodes; j++) {

		for (i = 0; i < number_of_nodes; i++) {

			A[i][j] = A[i][j] / number_of_outbound_links[j];
			A[i][j] = (-alpha)*A[i][j];

			if (i == j) {
				A[i][j] += 1.0;
			}
		}
	}

}

//Calculating the norm_1.
double calculate_error(int number_of_nodes) {

    int i;
    double error = 0.0;

    for (i = 0; i < number_of_nodes; i++) {

        error += fabs(page_rank[i] - page_rank_previous[i]);

    }

    return error;

}

//Actual PageRank calculation.
void gauss_seidel(int number_of_nodes, int max_number_of_iterations) {
    
	int i, j, k, l;
    double sum;

    for (k = 0; k < max_number_of_iterations; k++) {

        for (i = 0; i < number_of_nodes; i++) {

            sum = 0.0;
            for (j = 0; j < number_of_nodes; j++) {

                if (i != j) {

                    sum += A[i][j]*page_rank[j];
				
                }

            }

			//Storing the previous value of i-th page rank.
			page_rank_previous[i] = page_rank[i];

            page_rank[i] = (1/A[i][i]) * (b[i] - sum);

        }

		if (calculate_error(number_of_nodes) < desired_error) {
			return;
		}

    }

}

