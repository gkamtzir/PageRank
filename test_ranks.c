/**
* FILE: test_ranks.c
* THMMY, 7th semester, Parallel and Distributed Systems: 4th assignment
* This program tests if the result from the compressed algorithms are
* valid.
* Author:
*   Kamtziridis Georgios, 8542, gkamtzir@auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int get_total_number_of_nodes();

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

	int i, j, number_of_nodes, counter;
	double compressed, gauss;

	FILE *compressed_file, *gauss_file;

	compressed_file = fopen("./page_rank_compressed", "r");
	gauss_file = fopen("./page_rank_simple", "r");

	counter = 0;
	number_of_nodes = get_total_number_of_nodes();
	printf("%d \n", number_of_nodes);

	for (i = 0; i < number_of_nodes; i++) {
		
		if (fscanf(compressed_file,"%lf", &compressed) == 1 && fscanf(gauss_file, "%lf", &gauss)) {

			if (fabs(compressed - gauss) > 0.0000000001) {
				counter++;
				printf("Page %d (compressed): %f \n", i, compressed);
				printf("Page %d (simple): %f \n", i, gauss);
			}

		} else {

			printf("Could not read value\n");
			exit(1);

		}
		
	}

	printf("counter: %d \n", counter);


	fclose(compressed_file);
	fclose(gauss_file);

}

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
