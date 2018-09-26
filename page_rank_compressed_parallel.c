/**
* FILE: page_rank_compressed_parallel.c
* THMMY, 7th semester, Parallel and Distributed Systems: 4th assignment
* Parallel Page Rank algorithm using Gauss-Seidel approach.
* Author:
*   Kamtziridis Georgios, 8542, gkamtzir@auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void write_page_rank_to_file(double *page_rank, int number_of_nodes);
int get_total_number_of_nodes();
int get_number_of_links(int number_of_nodes, double *value_of_page, int *dangling_nodes);
void compressed_sparse_row(int number_of_nodes);
void prepare_values(double *values, int i, int number_of_nodes, int number_of_links, int parallel);
double get_element(int i, int j, int number_of_nodes, int number_of_links);
double calculate_error(int number_of_nodes);
void gauss_seidel(int number_of_nodes, int max_number_of_iterations, int number_of_links);

void graph_coloring(int number_of_nodes, int number_of_links);

typedef struct {
    int start;
    int end;
    int until;
    int color;
} Group;

//Global variables.
int *IA, *JA, *dangling_nodes;
double *value_of_page, *page_rank, *page_rank_previous, *b;
int number_of_threads;
double desired_error;
char folder_name[128];
Group *groups;
int number_of_groups;

int main(int argc, char **argv) {

    if (argc != 3) {
        printf("You have to provide a number that indicates which dataset to choose\n"
            "1: Gun Control (2955 nodes)\n2: Search Engines (11659 nodes)\n"
            "You, also, have to provide the number of threads as the second argument\n"
            "Example: comp_par 2 8\n");
        exit(1);
    }

    if (atoi(argv[1]) == 1) {

        sprintf(folder_name, "gun_control");

    } else {

        sprintf(folder_name, "search_engines");

    }

    number_of_threads = atoi(argv[2]);

    int i, j, number_of_nodes, number_of_links, max_number_of_iterations;
    double alpha, p;

    //Initializing.
    number_of_nodes = get_total_number_of_nodes();
    max_number_of_iterations = 50;
    desired_error = 0.000001;
    alpha = 0.85;
    p = 1.0 / number_of_nodes;

    //Allocating memory
    JA = (int *)malloc(sizeof(int) * number_of_nodes);
    value_of_page = (double *)malloc(sizeof(double) * number_of_nodes);
    dangling_nodes = (int *)malloc(sizeof(int) * number_of_nodes);
    page_rank = (double *)malloc(sizeof(double) * number_of_nodes);
    page_rank_previous = (double *)malloc(sizeof(double) * number_of_nodes);
    b = (double *)malloc(sizeof(double) * number_of_nodes);
    groups = (Group *)malloc(sizeof(Group) * number_of_nodes);

    if (JA == NULL || value_of_page == NULL || dangling_nodes == NULL
        || page_rank == NULL || page_rank_previous == NULL || b == NULL) {

        printf("Could not allocate memory. \n");
        exit(1);

    }

    /*
        Getting the total number of outbound (or inbound) links
        and populating the value_of_page and dangling_node vectors
        accordingly.
    */
    number_of_links = get_number_of_links(number_of_nodes, value_of_page, dangling_nodes);

    IA = (int *)malloc(sizeof(int) * number_of_links);

    if (IA == NULL) {

        printf("Could not allocate memory. \n");
        exit(1);

    }

    compressed_sparse_row(number_of_nodes);

    for (i = 0; i < number_of_nodes; i++) {

        //Calculating the (-a) * H vector.
        value_of_page[i] = value_of_page[i] *(-alpha);

        //Initializing the page_rank vector.
        page_rank[i] = 1.0 / number_of_nodes;
        page_rank[i] = 1.0 / number_of_nodes;

        //Initializing the 'b' vector.
    	b[i] = (1-alpha) * p;
    
    }

    graph_coloring(number_of_nodes, number_of_links);

    //Time the Page Rank algorithm.
    double start = omp_get_wtime();
    gauss_seidel(number_of_nodes, max_number_of_iterations, number_of_links);
    double end = omp_get_wtime();
    printf("%f \n", end - start);

    write_page_rank_to_file(page_rank, number_of_nodes);

    //Freeing memory.
    free(JA);
    free(value_of_page);
    free(dangling_nodes);
    free(page_rank);
    free(page_rank_previous);
    free(b);
    free(groups);

}

/*
    Writes the final Page Rank vector to a file.
*/
void write_page_rank_to_file(double *page_rank, int number_of_nodes) {

    FILE *file;
    int i;

    file= fopen("./page_rank_compressed", "w");

    for (i = 0; i < number_of_nodes; i++) {

        fprintf(file, "%f\n", page_rank[i]);
    
    }

    fclose(file);

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

/*
    Counts the total number of outbound (or inbound) links
    and populates the 'value_of_page' and 'dangling_nodes'
    vectors.
*/
int get_number_of_links(int number_of_nodes, double *value_of_page, int *dangling_nodes) {

    FILE *file;
	int i, j, number_of_non_zero_elements, counter;
    char filename[128];

    sprintf(filename, "./dataset/%s/adj_list", folder_name);

    number_of_non_zero_elements = 0;

	file = fopen(filename, "r");

	if (file == NULL) {

		printf("ERROR: could not open file 'adj_list'\n");
		exit(1);

	}

	for (i = 0; i < number_of_nodes; i++) {

        if (fscanf(file,"%*d: %d",&j) != 1) {
            printf("Could not read value\n");
            fclose(file);
            exit(1);
        }

        if (j == -1) {

            value_of_page[i] = 1.0 / number_of_nodes;
            dangling_nodes[i] = 1;
            continue;

        }

        counter = 0;

		while (j != -1) {

			counter++;
            number_of_non_zero_elements++;

            if (fscanf(file, "%d", &j) != 1) {
                printf("Could not read value\n");
                fclose(file);
                exit(1);
            }
		
        }

        dangling_nodes[i] = 0;
        value_of_page[i] = 1.0 / counter;

    }

    fclose(file);

    return number_of_non_zero_elements;

}

/*
    Compresses the given data using CSR algorithm.
    The exact algorithm is explained in the report
    and in reference [5].
*/
void compressed_sparse_row(int number_of_nodes) {

	FILE *file;
	int i, j, index, found;
    char filename[128];

    sprintf(filename, "./dataset/%s/inv_adj_list", folder_name);

    index = 0;

	file = fopen(filename, "r");

	if (file == NULL) {

		printf("ERROR: could not open file 'inv_adj_list'\n");
		exit(1);

	}

	for (i = 0; i < number_of_nodes; i++) {

        if (fscanf(file,"%*d: %d",&j) != 1) {
            printf("Could not read value\n");
            fclose(file);
            exit(1);
        }

        found = 0;
		
		while (j != -1) {

            IA[index] = j;
            if (!found) {

                JA[i] = index;
                found = 1;
            }

            index++;

            if (fscanf(file,"%d",&j) != 1) {
                printf("Could not read value\n");
                fclose(file);
                exit(1);
            }
		
        }

        if (!found) {

            JA[i] = -1;

        }

	}

	fclose(file);

}

/*
    Fetches the (i, j) element from the compressed matrix.
    The i and j correspond to the row and column of the element
    in the original matrix.
*/
double get_element(int i, int j, int number_of_nodes, int number_of_links) {

    int index, index_next_row;

    index = JA[i];

    if (index == -1) {

        return i == j ? 1.0 : 0.0;
    
    }

    if (i == number_of_nodes - 1) {

        index_next_row = number_of_links;

    } else {

        /*
            Searching the index value of IA's vector that indicates
            where the next row begins. If that index is -1 then
            we know that the next row is full of zeros. So, we keep 
            searching for an actual index. If we reach the end then 
            we set 'end' to equal to the overall length of the vector.
        */
        int end = 1;

        while (JA[i + end] == -1) {
            
            end++;
            if (i + end == number_of_nodes) {

                break;

            }

        }

        if (i + end == number_of_nodes) {

            index_next_row = number_of_links;
        
        } else {

            index_next_row = JA[i + end];
        
        }

    }

    while(1) {

        if ( IA[index] > j || index == index_next_row) {

            return i == j ? 1.0 : 0.0;

        } else if (IA[index] == j) {

            return i == j ? 1.0 + value_of_page[j] : value_of_page[j];

        }

        index++;
    
    }

}

/*
    Separates the nodes into groups of colors so they can
    be executed in parallel later. The exact algorithm is
    explained in the report.
*/
void graph_coloring(int number_of_nodes, int number_of_links) {

    int i, j, min, max;
    double value;

    int group = -1;

    for (j = 0; j < number_of_nodes; j++) {

        max = -1;
        min = number_of_nodes;

        /*
            If we have a dangling node then we consider it 
            a group on its own.
        */
        if (dangling_nodes[j] == 1) {
            
            //New color.
            group++;
            groups[group].color = group;
            groups[group].start = j;
            groups[group].end = j;
            groups[group].until = j;
            continue;
        }

        for (i = 0; i < number_of_nodes; i++) {

            value = get_element(i, j, number_of_nodes, number_of_links);

            if (fabs(value) > 0.00 && i < j) {

                if (i == j - 1) {

                    max = i;
                    break;

                }

                if (i > max) {

                    max = i;

                }

                if (i < min && i > j) {

                    min = i;

                }

            }
        
        }

        if (max == j - 1) {

            //New color.
            group++;
            groups[group].color = group;
            groups[group].start = j;
            groups[group].end = j;
            groups[group].until = min;

        } else if (max < j - 1) {

            //Merge with previous color if (j-1) color has not the same color with (max).
            if (max >= groups[group].start && max <= groups[group].end) {
                
                //New color.
                group++;
                groups[group].color = group;
                groups[group].start = j;
                groups[group].end = j;
                groups[group].until = min;

            } else {

                if (groups[group].until < j) {

                    //New color.
                    group++;
                    groups[group].color = group;
                    groups[group].start = j;
                    groups[group].end = j;
                    groups[group].until = min;

                } else {

                    groups[group].end = j;
                    if (groups[group].until > min) {
                        groups[group].until = min;
                    }

                }

            }

        } else {

            if (group > -1) {

                if (groups[group].until < j) {

                    //New color.
                    group++;
                    groups[group].color = group;
                    groups[group].start = j;
                    groups[group].end = j;
                    groups[group].until = min;

                } else {

                    groups[group].end = j;
                    if (groups[group].until > min) {
                        groups[group].until = min;
                    }

                }

            } else {

                //New color.
                group++;
                groups[group].color = group;
                groups[group].start = j;
                groups[group].end = j;
                groups[group].until = min;

            }

        }

    }

    number_of_groups = group + 1;

}

void prepare_values(double *values, int i, int number_of_nodes, int number_of_links, int parallel) {

    int j;

    #pragma omp parallel for private(j) num_threads(number_of_threads) if (parallel)
    for (j = 0; j < number_of_nodes; j++) {

        if (i != j) {

            //If the current node is dangling then we don't need to call 'get_element'.
            values[j] = (dangling_nodes[j] == 1) ? value_of_page[j] : get_element(i, j, number_of_nodes, number_of_links);

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

/*
    Calculates the Page Rank vector.
*/
void gauss_seidel(int number_of_nodes, int max_number_of_iterations, int number_of_links) {

    int i, j, k, l, group_size;
    double sum;

    for (k = 0; k < max_number_of_iterations; k++) {

        for (l = 0; l < number_of_groups; l++) {

            group_size = groups[l].end - groups[l].start + 1;

            /*
                If group_size is greater than 2 then we can use
                a parallel for. Otherwise, it's not efficient.
            */
            if (group_size > 2) {

                #pragma omp parallel for private(i, j, sum) num_threads(number_of_threads)
                for (i = groups[l].start; i <= groups[l].end; i++) {

                    sum = 0.0;
                    double *values = (double *)malloc(sizeof(double) * number_of_nodes);
                    prepare_values(values, i, number_of_nodes, max_number_of_iterations, 0);

                    for (j = 0; j < number_of_nodes; j++) {

                        if (i != j) {

                            sum += values[j] * page_rank[j];

                        }
                    
                    }

                    //Storing the previous value of i-th page rank.
                    page_rank_previous[i] = page_rank[i];

                    page_rank[i] = (dangling_nodes[i] == 1) ? (1/(value_of_page[i] + 1.0)) * (b[i] - sum) : ( 1 /get_element(i, i, number_of_nodes, number_of_links)) * (b[i] - sum);
                
                    free(values);
                
                }

            } else {
                
                for (i = groups[l].start; i <= groups[l].end; i++) {

                    sum = 0.0;
                    double *values = (double *)malloc(sizeof(double) * number_of_nodes);
                    prepare_values(values, i, number_of_nodes, max_number_of_iterations, 1);

                    for (j = 0; j < number_of_nodes; j++) {

                        if (i != j) {

                            sum += values[j] * page_rank[j];

                        }
                    
                    }

                    //Storing the previous value of i-th page rank.
                    page_rank_previous[i] = page_rank[i];

                    page_rank[i] = (dangling_nodes[i] == 1) ? (1/(value_of_page[i] + 1.0)) * (b[i] - sum) : ( 1 /get_element(i, i, number_of_nodes, number_of_links)) * (b[i] - sum);
                
                    free(values);
                
                }

            }

            if (calculate_error(number_of_nodes) < desired_error) {
                printf("Iterations: %d \n", k + 1);
                return;
            }

        }
    
    }

}