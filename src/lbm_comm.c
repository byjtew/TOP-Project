/********************  HEADERS  *********************/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include "lbm_comm.h"

// If OpenMP is not available, define OPENMP_AVAILABLE
#ifndef _OPENMP
#define OPENMP_AVAILABLE 0
#pragma message "OpenMP is not available"
#else
#define OPENMP_AVAILABLE 1
#pragma message "OpenMP is available"
#endif

/*******************  FUNCTION  *********************/
int lbm_helper_pgcd(int a, int b) {
	int c;
	while (b != 0) {
		c = a % b;
		a = b;
		b = c;
	}
	return a;
}

void lbm_comm_timers_start(lbm_comm_t *mesh_comm, int id) {
	if (id >= NB_TIMERS) {
		fprintf(stderr, "timer_id too high.\n");
		abort();
	}
	mesh_comm->timers[id][mesh_comm->current_timer[id]] = MPI_Wtime();
}

void lbm_comm_timers_stop(lbm_comm_t *mesh_comm, int id) {
	if (id >= NB_TIMERS) {
		fprintf(stderr, "timer_id too high.\n");
		abort();
	}
	mesh_comm->timers[id][mesh_comm->current_timer[id]] =
			MPI_Wtime() - mesh_comm->timers[id][mesh_comm->current_timer[id]];
	mesh_comm->current_timer[id]++;
}

/*******************  FUNCTION  *********************/
/**
 * Affiche la configuation du lbm_comm pour un rank donné
 * @param mesh_comm Configuration à afficher
**/
void lbm_comm_print(lbm_comm_t *mesh_comm) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf(" RANK %d ( LEFT %d RIGHT %d TOP %d BOTTOM %d CORNER %d, %d, %d, %d ) ( POSITION %d %d ) (WH %d %d ) \n", rank,
	       mesh_comm->left_id,
	       mesh_comm->right_id,
	       mesh_comm->top_id,
	       mesh_comm->bottom_id,
	       mesh_comm->corner_id[0],
	       mesh_comm->corner_id[1],
	       mesh_comm->corner_id[2],
	       mesh_comm->corner_id[3],
	       mesh_comm->x,
	       mesh_comm->y,
	       mesh_comm->width,
	       mesh_comm->height);
}

/*******************  FUNCTION  *********************/
int helper_get_rank_id(int nb_x, int nb_y, int rank_x, int rank_y) {
	if (rank_x < 0 || rank_x >= nb_x)
		return -1;
	else if (rank_y < 0 || rank_y >= nb_y)
		return -1;
	else
		return (rank_x + rank_y * nb_x);
}

/*******************  FUNCTION  *********************/
/**
 * Initialise un lbm_comm :
 * - Voisins
 * - Taille du maillage local
 * - Position relative
 * @param mesh_comm MeshComm à initialiser
 * @param rank Rank demandant l'initalisation
 * @param comm_size Taille totale du communicateur
 * @param width largeur du maillage
 * @param height hauteur du maillage
**/
void lbm_comm_init(lbm_comm_t *mesh_comm, int rank, int comm_size, int width, int height) {
	//vars
	int nb_x;
	int nb_y;
	int rank_x;
	int rank_y;

	if (rank == 0) {
		fprintf(stderr, "lbm_comm_init(): Hauteur du maillage %d\n", height);
		fprintf(stderr, "lbm_comm_init(): Largeur du maillage %d\n", width);
	}

	//compute splitting
	nb_y = lbm_helper_pgcd(comm_size, width);
	nb_x = comm_size / nb_y;
	//nb_x = lbm_helper_pgcd(comm_size,width);
	//nb_y = comm_size / nb_x;

	//check
	assert(nb_x * nb_y == comm_size);
	if (height % nb_y != 0 && rank == 0)
		fatal("Can't get a 2D cut for current problem size and number of processes.");

	//calc current rank position (ID)
	rank_x = rank % nb_x;
	rank_y = rank / nb_x;

	//setup nb
	mesh_comm->nb_x = nb_x;
	mesh_comm->nb_y = nb_y;

	//setup size (+2 for ghost cells on border)
	mesh_comm->width = width / nb_x + 2;
	mesh_comm->height = height / nb_y + 2;

	//setup position
	mesh_comm->x = rank_x * width / nb_x;
	mesh_comm->y = rank_y * height / nb_y;

	// Compute neighbour nodes id
	mesh_comm->left_id = helper_get_rank_id(nb_x, nb_y, rank_x - 1, rank_y);
	mesh_comm->right_id = helper_get_rank_id(nb_x, nb_y, rank_x + 1, rank_y);
	mesh_comm->top_id = helper_get_rank_id(nb_x, nb_y, rank_x, rank_y - 1);
	mesh_comm->bottom_id = helper_get_rank_id(nb_x, nb_y, rank_x, rank_y + 1);
	mesh_comm->corner_id[CORNER_TOP_LEFT] = helper_get_rank_id(nb_x, nb_y, rank_x - 1, rank_y - 1);
	mesh_comm->corner_id[CORNER_TOP_RIGHT] = helper_get_rank_id(nb_x, nb_y, rank_x + 1, rank_y - 1);
	mesh_comm->corner_id[CORNER_BOTTOM_LEFT] = helper_get_rank_id(nb_x, nb_y, rank_x - 1, rank_y + 1);
	mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] = helper_get_rank_id(nb_x, nb_y, rank_x + 1, rank_y + 1);

	mesh_comm->max_requests = width * 4 + 2 * 4 + 4;
	mesh_comm->requests = calloc(mesh_comm->max_requests, sizeof(MPI_Request));
	mesh_comm->statuses = calloc(mesh_comm->max_requests, sizeof(MPI_Request));

	//if more than 1 on y, need transmission buffer
	if (nb_y > 1) {
		mesh_comm->buffer = malloc(sizeof(double) * DIRECTIONS * width / nb_x);
	} else {
		mesh_comm->buffer = NULL;
	}

	for (int i = 0; i < NB_TIMERS; i++) {
		mesh_comm->current_timer[i] = 0;
		mesh_comm->timers[i] = calloc(ITERATIONS, sizeof(double));
	}


	//region Graph buffers

	int row_length = mesh_comm->width - 2, col_length = mesh_comm->height - 2;
	int borders_count = DIRECTIONS * (2 * (row_length + col_length) + 4);
	// [top-left, top, top-right, right, bot-right, bot, bot-left, left]
	mesh_comm->send_borders = (double *) calloc(borders_count, sizeof(double));
	mesh_comm->recv_borders = (double *) calloc(borders_count, sizeof(double));

	int nb_neighs = 8;
	int sources[8];
	sources[0] = mesh_comm->corner_id[CORNER_TOP_LEFT] >= 0 ? mesh_comm->corner_id[CORNER_TOP_LEFT] : rank;
	sources[1] = mesh_comm->top_id >= 0 ? mesh_comm->top_id : rank;
	sources[2] = mesh_comm->corner_id[CORNER_TOP_RIGHT] >= 0 ? mesh_comm->corner_id[CORNER_TOP_RIGHT] : rank;
	sources[3] = mesh_comm->right_id >= 0 ? mesh_comm->right_id : rank;
	sources[4] = mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] >= 0 ? mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] : rank;
	sources[5] = mesh_comm->bottom_id >= 0 ? mesh_comm->bottom_id : rank;
	sources[6] = mesh_comm->corner_id[CORNER_BOTTOM_LEFT] >= 0 ? mesh_comm->corner_id[CORNER_BOTTOM_LEFT] : rank;
	sources[7] = mesh_comm->left_id >= 0 ? mesh_comm->left_id : rank;


	mesh_comm->nb_per_neigh[0] = DIRECTIONS * (mesh_comm->corner_id[CORNER_TOP_LEFT] >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[1] = DIRECTIONS * row_length * (mesh_comm->top_id >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[2] = DIRECTIONS * (mesh_comm->corner_id[CORNER_TOP_RIGHT] >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[3] = DIRECTIONS * col_length * (mesh_comm->right_id >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[4] = DIRECTIONS * (mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[5] = DIRECTIONS * row_length * (mesh_comm->bottom_id >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[6] = DIRECTIONS * (mesh_comm->corner_id[CORNER_BOTTOM_LEFT] >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[7] = DIRECTIONS * col_length * (mesh_comm->left_id >= 0 ? 1 : 0);

	mesh_comm->displ_per_neigh[0] = 0;
	for (int i = 1; i < 8; i++)
		mesh_comm->displ_per_neigh[i] = mesh_comm->displ_per_neigh[i - 1] + mesh_comm->nb_per_neigh[i - 1];

	// endregion


	assert(nb_neighs >= 1);

	int weights[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0; i < 8; i++)
		weights[i] = mesh_comm->nb_per_neigh[i] > 0;

	MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, nb_neighs, sources, weights,
	                               nb_neighs, sources, weights,
	                               MPI_INFO_NULL, 1, &mesh_comm->comm_graph);
	if (rank == 0) fprintf(stderr, "Graph created.\n");
	/*
	usleep(rank*10);
	int source_neighs[8], dest_neighs[8];
	int out_nb_neighs = 0;
	MPI_Dist_graph_neighbors_count(mesh_comm->comm_graph, &nb_neighs, &out_nb_neighs, MPI_UNWEIGHTED);
	printf("[P%d]: Entry %d ; Out %d\n", rank, nb_neighs, out_nb_neighs);
	MPI_Dist_graph_neighbors(mesh_comm->comm_graph, 8, source_neighs, MPI_UNWEIGHTED, 8, dest_neighs, MPI_UNWEIGHTED);
	printf("[P%d]: SRC Graph_Neighbours{%d} : ", rank, nb_neighs);
	for (int i = 0; i < out_nb_neighs; i++)
		if (source_neighs[i] != rank)printf("%.02d ", source_neighs[i]);
		else printf("-- ");
	printf("\n");
	printf("[P%d]: DST Graph_Neighbours{%.02d}: ", rank, nb_neighs);
	for (int i = 0; i < nb_neighs; i++)
		if (dest_neighs[i] != rank)printf("%.02d ", dest_neighs[i]);
		else printf("-- ");
	printf("\n");


	MPI_Barrier(MPI_COMM_WORLD);
	 */
	//if debug print comm
#ifndef NDEBUG
	lbm_comm_print(mesh_comm);
#endif
}


/*******************  FUNCTION  *********************/
/**
 * Libere un lbm_comm
 * @param mesh_comm MeshComm à liberer
**/
void lbm_comm_release(lbm_comm_t *mesh_comm) {
	mesh_comm->x = 0;
	mesh_comm->y = 0;
	mesh_comm->width = 0;
	mesh_comm->height = 0;
	mesh_comm->right_id = -1;
	mesh_comm->left_id = -1;
	if (mesh_comm->buffer != NULL)
		free(mesh_comm->buffer);
	mesh_comm->buffer = NULL;
	free(mesh_comm->requests);
	free(mesh_comm->statuses);

	// Graph section
	free(mesh_comm->send_borders);
	free(mesh_comm->recv_borders);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == RANK_MASTER) {
		printf("Output timers.\n");
		// Open a file name timers.txt and write the timers
		char filename[64];
		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		for (int o = 0; o < NB_USED_TIMER; o++) {
			sprintf(filename, "timers-%d-t%d-n%d-i%d-o%d.txt", o, OPENMP_AVAILABLE, world_size, ITERATIONS,
			        WRITE_STEP_INTERVAL);
			FILE *fp = fopen(filename, "w");
			double mean = 0;
			for (int t = 0; t < mesh_comm->current_timer[o]; t++) mean += mesh_comm->timers[o][t];
			mean /= mesh_comm->current_timer[o];
			fprintf(fp, "0 %lf\n", mean);
			for (int t = 0; t < mesh_comm->current_timer[o]; t++)
				fprintf(fp, "%d %lf\n", t + 1, mesh_comm->timers[o][t]);
			fclose(fp);
		}
	}

	for (int i = 0; i < NB_TIMERS; i++)
		free(mesh_comm->timers[i]);
}

/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh Mesh a utiliser lors de l'échange des mailles fantomes
**/
void
lbm_comm_sync_ghosts_horizontal(lbm_comm_t *mesh_comm, Mesh *mesh, lbm_comm_type_t comm_type, int target_rank,
                                int x) {
	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type) {
		case COMM_SEND:
			MPI_Isend(Mesh_get_col(mesh, x), DIRECTIONS * (mesh_comm->height - 2), MPI_DOUBLE, target_rank, 0,
			          MPI_COMM_WORLD, &mesh_comm->requests[mesh_comm->current_request++]);
			break;
		case COMM_RECV:
			MPI_Irecv(Mesh_get_col(mesh, x), DIRECTIONS * (mesh_comm->height - 2), MPI_DOUBLE, target_rank, 0,
			          MPI_COMM_WORLD, &mesh_comm->requests[mesh_comm->current_request++]);
			break;
		default:
			fatal("Unknown type of communication.");
	}
}

/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh Mesh a utiliser lors de l'échange des mailles fantomes
**/
void lbm_comm_sync_ghosts_diagonal(lbm_comm_t *mesh_comm, Mesh *mesh, lbm_comm_type_t comm_type, int target_rank,
                                   int x, int y) {
	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type) {
		case COMM_SEND:
			MPI_Isend(Mesh_get_cell(mesh, x, y), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,
			          &mesh_comm->requests[mesh_comm->current_request++]);
			break;
		case COMM_RECV:
			MPI_Irecv(Mesh_get_cell(mesh, x, y), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,
			          &mesh_comm->requests[mesh_comm->current_request++]);
			break;
		default:
			fatal("Unknown type of communication.");
	}
}


/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh Mesh a utiliser lors de l'échange des mailles fantomes
**/
int lbm_comm_sync_ghosts_vertical(lbm_comm_t *mesh_comm, Mesh *mesh, lbm_comm_type_t comm_type, int target_rank,
                                  int y, double *buffer) {
	if (target_rank == -1) return 0;

	if (comm_type == COMM_SEND) {
		for (int x = 1; x < mesh->width - 2; x++)
			memcpy(buffer + (x - 1) * DIRECTIONS, Mesh_get_cell(mesh, x, y), sizeof(double) * DIRECTIONS);
		MPI_Isend(buffer, DIRECTIONS * (mesh_comm->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,
		          &mesh_comm->requests[mesh_comm->current_request++]);
	} else {
		MPI_Irecv(buffer, DIRECTIONS * (mesh_comm->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,
		          &mesh_comm->requests[mesh_comm->current_request++]);
		return 1;
	}
	return 0;
}

void checkMeshValid(Mesh *mesh) {
	for (int i = 0; i < mesh->width; i++) {
		for (int j = 0; j < mesh->height; j++) {
			lbm_mesh_cell_t cell = Mesh_get_cell(mesh, i, j);
			for (int k = 0; k < DIRECTIONS; k++) {
				double cell_k = cell[k];
				if (isnan(cell_k))
					fatal("NaN detected in mesh");
				else if (isinf(cell_k))
					fatal("Inf detected in mesh");
				else if (cell_k < 0 || cell_k > 8)
					fatal("Invalid value in mesh");
			}
		}
	}
}

// Todo: switch to mesh_comm->displs constants
// Todo: Thread in tasks and loops
void build_send_borders_from_mesh(Mesh *mesh, lbm_comm_t *mesh_comm) {
	int cur_idx = 0;
	int col_length = mesh->height - 2;

	memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_cell(mesh, 1, 1), sizeof(double) * DIRECTIONS);
	cur_idx += DIRECTIONS;

	for (int i = 1; i < mesh->width - 1; i++) {
		memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_cell(mesh, i, 1), sizeof(double) * DIRECTIONS);
		cur_idx += DIRECTIONS;
	}

	memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_cell(mesh, mesh->width - 2, 1), sizeof(double) * DIRECTIONS);
	cur_idx += DIRECTIONS;

	memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_col(mesh, mesh->width - 2),
	       sizeof(double) * DIRECTIONS * col_length);
	cur_idx += DIRECTIONS * col_length;

	memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_cell(mesh, mesh->width - 2, mesh->height - 2),
	       sizeof(double) * DIRECTIONS);
	cur_idx += DIRECTIONS;

	for (int i = 1; i < mesh->width - 1; i++) {
		memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_cell(mesh, i, mesh->height - 2), sizeof(double) * DIRECTIONS);
		cur_idx += DIRECTIONS;
	}

	memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_cell(mesh, 1, mesh->height - 2), sizeof(double) * DIRECTIONS);
	cur_idx += DIRECTIONS;

	memcpy(mesh_comm->send_borders + cur_idx, Mesh_get_col(mesh, 1), sizeof(double) * DIRECTIONS * col_length);
	cur_idx += DIRECTIONS * col_length;
}

// Todo: Thread the loop
void emplace_row_buffer(Mesh *mesh, double *buffer, int y) {
#pragma omp parallel for default(none) shared(mesh, buffer, y)
	for (int x = 1; x < mesh->width - 1; x++)
		memcpy(Mesh_get_cell(mesh, x, y), buffer + (x - 1) * DIRECTIONS, sizeof(double) * DIRECTIONS);
}

// Todo: unroll and thread in tasks
void emplace_recv_buffer(Mesh *mesh, lbm_comm_t *mesh_comm) {
	int *displs = mesh_comm->displ_per_neigh;
	double *buffer = mesh_comm->recv_borders;
#pragma omp parallel default(none) shared(mesh, buffer, mesh_comm, displs)
	{
#pragma omp single
		{
#pragma omp task default(none) shared(mesh_comm, mesh, buffer, displs) // 1 rown
			{
				if (mesh_comm->nb_per_neigh[1] > 0) // Top
					emplace_row_buffer(mesh, buffer + displs[1], 0);
			}
#pragma omp task default(none) shared(mesh_comm, mesh, buffer, displs) // 2 corners + 1 column
			{
				if (mesh_comm->nb_per_neigh[2] > 0) // Top-right
					memcpy(Mesh_get_cell(mesh, mesh->width - 1, 0), buffer + displs[2], DIRECTIONS);
				if (mesh_comm->nb_per_neigh[3] > 0) // Right
					memcpy(Mesh_get_col(mesh, mesh->width - 1), buffer + displs[3], DIRECTIONS * mesh->height - 2);
				if (mesh_comm->nb_per_neigh[4] > 0) // Bottom-right
					memcpy(Mesh_get_cell(mesh, mesh->width - 1, mesh->height - 1), buffer + displs[4], DIRECTIONS);
			}
#pragma omp task default(none) shared(mesh_comm, mesh, buffer, displs) // 1 row
			{
				if (mesh_comm->nb_per_neigh[5] > 0) // Bottom
					emplace_row_buffer(mesh, buffer + displs[5], mesh->height - 1);
			}
#pragma omp task default(none) shared(mesh_comm, mesh, buffer, displs) // 2 corners + 1 column
			{
				if (mesh_comm->nb_per_neigh[6] > 0) // Bottom-left
					memcpy(Mesh_get_cell(mesh, 0, mesh->height - 1), buffer + displs[6], DIRECTIONS);
				if (mesh_comm->nb_per_neigh[7] > 0) // Left
					memcpy(Mesh_get_col(mesh, 0), buffer + displs[7], DIRECTIONS * mesh->height - 2);
				if (mesh_comm->nb_per_neigh[0] > 0) // Top-left
					memcpy(Mesh_get_cell(mesh, 0, 0), buffer + displs[0], DIRECTIONS);
			}
#pragma omp taskwait
		}
	}

}

/*******************  FUNCTION  *********************/
void lbm_comm_ghost_exchange(lbm_comm_t *mesh_comm, Mesh *mesh, int rank) {
	/*
	// Build borders buffers (send & recv)
	fprintf(stderr, "Mesh size: %d x %d\n", mesh->width, mesh->height);
	fprintf(stderr, "Mesh_comm size: %d x%d\n", mesh_comm->width, mesh_comm->height);
	fprintf(stderr, "row_length: %d\n", row_length);
	fprintf(stderr, "col_length: %d\n", col_length);
	*/

	lbm_comm_timers_start(mesh_comm, TIMER_SEND_BUFFER_CREATE);
	build_send_borders_from_mesh(mesh, mesh_comm);
	lbm_comm_timers_stop(mesh_comm, TIMER_SEND_BUFFER_CREATE);


	/*fprintf(stderr, "send_counts: ");
	for(int i=0; i<8; i++) fprintf(stderr, "%d ", send_counts[i]);
	fprintf(stderr, "\n");
*/

	lbm_comm_timers_start(mesh_comm, TIMER_MESH_SYNC);
	MPI_Neighbor_alltoallv(mesh_comm->send_borders, mesh_comm->nb_per_neigh, mesh_comm->displ_per_neigh, MPI_DOUBLE,
	                       mesh_comm->recv_borders, mesh_comm->nb_per_neigh, mesh_comm->displ_per_neigh, MPI_DOUBLE,
	                       mesh_comm->comm_graph);
	lbm_comm_timers_stop(mesh_comm, TIMER_MESH_SYNC);
/*
	sleep(3*rank);
	printf("Process %d, recv:\n", rank);
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < send_counts[i]; j++)
			printf("%f ", recv_borders_buffer[sdispls[i] + j]);
		printf("\n");
		fflush(stdout);
	}
	printf("\n==========================================================\n\n");
*/

	lbm_comm_timers_start(mesh_comm, TIMER_RECV_BUFFER_EMPLACE);
	emplace_recv_buffer(mesh, mesh_comm);
	lbm_comm_timers_stop(mesh_comm, TIMER_RECV_BUFFER_EMPLACE);


#ifndef RELEASE_MODE
	checkMeshValid(mesh);
#endif
}

/*******************  FUNCTION  *********************/
/**
 * Rendu du mesh en effectuant une réduction a 0
 * @param mesh_comm MeshComm à utiliser
 * @param temp Mesh a utiliser pour stocker les segments
**/
void save_frame_all_domain(FILE *fp, Mesh *source_mesh, Mesh *temp, lbm_comm_t *mesh_comm) {
	// Todo: Switch to a Gather
	//vars
	int comm_size, rank;

	//get rank and comm size
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	lbm_comm_timers_start(mesh_comm, TIMER_OUTPUT_GATHER);
	/* If whe have more than one process */
	if (1 < comm_size) {
		MPI_Gather(source_mesh->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE,
		           temp->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE, RANK_MASTER,
		           MPI_COMM_WORLD);
		if (rank == RANK_MASTER) {
			/* Rank 0 receives & render other processes meshes */
			for (int i = 0; i < comm_size; i++) {
				Mesh tmp_mesh;
				tmp_mesh.width = source_mesh->width;
				tmp_mesh.height = source_mesh->height;
				tmp_mesh.cells = temp->cells + i * source_mesh->width * source_mesh->height * DIRECTIONS;
				save_frame(fp, &tmp_mesh);
			}
		}
	} else {
		/* Only 0 renders its local mesh */
		save_frame(fp, source_mesh);
	}
	lbm_comm_timers_stop(mesh_comm, TIMER_OUTPUT_GATHER);

}

