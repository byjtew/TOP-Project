/********************  HEADERS  *********************/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <pthread.h>
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
	printf(" RANK %d ( LEFT %d RIGHT %d) ( POSITION %d %d ) (WH %d %d ) \n", rank,
	       mesh_comm->left_id,
	       mesh_comm->right_id,
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
	nb_y = 1;
	nb_x = comm_size;

	//check
	assert(nb_x * nb_y == comm_size);
	assert(width % nb_x == 0); // If row-major
	if (height % nb_y != 0 && rank == 0)
		fatal("Can't get a 2D cut for current problem size and number of processes.");

	//calc current rank position (ID)
	rank_x = rank % nb_x;
	rank_y = rank / nb_x;
	printf("Rank %d is at %d %d\n", rank, rank_x, rank_y);

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

	mesh_comm->max_requests = width * 4 + 2 * 4 + 4;
	mesh_comm->requests = calloc(mesh_comm->max_requests, sizeof(MPI_Request));

	for (int i = 0; i < NB_TIMERS; i++) {
		mesh_comm->current_timer[i] = 0;
		mesh_comm->timers[i] = calloc(ITERATIONS, sizeof(double));
	}

	// Send and recv configuration (depending ROW/COL-MAJOR)
	int col_size = (mesh_comm->height - 2);
	mesh_comm->nb_per_neigh[0] = DIRECTIONS * col_size * (mesh_comm->left_id >= 0 ? 1 : 0);
	mesh_comm->nb_per_neigh[1] = DIRECTIONS * col_size * (mesh_comm->right_id >= 0 ? 1 : 0);
	// Top-left corner offset for the first row/col
	mesh_comm->recv_displ[0] = DIRECTIONS;
	// Bottom-left corner offset for the last row/col
	mesh_comm->recv_displ[1] = DIRECTIONS + mesh_comm->height * DIRECTIONS * (mesh_comm->width - 1);
	// Second row/col
	mesh_comm->send_displ[0] = mesh_comm->height * DIRECTIONS + DIRECTIONS;
	// Pre-last row/col
	mesh_comm->send_displ[1] = DIRECTIONS + mesh_comm->height * DIRECTIONS * (mesh_comm->width - 2);

#if MESH_SYNC_MODE == MESH_SYNC_GRAPH
#pragma message "Using graph sync"
	int nb_neighs = 2;
	int sources[2];
	sources[0] = mesh_comm->left_id >= 0 ? mesh_comm->left_id : rank;
	sources[1] = mesh_comm->right_id >= 0 ? mesh_comm->right_id : rank;

	int weights[2] = {1, 1};
	weights[0] = mesh_comm->nb_per_neigh[0] > 0;
	weights[1] = mesh_comm->nb_per_neigh[1] > 0;

	MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, nb_neighs, sources, weights,
	                               nb_neighs, sources, weights,
	                               MPI_INFO_NULL, 1, &mesh_comm->comm_graph);
	if (rank == 0) printf("Graph created.\n");
#else
#if MESH_SYNC_MODE == MESH_SYNC_CART
#pragma message "Using cartesian sync"
	// Cartesian buffers
	int dims[1];
	dims[0] = comm_size;
	int periods[1] = {0};
	int reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &mesh_comm->comm_graph);
	if (rank == 0) printf("Cartesian graph created.\n");
#endif
#endif

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
	if (mesh_comm->requests != NULL)
		free(mesh_comm->requests);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == RANK_MASTER) {
		char timers_names[NB_TIMERS][32] = {"sync_comm", "pre_sync", "post_sync", "io_gather_comm", "io_write",
		                                    "spec_cells", "collision", "propagation", "ghost_exchange_total"};
		printf("Output timers.\n");
		// Open a file name timers.txt and write the timers
		char filename[64];
		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		for (int o = 0; o < NB_USED_TIMER; o++) {
			printf("============================\nTIMER: %s\n", timers_names[o]);
			sprintf(filename, "timers-%.32s-n%d.txt", timers_names[o], world_size);
			FILE *fp = fopen(filename, "w");
			double mean = 0.0;
			for (int t = 0; t < mesh_comm->current_timer[o]; t++) mean += mesh_comm->timers[o][t];
			mean /= mesh_comm->current_timer[o];
			printf("Mean of %d: %lf\n", o, mean);

			double cumul = 0.0;
			for (int t = 0; t < mesh_comm->current_timer[o]; t++) cumul += mesh_comm->timers[o][t];
			printf("Cumul of %d: %lf\n", o, cumul);

			for (int t = 0; t < mesh_comm->current_timer[o]; t++)
				fprintf(fp, "%d %lf\n", t, mesh_comm->timers[o][t]);
			fclose(fp);
			printf("============================\n");
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


/*******************  FUNCTION  *********************/
void lbm_comm_ghost_exchange(lbm_comm_t *mesh_comm, Mesh *mesh, int rank) {
#if MESH_SYNC_MODE == MESH_SYNC_UNIT
	MPI_Request requests[4];
	int cur_request = 0;
#endif


#if MESH_SYNC_MODE == MESH_SYNC_UNIT
#pragma message "Using async unitary send & recv mesh synchronization"
	lbm_comm_timers_start(mesh_comm, TIMER_MESH_SYNC_COMM);
	if (mesh_comm->left_id >= 0) {
		MPI_Isend(Mesh_get_col(mesh, 1), mesh->height - 2, MPI_DOUBLE, mesh_comm->left_id, 0, MPI_COMM_WORLD,
							&requests[cur_request++]);
		MPI_Irecv(Mesh_get_col(mesh, 0), mesh->height - 2, MPI_DOUBLE, mesh_comm->left_id, 0, MPI_COMM_WORLD,
							&requests[cur_request++]);
	}
	if (mesh_comm->right_id >= 0) {
		MPI_Isend(Mesh_get_col(mesh, mesh->width - 2), mesh->height - 2, MPI_DOUBLE, mesh_comm->right_id, 0,
							MPI_COMM_WORLD,
							&requests[cur_request++]);
		MPI_Irecv(Mesh_get_col(mesh, mesh->width - 1), mesh->height - 2, MPI_DOUBLE, mesh_comm->right_id, 0,
							MPI_COMM_WORLD,
							&requests[cur_request++]);
	}
	MPI_Waitall(cur_request, requests, MPI_STATUSES_IGNORE);
	lbm_comm_timers_stop(mesh_comm, TIMER_MESH_SYNC_COMM);
#else
#if MESH_SYNC_MODE == MESH_SYNC_GRAPH || MESH_SYNC_MODE == MESH_SYNC_CART
	lbm_comm_timers_start(mesh_comm, TIMER_MESH_SYNC_COMM);
	MPI_Neighbor_alltoallv(Mesh_get_cell(mesh, 0, 0), mesh_comm->nb_per_neigh, mesh_comm->send_displ, MPI_DOUBLE,
	                       Mesh_get_cell(mesh, 0, 0), mesh_comm->nb_per_neigh, mesh_comm->recv_displ, MPI_DOUBLE,
	                       mesh_comm->comm_graph);
	lbm_comm_timers_stop(mesh_comm, TIMER_MESH_SYNC_COMM);
#endif
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

	/* If whe have more than one process */
	if (1 < comm_size) {
		lbm_comm_timers_start(mesh_comm, TIMER_IO_GATHER_COMM);
		MPI_Gather(Mesh_get_cell(source_mesh, 0, 0), source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE,
		           Mesh_get_cell(temp, 0, 0), source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE,
		           RANK_MASTER,
		           MPI_COMM_WORLD);
		lbm_comm_timers_stop(mesh_comm, TIMER_IO_GATHER_COMM);

		if (rank == RANK_MASTER) {
			lbm_comm_timers_start(mesh_comm, TIMER_IO_WRITE);
			for (int i = 0; i < comm_size; i++) {
				Mesh tmp_mesh;
				tmp_mesh.width = source_mesh->width;
				tmp_mesh.height = source_mesh->height;
				tmp_mesh.cells = temp->cells + i * source_mesh->width * source_mesh->height * DIRECTIONS;
				save_frame(fp, &tmp_mesh);
			}
			lbm_comm_timers_stop(mesh_comm, TIMER_IO_WRITE);
		}
	} else {
		/* Only 0 renders its local mesh */
		save_frame(fp, source_mesh);
	}

}

