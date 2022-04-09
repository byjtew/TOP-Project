/********************  HEADERS  *********************/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include "lbm_comm.h"

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


static inline double toMicroSeconds(double seconds) {
	return seconds * 1000000.0F;
}

double round_cell(const double *cell) {
	double sum = 0.0;
	for (int i = 0; i < DIRECTIONS; i++) sum += cell[i];
	return sum;
}

static inline void printInRed(const char *msg) {
	printf("\033[0;31m%s\033[0m", msg);
}

static inline void print_cell(lbm_mesh_cell_t cell) {
	double sum = round_cell(cell);
	if (isnan(sum))
		printInRed("   NaN    ");
	else if (isinf(sum))
		printInRed("   inf    ");
	else
		printf("%+.1e  ", round_cell(cell));
}

void print_mesh(Mesh *mesh) {
	int i, j;
	printf("\n");
	for (i = 0; i < 2; i++) {
		for (j = 0; j < mesh->width; j++) print_cell(Mesh_get_cell(mesh, j, i));
		printf("\n");
	}
	for (i = 2; i < mesh->height - 2; i++) {
		for (j = 0; j < 2; j++) print_cell(Mesh_get_cell(mesh, j, i));
		for (j = 2; j < mesh->width - 2; j++) printf("   ---    ");
		for (j = mesh->width - 2; j < mesh->width; j++) print_cell(Mesh_get_cell(mesh, j, i));
		printf("\n");
	}
	for (i = mesh->height - 2; i < mesh->height; i++) {
		for (j = 0; j < mesh->width; j++) print_cell(Mesh_get_cell(mesh, j, i));
		printf("\n");
	}
	printf("\n");
	fflush(stdout);
}

/**
 * @param borders | top-left corner | top row | top-right corner | right column | bottom-right corner | bottom row | bottom-left corner | left column
 * @param width mesh width
 * @param height mesh height
 */
void print_buffered_mesh(double **borders, int size) {
	printf("\n");
	for (int i = 0; i < size; i += DIRECTIONS)
		printf("%3.2lf ", round_cell(borders[i]));
	printf("\n");
}

double *top_send_buffer = NULL, *top_recv_buffer = NULL, *bot_send_buffer = NULL, *bot_recv_buffer = NULL;

void lbm_comm_ghost_exchange_init(lbm_comm_t *mesh_comm) {
	top_send_buffer = (double *) malloc(2 * sizeof(double) * DIRECTIONS * (mesh_comm->width - 2));
	top_recv_buffer = (double *) malloc(2 * sizeof(double) * DIRECTIONS * (mesh_comm->width - 2));
	bot_send_buffer = (double *) malloc(2 * sizeof(double) * DIRECTIONS * (mesh_comm->width - 2));
	bot_recv_buffer = (double *) malloc(2 * sizeof(double) * DIRECTIONS * (mesh_comm->width - 2));
}

void lbm_comm_ghost_exchange_release() {
	if (top_send_buffer != NULL)
		free(top_send_buffer);
	if (top_recv_buffer != NULL)
		free(top_recv_buffer);
	if (bot_send_buffer != NULL)
		free(bot_send_buffer);
	if (bot_recv_buffer != NULL)
		free(bot_recv_buffer);
}

/*******************  FUNCTION  *********************/
void lbm_comm_ghost_exchange(lbm_comm_t *mesh_comm, Mesh *mesh, int rank) {
	double timer;
	mesh_comm->current_request = 0;
	int used_top_recv = 0, used_bot_recv = 0;

	// region Horizontal
	if (rank == 0)
		timer = MPI_Wtime();
	lbm_comm_sync_ghosts_horizontal(mesh_comm, mesh, COMM_SEND, mesh_comm->right_id, mesh_comm->width - 2);
	lbm_comm_sync_ghosts_horizontal(mesh_comm, mesh, COMM_RECV, mesh_comm->right_id, mesh_comm->width - 1);
	lbm_comm_sync_ghosts_horizontal(mesh_comm, mesh, COMM_SEND, mesh_comm->left_id, 1);
	lbm_comm_sync_ghosts_horizontal(mesh_comm, mesh, COMM_RECV, mesh_comm->left_id, 0);
	if (rank == 0)
		fprintf(stderr, "Horizontal comms : %5.2lf\n", toMicroSeconds(MPI_Wtime() - timer));
	// endregion

	// region Vertical
	if (rank == 0)
		timer = MPI_Wtime();
	lbm_comm_sync_ghosts_vertical(mesh_comm, mesh, COMM_SEND, mesh_comm->bottom_id, mesh_comm->height - 2,
	                              bot_send_buffer);
	used_bot_recv = lbm_comm_sync_ghosts_vertical(mesh_comm, mesh, COMM_RECV, mesh_comm->bottom_id, mesh_comm->height - 1,
	                                              bot_recv_buffer);
	lbm_comm_sync_ghosts_vertical(mesh_comm, mesh, COMM_SEND, mesh_comm->top_id, 1, top_send_buffer);
	used_top_recv = lbm_comm_sync_ghosts_vertical(mesh_comm, mesh, COMM_RECV, mesh_comm->top_id, 0, top_recv_buffer);
	if (rank == 0)
		fprintf(stderr, "Vertical comms : %5.2lf\n", toMicroSeconds(MPI_Wtime() - timer));
	// endregion


	// region Diagonal
	if (rank == 0)
		timer = MPI_Wtime();
	//top left
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_SEND, mesh_comm->corner_id[CORNER_TOP_LEFT], 1, 1);
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_RECV, mesh_comm->corner_id[CORNER_BOTTOM_RIGHT],
	                              mesh_comm->width - 1, mesh_comm->height - 1);
	if (rank == 0)
		fprintf(stderr, "Top left comms : %5.2lf\n", toMicroSeconds(MPI_Wtime() - timer));


	if (rank == 0)
		timer = MPI_Wtime();
	//bottom left
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_SEND, mesh_comm->corner_id[CORNER_BOTTOM_LEFT], 1,
	                              mesh_comm->height - 2);
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_RECV, mesh_comm->corner_id[CORNER_TOP_RIGHT],
	                              mesh_comm->width - 1,
	                              0);
	if (rank == 0)
		fprintf(stderr, "Bottom left comms : %5.2lf\n", toMicroSeconds(MPI_Wtime() - timer));


	if (rank == 0)
		timer = MPI_Wtime();
	//top right
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_SEND, mesh_comm->corner_id[CORNER_TOP_RIGHT],
	                              mesh_comm->width - 2, 1);
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_RECV, mesh_comm->corner_id[CORNER_BOTTOM_LEFT], 0,
	                              mesh_comm->height - 1);
	if (rank == 0)
		fprintf(stderr, "Top right comms : %5.2lf\n", toMicroSeconds(MPI_Wtime() - timer));

	if (rank == 0)
		timer = MPI_Wtime();
	//bottom right
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_SEND, mesh_comm->corner_id[CORNER_BOTTOM_RIGHT],
	                              mesh_comm->width - 2, mesh_comm->height - 2);
	lbm_comm_sync_ghosts_diagonal(mesh_comm, mesh, COMM_RECV, mesh_comm->corner_id[CORNER_TOP_LEFT], 0, 0);
	if (rank == 0)
		fprintf(stderr, "Bottom right comms : %5.2lf\n", toMicroSeconds(MPI_Wtime() - timer));
	// endregion

	assert(mesh_comm->current_request <= mesh_comm->max_requests);
	MPI_Waitall(mesh_comm->current_request, mesh_comm->requests, mesh_comm->statuses);

	if (top_recv_buffer != NULL && used_top_recv > 0)
		for (int x = 1; x < mesh->width - 2; x++)
			memcpy(Mesh_get_cell(mesh, x, 0), top_recv_buffer + (x - 1) * DIRECTIONS, sizeof(double) * DIRECTIONS);
	if (bot_recv_buffer != NULL && used_bot_recv > 0)
		for (int x = 1; x < mesh->width - 2; x++)
			memcpy(Mesh_get_cell(mesh, x, mesh_comm->height - 1), bot_recv_buffer + (x - 1) * DIRECTIONS,
			       sizeof(double) * DIRECTIONS);
}

/*******************  FUNCTION  *********************/
/**
 * Rendu du mesh en effectuant une réduction a 0
 * @param mesh_comm MeshComm à utiliser
 * @param temp Mesh a utiliser pour stocker les segments
**/
void save_frame_all_domain(FILE *fp, Mesh *source_mesh, Mesh *temp) {
	// Todo: Switch to a Gather
	//vars
	int i = 0;
	int comm_size, rank;
	MPI_Status status;

	//get rank and comm size
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* If whe have more than one process */
	if (1 < comm_size) {
		if (rank == 0) {
			/* Rank 0 renders its local Mesh */
			save_frame(fp, source_mesh);
			/* Rank 0 receives & render other processes meshes */
			for (i = 1; i < comm_size; i++) {
				MPI_Recv(temp->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
				         &status);
				save_frame(fp, temp);
			}
		} else {
			/* All other ranks send their local mesh */
			MPI_Send(source_mesh->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE, 0, 0,
			         MPI_COMM_WORLD);
		}
	} else {
		/* Only 0 renders its local mesh */
		save_frame(fp, source_mesh);
	}

}

