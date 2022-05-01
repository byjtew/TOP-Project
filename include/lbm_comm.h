#ifndef MESH_COMM_H
#define MESH_COMM_H

/********************  HEADERS  *********************/
#include <mpi.h>
#include <stdlib.h>
#include "lbm_struct.h"

/*******************  DEFINITIONS  ******************/
/** Definition de l'ID du processus maître. **/
#define RANK_MASTER 0

/** Modes de communication. **/
#define MESH_SYNC_GRAPH 0
#define MESH_SYNC_UNIT_ASYNCHRONOUS 1
#define MESH_SYNC_CART 2
#define MESH_SYNC_UNIT_SYNCHRONOUS 3

#define MESH_SYNC_MODE MESH_SYNC_UNIT_ASYNCHRONOUS

/** Timers **/
#define TIMER_MESH_SYNC_COMM 0
#define TIMER_MESH_PRE_SYNC 1
#define TIMER_MESH_POST_SYNC 2
#define TIMER_IO_GATHER_COMM 3
#define TIMER_IO_WRITE 4
#define TIMER_SPECIAL_CELLS 5
#define TIMER_COLLISION 6
#define TIMER_PROPAGATION 7
#define TIMER_GHOST_EXCHANGE_TOTAL 8

#define NB_USED_TIMER 9

/********************  STRUCT  **********************/
/**
 * Structure utilisée pour stoquer les informations relatives aux communications.
**/
typedef struct lbm_comm_t_s {
		/** Position de la maille locale dans le maillage global (origine). **/
		int x;
		int y;
		/** Taille de la maille locale. **/
		int width;
		int height;
		int nb_x;
		int nb_y;
		/** Id du voisin de droite, -1 si aucun. **/
		/** Id du voisin de gauche, -1 si aucun. **/
		int left_id;
		int right_id;
		/** Requète asynchrone en cours
		 * @size width*2*2(rows) + 2*4(corners) + 2*2(columns)
		 * **/
		MPI_Request *requests;
		int max_requests;
		int current_request;

		// Measurements for the communication
#define NB_TIMERS 16
		double *timers[NB_TIMERS];
		int current_timer[NB_TIMERS];

		// Graph MPI communicator
		MPI_Comm comm_graph;
		int nb_per_neigh[2];
		int send_displ[2];
		int recv_displ[2];

		// IO synchronization
} lbm_comm_t;

/*******************  FUNCTION  *********************/
static inline int lbm_comm_width(lbm_comm_t *mc) {
	return mc->width;
}

/*******************  FUNCTION  *********************/
static inline int lbm_comm_height(lbm_comm_t *mc) {
	return mc->height;
}

/*******************  FUNCTION  *********************/
void lbm_comm_init(lbm_comm_t *mesh, int rank, int comm_size, int width, int height);

void lbm_comm_release(lbm_comm_t *mesh);

void lbm_comm_timers_start(lbm_comm_t *mesh_comm, int id);

void lbm_comm_timers_stop(lbm_comm_t *mesh_comm, int id);

void lbm_comm_print(lbm_comm_t *mesh);

/*******************  FUNCTION  *********************/
void lbm_comm_ghost_exchange(lbm_comm_t *mesh_comm, Mesh *mesh, int rank);

/*******************  FUNCTION  *********************/
void save_frame_all_domain(FILE *fp, Mesh *source_mesh, Mesh *temp, lbm_comm_t *mesh_comm);

#endif
