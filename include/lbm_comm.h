#ifndef MESH_COMM_H
#define MESH_COMM_H

/********************  HEADERS  *********************/
#include <mpi.h>
#include <stdlib.h>
#include "lbm_struct.h"

/*******************  DEFINITIONS  ******************/
/** Definition de l'ID du processus maître. **/
#define RANK_MASTER 0
#define TIMER_MESH_SYNC 0
#define TIMER_OUTPUT_GATHER 1
#define TIMER_SEND_BUFFER_CREATE 2
#define TIMER_RECV_BUFFER_EMPLACE 3

#define NB_USED_TIMER 4

/*********************  ENUM  ***********************/
/**
 * Definition des différents type de cellule pour savoir quel traitement y appliquer
 * lors du calcul.
**/
typedef enum lbm_corner_pos_e {
		CORNER_TOP_LEFT = 0,
		CORNER_TOP_RIGHT = 1,
		CORNER_BOTTOM_LEFT = 2,
		CORNER_BOTTOM_RIGHT = 3,
} lbm_corner_pos_t;

/*********************  ENUM  ***********************/
typedef enum lbm_comm_type_e {
		COMM_SEND,
		COMM_RECV
} lbm_comm_type_t;

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
		int right_id;
		/** Id du voisin de gauche, -1 si aucun. **/
		int left_id;
		int top_id;
		int bottom_id;
		int corner_id[4];
		/** Requète asynchrone en cours
		 * @size width*2*2(rows) + 2*4(corners) + 2*2(columns)
		 * **/
		MPI_Request *requests;
		int max_requests;
		int current_request;
		MPI_Status *statuses;
		lbm_mesh_cell_t buffer;

		// Measurements for the communication
#define NB_TIMERS 64
		double *timers[NB_TIMERS];
		int current_timer[NB_TIMERS];

		// Graph MPI communicator
		MPI_Comm comm_graph;
		int nb_per_neigh[8];
		int displ_per_neigh[8];
		double *send_borders;
		double *recv_borders;
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
void lbm_comm_sync_ghosts_wait(lbm_comm_t *mesh);

void lbm_comm_ghost_exchange_init(lbm_comm_t *mesh_comm);

void lbm_comm_ghost_exchange_release();

void lbm_comm_ghost_exchange(lbm_comm_t *mesh_comm, Mesh *mesh, int rank);

/*******************  FUNCTION  *********************/
void save_frame_all_domain(FILE *fp, Mesh *source_mesh, Mesh *temp, lbm_comm_t *mesh_comm);

#endif
