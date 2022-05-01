#ifndef LBM_PHYS_H
#define LBM_PHYS_H

/********************  HEADERS  *********************/
#include "lbm_struct.h"
#include "lbm_comm.h"

/********************** CONSTS **********************/
/**
 * Poids utilisé pour compenser les différentes de longueur des 9 vecteurs directions.
**/
#if DIRECTIONS == 9
static const double equil_weight[DIRECTIONS] = {
		4.0 / 9.0,
		1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
		1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};
//opposite directions, for bounce back implementation
static const int opposite_of[DIRECTIONS] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
/**
 * Definitions des 9 vecteurs de base utilisé pour discrétiser les directions sur chaque mailles.
**/
#if DIMENSIONS == 2
static const IntVector direction_matrix[DIRECTIONS] = {
		{+0, +0},
		{+1, +0},
		{+0, +1},
		{-1, +0},
		{+0, -1},
		{+1, +1},
		{-1, +1},
		{-1, -1},
		{+1, -1}
};
#else
#error Need to defined adapted direction matrix.
#endif
#else
#error Need to defined adapted equibirium distribution function
#endif


/*******************  FUNCTION  *********************/
//helper
double get_vect_norme_2_ividv(const IntVector vect1, const Vector vect2);

double get_vect_norme_2_dvdv(const Vector vect1, const Vector vect2);

double get_cell_density(const lbm_mesh_cell_t cell);

void get_cell_velocity(Vector v, lbm_mesh_cell_t cell, double cell_density);

double helper_compute_poiseuille(int i, int size);

/*******************  FUNCTION  *********************/
//collistion
double compute_equilibrium_profile(Vector velocity, double density, int direction);

void compute_cell_collision(lbm_mesh_cell_t cell_out, lbm_mesh_cell_t cell_in);

/*******************  FUNCTION  *********************/
//limit conditions
void compute_bounce_back(lbm_mesh_cell_t cell);

void compute_inflow_zou_he_poiseuille_distr(const Mesh *mesh, lbm_mesh_cell_t cell, int id_y);

void compute_outflow_zou_he_const_density(lbm_mesh_cell_t mesh);

/*******************  FUNCTION  *********************/
//main functions
void special_cells(Mesh *mesh, lbm_mesh_type_t *mesh_type, const lbm_comm_t *mesh_comm);

void collision(Mesh *mesh_out, const Mesh *mesh_in);

void propagation(Mesh *mesh_out, const Mesh *mesh_in);

#endif
