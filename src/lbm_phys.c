/********************  HEADERS  *********************/
#include <assert.h>
#include "lbm_config.h"
#include "lbm_struct.h"
#include "lbm_phys.h"
#include "lbm_comm.h"

/*******************  FUNCTION  *********************/
/**
 * Renvoie le résultat du produit des deux vecteurs passés en paramêtre.
**/
double get_vect_norme_2_ivdv(const IntVector vect1, const Vector vect2) {
	//vars
	int k;
	double res = 0.0;

	//loop on dimensions
	for (k = 0; k < DIMENSIONS; k++)
		res += vect1[k] * vect2[k];

	return res;
}

double get_vect_norme_2_dvdv(const Vector vect1, const Vector vect2) {
	//vars
	int k;
	double res = 0.0;

	//loop on dimensions
	for (k = 0; k < DIMENSIONS; k++)
		res += vect1[k] * vect2[k];

	return res;
}

/*******************  FUNCTION  *********************/
/**
 * Calcule la densité macroscopiques de la cellule en sommant ses DIRECTIONS
 * densités microscopiques.
**/
double get_cell_density(lbm_mesh_cell_t cell) {
	//vars
	int k;
	double res = 0.0;

	//errors
	assert(cell != NULL);

	//loop on directions
	for (k = 0; k < DIRECTIONS; k++)
		res += cell[k];

	//return res
	return res;
}

/*******************  FUNCTION  *********************/
/**
 * Calcule la vitesse macroscopiques de la cellule en sommant ses DIRECTIONS
 * densités microscopiques.
 * @param cell_density Densité macroscopique de la cellules.
**/
void get_cell_velocity(Vector v, lbm_mesh_cell_t const cell, double cell_density) {
	//vars
	int k, d;

	//errors
	assert(v != NULL);
	assert(cell != NULL);

	//loop on all dimensions
	for (d = 0; d < DIMENSIONS; d++) {
		//reset value
		v[d] = 0.0;

		//sum all directions
		for (k = 0; k < DIRECTIONS; k++)
			v[d] += cell[k] * direction_matrix[k][d];

		//normalize
		v[d] = v[d] / cell_density;
	}
}

/*******************  FUNCTION  *********************/
/**
 * Calcule le profile de densité microscopique à l'équilibre pour une direction donnée.
 * C'est ici qu'intervient un dérivé de navier-stokes.
 * @param velocity Vitesse macroscopique du fluide sur la maille.
 * @param density Densité macroscopique du fluide sur la maille.
 * @param direction Direction pour laquelle calculer l'état d'équilibre.
**/
double compute_equilibrium_profile(Vector velocity, double density, int direction) {
	//vars
	double v2;
	double p;
	double p2;
	double feq;

	//velocity norme 2 (v * v)
	v2 = get_vect_norme_2_dvdv(velocity, velocity);

	//calc e_i * v_i / c
	p = get_vect_norme_2_ivdv(direction_matrix[direction], velocity);
	p2 = p * p;

	//terms without density and direction weight
	feq = 1.0
	      + (3.0 * p)
	      + ((9.0 / 2.0) * p2)
	      - ((3.0 / 2.0) * v2);

	//mult all by density and direction weight
	feq *= equil_weight[direction] * density;

	return feq;
}

/*******************  FUNCTION  *********************/
/**
 * Calcule le vecteur de collision entre les fluides de chacune des directions.
**/
void compute_cell_collision(lbm_mesh_cell_t cell_out, lbm_mesh_cell_t const cell_in) {
	//vars
	int k;
	double density;
	Vector v;
	double feq;

	//compute macroscopic values
	density = get_cell_density(cell_in);
	get_cell_velocity(v, cell_in, density);

	//loop on microscopic directions
	for (k = 0; k < DIRECTIONS; k++) {
		//compute f at equilibr.
		feq = compute_equilibrium_profile(v, density, k);
		//compute fout
		cell_out[k] = cell_in[k] - RELAX_PARAMETER * (cell_in[k] - feq);
	}
}

/*******************  FUNCTION  *********************/
/**
 * Applique une reflexion sur les différentes directions pour simuler la présence d'un solide.
**/
void compute_bounce_back(lbm_mesh_cell_t cell) {
	/*//vars
	int k;
	double tmp[DIRECTIONS];*/

	//compute bounce back
/*	for (k = 0; k < DIRECTIONS; k++)
		tmp[k] = cell[opposite_of[k]];
	*/
	// nb: Equivalent to tmp = cell[0,3,4,1,2,7,8,5,6]
	//                         indx[0,1,2,3,4,5,6,7,8]
	double swap_tmp = cell[4];
	cell[4] = cell[2];
	cell[2] = swap_tmp;
	swap_tmp = cell[3];
	cell[3] = cell[1];
	cell[1] = swap_tmp;
	swap_tmp = cell[5];
	cell[5] = cell[7];
	cell[7] = swap_tmp;
	swap_tmp = cell[6];
	cell[6] = cell[8];
	cell[8] = swap_tmp;


	//compute bounce back
/*	for (k = 0; k < DIRECTIONS; k++)
		cell[k] = tmp[k];*/
}

/*******************  FUNCTION  *********************/
/**
 * Fournit la vitesse de poiseuille pour une position donnée en considérant un tube de taille donnée.
 * @param i Position pour laquelle on cherche la vitesse.
 * @param size diamètre du tube.
**/
double helper_compute_poiseuille(int i, int size) {
	double y = (double) (i - 1);
	double L = (double) (size - 1);
	return 4.0 * INFLOW_MAX_VELOCITY / (L * L) * (L * y - y * y);
}

/*******************  FUNCTION  *********************/
/**
 * Applique la méthode de Zou/He pour simler un fluidre entrant dans le domain de gauche vers la droite sur une
 * interface verticale. Le profile de vitesse du fluide entrant suit une distribution de poiseuille.
 * @param mesh Maillage considéré (surtout pour avoir la hauteur.)
 * @param cell Maille à mettre à jour.
 * @param id_y Position en y de la cellule pour savoir comment calculer la vitesse de poiseuille.
**/
void compute_inflow_zou_he_poiseuille_distr(const Mesh *mesh, lbm_mesh_cell_t cell, int id_y) {
	//vars
	double v;
	double density;

	//errors
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif

	//set macroscopic fluide info
	//poiseuille distr on X and null on Y
	//we just want the norm, so v = v_x
	v = helper_compute_poiseuille(id_y, mesh->height);

	//compute rho from u and inner flow on surface
	density = (cell[0] + cell[2] + cell[4] + 2 * (cell[3] + cell[6] + cell[7])) / (1.0 - v);

	//now compute unknown microscopic values
	cell[1] = cell[3];// + (2.0/3.0) * density * v_y <--- no velocity on Y so v_y = 0
	cell[5] = cell[7] - (1.0 / 2.0) * (cell[2] - cell[4])
	          + (1.0 / 6.0) * (density * v);
	//+ (1.0/2.0) * density * v_y    <--- no velocity on Y so v_y = 0
	cell[8] = cell[6] + (1.0 / 2.0) * (cell[2] - cell[4])
	          + (1.0 / 6.0) * (density * v);
	//- (1.0/2.0) * density * v_y    <--- no velocity on Y so v_y = 0

	//no need to copy already known one as the value will be "loss" in the wall at propagatation time
}

/*******************  FUNCTION  *********************/
/**
 * Applique la méthode de Zou/He pour simuler un fluide sortant du domaine de gauche vers la droite sur une
 * interface verticale. La condition appliquée pour construire l'équation est le maintient d'un gradiant de densité
 * nulle à l'interface.
 * @param mesh Maillage considéré (surtout pour avoir la hauteur.)
 * @param cell Maille à mettre à jour
 * @param id_y Position en y de la cellule pour savoir comment calculer la vitesse de poiseuille.
**/
void compute_outflow_zou_he_const_density(lbm_mesh_cell_t cell) {
	//vars
	const double density = 1.0;
	double v;

	//errors
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif

	//compute macroscopic v depeding on inner flow going onto the wall
	v = -1.0 + (1.0 / density) * (cell[0] + cell[2] + cell[4] + 2 * (cell[1] + cell[5] + cell[8]));

	//now can compute unknown microscopic values
	cell[3] = cell[1] - (2.0 / 3.0) * density * v;
	cell[7] = cell[5] + (1.0 / 2.0) * (cell[2] - cell[4])
	          //- (1.0/2.0) * (density * v_y)    <--- no velocity on Y so v_y = 0
	          - (1.0 / 6.0) * (density * v);
	cell[6] = cell[8] + (1.0 / 2.0) * (cell[4] - cell[2])
	          //+ (1.0/2.0) * (density * v_y)    <--- no velocity on Y so v_y = 0
	          - (1.0 / 6.0) * (density * v);
}

/*******************  FUNCTION  *********************/
/**
 * Applique les actions spéciales liées aux conditions de bords ou aux réflexions sur l'obstacle.
**/
void special_cells(Mesh *mesh, lbm_mesh_type_t *mesh_type, const lbm_comm_t *mesh_comm) {

// nb: do not parallelize this function under 10 processes
//#pragma omp parallel for default(none) shared(mesh, mesh_type, mesh_comm) schedule(runtime)
	for (int i = 1; i < mesh->width - 1; i++) {
		for (int j = 1; j < mesh->height - 1; j++) {
			switch (*(lbm_cell_type_t_get_cell(mesh_type, i, j))) {
				case CELL_FUILD:
					break;
				case CELL_BOUNCE_BACK:
					compute_bounce_back(Mesh_get_cell(mesh, i, j));
					break;
				case CELL_LEFT_IN:
					compute_inflow_zou_he_poiseuille_distr(mesh, Mesh_get_cell(mesh, i, j), j + mesh_comm->y);
					break;
				case CELL_RIGHT_OUT:
					compute_outflow_zou_he_const_density(Mesh_get_cell(mesh, i, j));
					break;
			}
		}
	}
}

/*******************  FUNCTION  *********************/
/**
 * Calcule les collisions sur chacune des cellules.
 * @param mesh Maillage sur lequel appliquer le calcul.
**/
void collision(Mesh *mesh_out, const Mesh *mesh_in) {
	//loop on all inner cells

#pragma omp parallel default(none) shared(mesh_in, mesh_out)
	{
#pragma omp for
		for (int i = 1; i < mesh_in->width - 1; i++)
			for (int j = 1; j < mesh_in->height - 1; j++)
				compute_cell_collision(Mesh_get_cell(mesh_out, i, j), Mesh_get_cell(mesh_in, i, j));
	}

}

/*******************  FUNCTION  *********************/
/**
 * Propagation des densités vers les mailles voisines.
 * @param mesh_out Maillage de sortie.
 * @param mesh_in Maillage d'entrée (ne doivent pas être les mêmes).
**/
void propagation(Mesh *mesh_out, const Mesh *mesh_in) {
	//loop on all cells
#pragma omp parallel for default(none) shared(mesh_in, mesh_out, direction_matrix)
	for (int i = 0; i < mesh_out->width; i++) {
		for (int j = 0; j < mesh_out->height; j++) {
			//for all direction
			for (int k = 0; k < DIRECTIONS; k++) {
				//compute destination point
				int ii = (i + direction_matrix[k][0]);
				int jj = (j + direction_matrix[k][1]);
				//propagate to neighboor nodes
				if ((ii >= 0 && ii < mesh_out->width) && (jj >= 0 && jj < mesh_out->height))
					Mesh_get_cell(mesh_out, ii, jj)[k] = Mesh_get_cell(mesh_in, i, j)[k];
			}
		}
	}
}
