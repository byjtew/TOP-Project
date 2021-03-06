#ifndef LBM_STRUCT_H
#define LBM_STRUCT_H

/********************  HEADERS  *********************/
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include "lbm_config.h"

/********************** TYPEDEF *********************/
/**
 * Une cellule est un tableau de DIRECTIONS doubles pour stoquer les
 * probabilités microscopiques (f_i).
**/
typedef double *lbm_mesh_cell_t;
/** Représentation d'un vecteur pour la manipulation des vitesses macroscopiques. **/
typedef double Vector[DIMENSIONS];
typedef int IntVector[DIMENSIONS];

/********************  STRUCT  **********************/
/**
 * Definit un maillage pour le domaine local. Ce maillage contient une bordure d'une cellule
 * contenant les mailles fantômes.
**/
typedef struct Mesh {
		/** Cellules du maillages (MESH_WIDTH * MESH_HEIGHT * DIRECTIONS). **/
		lbm_mesh_cell_t cells;
		/** Largeur du maillage local (mailles fantome comprises). **/
		int width;
		/** Largeur du maillage local (mailles fantome comprises). **/
		int height;
} Mesh;

/*********************  ENUM  ***********************/
/**
 * Definition des différents type de cellule pour savoir quel traitement y appliquer
 * lors du calcul.
**/
typedef enum lbm_cell_type_e {
		/** Cellule de fluide standard, uniquement application des collisions. **/
		CELL_FUILD,
		/** Cellules de l'obstacle ou des bordure supérieures et inférieurs. Application de réflexion. **/
		CELL_BOUNCE_BACK,
		/** Cellule de la paroie d'entrée. Application de Zou/He avec V fixé. **/
		CELL_LEFT_IN,
		/** Cellule de la paroie de sortie. Application de Zou/He avec gradiant de densité constant. **/
		CELL_RIGHT_OUT
} lbm_cell_type_t;

/********************  STRUCT  **********************/
/**
 * Tableau maitnenant les informations de type pour les cellules.
**/
typedef struct lbm_mesh_type_s {
		/** Type des cellules du maillages (MESH_WIDTH * MESH_HEIGHT). **/
		lbm_cell_type_t *types;
		/** Largeur du maillage local (mailles fantome comprises). **/
		int width;
		/** Largeur du maillage local (mailles fantome comprises). **/
		int height;
} lbm_mesh_type_t;

/********************  STRUCT  **********************/
/** Structure des en-têtes utilisée dans le fichier de sortie. **/
typedef struct lbm_file_header_s {
		/** Pour validation du format du fichier. **/
		uint32_t magick;
		/** Taille totale du maillage simulé (hors mailles fantômes). **/
		uint32_t mesh_width;
		/** Taille totale du maillage simulé (hors mailles fantômes). **/
		uint32_t mesh_height;
		/** Number of vertical lines. **/
		uint32_t lines;
} lbm_file_header_t;

/********************  STRUCT  **********************/
/** Une entrée du fichier, avec les deux grandeurs macroscopiques. **/
typedef struct lbm_file_entry_s {
		double v;
		double density;
} lbm_file_entry_t;

/********************  STRUCT  **********************/
/** Pour la lecture du fichier de sortie. **/
typedef struct lbm_data_file_s {
		FILE *fp;
		lbm_file_header_t header;
		lbm_file_entry_t *entries;
} lbm_data_file_t;


/*******************  FUNCTION  *********************/
void Mesh_init(Mesh *mesh, int width, int height);

void Mesh_release(Mesh *mesh);

/*******************  FUNCTION  *********************/
void lbm_mesh_type_t_init(lbm_mesh_type_t *mesh, int width, int height);

void lbm_mesh_type_t_release(lbm_mesh_type_t *mesh);

/*******************  FUNCTION  *********************/
void save_frame(FILE *fp, const Mesh *mesh);

/*******************  FUNCTION  *********************/
void fatal(const char *message);

/*******************  FUNCTION  *********************/
/** Mesh is in COL_MAJOR order.
 *          --- x --->
 * -------------------------------
 * |	|  |  |  |  |  |  |  |  |  |   |
 * |  |  |  |  |  |  |  |  |  |  |   |
 * |  |  |  |  |  |  |  |  |  |  |   y
 * |  |  |  |  |  |  |  |  |  |  |   |
 * |  |  |  |  |  |  |  |  |  |  |   |
 * |  |  |  |  |  |  |  |  |  |  |  \/
 * -------------------------------
 */

/** Mesh is in ROW_MAJOR order.
 *          --- width (x) --->
 *  -------------------------------
 *  |  |  |  |  |  |  |  |  |  |  |   |
 *  -------------------------------   |
 *  |  |  |  |  |  |  |  |  |  |  |   |
 *  -------------------------------   height (y)
 *  |  |  |  |  |  |  |  |  |  |  |   |
 *  ------------------------------- 	|
 *  |  |  |  |  |  |  |  |  |  |  |	 \/
 *  -------------------------------
 */

/**
 * Fonction à utiliser pour récupérer une cellule du maillage en fonction de ses coordonnées.
**/
static inline lbm_mesh_cell_t Mesh_get_cell(const Mesh *mesh, int x, int y) {
	// Col major: int idx = (x * mesh->height + y) * DIRECTIONS;
	int idx = (x * mesh->height + y) * DIRECTIONS;
	// Row major: int idx = (y * mesh->width + x) * DIRECTIONS;
	/*assert(mesh != NULL);
	assert(mesh->cells != NULL);
	assert(x >= 0);
	assert(x < mesh->width);
	assert(y >= 0);
	assert(y < mesh->height);
	assert(mesh->height * mesh->width * DIRECTIONS > idx);*/
	return &(mesh->cells[idx]);
}

/*******************  FUNCTION  *********************/
/**
 * Fonction à utiliser pour récupérer une ligne (suivant y, x = 1) du maillage en fonction de ses coordonnées.
**/
/*static inline lbm_mesh_cell_t Mesh_get_row(const Mesh *mesh, int y) {
	return Mesh_get_cell(mesh, 1, y);
}*/
static inline lbm_mesh_cell_t Mesh_get_col(const Mesh *mesh, int x) {
	return Mesh_get_cell(mesh, x, 1);
}

/*******************  FUNCTION  *********************/
/**
 * Fonction à utiliser pour récupérer un pointeur sur le type d'une cellule du maillage en fonction de ses coordonnées.
**/
static inline lbm_cell_type_t *lbm_cell_type_t_get_cell(const lbm_mesh_type_t *meshtype, int x, int y) {
	return &meshtype->types[x * meshtype->height + y];
}

#endif //LBM_STRUCT_H
