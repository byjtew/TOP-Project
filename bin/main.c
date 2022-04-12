/********************  HEADERS  *********************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <omp.h>
#include "lbm_config.h"
#include "lbm_struct.h"
#include "lbm_phys.h"
#include "lbm_init.h"
#include "lbm_comm.h"

#ifdef RELEASE_MODE
#pragma message "Release mode"
#else
#pragma message "Debug mode"
#endif

static int win_mutex = -1;
static int rank = -1;

/*******************  FUNCTION  *********************/


// Printf overload
void mpi_put(const char *format, ...) {
	int offset = 0;
	char buffer[1024];
	va_list args;
	if (rank == 0) // orange
		offset += sprintf(buffer + offset, "\033[0;33m");
	else // white
		offset += sprintf(buffer + offset, "\033[0;37m");

	offset += sprintf(buffer + offset, "[P%d]: ", rank);
	va_start(args, format);
	offset += vsprintf(buffer + offset, format, args);
	va_end(args);
	offset += sprintf(buffer + offset, "\n");
	// Remove color code
	offset += sprintf(buffer + offset, "\033[0m");
#ifndef RELEASE_MODE
	// Lock the mutex
	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_mutex);
	printf(buffer, args);
	// Unlock the mutex
	MPI_Win_unlock(0, win_mutex);
#endif
	fflush(stdout);
}

// Create a single mutex for all the processes, accessible by MPI_Win_lock
void create_mutex(void) {
	MPI_Win_create(&rank, sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_mutex);
	MPI_Win_fence(0, win_mutex);
	mpi_put("Creating mutex");
}

/**
 * Ecrit l'en-tête du fichier de sortie. Cet en-tête sert essentiellement à fournir les informations
 * de taille du maillage pour les chargements.
 * @param fp Descripteur de fichier à utiliser pour l'écriture.
**/
void write_file_header(FILE *fp, lbm_comm_t *mesh_comm) {
	//setup header values
	lbm_file_header_t header;
	header.magick = RESULT_MAGICK;
	header.mesh_height = MESH_HEIGHT;
	header.mesh_width = MESH_WIDTH;
	header.lines = mesh_comm->nb_y;

	// print header in stdout
	mpi_put("Writing header to file");
	mpi_put("\tMAGICK: %d", header.magick);
	mpi_put("\tMESH_HEIGHT: %d", header.mesh_height);
	mpi_put("\tMESH_WIDTH: %d", header.mesh_width);
	mpi_put("\tLINES: %d", header.lines);

	//write header
	size_t written = fwrite(&header, sizeof(lbm_file_header_t), 1, fp);

	if (written <= 0) {
		mpi_put("Error writing file header, written %d bytes, expected %d", written, sizeof(lbm_file_header_t));
		exit(1);
	}
}

/*******************  FUNCTION  *********************/
FILE *open_output_file(lbm_comm_t *mesh_comm) {
	//vars
	FILE *fp;

	//check if empty filename => so noout
	if (RESULT_FILENAME == NULL)
		return NULL;
	else
		mpi_put("Opening output file: %s", RESULT_FILENAME);

	//open result file
	fp = fopen(RESULT_FILENAME, "w");

	//errors
	if (fp == NULL) {
		fprintf(stderr, "Error opening output file %s\n", RESULT_FILENAME);
		perror(RESULT_FILENAME);
		abort();
	}

	//write header
	write_file_header(fp, mesh_comm);

	return fp;
}

void close_file(FILE *fp) {
	//close file
	fclose(fp);
}

/*******************  FUNCTION  *********************/
/**
 * Sauvegarde le résultat d'une étape de calcul. Cette fonction peu être appelée plusieurs fois
 * lors d'une sauvegarde en MPI sur plusieurs processus pour sauvegarder les un après-les autres
 * chacun des domaines.
 * Ne sont écrit que les vitesses et densités macroscopiques sous forme de flotant simple.
 * @param fp Descripteur de fichier à utiliser pour l'écriture.
 * @param mesh Domaine à sauvegarder.
**/
void save_frame(FILE *fp, const Mesh *mesh) {
	//write buffer to write float instead of double
	lbm_file_entry_t buffer[WRITE_BUFFER_ENTRIES];
	int i, j, cnt;
	double density;
	Vector v;
	double norm;

	//loop on all values
	cnt = 0;
	for (i = 1; i < mesh->width - 1; i++) {
		for (j = 1; j < mesh->height - 1; j++) {
			//compute macrospic values
			density = get_cell_density(Mesh_get_cell(mesh, i, j));
			get_cell_velocity(v, Mesh_get_cell(mesh, i, j), density);
			norm = sqrt(get_vect_norme_2(v, v));

			//fill buffer
			buffer[cnt].density = density;
			buffer[cnt].v = norm;
			cnt++;

			//errors
			assert(cnt <= WRITE_BUFFER_ENTRIES);

			//flush buffer if full
			if (cnt == WRITE_BUFFER_ENTRIES) {
				fwrite(buffer, sizeof(lbm_file_entry_t), cnt, fp);
				cnt = 0;
			}
		}
	}

	//final flush
	if (cnt != 0)
		fwrite(buffer, sizeof(lbm_file_entry_t), cnt, fp);
}

static inline double toMicroSeconds(double seconds) {
	return seconds * 1000000.0F;
}

/*******************  FUNCTION  *********************/
int main(int argc, char *argv[]) {
	//vars
	Mesh mesh;
	Mesh temp;
	// Mesh large enough to hold everyone's mesh (only for MASTER)
	Mesh temp_render;
	lbm_mesh_type_t mesh_type;
	lbm_comm_t mesh_comm;
	int i, comm_size;
	FILE *fp = NULL;
	const char *config_filename = NULL;

	// Init MPI and get current rank and communicator size.
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	create_mutex();

	//get config filename
	if (argc >= 2)
		config_filename = argv[1];
	else
		config_filename = "config.txt";

	//load config file and display it on master node
	load_config(config_filename);
	if (rank == RANK_MASTER)
		print_config();

	//init structures, allocate memory...
	lbm_comm_init(&mesh_comm, rank, comm_size, MESH_WIDTH, MESH_HEIGHT);
	Mesh_init(&mesh, lbm_comm_width(&mesh_comm), lbm_comm_height(&mesh_comm));
	Mesh_init(&temp, lbm_comm_width(&mesh_comm), lbm_comm_height(&mesh_comm));
	if (rank == RANK_MASTER) Mesh_init(&temp_render, comm_size * lbm_comm_width(&mesh_comm), lbm_comm_height(&mesh_comm));
	lbm_mesh_type_t_init(&mesh_type, lbm_comm_width(&mesh_comm), lbm_comm_height(&mesh_comm));

	//master open the output file
	if (rank == RANK_MASTER)
		fp = open_output_file(&mesh_comm);

	//setup initial conditions on mesh
	setup_init_state(&mesh, &mesh_type, &mesh_comm);
	setup_init_state(&temp, &mesh_type, &mesh_comm);


	/*sleep(rank);
	fprintf(stderr, "Rank: %d\n", rank);
	fprintf(stderr, "Comm size: %d x %d\n", mesh_comm.width, mesh_comm.height);
	fprintf(stderr, "Mesh size: %d x %d\n", mesh.width, mesh.height);
	MPI_Barrier(MPI_COMM_WORLD);*/

	//write initial condition in output file
	if (lbm_gbl_config.output_filename != NULL)
		save_frame_all_domain(fp, &mesh, &temp_render, &mesh_comm);

	omp_set_dynamic(1);

	omp_set_num_threads(4);

	//time steps
	double total_time = 0;
	double iteration_timer = 0;
	double sub_timer = 0;
	for (i = 1; i <= ITERATIONS; i++) {
		//print progress
		if (rank == RANK_MASTER)
			mpi_put("Progress [%5d / %5d]", i, ITERATIONS);

		if (rank == RANK_MASTER)
			iteration_timer = MPI_Wtime();


		if (rank == RANK_MASTER)
			sub_timer = MPI_Wtime();
		//compute special actions (border, obstacle...)
		special_cells(&mesh, &mesh_type, &mesh_comm);
		if (rank == RANK_MASTER)
			mpi_put("Special cells computed in %.2lf us", toMicroSeconds(MPI_Wtime() - sub_timer));

		if (rank == RANK_MASTER)
			sub_timer = MPI_Wtime();
		//compute collision term
		collision(&temp, &mesh);
		if (rank == RANK_MASTER)
			mpi_put("Collision computed in     \t%.2lf us", toMicroSeconds(MPI_Wtime() - sub_timer));

		if (rank == RANK_MASTER)
			sub_timer = MPI_Wtime();
		//propagate values from node to neighbors
		lbm_comm_ghost_exchange(&mesh_comm, &temp, rank);
		if (rank == RANK_MASTER)
			mpi_put("Ghost exchange computed in\t%.2lf us", toMicroSeconds(MPI_Wtime() - sub_timer));

		if (rank == RANK_MASTER)
			sub_timer = MPI_Wtime();
		propagation(&mesh, &temp);
		if (rank == RANK_MASTER)
			mpi_put("Propagation computed in   \t%.2lf us", toMicroSeconds(MPI_Wtime() - sub_timer));

		if (rank == RANK_MASTER) {
			double elapsed = MPI_Wtime() - iteration_timer;
			mpi_put("-- Iteration %.5d took    \t%.2lf us", i, toMicroSeconds(elapsed));
			total_time += MPI_Wtime() - iteration_timer;
		}

		//save step
		if (i % WRITE_STEP_INTERVAL == 0 && lbm_gbl_config.output_filename != NULL) {
			if (rank == RANK_MASTER) mpi_put("Saving step");
			save_frame_all_domain(fp, &mesh, &temp_render, &mesh_comm);
#ifdef RELEASE_MODE
			float percent = (float) i / (float) ITERATIONS * 100.0F;
			printf("\033[0;35m\r %f%% -- Iteration %.05d/%.05d\033[0m", percent, i,
			       ITERATIONS);
			fflush(stdout);
#endif
		}
	}
	if (rank == RANK_MASTER)
		printf("\n");

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == RANK_MASTER && fp != NULL) {
		mpi_put("Close output file");
		close_file(fp);
	}
	if (rank == RANK_MASTER)
		printf("\n-- Total time : %.2lf s ~ %.2lf µs / iter\n\n", total_time, toMicroSeconds(total_time) / ITERATIONS);

	//free memory
	if (rank == RANK_MASTER)
		mpi_put("Freeing memory...");
	lbm_comm_release(&mesh_comm);
	Mesh_release(&mesh);
	Mesh_release(&temp);
	if (rank == RANK_MASTER)Mesh_release(&temp_render);
	lbm_mesh_type_t_release(&mesh_type);

	//close MPI
	MPI_Finalize();
	return EXIT_SUCCESS;
}
