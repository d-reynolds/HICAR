/*
 * C wrapper for MPI_Comm_f2c conversion and AmgX resource creation
 * This allows Fortran code to convert MPI_Comm to a C pointer for AmgX
 */

#include <mpi.h>
#include <amgx_c.h>
#include <stdlib.h>
#include <stdio.h>

/*
 * Wrapper for AMGX_resources_create that accepts Fortran MPI communicator
 * Converts Fortran MPI_Comm to C MPI_Comm and calls AMGX_resources_create
 * 
 * Note: AMGX stores the pointer to the communicator, so we use a static variable
 * to ensure it persists for the lifetime of the program.
 */
int AMGX_resources_create_wrapper(AMGX_resources_handle *resources, 
                                   AMGX_config_handle config,
                                   int *fcomm,
                                   int num_devices,
                                   int *device_ids) {
    static MPI_Comm ccomm;  /* Static to persist for program lifetime */
    MPI_Fint fint_comm;
    int rc; //, comm_size, comm_rank;
    
    /* Dereference the pointer to get the Fortran integer handle */
    fint_comm = (MPI_Fint)(*fcomm);
    
    /* Convert Fortran MPI communicator to C MPI communicator */
    ccomm = MPI_Comm_f2c(fint_comm);
    
    /* Validate the communicator */
    if (ccomm == MPI_COMM_NULL) {
        fprintf(stderr, "ERROR: Converted MPI communicator is MPI_COMM_NULL\n");
        return -1;
    }
        
    /* Call the AMGX library function with the converted MPI communicator */
    /* Pass address of the static MPI_Comm variable */
    rc = AMGX_resources_create(resources, config, &ccomm, num_devices, device_ids);
        
    return rc;
}

/*
 * Device-pointer wrappers for AMGX functions.
 * When AMGX is configured in dDDI mode and passed device pointers
 * (via OpenACC host_data use_device), cudaMemcpy with cudaMemcpyDefault
 * performs D2D transfers automatically.
 */
int AMGX_vector_upload_device(AMGX_vector_handle vec, int n, int block_dim, const void* data) {
    return AMGX_vector_upload(vec, n, block_dim, data);
}

int AMGX_vector_download_device(AMGX_vector_handle vec, void* data) {
    return AMGX_vector_download(vec, data);
}

int AMGX_matrix_replace_coefficients_device(AMGX_matrix_handle mtx, int n, int nnz,
                                             const void* values, const void* diag) {
    return AMGX_matrix_replace_coefficients(mtx, n, nnz, values, diag);
}
