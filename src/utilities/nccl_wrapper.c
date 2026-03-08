/*
 * C wrapper for NCCL functions callable from Fortran via ISO_C_BINDING.
 * Follows the same pattern as mpi_comm_wrapper.c for AmgX.
 *
 * Compiled only when USE_NCCL is defined (set by CMake when NCCL is found).
 */

#include <nccl.h>
#include <cuda_runtime.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* Error checking macros */
#define NCCL_CHECK(cmd) do {                            \
    ncclResult_t r = cmd;                               \
    if (r != ncclSuccess) {                             \
        fprintf(stderr, "NCCL error %s:%d '%s'\n",     \
                __FILE__, __LINE__,                     \
                ncclGetErrorString(r));                 \
        return (int)r;                                  \
    }                                                   \
} while(0)

#define CUDA_CHECK(cmd) do {                            \
    cudaError_t e = cmd;                                \
    if (e != cudaSuccess) {                             \
        fprintf(stderr, "CUDA error %s:%d '%s'\n",     \
                __FILE__, __LINE__,                     \
                cudaGetErrorString(e));                 \
        return (int)e;                                  \
    }                                                   \
} while(0)

/*
 * Initialize NCCL communicator using MPI for unique ID broadcast.
 * comm  [out] - pointer to ncclComm_t handle
 * nranks [in] - total number of ranks in the communicator
 * rank   [in] - this rank's index (0-based)
 * f_comm [in] - Fortran MPI communicator handle (integer)
 * Returns 0 on success, non-zero on error.
 */
int nccl_comm_init(ncclComm_t *comm, int nranks, int rank, MPI_Fint f_comm) {
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    ncclUniqueId id;

    if (rank == 0) {
        ncclResult_t r = ncclGetUniqueId(&id);
        if (r != ncclSuccess) {
            fprintf(stderr, "NCCL error: ncclGetUniqueId failed: %s\n",
                    ncclGetErrorString(r));
            return (int)r;
        }
    }

    MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, c_comm);

    NCCL_CHECK(ncclCommInitRank(comm, nranks, id, rank));
    return 0;
}

/*
 * Destroy NCCL communicator.
 */
void nccl_comm_destroy(ncclComm_t comm) {
    ncclCommDestroy(comm);
}

/*
 * Create a CUDA stream for NCCL operations.
 * stream [out] - pointer to cudaStream_t handle
 * Returns 0 on success.
 */
int nccl_stream_create(cudaStream_t *stream) {
    CUDA_CHECK(cudaStreamCreateWithFlags(stream, cudaStreamNonBlocking));
    return 0;
}

/*
 * Synchronize a CUDA stream (blocks until all NCCL ops on this stream complete).
 * Returns 0 on success.
 */
int nccl_stream_synchronize(cudaStream_t stream) {
    CUDA_CHECK(cudaStreamSynchronize(stream));
    return 0;
}

/*
 * Destroy a CUDA stream.
 */
void nccl_stream_destroy(cudaStream_t stream) {
    cudaStreamDestroy(stream);
}

/*
 * Begin a NCCL group operation. All ncclSend/ncclRecv between
 * nccl_group_start and nccl_group_end are fused into a single operation.
 */
void nccl_group_start(void) {
    ncclGroupStart();
}

/*
 * End a NCCL group operation.
 */
void nccl_group_end(void) {
    ncclGroupEnd();
}

/*
 * NCCL point-to-point send (float/real).
 * buf    [in]  - device pointer to send buffer
 * count  [in]  - number of float elements to send
 * peer   [in]  - destination rank
 * comm   [in]  - NCCL communicator
 * stream [in]  - CUDA stream
 * Returns 0 on success.
 */
int nccl_send_float(const float *buf, int count, int peer,
                    ncclComm_t comm, cudaStream_t stream) {
    NCCL_CHECK(ncclSend(buf, (size_t)count, ncclFloat, peer, comm, stream));
    return 0;
}

/*
 * NCCL point-to-point recv (float/real).
 * buf    [out] - device pointer to receive buffer
 * count  [in]  - number of float elements to receive
 * peer   [in]  - source rank
 * comm   [in]  - NCCL communicator
 * stream [in]  - CUDA stream
 * Returns 0 on success.
 */
int nccl_recv_float(float *buf, int count, int peer,
                    ncclComm_t comm, cudaStream_t stream) {
    NCCL_CHECK(ncclRecv(buf, (size_t)count, ncclFloat, peer, comm, stream));
    return 0;
}
