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
 * comm      [out] - pointer to ncclComm_t handle
 * nranks    [in]  - total number of ranks in the communicator
 * rank      [in]  - this rank's index (0-based)
 * f_comm    [in]  - Fortran MPI communicator handle (integer)
 * device_id [in]  - CUDA device index to bind this rank to
 * Returns 0 on success, non-zero on error.
 */
int nccl_comm_init(ncclComm_t *comm, int nranks, int rank, MPI_Fint f_comm, int device_id) {
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    ncclUniqueId id;

    /* Ensure the correct CUDA device is set for this rank */
    CUDA_CHECK(cudaSetDevice(device_id));

    /* Force creation of the CUDA primary context.
     * OpenACC's acc_init() uses the CUDA driver API, which is invisible to
     * the runtime API that NCCL uses internally. cudaFree(0) is the standard
     * NVIDIA idiom to establish the runtime primary context. */
    CUDA_CHECK(cudaFree(0));

    /* Detect if all NCCL ranks share the same physical node.
     * If so, disable the external network plugin (libnccl-net.so) which can
     * segfault due to ABI incompatibilities on Cray/HPE systems with the
     * NVHPC-bundled NCCL. P2P/SHM transports suffice for intra-node comm.
     * overwrite=0 allows users to override via NCCL_NET_DISABLE=0. */
    MPI_Comm shm_comm;
    MPI_Comm_split_type(c_comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &shm_comm);
    int shm_size;
    MPI_Comm_size(shm_comm, &shm_size);
    MPI_Comm_free(&shm_comm);

    if (shm_size == nranks) {
        /* Empty string tells NCCL to skip loading any external net plugin.
         * Also disable InfiniBand as a belt-and-suspenders measure.
         * overwrite=0: user can override by setting these vars themselves. */
        setenv("NCCL_NET_PLUGIN", "", 0);
        setenv("NCCL_IB_DISABLE", "1", 0);
        if (rank == 0)
            fprintf(stderr, "[NCCL init] All %d ranks on same node, disabled network plugin\n", nranks);
    } else {
        if (rank == 0)
            fprintf(stderr, "[NCCL init] Multi-node detected (%d/%d local), network plugin enabled\n",
                    shm_size, nranks);
    }

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
 * Create a CUDA event with timing disabled (cheap; we only use it for
 * cross-stream dependency edges, never for elapsed-time queries).
 */
int cuda_event_create(cudaEvent_t *event) {
    CUDA_CHECK(cudaEventCreateWithFlags(event, cudaEventDisableTiming));
    return 0;
}

/*
 * Destroy a CUDA event.
 */
void cuda_event_destroy(cudaEvent_t event) {
    cudaEventDestroy(event);
}

/*
 * Record a CUDA event onto a stream. When the stream's preceding work
 * completes, the event fires. Non-blocking on the host.
 */
int cuda_event_record(cudaEvent_t event, cudaStream_t stream) {
    CUDA_CHECK(cudaEventRecord(event, stream));
    return 0;
}

/*
 * Make a CUDA stream wait for an event before launching subsequent work.
 * GPU-side dependency edge — non-blocking on the host.
 */
int cuda_stream_wait_event(cudaStream_t stream, cudaEvent_t event) {
    CUDA_CHECK(cudaStreamWaitEvent(stream, event, 0));
    return 0;
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

/*
 * NCCL point-to-point send (double).
 * Used by the native iterative wind solver for Krylov-vector halo
 * exchange — vectors are double precision.
 */
int nccl_send_double(const double *buf, int count, int peer,
                     ncclComm_t comm, cudaStream_t stream) {
    NCCL_CHECK(ncclSend(buf, (size_t)count, ncclDouble, peer, comm, stream));
    return 0;
}

/*
 * NCCL point-to-point recv (double).
 */
int nccl_recv_double(double *buf, int count, int peer,
                     ncclComm_t comm, cudaStream_t stream) {
    NCCL_CHECK(ncclRecv(buf, (size_t)count, ncclDouble, peer, comm, stream));
    return 0;
}

/*
 * NCCL all-reduce with sum op (double precision, in-place or out-of-place).
 * sendbuf [in]  - device pointer to local data
 * recvbuf [out] - device pointer for global result (may equal sendbuf for in-place)
 * count   [in]  - number of double elements
 * comm    [in]  - NCCL communicator
 * stream  [in]  - CUDA stream
 * Returns 0 on success.
 *
 * Used by the native iterative wind solver for the BiCGStab dot-product
 * reductions (sigma = <r̂,v> and the 5-value omega/rho/||r|| reduction).
 * Replaces MPI_Iallreduce + MPI_Wait so the reduction stays on-device,
 * stream-ordered with the OpenACC kernels that produced the partial sums.
 */
int nccl_allreduce_double_sum(const double *sendbuf, double *recvbuf,
                              int count, ncclComm_t comm, cudaStream_t stream) {
    NCCL_CHECK(ncclAllReduce(sendbuf, recvbuf, (size_t)count, ncclDouble,
                             ncclSum, comm, stream));
    return 0;
}
