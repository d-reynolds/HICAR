module mpi_utils_module
    use mpi_f08
    implicit none

contains

    function get_mpi_global_rank()
        implicit none
        integer :: get_mpi_global_rank

        ! Get the rank of this MPI process on the global communicator
        call MPI_Comm_rank(MPI_COMM_WORLD, get_mpi_global_rank)

    end function get_mpi_global_rank

end module mpi_utils_module