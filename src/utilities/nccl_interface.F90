!>------------------------------------------------------------
!!  Fortran interface module for NCCL C wrapper functions.
!!  Uses ISO_C_BINDING to call nccl_wrapper.c routines.
!!
!!  Only compiled/used when USE_NCCL preprocessor flag is set.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
module nccl_interface
    use iso_c_binding
    implicit none

    interface

        integer(c_int) function nccl_comm_init(comm, nranks, rank, f_comm, device_id) bind(C, name='nccl_comm_init')
            import :: c_ptr, c_int
            type(c_ptr), intent(out) :: comm
            integer(c_int), value :: nranks, rank, f_comm, device_id
        end function

        subroutine nccl_comm_destroy(comm) bind(C, name='nccl_comm_destroy')
            import :: c_ptr
            type(c_ptr), value :: comm
        end subroutine

        integer(c_int) function nccl_stream_create(stream) bind(C, name='nccl_stream_create')
            import :: c_ptr, c_int
            type(c_ptr), intent(out) :: stream
        end function

        integer(c_int) function nccl_stream_synchronize(stream) bind(C, name='nccl_stream_synchronize')
            import :: c_ptr, c_int
            type(c_ptr), value :: stream
        end function

        subroutine nccl_stream_destroy(stream) bind(C, name='nccl_stream_destroy')
            import :: c_ptr
            type(c_ptr), value :: stream
        end subroutine

        integer(c_int) function cuda_event_create(event) bind(C, name='cuda_event_create')
            import :: c_ptr, c_int
            type(c_ptr), intent(out) :: event
        end function

        subroutine cuda_event_destroy(event) bind(C, name='cuda_event_destroy')
            import :: c_ptr
            type(c_ptr), value :: event
        end subroutine

        integer(c_int) function cuda_event_record(event, stream) bind(C, name='cuda_event_record')
            import :: c_ptr, c_int
            type(c_ptr), value :: event, stream
        end function

        integer(c_int) function cuda_stream_wait_event(stream, event) bind(C, name='cuda_stream_wait_event')
            import :: c_ptr, c_int
            type(c_ptr), value :: stream, event
        end function

        subroutine nccl_group_start() bind(C, name='nccl_group_start')
        end subroutine

        subroutine nccl_group_end() bind(C, name='nccl_group_end')
        end subroutine

        integer(c_int) function nccl_send_float(buf, count, peer, comm, stream) bind(C, name='nccl_send_float')
            import :: c_ptr, c_int
            type(c_ptr), value :: buf
            integer(c_int), value :: count, peer
            type(c_ptr), value :: comm, stream
        end function

        integer(c_int) function nccl_recv_float(buf, count, peer, comm, stream) bind(C, name='nccl_recv_float')
            import :: c_ptr, c_int
            type(c_ptr), value :: buf
            integer(c_int), value :: count, peer
            type(c_ptr), value :: comm, stream
        end function

        ! Double precision send/recv for the native iterative wind solver
        integer(c_int) function nccl_send_double(buf, count, peer, comm, stream) bind(C, name='nccl_send_double')
            import :: c_ptr, c_int
            type(c_ptr), value :: buf
            integer(c_int), value :: count, peer
            type(c_ptr), value :: comm, stream
        end function

        integer(c_int) function nccl_recv_double(buf, count, peer, comm, stream) bind(C, name='nccl_recv_double')
            import :: c_ptr, c_int
            type(c_ptr), value :: buf
            integer(c_int), value :: count, peer
            type(c_ptr), value :: comm, stream
        end function

        ! All-reduce with sum op (double precision). sendbuf and recvbuf are device pointers.
        integer(c_int) function nccl_allreduce_double_sum(sendbuf, recvbuf, count, comm, stream) &
                bind(C, name='nccl_allreduce_double_sum')
            import :: c_ptr, c_int
            type(c_ptr), value :: sendbuf, recvbuf
            integer(c_int), value :: count
            type(c_ptr), value :: comm, stream
        end function

    end interface

end module nccl_interface
