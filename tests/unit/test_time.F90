!>------------------------------------------------------------
!! Test suite for HICAR's time machinery (time_object /
!! time_delta_object), salvaged from the ICAR-era test_calendar.F90
!! and test_time_obj.F90.
!!
!! 1. Calendar round-trip: for every supported calendar, MJD ->
!!    (y,m,d,h,m,s) -> MJD must be the identity to < 1 s over a
!!    ~2100-year span. A non-integer step samples all months/days and
!!    sub-day times; leap-rule bugs appear as a persistent offset, so
!!    the coarsened step (vs the original 0.1 d) loses no sensitivity.
!! 2. Accumulation precision: repeatedly adding a small time_delta_t
!!    must not drift against exact int64/real128 references.
!!------------------------------------------------------------
module test_time

    use iso_fortran_env,   only : real64, real128, int64
    use time_object,       only : time_type
    use time_delta_object, only : time_delta_t
    use testdrive,         only : new_unittest, unittest_type, error_type, check

    implicit none
    private

    public :: collect_time_suite

    real(real64), parameter :: MAX_RT_ERROR = 1e-5_real64  ! days (~0.86 s)

contains

    subroutine collect_time_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("calendar_gregorian", test_cal_gregorian), &
            new_unittest("calendar_standard",  test_cal_standard), &
            new_unittest("calendar_365day",    test_cal_365day), &
            new_unittest("calendar_noleap",    test_cal_noleap), &
            new_unittest("calendar_360day",    test_cal_360day), &
            new_unittest("set_from_string",    test_set_from_string), &
            new_unittest("delta_accumulation", test_delta_accumulation) &
          ]

    end subroutine collect_time_suite


    !> MJD -> date -> MJD round trip over ~2100 years for one calendar.
    subroutine calendar_roundtrip(error, calendar_name)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), intent(in) :: calendar_name

        type(time_type) :: time
        real(real64) :: mjd_in, mjd_out, max_err
        integer :: year, month, day, hour, minute, second
        character(len=128) :: msg

        call time%init(calendar_name)

        max_err = 0
        mjd_in = 365.0_real64
        do while (mjd_in <= 365.0_real64 * 2100.0_real64)
            call time%set(mjd_in)
            call time%date(year, month, day, hour, minute, second)
            mjd_out = time%date_to_mjd(year, month, day, hour, minute, second)

            if (day < 1 .or. day > 31 .or. month < 1 .or. month > 12) then
                write(msg, '(A,A,6I6)') trim(calendar_name), &
                    ": unphysical date ", year, month, day, hour, minute, second
                call check(error, .false., trim(msg))
                return
            endif

            max_err = max(max_err, abs(mjd_out - mjd_in))
            if (max_err > MAX_RT_ERROR) then
                write(msg, '(A,A,6I6,A,ES10.3)') trim(calendar_name), &
                    ": round-trip error at ", year, month, day, hour, minute, second, &
                    " err(days)=", max_err
                call check(error, .false., trim(msg))
                return
            endif

            ! non-integer step: samples all months/days and sub-day times
            mjd_in = mjd_in + 0.7321_real64
        end do

        call check(error, max_err <= MAX_RT_ERROR, &
                   trim(calendar_name)//": calendar round trip exceeded tolerance")
    end subroutine calendar_roundtrip

    subroutine test_cal_gregorian(error)
        type(error_type), allocatable, intent(out) :: error
        call calendar_roundtrip(error, "gregorian")
    end subroutine test_cal_gregorian

    subroutine test_cal_standard(error)
        type(error_type), allocatable, intent(out) :: error
        call calendar_roundtrip(error, "standard")
    end subroutine test_cal_standard

    subroutine test_cal_365day(error)
        type(error_type), allocatable, intent(out) :: error
        call calendar_roundtrip(error, "365-day")
    end subroutine test_cal_365day

    subroutine test_cal_noleap(error)
        type(error_type), allocatable, intent(out) :: error
        call calendar_roundtrip(error, "noleap")
    end subroutine test_cal_noleap

    subroutine test_cal_360day(error)
        type(error_type), allocatable, intent(out) :: error
        call calendar_roundtrip(error, "360-day")
    end subroutine test_cal_360day


    !> Date strings parse to the expected calendar date.
    subroutine test_set_from_string(error)
        type(error_type), allocatable, intent(out) :: error
        type(time_type) :: time
        integer :: year, month, day, hour, minute, second

        call time%init("gregorian")
        call time%set("2017-02-14 10:30:00")
        call time%date(year, month, day, hour, minute, second)

        call check(error, year == 2017 .and. month == 2 .and. day == 14 &
                          .and. hour == 10 .and. minute == 30 .and. second == 0, &
                   "set_from_string('2017-02-14 10:30:00') round trip failed")
    end subroutine test_set_from_string


    !> Advancing time 1e6 times by a 3600 s delta (~114 years of hourly
    !! steps) must track exact int64/real128 second counters without
    !! drifting. Uses the PRODUCTION increment pattern
    !! `set(mjd() + dt%days())` (flow_obj_t%increment_sim_time) — the
    !! original test's `t + dt` operator no longer exists on time_type.
    subroutine test_delta_accumulation(error)
        type(error_type), allocatable, intent(out) :: error

        type(time_type) :: t_obj
        type(time_delta_t) :: dt
        real(real128) :: t128
        integer(int64) :: i, i64
        integer(int64), parameter :: n = 1000000_int64

        call t_obj%init("gregorian")
        call t_obj%set(1980, 1, 1, 0, 0, 0)
        call dt%set(seconds=3600.0)
        t128 = t_obj%seconds()
        i64 = int(t_obj%seconds(), kind=int64)

        do i = 1, n
            call t_obj%set(t_obj%mjd() + dt%days())
            t128 = t128 + 3600.0_real128
            i64 = i64 + 3600_int64

            if (abs(t128 - t_obj%seconds()) > 1.0_real128) then
                call check(error, .false., &
                    "time accumulation drifted > 1 s from the real128 reference")
                return
            endif
            if (abs(i64 - int(t_obj%seconds(), kind=int64)) > 0_int64) then
                call check(error, .false., &
                    "time accumulation drifted from the int64 reference")
                return
            endif
        end do

        call check(error, .true., "unreachable")
    end subroutine test_delta_accumulation

end module test_time
