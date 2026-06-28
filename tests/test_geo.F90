!>------------------------------------------------------------
!! Test suite for the geo module (geo_reader.F90): the geometric
!! primitives behind the geographic look-up table (point_is_on_line,
!! point_in_poly — assertions salvaged from the ICAR-era
!! test_point_in_on.F90, including historical real-world failure
!! cases) and the geo_LUT / geo_interp(2d) interpolation itself on
!! synthetic grids (replacing the ICAR-era file-based test_geo.F90).
!!
!! Interpolation correctness uses two exact properties of bilinear
!! interpolation: a constant field is reproduced exactly (weights sum
!! to 1) and a field linear in lon/lat is reproduced to rounding.
!!------------------------------------------------------------
module test_geo

    use data_structures, only : interpolable_type
    use geo,             only : point_is_on_line, point_in_poly, geo_LUT, &
                                geo_interp, geo_interp2d
    use testdrive,       only : new_unittest, unittest_type, error_type, check

    implicit none
    private

    public :: collect_geo_suite

contains

    subroutine collect_geo_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("point_is_on_line",   test_point_is_on_line), &
            new_unittest("point_in_square",    test_point_in_square), &
            new_unittest("point_in_diamond",   test_point_in_diamond), &
            new_unittest("point_on_shared_edge", test_point_on_shared_edge), &
            new_unittest("point_in_concave",   test_point_in_concave), &
            new_unittest("point_in_triangles", test_point_in_triangles), &
            new_unittest("interp_constant",    test_interp_constant), &
            new_unittest("interp_linear",      test_interp_linear), &
            new_unittest("interp_3d",          test_interp_3d) &
          ]

    end subroutine collect_geo_suite


    subroutine test_point_is_on_line(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, point_is_on_line(0.5, 0.5, 0.0, 0.0, 1.0, 1.0), &
                   "(0.5,0.5) should be on the line (0,0)-(1,1)")
        if (allocated(error)) return
        call check(error, .not. point_is_on_line(0.6, 0.5, 0.0, 0.0, 1.0, 1.0), &
                   "(0.6,0.5) should be off the line (0,0)-(1,1)")
    end subroutine test_point_is_on_line


    subroutine test_point_in_square(error)
        type(error_type), allocatable, intent(out) :: error
        real :: poly(2,4)

        poly(1,:) = [0., 0., 1., 1.]
        poly(2,:) = [0., 1., 1., 0.]

        call check(error, point_in_poly(0.5, 0.5, poly), &
                   "(0.5,0.5) inside unit square")
        if (allocated(error)) return
        call check(error, .not. point_in_poly(1.01, 0.5, poly), &
                   "(1.01,0.5) outside unit square")
    end subroutine test_point_in_square


    !> Diamond polygon: vertex, edge, interior and exterior cases that
    !! historically tripped ray-casting implementations.
    subroutine test_point_in_diamond(error)
        type(error_type), allocatable, intent(out) :: error
        real :: poly(2,4)

        poly(1,:) = [0.5, 0.0, 0.5, 1.0]
        poly(2,:) = [0.0, 0.5, 1.0, 0.5]

        call check(error, point_in_poly(0.5, 0.5, poly),        "diamond: center")
        if (allocated(error)) return
        call check(error, .not. point_in_poly(-0.1, 0.5, poly), "diamond: left of left vertex")
        if (allocated(error)) return
        call check(error, .not. point_in_poly(0.1, 0.0, poly),  "diamond: outside, level with bottom vertex")
        if (allocated(error)) return
        call check(error, point_in_poly(0.51, 0.51, poly),      "diamond: just off-center")
        if (allocated(error)) return
        call check(error, point_in_poly(0.26, 0.25, poly),      "diamond: strictly inside near lower-left edge")
        if (allocated(error)) return
        call check(error, .not. point_in_poly(0.2, 0.25, poly), "diamond: just outside lower-left edge")
        if (allocated(error)) return
        call check(error, .not. point_in_poly(0.9, 0.25, poly), "diamond: outside lower-right edge")
        if (allocated(error)) return
        call check(error, .not. point_in_poly(0.2, 0.75, poly), "diamond: just outside upper-left edge")
        if (allocated(error)) return
        call check(error, point_in_poly(0.26, 0.75, poly),      "diamond: strictly inside near upper-left edge")
    end subroutine test_point_in_diamond


    !> Exact-on-edge semantics: the legacy point_in_poly included boundary
    !! points via a precision test; the current ray-cast leaves them
    !! unspecified (vertices excepted). What geo_LUT/find_surrounding
    !! actually requires is that a point on an edge SHARED by two
    !! polygons is claimed by at least one of them — otherwise the point
    !! is orphaned and the LUT search fails (the original ICAR bug).
    subroutine test_point_on_shared_edge(error)
        type(error_type), allocatable, intent(out) :: error
        real :: diamond(2,4), tri_left(2,3)

        ! (0.25,0.25) lies exactly on the diamond's lower-left edge,
        ! which it shares with the triangle (0.5,0)-(0,0.5)-(0,0)
        diamond(1,:) = [0.5, 0.0, 0.5, 1.0]
        diamond(2,:) = [0.0, 0.5, 1.0, 0.5]
        tri_left(1,:) = [0.5, 0.0, 0.0]
        tri_left(2,:) = [0.0, 0.5, 0.0]

        call check(error, point_in_poly(0.25, 0.25, diamond) .or. &
                          point_in_poly(0.25, 0.25, tri_left), &
                   "point on a shared edge must belong to at least one neighbor")
    end subroutine test_point_on_shared_edge


    !> Concave (U-shaped) polygons: the test point is inside the hull
    !! but the polygon boundary passes between it and parts of the hull.
    subroutine test_point_in_concave(error)
        type(error_type), allocatable, intent(out) :: error
        real :: poly(2,8)

        poly(1,:) = [0.0, 0.0, 1.0, 1.0, 0.75, 0.75, 0.25, 0.25]
        poly(2,:) = [0.0, 1.0, 1.0, 0.0, 0.00, 0.50, 0.50, 0.00]
        call check(error, point_in_poly(0.25, 0.5, poly), &
                   "concave U: point above the notch is inside")
        if (allocated(error)) return

        poly(1,:) = [0.0, 0.0, 1.0, 1.0,  0.75, 0.75, 0.25, 0.25]
        poly(2,:) = [0.0, 1.0, 1.0, 0.75, 0.75, 0.50, 0.50, 0.00]
        call check(error, point_in_poly(0.25, 0.5, poly), &
                   "concave S: point on the step is inside")
    end subroutine test_point_in_concave


    !> Real-world lat/lon triangles from historical find_surrounding
    !! failures in simulations (single-precision degenerate geometry).
    subroutine test_point_in_triangles(error)
        type(error_type), allocatable, intent(out) :: error
        real :: tri(2,3)
        real :: tri_mirror(2,3)

        ! Historical find_surrounding failure: the point lies ON edge B-C of
        ! its triangle to within single precision (cross product ~3e-7), so
        ! exact membership is unspecified — but it must belong to the
        ! triangle OR its mirror neighbor across that edge (no orphans).
        tri(1,:) = [35.4077759, 35.4736328, 35.4375458]
        tri(2,:) = [34.0487747, 34.0435371, 34.0188789]
        tri_mirror(1,:) = [35.4736328, 35.4375458, tri(1,2) + tri(1,3) - tri(1,1)]
        tri_mirror(2,:) = [34.0435371, 34.0188789, tri(2,2) + tri(2,3) - tri(2,1)]
        call check(error, point_in_poly(35.4495850, 34.0271034, tri) .or. &
                          point_in_poly(35.4495850, 34.0271034, tri_mirror), &
                   "historical problem point must be inside its triangle or the neighbor across the shared edge")
        if (allocated(error)) return

        tri(1,:) = [281.6548, 281.7756, 281.9876]
        tri(2,:) = [34.91632, 35.36251, 35.08936]
        call check(error, .not. point_in_poly(281.98, 35.0, tri), &
                   "point just outside the hypotenuse is outside")
        if (allocated(error)) return

        tri(1,:) = [281.6548, 282.1, 281.9876]
        tri(2,:) = [34.91632, 34.8, 35.08936]
        call check(error, point_in_poly(281.98, 35.0, tri), &
                   "point just inside the hypotenuse is inside")
        if (allocated(error)) return

        tri(1,:) = [290.0750, 289.8967, 289.7062]
        tri(2,:) = [38.79078, 38.35439, 38.64164]
        call check(error, point_in_poly(289.9, 38.72, tri), &
                   "point near the edge of a thin triangle is inside")
    end subroutine test_point_in_triangles


    !> Build a coarse regular lat/lon grid (lo) and a finer grid strictly
    !! inside it (hi), offset so the interpolation weights are fractional.
    subroutine make_grids(lo, hi)
        type(interpolable_type), intent(out) :: lo, hi
        integer, parameter :: nx_lo = 11, ny_lo = 11, nx_hi = 9, ny_hi = 9
        integer :: i, j

        allocate(lo%lat(nx_lo, ny_lo), lo%lon(nx_lo, ny_lo))
        do j = 1, ny_lo
            do i = 1, nx_lo
                lo%lon(i,j) = 0.0 + (i-1) * 1.0          ! 0 .. 10 deg
                lo%lat(i,j) = 40.0 + (j-1) * 1.0         ! 40 .. 50 deg
            enddo
        enddo

        allocate(hi%lat(nx_hi, ny_hi), hi%lon(nx_hi, ny_hi))
        do j = 1, ny_hi
            do i = 1, nx_hi
                hi%lon(i,j) = 2.3 + (i-1) * 0.6          ! 2.3 .. 7.1, off-node
                hi%lat(i,j) = 42.7 + (j-1) * 0.6         ! 42.7 .. 47.5
            enddo
        enddo
    end subroutine make_grids


    subroutine test_interp_constant(error)
        type(error_type), allocatable, intent(out) :: error
        type(interpolable_type) :: lo, hi
        real, allocatable :: fieldin(:,:), fieldout(:,:)
        real, parameter :: c0 = 273.15

        call make_grids(lo, hi)
        call geo_LUT(hi, lo)

        allocate(fieldin(size(lo%lon,1), size(lo%lon,2)), source=c0)
        allocate(fieldout(size(hi%lon,1), size(hi%lon,2)), source=0.0)

        call geo_interp2d(fieldout, fieldin, lo%geolut)

        ! interpolation weights must sum to 1: a constant field is exact
        call check(error, maxval(abs(fieldout - c0)) < 1e-4, &
                   "constant field not reproduced by geo_interp2d (weights do not sum to 1?)")
    end subroutine test_interp_constant


    subroutine test_interp_linear(error)
        type(error_type), allocatable, intent(out) :: error
        type(interpolable_type) :: lo, hi
        real, allocatable :: fieldin(:,:), fieldout(:,:), expected(:,:)
        real, parameter :: a = 5.0, b = 1.5, c = -0.75

        call make_grids(lo, hi)
        call geo_LUT(hi, lo)

        allocate(fieldin(size(lo%lon,1), size(lo%lon,2)))
        fieldin = a + b*lo%lon + c*lo%lat
        allocate(fieldout(size(hi%lon,1), size(hi%lon,2)), source=0.0)
        allocate(expected(size(hi%lon,1), size(hi%lon,2)))
        expected = a + b*hi%lon + c*hi%lat

        call geo_interp2d(fieldout, fieldin, lo%geolut)

        ! bilinear interpolation is exact for fields linear in lon/lat
        call check(error, maxval(abs(fieldout - expected)) < 1e-3, &
                   "linear field not reproduced by geo_interp2d within tolerance")
    end subroutine test_interp_linear


    subroutine test_interp_3d(error)
        type(error_type), allocatable, intent(out) :: error
        type(interpolable_type) :: lo, hi
        real, allocatable :: fieldin(:,:,:), fieldout(:,:,:)
        real, allocatable :: fieldin2d(:,:), fieldout2d(:,:)
        integer :: nx_lo, ny_lo, nx_hi, ny_hi

        call make_grids(lo, hi)
        call geo_LUT(hi, lo)

        nx_lo = size(lo%lon,1); ny_lo = size(lo%lon,2)
        nx_hi = size(hi%lon,1); ny_hi = size(hi%lon,2)

        ! two-level field, (i,k,j) ordering as used throughout HICAR
        allocate(fieldin(nx_lo, 2, ny_lo))
        fieldin(:,1,:) = 1.0 + 2.0*lo%lon + 0.5*lo%lat
        fieldin(:,2,:) = -3.0 + 0.25*lo%lon - 1.5*lo%lat
        allocate(fieldout(nx_hi, 2, ny_hi), source=0.0)

        call geo_interp(fieldout, fieldin, lo%geolut)

        ! each level of the 3D path must match the 2D path on the same data
        allocate(fieldin2d(nx_lo, ny_lo), fieldout2d(nx_hi, ny_hi))
        fieldin2d = fieldin(:,1,:)
        fieldout2d = 0.0
        call geo_interp2d(fieldout2d, fieldin2d, lo%geolut)
        call check(error, maxval(abs(fieldout(:,1,:) - fieldout2d)) < 1e-4, &
                   "geo_interp level 1 does not match geo_interp2d")
        if (allocated(error)) return

        fieldin2d = fieldin(:,2,:)
        fieldout2d = 0.0
        call geo_interp2d(fieldout2d, fieldin2d, lo%geolut)
        call check(error, maxval(abs(fieldout(:,2,:) - fieldout2d)) < 1e-4, &
                   "geo_interp level 2 does not match geo_interp2d")
    end subroutine test_interp_3d

end module test_geo
