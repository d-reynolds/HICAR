#!/usr/bin/env python
"""Generate the Witch of Agnesi idealized test case for the RANS wind solver.

Creates, inside tests/Test_Cases (the HICAR-Model/Test-Data checkout):
  domains/agnesi_250m.nc        -- hi-res domain: Agnesi ridge on a flat plain,
                                   same grid/format as flat_plane_250m.nc
  forcing_agnesi/agnesi_*.nc    -- idealized COSMO-format forcing: uniform
                                   westerly U0, dry, constant N stratification,
                                   flat forcing terrain (HSURF = 0)
  input/file_list_Agnesi.txt    -- absolute-path forcing file list

The case follows docs/rans_solver_math.md section 10 (validation step 3):
stationary mountain-wave response over h(x) = HM / (1 + ((x-xc)/A)^2),
uniform in y. Linear-theory scales for the default parameters:
  Nh/U  = 0.5   (weakly nonlinear)
  Na/U  = 2.5   (mostly hydrostatic waves)
  lambda_z = 2*pi*U/N ~ 6.3 km vertical wavelength
  w ~ U*HM/A = 2 m/s expected wave amplitude over the slopes

Run with a python that has netCDF4 (e.g. the wrf_py conda env).
"""

import os
import numpy as np
import netCDF4 as nc

# ------------------------- parameters -------------------------
HM = 100.0      # hill height [m] -- Nh/U = 0.1 keeps the wave linear so the
                # model can be compared quantitatively against linear theory.
                # hm=500 (Nh/U = 0.5, weakly nonlinear) also runs stably for
                # the full 4 h since the lateral absorber was added
                # (kRANS_LATERAL_TAU; before it, wave energy shoaled in the
                # relaxation taper and a ring-edge instability grew after
                # ~1-2 h): wave saturates at |w| ~ 2.6 m/s aloft, linear
                # correlation 0.91-0.97, with a slowly-creeping ~0.8 m/s
                # upstream residual at the inflow ring (physical columnar
                # upstream influence the ring must continuously absorb).
A = 2500.0      # hill half-width [m]
# y-envelope: the ridge is flat along y in the center and tapers smoothly to
# zero well before the N/S boundaries. Terrain must NOT run through the
# lateral Davies relaxation ring: the ring nudges winds toward the (no-wave)
# balanced forcing state, and a wave-bearing interior meeting a no-wave ring
# over terrain creates a shear line that destabilizes the prognostic solver
# (observed: 30 m/s ring-edge jets and deep +-12 m/s w columns within 1 h).
Y_FLAT = 8000.0    # half-length of the uniform ridge section [m]
Y_TAPER = 6000.0   # width of the cosine^2 taper to zero [m]
U0 = 10.0       # forcing wind speed, westerly [m/s]
N_BV = 0.01     # Brunt-Vaisala frequency [1/s]
TH0 = 280.0     # surface potential temperature [K]
P0 = 100000.0   # surface pressure [Pa]
QV0 = 1.0e-6    # specific humidity (essentially dry) [kg/kg]
DX = 250.0      # hi-res grid spacing [m]

NZF = 60        # forcing vertical levels
DZF0 = 50.0     # forcing first layer thickness [m]
RZF = 1.05      # forcing layer stretching ratio

HOURS = range(0, 5)   # forcing valid times: 2017-02-14 00..04 UTC

GRAV = 9.81
RD = 287.04
CP = 1004.5

TEST_CASES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "tests", "Test_Cases")


def theta_profile(z):
    """Constant-N potential temperature profile."""
    return TH0 * np.exp(N_BV**2 * z / GRAV)


def hydrostatic_p_t(z_levels):
    """Integrate hydrostatic pressure and temperature on a column of
    level-center heights z_levels (1d, ascending)."""
    th = theta_profile(z_levels)
    p = np.zeros_like(z_levels)
    t = np.zeros_like(z_levels)
    # surface to first level, then level to level, using layer-mean T
    z_prev, p_prev = 0.0, P0
    for k, z in enumerate(z_levels):
        # iterate twice on the layer-mean temperature
        pk = p_prev
        for _ in range(3):
            tk = th[k] * (pk / P0) ** (RD / CP)
            if k == 0:
                t_mean = 0.5 * (TH0 * 1.0 + tk)   # theta ~ T at the surface
            else:
                t_mean = 0.5 * (t[k - 1] + tk)
            pk = p_prev * np.exp(-GRAV * (z - z_prev) / (RD * t_mean))
        p[k] = pk
        t[k] = th[k] * (pk / P0) ** (RD / CP)
        z_prev, p_prev = z, pk
    return p, t


def make_domain():
    src = nc.Dataset(os.path.join(TEST_CASES, "domains", "flat_plane_250m.nc"))
    dst_path = os.path.join(TEST_CASES, "domains", "agnesi_250m.nc")
    if os.path.exists(dst_path):
        os.remove(dst_path)
    dst = nc.Dataset(dst_path, "w", format="NETCDF4")

    for name, dim in src.dimensions.items():
        dst.createDimension(name, len(dim))

    ny = len(src.dimensions["y"])
    nx = len(src.dimensions["x"])

    # Agnesi ridge: centered in x, uniform along the central y section and
    # tapered to zero before the N/S relaxation rings (see Y_FLAT/Y_TAPER)
    x = (np.arange(nx) - (nx - 1) / 2.0) * DX
    y = (np.arange(ny) - (ny - 1) / 2.0) * DX
    topo_1d = HM / (1.0 + (x / A) ** 2)
    yd = np.maximum(np.abs(y) - Y_FLAT, 0.0)
    envelope = np.where(yd < Y_TAPER,
                        np.cos(0.5 * np.pi * yd / Y_TAPER) ** 2, 0.0)
    topo = envelope[:, None] * topo_1d[None, :]

    dhdx = np.gradient(topo, DX, axis=1)
    dhdy = np.gradient(topo, DX, axis=0)
    slope_rad = np.arctan(np.hypot(dhdx, dhdy))
    aspect_rad = np.arctan2(-dhdy, -dhdx)   # downslope direction

    replace = {
        "topo": topo,
        "slope": np.degrees(slope_rad),
        "slope_rad": slope_rad,
        "aspect_rad": aspect_rad,
    }

    for name, var in src.variables.items():
        out = dst.createVariable(name, var.dtype, var.dimensions)
        out.setncatts({k: var.getncattr(k) for k in var.ncattrs()
                       if k != "_FillValue"})
        out[:] = replace.get(name, var[:])

    dst.description = ("Witch of Agnesi ridge (hm=%.0fm, a=%.0fm) on a flat "
                       "plain; generated by helpers/gen_agnesi_test.py" % (HM, A))
    dst.close()
    src.close()
    print("wrote", dst_path)


def make_forcing():
    # reuse the COSMO lat/lon grid so the hi-res domain is well inside it
    laf = nc.Dataset(os.path.join(TEST_CASES, "forcing", "laf2017021400.nc"))
    lat = laf["lat_1"][:].astype(np.float32)
    lon = laf["lon_1"][:].astype(np.float32)
    laf.close()
    ny, nx = lat.shape

    dzf = DZF0 * RZF ** np.arange(NZF)
    z_if = np.concatenate(([0.0], np.cumsum(dzf)))
    z_f = 0.5 * (z_if[:-1] + z_if[1:])          # level centers, ascending

    p_col, t_col = hydrostatic_p_t(z_f)

    out_dir = os.path.join(TEST_CASES, "forcing_agnesi")
    os.makedirs(out_dir, exist_ok=True)
    paths = []

    for hour in HOURS:
        path = os.path.join(out_dir, "agnesi_201702140%d.nc" % hour)
        if os.path.exists(path):
            os.remove(path)
        f = nc.Dataset(path, "w", format="NETCDF4")
        f.createDimension("time", None)
        f.createDimension("z_2", NZF)
        f.createDimension("y_1", ny)
        f.createDimension("x_1", nx)

        tv = f.createVariable("time", "f8", ("time",))
        tv.standard_name = "time"
        tv.units = "hours since 2017-02-14 %02d:00:00" % hour
        tv.calendar = "gregorian"
        tv[:] = [0.0]

        latv = f.createVariable("lat_1", "f4", ("y_1", "x_1"))
        latv.standard_name = "latitude"
        latv.units = "degrees_north"
        latv[:] = lat
        lonv = f.createVariable("lon_1", "f4", ("y_1", "x_1"))
        lonv.standard_name = "longitude"
        lonv.units = "degrees_east"
        lonv[:] = lon

        hsurf = f.createVariable("HSURF", "f4", ("y_1", "x_1"))
        hsurf.units = "m"
        hsurf.long_name = "Geometric height of the earth surface above msl"
        hsurf[:] = 0.0

        hfl = f.createVariable("HFL", "f4", ("z_2", "y_1", "x_1"))
        hfl.units = "m"
        hfl.long_name = "Geometric height of model full level above msl"
        hfl[:] = z_f[:, None, None] * np.ones((1, ny, nx), np.float32)

        def var3d(name, units, data_col, long_name):
            v = f.createVariable(name, "f4", ("time", "z_2", "y_1", "x_1"))
            v.units = units
            v.long_name = long_name
            v[:] = (np.asarray(data_col, np.float32)[None, :, None, None]
                    * np.ones((1, 1, ny, nx), np.float32))
            return v

        var3d("P", "Pa", p_col, "Air pressure")
        var3d("T", "K", t_col, "Air temperature")
        var3d("QV", "1", np.full(NZF, QV0), "Specific humidity")
        var3d("U", "m s-1", np.full(NZF, U0), "Grid eastward wind")
        var3d("V", "m s-1", np.zeros(NZF), "Grid northward wind")
        var3d("W", "m s-1", np.zeros(NZF), "Vertical velocity (geometric)")

        f.description = ("Idealized forcing for the Agnesi RANS test: uniform "
                         "U=%.1f m/s, N=%.3f 1/s, dry, flat terrain" % (U0, N_BV))
        f.close()
        paths.append(os.path.abspath(path))
        print("wrote", path)

    list_path = os.path.join(TEST_CASES, "input", "file_list_Agnesi.txt")
    with open(list_path, "w") as fl:
        for p in paths:
            fl.write('"%s"\n' % p)
    print("wrote", list_path)


if __name__ == "__main__":
    make_domain()
    make_forcing()
