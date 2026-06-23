#!/usr/bin/env python
"""Quick-look diagnostics for the Witch of Agnesi RANS validation case.

Usage: check_agnesi_rans.py <hicar_output_file.nc> [--plot out.png]

Reports, for each output time:
  * NaN / Inf counts and global min/max of u, v, w
  * centerline (mid-y) wave diagnostics: max |w|, height/x of the extrema,
    u min/max
  * 2-dx noise index of w at ~1.5 km AGL (ratio of grid-scale variance)
  * comparison against the linear-theory scale w ~ U*HM/A

With --plot, writes a centerline x-z contour figure of w and theta at the
last output time.
"""

import sys
import numpy as np
import netCDF4 as nc

HM, A, U0, N_BV = 100.0, 2500.0, 10.0, 0.01   # keep in sync with gen_agnesi_test.py
W_SCALE = U0 * HM / A   # linear-theory vertical-velocity scale [m/s]


def two_dx_noise(field):
    """Fraction of variance at the 2-dx scale along x (mid altitude row)."""
    f = field - field.mean()
    if f.size < 8 or np.allclose(f, 0):
        return 0.0
    spec = np.abs(np.fft.rfft(f)) ** 2
    return float(spec[-max(2, len(spec) // 8):].sum() / max(spec.sum(), 1e-30))


def linear_w_fourier(terr_row, z_cl, dx=250.0):
    """Exact 2-D linear stationary mountain-wave w (constant U, N) by
    Fourier synthesis of the actual centerline terrain row: for each
    along-x wavenumber k, w_hat = i*k*U*h_hat*exp(i*m*z) with
    m = sign(k)*sqrt(N^2/U^2 - k^2) for propagating modes (upward
    radiation) and i*sqrt(k^2 - N^2/U^2) for evanescent ones.
    Includes the anelastic sqrt(rho0/rho) ~ exp(z/2H) amplitude factor."""
    nx = terr_row.size
    k = 2 * np.pi * np.fft.rfftfreq(nx, d=dx)
    h_hat = np.fft.rfft(terr_row)
    l2 = (N_BV / U0) ** 2
    m = np.where(k ** 2 < l2,
                 np.sqrt(np.maximum(l2 - k ** 2, 0.0)),          # propagating
                 1j * np.sqrt(np.maximum(k ** 2 - l2, 0.0)))     # evanescent
    w = np.empty_like(z_cl)
    for kk in range(z_cl.shape[0]):
        z = z_cl[kk].mean()
        w_hat = 1j * k * U0 * h_hat * np.exp(1j * m * z)
        w[kk] = np.fft.irfft(w_hat, n=nx).real * np.exp(z / (2 * 9500.0))
    return w


def compare_linear(w_cl, z_cl, terr_cl, xm):
    """Correlation and amplitude ratio of the centerline model w against
    linear theory, in the window |x|<=10 km, 0.5 km <= z-h <= 5 km."""
    zagl = z_cl - terr_cl[None, :]
    wl = linear_w_fourier(terr_cl, z_cl)
    sel = (np.abs(np.broadcast_to(xm, z_cl.shape)) <= 10000.0) & \
          (zagl >= 500.0) & (zagl <= 5000.0)
    a, b = w_cl[sel], wl[sel]
    corr = float(np.corrcoef(a, b)[0, 1])
    ratio = float(np.sqrt((a ** 2).mean() / max((b ** 2).mean(), 1e-30)))
    return corr, ratio


def main():
    path = sys.argv[1]
    plot = "--plot" in sys.argv
    d = nc.Dataset(path)
    w = d["w"][:]            # (time, z, y, x) cell-centered real w
    u = d["u"][:]            # staggered in x
    v = d["v"][:]
    z = d["z"][:]            # (time?, z, y, x) heights
    if z.ndim == 4:
        z = z[0]
    nt, nz, ny, nx = w.shape
    jc = ny // 2
    terr_cl = d["terrain"][jc, :] if "terrain" in d.variables else z[0, jc, :]
    xm = (np.arange(nx) - (nx - 1) / 2.0) * 250.0

    print(f"file: {path}")
    print(f"dims: nt={nt} nz={nz} ny={ny} nx={nx}; linear w scale ~ {W_SCALE:.2f} m/s")
    times = d["time"][:] if "time" in d.variables else range(nt)

    for n in range(nt):
        bad = 0
        for name, fld in (("u", u[n]), ("v", v[n]), ("w", w[n])):
            nbad = np.count_nonzero(~np.isfinite(fld))
            bad += nbad
            if nbad:
                print(f"  !! {name} has {nbad} non-finite values")
        wc = w[n, :, jc, :]
        k15 = np.argmin(np.abs(z[:, jc, nx // 2] - z[0, jc, nx // 2] - 1500.0))
        noise = two_dx_noise(wc[k15])
        kmax, imax = np.unravel_index(np.argmax(np.abs(wc)), wc.shape)
        corr, ratio = compare_linear(wc, z[:, jc, :], np.asarray(terr_cl), xm)
        print(f"t={n}: |w|max={np.abs(wc).max():6.3f} m/s at k={kmax},i={imax} "
              f"(z~{z[kmax, jc, imax]:.0f}m)  u:[{u[n].min():6.2f},{u[n].max():6.2f}] "
              f"v:[{v[n].min():6.2f},{v[n].max():6.2f}] "
              f"2dx-noise(w@1.5km)={noise:.2f}  lin-corr={corr:5.2f} "
              f"amp-ratio={ratio:5.2f}  finite-ok={bad == 0}")

    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        out = sys.argv[sys.argv.index("--plot") + 1]
        wc = w[-1, :, jc, :]
        zc = z[:, jc, :]
        xkm = (np.arange(nx) - (nx - 1) / 2) * 0.25
        x2 = np.broadcast_to(xkm, zc.shape)
        fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
        lim = max(np.abs(wc).max(), 0.1)
        cf = axes[0].pcolormesh(x2, zc / 1000, wc, cmap="RdBu_r",
                                vmin=-lim, vmax=lim)
        fig.colorbar(cf, ax=axes[0], label="w [m/s]")
        axes[0].set_title(f"w centerline, last output (|w|max={np.abs(wc).max():.2f})")
        if "potential_temperature" in d.variables:
            th = d["potential_temperature"][-1, :, jc, :]
            cs = axes[1].contour(x2, zc / 1000, th, levels=20, colors="k",
                                 linewidths=0.5)
            axes[1].set_title("potential temperature (isentropes)")
        uc = d["u"][-1, :, jc, :]
        cf = axes[1].pcolormesh(x2[:, :], zc / 1000, 0.5 * (uc[:, 1:] + uc[:, :-1]),
                                cmap="viridis", alpha=0.6)
        fig.colorbar(cf, ax=axes[1], label="u [m/s]")
        for ax in axes:
            ax.plot(xkm, zc[0] / 1000, "k", lw=1.5)
            ax.set_ylabel("z [km]")
            ax.set_ylim(0, 8)
        axes[1].set_xlabel("x [km]")
        fig.tight_layout()
        fig.savefig(out, dpi=120)
        print("wrote", out)


if __name__ == "__main__":
    main()
