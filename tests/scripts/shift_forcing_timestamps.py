#!/usr/bin/env python3
"""Shift the timestamps of HICAR forcing files by a fixed number of hours.

The Test_Cases forcing shipped with the Test_Data repo (laf2017021400-03.nc)
covers 00-03 UTC on 2017-02-14 -- the middle of the night, so there is no
shortwave radiation to exercise. This utility copies each forcing file into a
destination directory while advancing BOTH

  * the filename timestamp   (lafYYYYMMDDHH.nc), and
  * the internal `time` coordinate (and its `units` reference date)

by <shift_hours> (default 10). The atmospheric state is byte-for-byte identical;
only the clock moves. With the default +10 h the four files land at 10-13 UTC
(~solar noon at the Gaudergrat domain), so the run exercises the shortwave path
in the snow models.

The originals are never modified -- shifted copies are written to <dst_dir>.

Usage:
  shift_forcing_timestamps.py <src_dir> <dst_dir> [shift_hours]

Example:
  shift_forcing_timestamps.py tests/Test_Cases/forcing \
                              tests/Test_Cases/forcing_midday 10
"""
import sys
import os
import re
import glob
import shutil
from datetime import timedelta

import netCDF4 as nc
import cftime

NAME_RE = re.compile(r"^laf(\d{4})(\d{2})(\d{2})(\d{2})\.nc$")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    src_dir = sys.argv[1]
    dst_dir = sys.argv[2]
    shift_h = int(sys.argv[3]) if len(sys.argv) > 3 else 10

    os.makedirs(dst_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(src_dir, "laf*.nc")))
    if not files:
        sys.exit(f"No laf*.nc forcing files found in {src_dir}")

    n = 0
    for f in files:
        base = os.path.basename(f)
        m = NAME_RE.match(base)
        if not m:
            print(f"  skip {base} (filename is not lafYYYYMMDDHH.nc)")
            continue
        yyyy, mm, dd, hh = map(int, m.groups())
        # filename-derived stamp, advanced by the shift
        old_name_dt = cftime.DatetimeGregorian(yyyy, mm, dd, hh)
        new_name_dt = old_name_dt + timedelta(hours=shift_h)
        new_name = (
            f"laf{new_name_dt.year:04d}{new_name_dt.month:02d}"
            f"{new_name_dt.day:02d}{new_name_dt.hour:02d}.nc"
        )
        out = os.path.join(dst_dir, new_name)
        shutil.copyfile(f, out)

        # Advance the internal time coordinate by the same amount. Each file's
        # stamp may be encoded in the `units` reference date (value == 0) with
        # a per-file unit ('hours since'/'days since') and calendar, so decode
        # to absolute dates, add the offset, and re-encode on a clean unit.
        ds = nc.Dataset(out, "r+")
        t = ds.variables["time"]
        units = t.units
        cal = getattr(t, "calendar", "standard")
        dates = cftime.num2date(t[:], units, cal)
        shifted = [d + timedelta(hours=shift_h) for d in dates]
        ref = shifted[0]
        new_units = "hours since " + ref.strftime("%Y-%m-%d %H:%M:%S")
        t[:] = cftime.date2num(shifted, new_units, cal)
        t.units = new_units
        ds.close()

        print(f"  {base} -> {new_name}   ({old_name_dt} -> {new_name_dt})")
        n += 1

    print(f"\nShifted {n} forcing file(s) by +{shift_h} h into {dst_dir}")


if __name__ == "__main__":
    main()
