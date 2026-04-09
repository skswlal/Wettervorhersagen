"""
DWD RADOLAN RV Nowcast — Download & Parse
==========================================
Fetches the latest RV composite from DWD Open Data, parses the full
2-hour, 25-step forecast, and extracts a point-level time series for
any lat/lon (default: Berlin).

Product facts
-------------
  Grid    : DE1200 — 1200 rows × 1100 cols, 1 km resolution
  Steps   : 25 frames × 5 min = 0 … 120 min ahead
  Unit    : mm / 5 min  (rain accumulation per interval)
  Updated : every 5 min; files appear ~2-3 min after nominal time
  URL     : https://opendata.dwd.de/weather/radar/radolan/rv/
  Filename: DE1200_RV{YYMMDDHHMM}.tar.bz2
  Inside  : 25 binary files, one per forecast step

Dependencies
------------
  pip install wradlib requests numpy pyproj
  (optional) pip install matplotlib   — for the plot
"""

import csv
import io
import re
import tarfile
from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import numpy as np
import requests
import wradlib as wrl

TZ_CEST = ZoneInfo("Europe/Berlin")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_URL   = "https://opendata.dwd.de/weather/radar/composite/rv/"
TARGET_LAT = 52.50 # Berlin
TARGET_LON =  13.44

# RADOLAN RV grid dimensions
NROWS, NCOLS = 1200, 1100


# ---------------------------------------------------------------------------
# 1. Discover the latest available RV file
# ---------------------------------------------------------------------------

def list_rv_files() -> list[str]:
    """Return all DE1200_RV*.tar.bz2 filenames from the DWD index."""
    resp = requests.get(BASE_URL, timeout=30)
    resp.raise_for_status()
    return re.findall(r'DE1200_RV\d{10}\.tar\.bz2', resp.text)


def latest_rv_url() -> str:
    files = sorted(list_rv_files())          # lexicographic sort == time sort
    if not files:
        raise RuntimeError("No RV files found at DWD Open Data")
    fname = files[-1]
    print(f"Latest RV file: {fname}")
    return BASE_URL + fname, fname


# ---------------------------------------------------------------------------
# 2. Download and unpack the tar.bz2 in memory
# ---------------------------------------------------------------------------

def download_rv(url: str) -> list[io.BytesIO]:
    """Download the RV tar.bz2 and return 25 BytesIO objects (one per step)."""
    print(f"Downloading … {url}")
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    raw = io.BytesIO(resp.content)

    buffers = []
    with tarfile.open(fileobj=raw, mode="r:bz2") as tf:
        members = sorted(tf.getnames())       # sorted = chronological
        for name in members:
            member_bytes = tf.extractfile(name).read()
            buf = io.BytesIO(member_bytes)
            buf.name = name                   # wradlib reads the name for metadata
            buffers.append(buf)

    print(f"  → {len(buffers)} forecast steps unpacked")
    return buffers


# ---------------------------------------------------------------------------
# 3. Parse each binary frame with wradlib
# ---------------------------------------------------------------------------

def parse_rv_frames(buffers: list[io.BytesIO]) -> tuple[list, list[datetime]]:
    """
    Return (frames, valid_times).

    frames      : list of 2-D numpy arrays, shape (NROWS, NCOLS), unit mm/5min
    valid_times : list of UTC datetimes for each frame
    """
    frames      = []
    valid_times = []

    for buf in buffers:
        data, attrs = wrl.io.read_radolan_composite(buf)

        # Mask nodata / clutter BEFORE scaling to avoid float comparison issues
        nodata = attrs.get("nodataflag", -9999)
        data = np.ma.masked_equal(data, nodata)

        # Apply precision factor (raw values are integers scaled by 1/precision)
        precision = attrs.get("precision", 0.01)
        data = data.astype(float) * precision

        # Clamp: rain rate cannot be negative
        data = np.ma.where(data < 0, 0.0, data)

        # Secondary (below-detection) pixels → 0 mm
        sec = attrs.get("secondary", np.array([], dtype=int))
        data.flat[sec] = 0.0

        frames.append(data)

        # Compute valid time from base time + prediction offset, convert to CEST
        base_time = attrs["datetime"].replace(tzinfo=timezone.utc)
        offset_min = attrs.get("predictiontime", 0)
        valid_times.append((base_time + timedelta(minutes=offset_min)).astimezone(TZ_CEST))

    return frames, valid_times


# ---------------------------------------------------------------------------
# 4. Build the 1200 × 1100 lat/lon grid  (RV uses the larger DE1200 grid)
# ---------------------------------------------------------------------------

def get_rv_latlon_grid() -> tuple[np.ndarray, np.ndarray]:
    """
    Return (lat, lon) arrays of shape (NROWS, NCOLS) in WGS-84.

    The RV product uses the extended German composite grid (DE1200).
    wradlib's get_radolan_grid() only covers the classic 900×900 grid,
    so we compute coordinates manually using pyproj.
    """
    from pyproj import Proj, Transformer

    # Origin (lower-left corner of pixel [0,0]) in RADOLAN stereographic km
    # From DWD Kompositformatbeschreibung for DE1200:
    x0, y0 = -673.4661, -5008.6445  # km

    # 1 km pixel spacing
    xs = x0 + np.arange(NCOLS) * 1.0   # (NCOLS,)
    ys = y0 + np.arange(NROWS) * 1.0   # (NROWS,)
    xx, yy = np.meshgrid(xs, ys)        # (NROWS, NCOLS)

    # RADOLAN polar-stereographic (in km, not m)
    proj_rad = Proj(
        proj="stere", lat_0=90, lat_ts=90, lon_0=10,
        k=0.93301270189,
        x_0=0, y_0=0,
        a=6370040, b=6370040,
        units="km"
    )
    proj_wgs = Proj(proj="latlong", datum="WGS84")

    transformer = Transformer.from_proj(proj_rad, proj_wgs, always_xy=True)
    lon, lat = transformer.transform(xx, yy)
    return lat, lon


# ---------------------------------------------------------------------------
# 5. Extract nearest grid point for a given lat/lon
# ---------------------------------------------------------------------------

def nearest_pixel(lat_grid, lon_grid, target_lat, target_lon) -> tuple[int, int]:
    dist = (lat_grid - target_lat) ** 2 + (lon_grid - target_lon) ** 2
    row, col = np.unravel_index(np.argmin(dist), dist.shape)
    return int(row), int(col)


def extract_point_series(
    frames: list,
    valid_times: list[datetime],
    row: int,
    col: int,
) -> list[dict]:
    """Return list of {time, mm_per_5min, mm_per_hour} for the target pixel."""
    results = []
    for frame, t in zip(frames, valid_times):
        raw = frame[row, col]
        val = 0.0 if np.ma.is_masked(raw) else float(raw)
        results.append({
            "valid_time":    t,
            "mm_per_5min":   val,
            "mm_per_hour":   val * 12,   # × 12 to convert 5-min accumulation
        })
    return results


# ---------------------------------------------------------------------------
# 6. Optional: pretty-print and plot
# ---------------------------------------------------------------------------

def print_forecast(series: list[dict], location_name: str):
    print(f"\n{'─'*58}")
    print(f"  RV nowcast for {location_name}")
    print(f"{'─'*58}")
    print(f"  {'Valid time (CEST)':<22}  {'mm/5min':>8}  {'mm/h':>8}")
    print(f"{'─'*58}")
    for row in series:
        mm_h = row["mm_per_hour"]
        bar = "" if np.isnan(mm_h) else "█" * min(int(mm_h * 2), 20)
        print(
            f"  {row['valid_time'].strftime('%Y-%m-%d %H:%M'):<22}"
            f"  {row['mm_per_5min']:>8.3f}"
            f"  {row['mm_per_hour']:>8.2f}  {bar}"
        )
    print(f"{'─'*58}\n")


def save_csv(series: list[dict], base_time_str: str):
    """Save the forecast series to a CSV file named by base time."""
    csv_path = f"rv_nowcast_{base_time_str}.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["valid_time_CEST", "mm_per_5min", "mm_per_hour"])
        for row in series:
            writer.writerow([
                row["valid_time"].strftime("%Y-%m-%d %H:%M"),
                f"{row['mm_per_5min']:.3f}",
                f"{row['mm_per_hour']:.2f}",
            ])
    print(f"CSV saved  → {csv_path}")


def plot_forecast(series: list[dict], location_name: str, base_time_str: str):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed — skipping plot")
        return

    times  = [r["valid_time"] for r in series]
    values = [0.0 if np.isnan(v) else v for v in (r["mm_per_hour"] for r in series)]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(times, values, width=0.003, color="#1D9E75", alpha=0.85, label="Rain rate")
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Valid time (CEST)")
    ax.set_ylabel("Rain rate  [mm h⁻¹]")
    ax.set_title(f"DWD RADOLAN RV — 2-hour nowcast\n{location_name}  "
                 f"({TARGET_LAT:.4f}°N, {TARGET_LON:.4f}°E)")
    ax.legend()
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plot_path = f"rv_nowcast_{base_time_str}.png"
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved → {plot_path}")
    plt.show()


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------

def main():
    # --- Download ---
    url, fname = latest_rv_url()
    buffers = download_rv(url)

    # --- Parse ---
    frames, valid_times = parse_rv_frames(buffers)
    print(f"\nForecast base time : {valid_times[0].strftime('%Y-%m-%d %H:%M CEST')}")
    print(f"Forecast end  time : {valid_times[-1].strftime('%Y-%m-%d %H:%M CEST')}")
    print(f"Steps              : {len(frames)}  (every 5 min)")
    print(f"Grid               : {frames[0].shape[0]} rows × {frames[0].shape[1]} cols")

    # --- Grid ---
    print("\nBuilding lat/lon grid …")
    lat_grid, lon_grid = get_rv_latlon_grid()

    # --- Locate target ---
    row, col = nearest_pixel(lat_grid, lon_grid, TARGET_LAT, TARGET_LON)
    actual_lat = lat_grid[row, col]
    actual_lon = lon_grid[row, col]
    print(f"\nTarget  : {TARGET_LAT:.4f}°N, {TARGET_LON:.4f}°E  (Berlin)")
    print(f"Nearest : row={row}, col={col}  →  {actual_lat:.4f}°N, {actual_lon:.4f}°E")

    # --- Extract & display ---
    series = extract_point_series(frames, valid_times, row, col)
    base_time_str = valid_times[0].strftime("%Y%m%d_%H%M")
    print_forecast(series, "Berlin")
    save_csv(series, base_time_str)
    plot_forecast(series, "Berlin", base_time_str)

    return series   # useful when calling from a notebook


if __name__ == "__main__":
    main()