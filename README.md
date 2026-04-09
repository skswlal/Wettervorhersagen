# Wettervorhersagen

## DWD RADOLAN RV Nowcast — Download & Parse

Fetches the latest RV composite from DWD Open Data, parses the full 2-hour, 25-step forecast, and extracts a point-level time series for any latitude/longitude (default: Berlin).

### Product Facts

- **Grid**: DE1200 — 1200 rows × 1100 cols, 1 km resolution
- **Steps**: 25 frames × 5 min = 0 … 120 min ahead
- **Unit**: mm / 5 min (rain accumulation per interval)
- **Updated**: Every 5 minutes; files appear ~2-3 min after nominal time
- **Source**: https://opendata.dwd.de/weather/radar/radolan/rv/
- **Filename**: `DE1200_RV{YYMMDDHHMM}.tar.bz2`

### Installation

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd Wettervorhersagen
   ```

2. **Create a Python virtual environment** (Python 3.9+):
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

### Usage

Run the main script to download the latest RV nowcast and extract a forecast for Berlin:

```bash
python dwd_rv_nowcast.py
```

The script will:
1. Download the latest RV composite from DWD Open Data
2. Parse 25 forecast steps
3. Extract point-level rain forecast for Berlin
4. **Output**:
   - Display forecast in the terminal
   - Save CSV: `rv_nowcast_YYYYMMDD_HHMM.csv`
   - Save PNG plot: `rv_nowcast_YYYYMMDD_HHMM.png` (if matplotlib is installed)

### Customization

Edit the configuration section in `dwd_rv_nowcast.py` to change the target location:

```python
TARGET_LAT = 52.50  # Berlin
TARGET_LON = 13.44
```

### Dependencies

- **wradlib** — RADOLAN file format parsing
- **requests** — HTTP downloads
- **numpy** — Array operations
- **pyproj** — Coordinate transformations (stereographic → WGS-84)
- **matplotlib** — Optional, for plotting
