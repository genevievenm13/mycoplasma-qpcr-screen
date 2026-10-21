"""
Mycoplasma contamination screen of cell-culture samples (Nt, CFL) by SYBR qPCR + melt-curve analysis.

Data: Bio-Rad CFX Maestro exports from a run on 2026-09-11 (protocol MycoSensor_Protocol.prcl).
Interpretation rules follow the Agilent MycoSensor QPCR Assay Kit manual (#302106):
  - Mycoplasma amplicon melts ~82 C, amplification-control (AC) amplicon ~85 C, primer-dimer ~74 C
    (values given for an Mx3000P; peak positions shift by ~1-2 C between instruments, so this script
    calibrates the AC peak on this instrument and expects the Mycoplasma peak ~3 C below it).

Usage:  python scripts/fix_xlsx.py data/raw data/clean   (once, repairs the CFX exports)
        python scripts/analyze_qpcr.py
Writes: results/calls.csv, figures/*.png
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[1]
DATA, FIG, RES = ROOT / "data" / "clean", ROOT / "figures", ROOT / "results"
FIG.mkdir(exist_ok=True); RES.mkdir(exist_ok=True)

# ---------- load ----------
def table(name):
    df = pd.read_excel(DATA / name, header=0)
    return df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]

def wide(name, xcol):
    df = table(name).set_index(xcol)
    df.columns = [c[0] + c[1:].zfill(2) for c in df.columns]  # A1 -> A01
    return df

cq    = table("Quantification_Cq_Results.xlsx")[["Well", "Sample", "Cq"]]
peaks = table("Melt_Curve_Peak_Results.xlsx")[["Well", "Melt Temperature", "Peak Height"]]
endp  = table("End_Point_Results.xlsx")[["Well", "End RFU"]]
amp   = wide("Quantification_Amplification_Results.xlsx", "Cycle")
deriv = wide("Melt_Curve_Derivative_Results.xlsx", "Temperature")

df = cq.merge(peaks, on="Well").merge(endp, on="Well")
df["Source"] = df["Sample"].str.split().str[0]
role_map = {"": "Test sample", "amp": "Sample + amp. control",
            "Myco": "Positive control", "Neg": "Negative control"}
df["Reaction"] = df["Sample"].str.split().str[1:].str.join(" ").map(role_map)

# ---------- melt-peak classification ----------
# Calibrate on this instrument: AC peak observed in both "+ amp" wells.
ac_tm = df.loc[df["Reaction"] == "Sample + amp. control", "Melt Temperature"].mean()
myco_tm = ac_tm - 3.0            # kit: Mycoplasma ~82 C vs AC ~85 C
def classify(tm):
    if tm >= ac_tm - 1.0:  return "AC product"
    if tm >= myco_tm - 1.5: return "Mycoplasma-size product"
    return "Primer-dimer"
df["Peak class"] = df["Melt Temperature"].apply(classify)

# The CFX summary reports only the dominant peak. Scan each derivative curve for every peak
# >= 15% of the well's tallest peak, so a minor target peak hiding behind primer-dimer is not missed.
def all_peaks(well):
    y = deriv[well].loc[68:]
    idx, _ = find_peaks(y.values, prominence=0.15 * y.max())
    return [(float(y.index[i]), float(y.iloc[i])) for i in idx]
df["All peaks"] = df["Well"].apply(lambda w: "; ".join(f"{t:.1f} C ({h:.0f})" for t, h in all_peaks(w)))
# A minor target product often shows only as a shoulder (no local maximum), so also measure the
# derivative signal inside the Mycoplasma window. Test wells with no target read ~0 there.
df["Myco-window signal"] = df["Well"].apply(lambda w: float(deriv[w].loc[myco_tm - 1.5: myco_tm + 0.5].max()))
SHOULDER_MIN = 50   # -dRFU/dT; test-sample wells read < 5 in this window
df["Myco-size peak"] = df["Myco-window signal"] >= SHOULDER_MIN

def call(row):
    r, p = row["Reaction"], row["Peak class"]
    if r == "Test sample":
        return "Mycoplasma detected" if row["Myco-size peak"] else "No Mycoplasma product"
    if r == "Sample + amp. control":
        return {"AC product": "AC amplified - no inhibition",
                "Mycoplasma-size product": "Heavy infection (AC suppressed)"}.get(p, "Inhibited")
    if r == "Positive control":
        if p == "Mycoplasma-size product": return "PASS"
        return "WEAK - target shoulder ~80 C, primer-dimer dominant" if row["Myco-size peak"] else "FAIL - no target product"
    if r == "Negative control":
        return "PASS (primer-dimer only)" if p == "Primer-dimer" else f"FAIL - contaminated ({p})"
df["Call"] = df.apply(call, axis=1)

order = ["Test sample", "Sample + amp. control", "Positive control", "Negative control"]
df["Reaction"] = pd.Categorical(df["Reaction"], order, ordered=True)
df = df.sort_values(["Source", "Reaction"], ascending=[False, True])
out = df[["Well", "Source", "Reaction", "Cq", "Melt Temperature", "Peak Height", "End RFU", "Peak class", "All peaks", "Myco-window signal", "Call"]]
out.round(2).to_csv(RES / "calls.csv", index=False)
print(f"Instrument-calibrated AC Tm = {ac_tm:.1f} C; expected Mycoplasma Tm ~ {myco_tm:.1f} C\n")
print(out.round(2).to_string(index=False))

# ---------- figures ----------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#8a8984", "axes.labelcolor": "#52514e",
                     "xtick.color": "#52514e", "ytick.color": "#52514e", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})
COLORS = dict(zip(order, ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]))
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"

def panel_lines(ax, data, wells, xlabel, ylabel, title):
    for _, row in df[df["Well"].isin(wells)].iterrows():
        ax.plot(data.index, data[row["Well"]], lw=2, color=COLORS[row["Reaction"]], label=row["Reaction"])
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title, loc="left", color=INK, fontsize=11, fontweight="bold")

sources = ["Nt", "CFL"]
wells = {s: df.loc[df["Source"] == s, "Well"].tolist() for s in sources}

# 1) amplification curves
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, s in zip(axes, sources):
    panel_lines(ax, amp, wells[s], "Cycle", "Fluorescence (RFU, baseline-subtracted)" if s == "Nt" else "", f"{s} sample set")
    for _, row in df[df["Source"] == s].iterrows():
        ax.plot(row["Cq"], np.interp(row["Cq"], amp.index, amp[row["Well"]]), "o", ms=7,
                color=COLORS[row["Reaction"]], mec="#fcfcfb", mew=2)
axes[1].legend(frameon=False, loc="upper left", fontsize=9)
fig.suptitle("Amplification curves (dots = Cq)", x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(); fig.savefig(FIG / "amplification_curves.png", dpi=200); plt.close(fig)

# 2) melt-curve derivatives with expected-product bands
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
bands = [(70, 76, "Primer-dimer"), (myco_tm - 1.5, ac_tm - 1.0, "Myco"), (ac_tm - 1.0, ac_tm + 1.5, "AC")]
BAND_FILL = {"Primer-dimer": "#efeee9", "Myco": "#e4ece6", "AC": "#efeee9"}
for ax, s in zip(axes, sources):
    for lo, hi, lab in bands:
        ax.axvspan(lo, hi, color=BAND_FILL[lab], zorder=0)
        ax.text((lo + hi) / 2, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, color=MUTED)
    panel_lines(ax, deriv.loc[68:], wells[s], "Temperature (°C)", "-d(RFU)/dT" if s == "Nt" else "", "")
    ax.set_title(f"{s} sample set", loc="left", color=INK, fontsize=11, fontweight="bold", pad=18)
axes[1].legend(frameon=False, loc="upper right", fontsize=9)
fig.suptitle("Melt-curve derivatives: which product formed in each well", x=0.01, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(); fig.savefig(FIG / "melt_curves.png", dpi=200); plt.close(fig)

# 3) Cq vs melt temperature map
fig, ax = plt.subplots(figsize=(9.5, 4.6))
for lo, hi, lab in bands:
    ax.axvspan(lo, hi, color=BAND_FILL[lab], zorder=0)
    ax.text((lo + hi) / 2, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, color=MUTED)
for _, row in df.iterrows():
    mk = "o" if row["Source"] == "Nt" else "s"
    ax.plot(row["Melt Temperature"], row["Cq"], mk, ms=10, color=COLORS[row["Reaction"]], mec="#fcfcfb", mew=2)
    off = {"Nt amp": (8, 6), "CFL amp": (8, -12), "Nt Myco": (8, 4), "CFL Myco": (-62, -4), "CFL": (-30, -4)}.get(row["Sample"], (8, -3))
    ax.annotate(row["Sample"], (row["Melt Temperature"], row["Cq"]), xytext=off, textcoords="offset points", fontsize=8.5, color=MUTED)
ax.invert_yaxis(); ax.set_xlim(70, 87)
ax.set_xlabel("Melt peak (°C)"); ax.set_ylabel("Cq  (higher on chart = earlier = more template)")
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="", color=COLORS[r], ms=8, label=r) for r in order]
h += [Line2D([], [], marker="o", ls="", color=MUTED, ms=8, label="Nt set"), Line2D([], [], marker="s", ls="", color=MUTED, ms=8, label="CFL set")]
ax.legend(handles=h, frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(1.01, 0.5))
ax.set_title("Every well: which product (x) and how early it amplified (y)", loc="left", color=INK, fontsize=11, fontweight="bold", pad=22)
fig.tight_layout(); fig.savefig(FIG / "cq_vs_tm.png", dpi=200); plt.close(fig)
print("\nFigures written to", FIG)
