from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

HERE = Path(__file__).resolve().parent

EXCEL_FILE = HERE / "master_data.xlsx"
OUTPUT_DIR = HERE / "Histogram_Same_Bins_Results"
OUTPUT_DIR.mkdir(exist_ok=True)

# Load data
df = pd.read_excel(EXCEL_FILE)

# Clean data
df["spectral_class"] = df["spectral_class"].astype(str).str.strip().str.upper()
df["period_d"] = pd.to_numeric(df["period_d"], errors="coerce")

df = df.dropna(subset=["period_d"])
df = df[df["period_d"] > 0]

# Select spectral classes and random sample
classes = ["B", "A", "F", "G", "K", "M", "RANDOM"]

df = df[df["spectral_class"].isin(classes)].copy()

# Calculate logarithmic periods
df["log_period"] = np.log10(df["period_d"])

# Separate samples
samples = {
    spectral_class: df[df["spectral_class"] == spectral_class]["log_period"]
    for spectral_class in classes
}

# Print sample sizes
print("\nNumber of periods:")
for spectral_class, periods in samples.items():
    print(f"{spectral_class}: {len(periods)}")

# Use identical bins for every histogram
all_periods = df["log_period"]

bins = np.linspace(
    all_periods.min(),
    all_periods.max(),
    9
)

# Colors
colors = {
    "B": "steelblue",
    "A": "darkorange",
    "F": "forestgreen",
    "G": "purple",
    "K": "firebrick",
    "M": "teal",
    "RANDOM": "gray"
}

# Individual histograms
for spectral_class, periods in samples.items():

    plt.figure(figsize=(8, 5))

    plt.hist(
        periods,
        bins=bins,
        color=colors[spectral_class],
        edgecolor="black"
    )

    plt.xlabel("log10(Orbital Period [days])")
    plt.ylabel("Number of Binary Systems")
    plt.title(f"Orbital Period Distribution - {spectral_class}")

    plt.tight_layout()

    plt.savefig(
        OUTPUT_DIR / f"Histogram_{spectral_class}.png",
        dpi=300
    )

    plt.close()

# Combined normalized histogram
plt.figure(figsize=(10, 6))

for spectral_class, periods in samples.items():

    plt.hist(
        periods,
        bins=bins,
        density=True,
        alpha=0.35,
        color=colors[spectral_class],
        edgecolor="black",
        label=spectral_class
    )

plt.xlabel("log10(Orbital Period [days])")
plt.ylabel("Probability Density")
plt.title("Orbital Period Distributions by Spectral Class")
plt.legend(title="Spectral Class")

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "Combined_Histogram.png",
    dpi=300
)

plt.close()

print("\nAll histograms successfully created!")
print(f"Graphs saved in: {OUTPUT_DIR}")