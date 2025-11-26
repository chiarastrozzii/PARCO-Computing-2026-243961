import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set(style="whitegrid")

INPUT = "merged_report.csv"
OUTDIR = "graphs"

os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(INPUT)

df["matrix_clean"] = df["matrix"].str.replace(r"[\\/ ]", "_", regex=True)
df["threads"] = pd.to_numeric(df["threads"], errors="coerce")
df["chunk"] = pd.to_numeric(df["chunk"], errors="coerce")
df["time_ms"] = pd.to_numeric(df["time_ms"], errors="coerce")

df["schedule"] = df["schedule"].fillna("none")
df["chunk"] = df["chunk"].fillna(0)

csr_seq = df[(df["format"]=="CSR") & (df["mode"]=="SEQ")].set_index("matrix")["time_ms"]
csb_seq = df[(df["format"]=="CSB") & (df["mode"].str.contains("SEQ"))].set_index("matrix")["time_ms"]

def compute_speedup(row):
    m = row["matrix"]
    if row["format"] == "CSR":
        if m in csr_seq.index:
            return csr_seq[m] / row["time_ms"]
    if row["format"] == "CSB":
        if m in csb_seq.index:
            return csb_seq[m] / row["time_ms"]
    return np.nan

df["speedup"] = df.apply(compute_speedup, axis=1)

comparison = df[df["mode"].str.contains("PAR")]

plt.figure(figsize=(12,6))
sns.barplot(data=comparison, x="matrix_clean", y="time_ms", hue="format")
plt.xticks(rotation=45)
plt.title("CSR vs CSB Parallel Time")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/csr_vs_csb_parallel_time.png")
plt.close()

plt.figure(figsize=(12,6))
sns.barplot(data=comparison, x="matrix_clean", y="speedup", hue="format")
plt.xticks(rotation=45)
plt.title("CSR vs CSB Speedup")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/csr_vs_csb_speedup.png")
plt.close()

matrices = df["matrix"].unique()

for matrix in matrices:
    sub = df[df["matrix"] == matrix]
    mat_clean = sub["matrix_clean"].iloc[0]

    plt.figure(figsize=(10,6))
    sns.lineplot(data=sub[sub["mode"].str.contains("PAR")],
                 x="threads", y="time_ms", hue="format", marker="o")
    plt.title(f"Time vs Threads — {matrix}")
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/{mat_clean}_time_threads.png")
    plt.close()

    plt.figure(figsize=(10,6))
    sns.lineplot(data=sub[sub["mode"].str.contains("PAR")],
                 x="threads", y="speedup", hue="format", marker="o")
    plt.title(f"Speedup vs Threads — {matrix}")
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/{mat_clean}_speedup_threads.png")
    plt.close()

    plt.figure(figsize=(10,6))
    sns.lineplot(
        data=sub[(sub["mode"].str.contains("PAR")) & (sub["schedule"]!="none")],
        x="threads", y="time_ms", hue="schedule", style="format", marker="o"
    )
    plt.title(f"Schedule Comparison — {matrix}")
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/{mat_clean}_schedule_comparison.png")
    plt.close()

    plt.figure(figsize=(10,6))
    sns.lineplot(
        data=sub[(sub["mode"].str.contains("PAR")) & (sub["chunk"]!=0)],
        x="chunk", y="time_ms", hue="format", style="threads", marker="o"
    )
    plt.title(f"Chunk Sensitivity — {matrix}")
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/{mat_clean}_chunk_sensitivity.png")
    plt.close()

    heat = sub[(sub["mode"].str.contains("PAR"))]
    if not heat.empty:
        pivot = heat.pivot_table(
            index="schedule", columns="chunk", values="time_ms"
        )
        plt.figure(figsize=(10,6))
        sns.heatmap(pivot, annot=True, cmap="viridis")
        plt.title(f"Schedule × Chunk Heatmap — {matrix}")
        plt.tight_layout()
        plt.savefig(f"{OUTDIR}/{mat_clean}_heatmap.png")
        plt.close()


TARGET_MATRIX = "matrix_mtk_bcsprw10.mtx"  
THREADS = 32                           
FIXED_CHUNK = 10                              

sub = df[
    (df["matrix"] == TARGET_MATRIX) &
    (df["threads"] == THREADS) &
    (df["chunk"] == FIXED_CHUNK) &
    (df["mode"].str.contains("PAR")) &
    (df["schedule"] != "none")
]

if not sub.empty:
    plt.figure(figsize=(10,6))
    sns.barplot(data=sub, x="schedule", y="time_ms", hue="format")
    plt.title(f"Scheduling Policies Comparison — {TARGET_MATRIX}\nthreads={THREADS}, chunk={FIXED_CHUNK}")
    plt.xlabel("Scheduling Policy")
    plt.ylabel("Runtime (ms)")
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/{sub['matrix_clean'].iloc[0]}_schedule_fixed_threads{THREADS}_chunk{FIXED_CHUNK}.png")
    plt.close()
else:
    print("No data available for selected matrix/chunk/threads combination.")

print("All graphs generated inside the 'graphs/' folder.")
