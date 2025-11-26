import pandas as pd

CSR_FILE = "report_CSR.csv"
CSB_FILE = "report_CSB.csv"
OUT_FILE = "merged_report.csv"

csr = pd.read_csv(CSR_FILE)
csb = pd.read_csv(CSB_FILE)

csr.columns = csr.columns.str.strip().str.lower()
csb.columns = csb.columns.str.strip().str.lower()

csr["format"] = "CSR"
csb["format"] = "CSB"

csr["mode"] = csr["mode"].str.upper()
csb["mode"] = csb["mode"].str.upper()

csr["threads"] = pd.to_numeric(csr["threads"], errors="coerce")
csr["chunk"] = pd.to_numeric(csr["chunk"], errors="coerce")
csb["threads"] = pd.to_numeric(csb["threads"], errors="coerce")
csb["chunk"] = pd.to_numeric(csb["chunk"], errors="coerce")

merged = pd.concat([csr, csb], ignore_index=True)

merged.to_csv(OUT_FILE, index=False)
print(f"Generated: {OUT_FILE}")
