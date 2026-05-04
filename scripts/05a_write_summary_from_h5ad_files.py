from pathlib import Path
import importlib.util

import anndata as ad
import pandas as pd

def main():
    REPO_ROOT = Path(".").resolve()
    QC_DIR = REPO_ROOT / "results" / "qc"

    spec = importlib.util.spec_from_file_location("qc", REPO_ROOT / "scripts" / "05_qc.py")
    qc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(qc)

    summary_rows = []

    for h5ad_path in sorted(QC_DIR.glob("*.h5ad")):
        adata = ad.read_h5ad(h5ad_path)
        summary_rows.append(qc.summarize_qc(h5ad_path.stem, adata))

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(QC_DIR / "qc_summary.tsv", sep="\t", index=False)

if __name__== "__main__":
    main()
