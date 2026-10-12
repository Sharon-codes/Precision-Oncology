# TCR benchmark construction audit

**Current manuscript:** *Batch-dependent negative sampling and marginal-preserving decoy construction in TCR–peptide benchmarks* (12 October 2026).

The completed revision is an executable computational audit of decoy construction using frozen public VDJdb records. It reports batch/order-dependent registry contradictions and antigen-frequency effects, a constrained closed-assignment correction, its coverage cost, and receptor-grouped metadata controls. It does not establish biological nonbinding or validate the original neural predictor's superiority claims.

## Paper and complete evidence

- [Main manuscript PDF](benchmark_audit/paper_netmahib/manuscript.pdf) · [editable LaTeX](benchmark_audit/paper_netmahib/manuscript.tex)
- [Supplementary methods and results](benchmark_audit/paper_netmahib/ESM_1.pdf)
- [Author review and evidence ledger](benchmark_audit/paper/author_only/review_report.pdf)
- [Complete reproduction instructions](benchmark_audit/README.md)
- [Three fresh reviewer reports](benchmark_audit/reviews/) and independent numerical checks
- [Vector figures and rendered figures](benchmark_audit/paper_netmahib/figures/)

`benchmark_audit/` contains all new code, original source data archive with license, row-index construction manifests, actual assignment witnesses, aggregate measurements, out-of-fold predictions, simulation controls and manuscript sources. Large prediction files are losslessly chunked for transport and automatically restored with checksum verification.

```bash
cd benchmark_audit
python -m venv .venv
source .venv/bin/activate
pip install -r requirements_reproducibility.txt
python run_all.py --verify-only
```

Run `python run_all.py --compile-pdfs` to regenerate the full research pipeline and PDFs after installing a TeX distribution. The numerical verification checks archived predictions; it is distinct from full retraining. The final independent review freshly refitted a selected subset, and separately checked all archived numerical metrics.

## Historical source and submission status

The original repository-root implementation is retained as historical source. The revised audit pins original commit `0ba2d538065e2d3a8b4642547eb77447ecc89a5b`, also supplied under `benchmark_audit/precision_repo/` with its original notices. Earlier manuscript performance claims are not endorsed by the new audit; the evidence ledger documents their disposition.

Prepared for consideration by *Network Modeling Analysis in Health Informatics and Bioinformatics* under the subscription route. The author-identified Springer manuscript includes Sharon Melhi and Rajat Bhardwaj of Amity School of Engineering and Technology, Amity University Bengaluru. See [journal packaging and author instructions](benchmark_audit/paper_netmahib/README.md).

No journal submission is represented by this update. Funding, competing interests and author contributions are finalized. Sharon Melhi is the corresponding author (sharonmelhi365@gmail.com). Coauthor approval and submission eligibility are confirmed by the corresponding author during journal upload. AI reviews are internal stress tests, not journal peer review. Scientific measurements are unchanged from research revision `a4162ded50ba24d76ab8d5134577eba6935f910b`; the updated manuscript uses author–year references and concise limitations with detailed implementation discussion in the supplement. The original anonymous manuscript remains in `benchmark_audit/paper/` for provenance.
