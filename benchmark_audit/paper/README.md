# LaTeX manuscript package

This is a substantially revised computational audit prepared for *Interdisciplinary Sciences: Computational Life Sciences*. It replaces unsupported original model-superiority claims with newly executed public-data construction experiments. Read `author_only/review_report.pdf` for the claim-by-claim disposition and remaining factual author checks.

## Documents

- `manuscript.tex` / `manuscript.pdf`: anonymous main manuscript in the Springer Nature class.
- `supplementary.tex` / `supplementary.pdf`: methods, sensitivity analyses and full retention tables.
- `references.bib`, `sn-basic.bst`, `sn-jnl.cls`: bibliography and publisher class.
- `figures/`: five main and one supplementary scientific figure, in vector PDF and EPS plus PNG previews.
- `tables/`: generated supplementary tables.
- `author_only/title_page.tex`: separate title page; author-specific declarations require factual completion.
- `author_only/review_report.pdf`: adversarial reviewer debate, withdrawn claims, evidence ledger and submission constraints.

## Compile

From this directory, using a TeX Live installation with `latexmk`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary.tex
```

The archive includes generated `.bbl` files as well as the editable `.bib` source. No neural model or external data download is needed to compile the documents. Plot and table scripts are included under `analysis/`, but regenerating them requires the separate full reproducibility package with its result folders.

## Submission handling

The journal uses double-anonymous review. Upload the main manuscript and title page in their respective categories. `author_only/` is not an anonymous reviewer supplement. The full reproducibility archive contains identifying source licenses and provenance; do not upload that author archive wholesale to anonymous reviewers. Preserve copyright notices when preparing any review supplement.

The factual author declarations and publication status are not invented. The source README contains an earlier statement about review elsewhere, which the authors must check. All authors need to approve the substantial AI-assisted analyses and revised manuscript. No journal submission has been made.

The publisher offers subscription publication without an APC; optional open access is chargeable. The scientific scope is benchmark construction, not experimentally verified biological nonbinding, clinical utility or guaranteed acceptance.
