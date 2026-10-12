# NetMAHIB manuscript package

Preferred journal: **Network Modeling Analysis in Health Informatics and Bioinformatics** (Springer Nature). Article type: original research. Choose **subscription publication**, for which the official journal page states that no APC applies. Open access is optional and chargeable.

## Files and compilation

- `manuscript.pdf` / `manuscript.tex`: author-identified scientific paper, Springer class with author–year references.
- `ESM_1.pdf` / `ESM_1.tex`: Online Resource 1, supplementary methods and results.
- `title_page.pdf` / `title_page.tex`: separate author details if requested by the portal.
- `figures/`: identical research figures in vector PDF/EPS and PNG; publisher-preferred filenames.
- `tables/`, `references.bib`, `sn-jnl.cls`, `sn-basic.bst`: editable tables and build dependencies.
- `AUTHOR_DECLARATIONS.md`: completed funding, competing-interest and contribution statements.
- `cover_letter.md`: editable unsent cover-letter draft.
- `JOURNAL_DECISION.md`: independent-agent debate, rationale and remaining fit risk.
- `FORMAT_PROVENANCE.json`: exact input and unchanged scientific-asset fingerprints.

Run `latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex`, then the same command for `ESM_1.tex` and `title_page.tex`.

The code and data remain in the complete research repository, pinned at revision `a4162ded50ba24d76ab8d5134577eba6935f910b`. Figure-generation scripts here require the repository's archived analysis inputs; `benchmark_audit/README.md` explains reproduction. This source package is not a substitute for the complete research archive.

## Current preparation status

The manuscript and supplement carry Sharon Melhi and Rajat Bhardwaj, with the supplied Amity affiliation. Sharon Melhi is the corresponding author at sharonmelhi365@gmail.com. The 218-word abstract, six keywords, author–year references and numbered sections follow the retrieved journal instructions. Figures were created programmatically; no generative image was used.

Funding, competing interests and author contributions have been completed from the author's instructions. The journal also requires contribution and competing-interest details in its submission interface. At upload, the corresponding author must confirm coauthor approval and that the work is not under consideration elsewhere. No submission has been made.

No journal-specific double-anonymous requirement was stated in the retrieved instructions; author details are requested. Follow any additional live portal requirement when submitting. No submission or contact with an editor has been made.

## Official sources checked 12 October 2026

- https://link.springer.com/journal/13721/aims-and-scope
- https://link.springer.com/journal/13721/how-to-publish-with-us
- https://link.springer.com/journal/13721/submission-guidelines

The scientific results and limitations are unchanged from the completed research revision. This update changes author information, reference format, access/provenance wording and journal packaging.
