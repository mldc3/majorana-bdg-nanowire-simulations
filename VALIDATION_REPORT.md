# Validation Report

This report documents the final quality-control pass for the GitHub-ready version of the Majorana BdG nanowire repository.

## Scope of the final package

The repository contains three progressively more complete models:

1. the one-dimensional Kitaev chain;
2. the self-consistent Bogoliubov–de Gennes nanowire;
3. the realistic proximitized semiconductor–superconductor nanowire with spatial proximity profiles and screened Hartree electrostatics.

For every model, the approved long-form documentation is written directly in GitHub-compatible Markdown. The mathematical expressions are native LaTeX source between Markdown math delimiters, not equation screenshots. The archived numerical graphs are embedded at the points where they are introduced and analysed. The supplied Python source files are retained separately and unchanged.

No PDF files are included. No `docs/` directory is used. The only images included in the main model documentation are numerical graphs.

## Automated checks

The following checks passed in the final validation run:

- all required root, model, theory, source-code, results and supplementary files are present and non-empty;
- all Markdown documents parse successfully with Pandoc's GitHub-Flavoured Markdown reader and TeX-math support;
- all display-math delimiters and fenced code blocks are balanced;
- no unresolved LaTeX document commands such as labels, equation references, citation commands, figure environments or document environments remain;
- all Markdown image paths resolve to existing PNG files;
- all 37 principal result graphs occur exactly once and in the approved order in their corresponding results documents;
- the two additional theory figures in the realistic model resolve to the same original graph files used in the results section;
- all principal graph files are byte-for-byte identical to the archived source figures;
- all three Python files compile successfully;
- all three Python files are byte-for-byte identical to their supplied source versions;
- the supplementary robustness study contains 20 parameter directories and exactly 100 PNG graphs;
- each supplementary parameter directory contains five graphs, `descripcion.txt` and `diagnostico_variantes.csv`;
- no file exceeds GitHub's 25 MiB browser-upload threshold;
- no PDF, non-PNG image or non-graph device schematic is present;
- the converted self-consistent and realistic theory documents retain every chapter, section and subsection heading from their approved LaTeX sources;
- the combined Kitaev theory and implementation document contains 174 numbered equations from the approved documentation.

## Source-code integrity

| File | SHA-256 |
|---|---|
| `models/01-kitaev-chain/code/kitaev_chain_1d.py` | `ec2ebf99ead84354b04bf9fa60e2843070c8cf289b24fb2985b6ff8964871cb4` |
| `models/02-self-consistent-bdg/code/self_consistent_bdg_nanowire.py` | `8b56cfa79d06015efcaa2f58900f3666d2b0f327ec66b3f82004b6801001939f` |
| `models/03-realistic-proximitized-nanowire/code/nanohilosimetricofinal.py` | `631cd5ba0cc6d8b73a5ac4b4275dfbe9d4900616c480f8e9fa59210bc8948825` |

## Figure inventory

| Model | Result graphs | Status |
|---|---:|---|
| Kitaev chain | 10 | complete, ordered and embedded once each |
| Self-consistent BdG nanowire | 13 | complete, ordered and embedded once each |
| Realistic proximitized nanowire | 14 | complete, ordered and embedded once each |
| Supplementary parameter study | 100 | complete across 20 parameter directories |

Every main result figure is followed by the approved explanation of what is calculated and why, and by the corresponding result and physical interpretation. No graph has been regenerated, edited, recoloured, relabelled or silently replaced.

## Documentation fidelity

The scientific prose, derivations, equation sequence, implementation discussion, result analysis and conclusions follow the approved documents. Changes are limited to the transformations required for direct GitHub rendering:

- page-oriented headings were converted to Markdown headings;
- equations were represented as inline or display math source;
- LaTeX-only document machinery was removed;
- graph references were converted to relative Markdown image links;
- page numbers, running headers and PDF pagination were not reproduced;
- unambiguous PDF text-extraction artefacts in the Kitaev material were restored to their intended mathematical notation.

No scientific section was deliberately shortened or replaced by a summary. The code remains in dedicated `.py` files, and the long implementation discussions remain in the combined theory documents in the same approved order.

## Final eight-criterion review

| Criterion | Score | Final assessment |
|---|---:|---|
| 1. Everything performed by the project is explained | 9.8/10 | The three models, algorithms, observables and model progression are documented in full. |
| 2. Explanations are developed, clear and understandable | 9.8/10 | Long-form theory and implementation material is preserved, with navigable Markdown headings. |
| 3. No structural, conversion or repository errors remain | 9.8/10 | Automated rendering, path, inventory and integrity checks report no remaining issues. |
| 4. Equations are correctly represented | 9.8/10 | Equations are native GitHub math, delimiters are balanced and numbered derivations are retained. |
| 5. Intermediate derivation steps are retained | 9.8/10 | The approved equation sequence and explanatory transitions remain in the Markdown documents. |
| 6. Each physical system is presented properly | 9.9/10 | Every model begins with its scope, assumptions, basis, Hamiltonian and relationship to the preceding model. |
| 7. Theory, figures, results and conclusions are fully developed | 9.9/10 | Every principal graph is present in order with calculation-purpose and result-interpretation text. |
| 8. Overall level of detail and repository readiness | 9.8/10 | The package is directly uploadable and requires no scientific rewriting by Copilot. |

**Final mean score: 9.83/10.** All eight criteria exceed 9.5 in the same final pass.

## Final status

The package is ready to upload to the root of the existing GitHub repository. Copilot should be used only to verify that the uploaded paths were preserved and to report any upload-induced flattening; it should not rewrite the Markdown, equations, source code, captions, analyses or conclusions.
