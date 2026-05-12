# GitHub Copilot Organization Prompt

Use this prompt after uploading all files from the flat ZIP to GitHub.

---

I uploaded a flat set of files for a polished scientific Python repository about Majorana zero modes in 1D superconducting systems. Please reorganize the repository into a clean professional structure without changing the scientific meaning of the text or the numerical logic of the scripts.

Required final structure:

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── src/
│   ├── kitaev_chain_1d.py
│   ├── self_consistent_bdg_nanowire.py
│   └── proximitized_bdg_nanowire.py
├── docs/
│   ├── theory_and_history.md
│   ├── results_and_discussion.md
│   └── implementation_notes.md
└── assets/
    ├── ASSETS_INDEX.md
    └── figures/
        ├── kitaev/
        ├── self_consistent/
        ├── realistic_nanowire/
        └── quasi_majoranas/
```

Move files as follows:

- Move all `*.py` files into `src/`.
- Move `theory_and_history.md`, `results_and_discussion.md`, and `implementation_notes.md` into `docs/`.
- Move `ASSETS_INDEX.md` into `assets/ASSETS_INDEX.md`.
- Move figures beginning with `fig_kitaev_` into `assets/figures/kitaev/`.
- Move figures beginning with `fig_selfconsistent_` into `assets/figures/self_consistent/`.
- Move figures beginning with `fig_realistic_` into `assets/figures/realistic_nanowire/`.
- Move figures beginning with `fig_quasi_` into `assets/figures/quasi_majoranas/`.

Then update every Markdown image link so the images render correctly after the move. For example, `fig_kitaev_bulk_bands.png` should become `assets/figures/kitaev/fig_kitaev_bulk_bands.png` from the README, and `../assets/figures/kitaev/fig_kitaev_bulk_bands.png` from files inside `docs/`.

Polish the README so it looks excellent on GitHub:

- Keep the title and project summary.
- Add a clean table of contents.
- Keep the visual highlights, but make sure image links work.
- Add a short “Scientific motivation” section.
- Add a short “How to run” section.
- Add a short “Limitations” section so the project does not overclaim experimental Majorana detection.
- Add links to the theory, results, and implementation docs.

Do not delete any selected figures. Do not replace the equations with plain text. Preserve the GitHub Markdown math blocks. Do not change the numerical methods unless the change is purely structural, such as fixing imports or guarding expensive execution behind `if __name__ == "__main__"`.

After reorganizing, verify:

1. all image links render,
2. all Python files compile with `python -m py_compile src/*.py`,
3. `README.md` links to the docs correctly,
4. no LaTeX build artifacts or macOS metadata files are present,
5. the repository looks professional enough to send to recruiters or technical interviewers.
