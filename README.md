# Training for Marita

Self-contained agent skills and training material, generated from scientific
books and papers with the [book-to-skill](https://github.com/zuwasi) converter
workflow for Amp / Claude Code style agents.

## Contents

| Folder | What it is |
|---|---|
| [`rozenson-marita/`](rozenson-marita/) | Knowledge base from **"The Non-Newtonian Approach to the Muco-Ciliary Transport"** by Marita Rozenson (PhD thesis, Weizmann Institute of Science, 1985). Oldroyd viscoelastic mucus modeling, ciliary penetration, Fourier solution of pulsatile two-layer flows. |
| [`schmied-mathematica-qm/`](schmied-mathematica-qm/) | Knowledge base from **"Using Mathematica for Quantum Mechanics: A Student's Manual"** by Roman Schmied (arXiv:1403.7050). Wolfram Language workflows for quantum mechanics: spin, coupled qubits, real-space propagation, Gross-Pitaevskii. |
| [`GroverSearchTraining/`](GroverSearchTraining/) | Lean 4 project skeleton for Grover-search training (Lakefile-based; `.lake/` build cache is not tracked). |

Each skill folder follows the same layout:

```
SKILL.md          core frameworks + chapter/topic index (load this first)
chapters/         per-chapter deep dives (loaded on demand)
glossary.md       key terms and definitions
patterns.md       solution techniques and methods
cheatsheet.md     quick-reference tables
```

## Usage

Copy a skill folder into your agent's skills directory, e.g. for Amp:

```
~/.config/agents/skills/rozenson-marita/
```

Then ask your agent about a topic (for example "penetration" or "Oldroyd
equation") and it will load the matching chapter file, or ask for a specific
chapter such as `ch04a`.

## License

MIT. See [LICENSE](LICENSE).

Note: source documents (such as `rozenson_marita.pdf`) are copyrighted by
their authors and are intentionally not distributed in this repository; the
skills contain only synthesized, structured knowledge bases.
