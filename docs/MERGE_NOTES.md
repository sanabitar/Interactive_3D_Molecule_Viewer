# Merge notes: uploaded line of work + audit line of work

**Sources.** (A) your uploaded project (`models.mjs`, Ramachandran regions/summary, flexible `peptide()`, `PM`/`E` library entries, `tdoc` tool descriptions, `viewer_parity` test). (B) the audit line (`structure`/`bonds` hardening, `session`, `ui`, `camera`, page fixes D1–D16, browser/Blender-mock/perf tests).

| File | Resolution |
|---|---|
| `src/core/analysis.mjs`, `build.mjs`, `models.mjs` | **A's version**; B's one-line stack-overflow fix (`maxDeviation` via `reduce`) re-applied on top of A's `analysis.mjs`. |
| `src/core/structure.mjs`, `bonds.mjs`, `properties.mjs` | **B's version** (strict supersets of A's: model/truncation/duplicate/element warnings, residue-template bond orders, no argument spreading). |
| `src/core/contacts|geometry|superpose.mjs` | identical on both sides. |
| `src/core/session|ui|camera.mjs` | only in B, kept. |
| `laura-viewer.html` | **B's page as the base** (all audit fixes), A's features ported in by bracket-matched extraction: `TD`/`tdoc` + 979 chars of `.td` CSS, the new Ramachandran handler, `SRCM`/`PM` helpers + 6 `PM` entries, 7 new `E` molecules, Ramachandran help text. Inlined `Core` regenerated from the merged modules (exposes `models`, `ramaSummary`, `RAMA_REGIONS`, `ramaRegion`, `transformAtoms`). |
| tests | A's `models`, `rama`, `viewer_parity` added; the rest were identical. `package.json` runs the union. |

**Deliberately NOT taken from A's page:** A's older copies of CSS rules that B had fixed (they would have re-introduced the 390 px overflow), and the pre-audit versions of code B corrected (Help handler override, missing dropdown handler, unescaped `innerHTML`, `Math.max(...spread)`).

**New work found by the merge:** ring perception (planar 5/6-rings of C/N with aromatic-range bond lengths become aromatic) so A's pyridine/imidazole entries draw consistently.

`laura-viewer.html` must stay at the project root: `tests/viewer_parity.test.mjs` reads it from there.

## Update: 16 viewer changes (2D view, labels, selected-only, rename/save imports, ...)
New core modules: `depict2d.mjs` (best-fit-plane projection), `labels.mjs` (screen-space declutter + label text), `ui.mjs#cleanName`. New tests: `ui2`, `tryit` (pins every number quoted in the tool "Try it" text). Browser suite: 49 tests (`--only 36,38` runs a subset).
Removed on purpose: the disabled "Cartoon & surface" style (not implemented, cannot be verified without real Three.js), the Teaching quiz.
