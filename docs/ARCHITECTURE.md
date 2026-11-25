# LAURA architecture

LAURA is plain files: `index.html` (markup only), `css/` and `js/`. There is no bundler and no build step for the app itself, so you can open `index.html` by double-click or serve the folder (`Run LAURA.bat`). The only generated code is in `js/generated/` (see below).

## Why classic scripts (not ES modules)
Browsers refuse to load ES modules from `file://`, and LAURA is meant to work when opened from disk. So every file is a classic `<script src>`. They share **one global scope**: a top-level `const`, `let` or `function` in one file is visible in every file loaded after it. That is how the code always worked (it used to be one big inline script); the split only changed where each part lives.

## Rules that keep it working
1. **Load order is part of the code.** A file may use anything defined in an *earlier* file at load time; anything may be used by functions *called later* (after start-up). The order is the order of the `<script>` tags in `index.html`, listed below. `tests/structure.test.mjs` checks the key constraints.
2. **Files only define things. `js/main.js` starts them.** The render loop (`renderLoop`) and the saved-imports restore (`restoreSavedImports`) are defined in `viewer/picking.js` and `app/imports.js` and started from `main.js`, the last app script. Do not start a loop, observer or timer at the top level of any other file.
3. **`js/generated/*` are built, not edited.** They come from `src/` (`hydrogens.mjs`, `layout2d.mjs`, and the labs with the cores they import) via `npm run build`. `npm test` fails if one is stale. Edit the file in `src/`, then rebuild.
4. **CSS order matters.** Several selectors (for example `.it` and `aside.r`) are defined more than once and later rules deliberately override earlier ones. The `<link>` order in `index.html` is the original rule order; do not reorder.
5. **No inline code in `index.html`** (no `<script>` body, `<style>` or `onclick=`); a test enforces it.

## Stylesheets (load order)
| # | File | Size | What it styles |
|---|---|---|---|
| 1 | `css/theme.css` | 0.8 KB | Colour variables for the light and dark themes (dark follows the OS unless the Dark/Light Mode button was used). |
| 2 | `css/base.css` | 1.9 KB | Reset, page grid, header/footer, buttons, form controls, headings and list items. |
| 3 | `css/viewer.css` | 1.1 KB | The 3D viewport: toolbar, loading/error overlay, atom tooltip, toast message and the presentation card. |
| 4 | `css/layout.css` | 2.1 KB | Panel layout: library | viewer | details grid, resize handles, details-panel toggle, collapsible sections, tool button grid, narrow-screen rules. |
| 5 | `css/controls.css` | 1.7 KB | Icon buttons, text-style buttons, Dark/Light Mode buttons and range sliders. |
| 6 | `css/dialogs.css` | 3.5 KB | The INFO card and the Help dialog. |
| 7 | `css/panels.css` | 5.4 KB | Section headers, the library panel (current molecule, search, filters, list), file inputs, analysis tool blocks, and the small-screen / dark-mode / reduced-motion overrides. |
| 8 | `css/overlays.css` | 2.4 KB | Presentation mode, the scrolling details panel, atom labels (element symbols and name labels) and the enlarged-diagram lightbox. |
| 9 | `css/import-and-plots.css` | 1.6 KB | The import section and dialog graphics (plots and images keep their aspect ratio). |

## Scripts (load order)
| # | File | Size | What it contains |
|---|---|---|---|
| 1 | `js/foundation/constants.js` | 0.9 KB | Shared constants and helpers: $ (querySelector), element colours and radii (CPK, VDW, COV), toast(), esc(). Loaded first. |
| 2 | `js/generated/hydrogens.js` | 6.4 KB | GENERATED from src/core/hydrogens.mjs. Do not edit. |
| 3 | `js/foundation/library.js` | 22.2 KB | The built-in molecule library (43 entries in LIB) and the `imported` list of user-added structures. |
| 4 | `js/foundation/core.js` | 40.0 KB | Core: parsers, bond perception, analysis, models and sessions, plus parseFile() and computeBonds(). |
| 5 | `js/viewer/scene.js` | 2.8 KB | three.js scene: renderer, camera, orbit controls, lighting, and the shared viewer state (mol, sel, vis, ...). |
| 6 | `js/viewer/display-options.js` | 1.8 KB | Toolbar toggles: the view button (ball and stick / chains) and the hydrogens switch (setH). |
| 7 | `js/viewer/molecule-render.js` | 3.3 KB | build(): turns the current molecule into spheres, bonds and (in Chains view) tubes, applying the display options. |
| 8 | `js/viewer/labels.js` | 6.8 KB | Label layer: element symbols on atom spheres (debounced so they do not flicker), atom-name labels, chain and contact tags, |
| 9 | `js/generated/layout2d.js` | 9.3 KB | GENERATED from src/core/layout2d.mjs. Do not edit. |
| 10 | `js/tools/view2d.js` | 8.2 KB | The 2D view dialog: generated structure diagram (uses Layout2D) or a flat projection of the 3D coordinates. |
| 11 | `js/viewer/chains-and-selection.js` | 2.8 KB | Chains view (tubes along the CA/P trace) and the selection highlight (drawSel0). |
| 12 | `js/viewer/picking.js` | 1.8 KB | Mouse picking and the hover tooltip, and the render loop. renderLoop() is only defined here; main.js starts it. |
| 13 | `js/app/loading.js` | 3.6 KB | Loading a molecule: the loading/error overlay, load(), the details panel (info) and the bond-order summary. |
| 14 | `js/app/library-list.js` | 3.0 KB | The molecule list and menu (search, sort, categories) and opening a local file (addFile). |
| 15 | `js/app/ligand-contacts.js` | 2.7 KB | Ligand contacts: list ligands, find the residues in contact, show, copy and export them. |
| 16 | `js/app/rcsb.js` | 3.0 KB | Fetching from RCSB: metadata, text search and downloading a structure by ID. |
| 17 | `js/ui/modal.js` | 1.1 KB | Modal dialog helper (modal(), focus trap) and need(), which asks for a structure to be loaded first. |
| 18 | `js/tools/tool-docs.js` | 4.4 KB | Title, icon and short description of each analysis tool (TD), and tdoc(), the header card its dialog shows. |
| 19 | `js/tools/analysis-tools.js` | 7.8 KB | Analysis dialogs: Ramachandran plot, contact map and sequence alignment (and the Labels checkbox). |
| 20 | `js/app/session.js` | 2.7 KB | Saving and restoring a session (structure, selection, camera, display settings). |
| 21 | `js/tools/teaching.js` | 1.0 KB | The guided "measure water" activity. |
| 22 | `js/tools/compare-and-rmsd.js` | 11.1 KB | Comparison overlays: Compare molecules (table and side-by-side) and Superposition & RMSD. |
| 23 | `js/tools/help.js` | 9.9 KB | Help and Analysis-guide text, and the tooltips of the toolbar and panel buttons. |
| 24 | `js/ui/panels.js` | 1.9 KB | Resizable library/details panels (widths are remembered) and the details-panel toggle (the panel starts closed). |
| 25 | `js/ui/header.js` | 2.3 KB | Header controls: the INFO card, fullscreen, and the Dark/Light Mode button. |
| 26 | `js/ui/lightbox.js` | 3.1 KB | Closing dialogs on outside click, the enlarge-a-diagram lightbox (zoom and pan), and the 2D view button. |
| 27 | `js/app/imports.js` | 3.3 KB | Imported structures: rename, save in this browser (IndexedDB), remove. restoreSavedImports() is started from main.js. |
| 28 | `js/app/export.js` | 5.3 KB | Downloads: XYZ and the Blender package (including the generated Blender import script). |
| 29 | `js/main.js` | 4.9 KB | Wiring and START-UP: theme, keyboard shortcuts and the remaining buttons. Loaded last; the only file that starts work at load. |
| 30 | `js/app/mirror-source.js` | 0.3 KB | Hands the current molecule to the Mirror lab (generated/labs.js reads window.lauraMirrorSource). |
| 31 | `js/generated/labs.js` | 61.8 KB | GENERATED from src/labs/*.mjs and the src/core modules they import. Do not edit. |

`foundation/core.js` is a hand-maintained copy of the pure modules in `src/core`; `tests/viewer_parity.test.mjs` fails if its behaviour drifts from `src/core`.

## Adding or moving a file
1. Create the file with a one-line header (`// js/<folder>/<name>.js: what it is for`; CSS: `/* css/<name>.css: ... */`).
2. Add its tag to `index.html` **at the right place in the order** (after everything it uses at load time, before `main.js`).
3. Run `npm test`: the structure test fails if a file is unreferenced, referenced twice, missing, has no header, or is out of order.

## Where did the old sections go?
| Old section of the single script | Now |
|---|---|
| helpers at the top (`$`, CPK, VDW, `toast`) | `foundation/constants.js` (`esc` moved here too) |
| hydrogens (inlined block) | `generated/hydrogens.js` |
| library | `foundation/library.js` |
| parsers & validation (the `Core` object) | `foundation/core.js` |
| measurements, three.js scene | `viewer/scene.js`, `viewer/display-options.js`, `viewer/molecule-render.js` (`build`) |
| label layer | `viewer/labels.js` |
| layout2d (inlined block) | `generated/layout2d.js` |
| 2D view (also held `tubes` and `drawSel0`) | `tools/view2d.js`, `viewer/chains-and-selection.js` |
| picking (and the render loop) | `viewer/picking.js` (loop started in `main.js`) |
| loading (also the molecule list) | `app/loading.js`, `app/library-list.js` |
| ligand contacts (also `esc`, RCSB fetch) | `app/ligand-contacts.js`, `app/rcsb.js` |
| analysis tools, session, teaching | `ui/modal.js`, `tools/tool-docs.js`, `tools/analysis-tools.js`, `app/session.js`, `tools/teaching.js`, `tools/compare-and-rmsd.js`, `tools/help.js` |
| resizable panels, header, outside-click/lightbox | `ui/panels.js`, `ui/header.js`, `ui/lightbox.js` |
| imports, export, wiring | `app/imports.js`, `app/export.js`, `main.js` |
| the small mirror script, the labs script | `app/mirror-source.js`, `generated/labs.js` |

## Not done (candidates for later)
- `foundation/core.js` duplicates `src/core` by hand; it could be generated like the other files once the few differences are reconciled.
- Some CSS rules are defined twice; they could be merged now that each file has one topic (this needs a visual check, which is why it was not part of the split).
- The scripts share globals. Moving to ES modules would give real imports but would end `file://` support.
