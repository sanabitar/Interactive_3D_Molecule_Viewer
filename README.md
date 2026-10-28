# LAURA tests

| Layer | Command | Needs |
|---|---|---|
| Unit + regression (scientific core, parsers, session, shortcuts, framing) | `npm test` | Node 22 |
| Real-browser regression (Playwright for Python, Chromium) | `python3 tests/browser/run_browser_tests.py` | `pip install playwright && playwright install chromium` |
| Blender-script logic (mock `bpy`, NOT Blender) | `python3 tests/blender_mock_check.py` | Playwright |
| Import timing / limits | `PERF=1,100,500,2000 python3 tests/perf_import.py` | Playwright |
| Fixtures | `node tools/make_fixtures.mjs` (synthetic); `tools/fetch_real_fixtures.sh` (real RCSB entries, needs network) | |

**Shim mode (default):** the browser tests replace Three.js with `tests/browser/three_shim.js` when no CDN is reachable. DOM/CSS/JS/logic are real; 3D rendering, picking and GPU memory are **not** exercised.
Run with `--real-cdn` on a networked machine for the real engine. **Still required externally:** run `import_molecule.py` inside Blender; test with real RCSB fixtures.

## v1.9 layout
```
laura-viewer.html    viewer (project root: viewer_parity reads it)
src/core/            pure logic (parsers, analysis, labs' cores, layout2d.mjs = 2D structure diagrams)
src/labs/            lab UIs: mirror, bond rotation (+ energy panel), sugar, chair, electron
tests/               *.test.mjs (npm test), fixtures/, browser/ (Playwright), blender_mock_check.py, perf_import.py
tools/               make_fixtures.mjs, fetch_real_fixtures.sh
docs/                VERSION-1.4 to 1.7, MERGE_NOTES, readme updates, archive/
```
The labs are ES modules, so open the page through a server (`Run LAURA.bat`). For the browser suite: run `python3 -m http.server 8799` in the project root, then `LAURA_URL=http://127.0.0.1:8799/laura-viewer.html python3 tests/browser/run_browser_tests.py`.

2D view: a generated structure diagram (rings as regular polygons, fused rings, zigzag chains, substituents outward; no stereochemistry; bridged/spiro systems and >80 heavy atoms fall back to the flat projection). The page's inlined copy of `layout2d.mjs` is checked by `tests/layout2d.test.mjs`.
