# LAURA v1.9: summary of changes since v1.8

v1.9 builds on v1.8 (Mirror lab, Bond rotation with Newman projection and ethane energy, Sugar lab, Chair lab, Electron lab). It improves the 2D view, selection, image viewing and the lab dialogs. The version number is for this summary only and does not appear anywhere in the app.

Run `Run LAURA.bat` (port 8771) and open the viewer. The labs are now built into `laura-viewer.html`, so the page also works when opened directly from a folder.

## 2D view
- **Fixed:** the 2D view now draws the actual molecule structure. Before, it was only a flat projection of the 3D model.
- **Structure diagram (new default):** rings are drawn as regular polygons and fused rings share an edge. Chains zigzag at 120 degrees, nitriles and CO2 are drawn straight, and substituents point outward. Double, triple and aromatic bonds are drawn (the second line sits inside rings). Heteroatoms are labelled and carbons are left unlabeled, as in a skeletal formula. The Hydrogens option writes H counts on labels (for example OH or NH2), and Element labels also labels carbons.
- **Flat projection (kept):** the original view is still available from a selector, and it is the automatic fallback for bridged or spiro ring systems and for molecules over 80 heavy atoms.
- No stereochemistry is drawn, and the caption says so.
- The image is drawn at twice the resolution, so zooming and saved PNGs stay sharp.

## Viewing and navigation
- **Zoom on images:** click any picture (the 2D view, Ramachandran plot, contact map and so on) to enlarge it, then zoom with the mouse wheel, pinch, the + and - buttons or double-click. Drag to pan, use Fit to reset, and press Esc or click outside to close. The keys +, - and 0 also work.
- **Centered view for "Selected only":** turning it on centers and zooms the view on the selected atoms. The yellow selection halo is hidden in this view. Clearing the selection, or deselecting the last atom, turns the mode off so the whole molecule is visible again.
- **Default rotation speed:** Rotate now starts slower, and the speed slider's default is lower.
- **Selection by chain:** in tube view, clicking a chain selects all of its atoms. Clicking it again deselects them. Other styles still select single atoms.

## Labs
- **Icons:** the lab buttons (Mirror lab, Bond rotation, Sugar lab, Chair lab, Electron lab) are now icons that match the other toolbar buttons. They sit right after Fullscreen, before Info and Help.
- **About section:** every lab has a short "About this lab" section that explains what it does, how to read it and its limits.
- **Sticky close button:** the x stays visible when you scroll down inside a lab. Labs close with the x only, and the old "Close" text is gone.
- **Close on outside click:** clicking outside a lab closes it, matching the other dialogs such as Ramachandran and Info.
- **Bond rotation error:** when the molecule has too few or too many atoms (outside 3 to 80), it shows the same alert as the Mirror lab.

## Help and interface cleanup
- **Help:** it now describes the labs, the structure diagram and projection views, enlarged-picture zoom, "Selected only", the slower rotation, and chain selection in tube view.
- **Teaching button removed:** the Teaching button and its Help entry are gone, and the section is now named "Session & export".
- **Theme button:** the white border around the Dark Mode and Light Mode button is gone.
- **Version text:** version numbers no longer appear in the page title or in the lab dialogs.

## Project and testing
- **Layout:** the project is organised into `src/core/` (logic), `src/labs/` (lab dialogs), `tests/` (tests, fixtures and browser tests), `tools/` and `docs/`.
- **New files:** `src/core/layout2d.mjs` is the 2D diagram engine. `tools/build_viewer.mjs` rebuilds the viewer after you edit a lab, by running `node tools/build_viewer.mjs` from the project root.
- **Tests:** 101 Node tests and 53 browser tests pass. The browser suite uses a Three.js stand-in, so 3D drawing and picking are not exercised there.

## Limits
- The structure diagram has no stereochemistry, and bridged or spiro ring systems fall back to the projection.
- Some built-in molecules have no hydrogens in their data, so the 2D diagram cannot show their H labels (for example phenol's OH).
- Drag-to-select is not implemented.
