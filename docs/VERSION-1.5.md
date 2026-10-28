# LAURA v1.5 — Ethane energy teaching model
Run Run LAURA.bat (port 8768). Select Ethane (staggered), open Bond rotation, and scroll down to the energy graph. Rotate with the slider, +/-60 degree controls, or click the curve. The marker and numeric energy stay synchronized with the molecule and Newman projection.

Model: E = 6[1 + cos(3 phi)] kJ/mol, with phi the projected front/rear C-H dihedral. The model uses the loaded geometry to establish its phase. The horizontal axis is rotation from that geometry, not absolute dihedral. Energies are relative to the staggered minimum. Barrier: approximately 12 kJ/mol.
Reference: https://openstax.org/books/organic-chemistry/pages/3-6-conformations-of-ethane

This approximate threefold teaching potential is only enabled for complete ethane connectivity with near-threefold projected substituent spacing. It is not a force-field or quantum calculation and does not include relaxation, bond stretching or angle bending. Unsupported molecules retain the rotation and projection features, but do not display numerical energy values.

Validation: node --test energy-core.test.mjs conformation-core.test.mjs mirror-core.test.mjs
14 automated tests passed. Browser validation covered staggered/eclipsed energies, clicking the curve to set rotation, and unsupported alanine messaging.

Remaining roadmap: sugar structures and configurations; Fischer and chair projections; electrons and lone pairs; additional validated energy models.
