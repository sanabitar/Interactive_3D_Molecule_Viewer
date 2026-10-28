# LAURA v1.6 — Sugar and Fischer lab
Run Run LAURA.bat (port 8769), then click Sugar lab.

Choose among 16 open-chain aldohexoses. Invert C2–C5 individually, mirror all centers, reset, or export the current Fischer projection as SVG. The lab identifies configurations and relationships (identical, epimers, other diastereomers, enantiomers) and highlights differences from the reference.

Examples: D-glucose inverted at C2 gives D-mannose; at C4 gives D-galactose. Mirroring every center gives L-glucose. The OH position at C5 determines the D/L series in the displayed canonical projection.

Scope: predefined open-chain aldohexose configurations. Not automatic conversion from arbitrary imported 3D coordinates. No cyclic forms, alpha/beta anomers, chair models, or solution populations yet. D/L is not a prediction of optical rotation sign.
References:
https://openstax.org/books/organic-chemistry/pages/25-4-configurations-of-the-aldoses
https://openstax.org/books/organic-chemistry/pages/25-3-d-l-sugars

Validation: 18 tests passed across sugar-core, energy-core, conformation-core and mirror-core suites. Browser checks confirmed C2 epimerization and the D/L glucose mirror pair. Main-viewer keyboard shortcuts are suppressed while a lab dialog is open.

Next stages: cyclic sugar/chair teaching views; electrons and lone pairs; stereochemistry-aware projections from imported 3D structures.
