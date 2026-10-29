We started working from v1.2, so I can’t confirm what changed between v1.0 and v1.2.
Version	Updates we made
v1.2	Updated the application title and credits. Added Run LAURA.bat for launching the app.
v1.3	Added Mirror lab: generate mirror images, rotate both views independently, check rigid superimposability, and overlay aligned structures.
v1.4	Added Bond rotation with angle controls and synchronized Newman projections, including staggered/eclipsed indicators.
v1.5	Added an interactive ethane conformational-energy curve, synchronized with rotation and clickable to change the angle.
v1.6	Added Sugar lab with 16 open-chain aldohexose configurations, Fischer projections, stereocenter inversion, mirror comparisons, epimer identification, and SVG export.
v1.7	Added Chair lab for glucose, mannose, and galactose: α/β controls, Haworth projections, rotatable 3D chairs, and axial/equatorial tracking during chair flips.
v1.8	Added Electron lab for water, ammonia, methane, and CO₂: bonding pairs, lone pairs, display toggles, electron counts, and molecular-geometry explanations.


Each new version preserves the previous version in a separate folder and includes its own launcher. The latest feature test suite has 27 passing tests, supplemented by browser checks.

Current limits: energy visualization covers an approximate ethane teaching model; sugar, chair, and electron labs use defined examples. General conversion of imported 3D structures into Fischer/chair projections and automatic lone-pair assignment remain unfinished.