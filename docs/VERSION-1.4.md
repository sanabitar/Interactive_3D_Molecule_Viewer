# LAURA v1.4
Launch Run LAURA.bat, then open Bond rotation from the top bar.
Start with Ethane (staggered). Move the slider or use +/-60 degrees to explore conformations alongside the synchronized Newman projection. Look down bond aligns the camera; Reverse viewing direction swaps the front and rear fragments and resets the rotation.

Supported: 3–80 atoms, non-ring single bonds between heavy atoms with substituents at both ends. Amide-like C–N bonds are excluded. Bond eligibility depends on inferred bond orders. No missing hydrogens are generated.

The lab preserves the loaded structure. It performs rigid geometric rotation only, without collision prevention, energy minimization, or an energy model. Newman labels identify adjacent atoms, not complete substituent groups. Axial substituents have no defined projected angle and are flagged rather than drawn in an arbitrary direction.

Validation: node --test conformation-core.test.mjs mirror-core.test.mjs
10 tests passed. Browser checks: ethane staggered/eclipsed, reverse direction, aligned camera, and cyclohexane exclusion.

Next planned stage: conformational energy visualization. Fischer/chair projections, sugar teaching models, and electron/lone-pair displays remain future stages.
