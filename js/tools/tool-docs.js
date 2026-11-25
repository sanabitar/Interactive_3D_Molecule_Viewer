// js/tools/tool-docs.js: Title, icon and short description of each analysis tool (TD), and tdoc(), the header card its dialog shows.
const TD={
ra:{t:'Ramachandran Plot',s:'Protein backbone geometry',i:'<circle cx="7" cy="16" r="1.5"/><circle cx="11" cy="9" r="1.5"/><circle cx="16" cy="14" r="1.5"/><circle cx="18" cy="6" r="1.5"/><path d="M3 3v18h18"/>',
p:'Examine backbone φ/ψ angles to understand protein conformation and identify unusual residues.',
o:'Each point represents an amino acid residue, plotted by its phi (x) and psi (y) dihedral angles. Shaded areas are schematic reference regions marking commonly favoured backbone conformations (α-helix, β-strand/PPII and left-handed helix). Glycine (triangles) and proline (squares) are marked separately because their allowed areas differ from other residues.',
h:'Pick “Right-handed α-helix (Ala12, idealised)” in the library and open Ramachandran: all 10 residues fall in the α-helix region. Then pick “β-strand (Ala12, idealised)”: all 10 move to the β/PPII region. Hover a point to see which residue it is, or click the plot to enlarge it.',
n:'Unusual angles are not necessarily structural errors. The shaded regions are schematic, not validation contours.'},
cmb:{t:'Contact Map',s:'Residue proximity within a chain',i:'<rect x="4" y="4" width="16" height="16" rx="1"/><path d="M4 9h16M4 14h16M9 4v16M14 4v16"/>',
p:'Reveal which residues lie close together in 3D to understand fold and packing.',
o:'A square matrix of the residues in one chain, in file order. A filled cell means the CA atoms of that residue pair are within the cutoff distance (default 8 Å).',
h:'Pick the α-helix model and open Contact map: each residue touches neighbours up to 4 positions away (a band 4 cells wide either side of the diagonal, 38 contacting pairs). Pick the β-strand model: the band narrows to 2 cells (21 pairs) because the chain is extended. Click the map to enlarge it.',
n:'CA–CA proximity is geometric and does not prove a physical interaction.'},
sq:{t:'Sequence Alignment',s:'Pairwise global comparison',i:'<path d="M4 6h7M13 6h7M4 12h4M10 12h10M4 18h10M16 18h4"/>',
p:'Compare two protein sequences to see how much of them matches.',
o:'Needleman–Wunsch global alignment (match +2, mismatch −1, gap −2, no substitution matrix). The result reports score, identity, coverage and gaps.',
h:'Paste ACGT and AGT and press Align: 3 of 4 columns are identical (75%) with one gap. Or load “Antiparallel helix pair (chains A and B)” (two chains) and press Align: 10 of 12 columns are identical.',
n:'Sequence identity alone does not establish shared function or structure.'},
rm:{t:'Superposition & RMSD',s:'Structural similarity in 3D',i:'<circle cx="9" cy="12" r="5"/><circle cx="15" cy="12" r="5"/>',
p:'Quantify how closely two structures match in three dimensions.',
o:'Structure B is moved onto the current structure A by the least-squares rotation and translation, and RMSD (Å) is reported over the paired atoms. In the overlay, cyan is A, magenta is B and grey lines join paired atoms.',
h:'Load “Antiparallel helix pair (chains A and B)”, keep chain A as the reference and choose chain B to move: the RMSD is about 0 Å over 12 CA pairs. Then load the α-helix model, choose the β-strand model as B: the RMSD is 6.7 Å over 12 pairs, a clear structural difference.',
n:'RMSD depends on which atoms are paired and is inflated by flexible regions; it is meaningful only for homologous proteins or the same molecule.'},
cmp:{t:'Compare Molecules',s:'Side-by-side properties',i:'<rect x="3" y="4" width="8" height="16" rx="1"/><rect x="13" y="4" width="8" height="16" rx="1"/>',
p:'Compare basic physical properties of two molecules side by side.',
o:'A table of formula, atom counts, molecular weight, bonds, radius of gyration and size for molecules A and B, with the B − A difference.',
h:'Compare Water with Methane: H2O weighs 18.015 g/mol and CH4 16.043 g/mol; the table also lists atom counts, bonds, radius of gyration and size, and you can export it as CSV.'
,n:'These are descriptive properties, not a chemical-similarity score.'}};
const tdoc=id=>{const d=TD[id];return `<details class="td"><summary title="Show or hide the description"><svg viewBox="0 0 24 24" aria-hidden="true">${d.i}</svg><span><span class="tt">${esc(d.t)}</span><span class="ts">${esc(d.s)}</span></span></summary><h4>Purpose</h4><p>${esc(d.p)}</p><h4>Overview</h4><p>${esc(d.o)}</p><h4>Try it</h4><p>${esc(d.h)}</p><div class="tn">ⓘ Scientific note: ${esc(d.n)}</div></details>`};
