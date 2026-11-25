// js/app/mirror-source.js: Hands the current molecule to the Mirror lab (generated/labs.js reads window.lauraMirrorSource).
window.lauraMirrorSource=()=>{if(!need())return null;return {atoms:mol.atoms,bonds:mol.bonds,orders:mol.orders,name:mol.entry.name};};
