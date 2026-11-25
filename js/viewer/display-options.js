// js/viewer/display-options.js: Toolbar toggles: the view button (ball and stick / chains) and the hydrogens switch (setH).
const ICON_BS='<svg viewBox="0 0 24 24" aria-hidden="true"><circle cx="7" cy="17" r="3"/><circle cx="17" cy="7" r="3"/><path d="M9.1 14.9l5.8-5.8"/></svg>',ICON_CH='<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M3 15c3-9 6 5 9-2s6 5 9-4"/><circle cx="3" cy="15" r="1.2"/><circle cx="21" cy="9" r="1.2"/></svg>';
// Two views: ball and stick, and chains (a smooth trace through each polymer chain's CA/P atoms). One small button toggles between them.
function syncStyBtn(){const tb=$('#sty').value==='tb',b=$('#stb');b.classList.toggle('on',tb);b.setAttribute('aria-pressed',tb);b.title=tb?'View: Chains. Click for Ball and stick.':'View: Ball and stick. Click for Chains.';b.setAttribute('aria-label',tb?'View: Chains (click for Ball and stick)':'View: Ball and stick (click for Chains)');if(b._s!==tb){b.innerHTML=tb?ICON_CH:ICON_BS;b._s=tb}}
function chainTraces(){const by={};vis.forEach(i=>{const a=mol.atoms[i];if(!a.hetero&&(a.name==='CA'||a.name==='P'))(by[a.chain]=by[a.chain]||[]).push(i)});return by}
function syncHBtn(){const on=$('#hh').checked,b=$('#hb');b.classList.toggle('on',on);b.setAttribute('aria-pressed',on);b.title=(on?'Hide':'Show')+' hydrogens (H)';b.setAttribute('aria-label',(on?'Hide':'Show')+' hydrogens')}
// One switch for hydrogens: the toolbar button, the H key, the Hydrogens checkbox and the checkbox in the 2D dialog all end up here.
let redraw2D=null;
function setH(on){$('#hh').checked=!!on;syncHBtn();if(mol){build();if(on&&!mol.atoms.some(a=>a.el==='H'))toast(mol.entry.noH||'This structure contains no hydrogen atoms, so there is nothing to show. Imported files only have the hydrogens they were saved with.',7000)}
 const d=$('#d2h');if(d)d.checked=!!on;if(redraw2D&&$('#d2c'))redraw2D()}
