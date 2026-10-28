import {ethanePhase,ethaneEnergy} from '../core/energy-core.mjs';
export function energyPanel(body,atoms,bonds,orders,onAngle){
 const section=document.createElement('section');section.innerHTML='<h3>Conformational energy · ethane teaching model</h3><p class="envalue" role="status" aria-live="polite"></p><svg class="engraph" viewBox="0 0 720 230" preserveAspectRatio="none" role="img" aria-label="Relative torsional energy versus rotation angle" style="height:230px;cursor:crosshair"></svg><p class="enhelp">Click the curve to set the rotation, or use the angle slider above. Energy is relative to the staggered minimum.</p><p><small>Approximate model: E = 6[1 + cos(3φ)] kJ/mol, with φ measured between a front and rear C–H bond. The curve uses the starting geometry, so 0° rotation means the loaded conformation. This is a teaching potential with a 12 kJ/mol barrier, not a force-field or quantum calculation. It omits bond stretching, angle bending, and relaxation. <a href="https://openstax.org/books/organic-chemistry/pages/3-6-conformations-of-ethane" target="_blank" rel="noopener noreferrer">Reference: OpenStax, Conformations of Ethane</a>.</small></p>';
 body.querySelector('.clstatus').after(section);const graph=section.querySelector('svg'),value=section.querySelector('.envalue');let phase=null;
 const add=(tag,attrs,text)=>{const el=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,v);if(text!==undefined)el.textContent=text;graph.append(el);return el;};
 const x=a=>60+(a+180)/360*630,y=e=>180-e/12*145;
 graph.addEventListener('pointerdown',e=>{if(phase===null)return;const r=graph.getBoundingClientRect(),px=(e.clientX-r.left)/r.width*720;onAngle(Math.round(Math.max(-180,Math.min(180,(px-60)/630*360-180))));});
 return (front,back,angle)=>{
  phase=ethanePhase(atoms,bonds,orders,front,back);graph.style.display=phase===null?'none':'block';section.querySelector('.enhelp').hidden=phase===null;
  if(phase===null){value.textContent='Energy model unavailable for this structure. Choose all-atom Ethane (staggered). Other molecules and substantially distorted geometries need a different validated model; their rotation and Newman views remain available.';return;}
  const e=ethaneEnergy(angle,phase);value.textContent=`Relative torsional energy: ${e.toFixed(2)} kJ/mol (${(e/4.184).toFixed(2)} kcal/mol) · rotation ${angle}°. `+(e<.05?'Near a staggered energy minimum.':e>11.95?'Near an eclipsed energy maximum.':'Between an energy minimum and maximum.');
  graph.replaceChildren();add('text',{x:60,y:18,fill:'#d4e7df','font-size':13},'Relative energy (kJ/mol)');
  for(const energy of [0,4,8,12]){add('line',{x1:60,x2:690,y1:y(energy),y2:y(energy),stroke:'#355347'});add('text',{x:47,y:y(energy)+4,fill:'#d4e7df','text-anchor':'end','font-size':12},energy);}
  for(const angle of [-180,-120,-60,0,60,120,180])add('text',{x:x(angle),y:200,fill:'#d4e7df','text-anchor':'middle','font-size':12},angle+'°');
  add('text',{x:375,y:223,fill:'#d4e7df','text-anchor':'middle','font-size':13},'Rotation from loaded conformation');
  const points=Array.from({length:361},(_,i)=>`${x(i-180)},${y(ethaneEnergy(i-180,phase))}`).join(' ');add('polyline',{points,fill:'none',stroke:'#6de5ee','stroke-width':3});
  add('line',{x1:x(angle),x2:x(angle),y1:30,y2:180,stroke:'#ffbd76','stroke-dasharray':'4 4'});add('circle',{cx:x(angle),cy:y(e),r:6,fill:'#ffbd76',stroke:'#10231e','stroke-width':2});
 };
}
