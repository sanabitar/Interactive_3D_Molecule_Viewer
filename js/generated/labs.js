// GENERATED FILE - do not edit. Source: src/labs/*.mjs and the src/core modules they import
// Rebuild with: node tools/build_viewer.mjs
(function(){try{
const __M={};
__M['geometry']=(()=>{
// Geometry. Inputs {x,y,z} in Angstrom; angles returned in degrees.
const sub=(a,b)=>[b.x-a.x,b.y-a.y,b.z-a.z],dot=(u,v)=>u[0]*v[0]+u[1]*v[1]+u[2]*v[2],cross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],len=u=>Math.hypot(...u);
const distance=(a,b)=>len(sub(a,b));
function angle(a,b,c){const u=sub(b,a),v=sub(b,c),m=len(u)*len(v);if(m===0)throw new RangeError('Angle undefined: coincident atoms.');return Math.acos(Math.max(-1,Math.min(1,dot(u,v)/m)))*180/Math.PI}
// IUPAC sign convention; range (-180, 180]
function dihedral(a,b,c,d){const b1=sub(a,b),b2=sub(b,c),b3=sub(c,d);if(!len(b2)||!len(cross(b1,b2))||!len(cross(b2,b3)))throw new RangeError('Dihedral undefined: collinear or coincident atoms.');
 return Math.atan2(len(b2)*dot(b1,cross(b2,b3)),dot(cross(b1,b2),cross(b2,b3)))*180/Math.PI}

return {distance,angle,dihedral};})();
__M['superpose']=(()=>{
// Least-squares rigid superposition (Horn 1987 quaternion method == Kabsch optimum). Equal weights, rotation+translation only (no reflection).
const {distance}=__M['geometry'];
function eig4(N){const A=N.map(r=>[...r]),V=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
 for(let sw=0;sw<80;sw++){let off=0;for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)off+=A[i][j]**2;if(off<1e-26)break;
  for(let p=0;p<3;p++)for(let q=p+1;q<4;q++){if(Math.abs(A[p][q])<1e-300)continue;const th=(A[q][q]-A[p][p])/(2*A[p][q]),t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),s=t*c;
   for(let k=0;k<4;k++){const x=A[k][p],y=A[k][q];A[k][p]=c*x-s*y;A[k][q]=s*x+c*y}
   for(let k=0;k<4;k++){const x=A[p][k],y=A[q][k];A[p][k]=c*x-s*y;A[q][k]=s*x+c*y}
   for(let k=0;k<4;k++){const x=V[k][p],y=V[k][q];V[k][p]=c*x-s*y;V[k][q]=s*x+c*y}}}
 return{values:[0,1,2,3].map(i=>A[i][i]),vectors:V}}
// Finds the rotation+translation that moves P (moving) onto Q (reference). P[i] pairs with Q[i].
function superpose(P,Q){const n=P.length;if(n!==Q.length)throw new RangeError('Coordinate sets differ in length.');if(n<3)throw new RangeError('At least 3 coordinate pairs are required.');
 const cen=X=>({x:X.reduce((s,a)=>s+a.x,0)/n,y:X.reduce((s,a)=>s+a.y,0)/n,z:X.reduce((s,a)=>s+a.z,0)/n}),cp=cen(P),cq=cen(Q);
 let Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0,G=0;
 for(let i=0;i<n;i++){const a=[P[i].x-cp.x,P[i].y-cp.y,P[i].z-cp.z],b=[Q[i].x-cq.x,Q[i].y-cq.y,Q[i].z-cq.z];
  Sxx+=a[0]*b[0];Sxy+=a[0]*b[1];Sxz+=a[0]*b[2];Syx+=a[1]*b[0];Syy+=a[1]*b[1];Syz+=a[1]*b[2];Szx+=a[2]*b[0];Szy+=a[2]*b[1];Szz+=a[2]*b[2];G+=a[0]**2+a[1]**2+a[2]**2+b[0]**2+b[1]**2+b[2]**2}
 const N=[[Sxx+Syy+Szz,Syz-Szy,Szx-Sxz,Sxy-Syx],[Syz-Szy,Sxx-Syy-Szz,Sxy+Syx,Szx+Sxz],[Szx-Sxz,Sxy+Syx,-Sxx+Syy-Szz,Syz+Szy],[Sxy-Syx,Szx+Sxz,Syz+Szy,-Sxx-Syy+Szz]];
 const{values,vectors}=eig4(N);let m=0;for(let k=1;k<4;k++)if(values[k]>values[m])m=k;
 let[w,x,y,z]=[0,1,2,3].map(r=>vectors[r][m]);const L=Math.hypot(w,x,y,z);w/=L;x/=L;y/=L;z/=L;
 const R=[[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],[2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],[2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]];
 const apply=p=>{const a=[p.x-cp.x,p.y-cp.y,p.z-cp.z];return{x:R[0][0]*a[0]+R[0][1]*a[1]+R[0][2]*a[2]+cq.x,y:R[1][0]*a[0]+R[1][1]*a[1]+R[1][2]*a[2]+cq.y,z:R[2][0]*a[0]+R[2][1]*a[1]+R[2][2]*a[2]+cq.z}};
 const deviations=P.map((p,i)=>distance(apply(p),Q[i]));
 return{rmsd:Math.sqrt(deviations.reduce((s,d)=>s+d*d,0)/n),rmsdFromEigenvalue:Math.sqrt(Math.max(0,(G-2*values[m])/n)),rotation:R,apply,deviations,n}}

return {superpose};})();
__M['mirror-core']=(()=>{
const {superpose}=__M['superpose'];
const reflect = atoms => atoms.map(a=>({...a,x:-a.x}));
// Enumerate graph automorphisms; only proper rotations are allowed in fitting.
function checkMirror(atoms,bonds,orders=[],{limit=50000,tolerance=0.02}={}) {
 if(atoms.length<3 || atoms.length>80) throw Error('Use a molecule with 3–80 atoms for this teaching tool.');
 if(atoms.some(a=>![a.x,a.y,a.z].every(Number.isFinite))) throw Error('Invalid coordinates.');
 const n=atoms.length,edge=Array.from({length:n},()=>Array(n).fill(0));
 bonds.forEach(([i,j],k)=>edge[i][j]=edge[j][i]=orders[k]||1);
 const el=a=>a.element||a.el;
 const sig=atoms.map((a,i)=>el(a)+'|'+edge[i].map((v,j)=>v?`${el(atoms[j])}:${v}`:'').filter(Boolean).sort().join(','));
 const candidates=sig.map(s=>sig.flatMap((t,j)=>s===t?[j]:[]));
 const sequence=Array.from({length:n},(_,i)=>i).sort((a,b)=>candidates[a].length-candidates[b].length);
 const mirror=reflect(atoms),map=Array(n).fill(-1),used=new Set();
 let best=null,visits=0,complete=true,matched=false;
 function visit(depth){
  if(matched)return;
  if(++visits>limit){complete=false;return;}
  if(depth===n){const fit=superpose(map.map(j=>mirror[j]),atoms);if(!best||fit.rmsd<best.rmsd)best={...fit,mapping:[...map]};if(fit.rmsd<=tolerance)matched=true;return;}
  const i=sequence[depth];
  for(const j of candidates[i]){if(used.has(j))continue;let ok=true;for(let k=0;k<n;k++)if(map[k]>=0&&edge[i][k]!==edge[j][map[k]]){ok=false;break;}if(!ok)continue;
   map[i]=j;used.add(j);visit(depth+1);used.delete(j);map[i]=-1;if(matched||!complete)return;
  }
 }
 visit(0);
 return {matched,complete,best,visits,tolerance,mirror};
}

return {reflect,checkMirror};})();
__M['conformation-core']=(()=>{
const xyz=a=>[a.x,a.y,a.z];
const sub=(a,b)=>a.map((v,i)=>v-b[i]);
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const unit=a=>{const n=Math.hypot(...a);if(n<1e-8)throw Error('Bond axis has zero length.');return a.map(v=>v/n);};
function movingSide(n,bonds,front,back){
 const adj=Array.from({length:n},()=>[]);for(const [a,b] of bonds){if((a===front&&b===back)||(a===back&&b===front))continue;adj[a].push(b);adj[b].push(a);}
 const seen=new Set([back]),todo=[back];while(todo.length){for(const j of adj[todo.pop()])if(!seen.has(j)){seen.add(j);todo.push(j);}}
 return seen.has(front)?null:[...seen];
}
function eligibleBonds(atoms,bonds,orders){
 const adj=Array.from({length:atoms.length},()=>[]);bonds.forEach(([a,b])=>{adj[a].push(b);adj[b].push(a);});
 return bonds.flatMap(([a,b],k)=>{const el=i=>atoms[i].element||atoms[i].el;
  if(orders[k]!==1||el(a)==='H'||el(b)==='H'||adj[a].length<2||adj[b].length<2)return [];
  // Exclude amide-like C–N bonds; their partial double-bond character restricts rotation.
  if(((el(a)==='C'&&el(b)==='N')||(el(b)==='C'&&el(a)==='N'))&&bonds.some(([u,v],j)=>orders[j]===2&&((u===(el(a)==='C'?a:b)&&['O','S'].includes(el(v)))||(v===(el(a)==='C'?a:b)&&['O','S'].includes(el(u))))))return [];
  const moving=movingSide(atoms.length,bonds,a,b);return moving?[{front:a,back:b,moving,index:k}]:[];
 });
}
function rotateBond(atoms,front,back,moving,degrees){
 const p=xyz(atoms[front]),u=unit(sub(xyz(atoms[back]),p)),angle=degrees*Math.PI/180,c=Math.cos(angle),s=Math.sin(angle),set=new Set(moving);
 return atoms.map((a,i)=>{if(!set.has(i)||i===front||i===back)return {...a};const v=sub(xyz(a),p),w=cross(u,v),d=dot(u,v),r=v.map((v,k)=>p[k]+v*c+w[k]*s+u[k]*d*(1-c));return {...a,x:r[0],y:r[1],z:r[2]};});
}
// Looking from front towards back. Basis remains anchored to the original coordinates.
function projectionBasis(atoms,front,back){const u=unit(sub(xyz(atoms[back]),xyz(atoms[front]))),helper=Math.abs(u[2])<.9?[0,0,1]:[0,1,0],x=unit(cross(u,helper)),y=cross(x,u);return {u,x,y};}
function newmanBranches(atoms,bonds,front,back,basis){
 return [front,back].flatMap((center,side)=>bonds.flatMap(([a,b])=>{const j=a===center?b:b===center?a:-1;if(j<0||j===front||j===back)return [];const v=sub(xyz(atoms[j]),xyz(atoms[center])),x=dot(v,basis.x),y=dot(v,basis.y);return [{index:j,side,x,y,angle:Math.atan2(y,x),axial:Math.hypot(x,y)<1e-7}];}));
}

return {movingSide,eligibleBonds,rotateBond,projectionBasis,newmanBranches};})();
__M['energy-core']=(()=>{
const {projectionBasis,newmanBranches}=__M['conformation-core'];
function ethanePhase(atoms,bonds,orders,front,back){
 const el=i=>atoms[i].element||atoms[i].el;
 if(atoms.length!==8||el(front)!=='C'||el(back)!=='C'||bonds.length!==7||orders.some(o=>o!==1))return null;
 const adj=atoms.map(()=>[]);for(const [a,b] of bonds){adj[a].push(b);adj[b].push(a);}
 if(!adj[front].includes(back)||adj[front].length!==4||adj[back].length!==4)return null;
 if(atoms.some((_,i)=>i!==front&&i!==back&&(el(i)!=='H'||adj[i].length!==1||![front,back].includes(adj[i][0]))))return null;
 const branches=newmanBranches(atoms,bonds,front,back,projectionBasis(atoms,front,back));
 if(branches.length!==6||branches.some(b=>b.axial))return null;
 for(const side of [0,1]){const a=branches.filter(b=>b.side===side).map(b=>b.angle).sort((a,b)=>a-b);for(let i=0;i<3;i++){const gap=(a[(i+1)%3]-a[i]+Math.PI*2)%(Math.PI*2);if(Math.abs(gap-2*Math.PI/3)>.08)return null;}}
 return (branches.find(b=>b.side===1).angle-branches.find(b=>b.side===0).angle)*180/Math.PI;
}
// Approximate threefold torsional teaching potential: 0 staggered, 12 kJ/mol eclipsed.
// Our projection basis has opposite angular handedness to the positive Rodrigues rotation.
function ethaneEnergy(rotation,phase){return 6*(1+Math.cos(3*(phase-rotation)*Math.PI/180));}

return {ethanePhase,ethaneEnergy};})();
__M['energy-panel']=(()=>{
const {ethanePhase,ethaneEnergy}=__M['energy-core'];
function energyPanel(body,atoms,bonds,orders,onAngle){
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

return {energyPanel};})();
__M['sugar-core']=(()=>{
// OH positions from C2 to C5 in the canonical open-chain Fischer projection.
// R/L here mean right/left on the page, NOT CIP stereodescriptors.
const ALDOHEXOSES={allose:'RRRR',altrose:'LRRR',glucose:'RLRR',mannose:'LLRR',gulose:'RRLR',idose:'LRLR',galactose:'RLLR',talose:'LLLR'};
function validate(pattern){if(!/^[RL]{4}$/.test(pattern))throw Error('An aldohexose needs four left/right OH positions.');return pattern;}
function mirrorSugar(pattern){return [...validate(pattern)].map(s=>s==='R'?'L':'R').join('');}
function flipSugar(pattern,carbon){validate(pattern);if(!Number.isInteger(carbon)||carbon<2||carbon>5)throw Error('Select stereocenter C2–C5.');return [...pattern].map((s,i)=>i===carbon-2?(s==='R'?'L':'R'):s).join('');}
function sugarName(pattern){validate(pattern);const series=pattern[3]==='R'?'D':'L',d=series==='D'?pattern:mirrorSugar(pattern);return series+'-'+Object.keys(ALDOHEXOSES).find(n=>ALDOHEXOSES[n]===d);}
function sugarRelationship(a,b){validate(a);validate(b);const differences=[...a].flatMap((s,i)=>s!==b[i]?[i+2]:[]);return {differences,type:differences.length===0?'Identical configuration':differences.length===4?'Enantiomers':differences.length===1?'Epimers (diastereomers)':'Diastereomers'};}

return {ALDOHEXOSES,validate,mirrorSugar,flipSugar,sugarName,sugarRelationship};})();
__M['chair-core']=(()=>{
const SUGARS={glucose:['down','up','down'],mannose:['up','up','down'],galactose:['down','up','up']};
function chairModel(sugar='glucose',anomer='beta',flipped=false){
 if(!SUGARS[sugar]||!['alpha','beta'].includes(anomer))throw Error('Unsupported sugar or anomer.');
 const signs=[anomer==='alpha'?-1:1,...SUGARS[sugar].map(v=>v==='up'?1:-1),1];
 const ring=Array.from({length:6},(_,i)=>{const t=i*Math.PI/3;return {x:1.45*Math.cos(t),y:1.45*Math.sin(t),z:(i%2?1:-1)*.25*(flipped?-1:1),label:i===5?'O5':'C'+(i+1)};});
 const groups=signs.map((sign,i)=>{const axialSign=(i%2?1:-1)*(flipped?-1:1),axial=sign===axialSign,t=i*Math.PI/3,p=ring[i],r=axial?0:1.2;return {carbon:i+1,label:i===4?'CH₂OH':'OH',direction:sign===1?'up':'down',position:axial?'axial':'equatorial',x:p.x+r*Math.cos(t),y:p.y+r*Math.sin(t),z:p.z+sign*(axial?1.25:.48)};});
 return {ring,groups,chair:flipped?'¹C₄':'⁴C₁',name:(anomer==='alpha'?'α':'β')+'-D-'+({glucose:'glucopyranose',mannose:'mannopyranose',galactose:'galactopyranose'}[sugar])};
}

return {SUGARS,chairModel};})();
__M['electron-core']=(()=>{
const atom=(el,x,y,z)=>({el,x,y,z});
function electronModel(id){
 const origin=el=>atom(el,0,0,0),r=Math.PI/180;let model;
 if(id==='water'){const a=52.25*r;model={name:'Water · H₂O',atoms:[origin('O'),atom('H',Math.sin(a)*.96,0,-Math.cos(a)*.96),atom('H',-Math.sin(a)*.96,0,-Math.cos(a)*.96)],bonds:[[0,1,1],[0,2,1]],lonePairs:[{atom:0,direction:[0,.82,.57]},{atom:0,direction:[0,-.82,.57]}],geometry:'Bent',electronGeometry:'Tetrahedral',angle:104.5};}
 else if(id==='ammonia'){const c=Math.sqrt((Math.cos(107*r)+.5)/1.5),s=Math.sqrt(1-c*c);model={name:'Ammonia · NH₃',atoms:[origin('N'),...Array.from({length:3},(_,i)=>atom('H',1.01*s*Math.cos(i*2*Math.PI/3),1.01*s*Math.sin(i*2*Math.PI/3),-1.01*c))],bonds:[[0,1,1],[0,2,1],[0,3,1]],lonePairs:[{atom:0,direction:[0,0,1]}],geometry:'Trigonal pyramidal',electronGeometry:'Tetrahedral',angle:107};}
 else if(id==='methane'){const s=1.09/Math.sqrt(3);model={name:'Methane · CH₄',atoms:[origin('C'),...[ [1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(p=>atom('H',...p.map(v=>v*s)))],bonds:[[0,1,1],[0,2,1],[0,3,1],[0,4,1]],lonePairs:[],geometry:'Tetrahedral',electronGeometry:'Tetrahedral',angle:109.47};}
 else if(id==='co2')model={name:'Carbon dioxide · CO₂',atoms:[origin('C'),atom('O',-1.16,0,0),atom('O',1.16,0,0)],bonds:[[0,1,2],[0,2,2]],lonePairs:[{atom:1,direction:[-.6,.8,0]},{atom:1,direction:[-.6,-.8,0]},{atom:2,direction:[.6,.8,0]},{atom:2,direction:[.6,-.8,0]}],geometry:'Linear',electronGeometry:'Linear',angle:180};
 else throw Error('Unknown electron teaching model.');
 return {...model,central:0,id};
}
function electronCounts(model){const bondingPairs=model.bonds.reduce((s,b)=>s+b[2],0),lonePairs=model.lonePairs.length,centralLonePairs=model.lonePairs.filter(p=>p.atom===model.central).length,centralBonds=model.bonds.filter(b=>b[0]===model.central||b[1]===model.central).length;return {bondingPairs,lonePairs,centralLonePairs,domains:centralBonds+centralLonePairs,valenceElectrons:2*(bondingPairs+lonePairs)};}

return {electronModel,electronCounts};})();
__M['mirror-lab']=(()=>{
const {reflect,checkMirror}=__M['mirror-core'];
const style=document.createElement('style');style.textContent=`#mirrorLab{background:var(--sf,#17231f);color:var(--tx,#fff);border:1px solid var(--bd,#777);border-radius:16px;width:min(1050px,94vw);max-height:90vh;overflow:auto;padding:22px}#mirrorLab::backdrop{background:#000a}#mirrorLab .mlviews{display:flex;gap:14px}#mirrorLab .mlpane{flex:1;min-width:0}#mirrorLab canvas{width:100%;height:310px;display:block;touch-action:none}#mirrorLab button{margin:5px}#mirrorLab p{line-height:1.5} @media(max-width:650px){#mirrorLab .mlviews{flex-direction:column}#mirrorLab canvas{height:220px}}`;document.head.append(style);
const button=document.createElement('button');button.className='ib';button.setAttribute('aria-label','Mirror lab');button.innerHTML='<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 3v18" stroke-dasharray="2 2.4"/><path d="M9 7 4 17h5zM15 7l5 10h-5z"/></svg>';button.title='Reflect a molecule, rotate both views, and test rigid superposition';document.querySelector('#inf').before(button);
button.onclick=()=>{
 const source=window.lauraMirrorSource();if(!source)return;
 const T=window.THREE;if(!T?.OrbitControls){alert('The 3D engine is not available. Check your internet connection and reload.');return;}
 const atoms=source.atoms.map(a=>({...a})),bonds=source.bonds.map(b=>[...b]),orders=[...source.orders];
 if(atoms.length>80||atoms.length<3){alert('Mirror lab supports molecules with 3–80 atoms. Choose a small molecule such as alanine, methane, or water.');return;}
 const dialog=document.createElement('dialog');dialog.id='mirrorLab';dialog.setAttribute('aria-label','Molecular mirror lab');
 dialog.innerHTML='<button class="mlclose" style="float:right;position:sticky;top:0;z-index:5" aria-label="Close mirror lab">×</button><h2>Chirality & mirror images</h2><details class="td"><summary><b>About this lab</b></summary><p><b>What:</b> reflects the loaded molecule through a plane to make its mirror image, lets you rotate both views independently, and tests whether the two can be superimposed by rotation and translation alone. <b>Read it:</b> a chiral molecule is not superimposable on its mirror image; an achiral one is. <b>Limits:</b> rigid geometry only (no conformational search), 3 to 80 atoms, and the answer depends on the loaded conformation.</p></details><p class="mlname"></p><p>Drag each view to rotate it independently; scroll to zoom. Reflection reverses the x coordinate. Rotating a view does not change the molecule.</p><div class="mlviews"><section class="mlpane"><h3>Original</h3><div class="mlhost"></div></section><section class="mlpane"><h3>Mirror image</h3><div class="mlhost"></div></section></div><button class="mlcheck">Check superimposability</button><button class="mlreset">Reset views</button><p class="mlresult" role="status" aria-live="polite">Ready to compare using rotation and translation only.</p><p><small>This tests the current rigid geometry and inferred bond graph, allowing equivalent atoms to exchange. A mismatch is not a definitive chirality assignment: missing atoms, coordinate noise, or different conformations can affect the result. No R/S assignment or conformational search is performed. Matching tolerance: 0.02 Å RMSD.</small></p>';
 dialog.querySelector('.mlname').textContent=source.name+' · '+atoms.length+' atoms';document.body.append(dialog);dialog.showModal();{let dn=false;const out=e=>{const r=dialog.getBoundingClientRect();return e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom};dialog.addEventListener('mousedown',e=>{dn=out(e)});dialog.addEventListener('click',e=>{if(dn&&out(e))dialog.close();dn=false});}
 const views=[],colors={H:0xffffff,C:0x929b9d,N:0x3050f8,O:0xff3030,F:0x50d050,Cl:0x30cc30,Br:0xa84025,I:0xa050d0,S:0xffd62a,P:0xff8000};
 const center=atoms.reduce((c,a)=>c.add(new T.Vector3(a.x,a.y,a.z)),new T.Vector3()).multiplyScalar(1/atoms.length);
 const radius=Math.max(1,...atoms.map(a=>new T.Vector3(a.x,a.y,a.z).distanceTo(center)));
 function build(A,tint,opacity=1){const group=new T.Group();const material=color=>new T.MeshStandardMaterial({color,transparent:opacity<1,opacity,roughness:.5});
  A.forEach(a=>{const sphere=new T.Mesh(new T.SphereGeometry((a.element||a.el)==='H'?.18:.3,20,14),material(tint??colors[a.element||a.el]??0xff80cc));sphere.position.set(a.x,a.y,a.z);group.add(sphere);});
  bonds.forEach(([i,j])=>{const a=new T.Vector3(A[i].x,A[i].y,A[i].z),b=new T.Vector3(A[j].x,A[j].y,A[j].z),d=b.clone().sub(a);if(d.length()<1e-8)return;const mesh=new T.Mesh(new T.CylinderGeometry(.075,.075,d.length(),12),material(tint??0xb5bfc0));mesh.position.copy(a.add(b).multiplyScalar(.5));mesh.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),d.normalize());group.add(mesh);});return group;}
 function clear(v){for(const g of v.groups){v.scene.remove(g);g.traverse(o=>{o.geometry?.dispose();o.material?.dispose();});}v.groups=[];}
 function display(v,A,overlay=false){clear(v);const c=A.reduce((s,a)=>s.add(new T.Vector3(a.x,a.y,a.z)),new T.Vector3()).multiplyScalar(1/A.length);const g=build(A);g.position.copy(c).negate();v.scene.add(g);v.groups.push(g);if(overlay){const ref=build(atoms,0x00ffff,.3);ref.position.copy(center).negate();v.scene.add(ref);v.groups.push(ref);}}
 try{dialog.querySelectorAll('.mlhost').forEach(host=>{const renderer=new T.WebGLRenderer({antialias:true});renderer.setPixelRatio(Math.min(devicePixelRatio,2));host.append(renderer.domElement);const scene=new T.Scene();scene.background=new T.Color(0x10231e);scene.add(new T.HemisphereLight(0xffffff,0x445544,1.8));const light=new T.DirectionalLight(0xffffff,1);light.position.set(4,7,8);scene.add(light);const camera=new T.PerspectiveCamera(42,1,.01,1000),controls=new T.OrbitControls(camera,renderer.domElement);const v={renderer,scene,camera,controls,groups:[]};views.push(v);});}catch(error){views.forEach(v=>{v.controls.dispose();v.renderer.dispose();});dialog.remove();alert('Unable to start mirror lab: '+error.message);return;}
 const reset=()=>{views.forEach((v,i)=>{display(v,i?reflect(atoms):atoms);v.camera.position.set(0,0,radius*3.8+2);v.controls.target.set(0,0,0);v.controls.update();});dialog.querySelector('.mlresult').textContent='Ready to compare using rotation and translation only.';};reset();
 let frame;function animate(){for(const v of views){const canvas=v.renderer.domElement,w=canvas.clientWidth,h=canvas.clientHeight;v.renderer.setSize(w,h,false);v.camera.aspect=w/h;v.camera.updateProjectionMatrix();v.controls.update();v.renderer.render(v.scene,v.camera);}frame=requestAnimationFrame(animate);}animate();
 dialog.querySelector('.mlcheck').onclick=()=>{const out=dialog.querySelector('.mlresult');try{const r=checkMirror(atoms,bonds,orders);out.textContent=r.matched?'Superimposable within tolerance for the atoms present.':r.complete?'Not superimposable as this rigid geometry.':'Inconclusive: the equivalent-atom search reached its limit.';if(r.best){out.textContent+=' Best tested RMSD: '+r.best.rmsd.toFixed(4)+' Å. Right: fitted mirror with original shown in transparent cyan. Bonds are shown as single sticks for clarity.';display(views[1],r.mirror.map(a=>({...a,...r.best.apply(a)})),true);views[1].camera.position.set(0,0,radius*3.8+2);views[1].controls.target.set(0,0,0);views[1].controls.update();}}catch(e){out.textContent=e.message;}};
 dialog.querySelector('.mlreset').onclick=reset;dialog.querySelector('.mlclose').onclick=()=>dialog.close();dialog.addEventListener('close',()=>{cancelAnimationFrame(frame);views.forEach(v=>{clear(v);v.controls.dispose();v.renderer.dispose();v.renderer.forceContextLoss();});dialog.remove();button.focus();},{once:true});
};

return {};})();
__M['conformation-lab']=(()=>{
const {energyPanel}=__M['energy-panel'];
const {eligibleBonds,movingSide,rotateBond,projectionBasis,newmanBranches}=__M['conformation-core'];
const style=document.createElement('style');style.textContent=`#confLab{background:var(--sf,#17231f);color:var(--tx,#fff);border:1px solid var(--bd,#777);border-radius:16px;width:min(1050px,94vw);max-height:90vh;overflow:auto;padding:22px}#confLab::backdrop{background:#000a}#confLab .clviews{display:flex;gap:16px}#confLab .clpane{flex:1;min-width:0}#confLab canvas,#confLab svg{width:100%;height:320px;display:block;background:#10231e;border-radius:8px}#confLab canvas{touch-action:none}#confLab button{margin:5px}#confLab select{max-width:100%;width:auto}#confLab p{line-height:1.5}#confLab .clslider{width:min(460px,90%)}#confLab .clcontrols{display:flex;gap:8px;align-items:center;flex-wrap:wrap}@media(max-width:650px){#confLab .clviews{flex-direction:column}#confLab canvas,#confLab svg{height:250px}}`;document.head.append(style);
const launch=document.createElement('button');launch.className='ib';launch.setAttribute('aria-label','Bond rotation');launch.innerHTML='<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M20 12a8 8 0 1 1-2.6-5.9"/><path d="M20 4v4.5h-4.5"/><circle cx="12" cy="12" r="1.6"/></svg>';launch.title='Explore bond conformations and Newman projections';document.querySelector('#inf').before(launch);
launch.onclick=()=>{
 const source=window.lauraMirrorSource();if(!source)return;
 const atoms=source.atoms.map(a=>({...a})),bonds=source.bonds.map(b=>[...b]),orders=[...source.orders],T=window.THREE;
 if(atoms.length>80||atoms.length<3){alert('Bond rotation supports molecules with 3–80 atoms. Choose a small molecule such as ethane, alanine, or butane.');return;}
 const dialog=document.createElement('dialog');dialog.id='confLab';dialog.setAttribute('aria-label','Bond rotation and Newman projection');dialog.innerHTML='<button class="clclose" style="float:right;position:sticky;top:0;z-index:5" aria-label="Close bond rotation">×</button><h2>Bond rotation & Newman projection</h2><details class="td"><summary><b>About this lab</b></summary><p><b>What:</b> rotates the two ends of a chosen single bond against each other and shows the matching Newman projection. For ethane an approximate threefold curve is shown: E = 6[1 + cos(3 phi)] kJ/mol, relative to the staggered form, with a barrier of about 12 kJ/mol. <b>Limits:</b> purely geometric rotation with no collision check, relaxation or force field; energies only for complete ethane; ring and amide-like C-N bonds are excluded; Newman labels name adjacent atoms, not whole groups. Reference: OpenStax Organic Chemistry, section 3.6.</p></details><p class="clname"></p><div class="clbody"></div>';dialog.querySelector('.clname').textContent=source.name;document.body.append(dialog);dialog.showModal();{let dn=false;const out=e=>{const r=dialog.getBoundingClientRect();return e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom};dialog.addEventListener('mousedown',e=>{dn=out(e)});dialog.addEventListener('click',e=>{if(dn&&out(e))dialog.close();dn=false});}
 let cleanup=()=>{};dialog.querySelector('.clclose').onclick=()=>dialog.close();dialog.addEventListener('close',()=>{cleanup();dialog.remove();launch.focus();},{once:true});
 const body=dialog.querySelector('.clbody');
 if(atoms.length>80||atoms.length<3){body.textContent='Choose a small molecule with 3–80 atoms. Ethane (staggered) is a good starting example.';return;}
 const choices=eligibleBonds(atoms,bonds,orders);
 if(!choices.length){body.textContent='No suitable bond found. Choose Ethane (staggered) to begin. This lab requires a non-ring single bond between non-hydrogen atoms with visible substituents at both ends; amide-like C–N bonds are excluded.';return;}
 if(!T?.OrbitControls){body.textContent='The 3D engine is unavailable. Check your connection and reload.';return;}
 body.innerHTML='<p>Rotate the rear fragment around a bond. The Newman projection looks from the front atom toward the rear atom and updates with the same coordinates.</p><div class="clcontrols"><label>Bond axis <select class="clbond" aria-label="Bond axis"></select></label><button class="clreverse">Reverse viewing direction</button><button class="clview">Look down bond</button></div><div class="clviews"><section class="clpane"><h3>3D conformation</h3><div class="clhost"></div></section><section class="clpane"><h3>Newman projection</h3><svg class="clsvg" viewBox="0 0 420 320" role="img" aria-label="Newman projection of selected bond"></svg></section></div><p class="cldirection"></p><label>Rotation from loaded conformation: <output class="clangle">0°</output><br><input class="clslider" aria-label="Bond rotation angle" type="range" min="-180" max="180" step="1" value="0"></label><div><button class="clminus">−60°</button><button class="clplus">+60°</button><button class="clreset">Reset rotation</button></div><p class="clstatus" role="status" aria-live="polite"></p><p><small>Cyan: front substituents; orange dashed lines: rear substituents. Labels identify the adjacent atom and its index, not an entire substituent group. Bonds appear as single sticks in 3D. Rotation is a rigid geometric operation without energy minimization or collision prevention. Eligibility uses inferred bond orders. Missing hydrogens are not added. Changes stay in this lab; the loaded structure is preserved.</small></p>';
 const q=s=>body.querySelector(s),label=i=>(atoms[i].element||atoms[i].el)+' #'+(i+1);
 choices.forEach((b,i)=>{const option=document.createElement('option');option.value=i;option.textContent=label(b.front)+' ↔ '+label(b.back);q('.clbond').append(option);});
 let front,back,moving,basis,coords=atoms,angle=0,reversed=false;
 const center=atoms.reduce((s,a)=>s.add(new T.Vector3(a.x,a.y,a.z)),new T.Vector3()).multiplyScalar(1/atoms.length),radius=Math.max(1,...atoms.map(a=>new T.Vector3(a.x,a.y,a.z).distanceTo(center)));
 let renderer;try{renderer=new T.WebGLRenderer({antialias:true});}catch(e){body.textContent='Could not initialize the 3D view: '+e.message;return;}
 renderer.setPixelRatio(Math.min(devicePixelRatio,2));q('.clhost').append(renderer.domElement);const scene=new T.Scene();scene.background=new T.Color(0x10231e);scene.add(new T.HemisphereLight(0xffffff,0x445544,1.5));const light=new T.DirectionalLight(0xffffff,1);light.position.set(4,7,8);scene.add(light);const camera=new T.PerspectiveCamera(42,1,.01,2000),controls=new T.OrbitControls(camera,renderer.domElement);camera.position.set(radius*1.5,radius,radius*3.8+2);
 const group=new T.Group();group.position.copy(center).negate();scene.add(group);const colors={H:0xffffff,C:0x929b9d,N:0x3050f8,O:0xff3030,F:0x50d050,Cl:0x30cc30,Br:0xa84025,I:0xa050d0,S:0xffd62a,P:0xff8000};
 const spheres=atoms.map(a=>{const mesh=new T.Mesh(new T.SphereGeometry((a.element||a.el)==='H'?.18:.3,18,12),new T.MeshStandardMaterial({color:colors[a.element||a.el]??0xff80cc,roughness:.5}));group.add(mesh);return mesh;});
 const sticks=bonds.map(()=>{const mesh=new T.Mesh(new T.CylinderGeometry(.07,.07,1,10),new T.MeshStandardMaterial({color:0xaabbb6,roughness:.5}));group.add(mesh);return mesh;});
 function svg(tag,attrs,text){const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;q('.clsvg').append(e);return e;}
 const updateEnergy=energyPanel(body,atoms,bonds,orders,a=>{angle=a;update();});
 function update(){updateEnergy(front,back,angle);coords=rotateBond(atoms,front,back,moving,angle);q('.clslider').value=angle;q('.clangle').textContent=angle+'°';spheres.forEach((m,i)=>m.position.set(coords[i].x,coords[i].y,coords[i].z));sticks.forEach((m,k)=>{const [i,j]=bonds[k],a=spheres[i].position,b=spheres[j].position,d=b.clone().sub(a);m.position.copy(a).add(b).multiplyScalar(.5);m.scale.y=d.length();m.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),d.normalize());m.material.color.setHex((i===front&&j===back)||(i===back&&j===front)?0xffcc55:0xaabbb6);});
  q('.clsvg').replaceChildren();svg('circle',{cx:210,cy:150,r:38,fill:'#18352c',stroke:'#ffbd76','stroke-width':3});const branches=newmanBranches(coords,bonds,front,back,basis);
  for(const side of [0,1])for(const b of branches.filter(b=>b.side===side&&!b.axial)){const r=side?38:0,c=Math.cos(b.angle),s=-Math.sin(b.angle),color=side?'#ffbd76':'#6de5ee';svg('line',{x1:210+c*r,y1:150+s*r,x2:210+c*104,y2:150+s*104,stroke:color,'stroke-width':side?2:3,...(side?{'stroke-dasharray':'7 4'}:{})});const lr=side?145:98;svg('text',{x:210+c*lr,y:150+s*lr+5,fill:color,'font-size':14,'text-anchor':'middle'},label(b.index));}
  svg('circle',{cx:210,cy:150,r:6,fill:'#6de5ee'});svg('text',{x:210,y:310,fill:'#d4e7df','text-anchor':'middle','font-size':13},'Front '+label(front)+' → Rear '+label(back));
  q('.cldirection').textContent='Looking '+label(front)+' → '+label(back)+'. The rear fragment rotates; the front fragment stays fixed. Drag the 3D view to inspect from another angle.';
  const f=branches.filter(b=>b.side===0&&!b.axial),r=branches.filter(b=>b.side===1&&!b.axial);let status='';if(f.length===3&&r.length===3){const differences=f.map(a=>Math.min(...r.map(b=>Math.abs(Math.atan2(Math.sin(a.angle-b.angle),Math.cos(a.angle-b.angle)))*180/Math.PI)));if(differences.every(v=>v<5))status='Near eclipsed alignment.';else if(differences.every(v=>Math.abs(v-60)<5))status='Near staggered alignment.';else status='Intermediate alignment.';}q('.clstatus').textContent=status+(branches.some(b=>b.axial)?' A substituent lies along the viewing axis; its projected direction is undefined.':'');
 }
 function choose(){const b=choices[Number(q('.clbond').value)];front=reversed?b.back:b.front;back=reversed?b.front:b.back;moving=movingSide(atoms.length,bonds,front,back);basis=projectionBasis(atoms,front,back);angle=0;update();}
 q('.clbond').onchange=()=>{reversed=false;choose();};q('.clreverse').onclick=()=>{reversed=!reversed;choose();};q('.clslider').oninput=e=>{angle=Number(e.target.value);update();};q('.clminus').onclick=()=>{angle=angle-60<-180?angle+300:angle-60;update();};q('.clplus').onclick=()=>{angle=angle+60>180?angle-300:angle+60;update();};q('.clreset').onclick=()=>{angle=0;update();};q('.clview').onclick=()=>{const p=new T.Vector3(atoms[front].x,atoms[front].y,atoms[front].z).sub(center),u=new T.Vector3(...basis.u);controls.target.copy(p);camera.up.set(...basis.y);camera.position.copy(p).addScaledVector(u,-(radius*3.8+2));controls.update();};choose();
 let frame;function animate(){const canvas=renderer.domElement,w=canvas.clientWidth,h=canvas.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();controls.update();renderer.render(scene,camera);frame=requestAnimationFrame(animate);}animate();cleanup=()=>{cancelAnimationFrame(frame);controls.dispose();group.traverse(m=>{m.geometry?.dispose();m.material?.dispose();});renderer.dispose();renderer.forceContextLoss();};
};

return {};})();
__M['sugar-lab']=(()=>{
const {ALDOHEXOSES,mirrorSugar,flipSugar,sugarName,sugarRelationship}=__M['sugar-core'];
const style=document.createElement('style');style.textContent=`#sugarLab{background:var(--sf,#17231f);color:var(--tx,#fff);border:1px solid var(--bd,#777);border-radius:16px;width:min(1000px,94vw);max-height:90vh;overflow:auto;padding:24px}#sugarLab::backdrop{background:#000a}#sugarLab .sgcols{display:flex;gap:20px}#sugarLab .sgcol{flex:1;min-width:0}#sugarLab svg{width:100%;height:350px;background:#10231e;border-radius:12px}#sugarLab button{margin:4px}#sugarLab select{width:auto;max-width:100%}#sugarLab p{line-height:1.5}#sugarLab .sgresult{padding:14px;border:1px solid var(--ac,#64ddd0);border-radius:8px}@media(max-width:650px){#sugarLab .sgcols{flex-direction:column}}`;document.head.append(style);
const launch=document.createElement('button');launch.className='ib';launch.setAttribute('aria-label','Sugar lab');launch.innerHTML='<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 3v18"/><path d="M7 7.5h10M7 12h10M7 16.5h10"/></svg>';launch.title='Fischer projections and aldohexose configurations';document.querySelector('#inf').before(launch);
launch.onclick=()=>{
 let reference=ALDOHEXOSES.glucose,current=reference;
 const dialog=document.createElement('dialog');dialog.id='sugarLab';dialog.setAttribute('aria-label','Sugar configurations and Fischer projections');dialog.innerHTML='<button class="sgclose" style="float:right;position:sticky;top:0;z-index:5" aria-label="Close sugar lab">×</button><h2>Sugars & Fischer projections</h2><details class="td"><summary><b>About this lab</b></summary><p><b>What:</b> the 16 open-chain aldohexoses as Fischer projections. Invert C2 to C5 one at a time or mirror every centre; the lab names the relationship to the reference (identical, epimers, other diastereomers, enantiomers). <b>Limits:</b> predefined configurations, not converted from imported 3D coordinates; no ring forms; D/L does not predict the sign of optical rotation. References: OpenStax Organic Chemistry, sections 25.3 and 25.4.</p></details><p>Explore open-chain aldohexoses</p><label>Reference sugar <select class="sgpreset" aria-label="Reference sugar"></select></label><div class="sgcols"><section class="sgcol"><h3 class="sgrefname"></h3><svg class="sgref" viewBox="0 0 340 350" role="img"></svg></section><section class="sgcol"><h3 class="sgname"></h3><svg class="sgcurrent" viewBox="0 0 340 350" role="img"></svg></section></div><div class="sgflips" aria-label="Invert a stereocenter"></div><div><button class="sgmirror">Mirror all centers</button><button class="sgreset">Reset to reference</button><button class="sgexport">Download Fischer SVG</button></div><p class="sgresult" role="status" aria-live="polite"></p><p class="sgseries"></p><p>In a Fischer projection, horizontal bonds point toward you and vertical bonds point away. The chain is vertical with CHO at the top. Swapping H and OH at one cross inverts that stereocenter. Highlighted crosses differ from the reference.</p><p><small>These are defined open-chain sugar configurations, not projections calculated from the molecule loaded in the main viewer. The lab does not model ring closure, α/β anomers, chair conformations, or solution populations. D/L describes configuration; it does not specify the sign of optical rotation. <a href="https://openstax.org/books/organic-chemistry/pages/25-4-configurations-of-the-aldoses" target="_blank" rel="noopener noreferrer">Fischer configurations</a> · <a href="https://openstax.org/books/organic-chemistry/pages/25-3-d-l-sugars" target="_blank" rel="noopener noreferrer">D/L convention</a>.</small></p>';
 document.body.append(dialog);const q=s=>dialog.querySelector(s);
 for(const [name,p] of Object.entries(ALDOHEXOSES))for(const series of ['D','L']){const option=document.createElement('option');option.value=series==='D'?p:mirrorSugar(p);option.textContent=series+'-'+name;q('.sgpreset').append(option);}q('.sgpreset').value=reference;
 for(let c=2;c<=5;c++){const button=document.createElement('button');button.textContent='Invert C'+c;button.onclick=()=>{current=flipSugar(current,c);render();};q('.sgflips').append(button);}
 function draw(target,p,differences){target.replaceChildren();target.setAttribute('aria-label',sugarName(p)+' Fischer projection; OH sides C2 to C5: '+[...p].map(s=>s==='R'?'right':'left').join(', '));const add=(tag,attrs,text)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;target.append(e);};
  add('rect',{width:340,height:350,fill:'#10231e'});add('line',{x1:170,y1:46,x2:170,y2:300,stroke:'#acc5bc','stroke-width':3});
  const text=(x,y,t,fill='#edf6f2')=>add('text',{x,y,fill,'text-anchor':'middle','font-size':19,'font-family':'Arial, sans-serif'},t);
  text(170,36,'CHO');text(285,36,'C1','#94b4a7');text(170,326,'CH₂OH');text(285,326,'C6','#94b4a7');
  [...p].forEach((s,i)=>{const y=85+i*60,changed=differences.includes(i+2);if(changed)add('rect',{x:50,y:y-23,width:245,height:43,rx:8,fill:'#413729','fill-opacity':.45,stroke:'#ffbd76'});add('line',{x1:105,x2:235,y1:y,y2:y,stroke:changed?'#ffbd76':'#6de5ee','stroke-width':3});add('circle',{cx:170,cy:y,r:4,fill:'#edf6f2'});text(77,y+6,s==='L'?'OH':'H');text(260,y+6,s==='R'?'OH':'H');text(310,y+6,'C'+(i+2),'#94b4a7');});
 }
 function render(){const r=sugarRelationship(reference,current);q('.sgrefname').textContent='Reference · '+sugarName(reference);q('.sgname').textContent='Your configuration · '+sugarName(current);draw(q('.sgref'),reference,[]);draw(q('.sgcurrent'),current,r.differences);q('.sgresult').textContent=r.type+(r.differences.length?' · Differ at '+r.differences.map(c=>'C'+c).join(', ')+'.':'.')+' '+(r.differences.length===4?'Every stereocenter is inverted.':r.differences.length===1?'Only one stereocenter is inverted.':'');q('.sgseries').textContent=(current[3]==='R'?'D':'L')+' series: the OH at C5, farthest from the carbonyl, is on the '+(current[3]==='R'?'right':'left')+'. OH positions C2 → C5: '+[...current].map(s=>s==='R'?'right':'left').join(' · ')+'.';}
 q('.sgpreset').onchange=e=>{reference=e.target.value;current=reference;render();};q('.sgmirror').onclick=()=>{current=mirrorSugar(current);render();};q('.sgreset').onclick=()=>{current=reference;render();};q('.sgexport').onclick=()=>{const copy=q('.sgcurrent').cloneNode(true);copy.setAttribute('xmlns','http://www.w3.org/2000/svg');const title=document.createElementNS('http://www.w3.org/2000/svg','title');title.textContent=sugarName(current)+' — open-chain Fischer projection';copy.prepend(title);const blob=new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=sugarName(current)+'-Fischer.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
 q('.sgclose').onclick=()=>dialog.close();dialog.addEventListener('close',()=>{dialog.remove();launch.focus();},{once:true});render();dialog.showModal();{let dn=false;const out=e=>{const r=dialog.getBoundingClientRect();return e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom};dialog.addEventListener('mousedown',e=>{dn=out(e)});dialog.addEventListener('click',e=>{if(dn&&out(e))dialog.close();dn=false});}
};

return {};})();
__M['chair-lab']=(()=>{
const {chairModel}=__M['chair-core'];
const style=document.createElement('style');style.textContent=`#chairLab{background:var(--sf,#17231f);color:var(--tx,#fff);border:1px solid var(--bd,#777);border-radius:16px;width:min(1050px,94vw);max-height:90vh;overflow:auto;padding:24px}#chairLab::backdrop{background:#000a}#chairLab .chviews{display:flex;gap:16px}#chairLab .chpane{flex:1;min-width:0}#chairLab canvas,#chairLab svg{width:100%;height:340px;display:block;background:#10231e;border-radius:10px}#chairLab canvas{touch-action:none}#chairLab button,#chairLab select{margin:5px;width:auto}#chairLab p{line-height:1.5}#chairLab table{width:100%;border-collapse:collapse}#chairLab td,#chairLab th{text-align:left;padding:6px;border-bottom:1px solid var(--bd,#444)}@media(max-width:650px){#chairLab .chviews{flex-direction:column}}`;document.head.append(style);
const launch=document.createElement('button');launch.className='ib';launch.setAttribute('aria-label','Chair lab');launch.innerHTML='<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M3 15l5-6 5 4 8-7"/><path d="M3 15l4 4 6-6M13 13l4 4 4-9"/></svg>';launch.title='Cyclic sugars, anomers and chair flips';document.querySelector('#inf').before(launch);
launch.onclick=()=>{
 const T=window.THREE;const dialog=document.createElement('dialog');dialog.id='chairLab';dialog.setAttribute('aria-label','Cyclic sugar chair lab');dialog.innerHTML='<button class="chclose" style="float:right;position:sticky;top:0;z-index:5" aria-label="Close chair lab">×</button><h2>Cyclic sugars & chair flips</h2><details class="td"><summary><b>About this lab</b></summary><p><b>What:</b> alpha/beta D-glucopyranose, D-mannopyranose and D-galactopyranose as a Haworth projection, an idealised rotatable 3D chair and a position table. Flipping 4C1 and 1C4 swaps axial and equatorial while up/down and the anomer stay fixed. <b>Limits:</b> a schematic with hydrogens omitted; no energies, ring-closure reaction or solution populations. Reference: OpenStax Organic Chemistry, section 25.5.</p></details><p>D-pyranose teaching models</p><label>Sugar <select class="chsugar" aria-label="Cyclic sugar"><option value="glucose">D-glucose</option><option value="mannose">D-mannose</option><option value="galactose">D-galactose</option></select></label><label>Anomer <select class="chanomer" aria-label="Anomer"><option value="beta">β (beta)</option><option value="alpha">α (alpha)</option></select></label><button class="chflip">Flip chair</button><button class="chreset">Reset camera</button><h3 class="chname"></h3><div class="chviews"><section class="chpane"><h3>Haworth projection</h3><svg class="chsvg" viewBox="0 0 420 340" role="img" aria-label="Haworth projection"></svg></section><section class="chpane"><h3 class="chchair"></h3><div class="chhost"></div></section></div><p class="chstatus" role="status" aria-live="polite"></p><table><thead><tr><th>Carbon</th><th>Group</th><th>Ring face</th><th>Chair position</th></tr></thead><tbody></tbody></table><p>Chair flipping exchanges axial and equatorial positions while retaining each group’s up/down orientation and the α/β configuration. For these D sugars, α has C1–OH down; β has C1–OH up. Changing anomer changes configuration at C1; it is not a chair flip.</p><p><small>The 3D view is an idealized schematic with labeled OH and CH₂OH groups, not an all-atom or energy-minimized structure. Hydrogens are omitted. Up/down means relative to the ring faces, not the screen after camera rotation. No ring-flip pathway, energy, or solution populations are calculated. <a href="https://openstax.org/books/organic-chemistry/pages/25-5-cyclic-structures-of-monosaccharides-anomers" target="_blank" rel="noopener noreferrer">Reference: OpenStax, cyclic sugars and anomers</a>.</small></p>';
 document.body.append(dialog);dialog.showModal();{let dn=false;const out=e=>{const r=dialog.getBoundingClientRect();return e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom};dialog.addEventListener('mousedown',e=>{dn=out(e)});dialog.addEventListener('click',e=>{if(dn&&out(e))dialog.close();dn=false});}const q=s=>dialog.querySelector(s);let flipped=false,renderer,frame,group;const scene=T?new T.Scene():null;let camera,controls;
 const disposeGroup=()=>{if(!group)return;scene.remove(group);group.traverse(o=>{o.geometry?.dispose();o.material?.map?.dispose();o.material?.dispose();});};
 q('.chclose').onclick=()=>dialog.close();dialog.addEventListener('close',()=>{cancelAnimationFrame(frame);disposeGroup();controls?.dispose();renderer?.dispose();renderer?.forceContextLoss();dialog.remove();launch.focus();},{once:true});
 if(T?.OrbitControls){try{renderer=new T.WebGLRenderer({antialias:true});renderer.setPixelRatio(Math.min(devicePixelRatio,2));q('.chhost').append(renderer.domElement);scene.background=new T.Color(0x10231e);scene.add(new T.HemisphereLight(0xffffff,0x445544,1.5));const light=new T.DirectionalLight(0xffffff,1);light.position.set(4,-5,8);scene.add(light);camera=new T.PerspectiveCamera(42,1,.01,100);camera.up.set(0,0,1);controls=new T.OrbitControls(camera,renderer.domElement);}catch(e){renderer?.dispose();renderer=null;q('.chhost').textContent='3D view unavailable. Haworth and the position table remain available.';}}else q('.chhost').textContent='3D engine unavailable. Haworth and the position table remain available.';
 function resetCamera(){if(!camera)return;camera.position.set(4,-6,4);controls.target.set(0,0,0);controls.update();}resetCamera();
 function draw3D(m){if(!renderer)return;disposeGroup();group=new T.Group();scene.add(group);
  const sphere=(p,color,r)=>{const mesh=new T.Mesh(new T.SphereGeometry(r,18,12),new T.MeshStandardMaterial({color,roughness:.5}));mesh.position.set(p.x,p.y,p.z);group.add(mesh);};
  const line=(p,q,color)=>{const a=new T.Vector3(p.x,p.y,p.z),b=new T.Vector3(q.x,q.y,q.z),d=b.clone().sub(a),mesh=new T.Mesh(new T.CylinderGeometry(.065,.065,d.length(),10),new T.MeshStandardMaterial({color}));mesh.position.copy(a.add(b).multiplyScalar(.5));mesh.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),d.normalize());group.add(mesh);};
  const label=(p,text,color)=>{const c=document.createElement('canvas');c.width=256;c.height=64;const ctx=c.getContext('2d');ctx.font='bold 34px Arial';ctx.fillStyle=color;ctx.textAlign='center';ctx.fillText(text,128,40);const texture=new T.CanvasTexture(c),sprite=new T.Sprite(new T.SpriteMaterial({map:texture,depthTest:false}));sprite.position.set(p.x,p.y,p.z+.23);sprite.scale.set(1.65,.42,1);group.add(sprite);};
  m.ring.forEach((p,i)=>{sphere(p,i===5?0xff4444:0x9aadaa,.18);line(p,m.ring[(i+1)%6],0xaabbb6);label(p,p.label,'#ffffff');});m.groups.forEach((p,i)=>{const color=p.position==='axial'?0xffbd76:0x6de5ee;line(m.ring[i],p,color);sphere(p,color,.15);label(p,p.label+' · '+p.position.slice(0,2),p.position==='axial'?'#ffbd76':'#6de5ee');});
 }
 function drawHaworth(m){const svg=q('.chsvg');svg.replaceChildren();svg.setAttribute('aria-label',m.name+' Haworth projection');const add=(tag,attrs,text)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;svg.append(e);};const pts=[[325,170],[270,230],[145,230],[80,170],[145,110],[270,110]];
  add('polygon',{points:pts.map(p=>p.join(',')).join(' '),fill:'none',stroke:'#c1d9cd','stroke-width':3});add('line',{x1:145,y1:230,x2:270,y2:230,stroke:'#c1d9cd','stroke-width':6});
  const text=(x,y,t,color='#edf6f2')=>add('text',{x,y,fill:color,'font-size':17,'font-family':'Arial','text-anchor':'middle'},t);
  pts.forEach(([x,y],i)=>{add('circle',{cx:x,cy:y,r:i===5?14:4,fill:'#10231e'});text(x+(i===5?0:20),y+(i===5?6:5),i===5?'O':'C'+(i+1),i===5?'#ff7777':'#edf6f2');});
  m.groups.forEach((g,i)=>{const [x,y]=pts[i],dy=g.direction==='up'?-50:50;add('line',{x1:x,y1:y,x2:x,y2:y+dy,stroke:'#6de5ee','stroke-width':2});text(x,y+dy+(dy<0?-8:20),g.label,'#6de5ee');});
  text(210,320,'Up/down configuration stays fixed during a flip.','#a6c7b8');
 }
 function update(){const m=chairModel(q('.chsugar').value,q('.chanomer').value,flipped);q('.chname').textContent=m.name;q('.chchair').textContent='3D chair · '+m.chair;q('.chflip').textContent='Flip to '+(flipped?'⁴C₁':'¹C₄');q('tbody').replaceChildren();m.groups.forEach(g=>{const row=document.createElement('tr');for(const v of ['C'+g.carbon,g.label,g.direction,g.position]){const cell=document.createElement('td');cell.textContent=v;row.append(cell);}q('tbody').append(row);});const axial=m.groups.filter(g=>g.position==='axial').length;q('.chstatus').textContent=m.chair+': '+axial+' axial and '+(5-axial)+' equatorial non-hydrogen substituents. Orange = axial; cyan = equatorial. Drag the 3D view to rotate.';drawHaworth(m);draw3D(m);}
 q('.chsugar').onchange=()=>{flipped=false;update();};q('.chanomer').onchange=update;q('.chflip').onclick=()=>{flipped=!flipped;update();};q('.chreset').onclick=resetCamera;update();if(renderer){function animate(){const c=renderer.domElement,w=c.clientWidth,h=c.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();controls.update();renderer.render(scene,camera);frame=requestAnimationFrame(animate);}animate();}
};

return {};})();
__M['electron-lab']=(()=>{
const {electronModel,electronCounts}=__M['electron-core'];
const style=document.createElement('style');style.textContent=`#electronLab{background:var(--sf,#17231f);color:var(--tx,#fff);border:1px solid var(--bd,#777);border-radius:16px;width:min(1000px,94vw);max-height:90vh;overflow:auto;padding:24px}#electronLab::backdrop{background:#000a}#electronLab .elviews{display:flex;gap:22px}#electronLab .elhost{flex:1.2;min-width:0}#electronLab .elinfo{flex:1;min-width:0}#electronLab canvas{width:100%;height:390px;display:block;background:#10231e;border-radius:12px;touch-action:none}#electronLab button,#electronLab select{margin:5px;width:auto}#electronLab .elcontrols{display:flex;gap:16px;flex-wrap:wrap;align-items:center;margin:12px 0}#electronLab p{line-height:1.5}#electronLab dl{display:grid;grid-template-columns:1.3fr 1fr;gap:12px}#electronLab dd{margin:0;font-weight:bold}#electronLab dt{color:var(--mu,#b8c9c2)}#electronLab .elstatus{padding:12px;border:1px solid var(--bd,#777);border-radius:8px}@media(max-width:650px){#electronLab .elviews{flex-direction:column}}`;document.head.append(style);
const launch=document.createElement('button');launch.className='ib';launch.setAttribute('aria-label','Electron lab');launch.innerHTML='<svg viewBox="0 0 24 24" aria-hidden="true"><circle cx="12" cy="12" r="1.7"/><ellipse cx="12" cy="12" rx="9" ry="3.6"/><ellipse cx="12" cy="12" rx="9" ry="3.6" transform="rotate(60 12 12)"/><ellipse cx="12" cy="12" rx="9" ry="3.6" transform="rotate(120 12 12)"/></svg>';launch.title='Bonding electron pairs, lone pairs and molecular geometry';document.querySelector('#inf').before(launch);
launch.onclick=()=>{
 const T=window.THREE,dialog=document.createElement('dialog');dialog.id='electronLab';dialog.setAttribute('aria-label','Electron pairs and molecular geometry');dialog.innerHTML='<button class="elclose" style="float:right;position:sticky;top:0;z-index:5" aria-label="Close electron lab">×</button><h2>Electrons & molecular geometry</h2><details class="td"><summary><b>About this lab</b></summary><p><b>What:</b> water, ammonia, methane and CO2 with bonding pairs, lone pairs and electron counts, explaining electron-region geometry versus molecular shape (VSEPR). <b>Limits:</b> Lewis-style dots and lobe guides, not orbitals or electron density; four neutral teaching models only; lone pairs are not assigned automatically for imported structures.</p></details><p>Neutral molecule teaching models</p><label>Molecule <select class="elpreset" aria-label="Electron model"><option value="water">Water (H₂O)</option><option value="ammonia">Ammonia (NH₃)</option><option value="methane">Methane (CH₄)</option><option value="co2">Carbon dioxide (CO₂)</option></select></label><button class="elreset">Reset camera</button><div class="elcontrols"><label><input type="checkbox" class="elbond" checked> Bonding pairs</label><label><input type="checkbox" class="ellone" checked> Lone pairs</label><label><input type="checkbox" class="elregions" checked> Lone-pair guides</label><label><input type="checkbox" class="ellabels" checked> Atom labels</label></div><h3 class="elname"></h3><div class="elviews"><div class="elhost"></div><section class="elinfo"><h3>Geometry at the central atom</h3><dl class="elfacts"></dl><p class="elexplain"></p></section></div><p class="elstatus" role="status" aria-live="polite"></p><p>Cyan dots = bonding electrons; gold dots = lone-pair electrons. Each pair has two dots. Drag to rotate and scroll to zoom. Turning a display off hides its markers without changing the molecule.</p><p><small>Dots are Lewis electron-counting symbols, not fixed electron positions or trajectories. Translucent lone-pair guides indicate schematic directions, not calculated orbitals or electron density. CO₂’s two pairs per double bond are counted separately; each double bond is one VSEPR region. Presets have defined bonding and charge; lone pairs are not inferred for imported molecules. <a href="https://openstax.org/books/chemistry-2e/pages/7-6-molecular-structure-and-polarity" target="_blank" rel="noopener noreferrer">Reference: OpenStax, molecular structure</a>.</small></p>';
 document.body.append(dialog);dialog.showModal();{let dn=false;const out=e=>{const r=dialog.getBoundingClientRect();return e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom};dialog.addEventListener('mousedown',e=>{dn=out(e)});dialog.addEventListener('click',e=>{if(dn&&out(e))dialog.close();dn=false});}const q=s=>dialog.querySelector(s);let model,renderer,camera,controls,frame,root,bondDots,loneDots,guides,labels;const scene=T?new T.Scene():null;
 function dispose(){if(!root)return;scene.remove(root);root.traverse(o=>{o.geometry?.dispose();o.material?.map?.dispose();o.material?.dispose();});root=null;}
 q('.elclose').onclick=()=>dialog.close();dialog.addEventListener('close',()=>{cancelAnimationFrame(frame);dispose();controls?.dispose();renderer?.dispose();renderer?.forceContextLoss();dialog.remove();launch.focus();},{once:true});
 if(T?.OrbitControls){try{renderer=new T.WebGLRenderer({antialias:true});renderer.setPixelRatio(Math.min(devicePixelRatio,2));q('.elhost').append(renderer.domElement);scene.background=new T.Color(0x10231e);scene.add(new T.HemisphereLight(0xffffff,0x445544,1.4));const light=new T.DirectionalLight(0xffffff,1);light.position.set(3,-5,6);scene.add(light);camera=new T.PerspectiveCamera(42,1,.01,100);camera.up.set(0,0,1);controls=new T.OrbitControls(camera,renderer.domElement);}catch(e){renderer?.dispose();renderer=null;}}
 if(!renderer)q('.elhost').textContent='The 3D engine is unavailable. Electron counts and geometry explanations remain available.';
 const reset=()=>{if(!camera)return;camera.position.set(2,-5,2.4);controls.target.set(0,0,0);controls.update();};q('.elreset').onclick=reset;
 function visibility(){const c=electronCounts(model);if(root){bondDots.visible=q('.elbond').checked;loneDots.visible=q('.ellone').checked;guides.visible=q('.elregions').checked&&q('.ellone').checked;labels.visible=q('.ellabels').checked;}q('.elstatus').textContent='Shown: '+(q('.elbond').checked?c.bondingPairs:0)+' bonding pairs and '+(q('.ellone').checked?c.lonePairs:0)+' lone pairs. Whole molecule: '+c.valenceElectrons+' valence electrons ('+c.bondingPairs+' bonding pairs + '+c.lonePairs+' lone pairs).';}
 for(const s of ['.elbond','.ellone','.elregions','.ellabels'])q(s).onchange=visibility;
 function draw(){if(!renderer)return;dispose();root=new T.Group();scene.add(root);bondDots=new T.Group();loneDots=new T.Group();guides=new T.Group();labels=new T.Group();root.add(bondDots,loneDots,guides,labels);
  const vec=a=>new T.Vector3(a.x,a.y,a.z),sphere=(g,p,r,color,unlit=false)=>{const m=new T.Mesh(new T.SphereGeometry(r,18,12),unlit?new T.MeshBasicMaterial({color}):new T.MeshStandardMaterial({color,roughness:.5}));m.position.copy(p);g.add(m);};
  const perp=d=>new T.Vector3().crossVectors(d,Math.abs(d.z)<.9?new T.Vector3(0,0,1):new T.Vector3(0,1,0)).normalize();
  const pair=(g,p,t,color)=>{sphere(g,p.clone().addScaledVector(t,.075),.048,color,true);sphere(g,p.clone().addScaledVector(t,-.075),.048,color,true);};
  model.atoms.forEach((a,i)=>{const p=vec(a);sphere(root,p,a.el==='H'?.16:.24,{H:0xffffff,C:0x9aadaa,N:0x466bff,O:0xff4545}[a.el]);const c=document.createElement('canvas');c.width=128;c.height=64;const ctx=c.getContext('2d');ctx.font='bold 32px Arial';ctx.fillStyle='#ffffff';ctx.textAlign='center';ctx.fillText(a.el+(i===0?' (center)':''),64,40,125);const sprite=new T.Sprite(new T.SpriteMaterial({map:new T.CanvasTexture(c),depthTest:false}));sprite.position.copy(p).add(new T.Vector3(0,0,.3));sprite.scale.set(.7,.35,1);labels.add(sprite);});
  model.bonds.forEach(([i,j,order])=>{const a=vec(model.atoms[i]),b=vec(model.atoms[j]),d=b.clone().sub(a),len=d.length(),u=d.clone().normalize(),t=perp(u);for(let k=0;k<order;k++){const offset=(k-(order-1)/2)*.16,p=a.clone().addScaledVector(t,offset),r=b.clone().addScaledVector(t,offset),mesh=new T.Mesh(new T.CylinderGeometry(.035,.035,len,10),new T.MeshStandardMaterial({color:0xaabbb6}));mesh.position.copy(p).add(r).multiplyScalar(.5);mesh.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),u);root.add(mesh);pair(bondDots,mesh.position.clone().addScaledVector(t,.09),u,0x5eeeff);}});
  model.lonePairs.forEach(lp=>{const p=vec(model.atoms[lp.atom]),u=new T.Vector3(...lp.direction).normalize(),end=p.clone().addScaledVector(u,.7);pair(loneDots,end,perp(u),0xffce63);const guide=new T.Mesh(new T.SphereGeometry(1,20,14),new T.MeshBasicMaterial({color:0xffce63,transparent:true,opacity:.16,depthWrite:false}));guide.position.copy(p).addScaledVector(u,.48);guide.scale.set(.23,.23,.48);guide.quaternion.setFromUnitVectors(new T.Vector3(0,0,1),u);guides.add(guide);});
 }
 function load(){model=electronModel(q('.elpreset').value);const c=electronCounts(model);q('.elname').textContent=model.name;q('.elfacts').replaceChildren();for(const [key,value] of [['Central atom',model.atoms[0].el],['Electron regions',c.domains],['Lone pairs at center',c.centralLonePairs],['Electron-region geometry',model.electronGeometry],['Molecular shape',model.geometry],['Bond angle (approx.)',model.angle.toFixed(1)+'°']]){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;q('.elfacts').append(dt,dd);}q('.elexplain').textContent={water:'Two bonds and two lone pairs surround oxygen. The electron regions are approximately tetrahedral, while the atoms form a bent shape.',ammonia:'Three bonds and one lone pair surround nitrogen. The atoms form a trigonal pyramid within an approximately tetrahedral arrangement of electron regions.',methane:'Four bonds surround carbon, with no lone pairs. Both electron-region geometry and molecular shape are tetrahedral.',co2:'Carbon has two double bonds and no lone pairs: two electron regions give a linear shape. Each oxygen has two lone pairs; these are included in the whole-molecule count.'}[model.id];draw();visibility();reset();}
 q('.elpreset').onchange=load;load();if(renderer){function animate(){const c=renderer.domElement,w=c.clientWidth,h=c.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();controls.update();renderer.render(scene,camera);frame=requestAnimationFrame(animate);}animate();}
};

return {};})();
}catch(err){console.error("LAURA labs failed to load",err)}})();
