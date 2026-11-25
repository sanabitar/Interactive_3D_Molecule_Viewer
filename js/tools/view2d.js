// js/tools/view2d.js: The 2D view dialog: generated structure diagram (uses Layout2D) or a flat projection of the 3D coordinates.
// ---------- 2D view ----------
function show2D(){if(!need())return;const A=mol.atoms,hOn=()=>$('#hh').checked;
 // Every atom of the view (hydrogens included); the Hydrogens option is applied while drawing so the dialog can switch it live.
 const all=sel.length?sel.slice():A.map((_,i)=>i).filter(i=>$('#ww').checked||A[i].res!=='HOH'),ids=()=>all.filter(i=>hOn()||A[i].el!=='H');
 if(all.length>5000){modal(`<b>2D view</b><p class="kv">${all.length} atoms are too many for a readable flat view. Select a smaller region (Select by… or click atoms) and try again.</p>`+CL);return}
 modal(`<b>2D view</b> <span class="kv">${sel.length?'selected atoms only ('+sel.length+')':'all atoms'}</span><div class="kv" id="d2n"></div><div class="row mt"><label>View <select id="d2m"><option value="dia">Structure diagram</option><option value="proj">Flat projection of 3D</option></select></label><label title="Draw hydrogen atoms as bonded atoms (same as the H button and the H key). Hidden: hydrogens are written on the heteroatom instead, e.g. OH or NH₂."><input type="checkbox" id="d2h"${hOn()?' checked':''}> Hydrogens</label><label><input type="checkbox" id="d2l" checked> Element labels</label><button id="d2s">Save PNG</button></div><canvas id="d2c" width="1520" height="1040" title="Click to enlarge"></canvas>`+CL);
 const drawProj=()=>{const cv=$('#d2c'),g=cv.getContext('2d'),W=cv.width/2,H=cv.height/2,css=getComputedStyle(document.documentElement),bg=css.getPropertyValue('--sf').trim()||'#fff',fg=css.getPropertyValue('--tx').trim()||'#111';
  const idl=ids(),pts=idl.map(i=>A[i]);g.setTransform(2,0,0,2,0,0);g.fillStyle=bg;g.fillRect(0,0,W,H);if(!pts.length)return;
  const P=Core.projectPlane(pts);let a0=1e9,a1=-1e9,b0=1e9,b1=-1e9;for(const p of P.xy){if(p[0]<a0)a0=p[0];if(p[0]>a1)a1=p[0];if(p[1]<b0)b0=p[1];if(p[1]>b1)b1=p[1]}
  const s=Math.min((W-90)/Math.max(a1-a0,1e-6),(H-90)/Math.max(b1-b0,1e-6),60),cx=W/2-(a0+a1)/2*s,cy=H/2+(b0+b1)/2*s,pos=new Map(idl.map((i,k)=>[i,[cx+P.xy[k][0]*s,cy-P.xy[k][1]*s,P.depth[k]]])),small=idl.length<=120;
  g.strokeStyle=fg;g.lineCap='round';g.lineWidth=small?2:1;
  mol.bonds.forEach(([i,j],k)=>{const a=pos.get(i),b=pos.get(j);if(!a||!b)return;const o=mol.orders[k]||1,dx=b[0]-a[0],dy=b[1]-a[1],L=Math.hypot(dx,dy)||1,nx=-dy/L,ny=dx/L,line=(d,dash)=>{g.setLineDash(dash?[4,3]:[]);g.beginPath();g.moveTo(a[0]+nx*d,a[1]+ny*d);g.lineTo(b[0]+nx*d,b[1]+ny*d);g.stroke()},w=small?3.2:1.6;
   if(o===2){line(-w);line(w)}else if(o===3){line(0);line(-w*1.4);line(w*1.4)}else if(o===1.5){line(-w*.6);line(w*.9,true)}else line(0)});g.setLineDash([]);
  [...idl].sort((p,q)=>pos.get(p)[2]-pos.get(q)[2]).forEach(i=>{const a=A[i],q=pos.get(i),r=small?10:2.6,hex='#'+cpk(a.el).toString(16).padStart(6,'0');g.beginPath();g.arc(q[0],q[1],r,0,7);g.fillStyle=hex;g.fill();g.lineWidth=1;g.strokeStyle=fg;g.stroke();
   if(small&&$('#d2l').checked){const c=cpk(a.el),lum=((c>>16&255)*.299+(c>>8&255)*.587+(c&255)*.114);g.fillStyle=lum>150?'#111':'#fff';g.font='bold 11px system-ui';g.textAlign='center';g.textBaseline='middle';g.fillText(a.el,q[0],q[1]+.5)}});
  g.fillStyle=fg;g.font='11px system-ui';g.textAlign='left';g.textBaseline='alphabetic';g.fillText(dn(mol.entry.name)+' · '+idl.length+' atoms · 2D projection (not a structure diagram)',10,H-10)};
 // Bond orders / charges as drawn: delocalised carboxylate shown as C=O + C–O⁻, 4-coordinate N as N⁺ (Layout2D.drawnOrdersAndCharges).
 const dr=Layout2D.drawnOrdersAndCharges(A,mol.bonds,mol.orders),lay=()=>Layout2D.layout2D(A,mol.bonds,dr.orders,ids(),{explicitH:true});
 let LY=lay();
 if(!LY.ok){$('#d2m').value='proj';$('#d2m').options[0].disabled=true}
 const drawDia=()=>{const cv=$('#d2c'),g=cv.getContext('2d'),W=cv.width/2,H=cv.height/2,css=getComputedStyle(document.documentElement),bg=css.getPropertyValue('--sf').trim()||'#fff',fg=css.getPropertyValue('--tx').trim()||'#111';
  g.setTransform(2,0,0,2,0,0);g.fillStyle=bg;g.fillRect(0,0,W,H);const pos=LY.pos;let a0=1e9,a1=-1e9,b0=1e9,b1=-1e9;for(const p of pos.values()){a0=Math.min(a0,p[0]);a1=Math.max(a1,p[0]);b0=Math.min(b0,p[1]);b1=Math.max(b1,p[1])}
  const s=Math.min((W-120)/Math.max(a1-a0,1e-6),(H-120)/Math.max(b1-b0,1e-6),52),cx=W/2-(a0+a1)/2*s,cy=H/2+(b0+b1)/2*s,T=i=>{const p=pos.get(i);return[cx+p[0]*s,cy-p[1]*s]},showC=$('#d2l').checked;
  g.font='bold '+Math.round(Math.max(13,Math.min(20,s*.4)))+'px system-ui';g.textBaseline='middle';g.textAlign='left';g.lineJoin='round';
  const nbr=new Map();for(const[a,b]of mol.bonds)if(pos.has(a)&&pos.has(b)){(nbr.get(a)||nbr.set(a,[]).get(a)).push(b);(nbr.get(b)||nbr.set(b,[]).get(b)).push(a)}
  // What is written at each atom: {pre,sym,post} or nothing (bare skeletal vertex). Hydrogens that are not drawn as atoms are written on the label (OH, NH₂, NH₃⁺).
  const lab=new Map();let nH=0,nHeavy=0;
  for(const i of pos.keys()){const e=A[i].el;if(e==='H'){nH++;lab.set(i,{pre:'',sym:'H',post:''});continue}nHeavy++;
   const ch=dr.charge.get(i)||0,ns=nbr.get(i)||[],heavyN=ns.filter(j=>A[j].el!=='H').length;if(e==='C'&&!ch&&!showC&&heavyN>0)continue;
   const n=(LY.hcount.get(i)||0)-(LY.hdrawn.get(i)||0),mx=ns.length?ns.reduce((t,j)=>t+pos.get(j)[0],0)/ns.length-pos.get(i)[0]:0; // bonds leaving to the right: write HO / H₂N so the symbol touches the bond
   lab.set(i,Layout2D.atomLabelParts(e,n,ch,mx>.15))}
  const r0=Math.max(8,s*.27),w=Math.min(8,Math.max(3,s*.12));g.lineCap='round';g.lineWidth=Math.max(1.6,s*.045);g.strokeStyle=fg;
  mol.bonds.forEach(([i,j],k)=>{if(!pos.has(i)||!pos.has(j))return;const o=dr.orders[k]||1,A0=T(i),B=T(j),dx=B[0]-A0[0],dy=B[1]-A0[1],L=Math.hypot(dx,dy)||1,ux=dx/L,uy=dy/L;let nx=-uy,ny=ux;
   const ti=lab.has(i)?r0:0,tj=lab.has(j)?r0:0,rg=LY.rings.filter(r=>r.atoms.includes(i)&&r.atoms.includes(j)).sort((p,q)=>p.atoms.length-q.atoms.length)[0];let inner=0;
   if(rg){const c=[cx+rg.center[0]*s,cy-rg.center[1]*s];inner=((c[0]-(A0[0]+B[0])/2)*nx+(c[1]-(A0[1]+B[1])/2)*ny)>0?1:-1}
   const line=(d,dash,sh=0)=>{const p0=[A0[0]+ux*(ti+sh)+nx*d,A0[1]+uy*(ti+sh)+ny*d],p1=[B[0]-ux*(tj+sh)+nx*d,B[1]-uy*(tj+sh)+ny*d];g.setLineDash(dash?[5,4]:[]);g.beginPath();g.moveTo(...p0);g.lineTo(...p1);g.stroke()};
   if(o===2){if(inner){line(0);line(inner*w*1.5,false,w*1.2)}else{line(-w/2);line(w/2)}}
   else if(o===3){line(0);line(-w);line(w)}
   else if(o===1.5){if(inner){line(0);line(inner*w*1.5,true,w*1.2)}else{line(-w/2);line(w/2,true)}}
   else line(0)});g.setLineDash([]);
  for(const[i,p]of lab){const q=T(i),e=A[i].el,col=(e==='C'||e==='H')?fg:'#'+cpk(e).toString(16).padStart(6,'0'),t=p.pre+p.sym+p.post,x0=q[0]-g.measureText(p.pre).width-g.measureText(p.sym).width/2; // the element symbol is centred on the atom
   g.lineWidth=6;g.strokeStyle=bg;g.strokeText(t,x0,q[1]+.5);g.fillStyle=col;g.fillText(t,x0,q[1]+.5)}
  g.fillStyle=fg;g.font='11px system-ui';g.textAlign='left';g.textBaseline='alphabetic';g.fillText(dn(mol.entry.name)+' · '+nHeavy+' heavy atoms'+(nH?' + '+nH+' H':hOn()&&!A.some(a=>a.el==='H')?' · this structure has no hydrogen atoms':'')+' · generated 2D diagram: idealised angles, no stereochemistry'+(LY.overlaps?' · '+LY.overlaps+' overlapping atom pair(s)':''),10,H-10)};
 const draw=()=>{const dia=$('#d2m').value==='dia';if(dia&&LY.ok)LY=lay();$('#d2n').innerHTML=dia?'A generated <b>structure diagram</b> with idealised bond angles: rings as regular polygons, chains zigzagged, substituents pointing outward. It does <b>not</b> show 3D geometry or stereochemistry. <b>Hydrogens</b> on: H atoms are drawn as real bonded atoms. Off: they are written on the heteroatom instead (OH, NH₂, NH₃⁺). <b>Element labels</b> also labels carbons (off by default, as in skeletal formulas).':(LY.ok?'':'<b>Structure diagram unavailable: '+LY.reason+'.</b> ')+'Orthographic projection of the 3D coordinates onto their best-fit plane. It is <b>not</b> a chemical structure diagram: planar molecules keep their exact shape, non-planar ones overlap. <b>Hydrogens</b> draws H atoms.';dia?drawDia():drawProj()};
 redraw2D=draw;$('#d2h').onchange=()=>setH($('#d2h').checked);$('#d2l').onchange=draw;$('#d2m').onchange=()=>{$('#d2l').checked=$('#d2m').value!=='dia';draw()};$('#d2s').onclick=()=>$('#d2c').toBlob(b=>b&&save(slug()+'_2D.png',b));$('#d2l').checked=$('#d2m').value!=='dia';draw()}
