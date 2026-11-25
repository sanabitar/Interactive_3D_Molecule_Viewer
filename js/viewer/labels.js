// js/viewer/labels.js: Label layer: element symbols on atom spheres (debounced so they do not flicker), atom-name labels, chain and contact tags,
// and the selection refresh (drawSel, refreshSel, refreshLabels, updateLabels).
// ---------- label layer (DOM, decluttered) ----------
const LBL=$('#lbl'),ELB=$('#elb');let lbItems=[],lbKey='',lbDirty=true;
// Element symbols on the atom spheres. An atom gets its symbol only if its sphere is big enough on screen to hold it (so zooming in reveals more),
// atoms are visited nearest first, an atom hidden behind a nearer sphere gets none, and symbols that would overlap an already placed one are dropped.
// Stability while rotating/zooming: every label is a DOM node KEYED TO ITS ATOM and positioned/scaled with `transform` only (no layout work per frame,
// smooth sub-pixel size instead of whole-pixel font steps); an atom that was labelled last frame keeps its label a little longer (smaller on/off threshold,
// higher priority, more tolerance for being covered), so labels no longer flicker on and off as the camera moves; labels fade out instead of vanishing.
const EL_ON_PX=9,EL_OFF_PX=7,EL_CAP=250,EL_MAX_ATOMS=8000,EL_BASE=14,EL_SHOW=4,EL_HIDE=12,elNodes=new Map(),elHit=new Map(),elMiss=new Map(),elCol={};let elShown=new Set(),elMol=null;
const elInk=el=>elCol[el]||(elCol[el]=(c=>((c>>16&255)*.299+(c>>8&255)*.587+(c&255)*.114)>128?['#111','0 0 3px rgba(255,255,255,.9)']:['#fff','0 0 3px rgba(0,0,0,.85)'])(cpk(el)));
function elementLabels(w,h){if(elMol!==mol){ELB.textContent='';elNodes.clear();elShown=new Set();elHit.clear();elMiss.clear();elMol=mol}
 const m=aMesh&&aMesh.userData.map,on=$('#el').checked&&mol&&ro&&ro.st==='bs'&&m&&m.length<=EL_MAX_ATOMS,now=performance.now(),items=new Map();let acc=[],P=[];
 if(on){const A=mol.atoms,v=new THREE.Vector3(),ppu=h/(2*Math.tan(cam.fov*Math.PI/360)),cand=[];
  for(const i of m){const a=A[i],was=elShown.has(i),dx=a.x-cen.x-cam.position.x,dy=a.y-cen.y-cam.position.y,dz=a.z-cen.z-cam.position.z,d=Math.hypot(dx,dy,dz)||1e-6,px=ro.R(a)*ppu/d;if(px<(was?EL_OFF_PX:EL_ON_PX))continue;
   v.set(a.x-cen.x,a.y-cen.y,a.z-cen.z);v.project(cam);if(!(v.z>-1&&v.z<1))continue;const x=(v.x*.5+.5)*w,y=(-v.y*.5+.5)*h;if(x<-px||x>w+px||y<-px||y>h+px)continue;cand.push({i,x,y,d,px,was})}
  cand.sort((p,q)=>p.d-q.d||p.i-q.i);const seen=[];
  for(const c of cand){const fs=Math.max(9,Math.min(26,c.px*1.3)),text=A[c.i].el;items.set(c.i,{id:c.i,x:c.x,y:c.y,w:text.length*fs*.66,h:fs*1.1,prio:(c.was?0:1e6)+c.d,fs,text,el:text}); // geometry of every candidate, wanted or not
   if(P.length>=EL_CAP*2)continue;const k=c.was?.55:.85;if(seen.some(s=>(s.x-c.x)**2+(s.y-c.y)**2<(s.px*k)**2)){seen.push(c);continue}seen.push(c); // centre hidden behind a nearer sphere
   P.push(items.get(c.i))} // sticky: labels shown last frame outrank new ones
  acc=Core.declutter(P,{width:w,height:h,cap:EL_CAP,pad:1})}
 // Debounce: an atom must be wanted for EL_SHOW frames in a row before its label appears, and unwanted for EL_HIDE frames before it goes (hard reasons -
 // off screen, sphere too small, label turned off - apply at once). Brief overlaps/occlusions while rotating therefore do not make labels twinkle.
 const want=new Set(acc.map(k=>P[k].id)),next=new Set();let pending=false;
 for(const id of want){if(elShown.has(id)){next.add(id);elMiss.delete(id)}else{const n=(elHit.get(id)||0)+1;if(n>=EL_SHOW){next.add(id);elHit.delete(id)}else{elHit.set(id,n);pending=true}}}
 for(const id of [...elHit.keys()])if(!want.has(id))elHit.delete(id);
 for(const id of elShown){if(want.has(id))continue;if(!items.has(id)){elMiss.delete(id);continue}const n=(elMiss.get(id)||0)+1;if(n<EL_HIDE){elMiss.set(id,n);next.add(id);pending=true}else elMiss.delete(id)}
 elShown=next;if(pending)lbDirty=true;
 for(const id of elShown){const q=items.get(id),ink=elInk(q.el);let n=elNodes.get(id);if(!n){n=document.createElement('div');n.className='el';n.textContent=q.text;ELB.appendChild(n);elNodes.set(id,n)}
  n.style.opacity='1';n._t=0;n.dataset.x=q.x.toFixed(1);n.dataset.y=q.y.toFixed(1);n.dataset.f=q.fs.toFixed(1);
  n.style.transform='translate('+q.x.toFixed(1)+'px,'+q.y.toFixed(1)+'px) translate(-50%,-50%) scale('+(q.fs/EL_BASE).toFixed(3)+')';if(n._c!==ink[0]){n.style.color=ink[0];n.style.textShadow=ink[1];n._c=ink[0]}}
 for(const[id,n]of elNodes){if(elShown.has(id))continue;if(!n._t){n._t=now;n.style.opacity='0'}else if(now-n._t>400){n.remove();elNodes.delete(id)}if(elNodes.has(id))lbDirty=true} // keep running until faded labels are gone
}
function drawSel(){drawSel0();refreshLabels()}
function refreshSel(){if($('#so').checked&&!sel.length){$('#so').checked=false;build();fit();return}if($('#so').checked)build();else drawSel()}
function refreshLabels(){const it=[];if(mol){if(ro.st==='tb')chainLab.forEach(c=>it.push({p:c.p,text:c.text,cls:'chain',prio:0}));cLab.forEach(c=>it.push({p:c.p,text:c.text,cls:'contact',prio:c.prio??1,off:c.off||0}));
 if($('#lb').checked){const n=sel.length,step=n>4000?Math.ceil(n/4000):1;for(let k=0;k<n;k+=step){const a=mol.atoms[sel[k]];it.push({p:[a.x-cen.x,a.y-cen.y,a.z-cen.z],text:Core.labelText(a),cls:'',prio:2})}}}
 cmpTags.forEach(c=>it.push({p:c.p,text:c.text,cls:'tag',prio:-3,off:c.off||0}));lbItems=it;lbDirty=true}
function updateLabels(){if(!ren)return;const w=vp.clientWidth,h=vp.clientHeight,T=ctl?ctl.target:{x:0,y:0,z:0},key=[cam.position.x,cam.position.y,cam.position.z,T.x,T.y,T.z,w,h,cam.fov].map(v=>typeof v==='number'?v.toFixed(3):v).join();
 if(key===lbKey&&!lbDirty)return;lbKey=key;lbDirty=false;elementLabels(w,h);const v=new THREE.Vector3(),P=[];
 for(const it of lbItems){v.set(it.p[0],it.p[1],it.p[2]);v.project(cam);if(!(v.z>-1&&v.z<1))continue;P.push({x:(v.x*.5+.5)*w,y:(-v.y*.5+.5)*h-(it.off||0),w:it.text.length*6.6+12,h:18,prio:it.prio,it})}
 const acc=Core.declutter(P,{width:w,height:h,cap:80,pad:3});while(LBL.children.length<acc.length){const d=document.createElement('div');d.className='lab';LBL.appendChild(d)}
 [...LBL.children].forEach((d,i)=>{if(i<acc.length){const q=P[acc[i]];d.style.display='block';d.textContent=q.it.text;d.className='lab '+q.it.cls;d.style.left=q.x+'px';d.style.top=q.y+'px'}else d.style.display='none'})}
function makeContactLabels(k){cLab=[];if(!lastC)return;const idx=new Map();mol.atoms.forEach((a,i)=>{if(!a.hetero){const key=Core.residueKey(a)+'|'+a.name;if(!idx.has(key))idx.set(key,i)}});
 lastC.contacts.forEach(c=>{const i=idx.get(c.residue+'|'+c.proteinAtom);if(i===undefined)return;const a=mol.atoms[i];cLab.push({p:[a.x-cen.x,a.y-cen.y,a.z-cen.z],text:`${c.resName}${c.residue.slice(c.residue.indexOf(':')+1)} · ${c.minDistance.toFixed(1)} Å`})});
 const lg=mol.atoms.filter(a=>a.hetero&&Core.residueKey(a)===k);if(lg.length){const m=lg.reduce((s,a)=>[s[0]+a.x,s[1]+a.y,s[2]+a.z],[0,0,0]).map(x=>x/lg.length);cLab.push({p:[m[0]-cen.x,m[1]-cen.y,m[2]-cen.z],text:`${lg[0].resName} (ligand)`,prio:-2,off:18})}}
