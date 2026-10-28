// Structure-diagram layout: turns a molecular graph into clean 2D coordinates (bond length = 1).
// Rings become regular polygons, fused rings share an edge, chains zigzag at 120 degrees, substituents point outward.
// Not supported (returns {ok:false}): bridged or spiro ring systems and more than maxHeavy heavy atoms. Stereochemistry is not drawn.
const TAU=Math.PI*2;
export function layout2D(atoms,bonds,orders,ids,opts={}){
 const maxHeavy=opts.maxHeavy??80,el=i=>atoms[i].el||atoms[i].element,set=new Set(ids.filter(i=>el(i)!=='H'));
 if(!set.size)return{ok:false,reason:'no heavy atoms'};
 if(set.size>maxHeavy)return{ok:false,reason:'more than '+maxHeavy+' heavy atoms'};
 const adj=new Map([...set].map(i=>[i,[]])),ord=new Map(),hc=new Map([...set].map(i=>[i,0])),ek=(a,b)=>a<b?a+'-'+b:b+'-'+a;
 bonds.forEach(([a,b],k)=>{const o=orders&&orders[k]||1;
  if(set.has(a)&&set.has(b)){adj.get(a).push(b);adj.get(b).push(a);ord.set(ek(a,b),o)}
  else if(set.has(a)&&el(b)==='H')hc.set(a,hc.get(a)+1);else if(set.has(b)&&el(a)==='H')hc.set(b,hc.get(b)+1)});
 const comps=[],seenA=new Set();
 for(const s of set){if(seenA.has(s))continue;const c=[s];seenA.add(s);for(let h=0;h<c.length;h++)for(const y of adj.get(c[h]))if(!seenA.has(y)){seenA.add(y);c.push(y)}comps.push(c)}
 const pos=new Map(),ringsOut=[];let offX=0;
 for(const nodes of comps){const r=layoutComponent(nodes);if(!r.ok)return r;
  let lo=1e9,hi=-1e9;for(const i of nodes){const x=r.pos.get(i)[0];lo=Math.min(lo,x);hi=Math.max(hi,x)}
  const dx=offX-lo;for(const i of nodes){const p=r.pos.get(i);pos.set(i,[p[0]+dx,p[1]])}
  for(const g of r.rings)ringsOut.push({atoms:g.atoms,center:[g.center[0]+dx,g.center[1]]});offX+=hi-lo+1.6}
 let overlaps=0;const ks=[...pos.keys()];for(let i=0;i<ks.length;i++)for(let j=i+1;j<ks.length;j++){const a=pos.get(ks[i]),b=pos.get(ks[j]);if(Math.hypot(a[0]-b[0],a[1]-b[1])<.5)overlaps++}
 return{ok:true,pos,rings:ringsOut,hcount:hc,overlaps,components:comps.length};
 function findRings(nodes){const seen=new Set(),all=[];
  for(const u of nodes)for(const v of adj.get(u)){if(u>v)continue;const prev=new Map([[u,-1]]),q=[u];let found=false;
   for(let h=0;h<q.length&&!found;h++){const x=q[h];for(const y of adj.get(x)){if(x===u&&y===v)continue;if(prev.has(y))continue;prev.set(y,x);if(y===v){found=true;break}q.push(y)}}
   if(!found)continue;const path=[];for(let x=v;x!==-1;x=prev.get(x))path.push(x);
   const key=[...path].sort((a,b)=>a-b).join(',');if(seen.has(key))continue;seen.add(key);all.push(path)}
  all.sort((a,b)=>a.length-b.length);const cov=new Set(),rings=[];
  for(const r of all){const es=r.map((x,i)=>ek(x,r[(i+1)%r.length]));if(es.some(e=>!cov.has(e))){rings.push(r);es.forEach(e=>cov.add(e))}}
  return rings}
 function layoutComponent(nodes){
  const rings=findRings(nodes),P=new Map(),centers=new Map(),sysOf=[],systems=[];
  const par=rings.map((_,i)=>i),find=i=>par[i]===i?i:(par[i]=find(par[i]));
  for(let i=0;i<rings.length;i++)for(let j=i+1;j<rings.length;j++){const sh=rings[i].filter(x=>rings[j].includes(x));
   if(!sh.length)continue;if(sh.length!==2||!adj.get(sh[0]).includes(sh[1]))return{ok:false,reason:'bridged or spiro ring system'};par[find(i)]=find(j)}
  rings.forEach((_,i)=>{const r=find(i);if(!systems[r])systems[r]={rings:[],placed:false};systems[r].rings.push(i);sysOf[i]=systems[r]});
  const atomSys=new Map();rings.forEach((r,i)=>r.forEach(a=>atomSys.set(a,sysOf[i])));
  const R=n=>1/(2*Math.sin(Math.PI/n)),queue=[];
  const put=(a,x,y)=>{if(!P.has(a)){P.set(a,[x,y]);queue.push(a)}};
  function placeSystem(sys,anchor,ax,ay,out){
   const r0=sys.rings.find(i=>rings[i].includes(anchor)),ring=rings[r0],n=ring.length,rad=R(n);let cx=0,cy=0,th=0;
   if(out!==null){cx=ax+rad*Math.cos(out);cy=ay+rad*Math.sin(out);th=out+Math.PI}
   const ia=ring.indexOf(anchor);for(let k=0;k<n;k++){const a=ring[(ia+k)%n],t=th+k*TAU/n;put(a,cx+rad*Math.cos(t),cy+rad*Math.sin(t))}
   centers.set(r0,[cx,cy]);sys.placed=true;let more=true;
   while(more){more=false;for(const i of sys.rings){if(centers.has(i))continue;const ri=rings[i],m=ri.length;
    for(const j of sys.rings){if(!centers.has(j))continue;const sh=ri.filter(x=>rings[j].includes(x));if(sh.length!==2)continue;
     const [a,b]=sh,pa=P.get(a),pb=P.get(b),pc=centers.get(j),mx=(pa[0]+pb[0])/2,my=(pa[1]+pb[1])/2;
     let nx=-(pb[1]-pa[1]),ny=pb[0]-pa[0];const L=Math.hypot(nx,ny);nx/=L;ny/=L;if(nx*(mx-pc[0])+ny*(my-pc[1])<0){nx=-nx;ny=-ny}
     const ap=1/(2*Math.tan(Math.PI/m)),c2=[mx+nx*ap,my+ny*ap],rr=R(m),sg=((pa[0]-c2[0])*(pb[1]-c2[1])-(pa[1]-c2[1])*(pb[0]-c2[0]))>0?1:-1;
     const ia2=ri.indexOf(a),ib2=ri.indexOf(b),s=((ia2+1)%m===ib2)?1:-1,ta=Math.atan2(pa[1]-c2[1],pa[0]-c2[0]);
     for(let k=0;k<m;k++){const at=ri[((ia2+s*k)%m+m)%m],t=ta+sg*k*TAU/m;put(at,c2[0]+rr*Math.cos(t),c2[1]+rr*Math.sin(t))}
     centers.set(i,c2);more=true;break}}}
  }
  const turn=new Map();
  if(rings.length)placeSystem(sysOf[0],rings[0][0],0,0,null);
  else{const far=(s)=>{const d=new Map([[s,0]]),q=[s];for(let h=0;h<q.length;h++)for(const y of adj.get(q[h]))if(!d.has(y)){d.set(y,d.get(q[h])+1);q.push(y)}return q[q.length-1]};put(far(far(nodes[0])),0,0)}
  for(let h=0;h<queue.length;h++){const a=queue[h],pa=P.get(a),U=adj.get(a).filter(x=>!P.has(x)),Pn=adj.get(a).filter(x=>P.has(x));if(!U.length)continue;
   const k=U.length,ang=Pn.map(x=>Math.atan2(P.get(x)[1]-pa[1],P.get(x)[0]-pa[0])).sort((p,q)=>p-q);let out=[],sgn=0;
   if(!Pn.length)out=k===1?[0]:k===2?[0,TAU/3]:U.map((_,i)=>i*TAU/k);
   else if(Pn.length===1){const f=ang[0],sum=adj.get(a).reduce((s,y)=>s+ord.get(ek(a,y)),0);
    if(k===1&&adj.get(a).length===2&&sum>=4)out=[f+Math.PI];
    else if(k===1){sgn=-(turn.get(a)||-1);out=[f+sgn*TAU/3]}
    else if(k===2)out=[f+TAU/3,f-TAU/3];else out=U.map((_,i)=>f+(i+1)*TAU/(k+1))}
   else{let bi=0,bg=-1;ang.forEach((t,i)=>{const g=((ang[(i+1)%ang.length]-t)%TAU+TAU)%TAU||TAU;if(g>bg){bg=g;bi=i}});out=U.map((_,i)=>ang[bi]+bg*(i+1)/(k+1))}
   U.forEach((u,i)=>{if(P.has(u))return;const t=out[i],x=pa[0]+Math.cos(t),y=pa[1]+Math.sin(t);
    if(atomSys.has(u)&&!atomSys.get(u).placed)placeSystem(atomSys.get(u),u,x,y,t);else{put(u,x,y);if(sgn)turn.set(u,sgn)}})}
  if(nodes.some(i=>!P.has(i)))return{ok:false,reason:'layout incomplete'};
  return{ok:true,pos:P,rings:rings.map((r,i)=>({atoms:r,center:centers.get(i)||[0,0]}))};
 }
}
