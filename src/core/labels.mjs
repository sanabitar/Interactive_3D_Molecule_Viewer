// Screen-space label placement. Greedy by priority; a label is dropped if it overlaps one already placed or lies off-screen.
export function declutter(items,{width=1000,height=700,pad=3,cap=80}={}){const order=items.map((it,i)=>i).sort((a,b)=>((items[a].prio??0)-(items[b].prio??0))||a-b),placed=[],grid=new Map(),CELL=80;
 const cells=r=>{const out=[];for(let gx=Math.floor(r[0]/CELL);gx<=Math.floor(r[2]/CELL);gx++)for(let gy=Math.floor(r[1]/CELL);gy<=Math.floor(r[3]/CELL);gy++)out.push(gx+','+gy);return out};
 for(const i of order){if(placed.length>=cap)break;const it=items[i];if(!(it.x>=0&&it.x<=width&&it.y>=0&&it.y<=height))continue;
  const r=[it.x-it.w/2-pad,it.y-it.h/2-pad,it.x+it.w/2+pad,it.y+it.h/2+pad],cs=cells(r);let hit=false;
  for(const k of cs){const L=grid.get(k);if(L&&L.some(q=>!(r[2]<=q[0]||r[0]>=q[2]||r[3]<=q[1]||r[1]>=q[3]))){hit=true;break}}
  if(hit)continue;placed.push(i);for(const k of cs){if(!grid.has(k))grid.set(k,[]);grid.get(k).push(r)}}
 return placed}
// Label text: real atom name (+ residue) for imported/protein atoms; element symbol for generated molecules whose "names" are just running numbers (C1, O2...).
export function labelText(a){const el=a.element||a.el||'';if(a.autoName)return el;const res=a.resName&&a.resSeq!==undefined&&a.resSeq!==0?` ${a.resName}${a.resSeq}${a.iCode||''}${a.chain?':'+a.chain:''}`:'';return (a.name||el)+res}
