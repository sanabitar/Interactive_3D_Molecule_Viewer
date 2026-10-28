const xyz=a=>[a.x,a.y,a.z];
const sub=(a,b)=>a.map((v,i)=>v-b[i]);
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const unit=a=>{const n=Math.hypot(...a);if(n<1e-8)throw Error('Bond axis has zero length.');return a.map(v=>v/n);};
export function movingSide(n,bonds,front,back){
 const adj=Array.from({length:n},()=>[]);for(const [a,b] of bonds){if((a===front&&b===back)||(a===back&&b===front))continue;adj[a].push(b);adj[b].push(a);}
 const seen=new Set([back]),todo=[back];while(todo.length){for(const j of adj[todo.pop()])if(!seen.has(j)){seen.add(j);todo.push(j);}}
 return seen.has(front)?null:[...seen];
}
export function eligibleBonds(atoms,bonds,orders){
 const adj=Array.from({length:atoms.length},()=>[]);bonds.forEach(([a,b])=>{adj[a].push(b);adj[b].push(a);});
 return bonds.flatMap(([a,b],k)=>{const el=i=>atoms[i].element||atoms[i].el;
  if(orders[k]!==1||el(a)==='H'||el(b)==='H'||adj[a].length<2||adj[b].length<2)return [];
  // Exclude amide-like C–N bonds; their partial double-bond character restricts rotation.
  if(((el(a)==='C'&&el(b)==='N')||(el(b)==='C'&&el(a)==='N'))&&bonds.some(([u,v],j)=>orders[j]===2&&((u===(el(a)==='C'?a:b)&&['O','S'].includes(el(v)))||(v===(el(a)==='C'?a:b)&&['O','S'].includes(el(u))))))return [];
  const moving=movingSide(atoms.length,bonds,a,b);return moving?[{front:a,back:b,moving,index:k}]:[];
 });
}
export function rotateBond(atoms,front,back,moving,degrees){
 const p=xyz(atoms[front]),u=unit(sub(xyz(atoms[back]),p)),angle=degrees*Math.PI/180,c=Math.cos(angle),s=Math.sin(angle),set=new Set(moving);
 return atoms.map((a,i)=>{if(!set.has(i)||i===front||i===back)return {...a};const v=sub(xyz(a),p),w=cross(u,v),d=dot(u,v),r=v.map((v,k)=>p[k]+v*c+w[k]*s+u[k]*d*(1-c));return {...a,x:r[0],y:r[1],z:r[2]};});
}
// Looking from front towards back. Basis remains anchored to the original coordinates.
export function projectionBasis(atoms,front,back){const u=unit(sub(xyz(atoms[back]),xyz(atoms[front]))),helper=Math.abs(u[2])<.9?[0,0,1]:[0,1,0],x=unit(cross(u,helper)),y=cross(x,u);return {u,x,y};}
export function newmanBranches(atoms,bonds,front,back,basis){
 return [front,back].flatMap((center,side)=>bonds.flatMap(([a,b])=>{const j=a===center?b:b===center?a:-1;if(j<0||j===front||j===back)return [];const v=sub(xyz(atoms[j]),xyz(atoms[center])),x=dot(v,basis.x),y=dot(v,basis.y);return [{index:j,side,x,y,angle:Math.atan2(y,x),axial:Math.hypot(x,y)<1e-7}];}));
}
