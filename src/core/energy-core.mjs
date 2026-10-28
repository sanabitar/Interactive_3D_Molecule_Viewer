import {projectionBasis,newmanBranches} from './conformation-core.mjs';
export function ethanePhase(atoms,bonds,orders,front,back){
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
export function ethaneEnergy(rotation,phase){return 6*(1+Math.cos(3*(phase-rotation)*Math.PI/180));}
