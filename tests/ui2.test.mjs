import test from 'node:test';import assert from 'node:assert/strict';
import {projectPlane} from '../src/core/depict2d.mjs';import {declutter,labelText} from '../src/core/labels.mjs';import {cleanName} from '../src/core/ui.mjs';
const rot=(p,k,t)=>{const L=Math.hypot(...k),u=k.map(v=>v/L),c=Math.cos(t),s=Math.sin(t),d=u[0]*p.x+u[1]*p.y+u[2]*p.z,cr=[u[1]*p.z-u[2]*p.y,u[2]*p.x-u[0]*p.z,u[0]*p.y-u[1]*p.x];return{x:p.x*c+cr[0]*s+u[0]*d*(1-c),y:p.y*c+cr[1]*s+u[1]*d*(1-c),z:p.z*c+cr[2]*s+u[2]*d*(1-c)}};
const D3=(a,b)=>Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z),D2=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1]);
test('2D view: a planar molecule rotated arbitrarily in 3D is flattened isometrically (all distances preserved)',()=>{
 const hex=Array.from({length:6},(_,i)=>({x:1.397*Math.cos(i*Math.PI/3),y:1.397*Math.sin(i*Math.PI/3),z:0})).concat([{x:2.5,y:0,z:0},{x:-2.5,y:0.3,z:0}]),P=hex.map(p=>{const r=rot(p,[1,2,3],.9);return{x:r.x+7,y:r.y-3,z:r.z+11}}),R=projectPlane(P);
 for(let i=0;i<P.length;i++){assert.ok(Math.abs(R.depth[i])<1e-9,'depth ~0 for planar');for(let j=i+1;j<P.length;j++)assert.ok(Math.abs(D3(P[i],P[j])-D2(R.xy[i],R.xy[j]))<1e-9,`pair ${i},${j}`)}assert.ok(R.variance[2]<1e-12)});
test('2D view: non-planar input never stretches distances; degenerate inputs are safe',()=>{const P=[{x:0,y:0,z:0},{x:2,y:0,z:0},{x:0,y:2,z:0},{x:1,y:1,z:1.5}],R=projectPlane(P);
 for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)assert.ok(D2(R.xy[i],R.xy[j])<=D3(P[i],P[j])+1e-9);assert.ok(R.variance[0]>=R.variance[1]&&R.variance[1]>=R.variance[2]);
 assert.deepEqual(projectPlane([{x:1,y:2,z:3}]).xy,[[0,0]]);assert.equal(projectPlane([]).xy.length,0);const line=projectPlane([0,1,2,3].map(i=>({x:i,y:2*i,z:-i})));assert.ok(line.xy.every(p=>Math.abs(p[1])<1e-9),'collinear points lie on one axis')});
test('label declutter: no overlaps, priority first, off-screen dropped, cap respected, identical positions keep only one',()=>{
 let s=7;const rnd=()=>(s=(s*1103515245+12345)&0x7fffffff)/0x7fffffff,items=Array.from({length:400},(_,i)=>({x:rnd()*1000,y:rnd()*700,w:40,h:18,prio:i%5}));
 const acc=declutter(items,{width:1000,height:700,cap:60,pad:3});assert.ok(acc.length>0&&acc.length<=60);
 for(let a=0;a<acc.length;a++)for(let b=a+1;b<acc.length;b++){const p=items[acc[a]],q=items[acc[b]];const ox=Math.abs(p.x-q.x)<p.w/2+q.w/2,oy=Math.abs(p.y-q.y)<p.h/2+q.h/2;assert.ok(!(ox&&oy),'accepted labels overlap')}
 const two=declutter([{x:100,y:100,w:30,h:14,prio:5},{x:100,y:100,w:30,h:14,prio:1}]);assert.deepEqual(two,[1],'lower prio value wins; the duplicate position is dropped');
 assert.deepEqual(declutter([{x:-20,y:50,w:10,h:10},{x:50,y:50,w:10,h:10}],{width:100,height:100}),[1],'off-screen label dropped');
 assert.equal(declutter(Array.from({length:50},(_,i)=>({x:10+i*20,y:50,w:10,h:10})),{width:2000,height:100,cap:7}).length,7,'cap');
 const dense=Array.from({length:5000},(_,i)=>({x:500+(i%7),y:300+(i%5),w:40,h:18}));assert.equal(declutter(dense,{width:1000,height:700}).length,1,'a pile of 5000 coincident labels collapses to one readable label')});
test('label text: atom name (+residue) for real structures, element symbol for generated molecules',()=>{
 assert.equal(labelText({element:'O',el:'O',name:'O2',autoName:1}),'O','generated "O2" is an id, not a name');assert.equal(labelText({element:'C',name:'CA',resName:'ALA',resSeq:5,iCode:'',chain:'A'}),'CA ALA5:A');
 assert.equal(labelText({element:'C',name:'C1',resName:'LIG',resSeq:101,iCode:'A',chain:''}),'C1 LIG101A');assert.equal(labelText({element:'C',name:'C1',resName:'',resSeq:0}),'C1')});
test('rename: names are sanitised',()=>{assert.equal(cleanName('  My   helix \n 2 '),'My helix 2');assert.equal(cleanName('a\u0000b\u0007c'),'a b c');assert.equal(cleanName('   '),null);assert.equal(cleanName(null),null);assert.equal(cleanName('x'.repeat(200)).length,80);assert.equal(cleanName('<b>x</b>'),'<b>x</b>','markup is kept as text; rendering escapes it')});
