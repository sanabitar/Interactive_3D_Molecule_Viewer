import test from 'node:test';import assert from 'node:assert/strict';
import {checkMirror,reflect} from '../src/core/mirror-core.mjs';
const tetra=els=>[[0,0,0],[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(([x,y,z],i)=>({x,y,z,el:els[i]}));
const bonds=[[0,1],[0,2],[0,3],[0,4]];
test('four distinct substituents cannot superimpose',()=>{const r=checkMirror(tetra(['C','H','F','Cl','Br']),bonds);assert.equal(r.matched,false);assert.equal(r.complete,true);assert.ok(r.best.rmsd>.1);});
test('equivalent atoms permit methane mirror superposition',()=>{const r=checkMirror(tetra(['C','H','H','H','H']),bonds);assert.equal(r.matched,true);assert.ok(r.best.rmsd<1e-8);});
test('two equivalent substituents remove tetrahedral chirality',()=>assert.equal(checkMirror(tetra(['C','H','H','Cl','Br']),bonds).matched,true));
test('planar geometry matches and reflection preserves distances',()=>{const a=[{x:0,y:0,z:0,el:'O'},{x:1,y:1,z:0,el:'H'},{x:-1,y:1,z:0,el:'H'}];assert.equal(checkMirror(a,[[0,1],[0,2]]).matched,true);assert.deepEqual(reflect(reflect(a)),a);});
test('bounded search reports incomplete and rejects invalid size',()=>{assert.equal(checkMirror(tetra(['C','H','H','H','H']),bonds,[],{limit:1}).complete,false);assert.throws(()=>checkMirror([],[]));});
