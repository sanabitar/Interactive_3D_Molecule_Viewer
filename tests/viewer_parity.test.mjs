import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';
import * as analysis from '../src/core/analysis.mjs';import {peptide} from '../src/core/build.mjs';import * as models from '../src/core/models.mjs';
// The viewer inlines the core modules. This guards against the inlined copy drifting from src/core.
const html=fs.readFileSync(new URL('../laura-viewer.html',import.meta.url),'utf8');
const i=html.indexOf('const Core=(()=>{'),j=html.indexOf('})();',html.indexOf('return{parsePDB,'))+5;
const Core=new Function(html.slice(i,j)+';return Core')();
const MAP=a=>({...a,element:a.el,resName:a.res,resSeq:a.resi,hetero:false,iCode:''});
test('viewer Core.peptide equals src peptide',()=>{const o={phi:[0,-60,-120,60],psi:[130,-40,130,40],seq:['ALA','GLY','SER','PRO']};assert.deepEqual(Core.peptide(4,o),peptide(4,o));assert.deepEqual(Core.peptide(12,{phi:-57,psi:-47}),peptide(12,{phi:-57,psi:-47}))});
test('viewer models equal src models (atom for atom)',()=>{for(const k of['helixLoopHelix','betaHairpin','ramaTestPeptide','antiparallelHelixPair'])assert.deepEqual(Core.models[k](),models[k](),k);assert.deepEqual(Core.models.MODEL_SEQ,models.MODEL_SEQ)});
test('viewer phiPsi / regions / summary equal src',()=>{const A=models.antiparallelHelixPair().map(MAP),B=models.ramaTestPeptide().map(MAP);
 for(const X of[A,B]){assert.deepEqual(Core.phiPsi(X),analysis.phiPsi(X));assert.deepEqual(Core.ramaSummary(Core.phiPsi(X)),analysis.ramaSummary(analysis.phiPsi(X)))}
 assert.deepEqual(Core.RAMA_REGIONS,analysis.RAMA_REGIONS);for(const[f,p,r]of[[-57,-47,'ALA'],[60,-120,'ALA'],[60,-120,'GLY'],[180,180,'ALA']])assert.equal(Core.ramaRegion(f,p,r),analysis.ramaRegion(f,p,r))});
