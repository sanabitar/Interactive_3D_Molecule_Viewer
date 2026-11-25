import test from 'node:test';import assert from 'node:assert/strict';
import {generate} from '../tools/build_viewer.mjs';import {read} from './helpers/viewer-source.mjs';
// js/generated/* are built from src/ by tools/build_viewer.mjs. If you edit src/core/hydrogens.mjs, src/core/lonepairs.mjs, src/core/layout2d.mjs or anything in src/labs,
// run `node tools/build_viewer.mjs`; this test fails until you do (so a stale generated file can never be committed unnoticed).
const files=generate();
test('four files are generated: hydrogens, lonepairs, layout2d, labs',()=>{assert.deepEqual(Object.keys(files).sort(),['js/generated/hydrogens.js','js/generated/labs.js','js/generated/layout2d.js','js/generated/lonepairs.js'])});
for(const[rel,text]of Object.entries(files))test(`${rel} is up to date with src/ (run: node tools/build_viewer.mjs)`,()=>{assert.equal(read(rel),text)});
test('generated files are marked as such and are classic scripts (no import/export), so index.html works from file://',()=>{
 for(const[rel,text]of Object.entries(files)){assert.match(text,/^\/\/ GENERATED FILE - do not edit/);assert.ok(!/^(import|export) /m.test(text),rel)}});
