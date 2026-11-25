import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import vm from 'node:vm';
import {indexHtml,localScripts,localStyles,read} from './helpers/viewer-source.mjs';
// The app is plain files: index.html + css/*.css + js/**/*.js, loaded as classic scripts in a fixed order (see docs/ARCHITECTURE.md).
const html=indexHtml(),scripts=localScripts(),styles=localStyles();
const walk=d=>fs.readdirSync(new URL('../'+d,import.meta.url),{withFileTypes:true}).flatMap(e=>e.isDirectory()?walk(d+'/'+e.name):[d+'/'+e.name]);
test('index.html is markup only: no inline <script> code, no <style>, no inline event handlers',()=>{
 assert.ok(!/<script(?![^>]*\ssrc=)[^>]*>/.test(html),'inline script');assert.ok(!/<style/i.test(html),'<style>');assert.ok(!/\son(click|change|input|load)=/i.test(html),'inline handler')});
test('every referenced file exists, is referenced once, and no file in css/ or js/ is left unreferenced',()=>{
 for(const f of [...scripts,...styles])assert.ok(fs.existsSync(new URL('../'+f,import.meta.url)),'missing '+f);
 assert.equal(new Set(scripts).size,scripts.length);assert.equal(new Set(styles).size,styles.length);
 assert.deepEqual(walk('js').sort(),[...scripts].sort(),'js/ files vs <script> tags');assert.deepEqual(walk('css').sort(),[...styles].sort(),'css/ files vs <link> tags')});
test('load order: constants first, hydrogens before the library, core before the viewer, main.js last of the app scripts, labs after main.js',()=>{
 const at=f=>{const i=scripts.indexOf('js/'+f);assert.ok(i>=0,f);return i};
 assert.equal(at('foundation/constants.js'),0);assert.ok(at('generated/hydrogens.js')<at('foundation/library.js'));assert.ok(at('foundation/core.js')<at('viewer/scene.js'));
 assert.ok(at('generated/layout2d.js')<at('tools/view2d.js'));assert.ok(at('app/mirror-source.js')>at('main.js'));assert.ok(at('generated/labs.js')>at('main.js'));
 const app=scripts.filter(s=>!['js/app/mirror-source.js','js/generated/labs.js'].includes(s));assert.equal(app.at(-1),'js/main.js')});
test('every script parses as a classic script (these files must also work from file://, so no ES modules)',()=>{
 for(const f of scripts){assert.ok(!/^(import|export) /m.test(read(f)),f+' uses import/export');new vm.Script(read(f),{filename:f})}});
test('start-up rule: only main.js starts work at load (renderLoop, restoreSavedImports); other files just define things',()=>{
 for(const f of scripts){const t=read(f);if(f==='js/main.js'){assert.ok(/\nrenderLoop\(\);restoreSavedImports\(\);/.test(t))}else assert.ok(!/^renderLoop\(\)|^restoreSavedImports\(\)/m.test(t),f)}});
test('stylesheets keep the original rule order: theme first, then base, and the later override files after them',()=>{
 assert.deepEqual(styles,['css/theme.css','css/base.css','css/viewer.css','css/layout.css','css/controls.css','css/dialogs.css','css/panels.css','css/overlays.css','css/import-and-plots.css'])});
test('every js/ and css/ file starts with a header saying what it is for',()=>{
 for(const f of scripts)if(!f.startsWith('js/generated/'))assert.match(read(f),new RegExp('^// '+f.replace(/[.*+?^${}()|[\]\\/]/g,'\\$&')+': \\S'),f);
 for(const f of styles)assert.match(read(f),new RegExp('^/\\* '+f.replace(/[.*+?^${}()|[\]\\/]/g,'\\$&')+': \\S'),f)});
