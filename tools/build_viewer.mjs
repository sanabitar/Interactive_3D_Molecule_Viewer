// Generates the three files in js/generated/ from src/ (src/ is the single source of truth for this code):
//   js/generated/hydrogens.js  <- src/core/hydrogens.mjs
//   js/generated/layout2d.js   <- src/core/layout2d.mjs
//   js/generated/labs.js       <- src/labs/*.mjs plus the src/core modules they import
// They are written as classic scripts (no ES modules), so index.html also works when it is opened by double-click (file://), not only from a server.
//
//   node tools/build_viewer.mjs          write the generated files
//   node tools/build_viewer.mjs --check  exit 1 if they are out of date (npm test runs the same check)
import fs from 'node:fs';import path from 'node:path';import {fileURLToPath,pathToFileURL} from 'node:url';
const defaultRoot=path.join(path.dirname(fileURLToPath(import.meta.url)),'..');
const LABS_ORDER=['geometry','superpose','mirror-core','conformation-core','energy-core','energy-panel','sugar-core','chair-core','electron-core','mirror-lab','conformation-lab','sugar-lab','chair-lab','electron-lab'];
const banner=src=>`// GENERATED FILE - do not edit. Source: ${src}\n// Rebuild with: node tools/build_viewer.mjs\n`;

// Returns {relativePath: fileContents} for every generated file.
export function generate(root=defaultRoot){
 const rd=p=>fs.readFileSync(path.join(root,p),'utf8').replace(/\r\n/g,'\n');
 const files={};
 // one src/core module wrapped as `const Name=(()=>{...return{exports}})();`
 for(const[name,file,out]of[['Hydrogens','src/core/hydrogens.mjs','js/generated/hydrogens.js'],['Layout2D','src/core/layout2d.mjs','js/generated/layout2d.js']]){
  const raw=rd(file),src=raw.replace(/^export /gm,''),names=[...raw.matchAll(/^export\s+(?:async\s+)?(?:function\*?|const|let|class)\s+([A-Za-z0-9_$]+)/gm)].map(m=>m[1]);
  files[out]=banner(file)+'const '+name+'=(()=>{\n'+src+'\nreturn{'+names.join(',')+'}})();\n'}
 // labs: every module becomes __M['name']=(()=>{...})(); imports between them are rewritten to lookups in __M
 let body='const __M={};\n';
 for(const n of LABS_ORDER){const f=fs.existsSync(path.join(root,'src/core',n+'.mjs'))?'src/core/'+n+'.mjs':'src/labs/'+n+'.mjs';let src=rd(f);const names=[];
  src=src.replace(/^import\s*\{([^}]*)\}\s*from\s*'(?:\.\/|\.\.\/core\/)([^']+)\.mjs';?/gm,(m,l,mod)=>{if(!LABS_ORDER.includes(mod))throw Error('unknown dependency '+mod+' in '+f);return `const {${l}}=__M['${mod}'];`});
  src=src.replace(/^export\s+((?:async\s+)?function\*?|const|let|class)\s+([A-Za-z0-9_$]+)/gm,(m,k,nm)=>{names.push(nm);return k+' '+nm});
  if(/^(import|export) /m.test(src))throw Error('unhandled import/export in '+f);
  body+=`__M['${n}']=(()=>{\n${src}\nreturn {${names.join(',')}};})();\n`}
 files['js/generated/labs.js']=banner('src/labs/*.mjs and the src/core modules they import')+'(function(){try{\n'+body+'}catch(err){console.error("LAURA labs failed to load",err)}})();\n';
 return files}

if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){
 const files=generate(),check=process.argv.includes('--check');let stale=0;
 for(const[rel,text]of Object.entries(files)){const p=path.join(defaultRoot,rel),cur=fs.existsSync(p)?fs.readFileSync(p,'utf8').replace(/\r\n/g,'\n'):null;
  if(check){if(cur!==text){stale++;console.error('OUT OF DATE: '+rel)}}else{fs.mkdirSync(path.dirname(p),{recursive:true});if(cur!==text)fs.writeFileSync(p,text);console.log((cur===text?'unchanged ':'wrote     ')+rel+' ('+text.length+' bytes)')}}
 if(check){if(stale){console.error('Run: node tools/build_viewer.mjs');process.exit(1)}console.log('generated files are up to date')}}
