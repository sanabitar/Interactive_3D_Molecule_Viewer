// Inlines src/labs + their src/core dependencies into laura-viewer.html (between the LABS markers),
// so the viewer works from a double-click (file://) as well as from a server. Run: node tools/build_viewer.mjs
import fs from 'node:fs';import path from 'node:path';import {fileURLToPath} from 'node:url';
const root=path.join(path.dirname(fileURLToPath(import.meta.url)),'..'),rd=p=>fs.readFileSync(path.join(root,p),'utf8').replace(/\r\n/g,'\n');
const order=['geometry','superpose','mirror-core','conformation-core','energy-core','energy-panel','sugar-core','chair-core','electron-core','mirror-lab','conformation-lab','sugar-lab','chair-lab','electron-lab'];
let out='const __M={};\n';
for(const n of order){const f=fs.existsSync(path.join(root,'src/core',n+'.mjs'))?'src/core/'+n+'.mjs':'src/labs/'+n+'.mjs';let src=rd(f);const names=[];
 src=src.replace(/^import\s*\{([^}]*)\}\s*from\s*'(?:\.\/|\.\.\/core\/)([^']+)\.mjs';?/gm,(m,l,mod)=>{if(!order.includes(mod))throw Error('unknown dependency '+mod+' in '+f);return `const {${l}}=__M['${mod}'];`});
 src=src.replace(/^export\s+((?:async\s+)?function\*?|const|let|class)\s+([A-Za-z0-9_$]+)/gm,(m,k,nm)=>{names.push(nm);return k+' '+nm});
 if(/^(import|export) /m.test(src))throw Error('unhandled import/export in '+f);out+=`__M['${n}']=(()=>{\n${src}\nreturn {${names.join(',')}};})();\n`}
let html=rd('laura-viewer.html');const a=html.indexOf('<!--LABS-BEGIN-->'),b=html.indexOf('<!--LABS-END-->');if(a<0||b<a)throw Error('LABS markers missing');
html=html.slice(0,a)+'<!--LABS-BEGIN-->\n<script>\n(function(){try{\n'+out+'}catch(err){console.error("LAURA labs failed to load",err)}})();\n</script>\n'+html.slice(b);
fs.writeFileSync(path.join(root,'laura-viewer.html'),html);console.log('labs inlined,',html.length,'bytes');
