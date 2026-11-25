// Reads the viewer's own source files for tests (the app is plain files: index.html + css/ + js/, no bundler).
import fs from 'node:fs';
const root=new URL('../../',import.meta.url);
export const read=rel=>fs.readFileSync(new URL(rel,root),'utf8').replace(/\r\n/g,'\n');
export const indexHtml=()=>read('index.html');
// local scripts and stylesheets in the order index.html loads them
export const localScripts=()=>[...indexHtml().matchAll(/<script src="(js\/[^"]+)"><\/script>/g)].map(m=>m[1]);
export const localStyles=()=>[...indexHtml().matchAll(/<link rel="stylesheet" href="(css\/[^"]+)">/g)].map(m=>m[1]);
export const viewerScript=()=>localScripts().map(read).join('\n');   // all application JS in load order
