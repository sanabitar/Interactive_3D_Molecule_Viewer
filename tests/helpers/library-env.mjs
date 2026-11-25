// Evaluates the viewer's library + Core (js/foundation/*, js/generated/hydrogens.js) in Node, without a DOM, so tests can inspect the 43 built-in molecules.
import {read} from './viewer-source.mjs';
export const env=new Function([read('js/foundation/constants.js'),read('js/generated/hydrogens.js'),read('js/foundation/library.js'),read('js/foundation/core.js')].join('\n')+'\nreturn {LIB,Core,computeBonds,Hydrogens};')();
// The atom records exactly as the viewer's load() builds them from a library entry (needed for residue-based bond perception).
export const viewerAtoms=entry=>entry.get().atoms.map(a=>({hetero:false,iCode:'',element:a.el,resName:a.res,resSeq:a.resi,...a}));
