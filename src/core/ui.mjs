// Pure UI decision helpers (testable without a DOM).
export function shortcutAllowed(e,{modalOpen=false}={}){if(e.ctrlKey||e.metaKey||e.altKey)return false;const k=String(e.key||'').toLowerCase();if(k==='escape')return true;if(modalOpen)return false;
 const t=e.target||{},tag=String(t.tagName||'').toUpperCase();if(t.isContentEditable||tag==='TEXTAREA'||tag==='SELECT')return false;if(tag==='INPUT'&&t.type!=='range')return false;
 if(k===' '&&(tag==='BUTTON'||tag==='SUMMARY'||tag==='A'||(t.getAttribute&&t.getAttribute('role')==='separator')))return false;return true}
// Display-name sanitising for renamed imports: trims, collapses whitespace, removes control characters, caps length. Returns null if nothing is left.
export function cleanName(s,max=80){const t=String(s??'').replace(/[\u0000-\u001f\u007f]/g,' ').replace(/\s+/g,' ').trim().slice(0,max).trim();return t.length?t:null}
