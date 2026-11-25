// js/ui/modal.js: Modal dialog helper (modal(), focus trap) and need(), which asks for a structure to be loaded first.
let lastFocus=null;
function modal(h){const m=$('#modal'),b=$('#mb');if(h){if(m.style.display!=='flex')lastFocus=document.activeElement;m.style.display='flex';b.innerHTML=h;b.setAttribute('role','dialog');b.setAttribute('aria-modal','true');b.setAttribute('aria-label','LAURA dialog');b.setAttribute('tabindex','-1');b.focus()}else{m.style.display='none';b.innerHTML='';if(lastFocus&&lastFocus.focus)lastFocus.focus();lastFocus=null}}
$('#modal').addEventListener('keydown',e=>{if(e.key!=='Tab')return;const f=[...$('#mb').querySelectorAll('button,input,select,textarea,summary,[tabindex]:not([tabindex="-1"])')].filter(x=>!x.disabled&&x.offsetParent!==null);if(!f.length){e.preventDefault();return}const first=f[0],last=f[f.length-1],cur=document.activeElement;
if(e.shiftKey&&(cur===first||cur===$('#mb'))){e.preventDefault();last.focus()}else if(!e.shiftKey&&cur===last){e.preventDefault();first.focus()}});const CL='<p><button onclick="modal()">Close</button></p>';
const need=()=>mol||(toast('Load a structure first.'),0);
