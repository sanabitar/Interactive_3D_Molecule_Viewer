// js/foundation/constants.js: Shared constants and helpers: $ (querySelector), element colours and radii (CPK, VDW, COV), toast(), esc(). Loaded first.
const $=s=>document.querySelector(s);
function maxOf(arr,f,floor=-Infinity){let m=floor;for(const x of arr){const v=f(x);if(v>m)m=v}return m} // no spread: Math.max(...bigArray) overflows the call stack above ~125k items
const CPK={H:0xffffff,C:0x8a8f98,N:0x3050f8,O:0xff2a2a,P:0xff8000,S:0xffd62a,X:0x1fa36b};
const VDW={H:1.2,C:1.7,N:1.55,O:1.52,P:1.8,S:1.8,X:1.5},COV={H:.31,C:.76,N:.71,O:.66,P:1.07,S:1.05,X:.7};
const vdw=e=>VDW[e]||1.7,cov=e=>COV[e]||.8,cpk=e=>CPK[e]??0xcc66cc;
function toast(m,ms=4500){const t=$('#toast');t.textContent=m;t.style.display='block';clearTimeout(toast.t);toast.t=setTimeout(()=>t.style.display='none',ms)}
const esc=t=>String(t).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
