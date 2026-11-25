// js/ui/panels.js: Resizable library/details panels (widths are remembered) and the details-panel toggle (the panel starts closed).
// ---------- resizable panels ----------
const M=$('main');
function setW(k,v){M.style.setProperty(k,v+'px')}
function saveLayout(){try{localStorage.setItem('laura-layout',JSON.stringify({lw:M.style.getPropertyValue('--lw'),rw:M.style.getPropertyValue('--rw')}))}catch(e){}} // widths only: the details panel always starts closed
function drag(h,onMove){h.addEventListener('pointerdown',e=>{e.preventDefault();h.setPointerCapture(e.pointerId);h.classList.add('on');const up=()=>{h.classList.remove('on');h.removeEventListener('pointermove',onMove);h.removeEventListener('pointerup',up);saveLayout()};h.addEventListener('pointermove',onMove);h.addEventListener('pointerup',up)})}
const CL2=(v,lo,hi)=>Math.max(lo,Math.min(hi,v));
drag($('#gl'),e=>setW('--lw',CL2(e.clientX-M.getBoundingClientRect().left,180,520)));
drag($('#gr'),e=>setW('--rw',CL2(M.getBoundingClientRect().right-e.clientX,220,560)));
$('#gl').ondblclick=()=>{M.style.removeProperty('--lw');saveLayout()};$('#gr').ondblclick=()=>{M.style.removeProperty('--rw');saveLayout()};
[['#gl','--lw',1,180,520,260],['#gr','--rw',-1,220,560,300]].forEach(([id,v,dir,lo,hi,def])=>$(id).addEventListener('keydown',e=>{if(e.key!=='ArrowLeft'&&e.key!=='ArrowRight')return;e.preventDefault();setW(v,CL2((parseFloat(M.style.getPropertyValue(v))||def)+(e.key==='ArrowRight'?1:-1)*dir*20,lo,hi));saveLayout()}));
function rightPanel(show){M.classList.toggle('nor',!show);const b=$('#rshow');b.classList.toggle('on',show);b.setAttribute('aria-expanded',show);b.title=show?'Close details panel':'Open details panel';b.setAttribute('aria-label',b.title)}
$('#rshow').onclick=()=>rightPanel(M.classList.contains('nor'));rightPanel(false);
try{const L=JSON.parse(localStorage.getItem('laura-layout')||'null');if(L){if(L.lw)M.style.setProperty('--lw',L.lw);if(L.rw)M.style.setProperty('--rw',L.rw)}}catch(e){}
