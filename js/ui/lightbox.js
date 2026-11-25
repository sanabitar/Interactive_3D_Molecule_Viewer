// js/ui/lightbox.js: Closing dialogs on outside click, the enlarge-a-diagram lightbox (zoom and pan), and the 2D view button.
// ---------- close-on-outside-click, enlarge diagrams, 2D button ----------
let bdDown=false;$('#modal').addEventListener('pointerdown',e=>{bdDown=e.target===$('#modal')});$('#modal').addEventListener('click',e=>{if(bdDown&&e.target===$('#modal'))modal();bdDown=false});
document.addEventListener('pointerdown',e=>{const t=$('#tch');if(t.style.display!=='block')return;if(t.contains(e.target)||e.target.closest('#te')||e.target.closest('#vp canvas')||e.target.closest('#modal'))return;tch()},true);
$('#mb').addEventListener('click',e=>{const c=e.target.closest&&e.target.closest('canvas');if(!c)return;try{$('#lbimg').src=c.toDataURL('image/png');lbReset();$('#lightbox').classList.add('open')}catch(x){toast('This picture cannot be enlarged.')}});
const LB=$('#lightbox'),LV=$('#lbview'),LI=$('#lbimg'),lbPts=new Map();let lbS=1,lbX=0,lbY=0,lbD=0;
function lbApply(){LI.style.transform='translate('+lbX+'px,'+lbY+'px) scale('+lbS+')';LV.classList.toggle('zm',lbS>1.001)}
function lbClamp(){const r=LV.getBoundingClientRect(),mx=r.width*(lbS-1)/2,my=r.height*(lbS-1)/2;lbX=Math.max(-mx,Math.min(mx,lbX));lbY=Math.max(-my,Math.min(my,lbY))}
function lbReset(){lbS=1;lbX=lbY=0;lbApply()}
function lbZoom(f,x,y){const r=LV.getBoundingClientRect(),cx=r.left+r.width/2+lbX,cy=r.top+r.height/2+lbY,s2=Math.max(1,Math.min(16,lbS*f));if(s2===lbS)return;const k=1-s2/lbS;lbX+=(x-cx)*k;lbY+=(y-cy)*k;lbS=s2;lbClamp();lbApply()}
const lbMid=()=>{const r=LV.getBoundingClientRect();return[r.left+r.width/2,r.top+r.height/2]},lbClose=()=>LB.classList.remove('open');
LV.addEventListener('wheel',e=>{e.preventDefault();lbZoom(Math.exp(-e.deltaY*(e.ctrlKey?.01:.0016)),e.clientX,e.clientY)},{passive:false});
LV.addEventListener('pointerdown',e=>{LV.setPointerCapture(e.pointerId);lbPts.set(e.pointerId,[e.clientX,e.clientY]);if(lbPts.size===2){const[a,b]=[...lbPts.values()];lbD=Math.hypot(a[0]-b[0],a[1]-b[1])}else if(lbS>1)LV.classList.add('pan')});
LV.addEventListener('pointermove',e=>{const p=lbPts.get(e.pointerId);if(!p)return;if(lbPts.size===1){if(lbS>1){lbX+=e.clientX-p[0];lbY+=e.clientY-p[1];lbClamp();lbApply()}lbPts.set(e.pointerId,[e.clientX,e.clientY])}else{lbPts.set(e.pointerId,[e.clientX,e.clientY]);const[a,b]=[...lbPts.values()],d=Math.hypot(a[0]-b[0],a[1]-b[1]);if(lbD>0)lbZoom(d/lbD,(a[0]+b[0])/2,(a[1]+b[1])/2);lbD=d}});
const lbUp=e=>{lbPts.delete(e.pointerId);lbD=0;LV.classList.remove('pan')};LV.addEventListener('pointerup',lbUp);LV.addEventListener('pointercancel',lbUp);
LV.addEventListener('dblclick',e=>lbZoom(lbS>1.5?1/lbS:3,e.clientX,e.clientY));
$('#lbp').onclick=()=>lbZoom(1.5,...lbMid());$('#lbm').onclick=()=>lbZoom(1/1.5,...lbMid());$('#lbr').onclick=lbReset;$('#lbx').onclick=lbClose;LB.addEventListener('click',e=>{if(e.target===LB)lbClose()});
addEventListener('keydown',e=>{if(!LB.classList.contains('open')||e.ctrlKey||e.metaKey||e.altKey)return;if(e.key==='+'||e.key==='='){e.preventDefault();lbZoom(1.5,...lbMid())}else if(e.key==='-'){e.preventDefault();lbZoom(1/1.5,...lbMid())}else if(e.key==='0'){e.preventDefault();lbReset()}});
$('#d2').onclick=show2D;
