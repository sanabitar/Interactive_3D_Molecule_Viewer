// js/main.js: Wiring and START-UP: theme, keyboard shortcuts and the remaining buttons. Loaded last; the only file that starts work at load.
// ---------- wiring ----------
// START-UP: the only place where work is started at load. Every other file just defines things, so a frame, observer or promise callback
// can never run before a later file has defined what it needs.
renderLoop();restoreSavedImports();
function theme(t){document.documentElement.dataset.theme=t;thBtn();try{localStorage.setItem('laura-theme',t)}catch(e){}requestAnimationFrame(()=>{setBg();build()})}
$('#th').onclick=()=>theme((document.documentElement.dataset.theme||(matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light'))==='dark'?'light':'dark');
$('#q').oninput=renderList;$('#mol').onchange=e=>{const m=[...LIB,...imported].find(x=>x.id===e.target.value);if(m)load(m)};$('#srt').onchange=()=>{try{localStorage.setItem('laura-sort',$('#srt').value)}catch(e){}renderList();molSel()};$('#cat').onchange=renderList;$('#stb').onclick=()=>{const s=$('#sty');s.dataset.u=1;s.value=s.value==='tb'?'bs':'tb';build()};$('#ww').onchange=$('#bo').onchange=build;$('#el').onchange=()=>{lbDirty=true};$('#hh').onchange=()=>setH($('#hh').checked);$('#hb').onclick=()=>setH(!$('#hh').checked);
$('#rs').onclick=fit;$('#fz').onclick=focusSel;$('#ar').onclick=()=>{if(ctl){ctl.autoRotate=!ctl.autoRotate;if(ctl.autoRotate){ctl.autoRotateSpeed=1;$('#sp').value=1}$('#spl').style.display=ctl.autoRotate?'flex':'none';$('#ar').classList.toggle('on',ctl.autoRotate);$('#ar').setAttribute('aria-pressed',ctl.autoRotate)}};$('#sp').oninput=e=>{if(ctl)ctl.autoRotateSpeed=+e.target.value};
function popSelBy(){const A=mol.atoms,u=k=>[...new Set(A.map(a=>a[k]).filter(v=>v!==''&&v!=null))].sort(),g=(l,k)=>{const v=u(k);return v.length?`<optgroup label="${l}">`+v.map(x=>`<option value="${l[0]}:${esc(x)}">${esc(x)}</option>`).join('')+'</optgroup>':''};
$('#selby').innerHTML='<option value="">Select by…</option>'+g('Element','el')+g('Chain','chain')+g('Residue','res')}
$('#selby').onchange=e=>{const v=e.target.value;if(!v||!mol)return;const f={E:a=>a.el,C:a=>a.chain,R:a=>a.res}[v[0]],val=v.slice(2);sel=visBase.filter(i=>f(mol.atoms[i])===val);e.target.value='';refreshSel()};
$('#sall').onclick=()=>{if(mol){sel=[...visBase];refreshSel()}};$('#sinv').onclick=()=>{if(mol){const S=new Set(sel);sel=visBase.filter(i=>!S.has(i));refreshSel()}};
$('#clr').onclick=()=>{sel=[];if(mol)refreshSel()};$('#so').onchange=()=>{if($('#so').checked&&!sel.length)toast('Select some atoms first; Selected only then hides everything else.');if(mol){build();if($('#so').checked&&sel.length)focusSel();else if(!$('#so').checked)fit()}};
let pmFS=false;
function pres(){document.body.classList.toggle('pm');const on=document.body.classList.contains('pm');if(on){const d=document.documentElement,warn=()=>toast('Your browser blocked full screen here (common when the page is embedded). The viewer now fills the whole page; use F11 or open the page in its own tab for true full screen.',8000);try{const p=d.requestFullscreen?d.requestFullscreen():null;if(p&&p.catch)p.catch(warn);else if(!p)warn()}catch(e){warn()}}else{pmFS=false;if(document.fullscreenElement)document.exitFullscreen()}requestAnimationFrame(()=>{size();fit();setTimeout(size,150)})}
document.addEventListener('fullscreenchange',()=>{const pm=document.body.classList.contains('pm');if(fsEl()&&pm)pmFS=true;else if(!fsEl()&&pm&&pmFS){pmFS=false;pres()}})
$('#pm').onclick=pres;$('#px').onclick=pres;
$('#bl').onclick=blender;$('#dl').onclick=()=>mol&&save(slug()+'.'+(mol.raw?mol.ext:'xyz'),mol.raw||xyz(mol.atoms));
$('#file').onchange=e=>{addFile(e.target.files[0]);e.target.value=''};
vp.ondragover=e=>e.preventDefault();vp.ondrop=e=>{e.preventDefault();addFile(e.dataTransfer.files[0])};
$('#fetch').onclick=()=>{const id=$('#pid').value.trim().toUpperCase();if(!/^[0-9][A-Z0-9]{3}$/.test(id))return toast('Enter a valid 4-character PDB ID.');fetchId(id)};
addEventListener('keydown',e=>{const mo=$('#modal').style.display==='flex'||IM.classList.contains('open')||$('#lightbox').classList.contains('open');if(!Core.shortcutAllowed(e,{modalOpen:mo}))return;const k=e.key.toLowerCase();
if(k==='escape'){if($('#lightbox').classList.contains('open'))$('#lightbox').classList.remove('open');else if(IM.classList.contains('open'))imOpen(false);else if($('#modal').style.display==='flex')modal();else if(document.body.classList.contains('pm'))pres();return}
if(k==='f')toggleFS();else if(k==='r')fit();else if(k==='z')focusSel();else if(k==='d')show2D();else if(k==='h')setH(!$('#hh').checked);else if(k===' '){e.preventDefault();$('#ar').click()}else if(k==='p')pres();else if(k==='t')$('#th').click()});
try{const t=localStorage.getItem('laura-theme');if(t)document.documentElement.dataset.theme=t}catch(e){}
try{const v=localStorage.getItem('laura-sort');if(v)$('#srt').value=v}catch(e){}
setBg();size();thBtn();cats();renderList();load(LIB[4]);
