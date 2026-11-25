// js/tools/teaching.js: The guided "measure water" activity.
function tch(h){const t=$('#tch');t.style.display=h?'block':'none';t.innerHTML=h||''}
function teach(){tch(`<b>Guided activity: measure water</b><ol><li><button id="t1">Load water</button></li><li>Click <b>H, then O, then H</b> in the viewer (O is the vertex). This card stays open while you click in the viewer; click anywhere else to close it.</li><li><button id="t2">Check my measurement</button></li></ol><div id="tm"></div><p><button id="tx">Close</button></p>`);
$('#t1').onclick=()=>load(LIB.find(l=>l.id==='water'));$('#tx').onclick=()=>tch();
$('#t2').onclick=()=>{const P=mol?sel.map(i=>mol.atoms[i]):[];$('#tm').textContent=!mol||mol.entry.id!=='water'?'Load water first.':P.length===3&&P[1].el==='O'?`Measured ${angle(P[0],P[1],P[2]).toFixed(2)}° (the model is built from the 104.5° reference, so this checks your selection technique).`:P.length===2?`Distance ${dist(P[0],P[1]).toFixed(3)} Å (model O–H 0.957 Å).`:'Select exactly H, O, H in that order.'}}

