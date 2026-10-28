const atom=(el,x,y,z)=>({el,x,y,z});
export function electronModel(id){
 const origin=el=>atom(el,0,0,0),r=Math.PI/180;let model;
 if(id==='water'){const a=52.25*r;model={name:'Water · H₂O',atoms:[origin('O'),atom('H',Math.sin(a)*.96,0,-Math.cos(a)*.96),atom('H',-Math.sin(a)*.96,0,-Math.cos(a)*.96)],bonds:[[0,1,1],[0,2,1]],lonePairs:[{atom:0,direction:[0,.82,.57]},{atom:0,direction:[0,-.82,.57]}],geometry:'Bent',electronGeometry:'Tetrahedral',angle:104.5};}
 else if(id==='ammonia'){const c=Math.sqrt((Math.cos(107*r)+.5)/1.5),s=Math.sqrt(1-c*c);model={name:'Ammonia · NH₃',atoms:[origin('N'),...Array.from({length:3},(_,i)=>atom('H',1.01*s*Math.cos(i*2*Math.PI/3),1.01*s*Math.sin(i*2*Math.PI/3),-1.01*c))],bonds:[[0,1,1],[0,2,1],[0,3,1]],lonePairs:[{atom:0,direction:[0,0,1]}],geometry:'Trigonal pyramidal',electronGeometry:'Tetrahedral',angle:107};}
 else if(id==='methane'){const s=1.09/Math.sqrt(3);model={name:'Methane · CH₄',atoms:[origin('C'),...[ [1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(p=>atom('H',...p.map(v=>v*s)))],bonds:[[0,1,1],[0,2,1],[0,3,1],[0,4,1]],lonePairs:[],geometry:'Tetrahedral',electronGeometry:'Tetrahedral',angle:109.47};}
 else if(id==='co2')model={name:'Carbon dioxide · CO₂',atoms:[origin('C'),atom('O',-1.16,0,0),atom('O',1.16,0,0)],bonds:[[0,1,2],[0,2,2]],lonePairs:[{atom:1,direction:[-.6,.8,0]},{atom:1,direction:[-.6,-.8,0]},{atom:2,direction:[.6,.8,0]},{atom:2,direction:[.6,-.8,0]}],geometry:'Linear',electronGeometry:'Linear',angle:180};
 else throw Error('Unknown electron teaching model.');
 return {...model,central:0,id};
}
export function electronCounts(model){const bondingPairs=model.bonds.reduce((s,b)=>s+b[2],0),lonePairs=model.lonePairs.length,centralLonePairs=model.lonePairs.filter(p=>p.atom===model.central).length,centralBonds=model.bonds.filter(b=>b[0]===model.central||b[1]===model.central).length;return {bondingPairs,lonePairs,centralLonePairs,domains:centralBonds+centralLonePairs,valenceElectrons:2*(bondingPairs+lonePairs)};}
