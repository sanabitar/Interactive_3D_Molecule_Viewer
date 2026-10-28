export const SUGARS={glucose:['down','up','down'],mannose:['up','up','down'],galactose:['down','up','up']};
export function chairModel(sugar='glucose',anomer='beta',flipped=false){
 if(!SUGARS[sugar]||!['alpha','beta'].includes(anomer))throw Error('Unsupported sugar or anomer.');
 const signs=[anomer==='alpha'?-1:1,...SUGARS[sugar].map(v=>v==='up'?1:-1),1];
 const ring=Array.from({length:6},(_,i)=>{const t=i*Math.PI/3;return {x:1.45*Math.cos(t),y:1.45*Math.sin(t),z:(i%2?1:-1)*.25*(flipped?-1:1),label:i===5?'O5':'C'+(i+1)};});
 const groups=signs.map((sign,i)=>{const axialSign=(i%2?1:-1)*(flipped?-1:1),axial=sign===axialSign,t=i*Math.PI/3,p=ring[i],r=axial?0:1.2;return {carbon:i+1,label:i===4?'CH₂OH':'OH',direction:sign===1?'up':'down',position:axial?'axial':'equatorial',x:p.x+r*Math.cos(t),y:p.y+r*Math.sin(t),z:p.z+sign*(axial?1.25:.48)};});
 return {ring,groups,chair:flipped?'¹C₄':'⁴C₁',name:(anomer==='alpha'?'α':'β')+'-D-'+({glucose:'glucopyranose',mannose:'mannopyranose',galactose:'galactopyranose'}[sugar])};
}
