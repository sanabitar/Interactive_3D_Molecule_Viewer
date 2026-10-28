// Camera framing: centre and distance that fit a set of points inside the narrower of the vertical/horizontal field of view.
export function framing(points,{fov=40,aspect=1,minRadius=1.5,margin=1.15}={}){if(!points||!points.length)return null;const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity],K=['x','y','z'];
 for(const p of points)K.forEach((k,i)=>{if(!Number.isFinite(p[k]))throw new RangeError('Non-finite coordinate in framing input.');lo[i]=Math.min(lo[i],p[k]);hi[i]=Math.max(hi[i],p[k])});
 const c={x:(lo[0]+hi[0])/2,y:(lo[1]+hi[1])/2,z:(lo[2]+hi[2])/2};let r=minRadius;for(const p of points)r=Math.max(r,Math.hypot(p.x-c.x,p.y-c.y,p.z-c.z));
 const v=fov*Math.PI/360,h=Math.atan(Math.tan(v)*aspect);return{center:c,radius:r,distance:margin*r/Math.sin(Math.min(v,h))}}
