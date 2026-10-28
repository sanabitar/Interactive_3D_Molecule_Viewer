// Geometry. Inputs {x,y,z} in Angstrom; angles returned in degrees.
const sub=(a,b)=>[b.x-a.x,b.y-a.y,b.z-a.z],dot=(u,v)=>u[0]*v[0]+u[1]*v[1]+u[2]*v[2],cross=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],len=u=>Math.hypot(...u);
export const distance=(a,b)=>len(sub(a,b));
export function angle(a,b,c){const u=sub(b,a),v=sub(b,c),m=len(u)*len(v);if(m===0)throw new RangeError('Angle undefined: coincident atoms.');return Math.acos(Math.max(-1,Math.min(1,dot(u,v)/m)))*180/Math.PI}
// IUPAC sign convention; range (-180, 180]
export function dihedral(a,b,c,d){const b1=sub(a,b),b2=sub(b,c),b3=sub(c,d);if(!len(b2)||!len(cross(b1,b2))||!len(cross(b2,b3)))throw new RangeError('Dihedral undefined: collinear or coincident atoms.');
 return Math.atan2(len(b2)*dot(b1,cross(b2,b3)),dot(cross(b1,b2),cross(b2,b3)))*180/Math.PI}
