// Least-squares rigid superposition (Horn 1987 quaternion method == Kabsch optimum). Equal weights, rotation+translation only (no reflection).
import {distance} from './geometry.mjs';
function eig4(N){const A=N.map(r=>[...r]),V=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
 for(let sw=0;sw<80;sw++){let off=0;for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)off+=A[i][j]**2;if(off<1e-26)break;
  for(let p=0;p<3;p++)for(let q=p+1;q<4;q++){if(Math.abs(A[p][q])<1e-300)continue;const th=(A[q][q]-A[p][p])/(2*A[p][q]),t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),s=t*c;
   for(let k=0;k<4;k++){const x=A[k][p],y=A[k][q];A[k][p]=c*x-s*y;A[k][q]=s*x+c*y}
   for(let k=0;k<4;k++){const x=A[p][k],y=A[q][k];A[p][k]=c*x-s*y;A[q][k]=s*x+c*y}
   for(let k=0;k<4;k++){const x=V[k][p],y=V[k][q];V[k][p]=c*x-s*y;V[k][q]=s*x+c*y}}}
 return{values:[0,1,2,3].map(i=>A[i][i]),vectors:V}}
// Finds the rotation+translation that moves P (moving) onto Q (reference). P[i] pairs with Q[i].
export function superpose(P,Q){const n=P.length;if(n!==Q.length)throw new RangeError('Coordinate sets differ in length.');if(n<3)throw new RangeError('At least 3 coordinate pairs are required.');
 const cen=X=>({x:X.reduce((s,a)=>s+a.x,0)/n,y:X.reduce((s,a)=>s+a.y,0)/n,z:X.reduce((s,a)=>s+a.z,0)/n}),cp=cen(P),cq=cen(Q);
 let Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0,G=0;
 for(let i=0;i<n;i++){const a=[P[i].x-cp.x,P[i].y-cp.y,P[i].z-cp.z],b=[Q[i].x-cq.x,Q[i].y-cq.y,Q[i].z-cq.z];
  Sxx+=a[0]*b[0];Sxy+=a[0]*b[1];Sxz+=a[0]*b[2];Syx+=a[1]*b[0];Syy+=a[1]*b[1];Syz+=a[1]*b[2];Szx+=a[2]*b[0];Szy+=a[2]*b[1];Szz+=a[2]*b[2];G+=a[0]**2+a[1]**2+a[2]**2+b[0]**2+b[1]**2+b[2]**2}
 const N=[[Sxx+Syy+Szz,Syz-Szy,Szx-Sxz,Sxy-Syx],[Syz-Szy,Sxx-Syy-Szz,Sxy+Syx,Szx+Sxz],[Szx-Sxz,Sxy+Syx,-Sxx+Syy-Szz,Syz+Szy],[Sxy-Syx,Szx+Sxz,Syz+Szy,-Sxx-Syy+Szz]];
 const{values,vectors}=eig4(N);let m=0;for(let k=1;k<4;k++)if(values[k]>values[m])m=k;
 let[w,x,y,z]=[0,1,2,3].map(r=>vectors[r][m]);const L=Math.hypot(w,x,y,z);w/=L;x/=L;y/=L;z/=L;
 const R=[[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],[2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],[2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]];
 const apply=p=>{const a=[p.x-cp.x,p.y-cp.y,p.z-cp.z];return{x:R[0][0]*a[0]+R[0][1]*a[1]+R[0][2]*a[2]+cq.x,y:R[1][0]*a[0]+R[1][1]*a[1]+R[1][2]*a[2]+cq.y,z:R[2][0]*a[0]+R[2][1]*a[1]+R[2][2]*a[2]+cq.z}};
 const deviations=P.map((p,i)=>distance(apply(p),Q[i]));
 return{rmsd:Math.sqrt(deviations.reduce((s,d)=>s+d*d,0)/n),rmsdFromEigenvalue:Math.sqrt(Math.max(0,(G-2*values[m])/n)),rotation:R,apply,deviations,n}}
