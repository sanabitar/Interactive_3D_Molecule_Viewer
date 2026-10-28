// Flat "2D view" of 3D coordinates: orthographic projection onto the best-fit plane (principal-component analysis).
// This is NOT a chemical structure diagram: it preserves distances only for planar molecules; non-planar atoms overlap.
function eigSymN(N){const n=N.length,A=N.map(r=>[...r]),V=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
 for(let sw=0;sw<80;sw++){let off=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)off+=A[i][j]**2;if(off<1e-28)break;
  for(let p=0;p<n-1;p++)for(let q=p+1;q<n;q++){if(Math.abs(A[p][q])<1e-300)continue;const th=(A[q][q]-A[p][p])/(2*A[p][q]),t=(th>=0?1:-1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),s=t*c;
   for(let k=0;k<n;k++){const x=A[k][p],y=A[k][q];A[k][p]=c*x-s*y;A[k][q]=s*x+c*y}for(let k=0;k<n;k++){const x=A[p][k],y=A[q][k];A[p][k]=c*x-s*y;A[q][k]=s*x+c*y}for(let k=0;k<n;k++){const x=V[k][p],y=V[k][q];V[k][p]=c*x-s*y;V[k][q]=s*x+c*y}}}
 return{values:A.map((r,i)=>r[i]),vectors:V}}
export function projectPlane(points){const n=points.length;if(!n)return{xy:[],depth:[],variance:[0,0,0]};const c=['x','y','z'].map(k=>points.reduce((s,p)=>s+p[k],0)/n),C=[[0,0,0],[0,0,0],[0,0,0]];
 for(const p of points){const d=[p.x-c[0],p.y-c[1],p.z-c[2]];for(let i=0;i<3;i++)for(let j=0;j<3;j++)C[i][j]+=d[i]*d[j]/n}
 const{values,vectors}=eigSymN(C),ord=[0,1,2].sort((a,b)=>values[b]-values[a]),ax=ord.map(k=>[vectors[0][k],vectors[1][k],vectors[2][k]]);
 const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],xy=[],depth=[];
 for(const p of points){const d=[p.x-c[0],p.y-c[1],p.z-c[2]];xy.push([dot(d,ax[0]),dot(d,ax[1])]);depth.push(dot(d,ax[2]))}
 return{xy,depth,variance:ord.map(k=>Math.max(0,values[k]))}}
