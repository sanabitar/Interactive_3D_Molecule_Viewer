// js/viewer/molecule-render.js: build(): turns the current molecule into spheres, bonds and (in Chains view) tubes, applying the display options.
function build(){syncHBtn();clearGrp();chainLab=[];if(!mol)return;const A=mol.atoms,hh=$('#hh').checked,ww=$('#ww').checked;let st=$('#sty').value==='tb'?'tb':'bs';
const ok=i=>{const a=A[i];return(hh||a.el!=='H')&&(ww||a.res!=='HOH')};vis=[];A.forEach((a,i)=>ok(i)&&vis.push(i));{const vs=new Set(vis);sel=sel.filter(i=>vs.has(i))}visBase=vis;if($('#so').checked&&sel.length)vis=sel.slice();
const hasChains=Object.values(chainTraces()).some(p=>p.length>=4);
if(st==='tb'&&!hasChains){toast('Chains view needs a polymer chain with at least 4 CA/P atoms; showing ball and stick.');st='bs'}
if(st==='bs'&&vis.length>40000&&hasChains){toast('Very large structure: switched to Chains view for performance.');st='tb'}
$('#sty').value=st;syncStyBtn();const sub=vis;
const inS=new Set(sub),useBO=$('#bo').checked;let B=mol.bonds.map(([i,j],k)=>[i,j,useBO?mol.orders[k]:1]).filter(([i,j])=>inS.has(i)&&inS.has(j));
const R={bs:a=>vdw(a.el)*.25,tb:()=>.4}[st],br={bs:.1}[st];ro={st,R};
const geo=new THREE.SphereGeometry(1,14,10),m=new THREE.MeshStandardMaterial({color:0xffffff,roughness:.45,visible:st!=='tb'});
aMesh=new THREE.InstancedMesh(geo,m,sub.length);const o=new THREE.Object3D(),c=new THREE.Color();aMesh.userData.map=sub;
sub.forEach((i,k)=>{const a=A[i],r=R(a);o.position.set(a.x-cen.x,a.y-cen.y,a.z-cen.z);o.scale.setScalar(r);o.quaternion.identity();o.updateMatrix();aMesh.setMatrixAt(k,o.matrix);aMesh.setColorAt(k,c.setHex(cpk(a.el)))});grp.add(aMesh);
if(st==='tb')tubes();else{const SH={1:[[0,1]],2:[[-.5,1],[.5,1]],3:[[-1,1],[0,1],[1,1]],1.5:[[-.5,1],[.5,.55]]},sp=(br||.1)*2.6,seg=[];
for(const[i,j,ob]of B){const sh=SH[ob]||SH[1],u=sh.length>1?offDir(i,j):[0,0,0];for(const[s0,rs]of sh)seg.push([i,j,[u[0]*s0*sp,u[1]*s0*sp,u[2]*s0*sp],rs])}
if(seg.length){const cm=new THREE.InstancedMesh(new THREE.CylinderGeometry(1,1,1,10),new THREE.MeshStandardMaterial({color:0xffffff,roughness:.45}),seg.length*2),up=new THREE.Vector3(0,1,0);let n=0;
for(const[i,j,f,rs]of seg){const a=A[i],b=A[j],p=new THREE.Vector3(a.x+f[0],a.y+f[1],a.z+f[2]).sub(cen),q=new THREE.Vector3(b.x+f[0],b.y+f[1],b.z+f[2]).sub(cen),d=q.clone().sub(p),L=d.length(),mid=p.clone().add(q).multiplyScalar(.5);o.quaternion.setFromUnitVectors(up,d.normalize());o.scale.set((br||.1)*rs,L/2,(br||.1)*rs);
o.position.copy(p.clone().add(mid).multiplyScalar(.5));o.updateMatrix();cm.setMatrixAt(n,o.matrix);cm.setColorAt(n++,c.setHex(cpk(a.el)));o.position.copy(mid.clone().add(q).multiplyScalar(.5));o.updateMatrix();cm.setMatrixAt(n,o.matrix);cm.setColorAt(n++,c.setHex(cpk(b.el)))}grp.add(cm)}}
drawSel()}
function offDir(i,j){return offDirG(mol.atoms,mol.adj,i,j)}
function offDirG(A,adj,i,j){const a=A[i],b=A[j],d=[b.x-a.x,b.y-a.y,b.z-a.z],L=Math.hypot(...d);d[0]/=L;d[1]/=L;d[2]/=L;
for(const[c,o]of[[i,a],[j,b]])for(const k of adj[c]){if(k===i||k===j)continue;const t=A[k],v=[t.x-o.x,t.y-o.y,t.z-o.z],n=[d[1]*v[2]-d[2]*v[1],d[2]*v[0]-d[0]*v[2],d[0]*v[1]-d[1]*v[0]];if(Math.hypot(...n)>.3){const w=[n[1]*d[2]-n[2]*d[1],n[2]*d[0]-n[0]*d[2],n[0]*d[1]-n[1]*d[0]],q=Math.hypot(...w);return[w[0]/q,w[1]/q,w[2]/q]}}
const r=Math.abs(d[0])<.9?[1,0,0]:[0,1,0],n=[d[1]*r[2]-d[2]*r[1],d[2]*r[0]-d[0]*r[2],d[0]*r[1]-d[1]*r[0]],q=Math.hypot(...n);return[n[0]/q,n[1]/q,n[2]/q]}
