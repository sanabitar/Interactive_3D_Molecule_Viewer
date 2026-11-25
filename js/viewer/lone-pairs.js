// js/viewer/lone-pairs.js: Lone electron pairs drawn as translucent purple lobes with two electron dots (not atoms: not in mol.atoms, not pickable, not counted).
// Counts and directions come from LonePairs (js/generated/lonepairs.js); drawLonePairs() is called from build() (Ball and stick view) when "Show lone pairs" is on.
const LP_COLOR=0xb388ff,LP_MAX_ATOMS=6000;let lpShown=0;   // lpShown: lone pairs drawn by the last build (-1: structure too large)
// Computed once per molecule, from the same bond orders and charges the 2D diagram uses (carboxylate, protonated acid, N+).
function lonePairData(){
 if(mol._lp)return mol._lp;
 if(mol.entry&&mol.entry.noH){mol._lp=[];return mol._lp}   // schematic models (B-DNA) have no real atoms
 const dr=Layout2D.drawnOrdersAndCharges(mol.atoms,mol.bonds,mol.orders);
 mol._lp=LonePairs.computeLonePairs(mol.atoms,mol.bonds,dr.orders,{charges:dr.charge});return mol._lp}
// inS: Set of the atom indices drawn; R(a): the sphere radius of an atom. The lobe starts at the sphere surface and points along the lone-pair direction.
function drawLonePairs(inS,R){lpShown=0;if(!$('#lp').checked||!mol)return;
 if(inS.size>LP_MAX_ATOMS){lpShown=-1;return}
 const A=mol.atoms,items=lonePairData().filter(l=>inS.has(l.atom));lpShown=items.length;if(!items.length)return;
 const lobes=new THREE.InstancedMesh(new THREE.SphereGeometry(1,18,12),new THREE.MeshStandardMaterial({color:LP_COLOR,transparent:true,opacity:.5,roughness:.35,depthWrite:false}),items.length),
  dots=new THREE.InstancedMesh(new THREE.SphereGeometry(1,8,6),new THREE.MeshBasicMaterial({color:0xffffff}),items.length*2),o=new THREE.Object3D(),zUp=new THREE.Vector3(0,0,1);let n=0;
 items.forEach((l,k)=>{const a=A[l.atom],d=l.dir,L=R(a)+.3,cx=a.x-cen.x+d[0]*L,cy=a.y-cen.y+d[1]*L,cz=a.z-cen.z+d[2]*L;
  o.quaternion.setFromUnitVectors(zUp,new THREE.Vector3(d[0],d[1],d[2]));o.position.set(cx,cy,cz);o.scale.set(.21,.21,.34);o.updateMatrix();lobes.setMatrixAt(k,o.matrix);
  let s=Math.abs(d[1])<.9?[-d[2],0,d[0]]:[0,d[2],-d[1]];const sl=Math.hypot(...s)||1;s=[s[0]/sl,s[1]/sl,s[2]/sl];   // across the lobe: the two electrons of the pair side by side
  for(const q of[-1,1]){o.quaternion.identity();o.position.set(cx+q*s[0]*.085,cy+q*s[1]*.085,cz+q*s[2]*.085);o.scale.setScalar(.065);o.updateMatrix();dots.setMatrixAt(n++,o.matrix)}});
 lobes.userData.lp=items.length;dots.userData.lpDots=items.length*2;grp.add(lobes);grp.add(dots)}
// Feedback when the user turns the option on and nothing can be shown.
function lonePairNotice(){if(!$('#lp').checked||!mol)return;
 if(mol.entry&&mol.entry.noH)return toast(mol.entry.noH,7000);
 if(ro&&ro.st!=='bs')return toast('Lone pairs are drawn in Ball and stick view; Chains view has no atoms to attach them to.',7000);
 if(lpShown===-1)return toast('Lone pairs are drawn for structures up to '+LP_MAX_ATOMS+' atoms.',7000);
 if(lpShown>0)return;
 if(!lonePairData().length)return toast(mol.atoms.some(a=>a.el==='H')?'This structure has no lone pairs.':'No lone pairs shown: this structure has no hydrogens, so N and O atoms that look short of bonds (N-H, O-H) are skipped rather than guessed.',8000);
 toast('None of the atoms shown has a lone pair.',6000)}
