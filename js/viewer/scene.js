// js/viewer/scene.js: three.js scene: renderer, camera, orbit controls, lighting, and the shared viewer state (mol, sel, vis, ...).
// Also resize handling and the fit / focus helpers.
// ---------- measurements ----------
const dist=Core.distance,angle=Core.angle,dihedral=Core.dihedral;
// ---------- three.js scene ----------
const vp=$('#vp'),scene=new THREE.Scene(),cam=new THREE.PerspectiveCamera(40,1,.1,5000);let ren;
try{ren=new THREE.WebGLRenderer({antialias:true});}catch(e){vp.insertAdjacentHTML('beforeend','<div style="padding:60px 20px">WebGL is not available in this browser. Try another browser or enable hardware acceleration.</div>')}
if(ren){vp.prepend(ren.domElement);ren.setPixelRatio(Math.min(devicePixelRatio,2))}
if(ren){ren.domElement.addEventListener('webglcontextlost',e=>{e.preventDefault();toast('The graphics context was lost (the browser reclaimed the GPU). LAURA will redraw when it is restored.',9000)});ren.domElement.addEventListener('webglcontextrestored',()=>{if(mol){build();drawSel()}toast('Graphics restored.')})}
scene.add(new THREE.AmbientLight(0xffffff,.65));const dl=new THREE.DirectionalLight(0xffffff,.8);dl.position.set(1,1.5,2);cam.add(dl);scene.add(cam);
const ctl=ren&&new THREE.OrbitControls(cam,ren.domElement);if(ctl){ctl.enableDamping=true;ctl.autoRotateSpeed=1}
const grp=new THREE.Group();scene.add(grp);const ovg=new THREE.Group();scene.add(ovg);
let mol=null,cen=new THREE.Vector3(),rad=5,vis=[],visBase=[],chainLab=[],cLab=[],cmpTags=[],aMesh=null,sel=[],selG=new THREE.Group(),tok=0,ro={},hl=[],lastC=null;grp.add(selG);
function setBg(){scene.background=new THREE.Color(getComputedStyle(vp).backgroundColor)}
function size(){if(!ren)return;const w=vp.clientWidth,h=vp.clientHeight;ren.setSize(w,h);cam.aspect=w/h;cam.updateProjectionMatrix()}
new ResizeObserver(size).observe(vp);
function focusSel(){if(!mol)return;if(!sel.length){toast('Select atoms first, then use Zoom to selection.');return}
const F=Core.framing(sel.map(i=>{const a=mol.atoms[i];return{x:a.x-cen.x,y:a.y-cen.y,z:a.z-cen.z}}),{fov:cam.fov,aspect:cam.aspect,minRadius:2.5}),T=ctl?ctl.target:new THREE.Vector3(),dir=cam.position.clone().sub(T);if(dir.length()<1e-6)dir.set(0,0,1);dir.normalize();
cam.position.set(F.center.x+dir.x*F.distance,F.center.y+dir.y*F.distance,F.center.z+dir.z*F.distance);cam.near=Math.max(.01,F.distance/100);cam.far=F.distance*50+rad*4;cam.updateProjectionMatrix();if(ctl){ctl.target.set(F.center.x,F.center.y,F.center.z);ctl.update()}}
function fit(){const d=rad/Math.sin(cam.fov*Math.PI/360)*1.15;cam.position.set(0,0,d);cam.near=d/100;cam.far=d*20;cam.updateProjectionMatrix();if(ctl){ctl.target.set(0,0,0);ctl.update()}}
function clearGrp(){for(const o of[...grp.children]){if(o===selG)continue;grp.remove(o);o.dispose&&o.dispose();o.geometry&&o.geometry.dispose();o.material&&o.material.dispose()}aMesh=null}
