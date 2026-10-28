// Recording stand-in for Three.js r128 (NO rendering). Used only because no Three.js build exists offline.
// Verifies the page's DOM/CSS/JS/logic in a real browser and counts geometry/material/texture lifecycle; it CANNOT verify GPU output.
(function(){const L={geo:0,mat:0,tex:0,inst:0,created:{geo:0,mat:0,tex:0,inst:0},disposed:{geo:0,mat:0,tex:0,inst:0},renders:0,insts:[]};window.__three=L;
class V{constructor(x=0,y=0,z=0){this.x=x;this.y=y;this.z=z}set(x,y,z){this.x=x;this.y=y;this.z=z;return this}copy(o){this.x=o.x;this.y=o.y;this.z=o.z;return this}add(o){this.x+=o.x;this.y+=o.y;this.z+=o.z;return this}sub(o){this.x-=o.x;this.y-=o.y;this.z-=o.z;return this}multiplyScalar(s){this.x*=s;this.y*=s;this.z*=s;return this}clone(){return new V(this.x,this.y,this.z)}length(){return Math.hypot(this.x,this.y,this.z)}normalize(){return this.multiplyScalar(1/(this.length()||1))}negate(){return this.multiplyScalar(-1)}toArray(){return[this.x,this.y,this.z]}fromArray(a){return this.set(a[0],a[1],a[2])}project(c){const d=Math.max(Math.hypot(c.position.x,c.position.y,c.position.z),1e-6),f=1/(d*Math.tan((c.fov||40)*Math.PI/360));this.x=this.x*f/(c.aspect||1);this.y=this.y*f;this.z=0;return this}}
const life=k=>{L.created[k]++;L[k]++};const dis=(o,k)=>{if(!o._d){o._d=1;L.disposed[k]++;L[k]--}};
class Geo{constructor(){life('geo');this.attributes={};this.boundingSphere=null}setAttribute(k,v){this.attributes[k]=v;return this}setFromPoints(p){this.attributes.position={count:p.length};return this}dispose(){dis(this,'geo')}}
class Mat{constructor(o={}){life('mat');Object.assign(this,{visible:true,map:null},o)}dispose(){dis(this,'mat')}}
class Tex{constructor(){life('tex')}dispose(){dis(this,'tex')}}
class Obj{constructor(){this.position=new V();this.scale={set(){},setScalar(){}};this.quaternion={identity(){},setFromUnitVectors(){}};this.matrix={};this.userData={};this.children=[];this.visible=true}updateMatrix(){this.matrix={p:this.position.toArray()}}add(...o){this.children.push(...o)}remove(o){this.children=this.children.filter(c=>c!==o)}}
class Mesh extends Obj{constructor(g,m){super();this.geometry=g;this.material=m}}
class Inst extends Mesh{constructor(g,m,n){super(g,m);this.count=n;life('inst');this.pos=[];L.insts.push(this);this.instanceMatrix={};this.instanceColor=null}setMatrixAt(i,m){this.pos[i]=m&&m.p}setColorAt(){this.instanceColor={}}dispose(){dis(this,'inst')}}
const sg=new Geo();class Sprite extends Obj{constructor(m){super();this.material=m;this.geometry=sg}}
class Color{constructor(){}setHex(){return this}setHSL(){return this}toArray(){return[1,1,1]}}
class Cam extends Obj{constructor(f){super();this.fov=f||40;this.aspect=1;this.near=.1;this.far=1e3}updateProjectionMatrix(){}}
class Light extends Obj{}
class Ctl{constructor(c,d){this.target=new V();this.autoRotate=false;this.autoRotateSpeed=2;this.enableDamping=false}update(){return this.autoRotate}}
class Ren{constructor(){this.domElement=document.createElement('canvas')}setPixelRatio(){}setSize(w,h){this.domElement.width=w;this.domElement.height=h}render(){L.renders++}}
window.THREE={Vector3:V,Vector2:V,Color,Scene:class extends Obj{},PerspectiveCamera:Cam,WebGLRenderer:Ren,OrbitControls:Ctl,Group:class extends Obj{},Object3D:Obj,AmbientLight:Light,DirectionalLight:Light,
InstancedMesh:Inst,Mesh,LineSegments:Mesh,Points:Mesh,Sprite,SphereGeometry:Geo,CylinderGeometry:Geo,BufferGeometry:Geo,TubeGeometry:Geo,CatmullRomCurve3:class{constructor(p){this.p=p}},Float32BufferAttribute:class{constructor(a,s){this.count=a.length/s}},
MeshStandardMaterial:Mat,MeshBasicMaterial:Mat,LineBasicMaterial:Mat,PointsMaterial:Mat,SpriteMaterial:Mat,CanvasTexture:Tex,Raycaster:class{setFromCamera(){}intersectObject(){return[]}}};})();
