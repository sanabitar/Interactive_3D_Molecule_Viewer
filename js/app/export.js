// js/app/export.js: Downloads: XYZ and the Blender package (including the generated Blender import script).
// ---------- export ----------
function save(name,data,type='text/plain'){const a=document.createElement('a'),u=URL.createObjectURL(data instanceof Blob?data:new Blob([data],{type}));a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),2000)}
const xyz=A=>A.length+'\nLAURA export (Angstrom)\n'+A.map(a=>`${a.el} ${a.x.toFixed(4)} ${a.y.toFixed(4)} ${a.z.toFixed(4)}`).join('\n')+'\n';
const slug=()=>mol.entry.name.replace(/[^\w]+/g,'_').replace(/^_|_$/g,'');
const PY=`import bpy, json
from mathutils import Vector
# LAURA Blender importer. Adds ONE new collection; does not touch other scene objects.
NAME="__NAME__"
DATA=json.loads(r'''__DATA__''')  # coordinates in Angstrom, centred on the centroid; 1 Angstrom = 1 Blender unit
COL={"H":(1,1,1),"C":(.55,.56,.6),"N":(.19,.31,.97),"O":(1,.16,.16),"P":(1,.5,0),"S":(1,.84,.16)}
VDW={"H":1.2,"C":1.7,"N":1.55,"O":1.52,"P":1.8,"S":1.8}
col=bpy.data.collections.new("LAURA_"+NAME); bpy.context.scene.collection.children.link(col)
def mat(n,c):
    m=bpy.data.materials.new("LAURA_"+n); m.diffuse_color=(*c,1); m.use_nodes=True
    m.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value=(*c,1); return m
bpy.ops.mesh.primitive_uv_sphere_add(radius=1,segments=24,ring_count=12); t=bpy.context.object; sm=t.data; bpy.data.objects.remove(t)
bpy.ops.mesh.primitive_cylinder_add(radius=1,depth=1,vertices=16); t=bpy.context.object; cm=t.data; bpy.data.objects.remove(t)
cm.materials.append(mat("bond",(.6,.6,.6)))
meshes={}
A=DATA["a"]
for i,(el,x,y,z) in enumerate(A):
    if el not in meshes:
        me=sm.copy(); me.name="LAURA_"+el; me.materials.append(mat(el,COL.get(el,(.8,.4,.8)))); meshes[el]=me
    o=bpy.data.objects.new("%s_%d"%(el,i+1),meshes[el]); o.location=(x,y,z); r=VDW.get(el,1.7)*0.25; o.scale=(r,r,r); col.objects.link(o)
adj={}
for i,j,od in DATA["b"]:
    adj.setdefault(i,[]).append(j); adj.setdefault(j,[]).append(i)
def off(i,j):  # in-plane perpendicular for parallel double/triple bond cylinders
    p=Vector(A[i][1:]); q=Vector(A[j][1:]); d=(q-p).normalized()
    for c,o in ((i,p),(j,q)):
        for k in adj.get(c,[]):
            if k in (i,j): continue
            n=d.cross(Vector(A[k][1:])-o)
            if n.length>0.3: return n.cross(d).normalized()
    ref=Vector((1,0,0)) if abs(d.x)<.9 else Vector((0,1,0))
    return d.cross(ref).normalized()
SH={1:[(0,1)],2:[(-.5,1),(.5,1)],3:[(-1,1),(0,1),(1,1)],1.5:[(-.5,1),(.5,.55)]}  # single, double, triple, aromatic/delocalised
for k,(i,j,od) in enumerate(DATA["b"]):
    p=Vector(A[i][1:]); q=Vector(A[j][1:]); d=q-p
    sh=SH.get(od,SH[1]); u=off(i,j) if len(sh)>1 else Vector((0,0,0))
    for n,(s0,rs) in enumerate(sh):
        o=bpy.data.objects.new("Bond_%d_%d"%(k,n),cm); o.location=(p+q)/2+u*s0*0.26; o.scale=(.1*rs,.1*rs,d.length)
        o.rotation_mode="QUATERNION"; o.rotation_quaternion=Vector((0,0,1)).rotation_difference(d); col.objects.link(o)
R=max(Vector(a[1:]).length for a in A)+2
cam=bpy.data.objects.new("LAURA_Camera",bpy.data.cameras.new("LAURA_Camera")); cam.location=(0,-R*3,R*.8); col.objects.link(cam)
tgt=bpy.data.objects.new("LAURA_Target",None); col.objects.link(tgt)
c=cam.constraints.new("TRACK_TO"); c.target=tgt; c.track_axis="TRACK_NEGATIVE_Z"; c.up_axis="UP_Y"
sun=bpy.data.objects.new("LAURA_Sun",bpy.data.lights.new("LAURA_Sun","SUN")); sun.location=(R,-R,R*2); sun.rotation_euler=(.8,.2,.6); col.objects.link(sun)
bpy.context.scene.camera=cam
print("LAURA: imported",len(A),"atoms and",len(DATA["b"]),"bonds")
`;
async function blender(){if(!mol)return;const idx=new Map(),a=[],b=[];vis.forEach(i=>{idx.set(i,a.length);const t=mol.atoms[i];a.push([t.el,+(t.x-cen.x).toFixed(4),+(t.y-cen.y).toFixed(4),+(t.z-cen.z).toFixed(4)])});
if(a.length>8000)return toast('Too large for the per-atom Blender script ('+a.length+' atoms; limit 8000). Download the structure file and use a protein add-on such as Molecular Nodes, or hide solvent/hydrogens first.',8000);
mol.bonds.forEach(([i,j],k)=>idx.has(i)&&idx.has(j)&&b.push([idx.get(i),idx.get(j),$('#bo').checked?mol.orders[k]:1]));const nm=slug(),py=PY.replace('__NAME__',nm).replace('__DATA__',JSON.stringify({a,b}));
const raw=mol.raw||xyz(mol.atoms),ext=mol.raw?mol.ext:'xyz',meta=JSON.stringify({name:mol.entry.name,category:mol.entry.cat,source:mol.entry.src,format:mol.entry.fmt,atoms:mol.atoms.length,exportedAtoms:a.length,bondsNote:'distance-derived'},null,1);
const rd=`LAURA Blender package: ${mol.entry.name}\n\n1. Open Blender (tested path not verified by LAURA authors; script uses standard bpy API).\n2. Scripting workspace > Open > import_molecule.py > Run Script.\n3. A new collection LAURA_${nm} with atoms, bonds, camera and sun is added; other objects are untouched.\n\nNotes: atoms are spheres at 0.25 x van der Waals radius; bonds are inferred from distances; bond orders (single/double/triple/aromatic) are inferred from bond lengths and drawn as parallel cylinders. The original structure file is included; the Blender scene is a visual model only.\n`;
if(window.JSZip){const z=new JSZip();z.file(nm+'.'+ext,raw);z.file('import_molecule.py',py);z.file('metadata.json',meta);z.file('README.txt',rd);save('LAURA_'+nm+'_blender.zip',await z.generateAsync({type:'blob'}))}else{save('import_molecule.py',py);save(nm+'.'+ext,raw);toast('ZIP library unavailable; downloaded files separately.')}}
