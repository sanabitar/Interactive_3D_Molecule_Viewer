#!/usr/bin/env python3
"""Runs LAURA's *generated* Blender script against a minimal MOCK of bpy/mathutils and checks the geometry it would create.
This is NOT a Blender run: it validates the script's own logic (object counts, positions, bond lengths, parallel-bond offsets). Real Blender validation is still required."""
import sys,types,math,json,re,zipfile,tempfile
from pathlib import Path
class Vector:
    def __init__(s,v=(0,0,0)):s.x,s.y,s.z=[float(a) for a in v]
    def _t(s):return (s.x,s.y,s.z)
    def __add__(s,o):return Vector((s.x+o.x,s.y+o.y,s.z+o.z))
    def __sub__(s,o):return Vector((s.x-o.x,s.y-o.y,s.z-o.z))
    def __mul__(s,k):return Vector((s.x*k,s.y*k,s.z*k))
    def __truediv__(s,k):return Vector((s.x/k,s.y/k,s.z/k))
    def cross(s,o):return Vector((s.y*o.z-s.z*o.y,s.z*o.x-s.x*o.z,s.x*o.y-s.y*o.x))
    def dot(s,o):return s.x*o.x+s.y*o.y+s.z*o.z
    @property
    def length(s):return math.sqrt(s.dot(s))
    def normalized(s):return s*(1/s.length)
    def rotation_difference(s,o):return ('quat',s._t(),o._t())
class Objs(list):
    def link(s,o):s.append(o)
class Coll:
    def __init__(s,n):s.name=n;s.objects=Objs();s.children=Objs()
class Obj:
    def __init__(s,name,data):s.name=name;s.data=data;s.location=(0,0,0);s.scale=(1,1,1);s.rotation_mode=None;s.rotation_quaternion=None;s.constraints=types.SimpleNamespace(new=lambda t:types.SimpleNamespace(type=t));s.rotation_euler=(0,0,0)
class Mesh:
    def __init__(s,n='m'):s.name=n;s.materials=[]
    def copy(s):return Mesh(s.name)
class Node:
    def __init__(s):s.inputs={'Base Color':types.SimpleNamespace(default_value=None)}
def build_mock():
    ctx=types.SimpleNamespace(object=None,scene=types.SimpleNamespace(collection=Coll('Scene'),camera=None))
    def add(**k):ctx.object=Obj('prim',Mesh())
    bpy=types.SimpleNamespace(context=ctx,data=types.SimpleNamespace(
        collections=types.SimpleNamespace(new=Coll),
        materials=types.SimpleNamespace(new=lambda n:types.SimpleNamespace(name=n,diffuse_color=None,use_nodes=False,node_tree=types.SimpleNamespace(nodes={'Principled BSDF':Node()}))),
        objects=types.SimpleNamespace(new=Obj,remove=lambda o:None),cameras=types.SimpleNamespace(new=lambda n:Mesh(n)),lights=types.SimpleNamespace(new=lambda n,t:Mesh(n))),
        ops=types.SimpleNamespace(mesh=types.SimpleNamespace(primitive_uv_sphere_add=add,primitive_cylinder_add=add)))
    sys.modules['bpy']=bpy;sys.modules['mathutils']=types.SimpleNamespace(Vector=Vector);return bpy
def check(label,py):
    bpy=build_mock();ns={};exec(compile(py,label,'exec'),ns);col=bpy.context.scene.collection.children[0];objs=col.objects;D=ns['DATA'];A=D['a'];B=D['b']
    names=[o.name for o in objs];atoms=[o for o in objs if not o.name.startswith(('Bond_','LAURA_'))];bonds=[o for o in objs if o.name.startswith('Bond_')]
    per={1:1,1.5:2,2:2,3:3};exp_bonds=sum(per.get(o,1) for _,_,o in B)
    assert len(atoms)==len(A),f'{label}: {len(atoms)} atom objects vs {len(A)}'
    assert len(bonds)==exp_bonds,f'{label}: {len(bonds)} bond cylinders vs {exp_bonds} expected from orders'
    assert {'LAURA_Camera','LAURA_Sun','LAURA_Target'}<=set(names),'camera/light/target missing'
    assert len(set(names))==len(names),'duplicate object names'
    cen=[sum(a[k] for a in A)/len(A) for k in(1,2,3)];assert max(abs(c) for c in cen)<1e-3,f'{label}: not centred {cen}'
    for o,a in zip(atoms,A):assert all(abs(o.location[i]-a[i+1])<1e-9 for i in range(3)),'atom position mismatch'
    # bond cylinders: length == interatomic distance; parallel members offset perpendicular to the bond
    k=0
    for i,j,od in B:
        p,q=Vector(A[i][1:]),Vector(A[j][1:]);d=q-p;n=per.get(od,1);grp=bonds[k:k+n];k+=n
        for c in grp:assert abs(c.scale[2]-d.length)<1e-9,'cylinder length != bond length'
        if n>1:
            mid=(p+q)/2
            for c in grp:off=(c.location if isinstance(c.location,Vector) else Vector(c.location))-mid;assert abs(off.dot(d.normalized()))<1e-9,'offset not perpendicular to bond'
            assert abs((grp[0].location-grp[1].location).length-0.26*(1 if od in(1.5,2) else 1))<1e-9 or n==3,'unexpected spacing'
    return f'{label}: {len(A)} atoms, {len(B)} bonds -> {len(atoms)} atom objects + {len(bonds)} bond cylinders; all geometry checks OK'
if __name__=='__main__':
    sys.argv=['x'];__file__=str(Path(__file__).parent/'browser'/'run_browser_tests.py');exec(open(__file__).read().split("@case('01")[0])
    import tempfile
    out=[]
    with sync_playwright() as pw:
        for mid in ['alanine_h','ethene','hcn','co2','benzene','pep_alpha','pyridine','imidazole','urea','pep_hlh','pep_pair']:
            pg=new_page(pw);sel_mol(pg,mid);sec(pg,'Teaching, session');z=Path(tempfile.mkdtemp())/'b.zip'
            with pg.expect_download() as d:pg.click('#bl')
            d.value.save_as(str(z));py=zipfile.ZipFile(z).read('import_molecule.py').decode();out.append(check(mid,py));pg._b.close()
        pg=new_page(pw);imp(pg,'alpha12.pdb');sec(pg,'Teaching, session');z=Path(tempfile.mkdtemp())/'b.zip'
        with pg.expect_download() as d:pg.click('#bl')
        d.value.save_as(str(z));out.append(check('imported alpha12.pdb',zipfile.ZipFile(z).read('import_molecule.py').decode()));pg._b.close()
    print('\n'.join(out));print('MOCK-BASED check passed for',len(out),'molecules (not a Blender run)')
