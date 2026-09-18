import json, time, numpy as np, sionna.rt as rt
from sionna.rt import load_scene, PlanarArray, Transmitter, Receiver, PathSolver
scene = load_scene(rt.scene.munich); scene.frequency=28e9
scene.tx_array = PlanarArray(num_rows=1, num_cols=1, vertical_spacing=0.5, horizontal_spacing=0.5, pattern="iso", polarization="V"); scene.rx_array=scene.tx_array
BS=[8.5,21.0,27.0]; tx=Transmitter(name="ue",position=[0,0,1.5]); rx=Receiver(name="bs",position=BS); scene.add(tx); scene.add(rx); solver=PathSolver()
res=[]; t0=time.time()
for x in np.arange(-120,121,8.0):
    for y in np.arange(-120,161,8.0):
        tx.position=[float(x),float(y),1.5]
        p=solver(scene,max_depth=0,los=True,specular_reflection=False,diffuse_reflection=False,refraction=False,synthetic_array=True,samples_per_src=5000,seed=1)
        a,_=p.cir(normalize_delays=False,out_type='numpy'); los=a.shape[-2]>0
        if los: res.append({'x':float(x),'y':float(y),'los':True}); continue
        p=solver(scene,max_depth=2,los=False,specular_reflection=True,diffuse_reflection=False,refraction=False,synthetic_array=True,samples_per_src=3*10**4,seed=1)
        a,_=p.cir(normalize_delays=False,out_type='numpy'); n=a.shape[-2]; pw=float(np.sum(np.abs(a)**2)) if n else 0.0
        res.append({'x':float(x),'y':float(y),'los':False,'npaths':int(n),'power':pw,'dist':float(np.hypot(x-BS[0],y-BS[1]))})
json.dump(res,open('t6_q12_rt_scan2.json','w')); print('scan2 done',round(time.time()-t0),'s; LOS',sum(r['los'] for r in res),'NLOS w/ paths',sum((not r['los']) and r['npaths']>0 for r in res))
nl=[r for r in res if (not r['los']) and r['npaths']>0]; nl.sort(key=lambda r:-r['power'])
for r in nl[:10]: print(r)
