"""Real converter parity and above-stack-limit regression; no OpenMX runs."""
import re,tempfile,subprocess,time,resource,unittest
from pathlib import Path
import numpy as np,h5py
project=Path(__file__).resolve().parents[1]
fixture=project/'tests/fixtures/kfold_coefficients'
name='kfold_g_vectors_direct'
old=project/'bin/openmx_postproc_standard_axes'
new=project/'bin/openmx_postproc_4d_plus'


class HeapConverterTests(unittest.TestCase):
    @unittest.skipUnless(old.is_file() and new.is_file() and (fixture/(name+'.dat')).is_file(),
                         "Requires local SPADExp executables and the four-point graphite fixture")
    def test_heap_conversion_matches_reference_above_stack_limit(self):
        base=(fixture/f'{name}.dat').read_text()
        mo=base.split('<MO.kpoint\n')[1].split('MO.kpoint>')[0]
        points=np.array([[float(x) for x in row.split()] for row in mo.strip().splitlines()])
        with tempfile.TemporaryDirectory(prefix='heap4d_') as d:
            root=Path(d)
            for suffix in ('.out','.v0.cube'):(root/(name+suffix)).symlink_to(fixture/(name+suffix))
            inp=root/'input.dat';inp.write_text(base)
            def run(binary):
                started=time.monotonic()
                with (root/'log').open('w') as log:
                    result=subprocess.run([binary,inp],cwd=root,stdout=log,stderr=subprocess.STDOUT,
                        preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_STACK,(8*1024*1024,resource.getrlimit(resource.RLIMIT_STACK)[1])))
                assert result.returncode==0,(root/'log').read_text()[-2000:]
                out=root/(name+'.hdf5');assert out.exists(),(root/'log').read_text()[-2000:]
                return out,time.monotonic()-started
            out,_=run(old)
            with h5py.File(out,'r') as f:
                paths=['Output/Band']
                f['Output/LCAO'].visititems(lambda k,v:paths.append('Output/LCAO/'+k) if isinstance(v,h5py.Dataset) else None)
                reference={p:f[p][:] for p in paths}
            out,_=run(new)
            with h5py.File(out,'r') as f:
                for p,a in reference.items():np.testing.assert_array_equal(f[p][:],a)
            print('4 points: private converter exactly matches original bands and all coefficients',flush=True)
            count=6912
            large=re.sub(r'(<SPADExp.range\n)[^\n]+',r'\g<1> 1 0 0 0 1 '+str(count),base)
            large+='\n<SPADExp.KFold.Map\n'
            for i in range(count):
                coords=' '.join(map(str,points[i%len(points)]))
                large+=f'{coords} {coords} 0 0 0 {i%len(points)}\n'
            large+='SPADExp.KFold.Map>\n';inp.write_text(large)
            out,elapsed=run(new)
            with h5py.File(out,'r') as f:
                assert f['Input/Kpath/Coordinates'].shape==(count,3)
                for p,a in reference.items():
                    np.testing.assert_array_equal(f[p][:],a[np.arange(count)%len(points)])
            print('6912 points: conversion completed, all repeated data exactly matched, seconds=',elapsed,
                  'stack bytes=',resource.getrlimit(resource.RLIMIT_STACK)[0],flush=True)
