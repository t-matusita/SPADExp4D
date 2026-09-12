"""4D regression tests; numerical jobs are simulated in isolated temporary files.

Run with QT_QPA_PLATFORM=offscreen in an environment containing PyQt6.
"""
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import h5py
import numpy as np
APP = Path(__file__).resolve().parents[1] / 'python'
sys.path.insert(0, str(APP))
import spadexp_4d_core as core
import _main as viewer


class MeshTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='spadexp4d-tests-')
        self.addCleanup(tmp.cleanup); self.root = Path(tmp.name)
        self.s = dict(viewer.DEFAULT_SETTINGS, N1=2, N2=3, N3=2)
        self.calls = []
        self.original_request = core.mesh_request
        def request(settings):
            return {'schema': 1, 'divisions': list(core.divisions(settings)), 'material': 'test graphite'}
        patches = patch.multiple(core, CACHE=self.root/'cache', RESTART=self.root/'restart',
                                 CHUNK_POINTS=5, mesh_request=request,
                                 ensure_postprocessor=lambda *args: None)
        patches.start(); self.addCleanup(patches.stop)

    def raw(self, folder, points, *args):
        self.calls.append(points.copy()); folder.mkdir(parents=True, exist_ok=True)
        p = folder/'raw.hdf5'; n = len(points)
        with h5py.File(p, 'w') as f:
            f.attrs['AxisConvention'] = 'a1,a2=in-plane; a3=surface-normal'
            for name in ('Atomic.Species', 'Atoms.SpeciesAndCoordinates', 'UnitCells'):
                f.create_group('Input/'+name)
            f['Input/Atomic.Species'].create_dataset('Labels', data=np.array([b'C']))
            f['Input/Atomic.Species'].create_dataset('Orbits', data=np.array([b's2p2d1']))
            f['Input/Atoms.SpeciesAndCoordinates'].create_dataset(
                'Labels', data=np.array([b'C', b'C', b'C', b'C']))
            f.create_dataset('Input/Kpath/Coordinates', data=points)
            f.create_dataset('Output/Band', data=np.tile(np.arange(16)*0.1, (n,1))+points[:,0,None]*0.01)
            f['Output'].attrs.update(EF_Eh=0.25, Spin='off',
                                     BandIndexStartZeroBased=0,
                                     BandIndexEndZeroBased=15)
            f.create_dataset('Output/Potential/V0', data=np.zeros((2,2,2)))
            for atom in range(1,5):
                for orb,m in (('s0',1),('s1',1),('p0',3),('p1',3),('d0',5)):
                    coeff = np.zeros((n,16,m,2),dtype='f8')
                    coeff[...,0] = np.pi + points[:,0,None,None]
                    coeff[...,1] = np.e + points[:,2,None,None]
                    f.create_dataset(f'Output/LCAO/{atom}_C/{orb}',data=coeff)
        return p, 1.0, 0.0

    def database(self):
        with patch.object(core,'calculate_chunk',side_effect=self.raw):
            return core.ensure_mesh(self.s)[0]

    def test_grid_unique_half_open_and_inside_first_bz(self):
        for ns in ((1,1,1),(2,3,4),(7,5,3)):
            points=core.grid_points(ns)
            self.assertEqual(points.shape,(np.prod(ns),3))
            residues=np.mod(points,1).round(10)
            self.assertEqual(len(np.unique(residues,axis=0)),len(points))
            d0=np.sum((points@core.G)**2,axis=1)
            for shift in core.SHIFTS:
                self.assertTrue(np.all(d0 <= np.sum(((points-shift)@core.G)**2,axis=1)+1e-10))
            np.testing.assert_allclose(core.grid_points(ns,0,min(3,len(points))),points[:3])

    def test_complete_compressed_cache_and_reuse(self):
        p=self.database()
        self.assertEqual([len(c) for c in self.calls],[12])
        with h5py.File(p,'r') as f:
            self.assertTrue(f.attrs['Complete']); self.assertTrue(f['Mesh/CompletedChunks'][:].all())
            coeff=f['Output/LCAO/1_C/p0']
            self.assertEqual(coeff.dtype,np.float32); self.assertEqual(coeff.compression,'gzip')
            expected=np.pi+core.grid_points((2,3,2))[:,0]
            np.testing.assert_allclose(coeff[:,0,0,0],expected,rtol=1e-7)
            self.assertEqual(f['Output/Band'].dtype,np.float64)
        with patch.object(core,'calculate_chunk',side_effect=AssertionError('OpenMX rerun')):
            p2,elapsed=core.ensure_mesh(dict(self.s,hv=800,eb_max=50,polarization='c+',a=1,dE=1))
        self.assertEqual(p,p2);self.assertEqual(elapsed,0)
        self.assertFalse(core.mesh_plan(dict(self.s,N3=4))['reuse'])

    def test_resume_skips_legacy_checkpoints_in_one_call(self):
        import json
        plan=core.mesh_plan(self.s);plan['folder'].mkdir(parents=True)
        points=core.grid_points((2,3,2)); raw,_,_=self.raw(self.root/'fixture',points)
        partial=plan['folder']/'eigenstates.partial.hdf5'
        with h5py.File(partial,'w') as f,h5py.File(raw,'r') as source:
            f.attrs.update(Request=json.dumps(plan['request'],sort_keys=True),Complete=False)
            mesh=f.create_group('Mesh');mesh.attrs.update(Divisions=[2,3,2],ReciprocalVectors_Ainv=core.G)
            mesh['FirstBZCoordinates']=points
            mesh['CompletedChunks']=[False,True,False]
            core.initialize_database(f,source,12)
            for key in ['Output/Band']+core.lcao_paths(source):f[key][5:10]=source[key][5:10]
        self.calls.clear()
        with patch.object(core,'calculate_chunk',side_effect=self.raw):
            result,_=core.ensure_mesh(self.s)
        self.assertEqual([len(c) for c in self.calls],[7])
        np.testing.assert_array_equal(self.calls[0],points[[0,1,2,3,4,10,11]])
        with h5py.File(result,'r') as f,h5py.File(raw,'r') as source:
            for key in ['Output/Band']+core.lcao_paths(source):
                np.testing.assert_allclose(f[key][:],source[key][:],rtol=1e-7)
        self.assertTrue(core.mesh_plan(self.s)['reuse'])

    def test_completed_legacy_openmx_output_is_imported_before_single_missing_job(self):
        import hashlib
        plan=core.mesh_plan(self.s);legacy=plan['folder']/'chunk_00000'
        legacy.mkdir(parents=True);points=core.grid_points((2,3,2),0,5)
        text=core.input_text(legacy.name,points)
        output=legacy/(legacy.name+'.out');output.write_text('completed legacy OpenMX output')
        (legacy/'openmx.log').write_text('The calculation was normally finished.')
        core.json_write(legacy/'openmx.json',{'input':hashlib.sha256(text.encode()).hexdigest(),
                                            'output':core.stamp(output)})
        jobs=[]
        def calculate(folder,points,*args):
            jobs.append(folder.name)
            raw,_,_=self.raw(folder,points)
            return raw,0 if folder==legacy else 1,0
        with patch.object(core,'calculate_chunk',side_effect=calculate):
            result,seconds=core.ensure_mesh(self.s)
        self.assertEqual(jobs[0],'chunk_00000')
        self.assertTrue(jobs[1].startswith('batch_'))
        self.assertEqual([len(c) for c in self.calls],[5,7])
        self.assertEqual(seconds,1)
        self.assertTrue(core.mesh_plan(self.s)['reuse'])

    def test_interrupted_import_reuses_same_single_job(self):
        jobs=[];checks=0
        def raw(folder,points,*args):
            jobs.append(folder);return self.raw(folder,points)
        def stop():
            nonlocal checks
            if jobs: checks+=1
            return checks>3
        with patch.object(core,'calculate_chunk',side_effect=raw),self.assertRaises(core.Cancelled):
            core.ensure_mesh(self.s,stop=stop)
        self.assertFalse(core.mesh_plan(self.s)['reuse'])
        with patch.object(core,'calculate_chunk',side_effect=raw):
            p,_=core.ensure_mesh(self.s)
        self.assertTrue(p.exists());self.assertEqual([len(c) for c in self.calls],[12,12])
        self.assertEqual(jobs[0],jobs[1])  # calculate_chunk reuses this job's successful OpenMX output.

    def test_invalid_database_not_reused(self):
        p=self.database()
        with h5py.File(p,'r+') as f: f.attrs['Complete']=False
        self.assertFalse(core.mesh_plan(self.s)['reuse'])
        p.write_bytes(b'broken')
        self.assertFalse(core.mesh_plan(self.s)['reuse'])

    def test_lookup_exact_periodic_and_nearest_in_cartesian_metric(self):
        p=self.database(); grid=core.EigenstateGrid(p);self.addCleanup(grid.close)
        points=core.grid_points((2,3,2));cart=points@core.G
        expected=np.arange(len(points))
        for shift in ((0,0,0),(2,-3,7),(-4,5,-6)):
            idx,error=grid.lookup(cart+np.array(shift)@core.G)
            np.testing.assert_array_equal(idx,expected); np.testing.assert_allclose(error,0,atol=1e-13)
        target=np.random.default_rng(4).uniform(-4,4,(15,3))@core.G
        idx,error=grid.lookup(target)
        for q,i,d in zip(target,idx,error):
            reduced=core.fold_first_bz((q@np.linalg.inv(core.G))[None,:])[0]@core.G
            distances=np.array([np.linalg.norm(cart+shift@core.G-reduced,axis=1) for shift in core.SHIFTS])
            self.assertAlmostEqual(d,distances.min(),places=10)
            self.assertAlmostEqual(d,distances[:,i].min(),places=10)

    def test_lookup_handles_highly_unequal_divisions(self):
        self.s.update(N1=35,N2=1,N3=1)
        p=self.database();grid=core.EigenstateGrid(p);self.addCleanup(grid.close)
        target=np.random.default_rng(5).uniform(-2,2,(20,3))@core.G
        _,error=grid.lookup(target)
        cart=grid.points@core.G
        for q,d in zip(target,error):
            reduced=core.fold_first_bz((q@np.linalg.inv(core.G))[None,:])[0]@core.G
            expected=min(np.linalg.norm(cart+shift@core.G-reduced,axis=1).min() for shift in core.SHIFTS)
            self.assertAlmostEqual(d,expected,places=10)

    def test_scalar_band_interpolation_is_periodic_and_does_not_touch_coefficients(self):
        self.s.update(N1=4,N2=5,N3=3)
        p=self.database()
        q=np.stack(np.unravel_index(np.arange(60),(4,5,3)),axis=-1)/np.array((4,5,3))
        with h5py.File(p,'r+') as f:
            bands=np.empty((60,16))
            bands[:]=2*q[:,0,None]+3*q[:,1,None]+4*q[:,2,None]+np.arange(16)*.01
            f['Output/Band'][:]=bands;ef=float(f['Output'].attrs['EF_Eh'])
        grid=core.EigenstateGrid(p);self.addCleanup(grid.close)
        target_frac=np.array([[.2,.3,.1],[.2+3,.3-2,.1+4]])
        bands,error=grid.interpolated_bands(target_frac@core.G)
        expected=(2*.2+3*.3+4*.1+np.arange(16)*.01-ef)*viewer.EH
        np.testing.assert_allclose(bands[0],expected,atol=1e-12)
        np.testing.assert_allclose(bands[1],expected,atol=1e-12)
        self.assertEqual(bands.shape,(2,16));self.assertEqual(error.shape,(2,))

    def test_adapter_preserves_band_coefficient_pair_and_requested_phase_coordinate(self):
        p=self.database();grid=core.EigenstateGrid(p);self.addCleanup(grid.close)
        requested=(grid.points[[4,1,4]]+np.array([2,-1,3]))@core.G
        dest=self.root/'path.hdf5';error=grid.write_path(dest,requested,{'Dimension':2})
        self.assertLess(error,1e-12)
        with h5py.File(dest,'r') as f:
            np.testing.assert_allclose(f['Input/Kpath/Coordinates'][:],requested@np.linalg.inv(core.G))
            np.testing.assert_array_equal(f['Input/Kpath/MeshIndex'][:],[4,1,4])
            for key in ['Output/Band']+core.lcao_paths(f):
                np.testing.assert_array_equal(f[key][:],grid.take(key,[4,1,4]))
            self.assertEqual(f['Output/LCAO/1_C/p0'].dtype,np.float32)

    def test_openmx_input_only_contains_folded_mo_points(self):
        points=core.grid_points((3,2,2));text=core.input_text('test',points)
        self.assertIn('MO.Nkpoint 12',text);self.assertIn('scf.restart on',text)
        rows=text.split('<MO.kpoint\n')[1].split('MO.kpoint>')[0]
        parsed=np.array([[float(v) for v in line.split()] for line in rows.splitlines()])
        np.testing.assert_allclose(parsed,points,atol=1e-11)
        self.assertNotIn('SPADExp.curved true',text)

    def test_silent_subprocess_stop(self):
        begin=time.monotonic()
        with self.assertRaises(core.Cancelled):
            core.run_command([viewer.sys.executable,'-c','import time; time.sleep(30)'],
                             self.root/'child.log',self.root,stop=lambda:time.monotonic()-begin>0.2)
        self.assertLess(time.monotonic()-begin,3)

    def test_actual_mesh_identity_ignores_all_pad_and_display_settings(self):
        original=self.original_request(self.s)
        changed=dict(self.s,hv=300,v0=10,wfnc=5,eb_min=2,eb_max=8,de=0.5,
                     direction='kz',a=0.5,b=-0.5,dk=0.1,dE=0.7,theta_max=45,
                     polarization='c+',alpha=10,beta=20,gamma=30,k_min=-1,k_max=1,
                     k_step=0.05,show_bz=False)
        self.assertEqual(original,self.original_request(changed))
        self.assertNotEqual(original,self.original_request(dict(changed,N1=3)))

    def test_failed_conversion_reuses_completed_openmx(self):
        points=core.grid_points((2,2,2));folder=self.root/'chunk_test'
        fail=True;calls=[]
        def command(args,log,cwd,stdin=None,expected=None,**kwargs):
            calls.append(args[0])
            if args[0]==core.OPENMX:
                Path(expected).write_text('test eigenstate output')
                Path(log).write_text('The calculation was normally finished.')
            else:
                if fail:raise RuntimeError('conversion failed')
                raw,_,_=self.raw(folder,points)
                raw.replace(expected)
            return 1.0
        with patch.object(core,'run_command',side_effect=command):
            with self.assertRaisesRegex(RuntimeError,'conversion failed'):
                core.calculate_chunk(folder,points)
            fail=False
            result,seconds,conversion_seconds=core.calculate_chunk(folder,points)
        self.assertEqual(calls.count(core.OPENMX),1)
        self.assertTrue(result.exists());self.assertEqual(seconds,0)
        self.assertGreaterEqual(conversion_seconds,0)

    def test_empty_spectral_window_skips_external_pad_and_caches_zero(self):
        path=self.database();grid=core.EigenstateGrid(path);self.addCleanup(grid.close)
        s=dict(self.s,hv=500,eb_min=100,eb_max=101,k_min=-1,k_max=1,k_step=1)
        with patch.multiple(viewer,INPUT=self.root/'input',OUTPUT=self.root/'output'):
            viewer.INPUT.mkdir();viewer.OUTPUT.mkdir()
            with patch.object(core,'run_command',side_effect=AssertionError('zero PAD is known')):
                result=viewer.pad_slice(grid,s,20,100,viewer.light_vectors(s)[0])
                again=viewer.pad_slice(grid,s,20,100,viewer.light_vectors(s)[0])
        np.testing.assert_array_equal(result['pad'],np.zeros((3,3)))
        np.testing.assert_array_equal(again['pad'],result['pad'])


class BatchTests(unittest.TestCase):
    def test_stop_preserves_closed_slices_in_memory_and_retry_completes(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); binary=root/'spadexp_4d';binary.write_text('private executable')
            mesh=root/'grid.h5';mesh.write_text('immutable grid')
            grid=SimpleNamespace(path=mesh)
            s=dict(viewer.DEFAULT_SETTINGS,k_min=-1,k_max=1,k_step=1)
            field=viewer.light_vectors(s)[0]; inputs=[]; results=[]
            def manifest(path,grid,s,v0,requested):
                inputs.append(list(requested))
            def command(args,log,cwd,**kwargs):
                output=kwargs['expected']
                for i,eb in enumerate(inputs[-1]):
                    with h5py.File(f'{output}.{i:06d}.h5','w') as f:
                        f.attrs.update(Complete=1,RawSpectrum=1,EB=eb,MaxSamplingDistance_Ainv=.1)
                        f['PAD']=np.full((3,3),eb+1);f['BandWeight']=np.full((3,3),eb+2)
                    if len(inputs)==1: raise core.Cancelled()
                with h5py.File(output,'w') as f:f.attrs.update(Complete=1,KernelSeconds=.01)
                return .1
            with patch.multiple(viewer,INPUT=root,OUTPUT=root), \
                 patch.object(viewer,'ensure_cpp_binary',return_value=binary), \
                 patch.object(viewer,'write_batch_input',side_effect=manifest), \
                 patch.object(core,'run_command',side_effect=command):
                with self.assertRaises(core.Cancelled):
                    viewer.pad_batch(grid,s,20,[0.,1.,2.],field,results.append)
                self.assertEqual([r['eb'] for r in results],[0.])
                results.clear()
                stats=viewer.pad_batch(grid,s,20,[0.,1.,2.],field,results.append)
                self.assertEqual(inputs,[[0.,1.,2.],[0.,1.,2.]])
                self.assertEqual([r['eb'] for r in results],[0.,1.,2.])
                self.assertEqual(stats['pad_computed_slices'],3)
                results.clear()
                stats=viewer.pad_batch(grid,dict(s,dk=.8,dE=1.5),20,[0.,1.,2.],field,results.append)
                self.assertEqual(len(inputs),3)
                self.assertEqual(stats['pad_computed_slices'],3)
                self.assertEqual(len(results),3)
                # The private executable has a separate PAD identity, so legacy
                # Python results cannot be mistaken for a completed C++ slice.
                self.assertNotEqual(viewer.pad_request(grid,s,20,0.,field),
                                    viewer.pad_request(grid,s,20,0.,field,binary))

    def test_lightweight_manifest_links_original_coefficient_grid(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);mesh=root/'grid.h5';source=root/'request.h5'
            with h5py.File(mesh,'w') as f:
                for name in ('Atomic.Species','Atoms.SpeciesAndCoordinates','UnitCells'):
                    f.create_group('Input/'+name)
                f['Output/LCAO/1_C/s0']=np.ones((12,16,1,2),dtype='f4')
                f.create_group('Mesh')
            s=dict(viewer.DEFAULT_SETTINGS,k_min=-1,k_max=1,k_step=1)
            viewer.write_batch_input(source,SimpleNamespace(path=mesh),s,20,[0.,1.])
            with h5py.File(source,'r') as f:
                self.assertIsInstance(f.get('Output',getlink=True),h5py.ExternalLink)
                self.assertEqual(f['Output/LCAO/1_C/s0'].shape[0],12)
                self.assertEqual(f['Input/Kpath/Coordinates'].shape,(9,3))
                np.testing.assert_array_equal(f['Batch/BindingEnergies'][:],[0.,1.])


def impulse_result():
    s=viewer.validate_settings(dict(viewer.DEFAULT_SETTINGS,hv=150,k_min=-3,k_max=3,k_step=1,
                                    eb_min=0,eb_max=6,de=1,dk=0,dE=0))
    raw={}
    for eb in viewer.energies(s):
        pad=np.zeros((7,7))
        if eb==3:pad[3,3]=1
        raw[eb]=dict(eb=eb,pad=pad,band_weight=pad.copy(),error=0,extent=(-3,3,-3,3))
    return dict(s=s,v0=20,raw_maps=raw,openmx_seconds=0,sampling_error=0)


class GaussianTests(unittest.TestCase):
    def test_cube_gaussian_axes_mass_and_repeat_changes_use_raw_data(self):
        result=impulse_result();s=result['s']
        zero=viewer.prepare_display(result,s)
        smooth=viewer.prepare_display(result,dict(s,dk=.5,dE=.5))
        cube=np.stack([smooth['maps'][eb]['pad'] for eb in viewer.energies(s)])
        self.assertAlmostEqual(cube.sum(),1,places=12)
        for index in ((4,3,3),(3,4,3),(3,3,4)):
            self.assertAlmostEqual(cube[index]/cube[3,3,3],np.exp(-2),places=12)
        np.testing.assert_array_equal(smooth['cuts'][3]['pad'],cube[3,3])
        restored=viewer.prepare_display(smooth,s)
        for eb in zero['maps']:
            np.testing.assert_array_equal(restored['maps'][eb]['pad'],zero['maps'][eb]['pad'])
            np.testing.assert_array_equal(result['raw_maps'][eb]['pad'],zero['maps'][eb]['pad'])

    def test_raw_identity_ignores_gaussian_but_tracks_energy_sampling(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'grid';path.write_text('grid')
            grid=SimpleNamespace(path=path);s=viewer.DEFAULT_SETTINGS;field=viewer.light_vectors(s)[0]
            request=viewer.pad_request(grid,s,20,0.,field)
            self.assertNotIn('dE',request);self.assertNotIn('dk',request)
            self.assertEqual(request,viewer.pad_request(grid,dict(s,dk=0,dE=2),20,0.,field))
            self.assertNotEqual(request,viewer.pad_request(grid,dict(s,de=.2),20,0.,field))


class FormattingTests(unittest.TestCase):
    def test_file_size_uses_iec_unit_and_exact_bytes(self):
        self.assertEqual(viewer.format_file_size(0),'0 B (0 bytes)')
        self.assertEqual(viewer.format_file_size(1536),'1.50 KiB (1,536 bytes)')
        self.assertEqual(viewer.format_file_size(3*1024**2),'3.00 MiB (3,145,728 bytes)')

class ViewerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app=viewer.QtWidgets.QApplication.instance() or viewer.QtWidgets.QApplication([])

    def setUp(self):
        # Do not let repository structures override the settings supplied by
        # individual GUI tests when Window.refresh_structures() runs.
        structures=patch.object(viewer,'structure_entries',return_value=[])
        structures.start();self.addCleanup(structures.stop)

    def test_gaussian_edits_update_from_memory_without_worker_or_commands(self):
        result=impulse_result();w=viewer.Window(result['s']);self.addCleanup(w.close)
        w.result=viewer.prepare_display(result,result['s'])
        with tempfile.TemporaryDirectory() as d, \
             patch.object(viewer,'INITIAL_SETTINGS_FILE',Path(d)/'settings.json'), \
             patch.object(core,'run_command',side_effect=AssertionError('external command')), \
             patch.object(core,'mesh_plan',side_effect=AssertionError('mesh planning')), \
             patch.object(viewer,'Worker',side_effect=AssertionError('PAD worker')), \
             patch.object(viewer.QtWidgets.QMessageBox,'warning',side_effect=AssertionError('unexpected GUI error')):
            w.dk.setValue(.5);w.dE.setValue(.5);w.calculate()
            self.assertIsNone(w.thread)
            saved=core.json_read(Path(d)/'settings.json')
            self.assertEqual((saved['dk'],saved['dE']),(.5,.5))
            self.assertLess(w.result['maps'][3]['pad'][3,3],1)
            self.assertIn('PAD計算なし',w.message.text())
            self.assertTrue(w.canvas.figure.axes[2].get_title().startswith('Band weight'))
            self.assertTrue(w.canvas.figure.axes[3].get_title().startswith('PAD'))
            w.dk.setValue(0);w.dE.setValue(0);w.refresh_gaussian()
            self.assertEqual(w.result['maps'][3]['pad'][3,3],1)
            self.assertEqual(w.result['raw_maps'][3]['pad'][3,3],1)

    def test_gaussian_edits_wait_for_calc_instead_of_starting_processing(self):
        w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
        with patch.object(w,'refresh_gaussian') as refresh:
            w.dk.setValue(.5); w.dE.setValue(.5); self.app.processEvents()
        refresh.assert_not_called()

    def test_current_eb_commits_only_after_enter(self):
        w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
        original=w.slider.value()
        w.eb_input.lineEdit().setText('2.0')
        self.app.processEvents()
        self.assertEqual(w.slider.value(),original)
        w.eb_input.lineEdit().returnPressed.emit()
        self.assertEqual(w.slider.value(),20)

    def test_eb_sweep_starts_at_begin_and_stops_after_last_frame(self):
        w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
        with patch.object(w,'redraw') as redraw:
            w.result={}
            w.slider.setRange(0,2);w.slider.setValue(1)
            w._slider_min=1.0;w._slider_step=.5
            redraw.reset_mock()
            w.start_eb_sweep()
            self.assertEqual(w.slider.value(),0)
            self.assertTrue(w.eb_sweep_timer.isActive())
            self.assertFalse(w.play_eb_button.isEnabled())
            self.assertTrue(w.stop_eb_button.isEnabled())
            w.advance_eb_sweep();self.assertEqual(w.slider.value(),1)
            w.advance_eb_sweep();self.assertEqual(w.slider.value(),2)
            self.assertFalse(w.eb_sweep_timer.isActive())
            self.assertTrue(w.play_eb_button.isEnabled())
            self.assertFalse(w.stop_eb_button.isEnabled())
            self.assertEqual(w.eb_input.value(),2.0)
            self.assertEqual(redraw.call_count,3)

    def test_manual_eb_entry_stops_sweep(self):
        w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
        w.result={}
        with patch.object(w,'redraw'):
            w.start_eb_sweep()
            self.assertTrue(w.eb_sweep_timer.isActive())
            w.eb_input.lineEdit().setText('2.0')
            w.eb_input.lineEdit().returnPressed.emit()
        self.assertFalse(w.eb_sweep_timer.isActive())
        self.assertTrue(w.play_eb_button.isEnabled())

    def test_ef_shift_keeps_zero_eb_display_editable(self):
        w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
        w.ef_shift.setValue(2.0)
        self.assertLessEqual(w.eb_input.minimum(),0.0)
        w.eb_input.lineEdit().setText('0.0')
        w.eb_input.lineEdit().returnPressed.emit()
        self.assertEqual(w.eb_input.value(),0.0)

    def test_startup_result_divisions_are_loaded_into_mesh_boxes(self):
        with tempfile.TemporaryDirectory() as directory:
            structure=Path(directory)/'structure'; record=structure/'results'/'saved'; record.mkdir(parents=True)
            settings=dict(viewer.DEFAULT_SETTINGS,N1=24,N2=24,N3=12)
            viewer.write_settings(record/'config.json',settings)
            with h5py.File(record/'eigenstates.hdf5','w') as file:
                file.create_group('Mesh').attrs['Divisions']=np.array([7,8,9],dtype='i8')
            (record/'manifest.json').write_text('{}')
            w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
            w.current_structure=structure
            w.refresh_results(record)
            self.assertEqual(tuple(box.value() for box in w.nboxes),(7,8,9))

    def test_photoelectron_trajectory_draws_bz_boundaries(self):
        result=impulse_result();w=viewer.Window(result['s']);self.addCleanup(w.close)
        w.result=viewer.prepare_display(result,result['s']);w.redraw()
        panel=next(ax for ax in w.canvas.figure.axes
                   if ax.get_title()=='kx-kz photoelectron trajectories')
        kx,kz=viewer.trajectory_bz_boundaries(panel.get_xlim(),panel.get_ylim()[1])
        dashed=[line for line in panel.lines if line.get_linestyle()=='--']
        self.assertTrue(kx);self.assertTrue(kz)
        self.assertEqual(len(dashed),len(kx)+len(kz))
        w.show_bz.setChecked(False)
        panel=next(ax for ax in w.canvas.figure.axes
                   if ax.get_title()=='kx-kz photoelectron trajectories')
        self.assertFalse([line for line in panel.lines if line.get_linestyle()=='--'])

    def test_mesh_widgets_steps_count_and_separate_paths(self):
        w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
        for box,value in zip(w.nboxes,(10,12,7)):box.setValue(value)
        self.assertIn('840',w.total_points.text())
        for i,label in enumerate(w.glabels):
            self.assertIn(f'{np.linalg.norm(core.G[i])/w.nboxes[i].value():.2f}',label.text())
        self.assertEqual(core.divisions(w.settings()),(10,12,7))
        self.assertEqual(viewer.INITIAL_SETTINGS_FILE.name,'spadexp_4d.ini.json')
        self.assertEqual(viewer.STRUCTURES,viewer.PROJECT/'structures')
        self.assertEqual(viewer.CACHE,viewer.DEFAULT_STRUCTURE/'cache')
        self.assertEqual(w.settings()['pad_engine'],'cpp')
        self.assertEqual(w.mesh_result.text(),'未取得')
        w.update_mesh_result({'mesh_bytes':1536,'mesh_seconds':2,'openmx_seconds':1,'mesh_path':'/tmp/test.h5'})
        self.assertIn('1.50 KiB (1,536 bytes)',w.mesh_result.text())
        self.assertIn('tmp/test.h5',w.mesh_result.text())
        self.assertIn('OpenMX 1.0秒',w.mesh_result.text())
        self.assertEqual(w.mesh_result.toolTip(),'/tmp/test.h5')
        w.pad_engine.setCurrentIndex(1)
        self.assertEqual(w.settings()['pad_engine'],'python')
        self.assertTrue(viewer.validate_settings(dict(viewer.DEFAULT_SETTINGS,use_first_bz_only=False))['use_first_bz_only'])

    def test_settings_integer_validation(self):
        for value in (0,-1,1.5,257,float('inf'),True,None):
            with self.subTest(value=value),self.assertRaises(ValueError):
                viewer.validate_settings(dict(viewer.DEFAULT_SETTINGS,N1=value))

    def test_cut_maps_axis_orientation_and_fixed_kz_does_not_change_pad(self):
        s=dict(viewer.DEFAULT_SETTINGS,k_min=-1,k_max=1,k_step=1,dk=0.0001,hv=150,a=0.5)
        axis=viewer.cartesian_axis(s);xx,yy=np.meshgrid(axis,axis)
        data=dict(eb=1,pad=10+xx+2*yy,band_weight=20+3*xx+yy,extent=(-1,1,-1,1))
        for direction in ('kx','ky'):
            cut=viewer.cut_slice(data,dict(s,direction=direction),20)
            x,y=(axis,0.5) if direction=='kx' else (0.5,axis)
            np.testing.assert_allclose(cut['pad'],10+x+2*y)
            np.testing.assert_array_equal(cut['pad'],viewer.cut_slice(data,dict(s,direction=direction,b=3),20)['pad'])
        cut=viewer.cut_slice(data,dict(s,direction='kz',a=0,b=0),20)
        self.assertEqual(cut['k'].shape,(1,))
        np.testing.assert_allclose(cut['k'],np.sqrt(169/viewer.H2M))

    def test_all_mesh_points_prepared_before_pad_and_success_emitted_once(self):
        s=dict(viewer.DEFAULT_SETTINGS,pad_engine="python",eb_max=0.2,k_min=-1,k_max=1,k_step=1)
        worker=viewer.Worker(s);done=[];errors=[];events=[]
        worker.done.connect(done.append);worker.error.connect(errors.append)
        handle,fake_name=tempfile.mkstemp(suffix='.hdf5');viewer.os.close(handle)
        fake=Path(fake_name);self.addCleanup(lambda:fake.unlink(missing_ok=True))
        def mesh(*args,**kwargs):events.append('mesh complete');return fake,0
        def pad(grid,s,v0,eb,*args):
            events.append(eb)
            return dict(eb=eb,pad=np.ones((3,3)),band_weight=np.ones((3,3)),error=0.01,extent=(-1,1,-1,1))
        from types import SimpleNamespace
        fake_grid=SimpleNamespace(path=fake,close=lambda:None)
        with patch.object(core,'ensure_mesh',side_effect=mesh), patch.object(core,'EigenstateGrid',return_value=fake_grid), \
             patch.object(viewer,'fixed_band',return_value={'error':0.02}), \
             patch.object(viewer,'pad_slice',side_effect=pad),patch.object(viewer,'mkdirs'):
            worker.run()
        self.assertEqual(errors,[]);self.assertEqual(len(done),1)
        self.assertEqual(events,['mesh complete',0,0.1,0.2])
        self.assertEqual(done[0]['sampling_error'],0.02)
        self.assertEqual(len(done[0]['cuts']),3)

    def test_closing_active_window_stops_worker_and_saves_local_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'spadexp_4d.ini.json'
            structure=Path(directory)/'structures'/'test';structure.mkdir(parents=True)
            w=viewer.Window(viewer.DEFAULT_SETTINGS);self.addCleanup(w.close)
            @viewer.QtCore.pyqtSlot()
            def run(worker):
                worker.stop_event.wait(3)
                worker.stopped.emit(dict(s=worker.s,maps={},cuts={}))
            with patch.object(viewer.Worker,'run',run),patch.object(viewer,'INITIAL_SETTINGS_FILE',path), \
                 patch.object(core,'mesh_plan',return_value={'reuse':True}):
                w.calculate();self.assertTrue(w.cancel_action.isEnabled())
                w.current_structure=structure;w.close()
                deadline=time.monotonic()+5
                while w.thread is not None and time.monotonic()<deadline:
                    self.app.processEvents();time.sleep(0.01)
                self.assertIsNone(w.thread)
            self.assertEqual(core.json_read(structure/'parameter.json')['N1'],viewer.DEFAULT_SETTINGS['N1'])

    def test_ef_shift_moves_spectra_inside_fixed_energy_window(self):
        result=impulse_result();s=result['s']
        result['line']={'k':np.array([-1.,1.]),'bands':np.array([[-1.],[-1.]])}
        w=viewer.Window(s);self.addCleanup(w.close)
        w.result=viewer.prepare_display(result,s);w.redraw()
        def panel(prefix):
            return next(ax for ax in w.canvas.figure.axes if ax.get_title().startswith(prefix))
        band0=panel('Fixed-coordinate band cut');cut0=panel('Band weight cut')
        band_y0=np.asarray(band0.lines[0].get_ydata()).copy()
        cut_y0=np.asarray(cut0.collections[0].get_coordinates())[...,1].copy()
        limits=band0.get_ylim()
        w.ef_shift.setValue(2.0)
        band1=panel('Fixed-coordinate band cut');cut1=panel('Band weight cut')
        np.testing.assert_allclose(band1.lines[0].get_ydata(),band_y0+2)
        np.testing.assert_allclose(cut1.collections[0].get_coordinates()[...,1],cut_y0+2)
        np.testing.assert_allclose(band1.get_ylim(),limits)
        np.testing.assert_allclose(cut1.get_ylim(),limits)
        self.assertIn('EF shift',band1.get_ylabel())


if __name__=='__main__':unittest.main()
