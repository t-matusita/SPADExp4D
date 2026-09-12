"""Native raw-bin regression using copied eigenstates; never runs OpenMX."""
from pathlib import Path
import sys
import shutil
import tempfile
import unittest
from unittest.mock import patch
import h5py
import numpy as np
APP = Path(__file__).resolve().parents[1] / 'python'
sys.path.insert(0, str(APP))
import _main as viewer

FIXTURE = Path(__file__).resolve().parent / 'fixtures/eigenstates.hdf5'


class RawSpectrumTests(unittest.TestCase):
    @unittest.skipUnless(FIXTURE.is_file() and viewer.SPADEXP_4D.is_file(),
                         'Requires local first-BZ fixture and private PAD executable')
    def test_half_open_bins_and_density_have_no_gaussian_tail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root / 'states.h5'; shutil.copyfile(FIXTURE, path)
            # Artificial energies isolate bin assignment; keep the real AO/geometry
            # and coefficients solely to exercise the complete native matrix path.
            with h5py.File(path, 'r+') as f:
                f['Output'].attrs['EF_Eh'] = 0.0
                bands = f['Output/Band'][:]; bands[:] = 10 / viewer.PAD_EH
                bands[:, 0] = -1.5 / viewer.PAD_EH
                bands[:, 1] = -0.5 / viewer.PAD_EH
                f['Output/Band'][:] = bands
            grid = viewer.core.EigenstateGrid(path)
            try:
                s = dict(viewer.DEFAULT_SETTINGS, hv=150, k_min=-1, k_max=1, k_step=1,
                         eb_min=0, eb_max=2, de=1, dE=5, dk=2)
                maps = []
                with patch.multiple(viewer, INPUT=root, OUTPUT=root), \
                     patch.object(viewer, 'ensure_cpp_binary', return_value=viewer.SPADEXP_4D):
                    viewer.pad_batch(grid, s, 20, [0., 1., 2.], viewer.light_vectors(s)[0], maps.append)
                for data, expected in zip(maps, (1, 1, 0)):
                    np.testing.assert_array_equal(data['band_weight'], np.full((3, 3), expected))
                self.assertGreater(sum(data['pad'].sum() for data in maps[:2]), 0)
                np.testing.assert_array_equal(maps[2]['pad'], np.zeros((3, 3)))
            finally:
                grid.close()
