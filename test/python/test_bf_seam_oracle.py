"""L4 square seam: gauge, declared subgroup and arbitrary-orbital projector."""
import copy
import json
from pathlib import Path
import unittest

import numpy as np

from bf_canonical_model import CanonicalModel
from bf_seam_oracle import directed_coefficients, gauge_error, direct_configuration_error, scaled_error
from bf_seam_phase_tool import (RectangularSeam, definition_rows, range_table,
                                validate_transform, validate_orbital_subgroup)


def square_model():
    root = Path(__file__).resolve().parent / "data/BackFlow_Seam_L4"
    _, rows, table = range_table(root / "rangebf.def")
    # Populate only the independent model inputs; no C output is invented.
    model = CanonicalModel.__new__(CanonicalModel)
    model.nsite, model.nsite2, model.nsize, model.ne = 16, 32, 16, 8
    model.nrange, model.nrangeidx, model.nprojbf = 5, 4, 10
    model.posbf = [[k for i, k, _, _ in rows if i == site] for site in range(16)]
    model.rangeidx = [[table.get((i,k), (-1,0))[0] for k in range(16)] for i in range(16)]
    model.seam_phase = [[table.get((i,k), (-1,0))[1] for k in range(16)] for i in range(16)]
    model.bfsubidx = [[0,1,2,3], [1,4,5,6], [2,5,7,8], [3,6,8,9]]
    _, raw = definition_rows(root / "orbitalidx.def", (2,3,4))
    model.orbital_idx = [[0]*16 for _ in range(16)]
    model.orbital_sign = [[1]*16 for _ in range(16)]
    for row in raw:
        if len(row) == 2:
            continue
        i,j,index = map(int,row[:3])
        model.orbital_idx[i][j] = index
        model.orbital_sign[i][j] = int(row[3]) if len(row) == 4 else 1
    model.nslater = max(map(max,model.orbital_idx))+1
    model.slater = [complex(np.sin(0.39*(i+1)), 0.23*np.cos(0.27*i)) for i in range(model.nslater)]
    model.projbf = [1.07]+[complex(0.03*i, -0.01*(i%3)) for i in range(1,10)]
    model.ele_idx = [0,1,3,5,6,10,12,15, 0,2,3,4,8,10,13,15]
    model.count = model._make_count(model.ele_idx)
    model.nmp = model.nqpfull = 4
    model.nsp = 1
    model.spin = [(0j,1+0j,0j)]
    model.weight = [0.25+0j]*4
    grid = RectangularSeam(4,4,"x","x-fast")
    transforms = [grid.translation(dx,dy) for dx,dy in ((0,0),(2,0),(0,2),(2,2))]
    model.transform,model.sign = map(list,zip(*transforms))
    for mapping,signs in transforms:
        validate_transform(table,mapping,signs)
        validate_orbital_subgroup(root / "orbitalidx.def",mapping,signs)
    return model


class SeamPhysics(unittest.TestCase):
    def test_square_contracts(self):
        model = square_model()
        errors = {"gauge":gauge_error(model), "projector":direct_configuration_error(model)}
        # The reference orbital is covariant only under its declared subgroup.
        base = np.asarray(model.build())[0]
        errors["subgroup"] = 0
        for mapping,signs in zip(model.transform,model.sign):
            translated = copy.deepcopy(model)
            translated.ele_idx = [mapping[i] for i in model.ele_idx]
            translated.count = translated._make_count(translated.ele_idx)
            index = np.asarray(mapping+[i+16 for i in mapping])
            g = np.asarray(signs+signs)
            actual = np.asarray(translated.build())[0][index][:,index]
            errors["subgroup"] = max(errors["subgroup"],scaled_error(actual,base*g[:,None]*g[None,:]))
        # General f is also covered: destroy the subgroup tying deliberately.
        arbitrary = copy.deepcopy(model)
        arbitrary.nslater = 256
        arbitrary.orbital_idx = np.arange(256).reshape(16,16).tolist()
        arbitrary.orbital_sign = [[1]*16 for _ in range(16)]
        arbitrary.slater = [complex(np.sin(i*i+0.3),np.cos(i+0.7)) for i in range(256)]
        errors["general_projector"] = direct_configuration_error(arbitrary)
        self.assertLessEqual(max(errors.values()),1e-12,errors)
        # Removing seam phases must be detected on this nonzero-Theta fixture.
        mutant = copy.deepcopy(model)
        mutant.seam_phase = [[1]*16 for _ in range(16)]
        self.assertGreater(scaled_error(model.build(),mutant.build()),1e-4)
        # Both oriented ProjBF contractions are nonzero, including seam terms.
        ia,ib,pa,qa,pb,qb = directed_coefficients(model,0,3,0)
        self.assertGreater(max(abs(pa)),1e-4)
        self.assertGreater(max(abs(pb)),1e-4)
        print(json.dumps(errors,sort_keys=True))

    def test_tool_rejects_ambiguity(self):
        for numbering in ("x-fast","y-fast"):
            grid = RectangularSeam(4,4,"xy",numbering)
            with self.assertRaisesRegex(ValueError,"antipodal"):
                grid.phase(grid.site(0,0),grid.site(2,0))
            with self.assertRaisesRegex(ValueError,"outside"):
                grid.phase(-1,0)
            self.assertEqual(grid.phase(grid.site(0,0),grid.site(3,0)),-1)
            self.assertEqual(grid.phase(grid.site(3,0),grid.site(0,0)),-1)
            self.assertEqual(grid.phase(grid.site(0,0),grid.site(3,3)),1)
        root = Path(__file__).resolve().parent / "data/BackFlow_Seam_L4"
        grid = RectangularSeam(4,4,"x","x-fast")
        with self.assertRaisesRegex(ValueError,"not covariant"):
            validate_orbital_subgroup(root / "orbitalidx.def",*grid.translation(1,0))


if __name__ == "__main__":
    unittest.main()
