import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from cardiac_model import Tissue2D, ETNNP


class TestETNNP:
    def test_initialization(self):
        cell = ETNNP()
        assert cell is not None
        assert cell.F_afterload == 0.0

    def test_main_function(self):
        cell = ETNNP()
        Y = np.zeros(26)
        Y[0] = 3.373e-5
        Y[1] = 0.9755
        Y[11] = -85.23
        
        result = cell.main(0.0, Y)
        assert result.shape == (26,)
        assert not np.isnan(result).any()


class TestTissue2D:
    def test_initialization_small(self):
        tissue = Tissue2D(nx=5, ny=5)
        assert tissue.nx == 5
        assert tissue.ny == 5
        assert tissue.n_cells == 25

    def test_neighbors(self):
        tissue = Tissue2D(nx=5, ny=5)
        
        center_idx = 2 * tissue.ny + 2
        neighbors = tissue._get_neighbors(center_idx)
        assert len(neighbors) == 8
        
        corner_idx = 0
        neighbors = tissue._get_neighbors(corner_idx)
        assert len(neighbors) == 3

    def test_node_coords(self):
        tissue = Tissue2D(nx=5, ny=5, dx=0.01, dy=0.01)
        
        x, y = tissue._get_node_coords(0)
        assert x == 0.0
        assert y == 0.0
        
        x, y = tissue._get_node_coords(6)
        assert x == 0.01
        assert y == 0.01

    def test_gap_junction_current(self):
        tissue = Tissue2D(nx=5, ny=5)
        tissue.gap_conductance = 1.0
        
        V = np.random.randn(tissue.n_cells)
        i_gap = tissue.compute_gap_junction_current(V, 12)
        assert isinstance(i_gap, float)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
