
import numpy as np
from numba import jit
from .constants import *
from .Ekb_mech import ETNNP
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time


class Tissue2D:
    def __init__(self, nx=50, ny=50, dx=0.01, dy=0.01):
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.n_cells = nx * ny
        
        self.D = 0.001
        self.gap_conductance = 0.0
        self.dt = 0.005
        
        self.cell = ETNNP()
        
        self.Y = np.zeros((self.n_cells, 26))
        self._initialize_cells()
        
        self._stabilize()
        
        self._build_mechanical_mesh()
    
    def _stabilize(self, n_steps=500, dt=0.005):
        for step in range(n_steps):
            t = step * dt
            dY = np.zeros((self.n_cells, 26))
            for i in range(self.n_cells):
                dY[i] = self.cell.main(t, self.Y[i])
            self.Y += dY * dt
            
            self.Y[:, 11] = np.clip(self.Y[:, 11], -150, 100)
            self.Y[:, 10] = np.clip(self.Y[:, 10], 0, 1)
            self.Y[:, 8] = np.clip(self.Y[:, 8], 0, 1)
            self.Y[:, 9] = np.clip(self.Y[:, 9], 0, 1)
            self.Y[:, 0] = np.clip(self.Y[:, 0], 0, 1)
            self.Y[:, 1] = np.clip(self.Y[:, 1], 0, 1)
            self.Y[:, 2] = np.clip(self.Y[:, 2], 0, 1)
            self.Y[:, 3] = np.clip(self.Y[:, 3], 0, 1)
            self.Y[:, 7] = np.clip(self.Y[:, 7], 0, 1)
            self.Y[:, 17] = np.clip(self.Y[:, 17], 0, 1)
            self.Y[:, 18] = np.clip(self.Y[:, 18], 0, 1)
        
    def _initialize_cells(self):
        Y0 = np.zeros(26)
        Y0[0] = 3.373e-5
        Y0[1] = 0.9755
        Y0[2] = 0.9953
        Y0[3] = 0.7888
        Y0[4] = 3.64
        Y0[5] = 0.000126
        Y0[6] = 0.00036
        Y0[7] = 0.9073
        Y0[8] = 0.7444
        Y0[9] = 0.7045
        Y0[10] = 0.00172
        Y0[11] = -85.23
        Y0[12] = 136.89
        Y0[13] = 0.00621
        Y0[14] = 0.4712
        Y0[15] = 0.0095
        Y0[16] = 8.604
        Y0[17] = 2.42e-8
        Y0[18] = 0.999998
        Y0[19] = 0.0
        Y0[20] = 0.0
        Y0[21] = 2.726318970e-6
        Y0[22] = 6.7e-5
        Y0[23] = 0.436321675
        Y0[24] = 0.436328344
        Y0[25] = 0.088805830
        
        for i in range(self.n_cells):
            self.Y[i] = Y0.copy()
    
    def _build_mechanical_mesh(self):
        self.n_dof = self.n_cells * 2
        
        n = self.nx * self.ny
        self.K_global = np.zeros((self.n_dof, self.n_dof))
        self.F_global = np.zeros(self.n_dof)
        
        self.elements = []
        for i in range(self.nx - 1):
            for j in range(self.ny - 1):
                n1 = i * self.ny + j
                n2 = n1 + 1
                n3 = n1 + self.ny
                n4 = n3 + 1
                self.elements.append([n1, n2, n3, n4])
        
        E = 100.0
        nu = 0.3
        D = E / (1 - nu**2)
        
        for elem in self.elements:
            x1, y1 = self._get_node_coords(elem[0])
            x2, y2 = self._get_node_coords(elem[1])
            x3, y3 = self._get_node_coords(elem[2])
            x4, y4 = self._get_node_coords(elem[3])
            
            Ke = self._quadrilateral_stiffness(x1, y1, x2, y2, x3, y3, x4, y4, D, nu)
            
            for ni, i in enumerate(elem):
                for nj, j in enumerate(elem):
                    for di in range(2):
                        for dj in range(2):
                            self.K_global[2*i+di, 2*j+dj] += Ke[2*ni+di, 2*nj+dj]
    
    def _get_node_coords(self, idx):
        i = idx // self.ny
        j = idx % self.ny
        return i * self.dx, j * self.dy
    
    def _quadrilateral_stiffness(self, x1, y1, x2, y2, x3, y3, x4, y4, D, nu):
        Ke = np.zeros((8, 8))
        
        ax = (x2 + x3 - x1 - x4) / 4
        ay = (y2 + y3 - y1 - y4) / 4
        bx = (x3 + x4 - x1 - x2) / 4
        by = (y4 + y1 - y2 - y3) / 4
        
        if abs(ax) < 1e-10 and abs(ay) < 1e-10:
            ax = (x2 - x1) / 2
            ay = (y2 - y1) / 2
            bx = (x4 - x1) / 2
            by = (y4 - y1) / 2
        
        J = ax * by - ay * bx
        
        B = np.array([
            [-1, 0, 1, 0, 1, 0, -1, 0],
            [0, -1, 0, -1, 0, 1, 0, 1],
            [-(by-ay), (bx-ax), -(by+ay), -(bx+ax), (by+ay), (bx+ax), (by-ay), -(bx-ax)]
        ]) / (4 * J)
        
        C = np.array([
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1-nu)/2]
        ]) * D
        
        Ke = J * B.T @ C @ B
        
        return Ke
    
    def _get_neighbors(self, idx):
        i = idx // self.ny
        j = idx % self.ny
        neighbors = []
        
        if i > 0:
            neighbors.append((idx - self.ny, 1.0))
        if i < self.nx - 1:
            neighbors.append((idx + self.ny, 1.0))
        if j > 0:
            neighbors.append((idx - 1, 1.0))
        if j < self.ny - 1:
            neighbors.append((idx + 1, 1.0))
        
        if i > 0 and j > 0:
            neighbors.append((idx - self.ny - 1, 0.707))
        if i > 0 and j < self.ny - 1:
            neighbors.append((idx - self.ny + 1, 0.707))
        if i < self.nx - 1 and j > 0:
            neighbors.append((idx + self.ny - 1, 0.707))
        if i < self.nx - 1 and j < self.ny - 1:
            neighbors.append((idx + self.ny + 1, 0.707))
        
        return neighbors
    
    def compute_gap_junction_current(self, V, idx):
        neighbors = self._get_neighbors(idx)
        i_gap = 0.0
        for n_idx, weight in neighbors:
            i_gap += self.gap_conductance * weight * (V[n_idx] - V[idx])
        return i_gap
    
    def step_mechanics(self, dt):
        forces = np.zeros(self.n_dof)
        
        for i in range(self.n_cells):
            cell = self.cell
            l_1 = self.Y[i, 23]
            l_2 = self.Y[i, 24]
            N = self.Y[i, 21]
            A = self.Y[i, 22]
            v = self.Y[i, 19]
            
            F_CE = llambda * cell.P_star(v) / cell.G_star(v) * N
            
            forces[2*i] = F_CE
            forces[2*i+1] = F_CE
        
        self.F_global = forces
        
        fixed_nodes = []
        for i in range(self.nx):
            fixed_nodes.append(2 * (i * self.ny) + 0)
            fixed_nodes.append(2 * (i * self.ny) + 1)
        
        free_dofs = np.setdiff1d(np.arange(self.n_dof), fixed_nodes)
        
        K_free = self.K_global[np.ix_(free_dofs, free_dofs)]
        F_free = self.F_global[free_dofs]
        
        try:
            u_free = np.linalg.solve(K_free, F_free)
            u = np.zeros(self.n_dof)
            u[free_dofs] = u_free
        except:
            u = np.zeros(self.n_dof)
        
        return u
    
    def simulate(self, t_max, dt, save_interval=10):
        n_steps = int(t_max / dt)
        n_save = n_steps // save_interval
        
        V_history = np.zeros((n_save, self.n_cells))
        Ca_history = np.zeros((n_save, self.n_cells))
        Force_history = np.zeros(n_save)
        Time_history = np.zeros(n_save)
        
        save_idx = 0
        
        print(f"Starting simulation: {n_steps} steps, dt={dt} ms")
        start_time = time.time()
        
        for step in range(n_steps):
            t = step * dt
            
            dY = np.zeros((self.n_cells, 26))
            
            for i in range(self.n_cells):
                dY[i] = self.cell.main(t, self.Y[i])
            
            V = self.Y[:, 11]
            
            for i in range(self.n_cells):
                i_gap = self.compute_gap_junction_current(V, i)
                dY[i, 11] += i_gap
            
            self.Y += dY * dt
            
            self.Y[:, 11] = np.clip(self.Y[:, 11], -150, 100)
            self.Y[:, 10] = np.clip(self.Y[:, 10], 0, 1)
            self.Y[:, 8] = np.clip(self.Y[:, 8], 0, 1)
            self.Y[:, 9] = np.clip(self.Y[:, 9], 0, 1)
            self.Y[:, 0] = np.clip(self.Y[:, 0], 0, 1)
            self.Y[:, 1] = np.clip(self.Y[:, 1], 0, 1)
            self.Y[:, 2] = np.clip(self.Y[:, 2], 0, 1)
            self.Y[:, 3] = np.clip(self.Y[:, 3], 0, 1)
            self.Y[:, 7] = np.clip(self.Y[:, 7], 0, 1)
            self.Y[:, 17] = np.clip(self.Y[:, 17], 0, 1)
            self.Y[:, 18] = np.clip(self.Y[:, 18], 0, 1)
            
            if step % 100 == 0:
                u = self.step_mechanics(dt)
                for i in range(self.n_cells):
                    self.Y[i, 23] += u[2*i] * 0.001
            
            if step % save_interval == 0:
                V_history[save_idx] = self.Y[:, 11]
                Ca_history[save_idx] = self.Y[:, 5]
                
                total_force = 0.0
                for i in range(self.n_cells):
                    cell = self.cell
                    v = self.Y[i, 19]
                    N = self.Y[i, 21]
                    F_CE = llambda * cell.P_star(v) / cell.G_star(v) * N
                    total_force += F_CE
                Force_history[save_idx] = total_force / self.n_cells
                
                Time_history[save_idx] = t
                save_idx += 1
            
            if step % 1000 == 0:
                elapsed = time.time() - start_time
                progress = (step / n_steps) * 100
                print(f"Progress: {progress:.1f}% | Time: {t:.1f} ms | Elapsed: {elapsed:.1f}s")
        
        print(f"Simulation complete in {time.time() - start_time:.1f}s")
        
        return {
            'V': V_history[:save_idx],
            'Ca': Ca_history[:save_idx],
            'Force': Force_history[:save_idx],
            'Time': Time_history[:save_idx]
        }
    
    def apply_stimulus(self, region='center', radius=5, amplitude=52.0):
        pass
    
    def set_stimulus_at_time(self, time, region='center', radius=5, amplitude=52.0):
        if region == 'center':
            cx = self.nx // 2
            cy = self.ny // 2
        elif isinstance(region, tuple):
            cx, cy = region
        else:
            cx = self.nx // 2
            cy = self.ny // 2
        
        cx = int(cx)
        cy = int(cy)
        
        for i in range(self.nx):
            for j in range(self.ny):
                dist = np.sqrt((i - cx)**2 + (j - cy)**2)
                if dist <= radius:
                    idx = i * self.ny + j
                    self.Y[idx, 11] = -85.0


def visualize_results(results, nx=50, ny=50):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    n_frames = results['V'].shape[0]
    frame = n_frames // 2
    
    V_frame = results['V'][frame].reshape(nx, ny)
    im1 = axes[0, 0].imshow(V_frame, cmap='RdBu', aspect='auto')
    axes[0, 0].set_title(f'Transmembrane Potential at t={results["Time"][frame]:.0f} ms')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('y')
    plt.colorbar(im1, ax=axes[0, 0], label='V (mV)')
    
    Ca_frame = results['Ca'][frame].reshape(nx, ny)
    im2 = axes[0, 1].imshow(Ca_frame, cmap='hot', aspect='auto')
    axes[0, 1].set_title(f'Intracellular Calcium at t={results["Time"][frame]:.0f} ms')
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('y')
    plt.colorbar(im2, ax=axes[0, 1], label='Ca (mM)')
    
    axes[1, 0].plot(results['Time'], results['Force'])
    axes[1, 0].set_xlabel('Time (ms)')
    axes[1, 0].set_ylabel('Average Force')
    axes[1, 0].set_title('Tissue Force')
    axes[1, 0].grid(True)
    
    center_idx = nx // 2 * ny + ny // 2
    V_center = results['V'][:, center_idx]
    axes[1, 1].plot(results['Time'], V_center)
    axes[1, 1].set_xlabel('Time (ms)')
    axes[1, 1].set_ylabel('V (mV)')
    axes[1, 1].set_title('Action Potential at Center')
    axes[1, 1].grid(True)
    
    plt.tight_layout()
    plt.savefig('tissue2D_results.png', dpi=150)
    plt.show()


if __name__ == "__main__":
    tissue = Tissue2D(nx=50, ny=50, dx=0.02, dy=0.02)
    
    tissue.apply_stimulus(region='center', radius=5)
    
    results = tissue.simulate(t_max=500, dt=0.05, save_interval=5)
    
    visualize_results(results)
