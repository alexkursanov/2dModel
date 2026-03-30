import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from cardiac_model import Tissue2D, visualize_results


def main():
    print("Initializing tissue...")
    tissue = Tissue2D(nx=20, ny=20, dx=0.02, dy=0.02)
    
    print("Applying stimulus...")
    tissue.apply_stimulus(region='center', radius=3)
    
    print("Running simulation...")
    results = tissue.simulate(t_max=200, dt=0.05, save_interval=10)
    
    print(f"Simulation complete!")
    print(f"Time points: {len(results['Time'])}")
    print(f"Max voltage: {results['V'].max():.2f} mV")
    print(f"Min voltage: {results['V'].min():.2f} mV")
    print(f"Max force: {results['Force'].max():.4f}")
    
    visualize_results(results, nx=20, ny=20)


if __name__ == "__main__":
    main()
