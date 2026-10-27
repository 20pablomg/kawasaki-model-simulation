"""Integration smoke tests for a compiled executable, run with KAWASAKI_EXE env var."""
import csv
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest

@pytest.mark.parametrize('periodic',[0,1])
def test_conservation_and_outputs(tmp_path,periodic):
    exe=os.environ.get('KAWASAKI_EXE')
    if not exe:pytest.skip('Set KAWASAKI_EXE to run C++ integration tests')
    stem=tmp_path/'run'
    subprocess.run([exe,'--L','8','--T','2.2','--m0','0.25','--sweeps','60','--burn','20','--stride','2','--seed','77','--periodic-y',str(periodic),'--output',str(stem)],check=True)
    with (tmp_path/'run_summary.csv').open() as f:row=next(csv.DictReader(f))
    assert abs(float(row['mean_magnetization'])-float(row['m0_realized']))<1e-12
    assert int(row['n_samples'])==20
    assert float(row['specific_heat'])>=0
    assert float(row['domain_response_proxy'])>=0
    profile=np.genfromtxt(tmp_path/'run_density.csv',delimiter=',',names=True)
    assert len(profile)==8
    assert abs(profile['mean_plus_density'].mean()-(1+float(row['m0_realized']))/2)<1e-12
    grid=np.loadtxt(tmp_path/'run_final.csv',delimiter=',')
    assert grid.shape==(8,8) and set(np.unique(grid))=={-1.,1.}
    assert abs(grid.mean()-float(row['m0_realized']))<1e-12
    if periodic==0:
        assert (grid[0,:]==-1).all() and (grid[-1,:]==1).all()
