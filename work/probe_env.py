#!/usr/bin/env python3
"""H800节点环境探测: nvcc/github/torch.cuda/openfold可装性."""
import subprocess, sys, os
def run(cmd):
    try:
        r=subprocess.run(cmd,shell=True,capture_output=True,text=True,timeout=60)
        return r.returncode, r.stdout[-800:], r.stderr[-800:]
    except Exception as e:
        return -1, "", str(e)
print("=== nvidia-smi ===",flush=True)
print(run("nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv")[1])
print("=== nvcc ===",flush=True)
rc,o,e=run("which nvcc && nvcc --version")
print("rc",rc,o,e)
print("=== torch.cuda ===",flush=True)
try:
    import torch; print("torch",torch.__version__,"cuda avail:",torch.cuda.is_available(),"device:",torch.cuda.get_device_name(0) if torch.cuda.is_available() else "none")
except Exception as ex: print("torch err",ex)
print("=== github 可达 ===",flush=True)
print(run("curl -sI --max-time 15 https://github.com | head -1")[1])
print("=== dl.fbaipublicfiles (ESMFold权重源) 可达 ===",flush=True)
print(run("curl -sI --max-time 15 https://dl.fbaipublicfiles.com/fair-esm/models/esmfold_v1.pt | head -3")[1])
print("=== 试装 openfold from git (--no-build-isolation, 仅探测,限时300s) ===",flush=True)
rc,o,e=run("pip install --no-build-isolation git+https://github.com/aqlaboratory/openfold.git 2>&1 | tail -25")
print("install rc:",rc); print(o[-1500:]); print("ERR:",e[-500:])
print("=== 验证 openfold import ===",flush=True)
try:
    from openfold.data.data_transforms import make_atom14_masks
    from openfold.model.structure_module import StructureModule
    print("openfold import OK -> ESMFold可跑!")
except Exception as ex:
    print("openfold import 失败:",ex, "-> 需备选OmegaFold")
print("PROBE_DONE",flush=True)
