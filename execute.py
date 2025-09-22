import os
import glob
import subprocess
from pathlib import Path
import signal
import sys

def run_amitex():
    # ---- Environment ----
    os.environ["LD_LIBRARY_PATH"] = (
        "/mnt/data/dg765/micro/FFT/amitex_fftp/libAmitex/lib:"
        + os.environ.get("LD_LIBRARY_PATH", "")
    )

    # ---- Constants ----
    AMITEX = "/mnt/data/dg765/micro/FFT/amitex_fftp/libAmitex/src/user_Miehe2/amitex_fftp"
    MPIRUN = ["mpirun", "-np", "18"]

    MATE = "FFT/mate/mate_PF.xml"
    LOAD = "FFT/load/char.xml"
    ALGO = "FFT/algo/algo_Miehe2.xml"

    # ---- Find VTK files ----
    vtk_files = glob.glob("FFT/micr/code_python/mesh/L*/h0.0002/iUC*.vtk")

    current_proc = None  # Track running process

    try:
        for vtk_file in vtk_files:
            vtk_path = Path(vtk_file)

            spacing_dir = vtk_path.parent.name        # e.g. h0.0002
            vf_dir = vtk_path.parent.parent.name      # e.g. L0.05_vf0.35
            vtk_basename = vtk_path.stem              # e.g. iUC100_vpMIN0.09

            resu_dir = Path("resu/dg765") / vf_dir / spacing_dir / vtk_basename
            resu_dir.mkdir(parents=True, exist_ok=True)

            print(f"\nRunning AMITEX for: {vtk_file}")
            print(f"Results in: {resu_dir}")

            output_prefix = resu_dir / "Load0.0"

            # Build command
            cmd = MPIRUN + [
                AMITEX,
                "-nm", str(vtk_file),
                "-m", MATE,
                "-c", LOAD,
                "-a", ALGO,
                "-s", str(output_prefix),
            ]

            # Start subprocess (instead of run)
            current_proc = subprocess.Popen(cmd)

            # Wait for it to finish
            current_proc.wait()
            current_proc = None  # Reset after success

    except KeyboardInterrupt:
        print("\nInterrupted. Killing child processes...")
        if current_proc is not None:
            current_proc.terminate()
            try:
                current_proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                current_proc.kill()
        sys.exit(1)

if __name__ == "__main__":
    run_amitex()

