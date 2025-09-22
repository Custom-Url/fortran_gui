# execute.py
import os
import glob
import subprocess
from pathlib import Path
import threading
import queue
import signal
import time

# ---------------------------
# Globals for process tracking
# ---------------------------
_current_proc = None
_current_proc_lock = threading.Lock()
_stop_event = threading.Event()
_log_queue = queue.Queue()
_hard_pause_state = False

# ---------------------------
# Logging helper
# ---------------------------
def log(msg: str):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    _log_queue.put(f"[{timestamp}] {msg}")

# ---------------------------
# Pause/Resume AMITEX
# ---------------------------
def pause_amitex():
    global _hard_pause_state
    with _current_proc_lock:
        proc = _current_proc
    if proc is None:
        log("No running process to pause.")
        return
    try:
        pgid = os.getpgid(proc.pid)
        os.killpg(pgid, signal.SIGSTOP)
        _hard_pause_state = True
        log("AMITEX paused (mid-simulation).")
    except Exception as e:
        log(f"Error pausing: {e}")

def resume_amitex():
    global _hard_pause_state
    with _current_proc_lock:
        proc = _current_proc
    if proc is None:
        log("No process to resume.")
        return
    try:
        pgid = os.getpgid(proc.pid)
        os.killpg(pgid, signal.SIGCONT)
        _hard_pause_state = False
        log("AMITEX resumed.")
    except Exception as e:
        log(f"Error resuming: {e}")


# ---------------------------
# Stop AMITEX
# ---------------------------
def stop_amitex():
    global _current_proc
    _stop_event.set()
    with _current_proc_lock:
        proc = _current_proc
    if proc is None:
        log("No AMITEX process running.")
        return

    try:
        pgid = os.getpgid(proc.pid)
        log(f"Stopping process group {pgid}...")
        os.killpg(pgid, signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            log("Process did not exit after SIGTERM; sending SIGKILL...")
            os.killpg(pgid, signal.SIGKILL)
            proc.wait()
        log("AMITEX process stopped.")
    except Exception as e:
        log(f"Error stopping process: {e}")
    finally:
        with _current_proc_lock:
            _current_proc = None

# ---------------------------
# Run AMITEX loop
# ---------------------------
def run_amitex_loop(
    vtk_glob_pattern="FFT/micr/code_python/mesh/L*/h0.0002/iUC*.vtk",
    mate="FFT/mate/mate_PF.xml",
    load="FFT/load/char.xml",
    algo="FFT/algo/algo_Miehe2.xml",
    amitex_exec="/mnt/data/dg765/micro/FFT/amitex_fftp/libAmitex/src/user_Miehe2/amitex_fftp",
    mpirun=["mpirun", "-np", "18"],
    resu_root="resu/dg765",
):
    global _current_proc
    os.environ["LD_LIBRARY_PATH"] = (
        "/mnt/data/dg765/micro/FFT/amitex_fftp/libAmitex/lib:"
        + os.environ.get("LD_LIBRARY_PATH", "")
    )
    vtk_files = sorted(glob.glob(vtk_glob_pattern))
    if not vtk_files:
        log(f"No VTK files found: {vtk_glob_pattern}")
        return
    log(f"Found {len(vtk_files)} VTK files.")

    try:
        for vtk_file in vtk_files:
            if _stop_event.is_set():
                log("Stop requested. Exiting loop.")
                break

            vtk_path = Path(vtk_file)
            spacing_dir = vtk_path.parent.name
            vf_dir = vtk_path.parent.parent.name
            vtk_basename = vtk_path.stem
            resu_dir = Path(resu_root) / vf_dir / spacing_dir / vtk_basename
            resu_dir.mkdir(parents=True, exist_ok=True)
            log(f"Running AMITEX for {vtk_file}")

            cmd = mpirun + [
                amitex_exec,
                "-nm", str(vtk_file),
                "-m", mate,
                "-c", load,
                "-a", algo,
                "-s", str(resu_dir / "Load0.0"),
            ]

            log("Launching: " + " ".join(cmd))
            with _current_proc_lock:
                _current_proc = subprocess.Popen(cmd, preexec_fn=os.setsid, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            proc = _current_proc
            if proc.stdout:
                for line in iter(proc.stdout.readline, b""):
                    log(line.decode(errors="ignore").rstrip())
            proc.wait()
            log(f"Process for {vtk_basename} finished with code {proc.returncode}")
            with _current_proc_lock:
                _current_proc = None

    finally:
        log("AMITEX run loop finished.")
        _stop_event.clear()

# ---------------------------
# Queue getter (for GUI polling)
# ---------------------------
def get_log_queue():
    return _log_queue

