"""Execute unchanged scalar PLM routines with a manufactured nonuniform grid.

This isolates reconstruction correctness and loop overhead, not full solver
scaling. Generated source, compiler products, and measurements stay here.
"""
from pathlib import Path
import ctypes
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
source = (ROOT / "src/modules/states/PLM_rec.f90").read_text()


def extract(start, end):
    a = source.index(start)
    b = source.index(end, a) + len(end)
    return source[a:b]


body = extract("real*8 function minmod_mc_slope", "end function minmod_mc_slope") + "\n" + extract("subroutine PLM_rec_scalar", "end subroutine PLM_rec_scalar")
program = """! Extracted production routines for direct comparison; bodies unchanged.
module scalar_reconstruction_probe
implicit none
integer :: N=500
integer, parameter :: Ng=2
real*8, allocatable :: r(:)
contains
""" + body + """
end module
program probe
use scalar_reconstruction_probe
use omp_lib
implicit none
real*8, allocatable :: q(:),ql(:),qr(:),left_reference(:),right_reference(:)
real*8 :: start_time,elapsed,difference,elapsed_best
integer :: j,k,rep,nt,team,case_id
integer, parameter :: repeats=3000
call omp_set_dynamic(.false.)
do case_id=1,2
N=500
if (case_id==2) N=8000
allocate(r(1-Ng:N+Ng),q(1-Ng:N+Ng),ql(1-Ng:N+Ng),qr(1-Ng:N+Ng), &
left_reference(1-Ng:N+Ng),right_reference(1-Ng:N+Ng))
do j=1-Ng,N+Ng
r(j)=exp(2d0*dble(j)/dble(N))
q(j)=1d0+0.2d0*sin(10d0*r(j))
enddo
call omp_set_num_threads(1)
call PLM_rec_scalar(q,left_reference,right_reference)
do k=0,4
nt=2**k
call omp_set_num_threads(nt)
team=0
!$omp parallel
!$omp single
team=omp_get_num_threads()
!$omp end single
!$omp end parallel
elapsed_best=huge(1d0)
difference=0d0
do rep=1,3
start_time=omp_get_wtime()
do j=1,repeats
call PLM_rec_scalar(q,ql,qr)
enddo
elapsed=omp_get_wtime()-start_time
elapsed_best=min(elapsed_best,elapsed)
difference=max(difference,maxval(abs(ql-left_reference)),maxval(abs(qr-right_reference)))
enddo
write(*,'(I6,1X,I3,1X,I3,1X,ES15.7,1X,ES15.7)') N,nt,team,elapsed_best/repeats,difference
enddo
deallocate(r,q,ql,qr,left_reference,right_reference)
enddo
end program
"""
(HERE / "scalar_plm_probe.f90").write_text(program)
subprocess.run(["gfortran", "-O3", "-fopenmp", "scalar_plm_probe.f90", "-o", "scalar_plm_probe"], cwd=HERE, check=True)
env = dict(os.environ, OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close", OMP_PLACES="cores", OPENBLAS_NUM_THREADS="1")
run = subprocess.run([str(HERE / "scalar_plm_probe")], cwd=HERE, env=env, check=True, text=True, capture_output=True, timeout=120)
measurements = []
for line in run.stdout.splitlines():
    n, requested, actual, seconds, error = line.split()
    measurements.append(dict(cells=int(n), requested_threads=int(requested), actual_threads=int(actual), seconds_each=float(seconds), max_abs_difference=float(error)))
library = ctypes.CDLL("/opt/miniconda3/lib/libopenblas.so.0")
library.openblas_get_config.restype = ctypes.c_char_p
library.openblas_get_num_threads.restype = ctypes.c_int
library.openblas_get_parallel.restype = ctypes.c_int
results = dict(measurements=measurements, stdout=run.stdout,
    openblas_config=library.openblas_get_config().decode(),
    openblas_default_threads=library.openblas_get_num_threads(),
    openblas_parallel_kind=library.openblas_get_parallel(),
    inherited_thread_environment={key: os.environ.get(key) for key in ["OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PROC_BIND", "OPENBLAS_NUM_THREADS"]},
    cpu_affinity_count=len(os.sched_getaffinity(0)),
    compiler=subprocess.check_output(["gfortran", "--version"],text=True).splitlines()[0],
    note="Best of three batches; each batch has 3000 calls. Timings are local microbenchmarks, not full EXHALE speedups. OpenBLAS query uses the parent environment, separate from the pinned probe environment.")
(HERE / "results.json").write_text(json.dumps(results, indent=2)+"\n")
print(json.dumps(results, indent=2))
