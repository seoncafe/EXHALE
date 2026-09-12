"""Focused audit: unchanged production routine bodies, isolated dependencies.

Generated fixtures and compiler products stay beside this script. This is not
an atmospheric integration test. Run with Python 3 and gfortran on PATH.
"""
from pathlib import Path
import importlib.util
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def routine(path, start, end):
    source = (ROOT / path).read_text()
    begin = source.lower().index(start.lower())
    finish = source.lower().index(end.lower(), begin) + len(end)
    return source[begin:finish]


sed_path = "src/modules/radiation/sed_read.f90"
chem_path = "src/modules/lower_atmosphere/diffusive_photochemistry.f90"
sed = routine(sed_path, "double precision function sed_band_integrated_flux", "end function sed_band_integrated_flux")
reader = routine(sed_path, "logical function sed_next_row", "end function sed_next_row")
chem = routine(chem_path, "subroutine equilibrate_chemistry_at_fixed_conserved_state(u,", "end subroutine equilibrate_chemistry_at_fixed_conserved_state")
driver = """module isolated
implicit none
integer, parameter :: dp=kind(1d0), N=1, Ng=1, n_species=1
integer :: eos_calls=0, scenario=1
character(len=256) :: sed_file
type ioniz_eq_ledger
integer :: n_nonfinite=0, n_offsimplex=0
end type
contains
""" + sed + "\n" + reader + "\n" + chem + """
subroutine pressure_and_temperature_at_fixed_conserved_state(u,f,p,t,ntot,ne)
real(dp), intent(in) :: u(3,0:2),f(0:2,1)
real(dp), intent(out) :: p(0:2),t(0:2),ntot(0:2),ne(0:2)
eos_calls=eos_calls+1
p=1d0; t=100d0; ntot=1d0; ne=0d0
if (scenario==1) t=100d0+dble(eos_calls)
end subroutine
subroutine ioniz_eq(t,rho,f,heat,cool,eta,ledger)
real(dp), intent(in) :: t(0:2),rho(0:2)
real(dp), intent(inout) :: f(0:2,1),heat(0:2),cool(0:2),eta(0:2)
type(ioniz_eq_ledger), intent(out) :: ledger
ledger=ioniz_eq_ledger()
if (scenario==2) ledger%n_offsimplex=1
end subroutine
end module
program probe
use isolated
implicit none
real(dp) :: u(3,0:2), f(0:2,1), p(0:2),t(0:2),h(0:2),c(0:2),eta(0:2)
integer :: cycles
logical :: ok
sed_file='linear_sed.txt'
write(*,'(A,ES24.16)') 'clipped_linear_flux=',sed_band_integrated_flux(950d0,1100d0)
sed_file='single_segment_sed.txt'
write(*,'(A,ES24.16)') 'single_segment_flux=',sed_band_integrated_flux(950d0,1100d0)
u=1d0; f=1d0; h=0d0; c=0d0; eta=0d0
call equilibrate_chemistry_at_fixed_conserved_state(u,f,p,t,h,c,eta,ok,cycles)
write(*,'(A,L1,A,I0,A,F10.6)') 'cycle_exhaustion_ok=',ok,', cycles=',cycles, &
 ', last_relative_increment=',1d0/t(1)
scenario=2; eos_calls=0
call equilibrate_chemistry_at_fixed_conserved_state(u,f,p,t,h,c,eta,ok,cycles)
write(*,'(A,L1,A,I0)') 'off_simplex_ok=',ok,', cycles=',cycles
end program
"""
(HERE / "isolated_production_routines.f90").write_text(driver)
(HERE / "linear_sed.txt").write_text("900 900\n1000 1000\n1100 1100\n1200 1200\n")
(HERE / "single_segment_sed.txt").write_text("900 900\n1200 1200\n")
subprocess.run(["gfortran", "-O0", "-fcheck=all", "-ffree-line-length-none", "isolated_production_routines.f90", "-o", "isolated_probe"], cwd=HERE, check=True)
result = subprocess.run([str(HERE / "isolated_probe")], cwd=HERE, check=True, text=True, capture_output=True)
results = {"isolated_fortran_stdout": result.stdout, "exact_linear_band_integral": (1100**2-950**2)/2}

mapper = ROOT / "src/utils/map_state_to_grid.py"
spec = importlib.util.spec_from_file_location("production_mapper", mapper)
module = importlib.util.module_from_spec(spec)
sys.dont_write_bytecode = True
spec.loader.exec_module(module)

# Two valid, positive, monotonic but mutually inconsistent input radius arrays.
src = HERE / "mismatched_source"
src.mkdir(exist_ok=True)
(src / "Hydro_ioniz.txt").write_text("# coupling mode=phys t_phys=123 certified=T\n# columns r rho v[km/s] T[K]\n1 10 1 100\n2 5 1 100\n")
(src / "Ion_species.txt").write_text("# columns r HI HII\n1.2 9 1\n2.2 2.5 2.5\n")
target = HERE / "target.txt"
target.write_text("1 0\n3 0\n")
mapped = HERE / "mapped"
run = subprocess.run([sys.executable, str(mapper), str(src), str(target), str(mapped)], text=True, capture_output=True)
results["mapper_exit_code"] = run.returncode
results["mapper_stdout"] = run.stdout
results["mapper_stderr"] = run.stderr
if run.returncode == 0:
    results["mapped_hydro"] = module.read(mapped / "Hydro_ioniz.txt")[1].tolist()
    results["mapped_species"] = module.read(mapped / "Ion_species.txt")[1].tolist()
results["limitations"] = ["Fortran chemistry dependencies are deterministic stubs; no atmospheric occurrence frequency is measured.", "Mapper fixtures test input validation, not a physically calibrated atmosphere."]
results["separate_mapper_cases"] = {}
for name, species, radii in [
    ("mismatched_radii_only", "1.2 9 1\n2.2 2.5 2.5\n", "1 0\n2 0\n"),
    ("out_of_range_only", "1 9 1\n2 2.5 2.5\n", "1 0\n3 0\n"),
    ("interior_remap_metadata", "1 9 1\n2 2.5 2.5\n", "1 0\n1.5 0\n"),
]:
    case_src = HERE / name / "source"
    case_src.mkdir(parents=True, exist_ok=True)
    (case_src / "Hydro_ioniz.txt").write_text((src / "Hydro_ioniz.txt").read_text())
    (case_src / "Ion_species.txt").write_text("# columns r HI HII\n" + species)
    case_target = HERE / name / "target.txt"
    case_target.write_text(radii)
    case_out = HERE / name / "output"
    case_run = subprocess.run([sys.executable, str(mapper), str(case_src), str(case_target), str(case_out)], text=True, capture_output=True)
    results["separate_mapper_cases"][name] = {"exit_code": case_run.returncode, "stderr": case_run.stderr}
    if case_run.returncode == 0:
        header, data = module.read(case_out / "Hydro_ioniz.txt")
        results["separate_mapper_cases"][name].update(header=header, hydro=data.tolist())
(HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
