#!/usr/bin/env bash
# Run from EXHALE_v1.00. No production build or output is modified.
set -euo pipefail
audit_dir=docs/audit_20260905
review_build_dir=$(mktemp -d /tmp/exhale-issues-review-20260909.XXXXXX)
{
  sed -n '1,$p' "$audit_dir/issues_review_prefix_20260909.f90"
  sed -n '/^      subroutine pgmres(/,/^      end subroutine pgmres/p' \
    src/modules/time_step/steady_newton.f90
  sed -n '/^      subroutine dogleg_step(/,/^      end subroutine dogleg_step/p' \
    src/modules/time_step/steady_newton.f90
  sed -n '/^[[:space:]]*pure real\*8 function enthalpy_flux_term_ratio(/,/^[[:space:]]*end function enthalpy_flux_term_ratio/p' \
    src/modules/post_process/post_process_adv.f90
  sed -n '1,$p' "$audit_dir/issues_review_suffix_20260909.f90"
} | gfortran -x f95 -ffree-form -ffree-line-length-none -O0 -fcheck=all \
    -J"$review_build_dir" - -o "$review_build_dir/review.x"
"$review_build_dir/review.x"
printf 'Isolated build retained at %s\n' "$review_build_dir"
