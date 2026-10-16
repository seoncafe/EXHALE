      module blas_thread_policy
      ! THE THREAD COUNT OF THE LAPACK LIBRARY THE BINARY LINKS. The GNU
      ! build links the OpenBLAS of the compiler's prefix (Makefile), whose
      ! own thread pool is independent of OpenMP: unset, it takes every
      ! logical CPU of the machine (72 here), which the OpenMP cap of
      ! init.f90 does not touch (review P3, 2026-09-12). The band
      ! factorizations of the stationary solve (dgbtrf on a 1500-row band)
      ! are far too small for such a team, so unless the run states
      ! OPENBLAS_NUM_THREADS itself the count is set to 1 here, at startup,
      ! through the library's own entry point resolved at run time
      ! (dlsym), so that a build against MKL (sequential by the Makefile)
      ! or the reference LAPACK, which has no such entry point, is left
      ! alone. The setup report states what was found and what was set.
      use iso_c_binding
      implicit none
      private
      public :: blas_threads_set_policy, blas_threads_report

      character(len=64), save :: blas_library_text = 'not resolved'
      character(len=96), save :: blas_policy_text  = 'not set'

      interface
         function c_dlopen(file, mode) bind(C, name='dlopen')            &
                  result(handle)
         import :: c_ptr, c_int, c_char
         character(kind=c_char), dimension(*) :: file
         integer(c_int), value :: mode
         type(c_ptr) :: handle
         end function c_dlopen
         function c_dlsym(handle, name) bind(C, name='dlsym')            &
                  result(addr)
         import :: c_ptr, c_funptr, c_char
         type(c_ptr), value :: handle
         character(kind=c_char), dimension(*) :: name
         type(c_funptr) :: addr
         end function c_dlsym
         subroutine openblas_set_num_threads_t(n) bind(C)
         import :: c_int
         integer(c_int), value :: n
         end subroutine openblas_set_num_threads_t
         function openblas_get_num_threads_t() bind(C) result(n)
         import :: c_int
         integer(c_int) :: n
         end function openblas_get_num_threads_t
      end interface

      contains

      subroutine blas_threads_set_policy()
      ! dlopen(NULL) hands back the running program, whose symbols include
      ! every shared library it was linked against; RTLD_NOW = 2 on glibc.
      type(c_ptr)    :: self
      type(c_funptr) :: fset, fget
      procedure(openblas_set_num_threads_t), pointer :: set_threads
      procedure(openblas_get_num_threads_t), pointer :: get_threads
      character(len=32) :: env
      integer :: n_before, n_after
      self = c_dlopen(c_null_char, 2_c_int)
      if (.not. c_associated(self)) then
         blas_library_text = 'dlopen of the program failed'
         return
      endif
      fset = c_dlsym(self, 'openblas_set_num_threads'//c_null_char)
      fget = c_dlsym(self, 'openblas_get_num_threads'//c_null_char)
      if (.not. c_associated(fset) .or. .not. c_associated(fget)) then
         blas_library_text = 'not OpenBLAS (no openblas_* entry points)'
         blas_policy_text  = 'left to the library'
         return
      endif
      call c_f_procpointer(fset, set_threads)
      call c_f_procpointer(fget, get_threads)
      blas_library_text = 'OpenBLAS'
      n_before = get_threads()
      call get_environment_variable('OPENBLAS_NUM_THREADS', env)
      if (len_trim(env) .gt. 0) then
         write(blas_policy_text,'(A,I0,A)') 'OPENBLAS_NUM_THREADS stated:', &
               n_before, ' thread(s), left as stated'
         return
      endif
      call set_threads(1_c_int)
      n_after = get_threads()
      write(blas_policy_text,'(A,I0,A,I0,A)') 'set to ', n_after,          &
            ' thread(s) (library default was ', n_before,                 &
            '; state OPENBLAS_NUM_THREADS to choose)'
      end subroutine blas_threads_set_policy

      subroutine blas_threads_report(unit)
      integer, intent(in) :: unit
      write(unit,'(A,A,A,A)') ' - LAPACK library: ', trim(blas_library_text), &
            '; BLAS threads ', trim(blas_policy_text)
      end subroutine blas_threads_report

      end module blas_thread_policy
