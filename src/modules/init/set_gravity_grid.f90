      module gravity_grid_construction
      ! Pre-evaluate gravity on cell centers and at interfaces
      
      use global_parameters
      use grav_func
      
      contains
      
      subroutine set_gravity_grid
      ! Subroutine to evaluate the gravitational potential on
      !    cell centers and edges
            
      ! Evaluate the gravitational potential at cell centers.
      Gphi_c = (/ ( phi(r(j)), j = 1-Ng,N+Ng) /)
      
      ! Evaluate the gravitational potential at cell interfaces.
      Gphi_i = (/ ( phi(r_edg(j)), j = 1-Ng,N+Ng) /)
      
      end subroutine set_gravity_grid
      
      end module gravity_grid_construction
