module fluid_scheme_cns
  use field, only : field_t
  use fluid_base, only : fluid_base_t
  

  !> Base type of compressible fluid formulations.
  type, abstract, extends(fluid_base_t) :: fluid_scheme_cns_t
     !> The momentum field
     type(field_t), pointer :: m_x => null()    !< x-component of Momentum
     type(field_t), pointer :: m_y => null()    !< y-component of Momentum
     type(field_t), pointer :: m_z => null()    !< z-component of Momentum
     type(field_t), pointer :: E => null()    !< Total energy

  end type fluid_scheme_cns_t
end module fluid_scheme_cns