module engpar
  use :: iso_c_binding
  public
  interface
  subroutine cengpar_initialize() bind(C, NAME='cengpar_initialize')
    use :: iso_c_binding
  end subroutine
  end interface
end module
