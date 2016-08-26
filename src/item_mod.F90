module item_mod

  implicit none
  private
  public:: item

  type item
    private
    character(len=*)    :: name
    class(*),   pointer :: value => null()
    type(item), pointer :: next  => null()
  contains
    procedure, non_overridable:: next_item
    procedure, non_overridable:: set_next_item
    procedure, non_overridable:: get_item
  end type

  interface item
    module procedure item_constructor
  end interface

  function item_constructor(name, value, next)
    character(len=*):: name
    class(*):: value
    class(item), pointer:: next
    class(item), pointer:: item_constructor

    allocate(item_constructor)
    item_constructor => next
    allocate(item_constructor%name,  source=name)
    allocate(item_constructor%value, source=value)
  end function

  function next_item(self)
    class(item):: self
    class(item), pointer:: next_item

    next_item => self%next
  end function

  subroutine set_next_item(self, next)
    class(item) :: self
    class(item), pointer :: next

    self%next => next
  end subroutine

  function get_item(self)
    class(item) :: self
    class(*), pointer :: get_item

    get_item => self%value
  end function getValue

end module
