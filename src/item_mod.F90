module item_mod
  implicit none
  private
  public:: item

  type item
    private
    class(*),pointer:: var => null()
    type(item),pointer:: next => null()
  contains
    procedure,non_overridable:: next_item
    procedure,non_overridable:: get_item
  end type

  interface item
    module procedure item_constructor
  end interface
contains
  function item_constructor(var,next)
    class(*):: var
    class(item),pointer:: next
    class(item),pointer:: item_constructor

    allocate(item_constructor)
    allocate(item_constructor%var,source=var)
    item_constructor%next => next
  end function

  function next_item(self)
    class(item):: self
    class(item),pointer:: next_item

    next_item => self%next
  end function

  function get_item(self)
    class(item):: self
    class(*),pointer:: get_item

    get_item => self%var
  end function
end module
