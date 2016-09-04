module list_mod
  use item_mod

  implicit none
  private
  public list

  type,abstract:: list
    private
    class(item), pointer:: first_item => null()
    class(item), pointer:: current_item => null()
  contains
    procedure,non_overridable:: add_item
    procedure,non_overridable:: get_item
    procedure,non_overridable:: next
    procedure,non_overridable:: moreitems
    procedure,non_overridable:: reset
  end type
contains
  subroutine add_item(self, var)
    class(list):: self
    class(*):: var
    class(item), pointer:: new_item

    if (.not.associated(self%first_item)) then
      self%first_item => item(var, null())
    else
      new_item => item(var,self%first_item)
      self%first_item => new_item
    end if
  end subroutine

  function get_item(self)
    class(list):: self
    class(*), pointer:: get_item

    get_item => self%current_item%get_item()
  end function

  subroutine next(self)
    class(list) :: self

    self%current_item => self%current_item%next_item()
  end subroutine

  function moreitems(self)
    class(list) :: self
    logical moreitems

    moreitems = associated(self%current_item)
  end function

  subroutine reset(self)
    class(list) :: self

    self%current_item => self%first_item
  end subroutine
end module
