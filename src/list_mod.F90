module list_mod
  use item_mod

  implicit none
  private
  public list

  type,abstract:: list
    private
    class(item), pointer:: first_item => null()
    class(item), pointer:: current_item => null()
    class(item), pointer:: last_item => null()
  contains
    procedure:: add_item
    procedure:: get_item
  end type

contains

  subroutine add_item(self, var)
    class(list):: self
    class(*):: var
    class(item), pointer:: new_item

    if (.not. associated(self%first_item)) then
      self%first_item => item(var, null())
      self%current_item => self%first_item
      self%last_item => self%first_item
    else
      new_item => item(var, null())
      call self%last_item%set_next_item(new_item)
      self%last_item => new_item
    end if
  end subroutine

  function get_item(self)
    class(list):: self
    class(*), pointer:: get_item

    get_item => self%current_item%get_item()
  end function
end module
