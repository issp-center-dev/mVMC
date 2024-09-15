module wrapper
  use, intrinsic :: iso_c_binding
  implicit none
contains
  function c_charptr_to_f_charptr(ccp) result(result)
    type(c_ptr),intent(in),value :: ccp
    character(:,c_char),pointer :: result
    interface
       function strlen(p) bind(c)
         import c_ptr, c_size_t
         type(c_ptr),value :: p
         integer(c_size_t) strlen
       end function strlen
    end interface
    result => convert_cptr(ccp,strlen(ccp))
  contains
    function convert_cptr(p, len)
      type(c_ptr),intent(in) :: p
      integer(c_size_t),intent(in) :: len
      character(len, c_char),pointer :: convert_cptr
      call c_f_pointer(p, convert_cptr)
    end function convert_cptr
  end function c_charptr_to_f_charptr

  integer function ilaenv_wrap(ispec, name, opts, n1, n2, n3, n4) bind(c, name="ilaenv_wrap")
    implicit none
    integer,intent(in),value :: ispec, n1, n2, n3, n4
    type(c_ptr),intent(in),value :: name, opts
    character(:,c_char),pointer :: namef, optsf
    integer ilaenv

    namef => c_charptr_to_f_charptr(name)
    optsf => c_charptr_to_f_charptr(opts)

    ilaenv_wrap = ILAENV(ispec, namef, optsf, n1, n2, n3, n4)
  end function ilaenv_wrap
end module wrapper
