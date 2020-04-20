module yomhook

  save

  logical :: lhook = .false.

contains

  subroutine dr_hook(proc_name, iswitch, proc_key)

    use parkind1, only : jprb

    character(len=*), intent(in)    :: proc_name
    integer,          intent(in)    :: iswitch
    real(jprb),       intent(inout) :: proc_key
    ! Do nothing!

  end subroutine dr_hook

end module yomhook
