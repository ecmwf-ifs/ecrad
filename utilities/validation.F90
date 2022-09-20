module validation
    use parkind1, only: jprb, jpim

    implicit none

    private
    public :: error_type

    type :: error_type
        real(kind=jprb) :: zminval  ! Minimum value of field
        real(kind=jprb) :: zmaxval  ! Maximum value of field
        real(kind=jprb) :: zmaxerr  ! Maximum error
        real(kind=jprb) :: zavgerr  ! Average error
        real(kind=jprb) :: zsumerr  ! Total error accumulated
        real(kind=jprb) :: zrelerr  ! Relative error
        logical :: pass
    contains
        procedure :: validate_r1
        procedure :: validate_r2
        procedure :: validate_r3
        generic :: validate => validate_r1, validate_r2, validate_r3
        procedure, private :: compute_relative_error
        procedure :: print => error_type_print
        procedure :: print_header => error_type_print_header
    end type error_type

contains

    subroutine error_type_print_header(self)
        class(error_type), intent(in) :: self

        print '(1x,a30,1x,a3,5(1x,a20))', &
            & 'Variable','Dim', 'MinValue','MaxValue','AbsMaxErr','AvgAbsErr/GP','MaxRelErr-%'
    end subroutine error_type_print_header

    subroutine error_type_print(self, name, dims)
        class(error_type), intent(in) :: self
        character(*), intent(in) :: name
        integer(kind=jpim), intent(in) :: dims
        character(len=5) :: clwarn
        real(kind=jprb) :: zrelerr

        !-- If you get 4 exclamation marks next to your error output,
        !   then it is likely that some uninitialized variables exists or
        !   some other screw-up -- watch out this !!!!
        if (self%pass) then
            clwarn = ' '
        else
            clwarn = ' !!!!'
        endif
        zrelerr = 100.0_JPRB * self%zrelerr

        print "(1x,a30,1x,i1,'D',5(1x,e20.13),a)", name, dims, &
         & self%zminval, self%zmaxval, self%zmaxerr, self%zavgerr, self%zrelerr, clwarn
    end subroutine error_type_print

    subroutine compute_relative_error(self, zsum)
        class(error_type), intent(inout) :: self
        real(kind=jprb), intent(in) :: zsum
        real(kind=jprb) :: zeps = epsilon(1.0_jprb)

        if (self%zsumerr < zeps) then
            self%zrelerr = 0.0_JPRB
        elseif (zsum < zeps) then
            self%zrelerr = self%zsumerr / (1.0_jprb + zsum)
        else
            self%zrelerr = self%zsumerr / zsum
        endif

        self%pass = self%zrelerr <= 10.0_JPRB * zeps
    end subroutine compute_relative_error

    subroutine validate_r1(self, name, ref, field, dim1, is_print)
        class(error_type), intent(inout) :: self
        character(*), intent(in) :: name
        real(kind=jprb), intent(in) :: ref(:), field(:)
        integer(kind=jpim), intent(in) :: dim1
        logical, intent(in), optional :: is_print
        integer(kind=jpim) :: jl
        real(kind=jprb) :: zdiff, zrefsum  ! Sum of all entries in the reference data

        self%zminval = minval(field(:))
        self%zmaxval = maxval(field(:))

        self%zmaxerr = 0.0_jprb
        self%zsumerr = 0.0_jprb
        zrefsum = 0.0_jprb

        do jl=1,dim1
            zdiff = abs(field(jl) - ref(jl))
            self%zmaxerr = max(self%zmaxerr, zdiff)
            self%zsumerr = self%zsumerr + zdiff
            zrefsum = zrefsum + abs(ref(jl))
        end do

        self%zavgerr = self%zsumerr / real(dim1,jprb)

        call self%compute_relative_error(zrefsum)

        if (present(is_print)) then
            if (is_print) call self%print(name, 1)
        end if
    end subroutine validate_r1

    subroutine validate_r2(self, name, ref, field, dim1, dim2, ncol, is_print)
        class(error_type), intent(inout) :: self
        character(*), intent(in) :: name
        real(kind=jprb), intent(in) :: ref(:,:), field(:,:)
        integer(kind=jpim), intent(in) :: dim1, dim2, ncol
        logical, intent(in), optional :: is_print
        integer(kind=jpim) :: jl, jk
        real(kind=jprb) :: zdiff, zrefsum  ! Sum of all entries in the reference data

        self%zminval = minval(field(:,:))
        self%zmaxval = maxval(field(:,:))

        self%zmaxerr = 0.0_jprb
        self%zsumerr = 0.0_jprb
        zrefsum = 0.0_jprb

        do jl=1,dim2
            do jk=1,dim1
                zdiff = abs(field(jk,jl) - ref(jk,jl))
                self%zmaxerr = max(self%zmaxerr, zdiff)
                self%zsumerr = self%zsumerr + zdiff
                zrefsum = zrefsum + abs(ref(jk,jl))
            end do
        end do

        self%zavgerr = self%zsumerr / real(ncol,jprb)

        call self%compute_relative_error(zrefsum)

        if (present(is_print)) then
            if (is_print) call self%print(name, 2)
        end if
    end subroutine validate_r2

    subroutine validate_r3(self, name, ref, field, dim1, dim2, dim3, ncol, is_print)
        class(error_type), intent(inout) :: self
        character(*), intent(in) :: name
        real(kind=jprb), intent(in) :: ref(:,:,:), field(:,:,:)
        integer(kind=jpim), intent(in) :: dim1, dim2, dim3, ncol
        logical, intent(in), optional :: is_print
        integer(kind=jpim) :: jl, jk, jm
        real(kind=jprb) :: zdiff, zrefsum  ! Sum of all entries in the reference data

        self%zminval = minval(field(:,:,:))
        self%zmaxval = maxval(field(:,:,:))

        self%zmaxerr = 0.0_jprb
        self%zsumerr = 0.0_jprb
        zrefsum = 0.0_jprb

        do jl=1,dim3
            do jk=1,dim2
                do jm=1,dim1
                    zdiff = abs(field(jm,jk,jl) - ref(jm,jk,jl))
                    self%zmaxerr = max(self%zmaxerr, zdiff)
                    self%zsumerr = self%zsumerr + zdiff
                    zrefsum = zrefsum + abs(ref(jm,jk,jl))
                end do
            end do
        end do

        self%zavgerr = self%zsumerr / real(ncol,jprb)

        call self%compute_relative_error(zrefsum)

        if (present(is_print)) then
            if (is_print) call self%print(name, 3)
        end if
    end subroutine validate_r3
end module validation
