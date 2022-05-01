program main
  use pde_setting
  implicit none

  real t, tp, u
  integer output_index
  integer i, j
  integer isub, jsub
  integer adjTcellNum
  integer, allocatable :: seed(:)
  integer size

  call read_xdata()

  call random_seed(size=size)
  allocate(seed(size))
  seed(:) = iseed
  call random_seed(put=seed)

  open (unit = 100, file='./out/logfile', action="write")

  call pde_init()

  t = 0.0
  tp = 0.0
  output_index = 0

  do while (t < tend)
     if (t .ge. tp) then
        !call cell_stat(t)
        call output_to_file(output_index)
        output_index = output_index + 1
        tp = tp + tpinc
        write(*, '(a20, f10.2)') 'simulation time = ', t
     end if

     call step()
     t = t + dt
  end do

  close(unit=100)

end program main
