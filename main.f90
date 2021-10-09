program main
  use setting
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

  call init_cell_pool()

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

     do i = 1, Lbox
        do j = 1, Lbox
           if ( cmat(i,j)%type == 0 ) then
              ! empty
              call random_number(u)
              if ( u < a(i,j)*dt ) then
                 cmat(i,j)%type = 2                                     
              end if
           else if ( cmat(i,j)%type .eq. 1 ) then
              ! M cell
              adjTcellNum = 0              
              do isub = i-Rc, i+Rc
                 if (isub > 0 .and. isub <= Lbox) then
                    do jsub = j-Rc, j+Rc
                       if (jsub > 0 .and. jsub <= Lbox) then
                          if ( (cmat(isub, jsub)%type .eq. 2).or.(cmat(isub, jsub)%type .eq. 3) ) then
                             adjTcellNum = adjTcellNum + 1
                          end if
                       end if
                    end do
                 end if
              end do
              if ( adjTcellNum > nc ) then
                 cmat(i,j)%type = 0
                 ! activate T cell
                 do isub = i-Rc, i+Rc
                    if (isub > 0 .and. isub <= Lbox) then
                       do jsub = j-Rc, j+Rc
                          if (jsub > 0 .and. jsub <= Lbox) then
                             if ( cmat(isub, jsub)%type == 2 ) then
                                cmat(isub, jsub)%type = 3
                             end if
                          end if
                       end do
                    end if
                 end do
              end if
           else if ( cmat(i,j)%type == 3 ) then
              ! activated T cell
              call random_number(u)
              if ( u < beta*dt ) then
                 cmat(i,j)%type = 0
              end if
           else if ( cmat(i,j)%type == 2 ) then
              !else if ( cmat(i,j)%type == 2 .or. cmat(i,j)%type == 3) then
              ! T cell
              call random_number(u)
              if ( u < beta*dt ) then
                 cmat(i,j)%type = 0
              end if
           else
              ! not suppose to happen
              !write(*, *) 'ERROR: unknown cell type'
              !stop
           end if
        end do
     end do

     call update_lambda()
     call update_phi(dt)
     call update_p()     
     call update_rate()
     t = t + dt
  end do

  close(unit=100)

end program main
