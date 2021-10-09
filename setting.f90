module setting
  integer Lbox
  integer :: iseed
  real :: tend, dt
  real :: alpha, beta, Amax
  integer :: Rc, nc
  real :: diff, lambda, gamma, phi_tiny
  real :: k_lambda
  integer :: model_type, read_lambda_from_file

  namelist /xdata/ Lbox, tend, dt, alpha, beta, Amax, diff, k_lambda, &
       iseed, tpinc, Rc, nc, lambda, gamma, phi_tiny, model_type, &
       read_lambda_from_file

  type cell
     integer type
     ! type == 0 : empty splot
     ! type == 1 : M cell
     ! type == 2 : T cell
     ! type == 3 : activated T cell
  end type cell
  type(cell), allocatable :: cmat(:,:)
  real, allocatable ::  phi(:,:)
  real, allocatable ::  phi_old(:,:)
  real, allocatable :: p(:,:)
  real, allocatable :: a(:,:)
  real, allocatable :: fb_lambda(:,:)
  real, allocatable :: lambda_field(:,:)  

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    write(*, *) 'Control parameters...'
    write(*, '(a20, i10)') 'Lbox = ', Lbox
    write(*, '(a20, i10)') 'iseed = ', iseed
    write(*, '(a20, f10.4)') 'tend = ', tend
    write(*, '(a20, f10.4)') 'dt = ', dt
    write(*, '(a20, f10.4)') 'alpha = ', alpha
    write(*, '(a20, f10.4)') 'beta = ', beta
    write(*, '(a20, f10.4)') 'Amax = ', Amax    
    write(*, '(a20, f10.4)') 'diff = ', diff
    write(*, '(a20, f10.4)') 'tpinc = ', tpinc
    write(*, '(a20, i10)') 'Rc = ', Rc
    write(*, '(a20, f10.4)') 'lambda = ', lambda
    write(*, '(a20, f10.4)') 'gamma = ', gamma
    write(*, '(a20, f10.4)') 'phi_tiny = ', phi_tiny
    write(*, '(a20, f10.4)') 'k_lambda = ', k_lambda
    write(*, '(a20, i10)') 'nc = ', nc
    write(*, '(a20, i10)') 'model_type = ', model_type
    write(*, '(a20, i10)') 'read_lambda_from_file = ', read_lambda_from_file

    open(9, file="out/control.csv")
    write(9, '(a20, a10)') 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)') 'Lbox,', Lbox
    write(9, '(a20, i10)') 'iseed,', iseed
    write(9, '(a20, f10.4)') 'tend,', tend
    write(9, '(a20, f10.4)') 'dt,', dt
    write(9, '(a20, f10.4)') 'alpha,', alpha
    write(9, '(a20, f10.4)') 'beta,', beta
    write(9, '(a20, f10.4)') 'Amax,', Amax    
    write(9, '(a20, f10.4)') 'diff,', diff
    write(9, '(a20, f10.4)') 'tpinc,', tpinc
    write(9, '(a20, i10)') 'Rc,', Rc
    write(9, '(a20, f10.4)') 'lambda,', lambda
    write(9, '(a20, f10.4)') 'gamma,', gamma
    write(9, '(a20, f10.4)') 'phi_tiny,', phi_tiny
    write(9, '(a20, f10.2)') 'k_lambda,', k_lambda    
    write(9, '(a20, i10)') 'nc,', nc
    write(9, '(a20, i10)') 'model_type,', model_type
    write(9, '(a20, i10)') 'read_lambda_from_file,', read_lambda_from_file

    close(8)
    close(9)
  end subroutine read_xdata

  subroutine init_cell_pool()
    implicit none
    real u1, u2
    integer i, j, ishift, jshift
    integer curb

    write(*, '(A)', advance='no') 'Initialize...'
    allocate(cmat(1:Lbox,1:Lbox))
    allocate(phi(1:Lbox, 1:Lbox))
    allocate(p(1:Lbox,1:Lbox))
    allocate(a(1:Lbox,1:Lbox))
    allocate(fb_lambda(1:Lbox,1:Lbox))
    allocate(lambda_field(1:Lbox,1:Lbox))
    allocate(phi_old(1:Lbox,1:Lbox))

    if ( read_lambda_from_file == 1 ) then
       ! read lambda_field
       open (unit = 81, file="matrix.txt", action="read")
       do i = 1, Lbox
          read(81, *) lambda_field(i,:)
          !pause
       end do
       close(81)
       lambda_field(i,j) = k_lambda*lambda_field(i,j)
    else if ( read_lambda_from_file == 0 ) then
       lambda_field = lambda
    else
       print *, "Parameter error: read_lambda_from_file can only be 0 or 1."
    end if

    ! maximum total amount of T cell is limited by Amax
    Amax = Amax * Lbox * Lbox * alpha
    
    do i = 1, Lbox
       do j = 1, Lbox
          cmat(i,j)%type = 0 ! No cell everywhere
          phi(i,j) = 0.0
          p(i,j) = 0.0
          a(i,j) = 0.0
          fb_lambda(i,j) = 0.0
       end do
    end do

    select case (model_type)
    case (1) ! central line
       cmat( (Lbox/4 - 5) : (Lbox/4 + 5), (Lbox/2 - 5) : (Lbox/2 + 5) )%type = 3
       cmat( (2*Lbox/4 - 5) : (2*Lbox/4 + 5), (Lbox/2 - 5) : (Lbox/2 + 5) )%type = 3
       cmat( (3*Lbox/4 - 5) : (3*Lbox/4 + 5), (Lbox/2 - 5) : (Lbox/2 + 5) )%type = 3          
       
    case (2) ! central line
       cmat( (Lbox/3 - 5) : (Lbox/3 + 5), (Lbox/8 - 5) : (Lbox/8 + 5) )%type = 3
       cmat( (Lbox/3 - 5) : (Lbox/3 + 5), (7*Lbox/8 - 5) : (7*Lbox/8 + 5) )%type = 3
       
    case (3) ! random
       
    end select

    ! cmat(170:180, 251:261)%type = 3
    ! cmat(320:330, 251:261)%type = 3
    ! cmat(470:480, 251:261)%type = 3
    
    ! ! for 512 two side
    ! cmat(320:330, 51:61)%type = 3
    ! cmat(320:330, 451:461)%type = 3    
    
    curb = 10
    ! randomly distribute M cells
    do i = curb, Lbox-curb, 4
       do j = curb, Lbox-curb, 4
          call random_number(u1)
          call random_number(u2)
          u1 = u1 - 0.5
          u2 = u2 - 0.5
          if ( u1>0.25 ) then
             ishift = i + 1
          else if ( u1<-0.25 ) then
             ishift = i - 1
          else
             ishift = i
          end if
          if ( u2>0.25 ) then
             jshift = j + 1
          else if ( u2<-0.25 ) then
             jshift = j - 1
          else
             jshift = j
          end if
          cmat(ishift,jshift)%type = 1
       end do
    end do
    write(*, *) 'Done.'
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i, j

    if ( index==0 ) then
       WRITE(filename,'(A7,I5.5,A4)') './out/f', index, '.dat'
       open (unit = 8, file=filename, action="write")
       do i = 1, Lbox
          do j = 1, Lbox-1
             write(8, '(F8.4, A2)', advance="no") lambda_field(i,j), ', '
          end do
          write(8, '(F8.4, A2)') lambda_field(i,j)
       end do
       close(8)
    end if
    
    WRITE(filename,'(A7,I5.5,A4)') './out/c', index, '.dat'
    open (unit = 11, file=filename, action="write")
    WRITE(filename,'(A12,I5.5,A4)') './out/lambda', index, '.dat'
    open (unit = 21, file=filename, action="write")
    WRITE(filename,'(A9,I5.5,A4)') './out/phi', index, '.dat'
    open (unit = 31, file=filename, action="write")
    WRITE(filename,'(A7,I5.5,A4)') './out/a', index, '.dat'
    open (unit = 41, file=filename, action="write")


    do i = 1, Lbox
       do j = 1, Lbox-1
          write(11, '(I5, A2)', advance="no") cmat(i,j)%type, ', '
          write(21, '(F8.4, A2)', advance="no") fb_lambda(i,j), ', '
          write(31, '(F8.4, A2)', advance="no") phi(i,j), ', '
          write(41, '(F8.4, A2)', advance="no") a(i,j), ', '
       end do
       write(11, '(I5)') cmat(i,Lbox)%type
       write(21, '(F8.4)') fb_lambda(i,j)
       write(31, '(F8.4)') phi(i,j)
       write(41, '(F8.4)') a(i,j)
    end do
    close(11)
    close(21)
    close(31)
    close(41)
  end subroutine output_to_file

  subroutine update_rate()
    implicit none
    integer i, j
    real temp
    real coeff
    temp = 0.0
    do i = 1, Lbox
       do j = 1, Lbox
          a(i,j) = alpha*p(i, j)
          if ( cmat(i,j)%type == 0 ) then
             temp = temp + a(i, j)
          end if
       end do
    end do    

    !print *, temp, Amax
    if ( temp > Amax ) then
       !print *, 'come here', temp, Amax
       coeff = Amax / temp
       do i = 1, Lbox
          do j = 1, Lbox
             a(i,j) = coeff*a(i,j)
          end do
       end do
    end if
  end subroutine update_rate

  subroutine update_lambda()
    implicit none
    integer i, j, isub, jsub
    fb_lambda = 0.0
    do i = 1, Lbox
       do j = 1, Lbox
          if ( cmat(i,j)%type == 3 ) then
             do isub = i-1, i+1
                if (isub > 0 .and. isub <= Lbox) then
                   do jsub = j-1, j+1
                      if (jsub > 0 .and. jsub <= Lbox) then
                         fb_lambda(isub,jsub) = lambda_field(isub, jsub)
                      end if
                   end do
                end if
             end do
          end if
       end do
    end do
  end subroutine update_lambda

  subroutine update_phi(march_time)
    implicit none
    real, intent(in) :: march_time
    integer i, j
    real dt_small
    integer istep, nstep
    nstep = 1
    
    dt_small = march_time/nstep
    do istep = 1, nstep
       ! source
       do i = 1, Lbox
          do j = 1, Lbox
             phi(i,j) = phi(i,j) + fb_lambda(i,j)*dt_small
          end do
       end do
       ! Neumann boundary condition (sort of)
       ! i = 1 or Lbox
       do j = 1, Lbox
          phi(1,j) = phi(2,j)
          phi(Lbox,j) = phi(Lbox-1,j)
       end do
       ! j = 1 or Lbox
       do i = 1, Lbox
          phi(i,1) = phi(i,2)
          phi(i,Lbox) = phi(i,Lbox-1)
       end do

       ! diffuse
       phi_old = phi
       do i = 2, Lbox-1
          do j = 2, Lbox-1
             phi(i, j) = phi_old(i,j) &
                  + dt_small*diff*(phi_old(i-1,j)+phi_old(i+1,j)+phi_old(i,j-1)+phi_old(i,j+1)-4.0*phi_old(i,j))
          end do
       end do

       ! decay
       phi = phi - dt_small*gamma*phi
    end do
  end subroutine update_phi

  subroutine update_p()
    implicit none
    integer i, j, isub, jsub
    !p = 1.0 + phi
    p = phi_tiny + phi
  end subroutine update_p

end module setting
