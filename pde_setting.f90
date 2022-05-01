module pde_setting
  integer Lbox
  integer :: iseed
  real :: tend, dt
  real :: alpha, beta, Amax
  integer :: Rc, nc
  real :: diff, lambda, gamma, phi_tiny
  real :: lambda_phi, lambda_eta, lambda_psi
  real :: k_eta, k_psi  
  integer :: model_type, read_lambda_from_file
  real :: h, h2
  
  namelist /xdata/ Lbox, tend, dt, alpha, beta, Amax, diff, &
       lambda_phi, lambda_eta, lambda_psi, k_eta, k_psi, &
       iseed, tpinc, Rc, nc, lambda, gamma, phi_tiny, model_type, &
       read_lambda_from_file

  real, allocatable :: phi(:,:)
  real, allocatable :: eta(:,:)
  real, allocatable :: psi(:,:)
  real, allocatable :: phi_old(:,:)  
  real, allocatable :: eta_old(:,:)
  real, allocatable :: psi_old(:,:)

contains
  subroutine read_xdata()
    implicit none
    open(8, file="pde_control.txt", status='OLD', recl=80, delim='APOSTROPHE')
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
    write(*, '(a20, f10.4)') 'k_eta = ', k_eta
    write(*, '(a20, f10.4)') 'k_psi = ', k_psi
    write(*, '(a20, f10.4)') 'lambda_psi = ', lambda_psi
    write(*, '(a20, f10.4)') 'lambda_phi = ', lambda_phi
    write(*, '(a20, f10.4)') 'lambda_eta = ', lambda_eta                
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
    write(9, '(a20, f10.2)') 'k_psi,', k_psi
    write(9, '(a20, f10.2)') 'k_eta,', k_eta
    write(9, '(a20, f10.2)') 'lambda_eta,', lambda_eta
    write(9, '(a20, f10.2)') 'lambda_psi,', lambda_psi
    write(9, '(a20, f10.2)') 'lambda_phi,', lambda_phi                
    write(9, '(a20, i10)') 'nc,', nc
    write(9, '(a20, i10)') 'model_type,', model_type
    write(9, '(a20, i10)') 'read_lambda_from_file,', read_lambda_from_file

    close(8)
    close(9)
  end subroutine read_xdata

  subroutine pde_init()
    implicit none
    real u1, u2
    integer i, j, ishift, jshift
    integer curb

    write(*, '(A)', advance='no') 'Initialize...'
    allocate(phi(1:Lbox, 1:Lbox))
    allocate(eta(1:Lbox,1:Lbox))
    allocate(psi(1:Lbox,1:Lbox))
    allocate(phi_old(1:Lbox, 1:Lbox))
    allocate(eta_old(1:Lbox,1:Lbox))
    allocate(psi_old(1:Lbox,1:Lbox))

    h = 1.0 / Lbox
    h2 = h * h
    do i = 1, Lbox
       do j = 1, Lbox
          phi(i,j) = 1.0
          eta(i,j) = 0.0
          psi(i,j) = 0.0
       end do
    end do

    select case (model_type)
    case (1) ! central line
       eta( (Lbox/2 - 5) : (Lbox/2 + 5), (Lbox/2 - 5) : (Lbox/2 + 5) ) = 10.0

    case (2) ! central line

    case (3) ! random

    end select
    write(*, *) 'Done.'
  end subroutine pde_init

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i, j

    WRITE(filename,'(A13,I5.5,A4)') './pde_out/phi', index, '.dat'
    open (unit = 11, file=filename, action="write")
    WRITE(filename,'(A13,I5.5,A4)') './pde_out/eta', index, '.dat'
    open (unit = 31, file=filename, action="write")
    WRITE(filename,'(A13,I5.5,A4)') './pde_out/psi', index, '.dat'
    open (unit = 41, file=filename, action="write")


    do i = 1, Lbox
       do j = 1, Lbox-1
          write(11, '(F8.4, A2)', advance="no") phi(i,j), ', '
          write(31, '(F8.4, A2)', advance="no") eta(i,j), ', '
          write(41, '(F8.4, A2)', advance="no") psi(i,j), ', '
       end do
       write(11, '(F8.4)') phi(i,j)
       write(31, '(F8.4)') eta(i,j)
       write(41, '(F8.4)') psi(i,j)
    end do
    close(11)
    close(31)
    close(41)
  end subroutine output_to_file

  function sigma(x)
    implicit none
    real, intent(in) :: x
    real sigma
    if ( x > 2.0 ) then
       sigma = 1.0
    else
       sigma = 0.0
    end if
    return
  end function sigma
    
  subroutine step()
    implicit none
    !real sigma
    integer i, j
    integer istep, nstep

    ! source and sink
    psi_old = psi
    eta_old = eta
    phi_old = phi
    do i = 1, Lbox
       do j = 1, Lbox
          psi(i,j) = psi_old(i,j) + dt*k_psi*sigma(eta_old(i,j))*phi_old(i,j) - dt*lambda_psi*psi_old(i,j)
          phi(i,j) = phi_old(i,j) - dt*lambda_phi*sigma(eta_old(i,j))*phi_old(i,j)
          eta(i,j) = eta_old(i,j) + dt*k_eta*psi_old(i,j) - dt*lambda_eta*eta_old(i,j)
          !write(*, *) i, j, eta_old(i, j), eta(i, j), psi_old(i, j), phi(i, j)
       end do
    end do
    
    ! diffuse
    ! 1. Neumann boundary condition (sort of)
    ! i = 1 or Lbox
    do j = 1, Lbox
       psi(1,j) = psi(2,j)
       psi(Lbox,j) = psi(Lbox-1,j)
       eta(1,j) = eta(2,j)
       eta(Lbox,j) = eta(Lbox-1,j)
    end do
    ! j = 1 or Lbox
    do i = 1, Lbox
       psi(i,1) = psi(i,2)
       psi(i,Lbox) = psi(i,Lbox-1)
       eta(i,1) = eta(i,2)
       eta(i,Lbox) = eta(i,Lbox-1)
    end do

    ! 2. Laplace operator
    psi_old = psi
    eta_old = eta
    phi_old = phi
    do i = 2, Lbox-1
       do j = 2, Lbox-1
          psi(i, j) = psi_old(i,j) &
               + dt*diff*(psi_old(i-1,j)+psi_old(i+1,j)+psi_old(i,j-1)+psi_old(i,j+1)-4.0*psi_old(i,j))/h2
          !eta(i, j) = eta_old(i,j) &
          !     + dt*diff*(eta_old(i-1,j)+eta_old(i+1,j)+eta_old(i,j-1)+eta_old(i,j+1)-4.0*eta_old(i,j))/h2
       end do
    end do

  end subroutine step

end module pde_setting
