    module streamtracer
    use omp_lib
    implicit none

    ! Connectivity tracer!

    ! ROT: Reason of termination
    ! If  ROT = 0: Still running
    !     ROT = 1: Out of steps
    !     ROT = 2: Out of domain
    !     ROT = -1: vmag = 0
    !     ROT = -2: NaN present

    integer :: ns, openmp_enabled = 0, debug = 0
    double precision :: ds
    double precision, dimension(3) :: xc

    contains

    subroutine thread_count(message)
    character(len=*), intent(in) :: message
    integer :: num_threads, max_threads

    num_threads = get_num_threads()
    max_threads = get_max_threads()

    !$ openmp_enabled = 1

    !$omp single
    open(unit=500, file='threads.txt', access='append')
    write(500,200) message, num_threads, max_threads
    close(500)
200 format(A, ': ', I4, ', ', I4)
    !$omp end single

    end subroutine thread_count

    subroutine set_num_threads(n_threads)
    integer, intent(in) :: n_threads

    !$ openmp_enabled = 1

    call omp_set_num_threads(n_threads)

    end subroutine set_num_threads

    integer function get_num_threads()

    !$ openmp_enabled = 1

    get_num_threads = omp_get_num_threads()

    end function get_num_threads

    integer function get_max_threads()

    !$ openmp_enabled = 1

    get_max_threads = omp_get_max_threads()

    end function get_max_threads

    subroutine streamline_array(x0, nlines, v, nx, ny, nz, d, dir, ns, cyclic, xs, vs, ROT, ns_out)
    ! INPUT:
    ! x0: seed points
    ! v: vector field
    ! nx, ny, nz: grid spacing in (x, y, z) directions
    ! d:
    ! dir: direction to trace (1=forwards, -1=backwards)
    ! ns: max steps
    !
    ! OUTPUT:
    ! xs: traced streamlines
    ! vs:
    ! ROT: reason of termination
    ! ns_out: number of steps taken
    double precision, dimension(nlines,3), intent(in) :: x0
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: ns, nx, ny, nz, dir, nlines
    integer, intent(in), dimension(3) :: cyclic
    double precision, dimension(nlines, ns, 3), intent(out) :: xs, vs
    integer, intent(out), dimension(nlines) :: ROT, ns_out
    double precision, dimension(3) :: x0_i
    double precision, dimension(ns, 3) :: xs_i, vs_i
    integer :: i, j

    !$ openmp_enabled = 1

    !$omp parallel default(firstprivate) shared(v, xs, vs, x0, ROT, ns_out)

    if(debug.gt.0) call thread_count('streamline_array')

    !$omp do schedule(dynamic)
    do i=1,nlines
        !DIR$ NOUNROLL
        x0_i = x0(i,:)
        call streamline(x0_i, v, nx, ny, nz, d, dir, ns, cyclic, xs_i, vs_i, ROT(i), ns_out(i))
        do j=1,ns_out(i)
            xs(i,j,:) = xs_i(j,:)
            vs(i,j,:) = vs_i(j,:)
        end do
    end do
  !$omp end do

    !$omp end parallel

    end subroutine streamline_array

    subroutine streamline(x0, v, nx, ny, nz, d, dir, ns_in, cyclic, xs, vs, ROT, ns_out)
      ! INPUT:
      ! x0: seed point
      ! v: vector field
      ! nx, ny, nz: grid spacing in (x, y, z) directions
      ! d:
      ! dir: direction to trace (1=forwards, -1=backwards)
      ! ns_in: max steps
      ! cyclic: if true, cyclic boundaries (instead of terminating)
      !
      ! OUTPUT:
      ! xs: traced streamline
      ! vs:
      ! ROT: reason of termination
      ! ns_out: number of steps taken
    implicit none
    double precision, dimension(3), intent(in) :: x0, d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: ns_in, nx, ny, nz, dir
    integer, intent(in), dimension(3) :: cyclic
    double precision, dimension(ns_in, 3), intent(out) :: xs, vs
    integer, intent(out) :: ROT, ns_out
    double precision, dimension(3) :: xi
    integer :: i

    !$ openmp_enabled = 1

    ns = ns_in

    ROT = 0

    ! Set all streamline points to (0, 0, 0) to start
    xs = 0
    ! Set all the output vectors to (1,1,1) to start
    vs = 1.

    ! Set the initial point to be the seed
    xs(1, :) = x0
    xi = x0

    call interpolate(xi, v, nx, ny, nz, d, vs(1 ,:))

    do i=2,ns

        ! Do a single step
        call RK4_update(xi, v, nx, ny, nz, d, dir)
        ! Check if we are out of bounds, and move if cyclic is on
        call check_bounds(xi, nx, ny, nz, d, cyclic, ROT)
        if(ROT.ne.0) exit
        ! Save the step value
        xs(i,:) = xi

        ! Calculate the local B vector
        call interpolate(xi, v, nx, ny, nz, d, vs(i,:))

    end do

    ns_out = i-1

    if (ROT.eq.0) then
    ROT = 1
    ns_out = ns
  end if

    end subroutine streamline


  double precision function vector_dot(v1, v2)
    double precision, dimension(3), intent(in) :: v1, v2

  vector_dot = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

  end function vector_dot

  double precision function vector_mag(v)
    double precision, dimension(3), intent(in) :: v

  vector_mag = sqrt(vector_dot(v, v))

  end function vector_mag

    subroutine RK4_update(xi, v, nx, ny, nz, d, dir)
    double precision, dimension(3), intent(inout) :: xi
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    double precision, dimension(3) :: xu
    integer, intent(in) :: nx, ny, nz, dir
    double precision, dimension(3) :: k1, k2, k3, k4

    !--- RK4 K parameters ---------------------------------------------------------------------
    call stream_function(xi, v, nx, ny, nz, d, dir, k1)

    !DIR$ NOUNROLL
    xu = xi+0.5*k1
    call stream_function(xu, v, nx, ny, nz, d, dir, k2)

    !DIR$ NOUNROLL
    xu = xi+0.5*k2
    call stream_function(xu, v, nx, ny, nz, d, dir, k3)

    !DIR$ NOUNROLL
    xu = xi+k3
    call stream_function(xu, v, nx, ny, nz, d, dir, k4)

    !--- Step ---------------------------------------------------------------------------------

    !DIR$ NOUNROLL
    xi = xi + (k1 + 2*k2 + 2*k3 + k4)/6.

    end subroutine RK4_update

    subroutine check_bounds(xi, nx, ny, nz, d, cyclic, ROT)
    ! INPUT
    ! xi: vector position
    ! nx, ny, nz: number of grid points in (x, y, z)
    ! d: vector of grid spacings in (x, y, z)
    ! cyclic: bool, if true then don't terminate at edge of box
    double precision, intent(in), dimension(3) :: d
    integer, intent(in) :: nx, ny, nz
    integer, intent(in), dimension(3) :: cyclic
    double precision, intent(out), dimension(3) :: xi
    integer, intent(out) :: ROT


    ROT = 0
    if ( isnan(xi(1)).or.isnan(xi(2)).or.isnan(xi(3)) ) then
        ROT = -2
    elseif(xi(1).lt.0.or.xi(1).gt.d(1)*nx) then
        if(cyclic(1).ne.0) then
          xi(1) = MOD(xi(1) + d(1)*nx, d(1)*nx)
        else
          ROT = 2
        end if
    elseif(xi(2).lt.0.or.xi(2).gt.d(2)*ny) then
      if(cyclic(2).ne.0) then
        xi(2) = MOD(xi(2) + d(2)*ny, d(2)*ny)
      else
        ROT = 2
      end if
    elseif(xi(3).lt.0.or.xi(3).gt.d(3)*nz) then
      if(cyclic(3).ne.0) then
        xi(3) = MOD(xi(3) + d(3)*nz, d(3)*nz)
      else
        ROT = 2
      end if
    end if

    end subroutine check_bounds

    subroutine stream_function(xI, v, nx, ny, nz, d, dir, f)
    implicit none
    double precision, dimension(3), intent(in) :: xI
    double precision, dimension(nx, ny, nz, 3), intent(in) :: v
    integer, intent(in) :: nx, ny, nz
    double precision, dimension(3), intent(in) :: d
    integer, intent(in) :: dir
    double precision, dimension(3), intent(out) :: f
    double precision :: vmag
    double precision, dimension(3) :: vI, distI
    integer, dimension(3) :: i0, i1
    double precision, dimension(2,2,2) :: cell

    !DIR$ NOUNROLL
    i0 = floor(xI/d)+1
    i0(1) = min(max(1,i0(1)), nx-1)
    i0(2) = min(max(1,i0(2)), ny-1)
    i0(3) = min(max(1,i0(3)), nz-1)

    i1 = i0+1

    !DIR$ NOUNROLL
    distI = xI/d+1-i0

    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 1)
    call interp_trilinear(distI, cell, vI(1))

    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 2)
    call interp_trilinear(distI, cell, vI(2))

    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 3)
    call interp_trilinear(distI, cell, vI(3))

    !call interp_trilinear2(distI, v(i0(1),i0(2),i0(3),1), v(i0(1),i0(2),i1(3),1), &
    !                              v(i0(1),i1(2),i0(3),1), v(i0(1),i1(2),i1(3),1), &
    !                              v(i1(1),i0(2),i0(3),1), v(i1(1),i0(2),i1(3),1), &
    !                              v(i1(1),i1(2),i0(3),1), v(i1(1),i1(2),i1(3),1), &
    !                              vI(1))
    !
    !call interp_trilinear2(distI, v(i0(1),i0(2),i0(3),2), v(i0(1),i0(2),i1(3),2), &
    !                              v(i0(1),i1(2),i0(3),2), v(i0(1),i1(2),i1(3),2), &
    !                              v(i1(1),i0(2),i0(3),2), v(i1(1),i0(2),i1(3),2), &
    !                              v(i1(1),i1(2),i0(3),2), v(i1(1),i1(2),i1(3),2), &
    !                              vI(2))
    !
    !call interp_trilinear2(distI, v(i0(1),i0(2),i0(3),3), v(i0(1),i0(2),i1(3),3), &
    !                              v(i0(1),i1(2),i0(3),3), v(i0(1),i1(2),i1(3),3), &
    !                              v(i1(1),i0(2),i0(3),3), v(i1(1),i0(2),i1(3),3), &
    !                              v(i1(1),i1(2),i0(3),3), v(i1(1),i1(2),i1(3),3), &
    !                              vI(3))

  vmag  = sqrt(vI(1)**2+vI(2)**2+vI(3)**2)
    !DIR$ NOUNROLL
    f = dir*vI/vmag*ds

    end subroutine stream_function

    subroutine interpolate(xI, v, nx, ny, nz, d, vI)
    implicit none
    double precision, dimension(3), intent(in) :: xI
    double precision, dimension(nx, ny, nz, 3), intent(in) :: v
    integer, intent(in) :: nx, ny, nz
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(3), intent(out) :: vI
    double precision, dimension(3) :: distI
    integer, dimension(3) :: i0, i1
    double precision, dimension(2,2,2) :: cell

    !DIR$ NOUNROLL
    i0 = floor(xI/d)+1
    i0(1) = min(max(1,i0(1)), nx-1)
    i0(2) = min(max(1,i0(2)), ny-1)
    i0(3) = min(max(1,i0(3)), nz-1)

    i1 = i0+1

    !DIR$ NOUNROLL
    distI = xI/d+1-i0

    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 1)
    call interp_trilinear(distI, cell, vI(1))

    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 2)
    call interp_trilinear(distI, cell, vI(2))

    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 3)
    call interp_trilinear(distI, cell, vI(3))

    end subroutine interpolate

    !--- Trilinear interpolation function -------------------------------------------------------------

    subroutine interp_trilinear(xd, f, fI)
    implicit none
    double precision, intent(in), dimension(3) :: xd
    double precision, intent(in), dimension(0:1,0:1,0:1) :: f
    double precision, intent(out) :: fI
    double precision, dimension(0:1,0:1) :: c
    double precision, dimension(3) :: m_xd
    double precision :: c0, c1

    !DIR$ NOUNROLL
    m_xd = 1-xd

    !--- Interpolate over x -----------------------------------------------------------------------

    c(0,0) = f(0,0,0)*m_xd(1) + f(1,0,0)*xd(1)
    c(1,0) = f(0,1,0)*m_xd(1) + f(1,1,0)*xd(1)
    c(0,1) = f(0,0,1)*m_xd(1) + f(1,0,1)*xd(1)
    c(1,1) = f(0,1,1)*m_xd(1) + f(1,1,1)*xd(1)

    !--- Interpolate over y -----------------------------------------------------------------------

    c0 = c(0,0)*m_xd(2) + c(1,0)*xd(2)
    c1 = c(0,1)*m_xd(2) + c(1,1)*xd(2)

    !--- Interpolate over z -----------------------------------------------------------------------

    fI = c0*m_xd(3) + c1*xd(3)

    end subroutine interp_trilinear

    subroutine interp_trilinear2(xd, f000, f001, f010, f011, f100, f101, f110, f111, fI)
    implicit none
    double precision, intent(in), dimension(3) :: xd
    double precision, intent(in) :: f000, f001, f010, f011, f100, f101, f110, f111
    double precision, intent(out) :: fI
    double precision :: c00, c10, c01, c11, c0, c1
    double precision, dimension(3) :: m_xd

    !DIR$ NOUNROLL
    m_xd = 1-xd

    !--- Interpolate over x -----------------------------------------------------------------------

    c00 = f000*m_xd(1) + f100*xd(1)
    c01 = f001*m_xd(1) + f101*xd(1)
    c10 = f010*m_xd(1) + f110*xd(1)
    c11 = f011*m_xd(1) + f111*xd(1)

    !--- Interpolate over y -----------------------------------------------------------------------

    c0 = c00*m_xd(2) + c10*xd(2)
    c1 = c01*m_xd(2) + c11*xd(2)

    !--- Interpolate over z -----------------------------------------------------------------------

    fI = c0*m_xd(3) + c1*xd(3)

    end subroutine interp_trilinear2

    end module streamtracer
