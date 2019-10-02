    module streamtracer
    use omp_lib
    implicit none
    
    ! Connectivity tracer!
    
    ! ROT: Reason of termination
    ! If  ROT = 0: Still running
    !     ROT = 1: Out of steps
    !     ROT = 2: Out of domain
    !     ROT = 3: In inner boundary
    !     ROT = -1: vmag = 0
    !     ROT = -2: NaN present
    ! Link: Connectivity in MS
    ! If  Link = 1: SW
    !     Link = 2: Closed
    !     Link = 3: North Open
    !     Link = 4: South Open
    !     Link = 5: SW Inc
    !     Link = 6: North Inc
    !     Link = 7: South Inc
    !     Link = 8: Inc Inc
    
    integer :: ns, openmp_enabled = 0, debug = 0
    double precision :: ds, r_IB=1.
    double precision, dimension(3) :: xc
    logical :: inner_boundary=.false.
    
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
    
    subroutine streamline_array(x0, nlines, v, nx, ny, nz, d, dir, ns, xs, vs, ROT, ns_out)
    double precision, dimension(nlines,3), intent(in) :: x0
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: ns, nx, ny, nz, dir, nlines
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
        call streamline(x0_i, v, nx, ny, nz, d, dir, ns, xs_i, vs_i, ROT(i), ns_out(i))
        do j=1,ns_out(i)
            xs(i,j,:) = xs_i(j,:)
            vs(i,j,:) = vs_i(j,:)
        end do
    end do
	!$omp end do
	
    !$omp end parallel
    
    end subroutine streamline_array

    subroutine streamline(x0, v, nx, ny, nz, d, dir, ns_in, xs, vs, ROT, ns_out)
    implicit none
    double precision, dimension(3), intent(in) :: x0, d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: ns_in, nx, ny, nz, dir
    double precision, dimension(ns_in, 3), intent(out) :: xs, vs
    integer, intent(out) :: ROT, ns_out
    double precision, dimension(3) :: xi
    integer :: i

    !$ openmp_enabled = 1
    
    ns = ns_in
    
    ROT = 0
    
    xs = 0
    xs(1, :) = x0
    xi = x0
    
    vs = 1.
    
    call interpolate(xi, v, nx, ny, nz, d, vs(1,:))
    
    do i=2,ns
        
        ROT = check_bounds(xi, nx, ny, nz, d)
        if(ROT.ne.0) exit
        
        call RK4_update(xi, v, nx, ny, nz, d, dir)
        
        xs(i,:) = xi
        
        call interpolate(xi, v, nx, ny, nz, d, vs(i,:))
        
    end do
    
    ns_out = i-1

    if (ROT.eq.0) then
		ROT = 1
		ns_out = ns
	end if
    
    end subroutine streamline
    
    subroutine connectivity_array(x0, nlines, v, nx, ny, nz, d, link, rotation)
    double precision, dimension(nlines,3), intent(in) :: x0
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, dimension(nlines), intent(out) :: link, rotation
    integer, intent(in) :: nx, ny, nz, nlines
    integer, dimension(nlines) :: ROT_f, ROT_r
    double precision, dimension(3) :: x0_i
    double precision :: rotation_i
    integer :: i

    !$ openmp_enabled = 1
    
    !$omp parallel default(firstprivate) shared(v, x0, ROT_f, ROT_r, link, rotation)
    
    if(debug.gt.0) call thread_count('connectivity_array')
    
    !$omp do schedule(dynamic)
    do i=1,nlines
        !DIR$ NOUNROLL
        x0_i = x0(i,:)
        call streamline_end(x0_i, v, nx, ny, nz, d, 1, ROT_f(i), rotation_i)
		rotation(i) = rotation_i
    end do
    !$omp end do
    
    !$omp do schedule(dynamic)
    do i=1,nlines
        !DIR$ NOUNROLL
        x0_i = x0(i,:)
        call streamline_end(x0_i, v, nx, ny, nz, d, -1, ROT_r(i), rotation_i)
		rotation(i) = rotation(i)+rotation_i
    end do
    !$omp end do
    
    !$omp do schedule(static)
    do i=1,nlines
        link(i) = categorise_end_pts(ROT_f(i), ROT_r(i))
    end do
    !$omp end do
    
    !$omp end parallel
    
    end subroutine connectivity_array
 
    subroutine streamline_end(x0, v, nx, ny, nz, d, dir, ROT, rotation)
    implicit none
    double precision, dimension(3), intent(in) :: x0, d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: nx, ny, nz, dir
    integer, intent(out) :: ROT
	double precision, intent(out) :: rotation
    double precision, dimension(3) :: xi, xm, vi, vm
	double precision :: ri, rm
	integer :: i
    
    ROT = 0
	rotation = 0
    
    xi = x0
	xm = x0
	
	call interpolate(xm, v, nx, ny, nz, d, vm)
	vm = vm/vector_mag(vm)

    do i=2,ns
        
		! Check for end of streamline
        ROT = check_bounds(xi, nx, ny, nz, d)
        if(ROT.ne.0) exit
		
		! Update streamline position
        call RK4_update(xi, v, nx, ny, nz, d, dir)
		
		! Calculate rotation
		call interpolate(xi, v, nx, ny, nz, d, vi)
		vi = vi/vector_mag(vi)
		
		rotation = rotation + acos(max(-1., min(1., vector_dot(vm, vi))))
		
		xm = xi
		vm = vi
        
    end do
 
    if (ROT.eq.0) ROT = 1
 
    end subroutine streamline_end
	
	double precision function vector_dot(v1, v2)
    double precision, dimension(3), intent(in) :: v1, v2
	
	vector_dot = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
	
	end function vector_dot
	
	double precision function vector_mag(v)
    double precision, dimension(3), intent(in) :: v
	
	vector_mag = sqrt(vector_dot(v, v))
	
	end function vector_mag	
    
    integer function categorise_end_pts(f, r)
    integer, intent(in) :: f, r
    integer :: link
    
    if(r==2 .and. f==2) then
		link = 1	! Solar Wind
	elseif(r==3 .and. f==3) then
		link = 2	! Closed
	elseif(r==3 .and. f==2) then
		link = 3	! North-Open
	elseif(r==2 .and. f==3) then
		link = 4	! South-Open
	elseif(r==2 .or. f==2) then
		link = 5	! SW-Inc
	elseif(r==3) then
		link = 6	! North-Inc
	elseif(f==3) then
		link = 7	! South-Inc
	else
		link = 8	! Inc-Inc
    end if  
            
    categorise_end_pts = link
    
    end function categorise_end_pts
    
    subroutine RK4_update(xi, v, nx, ny, nz, d, dir)
    double precision, dimension(3), intent(inout) :: xi
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    double precision, dimension(3) :: xu
    integer, intent(in) :: nx, ny, nz, dir
    double precision, dimension(3) :: k1, k2, k3, k4
    integer :: i

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
    
    double precision function check_bounds(xi, nx, ny, nz, d)
    double precision, intent(in), dimension(3) ::xi, d
    integer, intent(in) :: nx, ny, nz
    double precision :: ri
    
    check_bounds = 0
    
    if ( isnan(xi(1)).or.isnan(xi(2)).or.isnan(xi(3)) ) then
        check_bounds = -2
        return
    end if
    
    if(xi(1).lt.0.or.xi(1).gt.d(1)*nx) then
        check_bounds = 2
        return
    end if
    
    if(xi(2).lt.0.or.xi(2).gt.d(2)*ny) then
        check_bounds = 2
        return
    end if
    
    if(xi(3).lt.0.or.xi(3).gt.d(3)*nz) then
        check_bounds = 2
        return
    end if
    
    if(inner_boundary) then
        ri = sqrt((xi(1)-xc(1))**2+(xi(2)-xc(2))**2+(xi(3)-xc(3))**2)
        if(ri.le.r_IB) then
            check_bounds = 3
            return
        end if
    end if
    
    end function check_bounds
    
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