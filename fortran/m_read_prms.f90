!!! Read parameters for simulation !!!
! H. Saito 2023/4/2
! nsite: number of sites
module m_read_prms
  implicit none
  character(100), protected, save :: foobar
  character(10), allocatable, protected, save :: ele(:), kpiname(:), kpfname(:)
  integer(4), protected, save :: sw_spn, sw_shift_ef
  integer(4), allocatable ,protected, save :: nkpdiv(:), nkp(:)
  integer(4), protected, save :: nsite, nnsite, nele, countf, nsave, nwf
  integer(4), allocatable, protected, save :: natm(:)
  real(8), allocatable, protected, save :: kpi(:,:), kpf(:,:), site(:,:)
  complex(8), allocatable, protected, save :: hop_up(:,:,:), hop_dn(:,:,:)
  real(8), allocatable, protected, save :: pos(:,:)
  real(8), protected, save :: temp, aunit, plat(3,3), qlat(3,3), Vc
  real(8), parameter :: tol = 1.0d-6, ang2Bohr = 1.889726125, pi = 3.14159265358979323

contains
!!! *************************************************************************
  subroutine get_prms()
    !!! get parameters !!!
    character(100) :: sdummy
    integer(4) :: i, j, n, access, logic_key, idummy
    real(8) :: def_nwf, fdummy

    if ( .not. access( './INPUT', ' ' ) == 0 ) then
       write(*,*) 'Prepare "INPUT" for input parameters!'
       stop
    end if

    !!! default values
    logic_key=0

    call get_key('./INPUT', 2, 'nwf', fdummy, nwf, sdummy, logic_key)
    if ( .not. logic_key == 0 ) then
       print *, 'There is no number of Wanner orbitals (nwf) in "INPUT"!'
       stop
    end if
    call get_key('./INPUT', 2, 'sw_spn', fdummy, sw_spn, sdummy, logic_key)
    if ( .not. logic_key == 0 ) sw_spn = 0
    call get_key('./INPUT', 3, 'mat', fdummy, idummy, foobar, logic_key)
    if ( .not. logic_key == 0 ) then
       print *, 'Footbar of the file (mat) is not given in "INPUT"!'
       stop
    end if

    call get_POS('./POSCAR', nele, aunit, plat, qlat, ele, natm, nsite)

    call get_hop('./Hopping.up', hop_up)
    if ( sw_spn == 2 ) then
       call get_hop('./Hopping.dn', hop_dn)
    end if

    print *, '***** READ PARAMETERS *****'
    print *, 'Number of Wannier orbitals:', nwf
    print *, 'Unit of translation vectors:', aunit
    print *, 'Translation vectors:'
    do i = 1, 3
       print *, plat(:,i)
    end do
    print *, 'Reciprocal vectors:'
    do i = 1, 3
       print *, qlat(:,i)
    end do
    do i = 1, size(ele)
       print *, 'Element', i, ':', ele(i)
    end do
    print *, 'Number of atoms in the cell:', natm
    print *, 'Number of sites:', nsite
    print *, 'Positions:'
    do i = 1, size(pos(1,:))
       print *, i, pos(:,i)
    end do
  end subroutine get_prms

!!! *************************************************************************
  subroutine get_POS(infle, nele, aunit, plat, qlat, ele, natm, nsite)
    !!! get crystal structure from 'POSCAR' !!!
    implicit none 
    character(*), intent(in) :: infle
    character(*), intent(out), allocatable :: ele(:)
    character(100) :: line, tline
    integer(4) :: i, j, n, ios
    integer(4), intent(in) :: nele
    integer(4), intent(out) :: nsite
    integer(4), intent(out), allocatable :: natm(:)
    real(8), intent(out) :: aunit, plat(3,3), qlat(3,3)

    allocate(ele(nele))
    allocate(natm(nele))
    open(11,file=infle,action='read')
    n = 0
    do
       read(11,'(A)',iostat=ios) line
       tline = trim(line)
!       print *, tline
       if (ios < 0) exit
       if (tline == "") cycle

       if (n == 0) then
          n = n + 1
          cycle
       end if

       if (n == 1) then
          read(tline, *) aunit
       elseif (n > 1 .and. n < 5) then
          read(tline, *) plat(:,n-1)
       elseif (n == 5) then
          read(tline, *) ele(:)
       elseif (n == 6) then
          read(tline, *) natm(:)
          nsite = sum(natm(:))
          allocate(pos(nsite,3))
       elseif (n > 7) then
          read(tline, *) pos(:,n-7)
       end if
       n = n + 1
    end do
    close(11)
    
    Vc = ( (aunit*ang2Bohr)**3.0 ) &
         * dot_product( plat(:,1), ext_prod(plat(:,2), plat(:,3)) )
    qlat(:,1) = ext_prod(plat(:,2),plat(:,3))
    qlat(:,2) = ext_prod(plat(:,3),plat(:,1))
    qlat(:,3) = ext_prod(plat(:,1),plat(:,2))
    qlat(:,:) = (1.0d0/Vc) * qlat(:,:)
  end subroutine get_POS

!!! *************************************************************************
  subroutine get_syml()
    !!! get symmetry line from syml.foobar !!!
    integer(4) :: n, count, ios
    character(100) :: line, tline
    
    n = 1
    open(12, file='./syml.'//trim(foobar), status='old', action='read')
    do
       read(12,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0) exit
       if (tline(1:1) == '#' .or. line(1:1) == '!' .or. line(1:1) == '0' ) cycle
       n = n + 1
    end do
    close(12)

    allocate(nkpdiv(n), kpi(n,3), kpf(n,3), kpiname(n), kpfname(n))

    count = 1
    open(12, file='./syml.'//trim(foobar), status='old', action='read')
    do
       read(12,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0) exit
       if (tline(1:1) == '#' .or. line(1:1) == '!' .or. line(1:1) == '0' ) cycle
       read(tline, *) nkpdiv(count), kpi(count,:), kpf(count,:), kpiname(count), kpfname(count)
       count = count + 1
    end do
    close(12)
  end subroutine get_syml

!!! *************************************************************************
  subroutine get_hop(infle, hop)
    !!! get hopping of onsite and intersite !!!
    character(*), intent(in) :: infle
    character(10000) :: line, tline
    integer(4) :: i, j, k, ios, count
    integer(4) :: isite(1:nwf**4,3), ipair(1:nwf**4,2)
    real(8) :: re_hop(1:nwf**4), im_hop(1:nwf**4), sitev(1:nwf**4,3)
    complex(8) :: hopij(1:nwf**4), hopi(1:nwf*nwf), hopj(1:nwf)
    complex(8), intent(out) :: hop(1:nwf,1:nwf,1:nwf*nwf)

    count = 1
    open(11, file=infle, status='old', action='read')
    do
       read(11,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0) exit
       if (tline(1:1) == '#' .or. line(1:1) == '!') cycle
       read(tline, *) isite(count,:), sitev(count,:), ipair(count,:), &
            re_hop(count), im_hop(count)
       count = count + 1
    end do
    close(11)

    if ( .not. allocated(site) ) then
       allocate(site(nwf*nwf, 3))
    end if
    site(:,:) = sitev(1:nwf**4:nwf*nwf,:)

    hopij(:) = dcmplx( re_hop(:), im_hop(:) )

    do i = 1, nwf*nwf
       hopi(:) = hopij(i*nwf*nwf:(i+1)*nwf*nwf)
       do j = 1, nwf
          hopj(:) = hopi(j*nwf:(j+1)*nwf)
          do k = 1, nwf
             hop(j,k,i) = hopj(k)
          end do
       end do
    end do
  end subroutine get_hop

!!! *************************************************************************
  subroutine get_key(infle, sw, key, fval, ival, string, lgc_key)
    !!! get number or word next to the keyword !!!
    implicit none 
    character(*), intent(in) :: infle, key
    character(10000) :: line, tline
    character(100), intent(out) :: string
    integer(4) :: i, j, ios
    integer(4), intent(in) :: sw
    integer(4), intent(out) :: ival
    integer(4), intent(inout) :: lgc_key
    real(8), intent(out) :: fval
    
    fval = 0.0d0
    ival = 0
    string = ''
    lgc_key = -1

    open(11,file=infle,action='read')
    do
       read(11,'(A)',iostat=ios) line
       tline =trim(line)
       if (ios<0) exit
       if (tline(1:1) == '#' .or. line(1:1) == '!') then
          cycle
       else
          i = index(tline, trim(key))
          if (i==0) then
             cycle
          end if
          j = index(tline((i+len(trim(key))):), "=")
          if ( sw == 1 ) then
             read(tline(i+len(trim(key))+j:),*) fval
          else if ( sw == 2 ) then
             read(tline(i+len(trim(key))+j:),*) ival
          else if ( sw == 3 ) then
             read(tline(i+len(trim(key))+j:),*) string
          end if
          lgc_key = 0
          exit
       end if
    end do
    close(11)    
  end subroutine get_key
  
!!! *************************************************************************
  function ext_prod(a, b)
    !!! return external product !!!
    real(8), intent(in) :: a(3), b(3)
    real(8) :: ext_prod(3)
    
    ext_prod(1) = a(2)*b(3) - b(2)*a(3)
    ext_prod(2) = a(3)*b(1) - b(3)*a(1)
    ext_prod(3) = a(1)*b(2) - b(1)*a(2)
  end function ext_prod
end module m_read_prms
