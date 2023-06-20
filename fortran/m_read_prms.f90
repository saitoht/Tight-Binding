!!! Read parameters for simulation !!!
! H. Saito 2023/6/17
! sw_spn: consider spin down or not
! sw_shift_ef: consider the shift of Fermi energy or not
! nsite: number of sites
! nele: number of elements
! nwf: number of maximally localized Wannier functions
! natm: number of atoms
! Vc: volume of the cell (Ang^3)
! site: atomic sites
! hop_up/dn: hopping matrix of up/down spin 
! temp: temperature (K)
module m_read_prms
  implicit none
  character(100), protected, save :: foobar
  character(10), allocatable, protected, save :: ele(:), kpiname(:), kpfname(:)
  integer(4), protected, save :: sw_spn, sw_shift_ef, nal
  integer(4), allocatable ,protected, save :: nkpdiv(:), natm(:)
  integer(4), protected, save :: nsite, nele, countf, nsave, nwf
  real(8), allocatable, protected, save :: kpi(:,:), kpf(:,:), site(:,:)
  complex(8), allocatable, protected, save :: hop_up(:,:,:), hop_dn(:,:,:)
  real(8), allocatable, protected, save :: pos(:,:)
  real(8), protected, save :: temp, aunit, plat(3,3), qlat(3,3), Vc, VBZ
  real(8), parameter :: tol = 1.0d-6, ang2Bohr = 1.889726125, pi = 3.14159265358979323

contains
!!! *************************************************************************
  subroutine get_prms()
    !!! get parameters !!!
    implicit none
    character(100) :: sdummy
    integer(4) :: i, j, n, access, logic_key, idummy
    integer(4) :: def_nele
    real(8) :: def_nwf, fdummy

    if ( .not. access( './INPUT', ' ' ) == 0 ) then
       write(*,*) 'Prepare "INPUT" for input parameters!'
       stop
    end if

    !!! default values
    def_nele=1
    logic_key=0

    call get_key('./INPUT', 2, 'nwf', fdummy, nwf, sdummy, logic_key)
    if ( .not. logic_key == 0 ) then
       print *, 'There is no number of Wanner orbitals (nwf) in "INPUT"!'
       stop
    end if
    call get_key('./INPUT', 2, 'nele', fdummy, nele, sdummy, logic_key)
    if ( .not. logic_key == 0 ) nele = def_nele
    call get_key('./INPUT', 2, 'sw_spn', fdummy, sw_spn, sdummy, logic_key)
    if ( .not. logic_key == 0 ) sw_spn = 0
    call get_key('./INPUT', 3, 'mat', fdummy, idummy, foobar, logic_key)
    if ( .not. logic_key == 0 ) then
       print *, 'Footbar of the file (mat) is not given in "INPUT"!'
       stop
    end if

    !!! get the crystal structure information 
    call get_POS('./POSCAR')

    !!! get hopping matrices
    call get_hop('./Hopping.up', hop_up)
    if ( sw_spn == 2 ) then
       call get_hop('./Hopping.dn', hop_dn)
    end if
    
    call get_syml()

    !!! check input parameters
    print *, '***** READ PARAMETERS *****'
    print *, 'Number of Wannier orbitals:', nwf
    print *, 'Unit of translation vectors:', aunit
    print *, 'Translational vectors:'
    do i = 1, 3
       print *, plat(i,:)
    end do
    print *, 'Reciprocal vectors:'
    do i = 1, 3
       print *, qlat(i,:)
    end do
    do i = 1, size(ele)
       print *, 'Element', i, ':', ele(i)
    end do
    print *, 'Number of atoms in the cell:', natm
    print *, 'Number of sites:', nsite
    print *, 'Positions:'
    do i = 1, size(pos(:,1))
       print *, i, pos(i,:)
    end do
  end subroutine get_prms

!!! *************************************************************************
  subroutine get_POS(infle)
    !!! get crystal structure from 'POSCAR' !!!
    implicit none 
    character(*), intent(in) :: infle
    character(100) :: line, tline
    integer(4) :: i, j, n, ios

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
          read(tline, *) plat(n-1,:)
       elseif (n == 5) then
          read(tline, *) ele(:)
       elseif (n == 6) then
          read(tline, *) natm(:)
          nsite = sum(natm(:))
          allocate(pos(1:nsite,1:3))
       elseif (n > 7) then
          read(tline, *) pos(n-7,:)
       end if
       n = n + 1
    end do
    close(11)
    
    Vc = ( (aunit*ang2Bohr)**3.0 ) &
         * dot_product( plat(1,:), ext_prod(plat(2,:), plat(3,:)) )
    qlat(1,:) = ext_prod(plat(2,:),plat(3,:))
    qlat(2,:) = ext_prod(plat(3,:),plat(1,:))
    qlat(3,:) = ext_prod(plat(1,:),plat(2,:))
    qlat(:,:) = (1.0d0/Vc) * qlat(:,:)
    !!! volume of the Brillouin zone
    VBZ = dot_product( qlat(1,:), ext_prod(qlat(2,:), qlat(3,:)) )
  end subroutine get_POS

!!! *************************************************************************
  subroutine ck_CrsySYm()
    !!! check the crystal symmetry !!!


  end subroutine ck_CrsySYm

!!! *************************************************************************
  subroutine get_syml()
    !!! get symmetry line from syml.foobar !!!
    implicit none
    integer(4) :: n, count, ios
    character(150) :: line, tline
    
    n = 1
    open(12, file='./syml.'//trim(foobar), status='old', action='read')
    do
       read(12,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0) exit
       if (tline(1:1) == '#' .or. tline(1:1) == '!' .or. tline(1:1) == '0' ) cycle
       n = n + 1
    end do
    close(12)
    n = n - 1

    allocate(nkpdiv(1:n), kpi(1:n,1:3), kpf(1:n,1:3), kpiname(1:n), kpfname(1:n))

    count = 1
    open(12, file='./syml.'//trim(foobar), status='old', action='read')
    do
       read(12,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0 .or. count == n+1) exit
       if (tline(1:1) == '#' .or. tline(1:1) == '!' .or. tline(1:1) == '0' ) cycle
       read(tline, *) nkpdiv(count), kpi(count,1:3), kpf(count,1:3), kpiname(count), kpfname(count)
       count = count + 1
    end do
    close(12)
  end subroutine get_syml

!!! *************************************************************************
  subroutine get_hop(infle, hop)
    !!! get hopping of onsite and intersite !!!
    implicit none
    character(*), intent(in) :: infle
    character(10000) :: line, tline
    integer(4) :: i, j, k, ios, count, counti, countj
    integer(4), allocatable :: isite(:,:), ipair(:,:)
    real(8), allocatable :: re_hop(:,:,:), im_hop(:,:,:)
    complex(8), allocatable, intent(out) :: hop(:,:,:)

    count = 1
    open(11, file=infle, status='old', action='read')
    do
       read(11,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0) exit
       if (tline(1:1) == '#' .or. tline(1:1) == '!') cycle
       count = count + 1
    end do
    close(11)
    nal = count - 1

    print *, "nall: ", nal
    allocate(isite(1:nal/(nwf*nwf),1:3), ipair(1:nal/(nwf*nwf),1:2))
    allocate(re_hop(1:nwf,1:nwf,1:nal/(nwf*nwf)), im_hop(1:nwf,1:nwf,1:nal/(nwf*nwf)))
    allocate(site(1:nal/(nwf*nwf),1:3))
    allocate(hop(1:nwf,1:nwf,1:nal/(nwf*nwf)))

    count = 0
    counti = 0
    countj = 1
    open(111, file=infle, status='old', action='read')
    do
       read(111,'(A)',iostat=ios) line
       tline = trim(line)
       if (ios<0 .or. count == nal) exit
       if (tline(1:1) == '#' .or. tline(1:1) == '!') cycle
       read(tline, *) isite(countj,1:3), site(countj,1:3), ipair(countj,1:2), &
            re_hop(mod(counti,nwf)+1,mod(count,nwf)+1,countj), &
            im_hop(mod(counti,nwf)+1,mod(count,nwf)+1,countj)
       count = count + 1
       if ( mod(count, nwf) == 0 ) then
          counti = counti + 1
       end if
       if ( mod(count, nwf*nwf) == 0 ) then
          countj = countj + 1
       end if
    end do
    close(111)

    hop(:,:,:) = dcmplx( re_hop(:,:,:), im_hop(:,:,:) )
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
    implicit none
    real(8), intent(in) :: a(3), b(3)
    real(8) :: ext_prod(3)
    
    ext_prod(1) = a(2)*b(3) - b(2)*a(3)
    ext_prod(2) = a(3)*b(1) - b(3)*a(1)
    ext_prod(3) = a(1)*b(2) - b(1)*a(2)
  end function ext_prod
end module m_read_prms
