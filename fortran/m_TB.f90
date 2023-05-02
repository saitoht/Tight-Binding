!!! Tight Binding Modeling !!!
!!! H. Saito 2023/
module m_TB
  use m_read_prms, only: nwf, hop_up, hop_dn, sw_spn, sw_shift_ef, &
       nkp, site
  implicit none
  integer(4), protected, save :: nkpath, nkall
  real(8), protected, allocatable, save :: eigu(:,:), eigd(:,:), kline(:,:)
  real(8), protected, save :: chem_pot, chem_shift
  real(8), parameter :: pi=3.14159265358979, tol=1.0d-8
  complex(8), protected, allocatable, save :: hamk_up(:,:,:), hamk_dn(:,:,:)

contains
!!! *************************************************************************
  subroutine mk_kline( ki, kf, ndiv )
    !!! generate k-point line for band calculation !!!
    real(8), intent(in) :: ki(1:3), kf(1:3)
    real(8) :: kd(1:3), kpline(1:ndiv), kpline_vec(1:ndiv,1:3)
    real(8) :: kd_norm, ki_norm, kf_norm
    integer(4), intent(in) :: ndiv
    integer(4) :: i, j

    kpline = 0.0d0
    kpline_vec = 0.0d0

    kd(:) = kf(:) - ki(:)
    ki_norm = 2.0d0 * pi * sqrt( dot_product(ki, ki) )
    kf_norm = 2.0d0 * pi * sqrt( dot_product(kf, kf) )
    kd_norm = 2.0d0 * pi * sqrt( dot_product(kd, kd) )
    kf_norm = ki_norm + kd_norm
    do i = 0, ndiv
       do j = 1, 3
          kpline_vec(i+1,j) = 2.0d0*pi*ki(j) + 2.0d0*pi*((kf(j)-ki(j))*dble(i))/dble(ndiv)
       end do
       kpline(i+1) = ki_norm + ((kf_norm-ki_norm)*dble(i))/dble(ndiv)
    end do
  end subroutine mk_kline

!!! *************************************************************************
  subroutine find_nearest_pair()
    !!! find the nearest neighbor pairs !!!


  end subroutine find_nearest_pair

!!! *************************************************************************
  subroutine init_TBham(kp, hop, hamk)
    !!! generate Tight-Binding Hamiltonian for a specific k-point !!!
    complex(8) :: hamk_comp
    complex(8), intent(in) :: hop(1:nwf, 1:nwf, 1:nwf*nwf)
    complex(8), intent(out) :: hamk(1:nwf, 1:nwf)
    real(8), intent(in) :: kp(1:3)
    integer(4) :: i, j, k, kmd

    do i = 1, nwf
       do j = 1, nwf
          hamk_comp = ( 0.0d0, 0.0d0 )
          do k = 1, nwf*nwf
             hamk_comp = hamk_comp + hop(i,j,k) * ekr(kp(:), site(k,:)) !!! unit***
          end do
          hamk(i,j) = hamk_comp
       end do
    end do
  end subroutine init_TBham

!!! *************************************************************************
  subroutine init_band()
    !!! perform band calculation !!!
    character(4) :: string
    integer(4) :: i, j, k, count, counti, countj
    
    if ( .not. allocated(eigu) ) then
       allocate(eigu(1:nkall, 1:nwf*nwf))
       allocate(eigd(1:nkall, 1:nwf*nwf))
       allocate(hamk_up(1:nkall, 1:nwf, 1:nwf))
       allocate(hamk_dn(1:nkall, 1:nwf, 1:nwf))
    end if

    count = 0
    counti = 0
    print *, '********** Start Band Calculation **********'
    do i = 1, nkpath
       countj = 0
       do j = 1, nkp(i)
          if ( ( counti .ne. 0 ) .and. ( countj == 0 ) ) then
             eigu(count,:) = eigu(count-1,:)
             if ( sw_spn == 2 ) then
                eigd(count,:) = eigd(count-1,:)
             end if
             go to 10
          end if
          call init_TBham(kline(i,j), hop_up, hamk_up)
          call diag_ham(nwf, hamk_up, eigu(count,:))
          if ( sw_spn == 2 ) then ! up and dn
             call init_TBham(kline(i,j), hop_dn, hamk_dn)
             call diag_ham(nwf, hamk_dn, eigd(count,:))
             !!! include SOC mode in the future
          end if
10        count = count + 1
          countj = countj + 1
       end do
       counti = counti + 1
    end do
    print *, '********** Band Calculation Finished **********'

    !!! saving data !!!
    call init_DOS()
    call init_chemical_potential()
    eigu = eigu - chem_pot
    if ( sw_spn == 2 ) eigd = eigd - chem_pot
    
    count = 0
    do i = 1, nkpath
       write(string, '(i3)') i
       do j = 1, 3
          if ( string(1:1) == ' ' ) string(:) = string(2:)
       end do
       open(11, file='TBband_'//trim(string)//'.out', status='replace', action='write')
       if ( sw_spn == 2 ) then
          write(11,*) '### index   kpoint   eigu(k)(eV)   eigd(k)(eV) ###'
       else 
          write(11,*) '### index   kpoint   eigu(k)(eV) ###'
       end if
       do j = 1, nwf
          do k = 1, nkp(i)
             if ( sw_spn == 2 ) then
                write(11,*) j, kline(i,k), eigu(count+k,j), eigd(count+k,j)
             else
                write(11,*) j, kline(i,k), eigu(count+k,j)
             end if
          end do
          write(11,*) ''
       end do
       count = count + nkp(i)
       close(11)
    end do
    call save_DOS()
  end subroutine init_band

!!! *************************************************************************
  subroutine init_chemical_potential()
    !!! calculate chamical potential of the system !!!

    chem_pot = 0.0d0

    if ( sw_shift_ef == 0 ) chem_pot = chem_pot + chem_shift ! chemical potential shift
  end subroutine init_chemical_potential

!!! *************************************************************************
  subroutine init_DOS()
    !!! perform DOS calculation !!!


  end subroutine init_DOS

!!! *************************************************************************
  subroutine save_DOS()


  end subroutine save_DOS

!!! *************************************************************************
  subroutine diag_ham(n, A, ev)
    !!! diagonalization of Hamiltonian !!!
    implicit none    
    integer(4), intent(in) :: n
    complex(8), intent(inout) :: A(1:n,1:n)
    real(8), intent(out) :: ev(1:n)
    integer(4) :: lda, lwork, liwork, lrwork, info
    real(8), allocatable :: rwork(:)
    complex(8), allocatable :: work(:)
    complex(8) :: qw(1:3)
    real(8) :: qrw(1:3)
    integer(4), allocatable :: iwork(:)
    integer(4) :: qiw(1:3)
    character(1) :: jobz, uplo

    jobz='V'
    uplo='U'
    lda=n

    call zheevd(jobz, uplo, n, 0, lda, 0, qw, -1, qrw, -1, qiw, -1, info)
    if ( info .ne. 0 ) then
       write(6, '(A)') "   program stop @zheevd"
       write(6, '(A,i0)') "   info ==> ", info
       stop
    end if

    lrwork = int(qrw(1)) + 1
    lwork = int(dble(qw(1))) + 1
    liwork = qiw(1) + 1
    allocate( work(1:lwork), iwork(1:liwork), rwork(1:lrwork) )
    work(:) = 0.0d0
    iwork(:) = 0
    rwork(:) = 0.0d0

    call zheevd(jobz, uplo, n, A, lda, ev, work, lwork, rwork, lrwork, iwork, liwork, info)
    if ( info .ne. 0 ) then
       write(6, '(A)') "   program stop @zheevd"
       write(6, '(A,i0)') "   info --> ", info
       stop
    end if
    
    deallocate(work, iwork)
    return
  end subroutine diag_ham

!!! *************************************************************************
  subroutine linear_tetrahedron()
    !!! Linear Tetrahedron Method !!!
    implicit none
    


  end subroutine linear_tetrahedron
  
!!! *************************************************************************
  subroutine k_point_sampling()
    !!! Monkhorst & Pack (1976) Method !!!

!    kp(:) = cf(1)*qlat(1,:) + cf(2)*qlat(2,:) + cf(3)*qlat(3,:)

  end subroutine k_point_sampling

!!! *************************************************************************
  function ekr( k, r )
    !!! phase factor !!!
    implicit none
    real(8), intent(in) :: k(3), r(3)
    real(8) :: kdotr ! dot product of k and r
    complex(8) :: ekr ! phase factor (plane wave)

    kdotr = dot_product( k(:), r(:) )
    ekr = dcmplx( cos(kdotr), sin(kdotr) )
  end function ekr
end module m_TB
