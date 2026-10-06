! General dense GEMM with distinct inputs on every OpenMP thread. Compile twice
! using -DREAL_KIND=4 -Dgemm=sgemm and -DREAL_KIND=8 -Dgemm=dgemm.
program check_gemm_concurrent
  use omp_lib
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  integer, parameter :: rk=REAL_KIND
  integer, parameter :: shapes(3,8)=reshape([ &
    137,24,137, 137,1344,2847, 2847,1344,137, 1344,2847,137, &
    1280,137,1280, 137,1280,1280, 257,137,1344, 17,137,2847 ],[3,8])
  integer :: shape,tid,failures
  interface
    integer(c_int) function openblas_get_parallel() bind(C)
      import c_int
    end function
    integer(c_int) function openblas_get_num_threads() bind(C)
      import c_int
    end function
  end interface
  failures=0
  write(*,'(a,i0,a,i0)') 'General GEMM: precision=',rk,' OpenMP threads=',omp_get_max_threads()
  write(*,'(a,i0,a,i0)') 'BLAS internal parallel mode=',openblas_get_parallel(), &
    ' internal threads=',openblas_get_num_threads()
  if(openblas_get_parallel()/=0.or.openblas_get_num_threads()/=1) &
    error stop 'Expected single-threaded kernels with external-thread locking'
  do shape=1,size(shapes,2)
    write(*,'(a,3(i0,1x))') 'Testing m,n,k=',shapes(:,shape)
    flush(6)
!$omp parallel private(tid) reduction(+:failures)
    tid=omp_get_thread_num()
    failures=failures+test_shape(shapes(1,shape)+mod(tid,3), &
      shapes(2,shape)-mod(tid,5),shapes(3,shape)+tid,tid)
!$omp end parallel
    if(failures/=0) error stop 'General GEMM numerical or guard check failed'
  enddo
  write(*,'(a)') 'PASS: general GEMM references, finite values, transposes, beta, and guards'
contains
  real(rk) function avalue(i,j,thread)
    integer,intent(in)::i,j,thread
    avalue=real(mod(13*i+17*j+3*thread+i*mod(j,7),127)-63,rk)/128._rk
  end function
  real(rk) function bvalue(i,j,thread)
    integer,intent(in)::i,j,thread
    bvalue=real(mod(23*i+19*j+11*thread+j*mod(i,11),131)-65,rk)/128._rk
  end function
  real(rk) function cvalue(i,j,thread)
    integer,intent(in)::i,j,thread
    cvalue=real(mod(i+3*j+thread,17)-8,rk)/16._rk
  end function
  integer function test_shape(m,n,k,thread) result(bad)
    integer,intent(in)::m,n,k,thread
    real(rk),allocatable::a(:,:),b(:,:),c(:,:)
    real(8)::reference,asum,product,error,tolerance,max_error
    real(rk),parameter::alpha=1.25_rk,beta=-0.375_rk
    integer::mode,ar,ac,br,bc,i,j,l,sample,pad
    character::ta,tb
    external gemm
    bad=0
    pad=16+mod(thread,7)
    do mode=0,3
      ta='N';tb='N'
      if(btest(mode,0))ta='T'
      if(btest(mode,1))tb='T'
      ar=m;ac=k;br=k;bc=n
      if(ta=='T')then
        ar=k;ac=m
      endif
      if(tb=='T')then
        br=n;bc=k
      endif
      allocate(a(ar+pad,ac),b(br+pad,bc),c(m+pad,n))
      a=-11111._rk;b=-22222._rk;c=-33333._rk
      do j=1,ac
        do i=1,ar
          if(ta=='N')then
            a(i,j)=avalue(i,j,thread)
          else
            a(i,j)=avalue(j,i,thread)
          endif
        enddo
      enddo
      do j=1,bc
        do i=1,br
          if(tb=='N')then
            b(i,j)=bvalue(i,j,thread)
          else
            b(i,j)=bvalue(j,i,thread)
          endif
        enddo
      enddo
      do j=1,n
        do i=1,m
          c(i,j)=cvalue(i,j,thread)
        enddo
      enddo
      call gemm(ta,tb,m,n,k,alpha,a,ar+pad,b,br+pad,beta,c,m+pad)
      max_error=0.
      do sample=0,99
        i=1+mod(sample*7919,m);j=1+mod(sample*6971,n)
        if(sample==1)then
          i=m;j=n
        endif
        reference=0.;asum=0.
        do l=1,k
          product=real(avalue(i,l,thread),8)*real(bvalue(l,j,thread),8)
          reference=reference+product;asum=asum+abs(product)
        enddo
        reference=real(alpha,8)*reference+real(beta,8)*real(cvalue(i,j,thread),8)
        error=abs(real(c(i,j),8)-reference)
        max_error=max(error,max_error)
        tolerance=64._8*real(epsilon(1._rk),8)*max(1._8,asum)
        if(error>tolerance.or..not.ieee_is_finite(c(i,j)))bad=bad+1
      enddo
      if(.not.all(ieee_is_finite(c(1:m,:))).or.any(c(m+1:,:)/=-33333._rk).or. &
        any(a(ar+1:,:)/=-11111._rk).or.any(b(br+1:,:)/=-22222._rk))bad=bad+1
      if(bad/=0)then
!$omp critical
        write(*,'(a,i0,a,2a,a,es12.4)') 'FAIL thread=',thread,' transpose=',ta,tb, &
          ' max error=',max_error
!$omp end critical
      endif
      deallocate(a,b,c)
    enddo
  end function
end program
