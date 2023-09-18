
! GHK simulator for truncated normal (possibly singular). Truncation can
! be on an interval for each co-ordinate. Alternatively truncation can be
! specified as a < abs(X) < b, where a and b are positive. This gives a disjoint
! union of two intervals [a , b] and [-b , -a].  
!
! @param n: number of random samples to generate by GHK simulator
! @param d: dimension (d >= 1)
! @param means: mean vector of dimension d (d x 1)
! @param sroot: upper triangular square root of the covariance matrix (d x d)
! @param a: lower bounds (d x 1)
! @param b: upper bounds (d x 1)
! @param imod: Vector of indicators for each co-ordinate (1/0) of whether 
!              absolute-value (modulus) should be applied for truncation of
!              of that co-ordinate. (d x 1)
! @return return value X --> vector (n * d) --> can be coerced into a (n x d) matrix
! @return return value W --> vector (n) --> importance sampling weight of each sample.

subroutine rtmvnghk(n, d, means, sroot, a, b, imod, elen, epos, X, W)

IMPLICIT NONE

integer :: n, d, i, p, i1, j, j1, cj1, k, ind


double precision :: unifrnd, qnormr, pnormr, u, v, prob, Fa, Fb, mu_j, s2
double precision, dimension(n*d), INTENT(INOUT) :: X
double precision, dimension(n), INTENT(INOUT) :: W
double precision, dimension(d)       :: a, b, means, uu
double precision, dimension(2 * (d + 1))       :: cur_list, cur_list_u
double precision, dimension(2 * (d + 1))       :: cprob, iprob

double precision, dimension(d*d)    :: sroot
double precision, dimension(d, d)    :: C
double precision :: A1, B1, A2, B2, low1, high1, low2, high2

integer                               :: num, ncur, isdone
integer, dimension(d)                 :: imod, elen, epos

! initialise R random number generator
call rndstart()

ind = 0
do j=1,d
  do k=1,d
    ind = ind + 1
    C(k, j) = sroot(ind)
  end do
end do

ind = 0
!For all samples
i = 1
do while (i <= n)
  W(i) = 0
  isdone = 0
  cj1 = 0
  ! For all dimensions
  do p = 1,d
    uu(p)=9999999.0d0
    ncur = 0
    if (elen(p) == 0 ) then
      u = unifrnd()
      uu(p) = qnormr(u, 0.0d0, 1.0d0, 1, 0)
      isdone = 1
    else
      do j1 = 1,elen(p)
        cj1 = cj1 + 1
        !print *, p, " =p, cj1= ", cj1
        j = epos(cj1)
        s2 = 0.0d0
        do k = 1,(p-1)
          !print *, j, " ", p, " ", s2, " ", C(j, k), " ", uu(k)
          s2 = s2 + C(j, k) * uu(k)
        end do
        mu_j = means(j) + s2
        !print *, "a/b ", j, " ", p, " ", a(j), " ", b(j), " ", cj1, " ", C(j, p), " ", s2
        if (C(j, p) > 0) then
          A1  = (a(j)- mu_j) /C(j, p)
          B1 = (b(j)- mu_j) /C(j, p)
          A2 = (-b(j)- mu_j) /C(j, p)
          B2 = (-a(j)- mu_j) /C(j, p)
        else
          A1 = (b(j)- mu_j) /C(j, p)
          B1 = (a(j)- mu_j) /C(j, p)
          A2 = (-a(j)- mu_j) /C(j, p)
          B2 = (-b(j)- mu_j) /C(j, p)
        end if
        if(imod(j) == 0) then
          low1 = A1
          high1 = B1
          num = 1
        else
          call getlh(A1, B1, A2, B2, num, low1, high1, low2, high2)
        end if
        !print *, A1, " ", B1, " ", A2, " ", B2
        !print *, j, ", ", p , low1, ", ", high1, ", ", low2, ", ", high2, ", ", num
        call list_update(cur_list, cur_list_u, 2 * (d + 1), ncur, num, low1, high1, low2, high2)
        !do i1 = 1,ncur
          !print *, " CL ", cur_list(2*(i1 - 1) + 1), " - ", cur_list(2*(i1 - 1) + 2)
        !end do
        !print *, p, " ", ncur, " done"
      end do
      if (ncur == 0) then
        !print *, "Empty!!"
        isdone = 0
        !stop 0
	!rwarn("Empty")
        EXIT 
       endif
      !print *, i, ", ", ncur, " done done"
      do  i1 = 1,ncur
        Fa         = pnormr(cur_list(2*(i1 - 1)+ 1), 0.0d0, 1.0d0, 1, 0)
        if(Fa == 0.0d0) then
          Fa       = 1 - pnormr(cur_list(2*(i1 - 1)+ 1), 0.0d0, 1.0d0, 0, 0)
        endif
        Fb         = pnormr(cur_list(2*(i1 - 1)+ 2), 0.0d0, 1.0d0, 1, 0)
        if(Fb == 0.0d0) then
          Fb       = 1 - pnormr(cur_list(2*(i1 - 1)+ 2), 0.0d0, 1.0d0, 0, 0)
        endif
        !print *, i1, " l1 = ", cur_list(2*(i1 - 1)+ 1), " l2 = ", cur_list(2*(i1 - 1)+ 2), " Fa = ", Fa, "Fb = ", Fb
        iprob(i1)       = Fb - Fa
        if (i1 > 1) then
          cprob(i1)      = cprob(i1 - 1) + iprob(i1)
        else
          cprob(i1) = iprob(i1)
        end if
      end do
      if(cprob(ncur) /= cprob(ncur)) then
        !print *, "cprob= ", cprob(ncur)
        !stop 0
	!rexit("Unexpected NA in cprob")
      end if
      W(i) = W(i) + log(cprob(ncur))
      v = unifrnd() * cprob(ncur)
      do  i1 = 1,ncur
        if (v <= cprob(i1)) exit
      end do
      u = unifrnd()

      Fa         = pnormr(cur_list(2*(i1 - 1)+ 1), 0.0d0, 1.0d0, 1, 0)
      if(Fa == 0.0d0) then
        Fa       = 1 - pnormr(cur_list(2*(i1 - 1)+ 1), 0.0d0, 1.0d0, 0, 0)
      endif
      Fb         = pnormr(cur_list(2*(i1 - 1)+ 2), 0.0d0, 1.0d0, 1, 0)
      if(Fb == 0.0d0) then
        Fb       = 1 - pnormr(cur_list(2*(i1 - 1)+ 2), 0.0d0, 1.0d0, 0, 0)
      endif
      prob       = u * (Fb - Fa) + Fa
      uu(p)      = qnormr(prob, 0.0d0, 1.0d0, 1, 0)
      isdone = 1
    end if

    if (uu(p) > 1000.0) then
      !print *, u, " =u prob = ", prob, " uu= ", uu(p), " W = ", (W(i))
      uu(p) = 1000.0
      !stop 0
    else if (uu(p) < -1000.0) then
      !print *, u, " =u prob = ", prob, " uu= ", uu(p), " W = ", (W(i))
      uu(p) = -1000.0
      !stop 0
    else if(uu(p) /= uu(p)) then
      !print *, "NaN, uu = ", uu(p) 
      !stop 0
      !rexit("Unexpected NaN in uu(p)")
    end if

    if(isdone == 1) then
      s2 = 0.0d0
      do k = 1, p
        s2 = s2 + C(p, k) * uu(k)
      end do
      ind = (i-1) * d + p
      X(ind)    = means(p) + s2
    end if
  end do
  if(isdone == 0) then
    !stop 0
    CYCLE
  end if
  W(i) = exp(W(i))
  i = i + 1
end do

! reset R random number generator
call rndend()
end subroutine rtmvnghk

! 
! subroutine to calculate union of two possibly overlapping intervals
!
! if there is overlap result is a single interval (return value num = 1)
! if there is no overlap result is two disjoint intervals (return value num = 2) with lower interval first.
!
! @param A1 lower boundary of 1'st interval
! @param B1 upper boundary of 1'st interval
! @param A2 lower boundary of 2'nd interval
! @param B2 upper boundary of 2'nd interval
! @param num Result is an union of num intervals. num=1/2
! @param low1 lower boundary of 1'st return interval
! @param high1 upper boundary of 1'st return interval
! @param low2 lower boundary of 2'nd return interval
! @param high2 upper boundary of 2'nd return interval
subroutine getlh(A1, B1, A2, B2, num, low1, high1, low2, high2)

IMPLICIT NONE

double precision :: A1, B1, A2, B2, low1, high1, low2, high2
integer :: num
if(A1 > B1) then
  if(A2 > B2) then
    num = 0
  else
    low1 = A2
    high1 = B2
    num = 1
  endif
  return
end if
if(A2 > B2) then
  if(A1 > B1) then
    num = 0
  else
    low1 = A1
    high1 = B1
    num = 1
  endif
  return
end if
if (B1 < A2) then
  low1 =  A1
  high1 = B1
  low2 = A2
  high2 = B2
  num = 2
else if (B2 < A1) then
  low1 =  A2
  high1 = B2
  low2 = A1
  high2 = B1
  num = 2
else
  low1 = min(A1, A2)
  high1 = max(B1, B2)
  num = 1
end if
end subroutine getlh

! 
! subroutine to update the current "union of disjoint intervals" by intersecting with 1 or 2 intervals
!
! @param cur_list current list of boundary points of successive intervals
! @param nmax maximum dimension of list
! @param ncur current dimension of list
! @param num number of intervals whose disjoint union is to be updated. 1 or 2.
! @param low1 lower bound of first interval
! @param high1 upper bound of second interval
! @param low2 lower bound of first interval
! @param high2 upper bound of second interval
subroutine list_update(cur_list, cur_list_u, nmax, ncur, num, low1, high1, low2, high2)

IMPLICIT NONE

integer :: ncur, ncur_u, nmax, num, i, j, i0
double precision, dimension(nmax):: cur_list, cur_list_u
double precision :: low1, high1, low2, high2


if (ncur == 0) then
  cur_list(1) = low1
  cur_list(2) = high1
  ncur = 1
  if(num == 2) then
    cur_list(3) = low2
    cur_list(4) = high2
    ncur = 2
  end if
else
  i0 = 1
  ncur_u = 0
  do i = 1, 2*ncur - 1, 2
    if(low1 > cur_list(i + 1)) then
      CYCLE
    else if(high1 < cur_list(i + 1)) then
      i0 = i
      if(high1 < cur_list(i)) then
        exit
      else
        cur_list_u(2 * ncur_u + 1) = max(low1, cur_list(i))
        cur_list_u(2 * ncur_u + 2) = high1
        ncur_u = ncur_u + 1
        cur_list(i) = high1
      end if
      exit
    else
      cur_list_u(2 * ncur_u + 1) = max(low1, cur_list(i))
      cur_list_u(2 * ncur_u + 2) = cur_list(i+1)
      ncur_u = ncur_u + 1
      i0 = i + 2
      CYCLE
    end if
  end do
  if(num == 2) then
	!print *, "i0= ", i0
    do j = i0, 2*ncur - 1, 2
      if(low2 > cur_list(j + 1)) then
		!print *, j, " 000"
        CYCLE
      else if(high2 < cur_list(j + 1)) then
		!print *, j, " 111 ", high2, " ", cur_list(j), " ", cur_list(j+1)
        if(high2 < cur_list(j)) then
          exit
        else
          cur_list_u(2 * ncur_u + 1) = max(low2, cur_list(j))
          cur_list_u(2 * ncur_u + 2) = high2
          ncur_u = ncur_u + 1
          cur_list(j) = high2
          exit
        end if
      else
		!print *, j, " 222"
        cur_list_u(2 * ncur_u + 1) = max(low2, cur_list(j))
        cur_list_u(2 * ncur_u + 2) = cur_list(j+1)
        ncur_u = ncur_u + 1
        CYCLE
      end if
    end do
  end if
  do i = 1, 2*ncur_u, 1
    cur_list(i) = cur_list_u(i)
  end do
  ncur = ncur_u
end if

end subroutine list_update




