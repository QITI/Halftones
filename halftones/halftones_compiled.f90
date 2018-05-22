!
!
! fast_iterator.f90 - Example of fortran module to be called by Python
!
!
! Autor: Pedro Garcia Freitas [sawp@sawp.com.br]
! License: Creative Commons
!      <http://creativecommons.org/licenses/by-nc-nd/2.5/br/>
!
! see: http://www.sawp.com.br
!
! Sep 2011
!

module halftone
  implicit none
  real, dimension(0:9), parameter, private :: shiaufan_factors = &
    (/  0., 0., 0.53333333, 0., 0., 0.06666667, 0.13333333, 0.26666667, 0., 0./)
  real, dimension(0:11), parameter, private :: atkinson_factors = &
    (/ 0.125, 0.125, 0.0, 0.125, 0.125, 0.125, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0 /)
  real, dimension(0:6), parameter, private :: sierra2_factors = &
    (/ 0.25, 0.1875, 0.0625, 0.125, 0.1875, 0.125, 0.0625 /)
  real, dimension(0:11), parameter, private :: sierra3_factors = &
    (/ 0.15625, 0.09375, 0.0625, 0.125, 0.15625, 0.125, 0.0625, &
       0.0, 0.0625, 0.09375, 0.0625, 0.0 /)
  real, dimension(0:6), parameter, private :: burkes_factors = &
    (/ 0.25, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625 /)
  real, dimension(0:11), parameter, private :: stucki_factors = &
    (/ 0.19047619,  0.0952381, 0.04761905,  0.0952381 ,  0.19047619, &
       0.0952381 ,  0.04761905, 0.02380952,  0.04761905,  0.0952381, &
       0.04761905,  0.02380952 /)
  real, dimension(0:11), parameter, private :: jarvis_factors = &
    (/ 0.10417, 0.14583, 0.0625, 0.10417, 0.14583, 0.10417, 0.0625, &
    0.02083, 0.0625, 0.10417, 0.0625, 0.02083 /)
  real, dimension(0:3), parameter, private :: floyd_factors = &
    (/ 0.4374, 0.1875,  0.3125,  0.0625 /)
  real, dimension(0:3), parameter, private :: sierra24A_factors = &
    (/ 0.5, 0.25, 0.25, 0.0 /)
  real, parameter, private :: threshold = 127.5
  contains
    subroutine shiaufan_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, colm1, colm2, colp3, colm3
      real    :: error
      real, dimension(0:11) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 1
      cols = n - 6
      tmp = matrix

      do row = 0, rows - 1
        do col = 3, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:9) = (/ (error * shiaufan_factors(i), i=0, 9) /)
          colp3 = col + 3
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          colm1 = col - 1
          colm2 = col - 2
          colm3 = col - 3
          tmp(row, colp3) =  diff_error(0) + tmp(row, colp3)
          tmp(row, colp2) =  diff_error(1) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(2) + tmp(row, colp1)
          tmp(rowp1, colp3) =  diff_error(3) + tmp(rowp1, colp3)
          tmp(rowp1, colp2) =  diff_error(4) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(5) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(6) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(7) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(8) + tmp(rowp1, colm2)
          tmp(rowp1, colm3) =  diff_error(9) + tmp(rowp1, colm3)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine shiaufan_iterator


    subroutine atkinson_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, rowp2, colm1, colm2
      real    :: error
      real, dimension(0:11) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 2
      cols = n - 4
      tmp = matrix

      do row = 0, rows - 1
        do col = 2, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:11) = (/ (error * atkinson_factors(i), i=0, 11) /)
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          rowp2 = row + 2
          colm1 = col - 1
          colm2 = col - 2
          tmp(row, colp2) =  diff_error(0) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(1) + tmp(row, colp1)
          tmp(rowp1, colp2) =  diff_error(2) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(3) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(4) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(5) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(6) + tmp(rowp1, colm2)
          tmp(rowp2, colp2) =  diff_error(7) + tmp(rowp2, colp2)
          tmp(rowp2, colp1) =  diff_error(8) + tmp(rowp2, colp1)
          tmp(rowp2, col) =  diff_error(9) + tmp(rowp2, col)
          tmp(rowp2, colm1) =  diff_error(10) + tmp(rowp2, colm1)
          tmp(rowp2, colm2) =  diff_error(11) + tmp(rowp2, colm2)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine atkinson_iterator


    subroutine sierra24a_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-2, 0:n-3), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, rowp1, colm1
      real    :: error
      real, dimension(0:3) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 1
      cols = n - 2
      tmp = matrix
      do row = 0, rows - 1
        do col = 1, cols
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          !print *, row, col, n, m
          error = tmp(row, col) - dithered(row, col)
          diff_error(0:3) = (/ (error * sierra24A_factors(i), i=0, 3) /)
          colp1 = col + 1
          rowp1 = row + 1
          colm1 = col - 1
          tmp(row, colp1) =  diff_error(0) + tmp(row, colp1)
          tmp(rowp1, colp1) =  diff_error(1) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(2) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(3) + tmp(rowp1, colm1)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine sierra24a_iterator


    subroutine sierra2_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, colm1, colm2
      real    :: error
      real, dimension(0:6) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 1
      cols = n - 4
      tmp = matrix

      do row = 0, rows - 1
        do col = 2, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:6) = (/ (error * sierra2_factors(i), i=0, 6) /)
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          colm1 = col - 1
          colm2 = col - 2
          tmp(row, colp2) =  diff_error(0) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(1) + tmp(row, colp1)
          tmp(rowp1, colp2) =  diff_error(2) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(3) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(4) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(5) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(6) + tmp(rowp1, colm2)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine sierra2_iterator


    subroutine sierra3_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, rowp2, colm1, colm2
      real    :: error
      real, dimension(0:11) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 2
      cols = n - 4
      tmp = matrix

      do row = 0, rows - 1
        do col = 2, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:11) = (/ (error * sierra3_factors(i), i=0, 11) /)
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          rowp2 = row + 2
          colm1 = col - 1
          colm2 = col - 2
          tmp(row, colp2) =  diff_error(0) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(1) + tmp(row, colp1)
          tmp(rowp1, colp2) =  diff_error(2) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(3) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(4) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(5) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(6) + tmp(rowp1, colm2)
          tmp(rowp2, colp2) =  diff_error(7) + tmp(rowp2, colp2)
          tmp(rowp2, colp1) =  diff_error(8) + tmp(rowp2, colp1)
          tmp(rowp2, col) =  diff_error(9) + tmp(rowp2, col)
          tmp(rowp2, colm1) =  diff_error(10) + tmp(rowp2, colm1)
          tmp(rowp2, colm2) =  diff_error(11) + tmp(rowp2, colm2)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine sierra3_iterator


    subroutine burkes_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, colm1, colm2
      real    :: error
      real, dimension(0:6) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 1
      cols = n - 4
      tmp = matrix

      do row = 0, rows - 1
        do col = 2, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:6) = (/ (error * burkes_factors(i), i=0, 6) /)
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          colm1 = col - 1
          colm2 = col - 2
          tmp(row, colp2) =  diff_error(0) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(1) + tmp(row, colp1)
          tmp(rowp1, colp2) =  diff_error(2) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(3) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(4) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(5) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(6) + tmp(rowp1, colm2)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols+1)
    end subroutine burkes_iterator


    subroutine stucki_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, rowp2, colm1, colm2
      real    :: error
      real, dimension(0:11) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 2
      cols = n - 4
      tmp = matrix

      do row = 0, rows - 1
        do col = 2, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:11) = (/ (error * stucki_factors(i), i=0, 11) /)
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          rowp2 = row + 2
          colm1 = col - 1
          colm2 = col - 2
          tmp(row, colp2) =  diff_error(0) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(1) + tmp(row, colp1)
          tmp(rowp1, colp2) =  diff_error(2) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(3) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(4) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(5) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(6) + tmp(rowp1, colm2)
          tmp(rowp2, colp2) =  diff_error(7) + tmp(rowp2, colp2)
          tmp(rowp2, colp1) =  diff_error(8) + tmp(rowp2, colp1)
          tmp(rowp2, col) =  diff_error(9) + tmp(rowp2, col)
          tmp(rowp2, colm1) =  diff_error(10) + tmp(rowp2, colm1)
          tmp(rowp2, colm2) =  diff_error(11) + tmp(rowp2, colm2)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine stucki_iterator


    subroutine floyd_steinberg_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-2, 0:n-3), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, rowp1, colm1
      real    :: error
      real, dimension(0:3) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 1
      cols = n - 2
      tmp = matrix
      do row = 0, rows - 1
        do col = 1, cols
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          !print *, row, col, n, m
          error = tmp(row, col) - dithered(row, col)
          diff_error(0:3) = (/ (error * floyd_factors(i), i=0, 3) /)
          colp1 = col + 1
          rowp1 = row + 1
          colm1 = col - 1
          tmp(row, colp1) =  diff_error(0) + tmp(row, colp1)
          tmp(rowp1, colp1) =  diff_error(1) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(2) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(3) + tmp(rowp1, colm1)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine floyd_steinberg_iterator


    subroutine jarvis_iterator(matrix, m, n, returnParam)
      integer, intent(in) :: n, m
      real, dimension(0:m-1, 0:n-1), intent(in) :: matrix
      real, dimension(0:m-3, 0:n-5), intent(out) :: returnParam

      integer :: rows
      integer :: cols
      integer :: row, col, i
      integer colp1, colp2, rowp1, rowp2, colm1, colm2
      real    :: error
      real, dimension(0:11) :: diff_error
      integer, dimension(0:m-1, 0:n-1) :: dithered
      real, dimension(0:m-1, 0:n-1) :: tmp

      rows = m - 2
      cols = n - 4
      tmp = matrix

      do row = 0, rows - 1
        do col = 2, cols + 1
          if (tmp(row, col) >= threshold) then
            dithered(row,col) = 255
          else
            dithered(row,col) = 0
          end if
          error = - dithered(row, col) + tmp(row, col)
          diff_error(0:11) = (/ (error * jarvis_factors(i), i=0, 11) /)
          colp2 = col + 2
          colp1 = col + 1
          rowp1 = row + 1
          rowp2 = row + 2
          colm1 = col - 1
          colm2 = col - 2
          tmp(row, colp2) =  diff_error(0) + tmp(row, colp2)
          tmp(row, colp1) =  diff_error(1) + tmp(row, colp1)
          tmp(rowp1, colp2) =  diff_error(2) + tmp(rowp1, colp2)
          tmp(rowp1, colp1) =  diff_error(3) + tmp(rowp1, colp1)
          tmp(rowp1, col) =  diff_error(4) + tmp(rowp1, col)
          tmp(rowp1, colm1) =  diff_error(5) + tmp(rowp1, colm1)
          tmp(rowp1, colm2) =  diff_error(6) + tmp(rowp1, colm2)
          tmp(rowp2, colp2) =  diff_error(7) + tmp(rowp2, colp2)
          tmp(rowp2, colp1) =  diff_error(8) + tmp(rowp2, colp1)
          tmp(rowp2, col) =  diff_error(9) + tmp(rowp2, col)
          tmp(rowp2, colm1) =  diff_error(10) + tmp(rowp2, colm1)
          tmp(rowp2, colm2) =  diff_error(11) + tmp(rowp2, colm2)
        end do
      end do
      returnParam = dithered(0:rows-1, 2:cols + 1)
    end subroutine jarvis_iterator


    subroutine ordered_dithering3x3(quantized, h, l, returnParam)
      implicit none
      integer, intent(in) :: h, l
      integer, dimension(0:h-1, 0:l-1), intent(in) :: quantized
      integer, dimension(0:3*h-1, 0:3*l-1), intent(out) :: returnParam

      integer, dimension(3,3) :: selectedLevel
      integer :: i, j, k, r, xs, xf, ys, yf
      integer level

      k = 0
      r = 0

      do i = 0, h-1
        do j = 0, l-1
          level = quantized(i,j)
          select case (level)
            case (1)
              selectedLevel(1,:) = (/ 0, 0, 0 /)
              selectedLevel(2,:) = (/ 0, 0, 0 /)
              selectedLevel(3,:) = (/ 0, 0, 0 /)
            case (2)
              selectedLevel(1,:) = (/ 0, 0, 0 /)
              selectedLevel(2,:) = (/ 0, 1, 0 /)
              selectedLevel(3,:) = (/ 0, 0, 0 /)
            case(3)
              selectedLevel(1,:) = (/ 0, 0, 0 /)
              selectedLevel(2,:) = (/ 1, 1, 0 /)
              selectedLevel(3,:) = (/ 0, 0, 0 /)
            case (4)
              selectedLevel(1,:) = (/ 0, 0, 0 /)
              selectedLevel(2,:) = (/ 1, 1, 0 /)
              selectedLevel(3,:) = (/ 0, 1, 0 /)
            case (5)
              selectedLevel(1,:) = (/ 0, 0, 0 /)
              selectedLevel(2,:) = (/ 1, 1, 1 /)
              selectedLevel(3,:) = (/ 0, 1, 0 /)
            case(6)
              selectedLevel(1,:) = (/ 0, 0, 1 /)
              selectedLevel(2,:) = (/ 1, 1, 1 /)
              selectedLevel(3,:) = (/ 1, 0, 1 /)
            case (7)
              selectedLevel(1,:) = (/ 0, 0, 1 /)
              selectedLevel(2,:) = (/ 1, 1, 1 /)
              selectedLevel(3,:) = (/ 1, 1, 0 /)
            case (8)
              selectedLevel(1,:) = (/ 1, 0, 1 /)
              selectedLevel(2,:) = (/ 1, 1, 1 /)
              selectedLevel(3,:) = (/ 1, 1, 0 /)
            case (9)
              selectedLevel(1,:) = (/ 1, 0, 1 /)
              selectedLevel(2,:) = (/ 1, 1, 1 /)
              selectedLevel(3,:) = (/ 1, 1, 1 /)
            case default
              selectedLevel(1,:) = (/ 1, 1, 1 /)
              selectedLevel(2,:) = (/ 1, 1, 1 /)
              selectedLevel(3,:) = (/ 1, 1, 1 /)
          end select
          xs = i + k
          xf = xs + 3
          ys = j + r
          yf = ys + 3
          returnParam(xs:xf, ys:yf) = selectedLevel
          r = r + 2
        end do
          r = 0
          k = k + 2
      end do
    end subroutine ordered_dithering3x3


    subroutine ordered_comb4_iterator(masked, h, l, binary)
      implicit none
      integer, intent(in) :: h, l
      integer, dimension(0:h-1, 0:l-1), intent(in) :: masked
      integer, dimension(0:h-1, 0:4*l-1), intent(out) :: binary

      integer :: i, j, r
      integer, dimension(0:15, 4) :: level

      level(0, :) = (/ 0, 0, 0, 0 /)
      level(1, :) = (/ 0, 0, 0, 1 /)
      level(2, :) = (/ 0, 0, 1, 0 /)
      level(3, :) = (/ 0, 0, 1, 1 /)
      level(4, :) = (/ 0, 1, 0, 0 /)
      level(5, :) = (/ 0, 1, 0, 1 /)
      level(6, :) = (/ 0, 1, 1, 0 /)
      level(7, :) = (/ 0, 1, 1, 1 /)
      level(8, :) = (/ 1, 0, 0, 0 /)
      level(9, :) = (/ 1, 0, 0, 1 /)
      level(10, :) = (/ 1, 0, 1, 0 /)
      level(11, :) = (/ 1, 0, 1, 1 /)
      level(12, :) = (/ 1, 1, 0, 0 /)
      level(13, :) = (/ 1, 1, 0, 1 /)
      level(14, :) = (/ 1, 1, 1, 0 /)
      level(15, :) = (/ 1, 1, 1, 1 /)

      do i = 0, h - 1
        r = 0
        do j = 0, l - 1
          binary(i, j+r:j+r+3) = level(masked(i, j), :)
          r = r + 3
        end do
      end do
    end subroutine ordered_comb4_iterator


    subroutine ordered_comb2_iterator(masked, h, l, returnParam)
      implicit none
      integer, intent(in) :: h, l
      integer, dimension(0:h-1, 0:l-1), intent(in) :: masked
      integer, dimension(0:h-1, 0:2*l-1), intent(out) :: returnParam

      integer :: i, j, r
      integer, dimension(0:3, 2) :: level
      integer, dimension(0:h-1, 0:2*l-1) :: binary

      level(0, :) = (/ 0, 0 /)
      level(1, :) = (/ 0, 1 /)
      level(2, :) = (/ 1, 0 /)
      level(3, :) = (/ 1, 1 /)
      do i = 0, h - 1
        r = 0
        do j = 0, l - 1
          binary(i, j+r:j+r+1) = level(masked(i, j), :)
          r = r + 1
        end do
      end do
      returnParam = binary
    end subroutine ordered_comb2_iterator


    subroutine ordered_comb3_iterator(masked, h, l, returnParam)
      implicit none
      integer, intent(in) :: h, l
      integer, dimension(0:h-1, 0:l-1), intent(in) :: masked
      integer, dimension(0:h-1, 0:3*l-1), intent(out) :: returnParam

      integer :: i, j, r
      integer, dimension(0:7, 3) :: level
      integer, dimension(0:h-1, 0:3*l-1) :: binary

      level(0, :) = (/ 0, 0, 0 /)
      level(1, :) = (/ 0, 0, 1 /)
      level(2, :) = (/ 0, 1, 0 /)
      level(3, :) = (/ 0, 1, 1 /)
      level(4, :) = (/ 1, 0, 0 /)
      level(5, :) = (/ 1, 0, 1 /)
      level(6, :) = (/ 1, 1, 0 /)
      level(7, :) = (/ 1, 1, 1 /)
      do i = 0, h - 1
        r = 0
        do j = 0, l - 1
          binary(i, j+r:j+r+2) = level(masked(i, j), :)
          r = r + 2
        end do
      end do
      returnParam = binary
    end subroutine ordered_comb3_iterator
end module halftone


module inverse_halftone
  implicit none

  contains
    subroutine inverse_ordered_comb2_iterator(binary, h, l, returnParam)
      implicit none
      integer, intent(in) :: h, l
      integer, dimension(0:h-1, 0:l-1), intent(in) :: binary
      integer, dimension(0:h-1, 0:l/2-1), intent(out) :: returnParam

      integer :: i, j, m
      integer ::  selectedLevel, x
      integer, dimension(0:3, 0:1) :: level
      integer, dimension(0:h-1, 0:l/2-1) :: restored
      real :: top, base, color, intervalsize

      level(0, :) = (/ 0, 0/)
      level(1, :) = (/ 0, 1 /)
      level(2, :) = (/ 1, 0 /)
      level(3, :) = (/ 1, 1 /)
      intervalsize = 255.0 / 4.0
      do i = 0, h - 1
        do j = 0, l - 1, 2
          do x = 0, 3
            selectedLevel = 0
            do m = 0, 1
              if (binary(i, j + m) == level(x, m)) then
                selectedLevel = selectedLevel + 1
              end if
            end do
            if (selectedLevel == 2) then
              selectedLevel = x
              exit
            end if
          end do
          base = selectedLevel * intervalsize
          top = base + intervalsize
          color = base + (top - base) * rand()
          restored(i, j / 2) = int(color)
        end do
      end do
      returnParam = restored
    end subroutine inverse_ordered_comb2_iterator


    subroutine inverse_ordered_comb3_iterator(binary, h, l, returnParam)
      implicit none
      integer, intent(in) :: h, l
      integer, dimension(0:h-1, 0:l-1), intent(in) :: binary
      integer, dimension(0:h-1, 0:l/3-1), intent(out) :: returnParam

      integer :: i, j, m
      integer ::  selectedLevel, x
      integer, dimension(0:7, 0:2) :: level
      integer, dimension(0:h-1, 0:l/3-1) :: restored
      real :: top, base, color, intervalsize

      level(0, :) = (/ 0, 0, 0 /)
      level(1, :) = (/ 0, 0, 1 /)
      level(2, :) = (/ 0, 1, 0 /)
      level(3, :) = (/ 0, 1, 1 /)
      level(4, :) = (/ 1, 0, 0 /)
      level(5, :) = (/ 1, 0, 1 /)
      level(6, :) = (/ 1, 1, 0 /)
      level(7, :) = (/ 1, 1, 1 /)
      intervalsize = 255.0 / 8.0
      do i = 0, h - 1
        do j = 0, l - 1, 3
          do x = 0, 7
            selectedLevel = 0
            do m = 0, 2
              if (binary(i, j + m) == level(x, m)) then
                selectedLevel = selectedLevel + 1
              end if
            end do
            if (selectedLevel == 3) then
              selectedLevel = x
              exit
            end if
          end do
          base = selectedLevel * intervalsize
          top = base + intervalsize
          color = base + (top - base) * rand()
          restored(i, j / 3) = int(color)
        end do
      end do
      returnParam = restored
    end subroutine inverse_ordered_comb3_iterator


    subroutine inverse_ordered_dithering_iterator(bin, h, l, msize, returnp)
      implicit none
      integer, intent(in) :: h, l, msize
      integer, dimension(0:h-1, 0:l-1), intent(in) :: bin
      integer, dimension(0:h/msize, 0:l/msize), intent(out) :: returnp

      integer :: i=0, j=0, stepx=0, stepy=0
      integer :: chosedLevel
      integer, dimension(0:h+msize-1, 0:l+msize-1) :: tmp
      integer, dimension(msize, msize) :: mask
      real :: intervalSize, base, top, greyValue

      intervalsize = 255.0 / (msize * msize)
      tmp(0:h-1, 0:l-1) = bin
      do i = 0, msize - 1
        tmp(h + i, :) = 0
        tmp(:, l+i) = 0
      end do

      do i = 0, h/msize-1
        do j = 0, l/msize-1
          mask(:,:) = tmp(stepx+i:stepx+i+msize, stepy+j:stepy+j+msize)
          chosedLevel = sum(mask) - 1
          base = chosedLevel * intervalsize
          top = base + intervalsize
          greyValue = base + (top - base) * rand()
          returnp(i,j) = int(greyValue)
          stepy = stepy + msize - 1
        end do
        stepy = 0
        stepx = stepx + msize - 1
      end do
    end subroutine inverse_ordered_dithering_iterator
end module inverse_halftone


module fbih
  implicit none
  real, parameter, private :: DEBLUR_FACTOR = 1.5d0
  real, parameter, private :: NORM = 100000000
  real, parameter, private :: THRESHOLD_DIFF = 0
  real, parameter, private :: ERRTO = 0.25
  integer, parameter, private :: TRUE = 1
  integer, parameter, private :: FALSE = 0

  contains
    subroutine inverse_halftone(im, rows, cols, returnParam)
      implicit none
      integer, intent(in) :: rows, cols
      integer, dimension(0:rows-1, 0:cols-1), intent(in) :: im
      integer, dimension(0:rows-1, 0:cols-1), intent(out) :: returnParam

      integer :: i, j
      real, dimension(0:rows-1, 0:cols-1) :: inputImage
      real, dimension(0:rows-1, 0:cols-1) :: gaussianLvlOne
      real, dimension(0:rows-1, 0:cols-1) :: median
      real, dimension(0:rows-1, 0:cols-1) :: gaussianLvlTwo
      real, dimension(0:rows-1, 0:cols-1) :: gaussianLvlThree
      real, dimension(0:rows-1, 0:cols-1) :: quantiz

      do i = 0, rows-1
        do j = 0, cols - 1
          inputImage(i, j) = float(im(i, j))
        end do
      end do
      call GaussianFilter1(rows, cols, inputImage, gaussianLvlOne)
      call medianFilter3by3(rows, cols, gaussianLvlOne, median)
      call GaussianFilter2(rows, cols, median, gaussianLvlTwo)
      call GaussianFilter3(rows, cols, gaussianLvlTwo, gaussianLvlThree)
      call thresholdDiff(rows, cols, gaussianLvlTwo, gaussianLvlThree, quantiz)
      call normalise(rows, cols, quantiz, median, returnParam)
    end subroutine inverse_halftone


    subroutine medianFilter3by3(rows, cols, img, output)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: img
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: output

      integer :: row, col
      integer :: lastRow, lastCol, nextRow, nextCol
      real :: block(0:9)

      do row = 0, rows-1
        do col = 0, cols-1
          lastRow = abs(row - 1)
          lastCol = abs(col - 1)
          nextRow = abs(row + 1)
          nextCol = abs(col + 1)
          block(1) = img(lastRow, lastCol)
          block(2) = img(lastRow, col)
          if (nextCol >= cols) then
            block(3) = img(lastRow, 2 * cols - nextCol - 2)
            block(6) = img(row, 2 * cols - nextCol - 2)
          else
            block(3) = img(lastRow, nextCol)
            block(6) = img(row, nextCol)
          end if
          block(4) = img(row, lastCol)
          block(5) = img(row, col)
          if (nextRow >= rows) then
            block(7) = img(2 * rows - nextRow - 2, lastCol)
            block(8) = img(2 * rows - nextRow - 2, col)
          else
            block(7) = img(nextRow, nextCol)
            block(8) = img(nextRow, col)
          end if
          if (nextRow >= rows .and. nextCol >= cols) then
            block(9) = img(2 * rows - nextRow - 2, 2 * cols - nextCol - 2)
          end if
          if (nextRow >= rows .and. nextCol < cols) then
            block(9) = img(2 * rows - nextRow - 2, nextCol)
          end if
          if (nextRow < rows .and. nextCol >= cols) then
            block(9) = img(nextRow, 2 * cols - nextCol - 2)
          end if
          if (nextRow < rows .and. nextCol < cols) then
            block(9) = img(nextRow, nextCol)
          end if
          output(row, col) = getPixel(5, 9, block)
        end do
      end do
    end subroutine medianFilter3by3


    subroutine SWAP(a, b)
      implicit none
      real, intent(inout) :: a, b
      real :: temp
      temp = a
      a = b
      b = temp
    end subroutine SWAP


    real function getPixel(pos, blocksize, blk)
      implicit none
      integer, intent(in) :: pos
      integer, intent(in) :: blocksize
      real, dimension(blocksize), intent(in) :: blk

      integer :: i, j, l, mid, ir
      real :: temp, selected
      real, dimension(blocksize) :: block
      integer, save :: s = 0

      s = s + 1

      i = 0
      j = 0
      mid = 0
      temp = 0.0
      l = 1
      ir = blocksize
      block(:) = blk(:)
      do
        if (ir <= l + 1) then
          if (ir == l + 1 .and. block(ir) < block(l)) then
            call SWAP(block(ir), block(l))
          end if
          exit
        else
          i = l + 1
          j = ir
          mid = int((l + ir) / 2)
          call SWAP(block(l+1), block(mid))
          if (block(l) > block(ir)) then
            call SWAP(block(ir), block(l))
          end if
          if (block(l+1) > block(ir)) then
            call SWAP(block(ir), block(l+1))
          end if
          if (block(l) > block(l+1)) then
            call SWAP(block(l), block(l+1))
          end if
          selected = block(l+1)
          do
            do
              i = i + 1
              if (block(i) >= selected) exit
            end do
            do
              j = j - 1
              if (block(j) <= selected) exit
            end do
            if (j < i) exit
            call SWAP(block(j), block(i))
          end do
          block(l+1) = block(j)
          block(j) = selected
          if (j >= pos) then
            ir = j - 1
          end if
          if (j <= pos) then
            l = i
          end if
        end if
      end do
      getPixel = block(pos)
      return
    end function getPixel


    subroutine GaussianFilter1(m, n, s, output)
      implicit none
      integer, intent(in) :: m, n
      real, dimension(0:m-1, 0:n-1), intent(in) :: s
      real, dimension(0:m-1, 0:n-1), intent(out) :: output
      integer, dimension(9) :: filter = & 
        (/ 11, 135, 808, 2359, 3372, 2359, 808, 135, 11 /)
      call separable9x9FIRBinaryImage(m, n, s, filter, output)
    end subroutine GaussianFilter1


    subroutine separable9x9FIRBinaryImage(rows, cols, source, filter, dest)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: source
      integer, dimension(9), intent(in) :: filter
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: dest

      real, dimension(0:rows-1, 0:cols-1) :: ws

      call rowConvolutions9by9(rows, cols, source, filter, ws)
      call columnConvolutions9by9(rows, cols, ws, filter, dest)
    end subroutine separable9x9FIRBinaryImage


    subroutine rowConvolutions9by9(rows, cols, src, filter, ws)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: src
      integer, dimension(9), intent(in) :: filter
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: ws

      integer :: i, j, k
      integer, dimension(9) :: p
      real :: sums
      real, dimension(cols) :: sourceRow

      do i = 0, rows-1
        sourceRow(:) = src(i, :)
        do j = 0, cols-1
          sums = 0.0
          p(1) = j - 4
          p(2) = j - 3
          p(3) = j - 2
          p(4) = j - 1
          p(5) = j
          p(6) = j + 1
          p(7) = j + 2
          p(8) = j + 3
          p(9) = j + 4
          p(1) = abs(p(1))
          p(2) = abs(p(2))
          p(3) = abs(p(3))
          p(4) = abs(p(4))
          p(6) = reflect(p(6), cols)
          p(7) = reflect(p(7), cols)
          p(8) = reflect(p(8), cols)
          p(9) = reflect(p(9), cols)
          do k = 1, 9
            if (sourceRow(p(k)) /= 0) then
              sums = sums + filter(k)
            end if
          end do
          ws(i, j) = 255 * sums
        end do
      end do
    end subroutine rowConvolutions9by9


    subroutine columnConvolutions9by9(rows, cols, src, filter, dest)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: src
      integer, dimension(9), intent(in) :: filter
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: dest

      integer :: i, j
      integer, dimension(9) :: p
      real :: sums

      p = 0
      do j = 0, cols-1
        do i = 0, rows-1
          sums = 0.0
          p(1) = i - 4
          p(2) = i - 3
          p(3) = i - 2
          p(4) = i - 1
          p(5) = i
          p(6) = i + 1
          p(7) = i + 2
          p(8) = i + 3
          p(9) = i + 4
          p(1) = abs(p(1))
          p(2) = abs(p(2))
          p(3) = abs(p(3))
          p(4) = abs(p(4))
          p(6) = reflect(p(6), rows)
          p(7) = reflect(p(7), rows)
          p(8) = reflect(p(8), rows)
          p(9) = reflect(p(9), rows)
          sums = filter(5) * src(i, j)
          sums = sums + filter(1) * (src(p(1), j) + src(p(9), j))
          sums = sums + filter(2) * (src(p(2), j) + src(p(8), j))
          sums = sums + filter(3) * (src(p(3), j) + src(p(7), j))
          sums = sums + filter(4) * (src(p(4), j) + src(p(6), j))
          dest(i, j) = sums / NORM + ERRTO
        end do
      end do
    end subroutine columnConvolutions9by9


    subroutine GaussianFilter2(m, n, s, dest)
      implicit none
      integer, intent(in) :: m, n
      real, dimension(0:m-1, 0:n-1), intent(in) :: s
      real, dimension(0:m-1, 0:n-1), intent(out) :: dest
      integer, dimension(7) :: filter = & 
        (/ 44, 540, 2420, 3991, 2420, 540, 44 /)
      call separable7x7FIRGreyImage(m, n, s, filter, dest)
    end subroutine GaussianFilter2


    subroutine GaussianFilter3(m, n, s, dest)
      implicit none
      integer, intent(in) :: m, n
      real, dimension(0:m-1, 0:n-1), intent(in) :: s
      real, dimension(0:m-1, 0:n-1), intent(out) :: dest
      integer, dimension(7) :: filter = & 
        (/ 1, 103, 2075, 5641, 2075, 103, 1 /)
      call separable7x7FIRGreyImage(m, n, s, filter, dest)
    end subroutine GaussianFilter3


    subroutine separable7x7FIRGreyImage(rows, cols, source, filter, dest)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: source
      integer, dimension(7), intent(in) :: filter
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: dest

      integer :: i, j
      real, dimension(0:rows-1, 0:cols-1) :: ws
      real, dimension(0:cols-1) :: sourceRow
      integer :: p2, p3, p4, p6, p7, p8

      do i = 0, rows-1
        sourceRow = source(i,:)
        do j = 0, cols-1
          p2 = j - 3
          p3 = j - 2
          p4 = j - 1
          p6 = j + 1
          p7 = j + 2
          p8 = j + 3
          p2 = abs(p2)
          p3 = abs(p3)
          p4 = abs(p4)
          p6 = reflect(p6, cols)
          p7 = reflect(p7, cols)
          p8 = reflect(p8, cols)
          ws(i, j) = filter(4) * sourceRow(j)
          ws (i, j) = ws(i, j) + filter(1) * (sourceRow(p2) + sourceRow(p8))
          ws (i, j) = ws(i, j) + filter(2) * (sourceRow(p3) + sourceRow(p7))
          ws (i, j) = ws(i, j) + filter(3) * (sourceRow(p4) + sourceRow(p6))
        end do
      end do
      call columnConvolutions7by7(rows, cols, ws, filter, dest)
    end subroutine separable7x7FIRGreyImage


    subroutine columnConvolutions7by7(rows, cols, src, filter, dest)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: src
      integer, dimension(7), intent(in) :: filter
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: dest

      integer :: i, j
      integer :: p2, p3, p4, p6, p7, p8
      real :: sums

      do j = 0, cols-1
        do i = 0, rows-1
          p2 = i - 3
          p3 = i - 2
          p4 = i - 1
          p6 = i + 1
          p7 = i + 2
          p8 = i + 3
          p2 = abs(p2)
          p3 = abs(p3)
          p4 = abs(p4)
          p6 = reflect(p6, rows)
          p7 = reflect(p7, rows)
          p8 = reflect(p8, rows)
          sums = filter(4) * src(i,j)
          sums = sums + filter(1) * (src(p2, j) + src(p8, j))
          sums = sums + filter(2) * (src(p3, j) + src(p7, j))
          sums = sums + filter(3) * (src(p4, j) + src(p6, j))
          dest(i, j) = real(int(sums / NORM + ERRTO))
        end do
      end do
    end subroutine columnConvolutions7by7


    integer function reflect(x, y)
      implicit none
      integer, intent(in) :: x, y
      if (x >= y) then
        reflect = 2 * y - x - 2
      else
        reflect = x
      end if
    end function


    subroutine thresholdDiff(rows, cols, xi, zi, hie)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: xi, zi
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: hie

      real, dimension(0:rows-1, 0:cols-1) :: x, z
      real, dimension(0:rows-1, 0:cols-1) :: edgeMap
      real :: pixel
      integer :: i, j

      x(:,:) = xi(:,:)
      z(:,:) = zi(:,:)
      call binaryMedialFilter5by5(rows, cols, z, edgeMap)

      do i = 0, rows-1
        do j = 0, cols-1
          pixel = x(i, j) - z(i, j)
          if (pixel <= THRESHOLD_DIFF .and. pixel >= -THRESHOLD_DIFF) then
            z(i, j) = 0.0
          else
            z(i, j) = 1.0
          end if
          hie(i,j) = pixel
        end do
      end do
      do i = 0, rows-1
        do j = 0, cols-1
          hie(i, j) = hie(i, j) * z(i, j) * edgeMap(i, j)
        end do
      end do
    end subroutine thresholdDiff


    subroutine binaryMedialFilter5by5(rows, cols, img, filtered)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: img
      real, dimension(0:rows-1, 0:cols-1), intent(out) :: filtered

      integer :: row, col
      integer :: ones
      integer, parameter :: threshold = 12

      do row = 0, rows-1
        do col = 0, cols-1
          ones = bitBlockCounter(rows, cols, row, col, img)
          if (ones > threshold) then
            filtered(row, col) = 1.0
          else
            filtered(row, col) = 0.0
          end if
        end do
      end do
    end subroutine binaryMedialFilter5by5


    integer function bitBlockCounter(rows, cols, row, col, img)
      implicit none
      integer, intent(in) :: row, col, rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: img

      integer :: countOnes, blockX, blockY, r, c
      integer, parameter :: blocksize = 4
      integer, parameter :: threshold = 12

      countOnes = 0
      do blockX = 0, blocksize
        r = row - blockX
        do blockY = 0, blocksize
          c = col - blockY
          if (r >= 0 .and. c >= 0) then
            if (img(r, c) > 0) then
              countOnes = countOnes + 1
            end if
            if (countOnes > threshold) then
              exit
            end if
          end if
        end do
        if (countOnes > threshold) then
          exit
        end if
      end do
      bitBlockCounter = countOnes
    end function bitBlockCounter


    subroutine normalise(rows, cols, hie, y1, outputImage)
      implicit none
      integer, intent(in) :: rows, cols
      real, dimension(0:rows-1, 0:cols-1), intent(in) :: hie, y1
      integer, dimension(0:rows-1, 0:cols-1), intent(out) :: outputImage

      integer :: i, j
      real :: outputValue

      do i = 0, rows-1
        do j = 0, cols-1
          outputValue = ERRTO
          outputValue = outputValue + DEBLUR_FACTOR * (hie(i, j) + y1(i, j))
          outputImage(i, j) = int(outputValue)
        end do
      end do
    end subroutine normalise
end module fbih
