      subroutine getbasename(filename, indices)
            implicit none
            character(*), intent(in) :: filename
            integer, intent(out) :: indices(2)
            integer :: n, temp, slashtemp, extindex, slashindex

            n = len(filename)
            temp = index(filename, '/')
            slashtemp = 0
            slashindex = 0
            do while (temp.ne.0 .and. slashtemp.ne.n)
                slashtemp = slashtemp + temp
                slashindex = slashtemp
                temp = index(filename(slashtemp+1:), '/')
            end do

            extindex = index(filename(slashindex+1:), '.')
            if (extindex.eq.0) then
                extindex = n + 1
            else
                extindex = slashindex + extindex
            end if

            indices(1) = slashindex + 1
            indices(2) = extindex - 1
      end subroutine
