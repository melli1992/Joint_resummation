cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stringmod(input,length)
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc   Takes input (string), left-aligns it
ccc   keeping track of the length and adding
ccc!   a character (endchar) at the end.
ccc
cccccccccccccccccccccccccccccccccccccccccccccccc
      character input*(*)
      integer length
cccccccccccccccccccccccccccccccccccccccccccccccc

      input  = adjustl(input)
      length = len(input)

      do while( input(length:length) .eq. '')
         length = length - 1
      enddo

      return
      end
