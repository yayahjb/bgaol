      
      program randtester
      implicit none
      integer i
      real*8 v(6)
      integer iseed
      
      data iseed /19191/
      
      
      do 1000 i=1,1000
          call fgbrand(iseed, v, 6)
          write(*,"(6f14.11)") v
 1000 continue
 
      do 2000 i=1001,1000000
         call fgbrand(iseed, v, 6)
 2000 continue
 
      write(*,"(6f14.11)")
 
      do 3000 i=1000001,1000001+1000
          call fgbrand(iseed, v, 6)
          write(*,"(6f14.11)") v
 3000 continue
 
      end
