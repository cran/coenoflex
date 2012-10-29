      subroutine testdf(x,y,z,test)
c
      integer x
      integer y
      integer z
      integer test(x,y,z)
c
      do i=1,10
      write(6,'(3i3)') test(i,:,:)
      end do    
c
      return
c
      end
