      subroutine auteco(numspc,numgrd,argmnt,grdlst,numper,count,physio)
c
c* parameters
c
      integer numspc
      integer numgrd
c                                                                               
c* arglst                                                                
c                                                                               
      integer argmnt(10)    ! LIST OF COMMAND LINE ARGUMENTS               
      integer grdlst(10,10)                                    
      integer numper(10)                                      
      integer count            ! NUMBER OF ARGUMENTS PASSED                   
c                                                              
c* species
c
      double precision physio(numspc,numgrd+10)
c
c* passed out
c
      integer i
c
      do 10 i=1,count
      if (argmnt(i) .eq. 1) then
        call avephy(numspc,numgrd,numper,physio,grdlst,i)
      else if (argmnt(i) .eq. 2) then
        call minphy(numspc,numgrd,numper,physio,grdlst,i)
      else if (argmnt(i) .eq. 3) then
        call maxphy(numspc,numgrd,numper,physio,grdlst,i)
      else if (argmnt(i) .eq. 4) then
        call geophy(numspc,numgrd,numper,physio,grdlst,i)
      else if (argmnt(i) .eq. 5) then
        call irmphy(numspc,numgrd,numper,physio,grdlst,i)
      endif
   10 continue
c
      return
c
      end
c
c* coenoflex/avephy *************************************************
c
      subroutine avephy(numspc,numgrd,numper,physio,grdlst,stack)
c
c* parameters
c
      integer numspc
      integer numgrd
c
c* arglst                                                                
c                                                                               
      integer grdlst(10,10)                                                    
      integer numper(10)                                                      
c                                                                               
c* species
c
      double precision physio(numspc,numgrd+10)
c
c* passed in
c
      integer stack
c
      do 10 i=1,numspc
      tmp = 0.0
        do 11 j=1,numper(stack)
        tmp = tmp + physio(i,grdlst(stack,j)) 
   11   continue
      physio(i,stack+10) = tmp / numper(stack)
   10 continue
c
      return
c
      end
c
c* coenoflex/minphy ****************************************************
c
      subroutine minphy(numspc,numgrd,numper,physio,grdlst,stack)
c
c* parameters
c
      integer numspc
      integer numgrd
c
c* arglst                                                                
c                                                                               
      integer grdlst(10,10)                                                    
      integer numper(10)                                                      
c                                                                               
c* species
c
      double precision physio(numspc,numgrd+10)
c
c* passed
c
      integer stack
c
      do 10 i=1,numspc
      tmp = 1.0
        do 11 j=1,numper(stack)
        tmp = min(tmp, physio(i,grdlst(stack,j)))
   11   continue
      physio(i,stack+10) = tmp 
   10 continue
c
      return
c
      end
c
c* coenoflex/maxphy *****************************************************
c
      subroutine maxphy(numspc,numgrd,numper,physio,grdlst,stack)
c
c* common params
c
      integer numspc
      integer numgrd
c
c* arglst                                                                
c                                                                               
      integer grdlst(10,10)                                                    
      integer numper(10)                                                      
c                                                                               
c* species
c
      double precision physio(numspc,numgrd+10)
c
c* passed
c
      integer stack
c
      do 10 i=1,numspc
      tmp = 0.0
        do 11 j=1,numper(stack)
        tmp = max(tmp,physio(i,grdlst(stack,j)))
   11   continue
      physio(i,stack+10) = tmp 
   10 continue
c
      return
c
      end
c
c* coenoflex/geophy ***************************************************
c
      subroutine geophy(numspc,numgrd,numper,physio,grdlst,stack)
c
c* params
c
      integer numspc
      integer numgrd
c
c* arglst                                                                
c                                                                               
      integer grdlst(10,10)                                                    
      integer numper(10)                                                      
c                                                                               
c* species
c
      double precision physio(numspc,numgrd+10)
c
c* passed
c
      integer stack
c
      do 10 i=1,numspc
      tmp = 1.0
        do 11 j=1,numper(stack)
        tmp = tmp * physio(i,grdlst(stack,j)) 
   11   continue
      physio(i,stack+10) = tmp**(1.0/numper(stack))
   10 continue
c
      return
c
      end
c
c* coenoflex/irmphy ***************************************************
c
      subroutine irmphy(numspc,numgrd,numper,physio,grdlst,stack)
c
c* params
c
      integer numspc
      integer numgrd
c
c* arglst                                                                
c                                                                               
      integer grdlst(10,10)                                                    
      integer numper(10)                                                      
c                                                                               
c* species
c
      double precision physio(numspc,numgrd+10)
c
c* passed
c
      integer stack
c
      do 10 i=1,numspc
      tmp = 0.0
        do 11 j=1,numper(stack)
        if (physio(i,grdlst(stack,j)) .gt. 0.0) then
          tmp = tmp + (1.0/physio(i,grdlst(stack,j)))
        else
          physio(i,stack+10) = 0.0
          goto 10 
        endif
   11   continue
      physio(i,stack+10) = numper(stack) / tmp
   10 continue
c
      return
c
      end
