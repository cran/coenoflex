      subroutine totphy(numplt,numspc,numgrd,centrd,spcamp,physio)
c
c* parameters
c
      integer numplt
      integer numspc
      integer numgrd
c
c* species
c
      double precision spcamp(numspc,numgrd,5)
      double precision physio(numspc,numgrd+10)
c
c* plots
c
      double precision centrd(numplt,numgrd)
c
c* local
c
      double precision denom
      double precision sumsq
      double precision sum
c
c* coenoflex/phys ****************** one ***************************
c
      do 10 i=1,numplt
      numout = 0
      sum = 0.0
      sumsq = 0.0
        do 11 j=1,numspc
        denom = 0.0
          do 12 k=1,numgrd
          if (centrd(i,k) .le. spcamp(j,k,1)) then
            physio(j,k) = 0.0
          else if (centrd(i,k) .gt. spcamp(j,k,1) .and.
     +             centrd(i,k) .le. spcamp(j,k,2)) then
            physio(j,k) = 2 * ((centrd(i,k)-spcamp(j,k,1))/
     +             (spcamp(j,k,3)-spcamp(j,k,1)))**2
          else if (centrd(i,k) .gt. spcamp(j,k,2) .and.
     +             centrd(i,k) .le. spcamp(j,k,3)) then
            physio(j,k) = 1 - 2 * ((spcamp(j,k,3)-centrd(i,k))/
     +              (spcamp(j,k,3)-spcamp(j,k,1)))**2
          else if (centrd(i,k) .gt. spcamp(j,k,3) .and.
     +             centrd(i,k) .lt. spcamp(j,k,4)) then
            physio(j,k) = 1 - 2 * ((centrd(i,k)-spcamp(j,k,3))/
     +              (spcamp(j,k,5)-spcamp(j,k,3)))**2
          else if (centrd(i,k) .gt. spcamp(j,k,4) .and.
     +             centrd(i,k) .le. spcamp(j,k,5)) then
            physio(j,k) = 2 * ((spcamp(j,k,5) - centrd(i,k))/
     +               (spcamp(j,k,5)-spcamp(j,k,3)))**2
          else
            physio(j,k) = 0.0
          endif
   12     continue
   11   continue
c
   10 continue
c
      return
c
      end
