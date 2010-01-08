      subroutine fixplt(maxplt,numgrd,grdlth,grdprd,centrd,pltprd,
     +           numpts,index)
c
c* common params
c
      integer numgrd
      integer maxplt
c
c* grads
c
      double precision grdlth(numgrd)
      double precision grdprd(numgrd)
c
c* plots
c
      double precision centrd(maxplt,numgrd)
      double precision pltprd(maxplt)
c
c* scratch
c
      double precision size
      double precision expans
      double precision grdpos
      integer numpts(numgrd)
      integer totsam
      integer index(numgrd)
c
c* coenoflex/fixplt ******************** one ******************************
c
      totsam = 0
      size = 1.0
      do 11 i=1,numgrd
        totsam = totsam + 1
        size = size * grdlth(i)
   11 continue
c
      expans = (maxplt / size)**(1.0/totsam)
c
      numplt = 1
      do 12 i=1,numgrd
      numpts(i) = nint(expans * grdlth(i))
      numplt = numplt * numpts(i)
   12 continue
c
      totsam = 1
      do 13 i=1,numgrd
      index(i) = totsam 
      totsam = totsam * numpts(i)
   13 continue
c
      do 14 i=1,numplt
      pltprd(i) = 1.0
        do 15 j=1,numgrd
          centrd(i,j) = mod((i-1)/index(j),numpts(j)) *
     +                   (grdlth(j)/(numpts(j)-1))
        if (grdprd(j) .ne. 0.0) then
          grdpos = (centrd(i,j)-(grdlth(j)/2.0))/grdlth(j) *
     +               (grdprd(j)/100.0) + 1
          pltprd(i) =  grdpos * pltprd(i)
        endif
   15   continue
   14 continue
c
      maxplt = numplt
c
      return
c
      end
