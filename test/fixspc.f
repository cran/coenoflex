      subroutine fixspc(numspc,numgrd,spcamp,maxabu,grdlth,width,
     +                  variab,grdtyp,skew,hiecon,numpts,index)
c
c* parameters
c
      integer numspc
      integer numgrd
c
c* species
c
      double precision spcamp(numspc,numgrd,5)
      double precision maxabu(numspc)
c
c* grads
c
      double precision grdlth(numgrd)
      double precision width(numgrd)
      double precision variab(numgrd)
      integer grdtyp(numgrd)
c
c* passed
c
      double precision skew
      double precision hiecon
c
c* local
c
      double precision size
      double precision expans
      integer numpts(numgrd)
      integer spcpnt
      integer index(numgrd)
      double precision center
      double precision fudge
      double precision hcnadj
c
c***************************** one *********************************
c
      size = 1.0
      do 10 i=1,numgrd
      size = size * grdlth(i)
   10 continue
c
      expans = (numspc / size)**(1.0/numgrd)
c
      numspc = 1
      do 11 i=1,numgrd
      numpts(i) = nint(expans * grdlth(i))
      numspc = numspc * numpts(i)
   11 continue
c
      totsam = 1
      do 13 i=1,numgrd
      index(i) = totsam
      totsam = totsam * numpts(i)
   13 continue
c
      do 14 i=1,numspc
      if (skew .eq. 0) then
        maxabu(i) = 100.0
      else
        maxabu(i) = 0.0
        do 12 j=1,3
        maxabu(i) = maxabu(i) + rand()
   12   continue
        maxabu(i) = (maxabu(i)/3.0)**skew * 100.0
      endif
      hcnadj = 1.0 + ((maxabu(i)/100.0)-0.5) * hiecon

        do 15 j=1,numgrd
        range = grdlth(j) + width(j)
c       center = mod((i-1)/index(j),numpts(j)) *
c    +                   (grdlth(j)/(numpts(j)-1))
        center = mod((i-1)/index(j),numpts(j)) *
     +                (range/(numpts(j)-1)) - width(j)/2
        range = grdlth(j) + width(j)
        if (grdtyp(j) .eq. 1) then
c         spcamp(i,j,3) = center * range - (width(j)/2.0)
          spcamp(i,j,3) = center 
          fudge = (rand() - 0.5) * variab(j)/50.0 * width(j)
          spcamp(i,j,1) = spcamp(i,j,3) - width(j)*hcnadj + fudge
          fudge = (rand() - 0.5) * variab(j)/50.0 * width(j)
          spcamp(i,j,5) = spcamp(i,j,3) + width(j)*hcnadj + fudge
          spcamp(i,j,2) = (spcamp(i,j,1) + spcamp(i,j,3)) / 2.0
          spcamp(i,j,4) = (spcamp(i,j,3) + spcamp(i,j,5)) / 2.0
        else
          spcamp(i,j,2) = center * grdlth(j)
          spcamp(i,j,3) = grdlth(j)
          spcamp(i,j,1) = spcamp(i,j,2) - (spcamp(i,j,3)-spcamp(i,j,2))
          spcamp(i,j,4) = grdlth(j)
          spcamp(i,j,5) = grdlth(j)
          maxabu(i) = min(100.0,maxabu(i) * (1.5 - (1.0-center)))
        endif
   15   continue
   14 continue
c
      return
c
      end
