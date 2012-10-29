      subroutine rndspc(numspc,numgrd,spcamp,maxabu,grdlth,
     +                  alphad,width,variab,grdtyp,
     +                  skew,hiecon)
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
c* common grads
c
      double precision grdlth(numgrd)
      double precision alphad(numgrd)
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
      double precision fudge
      double precision hcnadj
c
c* coenoflex/rndspc *************** one *****************************
c

      do 10 i=1,numspc
      if (skew .eq. 0) then
        maxabu(i) = 100.0
      else
        maxabu(i) = 0.0
        do 11 j=1,3
        maxabu(i) = maxabu(i) + rand()
   11   continue
        maxabu(i) = (maxabu(i)/3.0)**skew * 100.0
      endif
      hcnadj = 1.0 + ((maxabu(i)/100.0)-0.5) * hiecon
        do 12 j=1,numgrd
        range = grdlth(j) + width(j)
        center = rand()**alphad(j)
        if (grdtyp(j) .eq. 1) then
          spcamp(i,j,3) = center * range - (width(j)/2.0)
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
   12   continue
c
   10 continue
c
      return
c
      end
