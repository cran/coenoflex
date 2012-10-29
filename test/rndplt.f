      subroutine rndplt(numplt,numgrd,centrd,grdlth,grdprd,pltprd)
c
      integer numplt
      integer numgrd
c
c* grads
c
      double precision grdlth(numgrd)
      double precision grdprd(numgrd)
c
c* plots
c
      double precision centrd(numplt,numgrd)
      double precision pltprd(numplt)
c
c* scratch
c
      double precision grdpos
c
c* coenoflex/rndplt ************** one ********************************
c
      do 10 i=1,numplt
        do 11 j=1,numgrd
        centrd(i,j) = rand() *grdlth(j)
   11   continue
c
      if (all(grdprd == 1.0)) then
        pltprd(i) = 1.0
      else
        pltprd(i) = 1.0
        do 12 j=1,numgrd
          if (grdprd(j) .ne. 0.0) then
            grdpos = (centrd(i,j)-(grdlth(j)/2.0))/grdlth(j) *
     +               (grdprd(j)/100.0) + 1
            pltprd(i) =  grdpos * pltprd(i)
          endif
   12   continue
      endif
   10 continue
c
      return
c
      end
