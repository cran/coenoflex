      subroutine coenoflex(numplt,numspc,numgrd,spcamp,physio,numabu,
     +           grdlth,grdprd,alphad,width,variab,grdsam,grdtyp,
     +           grdnam,centrd,pltprd,
     +           arglst,grdlst,numper)

c* common species
c
      double precision spcamp(numspc,numgrd,5)
      double precision physio(numspc,numgrd+10)
      double precision numabu(numspc)
c
c* common grads
c
      double precision grdlth(numgrd)
      double precision grdprd(numgrd)
      double precision alphad(numgrd)
      double precision width(numgrd)
      double precision variab(numgrd)
      integer grdsam(numgrd)
      integer grdtyp(numgrd)
c
c* common plots
c
      double precision centrd(numplt,numgrd)
      double precision pltprd(numplt)
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      integer srorf
      integer prorf
c
c* common arglst
c
      character*20 argmnt(10)    ! LIST OF COMMAND LINE ARGUMENTS
      integer grdlst(10,10)
      integer numper(10)
      integer count            ! NUMBER OF ARGUMENTS PASSED
c
c* passed
c
c     integer numxgd
      double precision cmpasy
      double precision cmpphy
      double precision maxtot
      double precision noise
      double precision slack
      double precision skew
      double precision hiecon
      character*255 line
      integer final
      integer seed
c
c* local
c
      logical*1 spcfil
      character*1 yorn
      character*20 respon
      character*20 prefix
      integer lenprf
c
c*********************** one ************************************
c
      call spclbl
      call pltlbl
      call totphy(cmpasy,cmpphy,maxtot,noise,slack,final)
c
      end
c
c* coenoflex **************** subroutine rndspc *******************
c
      subroutine rndspc(numspc,numgrd,spcamp,physio,maxabu,gradlth,
     +                  grdprd,alphad,width,variab,grdsam,grdtyp,
     +                  skew,hiecon)
c
      parameter (numspc=500)
      parameter (maxgrd=10)
c
c* common species
c
      double precision spcamp(numspc,maxgrd,5)
      double precision physio(numspc,maxgrd+10)
      double precision maxabu(numspc)
c
c* common grads
c
      double precision grdlth(maxgrd)
      double precision grdprd(maxgrd)
      double precision alphad(maxgrd)
      double precision width(maxgrd)
      double precision variab(maxgrd)
      integer grdsam(maxgrd)
      character*1 grdtyp(maxgrd)
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
c* passed
c
      double precision skew
      double precision hiecon
      integer seed
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
        if (grdtyp(j) .eq. 'e' .or. grdtyp(j) .eq. 'e') then
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
      write(2,100) i,maxabu(i)
        do 13 j=1,numgrd
        write(2,110) j,(spcamp(i,j,k),k=1,5)
   13   continue
   10 continue
c
      close (2)
c
      return
c
  100 format (i4,2x,f10.0)
  110 format (i4,2x,5f8.2)
c
      end
c
c********************** subroutine fixspc *************************
c
      subroutine fixspc(skew,hiecon)
c
      parameter (numspc=500)
      parameter (maxgrd=10)
c
c* common species
c
      double precision spcamp(numspc,maxgrd,5)
      double precision physio(numspc,maxgrd+10)
      double precision maxabu(numspc)
c
      common / species / spcamp,physio,maxabu
c
c* common grads
c
      double precision grdlth(maxgrd)
      double precision grdprd(maxgrd)
      double precision alphad(maxgrd)
      double precision width(maxgrd)
      double precision variab(maxgrd)
      integer grdsam(maxgrd)
      character*1 grdtyp(maxgrd)
      character*8 grdnam(maxgrd)
c
      common / grads / grdlth,grdprd,alphad,width,variab,
     +                 grdsam,grdtyp,grdnam
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
      common / params / numspc,numgrd,numplt,srorf,prorf
c
c* passed
c
      double precision skew
      double precision hiecon
      integer seed
c
c* local
c
      double precision size
      double precision expans
      integer numpts(maxgrd)
      integer spcpnt
      integer index(maxgrd)
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
        if (grdtyp(j) .eq. 'e' .or. grdtyp(j) .eq. 'e') then
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
      write(2,100) i,maxabu(i)
        do 17 j=1,numgrd
        write(2,110) j,(spcamp(i,j,k),k=1,5)
   17   continue
   14 continue
c
      return
c
c
  100 format (i4,2x,f10.0)
  110 format (i4,2x,5f8.2)
c
      end
c
c* coenoflex *********************** subroutine rndplt *******************
c
      subroutine rndplt()
c
      parameter (maxplt=5000)
      parameter (maxgrd=10)
c
c* common grads
c
      double precision grdlth(maxgrd)
      double precision grdprd(maxgrd)
      double precision alphad(maxgrd)
      double precision width(maxgrd)
      double precision variab(maxgrd)
      integer grdsam(maxgrd)
      character*1 grdtyp(maxgrd)
      character*8 grdnam(maxgrd)
c
      common / grads / grdlth,grdprd,alphad,width,variab,
     +                 grdsam,grdtyp,grdnam
c
c* common plots
c
      double precision centrd(maxplt,maxgrd)
      double precision pltprd(maxplt)
c
      common / plots / centrd,pltprd
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
      common / params / numspc,numgrd,numplt,srorf,prorf
c
c* passed
c
      integer seed
c
c* local
c
      character*8 pltlbl
c
c* coenoflex/rndplt ************** one ********************************
c
      write(3,100) (grdnam(i),i=1,numgrd)
      do 10 i=1,numplt
        do 11 j=1,numgrd
        centrd(i,j) = rand() *grdlth(j)
   11   continue
      if (i .lt. 10) then
        write(pltlbl,110) i
      else if (i .ge. 10 .and. i .lt. 100) then
        write(pltlbl,120) i
      else
        write(pltlbl,130) i
      endif
c
      pltprd(i) = 1.0
        do 12 j=1,numgrd
        if (grdprd(j) .ne. 0.0) then
          grdpos = (centrd(i,j)-(grdlth(j)/2.0))/grdlth(j) *
     +               (grdprd(j)/100.0) + 1
          pltprd(i) =  grdpos * pltprd(i)
        endif
   12   continue
      write(3,140) pltlbl,pltprd(i),(centrd(i,j),j=1,numgrd)
   10 continue
c
      return
c
  100 format (10x,'prod. ',10a8)
  110 format ('plt____',i1)
  120 format ('plt___',i2)
  130 format ('plt__',i3)
  140 format (a8,2x,10f8.3,10f8.1)
c
      end
c
c* coenoflex ******************** subroutine fixplt **********************
c
      subroutine fixplt()
c
      parameter (maxplt=5000)
      parameter(maxgrd=10)
c
c* common grads
c
      double precision grdlth(maxgrd)
      double precision grdprd(maxgrd)
      double precision alphad(maxgrd)
      double precision width(maxgrd)
      double precision variab(maxgrd)
      integer grdsam(maxgrd)
      character*1 grdtyp(maxgrd)
      character*8 grdnam(maxgrd)
c
      common / grads / grdlth,grdprd,alphad,width,variab,
     +                 grdsam,grdtyp,grdnam
c
c* common plots
c
      double precision centrd(maxplt,maxgrd)
      double precision pltprd(maxplt)
c
      common / plots / centrd,pltprd
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
      common / params / numspc,numgrd,numplt,srorf,prorf
c
c* passed
c
      integer seed
c
c* local
c
      double precision size
      double precision expans
      double precision grdpos
      integer numpts(maxgrd)
      integer totsam
      integer index(maxgrd)
      character*8 pltlbl(maxplt)
c
c* coenoflex/fixplt ******************** one ******************************
c
      totsam = 0
      size = 1.0
      do 11 i=1,numgrd
      if (grdsam(i) .eq. 1) then
        totsam = totsam + 1
        size = size * grdlth(i)
      endif
   11 continue
c
      expans = (numplt / size)**(1.0/totsam)
c
      numplt = 1
      do 12 i=1,numgrd
      if (grdsam(i) .eq. 1) then
        numpts(i) = nint(expans * grdlth(i))
        numplt = numplt * numpts(i)
      else
        numpts(i) = 0 
      endif
   12 continue
c
      do 10 i=1,numplt
      if (i .lt. 10) then
        write(pltlbl(i),110) i
      else if (i .ge. 10 .and. i .lt. 100) then
        write(pltlbl(i),120) i
      else
        write(pltlbl(i),130) i
      endif
   10 continue
c
      totsam = 1
      do 13 i=1,numgrd
      if (grdsam(i) .eq. 1) then
        index(i) = totsam 
         totsam = totsam * numpts(i)
      endif
   13 continue
c
c     index(1) = 1
c     do 13 i=2,numgrd
c     index(i) = index(i-1) * numpts(i-1)
c  13 continue
c
      write(3,100) (grdnam(i),i=1,numgrd)
      do 14 i=1,numplt
      pltprd(i) = 1.0
        do 15 j=1,numgrd
        if (grdsam(j) .eq. 1) then
          centrd(i,j) = mod((i-1)/index(j),numpts(j)) *
     +                   (grdlth(j)/(numpts(j)-1))
        else
          centrd(i,j) = rand() * grdlth(j)
        endif
        if (grdprd(j) .ne. 0.0) then
          grdpos = (centrd(i,j)-(grdlth(j)/2.0))/grdlth(j) *
     +               (grdprd(j)/100.0) + 1
          pltprd(i) =  grdpos * pltprd(i)
        endif
   15   continue
      write(3,140) pltlbl(i),pltprd(i),(centrd(i,j),j=1,numgrd)
   14 continue
c
      return
c
  100 format (10x,'prod. ',10a8) 
  110 format ('plt____',i1)
  120 format ('plt___',i2)
  130 format ('plt__',i3)
  140 format (a8,2x,f6.3,10f8.1)
c
      end
c
c* coenoflex ******************** subroutine totphy *************************
c
      subroutine totphy(cmpasy,cmpphy,maxtot,noise,slack,final)
c
      parameter (numspc=500)
      parameter (maxgrd=10)
      parameter (maxplt=5000)
c
c* common species
c
      double precision spcamp(numspc,maxgrd,5)
      double precision physio(numspc,maxgrd+10)
      double precision maxabu(numspc)
c
      common / species / spcamp,physio,maxabu
c
c* common plots
c
      double precision centrd(maxplt,maxgrd)
      double precision pltprd(maxplt)
c
      common / plots / centrd,pltprd
c
* common grads
c
      double precision grdlth(maxgrd)
      double precision grdprd(maxgrd)
      double precision alphad(maxgrd)
      double precision width(maxgrd)
      double precision variab(maxgrd)
      integer grdsam(maxgrd)
      character*1 grdtyp(maxgrd)
      character*8 grdnam(maxgrd)
c
      common / grads / grdlth,grdprd,alphad,width,variab,
     +                 grdsam,grdtyp,grdnam
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
      common / params / numspc,numgrd,numplt,srorf,prorf
c
c* common lists
c
      character*8 spclst(numspc)
      character*8 pltlst(maxplt)
c
      common / lists / spclst,pltlst
c
c* passed
c
      double precision cmpasy
      double precision cmpphy
      double precision maxtot
      double precision noise
      double precision slack
      integer final
      integer seed
c
c* local
c
      double precision abunda(numspc)
      double precision abuout(numspc)
      double precision denom
      double precision sumsq
      double precision sum
      integer spcout(numspc)
      integer numout
c
c* coenoflex/phys ****************** one ***************************
c
      write(2,150) (grdnam(i),i=1,numgrd),(spclst(i),i=1,numspc)
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
      call auteco
      call syneco(abunda,final,pltprd(i),noise,slack,
     +            maxtot,cmpasy,cmpphy)
c
        numout = 0
        do 13 j=1,numspc
        if (abunda(j) .gt. 0.1) then
          numout = numout + 1
          spcout(numout) = j
          abuout(numout) =  abunda(j)
        endif
   13   continue
c
        do 16 j=1,numout,8
        write(1,120) i,(spcout(k),abuout(k),k=j,min(numout,j+7))
   16   continue
      write(2,140) pltlst(i),(centrd(i,j),j=1,numgrd),
     +              (abunda(j),j=1,numspc)
   10 continue
c
      write(1,130)
c
      return
c
  120 format (i4,1x,8(i4,f5.1))
  130 format ('0')
  140 format (a10,510f6.1)
  150 format (' plot      ',510a10)
c
      end
c
c* coenoflex/spclbl************** subroutine spclbl ******************
c
      subroutine spclbl
c
      parameter (numspc=500)
      parameter (numplt=5000)
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
      common / params / numspc,numgrd,numplt,srorf,prorf
c
c* common lists
c
      character*8 spclst(numspc)
      character*8 pltlst(numplt)
c
      common / lists / spclst,pltlst
c
c* coenoflex/spclst ****************** one ******************************
c
      do 10 i=1,min(numspc,9)
      write(spclst(i),100) i
   10 continue
c
      do 11 i=10,min(99,numspc)
      write(spclst(i),110) i
   11 continue
c
      do 12 i=100,numspc
      write(spclst(i),120) i
   12 continue
c
      return
c
  100 format ('spc____',i1)
  110 format ('spc___',i2)
  120 format ('spc__',i3)
c
      end
c
c* coenoflex/pltlbl ************ subroutine pltlbl **********************
c
      subroutine pltlbl
c
      parameter (numplt=5000)
      parameter (numspc=500)
c
c* common params
c
      integer numspc
      integer numgrd
      integer numplt
      character*1 srorf
      character*1 prorf
c
      common / params / numspc,numgrd,numplt,srorf,prorf
c
c* common lists
c
      character*8 spclst(numspc)
      character*8 pltlst(numplt)
c
      common / lists / spclst,pltlst
c
      do 10 i=1,min(numplt,9)
      write(pltlst(i),100) i
   10 continue
c
      do 11 i=10,min(99,numplt)
      write(pltlst(i),110) i
   11 continue
c
      do 12 i=100,numplt
      write(pltlst(i),120) i
   12 continue
c
      return
c
  100 format ('plt____',i1)
  110 format ('plt___',i2)
  120 format ('plt__',i3)
c
      end
c
c* subroutine hlppfx *****************************************
c
      subroutine hlppfx
c
      write(6,100)
      return
c
  100 format (/,
     + ' Coenoflex writes five files for each simulation:',/,
     + '   (1) a vegetation data file in matrix format',/,
     + '   (2) a vegetation data file in Cornell condensed format',/,
     + '   (3) a site data file',/,
     + '   (4) a species amplitude file',/,
     + '   (5) a log of all parameters entered for a simulation',//,
     + ' The prefix is used as the name of these files, with',/,
     + ' extensions .veg, .crn, .sit, .spc, and .log respectively.',/,
     + ' The prefix must be suitable for a file name with your',/,
     + ' operating system',/)
c
      end
c
c* subroutine hlpgrd *****************************************
c
      subroutine hlpgrd
c
      write(6,100)
      return
c
  100 format (/,
     + ' Coenoflex simulates up to 10 gradients.  The',/,
     + ' number of gradients should be few compared to',/,
     + ' the number of plots, so that the distribution',/,
     + ' of plots is not too sparse.',//,
     + ' Enter a number between 1 and 10',/)
c
      end
c
c* subroutine hlpspc ********************************************
c
      subroutine hlpspc(numspc)
c
      integer numspc
c
      write(6,100) numspc
      return
c
  100 format (/,
     + ' Coenoflex simulates the distribution of species',/,
     + ' along gradients.  The number of species should be',/,
     + ' chosen to obtain a reasonable species richness ',/,
     + ' given the number of gradients.',//,
     + ' Enter a maximum of ',i3,/)
c
      end
c
c* subroutine hlpplt
c
      subroutine hlpplt(numplt)
c
      integer numplt
c
      write(6,100) numplt
      return
c
  100 format (/,
     + ' Coenoflex simulates the vegetation composition of ',/,
     + ' sample points along environmental gradients.  Enter',/,
     + ' a sufficient number of points to generate a suitable',/,
     + ' density given the number of gradients',//,
     + ' Enter a maximum of ',i3,/)
c
      end
c
c* subroutine hlpasy
c
      subroutine hlpasy
c
      write(6,100)
      return
c
  100 format (/,
     + ' Coenoflex has several mechanisms to simulate ',/,
     + ' competition.  The competition asymmetry coefficient',/,
     + ' simulates a competitive advantage to larger plants',/,
     + ' intended to simulate asymmetric competition such as',/,
     + ' for light.  The algorithm scales species abundance',/,
     + ' proportional to raw abundance raised to the competition',/,
     + ' asymmetry coefficient.  Accordingly, a value of 1.0 ',/,
     + ' confers no advantage, and increasingly larger values',/,
     + ' give large species more advantage.',//,
     + ' Enter a number between 1.0 and 3.0',/)
c
      end
c
c* subroutine hlpstd *************************************************
c
      subroutine hlpstd
c
      write(6,100)
      return
c
  100 format (/,
     + ' Coenoflex includes the option to standardize the mean',/,
     + ' total plot abundance to a specific value.  If the ',/,
     + ' gradients have productivity coefficients specified, this',/,
     + ' mean value will appear at the center of the simulated',/,
     + ' coenospace.',//,
     + ' Standardized plot totals of 100% is preferred by some',/,
     + ' software packages or algorithms.  However, it may mask',/,
     + ' unrealistic vegetation models (too sparse or too dense).',//,
     + ' Enter Y or N.  If Y, enter the standard cover (in %) at the',/,
     + ' prompt.',/)
c
      end
c
c* subroutine hlpnoi *************************************************
c
      subroutine hlpnoi
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex allows you to simulate noise in the species',/,
     +  ' response model.  Values are +/- in percent.  For example,',/,
     +  ' a response of 25 means plus or minus 25% of the actual ',/,
     +  ' value.',//,
     +  ' Enter a number between 0 and 100.',/)
c
      end
c
c* subroutine hlpskw ******************************************************
c
      subroutine hlpskw
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex simulates the maximum abundance of organisms',/,
     +  ' according an equation similar to a log-random ',/,
     +  ' distribution.  The skew coefficient controls how skewed',/,
     +  ' the distribution is toward low values.  A skew of 1.0',/,
     +  ' results in a normal distribution with a mean of 50% cover.',/,
     +  ' A skew of 2.0 results in a mean of about 28% and a Y ',/,
     +  ' intercept of about x. enter value between 1.0 and 4.0',//)
c
      end
c
c* subroutine hlpslk
c
      subroutine hlpslk
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex includes a form of noise known as slack.  In',/,
     +  ' contrast to noise, slack causes a species to be absent',/,
     +  ' in a favorable environment, due perhaps to poor dispersal',/,
     +  ' (or incorrect identification in the field).  The values ',/,
     +  ' are probabilities of absence in favorable conditions, for',/,
     +  ' example a value of 0.1 means a species will be absent in',/,
     +  ' 10% of favorble plots.',//,
     +  ' Enter a value between 0.00 and 1.00.',/)
c
      end
c
c* subroutne hlpxgd ****************************************************
c
      subroutine hlpxgd
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex simulates a number of missing gradients as',/,
     +  ' well as the explicitly represented gradients.  These',/,
     +  ' gradients do not appear in the site data file, and ',/,
     +  ' represent unmeasured environmental variability.  All',/,
     +  ' real vegetation data have such unmeasured variation',/,
     +  ' but you may choose to include or exclude this variation.',//,
     +  ' Enter a number betwen 0 and 10',/)
c
      end
c
c* subroutine hlpsed ****************************************************
c
      subroutine hlpsed
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex includes stochastic variation for random',/,
     +  ' species centroids and random plot locations.',/,
     +  ' Accordingly, it requires a random number seed.',//,
     +  ' Enter any 6-digit odd number.',/)
c
      end
c
c
c* subroutine hlpgnm
c
      subroutine hlpgnm
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex produces environment files for the simulated',/,
     +  ' environment in either fuzphy or cornell format.  in either',/,
     +  ' case, the environment variables are labeled with a name',/,
     +  ' up to 8 chracters (beginning with a / for fuzphy) which',/,
     +  ' you provide.  choose a name or acronym that is descriptive',/,
     +  ' within the 8 characters allowed.',/)
c
      end
c
c* subroutine hlpglt
c
      subroutine hlpglt
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex calculates gradient lengths in arbitrary units.',/,
     +  ' these units are used to compare gradient lengths in',/,
     +  ' relative values, so that for example a gradient of 200 ',/,
     +  ' units is twice as long as a gradient of 100 units.',//,
     +  ' in addition, species amplitudes are measured in the same',/,
     +  ' units so that a species with an amplitude of 50 would ',/,
     +  ' extend about one half the length of a gradient of 100 ',/,
     +  ' units.',//,
     +  ' thus, species mean amplitudes scaled according to gradient',/,
     +  ' lengths controls species turnover rates.',/)
c
      end
c
c* subroutine hlpgtp
c
      subroutine hlpgtp
c
      write(6,100)
      return
c
  100 format(/,
     +  ' Coenoflex uses two types of gradients, following austin',/,
     +  ' and smith (198x).  the first type of gradient is an ',/,
     +  ' environmental gradient which includes both direct ',/,
     +  ' environmental gradients and complex environmental',/,
     +  ' gradients sensu austin and smith.  species have modal',/,
     +  ' responses along these gradients with minimum, optimum, and',/,
     +  ' maximum values along the gradient (although not necessarily',/,
     +  ' within the range of gradient being sampled).  species have',/,
     +  ' independent optima chosen randomly from a sampling ',/,
     +  ' distribution depending on the alpha-diversity parameter ',/,
     +  ' for that gradient.',//,
     +  ' resource gradients differ in that all species have their',/,
     +  ' optimum at the high end of the gradient, with ',/,
     +  ' independently varying minima.  resource gradients should',/,
     +  ' have a productivity increase along then, specified with ',/,
     +  ' the gradient productivity coefficient.',/)
c
      end
c
c* subroutine hpald
c
      subroutine hlpald
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex can systematically vary the alpha-diversity along',/,
     +  ' gradients by modifying the location of species centroids.',/,
     +  ' values from 0.5 to 1.0 increase the alpha-diversity at the',/,
     +  ' high end of the gradient, and values from 1.0 to 2.0 ',/,
     +  ' increase alpha-diversity at the low end of the gradient.',/,
     +  ' values near 1.0 have small effects; the farther from 1.0',/,
     +  ' the more dramatic the effect.',/)
c
      end
c
c* subroutine hlpamp
c
      subroutine hlpamp
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex calculates the ecological amplitude of species in',/,
     +  ' gradient units.  thus, if the gradient length is 100, and',/,
     +  ' a species amplitude is 50, it will extend about one half ',/,
     +  ' the length of the gradient.  species amplitude, along with',/,
     +  ' the number of species, determines the number of species ',/,
     +  ' per plot and the rate of species turnover.',//,
     +  ' enter a number in gradient units.',/)
c
     +
      end
c
c* subroutine hlpprd
c
      subroutine hlpprd
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex can vary the productivity along gradients ',/,
     +  ' systematically.  to increase productivity along the ',/,
     +  ' gradient, enter the percent increase in productivity',/,
     +  ' from the low to high end of the gradient.  resource',/,
     +  ' gradients should have a significant increase associated',/,
     +  ' with them; environmental gradients may or may not,',//,
     +  ' enter a percent increase from 0 to 200%.')
c
      end
c
c* subroutine hlpvar
c
      subroutine hlpvar
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex varies individual species amplitudes according to ',/,
     +  ' the variance parameter.  previously you specified a mean',/,
     +  ' species amplitude along each gradient.  the variance',/,
     +  ' parameter modifies that value for each species by +/- %.',//,
     +  ' for example, if you specified a mean amplitude of 50 and a ',/,
     +  ' variance of 50%, species amplitudes would vary randomly ',/,
     +  ' from 25 units to 75 units.',//,
     +  ' enter a variance in percent',/)
c
      end
c
c* subroutine hlpprf
c
      subroutine hlpprf
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex can establish plot centers randomly or along',/,
     +  ' a grid.  if you select randomly, plot centers will have',/,
     +  ' random coordinates along each gradient with approximately',/,
     +  ' constant density throughout the sample space.',//,
     +  ' if you select a grid, Coenoflex will attempt to space plots',/,
     +  ' along the gradients in proportion to their length, again',/,
     +  ' attempting to maintain approximately constant density.  if',/,
     +  ' necessary, Coenoflex will change the number of plots',/,
     +  ' specified to achieve equal spacing.',//,
     +  ' enter a g (for grids) or r (for random).',/)
c
      end
c
c* subroutine hlpsrf
c
      subroutine hlpsrf
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex can establish species modes randomly or along',/,
     +  ' a grid.  if you select randomly, species modes will have',/,
     +  ' random coordinates along each gradient with approximately',/,
     +  ' constant density throughout the sample space.',//,
     +  ' if you select a grid, Coenoflex will attempt to space ',/,
     +  ' species modes along the gradients in proportion to their',/,
     +  ' length, again attempting to maintain approximately',/,
     +  ' constant density.  if necessary, Coenoflex will change',/,
     +  ' the number of species specified to achieve equal spacing.',//,
     +  ' enter a g (for grids) or r (for random).',/)
c
      end
c
c* subroutine hlpphy
c
      subroutine hlpphy
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex calculates physiological-based competition as ',/,
     +  ' determined by how optimal an environment is.  in other',/,
     +  ' words, species nearer their phsyiological optimum are',/,
     +  ' more competitive than species in less favorable',/,
     +  ' environments irrespective of abundance. the coefficient',/,
     +  ' represents a penalty, with higher values indicating',/,
     +  ' more severe competition effects.',//,
     +  ' enter a value from 1.0 to 3.0')
c
      end
c
c* subroutine hlpaut
c
      subroutine hlpaut
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex combines the physiological response of a species ',/,
     +  ' along each gradient into an autecological response using ',/,
     +  ' ',/,
     +  ' min(x,y,...) = minimum response of gradients x,y,... ',/,
     +  ' max(x,y,...) = maximum response of gradients x,y,... ',/,
     +  ' irm(x,y,...) = integrated rate methodology response ',/,
     +  '                                 of gradients x,y,... '/,
     +  ' geo(x,y,...) = geometric mean response of gradients ',/,
     +  '                                              x,y,... ',/,
     +  ' ave(x,y,...) = arithmetic mean rsponse of gradients ',/,
     +  '                                 of gradients x,y,... ',//,
     +  ' functions can be nested as e.g.  min(1,irm(2,3)) ',//,
     +  ' enter a function below ')
c
       end
c
c* subroutine hlphcn
c
      subroutine hlphcn
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex allows species abundances to correlated',/,
     +  ' with their ecological amplitude as specified in the ',/,
     +  ' hierarchical continuum model (Collins et al. 1993',/,
     +  ' J. Veg. Sci. x:xx-xx.) ',//,
     +  ' Values greater than 0.0 result in positive correlation',/,
     +  ' between amplitude and abundancem while 0.0 results',/,
     +  ' in independence ',//,
     +  ' enter [0.00-1.00 ')
c
      end
c
c* subroutine hlpsam
c
      subroutine hlpsam
c
      write(6,100)
      return
c
  100 format (/,
     +  ' Coenoflex can sample gradients independently from ',/,
     +  ' the actial gradients used in the simulation.',/,
     +  ' This way, hidden gradients are easily included ',/,
     +  ' realistically, and non-used demmy gradients can',/,
     +  ' also be used.  The list can only include gradients',/,
     +  ' up to the number specified earlier, however.',//,
     +  ' Enter as a comma-separated list. e.g.  1,2,3,4')
c
      end
c
* subroutine percen ***********************************************************
c
      subroutine percen(respon)
c
      character*20 respon
c
      do 10 i=1,20
      if (respon(i:i) .eq. '%') respon(i:i) = ' '
   10 continue
c
      return
c
      end
c
