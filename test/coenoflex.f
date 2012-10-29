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
