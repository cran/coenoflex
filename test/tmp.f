      subroutine coenoflex(numplt,numspc,numgrd,spcamp,physio,numabu,
     +           grdlth,grdprd,alphad,width,variab,grdsam,grdtyp,
     +           grdnam,centrd,pltprd,spclst,pltlst,
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
      character*1 grdtyp(numgrd)
      character*8 grdnam(numgrd)
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
      character*1 srorf
      character*1 prorf
c
c* common arglst
c
      character*20 argmnt(10)    ! LIST OF COMMAND LINE ARGUMENTS
      integer grdlst(10,10)
      integer numper(10)
      integer count            ! NUMBER OF ARGUMENTS PASSED
c
c* common lists
c
      character*8 spclst(numspc)
      character*8 pltlst(numplt)
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
      character*80 line
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
   10 write(6,'('' enter the file name prefix'',t60,'' : '')',
     +      advance='NO')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlppfx
        goto 10
      endif
      prefix = respon
      do 11 i=1,20
      if (prefix(i:i) .eq. ' ') then
        lenprf = i-1
        goto 12
      endif
   11 continue
      lenprf = 20
c
   12 open (unit=4,file=prefix(1:lenprf) // '.log',status='unknown')
      write(4,'('' enter the file name prefix'',t59,'' : '',a20)',
     +     advance='NO') prefix
c
      open (unit=1,file=prefix(1:lenprf) // '.veg',status='unknown')
      open (unit=2,file=prefix(1:lenprf) // '.spc',status='unknown')
      open (unit=3,file=prefix(1:lenprf) // '.sit',status='unknown')
c
   13 write(6,
     +   '(/,'' enter the number of gradients'',t54,''[1-10] : '')',
     +   advance='NO')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpgrd
        goto 13
      else
        read(respon,'(i2)',err=13) numgrd
        if (numgrd .gt. maxgrd) then
          call hlpgrd
          goto 13
        endif
      endif
      write(4,
     +  '(/,'' enter the number of gradients'',t54,''[1-10] : '',i4)',
     +  advance='NO') numgrd
  103 format(/,' enter the number of gradients',t54,'[1-10] : ',i4)
c
   14 write(6,
     + '(/,'' enter the number of species'',t53,''[1-500] : '')',
     +  advance='NO')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpspc(maxspc)
        goto 14
      else
        read(respon,'(i3)',err=14) numspc
        if (numspc .gt. maxspc) then
          call hlpspc(maxspc)
          goto 14
        endif
      endif
      write(4,
     +   '(/'' enter the number of species'',t53,''[1-500] : '',i4)')
     +   numspc
c
   15 write(6,'(/,'' enter the number of plots'',t53,''[1-500] : '')',
     +      advance='NO')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpplt(maxplt)
        goto 15
      else
        read(respon,'(i3)',err=15) numplt
        if (numplt .gt. maxplt) then
          call hlpplt(maxplt)
          goto 15
        endif
      endif
      write(4,'('' enter the number of plots'',t53,''[1-500] : '',i4)')
     +      numplt
c
c* coeno ************************* two **********************************
c
   20 write(6,
     +'(/,'' enter the competition asymmetry'',t49,'
     +  // '''[1.00-3.00] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpasy
        goto 20
      else
        read(respon,'(f4.2)',err=20) cmpasy
      endif
      write(4,
     +'(/,'' enter the competition asymmetry'',t49,''[1.00-3.00] : '','
     +   // 'f4.2)') cmpasy
c
   29 write(6,
     +'('' enter the physiological competition coeff.'','
     +      // 't49,''[1.00-3.00] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpphy
        goto 29
      else
        read(respon,'(f4.2)',err=20) cmpphy
      endif
      write(4,
     +'('' enter the physiological competition coeff.'','
     +      // 't49,''[1.00-3.00] : '',f5.2)') cmpphy
c
   24 write(6,'('' enter the skew'',t49,''[1.00-4.00] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpskw
        goto 24
      else
        read(respon,'(f4.2)',err=24) skew
      endif
      write(4,'('' enter the skew'',t49,''[1.00-4.00] : '',f4.2)')
     +       skew
c
   25 write(6,
     +'('' enter the amplitude/abundance correlation'','
     +    // 't49,''[0.00-1.00] : '',$)') 
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then 
        call hlphcn
        goto 25
      else
        read(respon,'(f5.3)') hiecon
      endif
      write(4,
     +'('' enter the amplitude/abundance correlation'','
     +    // 't49,''[0.00-1.00] : '',f5.3)') hiecon
c
   21 write(6,
     +'('' standardize mean total plot abundance? '','
     +    // 't52,''[y or n] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpstd
        goto 21
      else
        read(respon,'(a1)',err=21) yorn
        if (yorn .eq. 'y' .or. yorn .eq. 'y') then
          write(6,'('' enter the standard total percent cover.'','
     +    //  't51,''[10-500%] : '',$)')
          read(5,'(a20)') respon
          call percen(respon)
          read(respon,'(f4.0)',err=21) maxtot
        else if (yorn .eq. 'n' .or. yorn .eq. 'n') then
          maxtot = 0.0
        else
          call hlpstd
          goto 21
        endif
      endif
      write(4,
     +'('' standardize mean total plot abundance? '','
     +    // 't52,''[y or n] : '',3x,a1)') yorn
      if (yorn .eq. 'y' .or. yorn .eq. 'y') then
        write(4,'('' enter the standard total percent cover.'','
     +    //  't51,''[10-500%] : '',f4.0)') maxtot
      endif
c
   22 write(6,'('' enter the noise'',t52,''[0-100%] : '',$)')
      read(5,'(a20)') respon
      call percen(respon)
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpnoi
        goto 22
      else
        read(respon,'(f4.0)',err=22) noise
      endif
      write(4,'('' enter the noise'',t52,''[0-100%] : '',f4.0)')
     +      noise
c
   23 write(6,'('' enter the slack'',t49,''[0.00-1.00] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpslk
        goto 23
      else
        read(respon,'(f4.2)',err=23) slack
      endif
      write(4,'('' enter the slack'',t49,''[0.00-1.00] : '',f4.2)')
     +    slack
c
   26 write(6,'('' enter a random number seed'','
     +    // 't47,''[6 digit odd] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpsed
        goto 26
      else
        read(respon,'(i6)',err=26) seed
      endif
      write(4,'('' enter a random number seed'','
     +    // 't47,''[6 digit odd] : '',i6)') seed
      test = rand(seed)
c
c************************** three ******************************
c
      write(6,'(/,'' describe each gradient in turn. '')')
      write(4,'(/,'' describe each gradient in turn. '')')
c
      do 30 i=1,numgrd
   31 write(6,'(/,'' enter a name for gradient '',i3,t60,'' : '',$)') i
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpgnm
        goto 31
      else
        read(respon,'(a8)',err=31) grdnam(i)
      endif
      write(4,'(/,'' enter a name for gradient '',i3,t60,'' : '',a8)')
     +    i,grdnam(i)
c
   32 write(6,'('' enter gradient type for gradient '',i3)') i
      write(6,'('' enter e for environment or r for resource'')')
      write(6,'(t52,''[e or r] : '',$)') 
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpgtp
        goto 32
      else
        read(respon,'(a1)',err=32) grdtyp(i)
        if (grdtyp(i) .ne. 'e' .and. grdtyp(i) .ne. 'e' .and.
     +      grdtyp(i) .ne. 'r' .and. grdtyp(i) .ne. 'r') then
          call hlpgtp
          goto 32
        endif
      endif
      write(4,'('' enter gradient type for gradient '',i3)')
      write(4,'('' enter e for environment or r for resource'')')
      write(4,'(t52,''[e or r] : '',3x,a1)') grdtyp(i)
c
   33 write(6,'('' enter gradient length for gradient '','
     +  // 'i3,t51,''[10-1000] : '',$)') i
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpglt
        goto 33
      else
        read(respon,'(f4.0)',err=33) grdlth(i)
        if (grdlth(i) .eq. 0.0) then
          call hlpglt
          goto 33
        endif
      endif
      write(4,'('' enter gradient length for gradient '','
     +  // 'i3,t51,''[10-1000] : '',f4.0)') i,grdlth(i)
c
   34 write(6,'('' enter the mean width of species amplitudes'','
     +   //   't52,''[10-100] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpamp
        goto 34
      else
        read(respon,'(f4.0)',err=34) width(i)
        if (width(i) .eq. 0.0) then
          call hlpamp
          goto 34
        endif
      endif
      write(4,'('' enter the mean width of species amplitudes'','
     +   //   't52,''[10-100] : '',f4.0)') width(i)
c
   35 write(6,'('' enter the variability'',t52,''[0-50%] : '',$)')
      read(5,'(a20)') respon
      call percen(respon)
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpvar
        goto 35
      else
        read(respon,'(f4.0)',err=35) variab(i)
      endif
      write(4,
     +  '('' enter the variability'',t52,''[0-50%] : '',f4.0)')
     +  variab(i)
c
   36 write(6,'('' enter the productivity response'')')
      write(6,'('' along the gradient'',
     +  t52,''[0-200%] : '',$)')
      read(5,'(a20)') respon
      call percen(respon)
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpprd
        goto 36
      else
        read(respon,'(f4.0)',err=36) grdprd(i)
      endif
      write(4,'('' enter the productivity response'')')
      write(4,'('' along the gradient'',
     +     t52,''[0-200%] : '',f4.0)') grdprd(i)
c
   37 write(6,'('' enter the trend in alpha-diversity '','
     +     // 't51,''[0.5-2.0] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?' .or. respon .eq. ' ') then
        call hlpald
        goto 37
      else
        read(respon,'(f4.0)',err=37) alphad(i)
      endif
      write(4,'('' enter the trend in alpha-diversity '','
     +     // 't51,''[0.5-2.0] : '',f4.2)') alphad(i)
   30 continue
c
c************************** four *************************************
c
   43 write(6,'(/,'' random or grid species modes'')')
      write(6,'(t52,''[r or g] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?') then
        call hlpsrf
        goto 43
      else
        read(respon,'(a1)',err=41) srorf
        if (srorf .eq. 'R' .or. srorf .eq. 'r') then
          call rndspc(skew,hiecon)
        else if (srorf .eq. 'G' .or. srorf .eq. 'g') then
          call fixspc(skew,hiecon)
        else
          call hlpsrf
          goto 43
        endif
      endif
      write(4,'(/,'' random or grid species modes'')')
      write(4,'(t52,''[r or g] : '',3x,a1,$)') srorf
c
   40 write(6,'(/,'' enter the autecological function'','
     +    // '/,'' : '',$)') 
      read(5,'(a80)') line
      if (line(1:1) .eq. '?' .or. line(1:1) .eq. ' ') then
        call hlpaut
        goto 40
      else
        write(4,'(/,'' enter the autecological function'','
     +    // '/,'' : '',a80)') line
        call autpar(line)
        final = count + 10
      endif
c
   41 write(6,'('' enter the list of gradients to sample'','
     +   // '/'' : '',$)')
      read(5,'(a80)') line
      if (line(1:1) .eq. '?' .or. line(1:1) .eq. ' ') then
        call hlpsam
        goto 40
      else
        write(4,'(/,'' enter the list of gradients to sample'','
     +    // '/,'' : '',a80)') line
        call sampar(line,grdsam)
      endif
c
   42 write(6,'(/,'' random or grid plot locations'','
     +    //    't52,''[r or g] : '',$)')
      read(5,'(a20)') respon
      if (respon(1:1) .eq. '?') then
        call hlpprf
        goto 42
      else
        read(respon,'(a1)',err=41) prorf
        if (prorf .eq. 'R' .or. prorf .eq. 'r') then
          call rndplt
        else if (prorf .eq. 'G' .or. prorf .eq. 'g') then
          call fixplt
        else
          call hlpprf
          goto 42
        endif
      endif
      write(4,'(/,'' random or grid plot locations'','
     +    //    't52,''[r or g] : '',3x,a1)') prorf
c
      write(1,400) numgrd,numspc,srorf,numplt,prorf
      write(1,410)
c
      open (unit=2,file=prefix(1:lenprf) // '.dat',status='unknown')
c
      call spclbl
      call pltlbl
      call totphy(cmpasy,cmpphy,maxtot,noise,slack,final)
      write(1,420) (spclst(i),i=1,numspc)
      write(1,420) (pltlst(i),i=1,numplt)
c
  400 format (' gradients = ',i3,' species = ',i3,1x,a1,
     +        ' plots = ',i3,1x,a1)
  410 format ('(i4,1x,8(i4,f5.1))',/,t4,'8')
  420 format (10a8)
c
      end
