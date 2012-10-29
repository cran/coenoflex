      subroutine autpar(line,argmnt,grdlst,numper,count)
c
c* arglst
c
      integer argmnt(10)    ! LIST OF COMMAND LINE ARGUMENTS
      integer grdlst(10,10)
      integer numper(10)
      integer count            ! NUMBER OF ARGUMENTS PASSED
c
c* passed
c
      character*255 line          ! TEXT LINE CONTAINING ARGUMENTS
c
c* local
c
      integer stkpnt,numlft
      character*3 tmparg
c
c* autpar ********************* one *****************************
c
      maxnst = 0
      stkpnt = 0
      count = 0
c
      call tolower(line)
c
   10 numlft = 0
      maxnst = 0
      do 11 i=1,255
      if (line(i:i) .eq. '(') then
        numlft = numlft + 1
        maxnst = max(maxnst,numlft)
      else if (line(i:i) .eq. ')') then
        numlft = numlft - 1
      endif
   11 continue
c
c     if (numlft .ne. 0) write(6,*) ' unbalanced parentheses '
c
      numlft = 0
      do 12 i=1,255
      if (line(i:i) .eq. '(') then
        numlft = numlft + 1
        if (numlft .eq. maxnst) then
          count = count + 1
          tmparg = line(i-3:i-1)
          if (tmparg .eq. 'ave') then
            argmnt(count) = 1
          else if (tmparg .eq. 'min') then
            argmnt(count) = 2
          else if (tmparg .eq. 'max') then
            argmnt(count) = 3
          else if (tmparg .eq. 'geo') then
            argmnt(count) = 4
          else
            argmnt(count) = 5
          endif
          do 13 j=i,255
          if (line(j:j) .eq. ',' .or. line(j:j) .eq. ')') then
            numper(count) = numper(count) + 1
            if (line(j-2:j-2) .ne. '(' .and.
     +          line(j-2:j-2) .ne. ',') then
              read(line(j-2:j-1),'(i2)') grdlst(count,numper(count))
            else
              read(line(j-1:j-1),'(i1)') grdlst(count,numper(count))
            endif
          endif
          if (line(j:j) .eq. ')') then
            stkpnt = stkpnt + 1
            do 14 k=i-3,j
            line(k:k) = ' '
   14       continue
            write(line(i-3:i-1),'(i2)')  stkpnt + 10
            call collap(line)
            goto 10
          endif
   13     continue
        endif
      else if (line(i:i) .eq. ')') then
        numlft = numlft - 1
      endif
   12 continue
c
      return
c
      end
c
c* coenoflex/tolower ******************************************************
c
      subroutine tolower(line)
c
      character*255 line
c
      do 10 i=1,255
      if (ichar(line(i:i)) .ge. 65 .and.
     +    ichar(line(i:i)) .le. 90) then
        line(i:i) = char(ichar(line(i:i)) + 32)
      endif
   10 continue
c
      return
c
      end
c
c* coenoflex/collap **********************************************************
c
      subroutine collap(line)
c
      character*255 line
c
   10 do 11 i=255,1,-1
      if (line(i:i) .ne. ' ') then
        do 12 j=i,1,-1
        if (line(j:j) .eq. ' ') then
          do 13 k=j,i
          line(k:k) = line(k+1:k+1)
   13     continue
          goto 10
        endif
   12   continue
      endif
   11 continue
c
      return
c
      end
