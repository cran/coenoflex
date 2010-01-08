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
c
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
c
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
      double precision maxval
c
c* coenoflex/rndspc *************** one *****************************
c
      maxval = 0.0
      do 10 i=1,numspc
      if (skew .eq. 0) then
        maxabu(i) = 1.0
      else
        maxabu(i) = 0.0
        do 11 j=1,3
        maxabu(i) = maxabu(i) + rand()
   11   continue
        maxabu(i) = (maxabu(i)/3.0)**skew 
      endif
      maxval = max(maxval,maxabu(i))
   10 continue
c
      do 12 i=1,numspc
      maxabu(i) = maxabu(i) / maxval * 100.0
   12 continue
c
      do 13 i=1,numspc
      hcnadj = 1.0 + ((maxabu(i)/100.0)-0.5) * hiecon
        do 14 j=1,numgrd
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
          spcamp(i,j,1) = spcamp(i,j,2)-((spcamp(i,j,3)-spcamp(i,j,2)))
          spcamp(i,j,4) = grdlth(j)
          spcamp(i,j,5) = grdlth(j)
c         maxabu(i) = min(100.0,maxabu(i) * (1.5 - (1.0-center)))
          maxabu(i) = maxabu(i) * 
     +              (1-((grdlth(j)-spcamp(i,j,2))/grdlth(j)))
        endif
   14   continue
c
   13 continue
c
      return
c
      end
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
        if (grdtyp(j) .eq. 1) then
          center = mod((i-1)/index(j),numpts(j)) *
     +                (range/(numpts(j)-1)) - width(j)/2
c         spcamp(i,j,3) = center * range - (width(j)/2.0)
          spcamp(i,j,3) = center 
          fudge = (rand() - 0.5) * variab(j)/50.0 * width(j)
          spcamp(i,j,1) = spcamp(i,j,3) - width(j)*hcnadj + fudge
          fudge = (rand() - 0.5) * variab(j)/50.0 * width(j)
          spcamp(i,j,5) = spcamp(i,j,3) + width(j)*hcnadj + fudge
          spcamp(i,j,2) = (spcamp(i,j,1) + spcamp(i,j,3)) / 2.0
          spcamp(i,j,4) = (spcamp(i,j,3) + spcamp(i,j,5)) / 2.0
        else
          center = mod((i-1)/index(j),numpts(j)) *
     +             (grdlth(j)/(numpts(j)-1)) - width(j)/2
          spcamp(i,j,2) = center 
          spcamp(i,j,3) = grdlth(j)
          spcamp(i,j,1) = spcamp(i,j,2)-(spcamp(i,j,3)-spcamp(i,j,2))
          spcamp(i,j,4) = grdlth(j)
          spcamp(i,j,5) = grdlth(j)
c         maxabu(i) = min(100.0,maxabu(i) * (1.5 - (1.0-center)))
          maxabu(i) = maxabu(i) * 
     +              (1-((grdlth(j)-spcamp(i,j,2))/grdlth(j)))
        endif
   15   continue
   14 continue
c
      return
c
      end
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
      integer stkpnt,numlft,maxlin
      character*3 tmparg
c
c* autpar ********************* one *****************************
c
      maxnst = 0
      stkpnt = 0
      count = 0
      maxlin = 0
c
      call tolower(line)
c
   10 numlft = 0
      maxnst = 0
      do 11 i=1,80
      if (line(i:i) .eq. '(') then
        numlft = numlft + 1
        maxnst = max(maxnst,numlft)
      else if (line(i:i) .eq. ')') then
        numlft = numlft - 1
        if (numlft .eq. 0) maxlin = i
      endif
   11 continue
c
c     if (numlft .ne. 0) write(6,*) ' unbalanced parentheses '
c
      numlft = 0
      do 12 i=1,maxlin
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
c
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
c
      subroutine totphy(numplt,numspc,numgrd,centrd,spcamp,physio,
     +                  argmnt,grdlst,numper,count,
     +                  maxabu,abunda,pltprd,noise,slack,maxtot,
     +                  cmpasy,cmpphy,diff)
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
c* passed through
c
      integer argmnt(10)
      integer grdlst(10,10)
      integer numper(10)
      integer count
      double precision maxabu(numspc) 
      double precision abunda(numplt,numspc)
      double precision pltprd(numplt)
      double precision noise
      double precision slack
      double precision maxtot
      double precision cmpasy
      double precision cmpphy
      double precision diff(numspc)
      integer i
c
c* local
c
      double precision denom
c
c* coenoflex/phys ****************** one ***************************
c
      do 10 i=1,numplt
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
     +             centrd(i,k) .le. spcamp(j,k,4)) then
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
      call auteco(numspc,numgrd,argmnt,grdlst,numper,count,physio)
      call syneco(numplt,numspc,numgrd,physio,maxabu,
     +                 abunda,count+10,pltprd,noise,slack,
     +                 maxtot,cmpasy,cmpphy,diff,i)
   10 continue
c
      return
c
      end
c
      subroutine syneco(numplt,numspc,numgrd,physio,maxabu,
     +                 abunda,final,pltprd,noise,slack,
     +                 maxtot,cmpasy,cmpphy,diff,plot)
c
c* params 
c 
      integer numplt
      integer numspc
      integer numgrd
c
c* species
c
      double precision physio(numspc,numgrd+10)
      double precision maxabu(numspc)
c
c* passed
c
      double precision abunda(numplt,numspc)
      integer final
      double precision pltprd(numplt)
      double precision noise
      double precision slack
      double precision maxtot
      double precision cmpasy
      double precision cmpphy
      double precision diff(numspc)
      integer plot
c
c* local
c
      double precision extra
      double precision sumdif
c
c* coenoflex/syneco ************ one *******************************
c
      sum = 0.0
      wgtsum = 0.0
      sumdif = 0.0
c
      do 10 i=1,numspc
      if (physio(i,final) .le. 0.0) then
        abunda(plot,i) = 0.0
      else
        test = rand()
        if (test .lt. slack) then
          abunda(plot,i) = 0.0
          goto 10
        endif
        abunda(plot,i) = maxabu(i) * physio(i,final) * pltprd(plot)
        abunda(plot,i) = abunda(plot,i) +
     +        (rand()-0.5) * noise/50.0 * abunda(plot,i)
        sum = sum + abunda(plot,i)
        diff(i) = (1.0-physio(i,final))**cmpphy * abunda(plot,i)
        sumdif = sumdif + diff(i)
      endif
   10 continue
c
      if (maxtot .ne. 0.0) then
        if (sum .gt. maxtot*pltprd(plot) .and. sumdif .gt. 0.0) then
          extra = sum-(maxtot*pltprd(plot))
          do 11 i=1,numspc
          if (abunda(plot,i) .le. 0) goto 11
          abunda(plot,i) = abunda(plot,i) - diff(i)/sumdif * extra
   11     continue
        endif
c
        do 12 i=1,numspc
        if (abunda(plot,i) .gt. 0.0) then
          wgtsum = wgtsum + abunda(plot,i)**cmpasy
        endif
   12   continue
c
        do 13 i=1,numspc
        abunda(plot,i) = max(0.0,abunda(plot,i))
        if (abunda(plot,i) .gt. 0.0) then
          abunda(plot,i) = maxtot * 
     +             min(100.0,(abunda(plot,i)**cmpasy/wgtsum))
        endif
   13   continue
      endif
c
      return
c
      end

