      subroutine syneco(numplt,numspc,numgrd,physio,maxabu,
     +                 abunda,final,pltprd,noise,slack,
     +                 maxtot,cmpasy,cmpphy,diff)
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
      double precision abunda(numspc)
      integer final
      double precision pltprd
      double precision noise
      double precision slack
      double precision maxtot
      double precision cmpasy
      double precision cmpphy
c
c* local
c
      double precision extra
      double precision sumdif
      double precision diff(numspc)
c
c* coenoflex/syneco ************ one *******************************
c
      sum = 0.0
      wgtsum = 0.0
      sumdif = 0.0
c
      do 10 i=1,numspc
      if (physio(i,final) .le. 0.0) then
        abunda(i) = 0.0
      else
        test = rand()
        if (test .lt. slack) then
          abunda(i) = 0.0
          goto 10
        endif
        abunda(i) = maxabu(i) * physio(i,final) * pltprd
        abunda(i) = abunda(i) +
     +        (rand()-0.5) * noise/50.0 * abunda(i)
        sum = sum + abunda(i)
        diff(i) = (1.0-physio(i,final))**cmpphy * abunda(i)
        sumdif = sumdif + diff(i)
      endif
   10 continue
c
      if (maxtot .ne. 0.0) then
        if (sum .gt. maxtot*pltprd .and. sumdif .gt. 0.0) then
          extra = sum-(maxtot*pltprd)
          do 11 i=1,numspc
          if (abunda(i) .le. 0) goto 11
          abunda(i) = abunda(i) - diff(i)/sumdif * extra
   11     continue
        endif
c
        do 12 i=1,numspc
        if (abunda(i) .gt. 0.0) then
          wgtsum = wgtsum + abunda(i)**cmpasy
        endif
   12   continue
c
        do 13 i=1,numspc
        abunda(i) = max(0.0,abunda(i))
        if (abunda(i) .gt. 0.0) then
          abunda(i) = maxtot * min(100.0,(abunda(i)**cmpasy/wgtsum))
        endif
   13   continue
      endif
c
      return
c
      end

