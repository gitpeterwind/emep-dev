! <My_Reactions.inc - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 


    do iter = 1, toiter(k)
!
! The chemistry is iterated several times, more close to the ground
! than aloft.
! For some reason, it proved faster for some compilers to include files as given below
! with the if statements, than to use loops.
!if(k>=KCHEMTOP)then

   include 'My_FastReactions.inc'

!endif
!if(k>=6)then
!
!   include 'My_FastReactions.inc'
!
!endif
!
!if(k>=KEMISTOP)then
!
!   include 'My_FastReactions.inc'
!
!endif
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   end do !! End iterations
 ! Just before SO4
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  include 'My_SlowReactions.inc'

!ORIG!-> SO4 
!ORIG
!ORIG      P =   &
!ORIG      2e-12*aqrck(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
!ORIG       + aqrck(ICLRC1,K)*XNEW(H2O2) * xnew(SO2 )  &
!ORIG       + aqrck(ICLRC2,K)*XNEW(O3) * xnew(SO2 )  &
!ORIG       + aqrck(ICLRC3,K) * xnew(SO2 )  &
!ORIG      + 0.05*rcemis(QRCSO2,K) 
!ORIG      ! L = 0.0
!ORIG 
!ORIG
!ORIG       xnew(SO4)=  max(0.0,  xold(SO4) + dt2 * P)
!ORIG
!ORIG
!ORIG!-> pNO3
!ORIG
!ORIG      P =   &
!ORIG        rcmisc(10,k) * xnew(HNO3)
!ORIG      L = 0.0
!ORIG
!ORIG
!ORIG       xnew(pNO3)=  max(0.0,  xold(pNO3) + dt2 * P)
!ORIG
!ORIG
!ORIG!-> NH3 
!ORIG
!ORIG      P =   &
!ORIG     rcemis(QRCNH3,K) 
!ORIG      ! L = 0.0
!ORIG 
!ORIG
!ORIG       xnew(NH3)=  max(0.0,  xold(NH3) + dt2 * P)
!ORIG
!ORIG!-> PPM25
!ORIG
!ORIG      P =   &
!ORIG     rcemis(QRCPM25,K) 
!ORIG      ! L = 0.0
!ORIG       xnew(PPM25)=  max(0.0,  xold(PPM25) + dt2 * P)
!ORIG
!ORIG
!ORIG!-> PPMCO
!ORIG
!ORIG      P =   &
!ORIG     rcemis(QRCPMCO,K) 
!ORIG      ! L = 0.0
!ORIG       xnew(PPMCO)=  max(0.0,  xold(PPMCO) + dt2 * P)
!ORIG
!ORIG!-> Rn222
!ORIG
!ORIG      P =  rc_Rn222(k)
!ORIG      L = 1.0/(5.51614*24.0*3600.0)
!ORIG       xnew(Rn222)=  max(0.0, ( xold(Rn222) + dt2 * P))  &
!ORIG                     /(1.0 + dt2*L ) 
!ORIG
!ORIG!-> Pb210
!ORIG      P = 1.0/(5.51614*24.0*3600.0) * xnew(Rn222)
!ORIG      ! L = 0.0
!ORIG       xnew(Pb210)=  max(0.0,  xold(Pb210) + dt2 * P)
!ORIG
!ORIG!SeaS
!ORIG   if (SEASALT) then
!ORIG
!ORIG!-> Sea salt fine
!ORIG
!ORIG      P =   &
!ORIG     rcss(QSSFI,K) 
!ORIG      ! L = 0.0
!ORIG       xnew(SSfi)=  max(0.0,  xold(SSfi) + dt2 * P)
!ORIG
!ORIG!-> Sea salt coarse
!ORIG
!ORIG      P =   &
!ORIG     rcss(QSSCO,K) 
!ORIG      ! L = 0.0
!ORIG       xnew(SSco)=  max(0.0,  xold(SSco) + dt2 * P)
!ORIG   endif
!ORIG
!ORIG!-> AMSU 
!ORIG      ! P = 0.0
!ORIG 
!ORIG      ! L = 0.0
!ORIG 
!ORIG
!ORIG!       xnew(AMSU)=  max(0.0, xold(AMSU))
!ORIG
!ORIG
!ORIG!-> AMNI 
!ORIG      ! P = 0.0
!ORIG 
!ORIG      ! L = 0.0
!ORIG 
!ORIG
!ORIG!       xnew(AMNI)=  max(0.0, xold(AMNI))
!ORIG
