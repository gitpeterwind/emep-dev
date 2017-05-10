c-----------------------------------------------------------------------
c
c M1QN3 optimizer for EMEP 3D-var.
c
c Return status
c   The integer 'reverse' value is the first that should be check
c   in the calling program. Set to values other than '1' or '4'
c   to indicate an error, in the default code often '-1'.
c
c History
c   2015-02, Arjo Segers
c     Extended original code for use in collective call with
c     state vectors decomposed over processors.
c     Two pairs of vector sizes should be passed as arguments now:
c       n_glb, n_loc
c         The global (former 'n') and local (new) state vector size.
c         Global number is used for information and computation
c         of some thresholds. Local number defines the size of
c         state vector and gradient, and is also used for internal
c         storage of derived vectors.
c       ndz_glb, ndz_loc
c         The global (former 'ndz') and local (new) work array size.
c         Global number is used for computation of thresholds,
c         local number defines the actual array size.
c     Apart from introduction of local sizes, other changes were
c     necessary too:
c       - add collective calls for threshold computations,
c         for exampleto compute a maximum over all processors;
c       - moved collective calls (norms) out of 'impres' tests,
c         since processors might have a different logging level.
c       
c-----------------------------------------------------------------------

      subroutine m1qn3 (simul,prosca,ctonb,ctcab,
     &                  n_glb,n_loc,x,f,g,dxmin,df1,
     &                  epsg,normtype,impres,io,imode,omode,niter,nsim,
     &                  iz,dz,ndz_glb,ndz_loc,reverse,indic,izs,rzs,dzs)
c
c-----------------------------------------------------------------------
c
c     M1QN3, Version 3.3, October 2009
c
c     M1qn3 has two running modes: the SID (Scalar Initial Scaling) mode
c     and the DIS (Diagonal Initial Scaling) mode. Both do not require
c     the same amount of storage, the same subroutines, ...
c     In the description below, items that differ in the DIS mode with
c     respect to the SIS mode are given in brakets.
c
c     Use the following subroutines:
c         M1QN3A
c         DDD, DDDS
c         NLIS0 + DCUBE (Dec 88)
c         MUPDTS, DYSTBL.
c
c     The following routines are proposed to the user in case the
c     Euclidean scalar product is used:
c         DUCLID, DTONBE, DTCABE.
c
c     La sous-routine M1QN3 est une interface entre le programme
c     appelant et la sous-routine M1QN3A, le minimiseur proprement dit.
c
c     Le module PROSCA est sense realiser le produit scalaire de deux
c     vecteurs de Rn; le module DTONB est sense realiser le changement
c     de coordonnees correspondant au changement de bases: base
c     euclidienne -> base orthonormale (pour le produit scalaire
c     PROSCA); le module CTBAB fait la transformation inverse: base
c     orthonormale -> base euclidienne.
c
c     Iz is an integer working zone for M1QN3A, its dimension is 5.
c     It is formed of 5 scalars that are set by the optimizer:
c         - the dimension of the problem,
c         - an identifier of the scaling mode,
c         - the number of updates,
c         - two pointers.
c
c     Dz est la zone de travail pour M1QN3A, de dimension ndz.
c     Elle est subdivisee en
c         3 [ou 4] vecteurs de dimension n: d,gg,[diag,]aux
c         m scalaires: alpha
c         m vecteurs de dimension n: ybar
c         m vecteurs de dimension n: sbar
c
c     m est alors le plus grand entier tel que
c         m*(2*n+1)+3*n .le. ndz [m*(2*n+1)+4*n .le. ndz)]
c     soit m := (ndz-3*n) / (2*n+1) [m := (ndz-4*n) / (2*n+1)].
c     Il faut avoir m >= 1, donc ndz >= 5n+1 [ndz >= 6n+1].
c
c     A chaque iteration la metrique est formee a partir d un multiple
c     de l'identite [d'une matrice diagonale] D qui est mise a jour m
c     fois par la formule de BFGS en utilisant les m couples {y,s} les
c     plus recents.
c
c-----------------------------------------------------------------------
c
c     Authors: Jean Charles Gilbert, Claude Lemarechal, INRIA.
c
c     Copyright 2008, 2009, INRIA.
c
c     M1QN3 is distributed under the terms of the GNU General Public
c     License.
c
c-----------------------------------------------------------------------
c
c     This program is free software: you can redistribute it and/or
c     modify it under the terms of the GNU General Public License as
c     published by the Free Software Foundation, either version 3 of
c     the License, or (at your option) any later version.
c
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with this program.  If not, see
c     <http://www.gnu.org/licenses/>.
c
c-----------------------------------------------------------------------
c
      use MPIF90, only : MPI_COMM_WORLD
      use MPIF90, only : MPIF90_AllReduce, MPI_SUM
c
      implicit none
c
c     arguments
c
      external  :: simul, prosca, ctonb, ctcab

      ! problem size:
      integer, intent(in)               ::  n_glb
      integer, intent(in)               ::  n_loc
      ! current esitmate, next estimate
      double precision, intent(inout)   ::  x(n_loc)
      ! evaluation at x
      double precision, intent(inout)   ::  f
      ! gradient at x
      double precision, intent(inout)   ::  g(n_loc)
      ! thresholds:
      double precision, intent(in)      ::  dxmin
      double precision , intent(in)     ::  df1
      double precision, intent(inout)   ::  epsg
      ! which norm?
      character(len=3), intent(in)      ::  normtype
      ! messages:
      integer, intent(in)               ::  impres
      integer, intent(in)               ::  io
      ! input mode:
      !   (1) scaling modes: 0=DIS, 1=SIS
      !   (2) starting mode
      !   (3) when to call simulator
      integer, intent(in)               ::  imode(3)
      ! return status:
      integer, intent(out)              ::  omode
      ! counters:
      integer, intent(inout)            ::  niter
      integer, intent(inout)            ::  nsim
      ! settings:
      integer, intent(inout)            ::  iz(5)
      ! work array:
      integer, intent(in)               ::  ndz_glb
      integer, intent(in)               ::  ndz_loc
      double precision, intent(inout)   ::  dz(ndz_loc)
      ! communication mode:
      integer, intent(inout)            ::  reverse
      ! info:
      integer, intent(inout)            ::  indic
      ! user arrays:
      integer                           ::  izs(*)
      real                              ::  rzs(*)
      double precision                  ::  dzs(*)
c
c     parametres
c
      character(len=*), parameter  ::  rname = 'm1qn3/m1qn3'
c
c --- local variables
c
c     The variable 'reentry' is used to know where to jump when entering
c     the subroutine in reverse commnunication (it is better than using
c     'reverse', which is known outside m1qn3 and could be changed by
c     the user):
c     = 0: no jump, start at the begining of the subroutine
c     = 1: jump inside m1qn3a, where the simulator was called with
c          indic=1
c     = 2: jump inside m1qn3a/mlis3, where the simulator was called with
c          indic=4
c     The variable 'reentry' is set to 0 by the compiler and when m1qn3
c     has completed its job (see below); since the value is saved, it
c     has always the value 0 when entering the solver for the first
c     time when solving a new problem, even if it is in the same
c     program.
c
      logical inmemo,sscale
      integer ntravu,id,igg,idiag,iaux,ialpha,iybar,isbar,m,mmemo,
     &    reentry
      double precision gnorm,ps
      save inmemo,sscale,ntravu,id,igg,idiag,iaux,ialpha,iybar,isbar,m,
     &    mmemo,reentry,gnorm,ps
c
      integer               ::  status
      integer               ::  n_sum
c
c     data
c
      data reentry /0/
c
c --- function
c
      double precision ddot,dnrmi
c
c --- stop if reverse < 0 (m1qn3 should not be called with reverse < 0)
c
      if (reverse.lt.0) then
          write (io,'(/a,a,i0,a/)')
     &        " >>> m1qn3 should not be called with a negative reverse",
     &        " (=", reverse, ")"
          reverse=-1; return
      end if
c
c --- possible jumps
c     9999: go directly in m1qn3a
c
      if (reentry.gt.0) goto 9999
c
c---- license notice
c
      if (impres.ge.5) then
          write (io,934)
      end if
  934 format (1x,79("-")
     &    /1x,"M1QN3 Copyright (C) 2008, J. Ch. Gilbert, Cl. ",
     &        "Lemarechal."
     &    /1x,79("-")
     &    /1x,"This program comes with ABSOLUTELY NO WARRANTY. This is",
     &        " free software, and you"
     &    /1x,"are welcome to redistribute it under certain ",
     &        "conditions. See the file COPYING "
     &    /1x,"in the root directory of the M1QN3 distribution for ",
     &        "details."
     &    /1x,79("-"))
c
c---- impressions initiales et controle des arguments
c
      if (impres.ge.1) then
          write (io,900) n_glb,dxmin,df1,epsg,normtype,niter,nsim,impres
          if (reverse.gt.0) then
              write (io,'(5x,a)') "reverse communication"
          else
              write (io,'(5x,a)') "direct communication"
          end if
      end if
  900 format (/" M1QN3 (Version 3.3, October 2009): entry point"/
     &    5x,"dimension of the problem (n):",i14/
     &    5x,"absolute precision on x (dxmin):",9x,1pd9.2/
     &    5x,"expected decrease for f (df1):",11x,1pd9.2/
     &    5x,"relative precision on g (epsg):",10x,1pd9.2,
     &       " (",a3,"-norm)"/
     &    5x,"maximal number of iterations (niter):",i6/
     &    5x,"maximal number of simulations (nsim):",i6/
     &    5x,"printing level (impres):",15x,i4)
c
c---- check the arguments
c
      if (n_glb.le.0) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,'(/a)')
     &        " >>> m1qn3: n should be > 0"
          return
      end if
      if (niter.le.0) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,'(/a)')
     &        " >>> m1qn3: niter should be > 0"
          return
      end if
      if (nsim.le.0) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,'(/a)')
     &        " >>> m1qn3: nsim should be > 0"
          return
      end if
      if (dxmin.le.0.d0) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,'(/a)')
     &        " >>> m1qn3: dxmin should be > 0.d0"
          return
      end if
c     if (epsg.le.0.d0 .or. epsg.gt.1.d0) then
      if (epsg.le.0.d0) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,'(/a)')
     &        " >>> m1qn3: epsg should be > 0.d0"
          return
      end if
      if (epsg.ge.1.d0) then
          omode=1
          niter=0
          nsim=0
          epsg=1.d0
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,'(/a)')
     &        " >>> m1qn3: epsg is >= 1.d0, no need to make progress"
          goto 1000
      end if
      if ((normtype.ne.'two') .and.
     &    (normtype.ne.'sup') .and.
     &    (normtype.ne.'dfn')) then
          omode=2
          if (reverse.gt.0) reverse = -1
          write (io,'(/a,a,a/)') " >>> m1qn3: unknown norm type '",
     &        normtype, "'"
          return
      end if
      if (impres.lt.0) then
          omode=2
          if (reverse.gt.0) reverse = -1
          write (io,'(/a,i0/)')
     &        " >>> m1qn3: impres should be >= 0 and has the value ",
     &        impres
          return
      end if
      
      ! check consistency of state vector decomposition;
      ! the work array cannot be checked in this way
      ! since it is not decomposed but simply smaller:
      call MPIF90_AllReduce( n_loc, n_sum, MPI_SUM, 
     &                          MPI_COMM_WORLD, status )
      if ( status /= 0 ) then
        write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
        reverse=-1; return
      end if
      if ( n_sum /= n_glb ) then
        write (io,'("sum of local problem sizes is ",i0,/
     &              "while global size is ",i0," (local ",i0,")")')
     &                 n_sum, n_glb, n_loc
        write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
        reverse=-1; return
      end if
c
c---- what method
c
      if (imode(1).eq.0) then
          if (impres.ge.1) write (io,920)
  920     format (/" m1qn3: Diagonal Initial Scaling mode")
          sscale=.false.
      else
          if (impres.ge.1) write (io,921)
  921     format (/" m1qn3: Scalar Initial Scaling mode")
          sscale=.true.
      end if
c
      if ( (ndz_glb.lt.5*n_glb+1) .or. 
     &     ((.not.sscale).and.(ndz_glb.lt.6*n_glb+1)) ) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,922)
  922     format (/" >>> m1qn3: not enough memory allocated")
          return
      end if
c
c---- Compute m
c
      call mupdts (sscale,inmemo,n_glb,m,ndz_glb)
c
c     --- Check the value of m (if (y,s) pairs in core, m will be >= 1)
c
      if (m.lt.1) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,930)
  930     format (/" >>> m1qn3: m is set too small in mupdts")
          return
      end if
c
c     --- mmemo = number of (y,s) pairs in core memory
c
      mmemo=1
      if (inmemo) mmemo=m
c
      ntravu=2*(2+mmemo)*n_loc+m
      if (sscale) ntravu=ntravu-n_loc
      if (impres.ge.1) write (io,931) ndz_glb,ntravu,m
  931 format (/5x,"allocated memory (ndz) :",i9/
     &         5x,"used memory :           ",i9/
     &         5x,"number of updates :     ",i9)
      if (ndz_loc.lt.ntravu) then
          omode=2
          if (reverse.gt.0) reverse = -1
          if (impres.ge.1) write (io,922)
          return
      end if
c
      if (impres.ge.1) then
          if (inmemo) then
              write (io,932)
          else
              write (io,933)
          end if
      end if
  932 format (5x,"(y,s) pairs are stored in core memory")
  933 format (5x,"(y,s) pairs are stored by the user")
c
c---- cold start or warm restart ?
c     check iz: iz(1)=n, iz(2)=(0 if DIS, 1 if SIS),
c               iz(3)=m, iz(4)=jmin, iz(5)=jmax
c
      if (imode(2).eq.0) then
          if (impres.ge.1) write (io,940)
      else
          if ( iz(1) .ne. n_glb    .or.
     &         iz(2) .ne. imode(1) .or.
     &         iz(3) .ne. m        .or. 
     &         iz(4) .lt. 1        .or. 
     &         iz(5) .lt. 0        .or.
     &         iz(4) .gt. iz(3)    .or.
     &         iz(5) .gt. iz(3)          ) then
              omode=2
              if (reverse.gt.0) reverse = -1
              if (impres.ge.1) then
                  write (io,941)
                  if (iz(1).ne.n_glb) write (io,942)
                  if (iz(2).ne.imode(1)) write (io,943)
                  if (iz(3).ne.m) write (io,944)
                  if (iz(4).lt.1 .or. iz(5).lt.0 .or. iz(4).gt.iz(3)
     &                .or. iz(5).gt.iz(3)) write (io,945)
              end if
              return
          end if
          if (impres.ge.1) write (io,946)
      end if
  940 format (/" m1qn3: cold start"/1x)
  941 format (/" >>> m1qn3: inconsistent warm restart ")
  942 format (" >>> m1qn3: (the number of variables has changed)")
  943 format (" >>> m1qn3: (the scaling mode has changed)")
  944 format (" >>> m1qn3: (the number of updates has changed)")
  945 format (" >>> m1qn3: (wrong pointers)")
  946 format (/" m1qn3: warm restart"/1x)
      iz(1)=n_glb
      iz(2)=0
      if (sscale) iz(2)=1
      iz(3)=m
c
c---- split the working zone dz
c
      idiag  = 1
      iybar  = idiag + n_loc
      if (sscale) iybar=1
      isbar  = iybar + n_loc * mmemo
      id     = isbar + n_loc * mmemo
      igg    = id    + n_loc
      iaux   = igg   + n_loc
      ialpha = iaux  + n_loc
c
c---- call the optimization code
c
 9999 continue
      call m1qn3a (simul,prosca,ctonb,ctcab,
     &             n_glb, n_loc,x,f,g,dxmin,df1,epsg,
     &             normtype,impres,io,imode,omode,niter,nsim,
     &             inmemo,iz(3),iz(4),iz(5),
     &             dz(id),dz(igg),dz(idiag),dz(iaux),dz(ialpha),
     &             dz(iybar),dz(isbar),
     &             reverse,reentry,indic, izs,rzs,dzs)
      if ( reentry > 0 ) return
      if ( reentry < 0 ) then; reverse=-1; return; end if
c
c---- impressions finales
c
 1000 continue
      if (impres.ge.1) write (io,960) omode,niter,nsim,epsg
  960 format (/1x,79("-")/
     &        /" m1qn3: output mode is ",i2
     &        /5x,"number of iterations: ",i14
     &        /5x,"number of simulations: ",i13
     &        /5x,"realized relative precision on g: ",1pd9.2)
      if (normtype.eq.'two') then
          gnorm = sqrt(ddot(n_loc,g,1,g,1))
      elseif (normtype.eq.'sup') then
          gnorm = dnrmi(n_loc,g)
      elseif (normtype.eq.'dfn') then
          call prosca (n_loc,g,g,ps,izs,rzs,dzs)
          gnorm=dsqrt(ps)
      end if

      if (impres.ge.1) write (io,961) f,normtype,gnorm
  961 format (5x,"f             = ",1pd15.8
     &       /5x,a3,"-norm of g = ",1pd15.8)

      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine m1qn3a( simul, prosca, ctonb, ctcab,
     &                   n_glb, n, x, f, g, dxmin, df1,
     &                   epsg,normtype,impres,io,imode,omode,niter,nsim,
     &                   inmemo,m,jmin,jmax,
     &                   d,gg,diag,aux,alpha,ybar,sbar,
     &                   reverse,reentry,indic, izs,rzs,dzs)
c----
c
c     Code d optimisation proprement dit.
c
c     Return status: 
c       'reentry' is set to negative value in case of error
c
c----
c 
      use MPIF90, only : MPI_COMM_WORLD
      use MPIF90, only : MPIF90_AllReduce, MPI_MAX, MPI_SUM
c
      implicit none
c
c         arguments
c
      logical                         ::  inmemo
      character(len=3), intent(in)    ::  normtype
      integer, intent(in)             ::  n_glb
      integer, intent(in)             ::  n
      integer :: impres,io,imode(3),omode,niter,nsim,m,jmin,jmax,indic,
     &    reentry,izs(*)
      integer, intent(inout)  ::  reverse
      real                    ::  rzs(*)
      double precision x(n),f,g(n),dxmin,df1,epsg,d(n),gg(n),diag(n),
     &    aux(n),alpha(m),ybar(n,1),sbar(n,1),dzs(*)
      external simul,prosca,ctonb,ctcab
c
c         variables locales
c
      logical               ::  sscale,cold,warm,skip_update
      integer               ::  i,itmax,moderl,isim,jcour
      double precision      ::  d1
      double precision      ::  t
      double precision      ::  tmin, tmin_loc
      double precision      ::  tmax
      double precision      ::  gnorm
      double precision      ::  gnorms
      double precision      ::  eps1
      double precision      ::  ff
      double precision      ::  preco
      double precision      ::  precos
      double precision      ::  ys,den,dk,dk1,ps,ps2,hp0
      integer               ::  status
c
      save sscale,cold,warm,i,itmax,moderl,isim,jcour,d1,t,tmin,tmax,
     &    gnorm,gnorms,eps1,ff,preco,precos,ys,den,dk,dk1,ps,ps2,hp0
     
      double precision      ::  ps_loc, ps2_loc
c
c         parametres
c
      character(len=*), parameter  ::  rname = 'm1qn3/m1qn3a'
c    
      double precision rm1,rm2
      parameter (rm1=0.0001d+0,rm2=0.99d+0)
      double precision pi
      parameter (pi=3.1415927d+0)
      double precision rmin
c
c         function
c
      double precision ddot,dnrmi
c
c --- possible jumps
c     9998: call of the simulator in m1qn3a with indic = 1
c     9999: call of the simulator in mlis3 with indic = 4
c
      if (reentry.eq.1) goto 9998
      if (reentry.eq.2) goto 9999
c
c---- initialisation
c
      rmin=1.d-20
c
      sscale=.true.                      ! SIS mode
      if (imode(1).eq.0) sscale=.false.  ! DIS mode
c
      warm=.false.
      if (imode(2).eq.1) warm=.true.
      cold=.not.warm
c
      skip_update = .false.
c
      itmax=niter
      niter=0
      isim=1
      eps1=1.d+0
c
      call prosca (n,g,g,ps,izs,rzs,dzs)
      gnorm = dsqrt(ps)
      if (normtype.eq.'two') then
          gnorms = sqrt(ddot(n,g,1,g,1))
      elseif (normtype.eq.'sup') then
          gnorms = dnrmi(n,g)
      elseif (normtype.eq.'dfn') then
          gnorms = gnorm
      end if
      if (impres.ge.1) write (io,900) f,normtype,gnorms
  900 format (5x,"f             = ",1pd15.8
     &       /5x,a3,"-norm of g = ",1pd15.8)
      if (gnorms.lt.rmin) then
          omode=2
          if (impres.ge.1) write (io,901)
          goto 1000
      end if
  901 format (/" >>> m1qn3a: initial gradient is too small")
c
c     --- initialisation pour dd
c
      if (cold) then
          jmin=1
          jmax=0
      end if
      jcour=1
      if (inmemo) jcour=jmax
c
c     --- mise a l'echelle de la premiere direction de descente
c
      if (cold) then
c
c         --- use Fletcher's scaling and initialize diag to 1.
c
          precos=2.d+0*df1/gnorm**2
          do 10 i=1,n
              d(i)=-g(i)*precos
              diag(i)=1.d+0
   10     continue
          if (impres.ge.5) write(io,902) precos
  902     format (/" m1qn3a: descent direction -g: precon = ",d10.3)
      else
c
c         --- use the matrix stored in [diag and] the (y,s) pairs
c
          if (sscale) then
              call prosca (n,ybar(:,jcour),ybar(:,jcour),ps,izs,rzs,dzs)
              precos=1.d+0/ps
          end if
          do 11 i=1,n
              d(i)=-g(i)
  11      continue
          if (inmemo) then
              call dd (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax,
     &                 precos,diag,alpha,ybar,sbar,izs,rzs,dzs)
          else
              call dds (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax,
     &                  precos,diag,alpha,ybar,sbar,izs,rzs,dzs)
          end if
      end if
c
      if (impres.eq.3) write(io,903)
      if (impres.eq.4) write(io,903)
  903 format (/1x,79("-"))
  904 format (1x)
c
c     --- initialisation pour mlis3
c
      tmax=1.d+20
      call prosca (n,d,g,hp0,izs,rzs,dzs)
      if (hp0.ge.0.d+0) then
          omode=7
          if (impres.ge.1) write (io,905) niter,hp0
          goto 1000
      end if
  905 format (/" >>> m1qn3 (iteration ",i2,"): "
     &        /5x," the search direction d is not a ",
     &         "descent direction: (g,d) = ",d12.5)
c
c     --- compute the angle (-g,d)
c
      if (warm) then
          ! l2 norm of g:
          call prosca (n,g,g,ps,izs,rzs,dzs)
          ps=dsqrt(ps)
          ! l2 norm of d:
          call prosca (n,d,d,ps2,izs,rzs,dzs)
          ps2=dsqrt(ps2)
          ! compute direction:
          ps=hp0/ps/ps2
          ps=dmin1(-ps,1.d+0)
          ps=dacos(ps)
          d1=ps*180.d+0/pi
          ! info:
          if (impres.ge.5) write (io,906) sngl(d1)
      end if
  906 format (/" m1qn3: descent direction d: ",
     &        "angle(-g,d) = ",f5.1," degrees")
c
c---- Debut de l iteration. on cherche x(k+1) de la forme x(k) + t*d,
c     avec t > 0. On connait d.
c
c         Debut de la boucle: etiquette 100,
c         Sortie de la boucle: goto 1000.
c
100   niter=niter+1
      if (impres.ge.5) write(io,903)
      if (impres.ge.4) write(io,904)
      if (impres.ge.4) write (io,910) niter,isim,f,hp0
  910 format (" m1qn3: iter ",i0,", simul ",i0,
     &        ", f=",1pd15.8,", h'(0)=",d12.5)
c
c     --- free simulation if desired
c
      if (imode(3).gt.0) then
          if (mod(niter-1,imode(3)).eq.0) then
              indic=1
              if (reverse.gt.0) then
                  reentry = 1
                  return
              else
                  call simul(indic,n,x,f,g,izs,rzs,dzs)
              end if
          end if
      end if
 9998 continue
c
c     --- recherche lineaire et nouveau point x(k+1)
c
      do 101 i=1,n
          gg(i)=g(i)
101   continue
      ff=f
      if (impres.ge.5) write (io,911)
  911 format (/" m1qn3: line search")
c
c         --- calcul de tmin
c
      ! local maximum of |d|
      tmin_loc = 0.d+0
      do i = 1, n
        tmin_loc = max( tmin_loc, abs(d(i)) )
      end do
      ! global maximum:
      call MPIF90_AllReduce( tmin_loc, tmin, MPI_MAX, 
     &                          MPI_COMM_WORLD, status )
      if ( status /= 0 ) then
        write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
        reentry=-1; return
      end if

      ! reset using parameters:
      tmin = dxmin/tmin

      ! settings for mlis3:
      t = 1.d+0
      d1 = hp0
c
 9999 continue
      call mlis3(n,simul,prosca,x,f,d1,t,tmin,tmax,d,g,rm2,rm1,impres,
     &            io,moderl,isim,nsim,aux,reverse,reentry,indic,izs,rzs,
     &            dzs)
      if ( reentry > 0 ) return
      ! error ?
      if ( reentry < 0 ) return
c
c         --- mlis3 renvoie les nouvelles valeurs de x, f et g
c
      if (moderl.ne.0) then
          if (moderl.lt.0) then
c
c             --- calcul impossible
c                 t, g: ou les calculs sont impossibles
c                 x, f: ceux du t_gauche (donc f <= ff)
c
              omode=moderl
              goto 1000
          elseif (moderl.eq.1) then
c
c             --- descente bloquee sur tmax
c
              skip_update = .true.
c             omode=3
c             if (impres.ge.1) write(io,912) niter
c 912         format (/" >>> m1qn3 (iteration ",i0,
c    &                "): line search blocked on tmax: "/
c    &                " >>> possible reasons: bad scaling,",
c    &                " unbounded problem")
          elseif (moderl.eq.4) then
c
c             --- nsim atteint
c                 x, f: ceux du t_gauche (donc f <= ff)
c
              omode=5
              goto 1000
          elseif (moderl.eq.5) then
c
c             --- arret demande par l utilisateur (indic = 0)
c                 x, f: ceux en sortie du simulateur
c
              omode=0
              goto 1000
          elseif (moderl.eq.6) then
c
c             --- arret sur dxmin ou appel incoherent
c                 x, f: ceux du t_gauche (donc f <= ff)
c
              omode=6
              goto 1000
          end if
      else
          skip_update = .false.
      end if
c
c NOTE: stopping tests are now done after having updated the matrix, so
c that update information can be stored in case of a later warm restart
c
c     --- mise a jour de la matrice
c
      if (skip_update) then
          if (impres.ge.5) write(io,'(/a)')
     &        " m1qn3: matrix update is skipped"
      elseif (m.gt.0) then
c
c       --- mise a jour des pointeurs
c
        jmax=jmax+1
        if (jmax.gt.m) jmax=jmax-m
        if ((cold.and.niter.gt.m).or.(warm.and.jmin.eq.jmax)) then
            jmin=jmin+1
            if (jmin.gt.m) jmin=jmin-m
        end if
        if (inmemo) jcour=jmax
c
c       --- y, s et (y,s)
c
        do 400 i=1,n
            sbar(i,jcour)=t*d(i)
            ybar(i,jcour)=g(i)-gg(i)
400     continue

        ! norm of sbar(:,jcour), local might be empty:
        if ( n > 0 ) then
          call prosca(n,sbar(:,jcour),sbar(:,jcour),ps,izs,rzs,dzs)
        else
          call prosca(n,      (/0.0/),      (/0.0/),ps,izs,rzs,dzs)
        end if
        dk1 = dsqrt(ps)
        ! info ..
        if (impres.ge.5) then
          if (niter.gt.1) write (io,930) dk1/dk
  930     format(/" m1qn3: convergence rate, s(k)/s(k-1) = ",1pd12.5)
        end if
        ! store for next print:
        dk = dk1
        
        ! dot product between ybar(:,jcour) and sbar(:,jcour),
        ! local might be empty:
        if ( n > 0 ) then
          call prosca(n,ybar(:,jcour),sbar(:,jcour),ys,izs,rzs,dzs)
        else
          call prosca(n,      (/0.0/),      (/0.0/),ys,izs,rzs,dzs)
        end if
        if (ys.le.0.d+0) then
            omode=7
            if (impres.ge.1) write (io,931) niter,ys
  931       format (/" >>> m1qn3 (iteration ",i2,
     &              "): the scalar product (y,s) = ",d12.5
     &              /27x,"is not positive")
            goto 1000
        end if
c
c       --- ybar et sbar
c
        d1=dsqrt(1.d+0/ys)
        do 410 i=1,n
            sbar(i,jcour)=d1*sbar(i,jcour)
            ybar(i,jcour)=d1*ybar(i,jcour)
  410   continue
        if (.not.inmemo) call ystbl (.true.,ybar,sbar,n,jmax)
c
c       --- compute the scalar or diagonal preconditioner
c
        if (impres.ge.5) write(io,932)
  932   format (/" m1qn3: matrix update:")
c
c           --- Here is the Oren-Spedicato factor, for scalar scaling
c
        if (sscale) then
          if ( n > 0 ) then
            call prosca(n,ybar(:,jcour),ybar(:,jcour),ps,izs,rzs,dzs)
          else
            call prosca(n,      (/0.0/),      (/0.0/),ps,izs,rzs,dzs)
          end if
            
          precos=1.d+0/ps
c
          if (impres.ge.5) write (io,933) precos
  933     format (5x,"Oren-Spedicato factor = ",d10.3)
c
c         --- Scale the diagonal to Rayleigh s ellipsoid.
c             Initially (niter.eq.1) and for a cold start, this is
c             equivalent to an Oren-Spedicato scaling of the
c             identity matrix.
c
        else
          if ( n > 0 ) then
            call ctonb(n,ybar(:,jcour),aux,izs,rzs,dzs)
          else
            call ctonb(n,      (/0.0/),aux,izs,rzs,dzs)
          end if
          ! local sum:
          ps_loc = 0.0d0
          do i=1,n
            ps_loc = ps_loc + diag(i) * aux(i) * aux(i)
          end do
          ! global sum:
          call MPIF90_AllReduce( ps_loc, ps, MPI_SUM, 
     &                            MPI_COMM_WORLD, status )
          if ( status /= 0 ) then
            write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
            reentry=-1; return
          end if
          ! scale factor:
          d1 = 1.0d0 / ps
          ! info ..
          if (impres.ge.5) then
              write (io,934) d1
  934         format(5x,"fitting the ellipsoid: factor = ",1pd10.3)
          end if
          ! scale:
          do i=1,n
            diag(i)=diag(i)*d1
          end do
c
c         --- update the diagonal
c             (gg is used as an auxiliary vector)
c
          if ( n > 0 ) then
            call ctonb(n,sbar(:,jcour),gg,izs,rzs,dzs)
          else
            call ctonb(n,      (/0.0/),gg,izs,rzs,dzs)
          end if
          ! local sum:
          ps_loc = 0.0d0
          do i = 1, n
            ps_loc = ps_loc + gg(i) * gg(i) / diag(i)
          end do
          ! global sum:
          call MPIF90_AllReduce( ps_loc, ps, MPI_SUM, 
     &                            MPI_COMM_WORLD, status )
          if ( status /= 0 ) then
            write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
            reentry=-1; return
          end if
          ! store:
          den = ps
          do 431 i=1,n
              diag(i)=1.d0/
     &               (1.d0/diag(i)+aux(i)**2-(gg(i)/diag(i))**2/den)
              if (diag(i).le.0.d0) then
                  if (impres.ge.5) write (io,935) i,diag(i),rmin
                  diag(i)=rmin
              end if
  431     continue
  935     format (/" >>> m1qn3-WARNING: diagonal element ",i8,
     &             " is negative (",d10.3,"), reset to ",d10.3)
c
          ! >>> variables for info, 
          !     need to be computed by all procs ...

          ! local sum:
          ps_loc = 0.0d0
          do i = 1, n
            ps_loc = ps_loc + diag(i)
          end do
          ! global sum:
          call MPIF90_AllReduce( ps_loc, ps, MPI_SUM, 
     &                            MPI_COMM_WORLD, status )
          if ( status /= 0 ) then
            write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
            reentry=-1; return
          end if
          ! mean value:
          ps = ps / n_glb
          ! store:
          preco=ps
c
          ! local sum:
          ps2_loc = 0.d0
          do i = 1, n
            ps2_loc = ps2_loc + (diag(i)-ps)**2
          end do
          ! global sum:
          call MPIF90_AllReduce( ps2_loc, ps2, MPI_SUM, 
     &                            MPI_COMM_WORLD, status )
          if ( status /= 0 ) then
            write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
            reentry=-1; return
          end if
          ! std.dev.:
          ps2 = dsqrt(ps2/n)

          ! info ?
          if (impres.ge.5) then
           write (io,936) preco,ps2
  936      format (5x,"updated diagonal: average value = ",
     &             1pd10.3,", sqrt(variance) = ",d10.3)
          end if  ! impres >= 5
          
          ! <<<

        end if  ! SIS or DIS mode

      end if  ! matrix update ?
c
c     --- printings
c
c
c     --- tests d arret
c
      call prosca(n,g,g,ps,izs,rzs,dzs)
      if (normtype.eq.'two') then
          gnorm = sqrt(ddot(n,g,1,g,1))
      elseif (normtype.eq.'sup') then
          gnorm = dnrmi(n,g)
      elseif (normtype.eq.'dfn') then
          gnorm = dsqrt(ps)
      end if
      eps1 = gnorm/gnorms
c
      if (impres.eq.3) then
          if (mod(niter-1,50).eq.0) write(io,'(/a,a)')
     &        "  iter  simul  stepsize            f                |g|",
     &        "       |g|/|g0|"
          write(io,
     &        '(1x,i5,2x,i5,2x,1pd9.2,2x,d21.14,2x,d12.5,2x,d11.4)')
     &        niter, isim, t, f, gnorm, eps1
      end if
      if (impres.ge.5) write (io,940) eps1
  940 format (/" m1qn3: stopping criterion on g: ",1pd12.5)
      if (eps1.lt.epsg) then
          omode=1
          goto 1000
      end if
      if (niter.eq.itmax) then
          omode=4
          if (impres.ge.1) write (io,941) niter
  941     format (/" >>> m1qn3 (iteration ",i0,
     &            "): maximal number of iterations")
          goto 1000
      end if
      if (isim.gt.nsim) then
          omode=5
          if (impres.ge.1) write (io,942) niter,isim
  942     format (/" >>> m1qn3 (iteration ",i3,"): ",i6,
     &            " simulations (maximal number reached)")
          goto 1000
      end if
c
c     --- calcul de la nouvelle direction de descente d = - H.g
c
      if (m.eq.0) then
          preco=2.d0*(ff-f)/ps
          do 500 i=1,n
              d(i)=-g(i)*preco
  500     continue
      else
          do 510 i=1,n
              d(i)=-g(i)
  510     continue
          if (inmemo) then
              call dd (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax,
     &                  precos,diag,alpha,ybar,sbar,izs,rzs,dzs)
          else
              call dds (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax,
     &                   precos,diag,alpha,ybar,sbar,izs,rzs,dzs)
          end if
      end if
c
c         --- test: la direction d est-elle de descente ?
c             hp0 sera utilise par mlis3
c
      call prosca (n,d,g,hp0,izs,rzs,dzs)
      if (hp0.ge.0.d+0) then
          omode=7
          if (impres.ge.1) write (io,905) niter,hp0
          goto 1000
      end if
      ! l2 norm of g:
      call prosca (n,g,g,ps,izs,rzs,dzs)
      ps=dsqrt(ps)
      ! l2 norm of d:
      call prosca (n,d,d,ps2,izs,rzs,dzs)
      ps2=dsqrt(ps2)
      ! compute direction:
      ps=hp0/ps/ps2
      ps=dmin1(-ps,1.d+0)
      ps=dacos(ps)
      d1=ps
      d1=d1*180.d0/pi
      ! info ...
      if (impres.ge.5) write (io,906) sngl(d1)
c
c---- on poursuit les iterations
c
      goto 100
c
c --- n1qn3 has finished for ever
c
 1000 continue
      if (reverse.ne.0) reverse = -1
      reentry = 0
      nsim=isim
      epsg=eps1
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine dd (prosca,ctonb,ctcab,n,sscale,nm,depl,aux,jmin,jmax,
     &                precos,diag,alpha,ybar,sbar,izs,rzs,dzs)
c----
c
c     calcule le produit H.g ou
c         . H est une matrice construite par la formule de bfgs inverse
c           a nm memoires a partir de la matrice diagonale diag
c           dans un espace hilbertien dont le produit scalaire
c           est donne par prosca
c           (cf. J. Nocedal, Math. of Comp. 35/151 (1980) 773-782)
c         . g est un vecteur de dimension n (en general le gradient)
c
c     la matrice diag apparait donc comme un preconditionneur diagonal
c
c     depl = g (en entree), = H g (en sortie)
c
c     la matrice H est memorisee par les vecteurs des tableaux
c     ybar, sbar et les pointeurs jmin, jmax
c
c     alpha(nm) est une zone de travail
c
c     izs(*),rzs(*),dzs(*) sont des zones de travail pour prosca
c
c----
c
      implicit none
c
c         arguments
c
      logical sscale
      integer n,nm,jmin,jmax,izs(*)
      real rzs(*)
      double precision depl(n),precos,diag(n),alpha(nm),ybar(n,1),
     &    sbar(n,1),aux(n),dzs(*)
      external prosca,ctonb,ctcab
c
c         variables locales
c
      integer jfin,i,j,jp
      double precision r,ps
c
      jfin=jmax
      if (jfin.lt.jmin) jfin=jmax+nm
c
c         phase de descente
c
      do 100 j=jfin,jmin,-1
          jp=j
          if (jp.gt.nm) jp=jp-nm
          if ( n > 0 ) then
            call prosca (n,depl,sbar(:,jp),ps,izs,rzs,dzs)
          else
            call prosca (n,depl,   (/0.0/),ps,izs,rzs,dzs)
          end if
          r=ps
          alpha(jp)=r
          do 20 i=1,n
              depl(i)=depl(i)-r*ybar(i,jp)
20        continue
100   continue
c
c         preconditionnement
c
      if (sscale) then
          do 150 i=1,n
              depl(i)=depl(i)*precos
  150     continue
      else
          call ctonb(n,depl,aux,izs,rzs,dzs)
          do 151 i=1,n
              aux(i)=aux(i)*diag(i)
  151     continue
          call ctcab(n,aux,depl,izs,rzs,dzs)
      end if
c
c         remontee
c
      do 200 j=jmin,jfin
          jp=j
          if (jp.gt.nm) jp=jp-nm
          if ( n > 0 ) then
            call prosca(n,depl,ybar(:,jp),ps,izs,rzs,dzs)
          else
            call prosca(n,depl,   (/0.0/),ps,izs,rzs,dzs)
          end if
          r=alpha(jp)-ps
          do 120 i=1,n
              depl(i)=depl(i)+r*sbar(i,jp)
120       continue
200   continue
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine dds (prosca,ctonb,ctcab,n,sscale,nm,depl,aux,jmin,
     &                 jmax,precos,diag,alpha,ybar,sbar,izs,rzs,dzs)
c----
c
c     This subroutine has the same role as dd (computation of the
c     product H.g). It supposes however that the (y,s) pairs are not
c     stored in core memory, but on a devise chosen by the user.
c     The access to this devise is performed via the subroutine ystbl.
c
c----
c
      implicit none
c
c         arguments
c
      logical sscale
      integer n,nm,jmin,jmax,izs(*)
      real rzs(*)
      double precision depl(n),precos,diag(n),alpha(nm),ybar(n),sbar(n),
     &    aux(n),dzs(*)
      external prosca,ctonb,ctcab
c
c         variables locales
c
      integer jfin,i,j,jp
      double precision r,ps
c
      jfin=jmax
      if (jfin.lt.jmin) jfin=jmax+nm
c
c         phase de descente
c
      do 100 j=jfin,jmin,-1
          jp=j
          if (jp.gt.nm) jp=jp-nm
          call ystbl (.false.,ybar,sbar,n,jp)
          call prosca (n,depl,sbar,ps,izs,rzs,dzs)
          r=ps
          alpha(jp)=r
          do 20 i=1,n
              depl(i)=depl(i)-r*ybar(i)
20        continue
100   continue
c
c         preconditionnement
c
      if (sscale) then
          do 150 i=1,n
              depl(i)=depl(i)*precos
  150     continue
      else
          call ctonb (n,depl,aux,izs,rzs,dzs)
          do 151 i=1,n
              aux(i)=aux(i)*diag(i)
  151     continue
          call ctcab (n,aux,depl,izs,rzs,dzs)
      end if
c
c         remontee
c
      do 200 j=jmin,jfin
          jp=j
          if (jp.gt.nm) jp=jp-nm
          call ystbl (.false.,ybar,sbar,n,jp)
          call prosca (n,depl,ybar,ps,izs,rzs,dzs)
          r=alpha(jp)-ps
          do 120 i=1,n
              depl(i)=depl(i)+r*sbar(i)
120       continue
200   continue
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine mlis3 (nloc,simul,prosca,x,f,fpn,t,tmin,tmax,d,g,
     1                  amd,amf,imp,io,logic,nap,napmax,xn,
     1                  reverse,reentry,indic,izs,rzs,dzs)
c
c ----
c
c     mlis3 + minuscules + commentaires
c     + version amelioree (XII 88): interpolation cubique systematique
c       et anti-overflows
c     + declaration variables (II/89, JCG).
c     + barr is also progressively decreased (12/93, CL & JChG).
c       barmul is set to 5.
c
c     ----------------------------------------------------------------
c
c        en sortie logic =
c
c        0          descente serieuse
c        1          descente bloquee sur tmax
c        4          nap > napmax
c        5          arret demande par le simulateur
c        6          fonction et gradient pas d accord
c        < 0        contrainte implicite active
c
c     indic on entry (in reverse communication)
c
c       < 0: the simulator cannot compute f and g
c       = 0: the simulator wants to stop
c       > 0: the simulator has done its work
c
c     indic on return (in reverse communication)
c
c       = 4: the simulator is asked to compute f and g
c
c     reverse
c
c       = 0: direct communication
c       = 1: reverse communication
c
c     reentry on entry
c
c       = 2: reverse communication, return from a simulation, skip to
c            the place where the simulator was called (9999)
c
c     reentry return
c
c       < 0: some error
c       = 0: reverse communication, the linesearch has finished its job
c       = 2: reverse communication, a simulation is required
c
c ----
c
      use MPIF90, only : MPI_COMM_WORLD
      use MPIF90, only : MPIF90_AllReduce, MPI_LOR
c
      implicit none
c
c --- arguments
c
      integer, intent(in)   ::  nloc

      external simul,prosca

      double precision      ::  x(nloc)
      double precision      ::  f
      double precision      ::  fpn
      double precision      ::  t
      double precision      ::  tmin
      double precision      ::  tmax
      double precision      ::  d(nloc)
      double precision      ::  g(nloc)
      double precision      ::  amd
      double precision      ::  amf
      integer               ::  imp
      integer               ::  io
      integer               ::  logic
      integer               ::  nap
      integer               ::  napmax
      double precision      ::  xn(nloc)
      integer               ::  reverse
      integer               ::  reentry
      integer               ::  indic
      integer               ::  izs(*)
      real                  ::  rzs(*)
      double precision      ::  dzs(*)
c
c     parametres
c
      character(len=*), parameter  ::  rname = 'm1qn3/mlis3'
c
c --- variables locales
c
      logical t_increased
      integer i,indica,indicd
      double precision tesf,tesd,tg,fg,fpg,td,ta,fa,fpa,d2,fn,fp,ffn,fd,
     &     fpd,z,test,barmin,barmul,barmax,barr,gauche,droite,taa,ps
      save t_increased,i,indica,indicd,tesf,tesd,tg,fg,fpg,td,ta,fa,
     &    fpa,d2,fn,fp,ffn,fd,fpd,z,test,barmin,barmul,barmax,barr,
     &    gauche,droite,taa,ps
c
      logical               ::  changed
      logical               ::  any_changed
      integer               ::  status
c
c     begin
c
 1000 format (/4x," mlis3   ",4x,"fpn=",1pd10.3," d2=",d9.2,
     1 "  tmin=",d9.2," tmax=",d9.2)
 1001 format (/4x," mlis3",3x,"stop on tmin",5x,
     1   "stepsizes",11x,"functions",8x,"derivatives")
 1002 format (4x," mlis3",37x,1pd10.3,2d11.3)
 1003 format (4x," mlis3",1pd14.3,2d11.3)
 1004 format (4x," mlis3",37x,1pd10.3," indic=",i3)
 1005 format (4x," mlis3",14x,1pd18.8,2x,d21.14,2x,d11.4)
 1006 format (4x," mlis3",14x,1pd18.8,"      indic=",i3)
 1007 format (/4x," mlis3",10x,"tmin forced to tmax")
 1008 format (/4x," mlis3",10x,"inconsistent call")
c
c --- possible jump
c
      if (reentry.eq.2) goto 9999
c
c     - check if arguments are ok:
c      if (nloc.gt.0 .and. fpn.lt.0.d0 .and. t.gt.0.d0
c     1 .and. tmax.gt.0.d0 .and. amf.gt.0.d0
c     1 .and. amd.gt.amf .and. amd.lt.1.d0) go to 5
c     - new check: allow nloc==0, collective call over all procs
      if ( (nloc .ge. 0    ) .and.
     &     (fpn  .lt. 0.0d0) .and. 
     &     (t    .gt. 0.0d0) .and.
     &     (tmax .gt. 0.0d0) .and. 
     &     (amf  .gt. 0.0d0) .and.
     &     (amd  .gt. amf  ) .and.
     &     (amd  .lt. 1.0d0)       ) then
        ! input ok
        go to 5
      end if
      ! input not ok, leave with error status:
      logic = 6
      go to 999

      ! begin
    5 tesf=amf*fpn
      tesd=amd*fpn
      barmin=0.01d0
      barmul=5.d0
      barmax=0.3d0
      barr=barmin
      td=0.d0
      tg=0.d0
      fn=f
      fg=fn
      fpg=fpn
      ta=0.d0
      fa=fn
      fpa=fpn
      call prosca(nloc,d,d,ps,izs,rzs,dzs)
      d2=ps
c
c               elimination d un t initial ridiculement petit
c
c<
c     if (t.gt.tmin) go to 20
c     t=tmin
c     if (t.le.tmax) go to 20
c     if (imp.gt.0) write (io,1007)
c     tmin=tmax
c  20 if (fn+t*fpn.lt.fn+0.9d0*t*fpn) go to 30
c     t=2.d0*t
c     go to 20
c changed into
      if (t.lt.tmin) then
          t=tmin
          if (imp.ge.4) write (io,'(a)')
     &        ' mlis3: initial step-size increased to tmin'
          if (t.gt.tmax) then
              if (imp.gt.0) write (io,1007)
              tmin=tmax
          end if
      end if
      t_increased = .false.
      do while (fn+t*fpn.ge.fn+0.9d0*t*fpn)
          t_increased = .true.
          t = 2.d0*t
      end do
      if (t_increased .and. (imp.ge.4)) write (io,'(a,1pd10.3)')
     &    ' mlis3: initial step-size increased to ',t
c>
   30 indica=1
      logic=0
      if (t.gt.tmax) then
          t=tmax
          logic=1
      end if
      if (imp.ge.4) write (io,1000) fpn,d2,tmin,tmax
c
c     --- nouveau x
c         initialize xn to the current iterate
c         use x as the trial iterate
c
      do i = 1, nloc
          xn(i) = x(i)
          x(i)  = xn(i) + t * d(i)
      end do
c
c --- boucle
c
  100 nap=nap+1
      if(nap.gt.napmax) then
          logic=4
          fn=fg
          do i = 1, nloc
            xn(i) = xn(i) + tg * d(i)
          end do
          go to 999
      end if
c
c     --- appel simulateur
c
      indic=4
      if (reverse.gt.0) then
          reentry = 2
          return
      else
          call simul(indic,nloc,x,f,g,izs,rzs,dzs)
      end if
 9999 continue
      if (indic.eq.0) then
c
c         --- arret demande par l utilisateur
c
          logic = 5
          fn = f
          do i = 1, nloc
            xn(i) = x(i)
          end do
          go to 999
      end if
      if(indic.lt.0) then
c
c         --- les calculs n ont pas pu etre effectues par le simulateur
c
          td=t
          indicd=indic
          logic=0
          if (imp.ge.4) write (io,1004) t,indic
          t=tg+0.1d0*(td-tg)
          go to 905
      end if
c
c     --- les tests elementaires sont faits, on y va
c
      call prosca(nloc,d,g,ps,izs,rzs,dzs)
      fp=ps
c
c     --- premier test de Wolfe
c
      ffn=f-fn
      if(ffn.gt.t*tesf) then
          td=t
          fd=f
          fpd=fp
          indicd=indic
          logic=0
          if(imp.ge.4) write (io,1002) t,ffn,fp
          go to 500
      end if
c
c     --- test 1 ok, donc deuxieme test de Wolfe
c
      if(imp.ge.4) write (io,1003) t,ffn,fp
      if(fp.gt.tesd) then
          logic=0
          go to 320
      end if
      if (logic.eq.0) go to 350
c
c     --- test 2 ok, donc pas serieux, on sort
c
  320 fn=f
      do i = 1, nloc
        xn(i) = x(i)
      end do
      go to 999
c
c
c
  350 tg=t
      fg=f
      fpg=fp
      if(td.ne.0.d0) go to 500
c
c              extrapolation
c
      taa=t
      gauche=(1.d0+barmin)*t
      droite=10.d0*t
      call ecube (t,f,fp,ta,fa,fpa,gauche,droite)
      ta=taa
      if(t.lt.tmax) go to 900
      logic=1
      t=tmax
      go to 900
c
c              interpolation
c
  500 if ( indica .le. 0 ) then
          ta=t
          t=0.9d0*tg+0.1d0*td
          go to 900
      end if
      test=barr*(td-tg)
      gauche=tg+test
      droite=td-test
      taa=t
      call ecube (t,f,fp,ta,fa,fpa,gauche,droite)
      ta=taa
      if (t.gt.gauche .and. t.lt.droite) then
          barr=dmax1(barmin,barr/barmul)
c         barr=barmin
        else
          barr=dmin1(barmul*barr,barmax)
      end if
c
c --- fin de boucle
c     - t peut etre bloque sur tmax
c       (venant de l extrapolation avec logic=1)
c
  900 fa=f
      fpa=fp
  905 indica=indic
c
c --- faut-il continuer ?
c
      if (td.eq.0.d0) go to 950
      if (td-tg.lt.tmin) go to 920
c
c     --- limite de precision machine (arret de secours) ?
c
      ! init flag:
      changed = .false.
      ! loop over local elements:
      do i = 1, nloc
        ! new value:
        z = xn(i) + t * d(i)
        !!old:
        !if ( z.ne.xn(i) .and. z.ne.x(i) ) go to 950
        ! different? then leave:
        changed = (z /= xn(i)) .and. (z /= x(i))
        if ( changed ) exit
      end do
      ! any change on one of the processes?
      ! collective call with logical-or:
      call MPIF90_AllReduce( changed, any_changed, MPI_LOR, 
     &                          MPI_COMM_WORLD, status )
      if ( status /= 0 ) then
        write (io,*) 'ERROR in ', __FILE__, ' line ', __LINE__
        reentry=-1; return
      end if
      ! changed value was found ? then jump as in old code:
      if ( any_changed ) goto 950
      
      ! info ...
      if (imp.gt.3) write (io,'(5x,a)') "mlis3   no change in x"
c
c --- arret sur dxmin ou de secours
c
  920 logic=6
c
c     si indicd<0, derniers calculs non faits par simul
c
      if (indicd.lt.0) logic=indicd
c
c     si tg=0, xn = xn_depart,
c     sinon on prend xn=x_gauche qui fait decroitre f
c
      if (tg.eq.0.d0) go to 940
      fn=fg
      do i = 1, nloc
        xn(i) = xn(i) + tg * d(i)
      end do
  940 if (imp.le.3) go to 999
      write (io,1001)
      write (io,1005) tg,fg,fpg
      if (logic.eq.6) write (io,1005) td,fd,fpd
      if (logic.eq.7) write (io,1006) td,indicd
      go to 999
c
c               recopiage de x et boucle
c
  950 do i = 1, nloc
        x(i) = xn(i) + t * d(i)
      end do
      go to 100
c
c     --- linesearch finished, no skip at next entry in mlis3
  999 if (reverse.ne.0) reentry = 0
c
      return
c
      end subroutine mlis3
c
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine ecube (t,f,fp,ta,fa,fpa,tlower,tupper)
c
      implicit none
c
c --- arguments
c
      double precision sign,den,anum,t,f,fp,ta,fa,fpa,tlower,tupper
c
c --- variables locales
c
      double precision z1,b,discri
c
c           Using f and fp at t and ta, computes new t by cubic formula
c           safeguarded inside [tlower,tupper].
c
      z1=fp+fpa-3.d0*(fa-f)/(ta-t)
      b=z1+fp
c
c              first compute the discriminant (without overflow)
c
      if (dabs(z1).le.1.d0) then
          discri=z1*z1-fp*fpa
        else
          discri=fp/z1
          discri=discri*fpa
          discri=z1-discri
          if (z1.ge.0.d0 .and. discri.ge.0.d0) then
              discri=dsqrt(z1)*dsqrt(discri)
              go to 120
          end if
          if (z1.le.0.d0 .and. discri.le.0.d0) then
              discri=dsqrt(-z1)*dsqrt(-discri)
              go to 120
          end if
          discri=-1.d0
      end if
      if (discri.lt.0.d0) then
          if (fp.lt.0.d0) t=tupper
          if (fp.ge.0.d0) t=tlower
          go to 900
      end if
c
c  discriminant nonnegative, compute solution (without overflow)
c
      discri=dsqrt(discri)
  120 if (t-ta.lt.0.d0) discri=-discri
      sign=(t-ta)/dabs(t-ta)
      if (b*sign.gt.0.d+0) then
          t=t+fp*(ta-t)/(b+discri)
        else
          den=z1+b+fpa
          anum=b-discri
          if (dabs((t-ta)*anum).lt.(tupper-tlower)*dabs(den)) then
              t=t+anum*(ta-t)/den
            else
              t=tupper
          end if
      end if
  900 t=dmax1(t,tlower)
      t=dmin1(t,tupper)
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine mupdts (sscale,inmemo,n,m,nrz)
c
      implicit none
c
c         arguments
c
      logical sscale,inmemo
      integer n,m,nrz
c----
c
c     On entry:
c       sscale: .true. if scalar initial scaling,
c               .false. if diagonal initial scaling
c       n:      number of variables
c
c     This routine has to return:
c       m:      the number of updates to form the approximate Hessien H,
c       inmemo: .true., if the vectors y and s used to form H are stored
c                  in core memory,
c               .false. otherwise (storage of y and s on disk, for
c                  instance).
c     When inmemo=.false., the routine "ystbl", which stores and
c     restores (y,s) pairs, has to be rewritten.
c
c----
c
      if (sscale) then
          m=(nrz-3*n)/(2*n+1)
      else
          m=(nrz-4*n)/(2*n+1)
      end if
      inmemo=.true.
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine ystbl (store,ybar,sbar,n,j)
c----
c
c     This subroutine should store (if store = .true.) or restore
c     (if store = .false.) a pair (ybar,sbar) at or from position
c     j in memory. Be sure to have 1 <= j <= m, where m in the number
c     of updates specified by subroutine mupdts.
c
c     The subroutine is used only when the (y,s) pairs are not
c     stored in core memory in the arrays ybar(.,.) and sbar(.,.).
c     In this case, the subroutine has to be written by the user.
c
c----
c
      implicit none
c
c         arguments
c
      logical store
      integer n,j
      double precision ybar(n),sbar(n)
c
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine ctonbe (n,u,v,izs,rzs,dzs)
c
      implicit none
      integer n,izs(*)
      real rzs(*)
      double precision u(1),v(1),dzs(*)
c
      integer i
c
      do i=1,n
          v(i)=u(i)
      end do
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine ctcabe (n,u,v,izs,rzs,dzs)
c
      implicit none
      integer n,izs(*)
      real rzs(*)
      double precision u(1),v(1),dzs(*)
c
      integer i
c
      do i=1,n
          v(i)=u(i)
      end do
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine euclid (n,x,y,ps,izs,rzs,dzs)
c
      implicit none
      integer n,izs(*)
      real rzs(*)
      double precision x(n),y(n),ps,dzs(*)
c
      integer i
c
      ps=0.d0
      do i=1,n
          ps=ps+x(i)*y(i)
      end do
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      subroutine simul_rc (indic,n,x,f,g,izs,rzs,dzs)
c
      implicit none
      integer indic,n,izs(*)
      real rzs(*)
      double precision x(n),f,g(n),dzs(*)
c
      return
      end
c
c--------0---------0---------0---------0---------0---------0---------0--
c
      double precision function dnrmi (n,v)
c
      integer n
      double precision v(n)
c
c----
c
c     Commputes the infinty-norm of the vector v(n)
c
c----
c
c --- local variables
c
      integer i
      double precision norm
c
c --- compute
c
      norm = 0.d0
      if (n.gt.0) then
        do i=1,n
           norm = max(norm,abs(v(i)))
        end do
      end if
      dnrmi = norm
c
      return
      end
