c-----------------------------------------------------------------------
      subroutine iniproc_f(intracomm)
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      logical flag

       call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
       if ( mpi_is_initialized .eq. 0 ) then
          call mpi_init (ierr)
       endif

      ! create communicator
      call init_nek_comm(intracomm)
      np  = np_
      nid = nid_

      if(nid.eq.0) call printHeader

      ! check upper tag size limit
      call mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,nval,flag,ierr)
      if (nval.lt.(10000+max(lp,lelg))) then
         if(nid.eq.0) write(6,*) 'ABORT: MPI_TAG_UB too small!'
         call exitt
      endif

      IF (NP.GT.LP) THEN
         WRITE(6,*) 
     $   'ERROR: Code compiled for a max of',LP,' processors.'
         WRITE(6,*) 
     $   'Recompile with LP =',NP,' or run with fewer processors.'
         WRITE(6,*) 
     $   'Aborting in routine INIPROC.'
         call exitt
      endif

      ! set word size for REAL
      nekreal = mpi_real

      ifdblas = .false.

      ! set word size for INTEGER
      ! HARDCODED since there is no secure way to detect an int overflow
      isize = 4

      ! set word size for LOGICAL
      lsize = 4

      ! set word size for CHARACTER
      csize = 1
c
      PID = 0
      NULLPID=0
      NODE0=0
      NODE= NID+1

      if (nid.eq.0) then 
         write(6,*) 'Number of processors:',np
         WRITE(6,*) 'REAL    wdsize      :',WDSIZE
         WRITE(6,*) 'INTEGER wdsize      :',ISIZE
      endif

      call crystal_setup(cr_h,nekcomm,np)  ! set cr handle to new instance

      return
      end
c-----------------------------------------------------------------------
      subroutine init_nek_comm_f(intracomm)
      include 'mpif.h'
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
C
      call create_comm(intracomm) ! set up nekton specific communicator
c
      nid_  = mynode()
      np_   = numnodes()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gop_f( x, w, op, n)
c
c     Global vector commutative operation in single precision
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
      real*4 x(n), w(n)
      character*3 op
c
      call adelay
      if (op.eq.'+  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_prod,nekcomm,ierr)
      else
         write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
         call exitt
      endif

      call copy_f(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine igop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      integer x(n), w(n)
      character*3 op

      call adelay
      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
      endif

      call icopy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine i8gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      integer*8 x(n), w(n)
      character*3 op

      call adelay
      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
      endif

      call i8copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine csend(mtype,buf,len,jnid,jpid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call mpi_send (buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine crecv(mtype,buf,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call adelay
      call mpi_recv (buf,len,mpi_byte
     $              ,jnid,mtype,nekcomm,status,ierr)
c
      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine crecv3(mtype,buf,len,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call adelay
      call mpi_recv (buf,len,mpi_byte
     $            ,jnid,mtype,nekcomm,status,ierr)
      call mpi_get_count (status,mpi_byte,len,ierr)
c
      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      integer function numnodes()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_comm_size (nekcomm, numnodes , ierr)

      return
      end
c-----------------------------------------------------------------------
      integer function mynode()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer myid

      call mpi_comm_rank (nekcomm, myid, ierr)
      mynode = myid

      return
      end
c-----------------------------------------------------------------------
      real*8 function dnekclock()
      include 'mpif.h'
c
      dnekclock=mpi_wtime()
c
      return
      end
c-----------------------------------------------------------------------
      real*8 function dnekclock_sync()
      include 'mpif.h'
c
      call nekgsync()
      dnekclock_sync=mpi_wtime()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lbcast(ifif)
C
C     Broadcast logical variable to all processors.
C
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      logical ifif

      if (np.eq.1) return

      item=0
      if (ifif) item=1
      call adelay
      call bcast(item,isize)
      ifif=.false.
      if (item.eq.1) ifif=.true.

      return
      end
c-----------------------------------------------------------------------
      subroutine bcast(buf,len)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call adelay
      call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_comm(intracomm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

c      call mpi_comm_group (mpi_comm_world,itmp,ierr)
c      call mpi_comm_create (mpi_comm_world,itmp,icomm,ierr)
c      call mpi_group_free (itmp,ierr)

      call mpi_comm_dup(intracomm,nekcomm,ierr)

c     write(6,*) 'nekcomm:',nekcomm

      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
      call adelay
      call mpi_isend (x,len,mpi_byte,jnid,msgtag
     $       ,nekcomm,imsg,ierr)
      isend = imsg
c     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)
c
      return
      end
c-----------------------------------------------------------------------
      function irecv(msgtag,x,len)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
      call adelay
      call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,nekcomm,imsg,ierr)
      irecv = imsg
c     write(6,*) nid,' irecv:',imsg,msgtag,len
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine msgwait(imsg)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
c
c     write(6,*) nid,' msgwait:',imsg
c
      call adelay
      call mpi_wait (imsg,status,ierr)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine nekgsync()

      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call adelay
      call mpi_barrier(nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine exittr(stringi,rdata,idata)
      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      include 'SIZE'
      include 'TOTAL'

      call blank(stringo,132)
      call chcopy(stringo,stringi,132)
      len = indx1(stringo,'$',1)
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call chcopy(stringo(len),s25,25)

      if (nid.eq.0) write(6,1) (stringo(k),k=1,len+24)
    1 format('EXIT: ',132a1)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine exitti(stringi,idata)
      character*1 stringi(132)
      character*1 stringo(132)
      character*11 s11
      include 'SIZE'
      include 'TOTAL'

      call blank(stringo,132)
      call chcopy(stringo,stringi,132)
      len = indx1(stringo,'$',1)
      write(s11,11) idata
   11 format(1x,i10)
      call chcopy(stringo(len),s11,11)

      if (nid.eq.0) write(6,1) (stringo(k),k=1,len+10)
    1 format('EXIT: ',132a1)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine err_chk(ierr,string)
      character*1 string(132)
      character*1 ostring(132)
      character*10 s10
      include 'SIZE'
      include 'TOTAL'

      ierr = iglsum(ierr,1)
      if(ierr.eq.0) return 

      len = indx1(string,'$',1)
      call blank(ostring,132)
      write(s10,11) ierr
   11 format(1x,' ierr=',i3)

      call chcopy(ostring,string,len-1)
      call chcopy(ostring(len),s10,10)

      if (nid.eq.0) write(6,1) (ostring(k),k=1,len+10)
    1 format('ERROR: ',132a1)

      call exitt

      return
      end
c
c-----------------------------------------------------------------------
      subroutine exitt0
      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'

      real*4 papi_mflops
      integer*8 papi_flops

      if(nid.eq.0) write(6,*) 'Exitting....'

c     call print_stack()

      call mpi_finalize (ierr)
      call exit(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt
      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'

      real*4 papi_mflops
      integer*8 papi_flops
c
      call nekgsync()

#ifdef PAPI
      call nek_flops(papi_flops,papi_mflops)
#endif

      tstop  = dnekclock()
      ttotal = tstop-etimes
      nxyz   = nx1*ny1*nz1

      if (nid.eq.0) then 
         dtmp1 = 0
         dtmp2 = 0
         dtmp3 = 0
         if(istep.gt.0) then
           dtmp3 = 1.*papi_flops/1e6
         endif 
         write(6,*) ' '
         write(6,'(A)') 'call exitt: dying ...'
         write(6,*) ' '
c        call print_stack()
         write(6,*) ' '
         write(6,'(1(A,1p1e13.5,A,/))') 
     &       'total elapsed time             : ',ttotal, ' sec'
#ifdef PAPI
         write(6,'(2(A,1g13.5,/))') 
     &       'Gflops                         : ',dtmp3/1000.
     &      ,'Gflops/s                       : ',papi_mflops/1000.
#endif
      endif 

      nz1 = 1/(nx1-ny1)

      call mpi_finalize (ierr)
      call exit(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine printHeader


      return
      end
c-----------------------------------------------------------------------
      function igl_running_sum(in)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      integer x,w,r

      x = in  ! running sum
      w = in  ! working buff
      r = 0   ! recv buff

      call mpi_scan(x,r,1,mpi_integer,mpi_sum,nekcomm,ierr)
      igl_running_sum = r

      return
      end
c-----------------------------------------------------------------------
      subroutine platform_timer(ivb) ! mxm, ping-pong, and all_reduce timer

      include 'SIZE'
      include 'TOTAL'


      call mxm_test_all(nid,ivb)  ! measure mxm times
c     call exitti('done mxm_test_all$',ivb)

      call comm_test(ivb)         ! measure message-passing and all-reduce times

      return
      end
c-----------------------------------------------------------------------
      subroutine comm_test(ivb) ! measure message-passing and all-reduce times
                                ! ivb = 0 --> minimal verbosity
                                ! ivb = 1 --> fully verbose
                                ! ivb = 2 --> smaller sample set(shorter)

      include 'SIZE'
      include 'PARALLEL'

      call gop_test(ivb)   ! added, Jan. 8, 2008

      log_np=log2(np)
      np2 = 2**log_np
      if (np2.eq.np) call gp2_test(ivb)   ! added, Jan. 8, 2008

      io = 6
      n512 = min(512,np-1)

      do nodec=1,n512
         nodeb=nodec
         if (nodec.gt.256.and.np.gt.512)
     $      nodeb = 256+(np-256)*(nodec-256)/(n512-256) - 1
         nodeb = min(nodeb,np-1)

         call pingpong(alphas,betas,0,nodeb,.0005,io,ivb)
         if (nid.eq.0) write(6,1) nodeb,np,alphas,betas
    1    format(2i10,1p2e15.7,' alpha beta')
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pingpong(alphas,betas,nodea,nodeb,dt,io,ivb)

      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd1 = 3*lt,mwd2=100 000, mwd = max(mwd1,mwd2))
      common /scrns/ x(mwd),y(mwd)

      include 'mpif.h'
      integer status(mpi_status_size)

      character*10 fname

      if (nid.eq.nodea) then
         write(fname,3) np,nodeb
    3    format('t',i4.4,'.',i4.4)
         if (io.ne.6) open (unit=io,file=fname)
      endif

      call nekgsync
      call get_msg_vol(msg_vol,dt,nodea,nodeb) ! Est. msg vol for dt s

      nwds = 0
      if (nid.eq.nodea.and.ivb.gt.0) write(io,*)

      betas = 0  ! Reported inverse bandwidth
      count = 0

      do itest = 1,500

         nloop = msg_vol/(nwds+2)
         nloop = min(nloop,1000)
         nloop = max(nloop,1)

         len   = 8*nwds
     
         call ping_loop(t1,t0,len,nloop,nodea,nodeb,nid,x,y,x,y)

         if (nid.eq.nodea) then
            tmsg = (t1-t0)/(2*nloop)   ! 2*nloop--> Double Buffer
            tmsg = tmsg / 2.           ! one-way cost = 1/2 round-trip
            tpwd = tmsg                ! time-per-word
            if (nwds.gt.0) tpwd = tmsg/nwds
            if (ivb.gt.0) write(io,1) nodeb,np,nloop,nwds,tmsg,tpwd
    1       format(3i6,i12,1p2e16.8,' pg')

            if (nwds.eq.1) then
               alphas = tmsg
            elseif (nwds.gt.10000) then   ! "average" beta
               betas = (betas*count + tpwd)/(count+1)
               count = count + 1
            endif
         endif

         if (ivb.eq.2) then
            nwds = (nwds+1)*1.25
         else
            nwds = (nwds+1)*1.016
         endif
         if (nwds.gt.mwd) then
c        if (nwds.gt.1024) then
            if (nid.eq.nodea.and.io.ne.6) close(unit=io)
            call nekgsync
            return
         endif

      enddo

      if (nid.eq.nodea.and.io.ne.6) close(unit=io)
      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine pingpongo(alphas,betas,nodea,nodeb,dt,io,ivb)

      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd1 = 3*lt,mwd2=100 000, mwd = max(mwd1,mwd2))
      common /scrns/ x(mwd),y(mwd)

      include 'mpif.h'
      integer status(mpi_status_size)

      character*10 fname

      if (nid.eq.nodea) then
         write(fname,3) np,nodeb
    3    format('t',i4.4,'.',i4.4)
         if (io.ne.6) open (unit=io,file=fname)
      endif

      call nekgsync
      call get_msg_vol(msg_vol,dt,nodea,nodeb) ! Est. msg vol for dt s

      nwds = 0
      if (nid.eq.nodea.and.ivb.gt.0) write(io,*)

      betas = 0  ! Reported inverse bandwidth
      count = 0

      do itest = 1,500
         call nekgsync
         nloop = msg_vol/(nwds+2)
         nloop = min(nloop,1000)
         nloop = max(nloop,1)

         len   = 8*nwds
         jnid = mpi_any_source

         if (nid.eq.nodea) then

            msg  = irecv(itest,y,1)
            call csend(itest,x,1,nodeb,0)   ! Initiate send, to synch.
            call msgwait(msg)

            t0 = mpi_wtime ()
            do i=1,nloop
               call mpi_irecv(y,len,mpi_byte,mpi_any_source,i
     $                        ,nekcomm,msg,ierr)
               call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)
               call mpi_wait (msg,status,ierr)
            enddo
            t1 = mpi_wtime ()
            tmsg = (t1-t0)/nloop
            tmsg = tmsg / 2.       ! Round-trip message time = twice one-way
            tpwd = tmsg
            if (nwds.gt.0) tpwd = tmsg/nwds
            if (ivb.gt.0) write(io,1) nodeb,np,nloop,nwds,tmsg,tpwd
    1       format(3i6,i12,1p2e16.8,' pg')

            if (nwds.eq.1) then
               alphas = tmsg
            elseif (nwds.gt.10000) then
               betas = (betas*count + tpwd)/(count+1)
               count = count + 1
            endif

         elseif (nid.eq.nodeb) then

            call crecv(itest,y,1)           ! Initiate send, to synch.
            call csend(itest,x,1,nodea,0)

            t0 = dnekclock()
            do i=1,nloop
               call mpi_recv (y,len,mpi_byte
     $               ,jnid,i,nekcomm,status,ierr)
               call mpi_send (x,len,mpi_byte,nodea,i,nekcomm,ierr)
            enddo
            t1 = dnekclock()
            tmsg = (t1-t0)/nloop

         endif

         nwds = (nwds+1)*1.016
         if (nwds.gt.mwd) then
            if (nid.eq.nodea.and.io.ne.6) close(unit=io)
            call nekgsync
            return
         endif

      enddo

      if (nid.eq.nodea.and.io.ne.6) close(unit=io)
      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine get_msg_vol(msg_vol,dt,nodea,nodeb)
      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ x(3*lt),y(3*lt)
!
!     Est. msg vol for dt s
!
      msg_vol = 1000

      nwds  = min(1000,lt)
      nloop = 50
 
      tmsg = 0.
      call gop(tmsg,t1,'+  ',1)

      len = 8*nwds
      if (nid.eq.nodea) then

         msg  = irecv(1,y,1)
         call csend(1,x,1,nodeb,0)   ! Initiate send, to synch.
         call msgwait(msg)

         t0 = dnekclock()
         do i=1,nloop
            msg  = irecv(i,y,len)
            call csend(i,x,len,nodeb,0)
            call msgwait(msg)
         enddo
         t1   = dnekclock()
         tmsg = (t1-t0)/nloop
         tpwd = tmsg/nwds

      elseif (nid.eq.nodeb) then

         call crecv(1,y,1)           ! Initiate send, to synch.
         call csend(1,x,1,nodea,0)

         t0 = dnekclock()
         do i=1,nloop
            call crecv(i,y,len)
            call csend(i,x,len,nodea,0)
         enddo
         t1   = dnekclock()
         tmsg = (t1-t0)/nloop
         tmsg = 0.

      endif

      call gop(tmsg,t1,'+  ',1)
      msg_vol = nwds*(dt/tmsg)
c     if (nid.eq.nodea) write(6,*) nid,msg_vol,nwds,dt,tmsg,' msgvol'

      return
      end
c-----------------------------------------------------------------------
      subroutine gop_test(ivb)
      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      include 'mpif.h'
      integer status(mpi_status_size)

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd = 3*lt)
      common /scrns/ x(mwd),y(mwd)
      common /scruz/ times(2,500)
      common/scrcg/nwd(500)


      nwds = 1
      mtest = 0
      do itest = 1,500
         nwds = (nwds+1)*1.016
         if (nwds.gt.mwd) goto 100
         mtest = mtest+1
         nwd(mtest) = nwds
      enddo
  100 continue

      nwds = 1
      do itest = mtest,1,-1

         tiny = 1.e-27
         call cfill(x,tiny,mwd)
         nwds = nwd(itest)
         call nekgsync


         t0 = mpi_wtime ()
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         t1 = mpi_wtime ()

         tmsg = (t1-t0)/6 ! six calls
         tpwd = tmsg
         if (nwds.gt.0) tpwd = tmsg/nwds
         times(1,itest) = tmsg
         times(2,itest) = tpwd

      enddo
  101 continue


      if (nid.eq.0) then
         nwds = 1
         do itest=1,500
            if (ivb.gt.0.or.itest.eq.1) 
     $         write(6,1) np,nwds,(times(k,itest),k=1,2)
    1       format(i9,i12,1p2e16.8,' gop')
            nwds = (nwds+1)*1.016
            if (nwds.gt.mwd) goto 102
         enddo
  102    continue
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gp2_test(ivb)

      include 'SIZE'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd = 3*lt)
      common /scrns/ x(mwd),y(mwd)
      common /scruz/ times(2,500)

      call rzero(x,mwd)

      nwds = 1
      do itest = 1,500
         call gp2(x,y,'+  ',1,nid,np)

         t0 = mpi_wtime ()
         call gp2(x,y,'+  ',nwds,nid,np)
         call gp2(x,y,'+  ',nwds,nid,np)
         call gp2(x,y,'+  ',nwds,nid,np)
         call gp2(x,y,'+  ',nwds,nid,np)
         t1 = mpi_wtime ()

         tmsg = (t1-t0)/4 ! four calls
         tpwd = tmsg
         if (nwds.gt.0) tpwd = tmsg/nwds
         times(1,itest) = tmsg
         times(2,itest) = tpwd

         nwds = (nwds+1)*1.016
         if (nwds.gt.mwd) goto 101
      enddo
  101 continue


      if (nid.eq.0) then
         nwds = 1
         do itest=1,500
            if (ivb.gt.0.or.itest.eq.1) 
     $         write(6,1) np,nwds,(times(k,itest),k=1,2)
    1       format(i9,i12,1p2e16.8,' gp2')
            nwds = (nwds+1)*1.016
            if (nwds.gt.mwd) goto 102
         enddo
  102    continue
      endif

      return
      end
c-----------------------------------------------------------------------
      integer function xor(m,n)
c
c  If NOT running on a parallel processor, it is sufficient to
c  have this routine return a value of XOR=1.
c
c  Pick one of the following:
c
c  UNIX 4.2, f77:
       XOR = OR(M,N)-AND(M,N)
c
c  Intel FTN286:
c     XOR = M.NEQV.N
c
c  Ryan-McFarland Fortran
C      XOR = IEOR(M,N)
c
c     XOR = 0
c     IF(M.EQ.1 .OR.  N.EQ.1) XOR=1
c     IF(M.EQ.0 .AND. N.EQ.0) XOR=0
c     IF(M.EQ.1 .AND. N.EQ.1) XOR=0
c     IF(M.GT.1 .OR.N.GT.1 .OR.M.LT.0.OR.N.LT.0) THEN
c        PRINT*,'ERROR IN XOR'
c        STOP
c     ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine gp2( x, w, op, n, nid, np)
c
c     Global vector commutative operation using spanning tree.
c
c     Std. fan-in/fan-out

      real x(n), w(n)
      character*3 op

      integer bit, bytes, cnt, diff, spsize, i, 
     *   parent, troot, xor, root, lnp, log2
      logical ifgot

      integer type
      save    type
      data    type  /998/

      type  = type+100
      if (type.gt.9992) type=type-998
      typer = type-1
      bytes = 8*n

      root    = 0
      troot   = max0((nid/np)*np, root)
      diff    = xor(nid,troot)
      nullpid = 0

c     Accumulate contributions from children, if any
      level2=1
    5 continue
         level=level2
         level2=level+level
         if (mod(nid,level2).ne.0) goto 20
            call crecv(type,w,bytes)
            if (op.eq.'+  ') then
               do i=1,n
                  x(i) = x(i) + w(i)
               enddo
            elseif (op.eq.'*  ') then
               do i=1,n
                  x(i) = x(i) * w(i)
               enddo
            elseif (op.eq.'M  ') then
               do i=1,n
                  x(i) = max(x(i),w(i))
               enddo
            elseif (op.eq.'m  ') then
               do i=1,n
                  x(i) = min(x(i),w(i))
               enddo
            endif
         if (level2.lt.np) goto 5

c     Pass result back to parent
   20 parent = nid-level
      if (nid .ne. 0) call csend(type,x,bytes,parent,nullpid)

c     Await final answer from node 0 via log_2 fan out
      level=np/2
      ifgot=.false.
      if (nid.eq.root) ifgot=.true.

      lnp = log2(np)
      do i=1,lnp
        if (ifgot) then
           jnid=nid+level
           call csend(typer,x,bytes,jnid,nullpid)
        elseif (mod(nid,level).eq.0) then
           call crecv(typer,x,bytes)
           ifgot=.true.
        endif
        level=level/2
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ping_loop1(t1,t0,len,nloop,nodea,nodeb,nid,x,y)

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real x(1),y(1)

      include 'mpif.h'
      integer status(mpi_status_size)

      i=0
      if (nid.eq.nodea) then
         call nekgsync
         call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr)    ! 1b
         call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
c        call mpi_rsend(x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
         call msgwait(msg)                                          ! 1b

         t0 = mpi_wtime ()
         do i=1,nloop
            call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr) ! 2b
            call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
c           call mpi_rsend(x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
            call mpi_wait (msg,status,ierr)                         ! 2b
         enddo
         t1 = mpi_wtime ()

      elseif (nid.eq.nodeb) then

         call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr)    ! 1a
         call nekgsync
         call mpi_wait (msg,status,ierr)                            ! 1a

         j=i
         do i=1,nloop
            call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr) ! 2a
c           call mpi_rsend(x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_wait (msg,status,ierr)                         ! 2a
            j=i
         enddo
c        call mpi_rsend(x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb
         call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb

      else
         call nekgsync
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ping_loop2(t1,t0,len,nloop,nodea,nodeb,nid,x,y)

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real x(1),y(1)

      include 'mpif.h'
      integer status(mpi_status_size)

      i=0
      if (nid.eq.nodea) then
         call nekgsync
         call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr)    ! 1b
         call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
         call msgwait(msg)                                          ! 1b

         t0 = mpi_wtime ()
         do i=1,nloop
            call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
            call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr) ! 2b
            call mpi_wait (msg,status,ierr)                         ! 2b
         enddo
         t1 = mpi_wtime ()

      elseif (nid.eq.nodeb) then

         call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr)    ! 1a
         call nekgsync
         call mpi_wait (msg,status,ierr)                            ! 1a

         j=i
         do i=1,nloop
            call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr) ! 2a
            call mpi_wait (msg,status,ierr)                         ! 2a
            j=i
         enddo
         call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb

      else
         call nekgsync
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ping_loop(t1,t0,len,nloop,nodea,nodeb,nid,x1,y1,x2,y2)
c     Double Buffer : does 2*nloop timings

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real x1(1),y1(1),x2(1),y2(1)

      include 'mpif.h'
      integer status(mpi_status_size)

      itag=1
      if (nid.eq.nodea) then
         call mpi_irecv(y1,len,mpi_byte,nodeb,itag,nekcomm,msg1,ierr)   ! 1b 
         call nekgsync


         t0 = mpi_wtime ()
         do i=1,nloop
            call mpi_send (x1,len,mpi_byte,nodeb,itag,nekcomm,ierr)     ! 1a 
            call mpi_irecv(y2,len,mpi_byte,nodeb,itag,nekcomm,msg2,ierr)! 2b 
            call mpi_wait (msg1,status,ierr)                            ! 1b
            call mpi_send (x2,len,mpi_byte,nodeb,itag,nekcomm,ierr)     ! 2a 
            call mpi_irecv(y1,len,mpi_byte,nodeb,itag,nekcomm,msg1,ierr)! 3b 
            call mpi_wait (msg2,status,ierr)                            ! 2b
         enddo
         t1 = mpi_wtime ()
         call mpi_send (x1,len,mpi_byte,nodeb,itag,nekcomm,ierr)        ! nb
         call mpi_wait (msg1,status,ierr)                              ! nb

      elseif (nid.eq.nodeb) then

         call mpi_irecv(y1,len,mpi_byte,nodea,itag,nekcomm,msg1,ierr)   ! nb 
         call nekgsync


         do i=1,nloop
            call mpi_wait (msg1,status,ierr)                            ! 1a
            call mpi_send (x1,len,mpi_byte,nodea,itag,nekcomm,ierr)     ! 1b
            call mpi_irecv(y2,len,mpi_byte,nodea,itag,nekcomm,msg2,ierr)! 2a
            call mpi_wait (msg2,status,ierr)                            ! 2a 
            call mpi_send (x2,len,mpi_byte,nodea,itag,nekcomm,ierr)     ! 2b
            call mpi_irecv(y1,len,mpi_byte,nodea,itag,nekcomm,msg1,ierr)! 3a
         enddo
         call mpi_wait (msg1,status,ierr)                            ! 2a 
         call mpi_send (x1,len,mpi_byte,nodea,itag,nekcomm,ierr)        ! nb

      else
         call nekgsync
      endif

      return
      end


