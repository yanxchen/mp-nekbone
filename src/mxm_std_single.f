      subroutine mxmf2_f(A,N1,B,N2,C,N3)
c
c     unrolled loop version 
c
      real*4 a(n1,n2),b(n2,n3),c(n1,n3)

      if (n2.le.8) then
         if (n2.eq.1) then
            call mxf1_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxf2_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxf3_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxf4_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxf5_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxf6_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxf7_f(a,n1,b,n2,c,n3)
         else
            call mxf8_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxf9_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxf10_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxf11_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxf12_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxf13_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxf14_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxf15_f(a,n1,b,n2,c,n3)
         else
            call mxf16_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call mxf17_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call mxf18_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call mxf19_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call mxf20_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call mxf21_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call mxf22_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call mxf23_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call mxf24_f(a,n1,b,n2,c,n3)
         endif
      else
         call mxm44_0_f(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf1_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,1),b(1,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf2_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,2),b(2,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf3_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,3),b(3,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf4_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,4),b(4,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf5_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,5),b(5,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf6_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,6),b(6,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf7_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,7),b(7,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf8_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,8),b(8,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf9_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,9),b(9,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf10_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,10),b(10,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf11_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,11),b(11,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf12_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,12),b(12,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf13_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,13),b(13,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf14_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,14),b(14,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf15_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,15),b(15,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf16_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,16),b(16,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf17_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,17),b(17,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf18_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,18),b(18,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf19_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,19),b(19,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf20_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,20),b(20,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf21_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,21),b(21,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf22_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,22),b(22,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf23_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,23),b(23,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf24_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,24),b(24,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
     $             + a(i,24)*b(24,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxm44_0_f(a, m, b, k, c, n)
c
c matrix multiply with a 4x4 pencil 
c
      real*4 a(m,k), b(k,n), c(m,n)
      real*4 s11, s12, s13, s14, s21, s22, s23, s24
      real*4 s31, s32, s33, s34, s41, s42, s43, s44

      mresid = iand(m,3) 
      nresid = iand(n,3) 
      m1 = m - mresid + 1
      n1 = n - nresid + 1

      do i=1,m-mresid,4
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s41 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s42 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          s43 = 0.0d0
          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0
          s44 = 0.0d0
          do l=1,k
            s11 = s11 + a(i,l)*b(l,j)
            s12 = s12 + a(i,l)*b(l,j+1)
            s13 = s13 + a(i,l)*b(l,j+2)
            s14 = s14 + a(i,l)*b(l,j+3)

            s21 = s21 + a(i+1,l)*b(l,j)
            s22 = s22 + a(i+1,l)*b(l,j+1)
            s23 = s23 + a(i+1,l)*b(l,j+2)
            s24 = s24 + a(i+1,l)*b(l,j+3)

            s31 = s31 + a(i+2,l)*b(l,j)
            s32 = s32 + a(i+2,l)*b(l,j+1)
            s33 = s33 + a(i+2,l)*b(l,j+2)
            s34 = s34 + a(i+2,l)*b(l,j+3)

            s41 = s41 + a(i+3,l)*b(l,j)
            s42 = s42 + a(i+3,l)*b(l,j+1)
            s43 = s43 + a(i+3,l)*b(l,j+2)
            s44 = s44 + a(i+3,l)*b(l,j+3)
          enddo
          c(i,j)     = s11 
          c(i,j+1)   = s12 
          c(i,j+2)   = s13
          c(i,j+3)   = s14

          c(i+1,j)   = s21 
          c(i+2,j)   = s31 
          c(i+3,j)   = s41 

          c(i+1,j+1) = s22
          c(i+2,j+1) = s32
          c(i+3,j+1) = s42

          c(i+1,j+2) = s23
          c(i+2,j+2) = s33
          c(i+3,j+2) = s43

          c(i+1,j+3) = s24
          c(i+2,j+3) = s34
          c(i+3,j+3) = s44
        enddo
* Residual when n is not multiple of 4
        if (nresid .ne. 0) then
          if (nresid .eq. 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,n)
              s21 = s21 + a(i+1,l)*b(l,n)
              s31 = s31 + a(i+2,l)*b(l,n)
              s41 = s41 + a(i+3,l)*b(l,n)
            enddo
            c(i,n)     = s11 
            c(i+1,n)   = s21 
            c(i+2,n)   = s31 
            c(i+3,n)   = s41 
          elseif (nresid .eq. 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
          else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)
              s13 = s13 + a(i,l)*b(l,j+2)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)
              s23 = s23 + a(i+1,l)*b(l,j+2)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)
              s33 = s33 + a(i+2,l)*b(l,j+2)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
              s43 = s43 + a(i+3,l)*b(l,j+2)
            enddo
            c(i,j)     = s11 
            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 
            c(i,j+1)   = s12 
            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
            c(i,j+2)   = s13
            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43
          endif
        endif
      enddo

* Residual when m is not multiple of 4
      if (mresid .eq. 0) then
        return
      elseif (mresid .eq. 1) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,j)
            s12 = s12 + a(m,l)*b(l,j+1)
            s13 = s13 + a(m,l)*b(l,j+2)
            s14 = s14 + a(m,l)*b(l,j+3)
          enddo
          c(m,j)     = s11 
          c(m,j+1)   = s12 
          c(m,j+2)   = s13
          c(m,j+3)   = s14
        enddo
* mresid is 1, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n)
          enddo
          c(m,n) = s11
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s12 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-1)
            s12 = s12 + a(m,l)*b(l,n)
          enddo
          c(m,n-1) = s11
          c(m,n) = s12
          return
        else
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-2)
            s12 = s12 + a(m,l)*b(l,n-1)
            s13 = s13 + a(m,l)*b(l,n)
          enddo
          c(m,n-2) = s11
          c(m,n-1) = s12
          c(m,n) = s13
          return
        endif          
      elseif (mresid .eq. 2) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          s21 = 0.0d0
          s22 = 0.0d0
          s23 = 0.0d0
          s24 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,j)
            s12 = s12 + a(m-1,l)*b(l,j+1)
            s13 = s13 + a(m-1,l)*b(l,j+2)
            s14 = s14 + a(m-1,l)*b(l,j+3)

            s21 = s21 + a(m,l)*b(l,j)
            s22 = s22 + a(m,l)*b(l,j+1)
            s23 = s23 + a(m,l)*b(l,j+2)
            s24 = s24 + a(m,l)*b(l,j+3)
          enddo
          c(m-1,j)   = s11 
          c(m-1,j+1) = s12 
          c(m-1,j+2) = s13
          c(m-1,j+3) = s14
          c(m,j)     = s21
          c(m,j+1)   = s22 
          c(m,j+2)   = s23
          c(m,j+3)   = s24
        enddo
* mresid is 2, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n)
          enddo
          c(m-1,n) = s11
          c(m,n) = s21
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-1)
            s12 = s12 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-1)
            s22 = s22 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-1) = s11
          c(m-1,n)   = s12
          c(m,n-1)   = s21
          c(m,n)     = s22
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-2)
            s12 = s12 + a(m-1,l)*b(l,n-1)
            s13 = s13 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-2)
            s22 = s22 + a(m,l)*b(l,n-1)
            s23 = s23 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-2) = s11
          c(m-1,n-1) = s12
          c(m-1,n)   = s13
          c(m,n-2)   = s21
          c(m,n-1)   = s22
          c(m,n)     = s23
          return
        endif
      else
* mresid is 3
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0

          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0

          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0

          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0

          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,j)
            s12 = s12 + a(m-2,l)*b(l,j+1)
            s13 = s13 + a(m-2,l)*b(l,j+2)
            s14 = s14 + a(m-2,l)*b(l,j+3)

            s21 = s21 + a(m-1,l)*b(l,j)
            s22 = s22 + a(m-1,l)*b(l,j+1)
            s23 = s23 + a(m-1,l)*b(l,j+2)
            s24 = s24 + a(m-1,l)*b(l,j+3)

            s31 = s31 + a(m,l)*b(l,j)
            s32 = s32 + a(m,l)*b(l,j+1)
            s33 = s33 + a(m,l)*b(l,j+2)
            s34 = s34 + a(m,l)*b(l,j+3)
          enddo
          c(m-2,j)   = s11 
          c(m-2,j+1) = s12 
          c(m-2,j+2) = s13
          c(m-2,j+3) = s14

          c(m-1,j)   = s21 
          c(m-1,j+1) = s22
          c(m-1,j+2) = s23
          c(m-1,j+3) = s24

          c(m,j)     = s31 
          c(m,j+1)   = s32
          c(m,j+2)   = s33
          c(m,j+3)   = s34
        enddo
* mresid is 3, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n)
          enddo
          c(m-2,n) = s11
          c(m-1,n) = s21
          c(m,n) = s31
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-1)
            s12 = s12 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-1)
            s22 = s22 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-1)
            s32 = s32 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-1) = s11
          c(m-2,n)   = s12
          c(m-1,n-1) = s21
          c(m-1,n)   = s22
          c(m,n-1)   = s31
          c(m,n)     = s32
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-2)
            s12 = s12 + a(m-2,l)*b(l,n-1)
            s13 = s13 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-2)
            s22 = s22 + a(m-1,l)*b(l,n-1)
            s23 = s23 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-2)
            s32 = s32 + a(m,l)*b(l,n-1)
            s33 = s33 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-2) = s11
          c(m-2,n-1) = s12
          c(m-2,n)   = s13
          c(m-1,n-2) = s21
          c(m-1,n-1) = s22
          c(m-1,n)   = s23
          c(m,n-2)   = s31
          c(m,n-1)   = s32
          c(m,n)     = s33
          return
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mxm44_2_f(a, m, b, k, c, n)
      real*4 a(m,2), b(2,n), c(m,n)

      nresid = iand(n,3) 
      n1 = n - nresid + 1

      do j=1,n-nresid,4
         do i=1,m
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
            c(i,j+1) = a(i,1)*b(1,j+1)
     $             + a(i,2)*b(2,j+1)
            c(i,j+2) = a(i,1)*b(1,j+2)
     $             + a(i,2)*b(2,j+2)
            c(i,j+3) = a(i,1)*b(1,j+3)
     $             + a(i,2)*b(2,j+3)
         enddo
      enddo
      if (nresid .eq. 0) then
        return
      elseif (nresid .eq. 1) then
         do i=1,m
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      elseif (nresid .eq. 2) then
         do i=1,m
            c(i,n-1) = a(i,1)*b(1,n-1)
     $             + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      else
         do i=1,m
            c(i,n-2) = a(i,1)*b(1,n-2)
     $             + a(i,2)*b(2,n-2)
            c(i,n-1) = a(i,1)*b(1,n-1)
     $             + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mxm_test_all_f(nid,ivb)
c
c     Collect matrix-matrix product statistics
c
      external mxms,mxmur2,mxmur3,mxmd,mxmfb,mxmf3,mxmu4,mxmn2
      external mxmk2,mxmtr,mxmrg,madd,mxm,mxm44
c
      parameter (nn=24)
      parameter (nt=10)
      character*5 c(3,nt)
      real*4        s(nn,2,nt,3)
      real*4        a(nn,2,nt,3)

      call nekgsync

      do k=1,3   ! 3 tests:  N^2 x N, NxN, NxN^2
         call mxmtest_f(s(1,1, 1,k),nn,c(k, 1),mxm44 ,'mxm44',k,ivb)
         call mxmtest_f(s(1,1, 2,k),nn,c(k, 2),mxms  ,' std ',k,ivb)
         call mxmtest_f(s(1,1, 3,k),nn,c(k, 3),mxmur2,'mxmu2',k,ivb)
         call mxmtest_f(s(1,1, 4,k),nn,c(k, 4),mxmur3,'mxmu3',k,ivb)
         call mxmtest_f(s(1,1, 5,k),nn,c(k, 5),mxmd  ,'mxmd ',k,ivb)
         call mxmtest_f(s(1,1, 6,k),nn,c(k, 6),mxmfb ,'mxmfb',k,ivb)
         call mxmtest_f(s(1,1, 7,k),nn,c(k, 7),mxmu4 ,'mxmu4',k,ivb)
         call mxmtest_f(s(1,1, 8,k),nn,c(k, 8),mxmf3 ,'mxmf3',k,ivb)
         if (k.eq.2) ! Add works only for NxN case
     $   call mxmtest_f(s(1,1, 9,k),nn,c(k, 9),madd  ,'madd ',k,ivb)
         call mxmtest_f(s(1,1,10,k),nn,c(k,10),mxm   ,'mxm  ',k,ivb)
      enddo

      call nekgsync
      if (nid.eq.0) call mxm_analyze_f(s,a,nn,c,nt,ivb)
      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine initab_f(a,b,n)
      real*4 a(1),b(1)
      do i=1,n-1
         x  = i
         k = mod(i,19) + 2
         l = mod(i,17) + 5
         m = mod(i,31) + 3
         a(i) = -.25*(a(i)+a(i+1)) + (x*x + k + l)/(x*x+m)
         b(i) = -.25*(b(i)+b(i+1)) + (x*x + k + m)/(x*x+l)
      enddo
      a(n) = -.25*(a(n)+a(n)) + (x*x + k + l)/(x*x+m)
      b(n) = -.25*(b(n)+b(n)) + (x*x + k + m)/(x*x+l)
      return
      end
c-----------------------------------------------------------------------
      subroutine mxms_f(a,n1,b,n2,c,n3)
C----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C---------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
         N0=N1*N3
         DO 10 I=1,N0
            C(I,1)=0.
 10      CONTINUE
         DO 100 J=1,N3
         DO 100 K=1,N2
         BB=B(K,J)
         DO 100 I=1,N1
            C(I,J)=C(I,J)+A(I,K)*BB
 100     CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmu4_f(a,n1,b,n2,c,n3)
C----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C---------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
         N0=N1*N3
         DO 10 I=1,N0
            C(I,1)=0.
 10      CONTINUE
         i1 = n1 - mod(n1,4) + 1
            DO 100 J=1,N3
            DO 100 K=1,N2
            BB=B(K,J)
               DO I=1,N1-3,4
                  C(I  ,J)=C(I  ,J)+A(I  ,K)*BB
                  C(I+1,J)=C(I+1,J)+A(I+1,K)*BB
                  C(I+2,J)=C(I+2,J)+A(I+2,K)*BB
                  C(I+3,J)=C(I+3,J)+A(I+3,K)*BB
               enddo
               DO i=i1,N1
                  C(I  ,J)=C(I  ,J)+A(I  ,K)*BB
               enddo
 100        CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine madd_f (a,n1,b,n2,c,n3)
c
      real*4 a(n1,n2),b(n2,n3),c(n1,n3)
c
      do j=1,n3
      do i=1,n1
         c(i,j) = a(i,j)+b(i,j)
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmUR2_f(a,n1,b,n2,c,n3)
C----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C---------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
      if (n2.le.8) then
         if (n2.eq.1) then
            call mxmur2_1_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxmur2_2_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxmur2_3_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxmur2_4_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxmur2_5_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxmur2_6_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxmur2_7_f(a,n1,b,n2,c,n3)
         else
            call mxmur2_8_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxmur2_9_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxmur2_10_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxmur2_11_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxmur2_12_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxmur2_13_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxmur2_14_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxmur2_15_f(a,n1,b,n2,c,n3)
         else
            call mxmur2_16_f(a,n1,b,n2,c,n3)
         endif
      else
         N0=N1*N3
         DO 10 I=1,N0
            C(I,1)=0.
 10      CONTINUE
         DO 100 J=1,N3
         DO 100 K=1,N2
         BB=B(K,J)
         DO 100 I=1,N1
            C(I,J)=C(I,J)+A(I,K)*BB
 100     CONTINUE
      endif
      return
      end
c
      subroutine mxmur2_1_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,1),b(1,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_2_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,2),b(2,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_3_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,3),b(3,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_4_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,4),b(4,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_5_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,5),b(5,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_6_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,6),b(6,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_7_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,7),b(7,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_8_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,8),b(8,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_9_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,9),b(9,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_10_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,10),b(10,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_11_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,11),b(11,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_12_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,12),b(12,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_13_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,13),b(13,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_14_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,14),b(14,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_15_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,15),b(15,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
      subroutine mxmur2_16_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,16),b(16,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmUR3_f(a,n1,b,n2,c,n3)
C----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C---------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
      N0=N1*N3
      DO 10 I=1,N0
         C(I,1)=0.
 10   CONTINUE
      if (n3.le.8) then
         if (n3.eq.1) then
            call mxmur3_1_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.2) then
            call mxmur3_2_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.3) then
            call mxmur3_3_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.4) then
            call mxmur3_4_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.5) then
            call mxmur3_5_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.6) then
            call mxmur3_6_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.7) then
            call mxmur3_7_f(a,n1,b,n2,c,n3)
         else
            call mxmur3_8_f(a,n1,b,n2,c,n3)
         endif
      elseif (n3.le.16) then
         if (n3.eq.9) then
            call mxmur3_9_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.10) then
            call mxmur3_10_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.11) then
            call mxmur3_11_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.12) then
            call mxmur3_12_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.13) then
            call mxmur3_13_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.14) then
            call mxmur3_14_f(a,n1,b,n2,c,n3)
         elseif (n3.eq.15) then
            call mxmur3_15_f(a,n1,b,n2,c,n3)
         else
            call mxmur3_16_f(a,n1,b,n2,c,n3)
         endif
      else
         DO 100 J=1,N3
         DO 100 K=1,N2
         BB=B(K,J)
         DO 100 I=1,N1
            C(I,J)=C(I,J)+A(I,K)*BB
 100     CONTINUE
      endif
      return
      end
c
      subroutine mxmur3_16_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,16),c(n1,16)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         tmp11 =  b(k,11)
         tmp12 =  b(k,12)
         tmp13 =  b(k,13)
         tmp14 =  b(k,14)
         tmp15 =  b(k,15)
         tmp16 =  b(k,16)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
            c(i,11)  =  c(i,11) + a(i,k) * tmp11
            c(i,12)  =  c(i,12) + a(i,k) * tmp12
            c(i,13)  =  c(i,13) + a(i,k) * tmp13
            c(i,14)  =  c(i,14) + a(i,k) * tmp14
            c(i,15)  =  c(i,15) + a(i,k) * tmp15
            c(i,16)  =  c(i,16) + a(i,k) * tmp16
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_15_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,15),c(n1,15)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         tmp11 =  b(k,11)
         tmp12 =  b(k,12)
         tmp13 =  b(k,13)
         tmp14 =  b(k,14)
         tmp15 =  b(k,15)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
            c(i,11)  =  c(i,11) + a(i,k) * tmp11
            c(i,12)  =  c(i,12) + a(i,k) * tmp12
            c(i,13)  =  c(i,13) + a(i,k) * tmp13
            c(i,14)  =  c(i,14) + a(i,k) * tmp14
            c(i,15)  =  c(i,15) + a(i,k) * tmp15
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_14_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,14),c(n1,14)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         tmp11 =  b(k,11)
         tmp12 =  b(k,12)
         tmp13 =  b(k,13)
         tmp14 =  b(k,14)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
            c(i,11)  =  c(i,11) + a(i,k) * tmp11
            c(i,12)  =  c(i,12) + a(i,k) * tmp12
            c(i,13)  =  c(i,13) + a(i,k) * tmp13
            c(i,14)  =  c(i,14) + a(i,k) * tmp14
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_13_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,13),c(n1,13)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         tmp11 =  b(k,11)
         tmp12 =  b(k,12)
         tmp13 =  b(k,13)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
            c(i,11)  =  c(i,11) + a(i,k) * tmp11
            c(i,12)  =  c(i,12) + a(i,k) * tmp12
            c(i,13)  =  c(i,13) + a(i,k) * tmp13
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_12_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,12),c(n1,12)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         tmp11 =  b(k,11)
         tmp12 =  b(k,12)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
            c(i,11)  =  c(i,11) + a(i,k) * tmp11
            c(i,12)  =  c(i,12) + a(i,k) * tmp12
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_11_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,11),c(n1,11)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         tmp11 =  b(k,11)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
            c(i,11)  =  c(i,11) + a(i,k) * tmp11
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_10_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,10),c(n1,10)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         tmp10 =  b(k,10)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
            c(i,10)  =  c(i,10) + a(i,k) * tmp10
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_9_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,9),c(n1,9)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         tmp9  =  b(k, 9)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
            c(i, 9)  =  c(i, 9) + a(i,k) * tmp9
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_8_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,8),c(n1,8)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         tmp8  =  b(k, 8)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
            c(i, 8)  =  c(i, 8) + a(i,k) * tmp8
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_7_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,7),c(n1,7)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         tmp7  =  b(k, 7)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
            c(i, 7)  =  c(i, 7) + a(i,k) * tmp7
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_6_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,6),c(n1,6)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         tmp6  =  b(k, 6)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
            c(i, 6)  =  c(i, 6) + a(i,k) * tmp6
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_5_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,5),c(n1,5)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         tmp5  =  b(k, 5)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
            c(i, 5)  =  c(i, 5) + a(i,k) * tmp5
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_4_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,4),c(n1,4)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         tmp4  =  b(k, 4)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
            c(i, 4)  =  c(i, 4) + a(i,k) * tmp4
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_3_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,3),c(n1,3)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         tmp3  =  b(k, 3)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
            c(i, 3)  =  c(i, 3) + a(i,k) * tmp3
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_2_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,2),c(n1,2)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         tmp2  =  b(k, 2)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
            c(i, 2)  =  c(i, 2) + a(i,k) * tmp2
         enddo
c
      enddo
c
      return
      end
      subroutine mxmur3_1_f(a,n1,b,n2,c,n3)
      real*4 a(n1,n2),b(n2,1),c(n1,1)
c
      do k=1,n2
         tmp1  =  b(k, 1)
         do i=1,n1
            c(i, 1)  =  c(i, 1) + a(i,k) * tmp1
         enddo
      enddo
c
      return
      end
C----------------------------------------------------------------------
      subroutine mxmd_f(a,n1,b,n2,c,n3)
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C---------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
      REAL ONE,ZERO,EPS
C
C
C
      one=1.0
      zero=0.0
      call dgemm( 'N','N',n1,n3,n2,ONE,A,N1,B,N2,ZERO,C,N1)
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_f(a,n1,b,n2,c,n3)
C-----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C----------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
      integer wdsize
      save    wdsize
      data    wdsize/0/
c
c     First call: determine word size for dgemm/sgemm discrimination, below.
c
      if (wdsize.eq.0) then
         one = 1.0
         eps = 1.e-12
         wdsize = 8
         if (one+eps.eq.1.0) wdsize = 4
      endif
c
      if (n2.le.8) then
         if (n2.eq.1) then
            call mxmfb_1_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxmfb_2_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxmfb_3_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxmfb_4_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxmfb_5_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxmfb_6_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxmfb_7_f(a,n1,b,n2,c,n3)
         else
            call mxmfb_8_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxmfb_9_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxmfb_10_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxmfb_11_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxmfb_12_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxmfb_13_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxmfb_14_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxmfb_15_f(a,n1,b,n2,c,n3)
         else
            call mxmfb_16_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call mxmfb_17_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call mxmfb_18_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call mxmfb_19_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call mxmfb_20_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call mxmfb_21_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call mxmfb_22_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call mxmfb_23_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call mxmfb_24_f(a,n1,b,n2,c,n3)
         endif
      else
c
         one=1.0
         zero=0.0
         if (wdsize.eq.4) then
            call sgemm( 'N','N',n1,n3,n2,ONE,A,N1,B,N2,ZERO,C,N1)
         else
            call dgemm( 'N','N',n1,n3,n2,ONE,A,N1,B,N2,ZERO,C,N1)
         endif
 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_1_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,1),b(1,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_2_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,2),b(2,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_3_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,3),b(3,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_4_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,4),b(4,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_5_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,5),b(5,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_6_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,6),b(6,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_7_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,7),b(7,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_8_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,8),b(8,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_9_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,9),b(9,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_10_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,10),b(10,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_11_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,11),b(11,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_12_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,12),b(12,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_13_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,13),b(13,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_14_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,14),b(14,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_15_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,15),b(15,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_16_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,16),b(16,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_17_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,17),b(17,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_18_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,18),b(18,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_19_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,19),b(19,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_20_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,20),b(20,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_21_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,21),b(21,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_22_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,22),b(22,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_23_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,23),b(23,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmfb_24_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,24),b(24,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
     $             + a(i,24)*b(24,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_f(a,n1,b,n2,c,n3)
C-----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C----------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
      integer wdsize
      save    wdsize
      data    wdsize/0/
c
c     First call: determine word size for dgemm/sgemm discrimination, below.
c
      if (wdsize.eq.0) then
         one = 1.0
         eps = 1.e-12
         wdsize = 8
         if (one+eps.eq.1.0) wdsize = 4
      endif
c
      if (n2.le.8) then
         if (n2.eq.1) then
            call mxmf3_1_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxmf3_2_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxmf3_3_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxmf3_4_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxmf3_5_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxmf3_6_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxmf3_7_f(a,n1,b,n2,c,n3)
         else
            call mxmf3_8_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxmf3_9_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxmf3_10_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxmf3_11_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxmf3_12_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxmf3_13_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxmf3_14_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxmf3_15_f(a,n1,b,n2,c,n3)
         else
            call mxmf3_16_f(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call mxmf3_17_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call mxmf3_18_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call mxmf3_19_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call mxmf3_20_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call mxmf3_21_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call mxmf3_22_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call mxmf3_23_f(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call mxmf3_24_f(a,n1,b,n2,c,n3)
         endif
      else
c
         one=1.0
         zero=0.0
         if (wdsize.eq.4) then
            call sgemm( 'N','N',n1,n3,n2,ONE,A,N1,B,N2,ZERO,C,N1)
         else
            call dgemm( 'N','N',n1,n3,n2,ONE,A,N1,B,N2,ZERO,C,N1)
         endif
c
c        N0=N1*N3
c        DO 10 I=1,N0
c           C(I,1)=0.
c  10    CONTINUE
c        DO 100 J=1,N3
c        DO 100 K=1,N2
c        BB=B(K,J)
c        DO 100 I=1,N1
c           C(I,J)=C(I,J)+A(I,K)*BB
c 100    CONTINUE
 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_1_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,1),b(1,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_2_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,2),b(2,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_3_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,3),b(3,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_4_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,4),b(4,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_5_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,5),b(5,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_6_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,6),b(6,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_7_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,7),b(7,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_8_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,8),b(8,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_9_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,9),b(9,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_10_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,10),b(10,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_11_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,11),b(11,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_12_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,12),b(12,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_13_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,13),b(13,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_14_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,14),b(14,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_15_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,15),b(15,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_16_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,16),b(16,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_17_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,17),b(17,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_18_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,18),b(18,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_19_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,19),b(19,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_20_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,20),b(20,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_21_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,21),b(21,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_22_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,22),b(22,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_23_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,23),b(23,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmf3_24_f(a,n1,b,n2,c,n3)
c
      real*4 a(n1,24),b(24,n3),c(n1,n3)
c
      do i=1,n1
         do j=1,n3
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
     $             + a(i,24)*b(24,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxm44_f(a,n1,b,n2,c,n3)
C-----------------------------------------------------------------------
C
C     NOTE -- this code has been set up with the "mxmf3" routine
c             referenced in memtime.f.   On most machines, the f2
c             and f3 versions give the same performance (f2 is the
c             nekton standard).  On the t3e, f3 is noticeably faster.
c             pff  10/5/98
C
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C----------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
c
      if (n2.eq.1) then
         call mxm44_2_t_f(a,n1,b,2,c,n3)
      elseif (n2.eq.2) then
         call mxm44_2_t_f(a,n1,b,n2,c,n3)
      else
         call mxm44_0_t_f(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine mxm44_0_t_f(a, m, b, k, c, n)
*      subroutine matmul44_f(m, n, k, a, lda, b, ldb, c, ldc)
*      real*4*8 a(lda,k), b(ldb,n), c(ldc,n)
      real*4 a(m,k), b(k,n), c(m,n)
      real*4 s11, s12, s13, s14, s21, s22, s23, s24
      real*4 s31, s32, s33, s34, s41, s42, s43, s44
c
c matrix multiply with a 4x4 pencil 
c

      mresid = iand(m,3) 
      nresid = iand(n,3) 
      m1 = m - mresid + 1
      n1 = n - nresid + 1

      do i=1,m-mresid,4
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s41 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s42 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          s43 = 0.0d0
          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0
          s44 = 0.0d0
          do l=1,k
            s11 = s11 + a(i,l)*b(l,j)
            s12 = s12 + a(i,l)*b(l,j+1)
            s13 = s13 + a(i,l)*b(l,j+2)
            s14 = s14 + a(i,l)*b(l,j+3)

            s21 = s21 + a(i+1,l)*b(l,j)
            s22 = s22 + a(i+1,l)*b(l,j+1)
            s23 = s23 + a(i+1,l)*b(l,j+2)
            s24 = s24 + a(i+1,l)*b(l,j+3)

            s31 = s31 + a(i+2,l)*b(l,j)
            s32 = s32 + a(i+2,l)*b(l,j+1)
            s33 = s33 + a(i+2,l)*b(l,j+2)
            s34 = s34 + a(i+2,l)*b(l,j+3)

            s41 = s41 + a(i+3,l)*b(l,j)
            s42 = s42 + a(i+3,l)*b(l,j+1)
            s43 = s43 + a(i+3,l)*b(l,j+2)
            s44 = s44 + a(i+3,l)*b(l,j+3)
          enddo
          c(i,j)     = s11 
          c(i,j+1)   = s12 
          c(i,j+2)   = s13
          c(i,j+3)   = s14

          c(i+1,j)   = s21 
          c(i+2,j)   = s31 
          c(i+3,j)   = s41 

          c(i+1,j+1) = s22
          c(i+2,j+1) = s32
          c(i+3,j+1) = s42

          c(i+1,j+2) = s23
          c(i+2,j+2) = s33
          c(i+3,j+2) = s43

          c(i+1,j+3) = s24
          c(i+2,j+3) = s34
          c(i+3,j+3) = s44
        enddo
* Residual when n is not multiple of 4
        if (nresid .ne. 0) then
          if (nresid .eq. 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,n)
              s21 = s21 + a(i+1,l)*b(l,n)
              s31 = s31 + a(i+2,l)*b(l,n)
              s41 = s41 + a(i+3,l)*b(l,n)
            enddo
            c(i,n)     = s11 
            c(i+1,n)   = s21 
            c(i+2,n)   = s31 
            c(i+3,n)   = s41 
          elseif (nresid .eq. 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
          else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)
              s13 = s13 + a(i,l)*b(l,j+2)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)
              s23 = s23 + a(i+1,l)*b(l,j+2)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)
              s33 = s33 + a(i+2,l)*b(l,j+2)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
              s43 = s43 + a(i+3,l)*b(l,j+2)
            enddo
            c(i,j)     = s11 
            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 
            c(i,j+1)   = s12 
            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
            c(i,j+2)   = s13
            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43
          endif
        endif
      enddo

* Residual when m is not multiple of 4
      if (mresid .eq. 0) then
        return
      elseif (mresid .eq. 1) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,j)
            s12 = s12 + a(m,l)*b(l,j+1)
            s13 = s13 + a(m,l)*b(l,j+2)
            s14 = s14 + a(m,l)*b(l,j+3)
          enddo
          c(m,j)     = s11 
          c(m,j+1)   = s12 
          c(m,j+2)   = s13
          c(m,j+3)   = s14
        enddo
* mresid is 1, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n)
          enddo
          c(m,n) = s11
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s12 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-1)
            s12 = s12 + a(m,l)*b(l,n)
          enddo
          c(m,n-1) = s11
          c(m,n) = s12
          return
        else
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-2)
            s12 = s12 + a(m,l)*b(l,n-1)
            s13 = s13 + a(m,l)*b(l,n)
          enddo
          c(m,n-2) = s11
          c(m,n-1) = s12
          c(m,n) = s13
          return
        endif          
      elseif (mresid .eq. 2) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          s21 = 0.0d0
          s22 = 0.0d0
          s23 = 0.0d0
          s24 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,j)
            s12 = s12 + a(m-1,l)*b(l,j+1)
            s13 = s13 + a(m-1,l)*b(l,j+2)
            s14 = s14 + a(m-1,l)*b(l,j+3)

            s21 = s21 + a(m,l)*b(l,j)
            s22 = s22 + a(m,l)*b(l,j+1)
            s23 = s23 + a(m,l)*b(l,j+2)
            s24 = s24 + a(m,l)*b(l,j+3)
          enddo
          c(m-1,j)   = s11 
          c(m-1,j+1) = s12 
          c(m-1,j+2) = s13
          c(m-1,j+3) = s14
          c(m,j)     = s21
          c(m,j+1)   = s22 
          c(m,j+2)   = s23
          c(m,j+3)   = s24
        enddo
* mresid is 2, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n)
          enddo
          c(m-1,n) = s11
          c(m,n) = s21
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-1)
            s12 = s12 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-1)
            s22 = s22 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-1) = s11
          c(m-1,n)   = s12
          c(m,n-1)   = s21
          c(m,n)     = s22
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-2)
            s12 = s12 + a(m-1,l)*b(l,n-1)
            s13 = s13 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-2)
            s22 = s22 + a(m,l)*b(l,n-1)
            s23 = s23 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-2) = s11
          c(m-1,n-1) = s12
          c(m-1,n)   = s13
          c(m,n-2)   = s21
          c(m,n-1)   = s22
          c(m,n)     = s23
          return
        endif
      else
* mresid is 3
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0

          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0

          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0

          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0

          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,j)
            s12 = s12 + a(m-2,l)*b(l,j+1)
            s13 = s13 + a(m-2,l)*b(l,j+2)
            s14 = s14 + a(m-2,l)*b(l,j+3)

            s21 = s21 + a(m-1,l)*b(l,j)
            s22 = s22 + a(m-1,l)*b(l,j+1)
            s23 = s23 + a(m-1,l)*b(l,j+2)
            s24 = s24 + a(m-1,l)*b(l,j+3)

            s31 = s31 + a(m,l)*b(l,j)
            s32 = s32 + a(m,l)*b(l,j+1)
            s33 = s33 + a(m,l)*b(l,j+2)
            s34 = s34 + a(m,l)*b(l,j+3)
          enddo
          c(m-2,j)   = s11 
          c(m-2,j+1) = s12 
          c(m-2,j+2) = s13
          c(m-2,j+3) = s14

          c(m-1,j)   = s21 
          c(m-1,j+1) = s22
          c(m-1,j+2) = s23
          c(m-1,j+3) = s24

          c(m,j)     = s31 
          c(m,j+1)   = s32
          c(m,j+2)   = s33
          c(m,j+3)   = s34
        enddo
* mresid is 3, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n)
          enddo
          c(m-2,n) = s11
          c(m-1,n) = s21
          c(m,n) = s31
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-1)
            s12 = s12 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-1)
            s22 = s22 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-1)
            s32 = s32 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-1) = s11
          c(m-2,n)   = s12
          c(m-1,n-1) = s21
          c(m-1,n)   = s22
          c(m,n-1)   = s31
          c(m,n)     = s32
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-2)
            s12 = s12 + a(m-2,l)*b(l,n-1)
            s13 = s13 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-2)
            s22 = s22 + a(m-1,l)*b(l,n-1)
            s23 = s23 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-2)
            s32 = s32 + a(m,l)*b(l,n-1)
            s33 = s33 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-2) = s11
          c(m-2,n-1) = s12
          c(m-2,n)   = s13
          c(m-1,n-2) = s21
          c(m-1,n-1) = s22
          c(m-1,n)   = s23
          c(m,n-2)   = s31
          c(m,n-1)   = s32
          c(m,n)     = s33
          return
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mxm44_2_t_f(a, m, b, k, c, n)
      real*4 a(m,2), b(2,n), c(m,n)

      nresid = iand(n,3) 
      n1 = n - nresid + 1

      do j=1,n-nresid,4
         do i=1,m
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
            c(i,j+1) = a(i,1)*b(1,j+1)
     $             + a(i,2)*b(2,j+1)
            c(i,j+2) = a(i,1)*b(1,j+2)
     $             + a(i,2)*b(2,j+2)
            c(i,j+3) = a(i,1)*b(1,j+3)
     $             + a(i,2)*b(2,j+3)
         enddo
      enddo
      if (nresid .eq. 0) then
        return
      elseif (nresid .eq. 1) then
         do i=1,m
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      elseif (nresid .eq. 2) then
         do i=1,m
            c(i,n-1) = a(i,1)*b(1,n-1)
     $             + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      else
         do i=1,m
            c(i,n-2) = a(i,1)*b(1,n-2)
     $             + a(i,2)*b(2,n-2)
            c(i,n-1) = a(i,1)*b(1,n-1)
     $             + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mxmtest_f(s,nn,cn,mxmt,name,k,ivb)

      real*4        s(nn,2)   ! MFLOPS
      character*5 cn        ! name
      character*5 name
      external mxmt

      include 'SIZE'
      parameter (lt=4*lx1*ly1*lz1*lelt)
      common /scrns/ a(lt)
      common /scruz/ b(lt)
      common /scrmg/ c(lt)

      integer ll,icalld
      save    ll,icalld
      data    ll,icalld  /1,0/

      if (icalld.eq.0) then    !     Initialize matrices:
         icalld = icalld + 1
         time1 = dnekclock()
         call initab_f(a,b,lt)
         time2 = dnekclock()-time1
         if (nid.eq.0) write(6,*) 'mxm test init:',lt,time2,name
      endif


      cn = name

c     Rectangular matrix tests

      nn0 = 1
      nn1 = nn
      if (ivb.eq.0) then
         nn0 = lx1
         nn1 = lx1
      endif

      m = k
      do n=nn0,nn1
         n1 = n
         n2 = n
         n3 = n
         if (m.eq.1) n1 = n*n
         if (m.eq.3) n3 = n*n
         if (lt.gt.n1*n3) then
          n13 = max(n1,n3)
          loop = 250000/(n1*n2*n3) + 500
          if (name.eq.'madd ') loop = 200000/(n1*n3) + 5000

c-------------------------------------------------------
c         mem test
c-------------------------------------------------------

          t0    = dnekclock()
          overh = dnekclock()-t0
          time1 = dnekclock()
          do l=1,loop
            if (ll.ge.lt-n1*n3) ll = 1
            call mxmt(a(ll),n1,b(ll),n2,c(ll),n3)
            ll = ll+n1*n3
          enddo
          time2 = dnekclock()
          time = time2-time1 - overh
          iops=loop*n1*n3*(2*n2-1)
          if (name.eq.'madd ') iops = loop*n1*n3
c         write(6,*) loop,time,time2,time1,overh
          flops=iops/(1.0e6*time)
          s(n,1) = flops
c
          timel = time/loop
          if (nid.eq.0) write(6,199) n,n1,n2,n3,flops,timel,name
  199     format(i3,'m',1x,3i6,f10.4,e16.5,3x,a5,' mem')
c
c-------------------------------------------------------
c         fast test
c-------------------------------------------------------
c
          call mxmt(a,n1,b,n2,c,n3)
          t0    = dnekclock()
          overh = dnekclock()-t0
          time1 = dnekclock()
          do l=1,loop
            call mxmt(a,n1,b,n2,c,n3)
          enddo
          time2 = dnekclock()
          time = time2-time1 - overh
          iops=loop*n1*n3*(2*n2-1)
          if (name.eq.'madd ') iops = loop*n1*n3
          flops=iops/(1.0e6*time)
          s(n,2) = flops
          timel = time/loop
c
          if (nid.eq.0) write(6,198) n,n1,n2,n3,flops,timel,name
  198     format(i3,'f',1x,3i6,f10.4,e16.5,3x,a5,' fast')
c
        endif
       enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mxm_analyze_f(s,a,nn,c,nt,ivb)
      include 'SIZE'

      character*5 c(3,nt)
      real*4        s(nn,2,nt,3)  ! Measured Mflops, 3 cases
      real*4        a(nn,2,nt,3)
c                   ^  ^ ^  |__ N^2xN, NxN, NxN^2
c  matrix order N __|  | |__________which mxm
c                      |
c                      |__cached vs. noncached data
 

      integer itmax(200)

      nn0 = 1
      nn1 = nn
      if (ivb.eq.0) then
         nn0 = lx1
         nn1 = lx1
      endif

      do n = nn0,nn1
         fmax = 0.   ! Peak mflops
         do it=1,nt
            ai = 0.
            di = 0.
            do k=1,3
               if (s(n,1,it,k).gt.0) then  ! Take harmonic means of
                  ai = ai + 1./s(n,1,it,k) ! case I II and III for 
                  di = di + 1.             ! mem test, s(n,1...).
               endif
            enddo
            if (ai.gt.0) ai = di/ai
            a(n,1,it,1) = di/ai
            if (ai.gt.fmax.and.c(2,it).ne.'madd ') then
               fmax     = ai
               itmax(n) = it
            endif
         enddo
         it = itmax(n)
         if (nid.eq.0) write(6,3) n,it,c(2,it),(s(n,1,it,k),k=1,3),fmax
    3    format(i3,i2,1x,a5,4f12.0,' Peak harmonic')
      enddo
      call out_anal_f(s,a,nn,c,nt,itmax,'Harmonic',1,ivb)
c
c     Case by case
c
      do k=1,3
         do n = nn0,nn1
            fmax = 0.   ! Peak mflops
            do it=1,nt
               ai = s(n,1,it,k)
               if (ai.gt.fmax.and.c(2,it).ne.'madd ') then
                  fmax     = ai
                  itmax(n) = it
               endif
            enddo
         enddo
         if (k.eq.1) call out_anal_f(s,a,nn,c,nt,itmax,'Case N2N',k,ivb)
         if (k.eq.2) call out_anal_f(s,a,nn,c,nt,itmax,'Case NxN',k,ivb)
         if (k.eq.3) call out_anal_f(s,a,nn,c,nt,itmax,'Case NN2',k,ivb)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine out_anal_f(s,a,nn,c,nt,itmax,name8,k,ivb)
      include 'SIZE'

      character*5 c(3,nt)
      real*4        s(nn,2,nt,3)
      real*4        a(nn,2,nt,3)
      integer itmax(200)
      character*8 name8

      if (nid.ne.0) return

      nn0 = 1
      nn1 = nn
      if (ivb.eq.0) then
         nn0 = lx1
         nn1 = lx1
      endif


      do n=nn0,nn1
         it = itmax(n)
         write(6,1) n,s(n,1,it,k),c(2,it),name8
    1    format(i4,f14.0,4x,a5,4x,a8,'   MxM MFLOPS')
      enddo

      return
      end
c-----------------------------------------------------------------------
