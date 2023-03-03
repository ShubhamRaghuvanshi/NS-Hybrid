
            subroutine dermat_lag(intx,inty,intz,intxp1,intyp1,intzp1,xd,yd,zd)

        use mod_serial_fluid
        implicit none
        integer,intent(in)::intx,inty,intz,intxp1,intyp1,intzp1
        real*8,intent(in)::xd,yd,zd
        real*8,dimension(3,3)::gradu_p,gradu_p2,gradu_p3
        real*8::divl,Ql,Rl

!**********************Interpolating derivative matrix**************************

      call trilinear(xd,yd,zd,V1_x1(intx,inty,intz),V1_x1(intxp1,inty,intz),&
                      V1_x1(intx,intyp1,intz),V1_x1(intxp1,intyp1,intz),&
                      V1_x1(intx,inty,intzp1),V1_x1(intxp1,inty,intzp1),&
          V1_x1(intx,intyp1,intzp1),V1_x1(intxp1,intyp1,intzp1),gradu_p(1,1))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V1_x2(intx,inty,intz),V1_x2(intxp1,inty,intz),&
                      V1_x2(intx,intyp1,intz),V1_x2(intxp1,intyp1,intz),&
                      V1_x2(intx,inty,intzp1),V1_x2(intxp1,inty,intzp1),&
          V1_x2(intx,intyp1,intzp1),V1_x2(intxp1,intyp1,intzp1),gradu_p(2,1))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V1_x3(intx,inty,intz),V1_x3(intxp1,inty,intz),&
                      V1_x3(intx,intyp1,intz),V1_x3(intxp1,intyp1,intz),&
                      V1_x3(intx,inty,intzp1),V1_x3(intxp1,inty,intzp1),&
          V1_x3(intx,intyp1,intzp1),V1_x3(intxp1,intyp1,intzp1),gradu_p(3,1))

!-----------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V2_x1(intx,inty,intz),V2_x1(intxp1,inty,intz),&
                      V2_x1(intx,intyp1,intz),V2_x1(intxp1,intyp1,intz),&
                      V2_x1(intx,inty,intzp1),V2_x1(intxp1,inty,intzp1),&
          V2_x1(intx,intyp1,intzp1),V2_x1(intxp1,intyp1,intzp1),gradu_p(1,2))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V2_x2(intx,inty,intz),V2_x2(intxp1,inty,intz),&
                      V2_x2(intx,intyp1,intz),V2_x2(intxp1,intyp1,intz),&
                      V2_x2(intx,inty,intzp1),V2_x2(intxp1,inty,intzp1),&
          V2_x2(intx,intyp1,intzp1),V2_x2(intxp1,intyp1,intzp1),gradu_p(2,2))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V2_x3(intx,inty,intz),V2_x3(intxp1,inty,intz),&
                      V2_x3(intx,intyp1,intz),V2_x3(intxp1,intyp1,intz),&
                      V2_x3(intx,inty,intzp1),V2_x3(intxp1,inty,intzp1),&
          V2_x3(intx,intyp1,intzp1),V2_x3(intxp1,intyp1,intzp1),gradu_p(3,2))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V3_x1(intx,inty,intz),V3_x1(intxp1,inty,intz),&
                      V3_x1(intx,intyp1,intz),V3_x1(intxp1,intyp1,intz),&
                      V3_x1(intx,inty,intzp1),V3_x1(intxp1,inty,intzp1),&
          V3_x1(intx,intyp1,intzp1),V3_x1(intxp1,intyp1,intzp1),gradu_p(1,3))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V3_x2(intx,inty,intz),V3_x2(intxp1,inty,intz),&
                      V3_x2(intx,intyp1,intz),V3_x2(intxp1,intyp1,intz),&
                      V3_x2(intx,inty,intzp1),V3_x2(intxp1,inty,intzp1),&
          V3_x2(intx,intyp1,intzp1),V3_x2(intxp1,intyp1,intzp1),gradu_p(2,3))

!------------------------------------------------------------------------------

      call trilinear(xd,yd,zd,V3_x3(intx,inty,intz),V3_x3(intxp1,inty,intz),&
                      V3_x3(intx,intyp1,intz),V3_x3(intxp1,intyp1,intz),&
                      V3_x3(intx,inty,intzp1),V3_x3(intxp1,inty,intzp1),&
          V3_x3(intx,intyp1,intzp1),V3_x3(intxp1,intyp1,intzp1),gradu_p(3,3))

 
         gradu_p2 = matmul(gradu_p,gradu_p)
         gradu_p3 = matmul(gradu_p,gradu_p2)

         divl = gradu_p(1,1)+gradu_p(2,2)+gradu_p(3,3)
         Ql = -0.5d0*(gradu_p2(1,1)+gradu_p2(2,2)+gradu_p2(3,3))
         Rl = -(1.0d0/3.0d0)*(gradu_p3(1,1)+gradu_p3(2,2)+gradu_p3(3,3))

         write(14,*)divl,Ql,Rl
                         end subroutine dermat_lag
