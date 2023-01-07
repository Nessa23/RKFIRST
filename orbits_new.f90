program orb
implicit none
integer, parameter::sizek=20000
real(16) const(sizek),ini(sizek),dzeta(sizek),energ(sizek),ep1t(sizek),eps1t(sizek)
real(16) a
integer i,num
!integer, parameter :: dp = selected_real_kind(p = 15, r = 307)
integer, parameter::P=10100000 
real(16) ,dimension(P)::T,X1,X2,V1,V2,energy,delta,error
real(16)::Ti,Tf,X10,X20,V10,V20
integer::Nt,k
!real(dp)::array(1000)
!C(1000)=(/(i, i=0,20,0.1)/)
real(16):: omega=0.02829174,G=0.004535,alpha1=-0.001976995,m1,m2,mu
real(16)::ep1tt=10.44666776
!real(8)::eps1t
real(16)::e1(sizek)
real(16):: Xi,period
integer::n,d,j
m2=499.0
m1=1.0
mu=1+(m1/m2)
!eps1t=-3./2.*G*m2/(mu**2.*ep1t)!energy(0)
!e1=(0.45*1.5*G*m2)/(2.*mu**2.*ep1t)
a = 0.001 ! Initial value
const= (/((i*a),i=1,sizek)/)!!!CONSTANTS==========================

!ini=-10.44666776+const
!print*,'const=',ini
const = const + 1.5D0
do k=1,sizek
!==================================INITIAL CONDITIONS=======================
ini(k)=const(k)-10.44666776
X20=0
V10=0
ep1t(k)=const(k)-10.44666776

eps1t(k)=-3./2.*G*m2/(mu**2.*ep1tt)!energy(0)
e1(k)= 0.44D0 * eps1t(k)
V20=SQRT(2.*(e1(k)+(G*m2)/(mu**2.*(ABS(ini(k)))))-alpha1*ini(k)**2.)
!============================================================================
Nt=100000
Ti=0
Tf=300
if(Nt.gt.P) stop 'Nt>P'
call RK(T,X1,X2,V1,V2,Ti,Tf,ini(k),X20,V10,V20,Nt)

!-------------------------------------------------------------------------------------------------------------------
open(unit=11,STATUS='REPLACE', file='rk2v.dat')
open(unit=14,STATUS='REPLACE', file='error1.dat')
do i=1,Nt
 energy(i)=(V1(i)**2.+V2(i)**2.)/2.-G*m2/(mu**2.*((X1(i)**2.+X2(i)**2.)**(1./2)))+alpha1*X1(i)**2./2.
 error(i)=ABS((energy(i)-e1(k))/e1(k))
 
 !print* , 'V1(Ti)=',energy(i)
 
 write (11 ,*)T(i),X1(i) ,X2(i) ,V1(i),V2(i)
  write (14 ,*)T(i),error(i)

!==================CONDITION CHECK====================================================================
if (X2(i)*X2(i+1)<0 .and. V2(i)>=0) then
 !print* , 'period=',T(i)
 !print* , 'i=',i
 !print* , 'X2I-1=',X2(i) 
 !print* , 'X2I=',X2(i+1)

 
period=T(i)+(T(i+1)-T(i))/(X2(i+1)-X2(i))*(0.-X2(i))
 Xi=X1(i)+(X1(i+1)-X1(i))/(X2(i+1)-X2(i))*(0.-X2(i))!that's the interpolation

!open(unit=12, file='const5.dat')
!open(unit=13,file='error.dat')
 !dzeta(k)=DABS(Xi-ini(k))
 !print*,'delta=',dzeta
 !energ(i)=(V1(i)**2.+V2(i)**2.)/2.-G*m2/(mu**2.*((Xi**2.+X2(i)**2.)**(1./2.)))+(alpha1*Xi**2.)/2.!Great question: is there supposed to be Xi or X1(i)
 !eps1t=-3./2.*G*m2/(mu**2.*ep1t)!energy(0)
!e1=(0.223*1.5*G*m2)/(2.*mu**2.*ep1t)
 dzeta(k)=ABS((Xi-ini(k)) / ini(k))
 
 !print*,'error=',error(i)
 write (12 ,*)dzeta(k),const(k)
 
 if (dzeta(k) < 1.D-5) then
  print* , 'orbit found ',const(k)
  print* , 'period = ',period
  stop
 end if
 !That's the check up 
 !print* , 'warning large error for energy for ',const(k)
 ! end if
 end if
 
 !write(13,*)error(i),T(i)
!if (error(i)>=dzeta(k))then!That's the check up 
 !print* , 'const=',const(k)
  !end if
 !end if

enddo

!if (MAX1(error)>1.D-10) then
if (MAXVAL(error)>1.D-10) then!That's the check up 
 print* , 'warning: large error for energy for ', MAXVAL(error), const(k)
end if
close(11)
!close(12)
!close(13)
close(14)
enddo

end program orb
!---------------------------------------------------FUNCTIONS--DERIVIATIVES-------------------------------
real(16) function f1(t,x1,x2,v1,v2)
 implicit none
 real(16)::t,x1,x2,v1,v2
 f1=v1 !dx1/dt= v = x2
end function f1 
!----------------------------------------------------------------------------------
real(16) function f2(t,x1,x2,v1,v2)
 implicit none
 real(16)::t,x1,x2,v1,v2
 f2=v2 
end function f2
!----------------------------------------------------------------------------------
real(16) function f3(t,x1,x2,v1,v2)
implicit none
real:: omega=0.02829174,G=0.004535,alpha1=-0.001976995,m1,m2,mu
real(16)::t,x1,x2,v1,v2
m2=499.0
m1=1.0
mu=1+(m1/m2)
 f3=2.*omega*v2-x1*((G*m2)/(mu**2.*((x1**2.+x2**2.)**(3./2.)))+alpha1)
 !f3=2.*omega*v2-x1*((G*m2)/(mu**2.*(x1**2.+x2**2.)**(3./2.))+alpha1)
end function f3
!----------------------------------------------------------------------------------
real(16) function f4(t,x1,x2,v1,v2)
implicit none
real(16)::t,x1,x2,v1,v2
real:: omega=0.02829174,G=0.004535,alpha1=-0.001976995,m1,m2,mu
m2=499.0
m1=1.0
mu=1+(m1/m2)
 !f4=-2*omega*v1-x2*((G*m2)/(mu**2.*((x1**2.+x2**2.)**(3./2.))))
 f4=-2.*omega*v1-x2*G*m2/(mu**2.*sqrt(x1**2.+x2**2.)**3.)
end function f4
!−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
subroutine RK(T,X1,X2,V1,V2,Ti,Tf,X10,X20,V10,V20,Nt)
implicit none
integer :: Nt
real(16) ,dimension(Nt):: T,X1,X2,V1,V2
real(16) :: Ti,Tf,X10,X20,V10,V20
real(16) :: dt
real(16) :: TS,X1S,X2S,V1S,V2S

integer ::i
dt=(Tf-Ti)/(Nt-1)
T(1)=Ti
X1(1)=X10; X2(1)=X20
V1(1)=V10; V2(1)=V20
TS=Ti
X1S=X10; X2S= X20
V1S=V10; V2S= V20

do i=2,Nt
 call RKSTEP(TS,X1S,X2S,V1S,V2S,dt)
 T(i)=TS
 X1(i)=X1S; X2(i)=X2S
 V1(i)=V1S; V2(i)=V2S
enddo
end subroutine RK
! ======================================================== 

!======================================================================
subroutine RKSTEP(t,x1,x2,x3,x4,dt)
implicit none
real(16)::t,x1,x2,x3,x4,dt
real(16)::f1,f2,f3,f4
real(16)::k11,k12,k13,k14,k21,k22,k23,k24
real(16)::k31,k32,k33,k34,k41,k42,k43,k44
real(16)::h,h2,h6
h =dt !h =dt , integration step h2 =0.5D0*h !h2=h/2
h2=0.5*h
h6 =h/6.0 !h6=h/6
k11=f1(t,x1,x2,x3,x4)
k21=f2(t,x1,x2,x3,x4)
k31=f3(t,x1,x2,x3,x4)
k41=f4(t,x1,x2,x3,x4)
k12=f1 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41 )
k22=f2 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41 )
k32=f3 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41 )
k42=f4 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41 )
k13=f1 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42 )
k23=f2 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42 )
k33=f3 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42 )
k43=f4 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42 )
k14=f1 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43 )
k24=f2 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43 )
k34=f3 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43 )
k44=f4 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43 )
t =t+h
x1 =x1+h6*(k11+2.0D0*(k12+k13)+k14)
x2 =x2+h6*(k21+2.0D0*(k22+k23)+k24)
x3=x3+h6 * ( k31+2.0D0*( k32+k33 )+k34 )
x4=x4+h6 * ( k41+2.0D0 *( k42+k43 )+k44 )
end subroutine RKSTEP
!=====================================================================

