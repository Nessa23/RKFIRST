program orb
implicit none
integer, parameter::sizek=20000
real(16) const(sizek),ini(sizek),dzeta(sizek),energ(sizek),eps1t
real(16) a
integer i,num
!integer, parameter :: dp = selected_real_kind(p = 15, r = 307)
integer, parameter::P=10100000 
real(16) ,dimension(P)::T,X1,X2,V1,V2,energy,delta,error,DX1,DX2,DV1,DV2,epsilonn,R1
real(16)::Ti,Tf,X10,X20,V10,V20,DX10,DX20,DV10,DV20
integer::Nt,k
!real(dp)::array(1000)
!C(1000)=(/(i, i=0,20,0.1)/)
real(16):: omega=0.02829174,G=0.004535,alpha1=-0.001976995,m1,m2,mu
real(16)::ep1tt=10.44666776
!real(8)::eps1t
real(16)::e1,b
real(16):: Xi,period,INITIAL,ep1t,DIFFER,R0
integer::n,d,j
m2=499.0
m1=1.0
mu=1+(m1/m2)
!eps1t=-3./2.*G*m2/(mu**2.*ep1t)!energy(0)
!e1=(0.45*1.5*G*m2)/(2.*mu**2.*ep1t)
a = 0.001 ! Initial value
const= (/((i*a),i=1,sizek)/)!!!CONSTANTS==========================
b=G*m2/mu**2
!ini=-10.44666776+const
!print*,'const=',ini
const = const + 1.4D0
!do k=1,sizek
!==================================INITIAL CONDITIONS=======================
INITIAL=1.5340000253636390-10.44666776
X20=0
V10=0
ep1t=1.5340000253636390-10.44666776
DX10=1
DX20=0
DV10=0
DV20=0
eps1t=-3./2.*G*m2/(mu**2.*ep1tt)!energy(0)
e1= 0.44D0 * eps1t
V20=SQRT(2.*(e1+(G*m2)/(mu**2.*(ABS(INITIAL))))-alpha1*INITIAL**2.)
!============================================================================
Nt=100000
Ti=0
Tf=212.37421593286641

if(Nt.gt.P) stop 'Nt>P'
call RK(T,X1,X2,V1,V2,DX1,DX2,DV1,DV2,Ti,Tf,INITIAL,X20,V10,V20,DX10,DX20,DV10,DV20,Nt)

print*,DX1(Nt),DX2(Nt),DV1(Nt),DV2(Nt)
!-------------------------------------------------------------------------------------------------------------------
open(unit=11 ,file='FILE1.dat')
open(unit=12 ,file='ERROR1.dat')

do i=1,Nt
R1(i)=sqrt(X1(i)**2+X2(i)**2)
energy(i)=V1(i)*DV1(i)+V2(i)*DV2(i)+(alpha1+b/(R1(i)**3))*X1(i)*DX1(i)+b*X2(i)*DX2(i)/(R1(i)**3)
!print*,energy(i)
R0=SQRT(INITIAL**2+X20**2)
DIFFER=V10*DV10+V20*DV20+(alpha1+b/(R0**3))*INITIAL*DX10+b*X20*DX20/(R0**3)

epsilonn(i)=ABS((energy(i)-DIFFER)/DIFFER)

 write (11 ,*)T(i),X1(i) ,X2(i) ,V1(i),V2(i)
 write (12 ,*)T(i),epsilonn(i)
  
enddo
close(11)
close(12)
!--------------
!==================================INITIAL CONDITIONS=======================
INITIAL=1.5340000253636390-10.44666776
X20=0
V10=0
ep1t=1.5340000253636390-10.44666776
DX10=0
DX20=1
DV10=0
DV20=0
eps1t=-3./2.*G*m2/(mu**2.*ep1tt)!energy(0)
e1= 0.44D0 * eps1t
V20=SQRT(2.*(e1+(G*m2)/(mu**2.*(ABS(INITIAL))))-alpha1*INITIAL**2.)
!============================================================================
Nt=100000
Ti=0
Tf=212.37421593286641

if(Nt.gt.P) stop 'Nt>P'
call RK(T,X1,X2,V1,V2,DX1,DX2,DV1,DV2,Ti,Tf,INITIAL,X20,V10,V20,DX10,DX20,DV10,DV20,Nt)

print*,DX1(Nt),DX2(Nt),DV1(Nt),DV2(Nt)
!-------------------------------------------------------------------------------------------------------------------
open(unit=13 ,file='FILE2.dat')
open(unit=14 ,file='ERROR2.dat')

do i=1,Nt
R1(i)=sqrt(X1(i)**2+X2(i)**2)
energy(i)=V1(i)*DV1(i)+V2(i)*DV2(i)+(alpha1+b/(R1(i)**3))*X1(i)*DX1(i)+b*X2(i)*DX2(i)/(R1(i)**3)
!print*,energy(i)
R0=SQRT(INITIAL**2+X20**2)
DIFFER=V10*DV10+V20*DV20+(alpha1+b/(R0**3))*INITIAL*DX10+b*X20*DX20/(R0**3)

epsilonn(i)=ABS((energy(i)-DIFFER))

 write (13 ,*)T(i),X1(i) ,X2(i) ,V1(i),V2(i)
 write (14 ,*)T(i),epsilonn(i)
  
enddo
close(13)
close(14)
!==================================INITIAL CONDITIONS=======================
INITIAL=1.5340000253636390-10.44666776
X20=0
V10=0
ep1t=1.5340000253636390-10.44666776
DX10=0
DX20=0
DV10=1
DV20=0
eps1t=-3./2.*G*m2/(mu**2.*ep1tt)!energy(0)
e1= 0.44D0 * eps1t
V20=SQRT(2.*(e1+(G*m2)/(mu**2.*(ABS(INITIAL))))-alpha1*INITIAL**2.)
!============================================================================
Nt=100000
Ti=0
Tf=212.37421593286641

if(Nt.gt.P) stop 'Nt>P'
call RK(T,X1,X2,V1,V2,DX1,DX2,DV1,DV2,Ti,Tf,INITIAL,X20,V10,V20,DX10,DX20,DV10,DV20,Nt)

print*,DX1(Nt),DX2(Nt),DV1(Nt),DV2(Nt)
!-------------------------------------------------------------------------------------------------------------------
open(unit=15 ,file='FILE3.dat')
open(unit=16 ,file='ERROR3.dat')

do i=1,Nt
R1(i)=sqrt(X1(i)**2+X2(i)**2)
energy(i)=V1(i)*DV1(i)+V2(i)*DV2(i)+(alpha1+b/(R1(i)**3))*X1(i)*DX1(i)+b*X2(i)*DX2(i)/(R1(i)**3)
!print*,energy(i)
R0=SQRT(INITIAL**2+X20**2)
DIFFER=V10*DV10+V20*DV20+(alpha1+b/(R0**3))*INITIAL*DX10+b*X20*DX20/(R0**3)

epsilonn(i)=ABS((energy(i)-DIFFER))

 write (15 ,*)T(i),X1(i) ,X2(i) ,V1(i),V2(i)
 write (16 ,*)T(i),epsilonn(i)
  
enddo
close(15)
close(16)
!==================================INITIAL CONDITIONS=======================
INITIAL=1.5340000253636390-10.44666776
X20=0
V10=0
ep1t=1.5340000253636390-10.44666776
DX10=0
DX20=0
DV10=0
DV20=1
eps1t=-3./2.*G*m2/(mu**2.*ep1tt)!energy(0)
e1= 0.44D0 * eps1t
V20=SQRT(2.*(e1+(G*m2)/(mu**2.*(ABS(INITIAL))))-alpha1*INITIAL**2.)
!============================================================================
Nt=100000
Ti=0
Tf=212.37421593286641

if(Nt.gt.P) stop 'Nt>P'
call RK(T,X1,X2,V1,V2,DX1,DX2,DV1,DV2,Ti,Tf,INITIAL,X20,V10,V20,DX10,DX20,DV10,DV20,Nt)

print*,DX1(Nt),DX2(Nt),DV1(Nt),DV2(Nt)
!-------------------------------------------------------------------------------------------------------------------
open(unit=17 ,file='FILE4.dat')
open(unit=18 ,file='ERROR4.dat')

do i=1,Nt
R1(i)=sqrt(X1(i)**2+X2(i)**2)
energy(i)=V1(i)*DV1(i)+V2(i)*DV2(i)+(alpha1+b/(R1(i)**3))*X1(i)*DX1(i)+b*X2(i)*DX2(i)/(R1(i)**3)
!print*,energy(i)
R0=SQRT(INITIAL**2+X20**2)
DIFFER=V10*DV10+V20*DV20+(alpha1+b/(R0**3))*INITIAL*DX10+b*X20*DX20/(R0**3)

epsilonn(i)=ABS((energy(i)-DIFFER)/DIFFER)

 write (17 ,*)T(i),X1(i) ,X2(i) ,V1(i),V2(i)
 write (18 ,*)T(i),epsilonn(i)
  
enddo
close(11)
close(12)


end program orb
!---------------------------------------------------FUNCTIONS--DERIVIATIVES-------------------------------
real(16) function f1(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
 implicit none
 real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
 f1=v1 !dx1/dt= v = x2
end function f1 
!----------------------------------------------------------------------------------
real(16) function f2(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
 implicit none
 real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
 f2=v2 
end function f2

real(16) function f5(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
 implicit none
 real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
 f5=dv1
end function f5
!----------------------------------------------------------------------------------
real(16) function f6(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
 implicit none
 real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
 f6=dv2
end function f6
!----------------------------------------------------------------------------
real(16) function f3(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
implicit none
real:: omega=0.02829174,G=0.004535,alpha1=-0.001976995,m1,m2,mu
real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
m2=499.0
m1=1.0
mu=1+(m1/m2)
 f3=2.*omega*v2-x1*((G*m2)/(mu**2.*((x1**2.+x2**2.)**(3./2.)))+alpha1)
 !f3=2.*omega*v2-x1*((G*m2)/(mu**2.*(x1**2.+x2**2.)**(3./2.))+alpha1)
end function f3
!----------------------------------------------------------------------------------
real(16) function f4(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
implicit none
real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
real:: omega=0.02829174,G=0.004535,alpha1=-0.001976995,m1,m2,mu
m2=499.0
m1=1.0
mu=1+(m1/m2)
 !f4=-2*omega*v1-x2*((G*m2)/(mu**2.*((x1**2.+x2**2.)**(3./2.))))
 f4=-2.*omega*v1-x2*G*m2/(mu**2.*sqrt(x1**2.+x2**2.)**3.)
end function f4

real(16) function f7(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
implicit none
real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
real:: omega=0.02829174,G=0.004535,var,alpha1=-0.001976995,m1,m2,mu,beta,r1
m2=499.0
m1=1.0
mu=1+(m1/m2)
r1=sqrt(x1**2+x2**2)
beta=G*m2/(mu**2*r1**3)
 !f4=-2*omega*v1-x2*((G*m2)/(mu**2.*((x1**2.+x2**2.)**(3./2.))))
 !f7=-2.*omega*v1-x2*G*m2/(mu**2.*sqrt(x1**2.+x2**2.)**3.)
var=beta*(1-3*x1**2/r1**2)
f7=2*omega*dv2-(var+alpha1)*dx1+3*beta*x1*x2*dx2/r1**2
end function f7

real(16) function f8(t,x1,x2,v1,v2,dx1,dx2,dv1,dv2)
implicit none
real(16)::t,x1,x2,v1,v2,dx1,dx2,dv1,dv2
real:: omega=0.02829174,G=0.004535,var,alpha1=-0.001976995,m1,m2,mu,beta,r1
m2=499.0
m1=1.0
mu=1+(m1/m2)
r1=sqrt(x1**2+x2**2)
beta=G*m2/(mu**2*r1**3)
 !f4=-2*omega*v1-x2*((G*m2)/(mu**2.*((x1**2.+x2**2.)**(3./2.))))
 var=beta*(1-3*x2**2/r1**2)
 f8=-2*omega*dv1+3*beta*x1*x2*dx1/r1**2-var*dx2
end function f8
!−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
subroutine RK(T,X1,X2,V1,V2,DX1,DX2,DV1,DV2,Ti,Tf,X10,X20,V10,V20,DX10,DX20,DV10,DV20,Nt)
implicit none
integer :: Nt
real(16) ,dimension(Nt):: T,X1,X2,V1,V2,DX1,DX2,DV1,DV2
real(16) :: Ti,Tf,X10,X20,V10,V20,DX10,DX20,DV10,DV20
real(16) :: dt
real(16) :: TS,X1S,X2S,V1S,V2S,DX1S,DX2S,DV1S,DV2S

integer ::i
dt=(Tf-Ti)/(Nt-1)
T(1)=Ti
X1(1)=X10; X2(1)=X20
V1(1)=V10; V2(1)=V20
DX1(1)=DX10;DX2(1)=DX20
TS=Ti
X1S=X10; X2S= X20
V1S=V10; V2S= V20
DX1S=DX10;DX2S=DX20
DV1S=DV10;DV2S=DV20

do i=2,Nt
 call RKSTEP(TS,X1S,X2S,V1S,V2S,DX1S,DX2S,DV1S,DV2S,dt)
 T(i)=TS
 X1(i)=X1S; X2(i)=X2S
 V1(i)=V1S; V2(i)=V2S
 DX1(i)=DX1S;DX2(i)=DX2S
 DV1(i)=DV1S;DV2(i)=DV2S
enddo
end subroutine RK
! ======================================================== 

!======================================================================
subroutine RKSTEP(t,x1,x2,x3,x4,x5,x6,x7,x8,dt)
implicit none
real(16)::t,x1,x2,x3,x4,x5,x6,x7,x8,dt
real(16)::f1,f2,f3,f4,f5,f6,f7,f8
real(16)::k11,k12,k13,k14,k21,k22,k23,k24
real(16)::k31,k32,k33,k34,k41,k42,k43,k44,k51,k61,k52,k62,k63,k53,k64,k54,k71,k72,k73,k74,k81,k82,k83,k84
real(16)::h,h2,h6
h =dt !h =dt , integration step h2 =0.5D0*h !h2=h/2
h2=0.5D0*h
h6 =h/6.0 !h6=h/6
k11=f1(t,x1,x2,x3,x4,x5,x6,x7,x8)
k21=f2(t,x1,x2,x3,x4,x5,x6,x7,x8)
k31=f3(t,x1,x2,x3,x4,x5,x6,x7,x8)
k41=f4(t,x1,x2,x3,x4,x5,x6,x7,x8)
k51=f5(t,x1,x2,x3,x4,x5,x6,x7,x8)
k61=f6(t,x1,x2,x3,x4,x5,x6,x7,x8)
k71=f7(t,x1,x2,x3,x4,x5,x6,x7,x8)
k81=f8(t,x1,x2,x3,x4,x5,x6,x7,x8)


k12=f1 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41,x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81 )
k22=f2 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41,x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81 )
k32=f3 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41,x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81)
k42=f4 ( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41, x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81)
k52=f5( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41, x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81)
k62=f6( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41, x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81)
k72=f7( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41, x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81)
k82=f8( t+h2 , x1+h2 * k11 , x2+h2 * k21 , x3+h2 * k31 , x4+h2 * k41, x5+h2 * k51,x6+h2 * k61,x7+h2 * k71,x8+h2 * k81)

k13=f1 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42,x5+h2 * k52,x6+h2 * k62,x7+h2 * k72,x8+h2 * k82 )
k23=f2 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42,x5+h2 * k52,x6+h2 * k62 ,x7+h2 * k72,x8+h2 * k82 )
k33=f3 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42,x5+h2 * k52,x6+h2 * k62 ,x7+h2 * k72,x8+h2 * k82 )
k43=f4 ( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42,x5+h2 * k52,x6+h2 * k62 ,x7+h2 * k72,x8+h2 * k82 )
k53=f5( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42, x5+h2 * k52,x6+h2 * k62,x7+h2 * k72,x8+h2 * k82 )
k63=f6( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42, x5+h2 * k52,x6+h2 * k62,x7+h2 * k72,x8+h2 * k82 )
k73=f7( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42, x5+h2 * k52,x6+h2 * k62,x7+h2 * k72,x8+h2 * k82)
k83=f8( t+h2 , x1+h2 * k12 , x2+h2 * k22 , x3+h2 * k32 , x4+h2 * k42, x5+h2 * k52,x6+h2 * k62,x7+h2 * k72,x8+h2 * k82)

k14=f1 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43,x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83 )
k24=f2 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43,x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83)
k34=f3 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43,x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83 )
k44=f4 ( t+h , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43,x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83 )
k54=f5( t+h2 , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43, x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83)
k64=f6( t+h2 , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43, x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83)
k74=f7( t+h2 , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43, x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83)
k84=f8( t+h2 , x1+h * k13 , x2+h * k23 , x3+h * k33 , x4+h * k43, x5+h * k53,x6+h * k63,x7+h * k73,x8+h * k83)

t =t+h
x1 =x1+h6*(k11+2.0D0*(k12+k13)+k14)
x2 =x2+h6*(k21+2.0D0*(k22+k23)+k24)
x3=x3+h6 * ( k31+2.0D0*( k32+k33 )+k34 )
x4=x4+h6 * ( k41+2.0D0 *( k42+k43 )+k44 )
x5=x5+h6 * ( k51+2.0D0 *( k52+k53 )+k54 )
x6=x6+h6 * ( k61+2.0D0 *( k62+k63 )+k64 )
x7=x7+h6 * ( k71+2.0D0 *( k72+k73 )+k74 )
x8=x8+h6 * ( k81+2.0D0 *( k82+k83 )+k84 )
end subroutine RKSTEP
!=====================================================================

