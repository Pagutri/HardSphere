!factor de estructura de una mezcla
! 
program baxt
implicit none
Integer, Parameter      :: nesp= 2, nkmax = 2**12
Real(8), Parameter      :: pi= 3.1416d0, k0 = .01d0, dk = 0.01d0, alfa = 1.545d0, fact =1.d0/(6.d0*(pi**2))
Integer			:: nk, l, it,imolar,iphif, fip
Real(8), Dimension(nesp,nesp,nkmax)	:: Z, ZI, lamdaI, sqrtnh, csqrtn
Real(8)			:: k, a(nesp), b(nesp), psi(3), tol, large, tolF, ethatotal, x1, dx
Real(8)			::xch,xgr,phitotal,phi_mix
Real(8)			::ethamin, ethamax, errorF
Real(8)			:: phi_fluid, phi_glass, etha1, etha2
Real(8), Dimension(nesp)  	::  R, kc, suma, rho, error, etha, fm
Real(8), Dimension(nesp,nesp)   ::  IM, S, gama, Rhoa, RhoaI, one, oneI, two, twoI, three
Real(8), Dimension(nesp, 10000) :: gamalist
integer, dimension(nesp)	::flag
Real(8) :: rk1,rk2,rk3,rk4

R=(/1.d0, 0.5d0/) !! se definen los diametros
kc= (2.d0*pi*alfa)/R

!se define la matriz identidad
IM= 0.0d0
do l=1, nesp 

   IM(l,l)=1.0d0

end do


!!!!!!!!!!!!! fracciones de llenado de la mezcla binaria

etha=(/0.1d0,0.1d0/)
		
rho=(6.0d0*etha)/(pi*(R**3))


call sdkbaxter

				
call evalua_flag
				

	print*, "flag=",flag



contains

subroutine sdkbaxter



!se define las srqrtN de las especies
Rhoa= 0.0d0
RhoaI=0.0d0
do l=1, nesp 

   Rhoa(l,l)=sqrt(rho(l))
   RhoaI(l,l)=1.d0/sqrt(rho(l))

end do



open(unit = 1,file = 'sdk_HS2.dat', status = 'unknown')
lamdaI=0.0d0
do nk = 1, nkmax
	k = k0 + (nk-1)*dk

call baxter(k, S)

	Z(:,:,nk) = S

call inversion(Z(:,:,nk),ZI(:,:,nk))

 
	sqrtnh(:,:,nk) = matmul(Z(:,:,nk) - IM, RhoaI)
	csqrtn(:,:,nk) = matmul(RhoaI, IM - ZI(:,:,nk))

	do l=1, nesp
		lamdaI(l,l,nk)=(1.d0 + (k/kc(l))**2)
	end do
	
	write(1,*) sngl(k),sngl(Z(1,1,nk)),sngl(Z(1,2,nk)),sngl(Z(2,1,nk)),sngl(Z(2,2,nk))
end do


close(unit = 1)

end subroutine sdkbaxter



subroutine evalua_flag
integer	:: nk, it
real(8)	:: k


!definicion de la gama
gama=0.0d0
do l=1, nesp 
  gama(l,l)=1.d-10
end do

tol=1.d-6
error=1.d0
it=0
large=1.d+20

call sdkbaxter

do while(error(1) > tol .and. error(2) > tol .and. gama(1,1) < large .and. gama(2,2) < large)


	it=it + 1

	suma= 0.d0
	do nk=1, nkmax

		k = k0 + (nk-1)*dk
		
		one = IM + k**2*matmul(gama,lamdaI(:,:,nk))

		call inversion(one,oneI)

		two= IM + k**2*matmul(matmul(gama,lamdaI(:,:,nk)),ZI(:,:,nk))

		call inversion(two,twoI)

		three=matmul(matmul(csqrtn(:,:,nk),twoI),sqrtnh(:,:,nk))
		

		do l=1, nesp
			suma(l)= suma(l) + (k**4)*oneI(l,l)*three(l,l)
		end do
		
		!write(2,*) sngl(k), sngl((k**4)*oneI(1,1)*three(1,1)), sngl((k**4)*oneI(2,2)*three(2,2))
	end do

suma=fact*suma*dk

	do l=1, nesp
	
		error(l)=abs((1.d0/suma(l) - gama(l,l))/gama(l,l))
	
	end do



	do l=1, nesp

		gama(l,l)= 1.d0/suma(l)

	end do



	do l=1, nesp

		gamalist(l,it)= gama(l,l)

	end do



end do
print*,"it=",it ,"gama=",sngl(gama(1,1)), sngl(gama(2,2))

	do l=1, nesp

		flag(l)= sign(1.d0, gamalist(l,it) - 2.0d0*gamalist(l,it - 1) + gamalist(l,it - 2))
	
	end do


end subroutine evalua_flag





subroutine inversion(a,b)
	Real(8)			:: det
	Real(8), dimension(nesp,nesp), intent(in)	:: a
	Real(8), dimension(nesp,nesp), intent(out)	:: b

	det = a(1,1)*a(2,2)-a(1,2)*a(2,1)

	b(1,1) = a(2,2)/det
	b(1,2) = -a(1,2)/det
	b(2,1) = -a(2,1)/det
	b(2,2) = a(1,1)/det
return


end subroutine


subroutine baxter(k, S)
complex(8), dimension(nesp,nesp) :: QM, QMT
Real(8), Dimension(nesp, nesp)	:: S, SInv
Complex(8), Parameter	:: I = Cmplx(0.d0,1.d0)
integer 		:: l,j
Real(8)			:: k, a(nesp), b(nesp)





do j=1, 3

	psi(j)= (pi/6.d0)*dot_product(rho,R**j)

end do

do  l=1, nesp

	a(l)= (1.d0 - psi(3) + 3.d0*R(l)*psi(2))/(1.d0 - psi(3))**2

	b(l)=-((3.d0/2.d0)*((R(l))**2)*psi(2))/(1.d0 - psi(3))**2

end do 




do l=1, nesp
	do j=1,nesp

	QM(l,j)= IM(l,j) - (pi*(sqrt(rho(l)*rho(j)))/k**3)*(2.d0*exp(I*k*((R(l) + R(j))/2.d0))*		&
		(k*(b(l) + a(l)*((R(l) + R(j))/2.d0)) + a(l)*I) + exp(I*k*((R(l)- R(j))/2.d0))*		&
		(a(l)*(k**2)*I*(((R(l)-R(j))/2.d0)**2 - ((R(j) + R(l))/2.d0)**2)- 2.d0*I*a(l) - 2.d0*b(l)*I*(k**2)*(R(j)) - 2.d0*k*(a(l)*	&
		((R(l)-R(j))/2.d0) + b(l))))
	
	QMT(l,j)= IM(l,j) + (pi*(sqrt(rho(l)*rho(j)))/k**3)*(2.d0*exp(I*(-k)*((R(j) + R(l))/2.d0))*		&
		((-k)*(b(j) + a(j)*((R(j) + R(l))/2.d0)) + a(j)*I) + exp(I*(-k)*((R(j)- R(l))/2.d0))*		&
		(a(j)*(k**2)*I*(((R(j)-R(l))/2.d0)**2 - ((R(l)+ R(j))/2.d0)**2) - 2.d0*I*a(j) - 2.d0*b(j)*I*(k**2)*(R(l)) - 2.d0*(-k)*(a(j)*	&
		((R(j)-R(l))/2.d0) + b(j))))
	end do
end do


SInv = matmul(QMT,QM)


call inversion(SInv,S)


return

end subroutine baxter



end program baxt

