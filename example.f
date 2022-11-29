      program example
      implicit none
      real*8 futil
      external  FCN, futil 

      
      open(27, file="input.fit")  !! Initial value of the parameters
      open(28, file="fitline.out") !! Output of the data
      open(29, file="sminuit.dat") !! Results of the fitting

      call mintio(27,28,29)

!=============================================================================================
!CALL MINTIO (IREAD,IWRITE,ISAVE)
!Input parameters:
!IREAD -> Fortran unit number for reading (default 5).
!IWRITE -> Fortran unit number for writing (default 6).
!Isave -> Fortran unit number for saving (default 7).
!=============================================================================================
!!!

      call minuit(FCN,futil)

      close(27)
      close(28)
      close(29)

      end

     

      subroutine FCN(npar, grad, chisq, xval, iflag, futil)
!==========================================================================================
! CALL FCN(NPAR,GRAD,FVAL,XVAL,IFLAG,FUTIL)
!input parameters:
!NPAR -> number of currently variable parameters.
!XVAL -> vector of (constant and variable) parameters.
!IFLAG -> Indicates what is to be calculated (see example below).
!FUTIL -> Name of utilitary routine (if needed, it must be declared EXTERNAL and provided by the user).
!Output parameters:
!FVAL -> The calculated function value. -> usually it's chi**2
!GRAD -> The (optional) vector of first derivatives).
!==========================================================================================	

      implicit none
      real*8 xval(3), futil, grad, chisq, delchisq
      real*8 x,y,dy,kx,xx
      real*8 a,b,c
      common/data/x(1000),y(1000),dy(1000)
      integer iflag, npar, ndat, k
      external futil
      common/pars/a,b,c

      ndat = 100

      if(iflag.eq.1) then
!read input data, calculate any necessary constants, etc.
	open (15,file='data.dat')
      	do k=1,ndat
      		read(15,*) x(k), y(k), dy(k)
      	enddo
      endif 

      if(iflag .eq. 2)then
	!calculate GRAD, the first derivatives of FVAL
	!(this is optional)
      endif


!! Here the parameters of out function are identified with the values, of the array xval. The initial values of this array is given in 'input.fit', and minuit modifies them during the minimization.
      a = xval(1)
      b = xval(2)
      c= xval(3)

!Always calculate the value of the function, FVAL, which is usually a chisquare or log likelihood. Optionally, calculation of FVAL may involve
!FTHEO = FUTIL(....)
!It is responsability of user to pass any parameter values needed by FUTIL, either through arguments, or in a COMMON block

      chisq=0.0d0
      do k=1, ndat
       delchisq = (futil(x(k)) - y(k))**2/dy(k)**2
       chisq = chisq + delchisq
      end do

      if(iflag .eq. 3)then
!!will come here only after the fit is finished. Perform any final calculations, output fitted data, etc.

      write(28, *) "Final reasult of the fitting_ "
      write(28, *) "a = ",a
      write(28, *) "b = ",b
      write(28, *) "c = ",c
      write(28, *) "chisq = ",chisq
      write(28, *) "_____________________________ "

      write(*, *) "a = ",a
      write(*, *) "b = ",b
      write(*, *) "c = ",c
      write(*, *) "chisq = ",chisq

      open(unit=8,file='results.out')
      !!!!!! one can compute here the function with the 'freshly' fitted parameters.
      do k=0,1000
      xx = 20.d0 * dfloat(k)/1000.d0
      write(8,*) xx,futil(xx)
      enddo
      
      endif
      

      return
      end

      function futil(xx)
      implicit none
      real*8 futil,xx
      real*8 a,b,c
      common/pars/a,b,c
      
      futil = a*xx**b*EXP(-1.0*c*xx)
      
      return
      end






























