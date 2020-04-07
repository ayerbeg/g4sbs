C     program epc
C     real*8 ebeam
C     integer n1
C     real*8 z1
C     character*1 fermi1
C     character*1 MPI
C     character*3 PART
C     character*1 scal1
C     real*8 p
C     real*8 thp
C     real*8 rad_l
C     real*8 result
C     real*8 epc_func
C     
C     ebeam = 4000 ! in MeV
C     z1 = 6 ! AtomiC number
C     n1 = 6 ! number of neutrons


C     PART='P' ! particle
C     SCAL1 ='Y' 
C     P = 374
C     THP = 59.8
C     ral_l = 25.0
C     result = epc_func(EBEAM,Z1,N1,FERMI1,MPI,PART,SCAL1,P,THP,rad_l)
C     print *,'epc=',result
C     end
	


C       electron beam, atomiC number, number of neutrons,
C	particle (p, pi+-,pi0), momentum, angle
	
	real*8 function epc_func(ebeam, z1, n1, partID, p, thp)

C       electron beam, atomiC number, number of neutrons, fermi motion (y/n),
C 	multi pion (y/n), particle (p, pi+-), virtual or real photons (y/n),
C	momentum, angle
	
C	real*8 function epc_func(e1, z, n, fermi, mpi, scal, p,
C      $ thp)



C       vax version
C       electroproduction yields of nucleons and pions 
C       written by j.s. oconnell and j.w. lightbody, jr. 
C       national bureau of standards 
C       april 1988 
C       transverse scaling region added

C 	modified by oara to plot like older epc's
C       modified slightly by glen warren to compile under linux (g77) sept. 02

	implicit real*8 (a-h,o-z) 
	real*8 mp,mpi0,mn

	real*8 ebeam
	integer n1
        integer z1
	integer partID
        real*8 p
        real*8 thp


C	gaw	character part*3,topa*1,ans*1,mpi*1 !,scal*1,fermi*1
	character*1 topa,ans,mpi,scal,fermi,mort,repe,units
	character*3 part
	logical*1 mom

	common/sg/ia
	common/qd/qdf 
      	common/del/ip 
 	common /fer/ fermi,scal
      	common/m/z,n
	common/kf/e1
	common/pions/part
	
	dimension c(0:4,8),x_plt(125),s_qd(125),s_qf(125),s_del(125) 
	dimension a(0:4),s_tot(125),s_dsc(125)
	data am/938.28/,amp/139.6/,mpi0/135.9/,mn/939.56/

	pi=acos(-1.0d0)

	print *,"beam:", e1,"z:",z1, " n:", n1, "partID:", partID
        print *, "p:", p, "angle: ", thp

	
C	e1= 4000
C	z1=6
C	n1=6
C	partID=1
C	p=200
C	angle=45
	


	
	fermi = 'n'
	mpi   = 'n' 
	scal  = 'n' 

	
C 	print * 
C 	print 105 
C     105 	format(' enter electron energy [mev]  ')
C 	read *,e1 
C 	print *
C 	print 103 
C     103 	format(' enter z,n   ') 
C 	read *,z,n


	z=z1
	n=n1
	e1=ebeam
	
	ia=z+n
	
	if(ia.eq.2.or.ia.eq.3.or.ia.eq.4.or.ia.eq.12.or.ia.eq.16)then
	   print *
	   print 106,ia 
 106	   format(' (using a =',i3,' spectral function)')
	elseif(ia.eq.1)then 
	   print *, ' nucleon case' 
	else
	   print * 
	   print 107
 107	   format(' no specific spectral function for this nucleus')
	   if(ia.gt.12) print 108
 108	   format(' a = 16 spectral function  will be used') 
	endif 
	
c	loop over particles
	
 2	continue

C 	print *
C 	print *,'target nucleons at rest (return) or fermi motion ("y")?'
C 	read(*,101) fermi
C     
C 	print *
C 	print *,'virtual (return) or real photons ("y")?'
C 	read(*,101) scal
C     !	mpi = 'y'
C     !	if(scal.eq.'y') then
C 	print *
C 	print *,'include multi-pion production ("y") or not (return)?'
C 	read(*,101) mpi
C     !	end if
C     
C 	print * 
C 	part = '   '
C 	print 100 
C     100	format('enter particle type [p,n,pi+,pi-,pi0] '//
C     1         '(return to exit) > ')
C 	read(5,201) part
C     101 	format(a1)
C     201	format(a3)
C     
C       'an' is effective number of nucleons for one pion production 

	if(partID.eq.1)then
	   part = 'p'
	elseif(partID.eq.-1)then
	   part = 'n'
	elseif(partID.eq.2)then
	   part = 'pi+'
	elseif(partID.eq.-2)then
	   part = 'pi-'
	elseif(partID.eq.0)then
	   part = 'pi0'
	endif
	
	
	if(part.eq.'p')then 
	   an=n/3.+2.*z/3.
	   mp = am
	   ip=1 
	elseif(part.eq.'n')then 
	   an=z/3.+2.*n/3.
	   mp = mn
	   ip=-1
	elseif(part.eq.'pi+')then 
	   an=z/3.
	   mp = amp
	   ip=2 
	elseif(part.eq.'pi-')then 
	   an=n/3.
	   mp = amp
	   ip=2 
	elseif(part.eq.'pi0')then 
	   an=2.*(n+z)/3. 
	   mp = mpi0
	   ip=0 
	else
	   call exit
	endif 
	
	if(abs(ip).eq.1)then
	   if(ia.eq.1)then 
	      dlf = 0
	   else
	      if(ia.gt.1.and.ia.lt.5)then
		 dlf=ia
	      else 
		 dlf=7.
	      endif
	      
	      al=dlf 
	      qdf=al*n*z/float(ia) 
C     print * 
C     print 102
C     102	      format(' enter levinger factor or 0 for default value')
C     read *,al
	      if(al.eq.0.) al=dlf 
	      qdf=al*n*z/float(ia) 
	   end if
	endif 
	
c	select one p/theta or multiple values (only plot output)
C     
C 	if(e1.le.0.) stop 
C 	print *
C 	print *,'write a topdrawer output file? (y/n) '
C 	read(5,101)topa
C 	if((topa.eq.'y').or.(topa.eq.'y')) go to 800

c	single momentum,angle
C 
C  1	continue
C 	print *
C 	print *,' momentum ("m") or kinetiC energy (return)?'
C 	read(*,101) mort
C 	
C 	if(mort.eq.'m') then
C 	   
C  11	   continue
C 	   print *
C 	   print 104 
C  104	   format(' enter particle momentum [MeV/c], angle [deg] ') 
C 	   read *,p,thp
C 	   if(p.le.0.)go to 11
C 	   
C 	else
C 	   
C  12	   continue
C 	   print *
C 	   print 1040 
C  1040	   format(' enter particle kinetiC energy [MeV], angle [deg] ') 
C 	   read *,tke,thp
C 	   p = sqrt(tke**2 + 2*mp*tke)
C 	   print *,' particle momentum: ',p,' [MeV/c]'
C 	   if(tke.le.0.)go to 12
C 	   
C 	end if
C

	
  	print * 
 	if(scal.eq.'y') then
 	   write(6,1090) part 
 1090	   format(' (e,',a3,') cross sections in ub.c/(q.(GeV/c)^2.sr)') 
 	else
 	   write(6,109) part 
 109	   format(' (e,',a3,') cross sections in ub/((MeV/c).sr)') 
 	end if
	


	th=thp*pi/180.

	print * 
	if(abs(ip).eq.1)then
	   e=sqrt(p**2+am**2) 
	   tp=p**2/(e+am) 
	   aj=p/e		! converts cross section from 1/MeV to 1/MeV/c
	
	   if(ia.eq.1)then
	      d2qd=0. 
	      d2qf=0. 
	   elseif(ia.gt.1)then
	      call dep(e1,tp,th,ip,d2qd)
	      if(scal.ne.'y') then
		 d2qd=d2qd*aj
	      end if
	      print 110,d2qd
 110	      format('      quasi-deuteron =  ',e10.3)
	      if(scal.ne.'y') then
		 call ep(e1,tp,th,d2qf)

		 d2qf=d2qf*aj
		 print 112,d2qf
 112		 format('          quasi-free =  ',e10.3)
	      end if
	   endif
	elseif(abs(ip).eq.2.or.ip.eq.0)then 
	   e=sqrt(p**2+amp**2)
	   tp=p**2/(e+amp)
	   aj=p/e 
	   d2qd=0.
	   d2qf=0.
	else
	   call exit
	endif 


C CHECK THIS LINE!!! (CA)

	if(p.lt.500.)then 
!	if(p.lt.200.)then 
	   d2sc=0.
	   call delta(e1,tp,th,d2del) 
	   if(scal.eq.'y') then
	      d2del=an*d2del
	   else
	      d2del=an*d2del*aj
	   end if
	   print 111,d2del
 111	   format('               delta =  ',e10.3) 
	   print *
	else
	   d2del=0. 
	   if(mpi.eq.'y') then
	      if(abs(ip).eq.2.or.ip.eq.0)then
		 call s2pi(2,e1,tp,th,d2sc1) 
		 call s2pi(-2,e1,tp,th,d2sc2) 
		 if(part.eq.'pi+')then 
		    d2sc=z*d2sc1+n*d2sc2 
		 elseif(part.eq.'pi-')then 
		    d2sc=n*d2sc1+z*d2sc2 
		 elseif(part.eq.'pi0')then 
		    d2sc=ia*d2sc1
		 else
		    call exit
		 endif 
	      elseif(abs(ip).eq.1)then 
		 call s2pi(1,e1,tp,th,d2sc1) 
		 d2sc=ia*d2sc1 
	      endif
	   end if
	endif 


	print 117,d2sc
 117	format('                d2sC =  ',e10.3)
	total=d2qd+d2qf+d2del+d2sc
	print * 
	print 113,total 
 113	format(' total cross section (ub/MeV sr) =  ',e13.4e3)
!	print *,total
	print * 
	print 115,tp,p
 115	format(' kinetiC energy = ',f10.2,' MeV.  momentum = ',f10.2,
	1    ' MeV/c') 
!	1	' mev/c') 
	if(scal.eq.'y') then
	   totale=total*p*1.0d-6
	   print *, "HERE"
	   print 1160,totale
 1160	   format(' total in  ub/(q.MeV.sr) = ',g16.3)
	else
	   totale=total/aj 
	   print 116,totale,aj
 116	   format(' total in  ub/(MeV/c.sr) = ',g16.3,' p/e = ',g16.3)
	end if
        print * 
	print *

	
C 	print *,'repeat same a, e_e, particle ("y") or not (return)?'
C 	read(*,101) repe
C 	if(repe.eq.'y') then
C 	   go to 1
C 	else
C 	   go to 2
C 	end if


C	print *, scal, aj
	epc_func = totale

	print *,"coming back (1)"

	
C c	plot output
C 
C  800	continue
C 	print * 
C 	print 1000
C  1000	format('step in momentum or angle? (m/a)> ') 
C 	read(5,101)ans 
C 
C 	if((ans.eq.'m').or.(ans.eq.'m'))then
C 
C 	   mom = .true. 
C 
C 	   print * 
C 	   print 144 
C  144	   format('enter p_lo,p_hi,p_step [MeV/c], angle [deg] > ') 
C 	   read *,p_lo,p_hi,p_step,thp
C 
C 	   th=thp*pi/180. 
C 	   n_step = nint((p_hi-p_lo)/p_step) + 1 
C 
C 	   print 1144 
C  1144	   format('limits p_min, p_max for integral sigma_total > ')
C 	   read *,p_min, p_max
C 	else 
C 
C 	   mom = .false. 
C 
C 	   print *
C 	   print *,' momentum ("m") or kinetiC energy (return)?'
C 	   read(*,101) mort
C 
C 	   if(mort.eq.'m') then
C 
C  1110	      continue
C 	      print * 
C 	      print 1990
C  1990	      format('enter momentum (MeV/c),ang_lo,ang_hi,ang_step> ') 
C 	      read *,p,thp_lo,thp_hi,thp_step 
C 	      n_step = nint((thp_hi-thp_lo)/thp_step) + 1 
C 	      if(p.le.0.)go to 1110
C 
C 	   else
C 
C  1210	      continue
C 
C 	      print *
C 	      print 1041
C  1041	      format(' enter kinetiC energy [MeV], ang_lo,ang_hi,ang_step> ') 
C 	      read *,tke,thp_lo,thp_hi,thp_step 
C 	      n_step = nint((thp_hi-thp_lo)/thp_step) + 1 
C 	      p = sqrt(tke**2 + 2*mp*tke)
C 	      if(tke.le.0.)go to 1210
C 	      print *,' particle momentum: ',p,' [mev/c]'
C 
C 	   end if
C 
C 	endif 
C 
C !	open(unit=1,name='epc.top',status='unknown',form='formatted') 
C 	open(unit=1,file='epc.top',status='unknown',form='formatted') 
C 	write(1,10) '  set size 13 10' 
C 	write(1,10) ' (set font duplex' 
C ! 	write(1,10) '  set intensity 4' 
C 	write(1,10) '  title top ''(e,p) cross section''' 
C 	ucf = 1.0d0
C 
C 	if(scal.eq.'y') then
C 
C 	   print *
C 	   print *,' energy [1/MeV] ("y") or momentum [c/(GeV/c)^2] units?',
C 	1	'(return)?'
C 	   read(*,101) units
C 
C 	   if(units.eq.'y') then
C 	      write(1,10)'  title left ''d223S/dW0p1dp [Mb/(q.MeV.sr)]'
C 	      write(1,10)'  case       '' x xg  fx x    g''' 
C 
C 	      if(ans.eq.'m'.or.ans.eq.'M') then
C 		 write(1,10) '  title bottom ''t [MeV]''' 
C 	      else
C 		 write(1,10) '  title bottom ''Q2o3'''
C 		 write(1,10) '  case       ''gx x''' 
C 	      end if
C 
C 	   else
C 
C 	      write(1,10)'  title left ''d223S/dW0p1dp [Mb.c/(q.(GeV/c)223.sr)]'''
C !	1)]''' 
C !       1''' 
C 	      write(1,10)'  case       '' x xg  fx x    g m    m       x xm''' 
C 
C 	      if(ans.eq.'m'.or.ans.eq.'M') then
C !       write(1,10) '  title bottom ''p [GeV/c]''' 
C 		 write(1,10) '  title bottom ''p [MeV/c]''' 
C 	      else
C 		 write(1,10) '  title bottom ''Q2o3'''
C 		 write(1,10) '  case       ''gx x''' 
C 	      end if
C 	   end if
C 
C 	else
C 
C 	   write(1,10) '  title left ''d223S/dW0p1dp [nb/((MeV/c).sr)]''' 
C 	   write(1,10) '  case       '' x xg  fx x               m''' 
C 
C 	   if(ans.eq.'m'.or.ans.eq.'M') then
C 	      write(1,10) '  title bottom ''p [MeV/c]''' 
C 	   else
C !       write(1,10) '  title bottom ''Q d'''
C !	write(1,10) '  case       ''g m''' 
C 	      write(1,10) '  title bottom ''Q2o3'''
C 	      write(1,10) '  case       ''gx x''' 
C 	   end if
C 
C 	end if
C 
C 	write(1,'(''  title 8.5 9.1 size 1.5 '''' e001 ='',
C 	1f9.2,'' MeV  '',a3,'' '''' '')') e1,part
C !       #                  f9.2,'' mev  '',a3,'' ' '')') e1,part
C 	write(1,'(''  case ''''  x x '''' '')')
C 	if(ans.eq.'a') then
C 	   write(1,'(''  title 8.5 8.8  size 1.5 '''' a ='',i4,''; p ='',
C 	1f9.2,'' MeV/c'''' '')') ia,p
C 	else
C 	   write(1,'(''  title 8.5 8.8 size 1.5 '''' a ='',i4,''; theta ='',
C 	1f9.3,'' deg'''' '')') ia,thp
C 	end if
C 	write(1,10) '  title 8.5 8.4 size 1.5 ''-2-2-2-2-2-2-2-2-2- total''' 
C  	write(1,10) '  case                   '' u u u u u u u u u''' 
C 	write(1,10) '  title 8.5 8.1 size 1.5 ''- - - - qd''' 
C   	write(1,10) '  title 8.5 7.8  size 1.5 ''. . . . . . . Ps''' 
C 	write(1,10) '  case                    ''m m m m m m m g''' 
C 	write(1,10) '  title 8.5 7.5  size 1.5 ''-.-.-.-.- qf''' 
C 	write(1,10) '  case                    '' m m m m    '''
C 	write(1,10) '  title 8.5 7.3  size 1.5 ''-..-..-..- D''' 
C 	write(1,10) '  case                    '' mm mm mm  f'''
C 	write(1,10) '  set limits ymin 0.' 
C 	ind = 0 


C	loop over p or th
C   	sumsig =0.		! initialize
C   
C   
C   	print *,n_step	
C   	do j=1,n_step
C   	   ind = ind + 1 
C   	   if(mom)then 
C   	      p = p_lo + (j-1.)*p_step 
C   	      x_plt(j) = p
C   	   else 
C   	      thp = thp_lo + (j-1.)*thp_step 
C   	      th=thp*pi/180. 
C   	      x_plt(j) = thp
C   	   endif 
C   	   deltax = x_plt(2) - x_plt(1) ! for integral of s_tot
C   	   if(abs(ip).eq.1)then
C   	      e=sqrt(p**2+am**2) 
C   	      tp=p**2/(e+am) 
C   	      aj=p/e 
C   	      if(ia.eq.1)then
C   		 d2qd=0. 
C   		 d2qf=0. 
C   	      elseif(ia.gt.1)then
C   		 call dep(e1,tp,th,ip,d2qd)
C   		 if(scal.ne.'y') then
C   		    d2qd=d2qd*aj
C   		 end if
C   		 if(scal.ne.'y') then
C   		    call ep(e1,tp,th,d2qf)
C   		    d2qf=d2qf*aj
C   		 end if
C   	      endif
C   	   elseif(abs(ip).eq.2.or.ip.eq.0)then 
C   	      e=sqrt(p**2+amp**2)
C   	      tp=p**2/(e+amp)
C   	      aj=p/e 
C   	      d2qd=0.
C   	      d2qf=0.
C   	   else
C   	      call exit
C   	   endif 
C   !       if(p.lt.500.)then 
C   	   if(p.lt.200.)then 
C   	      d2sc=0.
C   	      call delta(e1,tp,th,d2del) 
C   	      if(scal.eq.'y') then
C   		 d2del=an*d2del
C   	      else
C   		 d2del=an*d2del*aj
C   	      end if
C   	   else
C   	      d2del=0. 
C   	      if(mpi.eq.'y') then
C   		 if(abs(ip).eq.2.or.ip.eq.0)then
C   		    call s2pi(2,e1,tp,th,d2sc1) 
C   		    call s2pi(-2,e1,tp,th,d2sc2)
C   		    if(part.eq.'pi+')then 
C   		       d2sc=z*d2sc1+n*d2sc2 
C   		    elseif(part.eq.'pi-')then 
C   		       d2sc=n*d2sc1+z*d2sc2 
C   		    elseif(part.eq.'pi0')then 
C   		       d2sc=ia*d2sc1
C   		    else
C   		       call exit
C   		    endif 
C   		 elseif(abs(ip).eq.1)then 
C   		    call s2pi(1,e1,tp,th,d2sc1) 
C   		    d2sc=ia*d2sc1 
C   		 endif
C   	      end if
C   	   endif 
C   	   
C   	   if(scal.eq.'y'.and.units.eq.'y') then
C   	      ucf = p*1.0d-9
C   	   else if(scal.eq.'y') then
C   	      ucf = 1.0d-3
C   	   else
C   	      ucf = 1.0d0
C   	   end if
C   
C   	   s_qd(j) = d2qd*1000.*ucf
C   	   s_qf(j) = d2qf*1000.*ucf
C   	   s_del(j) = d2del*1000.*ucf
C   	   s_dsc(j) = d2sc*1000.*ucf
C   	   s_tot(j) = s_qf(j) + s_qd(j) + s_del(j) + s_dsc(j)
C   	   if(ans.eq.'m'.or.ans.eq.'m') then
C   	      if(x_plt(j).ge.p_min.and.x_plt(j).le.p_max) then
C   !       sumsig = sumsig + s_tot(j)*deltax
C   		 sumsig = sumsig + s_tot(j)*p_step
C   	      end if
C   	   end if
C   	end do			! end loop over p or th

	
c	if(ans.eq.'m'.or.ans.eq.'m') then
c	   write(1,'(''   title 9.7 8.4  size 1.5 '''' IS ='',
c	1	g13.6,'' nb/sr  '')') sumsig
c	   write(1,'(''   case '''' mg '''' '')')
c	end if
c
c
c	
C       last point

C	x_plt(j) = x_plt(j-1)
C	x_plt(j+1) = x_plt(j-1) 
C	s_qd(j) = s_qd(j-1) 
C	s_qf(j) = s_qf(j-1)
C	s_del(j) = s_del(j-1)
C	s_dsc(j) = s_dsc(j-1)
C	s_tot(j) = s_qf(j) + s_qd(j) + s_del(j) + s_dsc(j)
C	s_qd(j+1) = s_qd(j-1) 
C	s_qf(j+1) = s_qf(j-1)
C	s_del(j+1) = s_del(j-1)
C	s_dsc(j+1) = s_dsc(j-1)
C	s_tot(j+1) = s_qf(j) + s_qd(j) + s_del(j) + s_dsc(j)

C       write the output topdrawer file 

C 	iflag = 0
C 	do i=1,ind-1,3
C 	   if(s_tot(i).ne.0.0.or.s_tot(i+1).ne.0.0.or.s_tot(i+2).ne.0) then
C 	      write(1,1010)x_plt(i),s_tot(i),x_plt(i+1),s_tot(i+1),
C 	1	   x_plt(i+2),s_tot(i+2) 
C 	      iflag = 1
C 	   end if
C 	enddo 
C 	if(iflag.ne.0) write(1,10) '  join 1 solid' 
C 
C 	iflag = 0
C 	do i=1,ind-1,3
C 	   if(s_qd(i).ne.0.0.or.s_qd(i+1).ne.0.0.or.s_qd(i+2).ne.0) then
C 	      write(1,1010)x_plt(i),s_qd(i),x_plt(i+1),s_qd(i+1),
C 	1	   x_plt(i+2),s_qd(i+2) 
C 	      iflag = 1
C 	   end if
C 	enddo 
C 	if(iflag.ne.0) write(1,10) ' join 1 dash' 
C 
C 	iflag = 0
C 	do i=1,ind-1,3
C 	   if(s_qf(i).ne.0.0.or.s_qf(i+1).ne.0.0.or.s_qf(i+2).ne.0) then
C 	      write(1,1010)x_plt(i),s_qf(i),x_plt(i+1),s_qf(i+1),
C 	1	   x_plt(i+2),s_qf(i+2) 
C 	      iflag = 1
C 	   end if
C 	enddo 
C 	if(iflag.ne.0) write(1,10) ' join 1 dotdash' 
C 
C 	iflag = 0
C 	do i=1,ind-1,3
C 	   if(s_dsc(i).ne.0.0.or.s_dsc(i+1).ne.0.0.or.s_dsc(i+2).ne.0) then
C 	      write(1,1010)x_plt(i),s_dsc(i),x_plt(i+1),s_dsc(i+1),
C 	1	   x_plt(i+2),s_dsc(i+2) 
C 	      iflag = 1
C 	   end if
C 	enddo 
C 	if(iflag.ne.0) write(1,10) '  join 1 dot' 
C 
C 	write(1,10) '  set pattern .01 .09 .01 .09 .01 .09 .1'
C 	iflag = 0
C 	do i=1,ind-1,3
C 	   if(s_del(i).ne.0.0.or.s_del(i+1).ne.0.0.or.s_del(i+2).ne.0) then
C 	      write(1,1010)x_plt(i),s_del(i),x_plt(i+1),s_del(i+1),
C 	1	   x_plt(i+2),s_del(i+2) 
C 	      iflag = 1
C 	   end if
C 	enddo 
C 	if(iflag.ne.0) write(1,10) '  join 1 pattern' 
C                  	close(unit=1)

 10	format(a) 
 1010	format(2x,3(2g12.5,'; '))


	print *,"coming back (END)"

	
C	go to 2 
	return ! because is a function (CA)
	end


	
*       vtp
	subroutine vtp(amt,am1,ei,w0,tp,th,gn)
C       tiator-wright virtual photon spectrum
C       phys. rev. c26,2349(1982) and nuc. phys. a379,407(1982)
	implicit real*8 (a-h,o-z) 
	data ame/.511/
	pi=acos(-1.0d0)
	ef0=ei-w0 
	aki=sqrt(ei**2-ame**2)
	akf0=sqrt(ef0**2-ame**2)
	akp=sqrt(tp**2+2.*am1*tp) 
	ep=tp+am1 
	ar=ei+amt-ep
	br=ef0*(akp*cos(th)-aki)/akf0 
C       brp=(akf0/ef0)**2*br
	a=ame**2-ei*ef0 
	b=aki*akf0
	d=-ame**2*br*(ei/ef0-1.)/ar 
	ap=a-d
	bp=b+d
	an1=1./137./2./pi*w0**2/aki**2
!       apb=-ame**2*(aki-akf0)**2/(ame**2+ei*ef0+aki*akf0)
!	print *,'apb',apb
	apb = a + b		! new
!	print *,'a+b',apb
	an1=an1*b/bp*(ar+br)/(ar-ap/bp*br)
	an2=1.-2.*a/w0**2 
	an4=((ap-bp)*(ar+br)/apb/(ar-br)) 
	if(an4.le.0.)go to 1
	an2=an2*log(an4)
	an3=-4.*b/w0**2 
	ane=an1*(an2+an3) 
	d0=amt+ei-ep+ef0/akf0*(akp*cos(th)-aki) 
	r=(amt+w0-ep/akp*w0*cos(th))/d0 
	gn=ane*r/w0 
	if(gn.lt.0.)gn=0. 
!	print *,'tp,w0,r,ane,gn',tp,w0,r,ane,gn
	return
 1	gn=0. 
	return
	end


	
*       dep
	subroutine dep(e1,tp,th,ip,d2qd)
C       quasi-deuteron cross section 
	implicit real*8 (a-h,o-z) 
	character*1 scal,fermi
	common/sg/ia
	common/qd/qdf 
 	common /fer/ fermi,scal
	data am/939./,amd/1876./
	if(ia.eq.1)goto 1 
	pn=sqrt(tp**2+2.*am*tp) 
	ep = tp + am
	call kine(amd,am,am,pn,th,w0,thc) 
	if(w0.ge.e1)go to 1 
	if(w0.le.0.)go to 1 
	w0g=w0/1000.
	call sigd(w0g,thc,ip,dsqd)
	call part(amd,am,am,pn,th,ajt,ajw)
	if(scal.ne.'y') then
	   call vtp(amd,am,e1,w0,tp,th,phi)
C       cross section in ub/mev-sr 
C	oara: ajw gives (w+-w-)/(p+-p-) so it is dw/dp, not dw/dt!
c	cross section is then in ub/((mev/c).sr) directly
!	print *,'qdf,phi,dsqd,ajw,ajt',qdf,phi,dsqd,ajw,ajt
	   d2qd=qdf*phi*dsqd*ajw*ajt*ep/pn
	else
	   brem=1./w0
C       cross section in ub/mev-sr-q
C	cross section was then in ub/(q.(mev/c).sr) directly
	   d2qd = qdf*brem*dsqd*ajw*ajt*ep/pn ! from gpc
!	print *,'qdf,brem,dsqd,ajw,ajt',qdf,brem,dsqd,ajw,ajt
!	d2qd = qdf*brem*dsqd*ajw*ajt	one step conversion: next line
!	d2qd = d2qd*ep*1.0d6/pn**2 ! convert to ub.c/(q.(gev/c)^2.sr)
	   d2qd = d2qd*1.0d6/pn	! convert to units ub.c/(q.(gev/c)^2.sr)
	end if
	return
 1	d2qd=0. 
	return
	end

	
*       sigd 
	subroutine sigd(e,th,ip,dsqd) 
C       deuteron cross section 
C       based on fit of thorlacius & fearing 
C       phys. rev. c33,1830(1986)
C       photon energy range 10 - 625 mev 
C       e[gev] in lab system 
C       th[rad] & dsqd[ub/sr] in center-of-momentum system 
	implicit real*8 (a-h,o-z) 
	dimension c0(8),c1(4),c2(4),c3(4),c4(4) 
	dimension a(0:4),b(4,4)
	data c0/2.61e2,-1.10e2,2.46e1,-1.71e1,5.76e0, 
	1    -2.05e0,2.67e-1,1.13e2/ 
	data c1/1.68e1,-4.66e1,2.56e0,-4.72e0/
	data c2/-2.03e2,-8.12e1,-4.05e0,-5.99e0/
	data c3/-1.77e1,-3.74e1,-5.07e-1,-5.40e0/ 
	data c4/-2.05e0,-7.05e0,9.40e-1,-2.05e0/
	x=cos(th) 
	if(e.le.625.)then 
C       test for neutron 
	   x=ip*x 
C       coeficients
	   a(0)=c0(1)*exp(c0(2)*e)+c0(3)*exp(c0(4)*e) 
	   a(0)=a(0)+(c0(5)+c0(6)*e)/(1.+c0(8)*(e-c0(7))**2)
	   dsqd=a(0)*pl(0,x)
	   do 2 l=1,4 
	      b(1,l)=c1(l) 
	      b(2,l)=c2(l) 
	      b(3,l)=c3(l) 
 2	      b(4,l)=c4(l) 
	      do 1 l=1,4 
		 a(l)=b(l,1)*exp(b(l,2)*e)+ b(l,3)*exp(b(l,4)*e)
 1		 dsqd=dsqd+a(l)*pl(l,x) 
	      elseif(e.lt..700)then 
		 dsqd=.3
	      elseif(e.lt..800)then 
		 dsqd=.15 
	      elseif(e.lt..900)then 
		 dsqd=.1
	      else
		 dsqd=55./(e-.350)
	      endif 
	      return
	      end

	
*       leg
	real*8 function pl(l,x) 
C       legendre polynomials 
	implicit real*8 (a-h,o-z) 
	if(l.eq.0)then
	   pl=1.
	elseif(l.eq.1)then
	   pl=x 
	elseif(l.eq.2)then
	   pl=.5*(3.*x**2-1.) 
	elseif(l.eq.3)then
	   pl=.5*(5.*x**3-3.*x) 
	elseif(l.eq.4)then
	   pl=1./8.*(35.*x**4-30.*x**2+3.)
	else
	   pl=0.
	endif 
	return
	end


	
*       delta
	subroutine delta(e1,tp,th,d2del)
C       photoproduction of nucleons and pions via delta
	implicit real*8 (a-h,o-z) 
	common/del/ip
	character*1 scal,fermi
 	common /fer/ fermi,scal
	data am/939./,amp/139./ 
	if(abs(ip).eq.1)then
	   am1=am 
	   am2=amp
	else
	   am1=amp
	   am2=am 
	endif 
	ep=tp+am1 
	pn=sqrt(ep**2-am1**2) 
!	print *
!	print *,'delta',pn
	if(fermi.ne.'y') then
	   call kine(am,am1,am2,pn,th,w,tc)
	   cpf = 1.
	else
	   call kindel(am,am1,am2,pn,th,w,tc,cpf)
	end if
!	print *,'pn,w,cpf',pn,w,cpf
	if(w.le.0.)go to 1
!       if(w.ge.e1)go to 1
	if(w.gt.e1)go to 1
	if(fermi.ne.'y') then
	   call part(am,am1,am2,pn,th,ajt,ajw) 
	else
	   call partdel(am,am1,am2,pn,th,ajt,ajw) 
	end if
	call sigma(w,tc,dsigg)
	if(scal.ne.'y') then
C       cross section in ub/mev-sr
	   call vtp(am,am1,e1,w,tp,th,phi) 
!	write(*,100) am,am1,e1,w,th,phi
 100	   format(2x,2f8.1,5g12.3)
	   d2del = phi*dsigg*ajt*cpf 
!	print *,'phi,ajt,dsigg,d2del,cpf',phi,ajt,dsigg,d2del,cpf
!	print *
	else
C       cross section in ub/mev-sr-q
	   brem = 1./w
	   dsig = brem*dsigg*ajt*ajw
!       d2del = dsig*cpf
!	d2del = d2del*ep*1.0d6/pn**2 ! convert to ub.c/(q.(gev/c)^2.sr)
	   d2del = dsig*ep/pn
	   d2del = d2del*1.0d6/pn ! convert to units ub.c/(q.(gev/c)^2.sr)
	   d2del = d2del*cpf
	end if
	return
 1	d2del=0.
	return
	end


	
*       part 
	subroutine part(amt,am1,am2,pn,tn,ajt,ajw)
C       partial derivatives
	implicit real*8 (a-h,o-z) 
	pi=acos(-1.0d0)
!       dt=pi/180.0d0
!       dp=10.0d0
	dt=pi/720.0d0
	dp=2.0d0
C       angle
	tnp=tn+dt 
	tnm=tn-dt 
	tn0=tn
	call kine(amt,am1,am2,pn,tnp,wp,tcp) 
	call kine(amt,am1,am2,pn,tnm,wm,tcm) 
	call kine(amt,am1,am2,pn,tn0,w,tc0) 
	den=cos(tnp)-cos(tnm) 
	den=abs(den)
!	print *,'den,tcp,tcm,tc0',den,tcp,tcm,tc0
	if(den.gt.1.0d-3.and.(w*wp*wm.gt.0.))then
	   ajt=(cos(tcp)-cos(tcm))/den
	   ajt=abs(ajt) 
	else
	   ajt=(cos(tc0)-cos(tcm))/(cos(tn0)-cos(tnm))
	   ajt=abs(ajt) 
	endif 
C       energy 
	pnp=pn+dp 
	pnm=pn-dp 
	call kine(amt,am1,am2,pnp,tn,wp,tc) 
	call kine(amt,am1,am2,pnm,tn,wm,tc) 
	am12 = am1**2
	tp = sqrt(pnp**2 + am12) - am1
	tm = sqrt(pnm**2 + am12) - am1
	ajw=(wp-wm)/(pnp-pnm) 
!	ajw=(wp-wm)/(tp-tm) 
	ajw=abs(ajw)
	return
	end


	
*       kine 
	subroutine kine(amt,am1,am2,pn,th,w,tc) 
C       computes cm variables from lab variables 
	implicit real*8 (a-h,o-z) 
	common /kf/ e1
	amt2 = amt**2
	am12 = am1**2
	ep=sqrt(pn**2+am12) 
	pnt=pn*sin(th)
	pnl=pn*cos(th)
!       anum=pn**2+am2**2-(amt-ep)**2 
	anum = am2**2 - am12 - amt2 +2*ep*amt
	den=2.*(pnl+amt-ep) 
	w=anum/den
	if(w.le.0.)w=0. 
C       invariant mass 
	ww2 = amt**2+2.*w*amt
	ww = sqrt(ww2)
!	print *,' p,w,e_thr ',pn,ww,w,anum/den
C       cm variables 
	pct=pnt 
	b=w/(amt+w) 
	g=(w+amt)/ww
	pcl=g*(pnl-b*ep)
	pcs=pcl**2+pct**2 
	pc=sqrt(pcs)
	cthc=pcl/pC 
	tc=acos(cthc) 
	return
	end


	
*       sigma
	subroutine sigma(e,thrcm,sigcm) 
C       real photon cross section in delta region
C       microbarns per steradian 
	implicit real*8 (a-h,o-z) 
	gam=100.
	pi=acos(-1.0d0)
	if(e.gt.420.)then 
	   sigcm=(1.+420./e)*90./4./pi 
	else
	   sigcm=360.*(5.-3.*cos(thrcm)**2)
	   sigcm=sigcm/16./pi/(1.+(e-320)**2/gam**2) 
	endif 
	return
	end


	
*       ep 
	subroutine ep(e1,tp,thp,dsep) 
C       electro proton production cross sections 
	implicit real*8 (a-h,o-z) 
	common/p/ph(10),wph(10) 
	data aml/.511/
	pi=acos(-1.0d0)
	call gausab(10,ph,wph,0.d0,2.*pi,pi)
	ak=sqrt(e1**2-aml**2) 
	call sep(ak,tp,thp,dsep)
	dsep=dsep*1.e4
C       cross section in ub/mev-sr 
	end


	
*       dot
	real*8 function dot(v,u)
	implicit real*8 (a-h,o-z) 
	dimension v(3),u(3) 
	dot=0.
	do 1 i=1,3
 1	   dot=dot+v(i)*u(i) 
	   return
	   end

	
*       cross
	subroutine cross(v,u,w) 
	implicit real*8 (a-h,o-z) 
	dimension v(3),u(3),w(3)
	w(1)=v(2)*u(3)-v(3)*u(2)
	w(2)=v(3)*u(1)-v(1)*u(3)
	w(3)=v(1)*u(2)-v(2)*u(1)
	return
	end

	
*       gausab 
	subroutine gausab(n,e,w,a,b,c)
	implicit real*8 (a-h,o-z) 
	dimension e(*),w(*) 
	data eps/1.d-16/
	if(a.ge.c.or.c.ge.b)stop
C       stops program if a, c, b are out of sequence
	pi=acos(-1.0d0)
	al=(c*(a+b)-2*a*b)/(b-a)
	be=(a+b-2*c)/(b-a)
	m=(n+1)/2 
	dn=n
	do 5 i=1,m
	   di=i 
	   x=pi*(4.d0*(dn-di)+3.d0)/(4.d0*dn+2.d0)
	   xn=(1.d0-(dn-1.d0)/(8.d0*dn*dn*dn))*cos(x) 
	   if(i.gt.n/2) xn=0
	   do 3 iter=1,10 
	      x=xn
	      y1=1.d0 
	      y=x 
	      if(n.lt.2) go to 2
	      do 1 j=2,n
		 dj=j 
		 y2=y1
		 y1=y 
 1		 y=((2.d0*dj-1.d0)*x*y1-(dj-1.d0)*y2)/dj
 2		 continue
		 ys=dn*(x*y-y1)/(x*x-1.d0) 
		 h=-y/ys 
		 xn=x+h
		 if(abs(h).lt.eps) go to 4 
 3	      continue
 4	      e(i)=(c+al*x)/(1.d0-be*x)
	      e(n-i+1)=(c-al*x)/(1.d0+be*x)
	      gew=2.d0/((1.d0-x*x)*ys*ys)
	      w(i)=gew*(al+be*c)/(1.d0-be*x)**2
	      w(n-i+1)=gew*(al+be*c)/(1.d0+be*x)**2
 5	   continue 
	   return
	   end


	
*       vect 
	subroutine vect(thp,the,phi,p,ak1,ak2)
C       cartesian components of electron and proton vectors
	implicit real*8 (a-h,o-z) 
	common/v/ak1v(3),ak2v(3),qv(3),pv(3),pp(3)
	pv(1)=p*sin(thp)
	pv(2)=0.
	pv(3)=p*cos(thp)
	ak1v(1)=0.
	ak1v(2)=0.
	ak1v(3)=ak1 
	ak2v(1)=ak2*sin(the)*cos(phi) 
	ak2v(2)=ak2*sin(the)*sin(phi) 
	ak2v(3)=ak2*cos(the)
	qv(1)=ak1v(1)-ak2v(1) 
	qv(2)=ak1v(2)-ak2v(2) 
	qv(3)=ak1v(3)-ak2v(3) 
	pp(1)=pv(1)-qv(1) 
	pp(2)=pv(2)-qv(2) 
	pp(3)=pv(3)-qv(3) 
	return
	end


	
*       amag 
	real*8 function amag(v) 
	implicit real*8 (a-h,o-z) 
	dimension v(3)
	amag=0. 
	do 1 i=1,3
 1	   amag=amag+v(i)**2 
	   amag=sqrt(amag) 
	   return
	   end


	
*       lept 
	subroutine lept(e1,e2,ak1,ak2,aml,qs,qus,the,v) 
C       lepton factors for coincidence cross section 
	implicit real*8 (a-h,o-z) 
	dimension v(5)
	v(1)=(qus/qs)**2*(e1*e2+ak1*ak2*cos(the)+aml**2)
	x=ak1*ak2*sin(the)
	v(2)=x**2/qs+qus/2. 
	v(3)=qus/qs*x/sqrt(qs)*(e1+e2)
	v(4)=x**2/qs
	v(5)=0. 
	return
	end


	
*       d4s
	subroutine d4s(ak1,ak2,the,p,pp,thqp,cphip,dsig)
C       fully differential cross section 
	implicit real*8 (a-h,o-z) 
	dimension v(5),w(5) 
	data am/939./,aml/.511/,a/855./ 
	qs=ak1**2+ak2**2-2.*ak1*ak2*cos(the)
	e1=sqrt(ak1**2+aml**2)
	e2=sqrt(ak2**2+aml**2)
	qus=2.*(e1*e2-ak1*ak2*cos(the)-aml**2)
	sm=2.*(1.44)**2/qus**2*ak2/ak1
	ep=sqrt(am**2+p**2) 
	ps=ep*p 
	fns=1./(1.+qus/a**2)**4 
	call lept(e1,e2,ak1,ak2,aml,qs,qus,the,v) 
	call form(qs,p,thqp,cphip,w)
	sum=0.
	do 1 i=1,5
 1	   sum=sum+v(i)*w(i) 
	   dsig=sm*ps*fns*sum*sgsl(pp) 
	   return
	   end


	
*       sthe 
	subroutine sthe(d2s)
C       integral over electron polar angle 
	implicit real*8 (a-h,o-z) 
	common/s/ ak1,ak2,the,p,thp 
	common/e/th1(12),wt1(12),th2(12),wt2(12)
	common/e1/th3(24),wt3(24) 
	d2s1=0. 
	do 1 i=1,12 
	   the=th1(i)
	   call sphi(d3s)
 1	   d2s1=d2s1+d3s*wt1(i)*sin(the) 
	   d2s2=0. 
	   do 2 i=1,12 
	      the=th2(i)
	      call sphi(d3s)
 2	      d2s2=d2s2+d3s*wt2(i)*sin(the) 
	      d2s3=0. 
	      do 3 i=1,24 
		 the=th3(i)
		 call sphi(d3s)
 3		 d2s3=d2s3+d3s*wt3(i)*sin(the) 
		 d2s=d2s1+d2s2+d2s3
		 return
		 end

	
*       sphi 
	subroutine sphi(d3s)
C       integrate over electron azimuthal angle
	implicit real*8 (a-h,o-z) 
	common/s/ ak1,ak2,the,p,thp 
	common/v/ak1v(3),ak2v(3),qv(3),pv(3),pp(3)
	common/p/ph(10),wph(10) 
	dimension qxp(3),ak1x2(3) 
	d3s=0.
	do 1 i=1,10 
	   phi=ph(i) 
	   call vect(thp,the,phi,p,ak1,ak2)
	   call cross(qv,pv,qxp) 
	   call cross(ak1v,ak2v,ak1x2) 
C       proton theta 
	   cthep=dot(pv,qv)/amag(pv)/amag(qv)
	   thqp=acos(cthep)
C       proton phi 
	   cphip=dot(qxp,ak1x2)
	   if (cphip.eq.0.) then 
	      cphip=1. 
	   else
	      cphip=cphip/amag(qxp)/amag(ak1x2)
	   endif 
	   ppm=amag(pp)
	   call d4s(ak1,ak2,the,p,ppm,thqp,cphip,dsig) 
 1	   d3s=d3s+dsig*wph(i) 
	   return
	   end

	
*       form 
	subroutine form(qs,p,thqp,cphip,w)
C       nuclear form factors 
	implicit real*8 (a-h,o-z) 
	common/m/zz,nn
	common/del/ip 
	dimension w(5)
	data am/939./,up/2.79/,un/-1.91/
	if(ip.eq.1)then 
	   z=zz 
	   n=0. 
	elseif(ip.eq.-1)then
	   z=0. 
	   n=nn 
	else
	   z=0. 
	   n=0. 
	endif 
	y=p/am*sin(thqp)
	w(1)=z
	w(2)=z*y**2 
	w(2)=w(2)+(z*up**2+n*un**2)*qs/2./am**2 
	w(3)=-2.*z*y*cphip
	w(4)=z*y**2*(2.*cphip**2-1.)
	w(5)=0. 
	return
	end


	
*       sep
	subroutine sep(ak,tp,thpp,d2s)
	implicit real*8 (a-h,o-z) 
	common/s/ak1,ak2,the,p,thp
	common/e/th1(12),wt1(12),th2(12),wt2(12)
	common/e1/th3(24),wt3(24) 
	data am/939./,aml/.511/,be/16./ 
	pi=acos(-1.0d0)
	thp=thpp
	ak1=ak
	ak2=ak1-tp-be 
C       gaussian points for the
	themax=aml*(ak1-ak2)/ak1/ak2
	call gausab(12,th1,wt1,0.d0,2.*themax,themax) 
	call gausab(12,th2,wt2,2.*themax,100.*themax,10.*themax)
	a3=100.*themax
	c3=a3+(pi-a3)/10. 
	call gausab(24,th3,wt3,a3,pi,c3)
	p=sqrt(tp**2+2.*tp*am)
        call sthe(d2s)
	if(ak2.le.0.)d2s=0. 
	return
	end 


*       sgsl 
	real*8 function sgsl(p) 
	implicit real*8 (a-h,o-z) 
C       p integral over sgsl normalized to 1/4pi 
	common/sg/ia
	common/qd/qdf 
	if(ia.eq.2)then 
C       begin 2-h
	   pp=p/197.3 
	   sgs=3.697-7.428*pp-2.257*pp**2 
	   sgs=sgs+3.618*pp**3-1.377*pp**4+.221*pp**5-.013*pp**6
	   if(sgs.lt.-293.)go to 1
	   sgs=exp(sgs) 
	   sgs=sgs/.18825/4./3.1416/(197.3)**3
	   sgsl=sgs/1.
	elseif(ia.eq.3)then 
C       begin 3-he 
	   if(-(p/33)**2.lt.-293.)go to 1 
	   sgs=2.4101e-6*exp(-p/33) 
	   sgs=sgs-1.4461e-6*exp(-(p/33)**2)
	   sgs=sgs+1.6871e-10*exp(-(p/493)**2)
	   sgsl=sgs/2.0d0
	elseif(ia.eq.4)then 
C       begin 4-he
	   if(-(p/113.24)**2.lt.-293.)go to 1 
	   sgs=1.39066e-6*exp(-(p/113.24)**2) 
	   sgs=sgs+3.96476e-9*exp(-(p/390.75)**2) 
	   sgsl=sgs/2.0d0
	   sgsl=sgsl/2.0d0/3.1416
	elseif(ia.gt.4.and.ia.lt.12)then
	   if(-(p/127)**2.lt.-293.)go to 1
	   sgs=1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2) 
	   sgs=sgs+1.7052e-9*exp(-(p/493)**2) 
	   sgsl=sgs/(float(ia)/2.0d0)
	elseif(ia.eq.12)then
C       begin 12-C 
	   if(-(p/127)**2.lt.-293.)go to 1
	   sgs=1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2) 
	   sgs=sgs+1.7052e-9*exp(-(p/493)**2) 
	   sgsl=sgs/6.
	else
C       begin 16-o 
	   if(-(p/120)**2.lt.-293.)go to 1
	   sgs=3.0124e-7*(1.+(p/120)**2)*exp(-(p/120)**2) 
	   sgs=sgs+1.1296e-9*exp(-(p/493)**2) 
	   sgsl=sgs/(float(ia)/2.0d0)
	endif 
	return
 1	sgsl=0. 
	return
	end


	
*       s2pi 
	subroutine s2pi(ip,e1,tp,th,d2sc) 
C       integral over scaling cross section
	implicit real*8 (a-h,o-z) 
	character*1 scal,fermi
	character*3 part
	
 	common /fer/ fermi,scal
	common /pions/ part	
	data am/939./,amp/139./ 
	
	if(part.eq.'pi0')then
	   afit = 11150.d0/12.d0 ! pi0 scaling fit a*exp(c*pt) a/12 = per nucleon
	   cfit = -8.67d0
	elseif(part.eq.'pi+')then
	   afit = 12146.d0/12.d0 ! pi+ scaling fit a*exp(c*pt) a/12 = per nucleon
	   cfit = -8.59d0
	elseif(part.eq.'pi-')then
	   afit = 10108.d0/12.d0 ! pi- scaling fit a*exp(c*pt) a/12 = per nucleon
	   cfit = -8.78d0
	endif             
	
	if(abs(ip).eq.1)then
C       one pion thr 
	   ap=am
	   am2=amp
	elseif(ip.eq.2.or.ip.eq.0)then
C       one pion thr 
	   ap=amp 
	   am2=am 
	elseif(ip.eq.-2)then
C       two pion thr 
	   ap=amp 
	   am2=am+amp 
	else
	   stop 
	endif 
	p=sqrt(tp**2+2.*ap*tp)
	e=tp+ap 
!	print *,'s2pi',p
	if(fermi.ne.'y') then
	   call kine(am,ap,am2,p,th,thr,tc)
	   cpf = 1.
	else
	   call kindel(am,ap,am2,p,th,thr,tc,cpf)
	end if
!       print *,'p,thr',p,thr
	if(thr.le.0.)goto 2 
!       if(e1.le.thr)goto 2 
	if(e1.lt.thr)goto 2 
	dw=(e1-thr)/20.0d0
	sum=0.0d0
	sumbr=0.0d0

!       do 1 i=1,20 
	do i=1,20 
	   w=thr+(float(i)-.5d0)*dw
	   if(w.lt.(e1-0.511)) then
	      call vtp(am,amp,e1,w,tp,th,gn)
!	write(*,100) i,am,amp,e1,w,dw,thr,gn
!       100	format(i4,2f8.1,5g12.3)

!	call wiser(w/1.e3,p/1.e3,th,f)	!test of scaling fit
	      pt = p/1.e3*sin(th) ! pion transverse pt
	      f = afit*exp(cfit*pt) ! from scaling fit to Yerevan 12C gamma data

!	print *,'p,th,pt,f,e1',p,th,pt,f,e1,afit,cfit
!	print *,'p,w,gn,f',p,w,gn,f

	      sum=sum+gn*f*dw 
	      sumbr = sumbr + f*dw/w

!	write(20,111)'p,th,sum,sumbr',i,p/1.e3,th,sum,sumbr
!	write( *,111)'p,th,sum,sumbr',i,p/1.e3,th,sum,sumbr
!111	format(a14,5g10.3)

!       1 continue
	   end if
	end do

	if(scal.eq.'y') then
	   d2sc=sumbr
	else
	   d2sc=sum*p**2/e*1.e-6
	end if
	d2sc = d2sc*cpf
!	print *,'d2pi',d2sc,cpf
	return
 2	d2sc=0. 
	return
	end


	
*       wiser
	subroutine wiser(w,p,th,f)
C       invariant inclusive cross section
C       units in gev 
	implicit real*8 (a-h,o-z) 
	common/del/ip 
	dimension a(7),b(7),c(7)
	data a/5.66e2,8.29e2,1.79,2.10,-5.49,-1.73,0./
	data b/4.86e2,1.15e2,1.77,2.18,-5.23,-1.82,0./
	data c/1.33e5,5.69e4,1.41,0.72,-6.77,1.90,-1.17e-2/ 
	data am/.939/,amp/.139/ 
	if(abs(ip).eq.1)then
	   ap=am
	elseif(abs(ip).eq.2.or.ip.eq.0)then 
	   ap=amp 
	else
	   stop 
	endif 
	e=sqrt(p**2+ap**2)
C       mandelstam variables 
	s=2.*w*am+am**2 
C       t=-2.*w*e+2.*w*p*cos(th) +ap**2 
	u=am**2-2.*am*e+ap**2 
C       fitting variables
	pt=p*sin(th)
	aml=sqrt(pt**2+ap**2) 
	call fxr(w,p,e,th,xr) 
C       fitted 
	if(ip.eq.2.or.ip.eq.0)then
	   x1=a(1)+a(2)/sqrt(s) 
	   x2=(1.-xr+a(3)**2/s)**a(4) 
	   x3=exp(a(5)*aml) 
	   x4=exp(a(6)*pt**2/e) 
	   f=x1*x2*x3*x4
	elseif(ip.eq.-2)then
	   x1=b(1)+b(2)/sqrt(s) 
	   x2=(1.-xr+b(3)**2/s)**b(4) 
	   x3=exp(b(5)*aml) 
	   x4=exp(b(6)*pt**2/e) 
	   f=x1*x2*x3*x4
	elseif(abs(ip).eq.1)then
	   x1=c(1)+c(2)/sqrt(s) 
	   x2=(1.-xr+c(3)**2/s)**c(4) 
	   x3=exp(c(5)*aml) 
	   x4=1./(1.+abs(u))**(c(6)+c(7)*s) 
	   f=x1*x2*x3*x4
	else
	   stop 
	endif 
	return
	end 
*       fxr
	subroutine fxr(w,p,e,th,xr) 
C       computes ratio of cm particle momentum to photon momentum
C       gev units
	implicit real*8 (a-h,o-z) 
	data am/.939/ 
	pt=p*sin(th)
	pl=p*cos(th)
C       lorentz transformation 
	b=w/(w+am)
	d=sqrt(2.*w*am+am**2) 
	g=(w+am)/d
	bg=b*g
C       cm variables
	wc=g*w-bg*w 
	plc=g*pl-bg*e 
	pc=sqrt(pt**2+plc**2) 
	xr=pc/wc
	return
	end

	
*       kindel
	subroutine kindel(amt,am1,am2,pn,th,w,tc,cpf) 
C       computes cm variables from lab variables 
	implicit real*8 (a-h,o-z) 
	common /kf/ e1
	data cosav/0.63662/
	amt2 = amt**2
	am12 = am1**2
	ep=sqrt(pn**2+am12) 
	pnt=pn*sin(th)
	pnl=pn*cos(th)
!       anum=pn**2+am2**2-(amt-ep)**2 
	anum = am2**2 - am12 - amt2 +2*ep*amt
	den=2.*(pnl+amt-ep) 
	w=anum/den
	if(w.le.0.)w=0. 
C       invariant mass 
	ww2 = amt**2+2.*w*amt
	ww = sqrt(ww2)
!	print *,' p,w,e_thr ',pn,ww,w,anum/den

c	section on fermi motion

	if(cpf.ge.0.) then	! cpf = -1. for jacobian computation
	   cpf = 1.0d0
	   if(w.ge.e1) then
	      dm2 = ww2 - amt2 
	      pfmin = e1*((dm2/2./e1)**2 - amt2)/dm2
	      cpf0 = sgsl(0.d0)
	      cpf = sgsl(pfmin)
	      if(cpf.gt.0.0d0) then
		 dpf = 10.0d0
		 pf = pfmin
		 dspf = cpf*pf
		 sumpf = dspf
		 sumcf = cpf
		 do while((dspf/sumpf).gt.1.0d-4)
		    pf = pf + dpf
		    cpf = sgsl(pf)
		    dspf = cpf*pf
		    sumpf = dspf + sumpf
		    sumcf = cpf + sumcf
		 end do
!	print *,'pf,dspf,sumpf',pf,dspf,sumpf
		 pfmean = sumpf/sumcf
		 pfmean = pfmean*cosav
		 cpf = sgsl(pfmean)/cpf0
		 cpf = 0.5*cpf
		 ef = sqrt(pfmean**2 + amt**2)
		 weff = dm2/(2.0d0*(ef + pfmean))
!	weff = dm2/(2.0d0*ef)
!	print *,'w,ww,pfmin,cpf0',w,ww,pfmin,cpf0
!	print *,'pfmean,cpf,weff',pfmean,cpf,weff
		 w = weff
	      else
		 w = 0.
	      end if
	   end if
	end if

C       cm variables 
	pct=pnt 
	b=w/(amt+w) 
	g=(w+amt)/ww
	pcl=g*(pnl-b*ep)
	pcs=pcl**2+pct**2 
	pc=sqrt(pcs)
	cthc=pcl/pc 
	tc=acos(cthc) 
	return
	end



	
*       partdel
	subroutine partdel(amt,am1,am2,pn,tn,ajt,ajw)
C       partial derivatives
	implicit real*8 (a-h,o-z) 
	pi=acos(-1.0d0)
!       dt=pi/180.0d0
!       dp=10.0d0
	dt=pi/720.0d0
	dp=2.0d0
C       angle
	tnp=tn+dt 
	tnm=tn-dt 
	tn0=tn
	cpf = -1.
	call kindel(amt,am1,am2,pn,tnp,wp,tcp,cpf)
	cpf = -1.
	call kindel(amt,am1,am2,pn,tnm,wm,tcm,cpf)
	cpf = -1.
	call kindel(amt,am1,am2,pn,tn0,w,tc0,cpf)
	den=cos(tnp)-cos(tnm) 
	den=abs(den)
!	print *,'den,tcp,tcm,tc0',den,tcp,tcm,tc0
	if(den.gt.1.0d-3.and.(w*wp*wm.gt.0.))then
	   ajt=(cos(tcp)-cos(tcm))/den
	   ajt=abs(ajt) 
	else
	   ajt=(cos(tc0)-cos(tcm))/(cos(tn0)-cos(tnm))
	   ajt=abs(ajt) 
	endif 
C       energy 
	pnp=pn+dp 
	pnm=pn-dp 
	call kine(amt,am1,am2,pnp,tn,wp,tc) 
	call kine(amt,am1,am2,pnm,tn,wm,tc) 
	am12 = am1**2
	tp = sqrt(pnp**2 + am12) - am1
	tm = sqrt(pnm**2 + am12) - am1
	ajw=(wp-wm)/(pnp-pnm) 
!	ajw=(wp-wm)/(tp-tm) 
	ajw=abs(ajw)
	return
	end


	
*       oara
	subroutine oara(w,p,th,dsdwdp)
	implicit real*8 (a-h,m-z)
	integer n
	common/sg/ia
	common/qd/qdf 
      	common/del/ip 
        character*1 scal,fermi
 	common /fer/ fermi,scal
      	common/m/z,n
	common/kf/e1
	
	data mpi/0.140/
	data b/7.0d0/,c/9.0d0/,stot/125.0d0/
	data cp/8.0d3/,t0/124.0d-3/

	if(abs(ip).eq.1) then
	   dsdwdp = stot*c*exp(-b*p**2)
	else if (abs(ip).eq.2.or.ip.eq.0) then
	   t = sqrt(p**2+mpi**2) - mpi
	   dsdwdp = cp*exp(-t/t0)
	end if
	return
	end

