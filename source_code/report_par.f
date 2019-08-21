      program report_parameters
      implicit none
      include 'par.for'
      real*8 beta_fs, Blambda, amp

      beta_fs = 0.0
      Blambda = ((U_1*lambda_z)/N_1)*sqrt((lambda_z/Radius))

      open(1, file = 'parameters.txt', status = 'unknown')

      write(1,*)"my_form=", my_form
      write(1,*)"type_cur=", type_cur
      write(1,*)"type_sec=", type_sec
      write(1,*)"U_1=", U_1
      write(1,*)"L_1=", L_1
      write(1,*)"N_1=", N_1
      write(1,*)"Radius=", Radius
      write(1,*)"lambda_z=", lambda_z
      write(1,*)"Re=", Re
      write(1,*)"Go2Re_base=", Go2Re_base
      write(1,*)"Go_base=", Go_base
      write(1,*)"Pr=", Pr
      write(1,*)"Blambda=", Blambda
      write(1,*)"beta=", beta
      write(1,*)"beta_fs=", beta_fs
      write(1,*)"amp_gv=", amp_gv
      write(1,*)"freq_gv=", fhz_gv
      write(1,*)"omega_gv=", omega_gv
      write(1,*)"freq_sec=", fhz_sec
      write(1,*)"omega=", omega
      write(1,*)"imax=", imax
      write(1,*)"jmax=", jmax
      write(1,*)"dx=", dx
      write(1,*)"dy=", dy
      write(1,*)"stf=", stf
      write(1,*)"dt=", dt
      write(1,*)"stpp=", stpp
      write(1,*)"tt=", tt
      write(1,*)"np=", np
      write(1,*)"i0=", i0
      write(1,*)"i1=", i1
      write(1,*)"i2=", i2
      write(1,*)"i3=", i3
      write(1,*)"i4=", i4
      write(1,*)"i5=", i5
      write(1,*)"kfour=", kfour
      write(1,*)"kphys=", kphys
      if (tt_base .EQ. 0) then
         write(1,*)"Não simula o baseflow2D"
      else
         write(1,*)"Simular o baseflow2D"
         write(1,*)"dt_base=", dt_base
         write(1,*)"tt_base=", tt_base
      end if
      stop
      endprogram
