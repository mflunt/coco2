Module mcmc_uncorr

contains

SUBROUTINE hbtdmcmc(beta, ad, ef, h_agg,y,n0, &
ad_param1, ad_param2, ef_param1, ef_param2,  sigma_model, sigma_measure, &
R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, &
para_temp, ad_pdf_all, ef_pdf_all, burn_in, &
stepsize_ad, stepsize_ef, ef_indices, nobs_gas1, &
nIt, nsub, nit_sub, &
nbeta, k, k_ef, nmeasure,  ydim1, ydim2, &
ad_out, ef_out, y_out, sigma_model_out, sigma_y_out, &
n0T_out, accept_ad, reject_ad, accept_ef, reject_ef, &
accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, &
tot_acc_ad, tot_acc_ef, tot_acc_sigma_y, &
accept_ad_all, reject_ad_all, accept_ef_all, reject_ef_all, &
accept_sigma_y_all, reject_sigma_y_all)

IMPLICIT NONE

! THIS SUBROUTINE IS TO PERFORM A HIERARCHICAL MCMC
! THIS CODE RELIES ON MEASUREMENTS BEING UNCORELLATED
! THIS ALLOWS MUCH FASTER CALCULATION
! THAT MEANS THIS SCRIPT IS GOOD FOR TESTING STUFF 

! REQUIRES:
! Inputs to be genrated in Python using script tdmcmc_template.py

! OUTPUTS:
! Will be passed back to tdmcmc_template.py

! PARALLEL TEMPERING:
! The inversion can be performed with parallel tempering or without
! The default is without, performed on a single processor

!!******************************************************************************
!! TO COMPILE: For Edinburgh:

!! 02/2022 openmp now should be -qopenmp
!!! f2py -L/geos/usrgeos/intel/Compiler/xe2016/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m test_nonlin --fcompiler=intelem --f90flags='-qopenmp' -liomp5 hbmcmc_nonlin.f90

!! INPUTS !!!!!!!!!!!!!!!!!!!
! Dimensions
INTEGER nbeta
INTEGER k
INTEGER k_ef
INTEGER nmeasure
INTEGER nit_sub
INTEGER nIt
INTEGER burn_in
INTEGER nsub
INTEGER ydim1
INTEGER ydim2
! Single Variables
INTEGER sigma_model_pdf
!INTEGER pdf_param1_pdf
!INTEGER pdf_param2_pdf
INTEGER para_temp
INTEGER nobs_gas1
! Input arrays
REAL stepsize_ad(k)
REAL stepsize_ef(k_ef)
INTEGER ad_pdf_all(k)
INTEGER ef_pdf_all(k_ef)
REAL beta(nbeta)
REAL ad(k, nbeta)  
REAL ef(k_ef, nbeta)               
REAL ad_param1(k)            
REAL ad_param2(k)
REAL ef_param1(k_ef)            
REAL ef_param2(k_ef)                
REAL h_agg(nmeasure,k)
REAL y(nmeasure) 
REAL n0(nmeasure, nbeta) 
REAL n0T(nbeta)
REAL sigma_y(nmeasure,nbeta)
REAL sigma_model(ydim2, nbeta)
REAL stepsize_sigma_y(ydim2)
INTEGER R_indices(ydim1,ydim2)
INTEGER ef_indices(k_ef)
REAL sigma_measure(nmeasure)
REAL error_structure(nmeasure) 
REAl sigma_model_hparam1(ydim2)
REAl sigma_model_hparam2(ydim2)
! Outputs
REAL ad_out(k,nit_sub) 
REAL ef_out(k_ef,nit_sub)  
REAL y_out(nmeasure,nit_sub)    
REAL sigma_y_out(nmeasure, nit_sub)
REAL sigma_model_out(ydim2,nit_sub)
REAL n0T_out(nit_sub)
REAL tot_acc_ad(k)
REAL tot_acc_ef(k_ef)
REAL tot_acc_sigma_y(ydim2)
INTEGER accept_ad(k)
INTEGER reject_ad(k)
INTEGER accept_ef(k_ef)
INTEGER reject_ef(k_ef)
INTEGER accept_swap, reject_swap
INTEGER accept_sigma_y, reject_sigma_y
INTEGER accept_ad_all(k,nbeta)
INTEGER reject_ad_all(k,nbeta)
INTEGER accept_ef_all(k_ef,nbeta)
INTEGER reject_ef_all(k_ef,nbeta)
INTEGER accept_sigma_y_all(nbeta), reject_sigma_y_all(nbeta)
! INTERMEDIATE VARIABLES
INTEGER it, ibeta, remain_it, pair1,pair2, ib, it_sub, remain     !remain_dim
INTEGER remain_swap,jj, kk
REAL u1,u2, randomu,pT_chain, beta1,beta2
REAL ad_it(k,nit_sub)   
REAL ef_it(k_ef,nit_sub)            
REAL y_it(nmeasure,nit_sub)                      
REAL sigma_model_ap(ydim2)           
REAL sigma_y_it(nmeasure,nit_sub), sigma_model_it(ydim2,nit_sub)     
REAL n0T_it(nit_sub)
REAL sigma_y_temp(nmeasure), y_error_temp(nmeasure)
INTEGER accept_batch_ad(k)
INTEGER reject_batch_ad(k)
INTEGER accept_batch_ef(k_ef)
INTEGER reject_batch_ef(k_ef)
INTEGER acc_y_batch(ydim2)
INTEGER rej_y_batch(ydim2)
REAL detval_temp
REAL C(nmeasure), sigma_yinv(nmeasure)
REAL n0_temp(nmeasure)
REAL n0T_temp
INTEGER ti
! SUBROUTINE INPUTS
REAL betaib, n0Tib
REAL adib(k), efib(k_ef), n0ib(nmeasure)           
REAL sigma_yib(nmeasure), sigma_modelib(ydim2)
! SUBROUTINE OUTPUTS
REAL n0Tib1                     
REAL adib1(k), efib1(k_ef), n0ib1(nmeasure)      
REAL sigma_yib1(nmeasure), sigma_modelib1(ydim2)
INTEGER reject_yib1, accept_yib1
INTEGER acceptadib1(k), rejectadib1(k), acc_badib1(k), rej_badib1(k)
INTEGER acceptefib1(k_ef), rejectefib1(k_ef), acc_befib1(k_ef), rej_befib1(k_ef)
INTEGER acc_byib1(ydim2), rej_byib1(ydim2)
REAL detval(nbeta)
REAL detvalib, detvalib1
REAL stepsize_sig_ib1(ydim2)
REAL stepsize_adib1(k)
REAL stepsize_efib1(k_ef)     

!! F2PY IN/OUT COMMANDS !!!!!!!
!f2py intent(in) beta,k, k_ef, ad, ef, h_agg,y,n0
!f2py intent(in) ad_param1, ad_param2, stepsize_ad
!f2py intent(in) ef_param1, ef_param2, ef_indices, stepsize_ef
!f2py intent(in) sigma_model, sigma_measure, para_temp, R_indices 
!f2py intent(in) sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y
!f2py intent(in) sigma_model_pdf, ad_pdf_all, ef_pdf_all, burn_in
!f2py intent(in) nIt,nsub,nit_sub, nobs_gas1 
!f2py intent(in) nbeta, nmeasure, ydim1, ydim2
!f2py intent(out) ad_out, ef_out, sigma_model_out, sigma_y_out,y_out
!f2py intent(out) n0T_out, reject_swap
!f2py intent(out) accept_ad, reject_ad, accept_swap
!f2py intent(out) accept_ef, reject_ef
!f2py intent(out) accept_sigma_y, reject_sigma_y
!f2py intent(out) tot_acc_sigma_y, tot_acc_ad, tot_acc_ef
!f2py intent(out) accept_ad_all, reject_ad_all
!f2py intent(out) accept_ef_all, reject_ef_all
!f2py intent(out) accept_sigma_y_all, reject_sigma_y_all

 
 !    call OMP_SET_NUM_THREADS(nbeta)     ! UNCOMMENT IF PARALLEL TEMPERING REQUIRED

  call init_random_seed()          ! Ensure random number generation starts from new point each time program is run
                                  ! Random seed only needs to be called once in a program.  
				    ! It is used to generate the first random number. 
				   ! Any random numbers generated after the first will use the previous random number as its seed.

! SET-UP INITIAL STARTING VALUES

accept_ad(:)=0
reject_ad(:)=0
accept_ef(:)=0
reject_ef(:)=0
accept_sigma_y=0
reject_sigma_y=0
accept_swap=0
reject_swap=0
it_sub=1

accept_ad_all(:,:)=0
reject_ad_all(:,:)=0
accept_ef_all(:,:)=0
reject_ef_all(:,:)=0
accept_sigma_y_all(:)=0
reject_sigma_y_all(:)=0

accept_batch_ad(:)=0
reject_batch_ad(:)=0
accept_batch_ef(:)=0
reject_batch_ef(:)=0
acc_y_batch(:)=0
rej_y_batch(:)=0

y_it(:,:)=0.

 do jj=1,ydim2   
        y_error_temp(R_indices(:,jj)) = sigma_model(jj,1)    ! Provided dim2 isn't too big then should be fine
 enddo  
   			   
 sigma_y_temp=sqrt(y_error_temp**2 + sigma_measure**2)      ! sigma_y is the combination of model and measurement error
 n0T(:)=sum((n0(:,1)/sigma_y_temp)**2)

 do ibeta = 1,nbeta
    sigma_y(:,ibeta) = sigma_y_temp
 enddo
 sigma_model_ap=sigma_model(:,1)
 detval(:) = sum(alog(sigma_y(:,1)))


! MCMC loop
!###############################################################
do it=1,(nIt+burn_in)
   
   !call random_number(u) 

 
       ! If doing fixed dimension MCMC:
       !remain_it = FLOOR(2*u) + 1    ! Choose random number between 1 and 2 - no reversible jump.
   remain_it = modulo(it,3)+1
   !remain_it = 1
  

!$OMP PARALLEL DO DEFAULT(SHARED) private(ibeta, betaib,adib,efib), &
!$OMP& private(n0ib,n0Tib,adib1,efib1,n0ib1,n0Tib1), &
!$OMP& private(acceptadib1, rejectadib1,acceptefib1, rejectefib1 ), &
!$OMP& private(sigma_yib, sigma_modelib, sigma_yib1, sigma_modelib1, accept_yib1, reject_yib1), &
!$OMP& private(detvalib, detvalib1),&
!$OMP& private(stepsize_adib1, acc_badib1, rej_badib1),&
!$OMP& private(stepsize_efib1, acc_befib1, rej_befib1),&
!$OMP& private(stepsize_sig_ib1, acc_byib1, rej_byib1),&
!$OMP& shared(ad,ef,n0,n0T, k, h_agg)
   do ibeta=1,nbeta


     if (para_temp .EQ. 1 .or. ibeta .EQ. 1) then 

       ! The following ib variables are necessary for PT
       betaib = beta(ibeta)
       efib  = ef(:,ibeta)
       adib  = ad(:,ibeta)
       !pdf_param1ib = pdf_param1(:,ibeta)
       !pdf_param2ib = pdf_param2(:,ibeta)

       n0ib = n0(:,ibeta)
       n0Tib = n0T(ibeta)

       sigma_yib = sigma_y(:,ibeta)
       sigma_modelib = sigma_model(:,ibeta)
       detvalib = detval(ibeta)

       
       if (remain_it .EQ. 1) then              ! EF UPDATE

            
            call ef_update(betaib,k, adib, efib, ef_param1,ef_param2, &
                         h_agg,n0ib,n0Tib,sigma_yib, stepsize_ef, ef_indices, k_ef, &
                         accept_batch_ef, reject_batch_ef, ef_pdf_all, it, burn_in,  nmeasure, nobs_gas1, &
                         efib1, n0ib1, n0Tib1, acceptefib1, rejectefib1, stepsize_efib1, acc_befib1, rej_befib1)


            ef(:,ibeta) = efib1
            n0(:,ibeta) = n0ib1
            n0T(ibeta) = n0Tib1 
            !pdf_param1(:,ibeta) = pdf_param1ib1
            !pdf_param2(:,ibeta) = pdf_param2ib1

            accept_ef_all(:,ibeta) = accept_ef_all(:,ibeta) + acceptefib1
            reject_ef_all(:,ibeta) = reject_ef_all(:,ibeta) + rejectefib1

            if (betaib .EQ. 1.) then 
               accept_ef(:) = accept_ef(:) + acceptefib1
               reject_ef(:) = reject_ef(:) + rejectefib1
               stepsize_ef=stepsize_efib1
               accept_batch_ef=acc_befib1
               reject_batch_ef=rej_befib1
            endif

          elseif (remain_it .EQ. 2) then  ! AD UPDATE

            call ad_update(betaib,k, adib, efib, ad_param1,ad_param2, &
                         h_agg,n0ib,n0Tib,sigma_yib, stepsize_ad, ef_indices, k_ef, &
                         accept_batch_ad, reject_batch_ad, ad_pdf_all, it, burn_in,  nmeasure, nobs_gas1,  &
                         adib1, n0ib1, n0Tib1, acceptadib1, rejectadib1, stepsize_adib1, acc_badib1, rej_badib1)


            ad(:,ibeta) = adib1
            n0(:,ibeta) = n0ib1
            n0T(ibeta) = n0Tib1 

            accept_ad_all(:,ibeta) = accept_ad_all(:,ibeta) + acceptadib1
            reject_ad_all(:,ibeta) = reject_ad_all(:,ibeta) + rejectadib1

            if (betaib .EQ. 1.) then 
               accept_ad(:) = accept_ad(:) + acceptadib1
               reject_ad(:) = reject_ad(:) + rejectadib1
               stepsize_ad=stepsize_adib1
               accept_batch_ad=acc_badib1
               reject_batch_ad=rej_badib1
            endif

          elseif (remain_it .EQ. 3) then  ! SIGMA_Y UPDATE
             
              call sigma_y_update(betaib, sigma_modelib, sigma_model_ap, sigma_measure, sigma_yib, detvalib, &
                 sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, R_indices, &
                 n0ib,n0Tib, acc_y_batch, rej_y_batch, it, burn_in, nmeasure, ydim1, ydim2, &
                 n0Tib1, accept_yib1, reject_yib1, sigma_yib1, sigma_modelib1, detvalib1, &
                 stepsize_sig_ib1, acc_byib1, rej_byib1) 

              sigma_y(:,ibeta) = sigma_yib1
              sigma_model(:,ibeta) = sigma_modelib1
              n0T(ibeta) = n0Tib1
              detval(ibeta) = detvalib1 

              accept_sigma_y_all(ibeta) = accept_sigma_y_all(ibeta) + accept_yib1
              reject_sigma_y_all(ibeta) = reject_sigma_y_all(ibeta) + reject_yib1

              if (betaib .EQ. 1.) then 
               accept_sigma_y = accept_sigma_y + accept_yib1
               reject_sigma_y = reject_sigma_y + reject_yib1
               stepsize_sigma_y=stepsize_sig_ib1
               acc_y_batch=acc_byib1
               rej_y_batch=rej_byib1
              endif

           endif     ! remain_it
           
      endif     !para_temp .EQ. 1 .or. ibeta .EQ. 1) 

   enddo    ! beta loop
!$OMP END PARALLEL DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Store xit and swap betas
  

   IF (it .GT. burn_in/2) THEN      ! Begin swaps after half of burn-in time
        remain_swap = modulo(it,2)
        IF (remain_swap .EQ. 1) THEN
          if (para_temp .EQ. 1) then
            call random_number(u1)   
            pair1 = FLOOR(nbeta*u1) + 1
            call random_number(u2)  
            pair2 = FLOOR(nbeta*u2) + 1
            if (pair1 .EQ. pair2) THEN
                if (pair2 .EQ. 1) then
                    pair2=pair1+1
                else 
                    pair2=pair1-1
                endif
            endif

            beta1=beta(pair1)*1.
            beta2=beta(pair2)*1.
            pT_chain = (beta2-beta1)*(n0T(pair2)/2.-n0T(pair1)/2.+detval(pair2)-detval(pair1))  ! detvals should be inverse determinants so signs this way round
            call random_number(randomu)
            if (alog(randomu) .LE. pT_chain) then           
                    beta(pair2)=beta1*1.         ! Uncomment for PT
                    beta(pair1)=beta2*1.         ! Uncomment for PT
                    accept_swap=accept_swap+1              
            else
                reject_swap=reject_swap+1
            endif      ! pT_chain if   
           else
                reject_swap=reject_swap+1
           endif   ! para_temp=1  
         ENDIF      ! reamin_swap =0 if
   ENDIF          ! it > burn_in/2
   


   IF (it .GT. burn_in) THEN     
        remain = modulo(it,nsub)          ! nsub typically = 100
        if (remain .EQ. 0) then

           do ib=1,nbeta                             ! uncomment for PT
               if (beta(ib) .EQ. 1.) then            ! uncomment for PT
          !        ib=1                   ! DEFAULT - AGAIN ib=1 comment and uncomment if statement if doing PT
                  ! STORE THE FOLLOWING VARIABLES AT THINNED nsub FREQUENCY
                  ad_it(:,it_sub)=ad(:,ib)
                  ef_it(:,it_sub)=ef(:,ib)
                  sigma_model_it(:,it_sub)=sigma_model(:,ib)
                  sigma_y_it(:,it_sub)=sigma_y(:,ib)
                  n0T_it(it_sub)=n0T(ib) 

                  do jj=1,k
                     if (ANY(ef_indices .EQ. jj)) then
                        y_it(1:nobs_gas1,it_sub) = y_it(1:nobs_gas1,it_sub) + &
                        (h_agg(1:nobs_gas1,jj)*ad(jj,ib))
                        do kk = 1,k_ef
                           if (ef_indices(kk) .EQ. jj) then
                              y_it(nobs_gas1+1:nmeasure,it_sub) = y_it(nobs_gas1+1:nmeasure,it_sub) + &
                              (h_agg(nobs_gas1+1:nmeasure,jj)*ad(jj,ib))*ef(kk,ib)
                              EXIT
                           endif
                        enddo
                     else
                         y_it(:,it_sub) = y_it(:,it_sub) + h_agg(:,jj)*ad(jj,ib)
                     endif
                  !   !!!!y_it(:,it_sub) = matmul(h_agg, ad(:,ib)*ef(:,ib))  
                  enddo
                  it_sub=it_sub+1
               endif                         ! uncomment for PT
            enddo                            ! uncomment for PT
        endif
   ENDIF           ! it >= burn_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo  ! It loop    ! END OF MCMC LOOP

! OUTPUTS
ad_out=ad_it
ef_out=ef_it
y_out = y_it
sigma_model_out=sigma_model_it
sigma_y_out=sigma_y_it
n0T_out=n0T_it
tot_acc_sigma_y=stepsize_sigma_y
tot_acc_ef=stepsize_ef
tot_acc_ad=stepsize_ad
END SUBROUTINE hbtdmcmc


SUBROUTINE ad_update(beta,k, ad, ef, ad_param1_all,ad_param2_all,  &
h_agg,n0,n0T,sigma_y, stepsize, ef_indices, k_ef, &
accept_batch, reject_batch, ad_pdf_all, it, burn_in, nmeasure, nobs_gas1,  &
ad_out, n0_out, n0T_out, accept_out, reject_out, stepsize_out, acc_batch_out, rej_batch_out) 

Implicit none 
INTEGER nmeasure, it, burn_in, k, k_ef
INTEGER nobs_gas1
REAL beta, n0T, n0T_out 
REAL av_acc
REAL ef(k_ef) 
REAL ad(k) 
REAL ad_out(k) 
REAL h_agg(nmeasure,k)   
REAL n0(nmeasure) 
REAL n0_out(nmeasure) 
REAL sigma_y(nmeasure)
REAL dy(nmeasure)
REAL n1(nmeasure) 
INTEGER ef_indices(k_ef)
INTEGER ad_pdf_all(k)
REAL ad_param1_all(k), ad_param2_all(k)
REAL pdf_param1        
REAL pdf_param2
REAL stepsize(k)
REAL accep_prob(k)
INTEGER accept(k), reject(k)
INTEGER accept_batch(k)
INTEGER reject_batch(k)
INTEGER accept_out(k), reject_out(k)
REAL stepsize_out(k)
REAL C(nmeasure)
INTEGER xi, ad_pdf, jj
REAL dx, n1T, pT, randomu, p0,p1
REAL ad_new(k)
REAL stepsize0
INTEGER acc_batch_out(k)
INTEGER rej_batch_out(k)
REAL aa,bb

 aa = 1
 bb = 0

accept=0
reject=0

! Propose new x values
! Currently this is done for all fixed x values and 5 elements of the variable regions
   ! CHANGE OF EMISSIONS VALUES

do xi=1,k     ! Loop through all elements of x
  
  ad_pdf = ad_pdf_all(xi)
  pdf_param1 = ad_param1_all(xi)
  pdf_param2 = ad_param2_all(xi)
  stepsize0 = stepsize(xi)

  dx = random_normal()*stepsize0

  ad_new = ad
  ad_new(xi) =ad(xi)+dx
  p0=0.
  p1=0.
  dy(:)=0.

  !! Need to find a way to loop through ef_indices
  ! Need to find what index of ef_indices it should be...
  
   if (ANY(ef_indices .EQ. xi)) then
      dy(1:nobs_gas1) = h_agg(1:nobs_gas1,xi)*dx 
      do jj = 1,k_ef
         if (ef_indices(jj) .EQ. xi) then
            dy(nobs_gas1+1:nmeasure)  = h_agg(nobs_gas1+1:nmeasure,xi)*ef(jj)*dx
            EXIT
         endif
      enddo
   else
      dy=h_agg(:,xi)*dx
   endif


   !if (ANY(ef_indices .EQ. xi)) then
   !  dy(1:nobs_gas1) = h_agg(1:nobs_gas1,xi)*dx 
   !  dy(nobs_gas1:nmeasure)  = h_agg(nobs_gas1:nmeasure,xi)*ef(xi-ef_indices(1)+1)*dx  ! Can't do this if not sequential
   !else
   !   dy=h_agg(:,xi)*dx  
   !endif

   n1=n0+dy

         n1T=sum((n1/sigma_y)**2)

         ! hyperparams are fixed below to be a single number - will only apply when x is a scaling of the prior
                  
         call calc_pdf(ad(xi),pdf_param1,pdf_param2,ad_pdf, p0)           ! Will apply whatever the PDF
         call calc_pdf(ad_new(xi),pdf_param1,pdf_param2,ad_pdf, p1)        
       
         pT = p1-p0 -0.5*(n1T - n0T)*beta ! All p1s already in log space

         !if (xi .EQ. 10) then
         !  !write(*,*) pT, n1T, n0T
         !  write(*,*) stepsize0
         !endif

         !if (xi .EQ. ef_indices(10)) then
         !  write(*,*) pT, n1T, n0T
         !endif

         if (ad_pdf .eq. 1) then
             if (ad_new(xi) .lt. pdf_param1) pT = -1.e20
             if (ad_new(xi) .gt. pdf_param2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN
             !ACCEPT
             ad(xi)=ad(xi)+dx
             n0=n1
             n0T=n1T
             if (it .GT. burn_in) then
                 accept(xi) = accept(xi) + 1
             endif
                                
         else
             !REJECT
             !if (beta .EQ. 1. .and. it .GT. burn_in) then 
             if (it .GT. burn_in) then 
                     reject(xi) = reject(xi) + 1
             endif
         endif   ! randomu condition

         if(beta .eq. 1. .and. it .le. burn_in) then
             if (alog(randomu) .LE. pT) THEN
                    !if (beta .EQ. 1. .and. it .GT. burn_in) accept_batch(xi) = accept_batch(xi) + 1
                    accept_batch(xi) = accept_batch(xi) + 1
                 
             else
                    reject_batch(xi) = reject_batch(xi) + 1
             endif

         endif

         if(beta .eq. 1 .and. it .le. burn_in .and. it .gt. 100 .and. modulo(it,500) .le. 1) then
         !if(beta .eq. 1 .and. it .le. burn_in .and. sum(accept_batch+reject_batch) .ge. 100*nIC1) then
                 !write(*,*) accept_batch
                 if(accept_batch(xi)+reject_batch(xi) .gt. 0) then
                     accep_prob(xi) = real(accept_batch(xi))/(accept_batch(xi)+reject_batch(xi))
                     !av_acc = min(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     if(accep_prob(xi) .lt. 0.2) stepsize(xi) = exp(alog(stepsize(xi)) - av_acc)
                     if(accep_prob(xi) .gt. 0.6) stepsize(xi) = exp(alog(stepsize(xi)) + av_acc)
                     accept_batch(xi) = 0
                     reject_batch(xi) = 0
                 endif
         endif

  enddo   ! xi loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ad_out=ad
n0_out=n0
n0T_out=n0T
accept_out=accept
reject_out=reject
stepsize_out=stepsize
acc_batch_out = accept_batch
rej_batch_out = reject_batch
END SUBROUTINE ad_update

SUBROUTINE ef_update(beta,k, ad, ef, ef_param1_all,ef_param2_all,  &
h_agg,n0,n0T,sigma_y, stepsize, ef_indices, k_ef, &
accept_batch, reject_batch, ef_pdf_all, it, burn_in, nmeasure, nobs_gas1,  &
ef_out, n0_out, n0T_out, accept_out, reject_out, stepsize_out, acc_batch_out, rej_batch_out) 

Implicit none 
INTEGER nmeasure, it, burn_in, k, k_ef
INTEGER nobs_gas1
REAL beta, n0T, n0T_out 
REAL av_acc
REAL ef(k_ef) 
REAL ad(k) 
REAL ef_out(k_ef) 
REAL h_agg(nmeasure,k)   
REAL n0(nmeasure) 
REAL n0_out(nmeasure) 
REAL sigma_y(nmeasure)
REAL dy(nmeasure)
REAL n1(nmeasure) 
INTEGER ef_indices(k_ef)
INTEGER ef_pdf_all(k_ef)
REAL ef_param1_all(k_ef), ef_param2_all(k_ef)
REAL pdf_param1        
REAL pdf_param2
REAL stepsize(k_ef)
REAL accep_prob(k_ef)
INTEGER accept(k_ef), reject(k_ef)
INTEGER accept_batch(k_ef)
INTEGER reject_batch(k_ef)
INTEGER accept_out(k_ef), reject_out(k_ef)
REAL stepsize_out(k_ef)
REAL C(nmeasure)
INTEGER xi, ef_pdf
REAL dx, n1T, pT, randomu, p0,p1
REAL ef_new(k_ef)
REAL stepsize0
INTEGER acc_batch_out(k_ef)
INTEGER rej_batch_out(k_ef)
REAL aa,bb

 aa = 1
 bb = 0

accept=0
reject=0

! Propose new x values
! Currently this is done for all fixed x values and 5 elements of the variable regions
   ! CHANGE OF EMISSIONS VALUES

do xi=1,k_ef     ! Loop through only elements of x where ef ratio updates apply
  
  ef_pdf = ef_pdf_all(xi)
  pdf_param1 = ef_param1_all(xi)
  pdf_param2 = ef_param2_all(xi)
  stepsize0 = stepsize(xi)

  dy(:)=0.
  dx = random_normal()*stepsize0

  ef_new = ef
  ef_new(xi) =ef(xi)+dx
  p0=0.
  p1=0.
         
         ! EF update changes should only apply to gas 2 
         dy(1:nobs_gas1)=0.
         dy(nobs_gas1:nmeasure)=h_agg(nobs_gas1:nmeasure,ef_indices(xi))*ad(ef_indices(xi))*dx 
         n1=n0+dy

         n1T=sum((n1/sigma_y)**2)


         ! hyperparams are fixed below to be a single number - will only apply when x is a scaling of the prior
                  
         call calc_pdf(ef(xi),pdf_param1,pdf_param2,ef_pdf, p0)           ! Will apply whatever the PDF
         call calc_pdf(ef_new(xi),pdf_param1,pdf_param2,ef_pdf, p1)        
       
         pT = p1-p0 -0.5*(n1T - n0T)*beta ! All p1s already in log space

         !if (xi .EQ. 10) then
         !  write(*,*) pT, n1T, n0T
         !endif

         if (ef_pdf .eq. 1) then
             if (ef_new(xi) .lt. pdf_param1) pT = -1.e20
             if (ef_new(xi) .gt. pdf_param2) pT = -1.e20
         endif

         call random_number(randomu)
         if (alog(randomu) .LE. pT) THEN
             !ACCEPT
             ef(xi)=ef(xi)+dx
             n0=n1
             n0T=n1T
             if (it .GT. burn_in) then
                 accept(xi) = accept(xi) + 1
             endif
                                
         else
             !REJECT
             !if (beta .EQ. 1. .and. it .GT. burn_in) then 
             if (it .GT. burn_in) then 
                     reject(xi) = reject(xi) + 1
             endif
         endif   ! randomu condition

         if(beta .eq. 1. .and. it .le. burn_in) then
             if (alog(randomu) .LE. pT) THEN
                    !if (beta .EQ. 1. .and. it .GT. burn_in) accept_batch(xi) = accept_batch(xi) + 1
                    accept_batch(xi) = accept_batch(xi) + 1
                 
             else
                    reject_batch(xi) = reject_batch(xi) + 1
             endif

         endif

         if(beta .eq. 1 .and. it .le. burn_in .and. it .gt. 100 .and. modulo(it,500) .le. 1) then
         !if(beta .eq. 1 .and. it .le. burn_in .and. sum(accept_batch+reject_batch) .ge. 100*nIC1) then
                 !write(*,*) accept_batch
                 if(accept_batch(xi)+reject_batch(xi) .gt. 0) then
                     accep_prob(xi) = real(accept_batch(xi))/(accept_batch(xi)+reject_batch(xi))
                     !av_acc = min(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     av_acc = max(0.01,1.0/sqrt(real(it/500))) !1.0/sqrt(real(accept(xi)+reject(xi)))
                     if(accep_prob(xi) .lt. 0.2) stepsize(xi) = exp(alog(stepsize(xi)) - av_acc)
                     if(accep_prob(xi) .gt. 0.6) stepsize(xi) = exp(alog(stepsize(xi)) + av_acc)
                     accept_batch(xi) = 0
                     reject_batch(xi) = 0
                 endif
         endif

  enddo   ! xi loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ef_out=ef
n0_out=n0
n0T_out=n0T
accept_out=accept
reject_out=reject
stepsize_out=stepsize
acc_batch_out = accept_batch
rej_batch_out = reject_batch
END SUBROUTINE ef_update

     
SUBROUTINE sigma_y_update(beta, sigma_model_current, sigma_model_ap, sigma_measure, sigma_y_current, &
detval_current,sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf, R_indices, &
n0,n0T, accept_batch, reject_batch, it, burn_in, nmeasure, dim1, dim2, &
n0T_out, accept_out, reject_out, sigma_y_out, sigma_model_out, detval_out, &
stepsize_sig_out, accept_batch_out, reject_batch_out) 

IMPLICIT NONE

! Dimensions
INTEGER nmeasure, accept, reject, it, burn_in, dim1, dim2
! Single Variables
REAL beta 
! Input arrays
REAL n0(nmeasure) 
REAL n0T
REAL detval_current
REAL sigma_measure(nmeasure)
REAL sigma_model_current(dim2)
REAL sigma_model_ap(dim2)
REAl sigma_y_current(nmeasure)
INTEGER R_indices(dim1,dim2)
REAL stepsize_sigma_y(dim2)
REAL sigma_model_hparam1(dim2)
REAL sigma_model_hparam2(dim2)
INTEGER sigma_model_pdf
! Outputs
REAL n0T_out, detval_out
REAL sigma_y_out(nmeasure)
REAL sigma_model_out(dim2)
INTEGER accept_out, reject_out
! Intermediate variables
INTEGER  yi, jj,ii
REAL randomu, dsigma_y, sigma_model_new, u, av_acc
REAL p0_sigma_y, p1_sigma_y, n1T, detval_new, pT
REAL y_error_new(nmeasure)
REAL sigma_y_new(nmeasure)
REAL stepsize_sig_out(dim2)
INTEGER accept_batch(dim2), reject_batch(dim2)
INTEGER accept_batch_out(dim2), reject_batch_out(dim2)
REAL accep_prob

accept=0
reject=0

 !do yi=1,dim2
    call random_number(u)   
    yi = FLOOR(dim2*u)+1
		
    ! Generate new value of sigma_y
    dsigma_y = random_normal()*stepsize_sigma_y(yi)*sigma_model_ap(yi)
    sigma_model_new = sigma_model_current(yi) + dsigma_y       

    call calc_pdf(sigma_model_current(yi), sigma_model_hparam1(yi), sigma_model_hparam2(yi), sigma_model_pdf, p0_sigma_y)
    call calc_pdf(sigma_model_new, sigma_model_hparam1(yi), sigma_model_hparam2(yi), sigma_model_pdf, p1_sigma_y)
	
    do jj=1,dim2   
        y_error_new(R_indices(:,jj)) = sigma_model_current(jj)    ! Provided dim2 isn't too big then should be fine
    enddo  

    ! Change one element of sigma_y
    y_error_new(R_indices(:,yi)) = sigma_model_new  ! R_indices = array of indices to which each sigma_y applies 
   			   
    sigma_y_new=sqrt(y_error_new**2 + sigma_measure**2)   

    n1T=sum((n0/sigma_y_new)**2)

    detval_new = sum(alog(sigma_y_new))           ! This is actually the inverse of the determinant

    ! Compute P1/P0 	
    pT= p1_sigma_y-p0_sigma_y - detval_new*beta + detval_current*beta - 0.5*(n1T - n0T)*beta       !*beta      ! -detval_new becasue it's inverse of true determinant
    if (sigma_model_pdf .eq. 1) then
       if (sigma_model_new .lt. sigma_model_hparam1(yi)) pT = -1.e20
       if (sigma_model_new .gt. sigma_model_hparam2(yi)) pT = -1.e20
    endif
  
    call random_number(randomu)     ! Generates uniformly distributed random number
    
    if(alog(randomu) .le. pT) then      
       !ACCEPT	
       sigma_model_current(yi) = sigma_model_new
       sigma_y_current = sigma_y_new
       detval_current = detval_new
       n0T=n1T
       !if(beta .eq. 1. .and. it .gt. burn_in) accept=accept + 1
       if(it .gt. burn_in) accept=accept + 1
       if(beta .eq. 1. .and. it .le. burn_in) accept_batch(yi)=accept_batch(yi) + 1
    else
       !;REJECT					
       !if(beta .eq. 1. .and. it .gt. burn_in) reject=reject + 1
       if(it .gt. burn_in) reject=reject + 1
       if(beta .eq. 1. .and. it .le. burn_in) reject_batch(yi)=reject_batch(yi) + 1
    endif

    if(beta .eq. 1 .and. it .le. burn_in .and. modulo(it,500) .le. 2) then
       if (it .gt. 1) then
    !if(beta .eq. 1 .and. it .le. burn_in .and. sum(accept_batch+reject_batch) .ge. 100) then
        av_acc =  max(0.01,1.0/sqrt(real(it)/500.))
        
        !write(*,*) pT, sigma_model_hparam1, sigma_model_current, sigma_model_new
        
        do ii=1,dim2 
          if(accept_batch(ii)+reject_batch(ii) .gt. 0) then          
             accep_prob = real(accept_batch(ii))/(accept_batch(ii)+reject_batch(ii))
             !if (ii .eq. 1) write(*,*) accep_prob
             if(accep_prob .lt. 0.2) stepsize_sigma_y(ii) = exp(alog(stepsize_sigma_y(ii)) - av_acc)
             if(accep_prob .gt. 0.6) stepsize_sigma_y(ii) = exp(alog(stepsize_sigma_y(ii)) + av_acc)
             accept_batch(ii) = 0
             reject_batch(ii) = 0
          endif
        enddo
       endif   ! it .gt. 1
    endif

n0T_out=n0T
sigma_y_out=sigma_y_current
sigma_model_out=sigma_model_current
accept_out=accept
reject_out=reject
detval_out=detval_current
stepsize_sig_out=stepsize_sigma_y
accept_batch_out=accept_batch
reject_batch_out=reject_batch
END SUBROUTINE sigma_y_update




FUNCTION random_normal() RESULT(fn_val)
IMPLICIT NONE
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, half =0.5, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal

subroutine calc_pdf(x,pdf_param1,pdf_param2,pdf,p1)     ! Subroutine for calculating P1 for specified PDF

 Implicit none
 ! Args
 integer,intent(in)   :: pdf
 real,intent(in)      :: x, pdf_param1, pdf_param2
 real,intent(out)     :: p1
 integer              :: positive
 real                 :: mu, sigma
 ! Locals
real,parameter        :: pi = 3.14159265        ! a fixed constant

!! ALL PDFS ARE IN LOG-SPACE
!  PDF .eq. 'UNIFORM'
  If(pdf .eq. 1) then                 ! Uniform PDF = 1	

     if (x .lt. pdf_param1 .OR. x .gt. pdf_param2) then
        p1 = -1.e20
     else
        p1= alog(1./(pdf_param2 - pdf_param1))
     endif
     positive=0
  Endif

!  PDF .eq. 'GAUSSIAN'
  If(pdf .eq. 2) then                ! Gaussian PDF = 2
     p1=alog(1./(pdf_param2*sqrt(2.*pi)))+(-1.*(x - pdf_param1)**2 / 2./(pdf_param2**2))
     positive=0
  Endif

!  PDF .eq. 'LOGNORMAL'
  If(pdf .eq. 3) then                    ! Lognormal PDF = 3

    mu = alog(pdf_param1) - 0.5*alog(1. + pdf_param2**2/pdf_param1**2)
    sigma=sqrt(alog((pdf_param2/pdf_param1)**2 + 1.))


     p1=alog(1./(x*sigma*sqrt(2.*pi)))+( -1.*(alog(x) - mu)**2 / 2./(sigma**2))
     positive=1
  Endif

!  PDF .eq. 'EXPONENTIAL'
  If(pdf .eq. 4) then                    ! Exponential PDF = 4
     p1=alog(pdf_param1)+(-1.*pdf_param1*x)
     positive=1
  Endif


  if((x .le. 0.) .and. (positive .eq. 1)) p1 = -1.e20       ! Since log(0) isn't defined need to make p1 a very large negative number to ensure proposal is rejected

end subroutine calc_pdf

subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s, getpid
            integer(8) :: count, tms
           
            call random_seed(size = n)
            allocate(seed(n))

            ! First try if the OS provides a random number generator
            open(unit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
end subroutine init_random_seed

End Module mcmc_uncorr

