SUBDIR =  HF jj2lsj90 jjgen90 rangular90  rbiotransform90  rci90  rcsfgenerate90  \
          rcsfinteract90 rcsfzerofirst90 rhfs90 rmcdhf90 rnucleus90 rtransition90  rwfnestimate90 sms90  \
	  rangular90_mpi rbiotransform90_mpi rci90_mpi rmcdhf90_mpi rtransition90_mpi

TARGETS	 = install

$(TARGETS):
	@for i in $(SUBDIR) ;\
	do \
		echo "Current directory: $$i" ; \
	                cd $$i ; \
			make; \
echo "....Leaving: $$i" ; echo; echo;\
		cd .. ; \
	done

clean :
	@for i in $(SUBDIR) ;\
	do \
	echo "Current directory: $$i" ; \
           cd $$i ; \
           make clean; \
           cd .. ; \
        done

