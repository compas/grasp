

SUBDIR = libmod lib9290 libmcp90 librang90 mpi90 libdvd90
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
clean:
	@for i in $(SUBDIR) ;\
        do \
                echo "Current directory: $$i" ; \
                        cd $$i ; \
                        $(MAKE) $@ ; \
                cd .. ; \
        done
