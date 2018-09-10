SUBDIR = lib appl tool
TARGETS	 = install
MKDIR = ${GRASP} ${GRASP}/bin ${GRASP}/lib

$(TARGETS):
	@for D in $(MKDIR) ; \
        do \
           test -d $$D || echo ... creating a new directory $$D... ; \
           test -d $$D || mkdir $$D ; \
        done
	@for i in $(SUBDIR) ;\
	do \
                cd $$i ; \
                echo "....Entering: $$i" ; \
                $(MAKE) $@ ; \
                cd .. ; \
                echo "....Leaving: $$i" ; echo; echo;\
        done
clean :
	@for i in $(SUBDIR) ;\
        do \
        echo "Current directory: $$i" ; \
           cd $$i ; \
           make clean; \
           cd .. ; \
        done

