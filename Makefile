export LIBMESH_DIR=../libmesh

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common

libmesh_CXXFLAGS += -fopenmp

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
muobjects	:= $(wildcard ../muparser/*.o)
###############################################################################



.PHONY: clean distclean

###############################################################################
# Target:
#
target 	   := ./shell-$(METHOD)


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) $(muobjects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)


echo:
	@echo "obj= $(objects) $(muobjects)"


# Useful rules.
clean:
	@rm -f $(objects) *~ .depend

distclean:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o

run: $(target)
	@$(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################
