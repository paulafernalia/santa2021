# Name of the final executable
TARGET              = release_exe
DEBUGTARGET         = debug_exe

# Decide whether the commands will be shwon or not
VERBOSE 			= FALSE

# Name the compiler
CCC 				= g++
CSTD 				= c++11
CPU 				= 64

# Compiler options 
CCFLAGS 			= -m$(CPU) -std=$(CSTD) -Wno-sign-compare
CCDEBUG				= -O0 -g -DDEBUG
CCRELEASE 			= -O3 -DNDEBUG
CPPLIB 				= -L$(GUROBI_HOME)/lib -lgurobi91 -lgurobi_g++4.2 -lm
INC                 = $(GUROBI_HOME)/include/

RMDIR 			    = rm -rf
ERRIGNORE 		    = 2>/dev/null

# Hide or not the calls depending of VERBOSE
ifeq ($(VERBOSE),TRUE)
	HIDE 			=  
else
	HIDE 			= @
endif

.PHONY: all clean debug help build

run: $(TARGET)
	./$(TARGET)

all: build
	./$(DEBUGTARGET)

build: $(DEBUGTARGET)

debug: build
	lldb $(DEBUGTARGET)

$(TARGET): santa2021.cpp
	$(CCC) $(CCFLAGS) $(CCRELEASE) -o $@ $< -I$(INC) $(CPPLIB) -lm

$(DEBUGTARGET): santa2021.cpp
	$(CCC) $(CCFLAGS) $(CCDEBUG) -o $@ $< -I$(INC) $(CPPLIB) -lm

clean:
	$(HIDE)$(RMDIR) *.o *_c++ *.class *.log *.lp *.bas *.dSYM $(ERRIGNORE);
	$(HIDE)$(RMDIR) $(DEBUGTARGET) $(TARGET) $(ERRIGNORE);
	$(HIDE)$(RMDIR) *.aux *.fls *.fdb_latexmk *.gz $(ERRIGNORE);
	@echo cleaning done!
