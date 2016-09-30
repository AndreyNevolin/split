TARGET = split
export BUILD_FLAGS =

SRCS = split.cpp

OUTDIR = .
OBJDIR = .objs
FULLTARGET = ${OUTDIR}/${TARGET}

OBJS = $(patsubst %.cpp, ${OBJDIR}/%.o, ${SRCS})

GCC = g++ -std=c++0x -Wall $(BUILD_FLAGS)

.PHONY: clean

default : BUILD_FLAGS += -s
default : ${FULLTARGET}

debug : BUILD_FLAGS += -DSPLIT_DEBUG -g -rdynamic
debug : ${FULLTARGET}

clean:
	-rm -f ${FULLTARGET} > /dev/null 2>&1
	-rm -f ${OBJS} > /dev/null 2>&1

${FULLTARGET}: ${OBJS}
	-mkdir -p ${OUTDIR} > /dev/null 2>&1
	${GCC} -o ${FULLTARGET} ${OBJS}

${OBJDIR}/%.o: %.cpp
	-mkdir -p ${OBJDIR} > /dev/null 2>&1
	${GCC} -c -o $@ $<
