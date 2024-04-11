PKGSRC  := $(shell basename `pwd`)
RM = rm -rf
RCMD = R --vanilla CMD
RSCRIPT = Rscript --vanilla

all: make-db clean

db:
	@ $(RSCRIPT) create-sqlite-db.R

clean:
	@ git clean -dfn