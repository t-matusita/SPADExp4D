.PHONY: all cpp check clean

all: cpp

cpp:
	$(MAKE) -C cpp all

check:
	$(MAKE) -C cpp check

clean:
	$(MAKE) -C cpp clean
