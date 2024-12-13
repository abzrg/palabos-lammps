all: build
	cmake --build build


init:
	cmake -S . -B build


build:
	$(MAKE) init


clean:
	cmake --build build --target clean
