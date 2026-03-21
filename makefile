PROJECT_NAME = nrgplusplus

.PHONY: default build clean

default: build 
	
build:
	mkdir -p build 
	cmake -G Ninja -Bbuild -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang .
	cmake -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=TRUE -Bbuild -DCMAKE_CXX_FLAGS=" -D_GLIBCXX_USE_TBB_PAR_BACKEND=0 -fsanitize=address,leak,undefined -Wno-narrowing -fsanitize-recover=all -g -O0" -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang .
	ninja -C build -j12
	


clean: 
	rm -rf build
	rm -rf CMakeCache.txt CMakeFiles Makefile cmake_install.cmake compile_commands.json
