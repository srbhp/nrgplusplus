projectName = twoChannelSiam

default: build 
	
build:
	rm -rf build 
	mkdir build
	# Clang for debug
	#cmake  -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=TRUE -Bbuild -DCMAKE_CXX_FLAGS="  -D_GLIBCXX_USE_TBB_PAR_BACKEND=0 -fsanitize=address,leak,undefined -Wno-narrowing -fsanitize-recover=all -g -O0"  -DCMAKE_CXX_COMPILER=clang++  .
	# GCC 
	#cmake  -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=TRUE -Bbuild -DCMAKE_CXX_FLAGS="  -D_GLIBCXX_USE_TBB_PAR_BACKEND=0 -fsanitize=address,leak,undefined -Wno-narrowing -fsanitize-recover=all -g -O0"  -DCMAKE_CXX_COMPILER=g++  .
	# Intel Compiler 
	cmake -Bbuild -DCMAKE_CXX_FLAGS=" -D_GLIBCXX_USE_TBB_PAR_BACKEND=0 -std=c++20 -qopenmp  -qmkl -fp-model precise " -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx  .
	# ninja -C build -j12
	make -C build -j12 all 
	# ninja -C build
	#make -C build/${projectName} -j12 all 
plot:
	#python plot.py
	#python3 spec.py
	#python3.10 rgflow.py 

clean: 
	rm -rf build 
