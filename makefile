V_cppstd=c++11
O_cppstd=-std=${V_cppstd}
O_optim=-O3
V_path_eigen=./eigen-3.3.9/
V_path_LZ_math=./LZ_cpp_lib-2.0-dev/

main.bin: ./src/main.cpp
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} \
		-D EIGEN_MPL2_ONLY ${O_optim} -o $@ $^
clean:
	rm -f main.bin
run:
	./main.bin
renew:
	make clean
	make
crun:
	make
	make run
debug:
	make clean
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} \
		-Wall -O0 -g -o main.bin ./src/main.cpp
#
genDesign:
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} \
		-D EIGEN_MPL2_ONLY -D LZ_DEF_compile_part_val_makefile=0 -o main_$@.bin ./src/main.cpp
mdljob:
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} \
		-D EIGEN_MPL2_ONLY -D LZ_DEF_compile_part_val_makefile=4 -o main_$@.bin ./src/main.cpp
variance:
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} \
		-D EIGEN_MPL2_ONLY -D LZ_DEF_compile_part_val_makefile=7 -o main_$@.bin ./src/main.cpp