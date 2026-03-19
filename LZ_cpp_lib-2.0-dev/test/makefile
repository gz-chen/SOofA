V_cppstd=c++11
O_cppstd=-std=${V_cppstd}
O_optim=-O0
#O_thread=-lpthread
O_thread=
V_path_eigen=/usr/local/include/eigen-3.3.9/
V_path_LZ_math=/home/liuze/LZ/LZ_Cpp/LZ_H/LZ_cpp_lib-2.0-dev/
main.bin: main.cpp
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} ${O_optim} ${O_thread} -o $@ $^
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
	${CXX} ${O_cppstd} -I ${V_path_LZ_math} -I ${V_path_eigen} ${O_optim} -g -o main.bin main.cpp