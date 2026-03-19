#This Linux shell script file generates the raw data for obtaining Table S10.
if [ ! -e "main_variance.bin" ]; then
	make -f makefile variance
fi
./main_variance.bin
