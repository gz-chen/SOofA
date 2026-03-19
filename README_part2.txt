DESCRIPTION
	These are computer programmes for producing some of the results
	in the article "Stratum order-of-addition designs"
	by Liushan Zhou, Ze Liu, Min-Qian Liu and Guanzhou Chen.
	
CONTENTS
	Fig1and2.sh
		A Linux shell script file which prints Figures 1 and 2 into "Rplots.pdf".
		
	TabS8.sh
		A Linux shell script file which generates the data in Table S8,
		with u = 1, 2, 4, 8, 16.
		
	TabS9.sh
		A Linux shell script file which generates the data in Table S9,
		with N_train = 504, 72, 36, 24, 18, 9.
		
	TabS10.sh
		A Linux shell script file which generates the raw data for obtaining Table S10.
		
	main_genDesign.bin
		A programme for producing some of the known stratum OofA designs.
		To compile the programme, input the following command in the command line:
			make -f makefile genDesign
		The text file "help-main_genDesign.txt" illustrates the usage of this programme.
		
DEPENDENCIES
	The C++ code depends on the Eigen library
	(https://eigen.tuxfamily.org/index.php?title=Main_Page),
	with version = 3.3.9 (higher version may or may not work!).
	For convenience, a copy of Eigen-3.3.9 is provided in this directory.
	
AUTHOR
	Ze Liu (stat.liuze@mail.nankai.edu.cn)
	