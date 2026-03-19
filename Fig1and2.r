source("./src/plot.r");
#Figure 1:
fun_plot_eg_motivation();
#Figure 2:
mat_des_0=fun_des_read("./des_fig/OofAdes_28_8_2_2_1.txt");
mat_des_1=fun_des_perminv(fun_des_read("./des_fig/OofAdes_28_8_0_partCOA2.txt"));
fun_plot_strat_pairall_compare(mat_des_0,mat_des_1,c(2,2,2,2),"",1);
