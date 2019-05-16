package PoolHap;

import java.io.PrintWriter;
import java.io.File;
import java.io.FileWriter;
import PoolHap.Parameters.GenParameters;

import java.time.format.DateTimeFormatter;  
import java.time.LocalDateTime;	

public class Main {

	public static void main(String[] args) throws Exception {
		
		GenParameters gp = new GenParameters(args[0]);
		String parameter_file = args[0];
		String prefix = args[1];
		String gs_var_pos = gp.inter_dir + prefix + "_vars.intra_freq.txt";
		String dc_out_file = gp.inter_dir + prefix + "_dc_plan.txt";
		int num_pools = Integer.parseInt(args[2]); 
		
		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");  

		System.out.println("PHX Initiation: " + dtf.format(LocalDateTime.now()));  

		HapConfig[] initial_local_haps = new HapConfig[num_pools];
		PrintWriter in_list = new PrintWriter(new FileWriter(new File(gp.inter_dir + prefix + "_p.in.list"))); 
		for (int p = 0; p < num_pools; p++) {
			GraphColoring pool_in = new GraphColoring(gp.inter_dir + prefix + "_p" + p + ".vef", gs_var_pos, gp.inter_dir + prefix + "_p" + p + ".in"); 
			System.out.println("Graph colouring for pool " + p + " is finished.");
			in_list.println(gp.inter_dir + prefix + "_p" + p + ".in"); 
			initial_local_haps[p] = pool_in.hapOut();
		}
		in_list.close();
		System.out.println("\nGC Finished: " + dtf.format(LocalDateTime.now()) + "\n");

		DivideConquer dc_maker = new DivideConquer(gs_var_pos, gp.inter_dir + prefix + "_p.in.list", parameter_file, dc_out_file);
		System.out.println("DC Finished: " + dtf.format(LocalDateTime.now()) + "\n");
		
		HapConfig final_global_haps; 
		if (dc_maker.region_solve) {
			HapConfig[] level_I_config = dc_maker.analyze_regions(dc_maker.regions_level_I, parameter_file, gp.inter_dir + prefix, 1);
			System.out.println("Level 1 solving Finished: " + dtf.format(LocalDateTime.now()) + "\n");
			HapConfig[] level_II_config = dc_maker.analyze_regions(dc_maker.regions_level_II, parameter_file, gp.inter_dir + prefix, 2);
			System.out.println("Level 2 solving Finished: " + dtf.format(LocalDateTime.now()) + "\n");  
	
		
			GraphColoring region_linker = new GraphColoring(level_I_config, level_II_config, gs_var_pos, gp.fragments);
			final_global_haps = region_linker.hapOut();
			final_global_haps.write_global_file_string(gp.out_dir + prefix + "_gc.inter_freq_vars.txt", false);
			System.out.println("\nGC Finished: " + dtf.format(LocalDateTime.now()) + "\n");
		} else {
			final_global_haps = new HapConfig(initial_local_haps); // new HapConfig(gp.inter_dir + prefix + "_p" + 0 + ".in",gs_var_pos,num_pools); // 
		}
		
		HapConfig[] final_local_haps = new HapConfig[num_pools]; 
		for (int p = 0; p < num_pools; p++) {
			HapLASSO inpool_lasso = new HapLASSO(p, gp.lambda, final_global_haps, gp.final_cutoff, "500000000", gp.inter_dir + prefix + "_p" + p);
			inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + prefix + "_p" + p + ".vef", null, gp.lasso_weights);
			double new_penalty = gp.lambda; 
			while (inpool_lasso.r2 == gp.min_r2) {
				new_penalty -= gp.lasso_penalty_step; 
				System.out.println("Local full-length LASSO has failed. Adjust the lambda penalty to " + new_penalty + " and trying again.");
				inpool_lasso.lambda = new_penalty; 
				inpool_lasso.estimate_frequencies_lasso(gp.inter_dir + prefix + "_p" + p + ".vef", null, gp.lasso_weights);
			}
			final_local_haps[p] = inpool_lasso.hapOut();

			System.out.println(final_local_haps[p].num_global_hap + " haplotypes were generated for pool " + p + ".");
		}
		System.out.println("\nLASSO Finished: " + dtf.format(LocalDateTime.now()) + "\n");

		HapConfig final_reconstruction = new HapConfig(final_local_haps); 
		System.out.println("Haplotype reconstruction finished. There are " + final_reconstruction.num_global_hap + " global haplotypes.");
		final_reconstruction.write2files(gp.out_dir + prefix + ".inter_freq_vars.txt", gp.out_dir + prefix + ".intra_freq.txt", "string", false);
	}
}
