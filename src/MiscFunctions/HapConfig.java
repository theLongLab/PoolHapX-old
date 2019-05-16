package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import MiscFunctions.LocusAnnotation;

public class HapConfig {
	
	/*
	 * @author  Quan Long. Oct 12, 2018
	 * 
	 * Recording the configuration of the haplotypes, both in-pool and global.
	 * 
	 * The objects can be a small region or a whole chromosome.  
	 * The object of this class will be input or output parameters of other relaying algorithms.
	 * 
	 * Corresponding to the data structure in this class, there are two files:
	 * (1) global_hap_file: 
	 * 		(1.1) Each row represents a locus; each column represents a haplotype. 
	 * 	   	(1.2) The first row is the header of the IDs of all haps (usually are indexes); 
	 * 			  the second row is the global frequency (or NaN if unknown).
	 * 		(1.3) The first column is the IDs of all loci.
	 * 				(1.3.1) chr_index;start_loc;end_loc;alleles (the alleles are separated by ":"). chr_index starts with zero. 
	 * 				(1.3.2) NOTE: if it is an indel, start_loc=end_loc; Here start_loc!=end_loc only if it is a region, instead of a primitive locus.
	 * 		(1.3) Use '\t' to separate columns.   
	 * (2) local_hap_file: local_haplotypes
	 * 		(2.1) Each row represents the in-pool frequencies of all haplotypes (many of them can be zero)
	 * 			the first row is the header of all hap IDs (in the global file)
	 * 			the first element is the pool ID 
	 * 			the rest elements are the frequencies of this haplotype in all pools: zero indicates the absence.
	 * 		(2.2) Each column represents a haplotype's frequencies in all pools.
	 * 			the first column is the IDs of all pools.
	 * 		(2.3) Use '\t' to separate columns. 
	 * 		(2.4) In general, it is recommended that the order of haplotypes in this file is the same as the order of 
	 * 			the haplotypes in the global hap file. But if it is not the case, the program can handle it.
	 * 		(Note that this matrix representation is convenient for analysis of relationship between pools, e.g., transmission or evolution.) 
	 * 
	 * Note that the information of how to encode alleles to numbers will be generated on-the-fly
	 * in the future, we may add a function to read them from a file.  
	 * 
	 */
	
	// Redundant variables for convenience
	public int num_global_hap;				// number of global haplotype in this region
	public int num_loci;					// number of loci in this region
	public int num_pools;					// number of pools under study	
	public HashMap<String, Integer> hapID2index;  // map the pool IDs to their indexes in this.hap_IDs. 
	// Main structures
	public String[] hap_IDs;				// #global_hap
	public String[][] global_haps_string; 	// #global_hap x #loci: haplotypes that ever show up in the total population. String coded 
	public double[][] global_haps; 			// #global_hap x #loci: haplotypes that ever show up in the total population. floating number coded. 
	public double[] global_haps_freq; 		// #global_hap
	public String[] pool_IDs;				// #pools 		
	public double[][] in_pool_haps_freq;	// #global_hap x #pools  
	public LocusAnnotation[] locusInfo; 	// #loci. Note that a locus can be a SNP or a region. 
	public double[][] inpool_site_freqs;	// #loci x #pools; Added by Quan Dec. 2018.
	public double[] prop_reconstruction;	// The proportion of the reference sequence length that was reconstructed at all. 
	// Summary statistics
	public double[] mu;						// The 'average' allele at each position i.e. average haplotype Average frequency of alternate allele.
	public double[][] sigma;				// The variance-covariance matrix for each position.
	public double logL; 					// The log-likellihood of this HapConfig explaining the observed variant data.
	public int est_ind_pool; 				// Estimated number of individuals per pool (for now, assumed same per pool but easily extended). Used to calculate logL.
	// Link to the solving method
	
	// Constructors
	
	/*
	 * Constructor, forming the object by the variables in the memory.
	 * This constructor is useful in the divide & conquer algorithm.
	 */
	public HapConfig(String[][] global_haps_string, double[] global_haps_freq, double[][] in_pool_haps_freq, double[][] inpool_site_freqs, LocusAnnotation[] locusInput, int num_pools, String[] hap_IDs, String[] pool_IDs, int est_ind_pool){
		this.num_global_hap=global_haps_freq.length;
		this.num_loci=locusInput.length;
		this.global_haps_freq=global_haps_freq.clone();
		Boolean setID = false;
		if(hap_IDs!=null){ // *** Changed from this.hap_IDs (global) to hap_IDs (parameter)  
			this.hap_IDs=hap_IDs.clone();
		}else{ // if no IDs assigned, use indexes.
			this.hap_IDs = new String[this.num_global_hap];
			setID = true;
		}
		this.global_haps_string=global_haps_string.clone();
		for(int k=0;k<this.num_global_hap;k++) {
			this.global_haps_string[k]=global_haps_string[k].clone();
			if (setID) {
				String id = ""; 
				for (int p = 0; p < this.num_loci; p++) id += global_haps_string[k][p]; 
				this.hap_IDs[k] = id;
			}
			// System.out.println(this.hap_IDs[k]);
		}
		if(in_pool_haps_freq!=null){ // *** Matched format as above. 
			this.in_pool_haps_freq=in_pool_haps_freq.clone();
			for(int j=0;j<this.num_global_hap;j++)
				for(int k=0;k<this.num_pools;k++)
					this.in_pool_haps_freq[j][k]=in_pool_haps_freq[j][k];
		}else{ // if no intra-pool frequencies provided, start at 0. 
			this.in_pool_haps_freq = new double[this.num_global_hap][this.num_pools]; 
		}
		if(inpool_site_freqs!=null){
			this.inpool_site_freqs=inpool_site_freqs.clone();
		}
		this.locusInfo=locusInput;  // Originally, clone(), which is not a deep-clone. But we assume that the locusInfo won't be changed in the algorithm.
		// System.out.print("\n" + locusInput.length + "\t" + this.num_loci + "\t" + this.locusInfo.length);
		/*
		/for (int v = 0; v < this.num_loci; v++) {
			System.out.println(v);
			System.out.println("HapConfig constructor input:\t" + locusInput[v].alleles_coding);
			this.locusInfo[v].alleles_coding = new HashMap<String, Float>();
			System.out.println("HapConfig constructor pre:\t" + this.locusInfo[v].alleles_coding.size());
			for (String k : locusInput[v].alleles_coding.keySet()) {
				System.out.println("setting:\t" + locusInput[v].alleles_coding.get(k));
				this.locusInfo[v].alleles_coding.put(k, locusInput[v].alleles_coding.get(k));
			}
			System.out.println("HapConfig constructor post:\t" + this.locusInfo[v].alleles_coding.size());
		} // Need to specifically copy the mappings over to the new locus storage object. Shallow copy alone doesn't work.
		*/
		this.construct_hapID2index_map();
		this.encoding_haps(); // initialize this.global_haps;
		this.num_pools=num_pools;
		if(pool_IDs!=null){ // *** Matched format as above. 
			this.pool_IDs=pool_IDs.clone();
		}else{ // if no IDs assigned, use indexes.
			this.pool_IDs = new String[this.num_pools];
			for(int k=0;k<this.num_pools;k++)
				this.pool_IDs[k]=k+"";
		}
		this.est_ind_pool = est_ind_pool; 
		this.update_sigma_mu_logL();
	}
	
	public HapConfig(HapConfig[] final_local_haps){
		boolean[] skip_pools = new boolean[final_local_haps.length]; 
		HashSet<Integer> all_pos = new HashSet<Integer>();
		ArrayList<HashMap<Integer,Integer>> pool2posind = new ArrayList<HashMap<Integer,Integer>>();  
		HashMap<Integer,LocusAnnotation> pos2loc_anno = new HashMap<Integer,LocusAnnotation>();
		for (int p = 0; p < final_local_haps.length; p++) {
			pool2posind.add(new HashMap<Integer,Integer>());
			if (final_local_haps[p].num_global_hap == 0) skip_pools[p] = true;
			// TODODONE Amass a list of positions here from the LocusAnnotations.
			else
				for (int i = 0; i < final_local_haps[p].num_loci; i++) {
					if (!all_pos.contains(final_local_haps[p].locusInfo[i].start_loc)) {
						all_pos.add(final_local_haps[p].locusInfo[i].start_loc);
						pos2loc_anno.put(final_local_haps[p].locusInfo[i].start_loc, final_local_haps[p].locusInfo[i]);
					}
					// System.out.println(p + "\t" + final_local_haps[p].locusInfo[i].start_loc + "\t" + i);
					pool2posind.get(p).put(final_local_haps[p].locusInfo[i].start_loc,i); 
				}
		}
		int[] all_pos_final = all_pos.stream().mapToInt(Integer::intValue).toArray();
		Arrays.sort(all_pos_final);
		this.num_loci = all_pos_final.length; 
		ArrayList<String> global_ids = new ArrayList<String>(); 
		ArrayList<ArrayList<String>> local_ids = new ArrayList<ArrayList<String>>(); 
		ArrayList<Double> tmp_recon = new ArrayList<Double>();
		for (int p = 0; p < final_local_haps.length; p++) {
			local_ids.add(new ArrayList<String>());
			for (int h = 0; h < final_local_haps[p].num_global_hap; h++) {
				if (skip_pools[p]) continue;	// If this pool failed to be reconstructed, skip it entirely.
				String curr_vc = "";
				// TODODONE Make accurate strings for each partial haplotype.
				for (int l = 0; l < this.num_loci; l++) {
					if (pool2posind.get(p).containsKey(all_pos_final[l])) curr_vc += final_local_haps[p].global_haps_string[h][pool2posind.get(p).get(all_pos_final[l])];
					else curr_vc += "*"; 
				}
				if (!global_ids.contains(curr_vc)) {
					global_ids.add(curr_vc);
					tmp_recon.add(final_local_haps[p].prop_reconstruction[h]); 
				}
				local_ids.get(p).add(curr_vc); 
			}
		}
		this.num_global_hap = global_ids.size();
		this.hap_IDs = new String[this.num_global_hap];
		this.num_pools = final_local_haps.length;
		this.prop_reconstruction = new double[this.num_global_hap]; 
		for (int h = 0; h < this.num_global_hap; h++) this.prop_reconstruction[h] = tmp_recon.get(h);
		this.global_haps_string = new String[this.num_global_hap][this.num_loci];
		this.in_pool_haps_freq = new double[this.num_global_hap][final_local_haps.length];
		 // TODODONE Manage this part properly.
		this.locusInfo = new LocusAnnotation[this.num_loci];
		double[] tmp_global_freq = new double[this.num_global_hap]; 
		for (int h = 0; h < this.num_global_hap; h++) {
			String[] hap_var_comp = global_ids.get(h).split(""); 
			for (int l = 0; l < this.num_loci; l++) {
				this.global_haps_string[h][l] = hap_var_comp[l];
				this.locusInfo[l] = pos2loc_anno.get(all_pos_final[l]);
			}
			for (int p = 0; p < this.num_pools; p++) {
				if (skip_pools[p]) continue;
				else 
					if (local_ids.get(p).contains(global_ids.get(h))) {
						int inpool_index = local_ids.get(p).indexOf(global_ids.get(h));
						this.in_pool_haps_freq[h][p] = final_local_haps[p].global_haps_freq[inpool_index];
						tmp_global_freq[h] += final_local_haps[p].global_haps_freq[inpool_index]; 
					}
			}
			this.hap_IDs[h] = h + "";
		}
		this.global_haps_freq = new double[this.num_global_hap];
		for (int h = 0; h < this.num_global_hap; h++) this.global_haps_freq[h] = tmp_global_freq[h] / this.num_pools; 
		try {
			this.inpool_site_freqs = final_local_haps[0].inpool_site_freqs.clone();
		} catch (NullPointerException e) {} 
		this.construct_hapID2index_map();
		this.encoding_haps();
		this.pool_IDs = new String[this.num_pools];
		for (int p = 0; p < this.num_pools; p++) this.pool_IDs[p] = p + "";
		this.est_ind_pool = 0; 
		// this.update_sigma_mu_logL();
	}

	/*
	 * Constructor that reads information from files. 
	 * This constructor is for the connection between modules (such as AEM, rjMCMC, etc.)
	 */
	public HapConfig(String global_hap_input_file, String in_pool_hap_input_file){
		try{
			// parse global_hap_input_file
			BufferedReader br=new BufferedReader(new FileReader(global_hap_input_file));
			String[] header_ids=br.readLine().split("\t");
			String[] header_freqs=br.readLine().split("\t");
			this.num_global_hap=header_ids.length-1;
			this.hap_IDs=new String[this.num_global_hap];
			this.global_haps_freq=new double[this.num_global_hap];
			for(int h=0;h<num_global_hap;h++){
				this.hap_IDs[h]=header_ids[h+1];
				this.global_haps_freq[h]=Double.parseDouble(header_freqs[h+1]);
			}
			String line=br.readLine();
			while(line!=null && !line.contains("Recombinate")){
				this.num_loci++;
				line=br.readLine();
			}br.close();
			// System.out.println(this.num_loci);
			this.global_haps_string=new String[this.num_global_hap][this.num_loci];
			this.locusInfo=new LocusAnnotation[this.num_loci];
			br=new BufferedReader(new FileReader(global_hap_input_file));
			line=br.readLine();line=br.readLine(); // skip two headers
			line=br.readLine();
			int loci_index=0;
			while(line!=null){
				String[] tmp=line.split("\t");
				if (tmp[0].contains("Recombinate")) {
					this.prop_reconstruction = new double[this.num_global_hap]; 
					for(int h_index=0;h_index<this.num_global_hap;h_index++) this.prop_reconstruction[h_index] = Double.parseDouble(tmp[h_index + 1]); 
				} else {
					this.locusInfo[loci_index]=new LocusAnnotation(tmp[0]); // the first column is the locus-info
					for(int h_index=0;h_index<this.num_global_hap;h_index++)
						this.global_haps_string[h_index][loci_index]=tmp[h_index+1];
					// System.out.println(this.locusInfo[loci_index].alleles_coding.get("1"));
					loci_index++;
					
				}
				line=br.readLine();
				// System.out.println(line);
			}br.close();
			this.construct_hapID2index_map();
			this.encoding_haps();  // initialize this.global_haps;
			// then read the local-hap file
			if (in_pool_hap_input_file == null) {
				this.num_pools = 20; 
				this.pool_IDs=new String[this.num_pools];
				for (int p = 0; p < this.num_pools; p++) this.pool_IDs[p] = p + "";
				this.in_pool_haps_freq=new double[this.num_global_hap][this.num_pools];
			} else {
				br=new BufferedReader(new FileReader(in_pool_hap_input_file));
				String[] header_hap_IDs=br.readLine().split("\t"); // the first row is the hap IDs.
				if (header_hap_IDs[0].equals("Hap_ID")) {	// A *_haps.intra_freq.txt (in-pool frequencies of haplotypes) has been provided.
					if(this.num_global_hap!=header_hap_IDs.length-1)
						System.out.println("WRONG: this.num_global_hap is not the same between two files!");
					line=br.readLine();
					while(line!=null){ // count how many pools
						this.num_pools++;
						line=br.readLine();
					}br.close();
					this.pool_IDs=new String[this.num_pools];
					this.in_pool_haps_freq=new double[this.num_global_hap][this.num_pools];
					br=new BufferedReader(new FileReader(in_pool_hap_input_file));
					line=br.readLine(); //read the file again, skip the header 
					line=br.readLine();
					int pool_index=0;
					while(line!=null){
						String[] tmp=line.split("\t");
						this.pool_IDs[pool_index]=tmp[0]; // assign pool IDs
						for(int h=1;h<header_hap_IDs.length;h++){
							// System.out.println(h + "\t" + tmp[h] + "\t" + this.hapID2index.get(tmp[h]));
							int h_index=this.hapID2index.get(header_hap_IDs[h]); // based on the hap-ID, find out the hap-index
							this.in_pool_haps_freq[h_index][pool_index]=Double.parseDouble(tmp[h]);
						}
						pool_index++;
						line=br.readLine();
					}
					br.close();
				} else { // A *_vars.intra_freq.txt (in-pool frequencies of alternate alleles) has been provided.
					this.num_pools = header_hap_IDs.length - 1; 
					this.pool_IDs=new String[this.num_pools];
					for (int p = 0; p < this.num_pools; p++) this.pool_IDs[p] = header_hap_IDs[p + 1];
					this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
					int locus_index = 0; 
					line = br.readLine();
					while(line != null){
						String tmp[] = line.split("\t"); 
						for (int p = 0; p < this.num_pools; p++) this.inpool_site_freqs[locus_index][p] = Double.parseDouble(tmp[p + 1]);
						locus_index++;
						line = br.readLine();
					} br.close();
				}
			}
		}catch(Exception e){e.printStackTrace(); this.num_global_hap = 0;}
	}
	
	public void construct_hapID2index_map(){
		this.hapID2index = new HashMap<String, Integer>();
		for(int h=0;h<this.num_global_hap;h++) {
			this.hapID2index.put(this.hap_IDs[h], h);
			// System.out.println(this.hap_IDs[h] + "\t" + h);
		}
	}
	
	/*
	 * Maps global_haps_string to global_haps (that is coded by integers). 
	 * 
	 * It will be done by invoking the function encoding_alleles in the class LocusAnnotation. 
	 * The key algorithm of how to design the encoding is implemented in the   
	 *    
	 */
	public void encoding_haps(){
		this.global_haps=new double[this.num_global_hap][this.num_loci];
		// System.out.println("\n" + this.num_global_hap  +"\t" + this.num_loci);
		// System.out.println("\nEncoded:" + this.locusInfo[0].alleles_coding.size());
		for(int h=0;h<this.num_global_hap;h++){
			for(int l=0;l<this.num_loci;l++){
				// if (h == this.num_global_hap - 1) System.out.println(h + "\t" + l + "\t" + this.global_haps_string[h][l] + "\t" + this.locusInfo[l].alleles_coding.get(this.global_haps_string[h][l]));
				try {this.global_haps[h][l]=
						this.locusInfo[l].
						alleles_coding.
						get(this.global_haps_string[h][l]);
				} catch (NullPointerException e) {this.global_haps[h][l]=-1;};
			}
		}
	}
	
	// Write (STDOUT and text file) functions 
	
	/*
	 * Prints intermediate haplotype (single-region) guesses to STDOUT. 
	 */
	public void write_global_stdout() {
		DecimalFormat df = new DecimalFormat("#.####");
		df.setRoundingMode(RoundingMode.CEILING);    
		System.out.println("Hap.\tVar. Comp.\tInter. Freq.");
		for(int h=0;h<this.num_global_hap;h++) {
			System.out.print(h + "\t");
			for(int l=0;l<this.num_loci;l++){
				System.out.print(this.global_haps_string[h][l]);  
			}
			System.out.println("\t" + df.format(this.global_haps_freq[h]));
		}
	}	
	
	/*
	 * Output the in global_haplotypes using string alleles. 
	 */
	public void write_global_file_string(String global_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(global_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+h);
			bw.write("\nFreq");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.global_haps_freq[h]);
			bw.write("\n");
			for(int l=0;l<this.num_loci;l++){
				bw.write(this.locusInfo[l].output2string());
				for(int h=0;h<this.num_global_hap;h++)
					bw.write("\t"+this.global_haps_string[h][l]);
				bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Output the in global_haplotypes using coded alleles. 
	 * 
	 * During the analysis, e.g., divide-and-conquer, when read by the program, these coded alleles 
	 * will become strings and then the next round of coding will be performed. 
	 */
	public void write_global_file_code(String global_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(global_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.hap_IDs[h]);
			bw.write("\nFreq");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+this.global_haps_freq[h]);
			bw.write("\n");
			for(int l=0;l<this.num_loci;l++){
				bw.write(this.locusInfo[l].output2string());
				for(int h=0;h<this.num_global_hap;h++)
					// the only difference between this method and "write_global_file_string" is the line below: 
					bw.write("\t"+this.global_haps[h][l]);  
				bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Output the in pool frequencies. 
	 */
	public void write_inpool(String inpool_hap_output_file, boolean append){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(inpool_hap_output_file, append));
			bw.write("Hap_ID");
			for(int h=0;h<this.num_global_hap;h++)
				bw.write("\t"+h);
			bw.write("\n");
			// System.out.println(this.in_pool_haps_freq.length + "\t" + this.in_pool_haps_freq[0].length);
			for(int p=0;p<this.num_pools;p++){
				bw.write(this.pool_IDs[p]);
				for(int h=0;h<this.num_global_hap;h++)	// TODO Report error! Formerly, this.num_pools.
					if (this.in_pool_haps_freq[h].length == 0) bw.write("\t0");
					else bw.write("\t"+this.in_pool_haps_freq[h][p]);
				bw.write("\n");
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void write_prop_recon(String prop_recon_file, boolean append){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(prop_recon_file, append));
			bw.write("Hap_ID\tProp_Recon\n");
			for(int h=0;h<this.num_global_hap;h++) bw.write(h + "\t" + this.prop_reconstruction[h] + "\n");
			bw.close();
		} catch(Exception e){e.printStackTrace();}
	}
	
	public void write2files(String global_hap_output_file, String in_pool_hap_output_file, String prop_recon_file, String type, boolean append){
		write_inpool(in_pool_hap_output_file, append);
		if(type.equals("string")){
			write_global_file_string(global_hap_output_file, append);
		}else if(type.equals("code")){
			write_global_file_code(global_hap_output_file, append);
		}else{
			System.out.println("Error: the type has to be string or code!");
		}
		write_prop_recon(prop_recon_file, append); 
	}

	/*
	 * Updates the composition of the main global variables when existing haplotypes are removed.
	 * @param Whether to remove the haplotype or not, and the total number to remove.  
	 */
	public void remHaps(boolean[] list_rem_haps, int num_rem_haps) {
		int tot_new_haps = this.num_global_hap - num_rem_haps; 
		// System.out.println(tot_new_haps + "\t" + this.num_global_hap + "\t" + num_rem_haps); 
		String[][] tmp_global_haps_string = new String[tot_new_haps][this.num_loci];
		String[] tmp_hap_IDs = new String[tot_new_haps];
		double tmp_tot_freq = 0; 
		int tmp_count = 0; 
		for(int h=0;h<this.num_global_hap;h++) {
			if (list_rem_haps[h]) {
				// System.out.print(h);
				continue;
			}
			tmp_global_haps_string[tmp_count]=this.global_haps_string[h].clone();
			tmp_hap_IDs[tmp_count]=this.hap_IDs[h];
			tmp_tot_freq += this.global_haps_freq[h];
			tmp_count++; 
		}
		this.global_haps_string = tmp_global_haps_string.clone(); 
		this.hap_IDs = tmp_hap_IDs.clone();
		double[] tmp_global_haps_freq = new double[tot_new_haps];
		int new_hap = 0;
		for(int h=0;h<this.num_global_hap;h++) {
			if (list_rem_haps[h]) continue;
			tmp_global_haps_freq[new_hap] = this.global_haps_freq[h] / tmp_tot_freq;
			new_hap++; 
		}
		this.global_haps_freq = tmp_global_haps_freq.clone(); 
		this.num_global_hap=this.global_haps_freq.length;
		this.construct_hapID2index_map();
		this.encoding_haps(); // initialize this.global_haps;
	}

	/*
	 * Calculates the average variant frequency, var.-covar. matrix for different primitive loci
	 * Calculates the log-likelihood of this HapConfig object explaining the data. 
	 * Updates the global variables mu, sigma, and logL respectively for this HapConfig object. 
	 */
	public void update_sigma_mu_logL(){
		this.mu=new double[this.num_loci];
	    // System.out.println("\nmu:"); 
		for(int l=0;l<this.num_loci;l++){
	    	for(int h=0;h<this.num_global_hap;h++){		    
	    		this.mu[l]=this.mu[l]+this.global_haps[h][l]*this.global_haps_freq[h];
	    	}
    		// System.out.print(this.mu[l] + "\t");
	    } 
	    double[][] eta=new double[this.num_loci][this.num_loci];
	    for(int l1=0;l1<this.num_loci;l1++){
	    	 for(int l2=0;l2<this.num_loci;l2++){
	    		 for(int h=0;h<this.num_global_hap;h++){
	 		    	eta[l1][l2]+=(this.global_haps[h][l1]*this.global_haps[h][l2]*this.global_haps_freq[h]);
	 		     }
			 }
	    }
	    // System.out.println("\n\nsigma:"); 
	    this.sigma=new double[this.num_loci][this.num_loci]; // TODO Update the sigma function to use LDx
	    for(int q1=0;q1<this.num_loci;q1++){
	    	 for(int q2=0;q2<this.num_loci;q2++){
	    		 this.sigma[q1][q2]=eta[q1][q2]-this.mu[q1]*this.mu[q2];
	    		 // System.out.print(this.sigma[q1][q2] + "\t");
	    	 }
	    	 // System.out.println();
	    }
		try {
		    this.logL=Algebra.logL_aems(this.sigma, this.mu, this.inpool_site_freqs);
		} catch (NullPointerException e) {} 
	    // Got rid of Algebra.times(this.sigma, this.est_ind_pool), Algebra.times(this.mu, this.est_ind_pool) because everything is in frequencies.
	}

	/*
	 * Returns haplotypes that are above a certain frequency cutoff.
	 * @param the frequency cutoff.
	 */
	public HapConfig clone(double freq_cutoff){
		if (freq_cutoff == 0) return new HapConfig(this.global_haps_string, this.global_haps_freq, this.in_pool_haps_freq, this.inpool_site_freqs, this.locusInfo, this.num_pools, this.hap_IDs, this.pool_IDs, this.est_ind_pool);
		boolean[] list_rem_haps = new boolean[this.num_global_hap];
		int num_rem_haps = 0; 
		for (int h = 0; h < this.num_global_hap; h++) {
			if (this.global_haps_freq[h] < freq_cutoff) {
				list_rem_haps[h] = true;
				num_rem_haps++;
			}
		}
		// for (int h = 0; h < this.num_global_hap; h++)
			// System.out.println(this.hap_IDs[h]); 
		remHaps(list_rem_haps, num_rem_haps);
		// for (int h = 0; h < this.num_global_hap; h++)
			// System.out.println(this.hap_IDs[h]); 
		return new HapConfig(this.global_haps_string, this.global_haps_freq, this.in_pool_haps_freq, this.inpool_site_freqs, this.locusInfo, this.num_pools, this.hap_IDs, this.pool_IDs, this.est_ind_pool);
	}

}