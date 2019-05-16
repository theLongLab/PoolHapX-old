package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

public class GraphColoring
{
	
	public Vector <Integer> readindex_arr_tmp;     // The index of the genotype when it was added to readinfo_arr_tmp. 
	public Vector <String> readinfo_arr_tmp;       // The list of all possible genotypes, ArrayList<site=allele;...site=allele;>
    public HashMap<String,Integer> output_ref_arr;
    public HashMap<String,String> conf_ref_arr;
    public int num_loci;
    public int num_pools;
	public LocusAnnotation[] locusInfo; 	// #loci. Note that a locus can be a SNP or a region. 
	public double[][] inpool_site_freqs;	// #loci x #pools; Added by Quan Dec. 2018.
	
	// new GraphColoring(gp.inter_dir + prefix + "_p" + p + ".vef", gs_var_pos, gp.inter_dir + prefix + "_p" + p + ".in"); 
	public GraphColoring(String vef, String gs_var_pos, String out_file) throws IOException {
		this.readindex_arr_tmp = new Vector<Integer>();
		this.readinfo_arr_tmp = new Vector<String>();
		int count = 0;
		HashMap<String,Integer>  geno_dict= new HashMap<String,Integer>();    // HashMap<pos=allele;, count>
		BufferedReader bufferedreader= new BufferedReader(new FileReader(vef));
		String line="";
		int max_num_geno=32878;   // The maximum number of times a genotype can be counted in a single pool. TODO Why does this exist?
		while ( (line =bufferedreader.readLine())!=null ){
			line= line.replace("\r","");
			String[] line_arr = line.split("\t");    // Read_Name    Pos=Allele;     //  Start   End
			if (line_arr[1].contains("=")) { // If the read contains a segregating site (i.e.: has a distinguishing genotype)...
				String tmp_geno=line_arr[1];    
				if (!geno_dict.containsKey(tmp_geno)) { // If this combination of alleles hasn't been recorded yet...
				    geno_dict.put(tmp_geno,1);  // ... add it.
					this.readindex_arr_tmp.add(count);  // The index of readinfo_arr_tmp that corresponds to this genotype.
					count=count+1;
					this.readinfo_arr_tmp.add(line_arr[1]);
				}else { // Almost no difference if it has already been recorded yet.
					int tmp_num = geno_dict.get(tmp_geno);
					geno_dict.remove(tmp_geno);
					geno_dict.put(tmp_geno, (tmp_num+1));
					if (geno_dict.get(tmp_geno)< max_num_geno) {             
		                this.readindex_arr_tmp.add(count);
						count=count+1;
						this.readinfo_arr_tmp.add(line_arr[1]);
					}
				}
			}
		}
		bufferedreader.close();
		this.solver(gs_var_pos);
		this.fileOut(out_file);
	}
	
	public GraphColoring(HapConfig[] level_1, HapConfig[] level_2, String gs_var_pos, int dosage) throws IOException {
		this.readindex_arr_tmp = new Vector<Integer>();
		this.readinfo_arr_tmp = new Vector<String>();
		int count = 0; 
        for (int r = 0; r < level_1.length; r++) {
	        for(int h = 0; h < level_1[r].num_global_hap; h++) {
	        	String curr_vc = "";
	            for (int l = 0; l < level_1[r].num_loci; l++) curr_vc += (level_1[r].locusInfo[l].start_loc + "=" + level_1[r].global_haps_string[h][l] + ";");
	            int hap_ct = (int) (level_1[r].global_haps_freq[h] * dosage); 
	            for (int c = 0; c < hap_ct; c++) {
	                this.readindex_arr_tmp.add(count);  // The index of readinfo_arr_tmp that corresponds to this genotype.
	                count++;
	                this.readinfo_arr_tmp.add(curr_vc);
	            }
	            // System.out.println(curr_vc);
	        }
        }
        for (int r = 0; r < level_2.length; r++) {
	        for(int h = 0; h < level_2[r].num_global_hap; h++) {
	        	String curr_vc = "";
	            for (int l = 0; l < level_2[r].num_loci; l++) curr_vc += (level_2[r].locusInfo[l].start_loc + "=" + level_2[r].global_haps_string[h][l] + ";");
	            int hap_ct = (int) (level_2[r].global_haps_freq[h] * dosage); 
	            for (int c = 0; c < hap_ct; c++) {
	                this.readindex_arr_tmp.add(count);  // The index of readinfo_arr_tmp that corresponds to this genotype.
	                count++;
	                this.readinfo_arr_tmp.add(curr_vc);
	            }
	            // System.out.println(curr_vc);
	        }
        }
        // System.out.println("Finished reading in all of the source files.");
        System.out.println("There are " + count + " individual fragments (reads or regional haplotypes) in the dataset.");
		this.solver(gs_var_pos);
	}
	
	public void solver(String gs_var_pos)    // Source_Path or Source_Path
    throws IOException
  {
		HashMap<Integer, Integer> pos_dict = new HashMap<Integer, Integer>(); // HashMap<Seg_Site, Index>
		Integer pos_index = 0; 
	    BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
	    String currLine = br.readLine(); // Skip header. 
	    currLine = br.readLine();
	    while (currLine != null) {
	    	pos_dict.put(Integer.parseInt(currLine.split(";")[1]), pos_index); 
	    	// System.out.println(currLine.split(";")[1] + "\t" +  pos_index);
	    	pos_index++; 
	    	currLine = br.readLine();
	    }
	    br.close();
		this.num_loci = pos_dict.size();
		
		br = new BufferedReader(new FileReader(gs_var_pos));
		currLine = br.readLine();
		currLine = br.readLine();
		int loci_index = 0;
		this.num_pools = currLine.split("\t").length - 1; 
		this.locusInfo = new LocusAnnotation[this.num_loci];
		this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
		while(currLine!=null){
			String[] tmp = currLine.split("\t");
			this.locusInfo[loci_index]=new LocusAnnotation(tmp[0]); // the first column is the locus-info
			for (int p = 0; p < this.num_pools; p++) this.inpool_site_freqs[loci_index][p] = Double.parseDouble(tmp[p + 1]);
			loci_index++;
			currLine = br.readLine();
		} br.close();
		// System.out.println("There are " + this.locusInfo.length + " positions.");
		
		Vector <Integer> readindex_arr = new Vector<Integer>();               // ArrayList<pos=allele;> in random order
		Vector <String> readinfo_arr = new Vector<String>();                // ArrayList<original index of pos=allele;> TODO Confirm!
        
		int[] index_arr_tmp = new int[this.readindex_arr_tmp.size()]; 
		for(int k = 0; k < index_arr_tmp.length; k++) 
		    index_arr_tmp[k] = k ; 
		
		ArrayList<Integer> list = new ArrayList<Integer>();
		for(int i = 0;i < index_arr_tmp.length;i++){
			list.add(index_arr_tmp[i]);
		}

		int[] index_arr = new int[this.readindex_arr_tmp.size()];  // A list of the indices corresponding to each read, in random order.
        // index_arr_tmp, list, and index_arr are copies of this.readindex_arr_tmp so far. 
		Collections.shuffle(list); // It shuffles a given list using the user provided source of randomness.
		// Now, index_arr_tmp, index_arr are copies of this.readindex_arr_tmp. list is now a randomized version of this.readindex_arr_tmp.
        Iterator<Integer> ite = list.iterator();
		int tmp_i=0;
        while(ite.hasNext()){  
//            System.out.println(ite.next().toString()+", ");  
        	index_arr[tmp_i]= ite.next();
        	tmp_i++;
        } 
        // Now, index_arr_tmp is a copy of this.readindex_arr_tmp. list, index_arr are now randomized versions of this.readindex_arr_tmp.
        
//        for(int i = 0;i < index_arr.length;i++){
//        	System.out.println(index_arr[i]);
//        }
        // index_arr is a randomized list of indices corresponding to the genotypes from the VEF file.
		  		
        for(int i = 0;i < this.readindex_arr_tmp.size();i++){
        	readinfo_arr.add(readinfo_arr_tmp.get(index_arr[i]));
        	readindex_arr.add(this.readindex_arr_tmp.get(index_arr[i]));
        }
        // readinfo_arr contains the randomized order list of the genotypes from the VEF file.
        // readindex_arr contains the original indices corresponding to the genotypes from the VEF file i.e.: the order they were read in.
        
//        for(int i = 0;i < readindex_arr.size();i++){
//        	System.out.println(readinfo_arr.get(i));
//        }
        
        
        ArrayList<ArrayList<Integer>> read_pos_2D_arr= new ArrayList<ArrayList<Integer>>(); // ArrayList<Genotype, ArrayList<Index_SS>>
        ArrayList<ArrayList<String>> read_geno_2D_arr= new ArrayList<ArrayList<String>>(); // ArrayList<Genotype, ArrayList<Allele_SS>>
        

//        System.out.println(readinfo_arr);
        for(int i = 0;i < readinfo_arr.size();i++){
        	String tmp_str = readinfo_arr.get(i);
        	String [] tmp_str_arr= tmp_str.split(";");
        	ArrayList<String>  i_geno_arr= new ArrayList<String>();
        	ArrayList<Integer>  i_pos_arr= new ArrayList<Integer>();
        	for (int k =0;k < tmp_str_arr.length;k++) {
        		String [] tmp2_arr = tmp_str_arr [k].split("=");
        		int pos = Integer.parseInt(tmp2_arr[0]);
        		i_pos_arr.add(pos_dict.get(pos));
        		i_geno_arr.add(tmp2_arr[1]);
        	}
        	read_pos_2D_arr.add(i_pos_arr);
        	read_geno_2D_arr.add(i_geno_arr);
        }
        
//        System.out.println(read_pos_2D_arr);
//        System.out.println(read_geno_2D_arr);
       
        
        ArrayList<ArrayList<Integer>> readadj_arr= new ArrayList<ArrayList<Integer>>(); // ArrayList<Genotype, ArrayList<Index_Genotype>>
        for (int i=0; i < read_pos_2D_arr.size();i++) { // For each possible genotype...
        	// if (i % 500 == 0) System.out.println(i + " fragments have been processed.");
        	ArrayList<Integer>  adj_arr = new ArrayList<Integer>();
        	for (int j=0; j < read_pos_2D_arr.size();j++) {    // ...comparing it to all other genotypes...
        		if (i!=j) {
        			boolean IsConnect= false;
        			for (int k=0;k < read_pos_2D_arr.get(i).size();k++) {
        				for (int l=0;l < read_pos_2D_arr.get(j).size();l++) {
        					if ((read_pos_2D_arr.get(i).get(k)== read_pos_2D_arr.get(j).get(l)) &&  
        					(!read_geno_2D_arr.get(i).get(k).equals(read_geno_2D_arr.get(j).get(l))) ) {
        						IsConnect=true;
        					} // If pos(geno_i,index_k) == pos(geno_j,index_l) and the alleles are the same, the two reads can connect.
        				}
        			}
        			if (IsConnect){
        				adj_arr.add(readindex_arr.get(j));  // Add index of geno_j to the list of possible connects for geno_i.
        			}
        		}
        	}
        	readadj_arr.add(adj_arr);
        }
        // System.out.println("Finished identifying all possible conflicts between genotype fragments.");
//        System.out.println(readadj_arr);
        
        int max_color= readindex_arr.size(); // Only if there are all 2^loci genotypes present.
        ArrayList<Integer>  read_color_arr= new ArrayList<Integer>();
        ArrayList<HashSet<Integer>> nb_color_arr  = new ArrayList<HashSet<Integer>>();
//		nb_color_arr.add(new HashSet());
		
        for (int i=0;i< readadj_arr.size();i++) {
        	nb_color_arr.add(new HashSet<Integer>());
        	read_color_arr.add(-1);
        }
        
        read_color_arr.set(0, 0);   // Set the first colour as the genotype at index 0.
        
        for (int i=0;i< readadj_arr.get(0).size();i++) {
        	int index = readadj_arr.get(0).get(i);
        	nb_color_arr.get(index).add(0);    // Add colour 0 as a possible colour to all genotypes that can connect with the genotype index 0.
        }
        
        int real_max_color=0;
        ArrayList<HashSet<String>> color_geno_set_arr = new ArrayList<HashSet<String>>();
        color_geno_set_arr.add(new HashSet<String>());
        color_geno_set_arr.get(0).add(readinfo_arr.get(0)); // This is the full genotype of a colour.
        
//        System.out.println(color_geno_set_arr);
        
        // Make the conflict graph. Basically, assign colours to any genotype fragments that can't go together. 
        while (true) {
        	int max_nb_color=-1;
        	int index=-1;
        	for (int i =0;i < readadj_arr.size();i++) {    // For each genotype...
        		 if  ((read_color_arr.get(i)==-1) && // If there hasn't been a colour assigned to that genotype...
        		 (nb_color_arr.get(i).size()> max_nb_color  ) ) { // ...and there are other genotypes that can go with it of the same colour...
        			 index=i;
        			 max_nb_color= nb_color_arr.get(i).size(); 
        		 }
        	}
        	if (index==-1) {   // If there are no colours left (?) end this step.
        		break;
        	}
        	int color=-1;
        	for (int i=0;i< max_color;i++) {   // For each genotype...
        		if (!nb_color_arr.get(index).contains(i)){    // If the possible colours list for that genotype doesn't contain colour i...
        			if (i< color_geno_set_arr.size()) {  // If colour i is smaller than the number of available colours.
        				if (!color_geno_set_arr.get(i).contains(readinfo_arr.get(index))) { // ...and 
        					color =i;
        					color_geno_set_arr.get(i).add(readinfo_arr.get(index));
        					break;
        				}
        			}else {
        				color_geno_set_arr.add(new HashSet<String>());
        				color_geno_set_arr.get(color_geno_set_arr.size()-1).add(readinfo_arr.get(index));
        				color= i;
        				break;
        			}
        		}
        	}
        	if (color> real_max_color) {
        		real_max_color= color;
        	}
        	read_color_arr.set(index, color);
        	for (int i=0;i< readadj_arr.get(index).size();i++) {
        		nb_color_arr.get(readadj_arr.get(index).get(i)).add(color);
        	}
        	// System.out.println(max_nb_color);
        }
        // System.out.println("Finished identifying potential full-genome genotypes i.e.: haplotypes.");
        // System.out.println(read_color_arr);
        
        String null_ref= "*";
        String conf_ref= "";
        for (int i=1;i<num_loci;i++ ){
        	null_ref= null_ref+"*";
        	conf_ref=conf_ref+"?";
        }
        
        String[] ref_arr = new String[real_max_color+1];
        String[] conf_arr = new String[real_max_color+1];
        for (int i=0;i<= real_max_color;i++) {
        	ref_arr[i] = null_ref;
        	conf_arr[i]= conf_ref;
        }
        
//        String ss= "1234567";
//        System.out.println(ss.substring(1, ss.length()));
        
        for (int i =0;i < read_color_arr.size();i++) {
        	int i_color= read_color_arr.get(i);
        	for (int j=0; j< read_pos_2D_arr.get(i).size();j++) {
//        		System.out.println(read_pos_2D_arr.get(i).get(j).toString());
        		int p = read_pos_2D_arr.get(i).get(j);
        		p++;
        		ref_arr[i_color]= ref_arr[i_color].substring(0, (p-1))+ read_geno_2D_arr.get(i).get(j).toString()
        				+ref_arr[i_color].substring(p, ref_arr[i_color].length());
        	}
        	if (read_pos_2D_arr.get(i).size()>1) {
        		for (int j=0; j< (read_pos_2D_arr.get(i).size()-1);j++) {
        			if (Math.abs(read_pos_2D_arr.get(i).get(j+1)
        					-read_pos_2D_arr.get(i).get(j)) ==1) {
        				int p= read_pos_2D_arr.get(i).get(j);
        				p++;
        				try {
            				conf_arr[i_color]= conf_arr[i_color].substring(0, (p-1))+ "-"
                    				+conf_arr[i_color].substring(p, conf_arr[i_color].length());
        				} catch (StringIndexOutOfBoundsException e){
        					System.out.println(conf_arr[i_color].substring(0, (p-1)) + "\t" + p + "\t" + conf_arr[i_color].length());
        				}
        			}
        		}
        	}
        	
        }
        // for (int i =0;i< conf_arr.length;i++)
        	// System.out.println( ref_arr[i]);
        
        this.output_ref_arr = new HashMap<String,Integer>();
        this.conf_ref_arr = new HashMap<String,String>();
        double completeness_cutoff = 0; 
        for (int i=0; i<= real_max_color;i++) {
        	if (count(ref_arr[i]) <= completeness_cutoff) {	// May implement this in the future.
        		// System.out.println(ref_arr[i]);
        		if (this.output_ref_arr.containsKey(ref_arr[i])) {
        			this.output_ref_arr.put(ref_arr[i],  this.output_ref_arr.get(ref_arr[i])+1);
        			conf_ref_arr.put(ref_arr[i], conf_arr[i]);
        		}else {
        			this.output_ref_arr.put(ref_arr[i], 1);
        			this.conf_ref_arr.put(ref_arr[i], conf_arr[i]);
            	}
        	} else {
        		// System.out.println(ref_arr[i] + "\t" + count(ref_arr[i]));
        	}
        }
        // if (this.conf_ref_arr.isEmpty()) {
        	
        // }
    }
    
	public void fileOut(String out_file) throws IOException {
        FileWriter mydata = new FileWriter(out_file,false);
		PrintWriter pw = new PrintWriter(mydata);
//        System.out.println(output_ref_arr );
        for (String entry : this.output_ref_arr.keySet()) {	// Changed iterator from Map<K,V> -> K because just getting the key requires less memory. 
        	String b= this.conf_ref_arr.get(entry);
        	String c= "";
        	for (int i=0;i< b.length();i++) {
        		c= c+ entry.substring(i,i+1)+b.substring(i, i+1);
        	}
            c= c+ entry.substring(entry.length()-1, entry.length());
            pw.write(c+"\t"+this.output_ref_arr.get(entry).toString()+"\n" );
//            System.out.println(c+"\t"+output_ref_arr.get(x).toString());
        }
        pw.flush();
		pw.close();
//        for (int i =0;i< real_max_color+1;i++) {
//        	System.out.println(ref_arr[i]);
//        }
		return;
	}
	
	public HapConfig hapOut() {
    	int num_global_hap = this.output_ref_arr.size();
    	String[][] global_haps_string = new String[num_global_hap][num_loci];
        int[] global_haps_ct = new int[num_global_hap]; 
        int tot_hap_ct = 0; 
    	int hap_index = 0;
    	for (String entry : this.output_ref_arr.keySet()) {
			String[] var_comp = entry.split(""); 
			for (int v = 0; v < num_loci; v++) global_haps_string[hap_index][v] = var_comp[v]; 
            global_haps_ct[hap_index] = this.output_ref_arr.get(entry); 
            tot_hap_ct += this.output_ref_arr.get(entry); 
			hap_index++;
    	}
    	double[] global_haps_freq = new double[num_global_hap];
    	System.out.println("There are " + num_global_hap + " haplotypes generated by GC.");
        for (int h = 0; h < num_global_hap; h++) global_haps_freq[h] = (double) global_haps_ct[h] / (double) tot_hap_ct;
    	return new HapConfig(global_haps_string, global_haps_freq, null, this.inpool_site_freqs, this.locusInfo, this.num_pools, null, null, 0);
	}
	
	int count(String gc_hap) {
		int unknown = 0; 
		String[] allele_comp = gc_hap.split(""); 
		for (String a : allele_comp) if (a.equals("*")) unknown++;
		return unknown;
	}
}
