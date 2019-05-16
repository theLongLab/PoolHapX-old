package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;

public class SiteInPoolFreqAnno {
	
	public int num_sites;
	public int num_pools;
	public String[] pool_IDs;
	public LocusAnnotation[] loci_annotations;  
	public double[][] inpool_freqs; 			//[this.num_sites][this.num_pools]
	
	/*
	 * Constructor: file format below: (first line: header; the rest are variant annotations (Column #1) and frequencies.)
	 *  
	 * 	Var_ID  0       1       2       3       4     
	 *	0;329;329;0:1   0.8054936896807721      0.8954512105649303      1.0     0.7920611798980335      0.807418
	 *	0;796;796;0:1   0.1968  0.10485592315901815     0.0     0.20496397117694154     0.1938584779706275      
	 *	0;833;833;0:1   0.20240320427236316     0.2144571885836223      0.0952253934382502      0.09605122732123
	 */
	public SiteInPoolFreqAnno(String input_file){
		try{
			// read the file first and figure out the length of related fields. 
			 BufferedReader br=new BufferedReader(new FileReader(input_file));
			 String line=br.readLine();
			 while(line.startsWith("#")){ // skip headers, if any
				 line=br.readLine();
			 }
			 String[] header=line.split("\t");
			 this.num_pools=header.length-1;
			 line=br.readLine();
			 while(line!=null){
				 this.num_sites++;
				 line=br.readLine();
			 }br.close();
			 // create the arrays based on the length. 
			 this.inpool_freqs=new double[this.num_sites][this.num_pools];
			 this.pool_IDs=new String[this.num_pools];
			 for(int p=0;p<this.num_pools;p++){
				 this.pool_IDs[p]=header[p+1];
			 }
			 this.loci_annotations=new LocusAnnotation[this.num_sites];
			 // read the file again to fill the fields.
			 br=new BufferedReader(new FileReader(input_file));
			 line=br.readLine();
			 while(line.startsWith("#")){ // skip headers, if any
				 line=br.readLine();
			 }
			 line=br.readLine(); // skip the header
			 int locus_index=0;
			 while(line!=null){
				 String[] tmp=line.split("\t");
				 this.loci_annotations[locus_index]=new LocusAnnotation(tmp[0]);
				 // System.out.println("ding\t" + this.loci_annotations[0].alleles_coding.size());
				 for(int p=0;p<this.num_pools;p++){
					 this.inpool_freqs[locus_index][p]=Double.parseDouble(tmp[p+1]);
				 }
				 locus_index++;
				 line=br.readLine();
			 }
		 }catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Constructor: copy by reference -- no clone! 
	 */
	public SiteInPoolFreqAnno(int num_sites, int num_pools, String[] pool_IDs, LocusAnnotation[] loci_annotations, double[][] inpool_freqs){
		this.num_sites=num_sites;
		this.num_pools=num_pools;
		this.pool_IDs=pool_IDs;
		this.loci_annotations=loci_annotations;  
		this.inpool_freqs=inpool_freqs;
	}
	
	/*
	 * extract the subset of loci and return a new object.
	 * To save memory, here the annotations and frequencies are copied by reference; 
	 */
	public SiteInPoolFreqAnno subset(int start_index, int end_index){
		int sub_num_sites=end_index-start_index+1;
		LocusAnnotation[] sub_loci_annotations=new LocusAnnotation[sub_num_sites];
		double[][] sub_inpool_freq=new double[sub_num_sites][];
		for(int i=start_index;i<=end_index;i++){
			sub_loci_annotations[i-start_index]=this.loci_annotations[i]; 
			sub_inpool_freq[i-start_index]=this.inpool_freqs[i];
		}
		return new SiteInPoolFreqAnno(sub_num_sites, this.num_pools, this.pool_IDs, sub_loci_annotations, sub_inpool_freq);
		
	}
}
