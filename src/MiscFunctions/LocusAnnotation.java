package MiscFunctions;

import java.util.HashMap;

public class LocusAnnotation {
	
	/* 
	 * @author Quan Long. Oct 12, 2018
	 * 
	 * Information and Annotation regarding an individual locus 
	 * 
	 * Multiple methods for encoding alleles into numbers will be implemented in this class,
	 * reflecting different understanding of the relationships between alleles.
	 */
	
	boolean is_region;	// whether this is a primitive locus (e.g., a SNP or indel) or a region that contains multiple primitive loci
	int chr_index; 	// index of the chromosome, starting from ZERO
	int start_loc; 
	int end_loc;	// For SNPs and indels, end_loc=start_loc; start_loc!=end_loc only if it is a region, instead of a primitive locus
	String[] alleles; // Note that the length of each allele are the same if this.is_region==false. 
	public HashMap<String, Double> alleles_coding;
	String annotation=null;	
	
	/*
	 *  constructor using multiple variables 
	 */
	public LocusAnnotation(boolean is_region, int chr_index, int start_loc, int end_loc, String[] alleles){
		this.is_region=is_region;
		this.chr_index=chr_index;
		this.start_loc=start_loc;
		this.end_loc =end_loc;
		this.alleles=alleles.clone();
		this.encode_alleles();
	}
	
	/*
	 * constructor using a string, formated as:
	 * chr_index;start_loc;end_loc;alleles (the alleles are separated by ":")
	 * chr_index starts with zero. 
	 * 
	 * NOTE: if it is an indel, start_loc=end_loc; Here start_loc!=end_loc only if it is a region, instead of a primitive locus.
	 * 
	 */
	public LocusAnnotation(String locus_info){
		String[] the_info=locus_info.split(";");
		this.chr_index=Integer.parseInt(the_info[0]);
		this.start_loc=Integer.parseInt(the_info[1]);;
		this.end_loc =Integer.parseInt(the_info[2]);;
		this.alleles=the_info[3].split(":");
		this.is_region=(this.start_loc!=this.end_loc);
		this.encode_alleles();
	}
	
	public String output2string(){
		String allels_str=this.alleles[0];
		for(int i=1;i<this.alleles.length;i++)
			allels_str=allels_str+":"+this.alleles[i];
		String out_str=this.chr_index+";"+this.start_loc+";"+this.end_loc+";"+allels_str;
		return out_str;
	}
	/*
	 * constructor with annotation string initialized.
	 */
	public LocusAnnotation(boolean is_region, int chr_index, int start_loc, int end_loc, String[] alleles, String annotation){
		this(is_region, chr_index, start_loc, end_loc, alleles);
		this.annotation=annotation;
	}
	
	/*
	 * constructor with annotation string initialized.
	 */
	public LocusAnnotation(String locus_info, String annotation){
		this(locus_info);
		this.annotation=annotation;
	}
	
	/*
	 * Map string-coded alleles to integers  
	 * 
	 * Please note that this map is not trivial when there are many alleles at the same locus. 
	 * This function serves as a core of how we understand the relationship between alleles.
	 * It impacts the calculation of the co-variance matrix between alleles (in AEM)  
	 * It may be rewritten for different applications.
	 */
	
	void encode_alleles(){
		// if it is not a region, then just code the alleles using 0,1,2,3,... in a random order.
		if(!is_region){
			this.alleles_coding=new HashMap<String, Double>();
			for(int i=0;i<this.alleles.length;i++)
				this.alleles_coding.put(this.alleles[i], (double) i);
		}else{  // it is a region. We assign the code based on the similarity of the alleles.
			// form a distance matrix 
			int allele_num=this.alleles.length;
			// first, all of the alleles should have the same string length (covering the same number of primitive loci in the region)! 
			int allele_string_length=this.alleles[0].length(); 
			for(int k=1;k<allele_num;k++){
				if(this.alleles[k].length()!=allele_string_length){
					System.out.println("ERROR: Incorrect allele_string_length!");
					return;
				}
			}
			double[][] the_distances=new double[allele_num][allele_num];
			for(int i=0;i<allele_num;i++){
				for(int j=0;j<allele_num;j++){
					for(int k=0;k<allele_string_length;k++){
						if(this.alleles[i].charAt(k)!=this.alleles[j].charAt(k)){
							the_distances[i][j]++;	// Raw sum of primitive locus differences between sub-haplotypes i and j. 
						}
					}
				}
			}			
			// project the distance matrix to a direction for single values. //TODO
			// double[] hap_codes= matrix_projection(the_distances); // TODO
		}		
	}
	
	/*
	 *  TODO: project the distance matrix to a vector.
	 */
	public static double[] matrix_projection(double[][] distance_matrix){
		return null;
	}
}
