package MiscFunctions; 

import java.util.ArrayList;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class Algebra {
	public static double[] minus(double[] a1, double[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out minus: a1.length!=a2.length");
			return null;
		}
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]-a2[k];
		}return result;
	}
	
	public static double[] minus(double[] a1, int[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out minus: a1.length!=a2.length");
			return null;
		}
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]-a2[k];
		}return result;
	}
	
	public static double[] minus(int[] a1, double[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out minus: a1.length!=a2.length");
			return null;
		}
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]-a2[k];
		}return result;
	}
	
	public static String intarray2string(int[] array){
		String the_str="";
		for(int k=0;k<array.length;k++)the_str=the_str+array[k];
		return the_str;
	}
	
	public static double sum(double[] a){
		double result=0;
		for(int k=0;k<a.length;k++)result+=a[k];
		return result;
	}
	
	public static double sum(ArrayList<Double> a){
		double result=0;
		for(int k=0;k<a.size();k++)result+=a.get(k);
		return result;
	}
	
	public static double[] add(double[] a1, double[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out add: a1.length!=a2.length");
			return null;
		}
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]+a2[k];
		}return result;
	}
	
	public static double[] add(double[] a1, double c){
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]+c;
		}return result;
	}
	
	public static double[] times(double[] a1, double[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out inner mutiplication: a1.length!=a2.length");
			return null;
		}
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]*a2[k];
		}return result;
	}
	
	public static double[] times(double[] a1, double c){
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]*c;
		}return result;
	}
	
	public static double[][] times(double[][] a1, double c){
		double[][] result=new double[a1.length][];
		for(int k=0;k<a1.length;k++){
			result[k]=times(a1[k], c);
		}return result;
	}
	
	public static double[] divide(double[] a1, double[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out inner divide: a1.length!=a2.length");
			return null;
		}
		double[] result=new double[a1.length];
		for(int k=0;k<a1.length;k++){
			result[k]=a1[k]/a2[k];
		}return result;
	}
	
	public static double inner_product(double[] a1, double[] a2){
		if(a1.length!=a2.length){
			System.out.println("Can't carry out inner product: a1.length!=a2.length");
			return Double.NaN;
		}double result=0;
		for(int k=0;k<a1.length;k++){
			result+=a1[k]*a2[k];
		}return result;
	}
	
	public static double mean(double[] p){
		return sum(p)/p.length;
		
	}
	public static void rmlow_and_normalize(double[] p, double cutoff){
		for(int k=0;k<p.length;k++){
			if(p[k]<cutoff)p[k]=0;
		}
		double sum=sum(p);
		for(int k=0;k<p.length;k++)p[k]=p[k]/sum;
	}
	
	/*
	 *  Given an array, divide all the entries by their sum to form a distribution 
	 */
	public static void normalize_ditribution(double[] p){
		double sum=sum(p);
		for(int k=0;k<p.length;k++)p[k]=p[k]/sum;
	}
	
	public static void normalize_ditribution(ArrayList<Double> p){
		double sum=sum(p);
		for(int k=0;k<p.size();k++){
			double x=p.get(k)/sum;
			p.set(k, x);
		}
	}
	
	/*
	 *  assigning NaN with an average. the overal sum of all NaN entries will the the "proportion" of others 
	 */
	public static void normalize_ditribution(ArrayList<Double> p, double proportion){
		double sum0=0;
		int total_NaN=0;
		for(int k=0;k<p.size();k++){
			if(!Double.isNaN(p.get(k)))sum0+=p.get(k);
			else total_NaN++;
		}
		if(sum0==0){ // all are NaN
			for(int k=0;k<p.size();k++)p.set(k, 1.0/p.size());
			return; 
		}
		double sum=sum0/(1-proportion);
		double nan_average=sum*proportion/total_NaN;
		for(int k=0;k<p.size();k++){
			if(!Double.isNaN(p.get(k))){
				double x=p.get(k)/sum;
				p.set(k, x);
			}else p.set(k, nan_average);
		}
	}
	
	public static String intarray2string(int[] a, String sep){
		String st="";
		for(int k=0;k<a.length-1;k++)st=st+a[k]+sep;
		st=st+a[a.length-1];
		return st;
	}
	
	public static double[] abs(double[] p){
		double[] abs=new double[p.length];
		for(int k=0;k<p.length;k++)abs[k]=Math.abs(p[k]);
		return abs;
	}
	
	public static double[] exp(double[] a){
		double[] exp=new double[a.length];
		for(int k=0;k<a.length;k++)exp[k]=Math.exp(a[k]);
		return exp;
	}
	
	public static double quadratic_form(double[] X, double[][] A){
		RealMatrix AA= MatrixUtils.createRealMatrix(A);
		double[] tXA=AA.preMultiply(X);
		return inner_product(tXA,X);
	}
	
	public static double[] diag(double[][] A){
		double[] diag=new double[A.length];
		for(int k=0;k<A.length;k++)diag[k]=A[k][k];
		return diag;
	}
	
	/*
	 *  density of normal distribution as likelihood 
	 */
	public static double logL_normal(double[][] sigma, double[] mu, double[] x){
		int n=mu.length;
		if(x.length!=n || sigma.length!=n || sigma[0].length!=n)System.out.println(x);
		double log_likelihood=Math.log(2*Math.PI)*(-n/2.0);
		SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(sigma));
	    double[] sv=svd.getSingularValues();
	    double log_p_det=0;
	    for(int k=0;k<sv.length;k++){
	    	if(sv[k]!=0)log_p_det=log_p_det+Math.log(sv[k]);
	    }
	    log_likelihood=log_likelihood-0.5*log_p_det;
	    RealMatrix si= svd.getSolver().getInverse();
	    double quadratic=quadratic_form(minus(x, mu), si.getData());
	    log_likelihood=log_likelihood + (-0.5*quadratic);
		return log_likelihood;
	}
	
	public static double distance(double[] x, double[] y){
		double dis=0;
		for(int k=0;k<x.length;k++)dis+=Math.abs(x[k]-y[k]);
		return dis;
	}
	
	public static boolean not_in_span(SingularValueDecomposition svd, double[][] sigma, double[] mu, double[][] x){
		double minimal=0.001;
		DecompositionSolver solver=svd.getSolver();
		for(int k=0;k<x.length;k++){
			double[] B=minus(mu, x[k]);
			double[] X=solver.solve(new ArrayRealVector(B)).toArray();	
			double[] Bprime=(new Array2DRowRealMatrix(sigma)).operate(X);
			double dis=distance(B, Bprime);
			if(dis>minimal)return true;
		}return false;	
	}
	
	/*
	 *  product of multiple density sharing the same normal distribution (\sigma and \mu) 
	 *  To save time of solving the matrix we can't run the above method k times. 
	 */
	public static double logL_aems(double[][] sigma, double[] mu, double[][] x){		// Formerly, logL_normal.
		// x is public double[][] inpool_site_freqs i.e.: #loci x #pools, so it should be x.length
		int n=mu.length;
		int dim_span=0;
		// System.out.println(x.length + "\t" + sigma.length + "\t" + sigma[0].length);
		if(x.length!=n || sigma.length!=n || sigma[0].length!=n) System.out.println("logL_aems: Length of arrays don't match!");	// If the number of variant positions don't match up between all three input parameters, report the error.
		SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(sigma));
	    double[] sv=svd.getSingularValues();	// Solve the covariate matrix using singular value decomposition to get its eigenvalues.
	    // if(not_in_span(svd, sigma, mu, x))return Double.NaN;	// (0-i) If the the haplotypes are not in the span return an NaN. 
	    double log_p_det=0;
	    for(int k=0;k<sv.length;k++){
	    	if(sv[k]!=0){	// If the singular value/eigenvalue at variant k is not 0, add log(SV) to log_p_det (starts at 0). Also, increment dim_span. 
	    		log_p_det=log_p_det+Math.log(sv[k]);
	    		dim_span++;
	    	}
	    }//System.out.println("dim_span: "+dim_span);
	    double constant=Math.log(2*Math.PI)*(-dim_span/2.0);	// Instead of using the number of variants (mu.length), Set the base value of the log-likelihood (LL) as...
	    constant=constant-0.5*log_p_det;	// ... log(2 * pi) / (1/2 * size_min_set_span), and subtract half of log_p_det from the log-likelihood.
	    									// Recall that log_p_det is the sum of the log(covariate eigenvalues).
	    RealMatrix si= svd.getSolver().getInverse();	// Get the inverse matrix of the covariate eigenvalues.
		int num_x=x[0].length; 
		double logL=0;
		for(int k=0;k<num_x;k++){	// For each pool in the double[][] data matrix...
			double[] data_inpool = new double[n]; 
			for (int v = 0; v < n; v++) data_inpool[v] = x[v][k];
			double quadratic=quadratic_form(minus(data_inpool, mu), si.getData());	// (0-ii) Get the inner product of (data - mu) and (inverse matrix of covariate eigenvalues). Add that base value to the running total of the log-likelihood.
		    double log_likelihood=constant + (-0.5*quadratic);	// ...subtract half of the quadratic form from the base value of the log-likelihood. 
		    logL=logL+log_likelihood;	// In other words, the base log-likelihood decreases in a different way for each pool, and total log-likelihood is a sum of all of the individual ones.
		}	
		//System.out.println("debug: "+logL);		
		return logL;
	}

	public static double logL_rjmcmc(double[][] sigma, double[] mu, double[] x){	// NOTE: This might be wayyyyyy too small to use. Will have to observe the values coming out of it.
		int n=mu.length;
		int dim_span=0;
		if(x.length!=n || sigma.length!=n || sigma[0].length!=n)	// Note that x is the DATA i.e.: the observed allele frequencies in pool p, as set by update_sigma_mu_logL.
			System.out.println("logL_rjmcmc: Length of arrays don't match!");	// If the number of variant positions don't match up between all three input parameters, report the error.
		SingularValueDecomposition svd=new SingularValueDecomposition(MatrixUtils.createRealMatrix(sigma));
	    double[] sv=svd.getSingularValues();	// Solve the covariate matrix using singular value decomposition to get its eigenvalues.
	    // if(not_in_span(svd, sigma, mu, x))return Double.NaN;	// (0-i) If the the haplotypes are not in the span return an NaN. 
	    double log_p_det=0;
	    for(int k=0;k<sv.length;k++){
	    	if(sv[k]!=0){	// If the singular value/eigenvalue at variant k is not 0, add log(SV) to log_p_det (starts at 0). Also, increment dim_span. 
	    		log_p_det=log_p_det+Math.log(sv[k]);
	    		dim_span++;
	    	}
	    }//System.out.println("dim_span: "+dim_span);
	    double constant=Math.log(2*Math.PI)*(-dim_span/2.0);	// Instead of using the number of variants (mu.length), Set the base value of the log-likelihood (LL) as...
	    constant=constant-0.5*log_p_det;	// ... log(2 * pi) / (1/2 * size_min_set_span), and subtract half of log_p_det from the log-likelihood.
	    									// Recall that log_p_det is the sum of the log(covariate eigenvalues).
	    RealMatrix si= svd.getSolver().getInverse();	// Get the inverse matrix of the covariate eigenvalues.
		double logL=0;
		double quadratic=quadratic_form(minus(x, mu), si.getData());	// (0-ii) Get the inner product of (data - mu) and (inverse matrix of covariate eigenvalues). Add that base value to the running total of the log-likelihood.
	    double log_likelihood=constant + (-0.5*quadratic);	// ...subtract half of the quadratic form from the base value of the log-likelihood. 
	    logL=logL+log_likelihood;	// In other words, the base log-likelihood decreases in a different way for each pool, and total log-likelihood is a sum of all of the individual ones.
		//System.out.println("debug: "+logL);		
		return logL;
	}
}
