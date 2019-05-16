package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class PairedReadLinker {
	
	public static void main(String[] args) throws IOException{	
		String prefix = args[0];
		HashMap<String,Integer> dict_read = new HashMap<String,Integer>();
		ArrayList<String>  readname_arr= new ArrayList<String >();
		ArrayList<String>  pos_arr= new ArrayList<String >();
		ArrayList<Integer>  readstart_arr= new ArrayList<Integer >();
		ArrayList<Integer>  readend_arr= new ArrayList<Integer >();
		ArrayList<String>  readrange_arr= new ArrayList<String >();
		BufferedReader bufferedreader= new BufferedReader(new FileReader(prefix + ".raw.vef"));
		String line= null;
		int index=0;
		
		while ( (line =bufferedreader.readLine())!=null ){
			line= line.replace("\r","");
			String[] line_arr = line.split("\t");
			String readname = line_arr[0];
			int readstart = Integer.parseInt(line_arr[3]);
			if (dict_read.containsKey(readname)) {
				String[] vars = line_arr[1].split(";|=");
				String to_add = ""; 
				for (int v = 0; v < vars.length; v += 2)
					if (!pos_arr.get(dict_read.get(readname)).contains(vars[v])) to_add = to_add + vars[v] + "=" + vars[v + 1] + ";"; 
				if 	(readstart>=readstart_arr.get(dict_read.get(readname))){
					readrange_arr.set(dict_read.get(readname), readrange_arr.get(dict_read.get(readname))+"\t"+line_arr[3]
							+ "\t"+line_arr[4]); 
					pos_arr.set(dict_read.get(readname), pos_arr.get(dict_read.get(readname))+to_add);
				}else {
					readrange_arr.set(dict_read.get(readname),line_arr[3]+ "\t"+line_arr[4]+"\t" 
									+readrange_arr.get(dict_read.get(readname))); 
					pos_arr.set(dict_read.get(readname), to_add+pos_arr.get(dict_read.get(readname)));		
				}
				
			}else {
				dict_read.put(readname, index);
				index ++;
				readname_arr.add(line_arr[0]);
				pos_arr.add(line_arr[1]);
				readstart_arr.add(Integer.parseInt(line_arr[3]));
				readend_arr.add(Integer.parseInt(line_arr[4]));
				readrange_arr.add(line_arr[3]+"\t"+line_arr[4]);
			}
		}
		bufferedreader.close();
		
		PrintWriter pw = new PrintWriter(new FileWriter(prefix + ".vef",false));		
		for (int i =0;i<readname_arr.size();i++ ) {
			String tmp_str = readname_arr.get(i)+"\t"+pos_arr.get(i)+"\t"+"//"+"\t"+readrange_arr.get(i)+"\n";
			pw.write(tmp_str);
//			System.out.println(tmp_str);
		}
        pw.flush();
		pw.close();
	}

}
