import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;

/************************ Reference ************************
 Q. Zhang, W. Fang, J. Sun, and Q. Wang, “Improved genetic algorithm
 for high-utility itemset mining,” IEEE Access, vol. 7, pp. 176 799–
 176 813, 2019.
 Experiments ran on a 64-bit Windows 10 system with an Intel Core i5
 CPU with 8 GB of RAM, 250SSD.
 **********************************************************
 Part of the code refers to SPMF:
 Fournier-Viger, P., Lin, C.W., Gomariz, A., Gueniche, T., Soltani,
 A., Deng, Z., Lam, H. T. (2016). The SPMF Open-Source Data Mining
 Library Version 2. Proc.  19th European Conference on Principles of
 Data Mining and Knowledge Discovery (PKDD 2016) Part III, Springer
 LNCS 9853, Â pp. 36-40.
 **********************************************************/



public class HUIM_Main {

	public static void main(String[] args) throws IOException{

		String input = fileToPath("contextHUIM.txt");
		int min_utility = 40;
		String output = "result.txt";

		Algorithm algorithm = new Algorithm();
		algorithm.runAlgorithm(input, output, min_utility);
		algorithm.printStats();
	}
	
	public static String fileToPath(String filename) throws UnsupportedEncodingException{
		URL url = Algorithm.class.getResource(filename);
		 return java.net.URLDecoder.decode(url.getPath(),"UTF-8");
	}

}
