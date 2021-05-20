import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

public class CoverageAnalyser {

	// private HashMap<String, ArrayList<Long>> posDepth = new HashMap<String,
	// ArrayList<Long>>();
	private int column = 2;
	private int numCategories = 2;
	
	private double normalDepth = 100.0;
	
	private int maxDepth = 1000;
	private double zetaK = 1.5;
	
	private int MAX_ITERATIONS = 100;

	private long[] count = null;

	private double[] freq = new double[3];

	boolean findNormal = false;

	// private ArrayList<Double> variance = new ArrayList<Double>();

	/*
	 * public void loadDepth(String fn) { try { BufferedReader br = null; if
	 * (fn.equals("-")) br = new BufferedReader(new
	 * InputStreamReader(System.in)); else br = new BufferedReader(new
	 * FileReader(fn));
	 * 
	 * do { ArrayList<String> row = Input.loadTableRow(br, "[\t ]"); if (row ==
	 * null) break; if (row.size() <= column) {
	 * System.err.println("Warning: skipping " + row); continue; }
	 * 
	 * String contig = row.get(0); ArrayList<Long> list = posDepth.get(contig);
	 * if (list == null) { list = new ArrayList<Long>(); posDepth.put(contig,
	 * list); } //only store last position with same depth... int lsize =
	 * list.size(); long cov = 0; if (lsize == 0) { cov =
	 * Long.parseLong(row.get(1)) list.add(cov); // genomic position
	 * list.add(Long.parseLong(row.get(column))); // depth } else { long lastcov
	 * = list.get(lsize - 1); long cov = Long.parseLong(row.get(column)); if
	 * (lastcov == cov) { list.set(lsize - 2, Long.parseLong(row.get(1))); //
	 * update position } else { list.add(Long.parseLong(row.get(1))); // genomic
	 * position list.add(cov); // depth } } } while (true); br.close(); } catch
	 * (Exception e) { e.printStackTrace(); System.err.println("Error in file "
	 * + fn); System.exit(-1); } }
	 */

	public void loadDepth(String fn) {
		count = new long[maxDepth + 2];
		try {
			BufferedReader br = null;
			if (fn.equals("-"))
				br = new BufferedReader(new InputStreamReader(System.in));
			else
				br = new BufferedReader(new FileReader(fn));

			do {
				ArrayList<String> row = Input.loadTableRow(br, "[\t ]");
				if (row == null)
					break;
				if (row.size() <= column) {
					System.err.println("Warning: skipping " + row);
					continue;
				}
				int depth = Integer.parseInt(row.get(column));
				if (depth > maxDepth)
					depth = maxDepth + 1;
				++count[depth];
			} while (true);
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in file " + fn);
			System.exit(-1);
		}
	}

	public static final double isqrtPI = 1.0 / Math.sqrt(Math.PI);
	public static final double isqrt2 = 1.0 / Math.sqrt(2.0);

	private class Gaussian {
		private double mean;
		private double stdev;
		private double istdev;
		private double istdev2sp;
		private double logIstdev2sp;

		// Misc.i2sqrtPI = 1.0 / Math.sqrt(2.0 * Math.PI);

		public Gaussian(double mean, double stdev) {
			istdev = isqrt2 / stdev;
			istdev2sp = istdev * isqrtPI;
			logIstdev2sp = Math.log(istdev2sp);
			this.mean = mean;
			this.stdev = stdev;
		}
		public double getMean()
		{
			return mean;
		}

		public double getStdev()
		{
			return stdev;
		}

		public double density(double x) {
			double tmp = (x - mean) * istdev;
			return istdev2sp * Math.exp(-tmp * tmp);
		}

		public double logDensity(double x) {
			double tmp = (x - mean) * istdev;
			return logIstdev2sp - tmp * tmp;
		}
	}
	
	//distribution zeta - 1  
	private class Zeta {
		//double k;
		int max;
		ArrayList<Double> d = new ArrayList<Double>();
		//ArrayList<Double> ld = new ArrayList<Double>();
		public Zeta(double k, int max) {
			//this.k = k;
			this.max = max;
			double sum = 0.0;
			for (int i = 0; i <= max; ++i) {
				double p = Math.pow(i + 1, -k);
				sum += p;
				d.add(p);
			}
			double iSum = 1.0 / sum;
			for (int i = 0; i <= max; ++i)
				d.set(i, d.get(i) * iSum);
		}
		public double density(int x) {
			if (x > max)
				return 0.0;
			return d.get(x);
		}
		//public double logDensity(int x) {
		//	return ld.get(x);
		//}
	}
	
	private double likelihood(boolean verbose) {
		Gaussian g[] = new Gaussian[numCategories];

		/*double rate1 = count[0] / (double) count[1];
		double rate2 = count[0] / (double) count[2];
		double rate3 = count[0] / (double) count[3];
		double k1 = Math.log(rate1) / Math.log(2.0);
		double k2 = Math.log(rate2) / Math.log(3.0);
		double k3 = Math.log(rate3) / Math.log(4.0);
		k = Math.min(k1, k2, k3);
		if (k > 1.0)
			;
		else k = 1.0;*/
		//double k = 1.2;
		
		Zeta z = new Zeta(zetaK, maxDepth);
		if (verbose)
			System.err.println("Zeta parameter " + zetaK);
		
		{ // Initialisation
			for (int j = 0; j < numCategories; ++j) {
				double mean = (normalDepth * (j + 1)) / (double) numCategories;
				g[j] = new Gaussian(mean, Math.sqrt(normalDepth) * Math.exp(2.0 * Math.random() - 1.0));
			}
			double sum = 0.0;
			for (int j = 0; j <= numCategories; ++j) {
				double r = Math.random(); 
				freq[j] = r;
				sum += r;
			}
			double iSum = 1.0 / sum;
			for (int j = 0; j <= numCategories; ++j)
				freq[j] *= iSum;
		}

		double d[] = new double[numCategories + 1];

		double likelihood = 0.0;
		for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
			double qFreq[] = new double[numCategories + 1];
			double qVar[] = new double[numCategories];
			double qSum[] = new double[numCategories];

			//System.err.println(freq[0] + "\t" + freq[1] + "\t" + freq[2]);
			likelihood = 0.0;
			for (int i = 0; i <= maxDepth; ++i) {

				double sum = 0.0;
				for (int j = 0; j <= numCategories; ++j) {
					double p = 0.0;
					if (j < numCategories)
						p = g[j].density(i) * freq[j];
					else 
						p = z.density(i) * freq[j];
					d[j] = p;
					sum += p;
				}
				//System.err.println("Sum " + sum);
				double iSum = 1.0 / sum;
				
				for (int j = 0; j <= numCategories; ++j) {
					double q = d[j] * iSum * count[i];
					qFreq[j] += q;
					if (j < numCategories) {
						double diff = (i - g[j].mean);
						qVar[j] += q * diff * diff;
						qSum[j] += q;
					}
				}
				likelihood += Math.log(sum) * count[i];
			}
			double sum = 0.0;
			for (int j = 0; j <= numCategories; ++j)
				sum += qFreq[j];
			//for (int j = 0; j <= numCategories; ++j)
			//	if (qFreq[j] < 0.1 * sum)
			//		qFreq[j] = 0.1 * sum;
			//sum = 0.0;
			//for (int j = 0; j <= numCategories; ++j)
			//	sum += qFreq[j];
			
			double iSum = 1.0 / sum;
			for (int j = 0; j <= numCategories; ++j)
				freq[j] = qFreq[j] * iSum;

			//variance in gaussians must be non-decreasing 
			for (int j = numCategories - 1; j >= 1; --j) {
				if (qVar[j] * qSum[j - 1] < qVar[j - 1] * qSum[j]) {
					double qs = qSum[j] + qSum[j - 1];
					double qv = qVar[j] + qVar[j - 1];
					qSum[j] = qs;
					qSum[j - 1] = qs;
					qVar[j] = qv;
					qVar[j - 1] = qv;
				}
			}
			for (int j = 0; j < numCategories; ++j) {
				double sdev = Math.sqrt(qVar[j] / qSum[j]);
				if (sdev < 1.0)
					sdev = 1.0;
				//System.err.println(freq[j] + "\t" + sdev);
				g[j] = new Gaussian((normalDepth * (j + 1)) / (double) numCategories, sdev);
				//System.err.println("mean=" + g[j].mean);
			}
			//System.err.println(freq[numCategories]);
			// System.err.println(Math.sqrt(iSum * qVar[0]) + "\t" +
			// Math.sqrt(iSum * qVar[1]) + "\t" + Math.sqrt(iSum * qVar[2]));
			if (verbose)
				System.err.println("logL=" + likelihood);
		}
		long totalCount = 0;
		for (int i = 0; i <= maxDepth; ++i) {
			totalCount += count[i];
		}
		if (verbose)
			System.err.println(freq[0] + "\t" + freq[1] + "\t" + g[0].getStdev() + "\t" + g[1].getStdev());

		for (int i = 0; i <= maxDepth; ++i) {
			
			double tmp[] = new double[numCategories + 2];
			for (int j = 0; j < numCategories; ++j)
				tmp[j] = g[j].density(i);

			//split zeta distribution to 0 (<=normalDepht/numCategories) and repeat (>=normalDepht)
			double zf = z.density(i);
			if (i <= g[0].mean)
				tmp[numCategories] = zf;
			else if (i >= g[numCategories - 1].mean)
				tmp[numCategories + 1] = zf;
			else {
				double t = (i - g[0].mean) / (g[numCategories - 1].mean - g[0].mean);
				tmp[numCategories] = zf * (1.0 - t);
				tmp[numCategories + 1] = zf * t;
			}

			double sum = 0.0;
			double max = 0.0;
			for (int j = 0; j <= numCategories + 1; ++j) {
				double p = tmp[j] * freq[Math.min(j, numCategories)];
				sum += p;
				max = Math.max(max, tmp[j]);
			}
			if (verbose) {
				System.out.print(i + "\t" + count[i] + "\t" + (totalCount * sum) );
			
				double iMax = 1.0 / max;
				System.out.print("\t" + tmp[numCategories] * iMax);
				for (int j = 0; j <= numCategories + 1; ++j)
					if (j != numCategories)
						System.out.print("\t" + tmp[j] * iMax);
				
				System.out.println();
			}
		}

		
		//System.err.println(num[0] + "\t" + num[1] + "\t" + num[2] + "\t" + num[3]);
		return likelihood;
	}

	
	public void analyse() {
		if (findNormal) {
			// Find normalDepth and other parameters...
			System.err.println("Estimating normal depth parameter...");
			double maxL = Double.NEGATIVE_INFINITY;
			int maxNd = 1;
			for (int nd = 1; nd + nd <= maxDepth; ++nd) {
				normalDepth = nd;
				double ll = likelihood(false);
				if (ll > maxL) {
					System.err.println(nd + "\t" + ll + "\t*");
					maxL = ll;
					maxNd = nd;
				} else
					System.err.println(nd + "\t" + ll);
			}
			normalDepth = maxNd;
			System.err.println("Value " + maxNd + " chosen for normal depth");
		}
		likelihood(true);
	}

	public void setNumCategories(int value) {
		numCategories = value;
		freq = new double[numCategories + 1];


		for (int i = 0; i <= numCategories; ++i)
			freq[i] = Math.random();
	}

	public void setNormalDepth(double value) {
		normalDepth = value;
		maxDepth = (int)(3 * value + 0.99);
		setNumCategories(numCategories);
	}

	public void setColumn(int value) {
		column = value;
	}

	public void setMaxDepth(int value) {
		maxDepth = value;
		findNormal = true;
	}

	public void setZeta(double value) {
		zetaK = value;
	}

	public static void usageInfo() {
//		pp.warning(new String[] { "depth", "numCategories", "normalDepth",
//				"column", "sample", "maxDepth", "zeta"});
		System.err.println("Usage: java CoverageAnalyser depth=depth.txt [parameters] ");
		System.err.println("       depth=file           output from samtools -a depth");
		System.err.println("       numCategories=NUM    the ploidy [2]");
		System.err.println("       normalDepth=NUM      the normal depth [not set, ML value searched]");
		System.err.println("       maxDepthDepth=NUM    search normalDepth between 1..NUM, applies only when normalDepth is not set [1000]");
        System.err.println("       column=NUM           Take column NUM from depth file [3]");
        System.err.println("       zeta=NUM             Zeta distribution parameter [1.5]");
		
		System.exit(0);
	}

	public static void main(String args[]) {
		ParameterParser pp = new ParameterParser();

		String extraParameters = "";
		for (int i = 0; i < args.length; ++i) {
			extraParameters += " " + args[i];
		}
		if (args.length == 0 || !pp.init(extraParameters)) {
			usageInfo();
		}
		System.out.println("#java CoverageAnalyser " + extraParameters);
		
		pp.warning(new String[] { "depth", "numCategories", "normalDepth",
				"column", "maxDepth", "zeta"});

		String depthFile = pp.getValueAsString("depth", null);
		if (depthFile == null)
			usageInfo();

		CoverageAnalyser ca = new CoverageAnalyser();
		int column = Integer.parseInt(pp.getValueAsString("column", "3"));
		ca.setColumn(column - 1);
		ca.setNumCategories(Integer.parseInt(pp.getValueAsString(
				"numCategories", "2")));

		String normalD = pp.getValueAsString("normalDepth", null);
		if (normalD != null) {
			ca.setNormalDepth(Double.parseDouble(normalD));
			//ca.setMaxDepth();
			String md = pp.getValueAsString("maxDepth", null);
			if (md != null)
				ca.setMaxDepth(Integer.parseInt(md));
		}
		else
			ca.setMaxDepth(Integer.parseInt(pp.getValueAsString("maxDepth",
					"1000")));

		ca.setZeta(Double.parseDouble(pp.getValueAsString("zeta", "1.5")));
		
		ca.loadDepth(depthFile);
		ca.analyse();

	}
}
