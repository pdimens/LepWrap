
public class LiftoverHaplotypes {
	private static void usageInfo()
	{
        System.err.println("usage: java LiftoverHaplotypes haplotypes=haplo.bed map=mapX.txt chain=all.chain >mapX_liftover.txt");
        System.err.println("       haplotypes=FILE      try liftover for markers in removed haplotypes listed in FILE (from findFullHaplotypes or from PlaceAndOrientContigs)");
        System.err.println("       map=FILE2            any file with contig and pos in the first columns, typically the map file");
        System.err.println("       chain=FILE3          the chain file");
	}
	
	public static void main(String[] args)
	{
		
        if (args.length == 0) {
        	usageInfo();
            System.exit(0);
        }
	    String extraParameters = "";
	    for (int i = 0; i < args.length; ++i) {
	            extraParameters +=  " " + args[i];
	    }
		ParameterParser pp = new ParameterParser();
		if (!pp.init(extraParameters)) {
        	usageInfo();
            System.exit(0);
		}
		pp.warning(new String[]{"map", "chain", "haplotypes"});
		
		String haplotypes = pp.getValueAsString("haplotypes", null);
		String chain = pp.getValueAsString("chain", null);
		String map = pp.getValueAsString("map", null);
		if (map == null || chain == null || haplotypes == null){
        	usageInfo();
            System.exit(0);
		}
		
		
		System.out.println("#java LiftoverHaplotypes" + extraParameters);
		
		PlaceAndOrientContigs poc = new PlaceAndOrientContigs();
		poc.liftover(haplotypes, chain, map);
	}

}
