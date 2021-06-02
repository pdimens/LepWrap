/**
    This file is part of Lep-MAP3.

    Lep-MAP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Lep-MAP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Lep-MAP.  If not, see <http://www.gnu.org/licenses/>.

	Copyright (C) 2013-2016 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki, University of Cambridge
	
*/
//Error message handling 
public class Error {
	public static void error(int code)
	{
		error(code, "");
	}
	
	public static void error(int code, String extrainfo)
	{
		//           options:
		String os = "         ";
		System.err.println("Error " + code); 
		switch (code) {
		case 3002:
			System.err.println("Error: could not load mapping file!");
		case 3001:
			System.err.println("usage: samtools mpileup -q 10 -Q10 -s $(cat sorted_bams)|java Pileup2Likelihoods [options] >post.txt");
			System.err.println("options:");
			System.err.println(os + "pileup=file          load pileup from a file [-]");
			System.err.println(os + "mapping=file         load individual names (same order as bams) from a file [mapping.txt]");
			System.err.println(os + "minCoverage=NUM      minimum coverage per individual [3]");
			System.err.println(os + "numLowerCoverage=NUM number (or proportion) individuals allowed with lower than minCoverage [0.3]");
			System.err.println(os + "minAlleleFreq=NUM    minimum number (or proportion) of an allele [0.1]");
			System.err.println(os + "minQuality=NUM       minimum quality value [0.001]");
			System.err.println(os + "minCoverageSum=NUM   minimum total (summed) coverage");
			System.err.println(os + "ploidy=NUM           ploidy [2]");
			System.err.println(os + "callIndels=1         call indels as well [not set]");
			break;
		case 2001:
			System.err.println("usage: java IBD data=file [options]");
			System.err.println("options:");
			System.err.println(os + "data=file           Loads input data in Lep-MAP3 format with pedigree");
			System.err.println(os + "posteriotFile=file  Loads input data in Lep-MAP3 posterior format");
			System.err.println(os + "vcfFile=file        Loads input data from a vcf");
			System.err.println(os + "MAFLimit=NUM [0.05]");
			System.err.println(os + "genotypeLimit=NUM [0.1]");
			System.err.println(os + "missingLimit=NUM");
			System.err.println(os + "missingLimit=NUM NUM2");
			System.err.println(os + "numThreads=NUM [4]");
			System.err.println(os + "parents=p1 [p2 ...]");
			System.err.println(os + "parents=file:parents.txt");
			System.err.println(os + "allParentChildPairs=1");
			break;
		case 1101:
			System.err.println("usage: java ScaffoldHMM2 scaffold+LG1.txt scaffold+LG2.txt");
			System.err.println("options:");
			System.err.println(os + "two input files give the two input maps with reference mapping information.");
			System.err.println(os + "Columns 1-3 in these files should contain scaffold_name(1), scaffold_pos(2), and LG_name(3)");
			break;
		case 1001:
			System.err.println("usage: java Filtering2 [options] data=file >file1_filtered.post");
			System.err.println("options:");
			System.err.println(os + "data=file          Loads input data in Lep-MAP3 posterior format");
			System.err.println(os + "dataTolerance=NUM  P-value limit for segregation distortion [0.001]");
			System.err.println(os + "removeNonInformative=1  Remove (after filtering) non-informative markers");
			System.err.println(os + "convert2Bialleleic=1    Convert data to biallelic");
			System.err.println(os + "outputHWE=1        Output segregation distortion for the markers");
			System.err.println(os + "MAFLimit=NUM       NUM>=1: Filter out markers with minimum allele frequency < NUM in each family [0]");
			System.err.println(os + "                   NUM<1:  Filter out markers with minimum allele rate < NUM in each family [0]");
			System.err.println(os + "missingLimit=NUM   NUM>=1: Filter out markers with > NUM missing individuals in each family [inf]");
			System.err.println(os + "                   NUM<1:  Filter out markers with missing rate > NUM in each family");
			System.err.println(os + "familyInformativeLimit=NUM Filter out markers with < NUM informative families [0]");
			System.err.println(os + "noSexFiltering=1   Do not filter sex markers, useful for distorted sex ratios [0]");
			System.err.println(os + "heterozygoteRate=NUM    Set heterozygote rate [0.5]");
			System.err.println(os + "                        Useful for selfing data (S2=0.25, S3=0.125, ...)");
			
			break;
		case 1002:
			System.err.println("Error: Cannot find the parent");
			break;
		case 1003:
			System.err.println("Error: Wrong number of columns in the input file");
			break;
		case 1004:
			System.err.println("Error: Unable to load input file or errors in the file");
			break;
		case 1005:
			System.err.println("Error: Internal error");
			break;
		case 1006:
			System.err.println("Error: Unknown sex code");
			break;
		case 1007:
			System.err.println("Error: Pedigree");
			break;
		case 1008:
			System.err.println("Error: Pedigree");
			break;
		case 1009:
			System.err.println("Error: Multiple parents");
			break;
		case 1010:
			System.err.println("Error: Only one lodLimit allowed");
			break;
		case 1011:
			System.err.println("Error: Unable to load map file");
			break;
		case 1012:
			System.err.println("Error: Only one lg paramater allowed");
			break;
		case 1013:
			System.err.println("Error: map file has different number of lines than there are markers");
			break;
		case 1014:
			System.err.println("Error: this lod3Mode (and theta combination) is not supported");
			break;
		case 901:
			System.err.println("usage: java OutputData [options] data=file map=map_file >posterior_clean.txt 2>prints.txt");
			System.err.println("options:");
			System.err.println(os + "map=map_file       LG map file. Typically generated by SeparateIdenticals or JoinSingles2.");
            System.err.println(os + "sizeLimit=NUM      Output data only with LGs with >= NUM markers [1]");
            System.err.println(os + "lod3Mode=NUM       Controls how LOD scores are computed between double informative markers [1]");
			break;
		case 1201:
			System.err.println("usage: java LMPlot map.txt");
			System.err.println("options:");
			System.err.println(os + "map.txt           Output from OrderMarkers2 with outputPhasedData=1 flag");
			System.err.println(os + "selfingPhase=1    calculate pattern distances by ignoring individuals heterozygote in the two markers");
			System.err.println(os + "outputAllDistances=1  output all distances...");
			
			//System.err.println(os + "informativeMask=STR    Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
            //System.err.println(os + "limit1=NUM             Maximum differences in the segregation prints to join markers [1]"); 
            //System.err.println(os + "limit2=NUM             Minimum number of segregative individuals to join markers [1]"); 
			break;			
		case 1301:
			System.err.println("usage: java QTL map.txt traits.txt [options]");
			System.err.println("options:");
			System.err.println(os + "map.txt          Output from OrderMarkers2 with outputPhasedData=1 (or 2) flag");
			System.err.println(os + "traits.txt       Traits, loaded as a single line, columns same order as offspring in OrderMarkers2");
			System.err.println(os + "scale=NUM        round the traits to integers after dividing with scale [0.1]");
			
			//System.err.println(os + "informativeMask=STR    Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
            //System.err.println(os + "limit1=NUM             Maximum differences in the segregation prints to join markers [1]"); 
            //System.err.println(os + "limit2=NUM             Minimum number of segregative individuals to join markers [1]"); 
			break;			
			
			
		case 801:
			System.err.println("usage: java ShortPath [options] prints.txt");
			System.err.println("options:");
			System.err.println(os + "prints.txt               Output (err-stream) from SeparateIdenticals, JoinSingles or OutputData");
			System.err.println(os + "informativeMask=STR      Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
            System.err.println(os + "sizeLimit=NUM            Use only markers that occur NUM times");
            System.err.println(os + "limit1=NUM               Maximum differences in the segregation prints to join markers [1]"); 
            System.err.println(os + "limit2=NUM               Minimum number of segregative individuals to join markers [1]"); 
            System.err.println(os + "begin=NUM                find paths starting from marker NUM");
            System.err.println(os + "end=NUM                  find paths ending to marker NUM");
			break;			
		case 601:
			System.err.println("usage: java JoinIdenticalLGs [options] map=map_file data=file >map.txt 2>prints.txt");
			System.err.println("Note: Joins identically segregating LGs together");
			System.err.println("options:");
			System.err.println(os + "map=map_file       LG map file. Typically generated by SeparateIdenticals or JoinSingles2.");
			System.err.println(os + "data=file          Loads input genotypes in Lep-MAP posterior format");
			System.err.println(os + "lodLimit=NUM       LOD score limit [10.0]");
			System.err.println(os + "informativeMask=STR      Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
			System.err.println(os + "theta=NUM          Fixed recombination fraction [0.0]");
            System.err.println(os + "(fe)maleTheta=NUM  Fixed recombination fraction separately for both sex [theta]");
            System.err.println(os + "betweenSameType=1  Only compute LOD scores between identically informative markers (applicable to single family data)");
            System.err.println(os + "                   Also one has to specify 3 LOD limits (paternal, maternal and both) in this case");
            
            System.err.println(os + "lod3Mode=NUM       Controls how LOD scores are computed between double informative markers [1]");
            System.err.println(os + "sizeLimit=NUM      Joing only LGs with >= NUM markers [1]");
			break;			
		case 701:
			System.err.println("usage: java JoinSingles2[All|Identicals] [options] map=map_file data=file >map.txt 2>prints.txt");
			System.err.println("options:");
			System.err.println(os + "map=map_file       Initial LG map file. Typically generated by SeparateChromosomes2, SeparateIdenticals or JoinSingles2*.");
			System.err.println(os + "lodLimit=NUM       LOD score limit [10.0]");
			System.err.println(os + "lodDifference=NUM  Required LOD difference [0.0]");
			System.err.println(os + "informativeMask=STR      Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
			System.err.println(os + "theta=NUM          Fixed recombination fraction [0.0 (Identicals)] [0.03 (All)]");
            System.err.println(os + "(fe)maleTheta=NUM  Fixed recombination fraction separately for both sex [theta]");
            System.err.println(os + "betweenSameType=1  Only compute LOD scores between identically informative markers (applicable to single family data and Join...Identicals)");
            System.err.println(os + "                   Also one has to specify 3 LOD limits (paternal, maternal and both) in this case");
            System.err.println(os + "lod3Mode=NUM       Controls how LOD scores are computed between double informative markers [1]");
            System.err.println(os + "numThreads=NUM     Use maximum of NUM threads (JoinSingles2All only) [1]");

            System.err.println(os + "maxDistance=NUM    Only calculate LOD scores between this many markers (JoinSingles2All only) [not set]");
            
            System.err.println(os + "distortionLod=1    Use segregation distortion aware LOD scores (JoinSingles2All only) [not set]");
            System.err.println(os + "iterate=1          Iterate single joining until no markers can be added (JoinSingles2All only) [not set]");
            System.err.println(os + "                   (iterating is much faster than running JoinSingles2All multiple times)");
			System.err.println(os + "mask=map_file2     Filter out markers not (in group) 1 in the map_file2");
			break;
		case 702:
			System.err.println("Error: BetweenSameType=1 only works with a single family and three lodLimits");
			break;
		case 703:
			System.err.println("Error: Only one lodLimit accepted");
			break;
		case 501:
			System.err.println("usage: java ParentCall2 [options] data=file");
			System.err.println("options:");
			System.err.println(os + "data=file          Loads genotype posteriors from a file (- for standard input)");
			System.err.println(os + "                   Column 1: contig, Column 2: pos");
			System.err.println(os + "                   Columns 3,...,N+2: pedigree information for lines 1-6, for lines > 6");  
			System.err.println(os + "                   columns 3,...,10*N+2: 10 posteriors for each individual and each genotype combination (AA,AC,..., GT,TT)"); 

			//System.err.println(os + "callLimit=NUM      Required log-odds difference to call a genotype [2.0]");
			System.err.println(os + "familyLimit=NUM    Required log-odds difference to call a SNP [2.0]");
			System.err.println(os + "ZLimit=NUM         Required log-odds difference to call a SNP with Z inheritance [inf]");
			System.err.println(os + "XLimit=NUM         Required log-odds difference to call a SNP with X inheritance [inf]");
			
			System.err.println(os + "removeNonInformative=1  Remove markers that are not informative");
			System.err.println(os + "ignoreParentOrder=1     Do not care about the order of parental genotypes [not set]");
			System.err.println(os + "outputParentPosterior=1  Outputs the genotype likelihoods (posteriors) for the parents as well[not set]");			
			System.err.println(os + "halfSibs=1              Look for identical parent names for half-sib parental genotype inference[not set]");			
			System.err.println(os + "vcfFile=file            Read genotype likelihoods (posteriors) from a vcf file");			
			System.err.println(os + "posteriorFile=file      Read genotype likelihoods (posteriors) from a text file");
			System.err.println(os + "outputRaw=1             Do not call parents, just output the raw data likelihoods (posteriors)");			
			break;
		case 502:
			System.err.println("Error: Cannot find the parent(s) (family " + extrainfo + ")");
			break;
		case 503:
			System.err.println("Error: Wrong number of columns in the input file");
			break;
		case 504:
			System.err.println("Error: Unable to load input file or errors in the file");
			break;
		case 505:
			System.err.println("Error: Internal error");
			break;
		case 506:
			System.err.println("Error: Unknown sex code (family:id " + extrainfo +  ")");
			break;
		case 507:
			System.err.println("Error: Pedigree");
			break;
		case 508:
			System.err.println("Error: Pedigree");
			break;
		case 509:
			System.err.println("Error: Multiple parents (family:id " + extrainfo + ")");
			break;
		case 511:
			System.err.println("Error: Sex of the parent does not match the pedigree (family:id " + extrainfo + ")");
			break;
		case 512:
			System.err.println("Error: Too many parents (family:id " + extrainfo + ")");
			break;
		case 513:
			System.err.println("Error: Grandparents do not match the pedigree or are not present (family:id " + extrainfo + ")");
			break;
		case 514:
			System.err.println("Error: Only one of vcfFile or posteriorFile parameters can be provided");
			break;
		case 515:
			System.err.println("Error: Same individual is twice in the data");
			break;
		case 517:
			System.err.println("Error: Individual names not found from the vcf file");
			break;
		case 518:
			System.err.println("Error: vcf file does not contain such field");
			break;
		case 519:
			System.err.println("Error: vcf file does not contain any of PL/GL/GT fields");
			break;
		case 520:
			System.err.println("Error: No parent(s) in a family " + extrainfo);
			break;
		case 1401:
			System.err.println("usage: java SeparateChromosomes2 [options] data=file >map.txt");
			System.err.println("options:");
			System.err.println(os + "data=file          Loads genotype posteriors from a file (- for standard input)");
			System.err.println(os + "                   Column 1: contig, Column 2: pos");
			System.err.println(os + "                   Columns 3,...,N+2: pedigree information for lines 1-6, for lines > 6");  
			System.err.println(os + "                   Columns 3,...,10*N+2: 10 posteriors for each individual and each genotype combination (AA,AC,..., GT,TT)"); 
			System.err.println(os + "lodLimit=NUM       LOD score limit [10.0]");
			System.err.println(os + "informativeMask=STR     Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
            System.err.println(os + "families=F1 [F2 ...]  Use only some families [not set]");

			System.err.println(os + "theta=NUM          Fixed recombination fraction [0.03]");
            System.err.println(os + "(fe)maleTheta=NUM  Fixed recombination fraction separately for both sex [theta]");
            System.err.println(os + "sizeLimit=NUM      Remove LGs with < NUM markers [1]");
            System.err.println(os + "numThreads=NUM     Use maximum of NUM threads [1]");
            System.err.println(os + "subsample=NUM      Use only a random NUM fraction of markers [1] (speedup is 1/NUM^2)");
            System.err.println(os + "samplePairs=NUM    Use only a random NUM fraction of marker pairs [1] (speedup is 1/NUM)");
            
            System.err.println(os + "phasedData=1       Data is phased [not set]");
            System.err.println(os + "grandparentPhase=1 Pphase data based on grandparents [not set]");
            
            System.err.println(os + "lod3Mode=NUM       Controls how LOD scores are computed between double informative markers [1]");
            System.err.println(os + "                   1: haplotypes match (4 alleles, max LOD = log(4^n)");
            System.err.println(os + "                   2: homozygotes match (3 alleles, max LOD = log(3^n))");
            System.err.println(os + "                   3: homozygotes or heterozygotes match (2 alleles, max LOD = log(2^n))");

            System.err.println(os + "distortionLod=1    Use segregation distortion aware LOD scores [not set]");

            System.err.println(os + "map=file           refine linkage group lg of map file");
            System.err.println(os + "lg=NUM             refine linkage group lg [1 if map is provided]");
            System.err.println(os + "renameLGs=0        do not rename linkage groups after refine");

            System.err.println(os + "minLod=NUM         minimum LOD value for each family (used for multi-family maps without parents) [not set]");
            break;			
		case 1402:
		case 1403:
			System.err.println("Error: Unable to load input file"); 
			break;			
		case 401:
			System.err.println("usage: java SeparateIdenticals [options] data=file >map.txt 2>prints.txt");
			System.err.println("options:");
			System.err.println(os + "data=file          Loads genotype posteriors from a file (- for standard input)");
			System.err.println(os + "                   Column 1: contig, Column 2: pos");
			System.err.println(os + "                   Columns 3,...,N+2: pedigree information for lines 1-6, for lines > 6");  
			System.err.println(os + "                   Columns 3,...,10*N+2: 10 posteriors for each individual and each genotype combination (AA,AC,..., GT,TT)"); 
			System.err.println(os + "lodLimit=NUM       LOD score limit [10.0]");
			System.err.println(os + "informativeMask=STR     Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
			System.err.println(os + "theta=NUM          Fixed recombination fraction [0.0]");
            System.err.println(os + "(fe)maleTheta=NUM  Fixed recombination fraction separately for both sex [0.0]");
            System.err.println(os + "sizeLimit=NUM      Remove LGs with < NUM markers [1]");
            System.err.println(os + "numThreads=NUM     Use maximum of NUM threads (most speedup if equals numParts) [1]");
            System.err.println(os + "numParts=NUM       Divide markers to NUM parts [1]");
            System.err.println(os + "removeSingles=0    Do not remove single markers (slower but does not miss markers)");
            System.err.println(os + "keepRate=NUM       Keep joined markers with this prob (1.0 = all pair-wise comparisons are done) [1.0]");
            System.err.println(os + "betweenSameType=1  Only compute LOD scores between identically informative markers (applicable to single family data)");
            System.err.println(os + "                   Also one has to specify 3 LOD limits (paternal, maternal and both) in this case");
            
            System.err.println(os + "lod3Mode=NUM       Controls how LOD scores are computed between double informative markers [1]");
            System.err.println(os + "                   1: haplotypes match (4 alleles, max LOD = log(4^n)");
            System.err.println(os + "                   2: homozygotes match (3 alleles, max LOD = log(3^n))");
            System.err.println(os + "                   3: homozygotes or heterozygotes match (2 alleles, max LOD = log(2^n))");
            break;
		case 402:
		case 403:
			System.err.println("Error: Unable to load input file"); 
			break;			
		case 301:
//			pp.warning(new String[]{});
			
			System.err.println("usage: java OrderMarkers2 [options] data=file.posterior");
			System.err.println(os + "data=file          Loads genotype posteriors from a file (- for standard input)");
			System.err.println(os + "                   Column 1: contig, Column 2: pos");
			System.err.println(os + "                   Columns 3,...,N+2: pedigree information for lines 1-6, for lines > 6");  
			System.err.println(os + "                   Columns 3,...,10*N+2: 10 posteriors for each individual and each genotype combination (AA,AC,..., GT,TT)"); 
			
			System.err.println(os + "map=chromosome_map_file LG map file. Typically generated by SeparateChromosomes2 or JoinSingles2.");
			System.err.println(os + "evaluateOrder=order.txt Load initial marker order (single chromosome) from a file");
			
			System.err.println(os + "informativeMask=STR     Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]");
            System.err.println(os + "useMorgan=1        Use Morgan (linear) mapping function");
            System.err.println(os + "useKosambi=1       Use Kosambi mapping function");
            System.err.println(os + "improveOrder=0     Do not improve the order (used to only (re)evaluate an order)");
            System.err.println(os + "numThreads=NUM     Use NUM threads [1]");
            System.err.println(os + "numMergeIterations=NUM  Run NUM iterations [6]");
            System.err.println(os + "chromosome=NUM     Order chromosome NUM only [all]");
            System.err.println(os + "families=F1 [F2 ...]  Use only some families [not set]");
            

            System.err.println(os + "scale=NUM NUM2     Scale posteriors by NUM (p -> p^NUM) with a maximum of NUM2 (>=NUM) times map 'end effect' correction [M/N 2]");
            System.err.println(os + "                   , where N is number of markers and M number of individuals (e.g. 3M/N 3, 100/N 3)");
            System.err.println(os + "scale=NUM          Scale posteriors \"NUM 2\"");
            System.err.println(os + "scale=NUM NUM2 NUM3  same as \"NUM NUM2\" but cap maximum scale to NUM3");

            System.err.println(os + "scaleMode=1        Use the old data scaling mode [2]");

            System.err.println(os + "minError=NUM       Set minimum posterior value [0.001]");
            System.err.println(os + "outputPhasedData=0 Do not output phased data");
            System.err.println(os + "outputPhasedData=1 Output phased data [1]");
            System.err.println(os + "outputPhasedData=2 Output phased data but mask uncertain haplotypes");
            System.err.println(os + "outputPhasedData=3 Output phased data + haplotype likelihoods");
            System.err.println(os + "outputPhasedData=4 Output phased data but mask uncertain haplotypes + haplotype likelihoods");
            System.err.println(os + "removeMarkers=m1 [ m2 m3 ...]  Remove markers");
            System.err.println(os + "sexAveraged=1      Calculate sex-averaged map distances");
            System.err.println(os + "phasedData=1       Input data is phased");
            System.err.println(os + "grandparentPhase=1 Use grandparents to phase data, removes markers that cannot be phased");

            System.err.println(os + "selfingPhase=1     Phase data so that homozygotes and heterozygotes are kept as in selfing crosses(experimental)");

            System.err.println(os + "recombination1=NUM Recombination rate for male [0.001]");
            System.err.println(os + "recombination2=NUM Recombination rate for female [0.001]");
            System.err.println(os + "interference1=NUM  Recombination interference for male [0.001]");
            System.err.println(os + "interference2=NUM  Recombination interference for female [0.001]");
            
            System.err.println(os + "identicalLimit=NUM Reduce the number of markers (conditional on the order)");
            System.err.println(os + "                   If the absolute probability difference between markers is < NUM they are collapsed [0.01]"); 
            System.err.println(os + "computeLODScores=file   Evaluate pair-wise LOD scores and store to file");
            System.err.println(os + "calculateIntervals=file [NUM=1] Evaluate and store to a file the interval(s) for each marker where it could be located");
            System.err.println(os + "                    within NUM likelihood treshold. Useful for matching physical and linkage postions of markers"); 
                        
            System.err.println(os + "randomPhase=1      Start the phasing algorithm from a random phase configuration");
            System.err.println(os + "                   Useful if the phasing does not converge properly with evaluateOrder"); 
            System.err.println(os + "hyperPhaser=1      Use 'hyper' (instead of super) phasing algorithm.");
            System.err.println(os + "                   Useful if the phasing does not converge properly"); 
            System.err.println(os + "phasingIterations=NUM   Run NUM phasing iterations [1].");
            System.err.println(os + "                   Useful if the phasing does not converge properly"); 

            System.err.println(os + "usePhysical=1 NUM  Use physical positions in the marker ordering, [not set, NUM=0.01]");
            System.err.println(os + "                   penalise adjacent markers in different contigs by NUM"); 
            
            System.err.println(os + "maskIgnoreParentOrder=1  Mask markers where the order of parental genotypes is not clear");
            System.err.println(os + "                         (ignoreParentOrder=1 in ParentCall2)");

            System.err.println(os + "refineParentOrder=1  Refine the unclear order of parental genotypes");
            System.err.println(os + "                     (requires evaluateOrder, ignoreParentOrder=1 in ParentCall2)");
            
			break;
		case 4406:
			System.err.println("Error: a chromosome map file must be provided (map=file)");
			break;
		case 302:
			System.err.println("Error: unable to load chromosome map file");
			break;
		case 303:
			System.err.println("Error: either evaluateOrder or map parameter must be provided");
			break;
		case 304:
			System.err.println("Error: unable to parse scale parameter");
			break;
		case 305:
			System.err.println("Error: do not use maskIgnoreParentOrder=1 with refineParentOrder=1 or grandParentPhase=1");
			break;
		case 101: 
			System.err.println("Error: Multiple IDs in the input file"); 
			break;
		case 102: 
			System.err.println("Error: Sex of the parent(s) is not specified"); 
			break;
		case 103: 
			System.err.println("Error: Unable to load input file"); 
			break;			
		case 298: 
		case 299: 
			System.err.println("Error: Trying to add multiple genotype data for a individual"); 
			break;			
		default:
			System.err.println("Error " + code); 
		}
		System.exit(-1);
	}
	public static void warning(int code)
	{
		switch (code) {
		case 3003:
			System.err.println("Warning: different length of pileup quality and alignment columns (add \"-q 1\" to samtools mpileup)!");
			System.err.println("cutting to min length, strings below:");
			break;
//			case 101: 
//				System.err.println("Warning 1"); 
//				break;
//			case 102: 
//				System.err.println("Warning 2"); 
//				break;
			default:
				System.err.println("Warning " + code); 
		}
	}
}
