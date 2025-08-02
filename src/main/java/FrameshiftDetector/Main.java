package FrameshiftDetector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import mapper.AlignmentParameters;
import mapper.AlignedBlock;
import mapper.Logger;
import mapper.QueryAlignment;

public class Main {
  private static void usage(String error) {
    System.out.println(error);
    System.out.println("Usage: java -jar ProteinFrameshiftDetector.jar [--dna <dna>] [--dnas <filepath> [--protein <protein>] [--proteins <filepath>]");
    System.out.println("Looks for frameshifts in the given proteins to make them look like the given DNA");
    System.exit(1);
  }

  public static void main(String[] args) throws IOException {
    List<String> dnas = new ArrayList<String>();
    List<String> dnaNames = new ArrayList<String>();
    List<String> proteins = new ArrayList<String>();
    List<String> proteinNames = new ArrayList<String>();
    for (int i = 0; i < args.length; i++) {
      String arg = args[i];
      if (arg.equals("--dna")) {
        i++;
        dnas.add(args[i]);
        dnaNames.add("dna" + dnas.size());
        continue;
      }
      if (arg.equals("--dnas")) {
        i++;
        String filepath = args[i];
        Map<String, String> sequences = FastaParser.parseFasta(new File(filepath));
        for (Map.Entry<String, String> entry: sequences.entrySet()) {
          dnas.add(entry.getValue());
          dnaNames.add(entry.getKey());
        }
        continue;
      }
      if (arg.equals("--protein")) {
        i++;
        proteins.add(args[i]);
        proteinNames.add("protein" + proteins.size());
        continue;
      }
      if (arg.equals("--proteins")) {
        i++;
        String filepath = args[i];
        Map<String, String> sequences = FastaParser.parseFasta(new File(filepath));
        for (Map.Entry<String, String> entry: sequences.entrySet()) {
          proteins.add(entry.getValue());
          proteinNames.add(entry.getKey());
        }

        continue;
      }

      usage("Unrecognized argument '" + arg + "'");
    }
    System.out.println("\nThis program isn't optimized but hopefully is helpful\n");
    for (int i = 0; i < proteins.size(); i++) {
      String protein = proteins.get(i);
      for (int j = 0; j < dnas.size(); j++) {
        String dna = dnas.get(j);
        System.out.println("");
        System.out.println("Comparing dnas[" + i + "] and proteins[" + j + "]");
        System.out.println("DNA = " + dna);
        System.out.println("Protein = " + protein);
        process(dna, dnaNames.get(j), protein, proteinNames.get(i));
      }
    }
  }

  private static void process(String dna, String dnaName, String protein, String proteinName) {
    List<String> equivalentDNAComponents = convertProteinToDNAs(dna, protein);
    String equivalentDNA = String.join("", equivalentDNAComponents);
    System.out.println("Converted protein " + proteinName + " to resemble dna " + dnaName);
    System.out.println("Original protein:" + protein);
    System.out.println("Protein to DNA: " + equivalentDNA);
    System.out.println("Original DNA: " + dna);
    Map<String, String> namedDNA = new HashMap<String, String>();
    namedDNA.put(dnaName, dna);
    mapper.ReferenceDatabase reference = mapper.Api.newDatabase(namedDNA, Logger.NoOpLogger);
    // split dna into pieces to align in case of breaks
    int numPieces = equivalentDNA.length() / 150 + 1;
    int cumulativeLength = 0;
    List<AlignedBlock> indels = new ArrayList<AlignedBlock>();
    for (String component: equivalentDNAComponents) {
      int startIndex = cumulativeLength;
      int endIndex = cumulativeLength + component.length();
      System.out.println("Trying to align " + proteinName + ".dna[" + startIndex + ":" + endIndex + "] = " + component);
      System.out.println("Translated protein[" + (startIndex / 3) + ":" + (endIndex/3) + "]: " + dnaToProtein(component));
      List<QueryAlignment> alignments = mapper.Api.align(component, reference, alignmentParameters(), Logger.NoOpLogger);
      if (alignments.size() < 1) {
        System.out.println("found no alignments");
      } else {
        QueryAlignment alignment = alignments.get(0);
        System.out.println(alignment.format());
        for (AlignedBlock block: alignment.getComponent(0).getSections()) {
          if (block.getIndelLength() != 0)
            indels.add(block);
        }
      }
      cumulativeLength = endIndex;
    }
    if (indels.size() < 1) {
      System.out.println("Summary for " + proteinName + " compared to " + dnaName + ": no indels");
    } else {
      System.out.println("Summary for " + proteinName + " compared to " + dnaName + ":");
      for (AlignedBlock block: indels) {
        if (block.getLengthA() > 0) {
          System.out.println("Insertion of length " + block.getLengthA() + "bp at " + block.getStartIndexB());
        } else {
          System.out.println("Deletion of length " + block.getLengthA() + "bp at " + block.getStartIndexB());
        }
      }
    }
  }

  private static AlignmentParameters alignmentParameters() {
    AlignmentParameters parameters = new AlignmentParameters();
    parameters.MutationPenalty = 5;
    parameters.InsertionStart_Penalty = 2;
    parameters.InsertionExtension_Penalty = 0.5;
    parameters.DeletionStart_Penalty = 2;
    parameters.DeletionExtension_Penalty = 0.5;
    parameters.MaxErrorRate = 0.3;
    parameters.AmbiguityPenalty = 0.1;
    parameters.UnalignedPenalty = 0.1;
    return parameters;
  }

  private static String extractNeighborhood(String text, int position, int radius) {
    int startIndex = Math.max(position - radius, 0);
    int endIndex = Math.min(position + radius, text.length());
    return text.substring(startIndex, endIndex);
  }


  // adds an insertion of <insertionLength> at position <index>
  // If <insertionLength> is negative, creates a deletion instead
  private static String addFrameshift(String dna, int index, int insertionLength) {
    StringBuilder builder = new StringBuilder();
    builder.append(dna.substring(0, index + insertionLength));
    builder.append(dna.substring(index));
    return builder.toString();
  }

  private static String dnaToProtein(String dna) {
    int numCodons = dna.length() / 3;
    return Codons.DNAToProtein(dna.substring(0, numCodons * 3));
  }

  private static int logRoundUp(int value, int radix) {
    int result = 0;
    while (value > 1) {
      result++;
      value /= radix;
    }
    return result;
  }

  private static int getKmerLength(String gene) {
    return logRoundUp(gene.length(), 20) + 2;
  }

  // converts a protein into DNA and breaks it at the positions where we're less confident of what it should be
  private static List<String> convertProteinToDNAs(String gene, String protein) {
    // generate proteins for each frameshift
    String gene0 = gene;
    String gene1 = addFrameshift(gene0, 0, 1);
    String gene2 = addFrameshift(gene0, 0, 2);
    String gene0Comp = reverseComplementDNA(gene0);
    String gene1Comp = addFrameshift(gene0Comp, 0, 1);
    String gene2Comp = addFrameshift(gene0Comp, 0, 2);

    // keep track of which protein kmers correspond to which DNA kmers
    int proteinKmerLength = getKmerLength(gene);
    Map<String, String> proteinToDNA = new HashMap<String, String>();
    collectProteinToDNA(gene0, proteinKmerLength, proteinToDNA);
    collectProteinToDNA(gene1, proteinKmerLength, proteinToDNA);
    collectProteinToDNA(gene2, proteinKmerLength, proteinToDNA);
    collectProteinToDNA(gene0Comp, proteinKmerLength, proteinToDNA);
    collectProteinToDNA(gene1Comp, proteinKmerLength, proteinToDNA);
    collectProteinToDNA(gene2Comp, proteinKmerLength, proteinToDNA);

    // attempt to convert protein kmers to DNA kmers
    List<String> components = new ArrayList<String>();
    StringBuilder equivalentDNA = new StringBuilder();
    String pendingDNAKmer = null;
    for (int proteinIndex = 0; proteinIndex < protein.length(); proteinIndex++) {
      String proteinKmer = null;
      String dnaKmer = null;
      if (proteinIndex + proteinKmerLength <= protein.length()) {
        proteinKmer = protein.substring(proteinIndex, proteinIndex + proteinKmerLength);
        dnaKmer = proteinToDNA.get(proteinKmer);
        //System.out.println("proteinKmer " + proteinKmer + " comes from dnaKmer " + dnaKmer);
      }
      if (dnaKmer != null) {
        // we got a unique match, so we can add another codon
        equivalentDNA.append(dnaKmer.substring(0, 3));
        pendingDNAKmer = dnaKmer.substring(3);
      } else {
        // We didn't get a unique match, so this is a good place to break
        if (equivalentDNA.length() >= 99) {
          components.add(equivalentDNA.toString());
          equivalentDNA = new StringBuilder();
        }

        if (pendingDNAKmer != null && pendingDNAKmer.length() >= 3) {
          // We didn't get a unique match but we can reuse the previous DNA kmer
          equivalentDNA.append(pendingDNAKmer.substring(0, 3));
          pendingDNAKmer = pendingDNAKmer.substring(3);
        } else {
          // We aren't sure what the DNA should be here
          equivalentDNA.append("NNN");
        }
      }
    }
    if (equivalentDNA.length() > 0) {
      components.add(equivalentDNA.toString());
    }
    return components;
  }

  private static String reverseComplementDNA(String text) {
    StringBuilder builder = new StringBuilder();
    for (int i = text.length() - 1; i >= 0; i--) {
      char here = text.charAt(i);
      char complement;
      if (here == 'A') {
        complement = 'T';
      } else {
        if (here == 'C') {
          complement = 'G';
        } else {
          if (here == 'G') {
            complement = 'C';
          } else {
            if (here == 'T') {
              complement = 'A';
            } else {
              complement = here;
            }
          }
        }
      }
      builder.append(complement);
    }
    return builder.toString();
  }

  private static void collectProteinToDNA(String gene, int proteinKmerLength, Map<String, String> results) {
    int dnaKmerLength = proteinKmerLength * 3;

    String protein = dnaToProtein(gene);
    int max = gene.length() - dnaKmerLength;

    for (int i = 0; i < max; i += 3) {
      String dnaKmer = gene.substring(i, i + dnaKmerLength);
      String proteinKmer = dnaToProtein(dnaKmer);
      if (results.containsKey(proteinKmer) && !dnaKmer.equals(results.get(proteinKmer))) {
        //System.out.println("proteinKmer " + proteinKmer + " comes from " + dnaKmer + ": conflict");
        results.put(proteinKmer, null); // multiple dna kmers for the same protein kmer
      } else {
        //System.out.println("proteinKmer " + proteinKmer + " comes from " + dnaKmer + ": unique so far");
        results.put(proteinKmer, dnaKmer);
      }
    }
  }

}
