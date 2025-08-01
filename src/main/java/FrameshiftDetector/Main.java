package FrameshiftDetector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Main {
  private static void usage(String error) {
    System.out.println(error);
    System.out.println("Usage: java -jar ProteinFrameshiftDetector.jar [--dna <dna>] [--dnas <filepath> [--protein <protein>] [--proteins <filepath>]");
    System.out.println("Looks for frameshifts in the given proteins to make them look like the given DNA");
    System.exit(1);
  }

  public static void main(String[] args) throws IOException {
    List<String> dnas = new ArrayList<String>();
    List<String> proteins = new ArrayList<String>();
    for (int i = 0; i < args.length; i++) {
      String arg = args[i];
      if (arg.equals("--dna")) {
        i++;
        dnas.add(args[i]);
        continue;
      }
      if (arg.equals("--dnas")) {
        i++;
        String filepath = args[i];
        Map<String, String> sequences = FastaParser.parseFasta(new File(filepath));
        for (String sequence: sequences.values()) {
          dnas.add(sequence);
        }
        continue;
      }
      if (arg.equals("--protein")) {
        i++;
        proteins.add(args[i]);
        continue;
      }
      if (arg.equals("--proteins")) {
        i++;
        String filepath = args[i];
        Map<String, String> sequences = FastaParser.parseFasta(new File(filepath));
        for (String sequence: sequences.values()) {
          proteins.add(sequence);
        }
        continue;
      }

      usage("Unrecognized argument '" + arg + "'");
    }
    System.out.println("\nThis program isn't optimized but hopefully is helpful\n");
    for (String protein: proteins) {
      for (String dna: dnas) {
        System.out.println("Comparing " + dna + " and " + protein);
        process(dna, protein);
      }
    }
  }

  private static void process(String dna, String protein) {
    int bestInsertionIndex = 0;
    int bestInsertionLength = 0;
    String bestMutatedDNA = "";
    double bestSimilarity = 0;
    for (int insertionIndex = 0; insertionIndex < dna.length() - 3; insertionIndex++) {
      for (int insertionLength = -2; insertionLength <= 2; insertionLength++) {
        if (insertionLength == 0 && insertionIndex > 0) {
          continue;
        }
        if (insertionIndex + insertionLength < 0) {
          continue;
        }
        if (insertionIndex + insertionLength >= dna.length()) {
          continue;
        }
        String shifted = addFrameshift(dna, insertionIndex, insertionLength);
        //System.out.println("Checking frameshift of " + insertionLength + " at " + insertionIndex + " = " + shifted);
        String translated = dnaToProtein(shifted);
        //System.out.println("Translated to " + translated);
        double similarity = compareProteins(protein, translated);
        //System.out.println("Similarity = " + similarity);
        if (similarity > bestSimilarity) {
          bestSimilarity = similarity;
          bestInsertionIndex = insertionIndex;
          bestInsertionLength = insertionLength;
          bestMutatedDNA = shifted;
        }
      }
    }
    System.out.println("");
    System.out.println("Most similar result is a frameshift of length " + bestInsertionLength + " at " + bestInsertionIndex + " with similarity of about " + bestSimilarity);
    System.out.println("Original DNA: " + dna);
    System.out.println("Shifted  DNA: " + bestMutatedDNA);
    System.out.println("Translated  : " + dnaToProtein(bestMutatedDNA));
    System.out.println("Protein     : " + protein);
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

  private static int getKmerLength(String protein) {
    return logRoundUp(protein.length(), 20) + 2;
  }

  private static double compareProteins(String a, String b) {
    int numMatches = 0;
    int numMismatches = 0;
    int kmerLength = getKmerLength(b);
    List<String> kmersA = extractProteinKmers(a, kmerLength);
    Set<String> kmersB = new HashSet<String>(extractProteinKmers(b, kmerLength));
    for (String kmer: kmersA) {
      if (kmersB.contains(kmer)) {
        numMatches++;
      } else {
        numMismatches++;
      }
    }
    double similarity = (double)numMatches / (double)(numMatches + numMismatches);
    //System.out.println("Compared " + a + " and " + b + " and got " + similarity);
    return similarity;
  }

  private static List<String> extractProteinKmers(String protein, int length) {
    List<String> kmers = new ArrayList<String>();
    int maxIndex = protein.length() - length + 1;
    for (int i = 0; i < maxIndex; i++) {
      String kmer = protein.substring(i, i + length);
      kmers.add(kmer);
    }
    //System.out.println("extracted " + kmers.size() + " kmers of length " + length);
    return kmers;
  }

}
