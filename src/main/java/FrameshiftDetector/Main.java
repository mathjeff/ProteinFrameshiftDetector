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
    for (int i = 0; i < proteins.size(); i++) {
      String protein = proteins.get(i);
      for (int j = 0; j < dnas.size(); j++) {
        String dna = dnas.get(j);
        System.out.println("");
        System.out.println("Comparing dnas[" + i + "] and proteins[" + j + "]");
        process(dna, j, protein, i);
      }
    }
  }

  private static void process(String dna, int dnaIndex, String protein, int proteinIndex) {
    int bestInsertionIndex = 0;
    int bestInsertionLength = 0;
    String bestMutatedDNA = "";
    double bestKmerMatchRate = 0;
    // check each possible frameshift
    for (int insertionLength = -2; insertionLength <= 2; insertionLength++) {
      int minIndex = Math.max(insertionLength, 2);
      int maxIndex = Math.min(dna.length() - insertionLength, dna.length());

      String shiftedLow = addFrameshift(dna, minIndex, insertionLength);
      String translatedLow = dnaToProtein(shiftedLow);
      double similarityLow = compareProteins(protein, translatedLow);

      String shiftedHigh = addFrameshift(dna, maxIndex, insertionLength);
      String translatedHigh = dnaToProtein(shiftedHigh);
      double similarityHigh = compareProteins(protein, translatedHigh);

      // binary search to find the best position to put this frameshift
      while (maxIndex > minIndex + 1) {
        int middleIndex = (maxIndex + minIndex) / 2;

        String shiftedMiddle = addFrameshift(dna, middleIndex, insertionLength);
        String translatedMiddle = dnaToProtein(shiftedMiddle);
        double similarityMiddle = compareProteins(protein, translatedMiddle);

        if (similarityLow < similarityHigh) {
          minIndex = middleIndex;
          shiftedLow = shiftedMiddle;
          translatedLow = translatedMiddle;
          similarityLow = similarityMiddle;
        } else {
          minIndex = middleIndex;
          shiftedHigh = shiftedMiddle;
          translatedHigh = translatedMiddle;
          similarityHigh = similarityMiddle;
        }

        if (similarityMiddle > bestKmerMatchRate) {
          bestKmerMatchRate = similarityMiddle;
          bestInsertionIndex = middleIndex;
          bestInsertionLength = insertionLength;
          bestMutatedDNA = shiftedMiddle;
        }
      }
    }
    // convert number of matching kmers to similarity
    double kmerMatchRate = bestKmerMatchRate;
    String translatedProtein = dnaToProtein(bestMutatedDNA);
    double bestSimilarity = Math.pow(kmerMatchRate, 1.0 / (double)getKmerLength(protein, translatedProtein));

    System.out.println("Most similar result for DNA " + dnaIndex + " and protein " + proteinIndex + " is a frameshift of length " + bestInsertionLength + " at " + bestInsertionIndex + " with similarity of about " + bestSimilarity);
    System.out.println("Original DNA: " + dna);
    System.out.println("Shifted  DNA: " + bestMutatedDNA);
    System.out.println("Translated  : " + translatedProtein);
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

  private static int getKmerLength(String protein1, String protein2) {
    return logRoundUp(Math.max(protein1.length(), protein2.length()), 20) + 2;
  }

  private static double compareProteins(String a, String b) {
    if (a.length() < b.length())
      return hashCompareProteins(a, b);
    else
      return hashCompareProteins(b, a);
  }

  private static double hashCompareProteins(String a, String b) {
    int numMatches = 0;
    int numMismatches = 0;
    int kmerLength = getKmerLength(a, b);
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
