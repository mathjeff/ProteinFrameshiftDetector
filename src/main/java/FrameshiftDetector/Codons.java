package FrameshiftDetector;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Codons
{
    private static List<Codon> all;
    public static List<Codon> getAll()
    {
        {
            if (all == null)
            {
                all = listAll();
            }
            return all;
        }
    }

    private static String stop = "Stop";

    private static Map<String, Codon> codonsByDNA;
    public static Map<String, Codon> getCodonsByDNA()
    {
        {
            if (codonsByDNA == null)
            {
                Map<String, Codon> dict = new HashMap<String, Codon>();
                for (Codon codon: getAll())
                {
                    if (dict.containsKey(codon.getDNA()))
                    {
                        throw new IllegalArgumentException("Multiple codons registered with the same dna " + codon.getDNA());
                    }
                    dict.put(codon.getDNA(), codon);
                }
                codonsByDNA = dict;
            }
            return codonsByDNA;
        }
    }
    public static Codon GetByDNA(String text)
    {
        text = text.replace("U", "T");
        Codon result = codonsByDNA.get(text.replace("U", "T"));
        return result;
    }

    private static Map<String, List<Codon>> codonsByProtein;
    public static Map<String, List<Codon>> getCodonsByProtein()
    {
        {
            if (codonsByProtein == null)
            {
                Map<String, List<Codon>> dict = new HashMap<String, List<Codon>>();
                for (Codon codon: getAll())
                {
                    if (!dict.containsKey(codon.getProtein()))
                    {
                        dict.put(codon.getProtein(), new ArrayList<Codon>());
                    }
                    dict.get(codon.getProtein()).add(codon);
                }
                codonsByProtein = dict;
            }
            return codonsByProtein;
        }
    }

    private static List<String> sortedValues;
    public static List<String> getSortedValues()
    {
        {
            if (sortedValues == null)
            {
                List<String> values = new ArrayList<String>();
                for (String item: codonsByProtein.keySet())
                {
                    values.add(item);
                }
                values.sort(null);
                sortedValues = values;
            }
            return sortedValues;
        }
    }

    public static int getNumUniqueValues()
    {
        {
            return getSortedValues().size();
        }
    }

    public static String DNAToProtein(String dna)
    {
        if (dna.length() % 3 != 0)
            throw new IllegalArgumentException("Attempted to convert dna to protein but its length " + dna.length() + " is not divisible by 3. Full content: " + dna);
        String result = "";
        for (int i = 0; i < dna.length(); i += 3)
        {
            String codonText = dna.substring(i, i + 3);
            Codon codon = Codons.getCodonsByDNA().get(codonText);
            if (codon == null)
              throw new IllegalArgumentException("Cannot find codon with text " + codonText + " from text " + dna);
            String amino = codon.getProtein();
            if (!amino.equals(stop))
            {
                result += amino;
            }
        }
        return result;
    }
    private static List<Codon> listAll()
    {
        List<Codon> result = new ArrayList<Codon>();
        result.add(new Codon("K", "AAA"));
        result.add(new Codon("N", "AAC"));
        result.add(new Codon("K", "AAG"));
        result.add(new Codon("N", "AAT"));

        result.add(new Codon("T", "ACA"));
        result.add(new Codon("T", "ACC"));
        result.add(new Codon("T", "ACG"));
        result.add(new Codon("T", "ACT"));

        result.add(new Codon("R", "AGA"));
        result.add(new Codon("S", "AGC"));
        result.add(new Codon("R", "AGG"));
        result.add(new Codon("S", "AGT"));

        result.add(new Codon("I", "ATA"));
        result.add(new Codon("I", "ATC"));
        result.add(new Codon("M", "ATG"));
        result.add(new Codon("I", "ATT"));

        result.add(new Codon("Q", "CAA"));
        result.add(new Codon("H", "CAC"));
        result.add(new Codon("Q", "CAG"));
        result.add(new Codon("H", "CAT"));

        result.add(new Codon("P", "CCA"));
        result.add(new Codon("P", "CCC"));
        result.add(new Codon("P", "CCG"));
        result.add(new Codon("P", "CCT"));

        result.add(new Codon("R", "CGA"));
        result.add(new Codon("R", "CGC"));
        result.add(new Codon("R", "CGG"));
        result.add(new Codon("R", "CGT"));

        result.add(new Codon("L", "CTA"));
        result.add(new Codon("L", "CTC"));
        result.add(new Codon("L", "CTG"));
        result.add(new Codon("L", "CTT"));

        result.add(new Codon("E", "GAA"));
        result.add(new Codon("D", "GAC"));
        result.add(new Codon("E", "GAG"));
        result.add(new Codon("D", "GAT"));

        result.add(new Codon("A", "GCA"));
        result.add(new Codon("A", "GCC"));
        result.add(new Codon("A", "GCG"));
        result.add(new Codon("A", "GCT"));

        result.add(new Codon("G", "GGA"));
        result.add(new Codon("G", "GGC"));
        result.add(new Codon("G", "GGG"));
        result.add(new Codon("G", "GGT"));

        result.add(new Codon("V", "GTA"));
        result.add(new Codon("V", "GTC"));
        result.add(new Codon("V", "GTG"));
        result.add(new Codon("V", "GTT"));

        result.add(new Codon(stop, "TAA"));
        result.add(new Codon("Y", "TAC"));
        result.add(new Codon(stop, "TAG"));
        result.add(new Codon("Y", "TAT"));

        result.add(new Codon("S", "TCA"));
        result.add(new Codon("S", "TCC"));
        result.add(new Codon("S", "TCG"));
        result.add(new Codon("S", "TCT"));

        result.add(new Codon(stop, "TGA"));
        result.add(new Codon("C", "TGC"));
        result.add(new Codon("W", "TGG"));
        result.add(new Codon("C", "TGT"));

        result.add(new Codon("L", "TTA"));
        result.add(new Codon("F", "TTC"));
        result.add(new Codon("L", "TTG"));
        result.add(new Codon("F", "TTT"));

        return result;
    }
}

