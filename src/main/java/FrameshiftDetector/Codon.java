
package FrameshiftDetector;

class Codon
{
    public Codon(String protein, String dna)
    {
        this.protein = protein;
        this.dna = dna;
        if (dna.length() != 3)
        {
            throw new IllegalArgumentException("codon dna must be of length 3, got '" + dna + "'");
        }
    }
    public String getProtein()
    {
        {
            return this.protein;
        }
    }
    public String getDNA()
    {
        {
            return this.dna;
        }
    }

    @Override public String toString()
    {
        return this.dna + "(" + this.protein + ")";
    }

    private String protein;
    private String dna;
}

