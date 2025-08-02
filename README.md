ProteinFrameshiftDetector is a prototype project for finding a single small frameshift that explains the difference between a DNA sequence and a protein sequence

To use it, do one of these:

```
java -jar ProteinFrameshiftDetector.jar --dna <DNASequence> --protein <ProteinSequence>
java -jar ProteinFrameshiftDetector.jar --dnas dnas.fasta --protein proteins.fasta
```

For example:

```
 $ java -jar build/libs/ProteinFrameshiftDetector.jar --dna ATGCGAGTGTTGAAGTTCGGCGGTACATCAGGGCAAATGCAGAACGTTTTCTGCGGGTT --protein MRVLKFGGTSVANAERFLRV

Comparing dnas[0] and proteins[0]
DNA = ATGCGAGTGTTGAAGTTCGGCGGTACATCAGGGCAAATGCAGAACGTTTTCTGCGGGTT
Protein = MRVLKFGGTSVANAERFLRV
Converted protein protein1 to resemble dna dna1
Original protein:MRVLKFGGTSVANAERFLRV
Protein to DNA: ATGCGAGTGTTGAAGTTCGGCGGTACATCANNNGCAAATGCAGAACGTTTTCTGCGGNNN
Original DNA: ATGCGAGTGTTGAAGTTCGGCGGTACATCAGGGCAAATGCAGAACGTTTTCTGCGGGTT
Trying to align protein1.dna[0:60] = ATGCGAGTGTTGAAGTTCGGCGGTACATCANNNGCAAATGCAGAACGTTTTCTGCGGNNN
Translated protein[0:20]: MRVLKFGGTS?ANAERFLR?
alignment at dna1 offset 0:
TGCGAGTGTTGAAGTTCGGCGGTACATCANNNGCAAATGCAGAACGTTTTCTGCGGNNN
TGCGAGTGTTGAAGTTCGGCGGTACATC-AGGGCAAATGCAGAACGTTTTCTGCGGGTT

Summary for protein1 compared to dna1:
Insertion of length 1bp at 29
```

To get it, do one of these:

* Download [a release](https://github.com/mathjeff/ProteinFrameshiftDetector/releases/) version for example [this one](https://github.com/mathjeff/ProteinFrameshiftDetector/releases/download/0.0.1/ProteinFrameshiftDetector.jar)

* Build it yourself

```
./gradlew shadowJar
java -jar build/libs/ProteinFrameshiftDetector.jar --dna <DNASequence> --protein <ProteinSequence>
```
