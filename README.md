ProteinFrameshiftDetector is a prototype project for finding a single small frameshift that explains the difference between a DNA sequence and a protein sequence

To use it, do this:

```
java -jar ProteinFrameshiftDetector.jar <DNASequence> <ProteinSequence>
```

For example:

```
$ java -jar build/libs/ProteinFrameshiftDetector.jar ATGCGAGTGTGGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAAGTTTTCTGCGGGTT MRVLKFGGTSVANAERFLRV

This program isn't optimized but hopefully is helpful


Most similar result is a frameshift of length 1 at 43 with similarity of about 0.6666666666666666
Original DNA: ATGCGAGTGTGGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAAGTTTTCTGCGGGTT
Shifted  DNA: ATGCGAGTGTGGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAAAGTTTTCTGCGGGTT
Translated  : MRVWKFGGTSVANAESFLRV
Protein     : MRVLKFGGTSVANAERFLRV
```

To get it, do one of these:

* Download [a release](https://github.com/mathjeff/ProteinFrameshiftDetector/releases/) version for example [this one](https://github.com/mathjeff/ProteinFrameshiftDetector/releases/download/0.0.1/ProteinFrameshiftDetector.jar)

* Build it yourself

```
./gradlew build
java -jar build/libs/ProteinFrameshiftDetector.jar <DNASequence> <ProteinSequence>
```
