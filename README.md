# Blockbuster

Once short read sequences are mapped to a reference genome, one will face the problem of dividing consecutive reads into blocks to detect specific expression patterns. Due to biological variability and sequencing inaccuracies, the read arrangement does not always show exact block boundaries. The blockbuster tool automatically assigns reads to blocks and gives a unique chance to actually see the different origins where the short reads come from.

The original program written in C, incorrectly manages memory allocations hence causing memory leaks. I fixed this bug in a way that the program writes only into a specific short part of the memory. Not only the memory leak was removed, but also the runtime performance increased 4-5 times compared to the original program. I also added a Python implementation, but it is much more slower than the compiled C code as expected. Open the link below for the usage.

Source: http://hoffmann.bioinf.uni-leipzig.de/LIFE/blockbuster.html
Langenberger D, Bermudez-Santana C, Hertel J, Hoffmann S, Khaitovitch P, Stadler PF: "Evidence for Human microRNA-Offset RNAs in Small RNA Sequencing Data", Bioinformatics (2009) vol. 25 (18) pp. 2298-301
