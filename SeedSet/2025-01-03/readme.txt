Change history:

2025-01-03
- Numerated the allele families initially based on IMGT, and KIMD.

2024-11-15
- Applied the new allele cluster similarity method to cluster the sequences.

2024-06-24
- Removed sequences that were not part of ony of the core reference sets, and also were not found in the WGS trios. These sequences
  entered the set during early analysis of tje trios, but now do not have any evidential foundation.

2024-05-10
- Revised the V gapping, based on a new manually-derived aa alignment

2024-05-06
- Updated IGK and IGL with changes made after WGS alignment of 7 children from trios
- Fixed a bug in the merge code which caused some alleles to be assigned to the wrong clusters. All clusters are now correct
  at the 905/95% similarity threshold. Across the sets, around 10% of alleles were reassigned as a result of this fix.
- Dropped some questionable sequences from the IGHV set following review by Mats Ohlin dicussed at the macaque WG in April 2024.
  (thse are sequences from the AIRR-C GLDB rhesus database which have not been observed in the WGS trios)

2024-03-20
- added sequence IGHV-WGXO (kimdb: IGHV3-5*01_S6059) - previously marekd as 'do not include' but there did not seem to be a justification for not including it

2024-03-05
- Updated IGH V and D following review of the WGS alignment

2024-02-28
- Update Ds with changes made after WGS alignment of 7 children from trios 

2024-02-19
- Re-cluster Ds at 90% similarity threshold

2024-02-18
- Fix issues with four PSVV IGHDD alleles that previously included the heptamer fragment CAGTGT at the 3' end

2024-01-15
- Removed 17 sequences from IGH D.fasta that were sub-sequences of other sequences in the file. No sequence names changed. No changes to other sets.

2023-10-29
- Updated to include V and J sequences from seven genomic mother/father/child trios (21 animals in all)

2023-09-14
- Created from all validated sequences in the AIRR-C GLDB rhesus database at 2023-04-12


Using the germline sets with IgBlast:

Please follow the procedure outlined here: https://ncbi.github.io/igblast/cook/How-to-set-up.html, using makeblastdb to generate the blast databases and 
then following the procedures to set up for a custom organism. 'aux' and 'ndm' files are provided for each locus. 