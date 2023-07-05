.. _biotypes:

=========
Biotypes
=========


BamQuery biotype classification is based on `Ensembl`_ and `Repeat_Masker`_ annotations.


.. _Ensembl: https://m.ensembl.org/info/genome/genebuild/biotypes.html


.. _Repeat_Masker: https://www.repeatmasker.org/



1. **Protein-coding Regions: peptide harbored in a transcript containing an open reading frame (ORF)**

	In this category the following transcripts types in gencode annotations are considered as **Protein-coding Regions**: |br|
	Protein_coding, IG genes and TR genes.

    * **5'UTR**: peptide harbored into the 5' UnTranslated Region of the transcript

    * **3'UTR**: peptide harbored into the 3' UnTranslated Region of the transcript

    * **In_frame**: peptide harbored in the CDS of the transcript, and for which at least 50% of its sequence is In Frame with the protein

    * **Frameshift**: peptide harbored in the CDS of the transcript, and for which at least 50% of its sequence is Frameshift with the protein

    * **CDS**: Peptide containing a SNP harbored in the CDS of the transcript. Peptide that is not In_Frame or Frameshift with the protein

    * **Junctions**: peptide harbored in the junctions (Intron-CDS, CDS-3'UTR, 5'UTR-CDS) of the transcript

    * **Other coding regions**: peptide harbored in exons into Immunoglobulin gene (IG_C_gene, IG_D_gene, IG_J_gene, IG_V_gene) or T cell receptor genes (TR_C_gene, TR_D_gene, TR_J_gene, TR_V_gene) 


2. **Non-coding RNAs: peptide harbored in a transcript that doesn't containing an open reading frame (ORF)**	

	In this category the following transcripts types in gencode annotations are considered as **Non-coding RNAs**: |br|
	IG_C_pseudogene, IG_J_pseudogene, IG_V_pseudogene, IG_pseudogene, Mt_rRNA, Mt_tRNA, TEC, TR_J_pseudogene, TR_V_pseudogene, lncRNA, miRNA, misc_RNA, non_stop_decay, nonsense_mediated_decay, polymorphic_pseudogene, processed_pseudogene, processed_transcript, pseudogene, rRNA, rRNA_pseudogene, retained_intron, ribozyme, sRNA, scRNA, scaRNA, snRNA, snoRNA, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, translated_unprocessed_pseudogene, unitary_pseudogene, unprocessed_pseudogene, vault_RNA.

	* **Non_coding Exons**: peptide harbored in the exons of the transcript

	* **Non_coding Junctions**: peptide harbored in the junctions (Exon-Intron, Intron-Exon) of the transcript


3. **Intergenic Regions: peptide harbored in a non-annotated region**

	* **Intergenic Regions**: peptide harbored in a non-annotated region either in Ensembl or in Repeat Masker


4. **Intronic Regions: peptide harbored in a intronic region of any type transcripts**

	* **Intronic Regions**: peptide harbored in an intronic region of any type of transcript, whether protein-coding or non-coding RNA

5. EREs: peptide harbored in a ERE (based on `Repeat Masker <S2>`_ annotations)
	* **LINE**: peptide harbored in a ERE class LINE

	* **LTR**: peptide harbored in a ERE class LTR

	* **SINE**: peptide harbored in a ERE class SINE

	* **Antisense_EREs**: peptide harbored in the antisense of a ERE of any class

	* **Other EREs**: peptide harbored in a ERE of other class (DNA, RC, RNA, Satellite, Simple_repeat, Unknown, Retroposon, Low_complexity, rRNA, scRNA, snRNA, srpRNA, tRNA)


.. |br| raw:: html

      <br>


