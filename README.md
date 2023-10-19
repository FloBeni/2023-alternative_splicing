[![DOI](https://zenodo.org/badge/575819362.svg)](https://zenodo.org/doi/10.5281/zenodo.7415114)

This archive includes all processed data that have been generated and used in the study 'Random genetic drift sets an upper limit on mRNA splicing accuracy in metazoans', as well as the scripts used to analyze the data and to generate the figures.

Developed by Florian Bénitière, Anamaria Necsulea and Laurent Duret. Université de Lyon, Université Lyon 1, CNRS, Laboratoire de Biométrie et Biologie Évolutive UMR 5558, F-69622 Villeurbanne, France.

### Directories content

-   The folder 'figure' contains the necessary materials to produce the figures of the manuscript:

    -   Intermediate pannels generated by R scripts to produce figures are stored in the directory 'pannels'.

    -   Stable images used in pannels are stored in 'images_library'.

    -   Directories labeled with the extension '\_generator' contain R scripts used to produce figures.

    -   The remaining directories are where the figures are stored.

-   Rmarkdown scripts located in the 'table_supp' directory are responsible for generating supplementary tables, which are also saved in the same directory.

-   The 'data' directory contains processed data in tab-separated text format, which is used for generating figures and conducting analyses.

    -   The 'per_species' directory contains two compressed tab-separated text format files for each of the 69 studied species:

        '*by\_*gene_analysis.tab' which contain per-gene data.

        <div>

        -   **gene_id**: ID of the gene found in the GFF.

        -   **gene_name**: Corresponds to the attribute 'gene=' in the GFF.

        -   **seq_id**: Refers to seqname in the GFF, corresponding to the name of the chromosome or scaffold.

        -   **start**: Start position.

        -   **end**: End position.

        -   **strand**: Defined as + (forward) or - (reverse).

        -   **type**: Feature type name, e.g. Gene, Variation, Similarity found in the GFF.

        -   **attributes**: Semicolon-separated list of tag-value pairs, providing additional information about each feature, extracted from the GFF.

        -   **weighted_fpkm**: Computed as the average FPKM across samples, weighted by the sequencing depth of each sample. The sequencing depth of a sample is the median *per*-base read coverage across BUSCO genes.

        -   **median_fpkm**: Median FPKM across all RNA-seq samples.

        -   **mean_fpkm**: Average FPKM across all RNA-seq samples.

        -   **std_fpkm**: Standard deviation FPKM across all RNA-seq samples.

        -   **exon_coverage**: Exonic read *per* bp (measured on all RNA-seq pooled).

        -   **busco_metazoa**: If the gene can be associated to a BUSCO genes without ambiguity (*i.e.* BUSCO genes were removed from the analysis if they were associated to more than one annotated gene or to an annotated gene that was associated to more than one BUSCO gene).

        </div>

        '*by\_*intron_analysis.tab' which contain per-intron data.

        <div>

        -   **gene_id**: ID of the gene found in the GFF.

        -   **seqname**: seqname in the GFF, corresponding to the name of the chromosome or scaffold

        -   **strand**: Defined as 1 (forward) or -1 (reverse).

        -   **splice5**: Corresponds to the 5' splice donor site.

        -   **splice3**: Corresponds to the 3' splice acceptor site.

        -   **n1**: Number of spliced reads corresponding to the precise excision of the focal intron ($N_s$ in the paper).

        -   **n2_spl5**: Number of reads corresponding to alternative splice variants relative to this intron (*i.e.* sharing 5' splice donor site) ($\mathrm{n2spl5+n2spl3=N_a}$ in the paper).

        -   **n2_spl3**: Number of reads corresponding to alternative splice variants relative to this intron (*i.e.* sharing 3' splice acceptor site) ($\mathrm{n2spl5+n2spl3=N_a}$ in the paper).

        -   **n3_spl3**: Number of unspliced reads, co-linear with the genomic sequence (at the 3' splice acceptor site) ($\mathrm{n3spl5+n3spl3=N_u}$ in the paper).

        -   **n3_spl5**: Number of unspliced reads, co-linear with the genomic sequence (at the 5' splice donor site). ($\mathrm{n3spl5+n3spl3=N_u}$ in the paper).

        -   **splice_variant_rate**: Proportion of reads alternatively spliced, $\mathrm{1-RAS=AS=\frac{N_a}{N_s~+~N_a}}$.

        -   **nonsplice_variant_rate**:Proportion of unspliced reads, $\mathrm{1-RANS=\frac{N_u}{2\times N_s~+~N_u}}$.

        -   **intron_class**: Three categories of introns: major introns, defined as those introns that have RANS $>$ 0.5 and RAS $>$ 0.5; minor introns, defined as those introns that have RANS $\leq$ 0.5 or RAS $\leq$ 0.5; unclassified introns, which do not satisfy the above conditions.

        -   **into_cds**: Check if intron is located within protein-coding regions. To do this, for each protein-coding gene, we extracted the start codons and the stop codons for all annotated isoforms. We then identified the minimum start codon and the maximum end codon positions and we excluded introns that were upstream or downstream of these extreme coordinates.

        -   **id**: Semicolon-separated list of tag-value seqname;gene_id;splice3;end:strand

        -   **fpkm**: Weighted_fpkm of the corresponding gene. Computed as the average FPKM across samples, weighted by the sequencing depth of each sample. The sequencing depth of a sample is the median *per*-base read coverage across BUSCO genes.

        -   **splicesite**: Identify the dinucletotides donor and acceptor. (ex: GT AG)

        -   **phase**: Extract the phase of the intron from the GTF file.

        -   **Annotation**: Check if the intron is annotated (*i.e.* present in the GFF)or not.

        -   **mira criptic_intron**: For minor introns sharing a boundary with a major intron: $\mathrm{Minor~intron~relative~abundance~MIRA_i=\frac{N_i^m}{N^M~+~N^m}}$

        -   **distance_from_major**: In bp the distance from the sharing major intron boundaries.

        -   **frame_shift**: The residual value resulting from the division of distance_from_major by 3.

        -   **have_abundant_sv**: For major intron, check if they have an habundant minor variants (*i.e.* MIRA \> 0.05*)*.

        </div>

    -   'Data1_supp.tab' contains summarize data and information per species such as the genome assembly, the clade, the list of the RNA-seq studied... Used in Fig. 1,3,4,6 and Supplementary Fig. 2,3,6,8,9.

    -   'Data2_supp.tab' contains the analyzes of the relation between the average MIRA and the frame preserving proportion (Fig. 4A).

    -   'Data3_supp.tab' contains the distribution of the relative abundance of the spliced isoform compared to other transcripts with alternative splice boundaries (RAS) or compared to unspliced transcripts (RANS) (divided into 5% bins), for protein-coding gene introns (Fig. 2).

    -   'Data4_supp.tab' is the analysis of the average AS rate between organs in different tissues (Fig. 3 and Supplementary Fig. 9).

    -   'Data5_supp.tab' is the analysis of the SNPs density in *Drosophila melanogaster* and *Homo sapiens* (Fig. 5 and Supplementary Fig. 4).

    -   'Data6_supp.tab' is the analysis of the relation between the FPKM and the average AS rate (Fig. 6 and Supplementary Fig. 5).

    -   'Data7_supp.tab' analysis of the coverage impacting the AS rate and the proportion of analyzable annotated introns (Supplemenary Fig. 1).

    -   'Data8_supp.tab' Proportion of frame preserving introns for every species in (Fig. 4 and Supplementary Fig. 7).

    -   'Data9_supp.pdf' contains sources and information regarding the lifespan and longevity of the species studied.

    -   'Data10_supp.tab' Table containing all the RNA-seq studied with information extracted from SRA runinfo table.

    -   'Data11_supp.tab' contains the proportion of major introns per reading frame class per species.

-   The 'pipelines' folder contains the bionformatics pipelines for three different purposes: to calculate the dN/dS ratio ('dNdS_pipeline'); to analyze various alternative splicing characteristics as part of our study ('AS_pipeline'); to generate data table located in the 'data' directory ('data_generator').

### 

Some packages that you may encounter in these R scripts should be install: - stringr - kableExtra - cowplot - ggplot2 - imager - RColorBrewer ...

### 

The script 'values_in_text_paper_generator.R' retrieves the values found in the main text of the paper from the data.
