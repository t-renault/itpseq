.. _whatis_itpseq:

================
What is iTP-Seq?
================

Uneven translation rates resulting from mRNA context, tRNA abundance, nascent amino acid sequence or various external factors play a key role in controlling the expression level, regulation and folding of the proteome.

Inverse toeprinting coupled to next-generation sequencing (iTP-Seq) [1]_ is a scalable *in vitro* method for characterizing bacterial translation landscapes, complementary to ribosome profiling (Ribo-Seq) [2]_ [3]_, a widely used method for determining transcriptome-wide protein synthesis rates *in vivo*.

In iTP-Seq, ribosome-protected mRNA fragments known as inverse toeprints are generated using RNase R, a highly processive 3' to 5' RNA exonuclease. Deep sequencing of these fragments reveals the position of the leading ribosome on each mRNA with codon resolution, as well as the full coding regions translated by these ribosomes.

.. figure:: /_static/whatis_itpseq_overview.png
  :width: 400
  :alt: iTP-Seq methodology overview
  :align: center

  **Figure 1.** Overview of the iTP-Seq methodology.

Unlike Ribo-Seq, in which short mRNA footprints are mapped to a sequenced genome, no *a priori* knowledge of the translated sequences is required for iTP-Seq, making it possible to work with fully customizable transcript libraries. As a standardized framework for inverse toeprint generation, amplification and sequencing, iTP-Seq can therefore be used in combination with different types of libraries, *in vitro* translation conditions and data analysis pipelines tailored to address a range of biological questions involving context-dependent translational arrest, as exemplified by our previous studies on ribosome-targeting antibiotics [4]_ [5]_.

A detailed experimental protocol for performing iTP-Seq can be found here and is outlined below.

.. figure:: /_static/whatis_itpseq_workflow.png
  :width: 500
  :alt: iTP-Seq experimental workflow
  :align: center

  **Figure 2.** Experimental workflow for iTP-Seq analysis.

Bioinformatic analysis of the resulting iTP-seq data can be performed using the ``itpseq`` software available on this website, as explained in our various :ref:`tutorials <tutorials>`.

References:

.. [1] Seip, B., Sacheau, G., Dupuy, D. & Innis, C. A. Ribosomal stalling landscapes revealed by high-throughput inverse toeprinting of mRNA libraries. Life Sci. Alliance 1, e201800148 (2018). `<https://doi.org/10.26508/lsa.201800148>`_
.. [2] Ingolia, N. T., Ghaemmaghami, S., Newman, J. R. S. & Weissman, J. S. Genome-Wide Analysis in Vivo of Translation with Nucleotide Resolution Using Ribosome Profiling. Science 324, 218–223 (2009). `<https://doi.org/10.1126/science.1168978>`_
.. [3] Mohammad, F., Green, R. & Buskirk, A. R. A systematically-revised ribosome profiling method for bacteria reveals pauses at single-codon resolution. eLife 8, e42591 (2019). `<https://doi.org/10.7554/elife.42591>`_
.. [4] Leroy, E. C., Perry, T. N., Renault, T. T. & Innis, C. A. Tetracenomycin X sequesters peptidyl-tRNA during translation of QK motifs. Nat. Chem. Biol. 19, 1091–1096 (2023). `<https://doi.org/10.1038/s41589-023-01343-0>`_
.. [5] Beckert, B. et al. Structural and mechanistic basis for translation inhibition by macrolide and ketolide antibiotics. Nat. Commun. 12, 4466 (2021). `<https://doi.org/10.1038/s41467-021-24674-9>`_
