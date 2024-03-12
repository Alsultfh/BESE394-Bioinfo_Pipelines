
# W5 Multiome
### Assignment 3:
Generate motifs for the top 5 differentially accessible peaks
per cell type. From JASPAR, describe their function, and
compare them between cell types.Are the motifs, revealing the underlying biology of the
different cell types?

## Workflow

####  1. Load `filtered_feature_bc_matrix` contains 
- scRNA data Gene Expression -> `🔵RNA counts`
- scATAC data peaks info -> `🟡ATAC data`

#### 2. Create seurat object with both assays

Format correctly the `🟡ATAC data`. Rownames in atac counts look like: chrX:162621438-162622286  
Need to create GRange for each row in counts ->   [1]     chr1 3093782-3094521      *
```
An object of class Seurat 
113165 features across 23990 samples within 2 assays 
Active assay: RNA (32285 features, 0 variable features)
 1 layer present: counts
 1 other assay present: ATAC
```

#### 3. Add info to the metadata

    For `🟡ATAC data` quality control:
    
    1. `NucleosomeSignal(mouse_brain)`: nucleosome_signal,nucleosome_percentile, nucleosome_group
    2. `TSSEnrichment`TSS.enrichment, TSS.percentile, high.tss
    3. `CountFragments`fragments
    4. `FRiP` FRiP


    For `🔵RNA counts` quality control:

    1. mitoPct
    2. RPSPct
    3. log10GenesPerUMI

#### 4. Visualize the data without any filtering
![picture alt](./content/imag/1stqc.png)

#### 5. Filter

- `nCount_ATAC` 1000 - 5*std(counts)=22594.5
- `nCount_RNA`  1000 - 5*std(counts)=35075.9
- `nFeature_ATAC`  < 5*std(featue)=10253.58
- `nFeature_RNA`  < 5*std(featue)= 7761.97
- `mitoPct` < 15


```
 mouse_brain <- subset(
  x = mouse_brain,
  subset = nCount_ATAC < 22594.5 &
    nCount_RNA < 35075.9 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nFeature_ATAC < 10253.58 &
    nFeature_RNA < 7761.97 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 & 
    FRiP > 0.2 &
    mitoPct < 15
)

```
![picture alt](./content/imag/2ndqc.png)

```
An object of class Seurat 
113165 features across 18099 samples within 2 assays 
Active assay: RNA (32285 features, 0 variable features)
 1 layer present: counts
 1 other assay present: ATAC
 ```

### 6. Process the find Clusters
For `🔵RNA counts`

1. Normalization 
2. Find Var Features `selection.method = "vst", nfeatures = 2000)`
3. Scaling
4. PCA 
5. FindNeigh FindClust `dims = 1:30`
6. FindClust `resolution = 0.4, algorithm = 3, cluster.name="RNA_clusters_03"`
7. UMAP `reduction.name = "rna_umap"`
7. Azimuth `reference = "mousecortexref"`

For `🟡ATAC data`:
1. FindTopFeatures `min.cutoff = 5`
2. TFIDF
3. SVD
4. UMAP `reduction = "lsi", reduction.name = "atac_umap"`

![picture alt](./content/imag/RNAATAC.png)

### 7. Integrate 

Integrate using pca(`🔵RNA`) and lsi (`🟡ATAC`) reductions

```mouse_brain <- FindMultiModalNeighbors(
  object = mouse_brain,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 1:40),
  verbose = TRUE
)
```
![picture alt](./content/imag/MULTIOME.png)


### 8. Differentially Accessible Peaks

We choose 1 inhibitory neuron LAMP5 and compared it against the 3 other subtypes of inhibitory neurons: SST (sst and sst chodl), PVALB, and VIP

```
FindMarkers(
  object = mouse_brain,
  ident.1 = c("Lamp5"), 
  ident.2 = c("Sst"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
```





### 9. Motif Analysis

![picture alt](./content/imag/MOTIFS.png)

### Biological Role of the Motifs
|            | **LAMP5 vs PVALB**                                                                                                                                            |
|------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **EGR1**   | The products of target genes it activates are required for differentitation and mitogenesis.                                                                  |
| **Wt1**    | Regulates the expression of numerous target genes, including EPO. Plays an essential role for development of the urogenital system. I                         |
| **ZBTB14** | Binds to trinucleotide repeats in promoter regions and acts as a repressor of the FMR1 gene. Transcriptional repressor of MYC and thymidine kinase promoters. |
| **EGR3**   | Transcription factor involved in muscle spindle development                                                                                                   |
| **KLF15**  | Transcriptional regulator that binds to the GA element of the CLCNKA promoter.                                                                                |
| **E2F6**   | Regulate a subset of E2F-dependent genes whose products are required for entry into the cell cycle but not for normal cell cycle progression                  |

---

|                 | **LAMP5 vs SST**                                                                                                                                                                                                                                                                                          |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **KLF15**       | Transcriptional regulator that binds to the GA element of the CLCNKA promoter.                                                                                                                                                                                                                            |
| **WT1**         | Regulates the expression of numerous target genes, including EPO. Plays an essential role for development of the urogenital system.                                                                                                                                                                      |
| **EGR1**        | The products of target genes it activates are required for differentitation and mitogenesis.                                                                                                                                                                                                              |
| **EGR3**        | Transcription factor involved in muscle spindle development                                                                                                                                                                                                                                               |
| **BHLHE22**     | As an area-specific transcription factor that regulates the postmitotic acquisition of area identities and elucidate the genetic hierarchy between progenitors and postmitotic neurons driving neocortical. May be required for the survival of a specific population of inhibitory neurons arealization. |
| **FOSL1::JUND** | DNA-binding transcription factor activity and RNA polymerase II transcription regulatory region sequence-specific DNA binding                                                                                                                                                                             |


---

|           | **LAMP5 vs  SST CHODL**                                                                                                                     |
|-----------|---------------------------------------------------------------------------------------------------------------------------------------------|
| **KLF15** | Transcriptional regulator that binds to the GA element of the CLCNKA promoter.                                                              |
| **EGR1**  | The products of target genes it activates are required for differentitation and mitogenesis.                                                |
| **EGR3**  | Transcription factor involved in muscle spindle development                                                                                 |
| **HES1**  | negative regulator of myogenesis by inhibiting the functions of MYOD1 and ASH1                                                              |
| **SP4**   | bind to the GC promoter region of a variety of genes, including those of the photoreceptor signal transduction system                       |
| **KLF14** | transcriptional co-repressor, and is induced by transforming growth factor-beta (TGF-beta) to repress TGF-beta receptor II gene expression. |

---

|             | **LAMP5 vs VIP**                                                                                                                                               |
|-------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **WT1**     | Regulates the expression of numerous target genes, including EPO. Plays an essential role for development of the urogenital system. I                          |
| **KLF15**   | Transcriptional regulator that binds to the GA element of the CLCNKA promoter.                                                                                 |
| **TCF12**   | Transcriptional regulator. Involved in the initiation of neuronal differentiation.                                                                             |
| **E2F6**    | Regulate a subset of E2F-dependent genes whose products are required for entry into the cell cycle but not for normal cell cycle progression                   |
| **BHLHE22** | area-specific transcription factor that regulates the postmitotic acquisition. May be required for the survival of a specific population of inhibitory neurons |
| **NHLH1**   | DNA-binding protein and may be involved in the control of cell-type determination, possibly within the developing nervous system.                              |

---


Most of the motifs identified are not revealing the underlying biology of the cell type. However is important to note the comparisson we are doing is from LAMP5 against (SST, PVALB, and VIP) maybe if we consider all other cell types identified (exhitatory and non-neuronal cells) we would get more menaingful results. Adittionally we only explored 6 motifs for each of the comparissons.

One motif worth mentioning is **BHLHE22** obtained from the comparisson of LAMP5 against SST and VIP. BHLHE22 is an area-specific transcription factor that regulates the postmitotic acquisition of area identities and elucidate the genetic hierarchy between progenitors and postmitotic neurons driving neocortical. May be required for the survival of a specific population of inhibitory neurons.
             |

## References
The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses (PMID: 27322403; Citations: 2,841) 
Stelzer G, Rosen R, Plaschkes I, Zimmerman S, Twik M, Fishilevich S, Iny Stein T, Nudel R, Lieder I, Mazor Y, Kaplan S, Dahary, D, Warshawsky D, Guan - Golan Y, Kohn A, Rappaport N, Safran M, and Lancet D
Current Protocols in Bioinformatics(2016), 54:1.30.1 - 1.30.33.doi: 10.1002 / cpbi.5 


## Extra

#### If we run with container

```
cd /ibex/scratch/gonzalac/BESE394a/W5_multiome
module load singularity 
singularity run bese_multiomeanalysis_latest.sif
```

 