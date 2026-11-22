# Deteksi Kanker Kolorektal dari Sekuens 16S rRNA (FASTQ) dengan Machine Learning

Pipeline Python yang mengolah **584 berkas FASTQ mentah** (16S rRNA region V4, Illumina MiSeq) dari sampel feses **292 orang** (120 kanker kolorektal, 172 normal) menjadi fitur machine learning, lalu membandingkan **Logistic Regression** dan **Random Forest** pada tiga representasi: **ESV**, **Genus**, dan **k-mer**.

- **Hasil utama:** ROC-AUC nested 5-fold CV **0,77** (Random Forest + k-mer) dan **0,765** (Random Forest + genus), jauh di atas baseline 0,50.
- **Bakteri penanda:** *Porphyromonas*, *Peptostreptococcus*, *Parvimonas*, *Fusobacterium*, dan *Gemella* lebih tinggi pada pasien kanker (signifikan setelah FDR).
- **Validasi:** hasil olahan dibandingkan dengan tabel OTU penulis paper (Shannon ρ = 0,81).

📄 Penjelasan lengkap (metode, hasil, keterbatasan): **[RINGKASAN_PROYEK.md](RINGKASAN_PROYEK.md)** · Notebook: **[Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb](Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb)** · Alur: [flowchart_alur_penelitian.png](flowchart_alur_penelitian.png)

## Data
| Data | Sumber | Di repo? |
|---|---|---|
| FASTQ mentah (584 berkas, 4,5 GB) | ENA/SRA [PRJNA290926](https://www.ebi.ac.uk/ena/browser/view/PRJNA290926), Baxter et al. (2016) | ❌ terlalu besar. Notebook **mengunduhnya otomatis** dari ENA dan mengecek MD5 |
| Metadata + label diagnosis | ENA Portal API + NCBI BioSample | ✅ `data/metadata_sampel.tsv` |
| Referensi 16S | [NCBI RefSeq 16S rRNA](https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/) | ✅ `data/reference/bacteria.16SrRNA.fna.gz` |
| Tabel OTU penulis (hanya untuk validasi) | [SchlossLab/Baxter_glne007Modeling_GenomeMed_2015](https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015) | ✅ `data/reference/penulis_otu_0.03.shared` |

## Cara menjalankan
```bash
pip install numpy pandas scipy scikit-learn statsmodels matplotlib seaborn joblib jupyter
jupyter notebook Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb   # lalu Run All
```
Saat pertama dijalankan, notebook mengunduh FASTQ ±4,5 GB ke `data/fastq/`. Lama jalan ±16 menit setelah data ada, dengan seed 42. Di Google Colab, aktifkan dua baris *mount Drive* di bagian 2.

## Isi repo
- `hasil_proyek2/results/`: tabel CSV, FASTA ESV, dan grafik 01–11
- `hasil_proyek2/models/`: model `.joblib` + skema fitur
- `Microbiota-based model improves the.pdf`: paper sumber data (open access, CC BY 4.0)

## Sitasi data
Baxter NT, Ruffin MT, Rogers MAM, Schloss PD (2016). *Microbiota-based model improves the sensitivity of fecal immunochemical test for detecting colonic lesions.* Genome Medicine 8:37. https://doi.org/10.1186/s13073-016-0290-3

> Proyek akademik, **bukan alat diagnosis klinis**.
