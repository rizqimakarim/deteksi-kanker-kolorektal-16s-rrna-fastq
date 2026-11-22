# Ringkasan Proyek 2: Deteksi Kanker Kolorektal dari Sekuens 16S rRNA (FASTQ) dengan Machine Learning

**Judul:** Deteksi Kanker Kolorektal Berdasarkan Data Sekuens 16S rRNA (FASTQ) Mikrobioma Feses Menggunakan Logistic Regression dan Random Forest

**Tema tugas:** Pengembangan Model Machine Learning untuk Deteksi Penyakit Berdasarkan Data Sekuens DNA/RNA

**Notebook:** [Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb](Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb) · **Flowchart:** [flowchart_alur_penelitian.png](flowchart_alur_penelitian.png)

---

## Hasil utama secara singkat

- **Data:** 584 berkas **FASTQ mentah** (R1 + R2, 4,48 GB) dari 292 orang: 120 kanker kolorektal dan 172 normal. Seluruh fitur dihitung langsung dari read sekuensing di notebook.
- **Performa model** (nested 5-fold CV, ROC-AUC):
  - Random Forest + fitur genus: **0,765**
  - Random Forest + fitur k-mer: **0,770**
  - Logistic Regression + fitur k-mer: **0,751**
  - Semua jauh di atas baseline acak (0,50).
  - Performa ini **sedang**, bukan sempurna. Hasil ini sejalan dengan literatur: mikrobioma saja belum cukup untuk diagnosis.
- **Bakteri penanda yang ditemukan:** *Porphyromonas*, *Peptostreptococcus*, *Parvimonas*, *Fusobacterium*, dan *Gemella* lebih banyak dan lebih sering ditemukan pada pasien kanker.
  - Kelima genus ini signifikan setelah koreksi FDR.
  - Genus yang sama juga menjadi kontributor terbesar model.
  - Temuan ini konsisten dengan paper sumber data dan studi independen.
- **Validasi:** hasil olahan FASTQ di notebook dibandingkan dengan tabel OTU yang dipublikasikan penulis paper. Keragaman Shannon berkorelasi **ρ = 0,81**, dan taksa dominan berkorelasi dengan median **0,88**.

---

## 1. Latar belakang dan tujuan

Kanker kolorektal termasuk kanker paling umum di dunia. Deteksi dini umumnya memerlukan kolonoskopi, yang invasif dan mahal. Komposisi bakteri usus (mikrobioma) diketahui berubah pada pasien kanker kolorektal. Karena itu, sekuensing gen 16S rRNA dari sampel feses berpotensi menjadi alat skrining **non-invasif**.

Tujuan proyek:
1. Membangun pipeline yang mengolah **FASTQ mentah** menjadi fitur machine learning.
2. Membandingkan **Logistic Regression** dan **Random Forest** untuk mendeteksi kanker kolorektal.
3. Membandingkan tiga representasi sekuens: **ESV**, **Genus**, dan **k-mer**.
4. Mengidentifikasi bakteri yang paling membedakan pasien kanker dari orang normal.
5. Mengevaluasi model secara jujur: tanpa kebocoran data, dengan validasi silang, dan dengan analisis sensitivitas.

## 2. Dataset

| Item | Keterangan |
|---|---|
| Sumber | ENA/SRA BioProject [PRJNA290926](https://www.ebi.ac.uk/ena/browser/view/PRJNA290926), Baxter et al. (2016), *Genome Medicine* |
| Jenis data | Amplikon 16S rRNA region V4, Illumina MiSeq 2 × 250 bp |
| Format | **FASTQ** terkompresi gzip (`*_1.fastq.gz` = R1, `*_2.fastq.gz` = R2) |
| Ukuran | 584 berkas, 4,48 GB, 16,1 juta pasang read, 8,09 Gbp |
| Label | Atribut `diagnosis` dari BioSample NCBI (dicek cocok 100% dengan metadata penulis) |
| Kohor | 292 orang, satu sampel feses per orang |
| Kelas | **cancer**: 120 · **normal**: 172 (122 Normal + 50 High Risk Normal = kolonoskopi normal) |
| Dikeluarkan | Adenoma (di luar fokus kanker vs normal) dan sampel *mock community* (kontrol teknis) |
| Referensi taksonomi | NCBI RefSeq 16S rRNA (format **FASTA**, 26.888 sekuens) |

Contoh satu record FASTQ dari dataset:

```
@SRR2143516.1 1/1                       ← header (ID read)
TACGTAGGGAGCAAGCGTTGTCCGGATTTACTGGG...    ← sekuens basa (250 bp)
+                                         ← pemisah
AAAAACB>3B22AFBEEAE2FFFGGGEAF55FGEH...    ← skor kualitas Phred per basa
```

---

## 3. Alur analisis

```
Unduh 584 FASTQ utuh dari ENA + verifikasi MD5
        ↓
Metadata & label diagnosis (BioSample) → kohor 292 orang
        ↓
QC read mentah R1/R2 (Phred per posisi, Q30, GC, N)
        ↓
Filter kualitas (200 bp, maxEE ≤ 2, tanpa N) → rarefaction 10.000 read
        ↓
Fitur: ESV · k-mer (k=6) · Genus (anotasi FASTA NCBI 16S)  + validasi silang data penulis
        ↓
Eksplorasi: genus terpilih → PCA → elbow/K-Means → Mann–Whitney + FDR → volcano & heatmap
        ↓
ML: split 80/20 → LR vs RF × 3 representasi → evaluasi test set + nested 5-fold CV
        ↓
Analisis sensitivitas (teknis, batch, usia) → kontribusi fitur → ekspor
```

### 3.1 Pengambilan data (bagian 2–4 notebook)
- Metadata run diambil dari ENA Portal API, lalu label `diagnosis` dari BioSample NCBI.
- Ada 49 sampel yang diurutkan lebih dari sekali. Untuk sampel seperti ini dipilih satu run dengan read terbanyak, sehingga setiap orang tetap diwakili satu sampel.
- Seluruh 584 berkas FASTQ diunduh utuh dan diverifikasi dengan **MD5 resmi ENA**. Berkas yang sudah ada dan cocok MD5-nya dipakai ulang; unduhan yang terputus dilanjutkan dari posisi terakhir.
- Catatan metadata: kolom instrumen di ENA sebagian tertulis "454 GS". Kolom ini salah tulis, karena format read dan skor kualitas jelas Illumina MiSeq.

### 3.2 QC read mentah (bagian 6)
Keunggulan data FASTQ dibanding tabel ekspresi olahan adalah **QC tingkat read** bisa dilakukan langsung.

| Metrik | Cancer | Normal |
|---|---|---|
| Rerata read R1 per sampel | 56.290 | 54.409 |
| Rerata skor Phred R1 / R2 | 34,8 / 32,9 | 34,6 / 32,7 |
| Basa ≥ Q30, R1 / R2 | 89,0% / 80,0% | 88,0% / 78,9% |
| GC content R1 | 52,5% | 52,3% |
| Basa N | ≈ 0% | ≈ 0% |
| Read lolos filter | 96,4% | 96,9% |

- Kualitas R2 lebih rendah dari R1, terutama di ujung read. Karena itu analisis memakai **R1 yang dipotong ke 200 bp**, dan R2 hanya dipakai untuk QC.
- Filter: buang read yang mengandung N atau yang *expected errors*-nya > 2 (Edgar & Flyvbjerg, 2015).
- **Rarefaction:** dari setiap sampel diambil acak 10.000 read yang lolos filter, sehingga kedalaman sekuensing setara. Tidak ada sampel yang dikeluarkan.

### 3.3 Ekstraksi fitur dan anotasi taksonomi (bagian 7–8)

| Representasi | Cara pembuatan | Ukuran tabel | Persentase nol |
|---|---|---|---|
| **ESV** | sekuens 200 bp identik dihitung per sampel (458.107 varian → 12.192 setelah filter) | 292 × 12.192 | 93,2% |
| **k-mer** | frekuensi semua kata 6 basa dari read | 292 × 4.096 | 9,4% |
| **Genus** | ESV dianotasi ke referensi FASTA NCBI 16S, lalu dijumlahkan per genus | 292 × 280 | 70,9% |

- **Anotasi taksonomi:**
  - Dari setiap sekuens referensi diambil 200 bp setelah primer 515F, yaitu posisi yang sama dengan awal read.
  - Tiap ESV dicocokkan ke referensi terdekat.
  - Genus ditetapkan bila identitasnya ≥ 97%. Dengan aturan ini, **78,8% ESV (84,8% read)** terklasifikasi sampai tingkat genus.
- **Keragaman:** pasien kanker memiliki keragaman mikrobioma sedikit lebih rendah.
  - ESV terdeteksi: rerata 788 vs 862.
  - Indeks Shannon: 4,40 vs 4,61.

### 3.4 Validasi silang dengan data olahan penulis (bagian 8b)
Hasil pipeline Python kita dibandingkan dengan tabel OTU 97% yang dipublikasikan penulis paper di GitHub (Schloss Lab). Pipeline penulis memakai mothur, jadi pipelinenya berbeda dari kita.
- Korelasi keragaman Shannon per sampel: **ρ = 0,81**
- Untuk 20 OTU teratas penulis, korelasi terbaik dengan ESV kita: **median 0,88**

Kedua angka ini membuktikan bahwa fitur yang kita hitung dari FASTQ mewakili sampel yang benar. Pemeriksaan ini sempat membongkar sebuah bug penting di versi awal notebook (lihat bagian 9). Data penulis **hanya dipakai sebagai kontrol kualitas**, tidak dipakai di model.

### 3.5 Eksplorasi tanpa label (bagian 9–11)
- **Genus terkait kanker kolorektal** (prevalensi = persentase orang yang membawa bakteri tersebut):

| Genus | Prevalensi pada cancer | Prevalensi pada normal |
|---|---|---|
| *Fusobacterium* | 42,5% | 19,2% |
| *Porphyromonas* | 39,2% | 11,6% |
| *Peptostreptococcus* | 31,7% | 2,3% |
| *Parvimonas* | 30,0% | 6,4% |
| *Gemella* | 34,2% | 19,8% |

- **PCA** (CLR, 1.878 ESV): PC1 menjelaskan 6,1% varians dan PC2 5,8%; 10 PC menjelaskan 27,8%. Varians tersebar luas, khas data mikrobioma.
- **K-Means (k = 2):** silhouette 0,28, tetapi ARI terhadap label hanya −0,02. Pengelompokan alami mikrobioma (misalnya tipe komunitas usus) **tidak sejalan** dengan status kanker, sehingga diperlukan pendekatan supervised.

### 3.6 Uji kelimpahan berbeda (bagian 12–13)
- **Uji:** Mann–Whitney U dua sisi pada kelimpahan relatif, dengan koreksi FDR Benjamini–Hochberg.
  - Uji tidak berpasangan tepat dipakai karena satu orang hanya menyumbang satu sampel, sehingga tidak ada pseudoreplikasi.
  - Uji dilakukan pada kelimpahan relatif, bukan nilai CLR. Pada tabel yang 93% nol, CLR menghasilkan signifikansi palsu (826 ESV "signifikan", kebanyakan semu).
- **Tingkat genus:** 169 genus diuji, **10 signifikan** (FDR < 0,05).
  - Lebih tinggi pada cancer (8): *Peptostreptococcus*, *Porphyromonas*, *Parvimonas*, *Fusobacterium*, *Evtepia*, *Gemella*, *Massiliprevotella*, *Klebsiella*.
  - Lebih tinggi pada normal (2): *Haemophilus*, *Agathobaculum/Butyricicoccus*.
- **Tingkat ESV:** 1.878 ESV diuji, 76 signifikan (12 lebih tinggi pada cancer, 64 lebih tinggi pada normal).
- **Visualisasi:** volcano-style plot dan heatmap 15 genus teratas.

### 3.7 Machine learning (bagian 14–16)
- **Pipeline:** filter prevalensi ≥ 10% → CLR → StandardScaler (khusus LR) → model.
  - Semua langkah dipelajari **hanya dari data training**.
  - `class_weight="balanced"` dipakai untuk menangani kelas yang tidak seimbang.
- **Evaluasi:**
  1. Split *stratified* 80/20: training 233 sampel, test 59 sampel. Tuning dilakukan dengan `GridSearchCV` 5-fold pada training saja.
  2. **Nested 5-fold CV** pada seluruh sampel, sebagai estimasi performa yang lebih stabil.

| Model | Hyperparameter yang di-tuning |
|---|---|
| Logistic Regression | C ∈ {0,001; 0,01; 0,1; 1} |
| Random Forest (500 pohon) | max_features ∈ {sqrt, 0,05}, min_samples_leaf ∈ {1, 3} |

---

## 4. Hasil evaluasi

### 4.1 Nested 5-fold cross-validation (estimasi utama)

| Fitur | Model | ROC-AUC (± sd) | F1 cancer | Accuracy | Balanced acc. |
|---|---|---|---|---|---|
| ESV | Logistic Regression | 0,714 ± 0,019 | 0,581 | 0,634 | 0,631 |
| ESV | Random Forest | 0,680 ± 0,031 | 0,492 | 0,634 | 0,604 |
| Genus | Logistic Regression | 0,745 ± 0,074 | 0,609 | 0,692 | 0,678 |
| **Genus** | **Random Forest** | **0,765 ± 0,090** | 0,573 | **0,729** | **0,688** |
| k-mer | Logistic Regression | 0,751 ± 0,050 | **0,614** | 0,702 | 0,683 |
| **k-mer** | **Random Forest** | **0,770 ± 0,065** | 0,554 | 0,722 | 0,679 |

### 4.2 Test set (59 sampel: 35 normal, 24 cancer; ambang 0,5)

| Fitur | Model | Accuracy | Recall cancer | Specificity | F1 cancer | ROC-AUC |
|---|---|---|---|---|---|---|
| k-mer | Logistic Regression | 0,695 | 0,625 | 0,743 | 0,625 | **0,782** |
| k-mer | Random Forest | 0,712 | 0,417 | 0,914 | 0,541 | 0,756 |
| Genus | Random Forest | 0,695 | 0,417 | 0,886 | 0,526 | 0,735 |
| Genus | Logistic Regression | 0,610 | 0,542 | 0,657 | 0,531 | 0,720 |
| ESV | Random Forest | 0,695 | 0,583 | 0,771 | 0,609 | 0,667 |
| ESV | Logistic Regression | 0,593 | 0,583 | 0,600 | 0,538 | 0,674 |
| — | Dummy (baseline) | 0,593 | 0,000 | 1,000 | 0,000 | 0,500 |

### 4.3 Interpretasi
1. **Semua model jauh di atas baseline.** Accuracy Dummy 0,593 terlihat "lumayan" hanya karena ia selalu menebak normal (kelas mayoritas), tetapi tidak mendeteksi satu kanker pun. Inilah alasan **ROC-AUC dan recall** lebih bermakna daripada accuracy pada data tidak seimbang.
2. **Genus dan k-mer lebih baik daripada ESV.** ESV sangat jarang terisi (93% nol), dan satu bakteri bisa terpecah menjadi banyak ESV akibat galat sekuensing. Menggabungkan ke genus mengatasi masalah ini, sedangkan k-mer bekerja tanpa referensi sama sekali.
3. **LR vs RF tidak berbeda bermakna.**
   - Selisih AUC-nya hanya 0,01–0,03, lebih kecil dari simpangan baku antar-fold.
   - Pada ambang 0,5, RF cenderung **lebih spesifik** (jarang salah menuduh kanker), sedangkan LR **lebih sensitif** (lebih banyak kanker terdeteksi).
4. **Recall kanker masih rendah (0,42–0,63).** Untuk skrining, ambang probabilitas perlu diturunkan agar lebih banyak kanker terdeteksi, dengan konsekuensi lebih banyak positif palsu.

### 4.4 Analisis sensitivitas (bagian 16b)

| Pemeriksaan | Hasil | Makna |
|---|---|---|
| Model dari metrik QC saja (kualitas, GC, % lolos filter) | AUC 0,622 | Ada sedikit perbedaan teknis antarkelompok, tetapi jauh di bawah model berbasis bakteri |
| Tag instrumen "MiSeq" di ENA (hanya pada 26 sampel cancer) | AUC tanpa sampel bertag: genus 0,76, k-mer 0,74–0,76 (hampir tidak berubah) | Tag ini bukan sumber performa |
| Dapatkah fitur menebak tag di antara sampel cancer? | AUC 0,39–0,57 (≈ acak) | Tag hanya kesalahan metadata, bukan batch teknis |
| **Usia** | median 64 th (cancer) vs 54 th (normal); usia saja: AUC 0,715 | **Confounder penting.** Sebagian perbedaan mikrobioma bisa terkait usia, bukan kanker |

---

## 5. Bakteri yang berkontribusi pada model (bagian 17, fitur Genus)

| Genus | RF importance | Koef. LR | Rerata % cancer | Rerata % normal | Arah |
|---|---|---|---|---|---|
| *Porphyromonas* | 0,044 | +0,157 | 0,583 | 0,0015 | ↑ cancer |
| *Peptostreptococcus* | 0,034 | +0,137 | 0,043 | 0,0003 | ↑ cancer |
| *Parvimonas* | 0,021 | +0,110 | 0,040 | 0,0008 | ↑ cancer |
| *Fusobacterium* | 0,020 | +0,109 | 0,332 | 0,0099 | ↑ cancer |
| *Evtepia* | 0,018 | +0,086 | 0,079 | 0,052 | ↑ cancer |
| *Merdibacter* | 0,017 | −0,099 | 0,002 | 0,005 | ↑ normal |
| *Agathobaculum/Butyricicoccus* | 0,014 | −0,101 | 0,085 | 0,140 | ↑ normal |
| *Gemmiger* | 0,013 | −0,121 | 0,536 | 0,806 | ↑ normal |
| *Collinsella* | 0,013 | +0,084 | 1,705 | 1,113 | ↑ cancer |
| *Anaerostipes* | 0,013 | −0,071 | 2,147 | 2,615 | ↑ normal |

- **Empat kontributor teratas** adalah bakteri oral-patogen yang paling sering dilaporkan pada kanker kolorektal, yaitu *Porphyromonas asaccharolytica*, *Peptostreptococcus stomatis*, *Parvimonas micra*, dan *Fusobacterium*. Spesies terdekat ini terlihat dari ESV kontributor teratas. Temuan ini konsisten dengan Baxter et al. (2016), Zeller et al. (2014), Kostic et al. (2012), dan Castellarin et al. (2012).
- **Genus yang lebih tinggi pada normal** antara lain *Agathobaculum/Butyricicoccus*, *Gemmiger*, dan *Anaerostipes*. Beberapa di antaranya (misalnya *Butyricicoccus* dan *Anaerostipes*) dikenal sebagai penghasil butirat. Butirat adalah asam lemak rantai pendek yang umumnya dikaitkan dengan usus sehat. Peran masing-masing genus perlu dicek ke literatur sebelum ditulis di laporan.
- Kontribusi model **bukan bukti sebab-akibat**: belum diketahui apakah bakteri memicu kanker atau kanker yang mengubah lingkungan usus.

---

## 6. Kesimpulan

1. Data FASTQ 16S rRNA mentah dapat diolah dari awal sampai menjadi model machine learning yang mendeteksi kanker kolorektal dengan **ROC-AUC ±0,75–0,78**.
   - Performa ini **sedang**: cukup untuk menunjukkan adanya sinyal biologis, tetapi belum cukup sebagai alat diagnosis tunggal.
   - Paper sumber data juga menggabungkan mikrobioma dengan tes FIT.
2. Representasi **Genus** dan **k-mer** lebih baik daripada **ESV**. **Logistic Regression** dan **Random Forest** berkinerja setara.
3. Model dan uji statistik secara independen menemukan bakteri yang sama (*Porphyromonas*, *Peptostreptococcus*, *Parvimonas*, *Fusobacterium*). Bakteri-bakteri ini dikenal dalam literatur sebagai penanda kanker kolorektal.
4. **Validasi dan kontrol kualitas sangat penting.** Pembandingan dengan data olahan penulis berhasil menemukan bug yang membuat fitur bakteri tertukar antarsampel. Tanpa pemeriksaan ini, kesimpulan proyek akan salah.

## 7. Keterbatasan

- Satu kohor (Amerika Utara, ±292 orang) tanpa validasi eksternal.
- **Usia** berbeda jauh antarkelompok (median 64 vs 54 tahun) dan belum dikoreksi di model.
- Metrik QC saja sudah sedikit memprediksi label (AUC 0,62); ada kemungkinan perbedaan teknis kecil antarkelompok.
- Pipeline sederhana: hanya R1 200 bp, tanpa *denoising* atau penghapusan chimera. Varian galat sekuensing masih ikut dalam tabel ESV.
- Region V4 tidak cukup untuk memastikan spesies; anotasi hanya bisa dipercaya sampai tingkat genus.
- Recall kanker pada ambang 0,5 masih rendah. Ambang perlu dioptimalkan bila tujuannya skrining.
- **Bukan alat diagnosis klinis.**

---

## 8. Cara mereproduksi

1. **Lokal (Jupyter):** buka [Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb](Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb) dari folder ini, lalu *Run all*. Berkas FASTQ di `data/fastq/` dipakai langsung setelah dicek MD5-nya.
2. **Google Colab:**
   - Unggah seluruh folder ke Google Drive.
   - Di bagian 2 notebook, aktifkan dua baris *mount Drive*, lalu *Run all*.
   - Bila folder `data/` tidak ada, notebook mengunduh FASTQ sendiri dari ENA (±4,5 GB).
3. Waktu jalan lokal ±16 menit; bagian terberat adalah nested cross-validation. Seed acak tetap = 42.

### Isi folder

| Berkas / folder | Isi |
|---|---|
| `Proyek2_Kanker_Kolorektal_16S_FASTQ.ipynb` | Notebook lengkap beserta output |
| `flowchart_alur_penelitian.png` | Flowchart metodologi |
| `RINGKASAN_PROYEK.md` | Dokumen ini |
| `data/fastq/` | 584 berkas FASTQ utuh (4,48 GB) |
| `data/metadata_sampel.tsv` | Metadata run + label diagnosis |
| `data/reference/` | Referensi FASTA NCBI 16S + tabel OTU penulis (untuk validasi) |
| `hasil_proyek2/results/` | Tabel CSV, FASTA ESV, dan grafik bernomor 01–11 untuk laporan/presentasi |
| `hasil_proyek2/models/` | Model `.joblib` + skema fitur |
| `hasil_proyek2.zip` | Arsip hasil (tanpa FASTQ) |

### Berkas hasil penting

| Berkas | Isi |
|---|---|
| `qc_fastq_per_sampel.csv` | QC per sampel (R1 & R2) |
| `esv_sequences.fasta` | Sekuens semua ESV (bisa dicek manual lewat BLAST) |
| `taksonomi_esv.csv` | Anotasi genus dan referensi terdekat setiap ESV |
| `mann_whitney_genus.csv` | Hasil uji statistik tingkat genus |
| `perbandingan_model_test.csv` | Metrik test set semua model |
| `kfold_nested_cv.csv` | Hasil nested 5-fold CV |
| `sensitivitas_batch.csv` | Analisis sensitivitas |
| `kontribusi_fitur_genus.csv` | RF importance + koefisien LR |

---

## 9. Riwayat pengerjaan

1. **Pemilihan dataset:** tugas mewajibkan data berformat FASTQ, sehingga dipilih data amplikon 16S rRNA dari Baxter et al. (2016) yang label diagnosisnya tersedia di metadata resmi.
2. **Unduhan:** 584 FASTQ utuh (R1 + R2, 4,48 GB), diverifikasi MD5.
3. **Prototipe:** pipeline diuji pada 6.000 read pertama per sampel.
4. **Bug ditemukan lewat validasi silang.**
   - Gejala: pada data utuh, fitur bakteri tampak acak (AUC ≈ 0,5) dan tidak berkorelasi dengan tabel OTU penulis (ρ ≈ 0), padahal jumlah read per sampel cocok (ρ = 0,95).
   - Penyebab: `pandas.DataFrame.from_dict(orient="index")` mengubah urutan baris, sehingga profil bakteri tertukar antarsampel.
   - Perbaikan: urutan baris dikunci, dan ditambahkan pemeriksaan otomatis yang menghitung ulang sampel acak.
   - Setelah diperbaiki, korelasi dengan data penulis naik menjadi ρ = 0,81 dan bakteri penanda kanker muncul jelas.
5. **Perbaikan metode uji:** uji Mann–Whitney dipindah dari nilai CLR ke kelimpahan relatif, karena CLR pada tabel yang sangat jarang menghasilkan signifikansi semu.
6. **Verifikasi akhir:** notebook dieksekusi penuh pada data utuh tanpa error.

---

## Referensi utama

- Baxter NT, Ruffin MT, Rogers MAM, Schloss PD (2016). *Microbiota-based model improves the sensitivity of fecal immunochemical test for detecting colonic lesions.* Genome Medicine 8:37. [doi:10.1186/s13073-016-0290-3](https://doi.org/10.1186/s13073-016-0290-3). Sumber dataset.
- Kozich JJ et al. (2013). *Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform.* Applied and Environmental Microbiology 79(17):5112–5120. [doi:10.1128/AEM.01043-13](https://doi.org/10.1128/AEM.01043-13)
- Edgar RC, Flyvbjerg H (2015). *Error filtering, pair assembly and error correction for next-generation sequencing reads.* Bioinformatics 31(21):3476–3482. [doi:10.1093/bioinformatics/btv401](https://doi.org/10.1093/bioinformatics/btv401)
- Callahan BJ, McMurdie PJ, Holmes SP (2017). *Exact sequence variants should replace operational taxonomic units in marker-gene data analysis.* The ISME Journal 11:2639–2643. [doi:10.1038/ismej.2017.119](https://doi.org/10.1038/ismej.2017.119)
- Gloor GB et al. (2017). *Microbiome datasets are compositional: and this is not optional.* Frontiers in Microbiology 8:2224. [doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)
- Zeller G et al. (2014). *Potential of fecal microbiota for early-stage detection of colorectal cancer.* Molecular Systems Biology 10:766. [doi:10.15252/msb.20145645](https://doi.org/10.15252/msb.20145645)
- Kostic AD et al. (2012). *Genomic analysis identifies association of Fusobacterium with colorectal carcinoma.* Genome Research 22:292–298. [doi:10.1101/gr.126573.111](https://doi.org/10.1101/gr.126573.111)
- Castellarin M et al. (2012). *Fusobacterium nucleatum infection is prevalent in human colorectal carcinoma.* Genome Research 22:299–306. [doi:10.1101/gr.126516.111](https://doi.org/10.1101/gr.126516.111)
- [NCBI RefSeq 16S rRNA](https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/) · [Repositori penulis (Schloss Lab)](https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015)
