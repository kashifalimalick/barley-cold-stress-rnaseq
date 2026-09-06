import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ---- Load count matrix ----
counts_file = "counts/barley_gene_counts_paired.txt"
df = pd.read_csv(counts_file, sep="\t", skiprows=1)

# featureCounts output: first 6 cols are annotation, rest are sample counts
sample_cols = df.columns[6:]
counts = df[["Geneid"] + list(sample_cols)].set_index("Geneid")
counts.columns = [c.split("/")[-1].replace(".sorted.bam", "") for c in counts.columns]
counts = counts.T  # samples as rows, genes as columns (required by pydeseq2)

print("Samples found:", list(counts.index))
print("Number of genes:", counts.shape[1])

# ---- Metadata: genotype/condition per sample ----
metadata = pd.DataFrame({
    "condition": ["control", "cold", "cold"]
}, index=["ERR2930754", "ERR2930757", "ERR2930761"])

metadata = metadata.loc[counts.index]  # ensure order matches

# ---- Filter low-count genes ----
counts = counts.loc[:, counts.sum(axis=0) >= 10]
print("Genes after filtering low counts:", counts.shape[1])

# ---- Run DESeq2 ----
dds = DeseqDataSet(counts=counts, metadata=metadata, design="~condition", refit_cooks=True)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "cold", "control"])
stat_res.summary()
results_df = stat_res.results_df
results_df.to_csv("counts/deseq2_results.csv")
print("\nSaved results to counts/deseq2_results.csv")

# ---- PCA plot ----
from sklearn.decomposition import PCA
norm_counts = np.log2(dds.layers["normed_counts"] + 1)
pca = PCA(n_components=2)
pcs = pca.fit_transform(norm_counts)
plt.figure(figsize=(6,5))
for i, sample in enumerate(counts.index):
    plt.scatter(pcs[i,0], pcs[i,1], s=100, label=f"{sample} ({metadata.loc[sample,'condition']})")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA of Barley Samples")
plt.legend()
plt.tight_layout()
plt.savefig("counts/pca_plot.png", dpi=150)
plt.close()

# ---- MA plot ----
plt.figure(figsize=(6,5))
plt.scatter(np.log10(results_df["baseMean"]+1), results_df["log2FoldChange"], s=5, alpha=0.4)
plt.axhline(0, color="red", lw=1)
plt.xlabel("log10(mean expression)")
plt.ylabel("log2 fold change")
plt.title("MA Plot: Cold vs Control")
plt.tight_layout()
plt.savefig("counts/ma_plot.png", dpi=150)
plt.close()

# ---- Volcano plot ----
results_df["neglog10padj"] = -np.log10(results_df["padj"].fillna(1))
plt.figure(figsize=(6,5))
sig = (results_df["padj"] < 0.05) & (results_df["log2FoldChange"].abs() > 1)
plt.scatter(results_df.loc[~sig,"log2FoldChange"], results_df.loc[~sig,"neglog10padj"], s=5, color="grey", alpha=0.4)
plt.scatter(results_df.loc[sig,"log2FoldChange"], results_df.loc[sig,"neglog10padj"], s=5, color="red")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(adjusted p-value)")
plt.title("Volcano Plot: Cold vs Control")
plt.tight_layout()
plt.savefig("counts/volcano_plot.png", dpi=150)
plt.close()

# ---- Heatmap of top 30 genes ----
top_genes = results_df.dropna(subset=["padj"]).sort_values("padj").head(30).index
heat_data = pd.DataFrame(norm_counts, index=counts.index, columns=counts.columns)[top_genes].T
plt.figure(figsize=(6,8))
sns.heatmap(heat_data, cmap="RdBu_r", center=heat_data.values.mean(), yticklabels=True)
plt.title("Top 30 Differentially Expressed Genes")
plt.tight_layout()
plt.savefig("counts/heatmap_top30.png", dpi=150)
plt.close()

print("\nAll plots saved in counts/: pca_plot.png, ma_plot.png, volcano_plot.png, heatmap_top30.png")
print("\nTop 10 significant genes:")
print(results_df.dropna(subset=["padj"]).sort_values("padj").head(10)[["baseMean","log2FoldChange","padj"]])
