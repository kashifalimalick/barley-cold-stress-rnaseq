import pandas as pd
from gprofiler import GProfiler

# Load DESeq2 results
results = pd.read_csv("counts/deseq2_results.csv", index_col=0)

# Get significant genes: padj < 0.05 and |log2FC| > 1
sig_genes = results[(results["padj"] < 0.05) & (results["log2FoldChange"].abs() > 1)]
gene_list = list(sig_genes.index)

print(f"Number of significant genes for enrichment: {len(gene_list)}")
print("Example gene IDs:", gene_list[:5])

# Query g:Profiler for barley (Hordeum vulgare)
gp = GProfiler(return_dataframe=True)
enrich_results = gp.profile(
    organism="hvulgare",
    query=gene_list,
    sources=["GO:BP", "GO:MF", "GO:CC", "KEGG"]
)

print(f"\nTotal enrichment terms found: {len(enrich_results)}")
enrich_results.to_csv("counts/go_kegg_enrichment.csv", index=False)
print("Saved to counts/go_kegg_enrichment.csv")

if len(enrich_results) > 0:
    print("\nTop 10 enriched terms:")
    print(enrich_results[["source", "name", "p_value", "intersection_size"]].head(10).to_string())
else:
    print("\nNo enrichment terms returned. This usually means the gene ID format doesn't match g:Profiler's barley database.")
