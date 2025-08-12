# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd

sys.path.append('scripts/Downstream_analysis/')
from GO_analysis import (
    set_plotting_style, 
    run_gsea_enrichment, 
    cluster_and_order_terms, 
    plot_go_dotplot
)

def main():
    set_plotting_style()
    
    # --- Configuration ---
    module_file_path = "results/downstream_analysis/05_polya_modules/polya_modules.csv"
    output_dir = "results/downstream_analysis/09_go_analysis_polya"

    gmt_file_path = "data/annotation/msigdb_v2024.1.Mm_GMTs/m5.all.v2024.1.Mm.symbols.gmt"
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading poly(A) modules from {module_file_path}...")
    modules_df = pd.read_csv(module_file_path)
    if 'gene_symbol' not in modules_df.columns:
        modules_df['gene_symbol'] = modules_df['Gene'].str.split('-', 1).str[0]

    enriched_df = run_gsea_enrichment(modules_df, gmt_file_path, module_col_name='Module')
    if enriched_df.empty: return

    final_df = cluster_and_order_terms(enriched_df)
    if final_df.empty: return

    final_df.to_csv(os.path.join(output_dir, "go_enrichment_results_polya.csv"), index=False)
    print(f"Final GO results table saved in {output_dir}")

    plot_go_dotplot(
        plot_df=final_df,
        output_path=os.path.join(output_dir, "go_enrichment_dotplot_polya.pdf"),
        title="GO Enrichment for Poly(A) Modules"
    )

if __name__ == "__main__":
    main()