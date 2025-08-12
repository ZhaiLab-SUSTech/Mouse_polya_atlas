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
    gene_sets_input_path = "results/downstream_analysis/12_fig4ab_scatters/fig4_gene_sets.csv"
    output_dir = "results/downstream_analysis/13_fig4c_go_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    gmt_file_path = "data/annotation/msigdb_v2024.1.Mm_GMTs/m2.cp.kegg.v2024.1.Mm.symbols.gmt"
    
    if not os.path.exists(gene_sets_input_path):
        print(f"Error: Gene set file not found!", file=sys.stderr)
        print(f"Please run '15_plot_fig4ab_scatter.py' first to generate it.", file=sys.stderr)
        sys.exit(1)
        
    modules_df = pd.read_csv(gene_sets_input_path)
    modules_df.rename(columns={'group': 'Module'}, inplace=True)

    enriched_df = run_gsea_enrichment(modules_df, gmt_file_path, module_col_name='Module')
    if enriched_df.empty: 
        print("No enrichment results found. Exiting.")
        return

    final_df = cluster_and_order_terms(enriched_df)
    if final_df.empty:
        print("No terms left after clustering. Exiting.")
        return
        
    final_df.to_csv(os.path.join(output_dir, "go_enrichment_results_fig4c.csv"), index=False)

    plot_go_dotplot(
        plot_df=final_df,
        output_path=os.path.join(output_dir, "fig4c_go_enrichment_dotplot.pdf"),
        title="GO Enrichment: Elongated vs Consistent",
        module_order=['elongated', 'consistent']
    )
    
    print(f"\nGO analysis for Fig 4c complete. Results saved to {output_dir}")

if __name__ == "__main__":
    main()