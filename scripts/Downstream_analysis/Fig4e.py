import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats import multitest as multi
from itertools import combinations
import os

def standardize_pair(g1, g2):
    if pd.isna(g1) or pd.isna(g2):
        return None
    return tuple(sorted([g1, g2]))

def get_gene_symbol_from_isoform(isoform_id):
    try:
        return isoform_id.split('.')[0][:-19]
    except (AttributeError, IndexError):
        return None

def load_and_prepare_data(base_path):
    # Load MouseNet v2 network data
    mousenet_path = os.path.join(base_path, "MouseNetV2_symbol.txt")
    mousenet_df = pd.read_csv(mousenet_path, sep='\t', names=['gene1', 'gene2', 'score'])
    mousenet_df['sorted_pair'] = mousenet_df.apply(
        lambda row: standardize_pair(row['gene1'], row['gene2']), axis=1
    )
    mousenet_interaction_set = set(mousenet_df['sorted_pair'].dropna())

    polya_module_path = os.path.join(base_path, "polya_modules.csv")
    polya_module = pd.read_csv(polya_module_path, sep=',')
    polya_module.columns = ['Gene', 'Module_polya', 'gene_symbol', 'gene_id', 'gene_id_short']

    expression_module_path = os.path.join(base_path, 'expression_modules.csv')
    expression_module = pd.read_csv(expression_module_path, sep=',')
    expression_module.columns = ['Gene', 'Module_expression']

    polya_module['correct_gene_symbol'] = polya_module['Gene'].apply(get_gene_symbol_from_isoform)

    merged_df = pd.merge(polya_module, expression_module, on='Gene', how='inner')
    merged_df.dropna(subset=['correct_gene_symbol', 'Module_polya', 'Module_expression'], inplace=True)
    
    return merged_df, mousenet_interaction_set

def classify_pair_by_isoform(row, id_to_modules_map):
    iso1_mods = id_to_modules_map.get(row['isoform1'])
    iso2_mods = id_to_modules_map.get(row['isoform2'])
    
    if not iso1_mods or not iso2_mods:
        return "Invalid"
    
    in_same_polya = iso1_mods['Module_polya'] == iso2_mods['Module_polya']
    in_same_expr = iso1_mods['Module_expression'] == iso2_mods['Module_expression']
    
    if in_same_polya and not in_same_expr:
        return "A_polyA_specific"
    elif not in_same_polya and in_same_expr:
        return "B_expr_specific"
    elif in_same_polya and in_same_expr:
        return "C_both_concordant"
    else:  # not in_same_polya and not in_same_expr
        return "D_both_discordant"

def analyze_enrichment(merged_df, mousenet_set):
    
    analyzed_isoforms = merged_df['Gene'].unique().tolist()
    isoform_to_data_map = merged_df.set_index('Gene').to_dict('index')

    universe_of_pairs = list(combinations(analyzed_isoforms, 2))
    universe_df = pd.DataFrame(universe_of_pairs, columns=['isoform1', 'isoform2'])
    universe_df['category'] = universe_df.apply(classify_pair_by_isoform, axis=1, id_to_modules_map=isoform_to_data_map)
    print(f"Classified {len(universe_df)} total isoform pairs.")

    isoform_to_gene_map = merged_df.set_index('Gene')['correct_gene_symbol'].to_dict()
    categories = ["A_polyA_specific", "B_expr_specific", "C_both_concordant", "D_both_discordant"]
    contingency_data = []

    for cat in categories:
        cat_isoform_df = universe_df[universe_df['category'] == cat]
        g1_symbols = cat_isoform_df['isoform1'].map(isoform_to_gene_map)
        g2_symbols = cat_isoform_df['isoform2'].map(isoform_to_gene_map)
        
        cat_gene_pairs = {
            standardize_pair(s1, s2) 
            for s1, s2 in zip(g1_symbols, g2_symbols) 
            if s1 is not None and s2 is not None and s1 != s2
        }
        
        in_mousenet = len(cat_gene_pairs.intersection(mousenet_set))
        not_in_mousenet = len(cat_gene_pairs) - in_mousenet
        contingency_data.append([in_mousenet, not_in_mousenet])

    table_df = pd.DataFrame(contingency_data,
                            index=categories,
                            columns=['In_MouseNet', 'Not_in_MouseNet'])
    table_df['Total_Gene_Pairs'] = table_df.sum(axis=1)
    table_df['Enrichment_Ratio'] = table_df['In_MouseNet'] / table_df['Total_Gene_Pairs']

    print("\n--- 4x2 Contingency Table (based on unique derived Gene pairs) ---")
    print(table_df)

    print("\n--- Overall Chi-squared Test ---")
    chi2, p_value, _, _ = chi2_contingency(table_df[['In_MouseNet', 'Not_in_MouseNet']])
    print(f"Chi-squared statistic: {chi2:.4f}")
    print(f"P-value: {p_value:.4e}")
    if p_value < 0.05:
        print("Result is significant: The proportion of MouseNet interactions differs across categories.")
    else:
        print("Result is not significant: No evidence that proportions differ across categories.")

    print("\n--- Post-hoc Pairwise Comparisons (Fisher's Exact Test) ---")
    comparisons = [
        ("A_polyA_specific", "B_expr_specific"),
        ("A_polyA_specific", "C_both_concordant"),
        ("A_polyA_specific", "D_both_discordant"),
        ("B_expr_specific", "C_both_concordant")
    ]
    p_values_pairwise = []
    print(f"{'Comparison':<43} | {'Odds Ratio':<12} | P-value")
    print("-" * 65)
    for cat1, cat2 in comparisons:
        sub_table = table_df.loc[[cat1, cat2], ['In_MouseNet', 'Not_in_MouseNet']]
        odds_ratio, p = fisher_exact(sub_table, alternative='two-sided')
        p_values_pairwise.append(p)
        print(f"{cat1:<20} vs {cat2:<20} | {odds_ratio:<12.4f} | {p:<.4e}")

    print("\n--- Corrected P-values (Bonferroni) ---")
    reject, pvals_corrected, _, _ = multi.multipletests(p_values_pairwise, alpha=0.05, method='bonferroni')
    for i, (cat1, cat2) in enumerate(comparisons):
        print(f"{cat1} vs {cat2}: Corrected p = {pvals_corrected[i]:.4e}, Significant: {reject[i]}")

if __name__ == "__main__":
    DATA_PATH = "../../data/annotation"
    
    merged_module_data, mousenet_interaction_set = load_and_prepare_data(DATA_PATH)
    analyze_enrichment(merged_module_data, mousenet_interaction_set)