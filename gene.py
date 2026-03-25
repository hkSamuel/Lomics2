from setting import *

def lomics_gene(input_question, var_llm, var_maxtoken, var_temp, var_iterate_gene, var_num_gene, var_max_concurrent_call, var_max_attempt, output_name, output_dir):
    ############################################################################################################
    # Functions
    def extract_genes_from_gmt(pathway, gmt_dict, ls_entrez_valid, var_num_gene):
        """
        Extract genes directly from GMT file for a given pathway.
        Returns up to var_num_gene genes, or all available if fewer exist.
        """
        # Get genes for this pathway from GMT
        genes = gmt_dict.get(pathway, [])
        
        if not genes:
            print(f"WARNING: No genes found for pathway '{pathway}' in GMT file.")
            return pd.DataFrame()
        
        # Filter to only valid Entrez IDs
        valid_genes = [g for g in genes if g and str(g) in ls_entrez_valid]
        
        if not valid_genes:
            print(f"WARNING: No valid Entrez IDs found for pathway '{pathway}'.")
            return pd.DataFrame()
        
        # Take up to var_num_gene genes (or all if fewer available)
        selected_genes = valid_genes[:var_num_gene] if len(valid_genes) > var_num_gene else valid_genes
        
        print(f"Pathway '{pathway}': Selected {len(selected_genes)} genes out of {len(genes)} total ({len(valid_genes)} valid)")
        
        # Create DataFrame with the format matching the old LLM-based approach
        df_genes = pd.DataFrame({
            'pathway': [pathway] * len(selected_genes),
            'gene': selected_genes,
            'succeeded_iterate': [0] * len(selected_genes),
            'var_iterate_gene': [var_iterate_gene] * len(selected_genes),
            'llm': ['gmt_direct'] * len(selected_genes),
            'temp': [0] * len(selected_genes),
            'num_pathway': [var_num_pathway] * len(selected_genes),
            'num_gene': [var_num_gene] * len(selected_genes),
            'entrez_valid': [True] * len(selected_genes)
        })
        
        return df_genes

    ############################################################################################################
    # Execution
    print("\n" + "="*80)
    print("GENE EXTRACTION: Using genes directly from GMT file (no LLM)")
    print("="*80)
    
    # Load HGNC data to get valid Entrez IDs
    with open(path_hgnc, 'r', encoding='utf-8') as file:
        json_hgnc = load(file)
    
    ls_entrez_valid = set([str(item['entrez_id']) for item in json_hgnc['response']['docs'] if 'entrez_id' in item])
    print(f"Loaded {len(ls_entrez_valid)} valid Entrez IDs from HGNC database")
    
    # Load pathway list
    df_pathway = pd.read_csv(os.path.join(output_dir, output_name + "_pathway.csv"))
    ls_pathway = df_pathway['pathway'].value_counts().nlargest(var_num_pathway).index.tolist()
    print(f"Processing {len(ls_pathway)} pathways")
    
    # Load GMT file
    gmt_dict = load_gmt(path_gmt)
    print(f"Loaded {len(gmt_dict)} pathways from GMT file\n")
    
    # Extract genes for each pathway
    df_output_list = []
    for pathway in ls_pathway:
        df_genes = extract_genes_from_gmt(pathway, gmt_dict, ls_entrez_valid, var_num_gene)
        if not df_genes.empty:
            df_output_list.append(df_genes)
    
    # Combine all results
    if df_output_list:
        df_output = pd.concat(df_output_list, ignore_index=True)
        df_output.to_csv(os.path.join(output_dir, output_name + "_gene.csv"), index=False)
        print(f"\n{'='*80}")
        print(f"SUMMARY: Extracted {len(df_output)} total gene entries across {len(ls_pathway)} pathways")
        print(f"Output saved to: {os.path.join(output_dir, output_name + '_gene.csv')}")
        print(f"{'='*80}\n")
    else:
        print("ERROR: No genes extracted for any pathway!")

