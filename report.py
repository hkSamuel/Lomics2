from setting import *
import setting
import time

def lomics_report(input_question, var_llm, var_maxtoken, var_temp, var_max_attempt, output_name, output_dir):
    ############################################################################################################
    # Functions
    async def prompt_explain(input_question, var_llm, var_maxtoken, var_temp, pathway, limit):
        prompt = f'''Instructions:
            - You are a bioinformatician analyzing transcriptomic data. 
            - You are tasked with explaining succinctly with 1 sentence why the following pathway: [ {pathway} ] is relevant for analysing the scientific question: [ {input_question} ].
            - You are requred to output a JSON object only.
            - You must not return any words before and after the JSON object only.
            - JSON object schema is specified in the following pydantic description:
                class pydantic_explain(BaseModel):
                    explain: str = Field(description="explain")'''
        try:
            response = await llm_call(prompt, var_llm, var_maxtoken, var_temp, limit)
            return response
        except Exception as e:
            setting.increment_error_count('report_prompt_error_count', f"prompt_explain failed for '{pathway}': {e}")
            print(f"prompt_explain failed.")

    async def prompt_df_explain(input_question, var_llm, var_maxtoken, var_temp, pathway, var_max_attempt, limit):
        for i in range(var_max_attempt):
            try:
                explain = await prompt_explain(input_question, var_llm, var_maxtoken, var_temp, pathway, limit)
                json_dict = loads(explain)
                json_dict['explain'] = re.sub(r',', '', json_dict['explain']) 
                df_explain = pd.DataFrame([{"pathway": pathway, "pathway_explain": json_dict['explain']}])
                df_explain["llm"] = var_llm
                df_explain["temp"] = var_temp
                return df_explain
            except Exception as e:
                setting.increment_json_error_count(f"report explain parse failed for '{pathway}': {e}")
                setting.increment_error_count('report_parse_error_count', f"report explain parse failed for '{pathway}': {e}")
                print(f"prompt_df_explain failed for pathway '{pathway}', attempt {i+1}/{var_max_attempt}: {e}")
        
        # If all attempts fail, return a DataFrame with a placeholder explanation
        print(f"WARNING: All attempts failed for pathway '{pathway}'. Using placeholder.")
        df_explain = pd.DataFrame([{"pathway": pathway, "pathway_explain": "Explanation unavailable"}])
        df_explain["llm"] = var_llm
        df_explain["temp"] = var_temp
        return df_explain

    async def batch_prompt_df_explain(input_question, var_llm, var_maxtoken, var_temp, var_max_attempt, output_name, output_dir):
        limit = asyncio.Semaphore(var_max_concurrent_call)
        tasks = [prompt_df_explain(input_question, var_llm, var_maxtoken, var_temp, pathway, var_max_attempt, limit) for pathway in ls_pathway]
        ls_df_explain = await asyncio.gather(*tasks)
        df_explain = pd.concat(ls_df_explain, ignore_index=True)
        
        # Add pathway repeat counts from the pathway CSV
        df_pathway_all = pd.read_csv(os.path.join(output_dir, output_name + "_pathway.csv"))
        repeat_counts = df_pathway_all['pathway'].value_counts().to_dict()
        df_explain['pathway_repeat_count'] = df_explain['pathway'].map(repeat_counts)
        
        # Prepend a title line with total tokens used and total execution time
        total_tokens = getattr(setting, 'total_tokens_used', None)
        start_time = getattr(setting, 'run_start_time', None)
        json_errors = getattr(setting, 'json_error_count', 0)
        llm_call_errors = getattr(setting, 'llm_call_error_count', 0)
        llm_jsondecode_expecting_value_errors = getattr(setting, 'llm_jsondecode_expecting_value_error_count', 0)
        llm_empty_content_errors = getattr(setting, 'llm_empty_content_error_count', 0)
        llm_finish_reason_length_count = getattr(setting, 'llm_finish_reason_length_count', 0)
        pathway_prompt_errors = getattr(setting, 'pathway_prompt_error_count', 0)
        pathway_schema_errors = getattr(setting, 'pathway_schema_error_count', 0)
        pathway_parse_errors = getattr(setting, 'pathway_parse_error_count', 0)
        report_prompt_errors = getattr(setting, 'report_prompt_error_count', 0)
        report_parse_errors = getattr(setting, 'report_parse_error_count', 0)
        max_concurrent_call = getattr(setting, 'var_max_concurrent_call', var_max_concurrent_call)
        generated_at = time.strftime('%Y-%m-%d %H:%M:%S')
        if start_time:
            elapsed = time.time() - start_time
        else:
            elapsed = None

        title_line = (
            f"# Report generated: generated_at={generated_at}, total_tokens={total_tokens if total_tokens is not None else 'unknown'}, "
            f"elapsed_seconds={elapsed:.2f}, json_errors_total={json_errors}, llm_call_errors={llm_call_errors}, llm_jsondecode_expecting_value_errors={llm_jsondecode_expecting_value_errors}, llm_empty_content_errors={llm_empty_content_errors}, llm_finish_reason_length_count={llm_finish_reason_length_count}, pathway_prompt_errors={pathway_prompt_errors}, pathway_schema_errors={pathway_schema_errors}, pathway_parse_errors={pathway_parse_errors}, report_prompt_errors={report_prompt_errors}, report_parse_errors={report_parse_errors}, maxtoken={var_maxtoken}, maxconcurrent_call={max_concurrent_call}"
            if elapsed is not None else
            f"# Report generated: generated_at={generated_at}, total_tokens={total_tokens if total_tokens is not None else 'unknown'}, "
            f"elapsed_seconds=unknown, json_errors_total={json_errors}, llm_call_errors={llm_call_errors}, llm_jsondecode_expecting_value_errors={llm_jsondecode_expecting_value_errors}, llm_empty_content_errors={llm_empty_content_errors}, llm_finish_reason_length_count={llm_finish_reason_length_count}, pathway_prompt_errors={pathway_prompt_errors}, pathway_schema_errors={pathway_schema_errors}, pathway_parse_errors={pathway_parse_errors}, report_prompt_errors={report_prompt_errors}, report_parse_errors={report_parse_errors}, maxtoken={var_maxtoken}, maxconcurrent_call={max_concurrent_call}"
        )

        report_path = os.path.join(output_dir, output_name + "_report.csv")
        # Write title line then append dataframe
        with open(report_path, 'w', encoding='utf-8') as fh:
            fh.write(title_line + "\n")
        df_explain.to_csv(report_path, mode='a', index=False)

    def output_gmx(df_gene, ls_pathway, gmx_output_path):
        pathway_dict = {}
        for pathway in ls_pathway:
            df_pathway = df_gene[df_gene["pathway"] == pathway]
            pathway_dict[pathway] = [str(int(g)) if isinstance(g, float) else str(g) for g in df_pathway["gene"].tolist()]
            with open(gmx_output_path, "w") as f:
                for pathway, genes in pathway_dict.items():
                    f.write(f"{pathway}\t")
                    f.write('\t'.join(genes))
                    f.write("\n")
        return True

    def output_gmt(df_gene, ls_pathway, gmt_output_path):
        pathway_dict = {}
        for pathway in ls_pathway:
            df_pathway = df_gene[df_gene["pathway"] == pathway]
            pathway_dict[pathway] = [str(int(g)) if isinstance(g, float) else str(g) for g in df_pathway["gene"].tolist()]
            with open(gmt_output_path, "w") as f:
                for pathway, genes in pathway_dict.items():
                    # Format: pathway_name    pathway_description   gene1   gene2   gene3...
                    f.write(f"{pathway}\t{pathway}\t")
                    f.write('\t'.join(genes))
                    f.write("\n")
        return True

    def transpose_file(gmx_output_path):
        with open(gmx_output_path, 'r') as input_file:
            reader = csv.reader(gmx_output_path, delimiter='\t')
            rows = [row for row in reader]
        max_length = max(len(row) for row in rows)
        transposed = [[row[i] if i < len(row) else '' for row in rows] for i in range(max_length)]
        with open(gmx_output_path, 'w', newline='') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            for row in transposed:
                writer.writerow(row)

    def transpose_file(gmx_output_path):
        with open(gmx_output_path, 'r') as input_file:
            reader = csv.reader(input_file, delimiter='\t')
            rows = [row for row in reader]
        max_length = max(len(row) for row in rows)
        transposed = [[row[i] if i < len(row) else '' for row in rows] for i in range(max_length)]
        with open(gmx_output_path, 'w', newline='') as output_file:
            writer = csv.writer(output_file, delimiter='\t')
            for row in transposed:
                writer.writerow(row)

    ############################################################################################################
    # Execution
    # Get pathways from pathway.csv (top 100 by frequency) instead of gene.csv
    # This ensures all selected pathways get explanations, even if some have no genes
    df_pathway = pd.read_csv(os.path.join(output_dir, output_name + "_pathway.csv"))
    ls_pathway = df_pathway['pathway'].value_counts().nlargest(var_num_pathway).index.tolist()
    print(f"\nGenerating explanations for {len(ls_pathway)} pathways...")
    
    asyncio.run(batch_prompt_df_explain(input_question, var_llm, var_maxtoken, var_temp, var_max_attempt, output_name, output_dir))
    
    # Now load gene data for GMT/GMX generation
    df_gene = pd.read_csv(os.path.join(output_dir, output_name + "_gene.csv"))
    ls_pathway_with_genes = df_gene['pathway'].unique()
    print(f"Found {len(ls_pathway_with_genes)} pathways with valid genes in gene.csv")
    
    if len(ls_pathway_with_genes) < len(ls_pathway):
        print(f"WARNING: {len(ls_pathway) - len(ls_pathway_with_genes)} pathways have no valid genes and will not appear in GMT/GMX files")
    
    print(f"\nBefore filtering - Total genes: {len(df_gene)}")
    print(f"Valid Entrez IDs: {len(df_gene[df_gene['entrez_valid'] == True])}")
    print(f"Invalid Entrez IDs: {len(df_gene[df_gene['entrez_valid'] == False])}")
    
    df_gene = df_gene[df_gene['entrez_valid'] == True]
    
    if len(df_gene) == 0:
        print("ERROR: No valid Entrez IDs found! Cannot generate GMT/GMX files.")
        print("Sample invalid genes:")
        df_all = pd.read_csv(os.path.join(output_dir, output_name + "_gene.csv"))
        print(df_all[df_all['entrez_valid'] == False]['gene'].head(20).tolist())
        return
    
    df_gene = df_gene.groupby(['pathway', 'gene']).size().reset_index(name='count')
    df_gene = df_gene.sort_values(['pathway', 'count'], ascending=[True, False])
    df_gene = df_gene.groupby('pathway').head(var_num_gene).reset_index(drop=True)
    
    # Use pathways that actually have genes for GMT/GMX output
    output_gmx(df_gene, ls_pathway_with_genes, os.path.join(output_dir, output_name + ".gmx"))
    output_gmt(df_gene, ls_pathway_with_genes, os.path.join(output_dir, output_name + ".gmt"))
    transpose_file(os.path.join(output_dir, output_name + ".gmx"))