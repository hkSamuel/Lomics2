from setting import *
import setting

def lomics_pathway(input_question, var_num_pathway, var_llm, var_maxtoken, var_temp, var_iterate_pathway, var_max_concurrent_call, var_max_attempt, output_name, output_dir):
    ############################################################################################################
    # Functions
    async def prompt_pathway(input_question, var_num_pathway, var_llm, var_maxtoken, var_temp, limit, var_max_attempt, gmt_pathways):
        pathways_list = list(gmt_pathways.keys()) if gmt_pathways else []
        pathways_str = ', '.join(pathways_list)
        prompt = f'''Instructions:
            - You are a bioinformatician analyzing transcriptomic data. 
            - You are tasked with selecting exactly[ {str(var_num_pathway)} ] pathways for analysis. 
            - Your scientific question is: [ {input_question} ]. 
            - Restrictively reference the following pathways from WikiPathways: {pathways_str}
            - You are requred to output a JSON object only.
            - You must not return any words before and after the JSON object.
            - JSON object schema is specified in the following pydantic description:
                class pydantic_pathway(BaseModel):'''
        for i in range(var_num_pathway):
            prompt += f'''
                    pathway{i+1}: str = Field(description="pathway{i+1}")'''
        for attempt in range(var_max_attempt):
            try:
                response = await llm_call(prompt, var_llm, var_maxtoken, var_temp, limit)
                return response
            except Exception as e:
                setting.increment_error_count('pathway_prompt_error_count', f"pathway prompt failed (attempt {attempt+1}/{var_max_attempt}): {e}")
                continue
        return None

    async def json_schema_pathway(llm_output, var_num_pathway):
        if llm_output is None:
            return None
        fields = {}
        for i in range(1, var_num_pathway+1):
            column_name = f"pathway{i}"
            fields[column_name] = (str, Field(..., title=column_name))
        pydantic_pathway_object = create_model("pathways", **fields, __base__=BaseModel)
        try:
            json_dict = loads(llm_output)
            validated_data = pydantic_pathway_object.parse_obj(json_dict)
            return llm_output
        except ValueError as e:
            setting.increment_json_error_count(f"pathway schema validation failed: {e}")
            setting.increment_error_count('pathway_schema_error_count', f"pathway schema validation failed: {e}")
            print("Pydantic validation error:", e)
            return None
        
    async def batch_prompt_pathway(input_question, var_num_pathway, var_llm, var_maxtoken, var_temp, var_iterate_pathway, var_max_concurrent_call, var_max_attempt, output_name, output_dir):
        gmt_dict = load_gmt(path_gmt)
        limit = asyncio.Semaphore(var_max_concurrent_call)
        target_count = var_num_pathway * var_iterate_pathway
        rows = []
        succeeded_iterate = 0
        gmt_set = set(gmt_dict.keys()) if gmt_dict else set()

        # Only iterate var_iterate_pathway times; retry within each iteration until valid
        for iteration in range(var_iterate_pathway):
            print(f"\n{'='*50}")
            print(f"Iteration {iteration + 1}/{var_iterate_pathway} | LLM: {var_llm}")
            print(f"{'='*50}")
            while True:
                llm_output = await prompt_pathway(
                    input_question,
                    var_num_pathway,
                    var_llm,
                    var_maxtoken,
                    var_temp,
                    limit,
                    var_max_attempt,
                    gmt_dict,
                )
                validated_llm_output = await json_schema_pathway(llm_output, var_num_pathway)
                if validated_llm_output is None:
                    continue

                try:
                    json_dict = loads(validated_llm_output)
                    json_ls = [value for key, value in json_dict.items()]
                    json_ls = [re.sub(r',', '', value).strip() for value in json_ls]

                    for pathway in json_ls:
                        rows.append({
                            "pathway": pathway,
                            "succeeded_iterate": succeeded_iterate,
                            "var_iterate_pathway": var_iterate_pathway,
                            "llm": var_llm,
                            "temp": var_temp,
                            "num_pathway": var_num_pathway,
                            "filled": False,
                        })
                    succeeded_iterate += 1
                    break
                except Exception as e:
                    setting.increment_json_error_count(f"pathway JSON parse after validation failed: {e}")
                    setting.increment_error_count('pathway_parse_error_count', f"pathway JSON parse after validation failed: {e}")
                    continue

        df_output = pd.DataFrame(
            rows[:target_count],
            columns=["pathway", "succeeded_iterate", "var_iterate_pathway", "llm", "temp", "num_pathway", "filled"],
        )
        repeat_counts = df_output["pathway"].value_counts().to_dict()
        df_output["pathway_repeat_count"] = df_output["pathway"].map(repeat_counts)
        df_output.to_csv(os.path.join(output_dir, output_name + "_pathway.csv"), index=False)

    ############################################################################################################
    # Execution
    asyncio.run(batch_prompt_pathway(input_question, var_num_pathway, var_llm, var_maxtoken, var_temp, var_iterate_pathway, var_max_concurrent_call, var_max_attempt, output_name, output_dir))
