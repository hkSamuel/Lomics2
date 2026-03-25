# Lomics - Large Language Models for Omics Sciences
# Version v1.1
# Date: January 5, 2025
# Developers: Chun-Ka WONG, Ali CHOO, Eugene C. C. CHENG, Wing-Chun SAN, Kelvin Chak-Kong CHENG, Hung-Fat TSE, Jason Wing-Hon WONG (The University of Hong Kong)
# http://www.lomics.ai

import os
from glob import glob
import pandas as pd
from litellm import acompletion
import asyncio
from json import loads, load
from pydantic import BaseModel, Field, create_model
import re
import csv
import argparse

os.environ["OPENAI_API_KEY"] = "api"
os.environ["OPENAI_API_BASE"] = "https://openrouter.ai/api/v1"
var_llm = "openrouter/meta-llama/llama-4-maverick"
var_maxtoken = 20000
var_temp = 0
var_max_concurrent_call = 20
var_num_pathway = 100
var_iterate_pathway = 10
var_num_gene = 100
var_iterate_gene = 3
var_max_attempt = 3
lomics = "v2.1"
current_dir = os.path.dirname(os.path.abspath(__file__))
path_hgnc = os.path.join(current_dir, "resources", "hgnc_20260206.json")
path_gmt = os.path.join(current_dir, "resources", "cleaned_wikipathways.gmt")

# Default output directory: user's Downloads folder
path_default_output = os.path.join(os.path.expanduser("~"), "Downloads")

# Runtime tracking
total_tokens_used = 0
run_start_time = None
json_error_count = 0
llm_call_error_count = 0
llm_jsondecode_expecting_value_error_count = 0
llm_empty_content_error_count = 0
llm_finish_reason_length_count = 0
pathway_prompt_error_count = 0
pathway_schema_error_count = 0
pathway_parse_error_count = 0
report_prompt_error_count = 0
report_parse_error_count = 0

def increment_json_error_count(context=None):
    global json_error_count
    json_error_count += 1
    if context:
        print(f"JSON error recorded ({json_error_count}): {context}")

def increment_error_count(counter_name, context=None):
    if counter_name == 'llm_call_error_count':
        global llm_call_error_count
        llm_call_error_count += 1
        count = llm_call_error_count
    elif counter_name == 'llm_jsondecode_expecting_value_error_count':
        global llm_jsondecode_expecting_value_error_count
        llm_jsondecode_expecting_value_error_count += 1
        count = llm_jsondecode_expecting_value_error_count
    elif counter_name == 'llm_empty_content_error_count':
        global llm_empty_content_error_count
        llm_empty_content_error_count += 1
        count = llm_empty_content_error_count
    elif counter_name == 'llm_finish_reason_length_count':
        global llm_finish_reason_length_count
        llm_finish_reason_length_count += 1
        count = llm_finish_reason_length_count
    elif counter_name == 'pathway_prompt_error_count':
        global pathway_prompt_error_count
        pathway_prompt_error_count += 1
        count = pathway_prompt_error_count
    elif counter_name == 'pathway_schema_error_count':
        global pathway_schema_error_count
        pathway_schema_error_count += 1
        count = pathway_schema_error_count
    elif counter_name == 'pathway_parse_error_count':
        global pathway_parse_error_count
        pathway_parse_error_count += 1
        count = pathway_parse_error_count
    elif counter_name == 'report_prompt_error_count':
        global report_prompt_error_count
        report_prompt_error_count += 1
        count = report_prompt_error_count
    elif counter_name == 'report_parse_error_count':
        global report_parse_error_count
        report_parse_error_count += 1
        count = report_parse_error_count
    else:
        return

    if context:
        print(f"Error recorded [{counter_name}] ({count}): {context}")

def load_gmt(path_gmt):
    gmt_dict = {}
    try:
        with open(path_gmt, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > 2:
                    pathway = parts[0]
                    genes = parts[2:]  # skip description
                    gmt_dict[pathway] = genes
    except FileNotFoundError:
        print(f"Warning: {path_gmt} not found. Proceeding without GMT data.")
    return gmt_dict

async def llm_call(prompt, var_llm, var_maxtoken, var_temp, limit):
    try:
        async with limit:
            response = await acompletion(
                model=var_llm,
                messages=[{
                    "role": "user",
                    "content": prompt,
                }],
                max_tokens=var_maxtoken,
                temperature=var_temp,
            )
            # Attempt to extract token usage from the response and accumulate
            try:
                usage = None
                if isinstance(response, dict):
                    usage = response.get('usage')
                else:
                    usage = getattr(response, 'usage', None)

                if usage is not None:
                    total = None
                    if isinstance(usage, dict):
                        total = usage.get('total_tokens') or (usage.get('prompt_tokens', 0) + usage.get('completion_tokens', 0))
                    else:
                        total = getattr(usage, 'total_tokens', None) or (getattr(usage, 'prompt_tokens', 0) + getattr(usage, 'completion_tokens', 0))

                    if total:
                        global total_tokens_used
                        total_tokens_used += int(total)
            except Exception:
                pass

            print("LLM response:")
            print(response)
            try:
                first_choice = response['choices'][0] if isinstance(response, dict) else getattr(response, 'choices')[0]
                finish_reason = first_choice.get('finish_reason') if isinstance(first_choice, dict) else getattr(first_choice, 'finish_reason', None)
                if str(finish_reason).lower() == 'length':
                    increment_error_count('llm_finish_reason_length_count', f"finish_reason=length for model {var_llm}")
            except Exception:
                pass
            # Extract message content as before
            try:
                content = response['choices'][0]['message']['content'] if isinstance(response, dict) else getattr(response, 'choices')[0].message.content
            except Exception:
                content = str(response)

            if content is None or (isinstance(content, str) and content.strip() == ""):
                increment_error_count('llm_empty_content_error_count', f"Empty/None content for model {var_llm}")

            print("LLM message:")
            print(content)
            # Try to extract JSON object from content
            try:
                content_json = content[content.find('{'):content.rfind('}')+1]
            except Exception:
                content_json = content
            print("LLM message updated:")
            print(content_json)
            return content_json
    except Exception as e:
        err_text = str(e)
        if ('JSONDecodeError: Expecting value' in err_text) or ('Expecting value: line' in err_text):
            increment_error_count('llm_jsondecode_expecting_value_error_count', err_text)
        increment_error_count('llm_call_error_count', str(e))
        print(e)
        return None