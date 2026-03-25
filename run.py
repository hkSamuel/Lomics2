from pathway import lomics_pathway
from gene import lomics_gene
from report import lomics_report
from setting import *
import setting
import time

def main():
    parser = argparse.ArgumentParser(description='Lomics: Large Language Models for Omic Studies v1.0', epilog='python run.py --question "scientific question" --outputname "output file name" --outputdir "output file directory"')
    parser.add_argument('--question', type=str, required=True, help='Scientific question')
    parser.add_argument('--outputname', type=str, required=True, help='Output file name')
    parser.add_argument('--outputdir', type=str, required=False, default=path_default_output, help='Output file directory (defaults to your Downloads folder)')
    args = parser.parse_args()
    input_question = args.question
    output_name = args.outputname
    output_dir = args.outputdir
    # record start time for overall execution
    setting.json_error_count = 0
    setting.llm_call_error_count = 0
    setting.llm_jsondecode_expecting_value_error_count = 0
    setting.llm_empty_content_error_count = 0
    setting.llm_finish_reason_length_count = 0
    setting.pathway_prompt_error_count = 0
    setting.pathway_schema_error_count = 0
    setting.pathway_parse_error_count = 0
    setting.report_prompt_error_count = 0
    setting.report_parse_error_count = 0
    setting.run_start_time = time.time()
    lomics_pathway(input_question, var_num_pathway, var_llm, var_maxtoken, var_temp, var_iterate_pathway, var_max_concurrent_call, var_max_attempt, output_name, output_dir)
    lomics_gene(input_question, var_llm, var_maxtoken, var_temp, var_iterate_gene, var_num_gene, var_max_concurrent_call, var_max_attempt, output_name, output_dir)
    lomics_report(input_question, var_llm, var_maxtoken, var_temp, var_max_attempt, output_name, output_dir)

if __name__ == "__main__":
    main()