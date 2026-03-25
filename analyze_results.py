#!/usr/bin/env python3
"""
Analyze LOMICS pathway results
- Calculate repeated pathways within same LLM iteration
- Calculate how many pathways match WikiPathways reference
- Generate comprehensive metrics CSV report
"""

import pandas as pd
import os
import argparse


def truncate_path(full_path):
    """Truncate path to start from 'Downloads' onwards."""
    if 'Downloads' in full_path:
        idx = full_path.index('Downloads')
        return full_path[idx:]
    return full_path


def discover_result_sets(search_root):
    """Find all <name>_pathway.csv files and their matching <name>_report.csv files."""
    pairs = []
    for root, _, files in os.walk(search_root):
        for filename in files:
            if not filename.endswith("_pathway.csv"):
                continue
            outputname = filename[:-len("_pathway.csv")]
            pathway_csv = os.path.join(root, filename)
            report_csv = os.path.join(root, f"{outputname}_report.csv")
            if os.path.exists(report_csv):
                pairs.append((outputname, pathway_csv, report_csv, os.path.join(root, f"{outputname}_analysis.csv")))
    return sorted(pairs, key=lambda x: x[0])


def load_wikipathways(gmt_path):
    """Load WikiPathways reference pathways from GMT file"""
    wikipathways = set()
    try:
        with open(gmt_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts:
                    wikipathways.add(parts[0])
        return wikipathways
    except Exception as e:
        print(f"Error loading WikiPathways: {e}")
        return set()


def parse_report_header(report_path):
    """Parse the header comment line from report.csv if it exists
    
    Handles multiple versions of report headers:
    - Older versions: json_errors, total_tokens, elapsed_seconds, etc.
    - Newer versions: generated_at, detailed error counters, etc.
    
    Returns a dict with all parsed key=value pairs, with graceful defaults for missing fields.
    """
    metadata = {}
    try:
        with open(report_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#'):
                # Remove "# Report generated: " prefix and parse remaining key=value pairs
                header_content = first_line.replace('# Report generated: ', '')
                
                # Split on ", " but be careful about values that might contain commas or spaces
                # Use a simple split approach: split on ', ' and handle each part
                parts = header_content.split(', ')
                
                for part in parts:
                    if '=' not in part:
                        continue
                    # Split only on first '=' to handle values with '=' in them
                    key, val = part.split('=', 1)
                    key = key.strip()
                    val = val.strip()
                    metadata[key] = val
                
                # For backward compatibility, handle old 'json_errors' as 'json_errors_total' if needed
                if 'json_errors' in metadata and 'json_errors_total' not in metadata:
                    metadata['json_errors_total'] = metadata['json_errors']
                    
    except Exception as e:
        print(f"Warning: Could not parse report header: {e}")
    
    return metadata


def analyze_pathways(pathway_csv, report_csv, wikipathways_gmt, output_csv):
    """
    Main analysis function
    """
    print("="*80)
    print("LOMICS Pathway Analysis")
    print("="*80)
    
    # Load data
    print(f"\nLoading pathway data from: {pathway_csv}")
    df_pathway = pd.read_csv(pathway_csv)
    
    print(f"Loading report data from: {report_csv}")
    # Read report, handling potential header comment
    with open(report_csv, 'r', encoding='utf-8') as f:
        first_line = f.readline()
        if first_line.startswith('#'):
            # Skip header line
            df_report = pd.read_csv(report_csv, skiprows=1)
            metadata = parse_report_header(report_csv)
        else:
            # No header, read normally
            df_report = pd.read_csv(report_csv)
            metadata = {}

    # Capture LLM model info from existing CSV columns when available
    llm_models = set()
    if 'llm' in df_pathway.columns:
        llm_models.update([str(x) for x in df_pathway['llm'].dropna().unique().tolist()])
    if 'llm' in df_report.columns:
        llm_models.update([str(x) for x in df_report['llm'].dropna().unique().tolist()])
    llm_models_str = ';'.join(sorted(llm_models)) if llm_models else 'unknown'
    
    print(f"Loading WikiPathways reference from: {wikipathways_gmt}")
    wikipathways_ref = load_wikipathways(wikipathways_gmt)
    
    print(f"\n{'='*80}")
    print("DATA SUMMARY")
    print(f"{'='*80}")
    print(f"Total pathway entries: {len(df_pathway)}")
    print(f"Unique pathways: {df_pathway['pathway'].nunique()}")
    print(f"WikiPathways reference count: {len(wikipathways_ref)}")
    print(f"Number of iterations: {df_pathway['succeeded_iterate'].max() + 1 if 'succeeded_iterate' in df_pathway.columns else 'N/A'}")
    print(f"LLM model(s): {llm_models_str}")
    
    # Analysis results container
    results = []
    
    # ========================================================================
    # ANALYSIS 1: Repeated pathways within same LLM iteration
    # ========================================================================
    print(f"\n{'='*80}")
    print("ANALYSIS 1: Repeated Pathways within Same LLM Call/Iteration")
    print(f"{'='*80}")
    
    if 'succeeded_iterate' in df_pathway.columns:
        # Group by iteration and count duplicates
        iterations = df_pathway['succeeded_iterate'].unique()
        
        total_within_iteration_duplicates = 0
        iteration_stats = []
        
        for iteration in sorted(iterations):
            df_iter = df_pathway[df_pathway['succeeded_iterate'] == iteration]
            pathway_counts = df_iter['pathway'].value_counts()
            duplicates = pathway_counts[pathway_counts > 1]
            
            num_duplicated_pathways = len(duplicates)
            num_duplicate_entries = (duplicates - 1).sum()  # Extra occurrences
            total_within_iteration_duplicates += num_duplicate_entries
            
            iteration_stats.append({
                'iteration': iteration,
                'total_pathways': len(df_iter),
                'unique_pathways': df_iter['pathway'].nunique(),
                'duplicated_pathway_names': num_duplicated_pathways,
                'total_duplicate_entries': num_duplicate_entries
            })
            
            if num_duplicated_pathways > 0:
                print(f"\nIteration {iteration}:")
                print(f"  Total pathways: {len(df_iter)}")
                print(f"  Unique pathways: {df_iter['pathway'].nunique()}")
                print(f"  Duplicated pathway names: {num_duplicated_pathways}")
                print(f"  Duplicated pathways:")
                for pathway, count in duplicates.items():
                    print(f"    - '{pathway}': {count} times")
        
        results.append({
            'metric': 'total_within_iteration_duplicate_entries',
            'value': total_within_iteration_duplicates,
            'description': 'Total number of duplicate pathway entries within the same LLM iteration'
        })
        
        results.append({
            'metric': 'avg_duplicates_per_iteration',
            'value': round(total_within_iteration_duplicates / len(iterations), 2),
            'description': 'Average duplicate entries per iteration'
        })
    else:
        print("No 'succeeded_iterate' column found - cannot analyze within-iteration duplicates")
    
    # ========================================================================
    # ANALYSIS 2: Overall pathway statistics
    # ========================================================================
    print(f"\n{'='*80}")
    print("ANALYSIS 2: Overall Pathway Statistics")
    print(f"{'='*80}")
    
    # Get top pathways
    pathway_counts = df_pathway['pathway'].value_counts()
    
    results.append({
        'metric': 'total_pathway_entries',
        'value': len(df_pathway),
        'description': 'Total pathway entries across all iterations'
    })
    
    results.append({
        'metric': 'unique_pathways',
        'value': df_pathway['pathway'].nunique(),
        'description': 'Number of unique pathway names'
    })

    results.append({
        'metric': 'llm_models',
        'value': llm_models_str,
        'description': 'LLM model name(s) extracted from llm column in pathway/report CSVs'
    })
    
    # Calculate how many of ALL generated pathways are in WikiPathways
    pathways_in_wiki = set(df_pathway['pathway'].unique()) & wikipathways_ref
    pathways_not_in_wiki = set(df_pathway['pathway'].unique()) - wikipathways_ref
    
    results.append({
        'metric': 'unique_pathways_in_wikipathways',
        'value': len(pathways_in_wiki),
        'description': 'Number of unique generated pathways found in WikiPathways reference'
    })
    
    results.append({
        'metric': 'unique_pathways_not_in_wikipathways',
        'value': len(pathways_not_in_wiki),
        'description': 'Number of unique generated pathways NOT found in WikiPathways reference'
    })
    
    results.append({
        'metric': 'pct_unique_pathways_in_wikipathways',
        'value': round(len(pathways_in_wiki) / df_pathway['pathway'].nunique() * 100, 2),
        'description': 'Percentage of unique pathways that match WikiPathways'
    })
    
    print(f"\nUnique pathways in WikiPathways: {len(pathways_in_wiki)}")
    print(f"Unique pathways NOT in WikiPathways: {len(pathways_not_in_wiki)}")
    print(f"Percentage match: {round(len(pathways_in_wiki) / df_pathway['pathway'].nunique() * 100, 2)}%")
    
    if pathways_not_in_wiki:
        print(f"\nPathways NOT in WikiPathways:")
        for p in sorted(pathways_not_in_wiki):
            print(f"  - {p}")
    
    # ========================================================================
    # ANALYSIS 3: Top pathways from report.csv (already preselected there)
    # ========================================================================
    print(f"\n{'='*80}")
    print("ANALYSIS 3: Top Pathways from report.csv")
    print(f"{'='*80}")

    # Prefer report.csv pathway list (already top pathways in pipeline output).
    # Keep first occurrence order, then take up to 100 rows.
    if 'pathway' in df_report.columns:
        top_candidates = [str(x).strip() for x in df_report['pathway'].dropna().tolist() if str(x).strip()]
        top_100_ordered = list(dict.fromkeys(top_candidates))[:100]
    else:
        print("Warning: 'pathway' column missing in report.csv, falling back to pathway.csv frequency ranking")
        top_100_ordered = pathway_counts.head(100).index.tolist()

    top_n = len(top_100_ordered)
    top_100_pathways = set(top_100_ordered)
    top_100_in_wiki = top_100_pathways & wikipathways_ref
    top_100_not_in_wiki = top_100_pathways - wikipathways_ref
    
    results.append({
        'metric': 'top_100_pathways_in_wikipathways',
        'value': len(top_100_in_wiki),
        'description': 'Number of top pathways listed in report.csv (up to 100) found in WikiPathways'
    })
    
    results.append({
        'metric': 'top_100_pathways_not_in_wikipathways',
        'value': len(top_100_not_in_wiki),
        'description': 'Number of top pathways listed in report.csv (up to 100) NOT found in WikiPathways'
    })
    
    print(f"\nTop pathways considered from report.csv: {top_n}")
    print(f"In WikiPathways: {len(top_100_in_wiki)}")
    print(f"NOT in WikiPathways: {len(top_100_not_in_wiki)}")
    print(f"Percentage match: {round((len(top_100_in_wiki) / top_n * 100), 2) if top_n > 0 else 0.0}%")
    
    if top_100_not_in_wiki:
        print(f"\nTop pathways from report.csv NOT in WikiPathways:")
        for p in top_100_ordered:
            if p in top_100_not_in_wiki:
                print(f"  - {p}")
    
    # ========================================================================
    # ANALYSIS 4: Repeat count statistics
    # ========================================================================
    print(f"\n{'='*80}")
    print("ANALYSIS 4: Pathway Repetition Statistics")
    print(f"{'='*80}")
    
    max_repeat = pathway_counts.max()
    avg_repeat = pathway_counts.mean()
    median_repeat = pathway_counts.median()
    
    results.append({
        'metric': 'max_pathway_repeat_count',
        'value': max_repeat,
        'description': 'Maximum times a single pathway was repeated'
    })
    
    results.append({
        'metric': 'avg_pathway_repeat_count',
        'value': round(avg_repeat, 2),
        'description': 'Average repeat count per pathway'
    })
    
    results.append({
        'metric': 'median_pathway_repeat_count',
        'value': median_repeat,
        'description': 'Median repeat count per pathway'
    })
    
    print(f"\nMax repeat count: {max_repeat}")
    print(f"Average repeat count: {round(avg_repeat, 2)}")
    print(f"Median repeat count: {median_repeat}")
    
    print(f"\nTop 10 most repeated pathways:")
    for i, (pathway, count) in enumerate(pathway_counts.head(10).items(), 1):
        in_wiki = "✓" if pathway in wikipathways_ref else "✗"
        print(f"  {i}. {pathway}: {count} times [{in_wiki} WikiPathways]")
    
    # ========================================================================
    # Add metadata from report if available
    # ========================================================================
    if metadata:
        for key, value in metadata.items():
            results.append({
                'metric': f'metadata_{key}',
                'value': value,
                'description': f'Metadata from report: {key}'
            })
    
    # ========================================================================
    # Output CSV listing pathways NOT in reference
    # ========================================================================
    not_in_ref_rows = []
    
    # Add unique generated pathways not in reference
    for p in sorted(pathways_not_in_wiki):
        repeat_count = pathway_counts.get(p, 0)
        not_in_ref_rows.append({
            'pathway': p,
            'type': 'unique_generated',
            'repeat_count': repeat_count
        })
    
    # Add top-100 pathways not in reference
    for p in top_100_ordered:
        if p in top_100_not_in_wiki:
            repeat_count = pathway_counts.get(p, 0)
            # Check if already added as unique_generated
            if not any(r['pathway'] == p and r['type'] == 'unique_generated' for r in not_in_ref_rows):
                not_in_ref_rows.append({
                    'pathway': p,
                    'type': 'top_100_from_report',
                    'repeat_count': repeat_count
                })
    
    if not_in_ref_rows:
        df_not_in_ref = pd.DataFrame(not_in_ref_rows)
        not_in_ref_csv = output_csv.replace('_analysis.csv', '_not_in_wikipathways.csv')
        df_not_in_ref.to_csv(not_in_ref_csv, index=False)
        print(f"Pathways not in reference saved to: {not_in_ref_csv}")
    
    # ========================================================================
    # Save results to CSV
    # ========================================================================
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_csv, index=False)
    
    print(f"\n{'='*80}")
    print(f"Analysis complete! Results saved to: {output_csv}")
    print(f"{'='*80}\n")
    
    return df_results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze LOMICS pathway results',
        epilog=(
            'Single set: python analyze_results.py --outputname myrun --outputdir outputs/my_folder\n'
            'Batch mode: python analyze_results.py --input-root /path/to/results'
        )
    )
    parser.add_argument(
        '--outputname',
        type=str,
        required=False,
        help='Output name (same as used in run.py)'
    )
    parser.add_argument(
        '--outputdir',
        type=str,
        default='outputs',
        help='Output directory (default: outputs)'
    )
    parser.add_argument(
        '--wikipathways',
        type=str,
        default='resources/cleaned_wikipathways.gmt',
        help='Path to WikiPathways GMT file (default: resources/cleaned_wikipathways.gmt)'
    )
    parser.add_argument(
        '--input-root',
        type=str,
        default=None,
        help='Batch mode: recursively find all *_pathway.csv + *_report.csv under this folder'
    )
    parser.add_argument(
        '--batch-summary',
        type=str,
        default=None,
        help='Optional path for combined summary CSV in batch mode'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.wikipathways):
        print(f"Error: WikiPathways file not found: {args.wikipathways}")
        return

    # Batch mode: recursively discover and process many runs
    if args.input_root:
        pairs = discover_result_sets(args.input_root)
        if not pairs:
            print(f"Error: No matched result sets found under: {args.input_root}")
            print("Expected pairs like: <name>_pathway.csv and <name>_report.csv in the same folder")
            return

        summary_rows = []
        print(f"Found {len(pairs)} result set(s). Starting batch analysis...\n")
        for outputname, pathway_csv, report_csv, output_csv in pairs:
            try:
                df = analyze_pathways(pathway_csv, report_csv, args.wikipathways, output_csv)
                metric_map = dict(zip(df['metric'], df['value']))
                summary_rows.append({
                    'outputname': outputname,
                    'pathway_csv': truncate_path(pathway_csv),
                    'report_csv': truncate_path(report_csv),
                    'analysis_csv': truncate_path(output_csv),
                    'llm_models': metric_map.get('llm_models'),
                    'unique_pathways': metric_map.get('unique_pathways'),
                    'unique_pathways_in_wikipathways': metric_map.get('unique_pathways_in_wikipathways'),
                    'pct_unique_pathways_in_wikipathways': metric_map.get('pct_unique_pathways_in_wikipathways'),
                    'top_100_pathways_in_wikipathways': metric_map.get('top_100_pathways_in_wikipathways'),
                    'total_within_iteration_duplicate_entries': metric_map.get('total_within_iteration_duplicate_entries'),
                })
            except Exception as e:
                print(f"ERROR processing {outputname}: {e}")

        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            summary_path = args.batch_summary or os.path.join(args.input_root, 'batch_analysis_summary.csv')
            summary_df.to_csv(summary_path, index=False)
            print(f"\nBatch summary saved to: {summary_path}")
        return

    # Single set mode
    if not args.outputname:
        print("Error: --outputname is required in single-set mode, or use --input-root for batch mode")
        return

    pathway_csv = os.path.join(args.outputdir, f"{args.outputname}_pathway.csv")
    report_csv = os.path.join(args.outputdir, f"{args.outputname}_report.csv")
    output_csv = os.path.join(args.outputdir, f"{args.outputname}_analysis.csv")

    if not os.path.exists(pathway_csv):
        print(f"Error: Pathway file not found: {pathway_csv}")
        return

    if not os.path.exists(report_csv):
        print(f"Error: Report file not found: {report_csv}")
        return

    analyze_pathways(pathway_csv, report_csv, args.wikipathways, output_csv)


if __name__ == "__main__":
    main()
