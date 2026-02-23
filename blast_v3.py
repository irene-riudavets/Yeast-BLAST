"""
BLAST XML Analysis Tool
A Python module for parsing and analyzing BLAST XML results.
Extracts key metrics from BLAST hits, classifies hit quality, identifies yeast
organisms, and generates detailed reports including TSV files and HTML visualizations.
Main Features:
- Parse BLAST XML output and extract alignment metrics
- Classify hit quality based on e-value, percent identity, and coverage
- Identify and track yeast species hits
- Generate detailed TSV reports with quality classifications
- Create HTML report with embedded visualizations
- Produce distribution plots for e-values, percent identity, and coverage
- Fallback to TSV output when matplotlib is unavailable
Key Functions:
- analyze(): Main entry point for BLAST result analysis
- parse_blast_results(): Parse BLAST XML and extract metrics
- extract_hit_info(): Extract relevant information from individual BLAST hits
- classify_hit_quality(): Classify hits as Good/OK/Poor based on thresholds
- is_yeast(): Identify if organism is a yeast species
- write_tsv_results(): Generate detailed results TSV
- write_summary(): Generate summary TSV with quality labels
- write_html_report(): Generate interactive HTML report
- create_plots(): Generate matplotlib visualizations or fallback TSV files
Dependencies:
- BioPython (Bio.Blast.NCBIXML)
- numpy (optional, for plotting)
- matplotlib (optional, for visualizations)
Usage:
    python blast_v3.py --xml <path_to_blast.xml> --out <output_dir> 
                       [--evalue <threshold>] [--top <n_species>]
Output Files:
- top_hits.tsv: Detailed results with all metrics
- top_hits_summary.tsv: Summary with quality classification
- blast_report.html: Interactive HTML report
- top_species.png: Bar chart of top species (if matplotlib available)
- distributions.png: Histograms of metrics (if matplotlib available)
- species_counts.tsv, evalues.tsv, percent_id.tsv, coverage.tsv: 
  Fallback data files (if matplotlib unavailable)
"""
import os
import re
import csv
import argparse
from collections import Counter
from Bio.Blast import NCBIXML


def extract_organism(title: str) -> str:
    """Extract organism name from BLAST hit title.

    For protein hits the organism is in square brackets, e.g. [Saccharomyces cerevisiae].
    For RNA/PDB hits there is no organism field; we return 'Unknown' in that case.
    """
    # Look for organism in square brackets (protein / standard format)
    m = re.search(r"\[([^\]]+)\]", title)
    if m:
        return m.group(1)
    # Alternative format with OS=
    m = re.search(r"OS=([^\[]+?)(?:\sGN=|$)", title)
    if m:
        return m.group(1).strip()
    return "Unknown"


def extract_pdb_id(hit_id: str, accession: str) -> str:
    """Extract PDB code from BLAST hit identifiers.

    Handles both blastp format  (pdb|4YBB|A)
    and blastn/gi format        (gi|2981298973|pdb|8ZGR|LA).
    Falls back to the accession field which BioPython stores as e.g. '8ZGR_LA'.
    """
    # Try to find pdb|XXXX| pattern in either id string
    for s in (hit_id or '', accession or ''):
        m = re.search(r'pdb\|([A-Za-z0-9]{4})\|', s, re.IGNORECASE)
        if m:
            return m.group(1).upper()
    # Fallback: accession is often "8ZGR_LA" — take the part before '_'
    if accession:
        return accession.split('_')[0].split('|')[0]
    return ''


def classify_hit_quality(evalue, pct_id, coverage, blast_type='blastp'):
    """Classify hit quality based on metrics.

    Uses different thresholds for blastn (RNA/DNA) vs blastp (protein).
    E-value of exactly 0 means better than float can represent — treated as best possible.
    """
    if pct_id is None or coverage is None or evalue is None:
        return 'Poor'

    # E-value 0 means vanishingly small — treat as best possible
    effective_evalue = evalue if evalue > 0 else 1e-300

    if blast_type == 'blastn':
        # RNA/nucleotide thresholds: identity is typically high; coverage can be low for rRNA
        if pct_id >= 85 and coverage >= 20 and effective_evalue < 1e-5:
            return 'Good'
        elif pct_id >= 70 and coverage >= 10 and effective_evalue < 1e-3:
            return 'OK'
        else:
            return 'Poor'
    else:
        # Protein thresholds (original)
        if pct_id >= 50 and coverage >= 70 and effective_evalue < 1e-5:
            return 'Good'
        elif pct_id >= 30 and coverage >= 50 and effective_evalue < 1e-3:
            return 'OK'
        else:
            return 'Poor'


def is_yeast(organism: str) -> bool:
    """Check if organism is a yeast species."""
    if not organism:
        return False
    yeast_keywords = [
        'saccharomyces', 'candida', 'pichia', 'komagataella', 'kluyveromyces',
        'schizosaccharomyces', 'yarrowia', 'debaryomyces', 'hansenula'
    ]
    org_lower = organism.lower()
    return any(keyword in org_lower for keyword in yeast_keywords)


def extract_hit_info(alignment, hsps, qlen, blast_type='blastp'):
    """Extract all relevant information from a BLAST hit.

    For blastn: aggregates coverage and identity across all HSPs (union of query
    spans for coverage; pooled identities/align-lengths for %identity).
    For blastp: uses the single best HSP.
    """
    title = alignment.title

    if blast_type == 'blastn':
        # Best (lowest) e-value across HSPs; treat 0 as best possible
        evalue = min((hsp.expect if hsp.expect > 0 else 1e-300) for hsp in hsps)
        # Highest bit-score
        bit_score = max(getattr(hsp, 'bits', 0) or 0 for hsp in hsps)

        # Union of query coordinate ranges to avoid double-counting overlaps
        ranges = sorted((min(hsp.query_start, hsp.query_end),
                         max(hsp.query_start, hsp.query_end)) for hsp in hsps)
        merged = []
        for start, end in ranges:
            if merged and start <= merged[-1][1]:
                merged[-1] = [merged[-1][0], max(merged[-1][1], end)]
            else:
                merged.append([start, end])
        union_span = sum(e - s for s, e in merged)

        total_identities = sum(getattr(hsp, 'identities', 0) or 0 for hsp in hsps)
        total_align_len  = sum(getattr(hsp, 'align_length', 0) or 0 for hsp in hsps)

        pct_id     = round(100.0 * total_identities / total_align_len, 2) if total_align_len else None
        coverage   = round(100.0 * union_span / qlen, 2) if qlen else None
        identities = total_identities
        align_len  = total_align_len

    else:
        # Single best HSP for protein BLAST
        hsp        = hsps[0]
        evalue     = hsp.expect
        bit_score  = getattr(hsp, 'bits', None)
        identities = getattr(hsp, 'identities', None)
        align_len  = getattr(hsp, 'align_length', None)
        pct_id     = round(100.0 * identities / align_len, 2) if (identities and align_len) else None
        coverage   = round(100.0 * align_len / qlen, 2)       if (align_len and qlen)       else None

    # --- PDB ID: use regex to extract 4-char code from hit_id or accession ---
    hit_id    = getattr(alignment, 'hit_id', '') or ''
    accession = getattr(alignment, 'accession', '') or ''
    pdb_id    = extract_pdb_id(hit_id, accession)

    organism = extract_organism(title)

    return {
        'title': title, 'evalue': evalue, 'bit_score': bit_score,
        'identities': identities, 'align_len': align_len, 'pct_id': pct_id,
        'coverage': coverage, 'organism': organism, 'pdb_id': pdb_id
    }


def parse_blast_results(xml_path: str, evalue_threshold: float = None):
    """Parse BLAST XML and extract key metrics.

    Auto-detects blastp vs blastn from <BlastOutput_program> and adjusts
    coverage/identity calculations accordingly.
    """
    rows = []
    species_counter = Counter()
    evalues, pct_ids, coverages = [], [], []

    with open(xml_path) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            # Detect BLAST program (BioPython stores it as blast_record.application)
            app = getattr(blast_record, 'application', 'blastp').lower()
            blast_type = 'blastn' if 'blastn' in app else 'blastp'

            query = blast_record.query
            qlen  = getattr(blast_record, 'query_length', None)

            best_alignment  = None
            best_hsps       = []
            best_evalue     = float('inf')

            best_yeast_alignment = None
            best_yeast_hsps      = []
            best_yeast_evalue    = float('inf')

            for alignment in blast_record.alignments:
                # Filter HSPs by e-value threshold if supplied
                valid_hsps = [hsp for hsp in alignment.hsps
                              if not evalue_threshold or hsp.expect <= evalue_threshold]
                if not valid_hsps:
                    continue

                # Representative e-value: minimum across all valid HSPs
                # (treat 0 as best possible for comparison)
                rep_evalue = min((h.expect if h.expect > 0 else 1e-300) for h in valid_hsps)

                if rep_evalue < best_evalue:
                    best_evalue    = rep_evalue
                    best_alignment = alignment
                    best_hsps      = valid_hsps

                organism = extract_organism(alignment.title)
                if is_yeast(organism) and rep_evalue < best_yeast_evalue:
                    best_yeast_evalue     = rep_evalue
                    best_yeast_alignment  = alignment
                    best_yeast_hsps       = valid_hsps

            # No hit
            if not best_alignment or not best_hsps:
                rows.append({
                    'query': query, 'title': None, 'evalue': None, 'bit_score': None,
                    'identities': None, 'align_len': None, 'pct_id': None,
                    'coverage': None, 'organism': None, 'pdb_id': None,
                    'yeast_hit': None, 'is_yeast': False, 'blast_type': blast_type
                })
                continue

            # For blastp pass only the single best HSP; for blastn pass all valid HSPs
            if blast_type == 'blastp':
                top_hsp = min(best_hsps, key=lambda h: h.expect)
                hsps_to_use = [top_hsp]
            else:
                hsps_to_use = best_hsps

            hit_info = extract_hit_info(best_alignment, hsps_to_use, qlen, blast_type)
            hit_info['query']      = query
            hit_info['blast_type'] = blast_type
            hit_info['is_yeast']   = is_yeast(hit_info['organism'])

            if not hit_info['is_yeast'] and best_yeast_alignment and best_yeast_hsps:
                if blast_type == 'blastp':
                    y_top = min(best_yeast_hsps, key=lambda h: h.expect)
                    y_hsps = [y_top]
                else:
                    y_hsps = best_yeast_hsps
                yeast_info = extract_hit_info(best_yeast_alignment, y_hsps, qlen, blast_type)
                hit_info['yeast_hit'] = yeast_info
            else:
                hit_info['yeast_hit'] = None

            species_counter[hit_info['organism']] += 1
            rows.append(hit_info)

            if hit_info['evalue']   is not None: evalues.append(hit_info['evalue'])
            if hit_info['pct_id']   is not None: pct_ids.append(hit_info['pct_id'])
            if hit_info['coverage'] is not None: coverages.append(hit_info['coverage'])

    return rows, species_counter, evalues, pct_ids, coverages


def write_tsv_results(rows, out_dir):
    """Write TSV results."""
    out_tsv = os.path.join(out_dir, 'top_hits.tsv')
    with open(out_tsv, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['Query', 'Top_hit_title', 'E_value', 'Bit_score', 'Identities', 
                        'Align_length', 'Percent_identity', 'Coverage_percent', 'Organism', 
                        'PDB_id', 'Is_Yeast'])
        for r in rows:
            writer.writerow([
                r['query'], r['title'], r['evalue'], r['bit_score'], 
                r['identities'], r['align_len'], r['pct_id'], r['coverage'], 
                r['organism'], r['pdb_id'], 'Yes' if r.get('is_yeast') else 'No'
            ])


def write_summary(rows, out_dir):
    """Write summary TSV with quality labels and yeast information."""
    summary_tsv = os.path.join(out_dir, 'top_hits_summary.tsv')
    with open(summary_tsv, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['Query', 'Top_hit', 'Percent_identity', 'Coverage_percent', 
                        'E_value', 'Bit_score', 'PDB_id', 'Organism', 'Quality', 'Is_Yeast',
                        'Best_Yeast_Hit', 'Yeast_Organism', 'Yeast_Pct_Id', 'Yeast_E_value'])
        for r in rows:
            quality = classify_hit_quality(r['evalue'], r['pct_id'], r['coverage'], r.get('blast_type', 'blastp'))
            yeast_hit = r.get('yeast_hit')
            
            # Yeast alternative info if available
            yeast_title = ''
            yeast_org = ''
            yeast_pct = ''
            yeast_eval = ''
            if yeast_hit:
                yeast_title = simplify_title(yeast_hit['title'])
                yeast_org = yeast_hit['organism']
                yeast_pct = yeast_hit['pct_id']
                yeast_eval = yeast_hit['evalue']
            
            writer.writerow([
                r['query'], r['title'], r['pct_id'], r['coverage'], 
                r['evalue'], r['bit_score'], r['pdb_id'], r['organism'], quality,
                'Yes' if r.get('is_yeast') else 'No',
                yeast_title, yeast_org, yeast_pct, yeast_eval
            ])


def simplify_title(title):
    """Extract description from a BLAST hit title.

    Handles both blastp format  ('pdb|4YBB|A Chain A, protein [Organism]')
    and blastn/gi format        ('gi|...|pdb|8ZGR|LA Chain LA, 25S rRNA (3393-MER)').
    Strips the identifier prefix and returns only the chain description from the
    first entry (before the first '>').
    """
    if not title:
        return ''
    # Take only the first hit entry (before '>gi|...' or '>pdb|...' alternatives)
    first_hit = title.split('>')[0].strip()
    # Remove leading pdb|XXXX|chain or gi|...|pdb|XXXX|chain identifier
    first_hit = re.sub(r'^(gi\|\d+\|)?pdb\|[^|]+\|\S+\s*', '', first_hit)
    return first_hit.strip()


def simplify_pdb_id(pdb_id):
    """Return the PDB ID as-is (already cleaned to 4-char code by extract_pdb_id)."""
    if not pdb_id:
        return ''
    # Strip any residual prefixes just in case
    clean_id = re.sub(r'^(pdb\|?|gi\|[^|]+\|)', '', pdb_id, flags=re.IGNORECASE)
    return clean_id.split('_')[0].split('|')[0].upper()


def write_html_report(rows, out_dir):
    """Generate HTML report with embedded visualizations."""
    report_html = os.path.join(out_dir, 'blast_report.html')
    
    with open(report_html, 'w') as fh:
        fh.write('''<!doctype html>
<html><head><meta charset="utf-8"><title>BLAST Report</title>
<style>
    body{font-family:Arial,sans-serif;margin:20px}
    table{border-collapse:collapse;width:100%;margin-top:20px}
    th,td{border:1px solid #ccc;padding:6px;text-align:left;font-size:14px}
    .Good{background:#d4f7d4}
    .OK{background:#fff4cc}
    .Poor{background:#ffd6d6}
    .NonYeast{background:#cbc3e3}
    img{max-width:100%;height:auto;margin:20px 0}
    .yeast-alt{font-size:0.9em;color:#666;margin-top:4px;padding-top:4px;border-top:1px dashed #999}
</style></head><body>
<h1>BLAST Top-Hits Summary</h1>
''')
        
        # Embed plots at the beginning - check if files exist
        plot_files = ['top_species.png', 'distributions.png']
        for img in plot_files:
            img_path = os.path.join(out_dir, img)
            if os.path.exists(img_path):
                fh.write(f'<div><img src="{img}" alt="{img}"></div>\n')
        
        # Write table
        fh.write('<h2>Top hits</h2>\n')
        fh.write('<p><strong>Note:</strong> Rows highlighted in light purple indicate non-yeast top hits. ')
        fh.write('Alternative yeast hits are shown below when available.</p>\n')
        fh.write('<table><thead><tr>')
        fh.write('<th>Query</th><th>Top hit</th><th>% id</th><th>Coverage</th>')
        fh.write('<th>E-value</th><th>Bit score</th><th>PDB</th><th>Organism</th><th>Quality</th>')
        fh.write('</tr></thead><tbody>\n')
        
        for r in rows:
            quality = classify_hit_quality(r['evalue'], r['pct_id'], r['coverage'], r.get('blast_type', 'blastp'))
            clean_title = simplify_title(r['title'])
            clean_pdb = simplify_pdb_id(r['pdb_id'])
            
            # Determine row class
            row_class = quality
            if not r['is_yeast']:
                row_class = 'NonYeast'
            
            fh.write(f'<tr class="{row_class}">')
            fh.write(f'<td>{r["query"]}</td>')
            
            # Top hit cell - may include yeast alternative
            hit_cell = clean_title
            if not r['is_yeast'] and r['yeast_hit']:
                yh = r['yeast_hit']
                yeast_title = simplify_title(yh['title'])
                yeast_pdb = simplify_pdb_id(yh['pdb_id'])
                yeast_quality = classify_hit_quality(yh['evalue'], yh['pct_id'], yh['coverage'], r.get('blast_type', 'blastp'))
                hit_cell += f'<div class="yeast-alt"><strong>Best yeast:</strong> {yeast_title} '
                hit_cell += f'({yh["organism"]}, {yh["pct_id"]}% id, E={yh["evalue"]}, PDB: {yeast_pdb})</div>'
            
            fh.write(f'<td>{hit_cell}</td>')
            fh.write(f'<td>{r["pct_id"] or ""}</td><td>{r["coverage"] or ""}</td>')
            fh.write(f'<td>{r["evalue"] or ""}</td><td>{r["bit_score"] or ""}</td>')
            fh.write(f'<td>{clean_pdb}</td><td>{r["organism"] or ""}</td>')
            fh.write(f'<td>{quality}</td></tr>\n')
        
        fh.write('</tbody></table>')
        fh.write('<p><strong>Legend:</strong> ')
        fh.write('<span class="Good" style="padding:3px 8px">Good</span> ')
        fh.write('<span class="OK" style="padding:3px 8px">OK</span> ')
        fh.write('<span class="Poor" style="padding:3px 8px">Poor</span> ')
        fh.write('<span class="NonYeast" style="padding:3px 8px">Non-yeast hit</span>')
        fh.write('</p>')
        fh.write('</body></html>')


def create_plots(species_counter, evalues, pct_ids, coverages, out_dir, top_n_species):
    """Create visualization plots."""
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        
        # Species bar chart
        most_common = species_counter.most_common(top_n_species)
        if most_common:
            labels, counts = zip(*most_common)
            plt.figure(figsize=(10, 6))
            plt.bar(range(len(labels)), counts, color='C0')
            plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
            plt.ylabel('Count')
            plt.title('Top species in top hits')
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, 'top_species.png'))
            plt.close()
        
        # Distribution plots
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        if evalues:
            ev = np.array(evalues, dtype=float)
            ev[ev <= 0] = 1e-200
            axes[0].hist(-np.log10(ev), bins=40, color='C1')
            axes[0].set_xlabel('-log10(E-value)')
            axes[0].set_ylabel('Count')
        
        if pct_ids:
            axes[1].hist(pct_ids, bins=40, color='C2')
            axes[1].set_xlabel('Percent identity')
        
        if coverages:
            axes[2].hist(coverages, bins=40, color='C3')
            axes[2].set_xlabel('Coverage (%)')
        
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'distributions.png'))
        plt.close()
        
    except ImportError:
        # Fallback: write data to TSV files
        _write_distribution_tsvs(species_counter, evalues, pct_ids, coverages, out_dir)


def _write_distribution_tsvs(species_counter, evalues, pct_ids, coverages, out_dir):
    """Write distribution data as TSV when plotting is unavailable."""
    # Species counts
    with open(os.path.join(out_dir, 'species_counts.tsv'), 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['Organism', 'Count'])
        for org, cnt in species_counter.most_common():
            writer.writerow([org, cnt])
    
    # E-values
    if evalues:
        with open(os.path.join(out_dir, 'evalues.tsv'), 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(['E_value'])
            for v in evalues:
                writer.writerow([v])
    
    # Percent identities
    if pct_ids:
        with open(os.path.join(out_dir, 'percent_id.tsv'), 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(['Percent_identity'])
            for v in pct_ids:
                writer.writerow([v])
    
    # Coverage
    if coverages:
        with open(os.path.join(out_dir, 'coverage.tsv'), 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(['Coverage'])
            for v in coverages:
                writer.writerow([v])


def analyze(xml_path: str, out_dir: str, evalue_threshold: float = None, top_n_species: int = 10):
    """Main analysis function."""
    os.makedirs(out_dir, exist_ok=True)
    
    # Parse BLAST results
    rows, species_counter, evalues, pct_ids, coverages = parse_blast_results(xml_path, evalue_threshold)
    
    # Create plots FIRST (before HTML report references them)
    create_plots(species_counter, evalues, pct_ids, coverages, out_dir, top_n_species)
    
    # Write outputs
    write_tsv_results(rows, out_dir)
    write_summary(rows, out_dir)
    write_html_report(rows, out_dir)
    
    print(f"Analysis complete. Results saved to {out_dir}")
    
    # Print summary stats
    yeast_count = sum(1 for r in rows if r.get('is_yeast'))
    non_yeast_count = sum(1 for r in rows if not r.get('is_yeast') and r['title'] is not None)
    print(f"  Total queries: {len(rows)}")
    print(f"  Yeast top hits: {yeast_count}")
    print(f"  Non-yeast top hits: {non_yeast_count}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze BLAST XML: produce top-hits TSV and species plot')
    parser.add_argument('--xml', default='/home/joaquin_arino/Documents/Ribosome/Blast/protein_blast.xml',
                       help='Path to BLAST XML file')
    parser.add_argument('--out', default='/home/joaquin_arino/Documents/Ribosome/Blast/Results',
                       help='Output directory')
    parser.add_argument('--evalue', type=float, default=None, 
                       help='Optional e-value threshold to filter HSPs')
    parser.add_argument('--top', type=int, default=20, 
                       help='Top N species to plot')
    
    args = parser.parse_args()
    analyze(args.xml, args.out, args.evalue, args.top)