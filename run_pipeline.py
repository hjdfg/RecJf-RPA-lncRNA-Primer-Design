"""
Main execution script for RecJf-RPA lncRNA Primer Design Pipeline

This script provides a complete workflow for:
1. Loading lncRNA target sequences
2. Designing primers and probes optimized for RPA
3. Validating designs for RecJf compatibility
4. Exporting results to CSV and detailed reports

Usage:
    python run_pipeline.py                          # Run all default targets
    python run_pipeline.py --target HOTAIR          # Run specific target
    python run_pipeline.py --target HOTAIR MALAT1   # Run multiple targets
    python run_pipeline.py --output results/        # Specify output directory
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

# Import pipeline modules
from primer_design import RPAPrimerDesigner, PrimerDesignResult
from sample_targets import get_lncrna_targets, get_target_info
from rpa_validator import validate_rpa_compatibility, print_validation_report


class RecJfRPAPipeline:
    """Main pipeline orchestrator for RecJf-RPA primer design."""
    
    def __init__(self, output_dir: str = "output"):
        """
        Initialize the pipeline.
        
        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.designer = RPAPrimerDesigner()
        self.results_summary = []
    
    def run_design_pipeline(self, target_names: list = None):
        """
        Run the complete design pipeline for specified targets.
        
        Args:
            target_names: List of target lncRNA names. If None, run all.
        """
        targets = get_lncrna_targets()
        target_info = get_target_info()
        
        # Determine which targets to process
        if target_names is None:
            targets_to_process = list(targets.keys())
        else:
            targets_to_process = [t for t in target_names if t in targets]
        
        if not targets_to_process:
            print("❌ No valid targets found!")
            return
        
        print(f"\n{'='*70}")
        print(f"  RecJf-RPA lncRNA Primer Design Pipeline")
        print(f"  {'='*70}")
        print(f"  Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Targets to process: {len(targets_to_process)}")
        print(f"  Output directory: {self.output_dir.absolute()}")
        print(f"{'='*70}\n")
        
        # Process each target
        for i, target_name in enumerate(targets_to_process, 1):
            print(f"\n[{i}/{len(targets_to_process)}] Processing: {target_name}")
            print("-" * 70)
            
            try:
                # Get sequence and info
                sequence = targets[target_name]
                info = target_info[target_name]
                
                # Design primers
                print(f"  ✓ Designing primers for {target_name}...")
                results = self.designer.design_primers(
                    sequence=sequence,
                    target_name=target_name,
                    primer_length_range=(17, 25),
                    probe_length=50,
                    gc_range=(0.40, 0.60),
                    max_tm_diff=5.0
                )
                
                # Print summary
                self.designer.print_summary(results)
                
                # Validate designs
                if results['primers'] and len(results['primers']) >= 2 and results['probe']:
                    fwd = results['primers'][0]['sequence']
                    rev = results['primers'][1]['sequence']
                    probe = results['probe']['sequence']
                    
                    print(f"  ✓ Validating RPA compatibility...")
                    validation = validate_rpa_compatibility(fwd, rev, probe)
                    print_validation_report(validation, target_name)
                    
                    results['validation'] = validation
                
                # Save results
                self._save_results(target_name, results, info)
                self.results_summary.append({
                    'target': target_name,
                    'status': 'SUCCESS',
                    'primer_count': len(results.get('primers', [])),
                    'probe_designed': results.get('probe') is not None,
                })
                
            except Exception as e:
                print(f"  ❌ Error processing {target_name}: {str(e)}")
                self.results_summary.append({
                    'target': target_name,
                    'status': 'FAILED',
                    'error': str(e),
                })
        
        # Print final summary
        self._print_pipeline_summary()
    
    def _save_results(self, target_name: str, results: dict, target_info: dict):
        """Save design results to files."""
        # Create target-specific directory
        target_dir = self.output_dir / target_name
        target_dir.mkdir(exist_ok=True)
        
        # Save CSV
        csv_file = target_dir / f"{target_name}_primers.csv"
        self.designer.save_results(results, str(csv_file))
        
        # Save detailed report
        report_file = target_dir / f"{target_name}_report.txt"
        self._write_detailed_report(report_file, target_name, results, target_info)
        
        # Save FASTA format
        fasta_file = target_dir / f"{target_name}_primers.fasta"
        self._write_fasta_format(fasta_file, target_name, results)
        
        print(f"  ✓ Results saved to {target_dir}/")
    
    def _write_detailed_report(self, filepath: Path, target_name: str, 
                               results: dict, target_info: dict):
        """Write a comprehensive design report."""
        with open(filepath, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"RecJf-RPA lncRNA Primer Design Report\n")
            f.write(f"Target: {target_name}\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")
            
            # Target information
            f.write("TARGET INFORMATION\n")
            f.write("-" * 80 + "\n")
            info = target_info
            f.write(f"Name: {info['name']}\n")
            f.write(f"Chromosome: {info['chromosome']}\n")
            f.write(f"Length: {info['length_bp']} bp\n")
            f.write(f"Function: {info['function']}\n")
            f.write(f"Cancer Association: {info['cancer_association']}\n")
            f.write(f"Expression Pattern: {info['expression_pattern']}\n\n")
            
            # Design parameters
            f.write("DESIGN PARAMETERS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Sequence Length: {results['sequence_length']} bp\n")
            f.write(f"Primer Length Range: 17-25 bp (RPA optimized)\n")
            f.write(f"Probe Length: 40-60 bp (RecJf protection zone)\n")
            f.write(f"GC Content Target: 40-60%\n")
            f.write(f"Max Tm Difference: 5.0°C\n\n")
            
            # Designed primers
            f.write("DESIGNED PRIMERS\n")
            f.write("-" * 80 + "\n")
            for primer in results.get('primers', []):
                f.write(f"\n{primer['type']} Primer:\n")
                f.write(f"  Sequence: 5'-{primer['sequence']}-3'\n")
                f.write(f"  Length: {primer['length']} bp\n")
                f.write(f"  GC Content: {primer['gc_content']*100:.1f}%\n")
                f.write(f"  Tm: {primer['tm_celsius']:.1f}°C\n")
                f.write(f"  Design Score: {primer['score']:.1f}/100\n")
                f.write(f"  Secondary Structure: {primer['secondary_structure']}\n")
                f.write(f"  Homopolymers: {primer['homopolymers']}\n")
            
            # Designed probe
            if results.get('probe'):
                f.write(f"\n\nProbe (RecJf Protected Region):\n")
                probe = results['probe']
                f.write(f"  Sequence: 5'-{probe['sequence']}-3'\n")
                f.write(f"  Length: {probe['length']} bp\n")
                f.write(f"  GC Content: {probe['gc_content']*100:.1f}%\n")
                f.write(f"  Tm: {probe['tm_celsius']:.1f}°C\n")
                f.write(f"  Note: {probe['note']}\n")
            
            # Validation
            if 'validation' in results:
                f.write(f"\n\nVALIDATION RESULTS\n")
                f.write("-" * 80 + "\n")
                validation = results['validation']
                f.write(f"Overall Compatible: {validation['overall_compatible']}\n")
                if validation['metrics']:
                    f.write(f"\nMetrics:\n")
                    for key, value in validation['metrics'].items():
                        if isinstance(value, float):
                            f.write(f"  {key}: {value:.1f}%\n")
                        else:
                            f.write(f"  {key}: {value}\n")
                if validation['warnings']:
                    f.write(f"\nWarnings:\n")
                    for warning in validation['warnings']:
                        f.write(f"  - {warning}\n")
                if validation['recommendations']:
                    f.write(f"\nRecommendations:\n")
                    for rec in validation['recommendations']:
                        f.write(f"  - {rec}\n")
            
            # Experimental notes
            f.write("\n\nEXPERIMENTAL NOTES FOR OPTIMIZATION\n")
            f.write("-" * 80 + "\n")
            f.write("""
1. REAGENT PREPARATION:
   - RecJf exonuclease: 1 unit/μL (NEB M0264)
   - RPA enzyme mix: Follow TwistAmp protocol
   - Mg2+ optimum: 3-5 mM for RPA reaction

2. REACTION CONDITIONS:
   - Temperature: 37°C (RPA optimal for RecJf-assisted detection)
   - Time: 15-30 minutes for amplification
   - Total volume: 50 μL standard reaction

3. PRIMER OPTIMIZATION:
   - Final concentration: 480 nM (240 nM each primer)
   - Temperature cycling not required
   - Isothermal amplification reduces background

4. PROBE OPTIMIZATION:
   - DNA probe concentration: 100-500 nM
   - RecJf digestion time: 5-10 minutes after RPA
   - Fluorescence readout: Real-time or endpoint

5. DETECTION:
   - Label-free detection using RecJf digestion
   - No RT step needed for lncRNA
   - Potential for integration with CRISPR-Cas systems

6. TROUBLESHOOTING:
   - If low signal: Increase probe concentration or RecJf time
   - If high background: Optimize Mg2+ concentration
   - If primer dimers: Reduce primer concentration or redesign
""")
            
            f.write("\n" + "=" * 80 + "\n")
    
    def _write_fasta_format(self, filepath: Path, target_name: str, results: dict):
        """Write primers and probe in FASTA format."""
        with open(filepath, 'w') as f:
            for i, primer in enumerate(results.get('primers', []), 1):
                ptype = primer['type']
                seq = primer['sequence']
                f.write(f">{target_name}__{ptype}_Primer\n")
                f.write(f"{seq}\n")
            
            if results.get('probe'):
                probe = results['probe']
                f.write(f">{target_name}__Probe\n")
                f.write(f"{probe['sequence']}\n")
    
    def _print_pipeline_summary(self):
        """Print final pipeline summary."""
        print("\n" + "=" * 70)
        print("  PIPELINE SUMMARY")
        print("=" * 70)
        
        successful = sum(1 for r in self.results_summary if r['status'] == 'SUCCESS')
        failed = sum(1 for r in self.results_summary if r['status'] == 'FAILED')
        
        print(f"\nTotal targets processed: {len(self.results_summary)}")
        print(f"✓ Successful: {successful}")
        print(f"✗ Failed: {failed}")
        
        print(f"\nResults saved to: {self.output_dir.absolute()}")
        print(f"\nGenerated files per target:")
        print(f"  - {{target}}_primers.csv          (primer sequences and properties)")
        print(f"  - {{target}}_primers.fasta        (FASTA format for synthesis)")
        print(f"  - {{target}}_report.txt           (detailed experimental report)")
        
        print("=" * 70 + "\n")


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="RecJf-RPA lncRNA Primer Design Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py                          # Run all default targets
  python run_pipeline.py --target HOTAIR          # Run specific target
  python run_pipeline.py --target HOTAIR MALAT1   # Run multiple targets
  python run_pipeline.py --output results/        # Specify output directory
        """
    )
    
    parser.add_argument(
        '--target',
        nargs='+',
        help='Specific lncRNA target(s) to design primers for'
    )
    parser.add_argument(
        '--output',
        default='output',
        help='Output directory for results (default: output)'
    )
    parser.add_argument(
        '--list-targets',
        action='store_true',
        help='List all available targets and exit'
    )
    
    args = parser.parse_args()
    
    # List available targets if requested
    if args.list_targets:
        targets = get_lncrna_targets()
        target_info = get_target_info()
        print("\nAvailable lncRNA Targets:")
        print("=" * 70)
        for name in sorted(targets.keys()):
            info = target_info[name]
            print(f"\n{name}")
            print(f"  Name: {info['name']}")
            print(f"  Chromosome: {info['chromosome']}")
            print(f"  Function: {info['function']}")
            print(f"  Cancer Assoc: {info['cancer_association']}")
        print("\n")
        return
    
    # Run pipeline
    pipeline = RecJfRPAPipeline(output_dir=args.output)
    pipeline.run_design_pipeline(target_names=args.target)


if __name__ == "__main__":
    main()
