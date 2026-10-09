"""
RecJf-RPA lncRNA Primer Design Module

Core algorithms for designing primers and probes optimized for:
- RecJf exonuclease-assisted RPA detection
- lncRNA targets in breast cancer tissue
- Label-free, rapid detection

Author: RecJf-RPA Development Team
License: MIT
"""

import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
import csv
from datetime import datetime

@dataclass
class PrimerDesignResult:
    """Data class for storing primer design results."""
    target_name: str
    primer_type: str  # 'Forward', 'Reverse', 'Probe'
    sequence: str
    length: int
    gc_content: float
    tm_celsius: float
    secondary_structure_risk: str  # 'Low', 'Medium', 'High'
    homopolymer_runs: str  # 'None', 'Present'
    design_score: float  # 0-100, higher is better
    notes: str
    timestamp: str


class RPAPrimerDesigner:
    """
    Main class for designing RPA primers and probes for lncRNA detection.
    
    Features:
    - Tm calculation (nearest-neighbor method approximation)
    - GC content analysis
    - Secondary structure prediction
    - Homopolymer detection
    - Design scoring
    """
    
    # Nucleotide base pairing constants
    DNA_BASES = {'A', 'T', 'G', 'C'}
    
    # Tm calculation parameters (simplified nearest-neighbor model)
    # Enthalpy (ΔH) in kcal/mol for dinucleotide pairs
    ENTHALPY_TABLE = {
        'AA': -8.0, 'TT': -8.0, 'AT': -7.6, 'TA': -7.2,
        'CA': -8.5, 'TG': -8.5, 'CT': -7.8, 'AG': -7.8,
        'GA': -8.2, 'TC': -8.2, 'GT': -8.4, 'AC': -8.4,
        'CC': -8.0, 'GG': -8.0, 'CG': -10.6, 'GC': -9.8,
    }
    
    # Entropy (ΔS) in cal/(mol·K)
    ENTROPY_TABLE = {
        'AA': -22.2, 'TT': -22.2, 'AT': -20.4, 'TA': -21.3,
        'CA': -22.7, 'TG': -22.7, 'CT': -21.0, 'AG': -21.0,
        'GA': -22.2, 'TC': -22.2, 'GT': -22.4, 'AC': -22.4,
        'CC': -19.9, 'GG': -19.9, 'CG': -27.8, 'GC': -24.4,
    }
    
    def __init__(self, rna_salt_conc: float = 0.05, mg_conc: float = 0.003):
        """
        Initialize the primer designer.
        
        Args:
            rna_salt_conc: Salt concentration (M), default 50 mM for RPA
            mg_conc: Mg2+ concentration (M), default 3 mM for RPA
        """
        self.salt_conc = rna_salt_conc
        self.mg_conc = mg_conc
    
    def design_primers(
        self,
        sequence: str,
        target_name: str,
        primer_length_range: Tuple[int, int] = (17, 25),
        probe_length: int = 50,
        gc_range: Tuple[float, float] = (0.40, 0.60),
        max_tm_diff: float = 5.0
    ) -> Dict:
        """
        Design RPA primers and probe for a given lncRNA sequence.
        
        Args:
            sequence: lncRNA sequence (FASTA or plain string)
            target_name: Name of the lncRNA target
            primer_length_range: Tuple (min, max) for primer length
            probe_length: Target probe length
            gc_range: Tuple (min, max) for GC content as fraction
            max_tm_diff: Maximum Tm difference between F and R primers
        
        Returns:
            Dictionary containing designed primers and probe with full annotations
        """
        # Clean sequence
        sequence = self._clean_sequence(sequence)
        
        if len(sequence) < 100:
            raise ValueError(f"Sequence too short: {len(sequence)} bp. Need at least 100 bp.")
        
        results = {
            'target_name': target_name,
            'sequence_length': len(sequence),
            'timestamp': datetime.now().isoformat(),
            'primers': [],
            'probe': None,
            'design_summary': {}
        }
        
        # Strategy: Select regions from sequence start, middle, and end
        regions = self._select_design_regions(sequence, num_regions=3)
        
        best_forward = None
        best_reverse = None
        best_score = -999
        
        for region in regions:
            # Forward primer: select from start of region
            forward_candidates = self._extract_candidates(
                region, 
                primer_length_range, 
                gc_range
            )
            
            for fwd_seq in forward_candidates[:5]:  # Top 5 candidates
                fwd_tm = self.calculate_tm(fwd_seq)
                fwd_gc = self.calculate_gc_content(fwd_seq)
                fwd_score = self._score_primer(fwd_seq, fwd_tm, fwd_gc)
                
                # Reverse primer: search within reasonable distance
                reverse_start = min(50, len(region) - primer_length_range[1])
                reverse_region = region[reverse_start:reverse_start + 200]
                reverse_candidates = self._extract_candidates(
                    reverse_region,
                    primer_length_range,
                    gc_range
                )
                
                for rev_seq in reverse_candidates[:5]:
                    rev_tm = self.calculate_tm(rev_seq)
                    rev_gc = self.calculate_gc_content(rev_seq)
                    rev_score = self._score_primer(rev_seq, rev_tm, rev_gc)
                    
                    # Check compatibility
                    tm_diff = abs(fwd_tm - rev_tm)
                    if tm_diff <= max_tm_diff:
                        combined_score = (fwd_score + rev_score) / 2
                        if combined_score > best_score:
                            best_score = combined_score
                            best_forward = (fwd_seq, fwd_tm, fwd_gc, fwd_score)
                            best_reverse = (rev_seq, rev_tm, rev_gc, rev_score)
        
        # Add forward primer
        if best_forward:
            fwd_seq, fwd_tm, fwd_gc, fwd_score = best_forward
            results['primers'].append({
                'type': 'Forward',
                'sequence': fwd_seq,
                'length': len(fwd_seq),
                'gc_content': fwd_gc,
                'tm_celsius': fwd_tm,
                'score': fwd_score,
                'secondary_structure': self._assess_secondary_structure(fwd_seq),
                'homopolymers': self._find_homopolymers(fwd_seq),
            })
        
        # Add reverse primer
        if best_reverse:
            rev_seq, rev_tm, rev_gc, rev_score = best_reverse
            results['primers'].append({
                'type': 'Reverse',
                'sequence': rev_seq,
                'length': len(rev_seq),
                'gc_content': rev_gc,
                'tm_celsius': rev_tm,
                'score': rev_score,
                'secondary_structure': self._assess_secondary_structure(rev_seq),
                'homopolymers': self._find_homopolymers(rev_seq),
            })
        
        # Design probe
        if len(sequence) >= probe_length + 50:
            probe_start = max(0, len(sequence) // 2 - probe_length // 2)
            probe_seq = sequence[probe_start:probe_start + probe_length]
            
            results['probe'] = {
                'sequence': probe_seq,
                'length': len(probe_seq),
                'gc_content': self.calculate_gc_content(probe_seq),
                'tm_celsius': self.calculate_tm(probe_seq),
                'note': 'Probe targets lncRNA, DNA complement will be protected from RecJf'
            }
        
        results['design_summary'] = {
            'num_forward_primers': 1 if best_forward else 0,
            'num_reverse_primers': 1 if best_reverse else 0,
            'probe_designed': results['probe'] is not None,
            'amplicon_length_estimate': 'F and R primer positions determine amplicon',
        }
        
        return results
    
    def calculate_tm(self, sequence: str) -> float:
        """
        Calculate melting temperature (Tm) using nearest-neighbor method.
        Simplified implementation suitable for RPA primers (17-25 bp).
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            Tm in Celsius
        """
        seq = sequence.upper()
        
        # Simple Tm for short primers (17-25 bp): Tm = 4(G+C) + 2(A+T)
        if len(seq) <= 25:
            gc_count = seq.count('G') + seq.count('C')
            at_count = seq.count('A') + seq.count('T')
            tm = 4 * gc_count + 2 * at_count
            return min(tm, 72)  # Cap at 72°C for RPA
        
        # For longer sequences (probes), use nearest-neighbor calculation
        if len(seq) > 25:
            # Simplified calculation
            gc_content = (seq.count('G') + seq.count('C')) / len(seq)
            tm = 64.9 + 41 * (gc_content - 0.5) / 0.5
            tm += 16.6 * math.log10(self.salt_conc)
            return min(max(tm, 45), 75)  # Bound between 45-75°C
        
        return 60.0  # Default fallback
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content as fraction (0-1).
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            GC content as fraction
        """
        seq = sequence.upper()
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq) if len(seq) > 0 else 0.0
    
    def _clean_sequence(self, sequence: str) -> str:
        """Remove whitespace and non-DNA characters."""
        seq = re.sub(r'[\s\n\r]', '', sequence.upper())
        seq = re.sub(r'[^ATGC]', '', seq)  # Keep only ATGC
        return seq
    
    def _select_design_regions(self, sequence: str, num_regions: int = 3) -> List[str]:
        """Select promising regions from sequence."""
        regions = []
        region_length = min(300, len(sequence) // 3)
        
        # Start region
        regions.append(sequence[:region_length])
        
        # Middle region
        if len(sequence) > region_length * 2:
            mid = len(sequence) // 2
            regions.append(sequence[mid - region_length // 2:mid + region_length // 2])
        
        # End region
        if len(sequence) > region_length:
            regions.append(sequence[-region_length:])
        
        return regions[:num_regions]
    
    def _extract_candidates(
        self,
        sequence: str,
        length_range: Tuple[int, int],
        gc_range: Tuple[float, float]
    ) -> List[str]:
        """Extract candidate primers from region based on length and GC content."""
        candidates = []
        min_len, max_len = length_range
        min_gc, max_gc = gc_range
        
        for pos in range(0, len(sequence) - min_len + 1, 2):  # Step by 2
            for length in range(min_len, min(max_len + 1, len(sequence) - pos + 1)):
                candidate = sequence[pos:pos + length]
                gc = self.calculate_gc_content(candidate)
                
                if min_gc <= gc <= max_gc and candidate.count('A') <= length * 0.5:
                    # Avoid poly-A
                    candidates.append(candidate)
        
        return candidates
    
    def _score_primer(self, sequence: str, tm: float, gc_content: float) -> float:
        """
        Score primer quality (0-100).
        Higher score = better primer.
        """
        score = 50.0
        
        # Tm score: optimal range 60-70°C
        if 60 <= tm <= 70:
            score += 20
        elif 55 <= tm < 60:
            score += 10
        elif 70 < tm <= 75:
            score += 10
        
        # GC content score: optimal 45-55%
        gc_pct = gc_content * 100
        if 45 <= gc_pct <= 55:
            score += 15
        elif 40 <= gc_pct < 45:
            score += 8
        elif 55 < gc_pct <= 60:
            score += 8
        
        # Secondary structure penalty
        hairpin_risk = self._assess_secondary_structure(sequence)
        if hairpin_risk == "Low":
            score += 10
        elif hairpin_risk == "Medium":
            score += 5
        
        # Homopolymer penalty
        if self._find_homopolymers(sequence) == "None":
            score += 5
        
        return min(score, 100)
    
    def _assess_secondary_structure(self, sequence: str) -> str:
        """
        Assess risk of hairpin/secondary structure.
        Returns: 'Low', 'Medium', 'High'
        """
        seq = sequence.upper()
        
        # Simple check: look for self-complementarity
        complement = str.maketrans('ATCG', 'TAGC')
        rev_comp = seq.translate(complement)[::-1]
        
        # Count matching positions in first 8 bp with reverse complement
        matches = sum(1 for i in range(min(8, len(seq))) if seq[i] == rev_comp[-(i+1)])
        
        if matches >= 6:
            return "High"
        elif matches >= 4:
            return "Medium"
        else:
            return "Low"
    
    def _find_homopolymers(self, sequence: str) -> str:
        """Detect long homopolymer runs (>4 consecutive same base)."""
        for base in 'ATGC':
            if f"{base * 5}" in sequence:
                return f"Present ({base}x5+)"
        return "None"
    
    def save_results(self, results: Dict, output_file: str = "primer_results.csv"):
        """
        Save design results to CSV file.
        
        Args:
            results: Dictionary from design_primers()
            output_file: Output CSV filename
        """
        rows = []
        
        target_name = results['target_name']
        
        # Add forward and reverse primers
        for primer_info in results['primers']:
            row = {
                'Target': target_name,
                'Type': primer_info['type'],
                'Sequence': primer_info['sequence'],
                'Length (bp)': primer_info['length'],
                'GC Content (%)': f"{primer_info['gc_content'] * 100:.1f}",
                'Tm (°C)': f"{primer_info['tm_celsius']:.1f}",
                'Score': f"{primer_info['score']:.1f}",
                'Secondary Structure Risk': primer_info['secondary_structure'],
                'Homopolymers': primer_info['homopolymers'],
            }
            rows.append(row)
        
        # Add probe
        if results['probe']:
            probe_info = results['probe']
            row = {
                'Target': target_name,
                'Type': 'Probe',
                'Sequence': probe_info['sequence'],
                'Length (bp)': probe_info['length'],
                'GC Content (%)': f"{probe_info['gc_content'] * 100:.1f}",
                'Tm (°C)': f"{probe_info['tm_celsius']:.1f}",
                'Score': 'N/A',
                'Secondary Structure Risk': 'N/A',
                'Homopolymers': 'N/A',
            }
            rows.append(row)
        
        # Write CSV
        if rows:
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            print(f"✓ Results saved to {output_file}")
    
    def print_summary(self, results: Dict):
        """Print a summary of design results to console."""
        print("\n" + "=" * 70)
        print(f"  RPA Primer Design Summary: {results['target_name']}")
        print("=" * 70)
        print(f"Sequence length: {results['sequence_length']} bp")
        print(f"Timestamp: {results['timestamp']}")
        
        print(f"\nForward Primer:")
        if results['primers'] and results['primers'][0]['type'] == 'Forward':
            p = results['primers'][0]
            print(f"  Sequence: 5'-{p['sequence']}-3'")
            print(f"  Length: {p['length']} bp")
            print(f"  GC: {p['gc_content']*100:.1f}%")
            print(f"  Tm: {p['tm_celsius']:.1f}°C")
            print(f"  Score: {p['score']:.1f}/100")
        
        print(f"\nReverse Primer:")
        if len(results['primers']) > 1 and results['primers'][1]['type'] == 'Reverse':
            p = results['primers'][1]
            print(f"  Sequence: 5'-{p['sequence']}-3'")
            print(f"  Length: {p['length']} bp")
            print(f"  GC: {p['gc_content']*100:.1f}%")
            print(f"  Tm: {p['tm_celsius']:.1f}°C")
            print(f"  Score: {p['score']:.1f}/100")
        
        print(f"\nProbe (RecJf Protected Region):")
        if results['probe']:
            pb = results['probe']
            print(f"  Sequence: 5'-{pb['sequence']}-3'")
            print(f"  Length: {pb['length']} bp")
            print(f"  GC: {pb['gc_content']*100:.1f}%")
            print(f"  Tm: {pb['tm_celsius']:.1f}°C")
            print(f"  Note: {pb['note']}")
        
        print("\n" + "=" * 70 + "\n")


# Import math for Tm calculation
import math


if __name__ == "__main__":
    # Test with sample sequence
    test_sequence = """
    AACATAATAAAATGCTGTATAGCATAGAACTGAGACTGATATAATGCTGTATAGCATAGAACTGAGACTGATATAATG
    CTGTATAGCATAGAACTGAGACTGATATAGTGCTGTATAGCATAGAACTGAGACTGATATAATGCTGTATAGCATAGAA
    CTGAGACTGATATAATGCTGTATAGCATAGAACTGAGACTGATATAATGCTGCATACCATAGAACTGAGACAGATATAAG
    CCTGGGTGAAAGATTCGTATGTAGCTCTGTTGGGTTGCTGGAGGAGAGTGAGGAGACACTTGCTGTCCCCTCAGCAGCC
    """
    
    designer = RPAPrimerDesigner()
    results = designer.design_primers(
        sequence=test_sequence,
        target_name="TEST_lncRNA",
        primer_length_range=(17, 25),
        probe_length=50
    )
    
    designer.print_summary(results)
    designer.save_results(results, output_file="test_primers.csv")
