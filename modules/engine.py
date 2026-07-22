""" App engine """

import os

import sqlite3
import gzip

from Bio import SeqIO, Entrez
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from BCBio import GFF

import modules.lib as lib

class Core:
    """ Calculus core"""
    def __init__(self, gui_cmd_buff, data_manager):
        self.data_manager = data_manager
        self.gui_command_buffer = gui_cmd_buff

    def run_comparing(self, patient_data_path, gene_id):
        """ Main comparing function """
        self.gui_command_buffer.put(("PROGRESS", [0, "Opening patient file..."]))

        # 1. Open file
        in_seq_handle = gzip.open(patient_data_path, "rt")
        patient_seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        if len(patient_seq_dict) != 1:
            lib.log(f"Expected one sequence in patient file, got {len(patient_seq_dict)}")
            return []

        # 2. Check & Download reference
        Entrez.email = lib.USER_EMAIL

        gene_chromosome_strand = self._resolve_gene_chromosome_strand(gene_id)

        search_handle = Entrez.esearch(
            db="nucleotide",
            term=f"{gene_id}[Gene Name] AND RefSeq[Keyword] AND Homo sapiens[Organism]"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            lib.log(f"Gene identifier {gene_id} not resolved on NCBI servers.")
            return []

        fetch_id = search_results["IdList"][0]

        ref_file_name = f"{gene_id}_reference.fasta.gz"
        ref_data_path = os.path.join(lib.GENES_CACHE_DIR, ref_file_name)

        self.gui_command_buffer.put(("PROGRESS", [1, "Loading a reference..."]))

        if not os.path.exists(ref_data_path):
            lib.log(f"Reference for {gene_id} not found locally. Fetching from NCBI...")

            try:
                fetch_handle = Entrez.efetch(db="nucleotide",
                                             id=fetch_id, rettype="fasta", retmode="text")
                fasta_data = fetch_handle.read()
                fetch_handle.close()

                with gzip.open(ref_data_path, "wt") as out_file:
                    out_file.write(fasta_data)

            except Exception as error:
                lib.log(f"Failed to download NCBI reference: {str(error)}")
                return []

        else: lib.log("Reference data file found. Skipping download.")

        # 3. Check & Download GFF3
        gff_file_name = f"{gene_id}_reference.gff3"
        ref_anno_path = os.path.join(lib.GENES_CACHE_DIR, gff_file_name)

        self.gui_command_buffer.put(("PROGRESS", [2, "Loading an annotation..."]))

        if not os.path.exists(ref_anno_path):
            lib.log(f"Annotation for {gene_id} not found locally. Fetching from NCBI...")
            try:
                fetch_handle = Entrez.efetch(db="nucleotide",
                                             id=fetch_id, rettype="gff3", retmode="text")
                gff_data = fetch_handle.read()
                fetch_handle.close()

                with open(ref_anno_path, "w", encoding="utf-8") as out_file:
                    out_file.write(gff_data)

            except Exception as error:
                lib.log(f"Failed to download NCBI annotation: {str(error)}")
                return []

        else: lib.log("Annotation file found. Skipping download.")

        in_seq_handle = gzip.open(ref_data_path, "rt")
        ref_seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        lib.dbg(f"Patient seq dict: {patient_seq_dict}, Ref seq dict: {ref_seq_dict}")

        patient_raw_record = list(patient_seq_dict.values())[0]
        ref_raw_record = list(ref_seq_dict.values())[0]

        # 4. Align sequences
        self.gui_command_buffer.put(("PROGRESS", [3, "Aligning..."]))

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")

        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        dna_alignment = aligner.align(ref_raw_record.seq, patient_raw_record.seq)[0]

        ref_blocks = dna_alignment.aligned[0]
        pat_blocks = dna_alignment.aligned[1]

        lib.dbg(f"DNA alignment score: {dna_alignment.score}")
        lib.dbg(f"Aligned blocks: {len(ref_blocks)}")

        # 5. Apply annotation
        self.gui_command_buffer.put(("PROGRESS", [4, "Applying annotation..."]))

        with open(ref_anno_path, "r", encoding="utf-8") as in_handle:
            gff_records = list(GFF.parse(in_handle))

        def is_mane_select(transcript_feature):
            tags = transcript_feature.qualifiers.get("tag", [])
            return any("MANE Select" in t or "MANE_Select" in t for t in tags)

        candidates = self._collect_transcripts(gff_records, gene_id)
        if not candidates:
            lib.log(f"No CDS-bearing transcript found for {gene_id}")
            return []

        mane_candidates = [c for c in candidates if is_mane_select(c[0])]
        if mane_candidates:
            chosen_transcript, ref_cds_feats = mane_candidates[0]
            lib.log(f"Using MANE Select transcript {chosen_transcript.id}")
        else:
            chosen_transcript, ref_cds_feats = max(
                candidates, key=lambda c: sum(len(f.location) for f in c[1])
            )
            lib.log(f"No MANE Select tag, using longest-CDS transcript {chosen_transcript.id}")

        ref_cds_locs = [f.location for f in ref_cds_feats]

        patient_cds_locs = [self._map_location_to_patient(loc, ref_blocks, pat_blocks)
                            for loc in ref_cds_locs]

        patient_cds_parts = [str(loc.extract(patient_raw_record.seq)) for loc in patient_cds_locs]
        ref_cds_parts = [str(loc.extract(ref_raw_record.seq)) for loc in ref_cds_locs]

        patient_dna_seq = Seq("".join(patient_cds_parts))
        ref_dna_seq = Seq("".join(ref_cds_parts))

        ref_dna_origin = self._build_nucleotide_origin_map(ref_cds_locs)
        patient_dna_origin = self._build_nucleotide_origin_map(patient_cds_locs)

        lib.dbg(f"ref_dna_origin length={len(ref_dna_origin)}, "
                f"ref_dna_seq length={len(ref_dna_seq)}")
        lib.dbg(f"patient_dna_origin length={len(patient_dna_origin)}, "
                f"patient_dna_seq length={len(patient_dna_seq)}")

        gene_strand = gene_chromosome_strand
        lib.dbg(f"Using TRUE chromosome strand for nucleotide reorientation: {gene_strand}")

        # 6. Turn sequence into aminoacids
        self.gui_command_buffer.put(("PROGRESS", [5, "Translating into aminoacids..."]))

        patient_amino_full = str(patient_dna_seq.translate())
        ref_amino_full = str(ref_dna_seq.translate())

        def truncate_at_stop(amino_str):
            stop_pos = amino_str.find("*")
            if stop_pos == -1:
                return amino_str, None
            return amino_str[:stop_pos], stop_pos

        patient_amino, patient_stop_pos = truncate_at_stop(patient_amino_full)
        ref_amino, ref_stop_pos = truncate_at_stop(ref_amino_full)

        patient_amino = Seq(patient_amino)
        ref_amino = Seq(ref_amino)

        lib.dbg(
            f"Patient CDS length: {len(patient_dna_seq)} bp, "
            f"frame remainder: {len(patient_dna_seq) % 3}"
        )
        lib.dbg(f"Ref stop at aa {ref_stop_pos}, patient stop at aa {patient_stop_pos}")

        protein_aligner = PairwiseAligner()
        protein_aligner.mode = "global"
        protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

        protein_aligner.open_gap_score = -10
        protein_aligner.extend_gap_score = -0.5

        protein_alignment = protein_aligner.align(ref_amino, patient_amino)[0]

        # 7. Extracting all mutations
        self.gui_command_buffer.put(("PROGRESS", [6, "Searching for mutations..."]))

        mutations = []

        # > Nonsense / Stop-loss
        if patient_stop_pos is not None and ref_stop_pos is not None:
            if patient_stop_pos < ref_stop_pos:
                wt_aa = ref_amino_full[patient_stop_pos] if patient_stop_pos < len(ref_amino_full) else "?"

                ref_codon, ref_origin = self._get_codon_info(
                    ref_dna_seq, ref_dna_origin, patient_stop_pos)
                pat_codon, pat_origin = self._get_codon_info(
                    patient_dna_seq, patient_dna_origin, patient_stop_pos)

                ref_codon, pat_codon, ref_origin, pat_origin = self._trim_common_flanks(
                    ref_codon, pat_codon, ref_origin, pat_origin
                )

                ref_codon, ref_local_pos = self._to_genomic_representation(
                    ref_codon, ref_origin, gene_strand)
                pat_codon, pat_local_pos = self._to_genomic_representation(
                    pat_codon, pat_origin, gene_strand)

                mutations.append({
                    "pos": patient_stop_pos + 1,
                    "type": "Nonsense",
                    "hgvs": f"p.{wt_aa}{patient_stop_pos + 1}*",
                    "ref": wt_aa,
                    "alt": "*",
                    "nt_ref": ref_codon,
                    "nt_alt": pat_codon,
                    "ref_local_pos": ref_local_pos,
                    "patient_local_pos": pat_local_pos,
                })
            elif patient_stop_pos > ref_stop_pos:
                new_aa = patient_amino_full[ref_stop_pos] if ref_stop_pos < len(patient_amino_full) else "?"

                ref_codon, ref_origin = self._get_codon_info(
                    ref_dna_seq, ref_dna_origin, ref_stop_pos)
                pat_codon, pat_origin = self._get_codon_info(
                    patient_dna_seq, patient_dna_origin, ref_stop_pos)

                ref_codon, pat_codon, ref_origin, pat_origin = self._trim_common_flanks(
                    ref_codon, pat_codon, ref_origin, pat_origin)

                ref_codon, ref_local_pos = self._to_genomic_representation(
                    ref_codon, ref_origin, gene_strand)
                pat_codon, pat_local_pos = self._to_genomic_representation(
                    pat_codon, pat_origin, gene_strand)

                mutations.append({
                    "pos": ref_stop_pos + 1,
                    "type": "Stop-loss",
                    "hgvs": f"p.*{ref_stop_pos + 1}{new_aa}ext*?",
                    "ref": "*",
                    "alt": new_aa,
                    "nt_ref": ref_codon,
                    "nt_alt": pat_codon,
                    "ref_local_pos": ref_local_pos,
                    "patient_local_pos": pat_local_pos,
                })
        elif patient_stop_pos is None:
            lib.log("No stop codon found in patient CDS translation.")

        # > Insertion / Deletion / Substuition
        for ref_range, pat_range in zip(protein_alignment.aligned[0], protein_alignment.aligned[1]):
            ref_start, ref_end = ref_range
            pat_start, pat_end = pat_range

            ref_slice = ref_amino[ref_start:ref_end]
            pat_slice = patient_amino[pat_start:pat_end]

            if ref_slice == pat_slice:
                continue

            # Substitution
            if len(ref_slice) == len(pat_slice):
                for index, (ref_char, pat_char) in enumerate(zip(ref_slice, pat_slice)):
                    if ref_char != pat_char:
                        pos = ref_start + index + 1

                        ref_codon, ref_origin = self._get_codon_info(ref_dna_seq, ref_dna_origin,
                                                                     ref_start + index)
                        pat_codon, pat_origin = self._get_codon_info(patient_dna_seq, patient_dna_origin,
                                                                     pat_start + index)

                        ref_codon, pat_codon, ref_origin, pat_origin = self._trim_common_flanks(
                            ref_codon, pat_codon, ref_origin, pat_origin
                        )

                        ref_codon, ref_local_pos = self._to_genomic_representation(
                            ref_codon, ref_origin, gene_strand)
                        pat_codon, pat_local_pos = self._to_genomic_representation(
                            pat_codon, pat_origin, gene_strand)

                        mutations.append({
                            "pos": int(pos),
                            "type": "Substitution",
                            "hgvs": f"p.{ref_char}{pos}{pat_char}",
                            "ref": str(ref_char),
                            "alt": str(pat_char),
                            "nt_ref": ref_codon,
                            "nt_alt": pat_codon,
                            "ref_local_pos": ref_local_pos,
                            "patient_local_pos": pat_local_pos,
                        })

            # Deletion
            elif len(ref_slice) > len(pat_slice):
                start_pos = ref_start + 1
                end_pos = ref_end

                if start_pos == end_pos:
                    hgvs = f"p.{ref_amino[ref_start]}{start_pos}del"
                else:
                    hgvs = (
                        f"p.{ref_amino[ref_start]}{start_pos}_"
                        f"{ref_amino[ref_end - 1]}{end_pos}del"
                    )

                nt_start = ref_start * 3
                nt_end = ref_end * 3

                if nt_end <= len(ref_dna_seq) and nt_end <= len(ref_dna_origin):
                    nt_ref_block = str(ref_dna_seq[nt_start:nt_end])
                    ref_origin_block = ref_dna_origin[nt_start:nt_end]
                else:
                    nt_ref_block, ref_origin_block = None, None

                nt_ref_block, ref_local_pos = self._to_genomic_representation(
                    nt_ref_block, ref_origin_block, gene_strand)

                mutations.append({
                    "pos": int(start_pos),
                    "type": "Deletion",
                    "hgvs": hgvs,
                    "ref": None,
                    "alt": None,
                    "nt_ref": nt_ref_block,
                    "nt_alt": "",
                    "ref_local_pos": ref_local_pos,
                    "patient_local_pos": None,
                })

            # Insertion
            elif len(ref_slice) < len(pat_slice):
                left_pos = ref_start
                right_pos = ref_start + 1

                left_aa = ref_amino[left_pos - 1] if left_pos - 1 >= 0 else None
                right_aa = ref_amino[right_pos - 1] if right_pos - 1 < len(ref_amino) else None

                if left_aa is not None and right_aa is not None:
                    hgvs = f"p.{left_aa}{left_pos}_{right_aa}{right_pos}ins{str(pat_slice)}"

                elif left_aa is not None and right_aa is None:
                    hgvs = f"p.{left_aa}{left_pos}_Terins{str(pat_slice)}"

                else: hgvs = f"p.0_{right_aa}{right_pos}ins{str(pat_slice)}"

                nt_start = pat_start * 3
                nt_end = pat_end * 3
                if nt_end <= len(patient_dna_seq) and nt_end <= len(patient_dna_origin):
                    nt_alt_block = str(patient_dna_seq[nt_start:nt_end])
                    patient_origin_block = patient_dna_origin[nt_start:nt_end]
                else:
                    nt_alt_block, patient_origin_block = None, None

                nt_alt_block, patient_local_pos = self._to_genomic_representation(
                    nt_alt_block, patient_origin_block, gene_strand)

                mutations.append({
                    "pos": int(left_pos),
                    "type": "Insertion",
                    "hgvs": hgvs,
                    "ref": None,
                    "alt": None,
                    "nt_ref": "",
                    "nt_alt": nt_alt_block,
                    "ref_local_pos": None,
                    "patient_local_pos": patient_local_pos,
                })

        resolved_nt_count = sum(1 for m in mutations
                                if m.get("nt_ref") is not None or m.get("nt_alt") is not None)

        lib.dbg(f"Nucleotide-level info resolved for "
                f"{resolved_nt_count}/{len(mutations)} mutation.")

        self.gui_command_buffer.put(("PROGRESS", [7, "Analysis completed."]))

        return mutations

    def find_mutations(self, anomalies, gene) -> list:
        """ Convert raw anomalies into mutation data via DB"""
        full_mutations_data = []
        lib.log(f"Parsing {len(anomalies)} mutations...")

        for ai, anomaly in enumerate(anomalies):
            try:
                local_pos, mutation_type, hgvs = anomaly["pos"], anomaly["type"], anomaly["hgvs"]
                ref_fetched, alt_fetched = anomaly["nt_ref"], anomaly["nt_alt"]
                #chromosome, position, ref_allele, alt_allele, clnvs, clinical_significance, disease_name
                mutation = self.data_manager.disease_database.find_mutation(hgvs, gene, alt_fetched)
                if mutation:
                    chrid, p, refer, alter, vs, sign, name = mutation

                    # Process values
                    pos : int = p if p is not None else local_pos
                    clnvs : str = vs if vs is not None else mutation_type

                    ref : str = refer if refer is not None else ref_fetched
                    alt : str = alter if alter is not None else alt_fetched

                    chrom : int = chrid if chrid is not None else "?"

                    full_mutations_data.append([
                        ai, f"chr{chrom}", gene, hgvs, pos, ref, alt, clnvs, sign, name])

            except IndexError as e:
                lib.log(f"Data format structure error: {e}")
                continue
            except sqlite3.Error as e:
                lib.log(f"Database error: {e}")
                continue
            except Exception as e:
                lib.log(f"Unexpected mutations fetcher error: {e}")
                continue

        return full_mutations_data

    # - Inner functions - #
    def _resolve_gene_chromosome_strand(self, gene_id) -> int:
        cache_path = os.path.join(lib.GENES_CACHE_DIR, f"{gene_id}_chrom_strand.txt")

        if os.path.exists(cache_path):
            with open(cache_path, "r", encoding="utf-8") as f_in:
                cached_value = f_in.read().strip()
            if cached_value in ("1", "-1"):
                lib.dbg(f"Using cached chromosome strand for {gene_id}: {cached_value}.")
                return int(cached_value)

        try:
            search_handle = Entrez.esearch(db="gene",
                                           term=f"{gene_id}[sym] AND Homo sapiens[Organism]")
            search_results = Entrez.read(search_handle)
            search_handle.close()

            if not search_results["IdList"]:
                lib.log(f"WARNING: could not resolve chromosome strand for {gene_id}.")
                return 1

            ncbi_gene_id = search_results["IdList"][0]

            summary_handle = Entrez.esummary(db="gene", id=ncbi_gene_id)
            summary_results = Entrez.read(summary_handle)
            summary_handle.close()

            doc = summary_results["DocumentSummarySet"]["DocumentSummary"][0]
            genomic_info_list = doc.get("GenomicInfo", [])

            if not genomic_info_list:
                lib.log(f"WARNING: no GenomicInfo for {gene_id} — defaulting to plus strand.")
                return 1

            info = genomic_info_list[0]
            raw_start = int(info["ChrStart"])
            raw_stop = int(info["ChrStop"])

            strand = -1 if raw_start > raw_stop else 1

            with open(cache_path, "w", encoding="utf-8") as f_out:
                f_out.write(str(strand))

            lib.dbg(f"Resolved chromosome strand for {gene_id}: {strand} "
                     f"(ChrStart={raw_start}, ChrStop={raw_stop})")

            return strand

        except Exception as e:
            lib.log(f"WARNING: failed to resolve chromosome strand for {gene_id}: {e}")
            return 1

    def _collect_transcripts(self, gff_records, gene_id) -> list:
        candidates = []
        seen_feature_keys = set()
        for rec in gff_records:
            for feature in rec.features:
                if feature.id == gene_id or gene_id in feature.qualifiers.get("Name", []):
                    key = (rec.id, feature.location.start, feature.location.end)
                    if key in seen_feature_keys:
                        continue
                    seen_feature_keys.add(key)

                    transcripts = [sf for sf in feature.sub_features
                                   if sf.type in ("mRNA", "transcript")]
                    if transcripts:
                        for t in transcripts:
                            cds_feats = sorted(
                                [sf for sf in t.sub_features if sf.type == "CDS"],
                                key=lambda x: int(x.location.start),
                                reverse=(feature.location.strand == -1),
                            )
                            if cds_feats:
                                candidates.append((t, cds_feats))
                    else:
                        cds_feats = sorted(
                            [sf for sf in feature.sub_features if sf.type == "CDS"],
                            key=lambda x: int(x.location.start),
                            reverse=(feature.location.strand == -1),
                        )
                        if cds_feats:
                            candidates.append((feature, cds_feats))

        return candidates

    def _map_ref_pos_to_patient(self, pos, ref_blocks, pat_blocks) -> int:
        for (r_start, r_end), (p_start, _) in zip(ref_blocks, pat_blocks):
            if r_start <= pos < r_end:
                return p_start + (pos - r_start)
        for (r_start, r_end), (p_start, _) in zip(ref_blocks, pat_blocks):
            if pos < r_start:
                lib.log(
                    f"CDS boundary at ref pos {pos} falls;"
                    f"Snapping to nearest block"
                )
                return p_start
        return pat_blocks[-1][1]

    def _map_location_to_patient(self, loc, ref_blocks, pat_blocks) -> int:
        p_start = int(self._map_ref_pos_to_patient(int(loc.start), ref_blocks, pat_blocks))
        p_end = int(self._map_ref_pos_to_patient(int(loc.end) - 1, ref_blocks, pat_blocks)) + 1

        return type(loc)(p_start, p_end, strand=loc.strand)

    def _build_nucleotide_origin_map(self, locs) -> list:
        origin = []
        for loc in locs:
            indices = list(range(int(loc.start), int(loc.end)))
            if loc.strand == -1:
                indices.reverse()
            origin.extend(indices)

        return origin

    def _get_codon_info(self, cds_seq, cds_origin, amino_pos_0based) -> tuple[str, str]:
        nt_start = amino_pos_0based * 3
        nt_end = nt_start + 3

        if amino_pos_0based < 0 or nt_end > len(cds_seq) or nt_end > len(cds_origin):
            return None, None

        codon = str(cds_seq[nt_start:nt_end])
        origin_slice = cds_origin[nt_start:nt_end]

        return codon, origin_slice

    def _trim_common_flanks(
            self, ref_str, alt_str, ref_origin_slice, alt_origin_slice) -> tuple[str, str, int, int]:
        if ref_str is None or alt_str is None:
            return ref_str, alt_str, ref_origin_slice, alt_origin_slice

        start = 0
        while start < len(ref_str) and start < len(alt_str) and ref_str[start] == alt_str[start]:
            start += 1

        end_ref = len(ref_str)
        end_alt = len(alt_str)
        while end_ref > start and end_alt > start and ref_str[end_ref - 1] == alt_str[end_alt - 1]:
            end_ref -= 1
            end_alt -= 1

        trimmed_ref = ref_str[start:end_ref]
        trimmed_alt = alt_str[start:end_alt]
        trimmed_ref_origin = ref_origin_slice[start:end_ref] if ref_origin_slice else ref_origin_slice
        trimmed_alt_origin = alt_origin_slice[start:end_alt] if alt_origin_slice else alt_origin_slice

        return trimmed_ref, trimmed_alt, trimmed_ref_origin, trimmed_alt_origin

    def _to_genomic_representation(self, nt_str, origin_slice, gene_strand) -> tuple[str, int]:
        if nt_str is None or not origin_slice:
            return nt_str, None

        genomic_pos = min(origin_slice)
        genomic_str = str(Seq(nt_str).reverse_complement()) if gene_strand == -1 else nt_str

        return genomic_str, genomic_pos
