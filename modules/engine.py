import modules.lib as lib

class Core:
    def __init__(self, gui_cmd_buff, dman):
        self.data_manager = dman
        self.gui_command_buffer = gui_cmd_buff

    def run(self, patient_file, reference_file, chromosome, position):
        results = []

        # --- Prepare data --- #
        # Patient DNA sector
        raw_lines = patient_file.readlines()
        clean_lines = [line.strip().upper() for line in raw_lines if not line.startswith(">")]
        patient_seq = "".join(clean_lines)

        seq_len = len(patient_seq)

        # Reference DNA sector
        reference_file.seek(0)
        
        letters_read = 0
        ref_buffer = []
        target_reached = False
        for line in reference_file:
            line = line.strip().upper()
            # Skip metadata and empty lines
            if not line or line.startswith(">"): continue
            line_len = len(line)
            
            if not target_reached and letters_read + line_len >= position:
                target_reached = True
                offset = position - letters_read - 1
                ref_buffer.append(line[offset:])
            elif target_reached: ref_buffer.append(line)
                
            letters_read += line_len
            
            if target_reached and sum(len(s) for s in ref_buffer) >= seq_len: break
                
        reference_sector = "".join(ref_buffer)[:seq_len]

        # --- Main check cycle --- #
        for n in range(seq_len):
            if patient_seq[n] != reference_sector[n]:
                results.append([position + n, reference_sector[n], patient_seq[n]])

        return results
    
    def find_mutations(self, dna_anomalies, chromosome_id):
        diseases = []
        
        return diseases