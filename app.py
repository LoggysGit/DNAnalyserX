import os

import queue
import threading

import modules.lib as lib
import modules.data_manager as dataManager
import modules.engine  as engine
import modules.interface as gui

gui_command_buffer = queue.Queue()
sys_command_buffer = queue.Queue()

data_manager = dataManager.DataManager(lib.DB_PATH, gui_command_buffer)
core = engine.Core(gui_command_buffer, data_manager)
app = gui.App(gui_command_buffer, sys_command_buffer, data_manager)

last_file_path = ""

def system_thread():
    while True:
        data_manager.handle_disease_db_update()

        command, payload = sys_command_buffer.get()

        match command:
            case "RUN":
                lib.log("Analysing started...")

                file_path, chrm, pos = payload

                # Open files
                if file_path.endswith(".gz"): selected_file_data = lib.gzip_open(file_path)
                else: selected_file_data = lib.open_file(file_path, "r")
                chr_data = data_manager.download_chromosome(chrm)

                last_file_path = file_path

                # Check
                if len(selected_file_data) > lib.MAX_NUCL_LENGTH:
                    lib.log(f"Selected file too long (>{lib.MAX_NUCL_LENGTH}). OverflowError.")
                    gui_command_buffer.put(("DONE", None))
                    break

                # Anaalyze & Send
                results = core.run_comparing(selected_file_data, chr_data, int(pos))

                # Purge temp
                data_manager.purge_temp()

                # Seek for diseases
                full_mutations_data = core.find_mutations(results, chrm)
                i = 0
                for d in full_mutations_data:
                    i += 1
                    d.insert(0, i)
                    gui_command_buffer.put(("MUTATION", d))

                gui_command_buffer.put(("DONE", None))
                lib.log("Analysing has ended.")

            case "EXPORT":
                mut_list, chr_ref = payload
                exp_path = os.path.join(lib.EXPORT_DIR, f"{last_file_path}_export.vcf")
                data_manager.save_mutations_to_vcf(exp_path, mut_list, chr_ref)

            case _: pass
        
        sys_command_buffer.task_done()

if __name__ == "__main__":
    sys_thread = threading.Thread(target=system_thread, daemon=True)
    sys_thread.start()

    app.mainloop()