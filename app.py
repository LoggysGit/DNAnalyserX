import queue
import threading

import modules.lib as lib
import modules.data_manager as dataManager
import modules.engine  as engine
import modules.interface as gui

gui_command_buffer = queue.Queue()
sys_command_buffer = queue.Queue()

data_manager = dataManager.DataManager()
core = engine.Core(gui_command_buffer, data_manager)
app = gui.App(gui_command_buffer, sys_command_buffer, data_manager)

def system_thread():
    while True:
        command, payload = sys_command_buffer.get()

        match command:
            case "RUN":
                lib.log("Analysing started...")

                file_path, chrm, pos = payload
                # Open files
                if file_path.endswith(".gz"): selected_file = lib.gzip_open(file_path)
                else: selected_file = lib.open_file(file_path, "r")
                chr_file = data_manager.download_chromosome(chrm)
                # Anaalyze & Send
                results = core.run(selected_file, chr_file, chrm, int(pos))
                for r in results: gui_command_buffer.put(("DNA_ANOMALY", r))
                # Close files
                selected_file.close()
                chr_file.close()
                # Purge temp
                data_manager.purge_temp()
                # Seek for diseases
                diseases = core.find_mutations(results, chrm)
                for d in diseases: gui_command_buffer.put(("DISEASE", d))

                gui_command_buffer.put(("DONE", None))
                lib.log("Analysing has ended.")

            case _: pass
        
        sys_command_buffer.task_done()

if __name__ == "__main__":
    sys_thread = threading.Thread(target=system_thread, daemon=True)
    sys_thread.start()

    app.mainloop()