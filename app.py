import queue
import threading

import modules.lib as lib
import modules.data_manager as dataManager
import modules.engine as engine
import modules.interface as gui

gui_command_buffer = queue.Queue()
sys_command_buffer = queue.Queue()

data_manager = dataManager.DataManager(lib.DB_PATH, gui_command_buffer)
core = engine.Core(gui_command_buffer, data_manager)
app = gui.App(gui_command_buffer, sys_command_buffer, data_manager)

def system_thread():
    while True:
        data_manager.handle_disease_db_update()

        command, payload = sys_command_buffer.get()

        match command:
            case "RUN":
                lib.log("Analysing started...")

                # Get payload
                data_file_path, gene_id = payload

                # Anaalyze
                results = core.run_comparing(data_file_path, gene_id)

                # Seek for diseases
                full_mutations_data = core.find_mutations(results, gene_id)
                full_mutations_data.sort(key=lambda r: r[0])

                # Send into interface
                for i in range(len(full_mutations_data)):
                    d = full_mutations_data[i]
                    gui_command_buffer.put(("MUTATION", d))

                gui_command_buffer.put(("DONE", None))
                lib.log("Analysing has ended.")

            case "EXPORT":
                mut_list, gene_ref, exp_dir = payload
                data_manager.save_mutations_to_vcf(exp_dir, mut_list, gene_ref)

            case _: pass
        
        sys_command_buffer.task_done()

if __name__ == "__main__":
    sys_thread = threading.Thread(target=system_thread, daemon=True)
    sys_thread.start()

    app.mainloop()