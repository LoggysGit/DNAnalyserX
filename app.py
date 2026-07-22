""" Main app file """

import queue
import threading

import modules.lib as lib
import modules.data_manager as data_manager
import modules.engine as engine
import modules.interface as gui

gui_command_buffer = queue.Queue()
sys_command_buffer = queue.Queue()

dm = data_manager.DataManager(lib.DB_PATH, gui_command_buffer)
core = engine.Core(gui_command_buffer, dm)
app = gui.App(gui_command_buffer, sys_command_buffer, dm)

def system_thread(swf):
    """ System thread"""
    while swf.is_set():
        dm.handle_disease_db_update()

        command, payload = sys_command_buffer.get()

        match command:
            case "RUN":
                lib.log("Analysing started...")

                # Get payload
                data_file_path, gene_id = payload

                # Analyse
                results = core.run_comparing(data_file_path, gene_id)

                # Seek for diseases
                full_mutations_data = core.find_mutations(results, gene_id)
                full_mutations_data.sort(key=lambda r: r[0])

                # Send into interface
                for di in full_mutations_data:
                    gui_command_buffer.put(("MUTATION", di))

                gui_command_buffer.put(("DONE", None))
                lib.log("Analysing has ended.")

            case "EXPORT":
                mut_list, gene_ref, exp_dir = payload
                dm.save_mutations_to_vcf(exp_dir, mut_list, gene_ref)

            case "CLOSE":
                swf.clear()

            case _: pass
        sys_command_buffer.task_done()

if __name__ == "__main__":
    # Set main working flag
    system_work_flag = threading.Event()
    system_work_flag.set()

    # Start system daemon
    sys_thread = threading.Thread(
    target=system_thread,
    args=(system_work_flag,),
    daemon=True
    )
    sys_thread.start()

    # Interface loop
    def on_closing():
        """ All-threads close function """
        lib.log("App closed.")
        system_work_flag.clear()
        app.destroy()

    app.protocol("WM_DELETE_WINDOW", on_closing)
    app.mainloop()
