""" Main app file """

import queue
import threading

from modules import lib, data_manager, engine, interface

gui_command_buffer = queue.Queue()
sys_command_buffer = queue.Queue()

dm = data_manager.DataManager(lib.DB_PATH, gui_command_buffer)
core = engine.Core(gui_command_buffer, dm)
app = interface.App(gui_command_buffer, sys_command_buffer, dm)

def system_thread():
    """ System thread """
    while True:
        item = sys_command_buffer.get()

        try:
            command, payload = item

            match command:
                # Analysis command #
                case "RUN":
                    # Check DB relevance first
                    dm.handle_disease_db_update()

                    lib.log("Analysing started...")

                    # Compare
                    data_file_path, gene_id = payload
                    results = core.run_comparing(data_file_path, gene_id)

                    # Extract mutations
                    full_mutations_data = core.find_mutations(results, gene_id)
                    full_mutations_data.sort(key=lambda r: r[0])

                    # Send results
                    for di in full_mutations_data:
                        gui_command_buffer.put(("MUTATION", di))

                    # Mark for done
                    gui_command_buffer.put(("DONE", None))
                    lib.log("Analysing has ended.")

                # Export file command #
                case "EXPORT":
                    mut_list, gene_ref, exp_dir = payload
                    dm.save_mutations_to_vcf(exp_dir, mut_list, gene_ref)

                # Daemon close (break) command #
                case "CLOSE":
                    break

                # Unknown command #
                case _:
                    lib.dbg(f"Unknown command: {command}")

        except Exception as e:
            lib.log(f"System thread error: {e}")

        finally:
            sys_command_buffer.task_done()

if __name__ == "__main__":
    # Start system daemon
    sys_thread = threading.Thread(target=system_thread, daemon=True)
    sys_thread.start()

    def on_closing():
        """ All-threads close function """
        lib.log("App closing...")
        sys_command_buffer.put(("CLOSE", None))
        app.destroy()

    app.protocol("WM_DELETE_WINDOW", on_closing)
    app.mainloop()
