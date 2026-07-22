""" Main app file """

import threading
from queue import Queue

from modules import lib, data_manager, engine, interface

def system_thread(
        sys_buffer: Queue,
        gui_buffer: Queue,
        data_manager_obj: data_manager.DataManager,
        core_engine: engine.Core ):
    """ System thread """

    while True:
        item = sys_buffer.get()

        try:
            command, payload = item

            match command:
                # Analysis command #
                case "RUN":
                    # Check DB relevance first
                    data_manager_obj.handle_disease_db_update()

                    lib.log("Analysing started...")

                    # Compare
                    data_file_path, gene_id = payload
                    results = core_engine.run_comparing(data_file_path, gene_id)

                    # Extract mutations
                    full_mutations_data = core_engine.find_mutations(results, gene_id)
                    full_mutations_data.sort(key=lambda r: r[0])

                    # Send results
                    for di in full_mutations_data:
                        gui_buffer.put(("MUTATION", di))

                    # Mark for done
                    gui_buffer.put(("DONE", None))
                    lib.log("Analysing has ended.")

                # Export file command #
                case "EXPORT":
                    mut_list, gene_ref, exp_dir = payload
                    data_manager_obj.save_mutations_to_vcf(exp_dir, mut_list, gene_ref)

                # Daemon close (break) command #
                case "CLOSE":
                    break

                # Unknown command #
                case _:
                    lib.dbg(f"Unknown command: {command}")

        except ValueError:
            lib.log(f"Wrong command structure (CMD, PAYLOAD): {item}")

        except Exception as e:
            lib.log(f"System thread error: {e}")

        finally:
            sys_buffer.task_done()

if __name__ == "__main__":
    # App command queues
    sys_command_buffer = Queue()
    gui_command_buffer = Queue()

    # App objects
    dm = data_manager.DataManager(lib.DB_PATH, gui_command_buffer)
    core = engine.Core(gui_command_buffer, dm)
    app = interface.App(gui_command_buffer, sys_command_buffer, dm)

    # Start system daemon
    sys_thread = threading.Thread(
        target=system_thread,
        args=(sys_command_buffer, gui_command_buffer, dm, core),
        daemon=True)
    sys_thread.start()

    # Start app loop
    def on_closing():
        """ All-threads close function """
        lib.log("App closing...")
        sys_command_buffer.put(("CLOSE", None))
        app.destroy()

    app.protocol("WM_DELETE_WINDOW", on_closing)
    app.mainloop()
