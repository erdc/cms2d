"""Process runner class."""
# 1. Standard python modules
import json
import os
import subprocess
import time
import shutil

# 2. Third party modules

# 3. Aquaveo modules

# 4. Local modules

__copyright__ = "(C) Copyright Aquaveo 2024"
__license__ = "All rights reserved"

# index for tuple in self.running_processes
PROCESS_IDX = 0
START_TIME_IDX = 1
FILE_NAME_IDX = 2


class ProcessRunner:
    """Class to run concurrent processes."""
    number_concurrent_processes = 1
    cms_exe = ''

    def __init__(self):
        """Constructor."""
        self.running_processes = []
        self.elapsed_times = []
        self.return_code = 0
        self.proc_num = 0

    def run_concurrent_processes(self, files):
        """Runs an executable concurrently on the list of files.

        Args:
            files (list(pathlib.Path)): list of files to process
        """
        for f in files:
            os.chdir(os.path.dirname(f))

            cmd = ProcessRunner.cms_exe
            out_file = os.path.join(os.getcwd(), f'out_file_{self.proc_num}.txt')
            self.proc_num += 1
            out_file_object = open(out_file, 'w')
            process_cmd = f'"{cmd}" "{os.path.basename(f)}"'
            print(process_cmd)

            start_time = time.time()
            self.running_processes.append((
                subprocess.Popen(args=process_cmd, stdout=out_file_object, cwd=os.getcwd()),
                start_time,
                f
            ))
            print(f'Current directory: {os.getcwd()}')
            print(f'Executing: {process_cmd}\n')
            while len(self.running_processes) == self.number_concurrent_processes:
                self._remove_finished_process()
                time.sleep(.2)

        while len(self.running_processes) > 0:
            self._remove_finished_process()
            time.sleep(.2)

    def _remove_finished_process(self):
        """Removes a finished process from list of processes."""
        for idx, process_data in enumerate(self.running_processes):
            process = process_data[PROCESS_IDX]
            if process.poll() is not None:
                self.running_processes.pop(idx)
                return
