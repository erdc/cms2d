"""Does the build and test for project."""
# 1. Standard python modules
import os
import shutil
import subprocess
import time

# 2. Third party modules

# 3. Aquaveo modules

# 4. Local modules

__copyright__ = "(C) Copyright 2024"
__license__ = "All rights reserved"


def shell(command):
    """Run a terminal command.

    Args:
        command (str): The command to run
    """
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print('cwd:', os.getcwd())
        raise RuntimeError('Failed to run command: ' + command)


def build(solution_path=None,
          build_config_names=None,
          project_name=None,
          exe_file_names=None,
          dev_env=None,
          dll_path=None,
          do_build=None,
          do_test=None,
          test_sub_directory=None,
          ):
    r"""Builds the project.

    Args:
        solution_path:       Example: '..\Intel_vs2022\CMS2D_V5.3.sln'
        build_config_names:  Example: ('Release')
        project_name:        Example: 'CMS2D_V5.3'
        exe_file_names:      Example: ('CMS2D_v5.3.exe')
        dev_env:             Example: "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\devenv.com"
        dll_path:            Example: "c:/Program Files (x86)/Common Files/Intel/Shared Libraries/redist/intel64/compiler" Used for finding the 'libiomp5.dll'
        do_build:            Example: 0 or 1 (0 to skip build)
        do_test:             Example: 0 or 1 (0 to skip test)
        test_sub_directory:  Example: 'testing'
    """
    flags = [solution_path, build_config_names, project_name, exe_file_names, dev_env, dll_path]
    if any(item is None for item in flags):
        return False

    start = time.time()

    directory_path = os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir))
    repo_dir = os.path.normpath(os.path.join(directory_path, '..', '..'))
    test_subdir = os.path.join(repo_dir, test_sub_directory)
    abs_path_solution = os.path.join(repo_dir, solution_path)

    # build executables
    if do_build:
        for configuration in build_config_names:
            shell(f'{dev_env} {abs_path_solution} /Clean {configuration}')
            shell(f'{dev_env} {abs_path_solution} /Build {configuration}')
    else:
        print('Skipping build step.')

    # run tests
    return_code = 0
    if do_test != 0:
        old_cwd = os.getcwd()
        os.chdir(test_subdir)

        import test_cms
        return_code1 = test_cms.main(os.path.dirname(abs_path_solution), dll_path)
        if return_code1 != 0:
            return_code = 1

        os.chdir(old_cwd)
    else:
        print('Skipping test step.')

    elapsed = time.time() - start
    print(f'Elapsed time: {elapsed} seconds.\n')
    return return_code
