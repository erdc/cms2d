"""Does the build and test for project."""
# 1. Standard python modules
import os
import shutil
import subprocess
import time

# 2. Third party modules

# 3. Aquaveo modules

# 4. Local modules

__copyright__ = "(C) Copyright Aquaveo 2020"
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


def replace_file(src, dst):
    """Overwrite a file on disk.

    Args:
        src (str): Path to the source file
        dst (str): Path to the destination file
    """
    if os.path.isfile(dst):
        os.remove(dst)
    shutil.copyfile(src, dst)


def build(solution_path=None,
          build_config_names=None,
          project_name=None,
          exe_file_names=None,
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
        do_build:            Example: 0 or 1 (0 to skip build)
        do_test:             Example: 0 or 1 (0 to skip test)
        test_sub_directory:  Example: 'testing'
    """
    if not solution_path or not build_config_names or not project_name or not exe_file_names:
        return False

    start = time.time()

    devenv = r'"C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\devenv.com"'
    this_file_directory_path = os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir))
    test_subdir = os.path.normpath(os.path.join(this_file_directory_path, test_sub_directory))

    # build executables
    if do_build:
        for configuration in build_config_names:
            shell(f'{devenv} {solution_path} /Clean {configuration}')
            shell(f'{devenv} {solution_path} /Build {configuration}')
    else:
        print('Skipping build step.')

    # run tests
    return_code = 0
    if do_test != 0:
        old_cwd = os.getcwd()
        os.chdir(test_subdir)

        import test_cms
        return_code1 = test_cms.main()
        if return_code1 != 0:
            return_code = 1

        os.chdir(old_cwd)
    else:
        print('Skipping test step.')

    elapsed = time.time() - start
    print(f'Elapsed time: {elapsed} seconds.\n')
    return return_code
