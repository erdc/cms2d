"""Build and test srh pre."""
# 1. Standard python modules
import argparse
import json
import os

# 2. Third party modules

# 3. Aquaveo modules

# 4. Local modules
import builder as builder
import process_runner as pr

__copyright__ = "(C) Copyright Aquaveo 2020"
__license__ = "All rights reserved"


def main():
    """Entry point for the build script."""
    directory_path = os.path.normpath(os.path.dirname(__file__))
    with open(os.path.join(directory_path,'build-config.json')) as f:
        config = json.load(f)

    parser = argparse.ArgumentParser(description="builder.")
    parser.add_argument('-b', '--build', help='Enter 0 to skip build', required=False, default=0)
    parser.add_argument('-t', '--test', help='Enter 0 to skip tests.', required=False, default=1)
    parser.add_argument('-tf', '--test_folder', help='Enter directory to search for .cmcards files.', required=False,
                        default='testing/tests/L05_Implicit_Test01')
    parser.add_argument('-n_proc', '--number_concurrent_processes', help='Enter number of concurrent processes',
                        required=False, default=2)
    args = parser.parse_args()

    pr.ProcessRunner.number_concurrent_processes = int(args.number_concurrent_processes)

    return_code = builder.build(
        solution_path=config['solution_path'],
        build_config_names=config['build_config_names'],
        project_name=config['project_name'],
        exe_file_names=config['exe_file_names'],
        dev_env=config['dev_env'],
        dll_path=config['dll_path'],
        do_build=int(args.build),
        do_test=int(args.test),
        test_sub_directory=args.test_folder,
    )
    return return_code


if __name__ == "__main__":
    exit(main())
