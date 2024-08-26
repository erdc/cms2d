import os
import numpy
import filecmp
from typing import Optional
from h5_file_to_text import write_h5_as_text_file
from pathlib import Path
import process_runner
import shutil


def compare_h5_files(base_file: str, out_file: str, to_exclude: Optional[list[str]] = None,
                     allow_close: bool = False) -> str:
    """Compare two H5 files for testing.

    Args:
        base_file: The path to the base coverage file.
        out_file: The path to the output coverage file.
        to_exclude: Optional list of dataset or attribute names to exclude from printing.
        allow_close: Allow floating point values to be close but not equal.
    """
    if to_exclude is None:
        to_exclude = []
    base_txt_file = os.path.splitext(base_file)[0] + '.h5diff.txt'
    out_txt_file = os.path.splitext(out_file)[0] + '.h5diff.txt'
    base_writer = write_h5_as_text_file(base_file, base_txt_file, to_exclude=to_exclude)
    out_writer = write_h5_as_text_file(out_file, out_txt_file, to_exclude=to_exclude)
    equal = filecmp.cmp(base_txt_file, out_txt_file)
    if not equal and allow_close:
        base_item_names = list(base_writer.compare_items.keys())
        out_item_names = list(out_writer.compare_items.keys())
        base_item_names.sort()
        out_item_names.sort()
        if base_item_names != out_item_names:
            return (f'files differ:\n  {base_txt_file}\n  {out_txt_file}\n'
                    'groups, datasets, and/or attributes differ')
        for name in base_item_names:
            base_data = base_writer.compare_items[name]
            out_data = out_writer.compare_items[name]
            base_float = isinstance(base_data, numpy.ndarray) and base_data.dtype.kind == 'f'
            out_float = isinstance(out_data, numpy.ndarray) and out_data.dtype.kind == 'f'
            if base_float and out_float:
                equal = numpy.allclose(base_data, out_data, equal_nan=True)
                if not equal:
                    message = (f'files differ:\n  {base_txt_file}\n  {out_txt_file}\n'
                               f'{name} not close')
                    return message
            elif base_data != out_data:  # pragma no cover - only hit for differing groups which should fail above
                message = (f'files differ:\n  {base_txt_file}\n  {out_txt_file}\n'
                           f'{name} not equal')
                return message
        equal = True
    if not equal:
        return f'files differ:\n  {base_txt_file}\n  {out_txt_file}'
    return ''


def main(vs_solution_path, dll_path):
    """Run cms tests.

    Args:
        vs_solution_path (string): Path to the visual studio solution file.
        dll_path (string): path to dll for openmp.
    """
    repo_path = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
    test_path = os.path.join(repo_path, 'testing', 'tests')
    files = [str(f.absolute()) for f in Path('.').rglob('*.cmcards')]
    skip_files = [str(f.absolute()) for f in Path('.').rglob('skip.txt')]
    skip_set = set()
    for sf in skip_files:
        print(f'Skipping {sf}')
        skip_set.add(os.path.dirname(sf))
    files = [f for f in files if os.path.dirname(f) not in skip_set]

    # Copy 'libiomp5md.dll' and executable to tests directory
    shutil.copyfile(os.path.join(repo_path, vs_solution_path, 'x64/Release/CMS2d_v5.3.exe'),
                    os.path.join(test_path, './CMS2D_v5.3.exe'))
    shutil.copyfile(os.path.join(dll_path, 'libiomp5md.dll'),
                    os.path.join(test_path, './libiomp5md.dll'))

    runner = process_runner.ProcessRunner()
    runner.run_concurrent_processes(files)
    # test CMS outputs
    return_code = _check_outputs(files)
    if os.path.isfile(os.path.join(test_path, 'CMS2d_v5.3.exe')):
        os.remove(os.path.join(test_path, 'CMS2d_v5.3.exe'))
        os.remove(os.path.join(test_path, 'libiomp5md.dll'))

    return return_code


def _check_outputs(files):
    """Check the test outputs.

    Args:
        files (list(str)): list of files being tested

    Returns:
        return_code (int): value for the return code (0-good, 1-bad)
    """
    return_code = 0
    for f in files:
        os.chdir(os.path.dirname(f))
        basefiles = [str(f.absolute()) for f in Path('./base').rglob('*.h5')]
        for basefile in basefiles:
            outfile = os.path.normpath(os.path.join(os.path.dirname(basefile), '..', os.path.basename(basefile)))
            if not os.path.isfile(outfile):
                print(f'TEST FAILED: {f}\n')
                return_code = 1
            val = compare_h5_files(basefile, outfile, allow_close=True)
            if val != '':
                print(f'{val}\n')
                return_code = 1

    return return_code
