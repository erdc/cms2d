"""Export an H5 file to text."""
# 1. Standard python modules
from dataclasses import dataclass, field
import io
import sys
from typing import Dict, List, Optional, TextIO, Union

# 2. Third party modules
import h5py
import numpy

# 3. Aquaveo modules

# 4. Local modules

__copyright__ = "(C) Copyright Aquaveo 2020"
__license__ = "All rights reserved"


@dataclass
class H5FileWriter:
    """Write an H5 file to text."""
    input_file: str
    output_file: TextIO
    to_exclude: List[str] = field(default_factory=list)
    compare_items: Dict[str, object] = field(default_factory=dict)
    h5_file: Optional[h5py.File] = None
    bool_write_single_attribute: bool = False

    def print_h5(self):
        """Print the H5 file."""
        if h5py.is_hdf5(self.input_file):
            self.h5_file = h5py.File(self.input_file)
            self._visit_h5_group(self.h5_file)
        else:
            self.output_file.write(f'Invalid or missing file: {self.input_file}\n')

    def _visit_h5_group(self, h5_group: h5py.Group):
        """Print an H5 group to a text file.

        Args:
            h5_group: The group to print.
        """
        self._visit_attributes(h5_group)
        self._visit_datasets(h5_group)
        self._visit_groups(h5_group)

    def _visit_attributes(self, h5_group: h5py.Group):
        """Print H5 group attributes to a text file.

        Args:
            h5_group: The group to print.
        """
        for attribute_name in h5_group.attrs.keys():
            if attribute_name in self.to_exclude:
                continue
            attribute = h5_group.attrs[attribute_name]
            group_label = f'{h5_group.name} attribute: "{attribute_name}"'
            compare_data = self._print_data(attribute, group_label)
            self.compare_items[group_label] = compare_data

    def _visit_datasets(self, h5_group: h5py.Group):
        """Print H5 group datasets to a text file.

        Args:
            h5_group (hdf5.Group): The group to print.
        """
        for dataset_name in h5_group.keys():
            if dataset_name in self.to_exclude:
                continue
            dataset = h5_group[dataset_name]
            if isinstance(dataset, h5py.Dataset):
                compare_data = self._print_data(dataset, dataset.name)
                self.compare_items[dataset.name] = compare_data

    def _visit_groups(self, h5_group: h5py.Group):
        """Print H5 groups within a group to a text file.

        Args:
            h5_group (hdf5.Group): The group to print.
        """
        for group_name in h5_group.keys():
            if group_name in self.to_exclude:
                continue
            group = h5_group[group_name]
            if isinstance(group, h5py.Group):
                self.compare_items[f'group: {group.name}'] = None
                self.output_file.write(f'{group.name}\n')
                self._visit_h5_group(group)

    def _print_data(self, data: Union[h5py.Dataset, numpy.ndarray], name: str) -> object:
        """Print H5 data.

        Args:
            data: The data to print.
            name: The name of the dataset.

        Returns:
            Data that can later be used for an almost equal comparison.
        """
        if data.shape == ():
            print_output = f'{name}: {data.dtype} = {str(data[()])}\n'
        elif len(data.shape) == 1:
            shape = f'({data.shape[0]})'
            if data.dtype.kind == 'S':
                value = data[0].decode('UTF-8')
                print_output = f"{name}: string{data.dtype.itemsize}{shape} = '{value}'\n"
            else:
                print_output = f'{name}: {data.dtype}{shape} = {data[:]}\n'
        else:
            file = io.StringIO()
            file.write(f'{name}: {data.dtype}{data.shape} =\n')
            array_string = numpy.array2string(data[:], threshold=sys.maxsize)
            for line in array_string.splitlines():
                if line != '':
                    file.write(line)
                    file.write('\n')
            print_output = file.getvalue()
        self.output_file.write(print_output)
        if data.dtype.kind == 'f' and data.shape != ():
            compare_data = data[:]
        else:
            compare_data = print_output
        return compare_data

    def write_single_attribute(self, attr):
        """Prints an attribute and everything below it in the h5 file.

        Args:
            attr: the attribute being checked
        """
        self.to_exclude = []
        self.h5_file = h5py.File(self.input_file)
        for group in self.h5_file.keys():
            if group == attr:
                continue
            else:
                self.to_exclude.append(group)
        self.print_h5()


def write_h5_as_text_file(input_h5_file: str, output_text_file: str, to_exclude: List[str]) -> H5FileWriter:
    """Write an H5 file as a text file.

    Args:
        input_h5_file: Path to the H5 file.
        output_text_file: Path to the output text file.
        to_exclude: List of dataset or attribute names to exclude from printing.
    """
    with open(output_text_file, 'w') as output_file:
        writer = H5FileWriter(input_h5_file, output_file, to_exclude)
        writer.print_h5()
        return writer


def write_single_attr_to_text_file(input_h5_file: str, output_text_file: str, single_attribute: str) -> H5FileWriter:
    """Write a single attribute from an H5 file as a text file.

    Args:
        input_h5_file: Path to the H5 file.
        output_text_file: Path to the output text file.
        single_attribute: The attribute that is being printed
    """
    with open(output_text_file, 'w') as output_file:
        writer = H5FileWriter(input_h5_file, output_file, bool_write_single_attribute=True)
        writer.write_single_attribute(single_attribute)
        return writer
