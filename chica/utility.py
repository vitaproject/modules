
from os import path


def get_example_data_path(relative_data_path):

    data_path = path.join(path.dirname(__file__), 'divertor_volume')
    data_path = path.join(data_path, relative_data_path)

    if not path.isfile(data_path):
        raise FileNotFoundError("The following example data file could not be located, "
                                "'{}'.".format(relative_data_path))

    return data_path
