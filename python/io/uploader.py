# calin/python/io/uploader.py -- Stephen Fegan -- 2021-08-02
#
# Classes to upload files to data to various systems
#
# Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris
#
# This file is part of "calin"
#
# "calin" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License version 2 or later, as published by
# the Free Software Foundation.
#
# "calin" is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

import io
import os

import matplotlib
import matplotlib.figure
import matplotlib.backends.backend_agg

class Uploader:
    def __init__(self):
        pass

    def upload_from_io(self, rel_filepath, mime_type, iostream):
        raise RuntimeError('upload_from_io: unimplemented in base class')

    def upload_png_from_figure(self, rel_filepath, figure):
        raise RuntimeError('upload_png_from_figure: unimplemented in base class')


class FilesystemUploader(Uploader):
    def __init__(self, root_directory):
        self.root_directory = os.path.normpath(os.path.expanduser(root_directory)) if root_directory else '.'
        if(not os.path.isdir(self.root_directory)):
            raise RuntimeError('Base path os not directory : '+self.root_directory)

    def make_path(self, rel_path):
        if(not rel_path):
            return ''
        rel_path = os.path.normpath(rel_path)
        abs_path = os.path.normpath(os.path.join(self.root_directory, rel_path))
        if((self.root_directory == '.' and (abs_path.startswith('../') or abs_path.startswith('/')))
                or (self.root_directory != '.' and not abs_path.startswith(self.root_directory))):
            raise RuntimeError('Cannot make path outside of base : '+rel_path)
        if(not os.path.isdir(abs_path)):
            (head, tail) = os.path.split(rel_path)
            self.make_path(head)
            # print("mkdir",abs_path)
            os.mkdir(abs_path)
        return abs_path

    def upload_from_io(self, rel_filepath, iostream):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        mode = 'wb' if iostream is io.StringIO else 'w'
        with open(abs_path, mode) as f:
            f.write(iostream.getvalue())

    def upload_png_from_figure(self, rel_filepath, figure):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        matplotlib.backends.backend_agg.FigureCanvasAgg(figure).print_png(abs_path)

# class GoogleDriveUploader(Uploader):
