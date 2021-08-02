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
import os.path

import matplotlib
import matplotlib.figure
import matplotlib.backends.backend_agg

import pickle
import os.path
import googleapiclient.http
import googleapiclient.discovery
# import google_auth_oauthlib.flow
import google.auth.transport.requests

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

    def upload_from_io(self, rel_filepath, mime_type, iostream):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        mode = 'wb' if iostream is io.StringIO else 'w'
        with open(abs_path, mode) as f:
            f.write(iostream.getvalue())

    def upload_png_from_figure(self, rel_filepath, figure):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        matplotlib.backends.backend_agg.FigureCanvasAgg(figure).print_png(abs_path)

class GoogleDriveUploader(Uploader):
    def __init__(self, token_file, root_folder_id, credentials_file = None):
        self.scopes = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
        self.root_folder_id = root_folder_id
        self.token_file = os.path.expanduser(token_file)
        self.credentials_file = os.path.expanduser(credentials_file)
        self.creds = None
        self.directory = {}
        self.auth()

    def auth(self):
        # The file token.pickle stores the user's access and refresh tokens, and is
        # created automatically when the authorization flow completes for the first
        # time.
        if os.path.exists(self.token_file):
            with open(self.token_file, 'rb') as token:
                self.creds = pickle.load(token)

        # If there are no (valid) credentials available, let the user log in.
        if not self.creds or not self.creds.valid:
            if self.creds and self.creds.expired and self.creds.refresh_token:
                self.creds.refresh(google.auth.transport.requests.Request())
            elif self.credentials_file and os.path.exists(self.credentials_file):
                flow = google_auth_oauthlib.flow.InstalledAppFlow.from_client_secrets_file(
                    self.credentials_file, self.scopes)
                creds = flow.run_local_server(port=0)
            else:
                # flow = google_auth_oauthlib.flow.InstalledAppFlow.from_client_secrets_file(
                #     'credentials.json', SCOPES)
                # creds = flow.run_local_server(port=0)
                raise RuntimeError('GoogleDriveUploader: could not find valid access token')

            # Save the credentials for the next run
            with open(self.token_file, 'wb') as token:
                pickle.dump(self.creds, token)

        self.drive_service = googleapiclient.discovery.build('drive', 'v3', credentials=self.creds)

    def make_path(self, rel_path):
        if(not rel_path):
            return self.root_folder_id
        rel_path = os.path.normpath(rel_path)
        if(rel_path.startswith('../')):
            raise RuntimeError('Cannot make path outside of base : '+rel_path)
        if(rel_path not in self.directory):
            (head, tail) = os.path.split(rel_path)
            parent = self.make_path(head)
            response = self.drive_service.files().list(\
                spaces='drive',
                fields='files(id, name)',
                q="name='%s' and '%s' in parents and mimeType='application/vnd.google-apps.folder'"%(tail,parent)).execute()
            files = response.get('files', [])
            if(files):
                self.directory[rel_path] = files[0].get('id')
            else:
                response = self.drive_service.files().create(\
                    body={ \
                        'name' : tail,
                        'mimeType' : 'application/vnd.google-apps.folder',
                        'parents' : [ parent ] },
                    fields='id').execute()
                self.directory[rel_path] = response.get('id')
        return self.directory[rel_path]

    def upload_from_io(self, rel_filepath, mime_type, iostream):
        (rel_path, filename) = os.path.split(rel_filepath)
        parent = self.make_path(rel_path)
        response = self.drive_service.files().list(\
            spaces='drive',
            fields='files(id, name)',
            q="name='%s' and '%s' in parents"%(filename,parent)).execute()
        files = response.get('files', [])
        media = googleapiclient.http.MediaIoBaseUpload(iostream, mimetype=mime_type)
        if(files):
            file_metadata = { \
                'mimeType' : mime_type }
            response = self.drive_service.files().update(\
                fileId=files[0].get('id'),
                body=file_metadata,
                media_body=media,
                fields='id').execute()
            pass
        else:
            file_metadata = { \
                'name': filename,
                'mimeType' : mime_type,
                'parents' : [ parent ]}
            response = self.drive_service.files().create(\
                body=file_metadata,
                media_body=media,
                fields='id').execute()
        return response.get('id')

    def upload_png_from_figure(self, rel_filepath, figure):
        canvas = matplotlib.backends.backend_agg.FigureCanvas(figure)
        output = io.BytesIO()
        canvas.print_png(output)
        return self.upload_from_io(rel_filepath, 'image/png', output)
