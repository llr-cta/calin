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
import sys
import os.path
import time
import fcntl

import matplotlib
import matplotlib.figure
import matplotlib.backends.backend_agg

import pickle
import os.path

import google.auth.transport.requests
import google.oauth2.credentials
import google_auth_oauthlib.flow
import googleapiclient.discovery
import googleapiclient.errors
import googleapiclient.http
import socket

class Uploader:
    def __init__(self, overwrite=True, loud=False):
        self.loud = loud
        self.overwrite = overwrite
        pass

    def do_single_upload_from_io(self, rel_filepaths, mime_type, iostream):
        raise RuntimeError('do_single_upload_from_io: unimplemented in base class')

    def upload_from_io(self, rel_filepaths, mime_type, iostream):
        if(type(rel_filepaths) is not list):
            rel_filepaths = [ rel_filepaths ]
        for rel_filepath in rel_filepaths:
            self.do_single_upload_from_io(rel_filepath, mime_type, iostream)

    def upload_png_from_figure(self, rel_filepaths, figure):
        canvas = matplotlib.backends.backend_agg.FigureCanvas(figure)
        output = io.BytesIO()
        canvas.print_png(output)
        return self.upload_from_io(rel_filepaths, 'image/png', output)

    def retrieve_sheet(self, sheet_id, row_start=0):
        raise RuntimeError('retrieve_sheet: unimplemented in base class')

    def append_row_to_sheet(self, sheet_id_and_tab_name, row, row_start=0):
        raise RuntimeError('retrieve_sheet: unimplemented in base class')

    def get_url(self, rel_filepath):
        raise RuntimeError('get_url: unimplemented in base class')

class FilesystemUploader(Uploader):
    def __init__(self, root_directory, overwrite=True, loud=False):
        self.root_directory = os.path.normpath(os.path.expanduser(root_directory)) if root_directory else '.'
        if(not os.path.isdir(self.root_directory)):
            raise RuntimeError('Base path os not directory : '+self.root_directory)
        super().__init__(overwrite=overwrite,loud=loud)

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

    def do_single_upload_from_io(self, rel_filepath, mime_type, iostream):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        mode = 'wb' if isinstance(iostream, io.BytesIO) else 'w'
        if(os.path.exists(abs_path)):
            if(self.overwrite):
                if(self.loud):
                    print("Updating:",rel_filepath)
            else:
                if(self.loud):
                    print("Skipping:",rel_filepath)
                return None
        else:
            if(self.loud):
                print("Uploading:",rel_filepath)
        with open(abs_path, mode) as f:
            f.write(iostream.getvalue())
        return abs_path

    def get_url(self, rel_filepath):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        if(os.path.exists(abs_path)):
            return "file://"+abs_path
        else:
            return ''

class GoogleDriveUploader(Uploader):
    def __init__(self, token_file, root_folder_id, credentials_file=None,
            overwrite=True, loud=False):
        self.ordinal = ["zeroth", "first", "second", "third", "fourth", "fifth",
            "sixth", "seventh", "eigth","ninth","tenth"]
        self.scopes = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
        self.root_folder_id = root_folder_id
        self.token_file = os.path.expanduser(token_file)
        self.credentials_file = os.path.expanduser(credentials_file)
        self.creds = None
        self.directory = {}
        self.drive_service = None
        self.sheets_service = None
        self.lockfile = open(self.token_file+".lock",'ab')
        self.lockcount = 0
        self.auth()
        super().__init__(overwrite=overwrite,loud=loud)

    def get_drive_service(self):
        return self.drive_service

    def get_sheets_service(self):
        return self.sheets_service

    def lock(self):
        if(self.lockcount == 0):
            fcntl.lockf(self.lockfile, fcntl.LOCK_EX)
        self.lockcount += 1

    def unlock(self):
        self.lockcount -= 1
        if(self.lockcount == 0):
            fcntl.lockf(self.lockfile, fcntl.LOCK_UN)

    def auth(self):
        self.lock()
        try:
            # The file token.json stores the user's access and refresh tokens, and is
            # created automatically when the authorization flow completes for the first
            # time.
            try:
                if os.path.exists(self.token_file):
                    self.creds = google.oauth2.credentials.Credentials.from_authorized_user_file(self.token_file, self.scopes)
            except:
                pass

            # If there are no (valid) credentials available, let the user log in.
            if not self.creds or not self.creds.valid:
                if self.creds and self.creds.expired and self.creds.refresh_token:
                    self.creds.refresh(google.auth.transport.requests.Request())
                elif self.credentials_file and os.path.exists(self.credentials_file):
                    flow = google_auth_oauthlib.flow.InstalledAppFlow.from_client_secrets_file(
                        self.credentials_file, self.scopes)
                    self.creds = flow.run_local_server(port=0)
                else:
                    raise RuntimeError('GoogleDriveUploader: could not find valid access token')

                # Save the credentials for the next run
                with open(self.token_file, 'w') as token:
                    token.write(self.creds.to_json())
        except:
            self.unlock()
            raise

        self.unlock()
        self.drive_service = googleapiclient.discovery.build('drive', 'v3', credentials=self.creds)
        self.sheets_service = googleapiclient.discovery.build('sheets', 'v4', credentials=self.creds)

    def make_path(self, rel_path, do_create = True):
        if(not rel_path):
            return self.root_folder_id
        rel_path = os.path.normpath(rel_path)
        if(rel_path.startswith('../')):
            raise RuntimeError('Cannot make path outside of base : '+rel_path)

        if(rel_path not in self.directory):
            if(do_create):
                self.lock()

            try:
                (head, tail) = os.path.split(rel_path)
                parent = self.make_path(head, do_create)
                if(not parent):
                    if(do_create == False):
                        return ''
                    else:
                        raise RuntimeError('Parent was not created')
                response = self.drive_service.files().list(\
                    spaces='drive',
                    fields='files(id, name)',
                    q="name='%s' and '%s' in parents and trashed=false and mimeType='application/vnd.google-apps.folder'"%(tail,parent)).execute()
                files = response.get('files', [])
                if(files):
                    self.directory[rel_path] = files[0].get('id')
                elif(do_create):
                    response = self.drive_service.files().create(\
                        body={ \
                            'name' : tail,
                            'mimeType' : 'application/vnd.google-apps.folder',
                            'parents' : [ parent ] },
                        fields='id').execute()
                    self.directory[rel_path] = response.get('id')
                else:
                    # do_create is false, so no need to unlock
                    return ''
            except:
                if(do_create):
                    self.unlock()
                raise()

            if(do_create):
                self.unlock()

        return self.directory[rel_path]

    def do_single_upload_from_io(self, rel_filepath, mime_type, iostream):
        (rel_path, filename) = os.path.split(rel_filepath)
        parent = self.make_path(rel_path)
        response = self.drive_service.files().list(\
            spaces='drive',
            fields='files(id, name)',
            q="name='%s' and '%s' in parents and trashed=false"%(filename,parent)).execute()
        files = response.get('files', [])
        media = googleapiclient.http.MediaIoBaseUpload(iostream, mimetype=mime_type)
        if(files):
            if(self.overwrite):
                if(self.loud):
                    print("Updating:",rel_filepath)
                file_metadata = { \
                    'mimeType' : mime_type }
                response = self.drive_service.files().update(\
                    fileId=files[0].get('id'),
                    body=file_metadata,
                    media_body=media,
                    fields='id').execute()
                return response.get('id')
            else:
                if(self.loud):
                    print("Skipping:",rel_filepath)
                return None
        else:
            if(self.loud):
                print("Uploading:",rel_filepath)
            file_metadata = { \
                'name': filename,
                'mimeType' : mime_type,
                'parents' : [ parent ]}
            response = self.drive_service.files().create(\
                body=file_metadata,
                media_body=media,
                fields='id').execute()
            return response.get('id')

    def upload_from_io(self, rel_filepaths, mime_type, iostream, max_try=5):
        if(type(rel_filepaths) is not list):
            rel_filepaths = [ rel_filepaths ]
        for rel_filepath in rel_filepaths:
            ntry = 0
            uploaded = False
            while(not uploaded):
                ntry += 1
                try:
                    self.do_single_upload_from_io(rel_filepath, mime_type, iostream)
                    uploaded = True
                except (googleapiclient.errors.HttpError, socket.timeout):
                    if(ntry<max_try):
                        if(ntry<len(self.ordinal)):
                            print("Upload failed on %s attempt, trying again"%self.ordinal[ntry], file=sys.stderr)
                        else:
                            print("Upload failed on attempt %d, trying again"%ntry, file=sys.stderr)
                        time.sleep(min(2**ntry,100))
                    else:
                        if(ntry<len(self.ordinal)):
                            print("Upload failed on %s and final attempt"%self.ordinal[ntry], file=sys.stderr)
                        else:
                            print("Upload failed on final attempt %d"%ntry, file=sys.stderr)
                        raise

    def get_sheet_id_and_tab_name(self, sheet_id_and_tab_name):
        bits = sheet_id_and_tab_name.split('#')
        if(len(bits) <= 1):
            return sheet_id_and_tab_name, ''
        elif(len(bits) == 2):
            return bits[0], bits[1]
        else:
            raise RuntimeError("Could not understand sheet and tab specification: "+sheet_id_and_tab_name)

    def retrieve_sheet(self, sheet_id_and_tab_name, row_start=0, max_try=2):
        sheet_id, range = self.get_sheet_id_and_tab_name(sheet_id_and_tab_name)
        if range:
            range = "'" + range + "'!"
        range += 'A%d:ZZZ'%(row_start+1)
        ntry = 0
        retrieved = False
        while(not retrieved):
            ntry += 1
            try:
                response = self.sheets_service.spreadsheets().values().get(
                    spreadsheetId=sheet_id,range=range).execute()
                retrieved = True
            except googleapiclient.errors.HttpError:
                if(ntry<max_try):
                    if(ntry<len(self.ordinal)):
                        print("Failed to retrieve sheet on %s attempt, trying again"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to retrieve sheet on attempt %d, trying again"%ntry, file=sys.stderr)
                    time.sleep(min(2**ntry,100))
                else:
                    if(ntry<len(self.ordinal)):
                        print("Failed to retrieve sheet on %s and final attempt"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to retrieve sheet on final attempt %d"%ntry, file=sys.stderr)
                    raise

        if response and 'values' in response:
            return response['values']
        else:
            return []

    def append_rows_to_sheet(self, sheet_id_and_tab_name, rows, row_start=0, max_try=2):
        sheet_id, range = self.get_sheet_id_and_tab_name(sheet_id_and_tab_name)
        if range:
            range = "'" + range + "'!"
        range += 'A%d:ZZZ'%(row_start+1)
        ntry = 0
        added = False
        while(not added):
            ntry += 1
            try:
                body = {
                    'values': rows
                }
                response = self.sheets_service.spreadsheets().values().append(
                    spreadsheetId=sheet_id, range=range,
                    valueInputOption='USER_ENTERED', body=body).execute()
                added = True
            except googleapiclient.errors.HttpError:
                if(ntry<max_try):
                    if(ntry<len(self.ordinal)):
                        print("Failed to append to sheet on %s attempt, trying again"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to append to sheet on attempt %d, trying again"%ntry, file=sys.stderr)
                    time.sleep(min(2**ntry,100))
                else:
                    if(ntry<len(self.ordinal)):
                        print("Failed to append to sheet on %s and final attempt"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to append to sheet on final attempt %d"%ntry, file=sys.stderr)
                    raise

        return response.get('updates').get('updatedCells')

    def append_row_to_sheet(self, sheet_id_and_tab_name, row, row_start=0, max_try=2):
        return self.append_rows_to_sheet(sheet_id_and_tab_name, [row], row_start, max_try)

    def get_url(self, rel_filepath):
        if(rel_filepath in self.directory):
            id = self.directory[rel_filepath]
            response = self.drive_service.files().get(fileId=id,
                fields='webViewLink').execute()
            return response.get('webViewLink')
        else:
            (rel_path, filename) = os.path.split(rel_filepath)
            parent = self.make_path(rel_path, do_create = False)
            if parent:
                response = self.drive_service.files().list(\
                    spaces='drive',
                    fields='files(webViewLink)',
                    q="name='%s' and '%s' in parents and trashed=false"%(filename,parent)).execute()
                files = response.get('files', [])
                for file in response.get('files', []):
                    return file.get('webViewLink')
        return ''

    def get_sheet_tab_dict(self, sheet_id, max_try=2):
        done = False
        tabs = dict()
        ntry = 0
        while(not done):
            ntry += 1
            try:
                sheet_metadata = self.sheets_service.spreadsheets().get(spreadsheetId=sheet_id,
                    fields='sheets(properties(title,sheetId))').execute()

                for sheet in sheet_metadata.get('sheets'):
                    tabs[sheet.get("properties").get('title')] = \
                        sheet.get("properties").get('sheetId')
                done = True
            except googleapiclient.errors.HttpError:
                if(ntry<max_try):
                    if(ntry<len(self.ordinal)):
                        print("Failed to get sheet tabs on %s attempt, trying again"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to get sheet tabs on attempt %d, trying again"%ntry, file=sys.stderr)
                    time.sleep(min(2**ntry,100))
                else:
                    if(ntry<len(self.ordinal)):
                        print("Failed to get sheet tabs on %s and final attempt"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to get sheet tabs on final attempt %d"%ntry, file=sys.stderr)
                    raise
        return tabs

    def get_sheet_tab_ids(self, sheet_id, max_try=2):
        done = False
        tabs = []
        ntry = 0
        while(not done):
            ntry += 1
            try:
                sheet_metadata = self.sheets_service.spreadsheets().get(spreadsheetId=sheet_id,
                    fields='sheets(properties(title,sheetId))').execute()

                for sheet in sheet_metadata.get('sheets'):
                    tabs.append(sheet.get("properties").get('sheetId'))
                done = True
            except googleapiclient.errors.HttpError:
                if(ntry<max_try):
                    if(ntry<len(self.ordinal)):
                        print("Failed to get sheet tabs on %s attempt, trying again"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to get sheet tabs on attempt %d, trying again"%ntry, file=sys.stderr)
                    time.sleep(min(2**ntry,100))
                else:
                    if(ntry<len(self.ordinal)):
                        print("Failed to get sheet tabs on %s and final attempt"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to get sheet tabs on final attempt %d"%ntry, file=sys.stderr)
                    raise
        return tabs

    def sort_sheet(self, sheet_id_and_tab_name, sort_column, ascending_order = True,
            row_start=0, max_try=2):
        sheet_id, tab_name = self.get_sheet_id_and_tab_name(sheet_id_and_tab_name)
        ntry = 0
        done = False
        tab_id = ''
        while(not done):
            ntry += 1
            try:
                if(tab_name):
                    tabs = self.get_sheet_tab_dict(sheet_id)
                    if(tab_name in tabs):
                        tab_id = tabs[tab_name]
                    else:
                        raise RuntimeError('Sheet not found ' + tab_name)
                else:
                    tab_id = self.get_sheet_tab_ids(sheet_id)[0]

                self.sheets_service.spreadsheets().batchUpdate(spreadsheetId=sheet_id,
                    body={
                        'requests' : [
                            {
                                'sortRange' : {
                                    'range' : {
                                        'sheetId' : tab_id,
                                        'startRowIndex' : row_start
                                    },
                                    'sortSpecs' : [
                                        {
                                            "sortOrder" : "ASCENDING" if ascending_order else "DESCENDING",
                                            "dimensionIndex": sort_column
                                        }
                                    ]
                                }
                            }
                        ]
                    }).execute()

                done = True
            except googleapiclient.errors.HttpError:
                if(ntry<max_try):
                    if(ntry<len(self.ordinal)):
                        print("Failed to sort sheet on %s attempt, trying again"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to sort sheet on attempt %d, trying again"%ntry, file=sys.stderr)
                    time.sleep(min(2**ntry,100))
                else:
                    if(ntry<len(self.ordinal)):
                        print("Failed to sort sheet on %s and final attempt"%self.ordinal[ntry], file=sys.stderr)
                    else:
                        print("Failed to sort sheet on final attempt %d"%ntry, file=sys.stderr)
                    raise
