import py7zr

def extract_7z(file_path, extract_path):
    with py7zr.SevenZipFile(file_path, 'r') as archive:
        archive.extractall(path=extract_path)