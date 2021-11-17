
import urllib.request
import zipfile
import pathlib
from html.parser import HTMLParser




class HMDB_Downloader:

    def get_file_name_from_url(self, url):
        url_path = pathlib.PurePosixPath(urllib.parse.urlparse(url).path)
        return  url_path.name.name

    def run(self,directory):

        URL_SPECTRA = 'http://specdb.wishartlab.com/downloads/exports/spectra_xml/hmdb_nmr_spectra.zip'
        URL_METABOLITES = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'

        spectra_name =self.get_file_name_from_url(URL_SPECTRA).stem
        target_directory =  pathlib.Path(directory,spectra_name)

        self.download_and_extract_xml_zip(URL_SPECTRA, target_directory)
        self.download_and_extract_xml_zip(URL_METABOLITES, target_directory)



    def download_and_extract_xml_zip(self, URL, directory):

        file = self.get_path_from_url(URL)
        zip_target = pathlib.Path(directory, file)

        target = pathlib.Path(directory)

        print(f'Beginning download of {file}...')
        print(f"    target directory is {target}")

        urllib.request.urlretrieve(URL, zip_target)

        zip_file = pathlib.Path(directory, file)
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(target)

        print(f"download and extract complete")
        print(f"    removing zip file  {zip_target}")

        zip_target.unlink


class BMRB_Downloader:
    class BMRB_Directory_Parser(HTMLParser):

        def __init__(self):
            super().__init__()
            self.files = []

        def handle_starttag(self, tag, attrs):
            if tag == 'a':
                href  = [elem for elem in attrs if elem[0] == 'href'][0][1]
                if href.startswith('bmst'):
                    file_name = f'{attrs[0][1][:-1]}.str'
                    href = f'{href}{file_name}'
                    self.files.append(href)

    def run(self,directory):
        URL = 'http://bmrb.io/ftp/pub/bmrb/metabolomics/entry_directories/'
        bmrb_directory = urllib.request.urlopen(URL).read().decode('utf-8')

        parser = BMRB_Downloader.BMRB_Directory_Parser()
        parser.feed(bmrb_directory)

        download_directory  = pathlib.Path(directory,'bmrb_nmr_spectra')
        download_directory.mkdir(exist_ok=True)

        for i,file in enumerate(parser.files):
            file_url = URL + file
            file =  pathlib.Path(file).parts[-1]
            target = pathlib.Path(download_directory, file)
            print (f'download {file} {i+1} of {len(parser.files)}')
            urllib.request.urlretrieve(file_url, target)

class MMCD_Downloader:
    # class BMRB_Directory_Parser(HTMLParser):
    #
    #     def __init__(self):
    #         super().__init__()
    #         self.files = []
    #
    #     def handle_starttag(self, tag, attrs):
    #         if tag == 'a':
    #             href  = [elem for elem in attrs if elem[0] == 'href'][0][1]
    #             if href.endswith('.str'):
    #                 self.files.append(href)

    def run(self,directory):
        URL = 'http://mmcd.nmrfam.wisc.edu/peaklist/'
        FILE_NAME_TEMPLATE = 'expnmr_%05i_3.txt'


        download_directory  = pathlib.Path(directory,'mmcd_nmr_spectra')
        download_directory.mkdir(exist_ok=True)

        for i in range(1,1000):
            file_name = FILE_NAME_TEMPLATE % i
            file_url = URL + file_name

            target = pathlib.Path(download_directory, file_name)
            print (f'download {file_name}')

            try:
                urllib.request.urlretrieve(file_url, target)
            except:
                print("  failed to download...")

if __name__ == '__main__':
    # downloader = HMDB_Downloader()
    # downloader.run('/Users/garythompson')

    downloader = BMRB_Downloader()
    downloader.run('.')

    # downloader = MMCD_Downloader()
    # downloader.run('/Users/garythompson')