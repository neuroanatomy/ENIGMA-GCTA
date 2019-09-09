#!/usr/bin/env python3

"""Download gcta and plink executables from the web"""

import os
import sys
import stat
import re
from io import BytesIO
from zipfile import ZipFile
from urllib.request import urlopen
from bs4 import BeautifulSoup
import config_dataset

def get_url_paths(url, ext=''):
    """List url paths conained in a HTML page"""
    soup = BeautifulSoup(urlopen(url), 'html.parser')
    return [node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)

    bin_dir = os.path.join(config.annex_dir, "bin")
    bin_subdirs = next(os.walk(bin_dir))[1]

    if sys.platform == "linux":
        plink_re = re.compile(r"plink_linux_x86_64_[0-9]")
        gcta_re = re.compile(r"gcta(?!.*_[a-z]*$)")
    elif sys.platform == "darwin":
        plink_re = re.compile(r"plink_mac_[0-9]")
        gcta_re = re.compile(r"gcta.*_mac")
    else:
        raise Exception("Unknown platform: " + sys.platform)

    plink_dirs = sorted(filter(plink_re.match, bin_subdirs))
    if not plink_dirs:
        print("Installing plink in {}".format(bin_dir))
        plink_dir_url = "http://s3.amazonaws.com/plink1-assets/"
        soup = BeautifulSoup(urlopen(plink_dir_url), 'lxml')
        url_paths = [node.text for node in soup.find_all('key') if plink_re.match(node.text)]
        if not url_paths:
            raise Exception("Unable to find plink url.")
        plink_url = plink_dir_url + url_paths[-1]
        print("Downloading {}...".format(plink_url))
        resp = urlopen(plink_url)
        zipfile = ZipFile(BytesIO(resp.read()))
        print("Done")
        zipfile.extractall(os.path.join(bin_dir, os.path.splitext(os.path.basename(plink_url))[0]))

    gcta_dirs = sorted(filter(gcta_re.match, bin_subdirs))
    if not gcta_dirs:
        print("Installing GCTA in {}".format(bin_dir))
        gcta_dir_url = "https://cnsgenomics.com/software/gcta/bin/"
        url_paths = sorted([url for url in get_url_paths(gcta_dir_url, 'zip')
                            if gcta_re.match(url.rpartition('.')[0])])
        if not url_paths:
            raise Exception("Unable to find GCTA url.")
        gcta_url = gcta_dir_url + url_paths[-1]
        print("Downloading {}...".format(gcta_url))
        resp = urlopen(gcta_url)
        zipfile = ZipFile(BytesIO(resp.read()))
        print("Done")
        zipfile.extractall(bin_dir)


    print("giving executable rights to downloaded executables")
    config = config_dataset.config_dataset(config_file)
    os.chmod(config.gcta, stat.S_IEXEC)
    os.chmod(config.plink, stat.S_IEXEC)


if __name__ == '__main__':
    main(sys.argv[1])
