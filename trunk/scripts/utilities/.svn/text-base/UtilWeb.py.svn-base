#!/usr/bin/env python

"""

"""

__author__ = "Ruben Acuna"
__copyright_ = "Copyright(c) 2011, ASU iGEM Team"

################################### IMPORTS ####################################
from ftplib import FTP
import urllib2

def fetchFileFTP(url, directory, url_filename, filename):
    server = FTP(url)
    server.login()
    server.cwd(directory)

    file = open(filename, "wb")
    server.retrbinary("RETR " + url_filename, file.write)
    file.close()

    server.close()
    file.close()

    return filename

def fetchFileHTTP(url, filename):
    server = urllib2.urlopen(url)

    file = open(filename,"wb")
    file.write(server.read())

    server.close()
    file.close()

    return filename
